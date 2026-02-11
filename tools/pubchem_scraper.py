"""
PubChem数据采集工具
用于从PubChem数据库批量采集化合物及其pKa值数据
"""

import os
import sys
import time
import json
import logging
import argparse
from pathlib import Path
from typing import List, Dict, Optional, Tuple, Set
from datetime import datetime

import requests
import pandas as pd
import yaml
from tqdm import tqdm
from rdkit import Chem

# 添加项目根目录到路径
PROJECT_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

# 配置日志
log_dir = PROJECT_ROOT / "data" / "raw"
log_dir.mkdir(parents=True, exist_ok=True)
log_file = log_dir / "scraper.log"

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler(str(log_file)),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)


class PubChemAPI:
    """PubChem REST API封装类"""
    
    def __init__(self, config: dict):
        self.base_url = config['api']['base_url']
        self.view_url = config['api'].get('view_url', "https://pubchem.ncbi.nlm.nih.gov/rest/pug_view")
        self.requests_per_second = config['api']['requests_per_second']
        self.timeout = config['api']['timeout']
        self.max_retries = config['api']['max_retries']
        self.retry_delay = config['api']['retry_delay']
        self.last_request_time = 0
        
    def _rate_limit(self):
        """实现速率限制"""
        min_interval = 1.0 / self.requests_per_second
        elapsed = time.time() - self.last_request_time
        if elapsed < min_interval:
            time.sleep(min_interval - elapsed)
        self.last_request_time = time.time()
    
    def _request(self, url: str, params: Optional[dict] = None, method: str = 'GET', data: Optional[dict] = None) -> Optional[dict]:
        """发送HTTP请求，带重试机制"""
        self._rate_limit()
        
        for attempt in range(self.max_retries):
            try:
                if method.upper() == 'POST':
                    response = requests.post(url, data=data, timeout=self.timeout)
                else:
                    response = requests.get(url, params=params, timeout=self.timeout)
                
                # 特殊处理404，不视为系统错误
                if response.status_code == 404:
                    return {}
                
                response.raise_for_status()
                return response.json() if response.text else {}
            except requests.exceptions.RequestException as e:
                # 400 错误可能是因为请求参数问题，记录但不重试过多
                if isinstance(e, requests.exceptions.HTTPError) and e.response.status_code == 400:
                     logger.warning(f"请求参数错误 (400): {url} - 可能该批次CID无数据")
                     return {}

                logger.warning(f"请求失败 (尝试 {attempt+1}/{self.max_retries}): {e}")
                if attempt < self.max_retries - 1:
                    time.sleep(self.retry_delay * (attempt + 1))
                else:
                    logger.error(f"请求最终失败: {url}")
                    return None
    
    def fetch_pka_annotations(self, heading: str = "Dissociation Constants", page: int = 1) -> List[Dict]:
        """
        直接获取拥有特定Heading数据的化合物列表
        使用 PUG View Annotations API
        
        API: https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/annotations/heading/{heading}/JSON
        
        Args:
            heading: 注解标题 (e.g. "Dissociation Constants")
            page: 页码
            
        Returns:
            包含CID和原始pKa文本的数据列表
            [
                {'cid': 123, 'pka_raw': '3.5 (at 25 C)', 'source': '...'},
                ...
            ]
        """
        # URL 编码 heading
        import urllib.parse
        heading_encoded = urllib.parse.quote(heading)
        
        # PUG View Annotations API
        url = f"{self.view_url}/annotations/heading/{heading_encoded}/JSON"
        params = {
            "heading_type": "Compound",
            "page": page
        }
        
        data = self._request(url, params=params)
        
        results = []
        if not data or 'Annotations' not in data:
            return results
            
        annotations = data['Annotations'].get('Annotation', [])
        
        for anno in annotations:
            try:
                # 获取CID
                cid = anno.get('LinkedRecords', {}).get('CID', [None])[0]
                if not cid:
                    continue
                
                # 获取数据值
                data_list = anno.get('Data', [])
                for item in data_list:
                    # 提取文本值
                    val_obj = item.get('Value', {})
                    string_val = ""
                    
                    if 'StringWithMarkup' in val_obj:
                        string_val = val_obj['StringWithMarkup'][0].get('String', '')
                    elif 'Number' in val_obj:
                        string_val = str(val_obj['Number'][0])
                    
                    if not string_val:
                        continue
                        
                    results.append({
                        'cid': cid,
                        'pka_raw': string_val,
                        'source': anno.get('SourceName', 'Unknown')
                    })
                    
            except Exception as e:
                logger.debug(f"解析注解失败: {e}")
                continue
                
        return results
    

    def extract_pka_values(self, raw_text: str) -> List[float]:
        """从原始文本中提取 pKa 数值（增强版）

        支持多种格式：
        - 单个值: "pKa = 4.5", "3.5 (at 25 C)"
        - 多值列表: "pKa: 2.1, 3.4; 5.0"
        - 范围: "3.5 - 3.7", "3.5 to 3.7", 使用平均值表示
        - 不同破折号/Unicode 连字符（– —）
        - 科学计数法
        会尝试排除显然不合理的数字（非常大的整数等）。最终返回去重后的浮点列表。
        """
        import re
        if not raw_text:
            return []

        text = raw_text
        # 规范化空白与破折号
        if '\u2013' in text or '\u2014' in text:
            text = text.replace('\u2013', '-').replace('\u2014', '-')
        text = text.replace('\u00A0', ' ')
        text = text.replace('\xb0', ' ').strip()

        pka_values: List[float] = []

        # 1) 提取带范围的表达: "3.5 - 3.7", "3.5 to 3.7"
        range_regex = re.compile(r'([-+]?\d+\.?\d*(?:[eE][-+]?\d+)?)\s*(?:-|–|—|to)\s*([-+]?\d+\.?\d*(?:[eE][-+]?\d+)?)', re.IGNORECASE)
        for m in range_regex.findall(text):
            try:
                a = float(m[0])
                b = float(m[1])
                avg = (a + b) / 2.0
                pka_values.append(avg)
            except Exception:
                pass

        # 2) 优先查找明确标注 pKa 的上下文数字
        # 捕获类似 "pKa = 4.5", "pKa: 4.5, 5.6", "pKa 4.5"
        pka_block_regex = re.compile(r'pKa[^0-9A-Za-z\-\+]*(?:[:=~])?\s*([\d\.,;\s\-to–—eE\+\-]+)', re.IGNORECASE)
        pka_blocks = pka_block_regex.findall(text)
        for block in pka_blocks:
            # split by common separators
            parts = re.split(r'[;,/\s]+', block)
            for part in parts:
                part = part.strip()
                if not part:
                    continue
                try:
                    # avoid stray non-numeric tokens
                    if re.match(r'^[-+]?\d+\.?\d*(?:[eE][-+]?\d+)?$', part):
                        pka_values.append(float(part))
                except:
                    pass

        # 3) 查找文中独立数字（当全文非常短或数字上下文像 pKa/dissociation）
        # 如果文本包含温度标注 (如 "25 C")，先移除以避免误识别
        text_no_temp = re.sub(r'\bat\s*\d+\b', '', text, flags=re.IGNORECASE)
        text_no_temp = re.sub(r'\b\d+\s*(?:°|º)?\s*[cC]\b', '', text_no_temp)

        # 如果文本很短（<=30 chars）并且是数字序列，直接接受
        if len(text_no_temp) <= 30 and re.match(r'^[-+]?\d+\.?\d*(?:[eE][-+]?\d+)?(?:[;,\s]+[-+]?\d+\.?\d*)*$', text_no_temp.strip()):
            for num in re.split(r'[;,\s]+', text_no_temp.strip()):
                try:
                    pka_values.append(float(num))
                except:
                    pass

        # 4) 作为兜底，从全文中抓取所有浮点数字，但谨慎过滤：
        if not pka_values:
            # 兜底提取时使用去温度后的文本
            candidates = re.findall(r'[-+]?\d+\.?\d*(?:[eE][-+]?\d+)?', text_no_temp)
            for c in candidates:
                try:
                    val = float(c)
                    # 过滤过大或过小的数值（非 pKa 范围）
                    if abs(val) < 100.0:
                        pka_values.append(val)
                except:
                    pass

        # 去重并按合理值筛选（-10 ~ 50 的宽阈值，后续 DataValidator 会更严格）
        cleaned = []
        for v in set(pka_values):
            if v is None:
                continue
            if v != v:
                continue
            if v < -10 or v > 50:
                continue
            cleaned.append(float(v))

        return sorted(cleaned)

    def get_compound_properties(self, cids: List[int]) -> Dict[int, Dict]:

        """
        批量获取化合物基本属性
        
        Args:
            cids: 化合物CID列表
            
        Returns:
            字典映射 {cid: properties}
        """
        properties = [
            'MolecularFormula',
            'MolecularWeight',
            'CanonicalSMILES',
            'IsomericSMILES',
            'ConnectivitySMILES',
            'IUPACName',
            'Title'
        ]
        
        # 将CID列表转换为字符串
        cid_str = ','.join(map(str, cids))
        url = f"{self.base_url}/compound/cid/{cid_str}/property/{','.join(properties)}/JSON"
        data = self._request(url)
        
        results = {}
        if not data or 'PropertyTable' not in data:
            return results
        
        try:
            props_list = data['PropertyTable']['Properties']
            for props in props_list:
                cid = props.get('CID')
                if not cid:
                    continue
                    
                results[cid] = {
                    'cid': cid,
                    'compound_name': props.get('Title', ''),
                    'molecular_formula': props.get('MolecularFormula', ''),
                    'molecular_weight': float(props.get('MolecularWeight', 0.0)),
                    'smiles': props.get('CanonicalSMILES') or props.get('IsomericSMILES') or props.get('ConnectivitySMILES', ''),
                    'iupac_name': props.get('IUPACName', '')
                }
        except Exception as e:
            logger.error(f"解析化合物属性失败: {e}")
            
        return results


class DataValidator:
    """数据验证器"""
    
    def __init__(self, config: dict):
        self.config = config['validation']
    
    def validate_pka(self, pka: float) -> bool:
        """验证pKa值是否在合理范围内"""
        return self.config['pka_min'] <= pka <= self.config['pka_max']
    
    def validate_smiles(self, smiles: str) -> bool:
        """验证SMILES字符串是否有效"""
        if not self.config['validate_smiles']:
            return True
        
        try:
            mol = Chem.MolFromSmiles(smiles)
            return mol is not None
        except:
            return False
    
    def validate_molecular_weight(self, mw: float) -> bool:
        """验证分子量"""
        return self.config['mw_min'] <= mw <= self.config['mw_max']
    
    def validate_compound(self, compound: dict) -> Tuple[bool, str]:
        """
        验证完整化合物数据
        
        Returns:
            (是否有效, 错误信息)
        """
        # 检查必需字段
        if not compound.get('smiles'):
            return False, "缺少SMILES"
        
        if not compound.get('pka'):
            return False, "缺少pKa值"
        
        # 验证SMILES
        if not self.validate_smiles(compound['smiles']):
            return False, "SMILES无效"
        
        # 验证pKa
        if not self.validate_pka(compound['pka']):
            return False, f"pKa值超出范围: {compound['pka']}"
        
        # 验证分子量
        if compound.get('molecular_weight'):
            if not self.validate_molecular_weight(compound['molecular_weight']):
                return False, f"分子量超出范围: {compound['molecular_weight']}"
        
        return True, ""


class PubChemScraper:
    """PubChem数据采集主类"""
    
    def __init__(self, config_path: str = "tools/config.yaml"):
        # 加载配置
        with open(config_path, 'r', encoding='utf-8') as f:
            self.config = yaml.safe_load(f)
        
        # 初始化组件
        self.api = PubChemAPI(self.config)
        self.validator = DataValidator(self.config)
        
        # 设置输出路径
        self.output_dir = Path(PROJECT_ROOT) / self.config['output']['output_dir']
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        self.output_file = self.output_dir / self.config['output']['output_file']
        self.log_file = self.output_dir / Path(self.config['output']['log_file']).name
        self.checkpoint_file = self.output_dir / Path(self.config['output']['checkpoint_file']).name
        
        # 已采集的化合物（用于去重）
        self.collected_cids: Set[int] = set()
        self.collected_smiles: Set[str] = set()
        self.compounds: List[Dict] = []
        
        # 加载断点
        self._load_checkpoint()
    
    def _load_checkpoint(self):
        """加载断点文件"""
        if self.checkpoint_file.exists():
            try:
                with open(self.checkpoint_file, 'r') as f:
                    checkpoint = json.load(f)
                    self.collected_cids = set(checkpoint.get('collected_cids', []))
                    logger.info(f"从断点恢复，已采集 {len(self.collected_cids)} 个化合物")
            except Exception as e:
                logger.warning(f"加载断点失败: {e}")
        
        # 加载已有数据
        if self.output_file.exists():
            try:
                df = pd.read_csv(self.output_file)
                self.collected_smiles = set(df['smiles'].dropna().unique())
                if 'cid' in df.columns:
                    self.collected_cids.update(df['cid'].dropna().astype(int).tolist())
                logger.info(f"加载已有数据: {len(df)} 条记录")
            except Exception as e:
                logger.warning(f"加载已有数据失败: {e}")
    
    def _save_checkpoint(self):
        """保存断点"""
        try:
            checkpoint = {
                'collected_cids': list(self.collected_cids)
            }
            with open(self.checkpoint_file, 'w') as f:
                json.dump(checkpoint, f)
        except Exception as e:
            logger.warning(f"保存断点失败: {e}")

    def _save_intermediate_results(self):
        """保存中间结果"""
        try:
            if not self.compounds:
                return

            df_new = pd.DataFrame(self.compounds)
            
            # 追加到现有文件或创建新文件
            if self.output_file.exists():
                df_existing = pd.read_csv(self.output_file)
                df = pd.concat([df_existing, df_new], ignore_index=True)
            else:
                df = df_new
            
            # 去重
            if self.config['advanced']['filter_duplicates']:
                df = df.drop_duplicates(subset=['smiles'], keep='first')
            
            df.to_csv(self.output_file, index=False)
            logger.info(f"保存 {len(self.compounds)} 条新数据到 {self.output_file}")
            
            # 清空缓存
            self.compounds = []
            
        except Exception as e:
            logger.error(f"保存中间结果失败: {e}")

    def collect_from_category(self, category: Optional[str] = None, max_compounds: int = 100) -> int:
        """
        基于 Annotation 的 pKa 采集（重构后的实现）

        逻辑：分页遍历 PubChem 注解（heading），解析 pKa 文本 -> 汇总 CID -> 批量获取属性 -> 验证并保存。

        Args:
            category: 占位参数（兼容旧接口），实际使用配置中的 `annotation_heading`。
            max_compounds: 本次采集的最大化合物数量

        Returns:
            成功采集的化合物数量
        """
        heading = self.config['collection'].get('annotation_heading', 'Dissociation Constants')
        page_size = self.config['collection'].get('page_size', 1000)
        batch_size = self.config['collection'].get('batch_size', 100)

        logger.info(f"开始基于注解采集: heading={heading}, 目标={max_compounds}")

        collected_count = 0
        page = 1

        # 按页遍历注解，直到达到目标或无更多数据
        while collected_count < max_compounds:
            annotations = self.api.fetch_pka_annotations(heading=heading, page=page)
            if not annotations:
                logger.info(f"第 {page} 页无注解，停止。")
                break

            # 收集本页所有有 pKa 的 CID，并保留原始文本
            cid_to_text: Dict[int, List[str]] = {}
            for anno in annotations:
                cid = anno.get('cid')
                if not cid or cid in self.collected_cids:
                    continue

                raw = anno.get('pka_raw', '')
                vals = self.api.extract_pka_values(raw)
                if not vals:
                    continue

                cid_to_text.setdefault(cid, []).append(raw)

            if not cid_to_text:
                page += 1
                continue

            cids_page = list(cid_to_text.keys())

            # 批量分块获取属性
            for i in range(0, len(cids_page), batch_size):
                batch_cids = cids_page[i:i+batch_size]
                props_map = self.api.get_compound_properties(batch_cids)

                for cid in batch_cids:
                    if cid not in props_map:
                        continue

                    props = props_map[cid]
                    smiles = props.get('smiles', '')
                    # 合并同一CID所有原始文本，再从每段文本提取 pKa 值
                    raw_texts = cid_to_text.get(cid, [])
                    pka_values = []
                    for rt in raw_texts:
                        pka_values.extend(self.api.extract_pka_values(rt))
                    pka_values = list(set(pka_values))

                    if not pka_values:
                        continue

                    # 全局SMILES去重
                    if self.config['advanced']['filter_duplicates'] and smiles in self.collected_smiles:
                        continue

                    pka_all_str = ";".join(map(str, pka_values))

                    for idx, pka in enumerate(pka_values):
                        compound = {
                            'id': f"pubchem_{cid}_{idx}",
                            'cid': cid,
                            'smiles': smiles,
                            'pka': pka,
                            'pka_all': pka_all_str,
                            'compound_name': props.get('compound_name', ''),
                            'molecular_formula': props.get('molecular_formula', ''),
                            'molecular_weight': props.get('molecular_weight', 0.0),
                            'iupac_name': props.get('iupac_name', ''),
                            'source': f"PubChem_Annotation",
                            'initial_charge': 0,
                            'uhf': 0
                        }

                        is_valid, error_msg = self.validator.validate_compound(compound)
                        if not is_valid:
                            continue

                        self.compounds.append(compound)

                    self.collected_cids.add(cid)
                    if smiles:
                        self.collected_smiles.add(smiles)
                    collected_count += 1

                    # 定期保存
                    if collected_count > 0 and collected_count % self.config['output']['checkpoint_frequency'] == 0:
                        self._save_intermediate_results()
                        self._save_checkpoint()

                    if collected_count >= max_compounds:
                        break

                if collected_count >= max_compounds:
                    break

            # 保存本页结果后继续下一页
            if self.compounds:
                self._save_intermediate_results()
                self._save_checkpoint()

            # 如果本页结果数量少于 page_size，说明已到末页
            if len(annotations) < page_size:
                break

            page += 1

        logger.info(f"注解采集完成，成功采集 {collected_count} 个化合物")
        return collected_count
    
    def run(self):
        """运行完整的数据采集流程"""
        logger.info("=" * 60)
        logger.info("开始PubChem数据采集")
        logger.info("=" * 60)
        
        target_count = self.config['collection']['target_count']
        total_collected = len(self.collected_cids)

        remaining = max(0, target_count - total_collected)
        if remaining <= 0:
            logger.info(f"已达到目标数量 {target_count}，停止采集")
        else:
            count = self.collect_from_category(category=None, max_compounds=remaining)
            total_collected += count
            logger.info(f"当前总计: {total_collected}/{target_count}")
        
        logger.info("=" * 60)
        logger.info(f"数据采集完成！总计采集 {total_collected} 个化合物")
        logger.info(f"数据保存至: {self.output_file}")
        logger.info("=" * 60)


def main():
    """主函数"""
    parser = argparse.ArgumentParser(description="PubChem化合物数据采集工具")
    parser.add_argument("--config", type=str, default="tools/config.yaml", help="配置文件路径")
    parser.add_argument("--target", type=int, help="目标采集数量（覆盖配置文件）")
    # `--category` 参数已移除；采集改为基于 PubChem 注解（Dissociation Constants）
    
    args = parser.parse_args()
    
    # 初始化采集器
    scraper = PubChemScraper(config_path=args.config)
    
    # 覆盖配置
    if args.target:
        scraper.config['collection']['target_count'] = args.target
    
    # 运行采集（基于注解）
    scraper.run()


if __name__ == "__main__":
    main()

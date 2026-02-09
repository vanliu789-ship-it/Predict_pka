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
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("tools/scraper.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)


class PubChemAPI:
    """PubChem REST API封装类"""
    
    def __init__(self, config: dict):
        self.base_url = config['api']['base_url']
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
    
    def _request(self, url: str, params: dict = None) -> Optional[dict]:
        """发送HTTP请求，带重试机制"""
        self._rate_limit()
        
        for attempt in range(self.max_retries):
            try:
                response = requests.get(url, params=params, timeout=self.timeout)
                response.raise_for_status()
                return response.json() if response.text else {}
            except requests.exceptions.RequestException as e:
                logger.warning(f"请求失败 (尝试 {attempt+1}/{self.max_retries}): {e}")
                if attempt < self.max_retries - 1:
                    time.sleep(self.retry_delay * (attempt + 1))
                else:
                    logger.error(f"请求最终失败: {url}")
                    return None
    
    def search_compounds(self, query: str, max_results: int = 100) -> List[int]:
        """
        搜索化合物，返回CID列表
        
        Args:
            query: 搜索关键词
            max_results: 最大结果数
            
        Returns:
            化合物CID列表
        """
        url = f"{self.base_url}/compound/name/{query}/cids/JSON"
        params = {'list_return': 'listkey'} if max_results > 100 else {}
        
        data = self._request(url, params)
        if not data:
            return []
        
        try:
            cids = data.get('IdentifierList', {}).get('CID', [])
            return cids[:max_results]
        except Exception as e:
            logger.error(f"解析搜索结果失败: {e}")
            return []
    
    def get_compound_properties(self, cid: int) -> Optional[Dict]:
        """
        获取化合物基本属性
        
        Args:
            cid: 化合物CID
            
        Returns:
            化合物属性字典
        """
        properties = [
            'MolecularFormula',
            'MolecularWeight',
            'CanonicalSMILES',
            'IsomericSMILES',
            'IUPACName',
            'Title'
        ]
        
        url = f"{self.base_url}/compound/cid/{cid}/property/{','.join(properties)}/JSON"
        data = self._request(url)
        
        if not data or 'PropertyTable' not in data:
            return None
        
        try:
            props = data['PropertyTable']['Properties'][0]
            return {
                'cid': props.get('CID'),
                'compound_name': props.get('Title', ''),
                'molecular_formula': props.get('MolecularFormula', ''),
                'molecular_weight': props.get('MolecularWeight', 0.0),
                'smiles': props.get('CanonicalSMILES', props.get('IsomericSMILES', '')),
                'iupac_name': props.get('IUPACName', '')
            }
        except Exception as e:
            logger.error(f"解析化合物属性失败 (CID: {cid}): {e}")
            return None
    
    def get_pka_data(self, cid: int) -> Optional[List[float]]:
        """
        获取化合物的pKa数据
        
        尝试多个数据源：
        1. PubChem实验数据
        2. 化合物描述中的pKa信息
        
        Args:
            cid: 化合物CID
            
        Returns:
            pKa值列表（可能有多个）
        """
        # 方法1: 从属性中获取
        url = f"{self.base_url}/compound/cid/{cid}/JSON"
        data = self._request(url)
        
        pka_values = []
        
        if data and 'Record' in data:
            try:
                record = data['Record']
                
                # 查找Section中的pKa信息
                sections = record.get('Section', [])
                for section in sections:
                    pka_values.extend(self._extract_pka_from_section(section))
                
            except Exception as e:
                logger.debug(f"从记录中提取pKa失败 (CID: {cid}): {e}")
        
        return pka_values if pka_values else None
    
    def _extract_pka_from_section(self, section: dict) -> List[float]:
        """从section中递归提取pKa值"""
        import re
        pka_values = []
        
        # 检查TOCHeading
        heading = section.get('TOCHeading', '').lower()
        if 'pka' in heading or 'dissociation' in heading:
            # 查找Information字段
            info_list = section.get('Information', [])
            for info in info_list:
                # 查找Value字段
                value = info.get('Value', {})
                if isinstance(value, dict):
                    string_val = value.get('StringWithMarkup', [{}])[0].get('String', '')
                    # 使用正则提取数字
                    matches = re.findall(r'pKa\s*[=:~]?\s*(-?\d+\.?\d*)', string_val, re.IGNORECASE)
                    for match in matches:
                        try:
                            pka_values.append(float(match))
                        except ValueError:
                            pass
        
        # 递归检查子section
        subsections = section.get('Section', [])
        for subsection in subsections:
            pka_values.extend(self._extract_pka_from_section(subsection))
        
        return pka_values


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
                'collected_cids': list(self.collected_cids),
                'timestamp': datetime.now().isoformat(),
                'total_collected': len(self.collected_cids)
            }
            with open(self.checkpoint_file, 'w') as f:
                json.dump(checkpoint, f)
        except Exception as e:
            logger.error(f"保存断点失败: {e}")
    
    def _save_intermediate_results(self):
        """保存中间结果"""
        if not self.compounds:
            return
        
        try:
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
    
    def collect_from_category(self, category: str, max_compounds: int = 100) -> int:
        """
        从特定化合物类别采集数据
        
        Args:
            category: 化合物类别
            max_compounds: 最大采集数量
            
        Returns:
            成功采集的化合物数量
        """
        logger.info(f"开始采集类别: {category}")
        
        # 搜索化合物
        cids = self.api.search_compounds(category, max_results=max_compounds * 3)  # 多搜索一些以备筛选
        logger.info(f"搜索到 {len(cids)} 个化合物CID")
        
        collected_count = 0
        
        for cid in tqdm(cids, desc=f"采集 {category}"):
            # 跳过已采集的
            if cid in self.collected_cids:
                continue
            
            # 获取化合物属性
            props = self.api.get_compound_properties(cid)
            if not props:
                continue
            
            # 获取pKa数据
            pka_values = self.api.get_pka_data(cid)
            if not pka_values:
                continue
            
            # 使用第一个pKa值（通常是最重要的）
            pka = pka_values[0]
            
            # 构建化合物数据
            compound = {
                'id': f"pubchem_{cid}",
                'cid': cid,
                'smiles': props['smiles'],
                'pka': pka,
                'compound_name': props['compound_name'],
                'molecular_formula': props['molecular_formula'],
                'molecular_weight': props['molecular_weight'],
                'iupac_name': props['iupac_name'],
                'source': f"PubChem_{category}",
                'initial_charge': 0,
                'uhf': 0
            }
            
            # 验证数据
            is_valid, error_msg = self.validator.validate_compound(compound)
            if not is_valid:
                logger.debug(f"化合物验证失败 (CID: {cid}): {error_msg}")
                continue
            
            # 去重检查
            if self.config['advanced']['filter_duplicates']:
                if compound['smiles'] in self.collected_smiles:
                    continue
            
            # 添加到结果
            self.compounds.append(compound)
            self.collected_cids.add(cid)
            self.collected_smiles.add(compound['smiles'])
            collected_count += 1
            
            # 定期保存
            if collected_count % self.config['output']['checkpoint_frequency'] == 0:
                self._save_intermediate_results()
                self._save_checkpoint()
            
            # 检查是否达到目标
            if collected_count >= max_compounds:
                break
        
        # 保存剩余数据
        if self.compounds:
            self._save_intermediate_results()
            self._save_checkpoint()
        
        logger.info(f"类别 {category} 采集完成，成功采集 {collected_count} 个化合物")
        return collected_count
    
    def run(self):
        """运行完整的数据采集流程"""
        logger.info("=" * 60)
        logger.info("开始PubChem数据采集")
        logger.info("=" * 60)
        
        target_count = self.config['collection']['target_count']
        categories = self.config['collection']['categories']
        
        total_collected = len(self.collected_cids)
        
        # 计算每个类别应采集的数量
        per_category = (target_count - total_collected) // len(categories) + 1
        
        for category in categories:
            if total_collected >= target_count:
                logger.info(f"已达到目标数量 {target_count}，停止采集")
                break
            
            count = self.collect_from_category(category, max_compounds=per_category)
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
    parser.add_argument("--category", type=str, help="指定单一类别进行采集")
    
    args = parser.parse_args()
    
    # 初始化采集器
    scraper = PubChemScraper(config_path=args.config)
    
    # 覆盖配置
    if args.target:
        scraper.config['collection']['target_count'] = args.target
    
    # 运行采集
    if args.category:
        scraper.collect_from_category(args.category, max_compounds=100)
    else:
        scraper.run()


if __name__ == "__main__":
    main()

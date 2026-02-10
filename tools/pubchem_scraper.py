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
    
    def _request(self, url: str, params: Optional[dict] = None) -> Optional[dict]:
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
        使用PubChem Compound搜索API
        
        Args:
            query: 搜索关键词
            max_results: 最大结果数
            
        Returns:
            化合物CID列表
        """
        # 使用compound/fastformula或name搜索
        # 对于化合物类别，使用不同的策略
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{query}/cids/JSON"
        
        data = self._request(url)
        if not data:
            logger.warning(f"搜索失败: {query}")
            return []
        
        try:
            cids = data.get('IdentifierList', {}).get('CID', [])
            return cids[:max_results]
        except Exception as e:
            logger.error(f"解析搜索结果失败: {e}")
            return []
    
    def get_cids_by_random_sampling(self, start_cid: int, count: int, max_attempts: Optional[int] = None) -> List[int]:
        """
        通过随机采样获取CID列表（更可靠的方法）
        
        Args:
            start_cid: 起始CID
            count: 需要的CID数量
            max_attempts: 最大尝试次数
            
        Returns:
            有效的CID列表
        """
        import random
        
        if max_attempts is None:
            max_attempts = count * 3
        
        cids = []
        attempts = 0
        current_cid = start_cid
        
        while len(cids) < count and attempts < max_attempts:
            # 随机跳跃，避免连续ID
            current_cid += random.randint(1, 100)
            
            # 检查CID是否有效
            url = f"{self.base_url}/compound/cid/{current_cid}/property/MolecularFormula/JSON"
            data = self._request(url)
            
            if data and 'PropertyTable' in data:
                cids.append(current_cid)
            
            attempts += 1
        
        return cids
    
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
        # 方法1: 从属性中获取 (使用 PUG View API)
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/{cid}/JSON"
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
        
        return list(set(pka_values)) if pka_values else None
    
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
                    
                    # 1. 尝试匹配范围 (e.g. "3.5 - 3.7")
                    range_matches = re.findall(r'pKa.*?(\d+\.?\d*)\s*-\s*(\d+\.?\d*)', string_val, re.IGNORECASE)
                    for start, end in range_matches:
                        try:
                            avg = (float(start) + float(end)) / 2
                            pka_values.append(avg)
                        except ValueError:
                            pass
                    
                    # 2. 尝试匹配单个值
                    # 使用正则提取数字
                    matches = re.findall(r'pKa\s*[=:~]?\s*(-?\d+\.?\d*)', string_val, re.IGNORECASE)
                    for match in matches:
                        try:
                            val = float(match)
                            pka_values.append(val)
                        except ValueError:
                            pass
        
        # 递归检查子section
        subsections = section.get('Section', [])
        for subsection in subsections:
            pka_values.extend(self._extract_pka_from_section(subsection))
        
        return list(set(pka_values))


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

    def _get_known_pka_compounds(self) -> List[int]:
        """
        返回已知含有pKa数据的化合物CID列表
        这些是常见的有机酸、碱和药物分子
        """
        known_cids = []
        
        # 1. 羧酸类 (pKa ~2-5)
        carboxylic_acids = [
            176,      # 乙酸 pKa 4.76
            264,      # 甲酸 pKa 3.75
            338,      # 丙酸 pKa 4.87
            1031,     # 丁酸 pKa 4.82
            7991,     # 戊酸 pKa 4.83
            8892,     # 己酸 pKa 4.85
            243,      # 苯甲酸 pKa 4.20
            445858,   # 水杨酸 pKa 2.97
            1060,     # 邻苯二甲酸 pKa 2.89
            10313,    # 草酸 pKa 1.25
            311,      # 琥珀酸 pKa 4.21
            1110,     # 马来酸 pKa 1.92
            444972,   # 富马酸 pKa 3.03
        ]
        known_cids.extend(carboxylic_acids)
        
        # 2. 酚类 (pKa ~8-11)
        phenols = [
            996,      # 苯酚 pKa 9.95
            135,      # 对甲酚 pKa 10.26
            7150,     # 间甲酚 pKa 10.09
            2879,     # 邻甲酚 pKa 10.28
            7342,     # 对硝基苯酚 pKa 7.15
            2394,     # 间硝基苯酚 pKa 8.35
            1970,     # 邻硝基苯酚 pKa 7.23
            6998,     # 儿茶酚 pKa 9.34
            289,      # 氢醌 pKa 9.85
            7054,     # 对氨基苯酚 pKa 10.30
        ]
        known_cids.extend(phenols)
        
        # 3. 胺类 (pKa ~9-11)
        amines = [
            6267,     # 乙胺 pKa 10.7
            6537,     # 二乙胺 pKa 10.9
            6115,     # 苯胺 pKa 4.6
            7515,     # 甲胺 pKa 10.6
            9270,     # 正丙胺 pKa 10.5
            7712,     # 吡啶 pKa 5.2
            1049,     # 吡啶 pKa 5.25
            8082,     # 咪唑 pKa 6.95
            795,      # 咪唑 pKa 7.0
            8248,     # 哌啶 pKa 11.1
            8252,     # 吗啉 pKa 8.5
        ]
        known_cids.extend(amines)
        
        # 4. 氨基酸类 (pKa ~2-10)
        amino_acids = [
            5950,     # L-丙氨酸 pKa 2.34, 9.69
            6137,     # 甘氨酸 pKa 2.34, 9.60
            6287,     # L-缬氨酸 pKa 2.32, 9.62
            6306,     # L-亮氨酸 pKa 2.36, 9.60
            6322,     # L-异亮氨酸 pKa 2.36, 9.68
            6274,     # L-苯丙氨酸 pKa 1.83, 9.13
            6305,     # L-色氨酸 pKa 2.38, 9.39
            6288,     # L-蛋氨酸 pKa 2.28, 9.21
            5960,     # L-脯氨酸 pKa 1.99, 10.60
            750,      # L-谷氨酸 pKa 2.19, 4.25, 9.67
            5961,     # L-天冬氨酸 pKa 1.88, 3.65, 9.60
            6106,     # L-赖氨酸 pKa 2.18, 8.95, 10.53
            6140,     # L-精氨酸 pKa 2.17, 9.04, 12.48
        ]
        known_cids.extend(amino_acids)
        
        # 5. 药物分子
        drugs = [
            2244,     # 阿司匹林 pKa 3.5
            3672,     # 布洛芬 pKa 4.91
            60823,    # 萘普生 pKa 4.15
            2519,     # 咖啡因 pKa 10.4
            3825,     # 尼古丁 pKa 8.0
            4409,     # 吗啡 pKa 8.0
            2585,     # 可待因 pKa 8.2
            3345,     # 扑热息痛 pKa 9.5
            4754,     # 苯巴比妥 pKa 7.4
            2157,     # 水杨酸钠 pKa 2.97
        ]
        known_cids.extend(drugs)
        
        # 6. 杂环化合物
        heterocycles = [
            1174,     # 吡咯 pKa 17.5
            9253,     # 吲哚 pKa 16.97
            1140,     # 呋喃
            8030,     # 噻吩
            1049,     # 吡啶 pKa 5.25
            9246,     # 喹啉 pKa 4.9
            8580,     # 异喹啉 pKa 5.4
        ]
        known_cids.extend(heterocycles)
        
        # 7. 其他重要化合物
        others = [
            962,      # 水 pKa 15.7
            1118,     # 硫酸 pKa -3
            313,      # 磷酸 pKa 2.15
            1032,     # 碳酸 pKa 6.35
            284,      # 氨 pKa 9.25
            1031,     # 硼酸 pKa 9.24
        ]
        known_cids.extend(others)
        
        return list(set(known_cids))

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
        
        # 1. 尝试使用搜索策略
        cids = self.api.search_compounds(category, max_results=max_compounds * 3)
        
        # 2. 如果搜索结果少，尝试使用随机采样策略 (针对已知类别)
        if len(cids) < max_compounds:
            category_ranges = {
                "carboxylic acids": (176, 10000),      # 从乙酸开始
                "phenols": (996, 20000),               # 从苯酚开始
                "amines": (6267, 30000),               # 从乙胺开始
                "amino acids": (5950, 15000),          # 从丙氨酸开始
                "pharmaceuticals": (2244, 50000),      # 从阿司匹林开始
                "benzoic acids": (243, 25000),         # 从苯甲酸开始
                "anilines": (6115, 35000),             # 从苯胺开始
                "pyridines": (1049, 40000),            # 从吡啶开始
                "imidazoles": (795, 45000),            # 从咪唑开始
                "thiols": (6970, 30000),               # 从乙硫醇开始
            }
            if category in category_ranges:
                start_cid, end_cid = category_ranges[category]
                logger.info(f"使用随机采样策略补充数据，范围: {start_cid} - {end_cid}")
                # 增加采样数量，确保能找到足够多的化合物
                random_cids = self.api.get_cids_by_random_sampling(start_cid, max_compounds * 10)
                cids.extend(random_cids)
        
        # 3. 始终添加已知化合物列表作为保底
        logger.info(f"添加已知化合物列表作为补充")
        cids.extend(self._get_known_pka_compounds())

        # 去重
        cids = list(set(cids))
        logger.info(f"准备处理 {len(cids)} 个候选化合物CID")
        
        collected_count = 0
        batch_size = 50
        
        with tqdm(total=len(cids), desc=f"采集 {category}") as pbar:
            for i in range(0, len(cids), batch_size):
                batch_cids = cids[i:i+batch_size]
                
                # 过滤已采集CID
                target_cids = [cid for cid in batch_cids if cid not in self.collected_cids]
                
                # 如果当前批次全部已采集，只更新进度条
                if not target_cids:
                    pbar.update(len(batch_cids))
                    continue
                
                # 批量获取属性
                props_map = self.api.get_compound_properties(target_cids)
                
                for cid in target_cids:
                    pbar.update(1)
                    
                    if cid not in props_map:
                        continue
                    
                    props = props_map[cid]
                    smiles = props['smiles']
                    
                    # 全局SMILES去重
                    if self.config['advanced']['filter_duplicates'] and smiles in self.collected_smiles:
                        continue
                    
                    # 获取pKa数据
                    pka_values = self.api.get_pka_data(cid)
                    if not pka_values:
                        continue
                    
                    # 记录所有pKa值供参考
                    pka_all_str = ";".join(map(str, pka_values))
                    
                    # 遍历每个pKa值，创建独立记录
                    for idx, pka in enumerate(pka_values):
                        # 构建化合物数据
                        compound = {
                            'id': f"pubchem_{cid}_{idx}",
                            'cid': cid,
                            'smiles': smiles,
                            'pka': pka,
                            'pka_all': pka_all_str,
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
                        
                        # 添加到结果
                        self.compounds.append(compound)
                    
                    # 标记为已采集
                    self.collected_cids.add(cid)
                    self.collected_smiles.add(smiles)
                    collected_count += 1
                    
                    # 定期保存
                    if collected_count % self.config['output']['checkpoint_frequency'] == 0:
                        self._save_intermediate_results()
                        self._save_checkpoint()
                    
                    # 检查是否达到目标
                    if collected_count >= max_compounds:
                        break
                
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

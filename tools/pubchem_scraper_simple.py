"""
PubChem数据采集工具 - 简化版
使用已知化合物CID列表直接采集，更可靠、更快速
"""

import os
import sys
import time
import json
import logging
import re
from pathlib import Path
from typing import List, Dict, Optional, Set
from datetime import datetime

import requests
import pandas as pd
from tqdm import tqdm
from rdkit import Chem

# 添加项目根目录到路径
PROJECT_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

# 配置日志
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


class SimplePubChemScraper:
    """简化版PubChem采集器 - 使用已知CID列表"""
    
    def __init__(self, output_dir: str = "data/raw"):
        self.base_url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
        self.output_dir = Path(PROJECT_ROOT) / output_dir
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        self.output_file = self.output_dir / "pubchem_compounds.csv"
        self.requests_per_second = 5
        self.last_request_time = 0
        
        # 已知的含pKa数据的化合物CID（扩展列表）
        self.known_pka_cids = self._build_cid_list()
        
        logger.info(f"初始化完成，已加载 {len(self.known_pka_cids)} 个已知化合物CID")
    
    def _build_cid_list(self) -> List[int]:
        """构建已知含pKa数据的化合物CID列表"""
        cids = []
        
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
        cids.extend(carboxylic_acids)
        
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
        cids.extend(phenols)
        
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
        cids.extend(amines)
        
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
        cids.extend(amino_acids)
        
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
        cids.extend(drugs)
        
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
        cids.extend(heterocycles)
        
        # 7. 其他重要化合物
        others = [
            962,      # 水 pKa 15.7
            1118,     # 硫酸 pKa -3
            313,      # 磷酸 pKa 2.15
            1032,     # 碳酸 pKa 6.35
            284,      # 氨 pKa 9.25
            1031,     # 硼酸 pKa 9.24
        ]
        cids.extend(others)
        
        return list(set(cids))  # 去重
    
    def _rate_limit(self):
        """速率限制"""
        min_interval = 1.0 / self.requests_per_second
        elapsed = time.time() - self.last_request_time
        if elapsed < min_interval:
            time.sleep(min_interval - elapsed)
        self.last_request_time = time.time()
    
    def _request(self, url: str, max_retries: int = 3) -> Optional[dict]:
        """发送HTTP请求"""
        self._rate_limit()
        
        for attempt in range(max_retries):
            try:
                response = requests.get(url, timeout=30)
                response.raise_for_status()
                return response.json() if response.text else {}
            except Exception as e:
                if attempt < max_retries - 1:
                    time.sleep(2 * (attempt + 1))
                else:
                    logger.debug(f"请求失败: {url}, {e}")
                    return None
    
    def get_compound_data(self, cid: int) -> Optional[Dict]:
        """获取单个化合物的完整数据"""
        # 1. 获取基本属性
        props_url = f"{self.base_url}/compound/cid/{cid}/property/MolecularFormula,MolecularWeight,CanonicalSMILES,IUPACName,Title/JSON"
        props_data = self._request(props_url)
        
        if not props_data or 'PropertyTable' not in props_data:
            return None
        
        try:
            props = props_data['PropertyTable']['Properties'][0]
        except (KeyError, IndexError):
            return None
        
        # 2. 获取pKa数据
        pka_values = self._extract_pka(cid)
        
        if not pka_values:
            return None
        
        # 3. 验证SMILES
        smiles = props.get('CanonicalSMILES', '')
        if not smiles or not self._validate_smiles(smiles):
            return None
        
        # 4. 构建结果
        compound = {
            'id': f"pubchem_{cid}",
            'cid': cid,
            'smiles': smiles,
            'pka': pka_values[0],  # 使用第一个pKa值
            'compound_name': props.get('Title', ''),
            'molecular_formula': props.get('MolecularFormula', ''),
            'molecular_weight': props.get('MolecularWeight', 0.0),
            'iupac_name': props.get('IUPACName', ''),
            'source': 'PubChem',
            'initial_charge': 0,
            'uhf': 0
        }
        
        return compound
    
    def _extract_pka(self, cid: int) -> List[float]:
        """从PubChem记录中提取pKa值"""
        # 获取完整记录
        record_url = f"{self.base_url}/compound/cid/{cid}/JSON"
        data = self._request(record_url)
        
        if not data or 'Record' not in data:
            return []
        
        pka_values = []
        
        try:
            record = data['Record']
            sections = record.get('Section', [])
            
            for section in sections:
                pka_values.extend(self._parse_section_for_pka(section))
        except Exception as e:
            logger.debug(f"提取pKa失败 (CID {cid}): {e}")
        
        return pka_values
    
    def _parse_section_for_pka(self, section: dict) -> List[float]:
        """递归解析section中的pKa数据"""
        pka_values = []
        
        # 检查标题
        heading = section.get('TOCHeading', '').lower()
        if 'pka' in heading or 'dissociation' in heading or 'ionization' in heading:
            # 查找数值
            info_list = section.get('Information', [])
            for info in info_list:
                value = info.get('Value', {})
                if isinstance(value, dict):
                    strings = value.get('StringWithMarkup', [])
                    for string_obj in strings:
                        text = string_obj.get('String', '')
                        # 提取pKa数值
                        matches = re.findall(r'pKa\s*[=:~]?\s*(-?\d+\.?\d*)', text, re.IGNORECASE)
                        for match in matches:
                            try:
                                pka = float(match)
                                if -5 < pka < 20:  # 合理范围
                                    pka_values.append(pka)
                            except ValueError:
                                pass
        
        # 递归检查子section
        subsections = section.get('Section', [])
        for subsection in subsections:
            pka_values.extend(self._parse_section_for_pka(subsection))
        
        return pka_values
    
    def _validate_smiles(self, smiles: str) -> bool:
        """验证SMILES有效性"""
        try:
            mol = Chem.MolFromSmiles(smiles)
            return mol is not None and mol.GetNumAtoms() > 0
        except:
            return False
    
    def collect(self, target_count: int = 100) -> pd.DataFrame:
        """采集数据"""
        logger.info("=" * 60)
        logger.info(f"开始采集数据，目标: {target_count} 个化合物")
        logger.info("=" * 60)
        
        results = []
        collected_smiles = set()
        
        # 随机打乱CID列表
        import random
        cids = self.known_pka_cids.copy()
        random.shuffle(cids)
        
        for cid in tqdm(cids, desc="采集进度"):
            if len(results) >= target_count:
                break
            
            compound = self.get_compound_data(cid)
            
            if compound:
                # 去重
                if compound['smiles'] not in collected_smiles:
                    results.append(compound)
                    collected_smiles.add(compound['smiles'])
        
        # 转换为DataFrame
        df = pd.DataFrame(results)
        
        # 保存
        if len(df) > 0:
            df.to_csv(self.output_file, index=False)
            logger.info(f"\n✅ 成功采集 {len(df)} 个化合物")
            logger.info(f"数据已保存至: {self.output_file}")
        else:
            logger.warning("未采集到任何数据")
        
        return df


def main():
    """主函数"""
    import argparse
    
    parser = argparse.ArgumentParser(description="PubChem数据采集工具（简化版）")
    parser.add_argument("--target", type=int, default=100, help="目标采集数量")
    parser.add_argument("--output", type=str, default="data/raw", help="输出目录")
    
    args = parser.parse_args()
    
    scraper = SimplePubChemScraper(output_dir=args.output)
    df = scraper.collect(target_count=args.target)
    
    if len(df) > 0:
        print("\n" + "=" * 60)
        print("📊 数据统计:")
        print(f"  总计: {len(df)} 个化合物")
        print(f"  pKa范围: {df['pka'].min():.2f} - {df['pka'].max():.2f}")
        print(f"  平均分子量: {df['molecular_weight'].mean():.2f}")
        print("=" * 60)


if __name__ == "__main__":
    main()

# 数据采集工具使用指南

## 🔧 问题修复说明

之前的版本遇到404错误是因为PubChem API不支持直接用化合物类别名称搜索。现在提供了两个版本：

### 版本1：简化版（推荐）⭐
**文件**: `pubchem_scraper_simple.py`

**特点**:
- ✅ 使用预定义的已知化合物CID列表
- ✅ 更快速、更可靠
- ✅ 无搜索API依赖
- ✅ 包含100+个常见含pKa化合物

**使用方法**:
```bash
# 快速测试（20个化合物）
python tools/quick_start.py

# 采集100个化合物（推荐）
python tools/pubchem_scraper_simple.py --target 100

# 采集更多
python tools/pubchem_scraper_simple.py --target 500
```

### 版本2：完整版（高级）
**文件**: `pubchem_scraper.py`

**特点**:
- 支持随机采样
- 更多配置选项
- 适合大规模采集

**使用方法**:
```bash
python tools/pubchem_scraper.py --target 1000
```

## 🚀 快速开始（3步完成）

### 步骤1: 更新环境
```bash
cd E:\code\pKa_projected\Predict_pka
conda env update -f environment.yml
conda activate pka_predictor
```

### 步骤2: 测试工具
```bash
python tools/quick_start.py
```

### 步骤3: 正式采集
```bash
# 采集100个化合物（约5分钟）
python tools/pubchem_scraper_simple.py --target 100
```

## 📊 预定义化合物类别

简化版包含以下类别的化合物（共100+个CID）：

| 类别 | 数量 | pKa范围 | 示例 |
|------|------|---------|------|
| 羧酸类 | ~15 | 2-5 | 乙酸(4.76), 苯甲酸(4.2) |
| 酚类 | ~10 | 8-11 | 苯酚(9.95), 对甲酚(10.26) |
| 胺类 | ~12 | 4-11 | 苯胺(4.6), 吡啶(5.2) |
| 氨基酸 | ~13 | 2-12 | 甘氨酸, 丙氨酸 |
| 药物分子 | ~10 | 3-10 | 阿司匹林(3.5), 布洛芬(4.91) |
| 杂环化合物 | ~7 | 5-17 | 咪唑(7.0), 吲哚(16.97) |

## 📈 采集效率

| 目标数量 | 预计时间 | 推荐用途 |
|---------|---------|----------|
| 20 | 1-2分钟 | 快速测试 |
| 100 | 5-8分钟 | 模型原型 |
| 500 | 25-40分钟 | 小规模训练 |
| 1000 | 50-80分钟 | 正式训练 |

**注意**: 实际时间取决于网络速度和PubChem响应速度

## 🎯 采集后的使用

### 1. 查看采集结果
```bash
# 查看CSV文件
Get-Content data/raw/pubchem_compounds.csv | Select-Object -First 10
```

### 2. 数据统计
```python
import pandas as pd
df = pd.read_csv('data/raw/pubchem_compounds.csv')
print(f"总计: {len(df)} 个化合物")
print(f"pKa范围: {df['pka'].min()} - {df['pka'].max()}")
print(df['pka'].describe())
```

### 3. 直接用于训练
```bash
python main.py --data data/raw/pubchem_compounds.csv --n_jobs 4
```

## ⚙️ 高级选项

### 自定义输出目录
```bash
python tools/pubchem_scraper_simple.py --target 100 --output data/my_compounds
```

### 添加更多CID
编辑 `pubchem_scraper_simple.py` 中的 `_build_cid_list()` 方法：

```python
def _build_cid_list(self) -> List[int]:
    cids = []
    
    # 添加你自己的CID
    my_compounds = [
        12345,  # 你的化合物1
        67890,  # 你的化合物2
    ]
    cids.extend(my_compounds)
    
    # ... 其余代码
```

## 🐛 故障排查

### 问题1: 网络错误
```
requests.exceptions.ConnectionError
```
**解决**: 
- 检查网络连接
- 确认可以访问 https://pubchem.ncbi.nlm.nih.gov
- 考虑使用VPN或代理

### 问题2: 采集数据很少
```
只采集到5个化合物
```
**原因**: 某些CID可能没有pKa数据或数据格式不标准

**解决**:
- 这是正常的，PubChem中含pKa数据的化合物有限
- 增加 `--target` 参数的值
- 考虑结合其他数据源（如ChEMBL）

### 问题3: RDKit验证失败
```
Invalid SMILES
```
**解决**: 工具会自动跳过无效SMILES，这是正常的过滤过程

## 🔍 数据质量检查

采集后建议检查数据质量：

```python
import pandas as pd
from rdkit import Chem

df = pd.read_csv('data/raw/pubchem_compounds.csv')

# 1. 检查SMILES有效性
valid_smiles = df['smiles'].apply(lambda x: Chem.MolFromSmiles(x) is not None)
print(f"有效SMILES: {valid_smiles.sum()}/{len(df)}")

# 2. 检查pKa分布
print(df['pka'].describe())

# 3. 检查重复
print(f"唯一SMILES: {df['smiles'].nunique()}/{len(df)}")

# 4. 可视化pKa分布
import matplotlib.pyplot as plt
df['pka'].hist(bins=20)
plt.xlabel('pKa')
plt.ylabel('Frequency')
plt.title('pKa Distribution')
plt.savefig('pka_distribution.png')
```

## 💡 数据增强建议

如果需要更多数据：

### 方法1: 多次运行（不会重复）
```bash
# 第一次
python tools/pubchem_scraper_simple.py --target 50

# 第二次（会追加）
python tools/pubchem_scraper_simple.py --target 50
```

### 方法2: 结合其他数据源
- **ChEMBL**: 大型生物活性数据库
- **DrugBank**: 药物数据库
- **文献数据**: 从论文中提取

### 方法3: 计算预测
对于没有实验pKa的化合物，可以使用其他工具预测：
- Marvin (ChemAxon)
- ACD/Labs
- 基于物理的方法（本项目）

## 📚 参考资源

- [PubChem REST API](https://pubchemdocs.ncbi.nlm.nih.gov/pug-rest)
- [RDKit文档](https://www.rdkit.org/docs/)
- [pKa数据库](https://www.chemicalize.com/)

## 🎯 下一步

采集完数据后：

1. ✅ 检查数据质量
2. ✅ 运行主程序训练模型
3. ✅ 评估模型性能
4. ✅ 根据需要调整采集策略

---

**最后更新**: 2026年2月10日
**状态**: ✅ 已修复404错误

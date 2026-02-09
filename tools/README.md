# PubChem数据采集工具使用指南

## 📖 简介

这是一个专为pKa预测项目设计的化合物数据采集工具，可以从PubChem数据库批量获取含有pKa值的化合物数据。

## 🚀 快速开始

### 1. 环境准备

首先更新conda环境（已添加所需依赖）：

```bash
cd E:\code\pKa_projected\Predict_pka
conda env update -f environment.yml
conda activate pka_predictor
```

### 2. 基础使用

**采集5000个化合物（默认配置）：**
```bash
python tools/pubchem_scraper.py
```

**指定目标数量：**
```bash
python tools/pubchem_scraper.py --target 1000
```

**仅采集特定类别：**
```bash
python tools/pubchem_scraper.py --category "carboxylic acids"
```

**使用自定义配置文件：**
```bash
python tools/pubchem_scraper.py --config my_config.yaml
```

### 3. 输出文件

采集完成后，数据将保存在：
- **主数据文件**: `data/raw/pubchem_compounds.csv`
- **日志文件**: `data/raw/collection_log.txt`
- **断点文件**: `data/raw/.checkpoint.json`（用于断点续传）

## ⚙️ 配置说明

配置文件位于 `tools/config.yaml`，主要参数：

### API设置
```yaml
api:
  requests_per_second: 5  # API请求速率（不要超过5）
  timeout: 30             # 请求超时时间
  max_retries: 3          # 失败重试次数
```

### 采集设置
```yaml
collection:
  target_count: 5000      # 目标化合物数量
  batch_size: 100         # 每批次大小
  categories:             # 采集的化合物类别
    - "carboxylic acids"
    - "phenols"
    - "amines"
    # ... 更多类别
```

### 验证规则
```yaml
validation:
  pka_min: -2.0          # pKa最小值
  pka_max: 20.0          # pKa最大值
  mw_min: 50.0           # 分子量最小值
  mw_max: 1000.0         # 分子量最大值
  validate_smiles: true  # 是否验证SMILES有效性
```

### 输出设置
```yaml
output:
  checkpoint_frequency: 50  # 每采集N个化合物保存一次
```

## 📊 数据格式

输出CSV文件包含以下字段：

| 字段 | 说明 | 示例 |
|------|------|------|
| `id` | 化合物ID | pubchem_12345 |
| `cid` | PubChem CID | 12345 |
| `smiles` | SMILES字符串 | CC(=O)O |
| `pka` | pKa值 | 4.76 |
| `compound_name` | 化合物名称 | Acetic acid |
| `molecular_formula` | 分子式 | C2H4O2 |
| `molecular_weight` | 分子量 | 60.05 |
| `iupac_name` | IUPAC名称 | ethanoic acid |
| `source` | 数据来源 | PubChem_carboxylic acids |
| `initial_charge` | 初始电荷 | 0 |
| `uhf` | 自旋多重度 | 0 |

## 🔄 断点续传

工具支持断点续传功能：

1. **自动保存进度**：每采集50个化合物自动保存
2. **意外中断恢复**：重新运行工具会自动从断点继续
3. **去重机制**：已采集的化合物不会重复采集

如需重新开始，删除以下文件：
```bash
Remove-Item data/raw/pubchem_compounds.csv
Remove-Item data/raw/.checkpoint.json
```

## 📈 性能优化

### 预计采集时间

由于PubChem API限制为5请求/秒：

| 化合物数量 | 预计时间 | 说明 |
|-----------|---------|------|
| 100 | ~5分钟 | 测试用 |
| 1000 | ~40分钟 | 小规模 |
| 5000 | ~3小时 | 推荐规模 |
| 10000 | ~6小时 | 大规模 |

**注意**：实际时间取决于：
- 网络速度
- pKa数据可用性（需要多次请求）
- 服务器响应速度

### 提高效率的技巧

1. **分批采集**：先采集1000个测试，验证效果后再扩展
2. **选择合适类别**：优先采集pKa数据丰富的类别
3. **后台运行**：使用`nohup`或`screen`进行长时间采集
4. **监控进度**：查看日志文件实时了解进度

## 🛠️ 高级用法

### 1. 添加自定义类别

编辑 `tools/config.yaml`：

```yaml
collection:
  categories:
    - "your_custom_category"
    - "another_category"
```

### 2. 调整验证规则

如果需要更宽松的验证：

```yaml
validation:
  pka_min: -5.0    # 扩大范围
  pka_max: 25.0
  validate_smiles: false  # 关闭SMILES验证（不推荐）
```

### 3. 集成到训练流程

采集完成后，直接用于训练：

```bash
# 1. 采集数据
python tools/pubchem_scraper.py --target 5000

# 2. 训练模型
python main.py --data data/raw/pubchem_compounds.csv --n_jobs 4
```

## 🐛 故障排查

### 问题1：请求频繁失败
**原因**：网络问题或API限制  
**解决**：
- 检查网络连接
- 增加 `retry_delay` 参数
- 降低 `requests_per_second`

### 问题2：pKa数据太少
**原因**：某些类别pKa数据稀缺  
**解决**：
- 调整 `categories` 列表
- 降低 `target_count`
- 考虑其他数据源（如ChEMBL）

### 问题3：SMILES验证失败
**原因**：RDKit无法解析某些SMILES  
**解决**：
- 检查RDKit安装
- 临时关闭验证（不推荐）
- 手动清理问题数据

### 问题4：采集中断
**原因**：程序崩溃或手动停止  
**解决**：
- 直接重新运行，会自动从断点继续
- 检查 `.checkpoint.json` 文件是否存在

## 📝 日志分析

查看采集日志：
```bash
Get-Content data/raw/collection_log.txt -Tail 50
```

日志包含信息：
- 采集进度
- 错误详情
- 验证失败原因
- API请求状态

## 🔐 注意事项

1. **遵守使用条款**：PubChem要求合理使用，不要过度请求
2. **数据许可**：确保使用符合PubChem数据使用政策
3. **网络稳定性**：长时间采集建议使用稳定网络
4. **存储空间**：5000个化合物约需5-10MB空间

## 📚 参考资源

- [PubChem REST API文档](https://pubchemdocs.ncbi.nlm.nih.gov/pug-rest)
- [PubChem使用政策](https://pubchemdocs.ncbi.nlm.nih.gov/programmatic-access)
- [RDKit文档](https://www.rdkit.org/docs/)

## 🤝 贡献

如需改进工具，请修改以下文件：
- `tools/pubchem_scraper.py` - 主程序逻辑
- `tools/config.yaml` - 默认配置
- `tools/README.md` - 本文档

## 📄 许可

本工具属于pKa预测项目的一部分，仅供学术研究使用。

---

**最后更新**: 2026年2月9日

# `preprocessor.py` 优化计划

根据您的建议，我将对 `src/preprocessor.py` 进行以下关键改进，以提高数据预处理的质量和效率。

## 1. 核心改进 (Implementation)

### A. 增加“脱盐”与“中和” (Desalting & Neutralization)
-   引入 `rdkit.Chem.SaltRemover` 去除无机盐（如 Na+, Cl-）。
-   引入 `rdkit.Chem.MolStandardize.rdMolStandardize.Uncharger`（或类似的标准化工具）将分子转化为中性形式，确保后续 pKa 预测的基准一致。
-   **代码变更**: 在 `MoleculePreprocessor` 类中新增 `_standardize_mol` 方法，并在 `process_smiles` 中首先调用它。

### B. 构象生成的鲁棒性 (Robust Conformer Generation)
-   优化 `process_smiles` 中的 Embed 逻辑。
-   如果 `AllChem.EmbedMolecule` 返回 -1（失败），则尝试使用 `useRandomCoords=True` 参数再次嵌入。
-   如果仍然失败，记录详细警告。

### C. 引入多进程并行 (Concurrency)
-   重构 `process_batch` 方法，使用 `concurrent.futures.ProcessPoolExecutor` 来并行处理分子列表。
-   这将显著加速大规模数据集的处理速度。
-   *注意*: 虽然 `main.py` 已经在更高级别实现了并行，但在 `preprocessor.py` 内部提供并行能力也是有益的，特别是当它被独立调用时。

## 2. 依赖更新
-   无需更新 `environment.yml`，因为 `rdkit` 已经包含所需的模块。

## 3. 验证 (Verification)
-   修改完成后，我将创建一个简单的测试脚本来验证：
    1.  带盐分子是否被正确清洗。
    2.  并行处理是否正常工作。

请确认是否执行此计划？

# `physics_engine.py` 与 `utils.py` 优化计划

根据您的建议，我将对物理引擎和工具模块进行以下深度改进，以提高特征提取的准确性和可靠性。

## 1. `src/physics_engine.py` 修改
-   **A. 从文件读取电荷 (Read Charges from File)**
    -   在运行完 xtb 后，在 `run_calculation` 方法中添加逻辑，读取生成的 `charges` 文件内容。
    -   将读取到的文件内容传递给解析函数，而不是仅依赖 stdout。
-   **B. 几何优化收敛性检查 (Convergence Check)**
    -   在解析 stdout 之前，先检查是否存在 "GEOMETRY OPTIMIZATION CONVERGED" 标志。如果未收敛，记录警告或视为失败（视策略而定，建议记录警告并标记）。

## 2. `src/utils.py` 修改
-   **A. 增强解析逻辑**
    -   更新 `parse_xtb_output` 函数签名，接受可选的 `charges_content` 参数。
    -   **优先**从 `charges_content` 解析电荷（通常更精确且包含所有原子）。
    -   如果文件不存在，回退到从 stdout 解析。
-   **C. 能量项的精细提取 (Detailed Energy Extraction)**
    -   添加对 `Solvation Free Energy` (广义 Born 溶剂化能) 的提取逻辑。
    -   确保 `Total Energy` 提取无误。

## 3. 验证 (Verification)
-   创建一个新的测试脚本 `test_physics_engine.py` (模拟 xtb 输出或实际运行小分子) 来验证：
    1.  是否正确读取了 `charges` 文件。
    2.  是否提取到了溶剂化能。
    3.  是否能正确识别未收敛状态。

请确认是否执行此计划？

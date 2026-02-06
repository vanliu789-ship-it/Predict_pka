# `utils.py` 优化计划

根据您的建议，我将对 `src/utils.py` 进行以下深度改进，以解决元素配对难题、优先使用 CM5 电荷并增强正则匹配的健壮性。

## 1. 电荷文件与元素的“配对” (Element-Charge Mapping)
-   **问题**: `charges` 文件仅包含坐标和电荷，缺乏元素符号；而 `stdout` 包含元素符号但精度可能略低。
-   **策略**: 结合两者。
    1.  从 `stdout` 中提取有序的元素列表（Element List）。
    2.  从 `charges` 文件中提取有序的电荷数值列表（Charge List）。
    3.  假设两者顺序一致（通常如此），将它们配对合并。
    4.  如果列表长度不一致，则回退到仅使用 `stdout` 的数据，以保证安全性。
-   **代码变更**: 修改 `parse_xtb_output` 和 `parse_charges_file`。

## 2. 优先使用 CM5 电荷 (Prefer CM5)
-   **策略**: 在解析 `stdout` 时，`xtb` 通常会输出 "Mulliken charges" 和 "Mulliken/CM5 charges" (或 "CM5 charges")。
-   **代码变更**: 修改 `extract_charges_from_output` 中的正则匹配优先级，明确寻找 `CM5 charges` 关键字。如果找到 CM5 表格，优先使用；否则回退到 Mulliken。

## 3. 正则表达式的健壮性 (Robust Regex)
-   **策略**: 放宽 `Gsolv` 的匹配规则，使其能同时匹配 `Gsolv` 和 `delta Gsolv`，以及不同的空白字符格式。
-   **代码变更**: 更新 `parse_xtb_output` 中的正则。

## 4. 验证 (Verification)
-   更新 `test_physics_parser.py` 以验证：
    1.  元素-电荷配对逻辑是否正确（模拟长度一致和不一致的情况）。
    2.  CM5 优先逻辑是否生效。
    3.  Gsolv 变体是否能被识别。

请确认是否执行此计划？

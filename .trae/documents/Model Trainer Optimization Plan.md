# 模型训练模块优化计划

根据您的专业建议，我将对 `src/model_trainer.py` 进行以下深度重构，以提升模型的性能、可解释性和可维护性。

## 1. 特征缩放 (Feature Scaling)
-   **问题**: 物理特征（如能量 -100 Eh）与指纹（0/1）量级差异巨大。
-   **策略**: 引入 `StandardScaler` 对物理特征进行标准化。
-   **实现**: 在 `prepare_data` 中，将物理特征部分进行 fit_transform，并将 Scaler 对象保存下来，以便预测时使用。

## 2. 模型持久化完整性 (Full Persistence)
-   **问题**: 仅保存模型对象，丢失了预处理参数（fp_radius, bits, scaler）。
-   **策略**: 创建一个自定义的 `PKaModel` 类（或字典结构），将 XGBoost 模型、Scaler、特征配置打包保存。
-   **实现**: 修改 `train_final_model`，保存一个包含所有元数据的字典。

## 3. 指纹生成效率 (Fingerprint Efficiency)
-   **问题**: `list(fp)` 转换慢。
-   **策略**: 使用 RDKit 的 `ConvertToNumpyArray` 直接转换，避免创建 Python list。
-   **实现**: 优化 `generate_fingerprints` 方法。

## 4. 特征重要性分析 (Feature Importance)
-   **问题**: 缺乏可解释性。
-   **策略**: 训练后提取 XGBoost 的 `feature_importances_`，并结合特征名称（指纹位点 + 物理特征名）输出 Top N 重要特征。
-   **实现**: 在 `train_final_model` 后增加分析步骤，并打印/返回重要性排行。

## 5. 新增物理特征
-   **补充**: 将之前新提取的 `gsolv` (溶剂化能) 加入特征列表。

请确认是否执行此计划？

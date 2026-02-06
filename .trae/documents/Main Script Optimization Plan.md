# `main.py` 优化计划

根据您的建议，我将对 `main.py` 进行以下关键改进，以提升性能、灵活性和可靠性。

## 1. 减少对象实例化开销 (Lazy Loading)
-   **问题**: `process_single_molecule` 每次调用都实例化 `MoleculePreprocessor` 和 `PhysicsEngine`，虽然开销不大但可避免。
-   **策略**: 使用全局变量（或 `multiprocessing` 的 `initializer`）在每个 worker 进程启动时仅初始化一次这些对象。
-   **实现**: 定义 `init_worker` 函数来初始化全局实例，并在 `process_single_molecule` 中复用它们。

## 2. 化学电荷的灵活性 (Flexible Charge)
-   **问题**: `charge=0` 和 `uhf=0` 是硬编码的。
-   **策略**:
    -   检查输入数据（DataFrame）中是否存在 `initial_charge` 和 `uhf` 列。
    -   如果存在，将其传递给 `physics_engine.run_calculation`。
    -   如果不存在，默认为 `0`。
-   **实现**: 修改 `process_single_molecule` 以从 `row` 中提取这些参数。

## 3. 中间状态的持久化 (Checkpointing)
-   **问题**: 长时间运行中途崩溃会导致数据丢失。
-   **策略**:
    -   **批次保存**: 在 `main` 循环中，不等待所有任务完成，而是每隔 N 个（例如 50-100）保存一次 CSV 到 `data/processed/features_partial.csv`。
    -   **断点续传**: 在启动时，检查是否存在已处理的特征文件（或部分文件），并过滤掉已完成的 `id`，只提交未完成的任务。
-   **实现**:
    -   在 `main` 中加载现有的 `processed_file`（如果存在），提取已处理的 ID 集合。
    -   过滤 `data_records`。
    -   在收集结果时，使用 `pool.imap_unordered` 配合计数器，定期追加写入文件（或简单地每批次写入）。为了简化，我们可以每完成一个 chunk 就 append 到 CSV。

请确认是否执行此计划？

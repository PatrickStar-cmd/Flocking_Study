# 时间序列验证分析

有效分析入口：`../../analyzeTimeSeriesValidation.m`。输出：`analysis-v2/`。Python 草稿因缺少 SciPy 未运行，不是本轮结果来源。

本目录包含 180 个 seed-condition MAT 文件（种子 101–130，6 条件），每个文件保存完整 `timeSeries`、配置、summary，以及 100/300 步 block 分析。分析以种子为单位，不能把 block 当独立重复。

本次结果是方法验证：严格阈值来自旧的时间均值分位数，当前 block-pass 只用于描述时间聚合规则的影响；不能据此宣称正式持续成功阈值或 95% 可靠下限。后续须在固定任务标准后用新种子确认。

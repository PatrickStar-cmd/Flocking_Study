# P2：同轨迹开环预测诊断（探索性开发集）

5个复用seed，3类片段S/T/V，2种遮挡长度，30个回放条件，0失败。10条物理输入（含5条新增转向记录），每个条件比较Hold、KF-CV、KF-GM，估计不反馈运动。

完整报告：[P2_REPORT.md](../../../experiments/kf-p2-20261102/P2_REPORT.md)。代码与复现：[README.md](../../../experiments/kf-p2-20261102/README.md)。运行前口径：[PROTOCOL.md](PROTOCOL.md)。

- `raw/input/`：完整2200步物理位置、航向和140条候选遮挡边的真值/观测MAT；其他8条可见边参与物理仿真，不在开环估计目标内。
- `raw/replay/`：30个原始预测NPZ及SHA清单，保留预200/遮挡D/恢复400步；全前缀用于滤波预热。
- `raw/matlab_reference.mat`：调用原版MATLAB生成的对照。
- `analysis/`：270行汇总、6300行逐年龄指标、逐seed配对差、bootstrap、校准及7图。
- `review/`：独立几何/统计/端口复核、AGY审核记录、3张补充图和分层误差表。
- `sources/`、`source_hashes.json`、`input_manifest.json`：冻结源码、上游来源及SHA；补充审核脚本来源另见最终发布审计。
- `protocol.json`、`design.json`、`cases.csv`、`selection_seed*.csv`：参数、真值驱动的选择规则及全部候选评分。
- `turn_pairing.csv`、`config_seed*.json`、`logs/`：转向配对、实际配置和保留的测试/运行日志。

MAT/NPZ通过Git LFS存储，取回方法见代码README。这里只是5种子的开发证据；P3–P5未执行，不能推断闭环GO收益、稳定性或安全保证。历史P1结果保持只读。

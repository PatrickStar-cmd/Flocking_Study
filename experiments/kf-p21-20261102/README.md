# P2.1 几何解耦与KF校准

本次实测已在全局Q选择门槛停止，未使用独立验证seed。正式结论与停止依据见[完整报告](../../effective-social-input/results/kf-p21-v1-20261102-r01/P21_FINAL_REPORT.md)及[判定说明](../../effective-social-input/results/kf-p21-v1-20261102-r01/DECISION_INTERPRETATION.md)，不把补充工作阈值混为原计划硬指标。

计划全文：PROTOCOL.md。运行前细则：RUN_RULES.md。基线Hao_L@81e0080。旧模型、P1/P2脚本和结果只读；不执行P3/P4。

执行顺序：baseline_audit.py → audit_p2_geometry.py → review_g0.py → test_p21.py → synthetic_bearings.py → review_g1.py → generateP21Inputs('development',2026111001:2026111012,'r01') → replay_q_grid.py → audit_independent.py → choose_q.py。若无可信全局q，冻结停止原因，不能强行选择最小训练误差的q。

只有FROZEN_DECISION.json及SHA写定、开发和源代码核验后，才生成冻结seed列表的validation参考轨迹，用replay_q_grid.py --frozen及validate_frozen.py做一次不调参的独立检验。验证n从开发配对差SD、0.025rad目标CI半宽计算，最少30。

Python3.11与P2的NumPy/SciPy/pandas/Matplotlib；MATLAB R2025a。仓库根执行。结果根effective-social-input/results/kf-p21-v1-20261102-r01。正式G0采用g0-r02；初版g0只保留QA历史（时间单位标签错误）。复现必须新run-id，不覆盖已有目录。generateP21Inputs的第四个可选参数runRoot可指定新结果根，先在该根通过G0/G1，再使用stage后缀r01；所有Python命令则显式传入对应新路径。无需修改冻结模型或覆盖旧文件。

## 原始预测的无损存储

raw/input保留2200步物理位置、航向和140条边真值/观测。每source/method的`*_fresh.npz`保存从第1步滤波的共享完整前缀；每D分支同时保留25个窗口的遮挡与恢复400步，用float64位模式相对共享序列的uint64 XOR压缩。engine.recover_branch(path,fresh_dict)可逐位恢复theta/omega/P/nis/age。预200步直接取共享fresh，初始检查点也在分支文件中。没有删去门限后预测或损失精度。

labels与GEOMETRY_FROZEN.json保存纯真值标签/连续速率；仅供离线诊断，未进入预测递推或在线门控。test_engine.py验证暖启动、全部Q的PSD、无污染和无损恢复；matlabQReference.m和review_matlab.py做原MATLAB6配置对照。audit_independent.py独立复算位置标签、完整SHA与非空窗口RMSE。

大型MAT/NPZ/压缩CSV使用局部Git LFS。完整表、原始预测、源码SHA、日志和失败版本保留。发布前逐文件审查，远端写入遵照用户授权；不推main，不强推，公共回退git revert。

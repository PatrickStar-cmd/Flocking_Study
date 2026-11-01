# 异中心神经场 × KF：首轮诊断 A 包

GitHub `Hao_L` 发布入口见 [GITHUB_PUBLICATION.md](GITHUB_PUBLICATION.md)。本次发布包含45个通过 Git LFS 管理的完整主实验 MAT；下文未提交状态记录的是实验完成时的情况。

本目录由 AGY 实现主体代码，Codex 主助手审查、修正并实际执行。基线为 `Hao_L@f0b797220db5d6f2b6ed883e3838614e67cebbe5`，工作分支为 `exp/haoli-kf-diagnostic-20261101`，目标分支为 `Hao_L`。没有修改原模型、历史结果、论文或 main，没有 commit/push/merge。

## 当前交付与证据

交付范围为用户诊断计划的 P0 + P1 A 包。P2 同轨迹开环比较、P3 神经场机制消融及正式独立验证没有执行。最终结果和证据见 [首轮诊断报告](KF_DIAGNOSTIC_A_REPORT.md)。

修订后的快速与全长八项 P0 实际通过；包括原模型精确回归、重复性、D0 四方法一致、遮挡前完整物理/神经/滤波/噪声/RNG 配对、真值隔离、输入质量、完整协方差半正定与分别触发的门限。证据在 `p0-reviewed-fast-v2.json`、`p0-reviewed-p1-v2.json` 和 `review/`。

本批属于**开发集探索性诊断**。五个 seed 是独立统计单位；P0 重复运行、参考轨迹、软件 fixture 和 Fresh 在两个 D 下的重复评价均不增加 n。图中区间为种子层逐点 95% t 区间，不能作同时置信带或确认性保证。

## 冻结配置

N20，M100，dt=.3 模型时间单位，场地 L1000，speedScale=.05，beta1000，连接指数 .5，等权原始强度 ht=.24，角核宽 .4 rad。固定循环索引选择 k15，遮挡其中7条入边，名称为 **index-fixed**，并非随机或真实几何遮挡。

种子 `2026110101–2026110105`；事件起点1401，总长2200，D20/120步（6/36模型时间单位）；预窗口1201:1400，遮挡窗口1401:1400+D，恢复窗口1401+D:1800+D，公共窗口1401:1920。公共 AUC 为 `dt * sum(520个逐步样本)`。预事件状态仍可能自然漂移，不将其称为已建立的稳态。

Fresh 每个 seed 跑一次 D0，保留相同 .01 rad 测量噪声；其余 none/hold/kf-mean/kalman 各跑两个 D，共45条唯一闭环轨迹。R=1e-4，q=1e-4，maxAge=120（模型时间单位），maxVariance=1 rad²；没有调参或受控转向。

## 关键定义

主比较是遮挡期 `meanPair(Hold) - meanPair(KFmean)`，正值为 KFmean 对距更小；GO 改进为 `GO(method)-GO(baseline)`。恢复期为次响应，KFcov 对 KFmean 用于检查协方差展宽。其余对照完整保留，报告探索性区间。

`uRecord(:,:,t)` 与 `positionsRecord(:,:,t)` 为推进第 t 步之前、即第 t−1 步之后的状态。公共检查点在推进 t0 之前保存，`stateAfterStep=t0-1`；绘图显示第 s 步后神经态时读取 `uRecord(:,:,s+1)`，最后一步读取 `finalNeuralState`。

误差均对各自闭环轨迹中的 `trueBearing` 评价，候选 KF 与实际注入分开。Hold 的候选取其实际保持值；Drop 不虚构隐藏边预测。P 与覆盖率只用于 KF，覆盖率为圆周高斯近似；不能代替 P2 的同轨迹开环因果比较。

## 代码与运行

- `generate_instrumented_simulator.py`：生成唯一命名的被动记录副本、数值路径差异和真实 SHA256。
- `runReviewedP0.m`：审查后的八项软件门槛；旧入口 `runP0AuditFast`、`runP0AuditFullSingleSeed` 转调此实现。
- `runReviewedBatch.m`：拒绝覆盖新目录，逐条保存 raw MAT、运行清单、源码和 SHA；失败写 FAILED 并抛错；每条验证共同前缀与实际生效边数。
- `analyze_kf_diagnostics.py`、`plot_kf_diagnostics.py`：流式分析、完整性检查、逐 seed 配对、固定种子聚类 bootstrap 和计划8图。
- `review/independent_review.py`、`review/verify_package.py`：主助手独立几何、神经相干度、统计和哈希复核。

在仓库根目录使用 MATLAB R2025a 和 Python3.11（numpy/scipy/pandas/matplotlib）：

```matlab
addpath('experiments/kf-diagnostic-20261101');
runP0AuditFast;              % 已有报告时拒绝覆盖
runP0AuditFullSingleSeed;    % 已有报告时拒绝覆盖
% 复现实验使用全新的 suffix，不创建新的独立验证证据
runKFDiagnosticBatch(2026110101:2026110105,'my-reproduction-001');
```

重新执行软件检查时，可通过 `runReviewedP0('fast',newReportPath)` 或 `runReviewedP0('p1',newReportPath)` 保留独立版本；批次入口使用目录内已审查的 v2 PASS 门槛。原始或数值代码变化需要重新审查与更新门槛，不沿用旧结果。

```powershell
py -3.11 experiments/kf-diagnostic-20261101/analyze_kf_diagnostics.py <run_dir> --output <new_analysis_dir>
py -3.11 experiments/kf-diagnostic-20261101/plot_kf_diagnostics.py <run_dir> --output <new_analysis_dir>
```

绘图读取同一分析目录的 contrasts.csv，拒绝覆盖图件。分析随机种子20261101，10000次整 seed 重采样，全部 D/方法/指标共用同一重采样索引矩阵；n<2不作区间推断。逐窗口长表为 summary.csv，逐 seed 对照为 contrasts.csv，估计误差为 estimator_diagnostics.csv。

## 输出和保护

正式 A 包目录：`effective-social-input/results/kf-diagnostic-v1-20261101-a01-reviewed/`。包含 config、audit、seed_ledger、manifest、sources、logs、raw 和 analysis。每条 MAT 包含实际 cfg、timeSeries、observationRecords、完整位置/航向/神经态、滤波诊断、公共检查点和传感噪声。

已有362个 tracked 文件的 SHA256 基线在 `review/protected_baseline.json`。原始 MAT 不进入普通 Git 提交，保留于本机；未来共享须单独安排数据传输。独立审核脚本与全部失败日志保留。

首版 AGY 代码中主响应窗口、分叉真值比较、PSD证明、检查点索引、完整性拒绝和绘图等问题已审查修正，清单及首版证据位于 `review/REQUIRED_FIXES.md`、`review/agy-p0-first.json`。早期 stock 参考的 MATLAB 元数据导出曾因不支持 fopen x 模式失败；轨迹在失败前已保存，独立复算后恢复了元数据，记录在 `review/stock_full_reference_recovered.json`；没有重写或隐藏失败。

回退采用停用本目录入口、保留原始结果；若未来提交并集成，使用对应源码提交的 git revert。当前无提交，无需合并或回退 main。

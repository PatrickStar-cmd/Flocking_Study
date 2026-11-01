# 主助手第一次代码审核：正式45条运行前必须修正

正确基线：`Hao_L@f0b7972`；实验分支 `exp/haoli-kf-diagnostic-20261101`。主助手已执行stock契约测试，并对仪器副本的四模式短程和2200步kalman逐字段精确回归，全部通过；`review/instrumentation-regression.json`有证据。不要再跑9条件benchmark浪费计算。

1. **P0-7目前没有证明PSD或对称性。** 两个对角元素非负不足以推出半正定；maxAsym/minEigen被定义但从未使用。新增真实rawP12/rawP21或每步minEig/maxAsym的被动记录，并实际检查P的对称性与2x2特征值/行列式。R匹配应首先断言cfg.bearingObservationVariance==cfg.bearingMeasurementStd^2，经验noise variance只是额外核验。同步generator和生成文件。
2. **P0前缀和D0的检查范围少于打印声称。** pre-check目前只比none和kalman，未比hold/kf-mean，也未比较完整tracker checkpoint状态；D0类似。对全部四方法比较物理、神经、preferred方向、sensorNoise、track x/P/age/bearing/weight，以及轨迹/timeSeries；完整prefix uRecord包含到t0之前的最终fresh更新，不能漏掉step1400。
3. **检查点标签有一个时间步偏差。** uRecord(:,:,t)和positionsRecord(:,:,t)是step t推进前（状态t-1）。目前checkpoint在t=t0-1时保存，实际是状态t0-2。改为在t=t0的推进前保存公共事件检查点，字段明确stateAfterStep=t0-1。图8若画step s之后的u，使用uRecord(:,:,s+1)，最后一步使用finalNeuralState，位置/航向同理。保留这种索引定义并写入README。
4. **P0-1和P0-6声称大于实际检查。** P0-1只检查存在与which，没有branch/commit/dirty/SHA。补真实版本哈希元数据（或引用主助手已核验的protected baseline和动态git，不能伪造）。P0-6目前只断言inputMass正数，需比较等权15→8的离散质量比例8/15，以及hold和KFcov在相同center/权重下卷积保质量的契约。原库testOcclusionPrediction覆盖后者，可调用并记录。
5. **失败不能显示COMPLETED。** batch与singleSeed捕获失败后仍把manifest.status写COMPLETED、single还打印All9 completed。失败须记录堆栈并状态COMPLETED_WITH_FAILURES或FAILED，最后显式报错；拒绝让失败输出进入统计。起跑前写config、branch/commit/dirty、source快照+SHA、seed ledger、P0验证引用，不能只在全部结束后才写audit。源代码和原始MAT按raw/、logs/分层或清楚记录实际布局；每条参数必须可追溯。输出suffix/seed检查防止路径逃出指定results新目录。seeds输入要求唯一合法整数。
6. **诊断定义要按实际注入区分。** 记录P截断原因须仅在使用P门限的KF模式适用；Hold忽略P，不能因为age与P同时超阈就把Hold计为P/联合截断。none的隐藏边不注入，不能把其后台KF轨迹称为Drop的有效预测。Kalman误差可同时记录candidate（即使截断）与actual injected，两者分开，不能静默丢掉坏预测。闭环各自轨迹误差不能称作同轨迹开环因果结论。

主助手独立复核：`review/independent_gate_verification.json`显示P11递推与rawVariance最大误差5.55e-17，第101个缺测步起关闭KF，D120共2800/16800=16.7%隐藏边时间被P截断，age截断0。首seed预窗口GO前100步均值0.2920、后100步0.3653，仍自然演化；保持预设事件时间作探索性诊断，明确不称稳态失效边界，不按GO筛种子。

已运行第一版runP0AuditFast并得到PASS日志，但上面的覆盖缺口意味着其ALL8 PASS尚不足以放行。原报告已另存review/agy-p0-first.json，日志review/agy-p0-first.*.log；修正后请输出带唯一版本名的新P0报告，保留旧证据，不覆盖。你不调用终端，修正后由主助手执行。

## 分析与绘图首版新增的必须修正项

7. **主响应被改错了窗口。** analyze_kf_diagnostics.py用pair_hold_rec-pair_mean_rec定义BenefitPair，但计划主响应明确是遮挡期occ；恢复期是次响应。请分别输出occ与rec的BenefitPair，不可将恢复期冒充主分析。同时输出各方法相对同seed Fresh的GO损失和对距增量，完整保留KFcov vs KFmean比较及GO/LO；摘要逐窗口长表要包含用户要求的commit/scenario/mask_rule/N/k/hidden_count/D_steps/dt/t0/window/GO_mean/LO_mean/pair_norm_mean/nearest_norm_mean/active_neighbors_mean/input_mass_mean/bump_concentration_mean/speed_mean/result_path。
8. **Hold误差目前用了另一分支真值。** 代码true_b来自Kalman rec，但hold_b来自Hold已分叉轨迹，拿它们相减是错误的。必须Hold估计对Hold自身rec.trueBearing，Kalman估计对Kalman自身truth，显式标签own-closed-loop diagnostics；不能混成P2同轨迹比较。不注入和被门限截断的candidate与实际注入分开。诊断CSV按seed×D×mode×age记录n_hidden、MAE、RMSE、P、近似coverage95、fresh NIS、age/P截断；统计量依seed汇总，不能把140条边当140个独立样本。
9. **统计随机性未冻结。** np.random.choice无种子，重复分析bootstrap会变。固定并记录analysis seed，按seed整组重采样同一矩阵；t临界值使用scipy.stats.t.ppf(0.975,n-1)，n<2不做推断，不硬编码非5时2.0。验证manifest全部45个唯一键且完成、5个seed无重复、所有窗口长度正确、Fresh仅5条；缺数据/重复不能静默覆盖字典或continue后当成功。
10. **缺少A包最关键的神经峰集中度图。** 首版8图用估计误差/P代替了计划的bump与其他核心图。请生成计划8张：P1_Fresh_Drop_Hold_KFmean_KFcov_GO.png、P1_PairDistance.png（相对各自同seed Fresh差）、P1_LocalOrder.png、P1_ActiveNeighbors.png（fresh和active）、P1_InputMass.png、P1_BumpConcentration.png、P1_GainScatter.png、P1_FailureCases.png（按最小seed固定选择，有位置/航向和u或positive activity，正确step索引）。估计器诊断图若保留可作为额外图，不替代这8张。种子区间用ddof=1与实际有效seed数，Fresh重复评价不增加n。
11. **全量加载内存风险。** 每条包含35MB uRecord+数十MB observations，当前分析器把45条全部常驻runs，画图再全量加载；可能数GB。请至少分析按seed流式、绘图仅存所需timeSeries和代表seed神经态；可以在每条处理后释放大字段，防止Python对象嵌套导致资源浪费。

# 遮挡与卡尔曼：实现记录

最新扩展结果（2026-09-15）：新增 `runPaperKFExtendedSimulation` 和
`plotPaperKFExtendedResults`，使用全新种子 201--230、固定身份 k=15、
遮挡隐藏数 {7,11,13}、T={3,6,10.2,18,30,60}，比较 fresh-drop、hold-last、
KF-mean、常角速度 KF-unc 与阻尼角速度 KF-GM-unc。共 2700 条配对记录，
无重复且每组 30 个种子。相对 fresh-drop，各类 KF 在全部扫描条件下均降低
遮挡期绝对 GO 偏差和遮挡后对距偏差，配对 bootstrap 区间不跨零；主要收益
随 T 增长而扩大。阻尼 KF 在 T<=30 将隐藏边平均方位误差相对 hold-last
降低约 11%--19%。结果见 `results/paper-kf-extended-v3/README.md`。
该实验仍是使用新种子的探索性检查点分叉，不是阈值化成功率的确认性验证，
也未验证碰撞安全。常角速度与阻尼角速度结果必须分开报告。

最新状态：一阶 Gauss-Markov 角速度候选的离线开发已经完成。它只改善了 T=10.2/30 的平均预测，未修复 T=3 退化，协方差一致性也未通过；见 `results/bearing-damped-model-v1/README.md`。因此没有替换现有线上常角速度实现，没有启动预测闭环 T 扫描。下一实现项改为完整状态检查点分叉，先比较无遮挡、删除和 hold。

完整状态检查点现已实现：仿真器保存并恢复 `finalPreferredDirections`，避免用当前航向错误重建 allocentric 神经元坐标。`testCheckpointFork` 验证 T=0、共同前缀和确定性重放。N20 两种子 pilot 已运行；预测分支仍未接入，结果边界见 `results/checkpoint-intervention-pilot-v1/README.md`。

最新修正：已依据初稿清单补齐匹配R的显式高斯方位测量噪声、跨分支共享噪声和滤波age/innovation/NIS诊断，见PREDICTION_CHECKLIST_REVIEW.md。下文关于理想无噪声观测的描述属于早期集成实验版本，不适用于新默认配置；历史结果保留，不重写。预测准确性和协方差校准仍未完成。

更新：已完成历史粗采样轨迹的预测筛查，结果见results/bearing-screen-v1/README.md。当前预测器在dt=30下无稳定优势，名义95%区间覆盖约50%–55%。完整检查点干预方案在OCCLUSION_INTERVENTION_PLAN.md；正式T扫描仍未启动。

## 已实现并测试

默认关闭的 scheduledOcclusion 分支已接入 simulateEffectiveSocialInput。旧等权路径逐项回归通过。第一版仅允许 allocentric、original 单边增益、ideal 观测；不启用已有FOV/几何遮挡，避免混淆。

固定身份集合：个体i选择循环编号之后的k个目标；遮挡固定隐藏前h条边，区间为 startStep <= step < startStep+durationSteps。这是有向固定拓扑测试，不等同于之前每步random选邻，也不是方位平衡或真实障碍物。k为固定目标数，新鲜输入数与激活输入数分开保存。

传感器 scheduledOcclusionInput 仅对当前可见目标计算观测距离和方位；隐藏边方位与权重为NaN。occlusionInputStep 不接收位置或真实速度。距离权重在隐藏期间保持最后观测；每条边age、weight、variance、fresh、active保存在observationRecords。

滤波状态为[方位,方位角速度]，F=[1 dt;0 1]，白噪声角加速度谱密度q，Q=q*[dt^3/3 dt^2/2;dt^2/2 dt]。观测矩阵[1 0]，创新atan2(sin(delta),cos(delta))，Joseph形式协方差更新。初始化角速度0、速度方差1。默认q=1e-4，R=1e-4，均为待校准设计参数；当前观测仍是理想真值观测，R是滤波器假设，并没有注入符合R的测量噪声。不能声称协方差已经校准。圆周单峰滤波不适用于任意大方差或多峰观测。

三对照在可见时都直接注入相同观测，后台滤波并行更新；隐藏后none删除、hold保持方位和权重、kalman传播状态并使用圆周热核卷积。FFT卷积对频率m乘exp(-Ptheta*m^2/2)，保留离散输入质量并降低峰；这是圆周离散近似，不是预测目标存在概率下降。默认maxAge=120模型时间、maxVariance=1 rad²，超限不注入。hold也受age门限约束。当前未实现存在概率和数据关联错误。

## 验证与结果

testOcclusionPrediction覆盖：旧等权回归、跨2pi创新、协方差PSD和无观测传播、隐藏目标位置改变不影响对应控制输入、遮挡起止边界、可见/激活数量、零时长对照一致、卷积质量与峰值、协方差截止。

runOcclusionIntegrationPilot完成六次测试：N8、Ns40、k4、隐藏2目标、种子136–137、600步。step301–400隐藏100步，即30模型时间；观测时刻从(step-1)*dt计，90开始、120结束。指数距离权重Rw=.5L。结果位于results/occlusion-integration-v1。

两种子新鲜数均2；none激活数2，hold/kalman为4。遮挡期平均GO分别约.33154和.72069，各模式相同到显示精度；输入质量none约.02973/.03680，hold/kalman约.05990/.06547。这是接线测试，不证明预测有益，也不是N20长时阈值实验。

## 尚未完成

正式阈值化T扫描、测量噪声校准、真实几何遮挡、近远/方位遮挡选择、数据关联、独立预测误差与协方差一致性验证仍未完成。`results/paper-kf-extended-v3` 是扩展后的探索性配对扫描。angularCoverage在新分支显式NaN；inputConcentration是实际场的一阶圆周矩，与旧未加权方位集中度定义不同，不能混合分析。不输出虚假的真值覆盖率。

下一项先做离线已知轨迹（匀角速/机动/跨界）的误差和协方差校准，再冻结N20的遮挡起点、比例、T网格及q/R/截止规则。完整扫描应包含零时长与无记忆、hold、kalman对照；不能在本轮功能测试后直接宣称预测降低了有效输入下限。

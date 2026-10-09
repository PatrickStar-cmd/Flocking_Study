# r3 数据字典

本包是协议 1.2 的 seed1 软件验收，不是科学收益验证。下列维度采用主配置：M=100、N=20、S=600；S 对应绝对步 6001…6600。MAT 每个模式文件只有一个顶层变量 `run`。

## 时间与采样

- `run.segment.state.position[2,N,S+1]`、`heading[N,S+1]`、`u[M,N,S+1]`：601 帧，绝对步末 6000…6600。第 j 步的步初状态为第 j 帧，步末为第 j+1 帧。
- `state.alpha[N,M]`：实际神经网格，allocentric 模式固定；不是统一的零起点网格。
- `data.step[S,1]`：绝对步；`localStep=step-6000`；`timeStart=(step-1)*dt`、`timeEnd=step*dt`，时间单位为模型时间。
- 输入、观测、真值及缓存诊断每步在运动更新前采样；GO/LO/距离、神经诊断在神经与运动更新后采样。
- `phase` 由事件起点及 duration 推导。d34 缺失 6201…6234，恢复 6235；d0 没有缺失 phase，也不生成 mask。

## 群体字段 [S,1]

`GO` 是实际速度单位向量的群体合矢量长度除以 N；静止速度贡献零。`LO` 是周期最短像距离下最近 5 邻居及自身的局部有序度均值。`pairDistance` 为所有无序对平均周期距离除以 L；`nearestDistance` 为个体最近距离均值除以 L。`driverIndex` 是本步控制缓存推进的绝对索引。无观测噪声，无逐步随机驱动；两级 checkpoint 及最终状态保存完整全局 RNG。

## 个体输入字段 [S,N]

| 字段 | 含义 |
|---|---|
| k_selected / k_fresh / k_active | 固定目标数 / 本步合法新鲜观测数 / 正权重实际激活边数 |
| A | 注入边总权重，sum(edgeWeight)，不是输入场积分 |
| Mass / Peak | (2π/M)sum(inputField) / 场最大值 |
| z1_real / z1_imag / z1_mod | 归一化加权方位一阶矩实部、虚部、模 |
| z2_real / z2_imag / z2_mod | 二阶矩，用于识别对向双峰等结构 |
| direction / direction_defined | 主方向（弧度，[-π,π]）及 A>0 且 R1>1e-8 的有效标记 |
| directionValid | direction_defined 的兼容别名 |
| weightEntropy | -sum(p log p)，自然对数，标签权重熵，不是方位熵 |
| weightEntropyNormalized | 旧 r2 数值，entropy/max(log(k_active),1)，明确另命名 |
| maxGap | 激活正权重方位的最大圆周间隙，弧度 |
| no_input | A=0；零场、A/Mass/Peak/k_active 和矩数值均为0 |
| normalization_unavailable | drop-normalized 在至少一条边缺失且无新鲜观测时为真 |

空输入时 `direction`、`maxGap`、`weightEntropy` 为 NaN，方向标记为 false；矩的0只是带 no_input 的占位值，不表示实测均匀分布。

## 逐边字段 [N,N,S]

第一个维度为 observer_id，第二个为 target_id。`selected/fresh/hidden/active` 是布尔标记。`observedBearing/observedWeight` 仅包含合法观测，隐藏及未选择边为 NaN。`lastBearing/lastWeight/age` 是本步缓存更新后的状态；最后合法方位来自步初观测。隐藏 hold 方位在缺失期不更新，年龄从 dt 递增；恢复首步年龄归0。未曾缓存的年龄为 Inf。

`injectedBearing/edgeWeight` 是实际注入的中心（弧度，[0,2π)）和权重；未激活中心为 NaN、权重为0。`trueBearing/trueDistance` 是独立离线评估通道，距离单位与 L 一致；它们不属于控制器参数。`controllerInputFromObservation` 公开入口只接收 alpha、合法缓存、obs、cfg、t，不能读取位置或真值。

`hiddenR1/hiddenR2/hiddenMaxGap`、`remaining*`、`all*` 均为 [S,N]，分别使用冻结隐藏集合、其余集合、全体目标在各分支本步真实几何中的结构；不用于反馈。空集合 R1/R2为0占位、gap为NaN。选择时刻的完整集合摘要、候选实际ID元组及随机流来源另保存在 checkpoint.maskDiagnostics 的20个cell中。

## 神经字段 [S,N]

`neuralActivityMass=sum(max(tanh(beta*u),0))` 为非负活动的离散总量，不是角积分；`bumpConcentration` 为神经活动归一化一阶合矢量模。`rawDecodedAngle` 是原始 atan2 解码（[0,2π)）；`decoded_direction_defined` 要求活动质量>0且 bumpConcentration>1e-8。无效时 `decodedAnalysisAngle=NaN`。分析标记不修改运动规则；实际航向保存在 state.heading，原运动阈值仍为 abs(cx)+abs(cy)>=1e-9。

## 两类 oracle

1. 每条闭环分支 `data.selfFreshField[M,N,S]`：使用该分支自身步初真实几何及实际 alpha 的 fresh 场。`selfFreshA/Mass/Peak`、`selfFreshError`（场最大绝对差）和 `selfFreshRelativeL2`（逐个体相对 L2）均为 [S,N]。L2 分母为0时为NaN。
2. `reference_same_state_replay.mat` 的 `replay`：固定 fresh 参考轨迹，不推进候选群体状态。保存四候选 field、A/Mass/Peak、矩/gap及 relativeL2；缓存由 fresh 轨迹共同前缀合法观测推进，隐藏期只使用 step6200-start 最后缓存。候选从各步 fresh 物理状态构造，不与闭环分支的自身 oracle 混用。

## 检查点与摘要

`checkpoint6000/6200` 是标量结构体，包含位置、航向、u、alpha、固定目标ID、冻结mask、20个诊断cell、缓存/年龄、完整rng、driverIndex、config。`checkpoints.mat` 提供独立 cp6000/cp6200。公开恢复入口从 cp.step+1 开始，不初始化或重置缓存。

`summary.csv` 共16行（4模式×4窗口）：缺失34步、恢复100步、恢复300步、固定事件后400步。GO/LO loss=fresh−branch；对距增量=branch−fresh；另存 normalized−drop、hold−normalized 同指标差。输入误差使用分支自身 oracle。所有摘要来自逐步 MAT，fresh 是共享参考，不作为重复独立样本。

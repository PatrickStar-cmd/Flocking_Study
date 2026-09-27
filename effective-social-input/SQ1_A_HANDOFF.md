# 科学问题一 A 交接说明：社会输入建模与局部机制 E1

## 1. 交接状态

- 负责人：Hao_L（A）
- 任务分支：`feature/sq1-e1-handoff`
- 代码基点：`YCwork@71583f9`
- 交接类型：E1 局部机制阶段性交付
- 证据等级：确定性、探索性开发证据
- 结果目录：`effective-social-input/results/sq1-local-v1/`

本交付可以作为 B 的 E2 无遮挡群体实验和 C 的 E3 缺失设计/统计分析的输入层基础。它不代表科学问题一的群体实验、独立验证或确认性结论已经完成。

## 2. 已完成工作

### 2.1 输入模型

单个 focal agent 的方位社会输入定义为

\[
I(\alpha_m,t)=\sum_j w_j(t)\exp\!\left[-\frac{d_c(\alpha_m,\theta_j(t))^2}{2\sigma^2}\right],
\]

其中 `theta_j` 是邻居方位，`w_j >= 0` 是输入权重，`d_c` 是圆周距离，`sigma` 是方位核宽度。

已明确区分：

- `k`：当前输入标签/邻居数；
- `A=sum(w)`：输入总权重；
- `p(theta)`：归一化方位结构；
- `k_fresh`：当前新鲜观测数；
- `k_active`：实际注入神经场的输入数。

当 `A=0` 时定义为零输入，不计算归一化方位形状。

### 2.2 输入指标

代码计算并保存：

- 离散输入积分 `inputMass`；
- 输入峰值 `fieldPeak`；
- 一阶圆周矩 `z1`；
- 二阶圆周矩 `z2`；
- 方位集中度、二阶集中度；
- 角度熵；
- 最大方位空隙和角度覆盖；
- 神经场解码方向；
- `bumpConcentration`。

### 2.3 E1 对照

在固定神经初态和预设方位序列下完成：

1. 强度缩放：总权重为 `0.5/1/2`；
2. 方位结构：集中、双峰、八方位分散；
3. 标签拆分等价性：同一方位拆成多个权重；
4. fresh bearing 与 hold-last bearing 的时间变化对照。

## 3. 文件清单

### 3.1 模型和代码

- `SQ1_LOCAL_INPUT_MODEL.md`：公式、假设、指标、局部神经场更新和证据边界。
- `localSocialInputMetrics.m`：单时刻输入场与指标计算。
- `replayLocalSocialInput.m`：固定神经初态的单体神经场回放。
- `runSq1LocalE1.m`：生成 E1 的 9 个确定性案例。
- `testLocalSocialInput.m`：E1 软件契约测试。
- `plotSq1LocalE1.m`：重建机制图。

### 3.2 已生成结果

`results/sq1-local-v1/` 包含：

- `summary.csv`：9 个案例的汇总指标；
- `e1-results.mat`：完整回放结果、逐步输入、神经状态和配置；
- `e1-mechanism.png`：输入积分、方位集中度、场峰值和 fresh/hold 解码方向图；
- `README.md`：结果说明；
- `audit.json`：参数、代码哈希和输出清单。

历史结果目录没有被覆盖。

## 4. 参数和复现命令

E1 参数：

- `NNeurons=100`；
- `NSteps=120`；
- `dt=0.3`；
- `beta=1000`；
- `globalInhibition=0`；
- `receptiveFieldWidth=0.4`；
- `connectivityExponent=0.5`；
- 初始神经状态为零；
- 预设方位序列；
- 无随机种子、无 KF、无距离、无遮挡。

MATLAB 运行：

```matlab
cd('effective-social-input')
testLocalSocialInput
runSq1LocalE1
plotSq1LocalE1
```

## 5. 已通过的验证

- 零输入边界测试通过；
- 离散输入积分复算通过；
- 同一预设输入重复回放结果完全一致；
- 标签拆分后输入场和神经状态最大差异小于 `1e-12`；
- 强度 `1.0/0.5` 的输入积分比为 `2.0`；
- 八方位分散案例的 `z1` 和 `z2` 均接近零；
- 核心群体仿真器最小 smoke test 通过；
- `git diff --check` 通过。

## 6. B 的使用方式

B 可以复用 `localSocialInputMetrics` 中的指标命名和圆周计算约定，将这些字段接入群体仿真器的逐步 `timeSeries`：

```text
k_selected, k_fresh, k_active, A, inputMass, fieldPeak,
z1_real, z1_imag, z2_real, z2_imag, maxGap,
GO, LO, pairDistance
```

B 负责 E2 的群体动力学、选邻算法、固定身份/动态重选方案、共同随机流、检查点恢复和批量运行。`fixed-total` 只用于机制拆分，不能直接作为生物控制律解释。

建议 B 在 E2 中继续保留：

- `original` 与 `fixed-total` 两种输入制度；
- 实际方位结构，而不是只记录策略标签；
- 输入字段与 GO、LO、凝聚距离的逐步时间序列；
- 配置快照、种子、代码提交号和运行日志。

## 7. C 的使用方式

C 可以使用 E1 的输入字段作为 E3 输入误差分析的统一口径，但需单独定义缺失干预：

- `fresh`：新鲜输入参考；
- `drop`：删除缺失邻居输入；
- `drop-normalized`：保持预先声明的缺失前总幅值的机制对照；
- `hold-last`：保持最后一次方位和权重。

遮挡期间不得读取隐藏目标真实方位、位置、速度或距离。离线 oracle 必须按受干预分支自身物理状态重建，只用于误差评价，不能直接作为控制器输入。

C 负责冻结缺失规则、T 网格、主响应、种子台账、配对统计和置信区间。E1 的 9 个确定性案例不是独立随机样本，不能用于群体层面样本量或确认性验证。

## 8. 尚未完成的工作

以下内容不属于本次 E1 交接的完成项：

1. E0 全套基线和完整状态恢复验收；
2. E2 无遮挡群体因素实验；
3. E3 连续缺失机制干预；
4. 独立种子验证和样本量冻结；
5. 群体 GO、LO、凝聚和网络连通性的统计结论；
6. 距离权重、预测补偿和 KF 相关实验。

## 9. 结论边界

E1 支持的结论是：在固定单体神经初态和给定输入序列下，输入总强度与方位结构分别能够改变局部神经场的幅值、峰形和解码响应；同一完整输入被拆成多个标签不会产生额外局部效应。

E1 不支持以下结论：某个邻居数是普适下限、某种方位结构一定改善群体运动、预测方法提高闭环稳定性、系统具有碰撞安全保证，或局部响应可以直接替代群体实验。

## 10. 交接后的建议顺序

1. B/C 审阅本文件和 `SQ1_LOCAL_INPUT_MODEL.md`，共同冻结字段和定义；
2. B 完成 E0/G0 基线与接口回归；
3. C 完成种子台账、缺失规则和统计终点冻结；
4. B 运行 E2，A 复核实际输入字段；
5. B/C 运行 E3，A 复核 `drop`、`drop-normalized`、`hold-last` 的输入属性；
6. 仅在规则冻结后，使用未参与方案选择的种子做独立验证。

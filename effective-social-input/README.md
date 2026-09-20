# Effective Social Input for Allocentric Flocking

本目录是在 Salahshour 与 Couzin (2025) 的公开 MATLAB 代码基础上建立的独立实验版本。原始代码位于相邻目录 `allocentric-flocking-original`，本目录不会修改原始文件。

## 当前研究范围

- 二维周期边界；
- 异中心环形神经场；
- 每个邻居只提供方位角高斯输入；
- 暂不使用距离调节、避撞、障碍遮挡或预测；
- 控制邻居数量和选邻策略；
- 对比固定单邻居输入和固定总输入两种制度；
- 输出原文使用的全局序参量 GO、拓扑局部序参量 LO、平均两两距离，以及输入强度、角度覆盖和神经场状态指标；在可见性扩展中同时记录新鲜可见邻居数与实际激活邻居数。

## 文件

- `defaultEffectiveInputConfig.m`：与作者演示代码一致的默认参数。
- `buildRingConnectivity.m`：构造环形神经场连接矩阵。
- `simulateEffectiveSocialInput.m`：仅方位社会输入的核心仿真器。
- `runEffectiveInputSweep.m`：运行可恢复的参数扫描并写出 CSV。
- `analyzeEffectiveInputSweep.m`：根据全邻居基准标定成功条件，计算成功率和邻居阈值。
- `runEffectiveInputSmokeTest.m`：快速验证核心不变量并生成小型结果。

## 快速开始

在 MATLAB 中运行：

```matlab
cd('D:\Work\Study\UESTC资料存档\课程\通关项目\coding\effective-social-input')
runEffectiveInputSmokeTest
```

运行小型扫描：

```matlab
results = runEffectiveInputSweep('quick');
analysis = analyzeEffectiveInputSweep(results);
```

运行论文级第一阶段扫描：

```matlab
results = runEffectiveInputSweep('full');
analysis = analyzeEffectiveInputSweep(results);
```

`full` 配置计算量很大，默认包含多个群体规模、30 个随机种子和 8000 个时间步。建议先检查 `quick` 输出，再按临界区域缩小参数范围。

运行原文参数下的长时基准/零输入配对对照：

```matlab
controls = runEffectiveInputControls(1:10);
```

结果写入 `results/effective-input-controls.csv` 和
`results/effective-input-controls-summary.csv`。该对照用于校准评价窗口和随机初态方差，不能替代邻居数扫描。

运行 `N=20` 的长时邻居阈值先导扫描：

```matlab
results = runEffectiveInputSweep('n20-long');
analysis = analyzeEffectiveInputSweep(results);
```

运行 `N=20`、30 个种子、全部整数 `k=0...19` 的正式基线扫描：

```matlab
results = runEffectiveInputSweep('n20-final');
analysis = analyzeEffectiveInputSweep(results);
```

该扫描使用 8 个并行 worker，并写入独立的 `effective-input-n20-final.*` 文件，不覆盖先导结果。

## 两种输入制度

`original`：对给定群体规模保持每个邻居幅值 `h0 = hTotal/N`。邻居减少时总输入同步减小，代表真实漏失。

`fixed-total`：只要至少有一个邻居，就把原始全连接参考强度 `(N-1)hTotal/N` 平均分配给当前选中的邻居。它只用于拆分“总强度”和“角度采样”效应，不建议直接解释为生物机制或无人机控制律。

## 选邻策略

- `all`：所有其他个体。
- `random`：每个时间步随机选择 `k` 个邻居。
- `nearest`：选择周期距离最近的 `k` 个邻居。
- `balanced`：先选最近邻，再贪心选择与已选方位集合最分散的邻居；距离用于选邻，但进入神经场的信号仍只有方位角。

`balanced` 是可复现的算法基线，不声称等同于动物视觉注意机制。

## 成功阈值

`analyzeEffectiveInputSweep` 从 `all + original` 基准运行中取：

- 平均全局极化的第 5 百分位作为极化下界；
- 归一化平均两两距离的第 95 百分位作为凝聚上界。

将两个指标相对于这些基准分位数和基准中位数归一化，取两者的较小值作为联合基准等价分数；该分数的下侧 5% 顺序统计量作为成功下界。这样保证有限样本中至少约 95% 的全邻居基准运行通过联合条件，避免线性插值分位数和两个边际分位数造成通过率偏低。随后在每个条件下寻找成功率至少为 0.95 的最小邻居数。该阈值是条件阈值，不是跨速度、规模和参数的普适常数。

## 仅方位可见性扩展（当前已实现）

将 `cfg.visibilityMode` 设为 `bearing-fov` 可启用相对于个体当前航向的视场门控和随机漏检；默认 `fieldOfView=2*pi`、`detectionProbability=1` 时与理想模式等价。该扩展仍不向控制器提供距离测量。`occlusionEnabled` 目前是二维圆盘几何的简化遮挡基线，需要给出正的 `agentRadius`；它只用于生成是否能看到目标的门控，不把距离作为社会输入幅值。

尚未加入方位测量噪声、预测记忆和障碍物；确认无预测失效边界后，再比较无记忆、最后值保持、方位-角速度预测和不确定性卷积输入。

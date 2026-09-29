# B-E2 无遮挡群体实验

本目录保存 B 负责的科学问题一 E2 探索性实验。实验把 A 已完成的局部社会输入模型接入完整群体仿真，研究邻居数量、输入总幅值、方位结构和选邻策略对群体行为的影响。

## 实验范围

本批实验使用：

- 群体规模：N=20；
- 神经元数：100；
- 仿真步数：8000；
- burn-in：6000；
- 记录步长：100；
- 随机种子：1–20；
- 理想无遮挡方位观测；
- 等权输入；
- 不使用遮挡、观测缺失、KF 或预测补偿。

这批结果是探索性证据，不能作为最终确认性结论，也不能推出普适的邻居数量阈值。

## 实验内容

### E1：输入强度制度

比较：

- `original`：每个选中邻居使用固定幅值，总输入随邻居数变化；
- `fixed-total`：把预设总输入幅值平均分配给当前选中的邻居。

使用 (k=2,8,19)，选邻策略为 `random`。

### E2：邻居数量

使用 (k=2,4,8,12,19)，分别比较 `original` 和 `fixed-total`。

### E3：方位结构

固定 (k=8) 和 `fixed-total`，比较：

- `angle-concentrated`；
- `random`；
- `angle-dispersed`。

该实验用于探索方位覆盖和方位集中程度与群体响应的关系。当前 `random` 是动态随机选邻，不能等同于严格预设的随机方位结构。

### E4：选邻策略

固定 (k=8) 和 `fixed-total`，比较：

- `random`；
- `nearest`；
- `balanced`；
- `angle-concentrated`；
- `angle-dispersed`。

提交前应复核 `balanced` 和 `angle-dispersed` 的实现是否确实不同；当前结果中两者出现完全相同的汇总值，不能直接当作两个独立策略的证据。

## 文件说明

- `run_B_E2_experiments.m)：统一运行 E1–E4；
- `simulateBGroupE2.m)：B 的群体仿真器副本；
- `B_E2_runs.csv`：480 次运行的逐种子汇总；
- `B_E2_grouped_summary.csv`：24 个条件的均值、样本标准差和正态近似 95% 置信区间；
- `B_E2_runs.mat`：MATLAB 汇总表；
- `raw/`：pilot 或运行器保存的逐步仿真结果；
- `audit.json`：参数、种子、文件哈希和证据边界；
- `B_E2_RESULT_EXPLANATION.md`：简要结果解释。

## 主要字段

- `meanGO)：群体整体有序性；
- `meanLO)：拓扑局部有序性；
- `meanPairDistanceNorm)：归一化平均两两距离；
- `meanInputAmplitude)：平均社会输入总幅值；
- `meanAngularCoverage`：方位覆盖；
- `meanInputConcentration`：方位集中度；
- `socialUnionLccFraction`：社会交互网络最大连通分量比例；
- `socialAlgebraicConnectivity`：社会交互网络代数连通度。

统计单位是独立随机种子，不是逐帧时间点。置信区间仅用于探索性描述。

## 复现命令

当前运行器仍需先完成个人路径移除，之后应在仓库的 `effective-social-input` 目录中运行：

```matlab
run_B_E2_experiments('pilot')
```

正式验证不能直接复用种子 1–20。应在实验规则、样本量和统计终点冻结后，使用未参与方案选择的独立种子。

## 证据边界

E2 只提供无遮挡基线。它不能证明某个 (k) 是普适阈值，不能证明某种方位结构在遮挡环境下最优，也不能证明 KF 预测补偿已经带来闭环收益。遮挡和预测实验应在后续阶段单独设计和验证。
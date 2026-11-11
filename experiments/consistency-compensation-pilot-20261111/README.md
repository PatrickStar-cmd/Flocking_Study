# 周期遮挡下的社会输入补偿：小规模一致性实验

## 实验状态与问题

**已完成**：6 个开发种子、24 个独立测试种子，测试集共 384 个“种子 × 保留边数 × 方法”条件单元。目的是先确认遮挡令一致性退化，再比较补偿方法是否能恢复一部分损失。本实验没有估计社会输入的通用阈值，也没有检验补偿收益与社会输入的相关关系。

采用 20 架飞机的异中心神经场模型，无领航者、目标或观测噪声。每个种子完整观测预运行 200 步；同一种子的各方法从相同检查点、相同遮挡掩码开始，但之后各自闭环演化。主条件是每架飞机保留 5/19 条邻居边，次条件保留 8/19 条。每个遮挡周期结束后，所有方法都有 50 步完整观测；预测器用**自身轨迹**的这段历史初始化下一周期。

主指标 `C_occ_tail`：最后 1000 步内，仅在遮挡步计算每步 `C(t)=min(GO(t),LO(t))`，再取平均；数值越高，一致性越好。另记录 GO、LO、全时段一致性、方位预测误差及同分支反事实神经场输入误差。场景标签 `S0` 按共同检查点的可见边计算：`mean_i sum_{j∈visible_i} 0.012 exp(cos(theta_ij−phi_i)−1)`；这是队友的三因素经验量，不是高斯神经场输入。

## 冻结设计

种子 `2026111101–2026111106` 只用于 Fresh−Drop 退化筛选和预测参数校准；`2026111107–2026111130` 才用于最终比较。开发集预先规定：主条件 r=5 的 Fresh−Drop 平均差至少 0.05，且至少 5/6 种子为正。两种候选周期都通过，按预定规则锁定较短的 `occ100_obs50`：**遮挡 100 步、完整观测 50 步，重复 40 次，总计 6000 步**。开发集 r=5 的 Fresh−Drop 均值为 +0.3773，6/6 种子为正。完整方案见 [PROTOCOL.md](PROTOCOL.md)，筛选记录见 [开发集报告](results/development/analysis/README.md)。

测试方法为 Fresh、Drop、Hold、Linear、DampedHolt、SG、KF、StrengthOnly。Fresh 是无缺失参考，Drop 丢弃隐藏边；StrengthOnly 仅恢复可见输入的名义总权重，不预测隐藏方位。DampedHolt 和 KF 的参数只用开发种子遮挡前的真实历史校准一次。

## 独立测试集结果

下表为 24 个测试种子的配对差值，括号为按**种子整体**重采样的逐项 95% bootstrap 区间（4000 次）。均未作多重比较校正。

| 条件 | 比较 | `C_occ_tail` 配对均值差 [95% 区间] | 正值种子 |
|---|---|---:|---:|
| r=5，实际 `S0` 均值 0.03040 | Fresh−Drop | +0.2677 [+0.1370, +0.3995] | 19/24 |
| r=5 | KF−Drop | +0.1038 [+0.0130, +0.1992] | 16/24 |
| r=5 | KF−Hold | +0.0497 [−0.0567, +0.1524] | 14/24 |
| r=5 | KF−DampedHolt | +0.0012 [−0.1050, +0.1064] | 12/24 |
| r=8，实际 `S0` 均值 0.04864 | Fresh−Drop | +0.1725 [+0.0437, +0.2933] | 19/24 |
| r=8 | KF−Drop | +0.0459 [−0.0559, +0.1373] | 15/24 |

主条件 r=5 在独立测试集满足预定退化判据。其 `C_occ_tail` 均值为 Fresh 0.7214、Drop 0.4537、KF 0.5575；KF 的均值提升约补回 Fresh−Drop 均值差的 **38.8%**。在同一条件下，KF−Drop 的 GO 差为 +0.1035 [0.0133, 0.1991]，LO 差为 +0.0938 [0.0352, 0.1537]。因此，在这个特定周期遮挡场景中，KF 相对 Drop 的改善同时反映在 GO 与 LO 上。

**结论边界**：本次数据支持主条件下 KF 相对 Drop 的一致性改善，尚不能证明 KF 优于 Hold、Linear 或 DampedHolt；这些直接比较的区间跨过零。r=8 的 KF−Drop 区间也跨过零。队友在永久遮挡设置下提出的 `S0≈0.06–0.09` 只是候选敏感区，不能当成本周期遮挡实验已验证的阈值。详细数字和图见 [测试集报告](results/test/analysis/README.md)、[配对汇总](results/test/analysis/paired_summary.csv) 与 [主指标图](results/test/analysis/paired_primary.png)。

## 文件与复核

| 文件 | 内容 |
|---|---|
| `config.json`、`PROTOCOL.md` | 运行前固定的条件、种子、指标与筛选规则 |
| `runConsistencyCompensation2.m`、`buildBearingForecasts.m`、`calibrateBearingParams.m`、`movementConsistency.m` | MATLAB 仿真、补偿算法、参数校准与指标代码 |
| `results/development/` | 两种候选时序的逐种子原始轨迹、CSV 和筛选结果 |
| `results/test/forecast_params.json` | 单次开发集校准参数 |
| `results/test/seed-*/` | 24 个测试种子的 `summary.csv` 与 `traces.mat` |
| `results/test/analysis/` | 图、逐种子配对差、条件均值、区间、审计和详细说明 |

`results/test/completion.json` 确认测试种子 `2026111107–2026111130` 共 24 个、每个周期设置运行 6000 步。最终审计读取了每个 MAT 文件内嵌的校准参数，全部与同一个 `forecast_params.json` 一致；并从原始轨迹复算了 384 个条件单元的 CSV 指标和遮挡时序。详见 [audit.json](results/test/analysis/audit.json)。最终文件的 SHA256：

| 文件 | SHA256 |
|---|---|
| `config.json` | `f833bb1e6f6ca5070fe811d0c413dfbe76508e3333e592276242fc1b3d22f0b4` |
| `runConsistencyCompensation2.m` | `e09ef21dcc8022cac50418bf5cd44c2f8e7721437ec88a1eb95f9a4a18c17dd9` |
| `results/test/forecast_params.json` | `33ae168e995ef0ee8c61148d6762e9f82ad7ab35c50ed3279755ee22684d8710` |

复核现有结果：在本目录运行 `python analyze_results.py test`，所需 Python 依赖见 `requirements.txt`。若要从零重跑，应先复制整个实验目录，在新副本中按 `plan → smoke → development → python analyze_results.py gate → calibrate → test → python analyze_results.py test` 顺序执行。**不要在已有的部分测试结果上重新运行 `calibrate`**；参数文件会被覆盖，分析脚本会阻止混用不同校准参数的结果。

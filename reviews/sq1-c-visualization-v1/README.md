# 科学问题一：探索性数据可视化

来源：B提交`63d56b5`，C规范视图12条件×20种子，以及配对bootstrap分析。仅绘制已有结果，不运行新仿真。

| 图 | 内容 | 读图说明 |
|---|---|---|
| [01_strength_effects.png](01_strength_effects.png) / [SVG](01_strength_effects.svg) | fixed-total减original，在k=2/4/8/12下的GO、LO和对距差 | 浅点是20个种子的差，彩色点和线是均值及逐项95%区间 |
| [02_neighbor_response.png](02_neighbor_response.png) / [SVG](02_neighbor_response.svg) | 固定总幅值下k=2/4/8/12/19的响应 | 浅线为种子，深线为均值；不是时间轨迹或置信带，不表示中间k已被测试 |
| [03_policy_effects.png](03_policy_effects.png) / [SVG](03_policy_effects.svg) | k=8时各策略减random | balanced与angle-dispersed视为同一算法，仅展示一次 |

GO/LO正差代表有序度更高；对距负差仅代表更紧凑，不等于更安全。区间来自20个种子整组bootstrap，20,000次重采样，计算种子20260929；没有多重比较校正，不能解释为同时95%覆盖或独立验证。

角度策略与random均动态选邻；分散和最近邻使用距离信息。图不能单独识别纯方位结构的因果作用，不能给出最优k或普适阈值。若需要输入峰演化、轨迹动画或遮挡响应可视化，应先取得对应正式逐步数据，不能用当前均值反推轨迹。

复现：在仓库根目录运行 `matlab -batch "run('reviews/sq1-c-visualization-v1/plot_exploratory.m')"`。

绘图脚本核查输入行数、每条件seed集合，并独立检查图中配对均值。当前仅新增本目录；未改源数据、原代码、main或历史结果，未提交推送。可单独撤回本目录，共享提交后使用git revert。

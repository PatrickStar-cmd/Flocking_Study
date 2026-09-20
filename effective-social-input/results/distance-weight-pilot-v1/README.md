# 距离权重 pilot v1

运行前设计：种子131–135，N20，8000步、burn-in6000，random选邻k7/15；四组合：equal/exponential × original/fixed-total。共40次。其他神经参数保持默认。采用指数a(r)=exp(-r/(0.5L))；0.5L是本轮探索性尺度，不是生物标定值。观测为无遮挡理想距离+方位，没有遮挡或预测，不将其称为仅方位传感器。

original保持h0=.24/N；fixed-total将权重归一化后乘.228，用于拆分总幅度与相对权重。周期最短距离用于观测，新增meanDistanceWeight与离散环积分meanInputMass时间序列，保留原邻接图为二元选择图，不能称其为幅值加权图。当前保存的是平均权重而不是每条边历史；可据轨迹研究几何，但逐边权重存储仍未实现。

每个MAT保留完整result，逐次保存可续跑。描述GO/LO/距离、输入幅值与输入质量，不从5个种子报告可靠阈值。固定增益下性能变差可能仅由总量下降造成，须看fixed-total对照。核心动力学不改；equal分支保留原浮点求和路径，testDistanceWeight验证旧版逐项相等以及周期距离、零输入与归一化。旧版回归快照是simulateEffectiveSocialInputBeforeDistance.m，不作为实验运行入口。

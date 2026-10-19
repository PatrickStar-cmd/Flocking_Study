# 数据字段与时间索引

## 输入与运行

RUN_MANIFEST.csv每条唯一运行一行，858行。batch为development/main/duration6/duration30/early/scale40/scale60；Fresh跨时长复用，由分析构造引用，不重复计物理运行。seed/N对应inputs/共享NPZ，RNG为PCG64。

HDF5维度：states/xy=(801,N,2)，states/headings=(801,N)，states/u=(801,100,N)，states/velocity=(800,N,2)，states/social_field=(800,100,N)。位置是原始周期区域坐标。头角为rad；静止时保留前航向，GO单位速度为零。

edges/维度均(800,15N)，edge_i/edge_j是观察者/目标的0起始编号。有向信息方向为j→i。true是评估器的无噪声方位；observed仅合法新鲜边有限，缺测为NaN。center是实际注入方位，inactive时为0但权重为0，不可把该0当作有效观测。weight为实际幅值，active为输入激活，fresh为本步观测可用。filter_theta/omega、P11/P12/P22、age记录影子KF和缓存；drop/hold不因影子KF改变注入规则。

agent/S、S_oracle、A、C为(800,N)：S使用经验方向核exp(cos(delta)-1)，不是神经场高斯核；A是输入幅值系数之和；C=S/A。社会场不包括任务提示。任务方向可从CSV command_deg和配置复建。

checkpoint/对应本场景遮挡开始前状态，初态索引start=200或140。final/为800步终态，包括KF、缓存和趋势历史。history=(10,N,N)，第0项最新，count记录连续合法观测数。

## 时间

CSV step为0—799；本行输入基于states[n]，运动指标基于states[n+1]。time=(n+1-start)*.3为步末时间标签。主窗口阴影[0,18]对应行200—259。恢复首行260标18.3。main最末标签180；early因起点不同最末标签198。

## 指标

GO：全部N对象单位速度合矢量模/N；GO_followers分母N−L；GO_leaders分母L。LO取每帧位置的最近5他者加自身，和固定观测图无关。

fixed5_distance：遮挡开始前为每架确定最近5个他者身份，之后固定，100条有向距离（N20）均值；遮挡前NaN。不同遮挡起点的Frozen集合不同；距离未除以区域边长。pair_distance为全部无序对的周期距离平均；minimum_distance为本帧最小距离。

follower_turn_error_deg沿用原协议命名，表示个体目标方向绝对角差均值，并非转向机制对照。follower_turn_angle_deg为圆周平均方向；follower_resultant是相应合矢量模，接近0时平均角不适定。

S_abs_error：先每跟随者计算abs(S_i−S_oracle_i)，再平均；oracle使用分支自己的步初无噪声几何和基础幅值.012。

bearing_error_deg只对隐藏且激活项的实际center评估，无活跃项NaN。shadow_kf_rmse_rad与shadow_kf_coverage对全部本步缺测项的影子KF评估，包括不激活项，不能混用两个口径。gated_age_edges/gated_variance_edges记录当前方法适用的门控原因，可能重叠，不应相加当独立条数。

## 分析文件

window_metrics_by_seed.csv：case/seed/mode/window的指标、相对Fresh有符号差与逐步绝对偏差。

paired_differences_by_seed.csv：method−baseline，GO正值好、目标误差正值差。统计独立单位种子。paired_statistics.csv：均值、SD、中位数、范围、配对t区间、bootstrap区间、p值与指定比较族Holm值。development只描述范围，不给确认性检验。

time_series_all_conditions.csv提供所有场景的完整800步含引用Fresh曲线。工程诊断不剔除种子。all_seed_trajectories包含原始周期轨迹，绘图仅用NaN断开跨边界线段，不替换原始坐标或额外展开。

# 数据读取说明

## 有效批次目录

`results/pilot_seed1to5_v1/`是五种子有效结果。audit.json的COMPLETE及actualTrajectoryCount=155必须与RUN_MANIFEST一致。

- `references/fresh_seedS_AXXX.mat`：变量fresh（轨迹）和cpbase（6200完整checkpoint）。A114/228/456表示A0乘1000。
- `conditions/seedS_hH_STRUCTURE_AXXX/checkpoint6200.mat`：变量cp、cond；冻结隐藏集合及完整合法缓存。
- 条件目录下`drop.mat`、`drop-normalized.mat`、`hold.mat`、`KF-mean.mat`：变量seg。
- `same_state_replay.mat`：变量replay，含五种方法候选input field及指标，candidate名字经matlab.lang.makeValidName处理。
- `prediction_diagnostics.mat`：变量diagnostic；R_6/R_4/R_2分别对应R=1e-6/1e-4/1e-2。

## MATLAB示例

```matlab
x=load('results/pilot_seed1to5_v1/conditions/seed1_h4_random_A228/KF-mean.mat');
t=x.seg.data.step;
plot(t,x.seg.data.pairDistance);
```

## seg/fresh轨迹结构

- config：实际参数；startStep是逐步记录前的状态步数；endStep是结束步数。
- data.step：T×1；data.timeStart/timeEnd为模型时间，不是秒。
- state.position：2×N×(T+1)；state.u：M×N×(T+1)；state.heading：N×(T+1)。
- state.alpha：N×M，固定allocentric世界神经角。
- 第j个输入对应state.position(:,:,j)，运动后为state.position(:,:,j+1)。
- data.inputField/selfFreshField：M×N×T。selfFresh为当前**分支自己的**真值几何参考，不等于独立fresh轨迹；闭环配对指标使用外部fresh。
- data.GO/LO：群体有序度；LO使用含自身及5个最近邻的局部方向合矢量。
- data.pairDistance：周期最短两两距离均值除以L；nearestDistance同样归一化。
- k_selected/k_fresh/k_active/A/Mass/Peak、z1/z2、direction/directionValid、selfFreshRelativeL2为T×N。
- edgeWeight/age/lastBearing/lastWeight/observedBearing/observedWeight/injectedBearing/trueBearing/trueDistance为N×N×T。
- selected/fresh/hidden/active为N×N×T布尔。
- predictionBearing/Velocity/Variance/P12/P22/NIS/Age/Error为N×N×T。未选择边为NaN；误差是圆周差弧度。Variance=P(theta,theta)，不是角速度方差。
- priorBearing/priorVariance/priorError为更新前预测，可用于观察测量修正。
- predictionNIS：更新步创新平方/(预测方差+R)；缺失步为NaN，不能把NaN变成0求平均。
- finalCheckpoint：可完整续跑；包含inputState.tracks及sensitivityTracks。

## CSV

- closed_loop_summary.csv：7条件×5seed×5方法×4窗口=700行；fresh有重复引用，但仅15条独立reference轨迹。
- same_state_summary.csv：7×5×5=175行，每行是遮挡34步×20个体的描述均值/最大值。
- prediction_diagnostics_summary.csv：7×5×3R=105行；RMSE先按每seed全部隐藏边/步计算，再按seed汇总；不是2720个独立实验。
- hidden_sets.csv：7×5×20=700行；hiddenIds为空格分隔的实际ID；maskBearingStep=6200。
- RUN_MANIFEST.csv：155个唯一闭环运行键；软件检查、same-state replay及R诊断不冒充额外科学轨迹。
- analysis/：Python独立复算、配对bootstrap、中文说明、静态图PNG/PDF。

缺失掩码在事件窗口同时施加于所有智能体，每个体隐藏h个邻居。方位结构规则可能因群体本已集中而相似，解释时应检查hidden_sets中的实际R1/R2/空隙。

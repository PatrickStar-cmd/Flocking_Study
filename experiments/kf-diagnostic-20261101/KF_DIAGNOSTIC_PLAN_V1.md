# Flocking 卡尔曼滤波（KF）诊断实验规范与实施计划书 (v1 - 修订版)

## 1. 实验定位与科学边界

本实验依据主助手指令与项目协议实施，专注于以下科学与工程边界：
1. **纯方位观测约束**：仅利用方位角观测 $\theta$，严禁假设或声称恢复了二维位置 $(x, y)$。
2. **被动诊断原则**：在不改变任何动力学、数值轨迹或随机数消耗的前提下，记录神经场状态 $u$、逐步观测与主要跟踪器诊断量、前遮挡完整检查点及边界截断原因。
3. **因果配对干预**：在同一随机种子、同一起始状态及同一观测噪声流下，对比无遮挡（Fresh Control）、删除（None）、保持（Hold）、卡尔曼均值（KF-mean）与卡尔曼协方差热核卷积（Kalman）。
4. **主响应与次响应时窗定义**：
   - **主响应窗口为遮挡期 (`occ`)**：$t \in [t_0, t_0+D-1]$；
   - **次响应窗口为恢复期 (`rec`)**：$t \in [t_0+D, t_0+D+399]$；
   - **公共比较窗口 (`common`)**：$t \in [1401, 1920]$（520 步），采用矩形离散数值积分 $\text{AUC} = dt \times \sum_{t=1401}^{1920} y(t)$。
5. **闭环估计器自评估准则**：
   - 估计误差必须基于**每个算法自身闭环运行轨迹**中的真实物理方位 (`rec.trueBearing`)；
   - 严禁将分叉后的某一算法估计与另一算法轨迹真值混合比较；
   - 候选估计 (candidate) 与实际注入估计 (actual injected) 分开记录。
6. **证据等级明确**：本阶段为软件工程审计（P0）与探索性配对实验（P1 A包，n=5），不作为普适收敛性或闭环安全保证。

---

## 2. 状态时序与检查点索引定义

1. **时序状态推进与记录时点**：
   - `positionsRecord(:,:,t)` 与 `uRecord(:,:,t)` 记录的是**进入 step t 动力学推进前**的瞬时物理与神经状态（即 step $t-1$ 结束时的后验状态）；
   - 因此，**step $s$ 执行完毕后的神经状态对应 `uRecord(:,:,s+1)`**；
   - 仿真第 2200 步执行完毕后的最终状态对应 `finalNeuralState`；
   - 公共前遮挡事件检查点定义在 $t = t_0$ 推进前保存，明确标记其代表 `stateAfterStep = t0 - 1`。

---

## 3. 仿真物理与传感器参数

- **智能体与神经场**：$N = 20$, $M = 100$, $dt = 0.3$, $\text{speedScale} = 0.05$, $L = 1000$, $\beta = 1000$, 神经连接指数 $0.5$。
- **社会输入**：allocentric = true, original 单边增益，等权（equal），$h_t = 0.24$, $\sigma = 0.4$。
- **选邻拓扑**：fixed-id 循环拓扑，邻居总数 $k = 15$。
- **遮挡机制**：scheduledOcclusion = true，每智能体隐藏前 $h = 7$ 条边（8条新鲜可见，7条隐藏），遮挡起始步 $t_0 = 1401$。
- **传感器与卡尔曼滤波**：
  - 测量噪声：$\text{measurementStd} = 0.01\text{ rad}$，观测方差 $R = 1\times 10^{-4}$。
  - 过程噪声：白噪声角加速度谱密度 $q = 1\times 10^{-4}$。
  - 状态转移矩阵 $F = \begin{bmatrix} 1 & dt \\ 0 & 1 \end{bmatrix}$，过程噪声协方差 $Q = q \begin{bmatrix} dt^3/3 & dt^2/2 \\ dt^2/2 & dt \end{bmatrix}$。
  - 门限设置：$\text{maxAge} = 120$（仿真时间单位，对应 $120/0.3 = 400$ 仿真步），$\text{maxVariance} = 1.0\text{ rad}^2$。
- **仿真总长与时窗**：
  - $NSteps = 2200$，$\text{burnIn} = 1000$，恢复评测 400 步。
  - 遮挡时长 $D \in \{20, 120\}$ 步（分别对应仿真时间 $6.0$ 与 $36.0$ 单位）。

---

## 4. 分析输出表结构规范

1. **`summary.csv`（逐窗口长表）**：
   - 包含列：`commit`, `scenario`, `mask_rule`, `N`, `k`, `hidden_count`, `D_steps`, `eval_D_steps`, `dt`, `t0`, `seed`, `mode`, `window`, `window_start_step`, `window_end_step`, `window_length_steps`, `GO_mean`, `LO_mean`, `pair_norm_mean`, `nearest_norm_mean`, `active_neighbors_mean`, `input_mass_mean`, `bump_concentration_mean`, `speed_mean`, `result_path`, `common_auc_GO_discrete`, `common_auc_pair_discrete`。
2. **`contrasts.csv`（逐 seed 逐 D 逐 window 配对表）**：
   - 主比较：`window = occ` 下 `baseline = hold`, `method = kf-mean`；
   - 次比较：`window = rec` 下各对比；
   - 展宽机制比较：`baseline = kf-mean`, `method = kalman`；
   - 相对无遮挡损失：`baseline = fresh`；
   - 列定义：`seed`, `D`, `window`, `baseline`, `method`, `is_primary_response`, `benefitPair`, `benefitGO`, `benefitLO`。
3. **`estimator_diagnostics.csv`（估计器自身闭环长表）**：
   - 按 `seed × D × mode × age` 细分；
   - 列定义：`seed`, `D`, `mode`, `step_in_occ`, `age`, `n_hidden`, `mae_candidate`, `rmse_candidate`, `mae_injected`, `rmse_injected`, `mean_P11`, `coverage95_candidate`, `fresh_nis_mean`, `age_cutoff_count`, `p_cutoff_count`。

---

## 5. 计划要求的 8 张核心图谱命名

1. `P1_Fresh_Drop_Hold_KFmean_KFcov_GO.png`
2. `P1_PairDistance.png`（相对同 seed Fresh 差值）
3. `P1_LocalOrder.png`
4. `P1_ActiveNeighbors.png`（区分 fresh 与 active 邻居数）
5. `P1_InputMass.png`
6. `P1_BumpConcentration.png`（神经场波包集中度）
7. `P1_GainScatter.png`（主响应 occ 窗口逐 seed 散点/条形图）
8. `P1_FailureCases.png`（固定规则最小 seed，正确的 post-step 索引与神经输出）

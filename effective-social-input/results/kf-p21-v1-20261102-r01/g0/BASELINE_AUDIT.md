# P2.1 G0 基线审计与几何分层报告

**审计状态**：`PASS`  
**核对时间**：2026-11-02  
**上游结果目录**：`D:\大二下专业课\卓越计划\Flocking_Study_KF_diagnostic_20261101\effective-social-input\results\kf-p2-v1-20261102-reviewed`  
**审计输出目录**：`D:\大二下专业课\卓越计划\Flocking_Study_KF_diagnostic_20261101\effective-social-input\results\kf-p21-v1-20261102-r01\g0`  

---

## 一、P2 原始结果复算与 SHA 校验

- **复算行数**：`270 / 270` 行（覆盖全部 15 cases × 2 durations × 3 methods × 3 windows）。
- **最大 RMSE 差异**：`1.110e-16 rad`（严格 `< 1e-10`）。
- **最大 NIS 差异**：`7.105e-15`（严格 `< 1e-10`）。
- **输入 MAT 文件 SHA 校验**：`全部通过 PASS`。
- **回放 NPZ 文件 SHA 校验**：`全部通过 PASS`。

## 二、G0 探索性描述分位点（仅作先导展示，严禁用于开发集调参）

- **G-clean 有效样本边数**：`3677`
- **连续角速度变化率三分位**：
  - `q1 (33.3%)`：`0.002764 rad/模型时间`
  - `q2 (66.7%)`：`0.007997 rad/模型时间`
- **声明**：本分位阈值仅用于历史 5 个种子描述性审计，后续正式开发集实验必须独立从 12 个新开发种子真值分布中冻结，绝不交叉复用。

## 三、交付产物检查清单

1. `audit.json`：精确机器可读校验结果。
2. `geometry_events.csv.gz`：逐 case/edge/time 几何标签及最近镜像分量压缩明细。
3. `geometry_counts.csv`：各几何分层边数与样本占比统计。
4. `error_by_geometry_seed.csv`：按种子和几何分层方法误差（RMSE/MAE/P90/P95/P99/极端误差）。
5. `nis_tail_by_geometry.csv`：新鲜观测时刻按突跳邻域及洁净状态分层的厚尾 NIS。
6. `coverage_by_age.csv`：遮挡年龄分段经验覆盖率与方差失配比。
7. `p_bins.csv`：P 方差分箱校准分析。
8. `figures/fig_G0_1_trajectories.png`：典型几何清洁与镜像跳变方位轨迹。
9. `figures/fig_G0_2_error_nis_tails.png`：几何分层绝对误差尾部 ECDF 及 NIS 厚尾分布。

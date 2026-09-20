# 方位遮挡预测模型：文献约束与选择

检索日期：2026-09-09。用途：约束连续遮挡时间 `T` 实验中的预测状态和证据边界，不把通用目标跟踪算法直接移植为已验证方案。

## 直接相关证据

1. Li 与 Jilkov（2003），[Survey of maneuvering target tracking. Part I: Dynamic Models](https://doi.org/10.1109/TAES.2003.1261132)。该综述强调运动模型各自的假设、关系和利弊；CV、加速度相关模型、协调转弯和多模型方法必须按可用状态与机动类型选择。它支持比较多个运动先验，但不支持在只有相对方位状态时直接宣称二维协调转弯 IMM 可观测。
2. Singer（1970），[Estimating Optimal Tracking Filter Performance for Manned Maneuvering Targets](https://doi.org/10.1109/TAES.1970.310128)。用有时间相关性的机动过程建立可实现的 Kalman 模型，并把性能与机动尺度、观测噪声和采样率联系起来。本项目借用“一阶 Gauss-Markov 相关过程”的思想，但将其施加于相对方位角速度；这是一项待数据检验的局部模型，不等同于 Singer 的 Cartesian 加速度模型。
3. Kirubarajan 与 Bar-Shalom（2001），[Bearings-only tracking of maneuvering targets using a batch-recursive estimator](https://doi.org/10.1109/7.953235)。论文明确指出标准递归估计器会因缺少初始目标距离而收敛差，并以 batch ML-PDA 初始化、IMM/PDA 维护和本舰机动增强可观测性。它说明二维 bearing-only IMM 需要处理距离初始化、平台机动和数据关联；当前环形神经场接口没有这些状态，因此不能只因 IMM 常用于机动目标就直接套用。
4. Kurz、Gilitschenski 与 Hanebeck（2016），[Recursive Bayesian filtering in circular state spaces](https://doi.org/10.1109/MAES.2016.150083)。论文用 wrapped normal 和 von Mises 分布处理圆周状态，并讨论非线性系统、测量函数及加性/非加性噪声。它支持显式处理角度周期性，也提醒当前 wrapped innovation KF 只是单峰局部近似；当方差接近整圈时不能把高覆盖率解释成有方向分辨力。
5. Sinopoli 等（2004），[Kalman Filtering With Intermittent Observations](https://doi.org/10.1109/TAC.2004.834121)。把观测到达建模为随机过程，分析间歇观测下误差协方差的统计行为及临界到达率。它支持遮挡时只做预测、协方差增长和有限 blackout 的条件性结论；不保证模型失配时预测均值优于最后值保持。

## 对当前模型的结论

当前传感器向控制接口提供带噪相对方位；距离只在距离权重机制实验中作为理想观测或最后观测权重使用。若改为 Cartesian 相对位置 EKF/IMM，需要新增距离或通过平台机动证明可观测性，还会改变研究问题。因此第一候选保持状态 `x=[theta, omega]`，将常角速度改为一阶 Gauss-Markov 角速度：

`d theta = omega dt,  d omega = -(1/tau) omega dt + sqrt(q) dW`。

离散均值转移为 `omega(k+1)=exp(-dt/tau) omega(k)`，旧角速度信息随遮挡时间衰减；`tau -> Inf` 精确退化为现有常角速度模型。过程协方差采用连续模型的精确离散化。角度更新仍使用 wrapped innovation，神经场接口仍使用预测方位和方位方差。

## 模型选择与停止条件

- 只用高频参考轨迹种子 138–139，在 `T={3,10.2,30}` 上按预测 MSE 选择 `q` 与 `tau`。`T={60,120}` 不参与选择，因为旧模型的 95% 区间已覆盖整圈，且在线分支会受 age/variance 截止。
- 种子 140–141 已被查看，只能作为开发评估，不能作为新确认性证据。
- 候选衰减模型必须消除 `T=3` 的明显退化，并在 `T=10.2,30` 相对 hold 和旧 CV 给出一致改善，才接入检查点闭环。
- 若通过开发评估，再生成从未使用的新种子作确认；若未通过，停止群体 `T` 扫描并把“预测不优于保持”作为模型结果报告。
- IMM 只在单模型仍出现可分离的静止/转动状态，且新鲜数据支持状态切换收益时再考虑；当前不实现。

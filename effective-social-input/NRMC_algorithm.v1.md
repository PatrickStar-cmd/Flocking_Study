# NRMC v1: Neural Response-Aware Multi-Model Compensation

状态：候选补偿算法定义；尚未实现闭环仿真或取得恢复效果数据。

本文件与 `bearingTrackStep.m`、`bearingTrackStepGaussMarkov.m`、`occlusionInputStep.m` 放在同一研究目录，供后续以相同观测历史实现和比较。NRMC 的输入仅包括当前个体的神经状态、可见方位、可见性、历史观测、原始邻居权重和已知模型参数。遮挡时不读取邻居的仿真真实位置、方位、速度或距离。`GO`、`LO` 只用于离线评价。

## 1. 符号与原始输入

个体为 `i`，邻居为 `j`，神经元偏好方位为 `alpha_l`，可见与遮挡的已选择邻居集合分别为 `V_i(t)` 和 `M_i(t)`。`a_ij` 是原模型中的单邻居输入幅值；等权对照采用与原模型一致的 `a_ij`。若传感器输出相对方位 `delta_ij`，则异中心方位 `y_ij = wrap(delta_ij + psi_i)`，其中 `psi_i` 为自身已知朝向，`wrap(x)=atan2(sin(x),cos(x))`。

原始核与可见输入为

```
K_sigma(phi) = exp[-wrap(phi)^2/(2 sigma^2)]
I_vis,i(alpha) = sum_{j in V_i} a_ij K_sigma(alpha-y_ij).
```

实现时核采用现有 `occlusionInputStep.m` 的周期角距离和神经元网格。当前可见观测使用实际测量值，保留原始输入幅值。

## 2. 三种边级方位模型

每条已观测边维护 H（角度随机游走）、C（常角速度）、D（阻尼角速度）三个预测器。遮挡时只预测；重新可见时以圆周创新更新。首次观测初始化角度为实测值，角速度为零，角度方差为 `R`，角速度方差为 `1`。

```
H: theta[t+1] = theta[t] + xi, Var(xi)=q_H dt.
C: x=[theta,omega]^T, F_C=[[1,dt],[0,1]],
   Q_C=q_C [[dt^3/3,dt^2/2],[dt^2/2,dt]].
D: rho=exp(-dt/tau), F_D=[[1,tau(1-rho)],[0,rho]],
   Q_D=q_D [[tau^2 dt-tau^3(2(1-rho)-(1-rho^2)/2),tau^2(1-rho)^2/2],
             [tau^2(1-rho)^2/2,tau(1-rho^2)/2]].
```

共同预测为 `x_m^- = F_m x_m`, `P_m^- = F_m P_m F_m^T + Q_m`；H 的状态为一维。可见时使用 `e_m=wrap(y-H_m x_m^-)`、`S_m=H_m P_m^- H_m^T+R`、`K_m=P_m^- H_m^T/S_m`、`x_m=x_m^-+K_m e_m`，协方差采用 Joseph 更新。角度结果映射回 `[0,2pi)`。D、C 的基本状态更新分别对应仓库中的 `bearingTrackStepGaussMarkov.m` 与 `bearingTrackStep.m`。

每条边维护模型概率 `pi_m`，初值均为 `1/3`。每次获得新观测时先设 `pi_m^-=(1-epsilon)pi_m+epsilon/3`，再以包裹正态预测似然 `WN(y;mu_m^-,S_m)` 更新并归一化。遮挡时冻结 `pi_m`，传播各模型状态。边级预测方位分布为

```
p_ij(theta) = sum_{m in {H,C,D}} pi_ij,m WN(theta;mu_ij,m,V_ij,m).
```

`V_ij,m` 是预测角度方差。包裹正态密度按 `sum_{k in Z} Normal(theta+2pi*k;mu,V)` 计算，数值实现截断到可忽略的尾部。

## 3. 候选神经场输入

对 `m in {H,C,D}` 计算单模型缺失输入：

```
B_i^m(alpha) = sum_{j in M_i} a_ij integral K_sigma(alpha-theta)
                WN(theta;mu_ij,m,V_ij,m) dtheta.
B_i^mix(alpha) = sum_{j in M_i} a_ij sum_m pi_ij,m integral
                   K_sigma(alpha-theta) WN(theta;mu_ij,m,V_ij,m) dtheta.
```

积分在圆周上进行，可复用现有输入核网格和圆周卷积。候选集合为 `A_i={gamma B_i^b : b in {H,C,D,mix}, gamma in {0,0.05,...,1}}`。零输入重复项去重，最多 81 个不同候选。

## 4. 神经响应目标与选择

使用当前仿真器的同一步神经场更新，`M` 为神经元数、`u_i` 为当前状态、`W` 为环连接矩阵、`h_b` 为全局抑制、`I_ext` 为该个体实际收到的外部刺激：

```
u_i^+(B) = u_i + dt[-u_i + W tanh(beta u_i)/M - h_b 1
                   + I_vis,i + I_ext,i + B].
g_i(B) = v_s sum_l [tanh(beta u_i,l^+(B))]_+
         [cos(alpha_l), sin(alpha_l)]^T.
```

从每条边的 `p_ij` 产生 `L` 组条件独立的可能方位 `theta_ij^(s)`；同一组中，每条遮挡边抽取一次。第 `s` 组的真实输入假设为 `B_i^(s)(alpha)=sum_{j in M_i} a_ij K_sigma(alpha-theta_ij^(s))`，对应响应为 `g_i(B_i^(s))`。计算 `vbar_i=L^(-1)sum_s g_i(B_i^(s))`，选择

```
B_i^* = argmin_{B in A_i} ||g_i(B)-vbar_i||_2^2,
I_NRMC,i = I_vis,i + B_i^*.
```

损失相同先取较小 `gamma`，再按 `H,C,D,mix` 排序。所有候选及样本共用同一当前 `u_i`、`I_vis,i` 和 `I_ext,i`。候选评价没有读取该时刻被遮挡的真实方位；真值仅供离线误差计算。样本生成使用单独记录的随机流。

该规则等价于在候选集合中最小化预测分布下的**一步神经运动响应平方误差**，因为响应平方误差的样本平均可分解为 `||g_i(B)-vbar_i||^2` 加一个与候选无关的常数。长期群体协同改善需要闭环实验测量。

## 5. 首版参数与实验接口

| 参数 | 开发初值 | 说明 |
| --- | ---: | --- |
| `dt,M,beta,sigma,W,v_s,a_ij` | 与原模型一致 | 对照组共享 |
| `q_H` | `1e-4` | 角度随机游走噪声强度 |
| `q_C` | `1e-4` | 常角速度过程噪声强度 |
| `q_D` | `1e-5` | 阻尼角速度过程噪声强度 |
| `tau` | `10` | 角速度衰减时间常数 |
| `R` | 与观测噪声匹配 | 数值下限 `1e-8` |
| `epsilon` | `0.02` | 模型概率混合率 |
| `L` | `64` | 每步响应场景数 |
| `Gamma` | `0:0.05:1` | 注入比例网格 |

这些数值是开发初值。开发集选择参数后冻结，在独立种子上评价。候选比较至少包括 Drop、Hold、CV-KF、阻尼 KF 加不确定性输入、IMM、NRMC 和完整观测。消融：`NRMC-Mix` 固定使用混合输入且 `gamma=1`；`NRMC-D` 仅允许阻尼模型的 21 个缩放候选。

主恢复指标沿用实验方案中 `C=min(GO,LO)>=0.8` 连续 100 步；同时报告恢复时间、恢复后达标时间和运行开销。弱化和持续退化分层报告。任何社会输入阈值均与种群规模、神经参数、选邻规则、观测模型和评价时长同时列出。

## 6. 验证顺序

1. 复核无遮挡时 NRMC 的新鲜输入与原模型逐步一致；验证遮挡接口没有向 NRMC 传递隐藏真值。
2. 用固定小例子复核三模型预测、概率归一化、核积分和 81 个候选的独立计数；用重复种子复核随机流。
3. 在同一初始检查点、遮挡计划和可用观测历史下与所有基准配对；各方法给出同等开发调参预算。
4. 冻结参数后在新种子上报告持续恢复与成本；用同一原始轨迹文件重算关键计数。

证据等级：算法设计。当前文件没有仿真结果、恢复收益或稳定性结论。该方法仅推断方位分布和神经输入；只有方位观测时，二维邻居位置仍需额外传感或可观测性条件。

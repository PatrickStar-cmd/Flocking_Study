# B → C 数据交付契约 v1

格式优先MAT（保留双精度数组）+CSV摘要+JSON审计。NA必须与真实0区分；复数保存实部/虚部。不能仅交汇总CSV。

交付阶段、明确数量、数组维度说明及文件包见 [B_DELIVERY_CHECKLIST.md](B_DELIVERY_CHECKLIST.md)。审计包须同时包括metadata.json与audit.json、字段字典、加载示例及实际运行状态表；计划清单不能作为已完成记录。

协议1.2的输入映射以 A_B_INPUT_INTERFACE.md 为准。每步个体另保存weightEntropy、神经活动质量、bumpConcentration、decoded_direction_defined及分析用解码角；原始神经解码值与模型运动航向单独保留。输入方向与神经方向的有效性不可混用。无输入时熵/方向/空隙为NA，零场和矩占位须配no_input标记。原A的angularEntropy/hasDirection若保留，必须标注legacy，不能作为正式分析字段。

| 层级 | 必须字段/结构 | 采样要求 |
|---|---|---|
| 批次 | protocol_version、协议/配置SHA、代码SHA、实际依赖路径/哈希、MATLAB版本、命令、worker、seed状态、错误记录 | 每批次 |
| 唯一运行键 | seed、checkpoint_id、mask_rule、duration_steps、mode、run_id；fresh有独立reference_id | 显式标注共享fresh，不复制为独立样本 |
| 配置 | N/M/dt/L、k/h/A0/w0、核sigma、神经/运动参数、无噪声、等权、固定ID规则、各时间窗口 | 不依赖读者从脚本默认值猜测 |
| 检查点 | 位置、航向、神经场u、偏好角alpha、全体目标ID、缓存/年龄、RNG/驱动索引、绝对步号 | 第6000步末和事件前第6200步末 |
| 掩码 | 每个体S_i/H_i/剩余ID、最后观测方位、mask生成器/seed、候选排序与tie规则、R1/R2/最大空隙 | 三mask均在共同前缀冻结，全部T共用 |
| 每步群体 | 局部/绝对步号、时间、phase、GO/LO、周期平均对距及最近距离（均注明归一化）、输入汇总 | 600个分支步骤每步保存 |
| 每步个体 | 位置、航向、k_selected/k_fresh/k_active、A、Mass、输入峰值、z1/z2实虚部/模、direction_defined、空隙、无输入/归一化不可用标记 | 每步每个体，不能只存群体平均 |
| 每步有向边 | observer_id/target_id、selected/hidden/fresh/active、observed_bearing、last_bearing、injected_bearing、weight、age | 可用N×N数组；隐藏observed_bearing=NaN |
| 每步输入场 | injected_field[M,N]、该分支自身fresh oracle[M,N]、神经网格alpha或固定网格引用 | 若空间大，用无损压缩，不通过降低采样率替代 |
| 参考同状态回放 | fresh参考状态下四种候选输入场、mask与缓存引用、各输入误差 | 区别于闭环自身oracle，不回流控制器 |
| 独立真值日志 | 分支自身真实方位/周期距离（仅评估）、几何时间戳 | 与控制器观测分开命名/存储，不作为控制输入 |
| 摘要 | 每seed/mask/T/mode/窗口的主次响应与reference_id | 从逐步数据自动生成，可由C独立复算 |

常规轨迹记录从分支初态j=0开始，共601帧；输入/输出指标j=1…600，共600点。位置/航向的状态时刻是步末，输入方位是步初，必须注明。输入字段的计算时点不可与步末几何混用。

每步保存完整神经场不是常规分析的强制项，但验收中的连续运行/检查点恢复测试必须比较逐步神经场；常规实验至少保留两级检查点和最终神经状态。运行失败时保留最后有效状态。

交付清单：README、metadata.json、配置快照、运行清单、原始MAT、summary.csv、日志、tests.json、文件SHA256清单。原目录存在则拒绝覆盖。C验收行数、键唯一性、seed完整性、phase边界和至少一个主终点的独立复算。

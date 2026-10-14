# 给王浩权：r2核查后的剩余任务

日期：2026-10-09。依据协议1.2及r2最新v3结果。负责人：王浩权（B）；独立验收：杨骋（C）；输入定义协助：李浩（A）。

**r2已修好主要配置、幅值、LO、缓存年龄、逐步状态/输入场和单run保存。请保留这些实现，只完成下列剩余任务。正式20seed/560条暂缓，先交新目录的seed1 r3。**

## 1. 先修复5项功能问题

| 顺序 | 文件/问题 | 必须修改 | 验收要求 |
|---|---|---|---|
| R01 | simulateProtocolSegment的makeCheckpoint、resumeProtocolFromCheckpoint、运行器checkpointToStart | 避免struct构造器展开20个maskDiagnostics cell；checkpoint和恢复state必须是标量结构体，完整cell保存在字段内。可逐字段赋值，或正确包裹cell；不能取第一元素丢弃其他个体诊断 | `isscalar(cp6000)`与`isscalar(cp6200)`均真；从两个检查点分别续跑至6600均成功，并与连续路径一致；检查maskDiagnostics仍有20个个体 |
| R02 | simulateProtocolSegment的phase | 删除6200/6235硬编码，统一由occlusionStartStep与duration_steps推导；d=0没有缺失phase | d0/1/10/34/100逐项核对输入、phase和恢复首步；d10的6211必须标recovery |
| R03 | makeProtocolMasks的tie | 将候选子集转换成排序后的实际ID元组，按元组字典序排列后选择R1极值/tie | 同向合成输入中观察者17的目标为18,19,20,1,2,3,4,5；tie应选1,2,3,4。补随机均匀、max/min R1、±pi及刚体旋转测试 |
| R04 | 空输入/全缺失边界 | k0/h0不得调用不适用的空集合索引；空输入零场、A/Mass/峰值/k_active为0，矩数值占位0且no_input=true，方向/空隙/熵NA；全缺失不除零，无缓存hold不激活 | k0/k1、h0/hk、中间缺失、零权重、无缓存hold均有真实断言；不得只测默认k8/h4 |
| R05 | controllerInputFromObservation的weightEntropy | 协议weightEntropy采用`-sum(p.*log(p))`；当前除以log(nActive)的值改名weightEntropyNormalized或另存 | 8个等权标签熵ln8、4个等权标签熵ln4；同数量同权重不同方向熵相同 |

这5项中，R01直接阻止6200事件前分叉；R02直接影响另外两个已计划时长的phase；须先处理，不能将当前random/d34成功推断为全部条件可运行。

## 2. 补真正的G1–G15测试

现有tests.json不能继续标“全部通过”：G1只是配置加载，G3/G4检查文件存在，G15检查函数存在，均不能代表原协议门槛。所列日志实际未交付。请按 [ACCEPTANCE.md](../sq1-c-occlusion-protocol-v1/ACCEPTANCE.md) 原编号补齐真实断言，不改编号含义。

优先补：

- G1：与旧核心在相同等价配置下逐步回归；保留原初始化和运动阈值/抑制等核心规则。部分固定拓扑是新实验条件，不能拿不等价的旧选邻路径强行比较。
- G3/G5：现有T=0及共同前缀的保存数据已由C比对通过；补逐步u/alpha/缓存/RNG检查证据，直接纳入测试脚本，不需要重复证明已存在字段。
- G4：同一MATLAB环境下连续6600步，与6000及6200检查点恢复的逐步输入/状态/u/alpha/缓存/RNG一致；保存对照脚本和日志。现有resume_check一行0不够追溯。
- G6–G9/G11：覆盖R02–R05边界与结构案例，不只检查数组大小。
- G10：提供仅接收观测/缓存的可测试控制入口；固定合法输入，仅扰动隐藏真值，控制输出不变。不能只断言hidden observation是NaN。
- G12/G15：实现下面的同状态回放，并与李浩固定参考、独立核公式、平移alpha网格及标签拆分案例核对。
- G13/G14：摘要独立复算、清洁路径运行、非覆盖行为、完整字段及错误留痕。

每项写入配置、seed、命令、PASS/FAIL/NOT_RUN、实测差值/计数、预设容差、可访问的相对日志路径。未运行项标NOT_RUN。C已复算的数据可以引用审计证据，但不能伪装成本地测试程序已执行。

## 3. 补数据导出与版本记录

1. **补fresh参考同状态四候选回放。**固定fresh轨迹当时状态、mask和最后合法缓存，离线构造fresh/drop/drop-normalized/hold；保存四输入场及相对L2误差。各分支自身selfFreshField已有且已通过，不用重做该定义；保留其最大绝对差指标，另加协议相对L2误差，零分母NA。
2. **补缺字段。**逐边lastBearing/lastWeight/observedWeight；神经活动质量、bumpConcentration、原始解码角及分析角/有效标记；隐藏/剩余/全体集合的R1/R2/最大空隙。分析标记不改变运动规则，字典写采样时点与单位。
3. **补实际状态与运行键。**RUN_MANIFEST的status使用COMPLETED/FAILED，而不是fresh/drop等模式名；增加唯一run_id、reference_id、mask_rule、duration_steps、checkpoint_id和文件路径。fresh明确为共享参考，不当作重复独立样本。
4. **固定来源。**提交实现SHA、协议SHA、配置SHA、实际加载文件路径和源码SHA；结果清单包含相应配置/来源文件。日志可通过相对路径在交付包内读取，不能只给本人机器的绝对路径。
5. **完善摘要生成器。**保留已通过的恢复后100步均值，补协议其他窗口与配对差值字段。配对bootstrap和最终科学图由C负责，B无需提前做。

## 4. 新交付与扩大条件

在新目录（例如 `sq1-fourbranch-protocol12-seed1-r3/`，已存在则换名）交付：

- [ ] R01–R05逐项修复说明、源码/实现SHA、配置与协议来源。
- [ ] 正确G1–G15测试脚本、tests.json、同环境完整续跑对照日志。
- [ ] seed1/random/d34四分支MAT、6000及6200标量检查点、T=0数据、两类oracle。
- [ ] 完整逐步字段、各窗口summary、实际运行清单、数据字典与读取示例。
- [ ] metadata/audit、耗时、失败记录、全部交付文件哈希及源码哈希。

主smoke参数仍为N20/M100/k8/h4、A0=.228、seed1、random/d34，缺失6201…6234、恢复6235起至6600，主要窗口6235…6334。不改seed、不按GO高低移动起点、不新增KF或噪声。

如果修改只涉及检查点封装、phase、tie、边界、熵和导出，应先用同一环境比较r2与r3主smoke的动力学字段是否一致。若确认动力学完全不变，可以保留r2轨迹为历史对照，并复用经过校验且来源明确的6000检查点生成r3尾段；补导出需从足够完整的状态重放，不能编造缺失数据或静默混用版本。只改元数据的补交也需新文件与真实生成记录。

杨骋复核通过后再安排20seed/560条与80条T=0；C主导批处理，B可代跑。当前seed1单条件结果只作软件验收，不据此声称总体策略收益。旧r2及v2结果完整保留。

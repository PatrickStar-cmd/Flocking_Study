# B：四分支缺失实验交接入口

发布分支：YCwork。协议科学设计版本：1.2。当前为实施与验收任务交接，尚未运行560条正式开发集实验。不要把计划清单视为已完成结果。

## 最新进度与当前任务（2026-10-10发布）

**基础r3软件已通过杨骋独立验收；你现在应开始批处理适配，无需重复修复旧验收中已解决的问题。** 固定基础提交为`feature/sq1-r3-software-package@54fb85d91151e0b28e8e570318b4a9937d6a0425`。

请先读新的 [批处理任务单](sq1-c-batch-freeze-v1/B_HANDOFF.md)、[配置包说明](sq1-c-batch-freeze-v1/README.md) 和 [批处理验收规则](sq1-c-batch-freeze-v1/BATCH_ACCEPTANCE.md)。r3证据见 [独立验收报告](sq1-b-r3-audit-v1/REVIEW.md)。

先提供多seed/mask/duration入口、动态键与汇总、共享基线/前缀/参考及完整离线回放；按新配置交seed1的28条主轨迹和4条T=0。杨骋核查后固定实际适配提交，再安排20seed全量执行。基础r3提交不能代替尚未交付的批处理实现SHA。

下列协议文件保留原字节，里面“待实现验收”和`execution.allowed=false`是当时的历史状态。当前配置允许适配检查，但`full_batch_allowed=false`，仍不能直接运行完整560条。

## 必读文件与顺序

| 文件 | 用途 |
|---|---|
| [B_DELIVERY_CHECKLIST.md](sq1-c-occlusion-protocol-v1/B_DELIVERY_CHECKLIST.md) | 先读：B做什么、先交什么、最终交付什么 |
| [PROTOCOL.md](sq1-c-occlusion-protocol-v1/PROTOCOL.md) | 研究问题、参数、时序、模式和响应定义 |
| [A_B_INPUT_INTERFACE.md](sq1-c-occlusion-protocol-v1/A_B_INPUT_INTERFACE.md) | A参考公式的适配、空输入及方向有效性 |
| [DATA_CONTRACT.md](sq1-c-occlusion-protocol-v1/DATA_CONTRACT.md) | 原始数据字段与采样要求 |
| [ACCEPTANCE.md](sq1-c-occlusion-protocol-v1/ACCEPTANCE.md) | G1–G15测试要求 |
| [EXECUTION_ORDER.md](sq1-c-occlusion-protocol-v1/EXECUTION_ORDER.md) | 实现、验收、批处理、统计顺序及负责人 |
| [protocol.json](sq1-c-occlusion-protocol-v1/protocol.json) | 机器可读设计参数；execution.allowed目前为false |
| [RUN_MANIFEST.csv](sq1-c-occlusion-protocol-v1/RUN_MANIFEST.csv) | 560个正式计划键，共享fresh参考 |
| [T0_TEST_MANIFEST.csv](sq1-c-occlusion-protocol-v1/T0_TEST_MANIFEST.csv) | 80个T=0软件测试计划键 |
| [check_design.py](sq1-c-occlusion-protocol-v1/check_design.py)、[design_check_v4.json](sq1-c-occlusion-protocol-v1/design_check_v4.json) | 当前设计算术与文档快照校验，不是仿真验收 |
| [A复核报告](sq1-a-recheck-992e71a-v1/REVIEW.md) | A已完成成果、可用范围及遗留问题 |

## 附带证据和复现文件

A复核目录同时附local_recheck.json、source_manifest.json、matlab.log、recheck_local.m、fetch_snapshot.py及source/固定源快照。源快照来自Hao_L@992e71a，只作参考/复现，不是B已集成的生产代码；禁止在该快照中修改历史结果。

协议目录的design_check.json、design_check_v2.json、design_check_v3.json保留为历史检查记录；当前规模为20seed、560条，不能从旧版112条记录读取执行规模。现有文档里“本地未提交/未推送”的表述是各次生成时的历史状态，本交接入口及实际Git历史反映发布状态。

在仓库根目录执行 `python reviews/sq1-c-occlusion-protocol-v1/check_design.py` 可检查计划数量、时间边界和参考文件；脚本会重新生成design_check_v4.json，重新审计请在新副本中运行并保留已发布记录。A的MATLAB复核也会写local_recheck.json，重新运行请使用独立目录。源快照已附，无需先下载才能阅读或复现；fetch_snapshot.py是原复核时的获取辅助脚本，来源以已保存manifest和固定SHA为准。

## B下一步

按上方最新任务完成批处理适配并提交完整条件包；G1–G15的基础软件验收已经通过，不再将旧任务清单当作当前未完成事项。适配通过后C主导正式批次，B可代跑并按契约交付。最终统计、bootstrap和科学解释由C负责。本次发布不修改A/B工作分支或main，不代表已获得总体科学收益。

如需撤回，使用git revert本次发布提交，不删除历史实验或改写公共历史。

# 杨骋最新验收与配置发布

发布日期：2026-10-10；发布分支：YCwork。内容截止2026-10-09验收与配置冻结。

## 发布内容

| 材料 | 结果与用途 |
|---|---|
| [首轮验收](../sq1-b-local-smoke-audit-v1/REVIEW.md)及整改清单 | 保留早期错误定位依据，属于历史软件审查 |
| [r2验收](../sq1-b-r2-audit-v1/REVIEW.md)及后续任务 | 保留修复进展与尚存问题，属于历史软件审查 |
| [r3验收](../sq1-b-r3-audit-v1/REVIEW.md)与脚本/日志/复算结果 | 固定54fb85d，96项原交付哈希匹配，G1–G15独立重跑通过 |
| [批量配置冻结包](../sq1-c-batch-freeze-v1/README.md) | 560主键、80个T=0键、逐项配置与窗口、源码/输入锁定 |
| [王浩权当前任务](../sq1-c-batch-freeze-v1/B_HANDOFF.md) | 批处理接口及seed1完整28+4条件包 |
| [杨骋后续验收规则](../sq1-c-batch-freeze-v1/BATCH_ACCEPTANCE.md) | 适配验收后才固定实际提交并放行全量实验 |

同步更新reviews主入口和B交接入口。已有9月29日探索分析图表继续保留，本次没有重新生成图表或运行仿真。

## 当前证据与历史记录

本次发布的是软件验收与执行配置准备，不是20seed缺失实验的科学结论。基础r3已通过，批处理实现尚未交付；冻结配置允许seed1适配检查，完整560条尚未放行。全量T=0计划80条，原r3只完成seed1的4条；不能把清单当作完成记录。

三轮审查报告的“未提交”“未推送”、当时HEAD和失败日志均保留为生成时历史。冻结包参考目录是原字节快照，不回写历史状态，不改哈希来掩盖来源。本发布说明与实际Git提交用于解释当前发布状态。

`sq1-b-r3-audit-v1/lfs_pending.json`是下载完成前的历史诊断；最终状态见同目录`source_audit.json`，8个真实MAT已获取且96项原清单匹配。G14故意触发的无效掩码报错属于预期负测试；其他早期失败日志同样保留。

## 下载与复现

阅读本次报告、配置和日志无需下载大文件。r3下载源码/原MAT、临时下载目录和C独立重跑原始结果由该审查目录`.gitignore`排除；没有在YCwork重复上传B的MAT/LFS内容。

在新克隆中校验配置：

```powershell
git fetch origin feature/sq1-r3-software-package
python reviews/sq1-c-batch-freeze-v1/validate_freeze.py
```

校验器依赖固定提交54fb85d及其祖先Git对象；不需要实际MAT。`prepare_freeze.py`是原构建脚本，不能对已存在的冻结包重复生成。

需要重新读取或生成r3验收数据时，在新的审查副本/新克隆中使用下列入口，并保留已发布证据：

```powershell
git lfs fetch origin 54fb85d91151e0b28e8e570318b4a9937d6a0425
python reviews/sq1-b-r3-audit-v1/fetch_package.py
```

之后按`audit_delivery.py`、`audit_saved_mat.m`和`reproduce_r3.m`核查。后两者需MATLAB；`reproduce_r3.m`生成新的软件结果，并不执行560条。部分原审查脚本会写同名JSON/日志，因此已有证据环境不要直接重跑覆盖。首轮/r2输入位于交付者提供的本地路径，未作为公共数据重复发布；其哈希与来源仍在审查记录中，不能假定新克隆已有这些旧MAT。

## 发布检查与回退

[verify_publication.py](verify_publication.py)检查暂存/指定提交的文件范围、原字节、读取入口链接、冻结清单哈希、运行计数和r3验收证据。[verification.json](verification.json)保存本次实测记录；检查脚本可以通过`--revision HEAD`对发布提交只读复验，不重写记录。

没有修改main、A/B分支、allocentric-flocking-original、历史实验结果或冻结协议。推荐将本发布作为协作入口使用，暂不据此自动合并main或执行完整科学批次。如需撤回，`git revert <本次发布提交>`；保留历史实验和公共提交，不强推。

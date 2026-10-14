# 实验入口

本目录提供实验索引；代码与历史结果保留原路径，避免破坏已有图表链接、源文件哈希和审计记录。

## 20对象固定3＋17角色的转向与连续方位缺测开发实验

- 实验标识：`informed-turn-N20-T18`；日期：2026-10-14。
- 证据等级：三个开发种子的探索性实验，不是确认性验证。
- [独立复刻协议，无结果版](../replication-protocols/20261014-informed-turn-n20-t18/实验复刻协议_无结果版.md)：包含参数、公式、随机数顺序、观测图、缺测时序、KF与指标定义。
- [实现与运行说明](../challenge2-informed-turn-v1/README.md)。
- [结果、图表与证据限制](../results/challenge2-informed-turn-v1-20261014/README.md)。
- [最终批次参数](../results/challenge2-informed-turn-v1-20261014/06-development-comparison/metadata.json)。
- [最终批次原始数据](../results/challenge2-informed-turn-v1-20261014/06-development-comparison/runs/)。
- [六面板指标图](../results/challenge2-informed-turn-v1-20261014/06-development-comparison/review-figures/occlusion-response-T-18.png)。
- [发布检查脚本](20261014-informed-turn-publication/verify_publication.py)。
- [发布检查记录](20261014-informed-turn-publication/publication_checks.json)。

结果包保留全部七个开发阶段，包括未通过设置、完整时序、长期个体诊断与运行日志。它采用固定3个外部提示对象、17个跟随者和k=15选边；不能当作无固定角色的20点全连接实验。

无结果复刻协议与结果目录相互独立。若要求复刻者事先不知结果，仅发送该Markdown及纯模型代码；本索引、实现README、结果文件和历史日志不属于盲复刻起始材料。

历史audit和README中的“未提交/未推送/main未修改”是实验完成时的记录，不被本次发布改写。当前发布状态以Git提交与PR为准。`.gitattributes`对本实验使用`-text`，使Git中的文件字节与原始SHA256审计一致。

## 随分支整合保留的既有研究材料

- [SQ1局部输入模型](../SQ1_LOCAL_INPUT_MODEL.md)。
- [既有审查、配对分析与时序分析索引](../../reviews/README.md)。

这些材料保留原有证据等级与历史路径。分支合并只整合版本历史，不改变实验参数或升级科研证据等级。

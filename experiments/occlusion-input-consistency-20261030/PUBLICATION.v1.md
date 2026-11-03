# 遮挡后初始社会输入与最终运动一致性：数据发布

本次运行 `run-v2-20261102` 已完成14000条轨迹：主实验500种子×20档遮挡，0步及5200步历史对照各100种子×20档。主实验遮挡前200步，遮挡后6000步。采用v2三因素S0与原始无补偿神经场。种子2026103001—2026103500，bootstrap种子2026103999、10000次按完整种子重采样。

## 阅读入口

- [中文报告](results/run-v2-20261102/delivery-v1/实验报告.v1.md)
- [主图及窗口对照](results/run-v2-20261102/delivery-v1/)
- [历史与总强度对照、逐轨迹汇总及复核](results/run-v2-20261102/analysis-v2/)
- [原始数据与配置](results/run-v2-20261102/)
- [运行、暂停和续跑日志及代码哈希](execution-run-v2-20261102/)
- [发布文件SHA256清单](publication_manifest.v1.json)

每个seed目录保存初态、历史检查点、各遮挡分支的完整位置/速度/逐步指标及summary.csv。MAT通过Git LFS发布，约38.9 GB。下载当前实验原始数据：

```powershell
git lfs install
git lfs pull origin YCwork --include="experiments/occlusion-input-consistency-20261030/results/run-v2-20261102/**/*.mat"
```

## 复现与验证

MATLAB执行入口为runOcclusionConsistency，使用新的run ID；resume只用于同配置的中断运行。v2/analyze_v2.py完成预定分析，deliver_results.py生成补充图及中文报告。执行和续跑脚本保留运行顺序。环境依赖见requirements.txt及execution.json。

此前静态检查、原模型model-smoke与续跑等价测试通过。9组原始轨迹独立复算GO/LO/C，最大C误差约1e-15；独立S0误差约4.86e-17。全部14000条记录完整且键唯一。图表提供PNG/PDF。

本次修复首条结构体记录赋值并增加已保存完整轨迹的续跑复用；模型方程、初态随机流、遮挡集合及种子保持冻结。代码来源提交06b1e07，实际执行版本及修复差异记录在运行快照和日志中。

证据类型：预设模拟设计下的探索性敏感区与平台分析。0.06—0.09为本条件下的候选敏感输入区域。原始模型、其他历史结果和论文保持原样。本次发布分支YCwork，main更新另行安排。公共回退采用本次发布提交的git revert，原始证据继续保留。

# GitHub 发布入口：KF 诊断 A 包

发布仓库：PatrickStar-cmd/Flocking_Study；目标分支：`Hao_L`。

本次发布包含已完成的 P0/P1 诊断代码、45条唯一轨迹的完整 MAT、配置、源码快照、逐种子统计、八张图和独立审核记录。结果仍为5种子开发集探索性证据；P2及后续阶段未执行。本发布不修改原始模型、历史结果、论文或 `main`。

## 阅读顺序

1. [完整诊断报告](KF_DIAGNOSTIC_A_REPORT.md)：结论、效应量、区间、证据限制与下一步。
2. [代码说明与复现命令](README.md)：MATLAB/Python入口、状态索引、运行参数与统计定义。
3. [完整结果目录](../../effective-social-input/results/kf-diagnostic-v1-20261101-a01-reviewed/)：`raw/`、`analysis/`、`sources/`、配置及运行记录。
4. [独立原始数据复核](review/independent-package-audit.json)及[估计器诊断复核](review/independent-diagnostics-audit.json)。

## 获取45个原始 MAT

原始数据共3,802,119,993字节，逐文件 SHA256 记录在结果目录的 `manifest.json` 和独立审核文件中。每个 MAT 保存实际配置、逐步指标、观测记录、物理与神经状态、检查点和传感噪声。

MAT 通过 Git LFS 发布；请安装 Git LFS 后克隆 `Hao_L`：

```powershell
git lfs install
git clone --branch Hao_L https://github.com/PatrickStar-cmd/Flocking_Study.git
cd Flocking_Study
git lfs pull --include="effective-social-input/results/kf-diagnostic-v1-20261101-a01-reviewed/raw/*.mat"
```

已克隆仓库时，先在适当工作区同步 `Hao_L`，再执行对应的 `git lfs pull`。完整 MAT 应约84 MB；若文件只有约130字节，里面包含 `version https://git-lfs.github.com/spec/v1`，则只是指针，尚未取回原始内容。GitHub网页源码压缩包是否包含LFS实物取决于仓库设置，复核时使用上述命令。

MATLAB使用 R2025a；Python使用3.11，分析依赖 NumPy、SciPy、pandas、Matplotlib，见 `requirements.txt`。运行新实验必须使用新的输出目录，已有原始结果保持只读。

## 发布范围与历史状态说明

- 代码及轻量结果使用普通 Git；45个主实验 MAT 使用局部 `raw/.gitattributes` 管理的 Git LFS。
- 软件调试时的重复 MAT、第一种子硬链接 fixture、缓存没有重复上传；它们不构成额外独立样本，相关软件日志与审核说明保留。
- 报告与审计中的“未提交、未推送”等文字记录的是实验完成时的状态。本次发布状态以当前 Git 提交和此入口为准，不改写实验时的来源记录。
- 你原先本地 `Hao_L` 的工作总结提交 `dfe5864` 保留在原工作区；本次只发布此实验，没有重写或删除该本地提交。
- 公共回退使用对应发布提交的 `git revert`，保留所有原始数据与审核记录。

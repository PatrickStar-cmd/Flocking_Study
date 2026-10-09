# SQ1 四分支 r3 软件验收包

目标：修复 r2 的 R01–R05，补齐协议1.2的真实 G1–G15、完整数据导出及来源记录，供杨骋（C）独立复核。主配置固定 seed1、N20/M100、k8/h4、A0=0.228、w0=0.0285、random/d34、无噪声、无KF。原 r2 和原论文代码保持不变。

最终有效结果目录及实测测试状态见根目录 `R3_STATUS.md`。同级 `_r3_build_history_20261009/` 内失败尝试和预检只作生成过程记录，不混入有效样本。

## C 的读取顺序

1. `R3_STATUS.md`、`FIXES.md`：本次变更、实际测试结果、与r2差异及范围。
2. 有效结果内 `metadata.json`、`audit.json`、`source_manifest.json`、`RUN_MANIFEST.csv`、`tests.json`。
3. `logs/G1.log`…`G15.log`、`resume_check.csv`、`T0_maxdiff.csv`、`weight_contract.csv`。
4. 四个模式MAT、`checkpoints.mat`、`T0_results.mat`、`reference_same_state_replay.mat`、`summary.csv`。
5. 根目录 `DATA_DICTIONARY.md`、`SHA256SUMS.csv`；可直接运行 `load_pilot_example.m`。

## 独立复现

MATLAB R2023a 是此次实际执行环境。进入本包后运行：

```matlab
restoredefaultpath;
addpath(pwd);
run_protocol12_seed1_r3('C_independent_recheck');
```

新名称会生成独立目录。已有名称会拒绝覆盖。运行器自动使用包内来源快照，源码不依赖原电脑绝对目录。完整运行结束会执行 G1–G15；查看 tests.json 的实际状态。要单独重跑门槛测试，可传入现有结果目录调用 run_tests，但这会更新该目录的测试日志，独立验收建议复制整包或创建新输出。

`sources/protocol` 是 C 固定协议的只读快照；`sources/A` 是 Hao_L@992e71a 的参考函数快照；`sources/legacy` 保留作者原核心和r2核心文件供版本比较。实际加载路径和哈希以 source_manifest 为准。没有把本地文件伪报为已提交Git版本；本地实现由逐文件SHA256固定。本任务没有提交、推送或更改main。

结果只能说明软件链路是否符合协议。seed1、mask、个体和时间步不是新增独立科学重复；未运行20seed/560轨迹、80条正式T0批次，也没有bootstrap或论文收益结论。C复核通过后再安排正式批次，A仍负责输入公式独立确认。

PowerShell中运行 `./hash_package.ps1 -Verify` 可核对交付清单的全部文件哈希和大小。SHA256SUMS.csv不包含自身；生成后修改文件会使校验失败。

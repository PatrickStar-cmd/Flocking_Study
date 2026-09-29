# 240 个规范运行与正式 MAT 的关联入口

日期：2026-09-29。任务第一项已完成。

## 结论

480 个正式 MAT 全部找到并通过 SHA256 核查，可映射到 240 个规范运行（12 条件 × seed 1–20）。额外 240 个标签副本与对应代表运行的全部已保存数据精确一致，支持沿用原去重口径进行后续时序分析。

这里的“240 个规范运行”不是 240 个彼此独立的统计重复。同一 seed 跨条件仍按配对设计处理；每个条件只有 20 个 seed。证据类型为软件/数据一致性核查，原实验仍属于探索性开发数据。

## 文件用途

| 文件 | 行数/用途 |
|---|---|
| `CANONICAL_RAW_INDEX.csv` | **后续分析的主入口，240 行**。每个规范运行恰好一个代表 MAT，含路径、SHA256、原始标签、seed、参数、字段及采样间隔 |
| `RAW_FILE_MAP.csv` | 480 行。保留每个标签文件及其规范键、代表文件、压缩包及包内路径、哈希、别名原因 |
| `alias_checks.json` | MATLAB 对全部文件的逐字段比较记录；包含 240 次非代表文件比较，另有 240 次代表文件自检 |
| `audit.json` | 输入、输出和核查脚本哈希，条件数和重复标签分布 |
| `matlab_alias_verification.log` | MATLAB R2024b 核查日志；没有运行任何仿真 |
| `build_raw_index.py`、`verify_aliases.m` | 可复现关联和等价核查的脚本 |

CSV 的 `raw_path` 和 `representative_raw_path` 均相对**仓库根目录**，不是相对本目录。绝对路径由仓库根目录与该字段拼接。MAT 中的数据变量名为 `rawResult`。源数据仍位于上轮审查目录，本轮不复制或修改它们。

代表文件继续沿用旧台账中最早的 CSV 行，不依据 GO 或任何结果选择。`canonical_policy` 是规范名，`selectionPolicy` 是实际原始配置标签；例如 C02 规范名为 balanced，而代表文件可保留 angle-dispersed 标签。C08 的规范名 all 表示 k=19 的全邻居等价条件，原始配置仍记录 random。

## 等价核查的范围与规则

- 固定源版本：`WHQ_work@0036487818d891c23d357f7af49030f323f4c9af`。
- 规范键沿用 `../sq1-c-condition-ledger-v1/` 的两个 CSV，未修改旧表。
- 对每份 MAT 比较 `timeSeries`、`trajectory`、`summary`、`finalPosition`、`finalHeadings`、`finalNeuralState`、`finalPreferredDirections`、`meanAdjacency`，同时核对顶层字段集合。
- 使用 MATLAB `isequaln` 精确数值比较，不使用宽松误差阈值掩盖差异；全部 14 项 timeSeries 均检查 8000 点且有限。
- 配置只有两项明确的等价转换：angle-dispersed → balanced；k=N−1 时 original/fixed-total → full-equivalent。其余配置字段必须精确一致。
- 全部 480 个原文件 SHA256 均匹配上轮来源清单。MAT 文件字节哈希可能不同，不能用文件哈希不同推断数值结果不同；反过来也不能只比较汇总均值便认定轨迹等价。
- 本次没有证明跨不同规范条件的随机驱动逐步一致，也没有恢复未保存的中间神经状态、观测缓存或 RNG 状态。

## 完整性与独立复核

| 检查 | 结果 |
|---|---|
| 原始标签/文件 | 480，无缺失、无重名 |
| 规范键/代表文件 | 240，每键恰好 1 个代表 |
| 条件覆盖 | C01–C12，各 20 个 seed，均为 1–20 |
| 别名比较 | 240/240 通过，全部已保存数据精确一致 |
| 配置比较 | 允许上述两种转换后，480/480 通过 |
| 每键标签数 | 100 个键有 1 个标签，80 个有 2 个，20 个有 3 个，40 个有 4 个 |
| 独立计数 | 用 PowerShell Import-Csv/Group-Object 独立得到相同条件数与标签分布：100+80+20+40=240；100+160+60+160=480 |

## 后续如何使用

后续分析只遍历 `CANONICAL_RAW_INDEX.csv`。按 `canonical_config_id` 和 `seed` 配对，读取 `rawResult.timeSeries`。需要追溯旧实验标签时查询 `RAW_FILE_MAP.csv`，不要把该表 480 行直接当统计样本。

时间序列每步保存：8000 点，dt=0.3，对应步骤 1–8000；原评价窗为 6001–8000。空间轨迹包含步骤 0，每 100 步一帧，共 81 帧，间隔 30 模型时间单位。不能将轨迹帧数与 timeSeries 长度混用。

下一项为按原 11 组对比制作时间过程与配对差值图。本轮没有提前运行该分析，也没有新增实验种子。

## 复现与保护

从仓库根目录运行 bundled Python 的 `build_raw_index.py`，再以 MATLAB 执行 `verify_aliases.m`，最后运行 `build_raw_index.py --finalize`。首次核查依赖前两轮审查目录中的来源文件与台账。

Python 对已存在输出仅接受字节完全一致的结果；MATLAB 拒绝覆盖已有 `alias_checks.json`。需要完整重做时复制脚本到新的同级审查目录，不删除旧证据。任一别名比较失败，最终 240 行入口表不会生成。

当前分支 `YCwork`，HEAD `71583f9`，尚未提交或推送。新增内容仅在本目录，原始代码、历史实验、旧台账和 `main` 未改动。已执行完整核查流程和独立计数，没有修改仿真器，因此不运行仿真 smoke test。建议作为 C 的分析入口保留；不需要合并其他成员的代码或自动合并 main。

回退仅需撤回本次新增目录，不涉及原始数据；若后续提交并共享，用 git revert 保留公共历史。

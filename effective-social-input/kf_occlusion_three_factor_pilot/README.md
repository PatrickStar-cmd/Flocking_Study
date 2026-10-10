# KF遮挡三因素开发实验

本文件夹是上传至 `WHQ_work` 的轻量实验包，包含可独立运行的MATLAB代码、来源快照、汇总CSV、PNG/PDF图表、日志，以及Python独立复核和作图脚本。输入是邻居方位，KF状态仅为`[方位角, 角速度]`。

**原始MAT及检查点未上传。** 270个MAT文件约5.87GB，保留在实验者电脑的 `10.10周/kf_occlusion_three_factor_pilot/`。其文件路径、大小和原始SHA256见 [RAW_DATA_MANIFEST.csv](RAW_DATA_MANIFEST.csv)。RUN_MANIFEST和历史审计记录中的MAT路径指向原始完整实验包，并不表示这些文件已在GitHub提供。

可先读 [RESULTS_说明.md](RESULTS_说明.md) 和 [三组扫描图](results/pilot_seed1to5_v1/analysis/figures/README_图表说明.md)。已上传结果对应seed1—5、7条件、155条闭环轨迹；中心条件KF平均输入误差低于hold，但seed2预测更差、协方差存在失配，仍属探索性结果。

## 先读哪些文件

1. `DESIGN.md`：本轮真正执行的实验条件与评价规则。
2. `RESULTS_说明.md`：批次完成后由实际结果生成的中文解释。
3. `results/pilot_seed1to5_v1/`：有效五种子开发pilot的汇总、图表及审计记录（不含MAT）。
4. `results/smoke_v2/`：早期seed1接通试验，不并入开发pilot统计。
5. `results/smoke_v1/`及`logs/initial_checks/`：保留的早期失败诊断，不计入科学样本。

根目录`tests/software_checks.*`与`software_smoke.log`也是最初调试记录，含已修正的测试索引失败；最终有效10项PASS请读`results/pilot_seed1to5_v1/logs/software_checks.json`。仿真过程中未剔除任何种子。Python复核结果见同一有效批次的`analysis/validation.json`。

## 独立复现

在MATLAB R2023a中，将当前文件夹切换到本实验目录，执行：

```matlab
restoredefaultpath;
addpath(pwd);
run_kf_occlusion_pilot('my_independent_pilot','pilot');
```

函数拒绝覆盖已有结果，请使用新名称。smoke只运行seed1的中心条件：

```matlab
run_kf_occlusion_pilot('my_smoke','smoke');
```

批次完成后可从命令行复核并生成图表：

```powershell
python analyze_results.py results/my_independent_pilot
```

以上分析命令需要MATLAB重新生成的原始MAT。仅下载GitHub轻量包，可校验已上传文件：

```powershell
python verify_uploaded_package.py
```

`verify_delivery.py`、`finalize_delivery.py`及原始DATA_DICTIONARY说明服务于本地完整包，不能直接用于缺少MAT的轻量包。历史 `analysis/validation.json` 是本地完整数据已执行复核的证据；轻量包哈希验证不等于重新复核原始轨迹。

Python依赖numpy、pandas、matplotlib、h5py。仿真只依赖本包MATLAB文件；无需给旧仓库addpath。

## 代码对应实验过程

- `run_kf_occlusion_pilot.m`：顺序调度，运行键、日志、参数、失败记录和实际轨迹登记。
- `simulateKFSegment.m`：神经场、Euler积分、运动、GO/LO、周期边界及完整逐步记录。
- `sensorObservation.m`：构造合法观测；隐藏方位与权重设为NaN。
- `controllerInputFromObservation.m`：五种补偿策略，控制器接口不接收真值或位置。
- `bearingTrackBankStep.m`：批量二维KF，圆周创新及Joseph协方差更新。
- `makeKFCheckpoint.m`、`checkpointToStart.m`：完整状态保存与续跑。
- `makeProtocolMasks.m`：固定隐藏集合；按实际ID枚举，并用独立随机流处理random规则。
- `run_kf_prediction_diagnostics.m`：fresh轨迹上的离线KF误差、覆盖率、R敏感性。
- `run_same_state_input_experiment.m`：同一fresh状态上的输入场比较。
- `run_closed_loop_experiment.m`、`summarizeClosedLoop.m`：6200步分叉及配对群体响应。
- `run_software_checks.m`：实际执行的10项接口和数值检查。
- `analyze_results.py`：读取原始MAT，独立复算、五种子配对汇总、图表与中文报告。

## 完整性与版本

运行前的MATLAB源码版本见结果内`source_manifest.json`，本次上传未修改仿真源码。`UPLOAD_SHA256SUMS.csv`记录此次上传文件的哈希；`verify_uploaded_package.py`可复核。来源快照是`sources/r3`及`sources/kf`，只用于版本对照。完整本地包的`SHA256SUMS.csv`未上传，原始MAT条目另存于`RAW_DATA_MANIFEST.csv`。

有效pilot采用主R=1e-4，另两个R都后台使用完整6200步相同观测历史。smoke_v2早期版本的替代R仅预热200步，不能与有效pilot混合。q沿用历史候选1e-6，本轮没有独立训练择参，也未依据闭环结果选参。本批是固定候选的探索性评估，不宣称完成独立KF准入或协方差校准。

本轻量包上传至GitHub的WHQ_work分支。本轮没有运行KF-uncertainty、真实几何遮挡或完整三因素交互组合。

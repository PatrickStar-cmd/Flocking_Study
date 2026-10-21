# 无遮挡方位观测实验：无噪声参考、直接观测、KF、IMM

源码交接日期：2026-10-21。目标分支：`Hao_L`；独立任务分支：`exp/no-occlusion-bearing-handoff-20261021`。包从 `Hao_L@992e71a` 新增，不改仓库已有实现、原论文代码和历史结果。

## 队友从这里开始

这套源码对应近期**无新增遮挡、无领航任务**的方位噪声实验，包含高斯／非高斯噪声、低噪声 σ=.01 rad／高噪声 σ=.10 rad，以及四种输入：

| 方法 | 实际进入神经环的方位 | 真值可否进入估计器 |
|---|---|---|
| `oracle` | 真实方位；无观测误差理想诊断组 | 不运行KF／IMM |
| `raw` | 带噪方位观测 | 否 |
| `selected-kf` | 同样带噪观测经冻结KF后的方位 | 否 |
| `imm` | 同样带噪观测经冻结IMM后的方位 | 否 |

这里“无噪声”指oracle诊断对照。其余三组使用对应噪声水平；不把R设成0来构造数值退化的KF。oracle用于诊断，不保证GO一定最大。

代码已包含全部运行依赖，**无需原作者电脑上的历史MAT、Python临时目录或绝对D盘路径**。原MATLAB动力学和滤波器采用正式实验的冻结版本，21个源码文件逐字节保留并有SHA256。新增可移植runner只负责配置、循环、配对、保存和计数。历史runner／分析器另存archive，用于追溯；它们可能需要原归档目录，不能代替下面的可移植入口。

## 环境与快速运行

MATLAB：本次实测R2025a，Base MATLAB即可运行本包入口，不需要Parallel Computing Toolbox。Python可选，用于独立复算和绘图：Python3.12实测通过，依赖见requirements.txt。正式历史分析使用过NumPy2.5.4／SciPy1.18.1／Matplotlib3.11.2。

先切到本文件所在目录，再启动MATLAB：

```matlab
restoredefaultpath;
addpath(fullfile(pwd,'matlab'));
testNoOcclusionHandoff(fullfile(pwd,'outputs','my-smoke-001'));
```

同样的命令行方式（在本目录执行）：

```powershell
matlab -singleCompThread -batch "restoredefaultpath; addpath(fullfile(pwd,'matlab')); testNoOcclusionHandoff(fullfile(pwd,'outputs','my-smoke-001'));"
```

每次使用**不存在的新输出目录**。已存在目录会拒绝覆盖；运行失败保留原目录，重试换编号。不要用addpath(genpath(...))把archive或仓库其他同名MATLAB实现加入路径。outputs目录在Git中忽略。

独立分析（同样在本目录执行）：

```powershell
python -m pip install -r requirements.txt
python tools/verify_sources.py
python analysis/analyze_handoff.py outputs/my-smoke-001
```

`verify_sources.py`只用Python标准库，检查21个冻结文件字节哈希。分析器从保存的周期位置差独立重建速度、GO和平均对距，核对summary；输出审计、MAT哈希、六对GO差值和均值图。样本量1的smoke不计算推断；n>1时每个噪声条件内使用seed配对t检验、Holm6、Bonferroni同时95%区间。它不复现旧分析的全部bootstrap／滤波器矩阵回放；完整旧分析源码见archive。

## 三个明确的运行配置

```matlab
% 仅短程软件检查；不是新的正式研究结果
runNoOcclusionHandoff(fullfile(pwd,'outputs','smoke-002'),'smoke');

% 重现旧6个种子：高斯/非高斯 × .01/.10 × 四组
runNoOcclusionHandoff(fullfile(pwd,'outputs','pilot-reproduction-001'),'pilot');

% 重现最新64种子：仅非高斯、sigma=.01、四组
runNoOcclusionHandoff(fullfile(pwd,'outputs','validation64-reproduction-001'),'validation64');
```

| profile | seed | 时间步／burn-in | 输出条件行 | 唯一轨迹 |
|---|---|---|---:|---:|
| smoke | 6299 | 80／40 | 16 | 13 |
| pilot | 6201–6206 | 8000／6000 | 96 | 78 |
| validation64 | 73001–73064 | 8000／6000 | 256 | 256 |

oracle对同一seed在噪声条件间复用，只存一份MAT。因此pilot中24行oracle对应6条唯一轨迹，smoke中4行oracle对应1条轨迹；不能把条件行当独立样本。每组summary的matFile指向实际文件。pilot／validation64均复用已做过的种子，**重跑是复现，不是新的独立验证集**。

本次只实际运行smoke及契约测试，没有为了交付重跑pilot／64seed全量实验。若要研究新问题，复制runner为新版本，先冻结问题、参数、种子和统计计划；不要改写现有冻结源文件。以上完整运行可能耗时数十分钟并产生较大MAT文件。

## 模型、噪声与配对契约

- N20、环形神经元M100、周期场地L1000、dt=.3、逐步随机选择k15、等权输入。无scheduledOcclusion、无几何遮挡、无turn cue、非openArena。
- 神经参数以`configs/base-config.json`为准：beta1000、globalInhibition0、speedScale.05、totalSocialAttraction.24、receptiveFieldWidth.4。沿用原神经场动力学，没有每步另加神经噪声。
- KF状态为`[视线方位theta,视线角速度omega]`，F=`[1,dt;0,1]`，H=`[1,0]`，Q=`q*[dt^3/3,dt^2/2;dt^2/2,dt]`，R=σ²；wrapped创新、Joseph协方差更新。冻结q=.1。
- IMM是相同F、不同Q的两个常角速度模型，q低1e-4／高.1；每模型步转移矩阵`[.98,.02;.10,.90]`。不是已经完成CV/CA/CT全集，也不是自动抗离群滤波器。
- mixture噪声：95%小高斯+5%标准差8倍的大高斯，除以sqrt(4.15)后乘目标σ；保证总体标准差与同级高斯一致。理论峰度35.8397。
- 四方法共享初态、每步选邻和噪声扰动；raw/KF/IMM同时核对噪声数组、随机流。轨迹一旦分化，真实方位和最后观测会不同，因此闭环不是“完全相同观测序列”比较。纯估计器比较使用`buildBearingReplaySequence.m`、`replayBearingModels.m`及archive中的IMM replay设计。
- 真值用于仿真几何、评估和oracle组；KF/IMM纯controller的接口只接收观测。软件测试改变附加真值／位置／离群标签后，controller输出不变。

## 已完成实验与证据范围

本包包含以下历史源码、参数、轻量摘要和溯源信息，大型原始MAT／NPZ不随源码提交：

| 历史系列 | 条件／规模 | 位置 |
|---|---|---|
| 初期高斯raw/KF | 5101–5106，σ.01/.10，q=1e-4，24运行 | archive/matlab/runNoisyBearingKFPilot.m；provenance/noisy-bearing-kf-pilot-20261014-v2 |
| 初期非高斯raw/KF | 同5101–5106与q=1e-4，24运行 | archive/matlab/runNonGaussianBearingKFPilot.m；对应provenance目录 |
| KF大RMSE诊断 | 复用低噪声轨迹，检查极端误差／周期几何跳变 | archive/analysis/diagnose_noisy_kf_rmse_20261014.m等 |
| IMM同观测回放 | 开发5101–5106；新6101–6106，24参考序列 | matlab/replayBearingModels.m；archive/matlab/runBearingIMMReplayPilot.m；provenance内frozen-selection.json |
| 三组闭环 | 6201–6206，2分布×2σ×3法，72轨迹 | archive/matlab/runBearingIMMClosedLoopPilot.m；对应provenance |
| 六种子oracle扩展 | 新6 oracle + 旧72轨迹，96行／78唯一轨迹 | archive/matlab/runOracleBearingDiagnostic.m；对应provenance |
| 最新64seed四组 | 73001–73064，mixture σ=.01，256轨迹 | archive/matlab/prepareNonGaussianFourGroupValidation.m、runNonGaussianFourGroupValidation.m；reports/validation64 |

**早期KF q=1e-4，当前冻结KF q=.1，两者不是同一算法配置。**便携pilot使用后者，重现的是后期四组条件；若重现初期24运行，按archive原runner和其metadata执行，并先修正输入／输出路径。归档脚本保留原貌，未保证所有旧路径在新机器存在。

最新64seed平均GO：oracle .6031、raw .6430、KF .6378、IMM .6293。六对比较同时区间均跨0、Holm校正p均1.0；旧6seed中oracle−raw +.2979的大差距没有复现。IMM自身闭环方位RMSE较raw低约11.7%，但未转成GO收益。不能据此声称普遍无效、已经等效或已证实有害，也不能否定其他有遮挡／领航任务的条件性收益。详细数字和限制见`reports/validation64/结论与研究决策.md`、`64新种子独立验证报告.md`及CSV／JSON。

reports中的历史报告按归档原文保留，其旧本地绝对链接／未携带图像可能无法在GitHub打开；报告正文、统计CSV/JSON和全部必需运行源码在包内。`provenance/validation64`的旧delivery-audit描述原本完整本机结果，不代表本次上传了其提到的原始MAT和全部图。

## 文件职责

- `matlab/`：可直接使用的冻结模型、KF、IMM、oracle、回放函数和测试；新增便携runNoOcclusionHandoff/testNoOcclusionHandoff。
- `configs/base-config.json`：最新冻结基本配置；噪声分布、σ、seed由profile逐条件设置。
- `analysis/analyze_handoff.py`：新输出结构的独立几何复算、六对比较和绘图。
- `archive/`：原实验runner与完整旧分析源码，可能依赖旧目录结构／数据。
- `provenance/`：冻结哈希、原metadata／轻量summary／方案、历史审计与本次交付测试记录。
- `reports/validation64/`：最新正式结果的可读结论与轻量统计。
- `tools/assemble_local_sources.py`：一次性本机来源导出脚本，含源机器路径，队友运行本包不需要执行。
- `outputs/`：运行时生成，Git忽略；与原仓库results完全分开。

## 已通过的交付测试与限制

本次在独立克隆、`restoredefaultpath`后实测：零噪声旧模型精确回归、KF真值隔离／wrapped创新／协方差、非高斯p=0退化与混合分布、IMM同模型退化／旋转等变／缺测转移／半正定、闭环配对、oracle=零噪raw全部通过。smoke13个MAT／16条件行独立重算GO、对距、速度，最大误差2.22e-14。21个冻结源码SHA一致。记录见`provenance/delivery-tests/`。

这些是**软件交付证据**。历史科学结果沿用各自证据等级；没有新增高／低噪声正式样本。预测状态只有方位与角速度，不能声称恢复二维邻居位置；平均GO不是持续稳定、安全或凝聚的保证。

## 队友接续建议

先跑smoke和独立复算，再决定重现pilot还是64seed。先区分“估计更准”和“群体更有序”，保留不利种子与oracle诊断；不要为获得正结果反复更换滤波器或筛选seed。本包默认保持无遮挡，若接续加入遮挡，另开实验版本并保留hold基线。

Git回退采用`git revert`撤销本次源码交接提交；不删除本地历史数据、不改main。该提交只新增本目录，队友可选择性cherry-pick或沿Hao_L继续工作。

# P2：固定轨迹开环方位预测诊断

本目录承接 P1 A 包，比较 Hold、KF-CV、KF-GM 对同一观测序列的预测。运行前设计见 [PROTOCOL.md](PROTOCOL.md)。这是复用 P1 五个种子的探索性开发评估，不是独立验证。估计不反馈轨迹，不能据此声称改善闭环 GO。

结果入口：`effective-social-input/results/kf-p2-v1-20261102-reviewed/`。S/V 是同一无外部转向轨迹中按真值变化率选出的慢/快片段；T 为匹配初态、随机驱动的受控旋转刺激。参考轨迹是带测量噪声的 Fresh 闭环，只有后续估计器回放是开环。物理轨迹没有缺失，人工遮挡只作用于回放观测。

## 获取原始数据

```powershell
git lfs install
git clone --branch Hao_L https://github.com/PatrickStar-cmd/Flocking_Study.git
cd Flocking_Study
git lfs pull --include="effective-social-input/results/kf-p2-v1-20261102-reviewed/raw/**"
```

完整 P2 输入 MAT 保留每步20个体位置/航向、140条被遮挡有向边的真值和带噪观测，以及GO/速度。NPZ保留每个条件的估计序列、方差、NIS、年龄与评价窗口，可重建全部统计和图；大文件使用 Git LFS。

## 从头复现

MATLAB R2025a；Python3.11。Python依赖与P1一致，见 `requirements.txt`。在仓库根目录运行。已有输出拒绝覆盖，每次使用新的suffix，不向旧结果目录写入。

1. 取回5条上游P1 Fresh MAT（其余40条不用于生成P2）：

```powershell
git lfs pull --include="effective-social-input/results/kf-diagnostic-v1-20261101-a01-reviewed/raw/*_mode_fresh_D0.mat"
```

2. MATLAB命令窗口：

```matlab
addpath('experiments/kf-p2-20261102');
out=generateP2Inputs('replicate1');
matlabParityReference(out);
```

3. Python命令（将结果路径替换为刚生成的目录）：

```powershell
py -3.11 experiments/kf-p2-20261102/test_replay.py
py -3.11 experiments/kf-p2-20261102/audit_inputs.py effective-social-input/results/kf-p2-v1-20261102-replicate1
py -3.11 experiments/kf-p2-20261102/review_parity.py effective-social-input/results/kf-p2-v1-20261102-replicate1
py -3.11 experiments/kf-p2-20261102/replay.py effective-social-input/results/kf-p2-v1-20261102-replicate1
py -3.11 experiments/kf-p2-20261102/analyze.py effective-social-input/results/kf-p2-v1-20261102-replicate1
py -3.11 experiments/kf-p2-20261102/audit_geometry.py effective-social-input/results/kf-p2-v1-20261102-replicate1
py -3.11 experiments/kf-p2-20261102/review_stats.py effective-social-input/results/kf-p2-v1-20261102-replicate1
py -3.11 experiments/kf-p2-20261102/supplement.py effective-social-input/results/kf-p2-v1-20261102-replicate1
```

`audit_inputs.py`独立核对位置→方位、噪声配对、片段选择和SHA；`review_parity.py`对照原版MATLAB滤波器。最终统计独立复核命令和实测结论见随结果发布的报告。MAT文件压缩封装可能受MATLAB版本/保存头时间影响；复现比较数值内容，不要求新MAT字节SHA等于历史文件。

## 解释约束

- 主响应是无门限预测的遮挡期圆周RMSE；门限有效子集单独报告，避免选择偏差。
- Hold使用最后一次原始带噪观测，协方差、覆盖率、NIS均不定义。
- CV与GM使用各自仓库默认q；这是两个默认配置的比较，不是单独操纵模型结构的因果实验。
- P2不做Q扫描，不优化参数，不进入P3、P4或新方法实验。后续校准必须另立版本，保留本轮负结果。
- 各层、各D均只有5个独立seed；边、帧和S/V复用片段不增加独立样本量。
- 上传目标为Hao_L；不改main。公共回退使用发布提交的git revert。

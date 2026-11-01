"""Create the review report from verified artifacts; never edits raw results."""
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[3]
EXP=Path(__file__).resolve().parents[1]
RUN=ROOT/'effective-social-input/results/kf-diagnostic-v1-20261101-a01-reviewed'
ANALYSIS=RUN/'analysis'


def main():
    audit=json.loads((EXP/'review/independent-package-audit.json').read_text(encoding='utf-8'))
    assert audit['status']=='PASS'
    groups=pd.read_csv(ANALYSIS/'contrasts_group_summary.csv')
    summary=pd.read_csv(ANALYSIS/'summary.csv')
    diag=pd.read_csv(ANALYSIS/'estimator_diagnostics.csv')
    manifest=json.loads((RUN/'manifest.json').read_text(encoding='utf-8'))
    p0=json.loads((EXP/'p0-reviewed-p1-v2.json').read_text(encoding='utf-8'))
    assert manifest['completed_count']==45 and manifest['failed_count']==0
    def select(duration,window,baseline,method,metric):
        row=groups[(groups.D==duration)&(groups.window==window)&(groups.baseline==baseline)&
                   (groups.method==method)&(groups.metric==metric)]
        assert len(row)==1
        return row.iloc[0]
    def interval(row):
        return f"{row['mean']:+.8f} [{row.ci95_t_lower:+.8f}, {row.ci95_t_upper:+.8f}]"
    rows=[]
    for duration in (20,120):
        for baseline,method,label in [('hold','kf-mean','KFmean−Hold'),('hold','kalman','KFcov−Hold'),('kf-mean','kalman','KFcov−KFmean')]:
            pair=select(duration,'occ',baseline,method,'benefitPair')
            go=select(duration,'occ',baseline,method,'benefitGO')
            rows.append(f"| {duration} | {label} | {interval(pair)} | {interval(go)} |")
    fresh=summary[(summary['mode']=='fresh')&(summary.window=='occ')].groupby('eval_D_steps').GO_mean.mean()
    drop=select(120,'occ','fresh','none','benefitGO')
    drop_pair=select(120,'occ','fresh','none','benefitPair')
    short=select(20,'occ','hold','kf-mean','benefitPair')
    long=select(120,'occ','hold','kf-mean','benefitPair')
    nis=diag[diag['mode']=='kf-mean'].groupby('D').fresh_nis_mean.mean()
    seed_rms=[]
    for (duration,mode,seed),group in diag.groupby(['D','mode','seed']):
        valid=group.rmse_candidate.notna()
        if valid.any():
            seed_rms.append({'D':duration,'mode':mode,'seed':seed,
                'rmse':float(np.sqrt(np.average(group.loc[valid,'rmse_candidate']**2,weights=group.loc[valid,'n_hidden'])))})
    rms=pd.DataFrame(seed_rms).groupby(['D','mode']).rmse.mean()
    max_geometry=max(max(record['geometry_errors'].values()) for record in audit['raw_reports'])
    max_bump=max(record['bump_max_error'] for record in audit['raw_reports'])
    report=f'''# KF 首轮诊断 A 包：完成情况与审核结论

## Material Passport

- 项目：PatrickStar-cmd/Flocking_Study；日期：2026-11-01。
- 代码基线：Hao_L@f0b797220db5d6f2b6ed883e3838614e67cebbe5。
- 工作分支：exp/haoli-kf-diagnostic-20261101；目标分支：Hao_L。
- 证据：软件测试 + 五种子开发集探索性诊断；不是确认性、稳态边界或安全保证。
- Verification Status：VERIFIED（确定性重复、原模型精确回归和本批原始数据独立复算）；不代表科学结论已独立外推验证。
- 协作：AGY 完成主体实现及分析/绘图重构，保持其默认模型和推理设置；Codex 主助手审查、纠错、执行和独立统计复核。

## 判断

本轮已完成计划的 P0 + P1 A 包：45条唯一闭环轨迹，0条失败，8张规定图，完整原始MAT、参数、噪声、神经状态、源码快照与哈希均保存。当前应进入 P2 同轨迹开环诊断，不应直接增加门控新方法或扩大多因素扫描。

长遮挡在本条件下确有可测的 GO 损失，Fresh 的平均 GO 仅为 {fresh.loc[20]:.4f}（D20评价窗口）和 {fresh.loc[120]:.4f}（D120窗口），因此不能把本轮无明显 KF GO 收益解释成“Fresh 已接近1的天花板”。D120 的 Drop 相对 Fresh 的平均 GO 损失为 {-drop['mean']:.6f}，逐种子95% t区间 [{-drop.ci95_t_upper:.6f}, {-drop.ci95_t_lower:.6f}]；平均对距增量为 {-drop_pair['mean']:.6f}L。

KFmean 相对 Hold 的遮挡期对距收益在D20为 {short['mean']:.8f}L，折算L=1000仅约 {1000*short['mean']:.6f} 个距离单位；这是很小的数值变化。本轮没有预先冻结最小实质收益阈值，不据此主张工程价值。D120平均收益为 {long['mean']:.8f}L，但区间跨0；两档D的GO增益区间也均跨0，尚不能声称KF稳定提高GO。不能把“某一档区间跨0、另一档不跨0”当作已检验的时长交互。

## 预定比较结果

以下均为遮挡期、n=5个配对seed；括号为逐项95% t区间，未作多重比较调整。对距收益定义为 baseline−method，GO收益定义为 method−baseline，正值均表示候选方法较好。主比较为 Hold 与 KFmean 的对距，其他为机制或次响应。

| D（步） | 比较 | 对距收益 / L：均值 [区间] | GO收益：均值 [区间] |
|---:|---|---:|---:|
{chr(10).join(rows)}

全部逐seed差值与恢复期、公共窗口结果见 [contrasts.csv](../../effective-social-input/results/kf-diagnostic-v1-20261101-a01-reviewed/analysis/contrasts.csv)，没有删去负结果。10000次bootstrap按整seed重采样，全部D/方法共用同一索引矩阵；小样本区间仅供开发诊断。

## 门限与预测层线索

方差门限的影响已由独立Riccati递推、stock原模型和新诊断记录三方核对：第101个缺测步开始 P11>1；D120 的最后20步关闭全部7条隐藏输入，共2800/16800=16.7%隐藏边时间被截断，全部5个seed的KFmean/KFcov一致。maxAge=120的单位是模型时间；120步只有36模型时间单位，年龄门限本轮没有触发。它不是新的滤波器缺陷，而是必须报告的实际策略变化。

按“先每seed聚合整个遮挡窗口、再对5个seed取均值”计算，候选方位RMSE在D20为 Hold {rms.loc[(20,'hold')]:.4f} rad、KFmean {rms.loc[(20,'kf-mean')]:.4f} rad；D120为 Hold {rms.loc[(120,'hold')]:.4f} rad、KFmean {rms.loc[(120,'kf-mean')]:.4f} rad。短期MAE与RMSE可能给出不同排序，不能只挑一个误差指标。新鲜可见边的KFmean平均NIS在D20约{nis.loc[20]:.2f}、D120约{nis.loc[120]:.2f}，提示创新分布/运动模型值得检查。

这些误差均针对各方法自己的闭环物理轨迹，不能替代同一真值与观测序列上的估计器比较，也不能证明预测精度到GO之间的中介机制。候选预测与实际注入分开保存，避免只统计未被门限关闭的容易样本。

## 软件与独立审核

- P0八项在短程和完整2200步实际通过，完整检查耗时{p0['seconds']:.1f}秒；物理、神经、滤波器、传感噪声与RNG的共同前缀均逐字段核对。
- 仪器副本保持原动力学与随机调用；四模式短程和完整原模型轨迹精确回归通过。
- 45个MAT及源码快照SHA全部核对；从位置独立重算GO/LO/周期对距/最近邻/速度，最大绝对误差{max_geometry:.3g}；从u独立重算神经峰集中度，最大误差{max_bump:.3g}。
- 200行窗口汇总、210行配对对照、126行分组区间和2800行估计器诊断完整；窗口均值与原MAT最大误差{audit['max_csv_error']:.3g}。配对符号与t区间全部独立复算。
- 首版存在主响应窗口错用恢复期、Hold对另一分支真值求误差、只看P对角就声称PSD、检查点提前一步、失败仍标完成等问题。主助手已修正并重新验证，首版记录和失败日志保留于review/。

独立证据：[independent-package-audit.json](review/independent-package-audit.json)、[P0完整审计](p0-reviewed-p1-v2.json)、[修正清单](review/REQUIRED_FIXES.md)。最早一次stock参考在轨迹保存后因MATLAB不支持fopen x模式而无法写元数据；已修正导出并独立核验保存的轨迹，保留原错误日志，不把该导出失败算成45条主批失败。

## 限制、未完成与下一步

本轮只有5个开发seed、N20、k15、index-fixed隐藏7条边、默认q与R、两个D和固定事件起点。预遮挡群体仍自然演化，未按高GO筛起点，也未声称已处稳态。126个探索性区间没有家族校正；均值与恢复过程高度依赖seed，不能外推为普适邻居阈值或稳定性保证。

P2、P3、正式P4及P5没有运行。下一步唯一优先项是P2：冻结同一物理轨迹、观测和缺测掩码，让Hold与KF面对相同输入，按年龄检验圆周MAE/RMSE、NIS与P；先解释短期RMSE代价及P截断，再决定是否进入神经场离线回放。本轮没有新方法收益、吸引域证明、二维位置恢复或安全结论。

统计谬误已检查11/11，详见 [statistical-review.json](review/statistical-review.json)。主要保留意见是小样本、探索性多重比较及闭环误差不能充当同轨迹因果比较。

## 交付与版本保护

原始结果目录：[kf-diagnostic-v1-20261101-a01-reviewed](../../effective-social-input/results/kf-diagnostic-v1-20261101-a01-reviewed/)。45个完整MAT在raw/；配置、来源、快照、SHA和日志均在同目录；图表与CSV在analysis/。主批总耗时约{manifest['total_elapsed_seconds']/60:.2f}分钟。

仅新增本诊断代码目录和新结果目录；已有362个tracked文件保持SHA一致，allocentric-flocking-original、历史results和论文未修改。分支HEAD仍为f0b7972，当前新增内容未提交、未推送、未合并；main未更改。建议先将经审查的源码及轻量结果整理到Hao_L，完整MAT单独保留/传输，不直接合入main。回退为停用新增入口；未来如提交后需撤回，用源码提交的git revert，保留全部实验结果。

## 规定八图

1. [GO](../../effective-social-input/results/kf-diagnostic-v1-20261101-a01-reviewed/analysis/P1_Fresh_Drop_Hold_KFmean_KFcov_GO.png)
2. [相对Fresh的对距](../../effective-social-input/results/kf-diagnostic-v1-20261101-a01-reviewed/analysis/P1_PairDistance.png)
3. [LO](../../effective-social-input/results/kf-diagnostic-v1-20261101-a01-reviewed/analysis/P1_LocalOrder.png)
4. [新鲜数与实际激活数](../../effective-social-input/results/kf-diagnostic-v1-20261101-a01-reviewed/analysis/P1_ActiveNeighbors.png)
5. [输入积分](../../effective-social-input/results/kf-diagnostic-v1-20261101-a01-reviewed/analysis/P1_InputMass.png)
6. [神经峰集中度](../../effective-social-input/results/kf-diagnostic-v1-20261101-a01-reviewed/analysis/P1_BumpConcentration.png)
7. [逐seed配对收益](../../effective-social-input/results/kf-diagnostic-v1-20261101-a01-reviewed/analysis/P1_GainScatter.png)
8. [预先固定的代表轨迹与神经态](../../effective-social-input/results/kf-diagnostic-v1-20261101-a01-reviewed/analysis/P1_FailureCases.png)
'''
    target=EXP/'KF_DIAGNOSTIC_A_REPORT.md'
    with target.open('x',encoding='utf-8') as stream:stream.write(report)
    (RUN/'README.md').write_text('# KF诊断 A包（探索性开发结果）\n\n45条唯一轨迹，5个seed，0失败。\n\n完整报告：[KF_DIAGNOSTIC_A_REPORT.md](../../../experiments/kf-diagnostic-20261101/KF_DIAGNOSTIC_A_REPORT.md)。\n\nraw/为完整MAT，analysis/为汇总、逐seed对照、估计器诊断和八图；sources/、audit.json、manifest.json和source_hashes.json为来源证据。P2及后续阶段未运行。\n',encoding='utf-8')
    pd.DataFrame(seed_rms).to_csv(EXP/'review/own-closed-loop-seed-rmse.csv',index=False)
    print(str(target))


if __name__=='__main__':main()

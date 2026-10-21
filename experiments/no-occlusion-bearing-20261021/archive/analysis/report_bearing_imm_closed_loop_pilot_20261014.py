"""Report fixed-parameter closed-loop results without tuning or seed exclusion."""
import argparse
import csv
import json
from pathlib import Path


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--run',type=Path,required=True)
    root=parser.parse_args().run.resolve();out=root/'independent-analysis'
    path=out/'IMM闭环集群测试报告.md';assert not path.exists()
    audit=json.loads((out/'audit.json').read_text(encoding='utf-8'))
    assert audit['status']=='PASS_INDEPENDENT_CLOSED_LOOP_RECONSTRUCTION'
    means=json.loads((out/'condition-means.json').read_text(encoding='utf-8'))
    stats=json.loads((out/'paired-statistics.json').read_text(encoding='utf-8'))
    frozen=json.loads((root/'upstream-frozen-selection.json').read_text(encoding='utf-8'))
    with (out/'independently-recomputed.csv').open(encoding='utf-8-sig',newline='') as f:seed_rows=list(csv.DictReader(f))
    keyed={(r['distribution'],r['measurementStd'],r['method']):r for r in means}
    conditions=[(d,s) for d in ['gaussian','mixture'] for s in [.01,.1]]
    names={'gaussian':'高斯','mixture':'非高斯混合'}
    positive=[r for r in stats['primary'] if r['CI95Low']>0 and r['pHolmFour']<.05]
    negative=[r for r in stats['primary'] if r['CI95High']<0 and r['pHolmFour']<.05]
    text=['# IMM闭环集群测试报告','',
      'Material Passport：academic-research-suite / experiment-agent；2026-10-14。Verification Status：VERIFIED（完整轨迹与指标独立重建、预定8条滤波器重放）。证据等级：开发选参后的新种子小规模闭环检验；非跨参数确认性研究。','',
      f'预设四个条件中，IMM相对于调参单Q KF有 {len(positive)}/4 项GO改善、{len(negative)}/4 项GO下降同时满足95%差值区间不含0及Holm校正p<.05。其他条件证据不足以确认方向，不把未显著写成等效。','',
      '## 本轮实际运行','',
      '新种子6201–6206，每种子Gaussian/同方差非Gaussian混合 × std .01/.10 rad × raw/selected-kf/imm。72闭环运行、24配对三元组完整纳入；没有筛除不利种子，也没有根据结果调整模型、参数、评价窗口或终点。','',
      'N20、M100、随机逐步k15、等权original h_t=.24、kernel .4、beta1000、globalInhibition0、speedScale .05、dt .3、周期L1000，8000步，评价6001:8000。无新增领航、遮挡、距离感知、神经场规则、协方差展宽或降权。每一步只将当前被选邻居的带噪观测或后验方位均值注入神经场。','',
      f'单Q KF固定q={frozen["bestSingleQ"]:g}；IMM低/高q={frozen["lowQ"]:g}/{frozen["immHighQ"]:g}、两状态[theta,omega]、每模型步转移[[.98,.02],[.10,.90]]，R=std²。参数来自上一阶段只按同观测方位RMSE选定的frozen-selection.json；不是按本轮GO选参，也不声称单KF已针对GO全局优化。','',
      '混合噪声95%小高斯/5%大高斯，大/小标准差比8、除以sqrt(4.15)匹配总体方差。Gaussian stream seed+15000、混合标签stream seed+25000。三个方法共享初始位置/航向/神经状态、全部选邻、噪声扰动与标签；实际观测角度在闭环几何分叉后不同。真值、位置、距离及标签只供传感器/评价，不提供给控制器。','',
      '## 实际GO与预设主要比较','',
      '下面为每条件6种子的平均GO；它不是持续有序比例或稳定性证明。主比较IMM−开发选定单Q KF，正值表示IMM更有序。paired t 95%区间为普通单项区间，p对四条件Holm校正，独立单位为种子。','',
      '| 噪声 | std | 直接观测GO | 调参KF GO | IMM GO | IMM−KF | 95%区间 | Holm p | 上升种子 |',
      '|---|---:|---:|---:|---:|---:|---|---:|---:|']
    for d,s in conditions:
        a,b,c=[keyed[d,s,m] for m in ['raw','selected-kf','imm']]
        r=next(r for r in stats['primary'] if r['distribution']==d and r['measurementStd']==s)
        text.append(f'| {names[d]} | {s:.2f} | {a["meanGO"]:.6f} | {b["meanGO"]:.6f} | {c["meanGO"]:.6f} | {r["difference"]:+.6f} | [{r["CI95Low"]:+.6f}, {r["CI95High"]:+.6f}] | {r["pHolmFour"]:.4g} | {r["positiveSeeds"]}/6 |')
    text+=['','## LO、凝聚与估计误差','',
      '这些是次要描述性指标；次要配对区间和未校正p保存在paired-statistics.json。对距更小只表示当前指标下更紧密，不作为安全或避碰证据。估计误差和真值输入误差都在各自闭环几何中评价，不是同观测序列上的纯降噪差值。','',
      '| 条件 | 方法 | LO | 平均对距/L | 方位RMSE(rad) | 输入对真值RMS |',
      '|---|---|---:|---:|---:|---:|']
    for d,s in conditions:
        for m,label in [('raw','直接观测'),('selected-kf','调参KF'),('imm','IMM')]:
            r=keyed[d,s,m]
            text.append(f'| {names[d]} {s:.2f} | {label} | {r["meanLO"]:.6f} | {r["meanPairDistanceNorm"]:.6f} | {r["usedBearingRMSE"]:.6f} | {r["inputTruthRMS"]:.6g} |')
    by_seed={(int(r['seed']),r['distribution'],float(r['measurementStd']),r['method']):r for r in seed_rows}
    opposite=[]
    for seed in range(6201,6207):
        for d,s in conditions:
            k=by_seed[seed,d,s,'selected-kf'];i=by_seed[seed,d,s,'imm']
            if float(i['usedBearingRMSE'])<float(k['usedBearingRMSE']) and float(i['meanGO'])<float(k['meanGO']):opposite.append((seed,d,s,k,i))
    text+=['',f'24个种子×条件配对中，有 {len(opposite)} 个同时出现“IMM方位RMSE更小、GO更低”。下面按种子/条件顺序列出前3个，仅用于说明指标排序不一致，不作典型轨迹或中介因果证明：','']
    for seed,d,s,k,i in opposite[:3]:
        text.append(f'- seed{seed}、{names[d]} std={s:.2f}：KF→IMM方位RMSE {float(k["usedBearingRMSE"]):.6f}→{float(i["usedBearingRMSE"]):.6f}，GO {float(k["meanGO"]):.6f}→{float(i["meanGO"]):.6f}。')
    text+=['','## 和前一阶段的关系','',
      '上一阶段在完全相同的观测序列上，用未参与选择的新种子6101–6106验证：IMM相对开发选定单Q KF的RMSE均值，高斯std .01/.10降低约24%/25%，非高斯std .01/.10降低约9%/17%，四个主要差值经Holm校正支持误差下降。它验证估计能力；本轮另用6201–6206检验闭环反馈效果。','',
      '估计更准与群体任务改善是不同证据层。闭环轨迹分叉后，运动几何又会改变角度变化率、周期图像切换、输入结构与后续误差。本轮没有单独操纵中介变量，因此不能仅凭RMSE与GO的同时变化证明某一种机制。','',
      '## 测试与完整独立复核','',
      'Gaussian/非Gaussian raw轨迹与此前实现逐步完全一致；共享初态、选邻、噪声与标签；控制器对附加真值、位置和标签的改变完全不敏感；相同模型退化的输入场最大差0。原IMM实现保持原样。','',
      f'全部72条轨迹独立重算GO、LO、归一化对距与速度，并从实际已用方位、观测和位置真值重算各自完整评价窗口600000个角误差点。全部72×2000步输入场也由记录的已用方位和真值重新构建，核对实际输入RMS；归档输入样本再次重建。最大指标差：{json.dumps(audit["maxErrors"],ensure_ascii=False)}。','',
      '72唯一键、24三元组、全部6新种子、每步300个被选边及配对随机流/扰动/标签完整通过。seed6201的四条件下KF与IMM各完整8000步用独立矩阵实现重放，共8条；未声称其余种子的滤波状态均第二实现重放。全部轨迹和统计均复核。','',
      '43,200,000个评价边×步误差并非43,200,000个独立样本；每条件的推断仍为n=6。p/区间依赖种子差值近似正态，样本很小时不能泛化为跨场景保证。流式复核只读取保存完毕且已发布三条CSV记录的MAT；最终统计在72条运行完成后生成。','',
      '## 边界与后续判断','',
      '结果限于当前N/k、神经参数、周期空间、两类观测分布与固定参数，没有证明持续稳定、安全、所有非高斯噪声鲁棒或二维位置恢复。只有方位/角速度估计；普通IMM仍可能追随大观测误差，角度融合使用局部Gaussian近似。','',
      '是否替换主算法应以预设闭环比较及相应适用条件判断。若GO区间跨0，继续保留IMM作为已检验的估计器候选；不能把上一阶段的RMSE收益写成flocking收益。下一研究步骤宜针对“估计—实际社会输入—群体运动”关系开展同状态的受控机制干预，另冻结方案与新种子，而非根据本轮有利种子反复调整算法。','',
      '## 文件、保护与回退','',
      '分支exp/bearing-imm-closed-loop-pilot-v1，HEAD f8b73cf，无commit/push/merge。新增bearingIMMClosedLoopControllerStep.m、bearingIMMClosedLoopInputStep.m、simulateBearingIMMClosedLoopPilot.m、testBearingIMMClosedLoop.m、runBearingIMMClosedLoopPilot.m；handoff中新增计划、独立复核和本报告生成脚本。','',
      '原始模型、旧IMM实现、此前所有结果、论文及四个既有tracked改动未修改。当前结果目录保存17项源码快照/哈希、上游冻结配置、分析和独立矩阵实现快照、完整MAT三元组、日志、软件/独立/交付审计和图表。','',
      '建议保留实验分支，代码审查后再合并研究实现；科学结论按上述范围保留。回退只撤销本阶段新增实现，历史及本次结果继续保留。','']
    path.write_text('\n'.join(text),encoding='utf-8')
    print(json.dumps({'report':str(path),'primaryGOImprovementConditions':len(positive),
                      'primaryGODeteriorationConditions':len(negative),'lowerRMSELowerGOPairs':len(opposite)},ensure_ascii=False,indent=2))


if __name__=='__main__':main()

"""Render a factual report from already independently checked artifacts."""
import argparse
import csv
import json
from pathlib import Path


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--run',type=Path,required=True)
    root=parser.parse_args().run.resolve()
    analysis=root/'independent-analysis'
    path=analysis/'IMM同观测序列测试报告.md'
    assert not path.exists()
    audit=json.loads((analysis/'audit.json').read_text(encoding='utf-8'))
    assert audit['status']=='PASS_INDEPENDENT_MATRIX_REPLAY_AND_STATISTICS'
    frozen=json.loads((root/'frozen-selection.json').read_text(encoding='utf-8'))
    means=json.loads((analysis/'condition-means.json').read_text(encoding='utf-8'))
    statistics=json.loads((analysis/'paired-statistics.json').read_text(encoding='utf-8'))
    with (analysis/'independently-recomputed.csv').open(encoding='utf-8-sig',newline='') as f:
        seed_rows=list(csv.DictReader(f))
    keyed={(r['distribution'],r['measurementStd'],r['method']):r for r in means}
    conditions=[(d,s) for d in ['gaussian','mixture'] for s in [.01,.1]]
    names={'gaussian':'高斯','mixture':'非高斯混合'}
    methods=['raw','original-KF','selected-single-KF','IMM']
    significant=sum(x['CI95High']<0 and x['pHolmFour']<.05 for x in statistics['primary'])
    text=['# IMM同观测序列测试报告','',
      'Material Passport：academic-research-suite / experiment-agent；2026-10-14。Verification Status：VERIFIED（软件与数值重放、统计复核）。科学证据为开发选参后的6个新种子小规模检验。','',
      f'预先指定的四个条件中，IMM相对于开发集选定的单Q KF有 {significant}/4 项RMSE差值区间完全低于0且Holm校正p<.05。具体条件和数值见下表；不能由此推出集群闭环收益或所有模型下的优势。','',
      '## 实际实验与公平比较','',
      '24条新参考轨迹由直接带噪观测控制产生，种子6101–6106，每种子高斯/同方差非高斯混合 × std .01/.10 rad。所有估计器读取同一条轨迹上的相同选邻、相同观测角度与缺测标记；估计器的输出不反馈到参考轨迹。没有IMM闭环仿真，不能给各离线方法分配GO差异。','',
      '条件沿用此前：N20、M100、随机逐步k15、等权original强度.24、神经核宽.4、beta1000、globalInhibition0、speedScale.05、dt.3、周期L1000，8000步，评价6001:8000。每序列600000个被选边×步，独立统计单位仍为种子，n=6。随机选邻造成同一边偶尔未更新，但没有额外人为遮挡。','',
      '非高斯分布沿用95%小高斯/5%大高斯，大/小标准差比8，除以sqrt(4.15)保持总体方差。混合标签、真实方位和距离均不提供给滤波器。R=std²，没有鲁棒门控、重置或自适应R。','',
      '## 冻结模型与选择','',
      '状态为邻居视线方位及其变化率[theta,omega]，F=[1 dt;0 1]，Q=q*[dt³/3 dt²/2;dt²/2 dt]。两个模型的状态方程相同、过程噪声不同；使用圆周均值、wrapped残差、切空间协方差及稳定log likelihood。它是局部Gaussian近似，不是非Gaussian测量下的精确后验。','',
      f'IMM低q=1e-4、高q={frozen["immHighQ"]:g}；最佳单Q KF q={frozen["bestSingleQ"]:g}。转移概率每模型步为[[.98,.02],[.10,.90]]，初始模型权重[.5,.5]。近反向混合的圆周均值不确定时选择概率更高分量并记录计数。','',
      '选择只使用已查看的5101–5106旧raw参考轨迹。高q候选[.01,.1,1,10]，单Q候选还含1e-4；按24个开发序列等权mean(log(RMSE/rawRMSE))分别选择IMM高q和最佳单Q。参数及分析实现哈希在首个新种子前冻结。全部候选分数保存，没有根据新种子重新选择。','',
      '运行器保留五个标签：raw、原KF、与IMM相同高Q的KF、开发选定单Q KF、IMM。选参后两个高Q单模型恰好相同，表中合并展示四个不同方法；120行汇总不代表120次独立重复。','',
      '## 新种子总体误差','',
      '数值为每条件六种子的RMSE算术均值，单位rad；降幅=1−IMM均值/调参KF均值。不是将360万个时间/边误差当作独立样本。','',
      '| 噪声 | std | 直接观测 | 原KF q=1e-4 | 调参单Q KF | IMM | 相对调参KF降幅 |',
      '|---|---:|---:|---:|---:|---:|---:|']
    for d,s in conditions:
        v=[keyed[d,s,m]['RMSE'] for m in methods]
        text.append(f'| {names[d]} | {s:.2f} | {v[0]:.6f} | {v[1]:.6f} | {v[2]:.6f} | {v[3]:.6f} | {(1-v[3]/v[2])*100:.1f}% |')
    text+=['','预设主要比较为IMM−调参单Q KF；负差值表示IMM误差更小。paired t 95%区间未作同时区间校正，p值对四项比较作Holm校正。n=6时区间依赖种子差值近似正态假设，不能外推为普遍保证。','',
      '| 噪声 | std | RMSE差值 | 95%区间 | Holm p | 更好种子数 |',
      '|---|---:|---:|---|---:|---:|']
    for r in statistics['primary']:
        text.append(f'| {names[r["distribution"]]} | {r["measurementStd"]:.2f} | {r["meanDifference"]:.6f} | [{r["CI95Low"]:.6f}, {r["CI95High"]:.6f}] | {r["pHolmFour"]:.4g} | {r["IMMImprovedSeeds"]}/6 |')
    by_seed={(int(r['seed']),r['distribution'],float(r['measurementStd']),r['method']):r for r in seed_rows}
    adverse=[]
    for r in seed_rows:
        if r['method']=='IMM':
            raw=by_seed[int(r['seed']),r['distribution'],float(r['measurementStd']),'raw']
            if float(r['RMSE'])>float(raw['RMSE']): adverse.append((r,raw))
    text += ['',f'24个新种子×条件序列中，仍有 {len(adverse)} 个序列的IMM总体RMSE高于直接观测；均保留，没有按优势筛除。']
    for r,raw in adverse:
        text.append(f'- 种子{r["seed"]}、{names[r["distribution"]]} std={float(r["measurementStd"]):.2f}：IMM {float(r["RMSE"]):.6f} rad，直接观测 {float(raw["RMSE"]):.6f} rad。这说明降低对调参KF的平均误差不等于每条轨迹都优于直接观测。')
    text+=['','## 降噪与突变代价分别看','',
      '以下是描述性机制指标，阈值在实验前固定：平缓点真实单步变化<.02 rad；突变点>.5 rad。类别只用于评价，没有输入估计器。平缓点可包含此前突变的恢复阶段，不等同于完全稳态段。','',
      '| 条件 | 调参KF平缓点RMSE | IMM平缓点RMSE | 调参KF突变当步RMSE | IMM突变当步RMSE | 调参KF P99 | IMM P99 |',
      '|---|---:|---:|---:|---:|---:|---:|']
    for d,s in conditions:
        a,b=keyed[d,s,'selected-single-KF'],keyed[d,s,'IMM']
        text.append(f'| {names[d]} {s:.2f} | {a["quietRMSE"]:.6f} | {b["quietRMSE"]:.6f} | {a["jumpRMSE"]:.6f} | {b["jumpRMSE"]:.6f} | {a["P99Abs"]:.6f} | {b["P99Abs"]:.6f} |')
    jump_worse=sum(keyed[d,s,'IMM']['jumpRMSE']>keyed[d,s,'selected-single-KF']['jumpRMSE'] for d,s in conditions)
    quiet_better=sum(keyed[d,s,'IMM']['quietRMSE']<keyed[d,s,'selected-single-KF']['quietRMSE'] for d,s in conditions)
    text += ['',f'四条件中，IMM在 {quiet_better}/4 项平缓点RMSE更小，但在 {jump_worse}/4 项突变当步RMSE更大。总体RMSE与瞬时响应需要分开判断；不能写成IMM在所有时刻更准确。P99按nearest-rank计算，突变很稀少时P99也不足以代表最大误差，完整尾部与事件保留。','',
      '恢复定义为突变后第三个连续模型步均fresh且绝对误差<=3std。最多观察30步=9模型时间单位，末尾或期限仍未恢复按右删失保存；事件/边/时间均不作为独立重复。KM恢复比例和“已恢复事件的中位延迟”仅作描述，不能用后者代表所有事件恢复时间。','',
      '| 条件 | 原KF 30步内恢复比例 | 调参KF | IMM |',
      '|---|---:|---:|---:|']
    for d,s in conditions:
        values=[keyed[d,s,m]['KMRecoveredBy30Steps'] for m in ['original-KF','selected-single-KF','IMM']]
        text.append(f'| {names[d]} {s:.2f} | {values[0]*100:.2f}% | {values[1]*100:.2f}% | {values[2]*100:.2f}% |')
    text+=['','## 软件与独立复核','',
      '相同两模型退化到原KF的最大差3.55e-15；角度旋转一致、缺测权重转移、初始化、极端残差稳定和协方差PSD通过。独立矩阵实现的500步合成fixture通过，覆盖带噪、缺测、环绕与近pi真跳变。','',
      f'24条序列的传感器重放与归档样本、逐步观测RMSE、最终随机流一致。Python从新参考位置独立复算真实方位、逐步选邻、周期切换及真实单步变化，再对全部24×8000步作第二KF/IMM矩阵实现重放。最大方位差 {audit["maxAngleDifference"]:.3g} rad，模型概率差 {audit["maxProbabilityDifference"]:.3g}，方差差 {audit["maxVarianceDifference"]:.3g}；独立指标最大差 {audit["maxMetricDifference"]:.3g}。120唯一键、每序列600000评价点与恢复事件计数均复核。','',
      '源码17项快照与哈希、开发输入哈希、新参考/观测/估计哈希、冻结参数、分析源码、日志与审计全部保存。验证期间预先已有四个tracked文件内容哈希保持不变。','',
      '## 能支持的结论与剩余工作','',
      '本轮支持在已列观测/几何条件和新种子上评估双Q角度IMM的估计收益与瞬时响应代价。普通IMM的高响应分量仍可能追随离群观测，不能称为自动鲁棒、安全保证、二维位置恢复或已改善群体运动。','',
      '尚未检验IMM输出接入神经场后的GO、LO、凝聚、连续遮挡预测及闭环适用边界；本轮选择也未覆盖所有转移矩阵、Q候选、噪声分布、N/k或非周期空间。后续闭环比较应继续使用冻结配置与新结果目录，保留调参单Q KF基线。','',
      '分支exp/bearing-imm-replay-pilot-v1，HEAD f8b73cf；没有commit/push/merge。只新增bearingIMMTrackStep.m、testBearingIMMTrack.m、buildBearingReplaySequence.m、replayBearingModels.m、runBearingIMMReplayPilot.m，以及handoff中的计划、fixture、分析和报告脚本。原始模型、此前Gaussian/非Gaussian结果、论文及已有四个源码改动未修改。','',
      '建议保留为实验分支，代码审查后可作为研究实现合并；当前证据不足以把IMM作为已确认的flocking改进替换主线。回退仅撤销新增实现，实验记录保留。','',
      '图表：paired-rmse.png、paired-p99.png为全部六种子的配对图；baseline-worst-event-example.png按固定条件seed6101中原KF最大误差选择，明确是失败示例，不能作为典型效果代表。','']
    path.write_text('\n'.join(text),encoding='utf-8')
    print(json.dumps({'report':str(path),'primarySupportedConditions':significant,'quietImprovedConditions':quiet_better,'jumpInstantWorseConditions':jump_worse},ensure_ascii=False,indent=2))


if __name__=='__main__':
    main()

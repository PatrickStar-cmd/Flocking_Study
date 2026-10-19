"""Finalize only after all scientific runs and independent verification pass."""
from pathlib import Path
import json,csv,hashlib
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parent
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(8*1024**2),b''):h.update(b)
    return h.hexdigest()

def finalize(out):
    out=Path(out);validation=json.loads((out/'checks/independent/validation.json').read_text())
    assert validation['status']=='PASS' and validation['runs']==858
    records=pd.read_csv(out/'RUN_MANIFEST.csv');assert len(records)==858
    agg=pd.read_csv(out/'analysis/paired_statistics.csv');raw=pd.read_csv(out/'analysis/window_metrics_by_seed.csv')
    pair=pd.read_csv(out/'analysis/paired_differences_by_seed.csv')
    figures=out/'analysis/figures'
    assert len(list((figures/'all_seed_time_series').glob('*.pdf')))==183
    assert len(list((figures/'all_seed_trajectories').glob('*.pdf')))==183
    assert np.isfinite(raw.GO).all() and np.isfinite(raw.GO_followers).all()
    core=agg[(agg.case=='main')&(agg.window=='occlusion')&(agg.method=='kf-mean')&(agg.baseline=='hold')&(agg.metric=='GO')].iloc[0]
    follower=agg[(agg.case=='main')&(agg.window=='occlusion')&(agg.method=='kf-mean')&(agg.baseline=='hold')&(agg.metric=='GO_followers')].iloc[0]
    task=agg[(agg.case=='main')&(agg.window=='occlusion')&(agg.method=='kf-mean')&(agg.baseline=='hold')&(agg.metric=='follower_turn_error_deg')].iloc[0]
    status='GO平均改善得到支持' if core.tCI_low>0 else ('GO平均下降' if core.tCI_high<0 else 'GO改善证据不足，区间跨0')
    text=f'''# 先看这里：领航者实验v2

858条科学运行与独立复核已完成；本目录代码和结果为完整实验交付，不是只有待运行代码。

## 这次研究什么

3个领航者接受目标方向提示，17个跟随者通过邻居方位输入响应；遮挡时17个跟随者暂时看不到3个领航者。原近似有序初态、噪声和神经场保持协议配置。主要问题是KF相对Hold是否改善性能，没有转向/无转向机制对照。

开发复刻18条（原3种子×6方法）；新评估840条（30种子、主7方法、时长/时机/40和60架规模扩展）。同一随机种子的各方法共享初态及完整噪声，全部保留，没有调参挑种子。

## 核心实际结果

**{status}。** 主场景遮挡期全群GO(KF-mean)−GO(Hold)为{core['mean']:+.6f}，95%配对t区间[{core.tCI_low:+.6f}, {core.tCI_high:+.6f}]；bootstrap区间[{core.bootstrap_low:+.6f}, {core.bootstrap_high:+.6f}]。

30个种子中{int(core.nPositive)}个正差、{int(core.nNegative)}个负差（仅数值符号）。GO差乘100是百分点，不是相对百分比。

跟随者GO差为{follower['mean']:+.6f}；跟随者目标方向绝对误差差为{task['mean']:+.6f}度，后者负值表示KF误差更小。扩展条件与其他方法结果请看RESULTS_说明.md和paired_statistics.csv，不能用一个主指标代替全部性能。

简单趋势外推获得了大部分相对Hold的收益，KF对它的额外GO增益约0.000136；KF-uncertainty主场景相对KF-mean的额外GO差区间跨0。较早遮挡条件有8/30种子下降，不能称各场景所有种子都改善。原工程诊断的Drop目标误差增量≥3°只在3/30评估种子达到，所有样本保留。

## 推荐阅读顺序

1. RESULTS_说明.md：实际结果与解释。
2. results/full_v2_run1/analysis/figures/main_six_panel.png：主场景七方法曲线。
3. results/full_v2_run1/analysis/figures/main_paired_seed_effects.png：全部30种子配对差。
4. results/full_v2_run1/analysis/figures/robustness_effects.png：全部扩展条件。
5. results/full_v2_run1/checks/independent/validation.json：全部858条独立复算PASS。
6. README.md：代码入口和复现命令；DATA_DICTIONARY.md：指标与时间标签。

原始数据在results/full_v2_run1/runs/，每条HDF5含完整800步输入诊断与801帧状态，另有800行CSV。所有种子的完整尾段图和轨迹在figures/all_seed_time_series及all_seed_trajectories，共各183页。共享初态和噪声在inputs/。原实验目录未改动；本次未上传GitHub。

结果适用于协议中构造的领航者任务与人工方位缺测，不是实体障碍仿真或真实飞行证明。主要报告效果量和区间，无实际安全阈值，也不声称所有其他指标无损。
'''
    (ROOT/'先看这里.md').write_text(text,encoding='utf-8')
    audit=json.loads((out/'audit.json').read_text());audit.update(status='COMPLETE',analysis='COMPLETE',independentValidation='PASS',
        figures=dict(seedTimeSeriesPDF=183,seedTrajectoryPDF=183),primaryMeanGO=core['mean'],primaryCI95=[core.tCI_low,core.tCI_high])
    (out/'audit.json').write_text(json.dumps(audit,indent=2),encoding='utf-8')
    # All sources, docs, inputs and outputs get a reviewable hash. Ignore regenerable bytecode.
    manifest=ROOT/'SHA256SUMS.csv';rows=[]
    for p in sorted(ROOT.rglob('*')):
        if p.is_file() and p!=manifest and '__pycache__' not in p.parts:
            rows.append(dict(file=p.relative_to(ROOT).as_posix(),bytes=p.stat().st_size,sha256=sha(p)))
    with manifest.open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=['file','bytes','sha256']);w.writeheader();w.writerows(rows)
    print(json.dumps(dict(status='COMPLETE',files=len(rows),bytes=sum(r['bytes'] for r in rows),primaryMeanGO=core['mean'])),flush=True)

if __name__=='__main__':
    import sys
    finalize(ROOT/(sys.argv[1] if len(sys.argv)>1 else 'results/full_v2_run1'))

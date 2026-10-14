"""Exploratory paired analysis of frozen B E2 canonical summaries; stdlib only."""
import csv
import hashlib
import json
import math
from pathlib import Path
import random
import statistics as st
from collections import Counter

ROOT=Path(__file__).resolve().parent
INPUT=ROOT.parent/'sq1-c-condition-ledger-v1/B_E2_CANONICAL_RUNS.csv'
NBOOT=20000
BOOT_SEED=20260929
def read(p):
    with p.open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
def write(name,rows):
    with (ROOT/name).open('w',encoding='utf-8-sig',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def quantile(x,p):
    pos=(len(x)-1)*p;i=int(pos);f=pos-i
    return x[i] if i==len(x)-1 else x[i]*(1-f)+x[i+1]*f

rows=read(INPUT)
expected=json.loads((INPUT.parent/'ledger_audit.json').read_text(encoding='utf-8'))['outputs'][INPUT.name]
assert hashlib.sha256(INPUT.read_bytes()).hexdigest()==expected
data={}
for r in rows:
    key=(r['canonical_config_id'],int(r['seed']))
    assert key not in data
    data[key]=r
seeds=list(range(1,21))
assert len(data)==240 and all((f'C{i:02d}',s) in data for i in range(1,13) for s in seeds)

contrasts=[
 ('strength_k2','输入制度','fixed-total − original，k=2','C04','C09'),
 ('strength_k4','输入制度','fixed-total − original，k=4','C05','C10'),
 ('strength_k8','输入制度','fixed-total − original，k=8','C06','C11'),
 ('strength_k12','输入制度','fixed-total − original，k=12','C07','C12'),
 ('count_4_2','fixed-total邻居数','k=4 − k=2','C05','C04'),
 ('count_8_4','fixed-total邻居数','k=8 − k=4','C06','C05'),
 ('count_12_8','fixed-total邻居数','k=12 − k=8','C07','C06'),
 ('count_19_12','fixed-total邻居数','k=19 − k=12','C08','C07'),
 ('policy_concentrated','k=8策略','angle-concentrated − random','C01','C06'),
 ('policy_balanced','k=8策略','balanced/angle-dispersed − random','C02','C06'),
 ('policy_nearest','k=8策略','nearest − random','C03','C06')]
metrics={'meanGO':'GO','meanLO':'LO','meanPairDistanceNorm':'归一化平均对距'}
# Same bootstrap index draws across all conditions preserve seed-level dependence.
rng=random.Random(BOOT_SEED)
weights=[list(Counter(rng.randrange(20) for _ in range(20)).items()) for _ in range(NBOOT)]
summary=[];differences=[]
for cid,family,label,treat,ref in contrasts:
    for metric,title in metrics.items():
        a=[float(data[(treat,s)][metric]) for s in seeds]
        b=[float(data[(ref,s)][metric]) for s in seeds]
        delta=[x-y for x,y in zip(a,b)]
        boot=sorted(sum(delta[i]*n for i,n in w)/20 for w in weights)
        mean=st.mean(delta)
        assert abs(mean-(st.mean(a)-st.mean(b)))<1e-12
        loo=[(sum(delta)-x)/19 for x in delta]
        directional=[-x if metric=='meanPairDistanceNorm' else x for x in delta]
        summary.append({'contrast_id':cid,'family':family,'contrast':label,
            'treatment_config':treat,'reference_config':ref,'metric':metric,'n_pairs':20,
            'treatment_mean':st.mean(a),'reference_mean':st.mean(b),
            'mean_difference':mean,'median_difference':st.median(delta),'sd_difference':st.stdev(delta),
            'bootstrap_ci95_low':quantile(boot,.025),'bootstrap_ci95_high':quantile(boot,.975),
            'favorable_direction_count':sum(x>1e-12 for x in directional),
            'unfavorable_direction_count':sum(x < -1e-12 for x in directional),
            'tie_count':sum(abs(x)<=1e-12 for x in directional),
            'leave_one_seed_out_mean_min':min(loo),'leave_one_seed_out_mean_max':max(loo),
            'interval_type':'pointwise percentile paired bootstrap; unadjusted; exploratory'})
        for s,x,y,d in zip(seeds,a,b,delta):
            differences.append({'contrast_id':cid,'metric':metric,'seed':s,
                                'treatment_value':x,'reference_value':y,'difference':d})
assert len(summary)==33 and len(differences)==660
assert all(r['favorable_direction_count']+r['unfavorable_direction_count']+r['tie_count']==20 for r in summary)
write('paired_effects.csv',summary)
write('paired_seed_differences.csv',differences)
audit={'B_commit':'63d56b51cb74bca509b0fb7dc29d1d5997d3a1ae','input':str(INPUT),
       'input_sha256':expected,'simulation_seeds':seeds,'canonical_rows':240,
       'contrasts':11,'metrics':list(metrics),'paired_differences':660,
       'bootstrap_replicates':NBOOT,'bootstrap_rng_seed':BOOT_SEED,
       'bootstrap_rng_seed_role':'computational resampling only, not a new simulation seed',
       'resampling_unit':'whole seed; shared draws across all contrasts and metrics',
       'interval':'linear-interpolated 2.5%/97.5% percentiles, pointwise, no multiplicity correction',
       'evidence':'post-hoc exploratory analysis of previously used development seeds; no new simulation',
       'checks':'PASS source hash, 12x20 unique pairs, difference of means identity, 33 summaries/660 deltas, direction counts',
       'outputs':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in ROOT.glob('*.csv')}}
(ROOT/'analysis_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf-8')
lines=['# B-E2：240行规范数据的探索性配对分析','',
'日期：2026-09-29；B源提交：`63d56b5`；本地分支：`YCwork@71583f9`。',
'', '## 1. 分析方法与解释边界','',
'输入为经核查的12个规范条件×20个种子，共240行。原始480行没有修改；跨章节重复、策略别名及全邻居等价制度不重复计样。',
'',
'本次在已查看数据的背景下选择11个描述性比较，分别分析GO、LO、归一化平均对距，共33项。所有差异为处理减参照；GO/LO正值表示更有序，对距负值表示更紧凑，不自动表示更安全或总体更优。',
'',
'按种子1–20整体重采样20,000次，各比较共用索引抽样，采用百分位法逐项95%区间；分析随机种子20260929只用于重采样。未做多重比较校正，不能将逐项区间解释为33项同时95%覆盖，也不报告确认性显著性或普适邻居阈值。',
'',
'B代码在相同seed下初始化，但不同策略的随机调用及闭环轨迹未在本轮逐步验收。因此配对是按已报告的共同seed设计，不等于已证明所有干预共享完整随机驱动。',
'',
'本数据均为无遮挡、动态选邻、等权的探索结果。balanced/angle-dispersed和nearest使用距离辅助；策略差异不能解释为纯方位因果效应。',
'',
'下面“同向种子”指GO/LO增加或对距减少的种子数，不是成功率，也不使用它宣称可靠性。']
lines += ['', '## 2. 当前可支持的探索性判断', '',
    '1. **输入制度是值得优先验证的因素。** k=2/4/8时，fixed-total相对original的GO差为+0.233/+0.365/+0.384，逐项区间均在零以上，对距差均为负。k=12的GO差约+0.104，其区间跨零，不能把前三个k的结果外推至全部k。',
    '2. **固定总幅值不产生清晰的邻居数单调规律。** GO条件均值按k=2/4/8/12/19约为0.453/0.607/0.642/0.618/0.737；相邻比较中只有4−2的GO逐项区间不跨零。8−4、12−8、19−12均跨零，不表示这些条件等效，也不能据此选择最优k。',
    '3. **角度分散不是当前数据中的自动优势。** balanced/angle-dispersed相对random的GO差约−0.150，区间约[−0.253,−0.037]；LO差也为负。两种策略仅3/20种子在GO方向上有利，对距差区间跨零。它们在动态重选、距离辅助和时序网络上不同，不能据此断言“角度分散本身有害”。',
    '4. **集中和最近邻策略的紧凑程度值得检查。** 相对random，集中策略对距差约+0.0755，最近邻约+0.0468，逐项区间均在零以上；但各自GO区间跨零。需要分别报告有序性与紧凑程度，不能只看单一分数。',
    '',
    '这些比较属于事后探索、未经多重校正。区间跨零不证明无效或等效；不跨零也不构成独立确认。更紧凑不等于碰撞更安全。',
    '',
    '**下一步实验建议：**先补正式轨迹和核心回归，复核动态随机选邻的实际网络时序；C再冻结代表性条件及缺失协议。不要根据本批高GO结果只保留成功seed或有利检查点。']
for family in ['输入制度','fixed-total邻居数','k=8策略']:
    lines+=['',f'## {family}','', '| 比较 | 指标 | 处理均值 | 参照均值 | 配对平均差 | 逐项95%区间 | 同向种子 |', '|---|---|---:|---:|---:|---|---:|']
    for r in summary:
        if r['family']!=family:continue
        lines.append(f"| {r['contrast']} | {metrics[r['metric']]} | {r['treatment_mean']:.5f} | {r['reference_mean']:.5f} | {r['mean_difference']:+.5f} | [{r['bootstrap_ci95_low']:+.5f}, {r['bootstrap_ci95_high']:+.5f}] | {r['favorable_direction_count']}/20 |")
lines+=['','## 文件和复现','',
'- [paired_effects.csv](paired_effects.csv)：33项配对效应、均值/中位数、标准差、区间及逐一删seed后的均值范围。',
'- [paired_seed_differences.csv](paired_seed_differences.csv)：660条逐seed差异，供独立复核。',
'- [analysis_audit.json](analysis_audit.json)：来源哈希、计算参数和输出哈希。',
'- [B_REQUIRED_DELIVERABLES.md](B_REQUIRED_DELIVERABLES.md)：B需要补交的具体数据/测试、理由和责任边界。',
'',
'在仓库根目录执行：`python reviews/sq1-c-paired-analysis-v1/analyze_paired.py`。无需运行MATLAB、无需第三方Python库。',
'',
'未运行新仿真、未更改B代码、未改历史结果或main。该目录为新增分析，可单独撤回；共享提交后使用git revert。']
(ROOT/'PAIRED_ANALYSIS.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
for r in summary:
    print(r['contrast_id'],r['metric'],round(r['mean_difference'],5),
          (round(r['bootstrap_ci95_low'],5),round(r['bootstrap_ci95_high'],5)),r['favorable_direction_count'])

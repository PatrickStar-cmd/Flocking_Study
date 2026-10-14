"""Paired seed bootstrap of predefined windows; no time-point resampling."""
import csv, hashlib, json, random, statistics as st
from collections import Counter
from pathlib import Path

ROOT=Path(__file__).resolve().parent
REPO=ROOT.parents[1]
def read(p):
    with p.open(encoding='utf-8-sig',newline='') as f: return list(csv.DictReader(f))
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def q(v,p):
    x=(len(v)-1)*p; i=int(x); f=x-i
    return v[i]*(1-f)+v[min(i+1,len(v)-1)]*f

idx=ROOT.parent/'sq1-c-raw-index-v1/CANONICAL_RAW_INDEX.csv'
previous=json.loads((idx.parent/'audit.json').read_text(encoding='utf-8'))
assert sha(idx)==previous['outputs'][idx.name]
for r in read(idx): assert sha(REPO/r['raw_path'])==r['raw_sha256']
raw=read(ROOT/'run_windows.csv'); pairs=read(ROOT/'paired_windows.csv')
assert len(raw)==2160 and len(pairs)==1980
old=read(ROOT.parent/'sq1-c-paired-analysis-v1/paired_effects.csv')
keys=[(r['contrast_id'],r['metric']) for r in old]
data={(r['contrast_id'],r['metric'],r['window'],int(r['seed'])):r for r in pairs}
assert len(data)==1980
rng=random.Random(20260929)
draws=[list(Counter(rng.randrange(20) for _ in range(20)).items()) for _ in range(20000)]
results=[]; max_identity=0; max_old=0
for cid,metric in keys:
    delta={}
    for window in ['main','early','late']:
        delta[window]=[]
        for s in range(1,21):
            r=data[cid,metric,window,s]
            v=float(r['treatment'])-float(r['reference'])
            assert abs(v-float(r['difference']))<1e-12
            delta[window].append(v)
    max_identity=max(max_identity,max(abs(m-(e+l)/2) for m,e,l in zip(delta['main'],delta['early'],delta['late'])))
    delta['late_minus_early']=[l-e for l,e in zip(delta['late'],delta['early'])]
    for window,v in delta.items():
        boot=sorted(sum(v[i]*n for i,n in draw)/20 for draw in draws)
        result=dict(contrast_id=cid,metric=metric,window=window,n_pairs=20,mean_difference=st.mean(v),
                    median_difference=st.median(v),sd_difference=st.stdev(v),ci95_low=q(boot,.025),ci95_high=q(boot,.975),
                    positive_count=sum(x>1e-12 for x in v),negative_count=sum(x< -1e-12 for x in v),
                    interval='pointwise paired percentile bootstrap; unadjusted; exploratory')
        results.append(result)
        if window=='main':
            o=next(r for r in old if (r['contrast_id'],r['metric'])==(cid,metric))
            max_old=max(max_old,abs(result['mean_difference']-float(o['mean_difference'])),
                        abs(result['ci95_low']-float(o['bootstrap_ci95_low'])),abs(result['ci95_high']-float(o['bootstrap_ci95_high'])))
assert max_identity<1e-12 and max_old<1e-12 and len(results)==132
out=ROOT/'window_effects.csv'; assert not out.exists()
with out.open('w',encoding='utf-8-sig',newline='') as f:
    writer=csv.DictWriter(f,fieldnames=list(results[0]));writer.writeheader();writer.writerows(results)
audit=dict(evidence='exploratory; previously used development seeds',bootstrap_replicates=20000,
           bootstrap_seed=20260929,resampling_unit='whole seed, shared draws across all comparisons/windows',
           max_full_window_identity_error=max_identity,max_previous_mean_and_ci_error=max_old,
           input_index_sha256=sha(idx),raw_hashes_verified=240,results=132,
           files={p.name:sha(p) for p in ROOT.iterdir() if p.suffix in ('.csv','.m','.py','.png','.svg')})
(ROOT/'analysis_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf-8')
print('PASS',json.dumps({k:v for k,v in audit.items() if k!='files'},ensure_ascii=False))
for cid,metric in keys:
    selected=[r for r in results if r['contrast_id']==cid and r['metric']==metric]
    print(cid,metric,' '.join(f"{r['window']}={r['mean_difference']:+.5f}[{r['ci95_low']:+.5f},{r['ci95_high']:+.5f}]" for r in selected))
lines=['# 时间过程与窗口配对分析','',
       '来源 WHQ_work@0036487；240规范运行、12条件、每条件20个已使用seed。探索性分析，不新增仿真。',
       '', '## 如何读图','',
       '每张图上排是处理条件（蓝）和参照条件（橙）的20种子平均曲线；下排灰线是各seed的配对差，蓝线是平均差。上排均值不表示每条轨迹均如此。下排GO/LO正值表示处理更有序，对距负值表示更紧凑。',
       '', '横轴为模拟步，乘0.3得到模型时间。6000步后是原评价窗，7000步把它分成等长两半。所有8000点均绘出，未平滑；未绘置信带，不能据曲线某一点宣布显著或持续成功。',
       '', '## 方法与核查','',
       '沿用原11组对比。主窗口6001–8000，前半6001–7000，后半7001–8000；另计算后半配对差减前半配对差，以描述效应是否随窗口变化。不是恢复时间、稳态检验或时间趋势的因果估计。',
       '', '按完整seed共用索引重采样20,000次，计算种子20260929，仅用于计算。线性插值百分位95%区间逐项未校正；132项结果与既有数据均为探索性。不从时间点数推导样本量。',
       '',f'主窗均值与旧CSV以及旧配对bootstrap结果一致，旧配对均值/区间最大差{max_old:.3g}；主窗等于两半平均的最大误差{max_identity:.3g}。240个MAT哈希全部通过。',
       '', '时间点自相关，图中均值可能掩盖种子差异；方位策略还混合邻居身份、距离辅助和动态网络。这里未验证跨策略共同随机驱动，不能据此分离纯方位因果效应。']
for j,cid in enumerate(dict.fromkeys(k[0] for k in keys),1):
    label=next(r['contrast'] for r in old if r['contrast_id']==cid)
    lines+=['',f'## {j}. {label}','',f'![{label}]({j:02d}_{cid}.png)','',
            '|指标|前半平均差|后半平均差|后半−前半 [逐项95%区间]|','|---|---:|---:|---|']
    for metric in ['meanGO','meanLO','meanPairDistanceNorm']:
        rr={r['window']:r for r in results if (r['contrast_id'],r['metric'])==(cid,metric)}
        d=rr['late_minus_early']; lines.append(f"|{metric}|{rr['early']['mean_difference']:+.5f}|{rr['late']['mean_difference']:+.5f}|{d['mean_difference']:+.5f} [{d['ci95_low']:+.5f}, {d['ci95_high']:+.5f}]|")
lines+=['','## 文件与边界','',
        '- `run_windows.csv`：2160行，各运行/指标/窗口的均值、时间标准差和极差。',
        '- `paired_windows.csv`：1980行，逐seed窗口配对值。',
        '- `window_effects.csv`：132项配对均值、区间、方向计数；方向计数不是成功率。',
        '- 11张PNG及同名SVG；DESIGN.md为本轮固定方案，审计JSON和MATLAB日志记录验证。',
        '', '当前数据没有观测缺失，不能据此证明缺失恢复、预测收益或安全；仍需B完成四分支接口及测试，C冻结缺失协议。',
        '', '工作分支YCwork@71583f9，新增内容仅在本目录，未提交或推送。原始代码、历史数据、此前分析和main未改动。建议保留为探索性材料；不作为确认性结论直接合并main。回退仅涉及本新增目录，共享提交后使用git revert。']
(ROOT/'REPORT.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')

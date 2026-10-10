"""Read raw MATLAB v7.3 data independently; aggregate by seed and plot.
Dependencies: numpy, pandas, h5py, matplotlib. No MATLAB summary code is used.
"""
from pathlib import Path
import sys, json
import numpy as np
import pandas as pd
import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parent
MODES=['fresh','drop','drop-normalized','hold','KF-mean']
COLORS=['#777777','#cc6677','#ddaa33','#4477aa','#228833']
WINDOWS={'occlusion':(6201,6234),'recovery100':(6235,6334),'recovery300':(6235,6534),'event400':(6201,6600)}
SCANS=[(['h2_random_A228','h4_random_A228','h6_random_A228'],['2 (25%)','4 (50%)','6 (75%)'],'Hidden count h'),
       (['h4_random_A228','h4_concentrated_A228','h4_dispersed_A228'],['Random','Concentrated','Dispersed'],'Hidden angular structure'),
       (['h4_random_A114','h4_random_A228','h4_random_A456'],['0.114','0.228','0.456'],'Full-observation weight A0')]
DRAW=np.random.default_rng(20261010).integers(0,5,size=(5000,5))

def array2(ds): return np.asarray(ds).T.squeeze()
def core(path,name='seg'):
    with h5py.File(path,'r') as f:
        d=f[name+'/data']
        return {k:array2(d[k]) for k in ['step','GO','LO','pairDistance','A','Mass','Peak','k_active','selfFreshRelativeL2']}

def recompute(result,expected):
    rows=[]; contracts=[]; cache={}
    for folder in sorted((result/'conditions').iterdir()):
        parts=folder.name.split('_'); seed=int(parts[0][4:]); cond='_'.join(parts[1:]); h=int(parts[1][1:]); A0=int(parts[3][1:])/1000
        refid=f'fresh_seed{seed}_A{round(A0*1000)}'
        if refid not in cache: cache[refid]=core(result/'references'/f'{refid}.mat','fresh')
        ref=cache[refid]
        for mode in MODES:
            d=ref if mode=='fresh' else core(folder/f'{mode}.mat')
            steps=d['step'].astype(int); idx=np.searchsorted(ref['step'],steps)
            assert np.array_equal(steps,ref['step'][idx])
            for window,(a,b) in WINDOWS.items():
                ix=(steps>=a)&(steps<=b); jj=idx[ix]
                row=dict(seed=seed,condition_id=cond,h=h,structure=parts[2],A0=A0,mode=mode,window=window,nSteps=int(ix.sum()))
                row.update(meanGO=float(d['GO'][ix].mean()),meanLO=float(d['LO'][ix].mean()),meanPairDistanceNorm=float(d['pairDistance'][ix].mean()),
                    pairDistanceIncrease=float((d['pairDistance'][ix]-ref['pairDistance'][jj]).mean()),
                    GOLoss=float((ref['GO'][jj]-d['GO'][ix]).mean()),LOLoss=float((ref['LO'][jj]-d['LO'][ix]).mean()),
                    meanActive=float(d['k_active'][ix].mean()),meanA=float(d['A'][ix].mean()),meanMass=float(d['Mass'][ix].mean()),meanPeak=float(d['Peak'][ix].mean()),
                    selfFreshL2Mean=float(np.nanmean(d['selfFreshRelativeL2'][ix])),selfFreshL2Max=float(np.nanmax(d['selfFreshRelativeL2'][ix])),
                    meanAbsolutePairDeviation=float(np.abs(d['pairDistance'][ix]-ref['pairDistance'][jj]).mean()))
                rows.append(row)
            if mode=='fresh': continue
            with h5py.File(folder/f'{mode}.mat','r') as f:
                dd=f['seg/data']; selected=np.asarray(dd['selected']).astype(bool); fresh=np.asarray(dd['fresh']).astype(bool)
                hidden=np.asarray(dd['hidden']).astype(bool); occ=(steps>=6201)&(steps<=6234)
                obs=np.asarray(dd['observedBearing']); weights=np.asarray(dd['observedWeight']); age=np.asarray(dd['age'])
                assert np.isnan(obs[hidden]).all() and np.isnan(weights[hidden]).all()
                assert np.isfinite(obs[fresh]).all() and np.all(age[fresh]==0)
                assert np.all(selected.sum(axis=1)==8) and np.all(hidden[occ].sum(axis=1)==h) and not hidden[~occ].any()
                expectedA=A0*(1-h/8) if mode=='drop' else A0
                assert np.max(np.abs(d['A'][occ]-expectedA))<1e-12
                assert np.all(d['k_active'][occ]==(8-h if mode in ['drop','drop-normalized'] else 8))
                assert np.max(np.abs(d['A'][~occ]-A0))<1e-12 and np.max(d['selfFreshRelativeL2'][steps==6235])<1e-12
                assert np.max(np.abs(age[steps==6234][hidden[steps==6234]]-10.2))<1e-12
                assert np.array_equal(array2(dd['driverIndex']).astype(int),steps)
                p11=np.asarray(dd['predictionVariance']); p12=np.asarray(dd['predictionP12']); p22=np.asarray(dd['predictionP22'])
                det=(p11*p22-p12**2)[selected]
                assert np.isfinite(p11[selected]).all() and np.isfinite(p22[selected]).all()
                assert p11[selected].min()>=0 and p22[selected].min()>=0 and det.min()>=-1e-12
                contracts.append(dict(run_id=f'{folder.name}_{mode}',status='PASS',nSteps=len(steps),hiddenEdgeSteps=int(hidden.sum()),minCovarianceDet=float(det.min())))
    raw=pd.DataFrame(rows); merged=raw.merge(expected,on=['seed','condition_id','mode','window'],suffixes=('_raw','_matlab'),validate='one_to_one'); assert len(merged)==700
    checks=[]
    for name in ['nSteps','meanGO','meanLO','meanPairDistanceNorm','pairDistanceIncrease','GOLoss','LOLoss','meanActive','meanA','meanMass','meanPeak','selfFreshL2Mean','selfFreshL2Max']:
        delta=float(np.max(np.abs(merged[name+'_raw']-merged[name+'_matlab'])))
        checks.append(dict(metric=name,maxAbsDifference=delta,status='PASS' if delta<1e-10 else 'FAIL'))
    assert all(x['status']=='PASS' for x in checks)
    return raw,pd.DataFrame(contracts),pd.DataFrame(checks)

def verify_replays(result,same,offline):
    ss=[]; ds=[]
    for folder in sorted((result/'conditions').iterdir()):
        seed=int(folder.name.split('_')[0][4:]); cond='_'.join(folder.name.split('_')[1:])
        with h5py.File(folder/'same_state_replay.mat','r') as f:
            rep=f['replay']; step=array2(rep['step']); occ=(step>=6201)&(step<=6234)
            candidates=rep['candidates']; ref=np.asarray(candidates['fresh/field'])[occ]
            for mode in MODES:
                field=np.asarray(candidates[mode.replace('-','_')+'/field'])[occ]
                rel=np.sqrt(((field-ref)**2).sum(axis=2)/(ref**2).sum(axis=2))
                row=same[(same.seed==seed)&(same.condition_id==cond)&(same['mode']==mode)].iloc[0]
                delta=max(abs(float(rel.mean())-row.relativeL2Mean),abs(float(rel.max())-row.relativeL2Max)); assert delta<1e-10
                ss.append(dict(seed=seed,condition_id=cond,mode=mode,maxDifference=delta,status='PASS'))
        with h5py.File(folder/'prediction_diagnostics.mat','r') as f:
            for name,R in [('R_6',1e-6),('R_4',1e-4),('R_2',1e-2)]:
                g=f['diagnostic/'+name]; e=np.asarray(g['error']); v=np.asarray(g['variance']); he=np.asarray(g['holdError']); valid=np.isfinite(e)
                rmse=float(np.sqrt((e[valid]**2).mean())); hold=float(np.sqrt((he[valid]**2).mean())); coverage=float((np.abs(e[valid])<=1.96*np.sqrt(v[valid])).mean())
                row=offline[(offline.seed==seed)&(offline.condition_id==cond)&np.isclose(offline.R,R,rtol=0,atol=1e-12)].iloc[0]
                delta=max(abs(rmse-row.rmseKF),abs(hold-row.rmseHold),abs(coverage-row.coverage95)); assert delta<1e-10
                ds.append(dict(seed=seed,condition_id=cond,R=R,maxDifference=delta,status='PASS'))
    return pd.DataFrame(ss),pd.DataFrame(ds)

def aggregate(df,groups,metrics):
    rows=[]
    for keys,g in df.groupby(groups,sort=False):
        if not isinstance(keys,tuple): keys=(keys,)
        g=g.sort_values('seed'); assert np.array_equal(g.seed.to_numpy(),np.arange(1,6))
        for metric in metrics:
            v=g[metric].to_numpy(); lo,hi=np.quantile(v[DRAW].mean(axis=1),[.025,.975])
            rows.append(dict(zip(groups,keys))|dict(metric=metric,nSeeds=5,mean=float(v.mean()),SD=float(v.std(ddof=1)),minimum=float(v.min()),maximum=float(v.max()),CI95_low=float(lo),CI95_high=float(hi)))
    return pd.DataFrame(rows)

def plots(result,out,raw,same,diag):
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})
    def save(fig,name):
        fig.savefig(out/'figures'/f'{name}.png',dpi=170,bbox_inches='tight'); fig.savefig(out/'figures'/f'{name}.pdf',bbox_inches='tight'); plt.close(fig)
    for src,metric,name,ylabel in [(same,'relativeL2Mean','same_state_L2','Same-state input relative L2 (lower)'),
        (raw[raw.window=='recovery100'],'pairDistanceIncrease','closed_loop_recovery100','Pair-distance change / L (signed)'),
        (raw[raw.window=='recovery100'],'meanAbsolutePairDeviation','absolute_pair_deviation100','Absolute pair-distance deviation / L')]:
        fig,axes=plt.subplots(1,3,figsize=(15,4),constrained_layout=True)
        for ax,(conds,labels,xlabel) in zip(axes,SCANS):
            for mode,color in zip(MODES[1:],COLORS[1:]):
                means=[]
                for i,c in enumerate(conds):
                    v=src[(src.condition_id==c)&(src['mode']==mode)].sort_values('seed')[metric].to_numpy(); means.append(v.mean())
                    ax.scatter(np.full(5,i)+(MODES.index(mode)-2.5)*.03,v,s=10,alpha=.3,color=color)
                ax.plot(range(3),means,'o-',label=mode,color=color,ms=4)
            ax.set_xticks(range(3),labels); ax.set_xlabel(xlabel); ax.axhline(0,color='#888888',lw=.7,ls='--')
        axes[0].set_ylabel(ylabel); axes[-1].legend(fontsize=8); save(fig,name)
    primary=diag[diag.primaryR==1]; fig,axes=plt.subplots(1,3,figsize=(15,4),constrained_layout=True)
    for ax,(conds,labels,xlabel) in zip(axes,SCANS):
        for metric,label,color in [('rmseKF','KF-mean',COLORS[4]),('rmseHold','hold',COLORS[3])]:
            ax.plot(range(3),[primary[primary.condition_id==c][metric].mean() for c in conds],'o-',label=label,color=color)
        ax.set_xticks(range(3),labels); ax.set_xlabel(xlabel)
    axes[0].set_ylabel('Circular bearing RMSE (rad; lower)'); axes[-1].legend(); save(fig,'bearing_prediction_RMSE')
    curves={m:[] for m in MODES[1:]}; inputs={m:[] for m in MODES[1:]}
    for seed in range(1,6):
        folder=result/'conditions'/f'seed{seed}_h4_random_A228'; ref=core(result/'references'/f'fresh_seed{seed}_A228.mat','fresh')
        for mode in MODES[1:]:
            d=core(folder/f'{mode}.mat'); curves[mode].append(d['pairDistance']-ref['pairDistance'][200:])
        with h5py.File(folder/'same_state_replay.mat','r') as f:
            for mode in MODES[1:]: inputs[mode].append(array2(f['replay/candidates/'+mode.replace('-','_')+'/relativeL2']).mean(axis=1))
    fig,axes=plt.subplots(1,2,figsize=(12,4),constrained_layout=True); x=(np.arange(400)+1)*.3
    for mode,color in zip(MODES[1:],COLORS[1:]):
        for ax,source in zip(axes,[inputs,curves]):
            y=np.array(source[mode]); ax.plot(x,y.mean(axis=0),color=color,label=mode); ax.fill_between(x,y.min(axis=0),y.max(axis=0),color=color,alpha=.08)
    for ax in axes: ax.axvspan(0,10.2,color='#999999',alpha=.15); ax.set_xlabel('Time since checkpoint (model units)')
    axes[0].set_ylabel('Same-state input relative L2'); axes[1].set_ylabel('Signed pair-distance change / L'); axes[1].legend(fontsize=8); save(fig,'center_condition_time_curves')
    masks=pd.read_csv(result/'hidden_sets.csv'); z=masks[(masks.A0==.228)&(masks.h==4)]
    fig,ax=plt.subplots(figsize=(8,4),constrained_layout=True)
    values=[z[z.structure==s].groupby('seed').R1.mean().to_numpy() for s in ['random','concentrated','dispersed']]
    ax.boxplot(values); ax.set_xticks([1,2,3],['Random','Concentrated','Dispersed']); ax.set_ylabel('Mean hidden R1 per seed'); save(fig,'realized_hidden_R1')
    fig,axes=plt.subplots(1,2,figsize=(11,4),constrained_layout=True)
    center=diag[diag.condition_id=='h4_random_A228']
    for R,marker in [(1e-6,'o'),(1e-4,'s'),(1e-2,'^')]:
        g=center[np.isclose(center.R,R,rtol=0,atol=1e-12)].sort_values('seed')
        axes[0].plot(g.seed,g.rmseKF,marker+'-',label=f'R={R:g}')
        axes[1].plot(g.seed,g.coverage95,marker+'-',label=f'R={R:g}')
    h=center[center.primaryR==1].sort_values('seed'); axes[0].plot(h.seed,h.rmseHold,'x--',color='#555555',label='hold')
    axes[1].axhline(.95,color='#555555',ls='--'); axes[1].set_ylim(0,1.04)
    for ax in axes: ax.set_xlabel('Development seed'); ax.set_xticks(range(1,6)); ax.legend(fontsize=8)
    axes[0].set_ylabel('Bearing RMSE (rad)'); axes[1].set_ylabel('Empirical 95% interval coverage'); save(fig,'R_sensitivity_by_seed')

def write_report(out,raw,same,diag,contrasts):
    primary=diag[diag.primaryR==1]; center='h4_random_A228'; c=primary[primary.condition_id==center]; si=same[same.condition_id==center]; cl=raw[(raw.condition_id==center)&(raw.window=='recovery100')]
    lines=['# 10.10周KF遮挡实验：实际结果与解释','',
      '本轮完成代码、实际仿真和独立复核。结果属于5个开发种子的单因素pilot，不能当作正式独立验证或完整三因素交互结论。','',
      '## 1. 实际完成了什么','',
      '- 种子1—5；20个体，每个体固定8个邻居，100个神经元；无噪声、无距离权重。',
      '- 第6201—6234步隐藏邻居方位，持续34步（10.2模型时间），第6235步恢复，运行至6600步。',
      '- 数量2/4/6、方位random/concentrated/dispersed、完整观测总权重A0=.114/.228/.456；三个扫描共享中心条件，共7个唯一条件。',
      '- fresh、drop、drop-normalized、hold、KF-mean五种方法；155条闭环轨迹、35组同状态回放、105项离线R诊断。',
      '- 10项MATLAB软件检查通过；Python从MAT复算700行群体窗口、175行同状态输入、105行离线诊断，并检查140条非fresh轨迹的数据契约。',
      '- KF状态为方位角和角速度；q=1e-6，主R=1e-4为假定方差/数值正则。另两个R仅做离线诊断。','',
      '## 2. 中心条件：h4/random/A0=.228','',
      f'五种子离线方位RMSE均值：KF={c.rmseKF.mean():.6g} rad，hold={c.rmseHold.mean():.6g} rad；KF在{int((c.rmseKF<c.rmseHold).sum())}/5个seed中更低。',
      f'主R的95%区间经验覆盖率平均{100*c.coverage95.mean():.2f}%，seed范围{100*c.coverage95.min():.2f}%—{100*c.coverage95.max():.2f}%；恢复首步平均NIS={c.recoveryNISmean.mean():.6g}。R是假定量，不把高覆盖率直接解释为传感器噪声已校准。','',
      '|方法|同状态输入L2均值↓|恢复100步对距增量/L（有符号）|恢复100步GO损失|','|---|---:|---:|---:|']
    for mode in MODES:
        s=si[si['mode']==mode]; g=cl[cl['mode']==mode]; lines.append(f'|{mode}|{s.relativeL2Mean.mean():.6g}|{g.pairDistanceIncrease.mean():.6g}|{g.GOLoss.mean():.6g}|')
    lines+=['','中心条件的逐seed预测诊断：','','|seed|KF RMSE|hold RMSE|95%覆盖率|恢复首步平均NIS|','|---|---:|---:|---:|---:|']
    for _,r in c.sort_values('seed').iterrows():
        lines.append(f'|{int(r.seed)}|{r.rmseKF:.6g}|{r.rmseHold:.6g}|{100*r.coverage95:.2f}%|{r.recoveryNISmean:.6g}|')
    bad=c[c.rmseKF>c.rmseHold]
    lines+=['',f"中心条件KF比hold预测更差的seed为：{', '.join(str(int(x)) for x in bad.seed) or '无'}。协方差诊断存在种子差异，应保留不利seed，不能据整体均值称KF普遍有效。"]
    lines+=['','## 3. 三条扫描结果','','下表均为五个seed的均值。有符号对距增量为branch−fresh；负值说明更紧凑，不能单凭负值断言控制更好。另保存绝对偏离作为补充，不替换预定主指标。','',
      '|条件|KF输入L2|hold输入L2|KF恢复100步对距增量/L|hold恢复100步对距增量/L|KF方位RMSE(rad)|','|---|---:|---:|---:|---:|---:|']
    conditions=list(dict.fromkeys(c for scan in SCANS for c in scan[0]))
    for cond in conditions:
        sk=same[(same.condition_id==cond)&(same['mode']=='KF-mean')]; sh=same[(same.condition_id==cond)&(same['mode']=='hold')]
        ck=raw[(raw.condition_id==cond)&(raw['mode']=='KF-mean')&(raw.window=='recovery100')]; ch=raw[(raw.condition_id==cond)&(raw['mode']=='hold')&(raw.window=='recovery100')]
        d=primary[primary.condition_id==cond]
        lines.append(f'|{cond}|{sk.relativeL2Mean.mean():.6g}|{sh.relativeL2Mean.mean():.6g}|{ck.pairDistanceIncrease.mean():.6g}|{ch.pairDistanceIncrease.mean():.6g}|{d.rmseKF.mean():.6g}|')
    p=contrasts[contrasts.measure=='pairDistanceIncrease']; a=contrasts[contrasts.measure=='meanAbsolutePairDeviation']
    lines+=['',f'35个seed×条件配对中，KF的恢复100步有符号对距增量小于hold的数量为{int((p.KFminusHold<0).sum())}/35；绝对对距偏离小于hold的数量为{int((a.KFminusHold<0).sum())}/35。它们共享5个seed，不能当作35个独立重复。','',
      '## 4. 限制与后续','','- 这是人工屏蔽固定邻居方位，尚未加入真实几何遮挡。',
      '- 总强度是完整观测权重A0，不是被隐藏边权重AH；等权时AH=h*A0/8。',
      '- 每个A0重新运行完整预热，因此强度变化包括前缀动力学变化。',
      '- 未实施KF-uncertainty，不能从KF均值结果推断不确定性传播收益。',
      '- 本轮沿用历史固定候选，没有训练择参；seed1—5是已使用开发种子，未完成独立确认性验证。',
      '- bootstrap以seed为单位，5000次、分析种子20261010，所有条件共享采样索引。5seed区间只能作探索性描述。',
      '- concentrated/dispersed按R1极值选择，不等于前方/后方遮挡；应查看实际隐藏R1。',
      '- 未定义成功率、失效/恢复时间，不能从时间均值虚构这些结论；也未检验完整三因素交互。','',
      '## 5. C的复核入口','','先读DESIGN.md、DATA_DICTIONARY.md和README.md。有效数据仅在results/pilot_seed1to5_v1。',
      '查看RUN_MANIFEST.csv、audit.json、logs/software_checks.json、analysis/validation.json与independent_recompute_check.csv；原始MAT位于references/和conditions/。',
      '图表在analysis/figures/，逐seed配对数值在paired_KF_contrasts_by_seed.csv，bootstrap汇总在paired_KF_contrast_aggregate.csv。正式扩展前冻结参数/协议并使用未参与开发的种子。','',
      '## 6. 组会可说的话','',
      f'“我接上了方位角和角速度的KF预测，完成5个种子、7种条件的探索性实验。依次检查离线方位预测、同状态输入场和检查点分叉闭环。中心条件下KF方位RMSE约{c.rmseKF.mean():.4g}弧度，hold约{c.rmseHold.mean():.4g}弧度。输入改善与闭环变化分别统计，目前还没有验证真实遮挡和完整三因素交互。”','',
      '## 7. 主张与证据','','|主张|证据|状态|','|---|---|---|',
      '|软件接口和续跑符合本次检查|10项实际测试、原始MAT数据契约复核|仅限已执行检查|',
      '|中心条件预测/输入表现|逐seed RMSE和同状态L2|探索性证据|',
      '|闭环偏离的改变|恢复100步逐seed配对和绝对偏离补充|按实际结果解释|',
      '|协方差已可靠校准|假定R的覆盖/NIS|不能宣称|','|真实遮挡和完整交互下有效|未运行|缺少证据|']
    text='\n'.join(lines)+'\n'; (out/'RESULTS_说明.md').write_text(text,encoding='utf-8'); (ROOT/'RESULTS_说明.md').write_text(text,encoding='utf-8')

def main(result):
    result=Path(result).resolve(); out=result/'analysis'; out.mkdir(exist_ok=True); (out/'figures').mkdir(exist_ok=True)
    audit=json.loads((result/'audit.json').read_text(encoding='utf-8')); assert audit['status']=='COMPLETE'
    manifest=pd.read_csv(result/'RUN_MANIFEST.csv'); closed=pd.read_csv(result/'closed_loop_summary.csv'); same=pd.read_csv(result/'same_state_summary.csv'); diag=pd.read_csv(result/'prediction_diagnostics_summary.csv')
    assert len(manifest)==155 and manifest.run_id.nunique()==155 and set(manifest.seed)==set(range(1,6))
    assert len(closed)==700 and len(same)==175 and len(diag)==105 and len(list((result/'conditions').iterdir()))==35
    assert (diag.numericalFailures==0).all() and (diag.warmup=='full_6200_steps').all()
    for p in manifest.file: assert (result/p).is_file()
    print('Recomputing 700 raw MAT windows and checking 140 branches...',flush=True)
    raw,contracts,checks=recompute(result,closed)
    raw.to_csv(out/'raw_recomputed_summary.csv',index=False); contracts.to_csv(out/'data_contract_checks.csv',index=False); checks.to_csv(out/'independent_recompute_check.csv',index=False)
    print('Recomputing 175 input fields and 105 prediction diagnostics...',flush=True)
    ss,ds=verify_replays(result,same,diag); ss.to_csv(out/'same_state_recompute_check.csv',index=False); ds.to_csv(out/'diagnostic_recompute_check.csv',index=False)
    aggregate(raw,['condition_id','mode','window'],['pairDistanceIncrease','GOLoss','LOLoss','meanAbsolutePairDeviation']).to_csv(out/'closed_loop_seed_aggregate.csv',index=False)
    aggregate(same,['condition_id','mode'],['relativeL2Mean','relativeL2Max']).to_csv(out/'same_state_seed_aggregate.csv',index=False)
    aggregate(diag,['condition_id','R'],['rmseKF','rmseHold','coverage95','recoveryNISmean']).to_csv(out/'prediction_seed_aggregate.csv',index=False)
    contrasts=[]
    for (seed,cond),g in raw[raw.window=='recovery100'].groupby(['seed','condition_id']):
        k=g[g['mode']=='KF-mean'].iloc[0]; h=g[g['mode']=='hold'].iloc[0]; d=g[g['mode']=='drop'].iloc[0]
        for metric in ['pairDistanceIncrease','GOLoss','LOLoss','meanAbsolutePairDeviation']:
            contrasts.append(dict(seed=seed,condition_id=cond,measure=metric,KFminusHold=k[metric]-h[metric],KFminusDrop=k[metric]-d[metric]))
    contrasts=pd.DataFrame(contrasts); contrasts.to_csv(out/'paired_KF_contrasts_by_seed.csv',index=False)
    aggregate(contrasts,['condition_id','measure'],['KFminusHold','KFminusDrop']).to_csv(out/'paired_KF_contrast_aggregate.csv',index=False)
    validation=dict(status='PASS',trajectoryCount=155,uniqueConditions=7,independentSeedCount=5,rawSummaryRows=700,branchDataContractChecks=140,sameStateChecks=175,offlineDiagnosticChecks=105,
        maxRecomputeDifference=float(checks.maxAbsDifference.max()),bootstrap=dict(unit='seed',draws=5000,seed=20261010,sharedDrawIndices=True),exclusions='none')
    plots(result,out,raw,same,diag); write_report(out,raw,same,diag,contrasts)
    (out/'validation.json').write_text(json.dumps(validation,indent=2),encoding='utf-8'); print(json.dumps(validation,indent=2),flush=True)
if __name__=='__main__':
    if len(sys.argv)!=2: raise SystemExit('python analyze_results.py results/pilot_seed1to5_v1')
    main(sys.argv[1])

"""Independent reconstruction from MAT trajectories; never reruns the simulator."""
from __future__ import annotations
import argparse, csv, hashlib, json, math, sys
from pathlib import Path

parser=argparse.ArgumentParser()
parser.add_argument('--run',required=True)
parser.add_argument('--deps',default='')
args=parser.parse_args()
if args.deps: sys.path.insert(0,args.deps)
import numpy as np
from scipy.io import loadmat
from scipy.stats import t as student_t
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

root=Path(args.run).resolve(); out=root/'independent-analysis'
if out.exists(): raise SystemExit('Refusing to overwrite analysis directory')
out.mkdir()
metadata=json.loads((root/'metadata.json').read_text(encoding='utf-8'))
assert metadata['expectedRuns']==24
rows=list(csv.DictReader((root/'summary.csv').open(encoding='utf-8-sig')))
assert len(rows)==24 and len({(r['seed'],r['noiseStd'],r['method']) for r in rows})==24
lookup={(int(r['seed']),float(r['noiseStd']),r['method']):r for r in rows}
recomputed=[]; errors={k:0.0 for k in ['GO','LO','distance','speed','sampleKernel','sampleInputError','sampleAngleError']}
series={}; hashes={}; pairing=[]
for path in sorted(root.glob('seed-*-pair.mat')):
    hashes[path.name]=hashlib.sha256(path.read_bytes()).hexdigest()
    pair=loadmat(path,squeeze_me=True,struct_as_record=False)['pair']
    a,b=pair.raw,pair.kf
    assert np.array_equal(a.initialPosition,b.initialPosition)
    assert np.array_equal(a.initialHeadings,b.initialHeadings)
    assert np.array_equal(a.initialNeuralStateUsed,b.initialNeuralStateUsed)
    assert np.array_equal(a.pilot.selectedHistory,b.pilot.selectedHistory)
    assert np.array_equal(a.pilot.finalNoiseStreamState,b.pilot.finalNoiseStreamState)
    for sa,sb in zip(np.atleast_1d(a.pilot.samples),np.atleast_1d(b.pilot.samples)):
        assert np.array_equal(sa.standardNoise,sb.standardNoise)
    pairing.append({'file':path.name,'initialization':True,'allStepSelection':True,'noiseStreamAndSamples':True})
    for label,r in [('raw',a),('kf-posterior',b)]:
        cfg=r.config; n=int(cfg.NAgents); L=float(cfg.arenaSize); dt=float(cfg.dt)
        burn=int(cfg.burnIn); seed=int(cfg.seed); sd=float(cfg.bearingMeasurementStd)
        assert int(cfg.NSteps)==8000 and burn==6000 and int(cfg.neighborCount)==15
        x=np.asarray(r.trajectory.x); y=np.asarray(r.trajectory.y)
        assert x.shape==(20,8001)
        dx=np.diff(x,axis=1); dy=np.diff(y,axis=1)
        dx-=L*np.round(dx/L); dy-=L*np.round(dy/L)
        vx=dx/dt; vy=dy/dt; speed=np.hypot(vx,vy)
        ux=np.divide(vx,speed,out=np.zeros_like(vx),where=speed>1e-12)
        uy=np.divide(vy,speed,out=np.zeros_like(vy),where=speed>1e-12)
        go=np.hypot(ux.sum(axis=0),uy.sum(axis=0))/n
        errors['GO']=max(errors['GO'],float(np.max(np.abs(go-np.asarray(r.timeSeries.globalOrder)))))
        mean_speed=speed.mean(axis=0)
        errors['speed']=max(errors['speed'],float(np.max(np.abs(mean_speed-np.asarray(r.timeSeries.meanSpeed)))))
        px=x[:,burn+1:].T; py=y[:,burn+1:].T
        ddx=px[:,:,None]-px[:,None,:]; ddy=py[:,:,None]-py[:,None,:]
        ddx-=L*np.round(ddx/L); ddy-=L*np.round(ddy/L)
        distances=np.hypot(ddx,ddy); ii,jj=np.triu_indices(n,1)
        pd=distances[:,ii,jj].mean(axis=1)/L
        neighbors=np.argsort(distances,axis=2)[:,:,:6]
        ux_win=ux[:,burn:].T; uy_win=uy[:,burn:].T
        lx=np.take_along_axis(ux_win[:,None,:],neighbors,axis=2).sum(axis=2)/6
        ly=np.take_along_axis(uy_win[:,None,:],neighbors,axis=2).sum(axis=2)/6
        lo=np.hypot(lx,ly).mean(axis=1)
        errors['LO']=max(errors['LO'],float(np.max(np.abs(lo-np.asarray(r.timeSeries.localOrder)[burn:]))))
        errors['distance']=max(errors['distance'],float(np.max(np.abs(pd-np.asarray(r.timeSeries.meanPairDistanceNorm)[burn:]))))
        for sample in np.atleast_1d(r.pilot.samples):
            step=int(sample.step)-1; ix=np.asarray(sample.linearIndex,dtype=int)-1
            obs=np.asarray(sample.observedBearing); truth=np.asarray(sample.trueBearing)
            used=np.asarray(sample.usedBearing); eps=np.asarray(sample.standardNoise)
            err=np.angle(np.exp(1j*(used-truth)))
            observed_err=np.angle(np.exp(1j*(obs-truth)))
            errors['sampleAngleError']=max(errors['sampleAngleError'],abs(float(np.sqrt(np.mean(err**2)))-float(r.timeSeries.usedBearingRMSE[step])))
            assert np.max(np.abs(np.angle(np.exp(1j*(obs-truth-sd*eps)))))<1e-12
            centers=np.full((n,n),np.nan)
            rec=ix%n; nei=ix//n; centers[rec,nei]=used
            true_centers=np.full((n,n),np.nan); true_centers[rec,nei]=truth
            ids=np.asarray(sample.orderedIds,dtype=int)-1; alpha=np.asarray(sample.alpha)
            fields=[]
            for c in (centers,true_centers):
                field=np.zeros_like(np.asarray(sample.input))
                for i in range(n):
                    delta=np.abs(alpha[i,:,None]-c[i,ids[i]][None,:])
                    delta=np.minimum(delta,2*np.pi-delta)
                    field[:,i]=float(cfg.totalSocialAttraction)/n*np.exp(-.5*delta**2/float(cfg.receptiveFieldWidth)**2).sum(axis=1)
                fields.append(field)
            errors['sampleKernel']=max(errors['sampleKernel'],float(np.max(np.abs(fields[0]-sample.input))),float(np.max(np.abs(fields[1]-sample.truthInput))))
            er=float(np.sqrt(np.mean((fields[0]-fields[1])**2)))
            errors['sampleInputError']=max(errors['sampleInputError'],abs(er-float(r.timeSeries.inputTruthRMS[step])))
            if label=='raw': assert np.array_equal(used,obs)
        row={'seed':seed,'noiseStd':sd,'method':label,'meanGO':float(go[burn:].mean()),
             'meanLO':float(lo.mean()),'meanPairDistanceNorm':float(pd.mean()),
             'meanSpeed':float(mean_speed[burn:].mean()),
             'observedBearingRMSE':float(np.sqrt(np.mean(np.asarray(r.timeSeries.observedBearingRMSE)[burn:]**2))),
             'usedBearingRMSE':float(np.sqrt(np.mean(np.asarray(r.timeSeries.usedBearingRMSE)[burn:]**2))),
             'inputTruthRMS':float(np.sqrt(np.mean(np.asarray(r.timeSeries.inputTruthRMS)[burn:]**2)))}
        archived=lookup[(seed,sd,label)]
        for key in row:
            if key not in ('seed','noiseStd','method'): assert abs(row[key]-float(archived[key]))<1e-9,(key,row[key],archived[key])
        recomputed.append(row); series[(seed,sd,label)]=go
assert len(pairing)==12 and len(recomputed)==24
assert max(errors.values())<1e-9,errors
fields=list(recomputed[0])
with (out/'raw-recomputed.csv').open('w',encoding='utf-8',newline='') as f:
    w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(recomputed)
stats=[]; seeds=list(range(5101,5107)); metrics=fields[3:]
for sd in (.01,.1):
    for metric in metrics:
        aa=np.array([next(r[metric] for r in recomputed if r['seed']==s and r['noiseStd']==sd and r['method']=='raw') for s in seeds])
        bb=np.array([next(r[metric] for r in recomputed if r['seed']==s and r['noiseStd']==sd and r['method']=='kf-posterior') for s in seeds])
        d=bb-aa; se=float(d.std(ddof=1)/math.sqrt(6)); mean=float(d.mean())
        tt=mean/se if se else 0.; p=float(2*student_t.sf(abs(tt),5)) if se else 1.
        margin=float(student_t.ppf(.975,5))*se
        stats.append({'noiseStd':sd,'metric':metric,'n':6,'rawMean':float(aa.mean()),'kfMean':float(bb.mean()),
                      'difference':mean,'ciLow':mean-margin,'ciHigh':mean+margin,'t':tt,'p':p,
                      'holmPrimaryP':None,'positiveSeeds':int((d>0).sum()),'negativeSeeds':int((d<0).sum()),
                      'seedDifferences':d.tolist()})
primary=sorted([s for s in stats if s['metric']=='meanGO'],key=lambda s:s['p'])
last=0.
for i,s in enumerate(primary): last=max(last,min(1.,(2-i)*s['p']));s['holmPrimaryP']=last
(out/'paired-statistics.json').write_text(json.dumps(stats,indent=2),encoding='utf-8')
with (out/'paired-statistics.csv').open('w',encoding='utf-8',newline='') as f:
    sf=[k for k in stats[0] if k!='seedDifferences'];w=csv.DictWriter(f,fieldnames=sf);w.writeheader()
    w.writerows({k:s[k] for k in sf} for s in stats)
plot_metrics=[('meanGO','Global order'),('meanLO','Local order'),('meanPairDistanceNorm','Periodic mean pair distance / L'),
              ('usedBearingRMSE','Used-bearing RMSE (rad)'),('inputTruthRMS','Input RMS vs truth in own geometry'),('meanSpeed','Mean speed (model units)')]
fig,axes=plt.subplots(2,3,figsize=(13,7),layout='constrained')
for ax,(metric,title) in zip(axes.flat,plot_metrics):
    for group,sd in enumerate((.01,.1)):
        xs=np.array([1,2])+3*group
        for seed in seeds:
            vals=[next(r[metric] for r in recomputed if r['seed']==seed and r['noiseStd']==sd and r['method']==m) for m in ('raw','kf-posterior')]
            ax.plot(xs,vals,'o-',color='#6a8499',alpha=.5,linewidth=.8,markersize=3)
        mean=[next(s for s in stats if s['noiseStd']==sd and s['metric']==metric)[k] for k in ('rawMean','kfMean')]
        ax.plot(xs,mean,'D-',color='#c2513c',linewidth=2,markersize=6)
    ax.set_xticks([1,2,4,5],['Raw\n0.01','KF\n0.01','Raw\n0.10','KF\n0.10']);ax.set_title(title,fontsize=10)
    ax.grid(axis='y',alpha=.2)
fig.suptitle('Continuous noisy bearings: paired pilot (6 seeds; k=15; steps 6001-8000)\nThin lines: seeds; red diamonds: means; noise labels are standard deviations in rad',fontsize=12)
fig.savefig(out/'paired-metrics.png',dpi=180);plt.close(fig)
fig,axes=plt.subplots(1,2,figsize=(12,4),layout='constrained')
for ax,sd in zip(axes,(.01,.1)):
    for method,color in [('raw','#376c96'),('kf-posterior','#c2513c')]:
        values=np.stack([series[(s,sd,method)] for s in seeds]); grouped=values.reshape(6,80,100).mean(axis=2)
        ax.plot(np.arange(1,81)*100,grouped.mean(axis=0),label=method,color=color)
    ax.axvspan(6000,8000,color='gray',alpha=.12);ax.set_ylim(0,1);ax.set_title(f'Noise std = {sd:.2f} rad')
    ax.set_xlabel('Model step (100-step block means)');ax.set_ylabel('Actual global order');ax.legend()
fig.savefig(out/'mean-GO-timeseries.png',dpi=180);plt.close(fig)
audit={'status':'PASS_INDEPENDENT_RAW_RECONSTRUCTION','runs':24,'uniqueKeys':24,'pairs':12,'seeds':seeds,
       'maxErrors':errors,'pairing':pairing,'inputHashes':hashes,'evidence':'exploratory, not independent validation',
       'limitations':['6 seeds','fixed q, no tuning','random k15 and original strength only','closed-loop geometry differs',
                      'posterior filtering, no occlusion','secondary intervals unadjusted','diagnostic angle/input checks use archived samples'],
       'fallacyChecksCovered':11}
(out/'audit.json').write_text(json.dumps(audit,indent=2),encoding='utf-8')
md=['# 连续带噪方位观测：直接输入与KF后验均值比较','',
    'Material Passport：academic-research-suite / experiment-agent；2026-10-14。Verification Status：ANALYZED（从原始MAT独立重建指标，不重跑科学仿真）。证据等级：6种子的探索性配对实验。','',
    'N20、k15随机选邻、原始等权增益、周期L1000、M100、8000步、评价6001:8000。无领航、无遮挡。q=1e-4，R=噪声方差。KF先时间预测、用当前带噪观测更新，再把后验均值用于与raw相同宽度和强度的输入核。','',
    '| 噪声std | 指标 | raw均值 | KF均值 | KF−raw | 配对95%区间 | 单项p |','|---|---|---:|---:|---:|---|---:|']
for s in stats:
    if s['metric']=='observedBearingRMSE':continue
    md.append(f"| {s['noiseStd']:.2f} | {s['metric']} | {s['rawMean']:.6g} | {s['kfMean']:.6g} | {s['difference']:+.6g} | [{s['ciLow']:+.6g}, {s['ciHigh']:+.6g}] | {s['p']:.4g} |")
md+=['','主比较族仅为两个噪声条件的实际GO。Holm校正p：'+', '.join(f"std={s['noiseStd']}: {s['holmPrimaryP']:.5g}" for s in primary)+
     '。次要指标未作多重校正。方位误差和输入误差在各自闭环几何下评价，不解释为固定轨迹的单独去噪收益。',
     '', '## 独立复核','',f'24条唯一记录、12对、全部6个种子齐全。最大原数组重构误差：`{json.dumps(errors)}`。',
     '配对初始位置/航向/神经状态、每步选邻顺序、噪声随机流终态和归档样本扰动全部一致。真值仅用于传感器生成和评价，纯控制函数不读取它。','',
     '## 解释限制','',
     '本轮结果限于固定q、连续带噪观测、随机k15与原强度条件；不代表所有KF、其他调参或遮挡预测的结论。平均距离更小不自动代表安全，GO均值不代表持续稳定。没有按结果调整q/R、筛选种子或重新定义终点。',
     '11项统计推断检查：按噪声分层；种子为独立单位；固定纳入全部6种子；不按处理后指标筛选；不宣称成功率；共同初态避免极端样本选择；无幸存者剔除；主族固定且次要比较标注；参数/窗口固定；闭环策略干预仅解释本模型；误差变化不当作已验证中介因果。','',
     '## 运行与回退','',
     '首个软件测试因浮点计算顺序差异失败，恢复旧核公式计算顺序后完全一致；首个运行目录因汇总首行结构体初始化错误中止，保留v1日志和快照，修复后在v2运行，科学参数不变。',
     '原始模型、既有实验结果和四个已有修改源码均未改变。本轮只新增隔离仿真与分析文件；未commit/push/merge。回退只撤销新增实现，结果继续作为历史记录保留。合并前需要审查新增代码，不把此探索性结果当成已确认优势。']
(out/'实验结果与独立复核.md').write_text('\n'.join(md)+'\n',encoding='utf-8')
print(json.dumps({'audit':audit['status'],'maxErrors':errors,'primary':primary,'selectedStats':[s for s in stats if s['metric'] in ('meanLO','meanPairDistanceNorm','usedBearingRMSE')]},ensure_ascii=False,indent=2))

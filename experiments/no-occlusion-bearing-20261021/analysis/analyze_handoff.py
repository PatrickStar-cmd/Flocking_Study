"""Portable independent geometry check and seed-paired GO statistics.

Usage: python analysis/analyze_handoff.py PATH_TO_NEW_RUN
Requires NumPy, SciPy, Matplotlib; refuses existing analysis output.
"""
from pathlib import Path
import argparse, csv, json, math, hashlib
import numpy as np
from scipy.io import loadmat
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

METHODS=['oracle','raw','selected-kf','imm']
CONTRASTS=[('oracle','raw'),('selected-kf','raw'),('imm','raw'),
           ('imm','selected-kf'),('oracle','selected-kf'),('oracle','imm')]

def main():
    parser=argparse.ArgumentParser(__doc__)
    parser.add_argument('run',type=Path)
    args=parser.parse_args(); run=args.run.resolve(); out=run/'independent-analysis'
    assert not out.exists(), 'Existing analysis is protected'
    completion=json.loads((run/'completion.json').read_text(encoding='utf-8'))
    assert completion['status']=='COMPLETED'
    with (run/'summary.csv').open(encoding='utf-8-sig',newline='') as f:
        rows=list(csv.DictReader(f))
    assert len(rows)==completion['summaryRows']
    keys=[(r['seed'],r['distribution'],r['measurementStd'],r['method']) for r in rows]
    assert len(keys)==len(set(keys)), 'Duplicate scientific keys'
    verified={}; errors=[]
    for row in rows:
        path=run/row['matFile']
        if row['matFile'] not in verified:
            r=loadmat(path,simplify_cells=True)['result']; cfg=r['config']
            x=np.asarray(r['trajectory']['x']); y=np.asarray(r['trajectory']['y'])
            assert x.shape[1]==int(cfg['NSteps'])+1, 'recordStride must equal one'
            L=float(cfg['arenaSize']); dt=float(cfg['dt']); burn=int(cfg['burnIn'])
            dx=(np.diff(x,axis=1)+L/2)%L-L/2
            dy=(np.diff(y,axis=1)+L/2)%L-L/2
            speed=np.hypot(dx,dy)/dt; moving=speed>1e-12
            ux=np.divide(dx/dt,speed,out=np.zeros_like(speed),where=moving)
            uy=np.divide(dy/dt,speed,out=np.zeros_like(speed),where=moving)
            go=np.hypot(ux.sum(axis=0),uy.sum(axis=0))/x.shape[0]
            # Independent distance computation at post-movement positions.
            xx=x[:,burn+1:]; yy=y[:,burn+1:]
            px=(xx[:,None,:]-xx[None,:,:]+L/2)%L-L/2
            py=(yy[:,None,:]-yy[None,:,:]+L/2)%L-L/2
            d=np.hypot(px,py); ii,jj=np.triu_indices(x.shape[0],1)
            values=dict(meanGO=float(go[burn:].mean()),meanPairDistanceNorm=float(d[ii,jj,:].mean()/L),
                        meanSpeed=float(speed[:,burn:].mean()))
            values['sha256']=hashlib.sha256(path.read_bytes()).hexdigest()
            verified[row['matFile']]=values
            del r,x,y,dx,dy,speed,px,py,d
        values=verified[row['matFile']]
        for metric in ['meanGO','meanPairDistanceNorm','meanSpeed']:
            error=abs(values[metric]-float(row[metric])); errors.append(error)
            assert error<2e-9,(row['matFile'],metric,error)
    assert len(verified)==completion['uniqueRuns']
    out.mkdir()
    groups={}
    for row in rows:
        groups.setdefault((row['distribution'],row['measurementStd']),{}).setdefault(int(row['seed']),{})[row['method']]=float(row['meanGO'])
    comparisons=[]
    for (dist,sd),seeds in sorted(groups.items()):
        assert all(set(v)==set(METHODS) for v in seeds.values())
        n=len(seeds); result=[]
        for a,b in CONTRASTS:
            d=np.array([v[a]-v[b] for _,v in sorted(seeds.items())]); mean=float(d.mean())
            se=float(d.std(ddof=1)/math.sqrt(n)) if n>1 else float('nan')
            half=float(stats.t.ppf(1-.05/(2*6),n-1)*se) if n>1 else float('nan')
            p=float(stats.ttest_1samp(d,0).pvalue) if n>1 and se>0 else (1.0 if n>1 and mean==0 else float('nan'))
            result.append(dict(distribution=dist,measurementStd=sd,contrast=f'{a}-{b}',n=n,mean=mean,
                               simultaneous95Low=mean-half,simultaneous95High=mean+half,p=p))
        if n>1:
            order=sorted(range(6),key=lambda i:result[i]['p']); previous=0
            for rank,i in enumerate(order):
                previous=max(previous,min(1,(6-rank)*result[i]['p'])); result[i]['holm6p']=previous
        else:
            for item in result: item['holm6p']=float('nan')
        comparisons.extend(result)
        fig,ax=plt.subplots(figsize=(7,4))
        means=[np.mean([v[m] for v in seeds.values()]) for m in METHODS]
        ax.bar(METHODS,means); ax.set(ylim=(0,1),ylabel='Mean GO',title=f'{dist}, sigma={sd}, n={n} paired seeds')
        fig.tight_layout(); fig.savefig(out/f'means-{dist}-{sd}.png',dpi=160); plt.close(fig)
    with (out/'contrasts.csv').open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(comparisons[0])); w.writeheader(); w.writerows(comparisons)
    audit=dict(status='PASS',evidence='independent geometry/software check; inference depends on run profile',
               uniqueMATFiles=len(verified),conditionRows=len(rows),maxIndependentMetricError=max(errors),
               hashesAndRecomputedMetrics=verified,
               inference='six GO comparisons per distribution/sigma, Bonferroni simultaneous95 intervals and Holm6; no inference at n1',
               warning='pilot/validation64 profiles reuse archived seeds; smoke is not scientific evidence')
    (out/'audit.json').write_text(json.dumps(audit,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in audit.items() if k!='hashesAndRecomputedMetrics'},indent=2))

if __name__=='__main__': main()

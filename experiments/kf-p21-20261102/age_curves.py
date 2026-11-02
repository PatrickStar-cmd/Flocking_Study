"""Complete raw branch age curves with the same seed/window weighting."""
import json,sys
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.io import loadmat
from replay_q_grid import group_masks


def run(stage):
    path=stage/'analysis/seed_age_curves.csv';assert not path.exists()
    inputs=json.loads((stage/'input_manifest.json').read_text(encoding='utf-8'))
    manifest=json.loads((stage/'raw/replay/manifest.json').read_text(encoding='utf-8'))
    cutfile=stage/('GEOMETRY_FROZEN.json' if 'development' in stage.name else 'GEOMETRY_APPLIED.json')
    cuts=json.loads(cutfile.read_text(encoding='utf-8'))['rate_cutpoints']
    methods=sorted(set(e['method'] for e in manifest['files']))
    totals={}
    for entry in inputs['inputs']:
        p=stage/entry['source'];stem=p.stem;source=loadmat(p,simplify_cells=True);seed=int(entry['seed'])
        for D in [20,120]:
            with np.load(stage/'labels'/f'{stem}_D{D}.npz') as z:
                labels=z['labels'];rates=z['rates'];starts=z['starts']
            ix=starts[None,:]-1+np.arange(D)[:,None];truth=source['true_bearing'][ix]
            for method in methods:
                with np.load(stage/'raw/replay'/f'{stem}_{method}_fresh.npz') as fresh:
                    with np.load(stage/'raw/replay'/f'{stem}_{method}_D{D}.npz') as branch:
                        ref=np.ascontiguousarray(fresh['theta'][ix,0]).view(np.uint64)
                        theta=np.bitwise_xor(branch['theta_xor'][:D],ref).view(np.float64)
                        cp=np.ascontiguousarray(fresh['P'][ix,0]).view(np.uint64)
                        P=np.bitwise_xor(branch['P_xor'][:D],cp).view(np.float64)[...,0,0]
                        age=branch['age'][:D]
                e=np.arctan2(np.sin(theta-truth),np.cos(theta-truth))
                for c in range(len(starts)):
                    masks=group_masks(labels[c],rates[c],cuts)
                    for geometry in ['all','clean-slow','clean-fast','G-switch','G-close']:
                        mask=masks[geometry]
                        if not mask.any():continue
                        key=(seed,D,method,geometry)
                        if key not in totals:totals[key]={'rmse':np.zeros(D),'coverage':np.zeros(D),'cutoff':np.zeros(D),'windows':0}
                        t=totals[key];t['windows']+=1
                        t['rmse']+=np.sqrt(np.mean(e[:,c,mask]**2,axis=1))
                        if method=='Hold':t['coverage'][:]=np.nan;reject=age>120
                        else:
                            t['coverage']+=(np.abs(e[:,c,mask])<=1.96*np.sqrt(P[:,c,None])).mean(axis=1)
                            reject=(age>120)|(P[:,c]>1)
                        t['cutoff']+=reject
    rows=[]
    for (seed,D,method,geometry),t in totals.items():
        for step in range(D):rows.append({'seed':seed,'D':D,'method':method,'geometry':geometry,'age_step':step+1,'age_time':(step+1)*.3,
            'mean_window_RMSE':t['rmse'][step]/t['windows'],'mean_window_coverage':t['coverage'][step]/t['windows'],
            'mean_window_cutoff_fraction':t['cutoff'][step]/t['windows'],'windows':t['windows']})
    frame=pd.DataFrame(rows);frame.to_csv(path,index=False)
    fig,axes=plt.subplots(2,2,figsize=(11,8))
    for row,geometry in enumerate(['clean-slow','clean-fast']):
        for col,D in enumerate([20,120]):
            ax=axes[row,col]
            for method in methods:
                f=frame[(frame.geometry==geometry)&(frame.D==D)&(frame.method==method)]
                if not len(f):continue
                grouped=f.groupby('age_step').mean_window_RMSE
                mean=grouped.mean();sd=grouped.std()
                ax.plot(mean.index,mean.values,label=method,lw=1.2)
                ax.fill_between(mean.index,np.maximum(0,mean-sd),mean+sd,alpha=.06)
            ax.set_title(f'{geometry}, D={D}');ax.set_xlabel('Missing age (steps; dt=.3 model time)');ax.set_ylabel('Mean-window circular RMSE (rad)')
            if row==0 and col==0:ax.legend(fontsize=7)
    title='Development' if 'development' in stage.name else 'Independent validation'
    fig.suptitle(f'{title}: posterior clean labels; {len(inputs["seeds"])} physical seeds\nLines: seed means; bands: +/-1 seed SD; no gate-based exclusions')
    fig.tight_layout(rect=(0,0,1,.93));fig.savefig(stage/'analysis/V2_age_curves.png',dpi=170);plt.close(fig)
    print('PASS age curves',stage,len(rows),flush=True)


if __name__=='__main__':run(Path(sys.argv[1]))

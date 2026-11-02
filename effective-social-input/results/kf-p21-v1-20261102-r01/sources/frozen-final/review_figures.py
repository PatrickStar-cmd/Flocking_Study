"""Final G0 tails use every seed/sample; repair G1 canvas clipping separately."""
import json,sys
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.io import loadmat
from geometry import geometry,classify_window,fresh_labels


def curve(values,tail=False):
    x=np.sort(np.concatenate(values));n=len(x)
    # Exact empirical probabilities at a thinned set of sorted full-data ranks.
    ranks=np.unique(np.r_[np.linspace(0,n-1,min(1500,n)).astype(int),
        n-1-np.geomspace(1,max(1,n),min(700,n)).astype(int),0,n-1])
    ranks=ranks[(ranks>=0)&(ranks<n)]
    probability=(ranks+1)/n
    if tail:probability=np.maximum(1/n,1-probability)
    return x[ranks],probability,n


def main(parent):
    dest=parent/'review/figures';assert not dest.exists();dest.mkdir()
    repo=Path(__file__).resolve().parent.parent.parent
    old=repo/'effective-social-input/results/kf-p2-v1-20261102-reviewed'
    d=json.loads((old/'design.json').read_text(encoding='utf-8'))
    errors={k:[] for k in ['G-clean','G-switch','G-close','G-ambiguous']}
    innovations={k:[] for k in ['G-clean','switch_event','switch_near2','G-close','G-ambiguous']}
    for case in d['cases']:
        g=geometry(loadmat(old/case['source'],simplify_cells=True));fresh=fresh_labels(g)
        for D in [20,120]:
            label,_=classify_window(g,case['t0'],D)
            with np.load(old/'raw/replay'/f"case_{case['case_id']}_D_{D}.npz") as z:
                e=z['theta_KF_CV'][200:200+D]-z['truth'][200:200+D]
                e=np.abs(np.arctan2(np.sin(e),np.cos(e)))
                for k in errors:
                    if (label==k).any():errors[k].append(e[:,label==k].ravel())
                fl=fresh[z['window_steps'].astype(int)-1];n=z['nis_KF_CV']
                for k in innovations:
                    a=n[fl==k];a=a[np.isfinite(a)]
                    if len(a):innovations[k].append(a)
    fig,axes=plt.subplots(1,2,figsize=(12,4.7))
    for k,parts in errors.items():
        if parts:
            x,y,n=curve(parts);axes[0].plot(x,y,label=f'{k}, n={n:,}')
    for k,parts in innovations.items():
        if parts:
            x,y,n=curve(parts,True);axes[1].loglog(np.maximum(x,1e-12),y,label=f'{k}, n={n:,}')
    axes[0].set_xlabel('Absolute circular occlusion error (rad)');axes[0].set_ylabel('Empirical CDF')
    axes[1].set_xlabel('Fresh NIS (log)');axes[1].set_ylabel('Empirical upper-tail probability (log)')
    for ax in axes:ax.grid(alpha=.2);ax.legend(fontsize=7)
    fig.suptitle('G0: all 5 reused physical seeds and all source samples\nPosterior geometry labels; overlapping windows descriptive; curves thinned only after full sorting',fontsize=10)
    fig.tight_layout(rect=(0,0,1,.88));fig.savefig(dest/'G0_full_tails.png',dpi=190,bbox_inches='tight');plt.close(fig)
    f=pd.read_csv(parent/'g1/metrics_by_age.csv');f=f[f.D==120]
    fig,axes=plt.subplots(1,5,figsize=(16,4))
    for ax,config in zip(axes,['A1','B2','C1','D1','E1']):
        t=f[f.config_id==config].groupby('age_time')[['err_Hold','err_KF_CV','err_KF_GM']].mean()
        for name in t:ax.plot(t.index,t[name],label=name.replace('err_',''),lw=1.3)
        ax.set_title(config);ax.set_xlabel('Age (model time)');ax.grid(alpha=.2)
    axes[0].set_ylabel('Mean absolute circular error (rad)');axes[0].legend(fontsize=7)
    fig.suptitle('G1 controlled mechanisms: 20 paired noise streams, D120\nSynthetic/software evidence only; not pooled with physical validation seeds',fontsize=10)
    fig.tight_layout(rect=(0,0,1,.84));fig.savefig(dest/'G1_mechanism_age.png',dpi=190,bbox_inches='tight');plt.close(fig)
    print('PASS final G0 full-data and G1 canvas QA figures')


if __name__=='__main__':main(Path(sys.argv[1]))

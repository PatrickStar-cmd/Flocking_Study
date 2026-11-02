"""Reviewer supplements: stratified contrasts, error distributions, fresh NIS."""
import json
from pathlib import Path
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def run(root):
    dest = root/'review/figures'
    assert not dest.exists()
    dest.mkdir()
    d = json.loads((root/'design.json').read_text(encoding='utf-8'))
    groups = pd.read_csv(root/'analysis/grouped.csv')
    colors = {'Hold':'#777777','KF_CV':'#0072b2','KF_GM':'#d55e00'}
    labels = {'Hold':'Hold','KF_CV':'KF-CV','KF_GM':'KF-GM'}
    fig, axes = plt.subplots(3,2,figsize=(11,10))
    for row, stratum in enumerate('STV'):
        for col, duration in enumerate([20,120]):
            ax = axes[row,col]
            subset = groups[(groups.stratum==stratum)&(groups.D==duration)&groups.metric.str.startswith('diff_')]
            for index, metric in enumerate(['diff_KF_CV_minus_Hold','diff_KF_GM_minus_Hold','diff_KF_GM_minus_KF_CV']):
                r = subset[subset.metric==metric].iloc[0]
                ax.errorbar(r['mean'],index,xerr=[[r['mean']-r.ci_lower],[r.ci_upper-r['mean']]],fmt='o',capsize=4,color=['#0072b2','#d55e00','#009e73'][index])
            ax.axvline(0,color='black',lw=.7)
            ax.set_yticks([0,1,2],['CV - Hold','GM - Hold','GM - CV'])
            ax.set_title(f'{stratum}, D={duration}; n=5 seeds')
            ax.set_xlabel('Unconditional circular RMSE difference (rad)')
            ax.grid(alpha=.25)
    fig.suptitle('Paired seed means and percentile bootstrap 95% intervals\nExploratory comparisons; intervals are not simultaneous')
    fig.tight_layout(rect=(0,0,1,.94)); fig.savefig(dest/'paired_intervals.png',dpi=180);plt.close(fig)
    fig, axes = plt.subplots(3,2,figsize=(11,10))
    for row, stratum in enumerate('STV'):
        for col, duration in enumerate([20,120]):
            ax=axes[row,col]
            cases=[c for c in d['cases'] if c['stratum']==stratum]
            errors={k:[] for k in labels}
            for case in cases:
                with np.load(root/'raw/replay'/f"case_{case['case_id']}_D_{duration}.npz") as raw:
                    for key in labels:
                        errors[key].append(raw['error_'+key][200:200+duration].ravel())
            for key in labels:
                ax.hist(np.concatenate(errors[key]),bins=np.linspace(-np.pi,np.pi,100),density=True,histtype='step',lw=1.3,color=colors[key],label=labels[key])
            ax.set_yscale('log'); ax.set_xlim(-np.pi,np.pi)
            ax.set_title(f'{stratum}, D={duration}; 5 reused seeds')
            ax.set_xlabel('Signed circular error (rad)'); ax.set_ylabel('Descriptive density (log)')
            if row==0 and col==0: ax.legend()
    fig.suptitle('Circular error distributions, including tails\nPooled edge-times are descriptive; no independent-sample inference')
    fig.tight_layout(rect=(0,0,1,.94));fig.savefig(dest/'circular_error_distributions.png',dpi=180);plt.close(fig)
    fig,axes=plt.subplots(3,2,figsize=(11,10))
    rows=[]
    for row,stratum in enumerate('STV'):
        cases=[c for c in d['cases'] if c['stratum']==stratum]
        for col,duration in enumerate([20,120]):
            ax=axes[row,col]
            for key in ['KF_CV','KF_GM']:
                seed_arrays=[]
                for case in cases:
                    with np.load(root/'raw/replay'/f"case_{case['case_id']}_D_{duration}.npz") as raw:
                        nis=raw['nis_'+key]
                        mean=np.full(len(nis),np.nan)
                        finite=np.isfinite(nis).any(axis=1)
                        mean[finite]=np.mean(nis[finite],axis=1)
                        seed_arrays.append(mean)
                arr=np.stack(seed_arrays); means=np.full(arr.shape[1],np.nan)
                valid=np.isfinite(arr).any(axis=0);means[valid]=np.nanmean(arr[:,valid],axis=0)
                steps=np.arange(len(means))-200
                ax.plot(steps*.3,means,color=colors[key],label=labels[key],lw=1)
                for t,value in zip(steps,means):
                    rows.append({'stratum':stratum,'D':duration,'method':labels[key],'relative_step':t,'mean_fresh_NIS':value,'n_seeds':5})
            ax.axvspan(0,duration*.3,color='gray',alpha=.12)
            ax.axhline(1,color='black',ls='--',lw=.8)
            ax.set_yscale('log');ax.set_title(f'{stratum}, D={duration}; 5-seed means')
            ax.set_xlabel('Time from occlusion onset (model units)');ax.set_ylabel('Fresh-observation NIS (log)')
            if row==0 and col==0:ax.legend()
    fig.suptitle('Innovation consistency before and after the missing interval\nMissing frames have no NIS; line at 1 is a Gaussian-model reference')
    fig.tight_layout(rect=(0,0,1,.94));fig.savefig(dest/'fresh_NIS_time.png',dpi=180);plt.close(fig)
    pd.DataFrame(rows).to_csv(root/'review/fresh_NIS_time.csv',index=False)
    print('Wrote 3 reviewer figures and fresh_NIS_time.csv')


if __name__=='__main__':
    run(Path(sys.argv[1]))

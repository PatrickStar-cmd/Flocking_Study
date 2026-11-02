"""Exploratory development paired intervals; explicitly not validation evidence."""
import sys
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from validate_frozen import describe


def main(stage):
    path=stage/'analysis/development_paired_exploratory.csv';assert not path.exists()
    data=pd.read_csv(stage/'analysis/seed_metrics.csv');data=data[data.scenario=='combined']
    rows=[];fig,axes=plt.subplots(1,2,figsize=(11,4.5))
    for ax,group,D in zip(axes,['clean-fast','clean-slow'],[120,20]):
        f=data[(data.geometry==group)&(data.D==D)]
        hold=f[f.method=='Hold'].set_index('seed')
        for i,method in enumerate(['CV_q0.0001','CV_q1e-05']):
            model=f[f.method==method].set_index('seed')
            ids=hold.index.intersection(model.index);delta=model.loc[ids].rmse-hold.loc[ids].rmse
            stats=describe(delta.to_numpy());rows.append({'evidence':'exploratory-development','geometry':group,'D':D,'method':method,**stats})
            ax.scatter(np.full(len(delta),i)+np.linspace(-.10,.10,len(delta)),delta,s=20)
            ax.errorbar(i,stats['mean'],yerr=[[stats['mean']-stats['ci_lower']],[stats['ci_upper']-stats['mean']]],fmt='ks',capsize=5)
        ax.axhline(0,color='black',ls='--',lw=.8);ax.set_xticks([0,1],['default CV q=1e-4','candidate q=1e-5'])
        ax.set_title(f'{group}, D={D}; n=12 physical seeds');ax.set_ylabel('CV - Hold window-mean RMSE difference (rad)');ax.grid(alpha=.2)
    fig.suptitle('DEVELOPMENT ONLY: paired seed points; mean and seed-bootstrap 95% intervals\nCandidate q=1e-5 did not pass the frozen calibration screen; no independent validation',fontsize=10)
    fig.tight_layout(rect=(0,0,1,.86));fig.savefig(stage/'analysis/development_paired_exploratory.png',dpi=180,bbox_inches='tight');plt.close(fig)
    pd.DataFrame(rows).to_csv(path,index=False);print(pd.DataFrame(rows).to_string(index=False))


if __name__=='__main__':main(Path(sys.argv[1]))

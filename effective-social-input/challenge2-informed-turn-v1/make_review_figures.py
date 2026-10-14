"""Additional review figures, exclusively new files, leaving prior figures intact."""
import argparse
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from analyze_results import read, COLORS, LABELS, shade, finish
from run_experiment import SEEDS, MODES


def make(stage):
    out=stage/'review-figures';out.mkdir(exist_ok=False)
    meta=json.loads((stage/'metadata.json').read_text(encoding='utf-8'));c=meta['config']
    T=18.;start=c['warm_steps']+c['lead_steps']
    raw={m:[read(stage/'runs'/f'seed-{s}_turn_{m}_T-{0 if m=="fresh" else T:g}.csv') for s in SEEDS] for m in MODES}
    metrics=[('GO','GO'),('LO','LO'),('fixed5_distance','Frozen 5-neighbor distance'),
             ('follower_turn_error_deg','Follower command error (deg)'),
             ('S','Empirical social input S'),('S_abs_error','Same-state social-input error |S - S oracle|')]
    fig,axes=plt.subplots(2,3,figsize=(15,8))
    for ax,(metric,label) in zip(axes.flat,metrics):
        visible=[]
        for mode in MODES:
            ds=raw[mode];t=ds[0]['time'];v=np.array([d[metric] for d in ds])
            mask=(t>=-4)&(t<=T+30)
            visible.extend(v[:,mask].ravel().tolist())
            ax.plot(t,v.mean(0),color=COLORS[mode],lw=1.5,label=LABELS[mode])
            ax.fill_between(t,v.min(0),v.max(0),color=COLORS[mode],alpha=.09)
        valid=np.array(visible);valid=valid[np.isfinite(valid)]
        lo=valid.min();hi=valid.max();pad=max((hi-lo)*.07,1e-5)
        ax.set_ylim(lo-pad,hi+pad);ax.set_xlim(-4,T+30)
        ax.set_ylabel(label);shade(ax,T)
    axes[0,0].legend(fontsize=8,ncol=2)
    fig.suptitle('N=20; T=18; 3 development seeds; 51 leader-observation edges masked\nMean lines; seed-range shading (not confidence intervals); no smoothing')
    fig.tight_layout();finish(fig,out/'occlusion-response-T-18.png')
    # Full-horizon trajectory view for every seed; all methods share axes.
    modes=['fresh','drop','hold','kf-mean','kf-uncertainty']
    for seed in SEEDS:
        fig,axes=plt.subplots(1,5,figsize=(20,5))
        allpoints=[]
        for ax,mode in zip(axes,modes):
            z=np.load(stage/'runs'/f'seed-{seed}_turn_{mode}_T-{0 if mode=="fresh" else T:g}.npz')
            xy=z['xy'];lo=c['warm_steps'];allpoints.append(xy[lo:].reshape(-1,2))
            for j in range(20):
                ax.plot(xy[lo:,j,0],xy[lo:,j,1],color='#009e73' if j<3 else '#687d9a',lw=1.2 if j<3 else .65,alpha=.9 if j<3 else .65)
            for ix,marker,label in [(lo,'o','Turn onset'),(start,'s','Blackout onset'),(start+60,'x','Restored')]:
                ax.scatter(xy[ix,:,0],xy[ix,:,1],s=13,marker=marker,c=['#009e73']*3+['#687d9a']*17,label=label)
            ax.set_title(LABELS[mode]);ax.set_xlabel('x (model units)');ax.set_aspect('equal',adjustable='box');ax.grid(alpha=.2)
        points=np.concatenate(allpoints);lower=points.min(0)-5;upper=points.max(0)+5
        for ax in axes:ax.set_xlim(lower[0],upper[0]);ax.set_ylim(lower[1],upper[1])
        axes[0].set_ylabel('y (model units)');axes[0].legend(fontsize=7)
        fig.suptitle(f'Full post-turn trajectories; development seed {seed}; T=18\nAll20 shown; green=3 leaders; blue=17 followers; common spatial limits')
        fig.tight_layout();finish(fig,out/f'full-trajectories-seed-{seed}-T-18.png')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--stage',type=Path,required=True)
    make(p.parse_args().stage)

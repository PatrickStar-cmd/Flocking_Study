"""Post hoc tail-risk inspection; descriptive thresholds do not replace gates."""
import argparse
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from run_experiment import SEEDS,MODES,write_csv
from analyze_results import LABELS,finish,shade


def inspect(stage):
    out=stage/'long-horizon-diagnostics';out.mkdir(exist_ok=False)
    c=json.loads((stage/'metadata.json').read_text())['config']
    start=c['warm_steps']+c['lead_steps'];rows=[];groups=[]
    for T in [6.,10.2,18.]:
        for seed in SEEDS:
            for mode in MODES:
                name=f'seed-{seed}_turn_{mode}_T-{0 if mode=="fresh" else T:g}'
                z=np.load(stage/'runs'/(name+'.npz'))
                e=np.abs((z['headings'][-30:,3:]-np.pi/4+np.pi)%(2*np.pi)-np.pi)*180/np.pi
                per=e.mean(0)
                for i,value in enumerate(per):
                    rows.append(dict(seed=seed,T=T,mode=mode,agent_id_1based=i+4,
                                     final30_mean_error_deg=value,over30deg=int(value>30),over90deg=int(value>90)))
                groups.append(dict(seed=seed,T=T,mode=mode,final30_group_error_deg=per.mean(),
                                   worst_agent_id_1based=int(per.argmax()+4),worst_agent_error_deg=per.max(),
                                   n_over30deg=int((per>30).sum()),n_over90deg=int((per>90).sum())))
    write_csv(out/'per-agent-final-errors.csv',rows);write_csv(out/'tail-summary.csv',groups)
    for seed in SEEDS:
        fig,axes=plt.subplots(1,5,figsize=(20,5),sharey=True)
        for ax,mode in zip(axes,['fresh','drop','hold','kf-mean','kf-uncertainty']):
            name=f'seed-{seed}_turn_{mode}_T-{0 if mode=="fresh" else 18:g}'
            z=np.load(stage/'runs'/(name+'.npz'))
            step=np.arange(800)+1
            command=np.minimum(np.pi/4,.025*np.maximum(0,(step-c['warm_steps'])*c['dt']))
            e=np.abs((z['headings'][1:,3:]-command[:,None]+np.pi)%(2*np.pi)-np.pi)*180/np.pi
            worst=e[-30:].mean(0).argmax();time=(step-start)*c['dt']
            for i in range(17):
                ax.plot(time,e[:,i],color='#d55e00' if i==worst else '#687d9a',alpha=.95 if i==worst else .35,lw=1.3 if i==worst else .6)
            ax.set_title(f'{LABELS[mode]}\nWorst endpoint agent: {worst+4}')
            ax.set_xlim(-4,time[-1]);ax.set_ylim(0,180);shade(ax,18)
        axes[0].set_ylabel('Individual follower command error (deg)')
        fig.suptitle(f'Post hoc individual tail inspection; development seed {seed}; T=18\nAll17 follower errors; orange=largest final-30-step mean error in each method; no exclusions')
        fig.tight_layout();finish(fig,out/f'agent-errors-seed-{seed}-T-18.png')
    print('Tail inspection saved. Descriptive 30/90deg counts do not change predeclared task gates.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--stage',type=Path,required=True)
    inspect(p.parse_args().stage)

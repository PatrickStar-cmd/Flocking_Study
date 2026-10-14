"""Rebuild unsmoothed development figures and summaries from saved raw files."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from run_experiment import MODES, SEEDS, write_csv

COLORS=dict(zip(MODES,['#202020','#d55e00','#cc79a7','#888888','#0072b2','#009e73']))
LABELS={'fresh':'Fresh','drop':'Drop','drop-normalized':'Drop normalized',
        'hold':'Hold','kf-mean':'KF mean','kf-uncertainty':'KF uncertainty'}


def read(path):
    with path.open(encoding='utf-8-sig',newline='') as f:
        rows=list(csv.DictReader(f))
    return {k:np.array([float(r[k]) for r in rows]) for k in rows[0]}


def shade(ax,T):
    ax.axvspan(0,T,color='#ffdca8',alpha=.3,zorder=-10)
    ax.axvline(0,color='#777',lw=.7)
    ax.axvline(T,color='#777',lw=.7,ls=':')
    ax.set_xlabel('Time from blackout onset (model units)')
    ax.grid(alpha=.18)


def finish(fig,path):
    fig.savefig(path,dpi=160,bbox_inches='tight')
    plt.close(fig)


def analyze(stage):
    meta=json.loads((stage/'metadata.json').read_text(encoding='utf-8'))
    c=meta['config']; dt=c['dt']; start=c['warm_steps']+c['lead_steps']
    with (stage/'manifest.csv').open(encoding='utf-8') as f:
        manifest=list(csv.DictReader(f))
    out=stage/'analysis'
    out.mkdir(exist_ok=False)
    raw={p.stem:read(p) for p in (stage/'runs').glob('*.csv')}
    def key(seed,mode,T,turn=1):
        return f'seed-{seed}_{"turn" if turn else "straight"}_{mode}_T-{0 if mode=="fresh" else T:g}'
    Ts=sorted({float(x['T']) for x in manifest})
    summaries=[]; gates=[]
    available=['fresh']+[m for m in MODES[1:] if any(r['mode']==m for r in manifest)]
    for T in Ts:
        event=slice(start,start+round(T/dt));post=slice(start+round(T/dt),start+round(T/dt)+100)
        for seed in SEEDS:
            fresh=raw[key(seed,'fresh',T)]
            drop=raw[key(seed,'drop',T)]
            sf=raw[key(seed,'fresh',T,0)];sd=raw[key(seed,'drop',T,0)]
            increase=np.mean(drop['follower_turn_error_deg'][event]-fresh['follower_turn_error_deg'][event])
            straight_increase=np.mean(sd['follower_turn_error_deg'][event]-sf['follower_turn_error_deg'][event])
            leader_error=fresh['leader_turn_error_deg'][-30:].mean()
            follower_error=fresh['follower_turn_error_deg'][-30:].mean()
            follower_turn=fresh['follower_turn_angle_deg'][-30:].mean()
            task=leader_error<=10 and follower_error<=15 and follower_turn>=30
            gate=dict(seed=seed,T=T,fresh_final_leader_error_deg=leader_error,
                      fresh_final_follower_error_deg=follower_error,
                      fresh_final_follower_turn_deg=follower_turn,
                      fresh_task_pass=int(task),drop_event_error_increase_deg=increase,
                      turn_specific_error_increase_deg=increase-straight_increase,
                      effective_mask_pass=int(task and increase>=3),
                      pre_event_follower_error_deg=float(fresh['follower_turn_error_deg'][start-1]),
                      event_hidden_bearing_change_deg=float(fresh['hidden_bearing_change_deg'][event].mean()))
            gates.append(gate)
            for mode in available:
                d=raw[key(seed,mode,T)]
                row=dict(seed=seed,T=T,mode=mode,
                         event_turn_error_deg=d['follower_turn_error_deg'][event].mean(),
                         event_turn_error_increase_deg=np.mean(d['follower_turn_error_deg'][event]-fresh['follower_turn_error_deg'][event]),
                         event_GO=d['GO'][event].mean(),event_LO=d['LO'][event].mean(),
                         event_GO_deviation=np.mean(np.abs(d['GO'][event]-fresh['GO'][event])),
                         event_LO_deviation=np.mean(np.abs(d['LO'][event]-fresh['LO'][event])),
                         event_fixed5_distance=d['fixed5_distance'][event].mean(),
                         post_fixed5_distance=d['fixed5_distance'][post].mean(),
                         post_fixed5_distance_increase=np.mean(d['fixed5_distance'][post]-fresh['fixed5_distance'][post]),
                         post_fixed5_distance_abs_deviation=np.mean(np.abs(d['fixed5_distance'][post]-fresh['fixed5_distance'][post])),
                         event_S=d['S'][event].mean(), event_A=d['A'][event].mean(),event_C=d['C'][event].mean(),
                         event_S_abs_error=d['S_abs_error'][event].mean(),
                         final_turn_error_deg=d['follower_turn_error_deg'][-30:].mean(),
                         event_mean_speed=d['mean_speed'][event].mean(),
                         event_bearing_error_deg=float(np.nanmean(d['bearing_error_deg'][event])) if mode not in ['fresh','drop','drop-normalized'] else np.nan,
                         event_shadow_coverage=float(np.nanmean(d['shadow_kf_coverage'][event])) if mode!='fresh' else np.nan,
                         event_active_hidden_edges=d['hidden_active_edges'][event].mean())
                summaries.append(row)
        # Mean and seed range; do not present n=3 seed ranges as confidence bands.
        metrics=[('GO','GO'),('LO','LO'),('fixed5_distance','Frozen 5-neighbor distance (model units)'),
                 ('follower_turn_error_deg','Follower command error (deg)'),
                 ('S','Empirical social input S (followers)'),('leader_follower_distance','Leader-follower distance (model units)')]
        fig,axes=plt.subplots(2,3,figsize=(16,8))
        for ax,(metric,label) in zip(axes.flat,metrics):
            for mode in available:
                ds=[raw[key(s,mode,T)] for s in SEEDS];time=ds[0]['time']
                values=np.array([d[metric] for d in ds])
                ax.plot(time,values.mean(0),color=COLORS[mode],lw=1.5,label=LABELS[mode])
                ax.fill_between(time,values.min(0),values.max(0),color=COLORS[mode],alpha=.09)
            ax.set_ylabel(label);shade(ax,T)
            ax.set_xlim(-c['lead_steps']*dt-3,time[-1])
        axes[0,0].legend(fontsize=8,ncol=2)
        fig.suptitle(f'N=20; 3 development seeds; T={T:g}; 51/300 directed edges masked\nLines: mean; shading: seed range; orange: missing-observation window')
        fig.tight_layout()
        fig.savefig(out/f'overview-T-{T:g}.png',dpi=160,bbox_inches='tight')
        for ax in axes.flat:ax.set_xlim(-c['lead_steps']*dt-3,T+30)
        fig.suptitle(f'Event and early-recovery zoom; N=20; n=3 development seeds; T={T:g}\nLines: mean; shading: seed range; orange: missing-observation window')
        finish(fig,out/f'overview-event-zoom-T-{T:g}.png')
        fig,axes=plt.subplots(2,3,figsize=(16,8))
        diagnostics=[('A','Social coefficient sum A'),('C','Empirical bearing factor C'),
                     ('S_abs_error','Same-state |S injected - S oracle|'),('S_memory','S from remembered hidden leaders'),
                     ('bearing_error_deg','Active hidden-edge bearing error (deg)'),('shadow_kf_coverage','Shadow KF nominal-95% coverage')]
        for ax,(metric,label) in zip(axes.flat,diagnostics):
            for mode in available:
                if metric=='shadow_kf_coverage' and mode not in ['kf-mean','kf-uncertainty']:continue
                if metric=='bearing_error_deg' and mode in ['fresh','drop','drop-normalized']:continue
                ds=[raw[key(s,mode,T)] for s in SEEDS];time=ds[0]['time']
                values=np.array([d[metric] for d in ds])
                ax.plot(time,values.mean(0),color=COLORS[mode],lw=1.4,label=LABELS[mode])
            if metric=='shadow_kf_coverage':ax.axhline(.95,color='black',ls=':',lw=.8)
            ax.set_ylabel(label);shade(ax,T);ax.set_xlim(-c['lead_steps']*dt-3,55)
        axes[0,0].legend(fontsize=8,ncol=2)
        fig.suptitle(f'Social-input and prediction diagnostics; n=3 development seeds; T={T:g}\nS/A/C exclude leader task cue; oracle quantities are evaluation only')
        fig.tight_layout();finish(fig,out/f'social-input-T-{T:g}.png')
        # Every seed shown, without selecting a visually successful example.
        for seed in SEEDS:
            fig,axes=plt.subplots(2,3,figsize=(16,8))
            for ax,(metric,label) in zip(axes.flat,metrics):
                for mode in available:
                    d=raw[key(seed,mode,T)]
                    ax.plot(d['time'],d[metric],color=COLORS[mode],lw=1.3,label=LABELS[mode])
                ax.set_ylabel(label);shade(ax,T);ax.set_xlim(-c['lead_steps']*dt-3,d['time'][-1])
            axes[0,0].legend(fontsize=8,ncol=2)
            fig.suptitle(f'Development seed {seed}; T={T:g}; raw unsmoothed curves')
            fig.tight_layout();finish(fig,out/f'seed-{seed}-curves-T-{T:g}.png')
    # Heading task gate and straight-control comparison at the longest horizon.
    T=max(Ts)
    fig,axes=plt.subplots(3,3,figsize=(16,11))
    for row,seed in enumerate(SEEDS):
        f=raw[key(seed,'fresh',T)];d=raw[key(seed,'drop',T)]
        axes[row,0].plot(f['time'],f['command_deg'],'k--',label='Command')
        axes[row,0].plot(f['time'],f['leader_turn_angle_deg'],color='#009e73',label='Fresh leaders')
        axes[row,0].plot(f['time'],f['follower_turn_angle_deg'],color=COLORS['fresh'],label='Fresh followers')
        axes[row,0].plot(d['time'],d['follower_turn_angle_deg'],color=COLORS['drop'],label='Drop followers')
        axes[row,0].set_ylabel(f'{seed}\nMean direction (deg)')
        for turn,ls in [(1,'-'),(0,'--')]:
            for mode in ['fresh','drop']:
                a=raw[key(seed,mode,T,turn)]
                axes[row,1].plot(a['time'],a['follower_turn_error_deg'],ls,color=COLORS[mode],label=f'{mode} {"turn" if turn else "straight"}')
        axes[row,1].set_ylabel('Follower command error (deg)')
        axes[row,2].plot(f['time'],f['hidden_bearing_change_deg'],color='#0072b2',label='Fresh leader-bearing change')
        axes[row,2].set_ylabel('Mean hidden-edge bearing change (deg)')
        for ax in axes[row]:shade(ax,T);ax.set_xlim(-c['lead_steps']*dt-3,f['time'][-1])
    for ax in axes[0]:ax.legend(fontsize=8)
    fig.suptitle('Task validity and information excitation; all 3 development seeds')
    fig.tight_layout();finish(fig,out/'task-validity.png')
    # Trajectories at every available T, using a common spatial extent per seed.
    for T in Ts:
        for seed in SEEDS:
            modes=[m for m in ['fresh','drop','hold','kf-mean','kf-uncertainty'] if m in available]
            fig,axes=plt.subplots(1,len(modes),figsize=(4.1*len(modes),4.8),squeeze=False)
            points=[]
            for ax,mode in zip(axes.flat,modes):
                z=np.load(stage/'runs'/(key(seed,mode,T)+'.npz'))
                xy=z['xy'];lo=start-c['lead_steps'];hi=min(start+round(T/dt)+100,len(xy)-1)
                points.append(xy[lo:hi+1].reshape(-1,2))
                for j in range(20):
                    color='#009e73' if j<3 else '#687d9a'
                    ax.plot(xy[lo:hi+1,j,0],xy[lo:hi+1,j,1],color=color,lw=1.4 if j<3 else .55,alpha=.9 if j<3 else .55)
                for ix,marker,label in [(lo,'o','Turn onset'),(start,'s','Blackout onset'),(start+round(T/dt),'x','Restored')]:
                    ax.scatter(xy[ix,:,0],xy[ix,:,1],s=17,marker=marker,c=['#009e73']*3+['#687d9a']*17,label=label)
                ax.set_title(LABELS[mode]);ax.set_xlabel('x (model units)');ax.set_aspect('equal',adjustable='box');ax.grid(alpha=.2)
            allxy=np.concatenate(points);lower=allxy.min(0)-3;upper=allxy.max(0)+3
            for ax in axes.flat:ax.set_xlim(lower[0],upper[0]);ax.set_ylim(lower[1],upper[1])
            axes[0,0].set_ylabel('y (model units)');axes[0,0].legend(fontsize=7,loc='best')
            fig.suptitle(f'Seed {seed}; T={T:g}; all20 trajectories; green=3 leaders; blue=17 followers\nFrom turn onset to 100 steps after restoration; same axes across methods')
            fig.tight_layout();finish(fig,out/f'seed-{seed}-trajectories-T-{T:g}.png')
    write_csv(out/'summary.csv',summaries);write_csv(out/'task-gates.csv',gates)
    # Paired seed bootstrap, explicitly exploratory at n=3, retain seed values.
    effects=[];rng=np.random.default_rng(2026101400)
    for T in Ts:
        for metric in ['event_turn_error_deg','event_GO_deviation','post_fixed5_distance_abs_deviation','event_S_abs_error','event_bearing_error_deg']:
            for method,base in [('drop','fresh'),('kf-mean','hold'),('kf-uncertainty','hold')]:
                if method not in available or base not in available:continue
                a=[next(r[metric] for r in summaries if r['seed']==s and r['T']==T and r['mode']==method) for s in SEEDS]
                b=[next(r[metric] for r in summaries if r['seed']==s and r['T']==T and r['mode']==base) for s in SEEDS]
                diff=np.array(a)-b
                if not np.isfinite(diff).all():continue
                boot=diff[rng.integers(0,3,(5000,3))].mean(1)
                effects.append(dict(T=T,metric=metric,method=method,baseline=base,n=3,
                                    mean_difference=diff.mean(),ci95_low=np.quantile(boot,.025),ci95_high=np.quantile(boot,.975),
                                    lower_count=int((diff<0).sum()),seed1=diff[0],seed2=diff[1],seed3=diff[2]))
    write_csv(out/'paired-effects.csv',effects)
    hashes={str(p.relative_to(stage)):hashlib.sha256(p.read_bytes()).hexdigest() for p in stage.rglob('*') if p.is_file() and 'analysis' not in p.parts}
    code_hashes={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in Path(__file__).parent.glob('*') if p.is_file()}
    audit=dict(n_seeds=3,gate_rows=len(gates),summary_rows=len(summaries),raw_runs=len(raw),
               prefix_and_T0_assertions='Passed during runner execution for all seeds and all methods',
               no_exclusions=True,matplotlib=matplotlib.__version__,hashes=hashes,code_hashes=code_hashes)
    (out/'audit.json').write_text(json.dumps(audit,indent=2),encoding='utf-8')
    print('GATES',json.dumps(gates,indent=2),flush=True)
    print('PAIRED',json.dumps(effects,indent=2),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--stage',type=Path,required=True)
    analyze(p.parse_args().stage)

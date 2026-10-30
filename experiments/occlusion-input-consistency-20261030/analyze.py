"""Analyze an explicitly executed run, preserving seed pairing across all S0."""
from pathlib import Path
import argparse,json,hashlib
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.io import loadmat

ROOT=Path(__file__).resolve().parent
CFG=json.loads((ROOT/'config.v1.json').read_text(encoding='utf-8'))

def validate_rows(df,profiles):
    required=['seed','profile','remaining','initialSocialInput','finalConsistency','finalGO','finalLO']
    assert set(required)<=set(df.columns)
    assert not df.duplicated(['seed','profile','remaining']).any()
    assert np.isfinite(df[required[2:]].to_numpy(dtype=float)).all()
    for p in profiles:
        z=df[df.profile==p['profile']]
        expected=set(range(p['seeds_first'],p['seeds_first']+p['seed_count']))
        assert set(z.seed)==expected
        assert len(z)==len(expected)*len(CFG['fresh_counts'])
        for seed,g in z.groupby('seed'):
            assert set(g.remaining)==set(CFG['fresh_counts'])
        assert np.allclose(z.initialSocialInput,CFG['model']['single_edge_weight']*z.remaining,atol=1e-12,rtol=0)
    assert set(df.profile)=={p['profile'] for p in profiles}
    for key in ['finalConsistency','finalGO','finalLO']:
        assert df[key].between(-1e-12,1+1e-12).all()
    assert (df.finalConsistency<=df.finalGO+1e-12).all()
    assert (df.finalConsistency<=df.finalLO+1e-12).all()

def paired_aggregate(df,profiles,repeats=10000):
    rng=np.random.default_rng(CFG['rng']['bootstrap_seed']); rows=[]
    metrics=['finalConsistency','finalGO','finalLO','finalConsistencyQ10',
             'finalJointHighFraction','finalMovingFraction','finalMeanSpeed','finalPairDistance','C_last500','C_last2000']
    for p in profiles:
        z=df[df.profile==p['profile']]; seeds=sorted(z.seed.unique())
        # One resampling index array applies to every r and every outcome.
        index=rng.integers(len(seeds),size=(repeats,len(seeds)))
        base={r:dict(profile=p['profile'],remaining=r,initialSocialInput=r*CFG['model']['single_edge_weight'],n=len(seeds)) for r in CFG['fresh_counts']}
        for metric in metrics:
            table=z.pivot(index='seed',columns='remaining',values=metric).loc[seeds,CFG['fresh_counts']].to_numpy(float)
            boot=np.empty((repeats,len(CFG['fresh_counts'])))
            for start in range(0,repeats,200):
                boot[start:start+200]=table[index[start:start+200]].mean(axis=1)
            lo,hi=np.quantile(boot,[.025,.975],axis=0)
            for col,r in enumerate(CFG['fresh_counts']):
                base[r].update({metric:table[:,col].mean(),metric+'_lo':lo[col],metric+'_hi':hi[col]})
        rows.extend(base.values())
    return pd.DataFrame(rows)

def verify_raw(run,df,profiles):
    checks=[]
    for p in profiles:
        seed=p['seeds_first']
        for r in [0,8,19]:
            file=run/f'seed-{seed}'/f'{p["profile"]}-r{r:02d}.mat'
            data=loadmat(file,simplify_cells=True); ts=data['series']; vel=data['velocity']; pos=data['position']
            # Reconstruct GO and LO independently from stored vectors and positions.
            go=[];lo=[]
            for t in range(CFG['post_steps']-CFG['terminal_steps'],CFG['post_steps']):
                v=vel[:,:,t]; speed=np.linalg.norm(v,axis=0); u=np.zeros_like(v); active=speed>1e-12
                u[:,active]=v[:,active]/speed[active]
                go.append(np.linalg.norm(u.sum(axis=1))/CFG['model']['N'])
                pp=pos[:,:,t+1]; delta=pp[:,:,None]-pp[:,None,:]; L=CFG['model']['L']; delta-=L*np.round(delta/L)
                distance=np.linalg.norm(delta,axis=0); local=[]
                for i in range(CFG['model']['N']):
                    others=[j for j in range(CFG['model']['N']) if j!=i]
                    nearest=sorted(others,key=lambda j:(distance[i,j],j))[:5]
                    local.append(np.linalg.norm(u[:,[i,*nearest]].sum(axis=1))/6)
                lo.append(np.mean(local))
            reconstructed=float(np.mean(np.minimum(go,lo)))
            row=df[(df.seed==seed)&(df.profile==p['profile'])&(df.remaining==r)].iloc[0]
            err=abs(reconstructed-row.finalConsistency)
            assert err<1e-10
            assert np.max(abs(np.asarray(go)-ts[-CFG['terminal_steps']:,0]))<1e-10
            assert np.max(abs(np.asarray(lo)-ts[-CFG['terminal_steps']:,1]))<1e-10
            checks.append(dict(seed=seed,profile=p['profile'],remaining=r,absolute_error=err))
    return checks

def main(run_id):
    assert run_id and '/' not in run_id and '\\' not in run_id
    run=ROOT/'results'/run_id
    assert json.loads((run/'config.json').read_text(encoding='utf-8'))==CFG,'Use the frozen run configuration.'
    execution=json.loads((run/'execution.json').read_text(encoding='utf-8')); assert execution['status']=='complete'
    assert not execution['seedOverride'],'Primary analysis requires the full frozen seed set.'
    profiles=[CFG['primary'],*CFG['history_controls']]
    files=sorted(run.glob('seed-*/summary.csv')); assert len(files)==CFG['primary']['seed_count']
    df=pd.concat([pd.read_csv(p) for p in files],ignore_index=True); validate_rows(df,profiles)
    out=run/'analysis-v1'; assert not out.exists(),'Use a fresh analysis directory/version.'; out.mkdir()
    df.to_csv(out/'by_seed.csv',index=False); a=paired_aggregate(df,profiles,CFG['statistics']['bootstrap_repeats'])
    a.to_csv(out/'aggregate.csv',index=False); checks=verify_raw(run,df,profiles)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(1,3,figsize=(15,4.8),sharex=True,sharey=True)
    for ax,p in zip(axes,profiles):
        z=a[a.profile==p['profile']].sort_values('initialSocialInput')
        for metric,label,color in [('finalConsistency','Joint C = mean min(GO, LO)','#146d90'),('finalGO','GO','#d88932'),('finalLO','LO','#6e9843')]:
            ax.plot(z.initialSocialInput,z[metric],'o-',label=label,color=color,ms=3)
            ax.fill_between(z.initialSocialInput,z[metric+'_lo'],z[metric+'_hi'],alpha=.14,color=color)
        ax.set_title(f"{p['profile']}: pre={p['pre_steps']}, n={p['seed_count']}")
        ax.set_xlabel('Initial social input S0 after occlusion'); ax.set_ylim(0,1.02); ax.grid(alpha=.2)
    axes[0].set_ylabel('Final movement consistency');axes[1].legend(fontsize=8)
    fig.suptitle('Uncompensated model: permanent occlusion; fixed terminal window; paired-seed 95% intervals')
    fig.tight_layout();fig.savefig(out/'input_vs_final_consistency.png',dpi=200);fig.savefig(out/'input_vs_final_consistency.pdf');plt.close(fig)
    fig,ax=plt.subplots(figsize=(8,5));p=CFG['primary'];z=df[df.profile==p['profile']]
    ax.scatter(z.initialSocialInput,z.finalConsistency,s=4,alpha=.035,color='#146d90')
    agg=a[a.profile==p['profile']].sort_values('initialSocialInput');ax.plot(agg.initialSocialInput,agg.finalConsistency,'o-',color='#146d90')
    ax.fill_between(agg.initialSocialInput,agg.finalConsistency_lo,agg.finalConsistency_hi,alpha=.25,color='#146d90')
    ax.set(xlabel='Initial social input S0 after occlusion',ylabel='Final joint movement consistency C',ylim=(0,1.02));ax.grid(alpha=.2)
    fig.tight_layout();fig.savefig(out/'primary_seed_distribution.png',dpi=200);fig.savefig(out/'primary_seed_distribution.pdf');plt.close(fig)
    audit=dict(status='PASS',rows=len(df),independent_seeds=500,raw_recompute=checks,
        source_summaries={str(p.relative_to(run)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files},
        confidence_unit='paired seed; 20 occlusion levels retained in each resample',
        finite_time_endpoint='terminal 1000 of 6000 post-event steps',model_changes='none: original neural dynamics, no prediction')
    (out/'audit.json').write_text(json.dumps(audit,indent=2),encoding='utf-8')

if __name__=='__main__':
    parser=argparse.ArgumentParser(); parser.add_argument('--run',required=True);args=parser.parse_args(); main(args.run)

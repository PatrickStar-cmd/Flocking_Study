"""Revised input axis over parent-protocol raw data; runs no model simulation."""
from pathlib import Path
import argparse,json,sys,hashlib
import numpy as np
import pandas as pd
from scipy.io import loadmat
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from three_factor_input import checkpoint_input

HERE=Path(__file__).resolve().parent;PARENT=HERE.parent
sys.path.insert(0,str(PARENT))
import analyze as base
CFG=json.loads((HERE/'config.v2.json').read_text(encoding='utf-8'))

def revised_rows(run,df,profiles):
    out=df.rename(columns={'initialSocialInput':'initialTotalStrength','initialVectorX':'initialLegacyVectorX','initialVectorY':'initialLegacyVectorY'}).copy()
    out['initialSocialInput']=np.nan;out['initialAngleWeight']=np.nan
    out['initialVectorX']=np.nan;out['initialVectorY']=np.nan
    for (seed,profile),g in out.groupby(['seed','profile']):
        cp=loadmat(run/f'seed-{seed}'/f'{profile}-checkpoint.mat',simplify_cells=True)['checkpoint']
        for r in base.CFG['fresh_counts']:
            z=checkpoint_input(cp,r,base.CFG['model'],CFG['kappa']);ix=g.index[g.remaining==r]
            assert len(ix)==1
            out.loc[ix,'initialSocialInput']=z['S'].mean();out.loc[ix,'initialAngleWeight']=z['angle_weight'].mean()
            out.loc[ix,'initialVectorX']=z['Vx'].mean();out.loc[ix,'initialVectorY']=z['Vy'].mean()
            assert abs(z['A'].mean()-out.loc[ix,'initialTotalStrength'].iloc[0])<1e-12
    assert np.isfinite(out.initialSocialInput).all()
    assert (out.initialSocialInput<=out.initialTotalStrength+1e-12).all()
    for _,g in out.groupby(['seed','profile']):
        assert (np.diff(g.sort_values('remaining').initialSocialInput)>=-1e-12).all()
    return out

def input_bins(df,repeats=10000):
    """Ratio of within-bin sums, with seed-cluster resampling for the curve."""
    rows=[];rng=np.random.default_rng(CFG['bootstrap_seed'])
    for profile,z in df.groupby('profile',sort=True):
        z=z.copy();x=z.initialSocialInput.to_numpy(float)
        width=CFG['bins']['positive_width'];maximum=CFG['bins']['maximum'];nbins=round(maximum/width)
        assert (x>=0).all() and (x<=maximum+1e-12).all()
        z['bin']=np.where(x==0,-1,np.minimum(np.floor(x/width).astype(int),nbins-1))
        seeds=sorted(z.seed.unique());bins=sorted(z['bin'].unique())
        grouped=z.groupby(['seed','bin'])
        def matrix(series):
            return series.unstack('bin',fill_value=0).reindex(index=seeds,columns=bins,fill_value=0).to_numpy(float)
        count=matrix(grouped.size());sx=matrix(grouped.initialSocialInput.sum());sy=matrix(grouped.finalConsistency.sum())
        boot=np.full((repeats,len(bins)),np.nan)
        for start in range(0,repeats,100):
            index=rng.integers(len(seeds),size=(min(100,repeats-start),len(seeds)))
            den=count[index].sum(axis=1);num=sy[index].sum(axis=1)
            boot[start:start+100]=np.divide(num,den,out=np.full_like(num,np.nan),where=den>0)
        for col,b in enumerate(bins):
            valid=boot[:,col][np.isfinite(boot[:,col])];assert len(valid)>0
            lo,hi=np.quantile(valid,[.025,.975]);n=count[:,col].sum()
            rows.append(dict(profile=profile,input_bin=int(b),input_mean=sx[:,col].sum()/n,
                consistency_mean=sy[:,col].sum()/n,lo=lo,hi=hi,cases=int(n),seeds=int(np.count_nonzero(count[:,col])),
                bootstrap_valid_fraction=len(valid)/repeats))
    return pd.DataFrame(rows)

def main(run_id):
    assert run_id and '/' not in run_id and '\\' not in run_id
    run=PARENT/'results'/run_id
    execution=json.loads((run/'execution.json').read_text(encoding='utf-8'))
    assert execution.get('status')=='complete' and not execution['seedOverride']
    assert json.loads((run/'config.json').read_text(encoding='utf-8'))==base.CFG
    profiles=[base.CFG['primary'],*base.CFG['history_controls']]
    files=sorted(run.glob('seed-*/summary.csv'));assert len(files)==500
    df=pd.concat([pd.read_csv(p) for p in files],ignore_index=True);base.validate_rows(df,profiles)
    out=run/'analysis-v2';assert not out.exists();out.mkdir()
    revised=revised_rows(run,df,profiles);revised.to_csv(out/'by_seed.csv',index=False)
    binned=input_bins(revised,CFG['bootstrap_repeats']);binned.to_csv(out/'input_bins.csv',index=False)
    # Keep the raw degree-conditioned summaries, but replace their x coordinate by actual S0.
    agg=base.paired_aggregate(df,profiles,CFG['bootstrap_repeats']).rename(columns={'initialSocialInput':'initialTotalStrength'})
    means=revised.groupby(['profile','remaining']).initialSocialInput.mean().rename('initialSocialInput').reset_index()
    agg=agg.merge(means,on=['profile','remaining'],validate='one_to_one');agg.to_csv(out/'condition_means.csv',index=False)
    raw_checks=base.verify_raw(run,df,profiles)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(1,3,figsize=(15,5),sharex=True,sharey=True)
    for ax,p in zip(axes,profiles):
        z=revised[revised.profile==p['profile']];b=binned[binned.profile==p['profile']].sort_values('input_mean')
        ax.scatter(z.initialSocialInput,z.finalConsistency,s=4,alpha=.035,color='#146d90')
        ax.plot(b.input_mean,b.consistency_mean,'o-',color='#146d90',ms=4,label='Prespecified S0 bins')
        ax.fill_between(b.input_mean,b.lo,b.hi,alpha=.2,color='#146d90')
        ax.set_title(f"{p['profile']}; n={p['seed_count']}");ax.set_xlabel('Initial three-factor social input S0');ax.set_ylim(0,1.02);ax.grid(alpha=.2)
    axes[0].set_ylabel('Final joint movement consistency');axes[1].legend(fontsize=8)
    fig.suptitle('S0 = mean sum w exp[cos(theta - event heading) - 1]; seed-cluster 95% intervals')
    fig.tight_layout();fig.savefig(out/'three_factor_input_vs_consistency.png',dpi=200);fig.savefig(out/'three_factor_input_vs_consistency.pdf');plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(11,4.7),sharey=True);primary=base.CFG['primary']['profile']
    z=revised[revised.profile==primary];a=agg[agg.profile==primary]
    for ax,key,label in zip(axes,['initialSocialInput','initialTotalStrength'],['Three-factor S0','Total strength A0']):
        ax.scatter(z[key],z.finalConsistency,s=4,alpha=.035,color='#146d90')
        ax.plot(a[key],a.finalConsistency,'o-',color='#146d90',ms=4)
        ax.set_xlabel(label);ax.set_ylim(0,1.02);ax.grid(alpha=.2)
    axes[0].set_ylabel('Final joint movement consistency');fig.suptitle('Same trajectories, two input summaries; lines connect remaining-count means')
    fig.tight_layout();fig.savefig(out/'social_input_and_strength_control.png',dpi=200);fig.savefig(out/'social_input_and_strength_control.pdf');plt.close(fig)
    audit=dict(status='PASS',measurement='v2 early-v3 three-factor input, kappa=1',rows=len(revised),
        raw_consistency_checks=raw_checks,source_model_run=run_id,models_executed_by_this_analyzer=0,
        source_summaries={str(p.relative_to(run)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files},
        binning='fixed 0.012-wide positive bins; exact-zero separate; case means with seed-cluster bootstrap')
    (out/'audit.json').write_text(json.dumps(audit,indent=2)+'\n',encoding='utf-8',newline='\n')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--run',required=True);args=p.parse_args();main(args.run)

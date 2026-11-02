"""Only development data enter the immutable global-Q and sample-size decision."""
import hashlib,json,math,subprocess,sys
from pathlib import Path
import numpy as np
import pandas as pd


def sha(path):return hashlib.file_digest(path.open('rb'),'sha256').hexdigest()


def choose(stage):
    assert 'development' in stage.name
    target=stage/'FROZEN_DECISION.json';assert not target.exists()
    reviewed=json.loads((stage/'review/independent_audit.json').read_text(encoding='utf-8'))
    assert reviewed['status']=='PASS' and reviewed['independent_seeds']==12
    frame=pd.read_csv(stage/'analysis/seed_metrics.csv');frame=frame[frame.scenario=='combined']
    coverage=pd.read_csv(stage/'analysis/coverage_by_age.csv.gz')
    cuts=json.loads((stage/'GEOMETRY_FROZEN.json').read_text(encoding='utf-8'))['rate_cutpoints']
    sweep=[];pairs={}
    for q in [1e-6,1e-5,1e-4,1e-3,1e-2]:
        method=f'CV_q{q:g}'
        summaries={}
        for group,D in [('clean-fast',120),('clean-slow',20)]:
            f=frame[(frame.geometry==group)&(frame.D==D)]
            hold=f[f.method=='Hold'].set_index('seed')
            model=f[f.method==method].set_index('seed')
            ids=hold.index.intersection(model.index)
            h=hold.loc[ids];m=model.loc[ids]
            support=(h.edge_windows>=5)&(m.edge_windows>=5)
            # Support failures are not deleted from published seed tables.
            enough=len(ids)>=10 and bool(support.all())
            summaries[group]={'n_seeds':len(ids),'enough':enough,'hold_mean':h.rmse.mean(),
                'model_mean':m.rmse.mean(),'delta':m.rmse-h.rmse,'ratio':m.rmse.mean()/h.rmse.mean()}
        cal=coverage[(coverage.method==method)&(((coverage.geometry=='clean-fast')&(coverage.D==120))|((coverage.geometry=='clean-slow')&(coverage.D==20)))]
        band_seed=cal.groupby(['geometry','D','age_band','seed']).coverage.mean()
        band=band_seed.groupby(['geometry','D','age_band']).mean()
        low=float(band.min());high=float(band.max())
        severe=not np.isfinite(low+high) or low<.80 or high>.995
        fast=summaries['clean-fast'];slow=summaries['clean-slow']
        admissible=fast['enough'] and slow['enough'] and slow['ratio']<=1.10 and not severe
        sweep.append({'q':q,'method':method,'objective_fast_D120_RMSE':fast['model_mean'],
            'fast_D120_hold_RMSE':fast['hold_mean'],'fast_D120_delta':fast['delta'].mean(),
            'slow_D20_ratio_to_Hold':slow['ratio'],'fast_n_seeds':fast['n_seeds'],'slow_n_seeds':slow['n_seeds'],
            'support_enough':fast['enough'] and slow['enough'],'age_coverage_min':low,'age_coverage_max':high,
            'coverage_all_reference_90_98':bool(((band>=.90)&(band<=.98)).all()),
            'severe_miscalibration':severe,'admissible':admissible})
        pairs[q]=fast['delta']
    table=pd.DataFrame(sweep);table.to_csv(stage/'analysis/q_sweep.csv',index=False)
    allowed=table[table.admissible]
    if len(allowed):
        best=allowed.objective_fast_D120_RMSE.min()
        equivalent=allowed[allowed.objective_fast_D120_RMSE<=1.01*best].copy()
        equivalent['distance']=np.abs(np.log10(equivalent.q/1e-4))
        selected=equivalent.sort_values(['distance','q']).iloc[0]
        q=float(selected.q);sd=float(pairs[q].std(ddof=1))
        n=max(30,int(math.ceil((1.96*sd/.025)**2)))
        validation_seeds=list(range(2026112001,2026112001+n))
        reason='minimum clean-fast D120 seed RMSE subject to frozen support/slow-protection/calibration rules'
    else:q=None;sd=None;n=0;validation_seeds=[];reason='No globally admissible q under the predeclared protection, support and severe-miscalibration constraints; validation gate remains closed.'
    here=Path(__file__).resolve().parent
    source_hashes={p.name:sha(p) for p in sorted(here.iterdir()) if p.suffix in {'.py','.m'} or p.name in {'PROTOCOL.md','RUN_RULES.md'}}
    repo=here.parent.parent
    upstream=list((repo/'effective-social-input').glob('*.m'))+list((repo/'experiments/kf-diagnostic-20261101').glob('*.m'))+[repo/'experiments/kf-p2-20261102/replay.py']
    upstream_hashes={p.relative_to(repo).as_posix():sha(p) for p in sorted(upstream)}
    decision={'status':'FROZEN' if q is not None else 'STOP_NO_ADMISSIBLE_Q','q_star':q,
        'reason':reason,'q_grid':[1e-6,1e-5,1e-4,1e-3,1e-2],'rate_cutpoints':cuts,
        'geometry':{'near_range':5,'half_box_margin':10,'aliasing_delta':.9*np.pi,'round':'MATLAB ties away from zero'},
        'source_commit':subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),
        'source_hashes':source_hashes,'upstream_source_hashes':upstream_hashes,'development_q_sweep_sha256':sha(stage/'analysis/q_sweep.csv'),
        'validation_n':n,'validation_seeds':validation_seeds,'development_paired_sd':sd,
        'desired_CI_halfwidth':.025,'minimum_scientific_difference_rad':.02,
        'primary':{'geometry':'clean-fast','D':120,'scenario':'combined','contrast':'CV(q*) - Hold, mean of within-seed window RMSE differences'},
        'protection':{'geometry':'clean-slow','D':20,'max_ratio':1.10},
        'bootstrap':{'seed':2026112100,'replicates':10000,'unit':'whole physical seed including all candidate windows and sources'},
        'validation_support':{'minimum_fraction_seeds':.8,'minimum_edge_windows_per_seed':5},
        'validation_not_generated_or_opened_before_this_file':True,
        'no_threshold_or_q_reselection_after_validation':True}
    target.write_text(json.dumps(decision,indent=2),encoding='utf-8')
    (stage/'FROZEN_DECISION.sha256').write_text(sha(target)+'\n',encoding='ascii')
    print(json.dumps({k:v for k,v in decision.items() if k not in ['source_hashes','validation_seeds']},indent=2))


if __name__=='__main__':choose(Path(sys.argv[1]))

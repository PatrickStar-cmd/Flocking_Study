"""Read-only frozen-decision checks and predeclared seed-cluster comparisons."""
import hashlib,json,sys
from pathlib import Path
import numpy as np
import pandas as pd


def sha(path):return hashlib.file_digest(path.open('rb'),'sha256').hexdigest()


def check(path):
    expected=(path.parent/'FROZEN_DECISION.sha256').read_text(encoding='ascii').strip()
    assert sha(path)==expected,'Frozen decision changed'
    frozen=json.loads(path.read_text(encoding='utf-8'))
    here=Path(__file__).resolve().parent
    for name,digest in frozen['source_hashes'].items():assert sha(here/name)==digest,('Frozen source changed',name)
    for name,digest in frozen['upstream_source_hashes'].items():assert sha(here.parent.parent/name)==digest,('Frozen upstream source changed',name)
    return frozen


def describe(values):
    values=np.asarray(values,dtype=float)
    finite=values[np.isfinite(values)]
    if not len(finite):return {'n_seeds':0,'mean':None,'median':None,'sd':None,'ci_lower':None,'ci_upper':None}
    rng=np.random.default_rng(2026112100)
    draws=values[rng.integers(0,len(values),size=(10000,len(values)))]
    count=np.isfinite(draws).sum(axis=1)
    bootstrap=np.divide(np.nansum(draws,axis=1),count,out=np.full(10000,np.nan),where=count>0)
    valid_bootstrap=bootstrap[np.isfinite(bootstrap)]
    return {'n_seeds':len(finite),'n_predeclared_seeds':len(values),'all_missing_bootstrap_draws':int((count==0).sum()),
        'mean':float(finite.mean()),'median':float(np.median(finite)),
        'sd':float(finite.std(ddof=1)) if len(finite)>1 else None,
        'ci_lower':float(np.quantile(valid_bootstrap,.025)),'ci_upper':float(np.quantile(valid_bootstrap,.975)),
        'fraction_KF_worse':float((finite>0).mean())}


def evaluate(stage,frozen_path):
    d=check(frozen_path)
    assert d['q_star'] is not None,'No globally admissible q: cannot claim independent calibrated validation'
    manifest=json.loads((stage/'input_manifest.json').read_text(encoding='utf-8'))
    assert sorted(map(int,manifest['seeds']))==d['validation_seeds']
    data=pd.read_csv(stage/'analysis/seed_metrics.csv')
    all_contrasts=[];grouped=[]
    for (geometry,D,scenario),f in data.groupby(['geometry','D','scenario']):
        hold=f[f.method=='Hold'].set_index('seed')
        for method in sorted(set(f.method)-{'Hold'}):
            model=f[f.method==method].set_index('seed');ids=hold.index.intersection(model.index)
            if not len(ids):continue
            delta=model.loc[ids].rmse-hold.loc[ids].rmse
            for seed,value in delta.items():all_contrasts.append({'seed':seed,'geometry':geometry,'D':D,'scenario':scenario,'method':method,'delta_RMSE':value,
                'Hold_RMSE':hold.loc[seed].rmse,'model_RMSE':model.loc[seed].rmse,'edge_windows':model.loc[seed].edge_windows,
                'relative_change':value/hold.loc[seed].rmse})
            grouped.append({'geometry':geometry,'D':D,'scenario':scenario,'method':method,
                'mean_relative_improvement':float((hold.loc[ids].rmse-model.loc[ids].rmse).mean()/hold.loc[ids].rmse.mean()),
                **describe(delta.reindex(d['validation_seeds']).to_numpy())})
    contrast=pd.DataFrame(all_contrasts);group=pd.DataFrame(grouped)
    out=stage/'analysis/validation_contrasts.csv';assert not out.exists();contrast.to_csv(out,index=False)
    group.to_csv(stage/'analysis/validation_summary.csv',index=False)
    method=f"CV_q{d['q_star']:g}"
    primary=group[(group.geometry=='clean-fast')&(group.D==120)&(group.scenario=='combined')&(group.method==method)]
    assert len(primary)==1
    primary=primary.iloc[0].to_dict()
    protection=contrast[(contrast.geometry=='clean-slow')&(contrast.D==20)&(contrast.scenario=='combined')&(contrast.method==method)]
    ratio=float(protection.model_RMSE.mean()/protection.Hold_RMSE.mean())
    support=contrast[(contrast.geometry=='clean-fast')&(contrast.D==120)&(contrast.scenario=='combined')&(contrast.method==method)]
    support_ok=len(support)>=np.ceil(.8*d['validation_n']) and bool((support.edge_windows>=5).all())
    report={'status':'PASS','frozen_decision_sha256':sha(frozen_path),'validation_n_predeclared':d['validation_n'],
        'primary':primary,'support_sufficient':support_ok,'slow_D20_ratio':ratio,
        'mean_reaches_002rad':primary['mean']<=-.02,'primary_CI_upper_below_zero':primary['ci_upper']<0,
        'short_term_protection_mean_ratio_pass':ratio<=1.10,
        'P3_not_executed':True,'interpretation':'Effect support, tail/coverage and generalization must be considered together; no GO or neural-field claim.'}
    path=stage/'analysis/validation_audit.json';assert not path.exists();path.write_text(json.dumps(report,indent=2),encoding='utf-8')
    print(json.dumps(report,indent=2))


if __name__=='__main__':
    if len(sys.argv)==2:print(json.dumps(check(Path(sys.argv[1])),indent=2))
    else:evaluate(Path(sys.argv[1]),Path(sys.argv[2]))

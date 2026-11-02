"""Independent raw verification of every calibration band used to stop/select Q."""
import json,sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.io import loadmat


def main(stage):
    target=stage/'review/calibration_audit.json';assert not target.exists()
    d=json.loads((stage/'FROZEN_DECISION.json').read_text(encoding='utf-8'))
    inputs=json.loads((stage/'input_manifest.json').read_text(encoding='utf-8'))
    table=pd.read_csv(stage/'analysis/coverage_by_age.csv.gz').set_index(['seed','scenario','method','D','t0','geometry','age_band'])
    methods=[f'CV_q{q:g}' for q in d['q_grid']]+['GM_default']
    verified=0;max_difference=0.;records=[]
    for item in inputs['inputs']:
        path=stage/item['source'];x=loadmat(path,simplify_cells=True);stem=path.stem
        for D,group in [(20,'clean-slow'),(120,'clean-fast')]:
            with np.load(stage/'labels'/f'{stem}_D{D}.npz') as label:
                labels=label['labels'];rates=label['rates'];starts=label['starts']
            ix=starts[None,:]-1+np.arange(D)[:,None];truth=x['true_bearing'][ix]
            for method in methods:
                with np.load(stage/'raw/replay'/f'{stem}_{method}_fresh.npz') as f:
                    with np.load(stage/'raw/replay'/f'{stem}_{method}_D{D}.npz') as b:
                        theta=np.bitwise_xor(b['theta_xor'][:D],np.ascontiguousarray(f['theta'][ix,0]).view(np.uint64)).view(np.float64)
                        P=np.bitwise_xor(b['P_xor'][:D],np.ascontiguousarray(f['P'][ix,0]).view(np.uint64)).view(np.float64)[...,0,0]
                err=np.arctan2(np.sin(theta-truth),np.cos(theta-truth))
                for c,t0 in enumerate(starts):
                    mask=(labels[c]=='G-clean')&((rates[c]<=d['rate_cutpoints'][0]) if D==20 else (rates[c]>d['rate_cutpoints'][1]))
                    if not mask.any():continue
                    for lo,hi in [(1,20),(21,60),(61,100),(101,120)]:
                        if lo>D:continue
                        e=err[lo-1:min(hi,D),c][:,mask];p=P[lo-1:min(hi,D),c]
                        coverage=float((np.abs(e)<=1.96*np.sqrt(p[:,None])).mean())
                        key=(int(item['seed']),item['scenario'],method,D,int(t0),group,f'{lo}-{hi}')
                        r=table.loc[key]
                        delta=abs(coverage-r.coverage);max_difference=max(max_difference,delta)
                        assert delta<2e-12
                        assert abs(np.sqrt(np.mean(e*e))-r.rmse)<2e-12
                        assert abs(p.mean()-r.mean_P)<2e-12
                        records.append({'seed':int(item['seed']),'method':method,'geometry':group,'age_band':f'{lo}-{hi}','coverage':coverage})
                        verified+=1
        print('CALIBRATION AUDIT PASS',stem,flush=True)
    raw=pd.DataFrame(records)
    folded=raw.groupby(['method','geometry','age_band','seed']).coverage.mean().groupby(['method','geometry','age_band']).mean()
    sweep=pd.read_csv(stage/'analysis/q_sweep.csv')
    for _,row in sweep.iterrows():
        cov=folded.xs(row.method,level='method')
        assert abs(cov.min()-row.age_coverage_min)<2e-12
        assert abs(cov.max()-row.age_coverage_max)<2e-12
    report={'status':'PASS','raw_calibration_band_rows_verified':verified,'max_coverage_difference':max_difference,
        'all_Q_screening_coverage_min_max_independently_verified':True,
        'no_criteria_changed_after_development_outcomes':True}
    target.write_text(json.dumps(report,indent=2),encoding='utf-8');print(json.dumps(report))


if __name__=='__main__':main(Path(sys.argv[1]))

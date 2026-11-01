"""Independent per-edge circular errors and estimator CSV verification."""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.io import loadmat


def main(run_dir,output):
    run=Path(run_dir); df=pd.read_csv(run/'analysis/estimator_diagnostics.csv')
    keys=['seed','D','mode','step_in_occ']; assert not df.duplicated(keys).any()
    df=df.set_index(keys); max_error=0.; checked=0; seed_stats=[]
    for path in sorted((run/'raw').glob('*.mat')):
        raw=loadmat(path,simplify_cells=True)['result']; cfg=raw['config']
        duration=int(cfg['occlusionDurationSteps']); seed=int(cfg['seed']); mode=cfg['predictionMode']
        if duration==0: continue
        all_errors=[]
        for age in range(1,duration+1):
            rec=raw['observationRecords'][1399+age]
            selected=np.isfinite(rec['trueBearing']); fresh=rec['fresh'].astype(bool)
            hidden=selected&~fresh
            truth=rec['trueBearing'][hidden]; injected=rec['injectedBearing'][hidden]
            candidate=(rec['rawX'][:,:,0][hidden] if mode in ('kalman','kf-mean') else
                       injected if mode=='hold' else np.full(len(truth),np.nan))
            error=np.arctan2(np.sin(candidate-truth),np.cos(candidate-truth))
            valid=np.isfinite(error); error=error[valid]
            if len(error): all_errors.extend(error.tolist())
            inj_error=np.arctan2(np.sin(injected-truth),np.cos(injected-truth))
            active=rec['active'].astype(bool)[hidden]&np.isfinite(inj_error)
            inj_error=inj_error[active]
            p=rec['rawVariance'][hidden]; codes=rec['cutoffReason'][hidden]
            values={'n_hidden':len(truth),
                    'mae_candidate':np.mean(np.abs(error)) if len(error) else np.nan,
                    'rmse_candidate':np.sqrt(np.mean(error**2)) if len(error) else np.nan,
                    'mae_injected':np.mean(np.abs(inj_error)) if len(inj_error) else np.nan,
                    'rmse_injected':np.sqrt(np.mean(inj_error**2)) if len(inj_error) else np.nan,
                    'mean_P11':np.mean(p) if mode in ('kalman','kf-mean') else np.nan,
                    'coverage95_candidate':np.mean(np.abs(error)<=1.96*np.sqrt(p[valid])) if mode in ('kalman','kf-mean') else np.nan,
                    'fresh_nis_mean':np.nanmean(rec['rawNis'][fresh]) if mode in ('kalman','kf-mean') else np.nan,
                    'P_trunc_fraction':np.isin(codes,[22,23]).mean(),
                    'age_trunc_fraction':np.isin(codes,[21,23]).mean(),
                    'n_active_hidden':active.sum()}
            row=df.loc[(seed,duration,mode,age)]
            for name,value in values.items():
                actual=float(row[name])
                if np.isnan(value):assert np.isnan(actual),(seed,duration,mode,age,name)
                else:
                    difference=abs(value-actual); max_error=max(max_error,difference)
                    assert np.isclose(value,actual,rtol=1e-10,atol=1e-12),(seed,duration,mode,age,name,value,actual)
            checked+=1
        if all_errors:
            seed_stats.append({'seed':seed,'D':duration,'mode':mode,
                'whole_occ_RMSE_rad':float(np.sqrt(np.mean(np.square(all_errors)))),
                'whole_occ_MAE_rad':float(np.mean(np.abs(all_errors)))})
        print(f'DIAGNOSTIC_AUDIT seed={seed} mode={mode} D={duration}',flush=True)
    assert checked==2800
    report={'status':'PASS','rows_checked':checked,'max_absolute_error':max_error,
            'note':'atan2(sin(error),cos(error)) independently reconstructs own-trajectory errors',
            'per_seed_whole_window_errors':seed_stats}
    with Path(output).open('x',encoding='utf-8') as stream:json.dump(report,stream,indent=2)
    print('PASS independent estimator diagnostics',flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('run_dir');parser.add_argument('output')
    args=parser.parse_args();main(args.run_dir,args.output)

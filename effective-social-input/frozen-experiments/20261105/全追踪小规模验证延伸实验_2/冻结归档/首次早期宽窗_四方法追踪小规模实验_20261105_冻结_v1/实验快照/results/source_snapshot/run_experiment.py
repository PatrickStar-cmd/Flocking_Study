"""Run all 65 fixed single-window trajectories, with hash-checked resume."""
import bootstrap
from pathlib import Path
import importlib.metadata
import os
import platform
import sys
import time
import numpy as np
import design as d
from common import ROOT,RESULTS,read_json,write_json,write_csv,digest,freeze_engine
from model import simulate
from metrics import series


def preflight():
    assert d.STEPS==120 and d.N==12 and len(d.jobs())==65
    assert d.PROFILES==('moving_moderate',) and d.SEEDS==(971001,971002,971003,971004,971005)
    assert d.STARTS==(20,40) and d.DURATIONS==(30,40)
    selected,hidden=d.topology()
    assert selected.sum()==77 and hidden.sum()==44 and not np.diag(selected).any()
    for start in d.STARTS:
        for width in d.DURATIONS:
            assert d.windows(width,start).sum()==width
            assert start>=20 and start+width+20<=d.STEPS
    provenance=read_json(ROOT/'SOURCE_PROVENANCE.json')
    assert digest(ROOT/'kf_reference_core.py')==provenance['unchanged_reference_core_sha256']
    source=Path(provenance['source_folder'])
    regression=[]
    if source.exists():
        for name,sha in provenance['source_files_sha256'].items():assert digest(source/name)==sha,name
        # Old seed only checks implementation regression, never enters new-seed outcomes.
        for method in d.METHODS:
            run=simulate(970001,'moving_moderate',method,0 if method=='fresh' else 40,0 if method=='fresh' else 40)
            old_id=f'970001_moving_moderate_{method}_D{0 if method=="fresh" else 40:02d}'
            with np.load(source/'results/trajectories'/f'{old_id}.npz') as old:
                for key in ('positions','headings'):
                    assert np.array_equal(run[key][:111],old[key][:111]),(method,key)
                for key in ('velocities','centers','active','filter_state','filter_covariance'):
                    assert np.array_equal(run[key][:110],old[key][:110]),(method,key)
                assert np.array_equal(run['noise'],old['noise'][:120])
            regression.append(method)
    write_json(RESULTS/'checks/PREFLIGHT.json',dict(status='PASS',planned_trajectories=65,source_folder_present=source.exists(),
        previous_experiment_unchanged=source.exists(),old_seed_prefix_regression_methods=regression,
        old_seed_is_not_an_experimental_sample=True,single_window_unclipped=True,reference_core_unchanged=True))


def main():
    RESULTS.mkdir(exist_ok=True)
    freeze_engine()
    preflight()
    for name in ('trajectories','timeseries'):(RESULTS/name).mkdir(exist_ok=True)
    write_json(RESULTS/'environment.json',dict(python=sys.version,executable=sys.executable,platform=platform.platform(),
        packages={p:importlib.metadata.version(p) for p in ('numpy','matplotlib')},
        threads={p:os.environ.get(p) for p in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS')}))
    manifest=RESULTS/'run_manifest.json'
    records=read_json(manifest) if manifest.exists() else []
    lookup={r['run_id']:r for r in records}
    jobs=d.jobs();assert len(set(d.run_id(*job) for job in jobs))==65
    begin=time.perf_counter()
    for seed,method,start,width in jobs:
        rid=d.run_id(seed,method,start,width)
        npzpath=RESULTS/'trajectories'/f'{rid}.npz';csvpath=RESULTS/'timeseries'/f'{rid}.csv'
        if rid in lookup:
            assert digest(npzpath)==lookup[rid]['npz_sha256'] and digest(csvpath)==lookup[rid]['csv_sha256'],rid
            continue
        run=simulate(seed,'moving_moderate',method,width,start)
        values=series(run)
        temporary=npzpath.with_suffix('.npz.partial')
        with temporary.open('wb') as f:np.savez_compressed(f,**run)
        os.replace(temporary,npzpath)
        write_csv(csvpath,[dict(step=k,time=(k+1)*d.DT,blocked=int(run['blocked'][k]),
            target_x=float(run['goals'][k+1,0]),target_y=float(run['goals'][k+1,1]),
            **{name:float(x[k]) for name,x in values.items()}) for k in range(d.STEPS)])
        rec=dict(run_id=rid,seed=seed,method=method,start=start,duration=width,npz_sha256=digest(npzpath),csv_sha256=digest(csvpath))
        records.append(rec);lookup[rid]=rec;write_json(manifest,records)
        write_json(RESULTS/'STATUS.json',dict(stage='simulating' if len(records)<65 else 'simulations_complete',completed=len(records),planned=65))
        print(f'[{len(records)}/65] {rid}, elapsed {time.perf_counter()-begin:.1f}s',flush=True)
    assert len(records)==65
    print('All 65 simulations complete.',flush=True)


if __name__=='__main__':main()

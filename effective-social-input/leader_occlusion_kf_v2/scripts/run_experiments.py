"""Exclusive-output, deterministic batch runner for the frozen v2 plan."""
from pathlib import Path
import os
for key in ['OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS']:os.environ[key]='1'
import argparse,csv,hashlib,json,time,sys,platform,shutil,traceback
from concurrent.futures import ProcessPoolExecutor,as_completed
from dataclasses import asdict
import numpy as np
import h5py
from threadpoolctl import threadpool_limits
import leader_model as model
from check_implementation import run_checks

ROOT=Path(__file__).resolve().parent

def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(8*1024**2),b''):h.update(b)
    return h.hexdigest()

def tasks():
    entries=[]
    def add(batch,seed,n,mode,T,start):
        rid=f'{batch}_seed{seed}_N{n}_{mode}_T{T}_start{start}'
        entries.append(dict(run_id=rid,batch=batch,seed=seed,N=n,mode=mode,T=T,start=start,
            durationSteps=round(T/.3),data=f'runs/{batch}/{rid}.h5',metrics=f'runs/{batch}/{rid}.csv'))
    for seed in model.DEVELOPMENT:
        for mode in model.ref.MODES:add('development',seed,20,mode,18,200)
    for seed in model.EVALUATION:
        for mode in model.MODES:add('main',seed,20,mode,18,200)
        for T in [6,30]:
            for mode in model.MODES[1:]:add(f'duration{T}',seed,20,mode,T,200)
        for mode in ['hold','kf-mean','kf-uncertainty']:add('early',seed,20,mode,18,140)
        for n in [40,60]:
            for mode in ['fresh','hold','kf-mean']:add(f'scale{n}',seed,n,mode,18,200)
    assert len(entries)==858 and len({r['run_id'] for r in entries})==858
    return entries

def dataset(group,key,a):
    a=np.asarray(a)
    if a.size==0:group.create_dataset(key,data=a);return
    opts=dict(compression='gzip',compression_opts=1,shuffle=True)
    if a.ndim>=2 and a.shape[0]>=800:opts['chunks']=(min(20,a.shape[0]),*a.shape[1:])
    group.create_dataset(key,data=a,**opts)

def save_dict(group,data):
    for k,a in data.items():dataset(group,k,a)

def worker(row,out):
    timer=time.perf_counter();out=Path(out);cfg=model.config(row['N'])
    inp=out/'inputs'/f'seed{row["seed"]}_N{row["N"]}.npz'
    with np.load(inp) as f:initial=tuple(f[k] for k in ['xy','heading','alpha','u','noise'])
    with threadpool_limits(limits=1):result=model.simulate(row['seed'],cfg,row['mode'],row['durationSteps'],row['start'],initial)
    hp=out/row['data'];cp=out/row['metrics'];hp.parent.mkdir(parents=True,exist_ok=True)
    with h5py.File(hp,'x') as f:
        f.attrs['config_json']=json.dumps(asdict(cfg));f.attrs['run_json']=json.dumps(row)
        f.attrs['input_path']=str(inp.relative_to(out));f.attrs['input_sha256']=row['input_sha256']
        for name in ['states','edges','agent','checkpoint','final']:save_dict(f.create_group(name),result[name])
        for name in ['alpha','selected','hidden','edge_i','edge_j','fixed_nearest']:dataset(f,name,result[name])
    with cp.open('x',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(result['rows'][0]));w.writeheader();w.writerows(result['rows'])
    return row|dict(status='COMPLETE',seconds=time.perf_counter()-timer,dataBytes=hp.stat().st_size,
        data_sha256=sha(hp),metrics_sha256=sha(cp))

def run(args):
    out=ROOT/'results'/args.name;out.mkdir(parents=True,exist_ok=False)
    (out/'inputs').mkdir();(out/'source_snapshot').mkdir();(out/'checks').mkdir()
    for p in ROOT.glob('*.py'):shutil.copy2(p,out/'source_snapshot'/p.name)
    for name in ['protocol_original.md','protocol_v2.md','EXECUTION_PROTOCOL.md']:
        if (ROOT/name).exists():shutil.copy2(ROOT/name,out/'source_snapshot'/name)
    manifest=tasks();source={p.name:sha(p) for p in (out/'source_snapshot').iterdir()}
    metadata=dict(status='PREPARED',sourceSHA256=source,main_config=asdict(model.config()),
        developmentSeeds=model.DEVELOPMENT,evaluationSeeds=model.EVALUATION,expectedRuns=858,workers=args.workers,
        RNG='NumPy Generator(PCG64)',python=sys.version,numpy=np.__version__,h5py=h5py.__version__,platform=platform.platform(),
        samplingOrder=['position[20,2]','heading[20]','neural[100,20]','noise[800,20,20]'],
        seedStatus='Reserved new evaluation IDs; absent from existing local scientific run manifests checked before launch',
        primaryComparison='main N20 T18 start200: occlusion GO kf-mean minus hold',
        exclusions='none; all scientific runs preserved',linearTrendWindow=10,parameterTuning='none; original fixed parameters',
        gateRules=dict(freshLeaderFinalErrorMaxDeg=10,freshFollowerFinalErrorMaxDeg=15,
            freshFollowerFinalAngleMinDeg=30,dropOcclusionErrorIncreaseMinDeg=3),
        mainEvidence='30 evaluation seeds; development results not pooled',BLASThreads=1)
    (out/'metadata.json').write_text(json.dumps(metadata,ensure_ascii=False,indent=2),encoding='utf-8')
    with threadpool_limits(limits=1):run_checks(out/'checks/preflight')
    inputs=[]
    for seed,n in sorted({(r['seed'],r['N']) for r in manifest}):
        cfg=model.config(n);initial=model.initialize(seed,cfg);path=out/'inputs'/f'seed{seed}_N{n}.npz'
        np.savez_compressed(path,**dict(zip(['xy','heading','alpha','u','noise'],initial)),
            selected=model.topology(cfg)[0],hidden=model.topology(cfg)[1])
        inputs.append(dict(seed=seed,N=n,path=str(path.relative_to(out)),sha256=sha(path)))
    ih={(r['seed'],r['N']):r['sha256'] for r in inputs}
    for r in manifest:r['input_sha256']=ih[(r['seed'],r['N'])]
    (out/'input_manifest.json').write_text(json.dumps(inputs,indent=2),encoding='utf-8')
    (out/'planned_manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
    records=[];failures=[];timer=time.perf_counter()
    with (out/'run.log').open('x',encoding='utf-8',buffering=1) as log:
        log.write('START; full v2 858-run plan; no straight mechanism controls\n')
        with ProcessPoolExecutor(max_workers=args.workers) as pool:
            futures={pool.submit(worker,r,str(out)):r for r in manifest}
            for future in as_completed(futures):
                row=futures[future]
                try:r=future.result();records.append(r)
                except Exception:
                    failure=row|dict(status='FAILED',error=traceback.format_exc());failures.append(failure);log.write(json.dumps(failure)+'\n')
                completed=len(records)+len(failures)
                if completed%10==0 or completed==858:
                    msg=f'{completed}/858 complete; failures={len(failures)}; elapsed={time.perf_counter()-timer:.1f}s'
                    print(msg,flush=True);log.write(msg+'\n')
                    (out/'progress.json').write_text(json.dumps(dict(completed=completed,successful=len(records),failed=len(failures),elapsed=time.perf_counter()-timer)),encoding='utf-8')
                with (out/'completed_manifest.json').open('w',encoding='utf-8') as f:json.dump(records,f,indent=2)
    if failures:
        (out/'failures.json').write_text(json.dumps(failures,indent=2),encoding='utf-8')
        raise RuntimeError(f'{len(failures)} scientific runs failed; partial data retained')
    records.sort(key=lambda r:r['run_id'])
    with (out/'RUN_MANIFEST.csv').open('x',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(records[0]));w.writeheader();w.writerows(records)
    audit=dict(status='SIMULATION_COMPLETE',expectedRuns=858,completedRuns=len(records),developmentRuns=18,evaluationRuns=840,
        elapsedSeconds=time.perf_counter()-timer,workers=args.workers,excludedSeeds=[],independentEvaluationSeeds=30,
        analysis='pending independent verification')
    (out/'audit.json').write_text(json.dumps(audit,indent=2),encoding='utf-8')
    print(json.dumps(audit),flush=True)
    return out

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--name',default='full_v2_run1');p.add_argument('--workers',type=int,default=6)
    args=p.parse_args();run(args)

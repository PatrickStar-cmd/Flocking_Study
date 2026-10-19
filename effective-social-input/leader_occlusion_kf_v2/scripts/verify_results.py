"""Independent position, metric, input and integrity checks of every scientific run."""
from pathlib import Path
import json,hashlib,time,csv
import numpy as np
import pandas as pd
import h5py
from scipy.special import ive

ROOT=Path(__file__).resolve().parent
def wrap(x):return np.angle(np.exp(1j*x))
def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(8*1024**2),b''):h.update(b)
    return h.hexdigest()

def verify(out):
    out=Path(out);dest=out/'checks/independent';dest.mkdir(parents=True,exist_ok=True)
    records=pd.read_csv(out/'RUN_MANIFEST.csv');assert len(records)==858
    assert not records.run_id.duplicated().any() and (records.status=='COMPLETE').all()
    metadata=json.loads((out/'metadata.json').read_text())
    for name,h in metadata['sourceSHA256'].items():assert sha(out/'source_snapshot'/name)==h
    for r in json.loads((out/'input_manifest.json').read_text()):assert sha(out/r['path'])==r['sha256']
    sources={p.name:sha(p) for p in ROOT.glob('*.py')}
    # Simulation versions must not drift after launch.
    for name in ['leader_model.py','reference_run_experiment.py','run_experiments.py','check_implementation.py']:
        assert sources[name]==metadata['sourceSHA256'][name],name
    checks=[];prefix_cache={};timer=time.perf_counter()
    order=records.copy();order['fresh_first']=(order['mode']!='fresh').astype(int)
    order=order.sort_values(['N','seed','fresh_first','batch','mode'])
    for k,r in enumerate(order.itertuples()):
        hp=out/r.data;csvpath=out/r.metrics
        assert sha(hp)==r.data_sha256 and sha(csvpath)==r.metrics_sha256,r.run_id
        d=pd.read_csv(csvpath);assert len(d)==800 and list(d.step)==list(range(800))
        start=int(r.start);duration=int(r.durationSteps);L=3 if r.N==20 else int(.15*r.N);F=r.N-L
        with h5py.File(hp) as f:
            p=f['states/xy'][:];heading=f['states/headings'][:];v=f['states/velocity'][:]
            assert p.shape==(801,r.N,2) and f['states/u'].shape==(801,100,r.N)
            assert np.isfinite(p).all() and np.isfinite(v).all() and np.all((p>=0)&(p<1000))
            displacement=np.diff(p,axis=0);displacement-=1000*np.round(displacement/1000)
            velocityerr=np.max(np.abs(displacement/.3-v));assert velocityerr<1e-10
            speed=np.linalg.norm(v,axis=2);unit=np.divide(v,speed[:,:,None],out=np.zeros_like(v),where=speed[:,:,None]>1e-12)
            go=np.linalg.norm(unit.mean(axis=1),axis=1);gof=np.linalg.norm(unit[:,L:].mean(axis=1),axis=1)
            E=np.degrees(np.abs(wrap(heading[1:,L:]-np.radians(d.command_deg.to_numpy())[:,None]))).mean(axis=1)
            maxerr=max(np.max(np.abs(go-d.GO)),np.max(np.abs(gof-d.GO_followers)),np.max(np.abs(E-d.follower_turn_error_deg)))
            assert maxerr<1e-10
            all_delta=p[1:,:,None,:]-p[1:,None,:,:];all_delta-=1000*np.round(all_delta/1000)
            dist=np.linalg.norm(all_delta,axis=3);pair=dist[:,np.triu_indices(r.N,1)[0],np.triu_indices(r.N,1)[1]].mean(axis=1)
            prefixkey=(r.seed,r.N)
            before=p[start];delta=before[:,None,:]-before[None,:,:];delta-=1000*np.round(delta/1000)
            initialdist=np.linalg.norm(delta,axis=2);np.fill_diagonal(initialdist,np.inf)
            ids=np.argsort(initialdist,axis=1)[:,:5];assert np.array_equal(ids,f['fixed_nearest'][:])
            fixed=dist[:,np.arange(r.N)[:,None],ids].mean(axis=(1,2))
            derr=max(np.max(np.abs(pair-d.pair_distance)),np.max(np.abs(fixed[start:]-d.fixed5_distance.to_numpy()[start:])))
            assert derr<1e-10 and d.fixed5_distance.iloc[:start].isna().all()
            selected=f['selected'][:].astype(bool);hidden=f['hidden'][:].astype(bool);ii=f['edge_i'][:];jj=f['edge_j'][:]
            assert selected.sum()==15*r.N and hidden.sum()==3*F and np.all(selected.sum(1)==15)
            observed=f['edges/observed'][:];fresh=f['edges/fresh'][:];active=f['edges/active'][:];weights=f['edges/weight'][:]
            masked=hidden[ii,jj];event=np.zeros(800,bool)
            if r.mode!='fresh':event[start:start+duration]=True
            missing=event[:,None]&masked[None]
            assert np.array_equal(~fresh,missing) and np.isnan(observed[missing]).all() and np.isfinite(observed[fresh]).all()
            assert np.array_equal(fresh.sum(axis=1),d.fresh_edges)
            assert np.array_equal((active&missing).sum(axis=1),d.hidden_active_edges)
            center=f['edges/center'][:];P11=f['edges/P11'][:];P12=f['edges/P12'][:];P22=f['edges/P22'][:]
            assert P11.min()>=0 and P22.min()>=0 and (P11*P22-P12*P12).min()>-1e-12
            hstart=heading[:-1,ii];factor=np.exp(np.cos(center-hstart)-1)
            if r.mode=='kf-uncertainty':
                cmask=missing&active;delta=center[cmask]-hstart[cmask];V=P11[cmask]
                harmonics=np.arange(1,24);factor[cmask]=ive(0,1)+2*np.sum(ive(harmonics,1)*np.exp(-.5*V[:,None]*harmonics**2)*np.cos(delta[:,None]*harmonics),axis=1)
            empirical=np.zeros((800,r.N))
            for i in range(r.N):empirical[:,i]=(weights[:,ii==i]*factor[:,ii==i]).sum(axis=1)
            truthdelta=p[:-1,jj]-p[:-1,ii];truthdelta-=1000*np.round(truthdelta/1000)
            truth=np.arctan2(truthdelta[:,:,1],truthdelta[:,:,0])%(2*np.pi)
            assert np.max(np.abs(wrap(truth-f['edges/true'][:])))<1e-12
            oracle=np.zeros_like(empirical)
            for i in range(r.N):oracle[:,i]=.012*np.exp(np.cos(truth[:,ii==i]-heading[:-1,i,None])-1).sum(axis=1)
            serr=max(np.max(np.abs(empirical[:,L:].mean(axis=1)-d.S)),
                np.max(np.abs(np.abs(empirical[:,L:]-oracle[:,L:]).mean(axis=1)-d.S_abs_error)))
            assert serr<1e-10
            # Rebuild social field at five predeclared rows without calling input().
            alpha=f['alpha'][:];field_err=0.
            for step in sorted(set([0,start,start+duration-1,start+duration,799])):
                if step>=800:continue
                g=np.exp(-.5*(wrap(alpha[ii]-center[step,:,None])/.4)**2)
                if r.mode=='kf-uncertainty':
                    cmask=missing[step]&active[step];freq=np.fft.fftfreq(100)*100
                    g[cmask]=np.maximum(0,np.fft.ifft(np.fft.fft(g[cmask],axis=1)*np.exp(-.5*P11[step,cmask,None]*freq**2),axis=1).real)
                rebuilt=np.zeros((100,r.N))
                for i in range(r.N):rebuilt[:,i]=(weights[step,ii==i,None]*g[ii==i]).sum(axis=0)
                field_err=max(field_err,float(np.max(np.abs(rebuilt-f['states/social_field'][step]))))
            assert field_err<1e-10
            if r.mode=='fresh':
                prefix_cache[prefixkey]=dict(p=p[:201],h=heading[:201],u=f['states/u'][:201],
                    checkpoint={name:f['checkpoint/'+name][:] for name in ['filter_x','filter_P','initialized','last_bearing','age']})
            else:
                ref=prefix_cache[prefixkey];assert np.array_equal(p[:start+1],ref['p'][:start+1])
                assert np.array_equal(heading[:start+1],ref['h'][:start+1]) and np.array_equal(f['states/u'][:start+1],ref['u'][:start+1])
                if start==200:
                    for name,a in ref['checkpoint'].items():assert np.array_equal(f['checkpoint/'+name][:],a,equal_nan=True)
            A=np.zeros((800,r.N))
            for i in range(r.N):A[:,i]=weights[:,ii==i].sum(axis=1)
            expA=.144 if r.mode=='drop' else .18
            if event.any():assert np.max(np.abs(A[event,L:]-expA))<1e-12
            assert np.max(np.abs(A[~event]-.18))<1e-12
            rec=start+duration
            if r.mode!='fresh':assert np.allclose(center[rec],observed[rec],atol=0,rtol=0) and np.all(fresh[rec])
        checks.append(dict(run_id=r.run_id,status='PASS',GOAndTaskMaxDifference=float(maxerr),distanceMaxDifference=float(derr),
            SMaxDifference=float(serr),fieldMaxDifference=float(field_err),velocityFromPositionMaxDifference=float(velocityerr)))
        if (k+1)%50==0:print(f'Independent checks {k+1}/858 elapsed {time.perf_counter()-timer:.1f}s',flush=True)
    pd.DataFrame(checks).to_csv(dest/'all_run_checks.csv',index=False,encoding='utf-8-sig')
    result=dict(status='PASS',runs=858,checksPerRun=['GO','GO_followers','task_error','pair_distance','fixed5_distance','S','S_abs_error',
        'social_field','shared_prefix','hidden_NaN','fresh_count','weights','KF_PSD','recovery','SHA256'],
        independentSeeds=30,developmentRuns=18,evaluationRuns=840,
        maximumMetricDifference=max(max(c['GOAndTaskMaxDifference'],c['distanceMaxDifference'],c['SMaxDifference'],c['fieldMaxDifference']) for c in checks),
        seconds=time.perf_counter()-timer,analysisSources=sources)
    (dest/'validation.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    print(json.dumps({k:v for k,v in result.items() if k!='analysisSources'}),flush=True)
    return result

if __name__=='__main__':
    import sys
    verify(ROOT/(sys.argv[1] if len(sys.argv)>1 else 'results/full_v2_run1'))

"""All predeclared candidates; truth labels never enter the filter recursion.

Raw branch sequences are exact float-bit XORs against shared full-fresh fits,
including all recovery400. Pre200 is the shared fresh sequence. Nothing is
removed because a deployment variance/age threshold rejects it.
"""
import argparse,gzip,hashlib,json,sys,time
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.io import loadmat
from engine import run,archive_branch
from geometry import geometry,classify_window,fresh_labels

STARTS=np.arange(1201,1682,20)
DURATIONS=[20,120]
GROUPS=['all','G-clean','G-switch','G-close','G-ambiguous','clean-slow','clean-mid','clean-fast']


def digest(path):return hashlib.file_digest(path.open('rb'),'sha256').hexdigest()


def jwrite(path,obj):
    assert not path.exists();path.write_text(json.dumps(obj,indent=2),encoding='utf-8')


def source(path):
    x=loadmat(path,simplify_cells=True)
    assert x['observed_bearing'].shape==x['true_bearing'].shape==(2200,140)
    assert len(set(zip(x['observer_id'],x['target_id'])))==140
    return x


def prepare(stage,frozen=None):
    manifest=json.loads((stage/'input_manifest.json').read_text(encoding='utf-8'))
    assert manifest['status']=='COMPLETED'
    dest=stage/'labels';assert not dest.exists();dest.mkdir()
    records=[];rates=[];cases=[]
    for entry in manifest['inputs']:
        path=stage/entry['source'];assert digest(path)==entry['sha256']
        x=source(path);g=geometry(x)
        stem=path.stem
        for D in DURATIONS:
            labels=[];speed=[]
            for t0 in STARTS:
                lb,rt=classify_window(g,int(t0),D)
                labels.append(lb);speed.append(rt)
                cases.append({'seed':int(entry['seed']),'scenario':entry['scenario'],'source':entry['source'],'t0':int(t0),'D':D})
            labels=np.array(labels,dtype='U16');speed=np.array(speed)
            rates.extend(speed[(labels=='G-clean')&np.isfinite(speed)].tolist())
            p=dest/f'{stem}_D{D}.npz'
            np.savez_compressed(p,labels=labels,rates=speed,starts=STARTS)
            records.append({'file':p.name,'sha256':digest(p),'seed':int(entry['seed']),'scenario':entry['scenario'],'D':D})
    if frozen is None:
        assert len(rates)>0,'No clean development samples'
        cuts=np.quantile(rates,[1/3,2/3]).tolist()
        decision={'status':'FROZEN_GEOMETRY','rate_cutpoints':cuts,'development_clean_edge_windows':len(rates),
                  'rule':'pooled development truth-only clean rates; slow<=q1,mid<=q2,fast>q2',
                  'labels':records,'starts':STARTS.tolist(),'durations':DURATIONS}
        jwrite(stage/'GEOMETRY_FROZEN.json',decision)
    else:
        cuts=frozen['rate_cutpoints']
        jwrite(stage/'GEOMETRY_APPLIED.json',{'frozen_decision_sha256':digest(frozen['_file']),'rate_cutpoints':cuts,'labels':records})
    pd.DataFrame(cases).to_csv(stage/'cases.csv',index=False)
    return manifest,cuts


def group_masks(labels,rates,cuts):
    clean=labels=='G-clean'
    return {'all':np.ones(140,bool),**{g:labels==g for g in GROUPS[1:5]},
        'clean-slow':clean&(rates<=cuts[0]),'clean-mid':clean&(rates>cuts[0])&(rates<=cuts[1]),
        'clean-fast':clean&(rates>cuts[1])}


def metrics(e,p,age,hold=False):
    absolute=np.abs(e);finite=np.isfinite(absolute)
    assert finite.all()
    out={'rmse':float(np.sqrt(np.mean(e*e))),'mae':float(absolute.mean()),
         'p90':float(np.quantile(absolute,.90)),'p95':float(np.quantile(absolute,.95)),
         'p99':float(np.quantile(absolute,.99)),'max_abs':float(absolute.max()),
         'extreme_fraction':float((absolute>1).mean()),'edge_times':e.size}
    if hold:
        out.update(mean_P=np.nan,coverage=np.nan,ratio_median=np.nan,ratio_p95=np.nan,ratio_p99=np.nan)
        reject=age>120
    else:
        ratio=e*e/p[:,None]
        out.update(mean_P=float(p.mean()),coverage=float((absolute<=1.96*np.sqrt(p[:,None])).mean()),
            ratio_median=float(np.median(ratio)),ratio_p95=float(np.quantile(ratio,.95)),ratio_p99=float(np.quantile(ratio,.99)))
        reject=(p>1)|(age>120)
    first=np.flatnonzero(reject)
    out['age_cutoff_fraction']=float((age>120).mean())
    out['P_cutoff_fraction']=float((p>1).mean()) if not hold else np.nan
    out['cutoff_fraction']=float(reject.mean());out['first_cutoff_step']=int(first[0]+1) if len(first) else np.nan
    longest=0;current=0
    for flag in reject:
        current=current+1 if flag else 0;longest=max(longest,current)
    out['longest_cutoff_steps']=longest
    return out


def execute(stage,validation=False,frozen_path=None):
    begin=time.time()
    frozen=None
    if validation:
        assert frozen_path is not None
        frozen=json.loads(frozen_path.read_text(encoding='utf-8'));frozen['_file']=frozen_path
        assert frozen['q_star'] is not None
    manifest,cuts=prepare(stage,frozen)
    raw=stage/'raw/replay';assert not raw.exists();raw.mkdir()
    analysis=stage/'analysis';assert not analysis.exists();analysis.mkdir()
    models=[('Hold',None,None),('CV_q1e-06',1e-6,None),('CV_q1e-05',1e-5,None),('CV_q0.0001',1e-4,None),('CV_q0.001',1e-3,None),('CV_q0.01',1e-2,None),('GM_default',1e-5,10)]
    if validation:
        models=[('Hold',None,None),('CV_q0.0001',1e-4,None),('GM_default',1e-5,10)]
        if frozen['q_star']!=1e-4:models.append((f"CV_q{frozen['q_star']:g}",frozen['q_star'],None))
    rows=[];coverage=[];nis_rows=[];bins=[];raw_manifest=[]
    for number,entry in enumerate(manifest['inputs'],1):
        path=stage/entry['source'];x=source(path);obs=x['observed_bearing'];truth=x['true_bearing']
        g=geometry(x);flabels=fresh_labels(g);stem=path.stem
        label_data={D:dict(np.load(stage/'labels'/f'{stem}_D{D}.npz',allow_pickle=False)) for D in DURATIONS}
        for method,q,tau in models:
            ishold=method=='Hold'
            fresh=run(obs,np.zeros(2200,bool),q=q or 0,tau=tau,hold=ishold)
            fresh_file=raw/f'{stem}_{method}_fresh.npz'
            np.savez_compressed(fresh_file,**fresh)
            raw_manifest.append({'file':fresh_file.name,'sha256':digest(fresh_file),'kind':'shared-fresh','source':entry['source'],'method':method})
            if not ishold:
                for layer in ['all','switch_event','switch_near2','G-clean','G-close','G-ambiguous']:
                    mask=np.ones((2200,140),bool) if layer=='all' else flabels==layer
                    values=fresh['nis'][:,0][mask];values=values[np.isfinite(values)]
                    if len(values):
                        nis_rows.append({'seed':int(entry['seed']),'scenario':entry['scenario'],'method':method,'geometry':layer,
                            'n':len(values),'mean':float(values.mean()),'median':float(np.median(values)),
                            **{f'p{v}':float(np.quantile(values,v/100)) for v in [90,95,99]}})
            for D in DURATIONS:
                K=D+400;indices=STARTS[None,:]-1+np.arange(K)[:,None]
                initial={key:fresh[key][STARTS-2,0] for key in ['theta','omega','P']}
                # Hold's last actual observation is already fresh['theta'].
                branch=run(obs[indices],np.arange(K)<D,q=q or 0,tau=tau,initial=initial,hold=ishold)
                f=raw/f'{stem}_{method}_D{D}.npz'
                archive_branch(f,branch,fresh,indices,STARTS,initial,D=D,method=method,source_sha256=entry['sha256'])
                raw_manifest.append({'file':f.name,'sha256':digest(f),'kind':'XOR-branch','fresh_file':fresh_file.name,'source':entry['source'],'method':method,'D':D})
                labels=label_data[D]['labels'];rates=label_data[D]['rates']
                for c,t0 in enumerate(STARTS):
                    error=branch['theta'][:D,c]-truth[indices[:D,c]]
                    error=np.arctan2(np.sin(error),np.cos(error))
                    masks=group_masks(labels[c],rates[c],cuts)
                    p=branch['P'][:D,c,0,0];age=branch['age'][:D]
                    meta={'seed':int(entry['seed']),'scenario':entry['scenario'],'method':method,'q':q,'D':D,'t0':int(t0)}
                    for group,mask in masks.items():
                        if not mask.any():
                            rows.append({**meta,'geometry':group,'n_edges':0,'rmse':np.nan,'edge_times':0})
                            continue
                        e=error[:,mask]
                        stat=metrics(e,p,age,ishold)
                        actual_rate=g['angular_rate_clean'][indices[:D,c]][:,mask]
                        estimated_rate=branch['omega'][:D,c][:,mask]
                        valid_rate=np.isfinite(actual_rate)&np.isfinite(estimated_rate)
                        rate_error=(estimated_rate-actual_rate)[valid_rate]
                        stat.update(omega_clean_count=int(valid_rate.sum()),
                            omega_rmse=float(np.sqrt(np.mean(rate_error**2))) if len(rate_error) else np.nan,
                            omega_mae=float(np.mean(np.abs(rate_error))) if len(rate_error) else np.nan)
                        rows.append({**meta,'geometry':group,'n_edges':int(mask.sum()),**stat})
                        if not ishold:
                            for lo,hi in [(1,20),(21,60),(61,100),(101,120)]:
                                stop=min(hi,D)
                                if lo>D:continue
                                ee=e[lo-1:stop];pp=p[lo-1:stop]
                                coverage.append({**meta,'geometry':group,'age_band':f'{lo}-{hi}','n_edges':int(mask.sum()),**metrics(ee,pp,age[lo-1:stop])})
                            for lo,hi in zip([0,1e-5,1e-4,.001,.01,.1,1],[1e-5,1e-4,.001,.01,.1,1,np.inf]):
                                valid=(p>=lo)&(p<hi)
                                if valid.any():bins.append({**meta,'geometry':group,'P_lo':lo,'P_hi':hi,'n':int(valid.sum()*mask.sum()),
                                    'mean_P':float(p[valid].mean()),'mean_error_squared':float(np.mean(e[valid]**2)),
                                    'coverage':float((np.abs(e[valid])<=1.96*np.sqrt(p[valid,None])).mean())})
                # Persist manifest after each independent source/model/D; no silent skip.
                manifest_path=raw/'manifest.progress.json'
                manifest_path.write_text(json.dumps({'status':'RUNNING','files':raw_manifest},indent=2),encoding='utf-8')
            print(f'PASS {number}/{len(manifest["inputs"])} {stem} {method}',flush=True)
        pd.DataFrame(rows).to_csv(analysis/'window_metrics.progress.csv.gz',index=False,compression='gzip')
    frame=pd.DataFrame(rows);frame.to_csv(analysis/'window_metrics.csv.gz',index=False,compression='gzip')
    pd.DataFrame(coverage).to_csv(analysis/'coverage_by_age.csv.gz',index=False,compression='gzip')
    pd.DataFrame(bins).to_csv(analysis/'p_bins.csv.gz',index=False,compression='gzip')
    pd.DataFrame(nis_rows).to_csv(analysis/'nis_tail_by_geometry.csv',index=False)
    valid=frame[frame.n_edges>0]
    def aggregate_rows(data,keys):
        return data.groupby(keys,as_index=False).agg(
            rmse=('rmse','mean'),mae=('mae','mean'),mean_window_p95=('p95','mean'),mean_window_p99=('p99','mean'),
            mean_window_coverage=('coverage','mean'),mean_cutoff_fraction=('cutoff_fraction','mean'),
            edge_windows=('n_edges','sum'),edge_times=('edge_times','sum'),n_windows=('t0','size'))
    combined=aggregate_rows(valid,['seed','geometry','D','method']);combined['scenario']='combined'
    separated=aggregate_rows(valid,['seed','scenario','geometry','D','method'])
    aggregate=pd.concat([combined,separated],ignore_index=True)
    aggregate.to_csv(analysis/'seed_metrics.csv',index=False)
    valid.groupby(['geometry','D','method'],as_index=False).agg(n_seeds=('seed','nunique'),n_windows=('t0','size'),edge_windows=('n_edges','sum'),edge_times=('edge_times','sum')).to_csv(analysis/'geometry_counts.csv',index=False)
    jwrite(raw/'manifest.json',{'status':'COMPLETED','files':raw_manifest,'model_count':len(models),'input_count':len(manifest['inputs']),
        'duration_seconds':time.time()-begin,'archive':'full fresh + exact XOR branch; recover_branch in engine.py; pre200 from fresh'})
    jwrite(analysis/'audit.json',{'status':'PASS','input_count':len(manifest['inputs']),'seed_count':len(manifest['seeds']),
        'candidate_starts':STARTS.tolist(),'model_count':len(models),'window_metric_rows':len(frame),'raw_files':len(raw_manifest),
        'rate_cutpoints':cuts,'all_geometry_retained':True,'thresholds_do_not_remove_errors':True})
    print('COMPLETE',stage,flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('stage',type=Path);parser.add_argument('--frozen',type=Path)
    args=parser.parse_args();execute(args.stage,args.frozen is not None,args.frozen)

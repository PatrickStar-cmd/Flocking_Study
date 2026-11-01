"""Independent raw-MAT, geometry, covariance, source-hash and CSV audit."""
import argparse
import hashlib
import json
import statistics
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.io import loadmat
from scipy.stats import t as student_t

from independent_review import verify_geometry


METRICS = {'GO_mean':'globalOrder', 'LO_mean':'localOrder',
           'pair_norm_mean':'meanPairDistanceNorm', 'nearest_norm_mean':'meanNearestDistanceNorm',
           'active_neighbors_mean':'meanActiveNeighbors', 'input_mass_mean':'meanInputMass',
           'bump_concentration_mean':'meanBumpConcentration', 'speed_mean':'meanSpeed'}


def sha(path):
    digest=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(1024*1024), b''):
            digest.update(block)
    return digest.hexdigest()


def verify(run_dir, analysis_dir, output):
    run_dir, analysis_dir, output=map(Path,[run_dir,analysis_dir,output])
    manifest=json.loads((run_dir/'manifest.json').read_text(encoding='utf-8'))
    assert manifest['status']=='COMPLETED' and manifest['failed_count']==0
    entries=list(manifest['runs'].values())
    assert len(entries)==45 and manifest['completed_count']==45
    assert len({entry['mat_file'] for entry in entries})==45
    for source in json.loads((run_dir/'source_hashes.json').read_text(encoding='utf-8')):
        assert sha(run_dir/source['relativePath'])==source['sha256']
    summary=pd.read_csv(analysis_dir/'summary.csv')
    contrasts=pd.read_csv(analysis_dir/'contrasts.csv')
    grouped=pd.read_csv(analysis_dir/'contrasts_group_summary.csv')
    diagnostics=pd.read_csv(analysis_dir/'estimator_diagnostics.csv')
    assert len(summary)==200 and len(contrasts)==210 and len(diagnostics)==2800
    means={}; raw_reports=[]; maximum_csv_error=0.; gate_counts=[]
    cfg_keys=('NAgents','NNeurons','NSteps','seed','occlusionDurationSteps')
    for number,entry in enumerate(entries,1):
        path=run_dir/entry['mat_file']; assert sha(path)==entry['sha256']
        content=loadmat(path,simplify_cells=True); result=content['result']; cfg=result['config']
        seed=int(cfg['seed']); duration=int(cfg['occlusionDurationSteps'])
        mode='fresh' if duration==0 else cfg['predictionMode']
        assert int(cfg['NAgents'])==20 and int(cfg['NNeurons'])==100 and int(cfg['NSteps'])==2200
        geometry=verify_geometry(path,content); assert geometry['pass']
        records=result['observationRecords']; codes=np.stack([r['cutoffReason'] for r in records])
        p_cut=int(np.isin(codes,[22,23]).sum()); age_cut=int(np.isin(codes,[21,23]).sum())
        expected_p=2800 if duration==120 and mode in ('kf-mean','kalman') else 0
        assert p_cut==expected_p and age_cut==0
        gate_counts.append({'seed':seed,'mode':mode,'D':duration,'P_cutoff_edge_steps':p_cut,'age_cutoff_edge_steps':age_cut})
        assert result['preOcclusionCheckpoint']['stateAfterStep']==1400
        # Reconstruct neural coherence from post-update u and fixed preferred angles.
        pre_u=result['uRecord']; alpha=result['finalPreferredDirections'].T
        post_u=np.concatenate([pre_u[:,:,1:],result['finalNeuralState'][:,:,None]],axis=2)
        positive=np.maximum(np.tanh(float(cfg['beta'])*post_u),0)
        cx=(positive*np.cos(alpha)[:,:,None]).sum(axis=0)
        cy=(positive*np.sin(alpha)[:,:,None]).sum(axis=0)
        bump=(np.hypot(cx,cy)/(positive.sum(axis=0)+np.finfo(float).eps)).mean(axis=0)
        bump_error=float(np.max(np.abs(bump-result['timeSeries']['meanBumpConcentration'])))
        assert bump_error<1e-12
        own_rows=summary[(summary.seed==seed)&(summary['mode']==mode)&(summary.D_steps==duration)]
        assert len(own_rows)==(8 if duration==0 else 4)
        for _,row in own_rows.iterrows():
            eval_d=int(row.eval_D_steps)
            bounds={'pre':(1201,1400),'occ':(1401,1400+eval_d),
                    'rec':(1401+eval_d,1800+eval_d),'common':(1401,1920)}
            start,end=bounds[row.window]
            assert int(row.window_start_step)==start and int(row.window_end_step)==end
            assert int(row.window_length_steps)==end-start+1
            for column,field in METRICS.items():
                value=float(np.mean(result['timeSeries'][field][start-1:end]))
                error=abs(value-float(row[column])); maximum_csv_error=max(maximum_csv_error,error)
                assert error<1e-12,(seed,mode,duration,row.window,column,error)
                means[(seed,mode,eval_d,row.window,column)]=value
        raw_reports.append({'seed':seed,'mode':mode,'D':duration,
                            'geometry_errors':geometry['max_absolute_errors'],
                            'bump_max_error':bump_error,'raw_sha256':entry['sha256']})
        print(f'RAW_AUDIT {number}/45 seed={seed} {mode} D={duration}',flush=True)
        del content,result,records,codes,pre_u,post_u,positive
    for _,row in contrasts.iterrows():
        key=(int(row.seed),int(row.D),row.window)
        for metric,column,sign in [('benefitPair','pair_norm_mean',-1),('benefitGO','GO_mean',1),('benefitLO','LO_mean',1)]:
            baseline=means[(key[0],row.baseline,key[1],key[2],column)]
            method=means[(key[0],row.method,key[1],key[2],column)]
            assert abs(sign*(method-baseline)-float(row[metric]))<1e-12
    for _,row in grouped.iterrows():
        subset=contrasts[(contrasts.D==row.D)&(contrasts.window==row.window)&
                         (contrasts.baseline==row.baseline)&(contrasts.method==row.method)]
        values=subset[row.metric].tolist(); assert len(values)==5
        mean=statistics.mean(values); sem=statistics.stdev(values)/np.sqrt(5)
        half=student_t.ppf(.975,4)*sem
        assert abs(mean-row['mean'])<1e-12
        assert abs(mean-half-row.ci95_t_lower)<1e-12
        assert abs(mean+half-row.ci95_t_upper)<1e-12
    report={'status':'PASS','raw_files':45,'summary_rows':len(summary),
            'contrast_rows':len(contrasts),'diagnostic_rows':len(diagnostics),
            'group_rows':len(grouped),'max_csv_error':maximum_csv_error,
            'raw_reports':raw_reports,'gate_counts':gate_counts,
            'scope':'Raw SHA/source SHA, periodic geometry GO/LO/pair/nearest/speed, neural bump, window means, paired contrasts, Student-t intervals'}
    with output.open('x',encoding='utf-8') as stream:
        json.dump(report,stream,ensure_ascii=False,indent=2)
    print('PASS independent package audit',flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(); parser.add_argument('run_dir'); parser.add_argument('analysis_dir'); parser.add_argument('output')
    args=parser.parse_args(); verify(args.run_dir,args.analysis_dir,args.output)

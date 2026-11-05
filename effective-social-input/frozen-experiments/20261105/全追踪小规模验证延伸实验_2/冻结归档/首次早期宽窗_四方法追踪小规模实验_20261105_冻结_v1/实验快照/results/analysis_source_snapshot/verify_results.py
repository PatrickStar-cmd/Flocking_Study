"""Independent raw motion, scalar KF, field replay and statistical audit."""
import bootstrap
import math
import statistics as st
from pathlib import Path
import numpy as np
from common import ROOT,RESULTS,read_json,read_csv,write_json,digest
from audit_math import scalar_filter,control_replay,circular,topology


def main():
    cfg=read_json(ROOT/'experiment_config.json')
    n,steps,dt=cfg['n'],cfg['steps'],cfg['dt']
    seeds,starts,widths,methods=cfg['pilot_seeds'],cfg['occlusion_starts_zero_based'],cfg['durations_steps'],cfg['methods']
    assert n==12 and steps==120 and dt==.3 and cfg['kf']['q']==cfg['kf']['R']==1e-4
    hashes=read_json(RESULTS/'frozen_source_hashes.json')
    for name,sha in hashes.items():assert digest(ROOT/name)==digest(RESULTS/'source_snapshot'/name)==sha,name
    provenance=read_json(ROOT/'SOURCE_PROVENANCE.json')
    source=Path(provenance['source_folder'])
    if source.exists():
        actual={p.relative_to(source).as_posix():digest(p) for p in source.rglob('*') if p.is_file() and '__pycache__' not in p.relative_to(source).parts}
        assert actual==provenance['source_files_sha256'],'Previous experiment modified'
    assert digest(ROOT/'kf_reference_core.py')==provenance['unchanged_reference_core_sha256']
    records=read_json(RESULTS/'run_manifest.json')
    expected={(seed,method,start,width) for seed in seeds for method,start,width in
              [('fresh',0,0)]+[(m,a,w) for a in starts for w in widths for m in methods[1:]]}
    assert len(records)==65 and {(r['seed'],r['method'],r['start'],r['duration']) for r in records}==expected
    selected,hidden=topology();scales=np.ones((12,12));scales[1:,0]=20.
    times=np.arange(steps+1)*dt
    goals=np.c_[235.+.6*times,200.+4*np.sin(2*np.pi*times/24.)]
    leader=np.empty((steps+1,2));leader[0]=[200.,200.]
    for k in range(steps):
        delta=goals[k]-leader[k];distance=math.hypot(*delta)
        leader[k+1]=leader[k]+dt*(delta/distance*min(1.5,distance/dt) if distance>1e-12 else np.zeros(2))
    raw={};fresh={};maximum_metric=maximum_motion=maximum_filter=maximum_cov=maximum_center=maximum_control=0.
    samples=0
    for number,rec in enumerate(records,1):
        rid=rec['run_id'];seed=rec['seed'];method=rec['method'];start=rec['start'];width=rec['duration']
        assert rid==f'{seed}_S{start:03d}_D{width:03d}_{method}'
        assert digest(RESULTS/'trajectories'/(rid+'.npz'))==rec['npz_sha256']
        assert digest(RESULTS/'timeseries'/(rid+'.csv'))==rec['csv_sha256']
        with np.load(RESULTS/'trajectories'/(rid+'.npz')) as z:r={k:z[k] for k in z.files}
        assert int(r['seed'])==seed and int(r['duration'])==width and int(r['occlusion_start'])==start
        assert str(r['method'])==method and str(r['profile'])=='moving_moderate'
        assert r['positions'].shape==(121,12,2) and np.isfinite(r['positions']).all()
        assert np.all(r['positions']>=0) and np.all(r['positions']<1000)
        assert np.array_equal(r['selected'],selected) and np.array_equal(r['hidden'],hidden)
        assert np.array_equal(r['edge_scales'],scales)
        assert np.max(np.abs(r['goals']-goals))<1e-12
        assert np.max(np.abs(r['positions'][:,0]-leader))<1e-9
        assert np.linalg.norm(r['velocities'][:,0],axis=1).max()<=1.5+1e-12
        assert np.array_equal(r['initial_u'][:,0],r['final_u'][:,0])
        blocked=np.array([start<=k<start+width and method!='fresh' for k in range(steps)])
        assert np.array_equal(r['blocked'],blocked) and int(blocked.sum())==width
        available=np.broadcast_to(selected,(steps,12,12)).copy();available[blocked]&=~hidden
        assert np.array_equal(available,r['fresh'])
        assert np.isnan(r['observations'][~available]).all()
        relative=r['positions'][:-1,None,:,:]-r['positions'][:-1,:,None,:]
        assert np.abs(relative).max()<500
        true=np.arctan2(relative[...,1],relative[...,0])
        expected_observed=(true+r['noise'])%(2*np.pi)
        assert np.max(np.abs(circular(expected_observed[available]-r['observations'][available])))<1e-10
        state,cov,last=scalar_filter(r['observations'],available,selected)
        maximum_filter=max(maximum_filter,float(np.abs(circular(state[...,0]-r['filter_state'][...,0])).max()),
                           float(np.abs(state[...,1]-r['filter_state'][...,1]).max()))
        maximum_cov=max(maximum_cov,float(np.abs(cov-r['filter_covariance']).max()))
        active=available if method in ('drop','fresh') else np.broadcast_to(selected,(steps,12,12)).copy()
        if method=='kf-mean':active&=cov[...,0,0]<=1.
        assert np.array_equal(r['active'],active)
        center=np.where(available,r['observations'],last if method=='hold' else state[...,0])
        center=np.where(active,center,0.)
        maximum_center=max(maximum_center,float(np.abs(circular(center-r['centers'])).max()))
        mass=(.012*scales*active).sum(2)
        assert np.max(np.abs(mass-r['input_mass']))<1e-12
        effective=(.012*scales*active*np.exp(np.cos(center-r['headings'][:-1,:,None])-1)).sum(2)
        assert np.max(np.abs(effective-r['effective_social_input']))<1e-12
        maximum_control=max(maximum_control,control_replay(r))
        maximum_motion=max(maximum_motion,float(np.abs(np.diff(r['positions'],axis=0)-dt*r['velocities']).max()))
        angle=[];distance=[];leader_mae=[]
        for k in range(steps):
            aa=[];dd=[]
            for agent in range(1,12):
                dx,dy=goals[k+1]-r['positions'][k+1,agent]
                vx,vy=r['velocities'][k,agent]
                dist=math.hypot(float(dx),float(dy));speed=math.hypot(float(vx),float(vy))
                assert dist>1e-9 and speed>1e-12
                aa.append(abs(math.atan2(float(vx*dy-vy*dx),float(vx*dx+vy*dy)))*180/math.pi)
                dd.append(dist);samples+=1
            angle.append(st.mean(aa));distance.append(st.mean(dd))
            leader_mae.append(st.mean(abs(math.atan2(math.sin(float(center[k,i,0]-true[k,i,0])),math.cos(float(center[k,i,0]-true[k,i,0]))))*180/math.pi for i in range(1,12)))
        raw[seed,method,start,width]={'angle_deg':angle,'distance':distance,'leader_mae':leader_mae,
                                     'hidden_coverage':float(active[blocked][:,hidden].mean()) if width else 1.,
                                     'maximum_hidden_angle_variance':float(cov[blocked][:,hidden,0,0].max()) if width else 0.}
        rows=read_csv(RESULTS/'timeseries'/(rid+'.csv'));assert len(rows)==steps
        for k,row in enumerate(rows):
            assert int(row['step'])==k and int(row['blocked'])==blocked[k]
            assert abs(float(row['time'])-(k+1)*dt)<1e-12
            assert abs(float(row['target_x'])-goals[k+1,0])<1e-12 and abs(float(row['target_y'])-goals[k+1,1])<1e-12
            for metric in ('angle_deg','distance'):
                maximum_metric=max(maximum_metric,abs(float(row[metric])-raw[seed,method,start,width][metric][k]))
            assert float(row['angle_valid_fraction'])==1. and float(row['stationary_fraction'])==0.
            assert abs(float(row['leader_target_distance'])-np.linalg.norm(goals[k+1]-leader[k+1]))<1e-10
        if method=='fresh':fresh[seed]={key:r[key] for key in ('positions','velocities','headings','noise','alpha','initial_u')}
        else:
            f=fresh[seed]
            for key in ('positions','headings'):assert np.array_equal(r[key][:start+1],f[key][:start+1])
            assert np.array_equal(r['velocities'][:start],f['velocities'][:start])
            for key in ('noise','alpha','initial_u'):assert np.array_equal(r[key],f[key])
        if number%10==0:print('Independent raw audit',number,'/65',flush=True)
    assert maximum_metric<1e-9 and maximum_motion<1e-10 and maximum_filter<1e-8 and maximum_cov<1e-9 and maximum_center<1e-8

    scopes=('pre20_mean','occlusion_mean','occlusion_last','occlusion_last3','post20_mean','fixed80_mean','full_mean')
    def ii(scope,a,w):
        if scope=='pre20_mean':return list(range(a-20,a))
        if scope=='occlusion_mean':return list(range(a,a+w))
        if scope=='occlusion_last':return [a+w-1]
        if scope=='occlusion_last3':return list(range(a+w-3,a+w))
        if scope=='post20_mean':return list(range(a+w,a+w+20))
        if scope=='fixed80_mean':return list(range(a,a+80))
        assert scope=='full_mean';return list(range(120))
    maximum_statistic=0.
    def close(a,b):
        nonlocal maximum_statistic
        error=abs(float(a)-float(b));assert math.isfinite(error) and error<1e-8,(a,b)
        maximum_statistic=max(maximum_statistic,error)
    ep={};endpoint_rows=read_csv(RESULTS/'analysis/endpoints_by_seed.csv')
    assert len(endpoint_rows)==1120
    for row in endpoint_rows:
        seed,start,width=int(row['seed']),int(row['start']),int(row['duration']);method=row['method'];scope=row['scope'];metric=row['metric']
        ts=raw[seed,method,0 if method=='fresh' else start,0 if method=='fresh' else width][metric]
        xx=[ts[k] for k in ii(scope,start,width)]
        close(row['mean'],st.mean(xx));close(row['sum_abs'],sum(xx));close(row['auc'],sum(xx)*dt)
        assert int(row['steps'])==len(xx)
        key=(seed,start,width,method,scope,metric);assert key not in ep;ep[key]=st.mean(xx)
    comparisons={'kf-mean - hold':('kf-mean','hold'),'kf-mean - drop':('kf-mean','drop'),
                 'kf-mean - fresh':('kf-mean','fresh'),'hold - fresh':('hold','fresh'),'drop - fresh':('drop','fresh')}
    pairs={};pair_rows=read_csv(RESULTS/'analysis/paired_by_seed.csv');assert len(pair_rows)==1400
    for row in pair_rows:
        seed,start,width=int(row['seed']),int(row['start']),int(row['duration']);scope,metric,contrast=row['scope'],row['metric'],row['contrast']
        a,b=comparisons[contrast];va,vb=ep[seed,start,width,a,scope,metric],ep[seed,start,width,b,scope,metric]
        for col,v in [('treatment',va),('reference',vb),('difference',va-vb),('relative_improvement_percent',100*(vb-va)/vb)]:close(row[col],v)
        pairs[seed,start,width,scope,metric,contrast]=va-vb
    summary_rows=read_csv(RESULTS/'analysis/paired_summary.csv');assert len(summary_rows)==280
    for row in summary_rows:
        a,w=int(row['start']),int(row['duration']);scope,metric,contrast=row['scope'],row['metric'],row['contrast']
        xx=[pairs[s,a,w,scope,metric,contrast] for s in seeds]
        for col,v in [('mean',st.mean(xx)),('sd',st.stdev(xx)),('minimum',min(xx)),('maximum',max(xx))]:close(row[col],v)
        assert int(row['n'])==5 and int(row['improved_seeds'])==sum(x<0 for x in xx) and int(row['worse_seeds'])==sum(x>0 for x in xx)
        ta,tb=comparisons[contrast]
        close(row['mean_relative_improvement_percent'],st.mean(100*(ep[s,a,w,tb,scope,metric]-ep[s,a,w,ta,scope,metric])/ep[s,a,w,tb,scope,metric] for s in seeds))
    performance_rows=read_csv(RESULTS/'analysis/performance_summary.csv');assert len(performance_rows)==224
    for row in performance_rows:
        xx=[ep[s,int(row['start']),int(row['duration']),row['method'],row['scope'],row['metric']] for s in seeds]
        for col,v in [('mean',st.mean(xx)),('sd',st.stdev(xx)),('minimum',min(xx)),('maximum',max(xx))]:close(row[col],v)
    diagnostics=read_csv(RESULTS/'analysis/estimator_diagnostics.csv');assert len(diagnostics)==60
    for row in diagnostics:
        a,w=int(row['start']),int(row['duration']);r=raw[int(row['seed']),row['method'],a,w]
        close(row['hidden_active_coverage'],r['hidden_coverage']);close(row['maximum_hidden_angle_variance'],r['maximum_hidden_angle_variance'])
        if row['method']!='drop':
            close(row['leader_bearing_mae_deg'],st.mean(r['leader_mae'][a:a+w]));close(row['leader_bearing_last_mae_deg'],r['leader_mae'][a+w-1])
        else:assert row['leader_bearing_mae_deg']==row['leader_bearing_last_mae_deg']==''
    primary=read_json(RESULTS/'analysis/PRIMARY_RESULTS.json')
    for row in primary['contrasts']:
        expected_row=next(r for r in summary_rows if r['start']=='40' and r['duration']=='40' and r['scope']=='occlusion_mean' and r['metric']==row['metric'] and r['contrast']==row['contrast'])
        close(row['mean'],expected_row['mean']);assert row['improved_seeds']==int(expected_row['improved_seeds'])
    write_json(RESULTS/'checks/INDEPENDENT_VERIFICATION.json',dict(status='PASS',unique_trajectories=65,
        tracking_aircraft_update_samples=samples,scalar_filter_edge_updates=65*77*120,raw_output_hashes_verified=130,
        maximum_tracking_metric_error=maximum_metric,maximum_displacement_error=maximum_motion,
        maximum_scalar_KF_error=maximum_filter,maximum_covariance_error=maximum_cov,maximum_input_center_error=maximum_center,
        maximum_neural_control_replay_error=maximum_control,maximum_statistic_error=maximum_statistic,
        paired_initial_histories_identical=True,hidden_observations_nan=True,all_primary_angle_samples_valid=True,
        previous_package_unchanged=source.exists(),endpoint_rows=1120,paired_rows=1400,summary_rows=280,performance_rows=224))
    print('INDEPENDENT VERIFICATION PASS',samples,'tracking samples;',maximum_metric,'maximum metric error',flush=True)


if __name__=='__main__':main()

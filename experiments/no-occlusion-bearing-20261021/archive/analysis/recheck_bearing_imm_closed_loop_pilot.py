"""Independent closed-loop trajectory, field and estimator reconstruction."""
import argparse
import csv
import hashlib
import json
import re
import sys
import time
from pathlib import Path

sys.path.insert(0,r'C:\Users\Administrator\AppData\Local\Temp\noisy-bearing-kf-pydeps')
import numpy as np
from scipy.io import loadmat
from scipy.stats import t as student_t
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import recheck_bearing_imm_replay_pilot as matrix_reference


def wrap(x):
    return np.arctan2(np.sin(x),np.cos(x))


def write_csv(path,rows):
    with path.open('w',encoding='utf-8-sig',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)


def ready_files(root,stream):
    """Only read a MAT after the producer closed it and published three CSV rows."""
    seen=set();deadline=time.monotonic()+7200
    while len(seen)<24:
        summary=root/('summary.csv' if (root/'summary.csv').exists() else 'progress-summary.csv')
        lookup={}
        if summary.exists():
            try:
                with summary.open(encoding='utf-8-sig',newline='') as f:rows=list(csv.DictReader(f))
                for r in rows:
                    if not r.get('elapsedSeconds'):continue
                    if r.get('method') not in ['raw','selected-kf','imm']:continue
                    lookup[int(r['seed']),r['distribution'],float(r['measurementStd']),r['method']]=r
            except (ValueError,KeyError):
                lookup={}
        found=False
        for path in sorted(root.glob('*-triple.mat')):
            if path.name in seen:continue
            match=re.fullmatch(r'(gaussian|mixture)-seed-(\d+)-std-(\d+)-triple',path.stem)
            assert match,path.name
            dist,seed,sd=match.group(1),int(match.group(2)),int(match.group(3))/1000
            if not all((seed,dist,sd,m) in lookup for m in ['raw','selected-kf','imm']):continue
            seen.add(path.name);found=True
            yield len(seen),path,lookup
        if len(seen)==24:return
        if not stream:raise AssertionError(f'Only {len(seen)} complete published triples')
        if time.monotonic()>deadline:raise TimeoutError('Stream audit reached its two-hour deadline; partial audit retained')
        if not found:time.sleep(10)


def check_one(r,label):
    cfg=r['config'];n=int(cfg['NAgents']);T=int(cfg['NSteps']);B=int(cfg['burnIn'])
    L=float(cfg['arenaSize']);dt=float(cfg['dt']);sigma=float(cfg['bearingMeasurementStd'])
    assert (n,T,B,int(cfg['neighborCount']))==(20,8000,6000,15)
    assert cfg['allocentric'] and not cfg['scheduledOcclusion'] and not cfg['turnExperiment']
    assert cfg['immLowQ']==1e-4 and cfg['immHighQ']==.1 and cfg['bearingProcessNoise']==.1
    assert np.array_equal(cfg['immTransition'],[[.98,.02],[.10,.90]])
    pilot=r['pilot'];ts=r['timeSeries'];x=r['trajectory']['x'];y=r['trajectory']['y']
    assert x.shape==y.shape==(20,8001)
    dx=np.diff(x,axis=1);dy=np.diff(y,axis=1)
    dx-=L*np.round(dx/L);dy-=L*np.round(dy/L)
    vx,vy=dx/dt,dy/dt;speed=np.hypot(vx,vy)
    ux=np.divide(vx,speed,out=np.zeros_like(vx),where=speed>1e-12)
    uy=np.divide(vy,speed,out=np.zeros_like(vy),where=speed>1e-12)
    go=np.hypot(ux.sum(axis=0),uy.sum(axis=0))/n
    mean_speed=speed.mean(axis=0)
    post_x=x[:,B+1:].T;post_y=y[:,B+1:].T
    px=post_x[:,:,None]-post_x[:,None,:];py=post_y[:,:,None]-post_y[:,None,:]
    px-=L*np.round(px/L);py-=L*np.round(py/L)
    distance=np.hypot(px,py);ii,jj=np.triu_indices(n,1)
    pd=distance[:,ii,jj].mean(axis=1)/L
    nearest=np.argsort(distance,axis=2)[:,:,:6]
    xx=np.take_along_axis(ux[:,B:].T[:,None,:],nearest,axis=2).sum(axis=2)/6
    yy=np.take_along_axis(uy[:,B:].T[:,None,:],nearest,axis=2).sum(axis=2)/6
    lo=np.hypot(xx,yy).mean(axis=1)
    errors={'GO':float(np.max(np.abs(go-ts['globalOrder']))),
            'LO':float(np.max(np.abs(lo-ts['localOrder'][B:]))),
            'distance':float(np.max(np.abs(pd-ts['meanPairDistanceNorm'][B:]))),
            'speed':float(np.max(np.abs(mean_speed-ts['meanSpeed'])))}
    # Inputs at model step t use the position at the start of that step.
    pre_x=x[:,:T].T;pre_y=y[:,:T].T
    rx=pre_x[:,None,:]-pre_x[:,:,None];ry=pre_y[:,None,:]-pre_y[:,:,None]
    rx-=L*np.round(rx/L);ry-=L*np.round(ry/L)
    edges=np.asarray(pilot['linearEdges']).astype(int)-1;rec=edges%n;nei=edges//n
    truth=np.arctan2(ry[:,rec,nei],rx[:,rec,nei])%(2*np.pi)
    ids=np.asarray(pilot['selectedHistory']).astype(int).transpose(2,0,1)-1
    chosen=np.zeros((T,n,n),dtype=bool)
    chosen[np.arange(T)[:,None,None],np.arange(n)[None,:,None],ids]=True
    fresh=chosen[:,rec,nei]
    obs=np.asarray(pilot['observedBearingHistory']);used=np.asarray(pilot['usedBearingEvaluation'])
    assert obs.shape==(T,380) and used.shape==(2000,380)
    assert np.array_equal(fresh,np.isfinite(obs))
    assert np.array_equal(fresh[B:],np.isfinite(used))
    assert np.all(fresh.sum(axis=1)==300)
    assert pilot['evalNoise']['count']==600000
    if label=='raw':assert np.array_equal(obs[B:],used,equal_nan=True)
    obs_error=wrap(obs[B:]-truth[B:]);use_error=wrap(used-truth[B:]);mask=fresh[B:]
    obs_step=np.sqrt(np.nansum(obs_error**2,axis=1)/300)
    use_step=np.sqrt(np.nansum(use_error**2,axis=1)/300)
    errors['observedAngleRMSE']=float(np.max(np.abs(obs_step-ts['observedBearingRMSE'][B:])))
    errors['usedAngleRMSE']=float(np.max(np.abs(use_step-ts['usedBearingRMSE'][B:])))
    alpha=np.asarray(r['finalPreferredDirections']);samples=pilot['samples']
    for sample in samples:assert np.array_equal(sample['alpha'],alpha)
    assert alpha.shape==(20,100)
    column_map=np.full((n,n),-1,dtype=int);column_map[rec,nei]=np.arange(380)
    total_sq=np.zeros(T-B);h=float(cfg['totalSocialAttraction'])/n;kernel=float(cfg['receptiveFieldWidth'])
    # Receiver-wise tensors keep memory bounded while checking ALL 2000 steps.
    for receiver in range(n):
        columns=column_map[receiver,ids[B:,receiver,:]]
        used_angles=np.take_along_axis(used,columns,axis=1)
        true_angles=np.take_along_axis(truth[B:],columns,axis=1)
        fields=[]
        for angle in [used_angles,true_angles]:
            delta=np.abs(alpha[receiver,None,:,None]-angle[:,None,:])
            delta=np.minimum(delta,2*np.pi-delta)
            fields.append(h*np.exp(-.5*delta**2/kernel**2).sum(axis=2))
        total_sq+=np.sum((fields[0]-fields[1])**2,axis=1)
    input_step=np.sqrt(total_sq/(100*n))
    errors['fullEvaluationInputRMS']=float(np.max(np.abs(input_step-ts['inputTruthRMS'][B:])))
    sample_error=0.
    for sample in samples:
        step=int(sample['step'])-1;ix=np.asarray(sample['linearIndex']).astype(int)-1
        assert np.array_equal(ix,np.flatnonzero(chosen[step].ravel(order='F')))
        re=ix%n;ne=ix//n
        assert np.max(np.abs(wrap(sample['trueBearing']-np.arctan2(ry[step,re,ne],rx[step,re,ne]))))<1e-12
        assert np.max(np.abs(wrap(sample['observedBearing']-sample['trueBearing']-sigma*sample['standardNoise'])))<1e-12
        fields=[]
        for angle in [sample['usedBearing'],sample['trueBearing']]:
            centers=np.full((n,n),np.nan);centers[re,ne]=angle
            field=np.zeros((100,n))
            order=np.asarray(sample['orderedIds']).astype(int)-1
            for receiver in range(n):
                delta=np.abs(alpha[receiver,:,None]-centers[receiver,order[receiver]][None,:])
                delta=np.minimum(delta,2*np.pi-delta)
                field[:,receiver]=h*np.exp(-.5*delta**2/kernel**2).sum(axis=1)
            fields.append(field)
        sample_error=max(sample_error,float(np.max(np.abs(fields[0]-sample['input']))),float(np.max(np.abs(fields[1]-sample['truthInput']))))
    errors['archivedFieldSample']=sample_error
    assert max(errors.values())<1e-9,errors
    row={'seed':int(cfg['seed']),'distribution':cfg['noiseDistribution'],'measurementStd':sigma,
         'method':cfg['noisyBearingMode'],'meanGO':float(go[B:].mean()),'meanLO':float(lo.mean()),
         'meanPairDistanceNorm':float(pd.mean()),'meanSpeed':float(mean_speed[B:].mean()),
         'observedBearingRMSE':float(np.sqrt(np.mean(obs_error[mask]**2))),
         'usedBearingRMSE':float(np.sqrt(np.mean(use_error[mask]**2))),
         'inputTruthRMS':float(np.sqrt(np.mean(input_step**2))),'evaluationPointCount':int(mask.sum())}
    filter_check=None
    if int(cfg['seed'])==6201 and label!='raw':
        state=None;max_angle=0.;max_probability=0.
        for step in range(T):
            if label=='kf':
                state,_,_=matrix_reference.kf_matrix(state,obs[step],fresh[step],dt,.1,sigma**2)
                estimate=state['x'][:,0]
            else:
                state,estimate,covariance=matrix_reference.imm_matrix(state,obs[step],fresh[step],dt,.1,sigma**2,np.array([[.98,.02],[.10,.90]]))
                estimate=estimate[:,0]
                if step>=B:
                    max_probability=max(max_probability,float(np.max(np.abs(state['mu'][:,1]-pilot['highProbabilityEvaluation'][step-B]))))
            if step>=B:
                f=fresh[step]
                max_angle=max(max_angle,float(np.max(np.abs(wrap(estimate[f]-used[step-B,f])))))
        assert max_angle<1e-8 and max_probability<1e-8
        filter_check={'seed':6201,'distribution':cfg['noiseDistribution'],'measurementStd':sigma,
                      'method':cfg['noisyBearingMode'],'steps':8000,'maxAngleDifference':max_angle,
                      'maxModelProbabilityDifference':max_probability}
    return row,errors,go,filter_check


def comparison(rows,distribution,sd,metric,baseline):
    group=[r for r in rows if r['distribution']==distribution and r['measurementStd']==sd]
    key={(r['seed'],r['method']):r for r in group};seeds=list(range(6201,6207))
    a=np.array([key[s,baseline][metric] for s in seeds]);b=np.array([key[s,'imm'][metric] for s in seeds])
    d=b-a;se=float(d.std(ddof=1)/np.sqrt(6));mean=float(d.mean());margin=float(student_t.ppf(.975,5))*se
    p=float(2*student_t.sf(abs(mean/se),5)) if se else (0. if mean else 1.)
    return {'distribution':distribution,'measurementStd':sd,'metric':metric,'baseline':baseline,
            'nSeeds':6,'baselineMean':float(a.mean()),'immMean':float(b.mean()),'difference':mean,
            'CI95Low':mean-margin,'CI95High':mean+margin,'pUnadjusted':p,
            'positiveSeeds':int(np.sum(d>0)),'negativeSeeds':int(np.sum(d<0)),
            'differences':d.tolist()}


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--run',type=Path,required=True)
    parser.add_argument('--stream',action='store_true')
    args=parser.parse_args();root=args.run.resolve();out=root/'independent-analysis'
    assert not out.exists(),'Existing analysis is protected';out.mkdir()
    meta=json.loads((root/'metadata.json').read_text(encoding='utf-8'))
    if not args.stream:
        complete=json.loads((root/'completion.json').read_text(encoding='utf-8'))
        assert complete['status']=='COMPLETED' and complete['runs']==72 and complete['triples']==24
    recomputed=[];checks=[];pairing=[];hashes={};series={};max_errors={};filter_checks=[]
    for fi,path,lookup in ready_files(root,args.stream):
        hashes[path.name]=hashlib.sha256(path.read_bytes()).hexdigest()
        triple=loadmat(path,simplify_cells=True)['triple'];a=triple['raw']
        for label in ['kf','imm']:
            b=triple[label]
            for field in ['initialPosition','initialHeadings','initialNeuralStateUsed']:
                assert np.array_equal(a[field],b[field])
            for field in ['selectedHistory','finalNoiseStreamState','finalMixtureStreamState']:
                assert np.array_equal(a['pilot'][field],b['pilot'][field])
            assert a['pilot']['evalNoise']==b['pilot']['evalNoise']
            for sa,sb in zip(a['pilot']['samples'],b['pilot']['samples']):
                assert np.array_equal(sa['standardNoise'],sb['standardNoise'])
                assert np.array_equal(sa['standardNormal'],sb['standardNormal'])
                assert np.array_equal(sa['outlierMask'],sb['outlierMask'])
            assert np.array_equal(np.isnan(a['pilot']['observedBearingHistory']),np.isnan(b['pilot']['observedBearingHistory']))
        pairing.append({'file':path.name,'initialization':True,'allStepSelection':True,
                        'latentNoiseAndLabels':True,'sameFreshMask':True})
        for label in ['raw','kf','imm']:
            r=triple[label];row,errors,go,filter_check=check_one(r,label)
            old=lookup[row['seed'],row['distribution'],row['measurementStd'],row['method']]
            for name,value in row.items():
                if name=='evaluationPointCount':continue
                if name in ['distribution','method']:assert old[name]==value
                else:assert abs(float(old[name])-value)<1e-9,(path.name,label,name,old[name],value)
            recomputed.append(row);series[row['seed'],row['distribution'],row['measurementStd'],row['method']]=go
            for name,value in errors.items():max_errors[name]=max(max_errors.get(name,0.),value)
            checks.append({'file':path.name,'method':row['method'],'errors':errors,'evaluationPoints':row['evaluationPointCount']})
            if filter_check is not None:filter_checks.append(filter_check)
            print(f'RECHECK {fi}/24 {path.name} {label}: GO={row["meanGO"]:.6f} angle={row["usedBearingRMSE"]:.6f} maxMetricError={max(errors.values()):.3g}',flush=True)
        del triple,a,b,r
    for _ in range(30):
        if (root/'completion.json').exists():break
        time.sleep(2)
    complete=json.loads((root/'completion.json').read_text(encoding='utf-8'))
    assert complete['status']=='COMPLETED' and complete['runs']==72 and complete['triples']==24
    with (root/'summary.csv').open(encoding='utf-8-sig',newline='') as f:archived=list(csv.DictReader(f))
    assert len(archived)==len({(r['seed'],r['distribution'],r['measurementStd'],r['method']) for r in archived})==72
    assert len(recomputed)==72 and len(filter_checks)==8
    assert {r['seed'] for r in recomputed}==set(range(6201,6207))
    write_csv(out/'independently-recomputed.csv',recomputed)
    conditions=[(d,s) for d in ['gaussian','mixture'] for s in [.01,.1]]
    primary=[comparison(recomputed,d,s,'meanGO','selected-kf') for d,s in conditions]
    order=np.argsort([r['pUnadjusted'] for r in primary]);previous=0.
    for rank,index in enumerate(order):
        adjusted=min(1.,max(previous,(4-rank)*primary[index]['pUnadjusted']))
        primary[index]['pHolmFour']=adjusted;previous=adjusted
    secondary=[comparison(recomputed,d,s,m,b) for d,s in conditions
               for b in ['raw','selected-kf']
               for m in ['meanGO','meanLO','meanPairDistanceNorm','meanSpeed','usedBearingRMSE','inputTruthRMS']
               if not (m=='meanGO' and b=='selected-kf')]
    statistics={'primary':primary,'secondaryUnadjustedExploratory':secondary}
    (out/'paired-statistics.json').write_text(json.dumps(statistics,indent=2),encoding='utf-8')
    write_csv(out/'primary-paired-statistics.csv',[{k:v for k,v in r.items() if k!='differences'} for r in primary])
    aggregate=[]
    for d,s in conditions:
        for method in ['raw','selected-kf','imm']:
            group=[r for r in recomputed if r['distribution']==d and r['measurementStd']==s and r['method']==method]
            aggregate.append({'distribution':d,'measurementStd':s,'method':method,'nSeeds':6,
                              **{m:float(np.mean([r[m] for r in group])) for m in ['meanGO','meanLO','meanPairDistanceNorm','meanSpeed','usedBearingRMSE','inputTruthRMS']}})
    write_csv(out/'condition-means.csv',aggregate)
    (out/'condition-means.json').write_text(json.dumps(aggregate,indent=2),encoding='utf-8')
    colors=['#5c7386','#a16a47','#4269a6'];methods=['raw','selected-kf','imm'];labels=['Raw','Selected KF','IMM']
    for metric,filename,title in [('meanGO','paired-go.png','Actual closed-loop global order'),
                                   ('meanLO','paired-lo.png','Actual closed-loop local order'),
                                   ('meanPairDistanceNorm','paired-distance.png','Actual periodic mean pair distance / L'),
                                   ('usedBearingRMSE','paired-angle-rmse.png','Bearing RMSE in each strategy own closed-loop geometry')]:
        fig,axes=plt.subplots(2,2,figsize=(10,7),layout='constrained')
        for ax,(d,s) in zip(axes.flat,conditions):
            key={(r['seed'],r['method']):r for r in recomputed if r['distribution']==d and r['measurementStd']==s}
            values=np.array([[key[seed,m][metric] for m in methods] for seed in range(6201,6207)])
            for value in values:ax.plot(range(3),value,'o-',color='#b3bac3',alpha=.6,lw=.8,ms=3)
            for x,mean,color in zip(range(3),values.mean(axis=0),colors):ax.plot(x,mean,'o',color=color,ms=7,zorder=4)
            ax.set_xticks(range(3),labels);ax.set_title(f'{d}, sigma={s:.2f}');ax.grid(axis='y',alpha=.2)
            ax.set_ylabel(metric)
            if metric in ['meanGO','meanLO']:ax.set_ylim(0,1.03)
            if metric=='usedBearingRMSE':ax.set_yscale('log')
        fig.suptitle(title+'\nSix new paired seeds; fixed upstream parameters. Lines: seeds, colored points: means.',fontsize=12)
        fig.savefig(out/filename,dpi=170);plt.close(fig)
    fig,axes=plt.subplots(2,2,figsize=(11,7),layout='constrained')
    for ax,(d,s) in zip(axes.flat,conditions):
        for method,label,color in zip(methods,labels,colors):
            a=np.stack([series[seed,d,s,method] for seed in range(6201,6207)]).reshape(6,80,100).mean(axis=2)
            ax.plot(np.arange(1,81)*100,a.mean(axis=0),label=label,color=color)
        ax.axvspan(6000,8000,color='gray',alpha=.12);ax.set_ylim(0,1);ax.set_title(f'{d}, sigma={s:.2f}')
        ax.set_xlabel('Model step (100-step block means)');ax.set_ylabel('Actual GO');ax.legend(fontsize=9)
    fig.suptitle('Closed-loop GO time courses; evaluation window fixed at steps 6001-8000',fontsize=12)
    fig.savefig(out/'mean-go-timeseries.png',dpi=170);plt.close(fig)
    audit={'status':'PASS_INDEPENDENT_CLOSED_LOOP_RECONSTRUCTION','runs':72,'uniqueKeys':72,'triples':24,
           'seeds':list(range(6201,6207)),'independentUnit':'seed','evaluationPointsPerRun':600000,
           'maxErrors':max_errors,'pairing':pairing,'checks':checks,'fileHashes':hashes,
           'filterSecondImplementationChecks':filter_checks,'filterReplayScope':'seed6201, both distributions, both noise levels, KF/IMM; 8 complete8000-step replays',
           'inputFieldReplayScope':'all2000 evaluation steps of all72 runs plus all archived samples',
           'evidence':'new-seed small-sample closed-loop test; fixed upstream parameters; no tuning against these seeds',
           'limitations':['6seeds percondition','fixedN20/k15/neuralparameters/periodicgeometry',
                          'closed-loop bearings differ after geometry diverges','secondarycomparisons exploratory',
                          'means are not sustained stability or safety guarantees']}
    (out/'audit.json').write_text(json.dumps(audit,indent=2),encoding='utf-8')
    print(json.dumps({'status':audit['status'],'maxErrors':max_errors,'independentFullFilterReplays':8,'primary':primary,'conditionMeans':aggregate},indent=2),flush=True)


if __name__=='__main__':main()

"""All fixed cells, both tracking metrics, paired by independent seed."""
import bootstrap
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import design as d
from common import RESULTS,read_csv,read_json,write_csv,write_json

SCOPES=('pre20_mean','occlusion_mean','occlusion_last','occlusion_last3','post20_mean','fixed80_mean','full_mean')
METRICS=('angle_deg','distance')
CONTRASTS=(('kf-mean','hold'),('kf-mean','drop'),('kf-mean','fresh'),('hold','fresh'),('drop','fresh'))
COLORS={'fresh':'#222222','drop':'#D55E00','hold':'#009E73','kf-mean':'#0072B2'}
LABELS={'fresh':'Fresh','drop':'Drop','hold':'Hold','kf-mean':'KF'}


def indices(scope,start,width):
    end=start+width
    a,b={'pre20_mean':(start-20,start),'occlusion_mean':(start,end),
         'occlusion_last':(end-1,end),'occlusion_last3':(end-3,end),
         'post20_mean':(end,end+20),'fixed80_mean':(start,start+80),'full_mean':(0,d.STEPS)}[scope]
    assert 0<=a<b<=d.STEPS
    return np.arange(a,b)


def main():
    out=RESULTS/'analysis';out.mkdir(exist_ok=True)
    figs=RESULTS/'figures';figs.mkdir(exist_ok=True)
    manifest=read_json(RESULTS/'run_manifest.json');assert len(manifest)==65
    data={}
    for r in manifest:
        rows=read_csv(RESULTS/'timeseries'/(r['run_id']+'.csv'))
        data[r['seed'],r['method'],r['start'],r['duration']]={name:np.array([float(row[name]) for row in rows])
            for name in ('angle_deg','distance','angle_valid_fraction','stationary_fraction','follower_GO','leader_target_distance')}
    def get(seed,method,start,width):return data[seed,method,0 if method=='fresh' else start,0 if method=='fresh' else width]
    endpoint_rows=[];ep={}
    for seed in d.SEEDS:
        for start in d.STARTS:
            for width in d.DURATIONS:
                for method in d.METHODS:
                    for scope in SCOPES:
                        ii=indices(scope,start,width)
                        for metric in METRICS:
                            values=get(seed,method,start,width)[metric][ii]
                            assert np.isfinite(values).all()
                            value=float(values.mean());ep[seed,start,width,method,scope,metric]=value
                            endpoint_rows.append(dict(seed=seed,start=start,duration=width,method=method,scope=scope,metric=metric,
                                steps=len(ii),mean=value,sum_abs=float(values.sum()),auc=float(values.sum()*d.DT)))
    write_csv(out/'endpoints_by_seed.csv',endpoint_rows)
    paired=[];pair_index={};summary=[];performance=[]
    for start in d.STARTS:
        for width in d.DURATIONS:
            for scope in SCOPES:
                for metric in METRICS:
                    for method in d.METHODS:
                        x=np.array([ep[s,start,width,method,scope,metric] for s in d.SEEDS])
                        performance.append(dict(start=start,duration=width,method=method,scope=scope,metric=metric,n=len(x),
                                                mean=float(x.mean()),sd=float(x.std(ddof=1)),minimum=float(x.min()),maximum=float(x.max())))
                    for a,b in CONTRASTS:
                        contrast=a+' - '+b
                        differences=[];relative=[]
                        for seed in d.SEEDS:
                            va,vb=ep[seed,start,width,a,scope,metric],ep[seed,start,width,b,scope,metric]
                            delta=va-vb;improvement=100*(vb-va)/vb
                            row=dict(seed=seed,start=start,duration=width,scope=scope,metric=metric,contrast=contrast,
                                     treatment=va,reference=vb,difference=delta,relative_improvement_percent=improvement)
                            paired.append(row);pair_index[seed,start,width,scope,metric,contrast]=delta
                            differences.append(delta);relative.append(improvement)
                        x=np.array(differences)
                        summary.append(dict(start=start,duration=width,scope=scope,metric=metric,contrast=contrast,n=len(x),
                            mean=float(x.mean()),sd=float(x.std(ddof=1)),minimum=float(x.min()),maximum=float(x.max()),
                            improved_seeds=int((x<0).sum()),worse_seeds=int((x>0).sum()),ties=int((x==0).sum()),
                            mean_relative_improvement_percent=float(np.mean(relative))))
    write_csv(out/'paired_by_seed.csv',paired);write_csv(out/'paired_summary.csv',summary)
    write_csv(out/'performance_summary.csv',performance)
    diagnostics=[]
    for r in manifest:
        if r['method']=='fresh':continue
        with np.load(RESULTS/'trajectories'/(r['run_id']+'.npz')) as z:
            a,w=r['start'],r['duration'];ii=np.arange(a,a+w)
            xy=z['positions'][ii]
            bearing=np.arctan2(xy[:,None,:,1]-xy[:,:,None,1],xy[:,None,:,0]-xy[:,:,None,0])
            delta=z['centers'][ii]-bearing
            error=np.degrees(np.abs(np.arctan2(np.sin(delta),np.cos(delta))))
            active=z['active'][ii]
            leader_error=error[:,1:,0]
            if r['method']=='drop':
                # Missing inputs have no estimate in Drop. Do not score its zero placeholder as a prediction.
                diagnostics.append(dict(seed=r['seed'],start=a,duration=w,method=r['method'],
                    interpretation='no prediction; missing inputs removed',leader_bearing_mae_deg='',leader_bearing_last_mae_deg='',
                    hidden_active_coverage=float(active[:,z['hidden']].mean()),maximum_hidden_angle_variance=float(z['filter_covariance'][ii][:,z['hidden'],0,0].max())))
            else:
                assert active[:,z['hidden']].all()
                diagnostics.append(dict(seed=r['seed'],start=a,duration=w,method=r['method'],
                    interpretation='estimation error on this method own closed-loop trajectory, not common geometry',
                    leader_bearing_mae_deg=float(leader_error.mean()),leader_bearing_last_mae_deg=float(leader_error[-1].mean()),
                    hidden_active_coverage=float(active[:,z['hidden']].mean()),maximum_hidden_angle_variance=float(z['filter_covariance'][ii][:,z['hidden'],0,0].max())))
    write_csv(out/'estimator_diagnostics.csv',diagnostics)
    primary=[]
    for metric in METRICS:
        for contrast in ('kf-mean - hold','kf-mean - drop'):
            primary.append(next(r for r in summary if r['start']==40 and r['duration']==40 and r['scope']=='occlusion_mean'
                                and r['metric']==metric and r['contrast']==contrast))
    write_json(out/'PRIMARY_RESULTS.json',dict(primary_start=40,primary_duration=40,scope='occlusion_mean',n=5,
        all_four_mean_contrasts_negative=all(r['mean']<0 for r in primary),contrasts=primary,
        result_scope='exploratory early single-window tracking; not a causal test of initial randomness'))

    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    cells=[(s,w) for s in d.STARTS for w in d.DURATIONS]
    def save(fig,name):
        fig.tight_layout(rect=(0,0,1,.94));fig.savefig(figs/name,dpi=170);plt.close(fig)
    def perf(start,width,scope,metric,method):
        return next(r for r in performance if r['start']==start and r['duration']==width and r['scope']==scope and r['metric']==metric and r['method']==method)
    fig,axes=plt.subplots(2,4,figsize=(15,7.8))
    for col,(start,width) in enumerate(cells):
        for row,metric in enumerate(METRICS):
            ax=axes[row,col]
            rows=[perf(start,width,'occlusion_mean',metric,m) for m in d.METHODS]
            ax.bar(np.arange(4),[r['mean'] for r in rows],yerr=[r['sd'] for r in rows],capsize=3,
                   color=[COLORS[m] for m in d.METHODS],alpha=.85)
            ax.set(xticks=np.arange(4),xticklabels=[LABELS[m] for m in d.METHODS],title=f'Start={start}, width={width}',
                   ylabel='Angle error (deg)' if row==0 else 'Target distance (model units)')
            ax.set_ylim(bottom=0);ax.grid(axis='y',alpha=.2)
    fig.suptitle('Single early occlusion | mean occluded tracking error | n=5 seeds | error bars: sample SD, not CI')
    save(fig,'01_four_method_tracking.png')
    fig,axes=plt.subplots(2,2,figsize=(12,8))
    xx=np.arange(4);offsets=np.linspace(-.16,.16,len(d.SEEDS))
    for col,contrast in enumerate(('kf-mean - hold','kf-mean - drop')):
        for row,metric in enumerate(METRICS):
            ax=axes[row,col]
            for offset,seed in zip(offsets,d.SEEDS):
                y=[pair_index[seed,s,w,'occlusion_mean',metric,contrast] for s,w in cells]
                ax.scatter(xx+offset,y,s=25,alpha=.7,label=str(seed))
            rows=[next(r for r in summary if r['start']==s and r['duration']==w and r['scope']=='occlusion_mean' and r['metric']==metric and r['contrast']==contrast) for s,w in cells]
            ax.errorbar(xx,[r['mean'] for r in rows],yerr=[r['sd'] for r in rows],fmt='kD',capsize=4,label='Mean +/- SD')
            ax.axhline(0,color='gray',linestyle='--')
            ax.set(xticks=xx,xticklabels=[f'S{s}/D{w}' for s,w in cells],title=contrast.replace('kf-mean','KF'),
                   ylabel='Angle difference (deg)' if row==0 else 'Distance difference (model units)')
            ax.grid(alpha=.2)
    axes[0,0].legend(fontsize=8,ncol=2)
    fig.suptitle('Paired occlusion-mean differences | negative = KF better | each dot is one independent seed')
    save(fig,'02_KF_paired_comparisons.png')
    fig,axes=plt.subplots(2,1,figsize=(12,8),sharex=True)
    t=(np.arange(d.STEPS)+1)*d.DT
    for ax,metric in zip(axes,METRICS):
        for method in d.METHODS:
            arr=np.array([get(s,method,40,40)[metric] for s in d.SEEDS])
            mean,sd=arr.mean(0),arr.std(0,ddof=1)
            ax.plot(t,mean,color=COLORS[method],label=LABELS[method],linewidth=1.8)
            ax.fill_between(t,mean-sd,mean+sd,color=COLORS[method],alpha=.08)
        ax.axvspan(40*d.DT,80*d.DT,color='gray',alpha=.17)
        ax.set_ylabel('Angle error (deg)' if metric=='angle_deg' else 'Target distance (model units)')
        ax.grid(alpha=.2)
    axes[0].legend(ncol=4);axes[-1].set_xlabel('Model time (not calibrated seconds)')
    fig.suptitle('Primary condition S40/D40 | mean +/- sample SD | gray: single occlusion interval')
    save(fig,'03_primary_error_timeseries.png')
    runs={}
    for method in d.METHODS:
        rid=d.run_id(971001,method,0 if method=='fresh' else 40,0 if method=='fresh' else 40)
        with np.load(RESULTS/'trajectories'/(rid+'.npz')) as z:runs[method]={k:z[k] for k in ('positions','goals')}
    points=np.concatenate([r['positions'].reshape(-1,2) for r in runs.values()]+[runs['fresh']['goals']])
    low,high=points.min(0)-3,points.max(0)+3
    fig,axes=plt.subplots(2,2,figsize=(12,8))
    for ax,method in zip(axes.flat,d.METHODS):
        r=runs[method]
        for agent in d.FOLLOWERS:
            ax.plot(r['positions'][:,agent,0],r['positions'][:,agent,1],color=COLORS[method],alpha=.5,linewidth=.9,
                    label='Followers' if agent==1 else None)
        ax.plot(r['positions'][:,0,0],r['positions'][:,0,1],color='#E69F00',linewidth=2,label='Leader')
        ax.plot(r['goals'][:,0],r['goals'][:,1],'k--',linewidth=1.5,label='Moving target')
        ax.scatter(r['positions'][0,1:,0],r['positions'][0,1:,1],s=15,color=COLORS[method])
        ax.scatter(r['positions'][-1,1:,0],r['positions'][-1,1:,1],marker='x',s=25,color=COLORS[method])
        ax.set(title=LABELS[method],xlabel='x (model units)',ylabel='y (model units)',xlim=(low[0],high[0]),ylim=(low[1],high[1]))
        ax.set_aspect('equal',adjustable='box');ax.grid(alpha=.2)
    axes[0,0].legend(fontsize=8)
    fig.suptitle('Prespecified illustration: seed 971001, S40/D40 | circles=start; crosses=end | same axes')
    save(fig,'04_primary_trajectories.png')
    for r in primary:print('PRIMARY',r,flush=True)
    print('All four conditions and seven scopes retained.',flush=True)


if __name__=='__main__':main()

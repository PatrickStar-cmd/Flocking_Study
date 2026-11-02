"""Independent position audit and complete G0 supports/calibration tables."""
import json,sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.io import loadmat
from audit_independent import position_labels


def main(old,out):
    assert json.loads((out/'audit.json').read_text(encoding='utf-8'))['status']=='PASS'
    assert not (out/'independent_labels.json').exists()
    design=json.loads((old/'design.json').read_text(encoding='utf-8'))
    cache={};labels={};rates={};pool=[];count=pd.read_csv(out/'geometry_counts.csv')
    for case in design['cases']:
        source=case['source']
        if source not in cache:
            x=loadmat(old/source,simplify_cells=True)
            pos=x['positions'];i=x['observer_id'].astype(int)-1;j=x['target_id'].astype(int)-1
            raw=(pos[:,j,:]-pos[:,i,:]).transpose(2,1,0)
            image=np.copysign(np.floor(np.abs(raw/1000)+.5),raw)
            delta=raw-1000*image;distance=np.linalg.norm(delta,axis=-1)
            switch=np.zeros((2200,140),bool);switch[1:]=(np.abs(delta[1:]-delta[:-1])>500).any(axis=-1)
            truth=x['true_bearing'];change=np.full_like(truth,np.nan)
            change[1:]=np.arctan2(np.sin(truth[1:]-truth[:-1]),np.cos(truth[1:]-truth[:-1]))
            near=distance<5;alias=np.abs(change)>=.9*np.pi
            rate=change/.3
            rate[1:][switch[1:]|near[1:]|near[:-1]]=np.nan
            cache[source]={'source':x,'delta':delta,'image':image,'distance':distance,'switch':switch,'near':near,
                'half':np.min(500-np.abs(delta),axis=-1)<10,'alias':alias,'rate':rate}
        g=cache[source];start=int(case['t0'])-1
        for D in [20,120]:
            lb=position_labels(g['source'],np.array([case['t0']]),D)[0]
            r=np.sqrt(np.nanmean(g['rate'][start:start+D]**2,axis=0))
            labels[(case['case_id'],D)]=lb;rates[(case['case_id'],D)]=r
            pool.extend(r[lb=='G-clean'].tolist())
            for group in ['G-clean','G-switch','G-close','G-ambiguous']:
                row=count[(count.case_id==case['case_id'])&(count.D==D)&(count.geometry_class==group)]
                assert len(row)==1 and row.iloc[0].count_edges==(lb==group).sum()
    checked=0
    for chunk in pd.read_csv(out/'geometry_events.csv.gz',chunksize=150000):
        for source,frame in chunk.groupby('source'):
            g=cache[source];step=frame.step_global.to_numpy(int);edge=frame.edge.to_numpy(int)
            assert np.max(np.abs(g['delta'][step,edge,0]-frame.dx_min))<=.000051
            assert np.max(np.abs(g['delta'][step,edge,1]-frame.dy_min))<=.000051
            assert np.array_equal(g['image'][step,edge,0],frame.image_index_x)
            assert np.array_equal(g['image'][step,edge,1],frame.image_index_y)
            for key,column in [('switch','image_switch'),('near','near_range'),('half','near_half_box'),('alias','near_aliasing')]:
                assert np.array_equal(g[key][step,edge],frame[column].to_numpy(bool))
            checked+=len(frame)
    cuts=np.quantile(pool,[1/3,2/3]);rows=[];coverage=[];support=[]
    for case in design['cases']:
        for D in [20,120]:
            lb=labels[(case['case_id'],D)];r=rates[(case['case_id'],D)];clean=lb=='G-clean'
            groups={'all':np.ones(140,bool),**{g:lb==g for g in ['G-clean','G-switch','G-close','G-ambiguous']},
                'clean-slow':clean&(r<=cuts[0]),'clean-mid':clean&(r>cuts[0])&(r<=cuts[1]),'clean-fast':clean&(r>cuts[1])}
            with np.load(old/'raw/replay'/f"case_{case['case_id']}_D_{D}.npz") as raw:
                for group,mask in groups.items():
                    support.append({'seed':int(case['seed']),'case_id':case['case_id'],'D':D,'geometry':group,'n_edges':int(mask.sum()),'edge_times':int(mask.sum()*D)})
                    for method,key in [('Hold','Hold'),('KF-CV','KF_CV'),('KF-GM','KF_GM')]:
                        e=raw['theta_'+key][200:200+D][:,mask]-raw['truth'][200:200+D][:,mask]
                        e=np.arctan2(np.sin(e),np.cos(e));a=np.abs(e)
                        row={'seed':int(case['seed']),'case_id':case['case_id'],'D':D,'geometry':group,'method':method,'n_edges':int(mask.sum()),'edge_times':e.size}
                        if e.size:row.update(rmse=float(np.sqrt(np.mean(e*e))),mae=float(a.mean()),p90=float(np.quantile(a,.90)),p95=float(np.quantile(a,.95)),p99=float(np.quantile(a,.99)),max_abs=float(a.max()))
                        else:row.update(rmse=np.nan,mae=np.nan,p90=np.nan,p95=np.nan,p99=np.nan,max_abs=np.nan)
                        rows.append(row)
                        if method=='Hold':continue
                        for lo,hi in [(1,20),(21,60),(61,100),(101,120)]:
                            if lo>D:continue
                            ee=e[lo-1:min(hi,D)];p=raw['p11_'+key][200+lo-1:200+min(hi,D)][:,mask]
                            c={**row,'age_band':f'{lo}-{hi}','n_samples':ee.size}
                            c.update(mean_squared_error=float(np.mean(ee**2)) if ee.size else np.nan,
                                mean_P=float(p.mean()) if ee.size else np.nan,
                                coverage=float((np.abs(ee)<=1.96*np.sqrt(p)).mean()) if ee.size else np.nan)
                            coverage.append(c)
    pd.DataFrame(rows).to_csv(out/'error_by_geometry_complete.csv',index=False)
    pd.DataFrame(coverage).to_csv(out/'coverage_by_geometry_age.csv',index=False)
    f=pd.DataFrame(support);f['has_data']=f.n_edges>0
    f.groupby(['geometry','D'],as_index=False).agg(n_seeds=('seed',lambda x:x[f.loc[x.index,'has_data']].nunique()),
        edge_windows=('n_edges','sum'),edge_times=('edge_times','sum'),nonempty_cases=('has_data','sum')).to_csv(out/'geometry_counts_complete.csv',index=False)
    report={'status':'PASS','position_event_rows_independently_verified':checked,'all_case_edge_labels_verified':True,
        'coordinate_output_rounding_tolerance':.000051,'pilot_only_rate_cutpoints':cuts.tolist(),
        'missing_support_explicit_NaN':True,'development_thresholds_not_selected_from_old_seeds':True}
    (out/'independent_labels.json').write_text(json.dumps(report,indent=2),encoding='utf-8');print(json.dumps(report))


if __name__=='__main__':main(Path(sys.argv[1]),Path(sys.argv[2]))

"""Complete G0 NIS supports and quantify descriptive innovation-tail mass."""
import json,sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.io import loadmat
from geometry import geometry,fresh_labels,classify_window


def run(old,out):
    target=out/'nis_tail_complete.csv';assert not target.exists()
    design=json.loads((old/'design.json').read_text(encoding='utf-8'))
    rows=[];bins=[]
    for case in design['cases']:
        g=geometry(loadmat(old/case['source'],simplify_cells=True));fl=fresh_labels(g)
        for D in [20,120]:
            labels,_=classify_window(g,case['t0'],D)
            with np.load(old/'raw/replay'/f"case_{case['case_id']}_D_{D}.npz") as raw:
                indices=raw['window_steps'].astype(int)-1;point=fl[indices]
                missing=raw['missing'].astype(bool)
                for method,key in [('KF-CV','KF_CV'),('KF-GM','KF_GM')]:
                    nis=raw['nis_'+key]
                    total=float(np.nansum(nis[~missing]))
                    for layer in ['all','switch_event','switch_near2','G-clean','G-close','G-ambiguous']:
                        mask=np.ones(point.shape,bool) if layer=='all' else point==layer
                        mask&=~missing[:,None]
                        values=nis[mask];values=values[np.isfinite(values)]
                        row={'case_id':case['case_id'],'seed':int(case['seed']),'D':D,'method':method,'geometry':layer,'n_edge_times':len(values)}
                        if len(values):row.update(mean=float(values.mean()),median=float(np.median(values)),p90=float(np.quantile(values,.9)),p95=float(np.quantile(values,.95)),p99=float(np.quantile(values,.99)),sum_NIS=float(values.sum()),fraction_total_NIS=float(values.sum()/total))
                        else:row.update(mean=np.nan,median=np.nan,p90=np.nan,p95=np.nan,p99=np.nan,sum_NIS=np.nan,fraction_total_NIS=np.nan)
                        rows.append(row)
                    e=raw['theta_'+key][200:200+D]-raw['truth'][200:200+D];e=np.arctan2(np.sin(e),np.cos(e))
                    p=raw['p11_'+key][200:200+D]
                    for layer in ['all','G-clean','G-switch','G-close','G-ambiguous']:
                        mask=np.ones(140,bool) if layer=='all' else labels==layer
                        for lo,hi in zip([0,1e-5,1e-4,.001,.01,.1,1],[1e-5,1e-4,.001,.01,.1,1,np.inf]):
                            ok=(p[:,mask]>=lo)&(p[:,mask]<hi)
                            pp=p[:,mask][ok];ee=e[:,mask][ok]
                            if len(pp):bins.append({'case_id':case['case_id'],'seed':int(case['seed']),'D':D,'method':method,'geometry':layer,
                                'P_lo':lo,'P_hi':hi,'n':len(pp),'mean_P':float(pp.mean()),'mean_squared_error':float(np.mean(ee*ee)),
                                'coverage':float((np.abs(ee)<=1.96*np.sqrt(pp)).mean())})
    frame=pd.DataFrame(rows);frame.to_csv(target,index=False)
    pd.DataFrame(bins).to_csv(out/'p_bins_by_geometry.csv',index=False)
    summary=frame.groupby(['method','geometry'],as_index=False).agg(n_edge_times=('n_edge_times','sum'),sum_NIS=('sum_NIS','sum'))
    totals=summary[summary.geometry=='all'].set_index('method')
    summary['fraction_edge_times']=[r.n_edge_times/totals.loc[r.method].n_edge_times for _,r in summary.iterrows()]
    summary['fraction_NIS_sum']=[r.sum_NIS/totals.loc[r.method].sum_NIS for _,r in summary.iterrows()]
    summary.to_csv(out/'nis_mass_contribution.csv',index=False)
    print(summary.to_string(index=False))


if __name__=='__main__':run(Path(sys.argv[1]),Path(sys.argv[2]))

"""Independent position algebra, source pairing, raw hashes and RMSE audit."""
import hashlib,json,sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.io import loadmat


def sha(path):return hashlib.file_digest(path.open('rb'),'sha256').hexdigest()


def position_labels(x,starts,D):
    p=x['positions'];i=x['observer_id'].astype(int)-1;j=x['target_id'].astype(int)-1
    dx=(p[:,j,:]-p[:,i,:]).transpose(2,1,0)
    rounded=np.copysign(np.floor(np.abs(dx/1000)+.5),dx)
    delta=dx-1000*rounded
    distance=np.sqrt(np.sum(delta**2,axis=-1))
    jump=np.zeros((2200,140),bool);jump[1:]=(np.abs(delta[1:]-delta[:-1])>500).any(axis=2)
    angle=np.arctan2(delta[...,1],delta[...,0]);change=np.zeros_like(angle)
    change[1:]=np.arctan2(np.sin(angle[1:]-angle[:-1]),np.cos(angle[1:]-angle[:-1]))
    near=distance<5;boundary=np.min(500-np.abs(delta),axis=-1)<10
    alias=np.abs(change)>=.9*np.pi
    label=[]
    for t0 in starts:
        s=int(t0)-1;z=np.full(140,'G-ambiguous',dtype='U16')
        switch=jump[s:s+D].any(axis=0);close=near[s:s+D].any(axis=0)
        clean=~jump[s-200:s+D].any(axis=0)&~alias[s-200:s+D].any(axis=0)&~boundary[s:s+D].any(axis=0)&~close
        z[clean]='G-clean';z[close]='G-close';z[switch]='G-switch';label.append(z)
    actual=x['true_bearing'];error=np.arctan2(np.sin(angle-actual),np.cos(angle-actual))
    assert np.abs(error).max()<2e-12
    return np.array(label)


def audit(stage):
    out=stage/'review';assert not out.exists();out.mkdir()
    inputs=json.loads((stage/'input_manifest.json').read_text(encoding='utf-8'))
    raw_manifest=json.loads((stage/'raw/replay/manifest.json').read_text(encoding='utf-8'))
    assert inputs['status']=='COMPLETED' and raw_manifest['status']=='COMPLETED'
    if 'development' in stage.name:assert set(map(int,inputs['seeds']))==set(range(2026111001,2026111013))
    else:assert not set(map(int,inputs['seeds']))&set(range(2026111001,2026111013))
    assert len(inputs['inputs'])==2*len(inputs['seeds'])
    for item in raw_manifest['files']:assert sha(stage/'raw/replay'/item['file'])==item['sha256']
    table=pd.read_csv(stage/'analysis/window_metrics.csv.gz')
    table=table.set_index(['seed','scenario','method','D','t0','geometry'])
    definitions=stage/('GEOMETRY_FROZEN.json' if 'development' in stage.name else 'GEOMETRY_APPLIED.json')
    cuts=json.loads(definitions.read_text(encoding='utf-8'))['rate_cutpoints']
    methods=sorted({a['method'] for a in raw_manifest['files']})
    n_rows=0;maximum_error=0.;min_eigen=np.inf
    for item in inputs['inputs']:
        path=stage/item['source'];assert sha(path)==item['sha256']
        x=loadmat(path,simplify_cells=True);seed=int(item['seed']);scenario=item['scenario'];stem=path.stem
        assert x['observed_bearing'].shape==(2200,140)
        assert np.array_equal(x['observer_id'],np.repeat(np.arange(1,21),7))
        if scenario=='turn':
            base=loadmat(stage/'raw/input'/f'seed{seed}_base.mat',simplify_cells=True)
            assert np.array_equal(base['positions'][:,:,:1401],x['positions'][:,:,:1401])
            assert np.array_equal(base['headings'][:,:1401],x['headings'][:,:1401])
            noise_a=np.arctan2(np.sin(base['observed_bearing']-base['true_bearing']),np.cos(base['observed_bearing']-base['true_bearing']))
            noise_b=np.arctan2(np.sin(x['observed_bearing']-x['true_bearing']),np.cos(x['observed_bearing']-x['true_bearing']))
            assert np.abs(noise_a-noise_b).max()<3e-15
        for D in [20,120]:
            with np.load(stage/'labels'/f'{stem}_D{D}.npz') as label_file:
                labels=label_file['labels'];rates=label_file['rates'];starts=label_file['starts']
            assert np.array_equal(starts,np.arange(1201,1682,20))
            assert np.array_equal(labels,position_labels(x,starts,D))
            truth=x['true_bearing'][starts[None,:]-1+np.arange(D)[:,None]]
            for method in methods:
                with np.load(stage/'raw/replay'/f'{stem}_{method}_fresh.npz') as full:
                    with np.load(stage/'raw/replay'/f'{stem}_{method}_D{D}.npz') as branch:
                        indices=starts[None,:]-1+np.arange(branch['theta_xor'].shape[0])[:,None]
                        reference=np.ascontiguousarray(full['theta'][indices,0]).view(np.uint64)
                        theta=np.bitwise_xor(branch['theta_xor'],reference).view(np.float64)[:D]
                        pref=np.ascontiguousarray(full['P'][indices,0]).view(np.uint64)
                        covariance=np.bitwise_xor(branch['P_xor'],pref).view(np.float64)
                        if method!='Hold':
                            eigen=float(np.linalg.eigvalsh(covariance).min());min_eigen=min(min_eigen,eigen);assert eigen>=-1e-12
                err=np.arctan2(np.sin(theta-truth),np.cos(theta-truth))
                for c,t0 in enumerate(starts):
                    lb=labels[c];rt=rates[c];clean=lb=='G-clean'
                    masks={'all':np.ones(140,bool),**{g:lb==g for g in ['G-clean','G-switch','G-close','G-ambiguous']},
                        'clean-slow':clean&(rt<=cuts[0]),'clean-mid':clean&(rt>cuts[0])&(rt<=cuts[1]),'clean-fast':clean&(rt>cuts[1])}
                    for group,mask in masks.items():
                        row=table.loc[(seed,scenario,method,D,int(t0),group)]
                        assert row.n_edges==mask.sum()
                        if not mask.any():assert np.isnan(row.rmse);continue
                        expected=float(np.sqrt(np.mean(err[:,c,mask]**2)))
                        difference=abs(expected-row.rmse);maximum_error=max(maximum_error,difference)
                        assert difference<1e-12
                        n_rows+=1
        print('AUDIT PASS',stem,flush=True)
    report={'status':'PASS','all_raw_sha_verified':len(raw_manifest['files']),'physical_inputs':len(inputs['inputs']),
        'independent_seeds':len(inputs['seeds']),'nonempty_window_RMSE_rows_independently_recomputed':n_rows,
        'max_RMSE_difference':maximum_error,'minimum_covariance_eigenvalue':min_eigen,
        'position_labels_independently_recomputed':True,'base_turn_prefix_and_noise_matched':True}
    (out/'independent_audit.json').write_text(json.dumps(report,indent=2),encoding='utf-8');print(json.dumps(report))


if __name__=='__main__':audit(Path(sys.argv[1]))

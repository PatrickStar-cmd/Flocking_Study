"""Batched warm starts and exact compressed archives vs untouched P2 filters."""
import json
from pathlib import Path
import sys
import numpy as np
from engine import run,archive_branch,recover_branch

sys.path.insert(0,str(Path(__file__).resolve().parent.parent/'kf-p2-20261102'))
from replay import replay


def main(out):
    assert not out.exists();out.mkdir(parents=True)
    rng=np.random.default_rng(2718);obs=np.mod(np.cumsum(rng.normal(0,.02,(500,4)),axis=0)+1.,2*np.pi)
    max_error=0.
    for index,(q,tau) in enumerate([(x,None) for x in [1e-6,1e-5,1e-4,1e-3,1e-2]]+[(1e-5,10)]):
        fresh=run(obs,np.zeros(500,bool),q=q,tau=tau)
        starts=np.array([201,251]);indices=starts[None,:]-1+np.arange(220)[:,None]
        initial={k:fresh[k][starts-2,0] for k in ['theta','omega','P']}
        missing=np.arange(220)<120
        branch=run(obs[indices],missing,q=q,tau=tau,initial=initial)
        for b,start in enumerate(starts):
            mask=np.zeros(500,bool);mask[start-1:start+119]=True
            old=replay(obs,mask,dt=.3,q=q,R=1e-4,tau=tau)
            for key,newname,col in [('theta','theta',None),('p11','P',(0,0)),('p12','P',(0,1)),('p22','P',(1,1)),('nis','nis',None)]:
                new=branch[newname][:,b]
                if col is not None:new=np.broadcast_to(new[:,col[0],col[1],None],(220,4))
                expected=old[key][indices[:,b]]
                assert np.allclose(new,expected,atol=1e-10,rtol=1e-10,equal_nan=True),(q,tau,key)
                max_error=max(max_error,float(np.nanmax(np.abs(new-expected))))
        path=out/f'model{index}.npz'
        archive_branch(path,branch,fresh,indices,starts,initial)
        restored=recover_branch(path,fresh)
        for key in branch:np.testing.assert_array_equal(restored[key],branch[key])
        covariance=branch['P'];assert np.linalg.eigvalsh(covariance).min()>-1e-12
        changed=obs[indices].copy();changed[missing]=-9999.
        poisoned=run(changed,missing,q=q,tau=tau,initial=initial)
        for key in branch:np.testing.assert_array_equal(poisoned[key],branch[key])
    report={'status':'PASS','models_checked':6,'warm_start_cases_per_model':2,
        'all_q_Joseph_PSD':True,'missing_observation_invariance':True,
        'lossless_XOR_roundtrip':True,'max_difference_from_untouched_P2':max_error}
    (out/'audit.json').write_text(json.dumps(report,indent=2),encoding='utf-8');print(json.dumps(report))


if __name__=='__main__':main(Path(sys.argv[1]))

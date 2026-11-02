"""Independent G1 raw-data and mechanism audit before new physical seeds."""
import hashlib,json,sys
from pathlib import Path
import numpy as np
import pandas as pd


def main(root):
    path=root/'audit.json';assert not path.exists()
    manifest=json.loads((root/'manifest.json').read_text(encoding='utf-8'))
    assert manifest['status']=='COMPLETED' and manifest['total_conditions']==320
    summary=pd.read_csv(root/'synthetic_summary.csv')
    assert len(summary)==960
    noises={};count=0
    for entry in manifest['runs']:
        f=root/'raw'/entry['npz_file']
        assert hashlib.file_digest(f.open('rb'),'sha256').hexdigest()==entry['npz_sha256']
        with np.load(f,allow_pickle=False) as raw:
            D=int(raw['D']);seed=int(raw['seed']);config=str(raw['config_id'])
            assert raw['truth'].shape==(720,1)
            missing=raw['missing'].ravel();assert missing.sum()==D and missing[200:200+D].all()
            noise=np.arctan2(np.sin(raw['observations']-raw['truth']),np.cos(raw['observations']-raw['truth']))
            key=(seed,config)
            if key in noises:assert np.allclose(noises[key],noise,atol=2e-15,rtol=0)
            else:noises[key]=noise
            if config=='D1':assert raw['near_range'][200:200+D].any() and not raw['image_switch'].any() and raw['distance'].min()>=2
            if config=='E1':assert raw['image_switch'][200:200+D].sum()==1 and not raw['near_range'].any()
            for method,name in [('Hold','Hold'),('KF-CV','KF_CV'),('KF-GM','KF_GM')]:
                e=raw['theta_'+name][200:200+D]-raw['truth'][200:200+D]
                e=np.arctan2(np.sin(e),np.cos(e))
                r=summary[(summary.seed==seed)&(summary.config_id==config)&(summary.D==D)&(summary.method==method)].iloc[0]
                assert abs(np.sqrt(np.mean(e*e))-r.circular_rmse)<1e-12
                if method=='Hold':assert np.isnan(raw['p11_Hold']).all()
                count+=1
    report={'status':'PASS','raw_conditions':320,'summary_rows_recomputed':count,'noise_matched_across_D':True,
        'D_close_passage_verified':True,'E_single_image_switch_verified':True,
        'evidence':'software and synthetic mechanisms, not physical group validation'}
    path.write_text(json.dumps(report,indent=2),encoding='utf-8');print(json.dumps(report))


if __name__=='__main__':main(Path(sys.argv[1]))

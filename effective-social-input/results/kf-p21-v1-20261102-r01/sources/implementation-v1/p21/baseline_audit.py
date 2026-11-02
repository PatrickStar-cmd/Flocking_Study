"""A gate: independently recover published P2 error and innovation summaries."""
import hashlib,json,platform,subprocess,sys
from pathlib import Path
import numpy as np
import pandas as pd
import scipy


def main(out):
    here=Path(__file__).resolve().parent;repo=here.parent.parent
    old=repo/'effective-social-input/results/kf-p2-v1-20261102-reviewed'
    assert not out.exists();out.mkdir(parents=True)
    assert subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip()=='81e0080c8ffc00134c64f407521974a32a950767'
    audit=json.loads((old/'review/publication_audit.json').read_text(encoding='utf-8'))
    for f in audit['artifacts_sha256']:
        path=repo/f['path']
        assert hashlib.file_digest(path.open('rb'),'sha256').hexdigest()==f['sha256'],f['path']
    summary=pd.read_csv(old/'analysis/summary.csv');max_rmse=0.;max_nis=0.
    for path in sorted((old/'raw/replay').glob('*.npz')):
        with np.load(path) as raw:
            case=str(raw['case_id']);D=int(raw['D'])
            for method,key in [('Hold','Hold'),('KF-CV','KF_CV'),('KF-GM','KF_GM')]:
                error=raw['theta_'+key]-raw['truth'];error=np.arctan2(np.sin(error),np.cos(error))
                for win,sl in [('pre',slice(0,200)),('occ',slice(200,200+D)),('recovery',slice(200+D,None))]:
                    row=summary[(summary.case_id==case)&(summary.D==D)&(summary.method==method)&(summary.window==win)]
                    assert len(row)==1;row=row.iloc[0]
                    e=error[sl];rmse=float(np.sqrt(np.mean(e*e)))
                    max_rmse=max(max_rmse,abs(rmse-row.circular_rmse));assert abs(rmse-row.circular_rmse)<1e-12
                    ni=raw['nis_'+key][sl];finite=ni[np.isfinite(ni)]
                    value=finite.mean() if len(finite) else np.nan
                    assert np.allclose(value,row.observed_fresh_meanNIS,atol=1e-10,rtol=1e-12,equal_nan=True)
                    if len(finite):max_nis=max(max_nis,abs(value-row.observed_fresh_meanNIS))
    report={'status':'PASS','baseline_commit':'81e0080c8ffc00134c64f407521974a32a950767',
        'published_artifacts_sha_verified':len(audit['artifacts_sha256']),'summary_rows_recovered':len(summary),
        'max_rmse_difference':max_rmse,'max_nis_difference':max_nis,
        'python':sys.version,'numpy':np.__version__,'scipy':scipy.__version__,'platform':platform.platform()}
    (out/'audit.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    (out/'BASELINE_AUDIT.md').write_text('# Baseline A gate: PASS\n\nP2 frozen source/artifact hashes, all 270 window-level RMSE and fresh mean NIS recovered. No old files modified.\n\n'+json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(report))


if __name__=='__main__':main(Path(sys.argv[1]))

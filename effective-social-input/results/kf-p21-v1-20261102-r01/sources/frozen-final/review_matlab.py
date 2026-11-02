import json,sys
from pathlib import Path
import numpy as np
from scipy.io import loadmat
from engine import run


def main(root):
    reference=loadmat(root/'matlab_q_reference.mat',simplify_cells=True)
    maxima={key:0. for key in ['theta','omega','P','nis','age']}
    for ref in reference['refs']:
        actual=run(reference['obs'],reference['missing'].astype(bool),q=ref['q'],tau=ref['tau'])
        for key in maxima:
            expected=ref[key]
            observed=actual[key][:,0] if key!='age' else actual[key]
            if key=='P':expected=expected[:,0]
            if key=='age':observed=actual[key][:,None];observed=np.broadcast_to(observed,expected.shape)
            assert np.allclose(observed,expected,atol=2e-10,rtol=1e-10,equal_nan=True),(key,ref['q'])
            maxima[key]=max(maxima[key],float(np.nanmax(np.abs(observed-expected))))
    report={'status':'PASS','models':6,'steps':2200,'edges':3,'max_absolute_differences':maxima}
    path=root/'matlab_audit.json';assert not path.exists();path.write_text(json.dumps(report,indent=2),encoding='utf-8');print(json.dumps(report))


if __name__=='__main__':main(Path(sys.argv[1]))

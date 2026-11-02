"""Flag Gaussian angle intervals wider than a half-circle, without dropping data."""
import json,sys
from pathlib import Path
import numpy as np
import pandas as pd


def main(stage):
    output=stage/'analysis/uncertainty_domain.csv';assert not output.exists()
    manifest=json.loads((stage/'raw/replay/manifest.json').read_text(encoding='utf-8'))
    rows=[]
    for record in manifest['files']:
        if record['kind']!='XOR-branch' or record['method']=='Hold':continue
        D=record['D']
        with np.load(stage/'raw/replay'/record['fresh_file']) as fresh:
            with np.load(stage/'raw/replay'/record['file']) as branch:
                starts=branch['starts'].astype(int);ix=starts[None,:]-1+np.arange(D)[:,None]
                reference=np.ascontiguousarray(fresh['P'][ix,0]).view(np.uint64)
                P=np.bitwise_xor(branch['P_xor'][:D],reference).view(np.float64)[...,0,0]
        for lo,hi in [(1,20),(21,60),(61,100),(101,120)]:
            if lo>D:continue
            block=P[lo-1:min(hi,D)]
            rows.append({'source':record['source'],'method':record['method'],'D':D,'age_band':f'{lo}-{hi}',
                'candidate_window_steps':block.size,'gaussian_halfwidth_ge_pi_fraction':float((1.96*np.sqrt(block)>=np.pi).mean()),
                'max_gaussian_halfwidth_rad':float(1.96*np.sqrt(block.max())),'min_P':float(block.min()),'max_P':float(block.max())})
    pd.DataFrame(rows).to_csv(output,index=False)
    print('Wrote',output,len(rows),flush=True)


if __name__=='__main__':main(Path(sys.argv[1]))

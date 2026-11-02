"""NIS tails for the actual masked branches, separately pre and recovery.

Candidate windows overlap and are not additional independent observations.
Each output row is clustered by the physical seed; missing frames have no NIS.
"""
import json,sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.io import loadmat
from geometry import geometry,fresh_labels


def main(stage):
    target=stage/'analysis/branch_nis_tail_by_geometry_seed.csv';assert not target.exists()
    manifest=json.loads((stage/'input_manifest.json').read_text(encoding='utf-8'))
    rawmanifest=json.loads((stage/'raw/replay/manifest.json').read_text(encoding='utf-8'))
    methods=sorted({e['method'] for e in rawmanifest['files']}-{'Hold'})
    rows=[]
    layers=['all','switch_event','switch_near2','G-clean','G-close','G-ambiguous']
    for seed in manifest['seeds']:
        entries=[e for e in manifest['inputs'] if e['seed']==seed]
        sources=[]
        for e in entries:
            p=stage/e['source'];x=loadmat(p,simplify_cells=True)
            sources.append((p.stem,fresh_labels(geometry(x))))
        for method in methods:
            for D in [20,120]:
                pools={(phase,g):[] for phase in ['pre','recovery'] for g in layers}
                for stem,labels in sources:
                    with np.load(stage/'raw/replay'/f'{stem}_{method}_fresh.npz') as fresh:
                        with np.load(stage/'raw/replay'/f'{stem}_{method}_D{D}.npz') as branch:
                            starts=branch['starts'].astype(int)
                            indices=starts[None,:]-1+np.arange(D+400)[:,None]
                            reference=np.ascontiguousarray(fresh['nis'][indices,0]).view(np.uint64)
                            nis=np.bitwise_xor(branch['nis_xor'],reference).view(np.float64)
                            before=starts[None,:]-201+np.arange(200)[:,None]
                            pair={'pre':(fresh['nis'][before,0],labels[before]),'recovery':(nis[D:],labels[indices[D:]])}
                            for phase,(values,fl) in pair.items():
                                for layer in layers:
                                    a=values.ravel() if layer=='all' else values[fl==layer]
                                    a=a[np.isfinite(a)]
                                    if len(a):pools[(phase,layer)].append(a)
                for (phase,layer),parts in pools.items():
                    values=np.concatenate(parts) if parts else np.array([])
                    r={'seed':int(seed),'D':D,'method':method,'phase':phase,'geometry':layer,'n_window_edge_times':len(values)}
                    if len(values):r.update(mean=float(values.mean()),median=float(np.median(values)),p90=float(np.quantile(values,.9)),p95=float(np.quantile(values,.95)),p99=float(np.quantile(values,.99)),max_NIS=float(values.max()))
                    else:r.update(mean=np.nan,median=np.nan,p90=np.nan,p95=np.nan,p99=np.nan,max_NIS=np.nan)
                    rows.append(r)
        print('PASS actual branch NIS',int(seed),flush=True)
    pd.DataFrame(rows).to_csv(target,index=False)
    print('Wrote',target,len(rows),flush=True)


if __name__=='__main__':main(Path(sys.argv[1]))

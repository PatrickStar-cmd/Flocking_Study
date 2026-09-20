"""Aggregate block metrics by seed and run endpoint sanity tests."""
import csv, json
from pathlib import Path

ROOT=Path(__file__).resolve().parent
SRC=ROOT/'results/block-metric-pilot/block-metrics.csv'
OUT=ROOT/'results/block-metric-pilot/block-summary.json'

def weighted(rows, key, target):
    total=0; passed=0; blocks=[]
    for r in rows:
        n=int(r['blockSize']); ok=float(r[key])
        total += n; passed += n*ok; blocks.append(ok)
    return {'weightedPassFraction': passed/total if total else 0,
            'blocks':len(blocks),'steps':total,'passBlocks':sum(blocks)}

def synthetic():
    # A block passes if its mean metric is successful; final endpoint is the
    # duration-weighted fraction of successful blocks/steps.
    cases={
      'complete_success':[True]*20,
      'complete_failure':[False]*20,
      'short_failure':[True]*9+[False]+[True]*10,
      'terminal_failure':[True]*19+[False],
    }
    return {name:{'weightedPassFraction':sum(x)/len(x),'expected':exp}
            for (name,x),exp in zip(cases.items(),[1,0,.95,.95])}

def main():
    with SRC.open(encoding='utf-8-sig',newline='') as f: raw=list(csv.DictReader(f))
    groups={}
    for r in raw:
        key=(r['seed'],r['policy'],r['strengthMode'],r['k'],r['blockSize'])
        groups.setdefault(key,[]).append(r)
    summary=[]
    for key,rows in sorted(groups.items()):
        q=weighted(rows,'passBlocks',.9)
        summary.append(dict(seed=int(key[0]),policy=key[1],strengthMode=key[2],k=int(key[3]),blockSize=int(key[4]),**q))
    result={'source':str(SRC),'rows':len(raw),'seedConditionBlockGroups':len(summary),
            'summary':summary,'syntheticTests':synthetic(),
            'limitations':['The source pilot thresholds are exploratory and were not independently calibrated for these block means.',
                           'Block rows are aggregated by seed; no block is treated as an independent replicate.',
                           'A final endpoint must specify task-based thresholds and a scientifically justified block duration.']}
    OUT.write_text(json.dumps(result,indent=2),encoding='utf-8')
    print(json.dumps(result['syntheticTests'],indent=2))
    for r in summary: print(r)
if __name__=='__main__': main()

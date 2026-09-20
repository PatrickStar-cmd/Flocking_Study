"""Audit the N20 sweep without modifying simulations or source results."""
import csv
import hashlib
import json
import math
import random
from collections import defaultdict
from pathlib import Path
from statistics import mean, median

ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / 'results/effective-input-n20-100.csv'
OUT = ROOT / 'results/n20-100-statistics'
FIELDS = ['meanGO', 'meanLO', 'meanPairDistanceNorm']


def quantile(x, p):
    x = sorted(x)
    t = (len(x)-1)*p
    i = int(t)
    return x[i] + (x[min(i+1, len(x)-1)]-x[i])*(t-i)


def calibrate(base):
    cols = list(zip(*base))
    limits = [quantile(cols[0], .05), quantile(cols[1], .05), quantile(cols[2], .95)]
    centers = [median(c) for c in cols]
    scales = [max(abs(c-t), 2.220446049250313e-16) for c,t in zip(centers, limits)]
    def score(v):
        return min((v[0]-limits[0])/scales[0], (v[1]-limits[1])/scales[1], (limits[2]-v[2])/scales[2])
    cutoff = sorted(map(score, base))[max(0, math.ceil(len(base)*.05)-1)]
    return limits, scales, cutoff, score


def wilson(s, n):
    z = 1.96
    p = s/n
    den = 1+z*z/n
    c = (p+z*z/(2*n))/den
    h = z*math.sqrt(p*(1-p)/n+z*z/(4*n*n))/den
    return max(0,c-h), min(1,c+h)


def isotonic(values):
    blocks = []
    for v in values:
        blocks.append([v, 1])
        while len(blocks)>1 and blocks[-2][0]>blocks[-1][0]:
            b, a = blocks.pop(), blocks.pop()
            blocks.append([(a[0]*a[1]+b[0]*b[1])/(a[1]+b[1]), a[1]+b[1]])
    return [v for v,n in blocks for _ in range(n)]


def first(ks, rates):
    return next((k for k,p in zip(ks,rates) if p>=.95), None)


def write(name, rows):
    with (OUT/name).open('w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def main():
    with SOURCE.open(encoding='utf-8-sig', newline='') as f:
        raw = list(csv.DictReader(f))
    unique = {}
    duplicates = 0
    for row in raw:
        key = tuple(row[c] for c in ['N','seed','policy','strengthMode','k'])
        assert all(math.isfinite(float(row[c])) for c in FIELDS)
        if key in unique:
            assert row == unique[key], ('Nonidentical duplicate', key)
            assert row['k']=='0'
            duplicates += 1
        else:
            unique[key] = row
    groups = defaultdict(dict)
    for row in unique.values():
        assert (row['N'],row['nNeurons'],row['nSteps']) == ('20','100','8000')
        groups[(row['policy'],row['strengthMode'],int(row['k']))][int(row['seed'])] = tuple(float(row[c]) for c in FIELDS)
    assert len(groups)==121
    assert all(set(g)==set(range(1,101)) for g in groups.values())
    base = groups[('all','original',19)]
    limits, scales, cutoff, score = calibrate(list(base.values()))
    summary = []
    for key,g in sorted(groups.items()):
        vals = list(g.values())
        s = sum(score(v)>=cutoff for v in vals)
        lo,hi = wilson(s,100)
        row = dict(policy=key[0],strengthMode=key[1],k=key[2],trials=100,successes=s,successRate=s/100,wilsonLower95=lo,wilsonUpper95=hi)
        for j,name in enumerate(FIELDS):
            row[name] = mean(v[j] for v in vals)
            single = sum(v[j]>=limits[j] if j<2 else v[j]<=limits[j] for v in vals)
            row[name+'Passes'] = single
            row[name+'WilsonLower95'],row[name+'WilsonUpper95'] = wilson(single,100)
            row[name+'PairedDifference'] = mean(g[seed][j]-base[seed][j] for seed in base)
        row['strictJointSuccesses'] = sum(v[0]>=limits[0] and v[1]>=limits[1] and v[2]<=limits[2] for v in vals)
        summary.append(row)
    curves = defaultdict(list)
    for row in summary:
        curves[(row['policy'],row['strengthMode'])].append(row)
    thresholds = []
    for key,rows in curves.items():
        ks = [r['k'] for r in rows]
        fitted = isotonic([r['successRate'] for r in rows])
        for r,p in zip(rows,fitted):
            r['isotonicSuccessRate'] = p
        thresholds.append(dict(policy=key[0],strengthMode=key[1],empirical=first(ks,[r['successRate'] for r in rows]),isotonic=first(ks,fitted),wilsonLower95=first(ks,[r['wilsonLower95'] for r in rows])))
    # Resample entire seeds together, including baseline calibration and all k.
    rng = random.Random(20260907)
    bootstrap = defaultdict(list)
    for _ in range(1000):
        seeds = rng.choices(range(1,101),k=100)
        _,_,c,sc = calibrate([base[s] for s in seeds])
        for key,rows in curves.items():
            ks = [r['k'] for r in rows]
            rates = [sum(sc(groups[(*key,k)][s])>=c for s in seeds)/100 for k in ks]
            bootstrap[key].append(first(ks,isotonic(rates)))
    for row in thresholds:
        samples = bootstrap[(row['policy'],row['strengthMode'])]
        found = [v for v in samples if v is not None]
        row['bootstrapReplicates'] = len(samples)
        row['bootstrapNoCrossingFraction'] = 1-len(found)/len(samples)
        row['conditionalBootstrapP025'] = quantile(found,.025) if found else None
        row['conditionalBootstrapP975'] = quantile(found,.975) if found else None
    OUT.mkdir(exist_ok=True)
    write('group-success.csv',summary)
    write('thresholds.csv',thresholds)
    audit = dict(source=str(SOURCE),sha256=hashlib.sha256(SOURCE.read_bytes()).hexdigest(),rawRows=len(raw),identicalDuplicatesRemoved=duplicates,uniqueRows=len(unique),groups=len(groups),seedsPerGroup=100,baselineMarginalLimits=limits,scales=scales,jointScoreCutoff=cutoff,effectiveJointLimits=[limits[0]+scales[0]*cutoff,limits[1]+scales[1]*cutoff,limits[2]-scales[2]*cutoff],bootstrapSeed=20260907,limitations=['Wilson intervals condition on a threshold estimated from these same seeds; they do not include calibration uncertainty or multiple comparisons.', 'Bootstrap resamples seeds and recalibrates; conditional quantiles omit replicates with no crossing and are not a simultaneous confidence set.', 'Isotonic fits assume monotonic success, which is not guaranteed by neural dynamics.', 'Time-averaged CSV cannot establish sustained success, collisions or failure times.'])
    (OUT/'audit.json').write_text(json.dumps(audit,indent=2),encoding='utf-8')
    print(json.dumps(audit,indent=2))
    print(json.dumps(thresholds,indent=2))
    print('BASELINE',summary[0])


if __name__ == '__main__':
    main()

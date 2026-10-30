"""Constructed-table tests only. This file never calls a model simulator."""
from pathlib import Path
import ast,json,math
import numpy as np
import pandas as pd
import analyze

ROOT=Path(__file__).resolve().parent; cfg=analyze.CFG
profiles=[cfg['primary'],*cfg['history_controls']]
rows=[]
for p in profiles:
    for seed in range(p['seeds_first'],p['seeds_first']+p['seed_count']):
        for r in cfg['fresh_counts']:
            # Artificial fixture, an algebraic sequence rather than a model response.
            c=.8*r/19
            rows.append(dict(seed=seed,profile=p['profile'],remaining=r,initialSocialInput=.012*r,
                finalConsistency=c,finalGO=r/19,finalLO=c,finalConsistencyQ10=c,
                finalJointHighFraction=float(c>=.9),finalMovingFraction=1,finalMeanSpeed=1,
                finalPairDistance=.2,C_last500=c,C_last2000=c))
df=pd.DataFrame(rows); analyze.validate_rows(df,profiles); assert len(df)==14000
a=analyze.paired_aggregate(df,profiles,repeats=50)
assert len(a)==60
assert np.max(abs(a.finalConsistency-.8*a.remaining/19))<1e-12
assert np.max(abs(a.finalConsistency_lo-a.finalConsistency_hi))<1e-12
independent_mean=math.fsum(df[df.profile=='short_history'].finalConsistency)/10000
assert abs(independent_mean-.4)<1e-14
rejected=0
for bad in [pd.concat([df,df.iloc[:1]]),df.iloc[1:],df.assign(initialSocialInput=0)]:
    try:analyze.validate_rows(bad,profiles)
    except AssertionError:rejected+=1
assert rejected==3
for name in ['analyze.py','test_preparation.py']:
    ast.parse((ROOT/name).read_text(encoding='utf-8'))
assert not (ROOT/'results').exists(),'Preparation should create zero model outputs.'
check=ROOT/'checks';check.mkdir(exist_ok=True)
report=dict(status='PASS',evidence='constructed metric and aggregation fixtures only',
    model_steps_executed=0,model_trajectories_executed=0,artificial_rows=14000,
    planned_primary_seeds=500,planned_primary_trajectories=10000,planned_control_trajectories=4000,
    aggregate_conditions=60,rejected_invalid_tables=rejected,independent_fixture_mean=independent_mean,
    input_axis_endpoints=[0,.228],coherence_endpoint='mean of per-time min(GO,LO)',
    Python_source_parse='PASS')
(check/'python-preparation.json').write_text(json.dumps(report,indent=2),encoding='utf-8',newline='\n')
print(json.dumps(report,indent=2))

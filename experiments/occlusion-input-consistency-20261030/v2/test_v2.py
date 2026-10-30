"""Algebraic and constructed-data checks. No neural model is executed."""
from pathlib import Path
import ast,json,math,tempfile
import numpy as np
import pandas as pd
from scipy.io import savemat
from three_factor_input import social_input,checkpoint_input
import analyze_v2 as analysis

HERE=Path(__file__).resolve().parent
w=np.full((3,4),.012);theta=np.repeat(np.array([0,np.pi/2,np.pi])[:,None],4,axis=1);phi=np.zeros(3)
z=social_input(w,theta,phi)
expected=.048*np.exp(np.array([0,-1,-2]))
assert np.allclose(z['S'],expected,atol=1e-15,rtol=0)
assert np.allclose(z['S'],z['n']*z['mean_weight']*z['angle_weight'])
assert np.allclose(social_input(2*w,theta,phi)['S'],2*z['S'])
assert np.allclose(social_input(np.tile(w,(1,2)),np.tile(theta,(1,2)),phi)['S'],2*z['S'])
assert np.allclose(social_input(w,theta+.7,phi+.7)['S'],z['S'])
assert np.allclose(social_input(w,-theta,-phi)['S'],z['S'])
empty=social_input(np.zeros((1,4)),np.full((1,4),np.nan),np.zeros(1))
assert empty['S'][0]==empty['A'][0]==empty['angle_weight'][0]==0
try:social_input(w,theta,phi,kappa=0)
except AssertionError:pass
else:raise AssertionError('The three-factor form requires nonconstant angular weighting.')
model=analysis.base.CFG['model'];n=model['N'];k=model['full_neighbors']
cp=dict(position=np.stack([np.arange(n,dtype=float),np.zeros(n)]),heading=np.zeros(n),
    referenceAxis=np.full(n,np.pi/2),maskOrder=np.tile(np.arange(1,k+1),(n,1)))
measure=[checkpoint_input(cp,r,model)['S'].mean() for r in range(k+1)]
assert np.all(np.diff(measure)>0) and measure[0]==0
assert measure[-1]<k*model['single_edge_weight']
rows=[dict(seed=7,profile='constructed',remaining=r,initialSocialInput=.012*r,
    initialVectorX=0,initialVectorY=0,finalConsistency=.4) for r in range(k+1)]
with tempfile.TemporaryDirectory(prefix='three-factor-fixture-') as folder:
    folder=Path(folder);(folder/'seed-7').mkdir()
    savemat(folder/'seed-7/constructed-checkpoint.mat',{'checkpoint':cp})
    derived=analysis.revised_rows(folder,pd.DataFrame(rows),[])
    assert np.allclose(derived.initialSocialInput,measure)
    assert np.allclose(derived.initialTotalStrength,.012*derived.remaining)
fixture=[]
for profile in ['short_history','formation','long_history']:
    for seed in range(8):
        for r in range(20):
            fixture.append(dict(seed=seed,profile=profile,remaining=r,
                initialSocialInput=.012*r*(.2+.03*seed),finalConsistency=.4))
df=pd.DataFrame(fixture);bins=analysis.input_bins(df,repeats=100)
assert bins.cases.sum()==len(df) and np.allclose(bins.consistency_mean,.4)
assert np.allclose(bins.lo,.4) and np.allclose(bins.hi,.4)
for name in ['three_factor_input.py','analyze_v2.py','test_v2.py']:
    ast.parse((HERE/name).read_text(encoding='utf-8'))
assert not (HERE.parent/'results').exists()
out=HERE/'checks';out.mkdir(exist_ok=True)
report=dict(status='PASS',model_steps_executed=0,model_trajectories_executed=0,
    fixture_type='analytic bearings and constructed summary/checkpoint tables',
    four_neighbor_input_front_side_back=list(map(float,z['S'])),
    factorization=True,strength_homogeneity=True,count_additivity=True,frame_rotation_invariance=True,
    zero_input_exact=True,mask_monotonicity=True,event_heading_used=True,
    legacy_strength_retained=True,constructed_rows=len(df),bin_cases_retained=int(bins.cases.sum()),
    independent_fixture_mean=math.fsum(df.finalConsistency)/len(df),static_parse='PASS')
(out/'preparation.v2.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8',newline='\n')
print(json.dumps(report,indent=2))

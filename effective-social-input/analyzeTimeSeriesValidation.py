import json, glob
from pathlib import Path
from statistics import mean, median
import scipy.io as sio

ROOT=Path(__file__).resolve().parent
OUT=ROOT/'results/time-series-validation-20260909'
files=sorted(glob.glob(str(OUT/'seed-*-condition-*.mat')))
rows=[]
for f in files:
    x=sio.loadmat(f,squeeze_me=True,struct_as_record=False)
    s1=x['summary100']; s3=x['summary300']; r=x['result']; cfg=r.config
    def get(o,n): return float(getattr(o,n))
    rows.append({'file':Path(f).name,'seed':int(get(cfg,'seed')),'condition':int(Path(f).stem.split('-')[-1]),'pointFraction':get(s1,'pointFraction'),'weightedPass100':get(s1,'durationWeightedPass'),'weightedPass300':get(s3,'durationWeightedPass'),'meanGO':get(r.summary,'globalOrder'),'meanLO':get(r.summary,'localOrder'),'meanPairDistance':get(r.summary,'meanPairDistanceNorm'),'longestFailureTime':get(s1,'longestFailureTime'),'failureRuns':get(s1,'failureRuns'),'recoveries':get(s1,'recoveries')})
conditions={1:'all original k=19',2:'random original k=0',3:'random fixed-total k=7',4:'random fixed-total k=15',5:'balanced fixed-total k=15',6:'random original k=15'}
summary=[]
for c,n in conditions.items():
 q=[x for x in rows if x['condition']==c]
 summary.append({'condition':n,'n':len(q),'meanGO':mean(x['meanGO'] for x in q),'meanLO':mean(x['meanLO'] for x in q),'meanPairDistance':mean(x['meanPairDistance'] for x in q),'meanPointFraction':mean(x['pointFraction'] for x in q),'meanWeightedPass100':mean(x['weightedPass100'] for x in q),'meanWeightedPass300':mean(x['weightedPass300'] for x in q),'medianLongestFailureTime':median(x['longestFailureTime'] for x in q),'meanFailureRuns':mean(x['failureRuns'] for x in q),'meanRecoveries':mean(x['recoveries'] for x in q)})
out={'files':len(files),'expected':180,'conditions':summary,'rows':rows,'limitations':['Block results are descriptive and aggregated by seed; blocks are not independent replicates.','Thresholds were frozen from prior data but were originally defined on time means, so block-pass results are exploratory.','No collision or task-specific success endpoint is included.']}
(OUT/'analysis-summary.json').write_text(json.dumps(out,indent=2),encoding='utf-8')
for x in summary: print(x)

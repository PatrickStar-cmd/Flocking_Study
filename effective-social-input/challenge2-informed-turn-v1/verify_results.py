"""Independent checks: raw states -> math definitions -> saved CSV, no simulator imports."""
import argparse
import csv
import json
import math
from pathlib import Path
import numpy as np


def read(p):
    with p.open(encoding='utf-8-sig',newline='') as f:
        return list(csv.DictReader(f))


def verify(stage):
    meta=json.loads((stage/'metadata.json').read_text(encoding='utf-8'))
    c=meta['config'];start=c['warm_steps']+c['lead_steps'];L=c['arena']
    errors={k:0. for k in ['GO','LO','fixed5_distance','follower_turn_error_deg','S_oracle']}
    checked=0
    def distance(a,b):
        dx=b[0]-a[0];dy=b[1]-a[1]
        dx-=L*round(dx/L);dy-=L*round(dy/L)
        return math.hypot(dx,dy)
    def angle_error(a,b):return abs(math.atan2(math.sin(a-b),math.cos(a-b)))
    files=list((stage/'runs').glob('*.npz'))
    for p in files:
        z=np.load(p); rows=read(p.with_suffix('.csv'))
        xy=z['xy'];velocity=z['velocity'];head=z['headings'];selected=z['selected'];hidden=z['hidden']
        assert selected.sum()==300 and hidden.sum()==51 and np.all(hidden[3:,:3])
        assert np.all(selected.sum(1)==15) and not np.diag(selected).any()
        fixed=[sorted([j for j in range(20) if j!=i],key=lambda j:distance(xy[start,i],xy[start,j]))[:5] for i in range(20)]
        event=[r for r in rows if int(r['blackout'])==1]
        for r in event:
            assert int(r['fresh_edges'])==249
            mode=p.stem.split('_')[2]
            if mode in ['drop','drop-normalized']:assert int(r['active_edges'])==249
            else:assert int(r['active_edges'])==300  # All pilot horizons below configured cutoff.
        for step in sorted(set(list(range(start,len(rows),10))+[start,len(rows)-1])):
            r=rows[step];pos=xy[step+1]
            unit=[]
            for v in velocity[step]:
                speed=math.hypot(*v)
                unit.append((v[0]/speed,v[1]/speed) if speed>1e-12 else (0.,0.))
            GO=math.hypot(sum(v[0] for v in unit),sum(v[1] for v in unit))/20
            local=[]
            for i in range(20):
                neighbors=sorted(range(20),key=lambda j:distance(pos[i],pos[j]))[:6]
                local.append(math.hypot(sum(unit[j][0] for j in neighbors),sum(unit[j][1] for j in neighbors))/6)
            D=sum(distance(pos[i],pos[j]) for i in range(20) for j in fixed[i])/100
            command=math.radians(float(r['command_deg']))
            E=sum(math.degrees(angle_error(head[step+1,i],command)) for i in range(3,20))/17
            S=0.
            for i in range(3,20):
                for j in range(20):
                    if selected[i,j]:
                        dx=xy[step,j,0]-xy[step,i,0];dy=xy[step,j,1]-xy[step,i,1]
                        dx-=L*round(dx/L);dy-=L*round(dy/L)
                        theta=math.atan2(dy,dx)
                        S+=c['h0']*math.exp(math.cos(theta-head[step,i])-1.)
            vals={'GO':GO,'LO':sum(local)/20,'fixed5_distance':D,'follower_turn_error_deg':E,'S_oracle':S/17}
            for k,val in vals.items():
                errors[k]=max(errors[k],abs(val-float(r[k])))
            checked+=1
    assert all(v<1e-9 for v in errors.values()),errors
    manifest=read(stage/'manifest.csv')
    prefix_checks=0
    for r in manifest:
        a=np.load(stage/'runs'/(r['name']+'.npz'));b=np.load(stage/'runs'/(r['control']+'.npz'))
        assert np.array_equal(a['xy'][:start+1],b['xy'][:start+1])
        assert np.array_equal(a['noise'],b['noise'])
        for k in [k for k in a.files if k.startswith('checkpoint_')]:
            assert np.array_equal(a[k],b[k],equal_nan=True)
        prefix_checks+=1
    summary=read(stage/'analysis'/'summary.csv');gates=read(stage/'analysis'/'task-gates.csv')
    keys=[(r['seed'],r['T'],r['mode']) for r in summary]
    assert len(keys)==len(set(keys))
    mean_checks=0
    for r in gates:
        seed=r['seed'];T=float(r['T']);end=start+round(T/c['dt'])
        fresh=read(stage/'runs'/f'seed-{seed}_turn_fresh_T-0.csv')
        drop=read(stage/'runs'/f'seed-{seed}_turn_drop_T-{T:g}.csv')
        delta=sum(float(drop[t]['follower_turn_error_deg'])-float(fresh[t]['follower_turn_error_deg']) for t in range(start,end))/(end-start)
        assert abs(delta-float(r['drop_event_error_increase_deg']))<1e-10
        mean_checks+=1
    result=dict(passed=True,raw_runs=len(files),state_samples_checked=checked,
                prefix_noise_checks=prefix_checks,independent_paired_means=mean_checks,
                summary_unique_rows=len(keys),max_abs_errors=errors)
    target=stage/'analysis'/'independent-verification.json'
    assert not target.exists()
    target.write_text(json.dumps(result,indent=2),encoding='utf-8')
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--stage',type=Path,required=True)
    verify(p.parse_args().stage)

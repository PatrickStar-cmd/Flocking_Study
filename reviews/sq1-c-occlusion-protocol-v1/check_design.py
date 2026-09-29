"""Check design arithmetic and document references, not simulation correctness."""
import hashlib
import json
import math
import csv
import io
from pathlib import Path

root=Path(__file__).resolve().parent
doc=json.loads((root/'protocol.json').read_text(encoding='utf-8'))
p=doc['proposal']
assert doc['user_confirmed']['fixed_target_identity']
assert doc['user_confirmed']['bearing_measurement_std']==0
assert not doc['execution']['allowed']
assert doc['revision']=='1.2'
assert doc['a_reference']['commit']=='992e71aba3b26896335b52a3182cf01cfd542a9d'
assert doc['a_reference']['requires_B_adapter_acceptance']
assert (root/doc['a_reference']['interface_contract']).is_file()
assert (root/doc['a_reference']['audit']).is_file()
assert doc['acceptance_gate_count']==15
assert (root/'B_DELIVERY_CHECKLIST.md').is_file()
assert sum(line.startswith('| G') for line in (root/'ACCEPTANCE.md').read_text(encoding='utf-8').splitlines())==15
assert p['event_start_local_step']==p['common_prefix_steps']+1
assert math.isclose(p['baseline_total_weight'],(p['N']-1)*p['ht']/p['N'])
assert 0<p['hidden_count']<p['k']<p['N']
intervals=[]
for d,t in zip(p['duration_steps'],p['duration_time']):
    assert math.isclose(d*p['dt'],t)
    start=p['event_start_local_step']; last=start+d-1; restore=last+1
    assert last-start+1==d
    assert restore+p['extended_post_steps']-1<=p['branch_steps']
    intervals.append(dict(d=d,T=t,event=[start,last],restore=restore,
                          primary_post=[restore,restore+p['primary_post_steps']-1],
                          extended_post=[restore,restore+p['extended_post_steps']-1]))
n=len(p['simulation_seeds']); m=len(p['mask_rules']); t=len(p['duration_steps'])
assert p['simulation_seeds']==list(range(1,21))
rows=[]; tests=[]
for seed in p['simulation_seeds']:
    ref=f'S{seed:03d}-fresh'
    common=dict(seed=seed,checkpoint_id=f'S{seed:03d}-step6000',event_checkpoint_id=f'S{seed:03d}-step6200',reference_id=ref,branch_steps=p['branch_steps'])
    rows.append(dict(run_id=ref,**common,mask_rule='not_applicable',duration_steps=0,mode='fresh'))
    for mask in p['mask_rules']:
        for d in p['duration_steps']:
            for mode in p['modes'][1:]:
                rows.append(dict(run_id=f'S{seed:03d}-{mask}-d{d}-{mode}',**common,mask_rule=mask,duration_steps=d,mode=mode))
    for mode in p['modes']:
        tests.append(dict(run_id=f'TEST-S{seed:03d}-d0-{mode}',**common,mask_rule='not_applicable',duration_steps=0,mode=mode))
assert len(rows)==560 and len(tests)==80
assert len({r['run_id'] for r in rows})==560
assert all(sum(r['seed']==s for r in rows)==28 for s in p['simulation_seeds'])
for name,records in [('RUN_MANIFEST.csv',rows),('T0_TEST_MANIFEST.csv',tests)]:
    buf=io.StringIO(newline=''); writer=csv.DictWriter(buf,fieldnames=list(records[0]));writer.writeheader();writer.writerows(records)
    content=buf.getvalue().encode('utf-8-sig'); path=root/name
    if path.exists(): assert path.read_bytes()==content,'Refusing changed run manifest'
    else: path.write_bytes(content)
a=p['baseline_total_weight']; k=p['k']; h=p['hidden_count']
report=dict(status='PASS design arithmetic only; no simulator test executed',
            proposal_not_execution=n>0,per_seed_trajectories=1+m*t*3,
            pilot_trajectories=n*(1+m*t*3),minimum_t0_test_trajectories=n*4,
            candidate_hidden_subsets_per_agent=math.comb(k,h),time_intervals=intervals,
            expected_amplitudes=dict(fresh=a,drop=a*(k-h)/k,drop_normalized=a,hold=a),
            source_hashes={x.name:hashlib.sha256(x.read_bytes()).hexdigest() for x in root.iterdir() if x.suffix in ('.md','.json','.py','.csv') and not x.name.startswith('design_check')})
(root/'design_check_v4.json').write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
print(json.dumps({x:y for x,y in report.items() if x!='source_hashes'},ensure_ascii=False,indent=2))

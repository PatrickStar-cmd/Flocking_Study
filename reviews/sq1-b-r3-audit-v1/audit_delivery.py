"""Independent provenance, manifests and test-log checks after LFS hydration."""
import csv, hashlib, json, subprocess
from pathlib import Path
r=Path(__file__).resolve().parent;s=r/'source';p=s/'results/sq1-fourbranch-protocol12-seed1-r3-results-v4';repo=r.parents[1]
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def git(*args):return subprocess.check_output(['git',*args],cwd=repo)
manifest=json.loads((p/'source_manifest.json').read_text(encoding='utf-8'))
sources=[dict(path=x['relativePath'],match=sha(s/x['relativePath'])==x['SHA256'].lower()) for x in manifest['sourceFiles']]
snapshots=[]
for path in (s/'sources/A').glob('*'):
    snapshots.append(dict(path=path.name,type='A',match=path.read_bytes()==git('show','992e71a:effective-social-input/'+path.name)))
for path in (s/'sources/protocol').glob('*'):
    snapshots.append(dict(path=path.name,type='protocol',match=path.read_bytes()==git('show','1e11d64:reviews/sq1-c-occlusion-protocol-v1/'+path.name)))
legacy=s/'sources/legacy/simulateCollectiveFlockingNeuralField.m'
snapshots.append(dict(path=legacy.name,type='original',match=legacy.read_bytes()==git('show','origin/main:allocentric-flocking-original/simulateCollectiveFlockingNeuralField.m')))
tests=json.loads((p/'tests.json').read_text())
logs=[dict(test_id=t['test_id'],status=t['status'],log_exists=(p/t['logPath']).is_file(),reports_pass=(t['test_id']+' PASS') in (p/t['logPath']).read_text(errors='replace')) for t in tests]
t0=list(csv.DictReader((p/'T0_TEST_MANIFEST.csv').open(encoding='utf-8-sig')))
runs=list(csv.DictReader((p/'RUN_MANIFEST.csv').open(encoding='utf-8-sig')))
out=dict(source_checks=sources,snapshot_checks=snapshots,test_log_checks=logs,T0=dict(planned=len(t0),completed=sum(t['status']=='COMPLETED' for t in t0),not_run=sum(t['status']=='NOT_RUN' for t in t0)),runs=runs,protocol_hash_matches=sha(s/'sources/protocol/PROTOCOL.md')==manifest['protocolSHA'].lower(),config_hash_matches=sha(s/'defaultFourBranchConfig.m')==manifest['configSHA'].lower())
(r/'delivery_audit.json').write_text(json.dumps(out,ensure_ascii=False,indent=2),encoding='utf-8')
assert all(x['match'] for x in sources+snapshots)
assert all(x['log_exists'] and x['reports_pass'] for x in logs)
assert len(runs)==4 and len({x['run_id'] for x in runs})==4
print(json.dumps(dict(source_hash_matches=len(sources),frozen_snapshot_matches=len(snapshots),PASS_logs=len(logs),T0=out['T0']),ensure_ascii=False))

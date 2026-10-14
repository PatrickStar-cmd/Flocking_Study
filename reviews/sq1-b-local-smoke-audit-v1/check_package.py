"""Read-only inventory/hash verification of B's local package."""
import csv, hashlib, json
from pathlib import Path
root=Path(__file__).resolve().parent
src=Path('D:/Work/Study/UESTC资料存档/课程/通关项目/paper/experiment/sq1_fourbranch_pilot_20260929')
pkg=src/'results/sq1-fourbranch-pilot-v1-seed1-retry'
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(8*1024*1024),b''): h.update(b)
    return h.hexdigest()
checks=[]
for row in csv.DictReader((pkg/'SHA256SUMS.csv').open(encoding='utf-8-sig')):
    p=pkg/row['File']; actual=sha(p)
    checks.append(dict(file=row['File'],size=p.stat().st_size,sha256=actual,match=actual==row['SHA256'].lower(),size_match=p.stat().st_size==int(row['Bytes'])))
sources=[]
for row in json.loads((src/'source_manifest.json').read_text(encoding='utf-8-sig')):
    actual=sha(src/row['file']);sources.append(dict(file=row['file'],sha256=actual,match=actual==row['sha256'].lower()))
report=dict(source=str(src),package=str(pkg),file_checks=checks,source_checks=sources,
    run_manifest=list(csv.DictReader((pkg/'RUN_MANIFEST.csv').open(encoding='utf-8-sig'))),
    t0_manifest=list(csv.DictReader((pkg/'T0_TEST_MANIFEST.csv').open(encoding='utf-8-sig'))),
    tests=json.loads((pkg/'tests.json').read_text(encoding='utf-8-sig')))
(root/'package_audit.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(dict(files=len(checks),file_hash_matches=sum(x['match'] for x in checks),source_hash_matches=sum(x['match'] for x in sources),sources=len(sources),run_rows=len(report['run_manifest']),t0_rows=len(report['t0_manifest']))))

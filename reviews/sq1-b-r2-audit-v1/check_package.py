"""Inventory B r2 without changing its files."""
import csv, hashlib, json
from pathlib import Path
root=Path(__file__).resolve().parent
src=Path('D:/Work/Study/UESTC资料存档/课程/通关项目/WHQ/sq1-fourbranch-protocol12-seed1-r2')
pkg=src/'results/sq1-fourbranch-protocol12-seed1-r2-results-v3'
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
    return h.hexdigest()
rows=list(csv.DictReader((pkg/'SHA256SUMS.csv').open(encoding='utf-8-sig')))
checks=[dict(file=r['File'],sha256=sha(pkg/r['File']),match=sha(pkg/r['File'])==r['SHA256'].lower(),size_match=(pkg/r['File']).stat().st_size==int(r['Bytes'])) for r in rows]
tests=json.loads((pkg/'tests.json').read_text(encoding='utf-8-sig'))
sources=[dict(file=p.name,sha256=sha(p)) for p in sorted(src.glob('*.m'))]
report=dict(source=str(src),package=str(pkg),file_checks=checks,source_hashes=sources,
    test_count=len(tests),declared_pass=sum(t['status']=='PASS' for t in tests),
    local_test_logs=list(map(str,pkg.glob('logs/*.log'))),
    run_manifest=list(csv.DictReader((pkg/'RUN_MANIFEST.csv').open(encoding='utf-8-sig'))),
    t0_manifest=list(csv.DictReader((pkg/'T0_TEST_MANIFEST.csv').open(encoding='utf-8-sig'))),
    unhashed_files=[p.name for p in pkg.iterdir() if p.is_file() and p.name not in {r['File'] for r in rows} and p.name!='SHA256SUMS.csv'])
(root/'package_audit.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(dict(hash_matches=sum(r['match'] for r in checks),hash_count=len(checks),source_files=len(sources),test_logs=len(report['local_test_logs']),runs=len(report['run_manifest']),T0=len(report['t0_manifest']))))

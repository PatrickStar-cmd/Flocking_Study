"""Pinned, read-only source extraction and CSV/ZIP integrity audit."""
import csv, hashlib, io, json, subprocess, zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parents[1]
NEW = '0036487818d891c23d357f7af49030f323f4c9af'
OLD = '63d56b51cb74bca509b0fb7dc29d1d5997d3a1ae'

def git(*args):
    return subprocess.check_output(['git', *args], cwd=REPO)

def save(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        assert path.read_bytes() == data, f'Will not overwrite: {path}'
    else:
        path.write_bytes(data)

paths = git('diff', '--name-only', '-z', OLD, NEW).decode().split('\0')
manifest, members = [], []
for path in filter(None, paths):
    data = git('show', f'{NEW}:{path}')
    save(ROOT / 'source' / path, data)
    manifest.append(dict(path=path, size=len(data), sha256=hashlib.sha256(data).hexdigest()))
    if path.endswith('.zip'):
        with zipfile.ZipFile(io.BytesIO(data)) as z:
            assert z.testzip() is None
            for info in z.infolist():
                if info.is_dir():
                    continue
                body = z.read(info)
                assert info.filename.endswith('.mat'), info.filename
                name = Path(info.filename).name
                assert name not in [m['file'] for m in members], name
                save(ROOT / 'raw' / name, body)
                members.append(dict(archive=path, member=info.filename, file=name,
                                    size=len(body), sha256=hashlib.sha256(body).hexdigest()))

old_bytes = git('show', f'{OLD}:effective-social-input/run_B_E2/run_B_E2_results/B_E2_runs.csv')
new_bytes = (ROOT/'source/run_B_1_03_summary_results/B_E2_runs.csv').read_bytes()
old = list(csv.DictReader(io.StringIO(old_bytes.decode('utf-8-sig'))))
new = list(csv.DictReader(io.StringIO(new_bytes.decode('utf-8-sig'))))
key = lambda r: tuple(r[c] for c in ('experiment','strengthMode','selectionPolicy','k','seed'))
od, nd = {key(r):r for r in old}, {key(r):r for r in new}
assert len(old)==len(new)==len(od)==len(nd)==480
assert od.keys()==nd.keys()
metrics = [c for c in new[0] if c not in ('experiment','strengthMode','selectionPolicy','k','seed')]
errors = {c:max(abs(float(nd[k][c])-float(od[k][c])) for k in nd) for c in metrics}
expected = {f"{r['experiment']}_{r['strengthMode']}_{r['selectionPolicy']}_k{r['k']}_seed{r['seed']}.mat" for r in new}
assert {m['file'] for m in members}==expected
canonical = set()
for r in new:
    mode, policy, k = r['strengthMode'], r['selectionPolicy'], int(r['k'])
    if policy=='angle-dispersed': policy='balanced'
    if k==19: mode,policy='full-equivalent','all'
    canonical.add((mode,policy,k,int(r['seed'])))
report = dict(old_commit=OLD,new_commit=NEW,changed_files=len(manifest),
              zip_count=sum(m['path'].endswith('.zip') for m in manifest),
              archive_bytes=sum(m['size'] for m in manifest if m['path'].endswith('.zip')),
              raw_count=len(members),csv_rows=len(new),canonical_run_count=len(canonical),
              csv_byte_identical=old_bytes==new_bytes,csv_max_abs_differences=errors,
              seed_values=sorted({int(r['seed']) for r in new}),
              source_manifest=manifest,raw_manifest=members)
save(ROOT/'source_audit.json', (json.dumps(report,ensure_ascii=False,indent=2)+'\n').encode())
print(json.dumps({k:v for k,v in report.items() if k not in ('source_manifest','raw_manifest')},indent=2))

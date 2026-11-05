"""Read-only complete-package verification using the Python standard library."""
from pathlib import Path
import csv
import hashlib
import json

ROOT=Path(__file__).resolve().parent


def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
    return h.hexdigest()


def main():
    with (ROOT/'DELIVERY_MANIFEST.csv').open(encoding='utf-8-sig',newline='') as f:rows=list(csv.DictReader(f))
    actual={p.relative_to(ROOT).as_posix() for p in ROOT.rglob('*') if p.is_file()
            and '__pycache__' not in p.relative_to(ROOT).parts and p!=ROOT/'DELIVERY_MANIFEST.csv'}
    assert actual=={r['path'] for r in rows},'File-set mismatch'
    for r in rows:
        p=ROOT/r['path'];assert p.stat().st_size==int(r['bytes']) and sha(p)==r['sha256'],r['path']
    for manifest,snapshot in [('frozen_source_hashes.json','source_snapshot'),('analysis_source_hashes.json','analysis_source_snapshot')]:
        hashes=json.loads((ROOT/'results'/manifest).read_text(encoding='utf-8'))
        for name,digest in hashes.items():assert sha(ROOT/name)==sha(ROOT/'results'/snapshot/name)==digest,name
    audit=json.loads((ROOT/'results/checks/INDEPENDENT_VERIFICATION.json').read_text(encoding='utf-8'))
    assert audit['status']=='PASS' and audit['unique_trajectories']==65
    visual=json.loads((ROOT/'results/checks/VISUAL_QA.json').read_text(encoding='utf-8'))
    assert visual['status']=='PASS'
    for name,digest in visual['figure_sha256'].items():assert sha(ROOT/'results/figures'/name)==digest,name
    status=json.loads((ROOT/'results/STATUS.json').read_text(encoding='utf-8'))
    assert status['stage']=='COMPLETE' and status['completed']==status['planned']==65
    print('DELIVERY PASS:',len(rows),'files;',sum(int(r['bytes']) for r in rows),'bytes',flush=True)


if __name__=='__main__':main()

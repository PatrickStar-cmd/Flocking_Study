"""Read-only frozen-package check; does not run simulations or alter the snapshot."""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import stat
import subprocess
import sys
import zipfile

ROOT=Path(__file__).resolve().parent

def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
    return h.hexdigest()

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--zip',type=Path,help='Optionally verify the adjacent ZIP and its .sha256 file')
    args=parser.parse_args()
    with (ROOT/'FROZEN_MANIFEST.csv').open(encoding='utf-8-sig',newline='') as f:rows=list(csv.DictReader(f))
    assert len({r['path'] for r in rows})==len(rows)
    actual={p.relative_to(ROOT).as_posix() for p in ROOT.rglob('*') if p.is_file()
            and '__pycache__' not in p.relative_to(ROOT).parts and p!=ROOT/'FROZEN_MANIFEST.csv'}
    assert actual=={r['path'] for r in rows},'Frozen file-set changed'
    for r in rows:
        p=ROOT/r['path']
        assert p.stat().st_size==int(r['bytes']) and sha(p)==r['sha256'],r['path']
    meta=json.loads((ROOT/'冻结索引.json').read_text(encoding='utf-8'))
    snapshot=ROOT/'实验快照'
    original=meta['source_file_sha256']
    snapshot_files={p.relative_to(snapshot).as_posix():sha(p) for p in snapshot.rglob('*')
                    if p.is_file() and '__pycache__' not in p.relative_to(snapshot).parts}
    assert snapshot_files==original,'Snapshot differs from the verified source at freeze time'
    assert sha(snapshot/'experiment_config.json')==meta['config_sha256']
    assert sha(snapshot/'results/analysis/PRIMARY_RESULTS.json')==meta['primary_results_sha256']
    subprocess.run([sys.executable,str(snapshot/'verify_delivery.py')],check=True)
    readonly=sum(not bool(p.stat().st_mode & stat.S_IWRITE) for p in ROOT.rglob('*') if p.is_file())
    zip_checked=False
    if args.zip is not None:
        archive=args.zip.resolve()
        expected=archive.with_name(archive.name+'.sha256').read_text(encoding='utf-8').split()[0]
        assert sha(archive)==expected,'ZIP SHA256 mismatch'
        paths=[p for p in ROOT.rglob('*') if p.is_file() and '__pycache__' not in p.relative_to(ROOT).parts]
        with zipfile.ZipFile(archive) as z:
            assert z.testzip() is None,'ZIP CRC check failed'
            expected_names={ROOT.name+'/'+p.relative_to(ROOT).as_posix() for p in paths}
            assert set(z.namelist())==expected_names and len(z.namelist())==len(expected_names)
            for p in paths:
                name=ROOT.name+'/'+p.relative_to(ROOT).as_posix()
                h=hashlib.sha256()
                with z.open(name) as f:
                    for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
                assert h.hexdigest()==sha(p),name
        zip_checked=True
    print('FROZEN PASS:',meta['freeze_id'],';',len(snapshot_files),'source files identical;',len(rows)+1,
          'frozen files;',readonly,'read-only files; ZIP checked:',zip_checked)

if __name__=='__main__':main()

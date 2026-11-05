"""Atomic output, SHA256 and frozen scientific-source snapshots."""
from pathlib import Path
import csv
import hashlib
import json
import os

ROOT=Path(__file__).resolve().parent
RESULTS=ROOT/'results'
ENGINE_FILES=('bootstrap.py','common.py','design.py','model.py','metrics.py','kf_reference_core.py',
              'run_experiment.py','experiment_config.json','EXPERIMENT_PROTOCOL.md','SOURCE_PROVENANCE.json','requirements.lock.txt')


def digest(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(1024*1024),b''): h.update(block)
    return h.hexdigest()


def read_json(path):
    return json.loads(Path(path).read_text(encoding='utf-8-sig'))


def write_json(path,value):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    temporary=path.with_name(path.name+'.partial')
    temporary.write_text(json.dumps(value,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    os.replace(temporary,path)


def write_csv(path,rows):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    assert rows, 'Empty result table'
    temporary=path.with_name(path.name+'.partial')
    with temporary.open('w',encoding='utf-8-sig',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    os.replace(temporary,path)


def read_csv(path):
    with Path(path).open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))


def freeze_engine():
    hashes={name:digest(ROOT/name) for name in ENGINE_FILES}
    path=RESULTS/'frozen_source_hashes.json'
    if path.exists():
        assert read_json(path)==hashes,'Scientific sources changed; use a new version.'
    else:
        directory=RESULTS/'source_snapshot';directory.mkdir(parents=True,exist_ok=True)
        for name in hashes:(directory/name).write_bytes((ROOT/name).read_bytes())
        write_json(path,hashes)
    return hashes

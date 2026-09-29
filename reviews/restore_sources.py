"""Restore pinned source snapshots without overwriting any differing file.

Run from any directory after fetching the WHQ_work branch. Requires Python 3.
Uses only git and Python's standard library, not the simulation runners.
"""
import hashlib
import io
import json
import subprocess
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
count = 0

def restore(path, body, expected):
    global count
    assert hashlib.sha256(body).hexdigest() == expected, f'Hash mismatch: {path}'
    path = path.resolve()
    assert path.is_relative_to(ROOT.resolve())
    if path.exists():
        assert path.read_bytes() == body, f'Existing file differs; preserved: {path}'
    else:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(body)
    count += 1

def blob(commit, path):
    return subprocess.check_output(['git', 'show', f'{commit}:{path}'], cwd=REPO)

old_root = ROOT/'sq1-b-e2-63d56b5'
old = json.loads((old_root/'source_manifest.json').read_text(encoding='utf-8'))
for entry in old['files']:
    restore(old_root/'source'/entry['path'], blob(old['commit'], entry['path']), entry['sha256'])

new_root = ROOT/'sq1-whq-update-0036487'
new = json.loads((new_root/'source_audit.json').read_text(encoding='utf-8'))
raw = {(r['archive'], r['member']): r for r in new['raw_manifest']}
restored_raw = set()
for entry in new['source_manifest']:
    body = blob(new['new_commit'], entry['path'])
    restore(new_root/'source'/entry['path'], body, entry['sha256'])
    if entry['path'].endswith('.zip'):
        with zipfile.ZipFile(io.BytesIO(body)) as z:
            assert z.testzip() is None
            for member in z.infolist():
                if member.is_dir(): continue
                key = (entry['path'], member.filename)
                r = raw[key]
                restore(new_root/'raw'/r['file'], z.read(member), r['sha256'])
                assert key not in restored_raw
                restored_raw.add(key)
assert len(restored_raw) == len(raw) == 480
print(f'PASS: {count} pinned source/raw files restored or verified, including 480 formal MAT files.')

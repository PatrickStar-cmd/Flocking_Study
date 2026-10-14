"""Verify staged publication bytes, counts, links, and recorded output hashes."""
import csv
import hashlib
import json
import re
import subprocess
from pathlib import Path

ROOT=Path(__file__).resolve().parent
REPO=ROOT.parent
def git(*args): return subprocess.check_output(['git',*args],cwd=REPO)
paths=[p for p in git('diff','--cached','--name-only','-z').decode().split('\0') if p]
assert paths
assert all(p.startswith('reviews/') or p=='科学问题一_三人分工与研究计划.md' for p in paths)
assert not any('/source/' in p or '/raw/' in p for p in paths)
for p in paths:
    if p.startswith('reviews/'):
        assert git('show',f':{p}')==(REPO/p).read_bytes(),p
assert sum(p.endswith('.png') for p in paths)==14
assert sum(p.endswith('.svg') for p in paths)==14
for filename,n in [('sq1-c-raw-index-v1/CANONICAL_RAW_INDEX.csv',240),
                   ('sq1-c-time-analysis-v1/window_effects.csv',132),
                   ('sq1-c-time-analysis-v1/paired_windows.csv',1980),
                   ('sq1-c-paired-analysis-v1/paired_effects.csv',33)]:
    with (ROOT/filename).open(encoding='utf-8-sig',newline='') as f:
        assert len(list(csv.DictReader(f)))==n
for file,key in [('sq1-c-time-analysis-v1/analysis_audit.json','files'),
                 ('sq1-c-raw-index-v1/audit.json','outputs'),
                 ('sq1-c-paired-analysis-v1/analysis_audit.json','outputs')]:
    audit_path=ROOT/file
    for name,h in json.loads(audit_path.read_text(encoding='utf-8'))[key].items():
        assert hashlib.sha256((audit_path.parent/name).read_bytes()).hexdigest()==h,name
links=0
for file in ['README.md','sq1-c-time-analysis-v1/README.md','sq1-c-time-analysis-v1/REPORT.md',
             'sq1-c-visualization-v1/README.md']:
    p=ROOT/file
    for target in re.findall(r'\]\(([^)]+)\)',p.read_text(encoding='utf-8')):
        if '://' in target or target.startswith('#'):continue
        dest=(p.parent/target.split('#')[0]).resolve()
        assert dest.is_file(),(file,target)
        assert dest.relative_to(REPO).as_posix() in paths,(file,target)
        links+=1
print(json.dumps(dict(status='PASS',staged_files=len(paths),png=14,svg=14,
                      verified_reader_links=links,staged_bytes=sum((REPO/p).stat().st_size for p in paths),
                      staged_bytes_match=True,analysis_hashes_match=True),indent=2))

"""Verify uploaded files, frozen MATLAB source, and summary completeness.
This does not repeat raw-MAT scientific validation: MAT files are local only.
"""
from pathlib import Path
import csv, hashlib, json

ROOT=Path(__file__).resolve().parent
def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()
def rows(path):
    with path.open(encoding='utf-8-sig',newline='') as f:
        return list(csv.DictReader(f))
def main():
    manifest=rows(ROOT/'UPLOAD_SHA256SUMS.csv')
    for row in manifest:
        p=ROOT/row['path']
        assert p.is_file() and p.stat().st_size==int(row['bytes']),row['path']
        assert sha(p)==row['sha256'],row['path']
    result=ROOT/'results'/'pilot_seed1to5_v1'
    for row in json.loads((result/'source_manifest.json').read_text(encoding='utf-8')):
        assert sha(ROOT/row['file'])==row['sha256'],row['file']
    runs=rows(result/'RUN_MANIFEST.csv')
    assert len(runs)==155 and len({r['run_id'] for r in runs})==155
    assert {int(r['seed']) for r in runs}=={1,2,3,4,5}
    for filename,count in [('closed_loop_summary.csv',700),('same_state_summary.csv',175),('prediction_diagnostics_summary.csv',105)]:
        assert len(rows(result/filename))==count,filename
    raw=rows(ROOT/'RAW_DATA_MANIFEST.csv')
    assert len(raw)==270 and all(x['path'].endswith('.mat') for x in raw)
    print(json.dumps(dict(status='PASS',uploadedFilesVerified=len(manifest),recordedTrajectoryCount=155,
        excludedMATFiles=len(raw),rawMATValidation='not repeated; local-only data'),indent=2))
if __name__=='__main__': main()

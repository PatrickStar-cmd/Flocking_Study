"""Freeze final delivery hashes after all simulations and analyses complete."""
from pathlib import Path
import json,csv,hashlib,sys
ROOT=Path(__file__).resolve().parent
def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda:f.read(4*1024**2),b''): h.update(b)
    return h.hexdigest()
def main():
    result=ROOT/'results'/'pilot_seed1to5_v1'
    audit=json.loads((result/'audit.json').read_text(encoding='utf-8'))
    validation=json.loads((result/'analysis'/'validation.json').read_text(encoding='utf-8'))
    assert audit['status']=='COMPLETE' and validation['status']=='PASS'
    manifest=json.loads((result/'source_manifest.json').read_text(encoding='utf-8'))
    assert all(sha(ROOT/x['file'])==x['sha256'] for x in manifest),'MATLAB simulation source changed since run began'
    source_check={'status':'PASS','frozenMATLABFiles':len(manifest),'source':'results/pilot_seed1to5_v1/source_manifest.json'}
    (result/'analysis'/'source_integrity.json').write_text(json.dumps(source_check,indent=2),encoding='utf-8')
    paths=[p for p in ROOT.rglob('*') if p.is_file() and '__pycache__' not in p.parts and p.name not in ['SHA256SUMS.csv','delivery_verification.log']]
    with (ROOT/'SHA256SUMS.csv').open('w',encoding='utf-8-sig',newline='') as f:
        w=csv.DictWriter(f,fieldnames=['path','bytes','sha256']); w.writeheader()
        for p in sorted(paths): w.writerow(dict(path=p.relative_to(ROOT).as_posix(),bytes=p.stat().st_size,sha256=sha(p)))
    print(json.dumps(dict(status='PASS',files=len(paths),totalBytes=sum(p.stat().st_size for p in paths),sourceCheck=source_check),indent=2))
if __name__=='__main__': main()

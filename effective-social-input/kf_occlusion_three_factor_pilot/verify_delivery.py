"""交付包结构和哈希检查；不修改结果。"""
from pathlib import Path
import hashlib, json, sys, csv
def sha(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''): h.update(b)
    return h.hexdigest()
def main(root):
    root=Path(root); result=root/'results'/'pilot_seed1to5_v1'
    assert (result/'audit.json').exists(), 'pilot audit missing'
    audit=json.loads((result/'audit.json').read_text(encoding='utf-8'))
    man=json.loads((result/'RUN_MANIFEST.json').read_text(encoding='utf-8'))
    assert audit['status']=='COMPLETE' and len(man)==155
    assert len(list((result/'conditions').glob('*/KF-mean.mat')))==35
    assert len(list((result/'conditions').glob('*/checkpoint6200.mat')))==35
    assert (result/'closed_loop_summary.csv').exists()
    rows=list(csv.DictReader((root/'SHA256SUMS.csv').open(encoding='utf-8-sig')))
    for row in rows:
        p=root/row['path']; assert p.is_file(),row['path']
        assert p.stat().st_size==int(row['bytes']) and sha(p)==row['sha256'],row['path']
    print(json.dumps({'status':'PASS','actualTrajectoryCount':len(man),'conditionFolders':35,'filesVerified':len(rows)},ensure_ascii=False,indent=2))
if __name__=='__main__': main(Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).resolve().parent)

"""Fetch public GitHub LFS media with SHA256 verification; no credentials."""
import concurrent.futures, hashlib, json, urllib.request
from pathlib import Path
root=Path(__file__).resolve().parent
pending=json.loads((root/'lfs_pending.json').read_text())
def fetch(item):
    oid=item['oid'];obj=root.parents[1]/'.git/lfs/objects'/oid[:2]/oid[2:4]/oid
    if obj.exists():return dict(oid=oid,status='already_present')
    url='https://media.githubusercontent.com/media/PatrickStar-cmd/Flocking_Study/54fb85d91151e0b28e8e570318b4a9937d6a0425/'+item['path']
    dst=root/'downloads'/oid;dst.parent.mkdir(exist_ok=True);h=hashlib.sha256();size=0
    try:
        with urllib.request.urlopen(url,timeout=30) as res,dst.open('wb') as out:
            while True:
                b=res.read(4*1024*1024)
                if not b:break
                out.write(b);h.update(b);size+=len(b)
        assert h.hexdigest()==oid and size==item['size']
        obj.parent.mkdir(parents=True,exist_ok=True)
        if not obj.exists():dst.rename(obj)
        return dict(oid=oid,status='verified',size=size)
    except Exception as e:return dict(oid=oid,status='error',error=str(e))
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
    records=[]
    for record in pool.map(fetch,pending):records.append(record);print(json.dumps(record),flush=True)
(root/'media_download.json').write_text(json.dumps(records,indent=2))

"""Extract a pinned Git tree and verified LFS objects without checking out B's branch."""
import csv, hashlib, json, subprocess
from pathlib import Path
root=Path(__file__).resolve().parent; repo=root.parents[1]
commit='54fb85d91151e0b28e8e570318b4a9937d6a0425'
prefix='effective-social-input/results/sq1-fourbranch-protocol12-seed1-r3/'
def git(*args):return subprocess.check_output(['git',*args],cwd=repo)
paths=git('ls-tree','-r','--name-only','-z',commit,'--',prefix).decode().strip('\0').split('\0')
records=[];pending=[]
for name in paths:
    raw=git('show',f'{commit}:{name}'); content=raw; lfs=None
    if raw.startswith(b'version https://git-lfs.github.com/spec/v1\n'):
        lines=raw.decode().splitlines();oid=lines[1].split('sha256:')[1];size=int(lines[2].split()[1]);lfs=dict(oid=oid,size=size)
        obj=repo/'.git/lfs/objects'/oid[:2]/oid[2:4]/oid
        if not obj.exists():pending.append(dict(path=name,oid=oid,size=size));continue
        content=obj.read_bytes();assert len(content)==size and hashlib.sha256(content).hexdigest()==oid
    dest=root/'source'/name[len(prefix):];dest.parent.mkdir(parents=True,exist_ok=True)
    if dest.exists():assert dest.read_bytes()==content
    else:dest.write_bytes(content)
    records.append(dict(path=name[len(prefix):],sha256=hashlib.sha256(content).hexdigest(),size=len(content),lfs=lfs))
if pending:
    (root/'lfs_pending.json').write_text(json.dumps(pending,indent=2),encoding='utf-8')
    print(json.dumps(dict(extracted=len(records),pending_lfs=len(pending))));raise SystemExit(0)
manifest=list(csv.DictReader((root/'source/SHA256SUMS.csv').open(encoding='utf-8-sig')))
checks=[]
for row in manifest:
    key=next(k for k in row if k.lower() in ('file','path','relativepath','relative_path'));p=root/'source'/row[key].replace('\\','/')
    hashkey=next(k for k in row if k.lower()=='sha256');sizekey=next(k for k in row if k.lower() in ('bytes','size'))
    actual=hashlib.sha256(p.read_bytes()).hexdigest();checks.append(dict(file=row[key],hash_match=actual==row[hashkey].lower(),size_match=p.stat().st_size==int(row[sizekey])))
out=dict(commit=commit,files=records,delivery_hash_checks=checks)
(root/'source_audit.json').write_text(json.dumps(out,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(dict(files=len(records),LFS=sum(r['lfs'] is not None for r in records),delivery_entries=len(checks),hash_matches=sum(c['hash_match'] for c in checks))))

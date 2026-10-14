"""Extract the pinned A handoff and exact dependency for read-only testing."""
import hashlib,json,subprocess
from pathlib import Path
root=Path(__file__).resolve().parent; repo=root.parents[1]
commit='992e71aba3b26896335b52a3182cf01cfd542a9d'
def git(*args):return subprocess.check_output(['git',*args],cwd=repo)
paths=git('diff','--name-only','-z','origin/main',commit,'--','effective-social-input').decode().strip('\0').split('\0')
paths+=['effective-social-input/buildRingConnectivity.m']
manifest=[]
for p in paths:
    data=git('show',f'{commit}:{p}'); dst=root/'source'/p
    dst.parent.mkdir(parents=True,exist_ok=True)
    if dst.exists():assert dst.read_bytes()==data
    else:dst.write_bytes(data)
    manifest.append(dict(path=p,sha256=hashlib.sha256(data).hexdigest(),bytes=len(data)))
audit=json.loads((root/'source/effective-social-input/results/sq1-local-v1/audit.json').read_text())
hash_checks={}
for name,h in audit['sourceSHA256'].items():
    data=(root/'source/effective-social-input'/name).read_bytes()
    actual=hashlib.sha256(data).hexdigest()
    hash_checks[name]=dict(exact=actual==h.lower(),lf_normalized=hashlib.sha256(data.replace(b'\r\n',b'\n')).hexdigest()==h.lower())
out=dict(commit=commit,files=manifest,published_code_hash_checks=hash_checks,
         published_audit_commit=audit['gitCommit'],published_audit_branch=audit['branch'])
(root/'source_manifest.json').write_text(json.dumps(out,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
print(json.dumps({k:v for k,v in out.items() if k!='files'},indent=2))

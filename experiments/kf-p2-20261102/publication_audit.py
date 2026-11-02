"""Final immutable artifact manifest and protected-file check before staging."""
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys


def sha(path):
    return hashlib.file_digest(path.open('rb'),'sha256').hexdigest()


def run(root):
    root=root.resolve()
    here=Path(__file__).resolve().parent
    repo=here.parent.parent
    path=root/'review/publication_audit.json'
    assert not path.exists()
    old=json.loads((root/'source_hashes.json').read_text(encoding='utf-8'))
    for entry in old['sources']:
        assert sha(root/entry['path']) == entry['sha256']
        if entry['path'].startswith('sources/p2/'):
            assert sha(here/Path(entry['path']).name) == entry['sha256']
    diff=subprocess.check_output(['git','diff','--name-only','c0ab184','--','.',
        ':!experiments/kf-p2-20261102',':!effective-social-input/results/kf-p2-v1-20261102-reviewed'],text=True)
    assert not diff.strip(),diff
    for name in ['input_audit.json','matlab_parity.json','statistical_audit.json']:
        assert json.loads((root/'review'/name).read_text(encoding='utf-8'))['status']=='PASS'
    reviewer=root/'sources/publication'; assert not reviewer.exists();reviewer.mkdir()
    for name in ['supplement.py','publication_audit.py']:
        shutil.copyfile(here/name,reviewer/name)
    entries=[]
    for parent in [here,root]:
        for f in sorted(parent.rglob('*')):
            if not f.is_file() or '__pycache__' in f.parts or f.suffix=='.pyc':continue
            entries.append({'path':f.relative_to(repo).as_posix(),'bytes':f.stat().st_size,'sha256':sha(f)})
    binaries=[e for e in entries if Path(e['path']).suffix in {'.mat','.npz'}]
    assert len(binaries)==41
    record={'status':'PASS','source_commit':old['base_commit'],
        'branch':old['branch'],'target_branch':'Hao_L','evidence':'exploratory-development',
        'protected_tracked_diff_before_staging':diff,
        'raw_binary_files':len(binaries),'raw_binary_bytes':sum(e['bytes'] for e in binaries),
        'independent_seeds':5,'open_loop_conditions':30,'main_figures':7,'reviewer_figures':3,
        'artifacts_sha256':entries,
        'manifest_scope':'All publication files before this manifest; this manifest cannot hash itself. Git commit records final identities.'}
    path.write_text(json.dumps(record,indent=2,ensure_ascii=False),encoding='utf-8')
    print(json.dumps({k:v for k,v in record.items() if k!='artifacts_sha256'},ensure_ascii=True))


if __name__=='__main__':
    run(Path(sys.argv[1]))

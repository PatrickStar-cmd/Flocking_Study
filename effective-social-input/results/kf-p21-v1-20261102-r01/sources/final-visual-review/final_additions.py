"""Additive final visual review audit; all earlier hashed artifacts stay intact."""
import hashlib,json,shutil,sys
from pathlib import Path
from validate_frozen import check


def sha(path):return hashlib.file_digest(path.open('rb'),'sha256').hexdigest()


def main(root):
    root=root.resolve();here=Path(__file__).resolve().parent;repo=here.parent.parent
    target=root/'FINAL_REVIEW_AUDIT.json';assert not target.exists()
    first=root/'FINAL_AUDIT.json';previous=json.loads(first.read_text(encoding='utf-8'))
    frozen=check(root/'development-r01/FROZEN_DECISION.json');assert frozen['q_star'] is None
    for f in previous['artifacts_sha256']:assert sha(repo/f['path'])==f['sha256'],f['path']
    destination=root/'sources/final-visual-review';assert not destination.exists();destination.mkdir()
    for name in ['review_Q2.py','final_additions.py']:shutil.copyfile(here/name,destination/name)
    original={e['path'] for e in previous['artifacts_sha256']}
    current=[];additions=[]
    for parent in [here,root]:
        for p in sorted(parent.rglob('*')):
            if not p.is_file() or '__pycache__' in p.parts:continue
            path=p.relative_to(repo).as_posix();entry={'path':path,'bytes':p.stat().st_size,'sha256':sha(p)}
            current.append(entry)
            if path not in original:additions.append(entry)
    result={'status':'PASS','previous_manifest_sha256':sha(first),'previous_artifacts_all_unchanged':True,
        'frozen_decision_sha256':sha(root/'development-r01/FROZEN_DECISION.json'),
        'only_additive_visual_QA_after_first_manifest':True,'official_entry':'INDEX.md','official_report':'P21_FINAL_REPORT_v2.md',
        'official_Q2':'review/figures/Q2_chronological.png','validation_not_run':True,
        'files_before_this_manifest':len(current),'total_bytes':sum(f['bytes'] for f in current),
        'artifacts_sha256':current,'added_artifacts':additions,
        'manifest_scope':'Complete final package before this manifest; excludes Python cache and cannot hash itself; Git identifies this manifest.'}
    target.write_text(json.dumps(result,indent=2,ensure_ascii=False),encoding='utf-8')
    print(json.dumps({k:v for k,v in result.items() if k not in ['artifacts_sha256','added_artifacts']},ensure_ascii=True,indent=2))


if __name__=='__main__':main(Path(sys.argv[1]))

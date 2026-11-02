"""Final immutable local package audit; never publishes or rewrites a result."""
import hashlib,json,re,shutil,subprocess,sys
from pathlib import Path
from validate_frozen import check


def sha(path):return hashlib.file_digest(path.open('rb'),'sha256').hexdigest()


def main(root):
    root=root.resolve();here=Path(__file__).resolve().parent;repo=here.parent.parent
    target=root/'FINAL_AUDIT.json';assert not target.exists()
    frozen_path=root/'development-r01/FROZEN_DECISION.json';frozen=check(frozen_path)
    assert frozen['status']=='STOP_NO_ADMISSIBLE_Q' and frozen['q_star'] is None
    assert frozen['validation_n']==0 and frozen['validation_seeds']==[]
    assert not list(root.glob('validation*')),'Unexpected validation artifacts'
    for path in ['baseline/audit.json','g0-r02/audit.json','g0-r02/independent_labels.json','g1/audit.json',
                 'engine-tests-r01/audit.json','engine-tests-r01/matlab_audit.json',
                 'development-r01/review/independent_audit.json','development-r01/review/calibration_audit.json']:
        assert json.loads((root/path).read_text(encoding='utf-8'))['status']=='PASS',path
    baseline='81e0080c8ffc00134c64f407521974a32a950767'
    protected=subprocess.check_output(['git','diff','--name-only',baseline,'--','.',
        ':!experiments/kf-p21-20261102',':!effective-social-input/results/kf-p21-v1-20261102-r01'],text=True)
    assert not protected.strip(),protected
    protocol=json.loads((root/'protocol.json').read_text(encoding='utf-8'))
    assert sha(here/'PROTOCOL.md')==protocol['original_protocol_sha256']
    original=Path('C:/Users/Administrator/Downloads/P2.1_几何解耦与KF校准_详细实验计划书.md')
    assert sha(original)==sha(here/'PROTOCOL.md'),'Original provided protocol bytes changed'
    earlier=json.loads((root/'source_hashes.json').read_text(encoding='utf-8'))
    for source in earlier['files']:assert sha(root/source['path'])==source['sha256']
    destination=root/'sources/frozen-final';assert not destination.exists();destination.mkdir()
    for p in sorted(here.iterdir()):
        if p.is_file() and (p.suffix in {'.py','.m','.md'} or p.name in {'.gitattributes','.gitignore'}):
            shutil.copyfile(p,destination/p.name)
    artifacts=[];sensitive=[]
    for parent in [here,root]:
        for p in sorted(parent.rglob('*')):
            if not p.is_file() or '__pycache__' in p.parts:continue
            digest=sha(p)
            artifacts.append({'path':p.relative_to(repo).as_posix(),'bytes':p.stat().st_size,'sha256':digest})
            if p.suffix in {'.py','.m','.md','.json','.log','.txt'}:
                if re.search(rb'(gh[pousr]_[A-Za-z0-9]{25,}|(?<![A-Za-z0-9])sk-[A-Za-z0-9]{35,})',p.read_bytes()):sensitive.append(p.relative_to(repo).as_posix())
    assert not sensitive,sensitive
    lfs=[e for e in artifacts if Path(e['path']).suffix in {'.mat','.npz','.gz'}]
    report={'status':'PASS','completion':'STOPPED_AT_PREREGISTERED_DEVELOPMENT_GATE',
        'branch':subprocess.check_output(['git','branch','--show-current'],text=True).strip(),
        'baseline_commit':baseline,'frozen_source_commit':frozen['source_commit'],
        'frozen_decision_sha256':sha(frozen_path),'original_protocol_sha256':sha(original),
        'old_model_and_P1_P2_tracked_diff':protected,'development_physical_seeds':12,'development_physical_inputs':24,
        'synthetic_conditions':320,'candidate_source_windows':1200,'method_window_conditions':8400,
        'independently_verified_nonempty_RMSE_rows':58555,'independently_verified_calibration_band_rows':18000,
        'validation_seed_count_used':0,'reserved_validation_seeds_untouched':list(range(2026112001,2026112031)),
        'P3_P4_not_run':True,'remote_publication_not_performed':True,'sensitive_token_scan_paths':sensitive,
        'files_before_this_manifest':len(artifacts),'bytes_before_this_manifest':sum(e['bytes'] for e in artifacts),
        'LFS_files':len(lfs),'LFS_bytes':sum(e['bytes'] for e in lfs),'artifacts_sha256':artifacts,
        'manifest_scope':'All local package files before this manifest; excludes Python cache and cannot hash itself. Final Git commit identifies the manifest.'}
    target.write_text(json.dumps(report,indent=2,ensure_ascii=False),encoding='utf-8')
    print(json.dumps({k:v for k,v in report.items() if k not in ['artifacts_sha256','reserved_validation_seeds_untouched']},ensure_ascii=True,indent=2))


if __name__=='__main__':main(Path(sys.argv[1]))

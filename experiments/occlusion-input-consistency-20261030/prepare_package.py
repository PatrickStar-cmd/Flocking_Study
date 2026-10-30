"""Freeze publication files and protected-worktree fingerprints; zero model calls."""
from pathlib import Path
import hashlib,json,subprocess,platform

HERE=Path(__file__).resolve().parent; REPO=HERE.parents[1]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def git(*args):return subprocess.check_output(['git',*args],cwd=REPO,stderr=subprocess.DEVNULL).decode('utf-8').strip()
def files():
    return sorted(p for p in HERE.rglob('*') if p.is_file() and not any(x in p.relative_to(HERE).parts for x in ['results','__pycache__']) and p.suffix!='.log' and p.name not in ['manifest.v1.json','publication.v1.json'])

cfg=json.loads((HERE/'config.v1.json').read_text(encoding='utf-8'))
for name in ['python-preparation.json','matlab-preparation.json']:
    check=json.loads((HERE/'checks'/name).read_text(encoding='utf-8'))
    assert check['status']=='PASS' and check['model_steps_executed']==check['model_trajectories_executed']==0
assert not (HERE/'results').exists()
assert cfg['primary']['seed_count']*len(cfg['fresh_counts'])==10000
assert sum(p['seed_count'] for p in cfg['history_controls'])*len(cfg['fresh_counts'])==4000
paper=REPO/'paper/论文初稿_KF预测补偿.v7.tex'
paper_text=paper.read_text(encoding='utf-8')
assert 'q^*(\\delta)=1' in paper_text and 'S_i^*=A_i=n_i\\bar w_i' in paper_text
protected=[p for p in (REPO/'allocentric-flocking-original').rglob('*') if p.is_file()]
protected.extend([paper,REPO/'paper/论文草稿_KF预测补偿.v7.pdf'])
manifest=dict(status='prepared_not_run',date='2026-10-30',repository='PatrickStar-cmd/Flocking_Study',branch='YCwork',
 base_commit=git('rev-parse','HEAD'),base_tree=git('rev-parse','HEAD^{tree}'),local_main=git('rev-parse','main'),
 python=platform.python_version(),model_steps_executed=0,model_trajectories_executed=0,
 input_definition='paper v7 eq:chosen_empirical_formula, q*=1, S_i*=sum_j w_ij; population mean immediately after masking',
 model_parameters=cfg['model'],prospective_trajectories=cfg['trajectory_counts'],
 terminal_measure='mean of per-time min(GO,LO), last1000 of6000 steps',
 working_tracked_diff_sha256=hashlib.sha256(git('diff','--binary','--','.',':(exclude)experiments/occlusion-input-consistency-20261030').encode('utf-8')).hexdigest(),
 protected_sha256={str(p.relative_to(REPO)).replace('\\','/'):sha(p) for p in protected},
 package_files=[dict(path=str(p.relative_to(REPO)).replace('\\','/'),sha256=sha(p),bytes=p.stat().st_size) for p in files()])
assert manifest['branch']=='YCwork'
(HERE/'manifest.v1.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2)+'\n',encoding='utf-8',newline='\n')
print(json.dumps(dict(status=manifest['status'],files=len(manifest['package_files'])+1,model_steps=0,prospective=manifest['prospective_trajectories']),ensure_ascii=False))

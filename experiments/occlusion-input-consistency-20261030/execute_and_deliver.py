"""Authorized fixed-protocol MATLAB execution followed by v2 analysis."""
import argparse,datetime,hashlib,json,os,shutil,subprocess,sys,traceback
from pathlib import Path
ROOT=Path(__file__).resolve().parent;REPO=ROOT.parents[1]
def main(run_id):
    run=ROOT/'results'/run_id;assert not run.exists()
    job=ROOT/('execution-'+run_id);job.mkdir(exist_ok=False)
    state={'status':'starting','run_id':run_id,'started':datetime.datetime.now().isoformat()}
    def status(value):
        state['status']=value;state['updated']=datetime.datetime.now().isoformat()
        (job/'status.json').write_text(json.dumps(state,ensure_ascii=False,indent=2),encoding='utf-8')
    try:
        snap=job/'sources';snap.mkdir()
        for p in ROOT.glob('*.m'):shutil.copy2(p,snap/p.name)
        for p in [ROOT/'config.v1.json',ROOT/'execute_and_deliver.py',ROOT/'deliver_results.py']:shutil.copy2(p,snap/p.name)
        shutil.copytree(ROOT/'v2',snap/'v2',ignore=shutil.ignore_patterns('__pycache__'))
        state['head']=subprocess.check_output(['git','rev-parse','HEAD'],cwd=REPO,text=True).strip()
        state['origin_main']=subprocess.check_output(['git','rev-parse','origin/main'],cwd=REPO,text=True).strip()
        state['sources']={str(p.relative_to(snap)):hashlib.sha256(p.read_bytes()).hexdigest() for p in snap.rglob('*') if p.is_file()}
        (job/'runner-fix.diff').write_bytes(subprocess.check_output(['git','diff','--','experiments/occlusion-input-consistency-20261030/runOcclusionConsistency.m'],cwd=REPO))
        status('model_running')
        cmd=['D:/Program/MatLab/MATLAB R2024b/bin/matlab.exe','-wait','-singleCompThread','-batch',f"addpath('experiments/occlusion-input-consistency-20261030');runOcclusionConsistency('execute','{run_id}');delete(gcp('nocreate'));exit(0,'force');",'-logfile',str(job/'matlab.log')]
        with (job/'model-process.log').open('w',encoding='utf-8') as log:subprocess.run(cmd,cwd=REPO,stdout=log,stderr=subprocess.STDOUT,check=True)
        status('v2_analysis_running')
        env=os.environ.copy();env['PYTHONUTF8']='1'
        with (job/'analysis.log').open('w',encoding='utf-8') as log:
            subprocess.run([sys.executable,str(ROOT/'v2/analyze_v2.py'),'--run',run_id],cwd=REPO,stdout=log,stderr=subprocess.STDOUT,check=True,env=env)
            status('report_and_additional_figures_running')
            subprocess.run([sys.executable,str(ROOT/'deliver_results.py'),'--run',run_id],cwd=REPO,stdout=log,stderr=subprocess.STDOUT,check=True,env=env)
        state['report']=str(run/'delivery-v1/实验报告.v1.md');status('complete')
    except BaseException:
        state['error']=traceback.format_exc();status('failed');raise
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--run',required=True);main(p.parse_args().run)

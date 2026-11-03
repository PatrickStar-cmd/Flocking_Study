"""Resume saved trajectories after an explicit user restart."""
import argparse,datetime,json,os,re,subprocess,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parent;REPO=ROOT.parents[1]
def main(run_id):
    assert re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_-]*',run_id)
    run=ROOT/'results'/run_id;job=ROOT/('execution-'+run_id)
    state=json.loads((job/'status.json').read_text(encoding='utf-8'));assert state['status']=='paused'
    tag=datetime.datetime.now().strftime('%Y%m%d-%H%M%S')
    def status(value):
        state.update(status=value,updated=datetime.datetime.now().isoformat())
        (job/'status.json').write_text(json.dumps(state,ensure_ascii=False,indent=2),encoding='utf-8')
    try:
        status('model_running')
        cmd=['D:/Program/MatLab/MATLAB R2024b/bin/matlab.exe','-wait','-singleCompThread','-batch',
             f"addpath('experiments/occlusion-input-consistency-20261030');runOcclusionConsistency('resume','{run_id}');delete(gcp('nocreate'));exit(0,'force');",'-logfile',str(job/f'matlab-resume-{tag}.log')]
        with (job/f'resume-process-{tag}.log').open('w',encoding='utf-8') as log:
            subprocess.run(cmd,cwd=REPO,stdout=log,stderr=subprocess.STDOUT,check=True)
        status('v2_analysis_running');env=os.environ.copy();env['PYTHONUTF8']='1'
        with (job/f'analysis-resume-{tag}.log').open('w',encoding='utf-8') as log:
            subprocess.run([sys.executable,str(ROOT/'v2/analyze_v2.py'),'--run',run_id],cwd=REPO,stdout=log,stderr=subprocess.STDOUT,check=True,env=env)
            status('report_and_additional_figures_running')
            subprocess.run([sys.executable,str(ROOT/'deliver_results.py'),'--run',run_id],cwd=REPO,stdout=log,stderr=subprocess.STDOUT,check=True,env=env)
        status('complete')
    except BaseException as error:
        state['error']=repr(error);status('failed');raise
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--run',default='run-v2-20261102');main(p.parse_args().run)

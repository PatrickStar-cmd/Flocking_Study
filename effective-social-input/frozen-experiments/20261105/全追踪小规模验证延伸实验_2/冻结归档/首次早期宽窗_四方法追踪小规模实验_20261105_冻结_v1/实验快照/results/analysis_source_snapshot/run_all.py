"""Run or resume simulations, analyses, audits and final deliverable checks."""
from pathlib import Path
import subprocess
import sys

ROOT=Path(__file__).resolve().parent


def main():
    for name in ('run_experiment.py','analyze_results.py','verify_results.py','build_report.py','finalize_delivery.py','verify_delivery.py'):
        print('RUN',name,flush=True)
        subprocess.run([sys.executable,'-u',str(ROOT/name)],cwd=ROOT,check=True)


if __name__=='__main__':main()

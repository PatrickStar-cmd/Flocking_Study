"""One-time provenance-preserving export, not needed to run the delivered package."""
from pathlib import Path
import json, shutil, hashlib

SOURCE = Path('D:/大三上专业课/卓越计划')
PACKAGE = Path(__file__).resolve().parents[1]
RESULTS = SOURCE / 'effective-social-input/results'
LATEST = RESULTS / 'task3-kf-boundary-v1/non-gaussian-four-group-validation/20261021-v1'

def copy(src, dst):
    assert src.is_file(), src
    assert not dst.exists(), dst
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(src, dst)

def main():
    old_runners = {'prepareNonGaussianFourGroupValidation.m', 'runNonGaussianFourGroupValidation.m', 'runBearingIMMClosedLoopPilot.m'}
    for src in (LATEST / 'source-snapshot').glob('*.m'):
        folder = 'archive/matlab' if src.name in old_runners else 'matlab'
        copy(src, PACKAGE / folder / src.name)
    extras = ['testBearingIMMTrack.m', 'testNoisyBearingKFPilot.m', 'testNonGaussianBearingKFPilot.m',
              'buildBearingReplaySequence.m', 'replayBearingModels.m']
    for name in extras:
        copy(SOURCE/'effective-social-input'/name, PACKAGE/'matlab'/name)
    runners = ['runNoisyBearingKFPilot.m', 'runNonGaussianBearingKFPilot.m', 'runBearingIMMReplayPilot.m', 'runOracleBearingDiagnostic.m']
    for name in runners:
        copy(SOURCE/'effective-social-input'/name, PACKAGE/'archive/matlab'/name)
    scripts = ['analyzeNonGaussianFourGroupValidation.py', 'plotNonGaussianValidationSummary.py',
               'analyze_oracle_bearing_diagnostic_20261021.py', 'recheck_noisy_bearing_kf_pilot.py',
               'diagnose_noisy_kf_rmse_20261014.m', 'recheck_noisy_kf_rmse_diagnosis.py',
               'recheck_non_gaussian_bearing_kf_pilot.py', 'recheck_bearing_imm_replay_pilot.py',
               'report_bearing_imm_replay_pilot_20261014.py', 'recheck_bearing_imm_closed_loop_pilot.py',
               'report_bearing_imm_closed_loop_pilot_20261014.py']
    for name in scripts:
        frozen = LATEST/'analysis-scripts'/name
        copy(frozen if frozen.exists() else SOURCE/'handoff'/name, PACKAGE/'archive/analysis'/name)
    for name in ['metadata.json','pre-run-plan.md','software-test-audit.json','methodology-audit.json','delivery-audit.json']:
        copy(LATEST/name, PACKAGE/'provenance/validation64'/name)
    for name in ['statistics.json','contrasts.csv','method-means.csv','per-seed-differences.csv',
                 '64新种子独立验证报告.md','结论与研究决策.md']:
        copy(LATEST/'independent-analysis'/name, PACKAGE/'reports/validation64'/name)
    for folder in ['noisy-bearing-kf-pilot-20261014-v2','non-gaussian-bearing-kf-pilot-20261014-v1',
                   'bearing-imm-replay-pilot-20261014-v1','bearing-imm-closed-loop-pilot-20261014-v1',
                   'oracle-bearing-diagnostic-20261021-v1']:
        for name in ['metadata.json','completion.json','summary.csv','frozen-selection.json']:
            src=RESULTS/folder/name
            if src.exists(): copy(src,PACKAGE/'provenance'/folder/name)
    meta=json.loads((LATEST/'metadata.json').read_text(encoding='utf-8'))
    (PACKAGE/'configs').mkdir(exist_ok=True)
    (PACKAGE/'configs/base-config.json').write_text(json.dumps(meta['baseConfig'],ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    manifest=[]
    for src in (LATEST/'source-snapshot').glob('*.m'):
        dst=PACKAGE/('archive/matlab' if src.name in old_runners else 'matlab')/src.name
        h=hashlib.sha256(dst.read_bytes()).hexdigest()
        assert h==meta['sourceHashes'][src.name.replace('.','_')]
        manifest.append(dict(file=dst.relative_to(PACKAGE).as_posix(),sha256=h,source='validation64/source-snapshot/'+src.name))
    (PACKAGE/'provenance/frozen-source-manifest.json').write_text(json.dumps(manifest,indent=2)+'\n',encoding='utf-8')
    print(f'Copied {len(manifest)} exact frozen MATLAB files; no MAT/NPZ trajectories exported.')

if __name__=='__main__': main()

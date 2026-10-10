"""Build a new C execution-design snapshot; never run MATLAB or overwrite output."""
from pathlib import Path
import csv
import hashlib
import json
import subprocess

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
DESIGN = ROOT / 'reviews/sq1-c-occlusion-protocol-v1'
AUDIT = ROOT / 'reviews/sq1-b-r3-audit-v1'
SOURCE = AUDIT / 'source'
B_COMMIT = '54fb85d91151e0b28e8e570318b4a9937d6a0425'
B_PREFIX = 'effective-social-input/results/sq1-fourbranch-protocol12-seed1-r3/'


def sha(data):
    return hashlib.sha256(data).hexdigest()


def git(*args):
    return subprocess.check_output(['git', *args], cwd=ROOT)


def save(relative, data):
    p = OUT / relative
    p.parent.mkdir(parents=True, exist_ok=True)
    with p.open('xb') as f:
        f.write(data)


def save_json(relative, value):
    save(relative, (json.dumps(value, ensure_ascii=False, indent=2) + '\n').encode('utf-8'))


def windows(d):
    return {
        'occlusion': [6201, 6200+d],
        'recovery100': [6201+d, 6300+d],
        'recovery300': [6201+d, 6500+d],
        'event400': [6201, 6600],
    }


def main():
    assert git('branch', '--show-current').decode().strip() == 'YCwork'
    protocol_commit = git('rev-parse', '1e11d64').decode().strip()
    proposal = json.loads((DESIGN/'protocol.json').read_text(encoding='utf-8-sig'))['proposal']
    tests = json.loads((AUDIT/'C_tests.json').read_text(encoding='utf-8-sig'))
    assert len(tests) == 15 and all(t['status'] == 'PASS' for t in tests)
    generated = ['execution_config.json', 'SOURCE_LOCK.json', 'RUN_BINDINGS.json',
                 'T0_BINDINGS.json', 'ADAPTER_CHECK.json', 'SHA256SUMS.json']
    assert not any((OUT/f).exists() for f in generated), 'New output required; do not overwrite a freeze.'
    frozen_sources = []
    for name in ['protocol.json', 'RUN_MANIFEST.csv', 'T0_TEST_MANIFEST.csv',
                 'PROTOCOL.md', 'ACCEPTANCE.md', 'DATA_CONTRACT.md',
                 'A_B_INPUT_INTERFACE.md', 'B_DELIVERY_CHECKLIST.md', 'EXECUTION_ORDER.md']:
        p = DESIGN/name
        data = p.read_bytes()
        assert data == git('show', f'{protocol_commit}:reviews/sq1-c-occlusion-protocol-v1/{name}')
        target = 'references/protocol/'+name
        save(target, data)
        frozen_sources.append({'path': target, 'sha256': sha(data), 'git_commit': protocol_commit,
                               'git_path': f'reviews/sq1-c-occlusion-protocol-v1/{name}'})
    evidence = ['REVIEW.md', 'B_NEXT_TASKS.md', 'C_tests.json', 'C_matlab.log', 'C_audit.json',
                'saved_mat_audit.json', 'delivery_audit.json', 'source_audit.json']
    for name in evidence:
        data = (AUDIT/name).read_bytes()
        target = 'references/c_acceptance/'+name
        save(target, data)
        frozen_sources.append({'path': target, 'sha256': sha(data), 'git_commit': None,
                               'local_origin': f'reviews/sq1-b-r3-audit-v1/{name}'})
    manifest_path = SOURCE/'results/sq1-fourbranch-protocol12-seed1-r3-results-v4/source_manifest.json'
    source_manifest = json.loads(manifest_path.read_text(encoding='utf-8-sig'))
    core = {
        'angularSetMetrics.m', 'buildFreshOracle.m', 'buildRingConnectivity.m',
        'checkpointToStart.m', 'controllerInputFromObservation.m', 'decodeNeuralActivity.m',
        'defaultFourBranchConfig.m', 'fixedTargetIds.m', 'fourBranchInputStep.m',
        'initializeCoreState.m', 'initializeInputState.m', 'makeProtocolMasks.m',
        'resumeProtocolFromCheckpoint.m', 'sensorObservation.m', 'simulateProtocolSegment.m',
    }
    source_files = []
    for r in source_manifest['sourceFiles']:
        relative = r['relativePath']
        raw = git('show', f'{B_COMMIT}:{B_PREFIX}{relative}')
        assert sha(raw) == r['SHA256'].lower()
        source_files.append({'path': relative, 'sha256': sha(raw),
                             'role': 'frozen_core' if relative in core else 'reference_or_adapter'})
    assert len([r for r in source_files if r['role'] == 'frozen_core']) == len(core)
    lock = {'schema': 'sq1-source-lock-1', 'accepted_branch': 'feature/sq1-r3-software-package',
            'accepted_commit': B_COMMIT, 'package_git_path': B_PREFIX.rstrip('/'),
            'protocol_commit': protocol_commit, 'A_commit': '992e71aba3b26896335b52a3182cf01cfd542a9d',
            'source_files': source_files, 'frozen_inputs': frozen_sources,
            'batch_adapter_commit': None,
            'rule': 'Bind the future adapter SHA in a new release record. Core changes invalidate automatic reuse of seed1 acceptance.'}
    save_json('SOURCE_LOCK.json', lock)
    base = {
        'NAgents': 20, 'NNeurons': 100, 'NSteps': 6600, 'dt': 0.3,
        'arenaSize': 1000, 'speedScale': 0.05, 'beta': 1000, 'globalInhibition': 0,
        'connectivityExponent': 0.5, 'receptiveFieldWidth': 0.4, 'allocentric': True,
        'k': 8, 'neighborCount': 8, 'h': 4, 'A0': 0.228, 'w0': 0.0285,
        'checkpointSteps': [6000, 6200], 'maskReadyStep': 6200,
        'occlusionStartStep': 6201, 'recordStride': 1, 'recordFromStep': 6000,
        'protocolVersion': '1.2', 'headingThreshold': 1e-9,
        'initialXY': [], 'initialHeadings': [], 'initialNeuralState': [],
    }
    assert proposal['simulation_seeds'] == list(range(1, 21))
    config = {
        'schema': 'sq1-batch-execution-freeze-1', 'freeze_id': 'sq1-c-batch-freeze-v1',
        'date': '2026-10-09', 'owner': '杨骋', 'implementation_owner': '王浩权',
        'protocol_revision': '1.2', 'evidence_type': 'development_exploratory',
        'state': 'DESIGN_FROZEN_BATCH_ADAPTER_PENDING',
        'accepted_implementation_commit': B_COMMIT,
        'batch_adapter_commit': None,
        'authorization': {'seed1_software_accepted': True, 'batch_adapter_check_allowed': True,
                          'full_batch_allowed': False,
                          'release_requires': ['adapter_commit_pinned', 'source_hashes_checked',
                                               'seed1_28_main_and_4_T0_accepted', 'environment_and_resources_recorded'],
                          'reason': 'Base software accepted; multi-condition adapter has not yet been delivered.'},
        'base_matlab_config': base,
        'sensing': {'bearing_noise_std': 0, 'distance_weight': 'equal', 'turn_cue': False,
                    'hold_age_limit': None, 'predictor': None},
        'design': {'seeds': list(range(1, 21)), 'mask_rules': proposal['mask_rules'],
                   'mask_display_names': {'random': 'random', 'max_hidden_R1': 'concentrated', 'min_hidden_R1': 'dispersed'},
                   'duration_steps': [10, 34, 100], 'modes': proposal['modes'],
                   'mask_stream': proposal['mask_stream'],
                   'target_rule': proposal['target_rule'],
                   'mask_geometry_step': 6200, 'mask_geometry_time': 'step_start_last_legal_observation',
                   'mask_geometry_source': 'checkpoint6200.lastBearing or saved start-of-step6200 position',
                   'mask_geometry_forbidden': 'checkpoint6200 end-of-step position used as a new observation'},
        'counts': {'seeds': 20, 'main_runs': 560, 'fresh_references': 20, 'nonfresh_runs': 540,
                   'T0_runs': 80, 'main_runs_per_seed': 28, 'T0_runs_per_seed': 4,
                   'logical_checkpoints': 40, 'mask_sets': 60, 'same_state_replay_contexts': 180,
                   'required_nonfresh_window_rows': 2160, 'primary_nonfresh_rows': 540},
        'files': {'main_manifest': 'references/protocol/RUN_MANIFEST.csv',
                  'T0_manifest': 'references/protocol/T0_TEST_MANIFEST.csv',
                  'main_bindings': 'RUN_BINDINGS.json', 'T0_bindings': 'T0_BINDINGS.json',
                  'source_lock': 'SOURCE_LOCK.json', 'adapter_check': 'ADAPTER_CHECK.json'},
        'summary': {'key': ['run_id', 'comparison_mask_rule', 'comparison_duration_steps', 'window'],
                    'required_modes': ['drop', 'drop-normalized', 'hold'],
                    'fresh_reference_window': 'Use the comparison duration, never fresh duration=0, to align windows.',
                    'optional_fresh_context_rows': 720,
                    'optional_fresh_rows_are_independent_runs': False,
                    'primary_endpoint': 'pairDistanceIncrease in recovery100',
                    'pairDistanceIncrease': 'mean(D_branch-D_fresh); normalized by arenaSize',
                    'GOLoss': 'mean(GO_fresh-GO_branch)', 'LOLoss': 'mean(LO_fresh-LO_branch)',
                    'contrasts': ['drop-normalized minus drop', 'hold minus drop-normalized']},
        'statistics': {'resampling_unit': 'seed_complete_paired_block', 'bootstrap_replicates': 20000,
                       'bootstrap_seed': 20260929, 'shared_resample_indices': True,
                       'interval': 'linear_percentile_95', 'multiple_comparison_adjustment': False,
                       'report': ['per-seed', 'mean', 'median', 'range', 'paired_bootstrap_interval'],
                       'scientific_success_threshold': None},
        'scheduler': {'unit': 'seed', 'order': list(range(1, 21)), 'execution': 'serial', 'workerCount': 0,
                      'parallel_release': 'Only after measured memory and schedule/RNG equivalence are recorded.',
                      'T0_before_nonfresh_forks': True, 'fresh_runs_once_per_seed': True,
                      'prefix_steps_reused': 200, 'fork_steps': 400, 'assemble_steps': 600,
                      'technical_failure_action': 'stop and preserve logs/partial data; no silent exclusions',
                      'bad_scientific_effect_is_technical_failure': False},
        'output_policy': {'new_directory_required': True, 'root': 'effective-social-input/results/',
                          'adapter_name_pattern': 'sq1-fourbranch-batch-adapter-v1-<run_id>',
                          'main_name_pattern': 'sq1-fourbranch-pilot-v1-<run_id>',
                          'overwrite_allowed': False, 'reuse_r3_smoke_by_default': False,
                          'resume': 'Only complete tasks with matching release/config/source/checkpoint/data hashes; append a new execution log.',
                          'state_values': ['PLANNED', 'RUNNING', 'GENERATED_PENDING_REVIEW', 'ACCEPTED', 'FAILED']},
    }
    save_json('execution_config.json', config)
    bindings = {}
    for manifest, result in [('RUN_MANIFEST.csv', 'RUN_BINDINGS.json'), ('T0_TEST_MANIFEST.csv', 'T0_BINDINGS.json')]:
        with (DESIGN/manifest).open(encoding='utf-8-sig', newline='') as f:
            rows = list(csv.DictReader(f))
        records = []
        for r in rows:
            seed, d = int(r['seed']), int(r['duration_steps'])
            cfg = dict(base, seed=seed, mode=r['mode'], maskRule=r['mask_rule'],
                       occlusionDurationSteps=d, recoveryStartStep=6201+d)
            records.append(dict(r, seed=seed, duration_steps=d, branch_steps=int(r['branch_steps']),
                                status='PLANNED', matlab_config=cfg,
                                windows=windows(d) if d else {},
                                output_relative_path=('runs/' if manifest == 'RUN_MANIFEST.csv' else 'T0/')+r['run_id']+'.mat'))
        bindings[result] = records
        save_json(result, {'schema': 'sq1-run-bindings-1', 'manifest': 'references/protocol/'+manifest,
                           'accepted_base_commit': B_COMMIT, 'execution_started': False, 'runs': records})
    check = {'schema': 'sq1-batch-adapter-check-1', 'seed': 1, 'evidence_type': 'software_test',
             'status': 'PLANNED', 'selection_reason': 'Prespecified seed1; all conditions retained regardless of outcomes.',
             'main_run_ids': [r['run_id'] for r in bindings['RUN_BINDINGS.json'] if r['seed'] == 1],
             'T0_run_ids': [r['run_id'] for r in bindings['T0_BINDINGS.json'] if r['seed'] == 1],
             'required_mask_sets': 3, 'same_state_replay_contexts': 9,
             'required_nonfresh_window_rows': 108, 'primary_nonfresh_rows': 27,
             'optional_fresh_context_rows': 36,
             'release_record': 'A new record must pin actual adapter SHA, source/config hashes, environment, C acceptance evidence and full_batch_allowed.'}
    save_json('ADAPTER_CHECK.json', check)
    inventory = []
    for p in sorted(OUT.rglob('*')):
        if p.is_file() and p.name not in ['SHA256SUMS.json', 'VALIDATION.json'] and '__pycache__' not in p.parts:
            inventory.append({'path': p.relative_to(OUT).as_posix(), 'sha256': sha(p.read_bytes()), 'bytes': p.stat().st_size})
    save_json('SHA256SUMS.json', {'schema': 'sq1-freeze-hashes-1', 'files': inventory,
                                'excluded': ['SHA256SUMS.json', 'VALIDATION.json']})
    print(f'Created execution freeze: {len(bindings["RUN_BINDINGS.json"])} main + {len(bindings["T0_BINDINGS.json"])} T0 planned keys; no simulation started.')


if __name__ == '__main__':
    main()

"""Read-only verification of frozen design; --record creates a NEW validation log."""
from pathlib import Path
from collections import Counter
from itertools import product
import argparse
import copy
import csv
import hashlib
import json
import subprocess

BASE = Path(__file__).resolve().parent
ROOT = BASE.parents[1]
PIN = '54fb85d91151e0b28e8e570318b4a9937d6a0425'
MASKS = ['random', 'max_hidden_R1', 'min_hidden_R1']
MODES = ['drop', 'drop-normalized', 'hold']


def read_json(path):
    return json.loads((BASE/path).read_text(encoding='utf-8-sig'))


def require(condition, reason):
    if not condition:
        raise ValueError(reason)


def validate_runs(rows, t0=False):
    expected = {}
    for s in range(1, 21):
        tag = f'S{s:03d}'
        if t0:
            for mode in ['fresh']+MODES:
                expected[f'TEST-{tag}-d0-{mode}'] = (s, 'not_applicable', 0, mode)
        else:
            expected[f'{tag}-fresh'] = (s, 'not_applicable', 0, 'fresh')
            for mask, d, mode in product(MASKS, [10, 34, 100], MODES):
                expected[f'{tag}-{mask}-d{d}-{mode}'] = (s, mask, d, mode)
    keys = [r['run_id'] for r in rows]
    require(len(keys) == len(set(keys)), 'duplicate run_id')
    require(set(keys) == set(expected), 'missing or extra run keys')
    paths = [r['output_relative_path'] for r in rows]
    require(len(paths) == len(set(paths)), 'duplicate output path')
    constants = {'NAgents': 20, 'NNeurons': 100, 'NSteps': 6600, 'dt': .3,
                 'arenaSize': 1000, 'speedScale': .05, 'beta': 1000, 'globalInhibition': 0,
                 'connectivityExponent': .5, 'receptiveFieldWidth': .4, 'allocentric': True,
                 'k': 8, 'neighborCount': 8, 'h': 4, 'A0': .228, 'w0': .0285,
                 'checkpointSteps': [6000, 6200], 'maskReadyStep': 6200,
                 'occlusionStartStep': 6201, 'recordStride': 1, 'recordFromStep': 6000,
                 'protocolVersion': '1.2', 'headingThreshold': 1e-9,
                 'initialXY': [], 'initialHeadings': [], 'initialNeuralState': []}
    for r in rows:
        s, mask, d, mode = expected[r['run_id']]
        require((r['seed'], r['mask_rule'], r['duration_steps'], r['mode']) == (s, mask, d, mode), 'key/config mismatch')
        require(r['checkpoint_id'] == f'S{s:03d}-step6000', 'checkpoint mismatch')
        require(r['event_checkpoint_id'] == f'S{s:03d}-step6200', 'event checkpoint mismatch')
        require(r['reference_id'] == f'S{s:03d}-fresh', 'reference mismatch')
        require(r['status'] == 'PLANNED' and r['branch_steps'] == 600, 'incorrect planned status/steps')
        c = r['matlab_config']
        require(all(c[k] == v for k, v in constants.items()), 'changed static scientific config')
        require((c['seed'], c['maskRule'], c['occlusionDurationSteps'], c['mode']) == (s, mask, d, mode), 'MATLAB dynamic config mismatch')
        require(c['recoveryStartStep'] == 6201+d, 'wrong recovery boundary')
        p = Path(r['output_relative_path'])
        require(not p.is_absolute() and '..' not in p.parts, 'unsafe output path')
        if d == 0:
            require(r['windows'] == {}, 'fresh/T0 cannot define a positive-duration event')
        else:
            # Check inclusive bounds and lengths independently from builder function.
            w = r['windows']
            require(set(w) == {'occlusion', 'recovery100', 'recovery300', 'event400'}, 'missing window')
            for name, length in [('occlusion', d), ('recovery100', 100), ('recovery300', 300), ('event400', 400)]:
                left, right = w[name]
                require(right-left+1 == length and 6001 <= left <= right <= 6600, 'window length/range mismatch')
            require(w['occlusion'] == [6201, 6200+d], 'wrong event window')
            require(w['event400'] == [6201, 6600], 'wrong event400 window')
            require(w['recovery100'][0] == w['occlusion'][1]+1, 'recovery not immediate')
            require(w['recovery300'][0] == w['recovery100'][0], 'inconsistent recovery windows')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--record', action='store_true')
    args = parser.parse_args()
    config = read_json('execution_config.json')
    lock = read_json('SOURCE_LOCK.json')
    main_runs = read_json('RUN_BINDINGS.json')['runs']
    t0_runs = read_json('T0_BINDINGS.json')['runs']
    validate_runs(main_runs)
    validate_runs(t0_runs, True)
    for manifest, rows in [(config['files']['main_manifest'], main_runs), (config['files']['T0_manifest'], t0_runs)]:
        with (BASE/manifest).open(encoding='utf-8-sig', newline='') as f:
            original = list(csv.DictReader(f))
        require(len(original) == len(rows), 'manifest count differs')
        for original_row, row in zip(original, rows):
            require(all(str(row[k]) == v for k, v in original_row.items()), 'snapshot/binding drift')
    require(config['accepted_implementation_commit'] == lock['accepted_commit'] == PIN, 'base SHA changed')
    require(config['batch_adapter_commit'] is None and lock['batch_adapter_commit'] is None, 'adapter SHA must be pending')
    require(config['authorization']['seed1_software_accepted'] is True, 'software acceptance missing')
    require(config['authorization']['full_batch_allowed'] is False, 'premature full-batch release')
    require(config['authorization']['batch_adapter_check_allowed'] is True, 'adapter check incorrectly blocked')
    require(config['scheduler']['workerCount'] == 0 and config['scheduler']['execution'] == 'serial', 'unmeasured parallelism')
    require(config['design']['mask_rules'] == MASKS, 'mask aliases changed')
    require(config['statistics']['resampling_unit'] == 'seed_complete_paired_block' and config['statistics']['bootstrap_seed'] == 20260929 and config['statistics']['bootstrap_replicates'] == 20000, 'statistics drift')
    require(config['design']['mask_geometry_time'] == 'step_start_last_legal_observation', 'incorrect mask geometry timing')
    c_tests = read_json('references/c_acceptance/C_tests.json')
    require(len(c_tests) == 15 and all(t['status'] == 'PASS' for t in c_tests), 'acceptance evidence not 15 PASS')
    hashes = read_json('SHA256SUMS.json')['files']
    require(len({r['path'] for r in hashes}) == len(hashes), 'duplicate hash entries')
    for r in hashes:
        raw = (BASE/r['path']).read_bytes()
        require(hashlib.sha256(raw).hexdigest() == r['sha256'] and len(raw) == r['bytes'], f'hash mismatch: {r["path"]}')
    for r in lock['frozen_inputs']:
        raw = (BASE/r['path']).read_bytes()
        require(hashlib.sha256(raw).hexdigest() == r['sha256'], 'input lock hash mismatch')
        if r['git_commit']:
            original = subprocess.check_output(['git', 'show', f'{r["git_commit"]}:{r["git_path"]}'], cwd=ROOT)
            require(original == raw, 'protocol Git snapshot mismatch')
    require(len(lock['source_files']) == 43, 'source lock incomplete')
    for r in lock['source_files']:
        raw = subprocess.check_output(['git', 'show', f'{PIN}:{lock["package_git_path"]}/{r["path"]}'], cwd=ROOT)
        require(hashlib.sha256(raw).hexdigest() == r['sha256'], 'source Git pin mismatch')
    fresh = [r for r in main_runs if r['mode'] == 'fresh']
    nonfresh = [r for r in main_runs if r['mode'] != 'fresh']
    counts = {'main_runs': len(main_runs), 'T0_runs': len(t0_runs), 'fresh_references': len(fresh),
              'nonfresh_runs': len(nonfresh), 'required_nonfresh_window_rows': sum(len(r['windows']) for r in nonfresh),
              'primary_nonfresh_rows': sum('recovery100' in r['windows'] for r in nonfresh),
              'logical_checkpoints': len({r['checkpoint_id'] for r in main_runs})+len({r['event_checkpoint_id'] for r in main_runs}),
              'mask_sets': len({(r['seed'], r['mask_rule']) for r in nonfresh}),
              'same_state_replay_contexts': len({(r['seed'], r['mask_rule'], r['duration_steps']) for r in nonfresh})}
    require(all(config['counts'][k] == v for k, v in counts.items()), 'declared count mismatch')
    require(Counter(r['seed'] for r in main_runs) == Counter({s: 28 for s in range(1,21)}), 'incomplete per-seed pairing')
    require(Counter(r['seed'] for r in t0_runs) == Counter({s: 4 for s in range(1,21)}), 'incomplete per-seed T0')
    check = read_json('ADAPTER_CHECK.json')
    require(set(check['main_run_ids']) == {r['run_id'] for r in main_runs if r['seed'] == 1}, 'adapter main subset wrong')
    require(set(check['T0_run_ids']) == {r['run_id'] for r in t0_runs if r['seed'] == 1}, 'adapter T0 subset wrong')
    # Negative controls exercise meaningful faults without changing any saved file.
    mutations = []
    a = copy.deepcopy(main_runs); a[1]['run_id'] = a[0]['run_id']; mutations.append(('duplicate_key', a))
    a = copy.deepcopy(main_runs); a[1]['reference_id'] = 'S002-fresh'; mutations.append(('wrong_pairing', a))
    a = copy.deepcopy(main_runs); a[1]['windows']['recovery100'][0] += 1; mutations.append(('off_by_one', a))
    a = copy.deepcopy(main_runs); a[1]['matlab_config']['A0'] = .24; mutations.append(('changed_amplitude', a))
    a = copy.deepcopy(main_runs); a.pop(); mutations.append(('missing_condition', a))
    negatives = {}
    for label, rows in mutations:
        try:
            validate_runs(rows)
        except ValueError:
            negatives[label] = 'EXPECTED_REJECTION'
        else:
            raise ValueError('Negative control was accepted: '+label)
    result = {'status': 'PASS', 'date': '2026-10-09', 'evidence_type': 'execution_design_validation',
              'independent_counts': counts, 'hashed_files_checked': len(hashes),
              'pinned_sources_checked': len(lock['source_files']), 'core_files_locked': sum(r['role']=='frozen_core' for r in lock['source_files']),
              'seed1_software_gates': '15 PASS', 'negative_controls': negatives,
              'simulations_started': False, 'full_batch_allowed': False}
    if args.record:
        with (BASE/'VALIDATION.json').open('x', encoding='utf-8', newline='\n') as f:
            json.dump(result, f, ensure_ascii=False, indent=2); f.write('\n')
    print(json.dumps(result, ensure_ascii=False))


if __name__ == '__main__':
    main()

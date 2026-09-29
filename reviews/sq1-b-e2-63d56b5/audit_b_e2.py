"""Read pinned Git blobs and audit B's E2 tables; never run the simulator."""
import csv
import hashlib
import io
import json
import math
from pathlib import Path
import statistics
import subprocess
from collections import defaultdict, Counter

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parents[1]
SHA = '63d56b51cb74bca509b0fb7dc29d1d5997d3a1ae'
PREFIX = 'effective-social-input/run_B_E2/'

def git(*args):
    return subprocess.check_output(['git', '-C', str(REPO), *args])

def sha256(data):
    return hashlib.sha256(data).hexdigest()

def table(path):
    return list(csv.DictReader(path.open(encoding='utf-8-sig', newline='')))

def text(path):
    data = path.read_bytes()
    return data.decode('utf-16') if data[:2] in (b'\xff\xfe', b'\xfe\xff') else data.decode('utf-8-sig')

files = git('ls-tree', '-r', '--name-only', SHA).decode().splitlines()
selected = [p for p in files if p.startswith(PREFIX) or
            (p.startswith('effective-social-input/') and p.count('/') == 1 and
             (p.endswith('.m') or p.rsplit('/', 1)[-1] in
              ('AGENTS.md', 'AGENT_WORK_PROTOCOL.md', 'RESEARCH_ROADMAP.md')))]
manifest = []
for p in selected:
    raw = git('show', f'{SHA}:{p}')
    dst = ROOT / 'source' / p
    dst.parent.mkdir(parents=True, exist_ok=True)
    if dst.exists():
        assert dst.read_bytes() == raw, f'Refuse to overwrite changed snapshot: {dst}'
    else:
        dst.write_bytes(raw)
    manifest.append({'path': p, 'bytes': len(raw), 'sha256': sha256(raw)})
(ROOT / 'source_manifest.json').write_text(json.dumps({
    'remote': 'https://github.com/PatrickStar-cmd/Flocking_Study.git',
    'branch': 'WHQ_work', 'commit': SHA,
    'scope': 'Complete run_B_E2 directory, top-level MATLAB sources and research instructions',
    'files': manifest}, ensure_ascii=False, indent=2), encoding='utf-8')

base = ROOT / 'source' / PREFIX
rows = table(base / 'run_B_E2_results/B_E2_runs.csv')
grouped = table(base / 'run_B_E2_results/B_E2_grouped_summary.csv')
keys = ['strengthMode', 'selectionPolicy', 'k', 'seed']
metrics = [k for k in rows[0] if k not in ['experiment', *keys]]
by_config = defaultdict(list)
by_group = defaultdict(list)
for r in rows:
    by_config[tuple(r[k] for k in keys)].append(r)
    by_group[tuple(r[k] for k in ['experiment', *keys[:-1]])].append(r)

duplicates = []
for key, rr in by_config.items():
    if len(rr) > 1:
        same = all(all(r[m] == rr[0][m] for m in metrics) for r in rr[1:])
        duplicates.append({'key': key, 'occurrences': len(rr),
                           'experiments': [r['experiment'] for r in rr], 'metrics_identical': same})

max_error = defaultdict(float)
for g in grouped:
    rr = by_group[tuple(g[k] for k in ['experiment', *keys[:-1]])]
    assert len(rr) == int(g['nRuns'])
    assert len(set(r['seed'] for r in rr)) == len(rr)
    for m in metrics:
        vals = [float(r[m]) for r in rr]
        avg, sd = statistics.mean(vals), statistics.stdev(vals)
        expected = {'mean': avg, 'sd': sd,
                    'ci95_low': avg - 1.96 * sd / math.sqrt(len(rr)),
                    'ci95_high': avg + 1.96 * sd / math.sqrt(len(rr))}
        for suffix, value in expected.items():
            max_error[suffix] = max(max_error[suffix], abs(value - float(g[m + '_' + suffix])))

equivalent = defaultdict(list)
for r in rows:
    policy = 'balanced' if r['selectionPolicy'] == 'angle-dispersed' else r['selectionPolicy']
    mode = 'full-equivalent' if r['k'] == '19' else r['strengthMode']
    equivalent[(mode, policy, r['k'], r['seed'])].append(r)

amplitude_errors = []
for r in rows:
    expected = .228 if r['strengthMode'] == 'fixed-total' else .012 * int(r['k'])
    amplitude_errors.append(abs(float(r['meanInputAmplitude']) - expected))

audit = json.loads((base / 'run_B_E2_audit.json').read_text(encoding='utf-8-sig'))
hash_checks = []
for entry in audit['files']:
    old = entry['path']
    mapped = 'run_B_E2_README.md' if old == 'README.md' else old.replace('results/', 'run_B_E2_results/', 1)
    p = base / mapped
    raw = p.read_bytes()
    normalized = raw.replace(b'\r\n', b'\n')
    variants = {'exact': raw, 'LF': normalized, 'CRLF': normalized.replace(b'\n', b'\r\n')}
    match = [k for k, v in variants.items() if sha256(v) == entry['sha256'].lower()]
    hash_checks.append({'declared_path': old, 'actual_path': mapped, 'match': match,
                        'declared_bytes': entry['bytes'], 'actual_bytes': len(raw)})

sim = text(base / 'simulateBGroupE2.m')
balanced = sim.split("case 'balanced'", 1)[1].split("case 'angle-dispersed'", 1)[0].strip()
dispersed = sim.split("case 'angle-dispersed'", 1)[1].split("case 'angle-concentrated'", 1)[0].strip()
runner = text(base / 'run_B_E2_experiments.m')
raw_paths = sorted((base / 'run_B_E2_results/raw').glob('*.mat'))
result = {
    'source_commit': SHA, 'downloaded_files': len(manifest),
    'B_directory_files': sum(p.startswith(PREFIX) for p in selected),
    'rows': len(rows), 'groups_with_experiment_labels': len(by_group),
    'seeds': sorted(set(int(r['seed']) for r in rows)),
    'unique_keys_without_experiment': len(by_config),
    'unique_keys_after_known_algorithm_equivalences': len(equivalent),
    'duplicate_config_groups': duplicates,
    'all_repeated_config_metrics_identical': all(d['metrics_identical'] for d in duplicates),
    'all_algorithm_equivalent_metrics_identical': all(
        all(all(r[m] == rr[0][m] for m in metrics) for r in rr[1:]) for rr in equivalent.values()),
    'nonfinite_numeric_cells': sum(not math.isfinite(float(r[m])) for r in rows for m in metrics),
    'group_statistic_max_absolute_errors': dict(max_error),
    'amplitude_formula_max_error': max(amplitude_errors),
    'active_count_max_error': max(abs(float(r['meanActiveNeighbors']) - int(r['k'])) for r in rows),
    'balanced_and_dispersed_code_identical': balanced == dispersed,
    'raw_mat_count': len(raw_paths), 'raw_seeds_by_filename': dict(Counter(p.stem.rsplit('seed', 1)[-1] for p in raw_paths)),
    'runner_utf16': (base / 'run_B_E2_experiments.m').read_bytes()[:2] == b'\xff\xfe',
    'runner_contains_NUL': '\x00' in runner,
    'original_audit_hash_checks': hash_checks,
    'limitations': 'CSV and static-code audit; MATLAB inspection reported separately; no population simulation or significance claim.'
}
assert len(rows) == 480 and len(by_group) == 24
assert all(len(rr) == 20 and {int(r['seed']) for r in rr} == set(range(1, 21)) for rr in by_group.values())
assert result['nonfinite_numeric_cells'] == 0
assert all(0 <= float(r[m]) <= 1 for r in rows for m in ['meanGO', 'meanLO', 'meanAngularCoverage', 'meanInputConcentration'])
assert result['all_repeated_config_metrics_identical']
assert result['all_algorithm_equivalent_metrics_identical']
assert max(max_error.values()) < 2e-8
assert result['amplitude_formula_max_error'] < 1e-12
assert result['active_count_max_error'] == 0
assert all('exact' in c['match'] for c in hash_checks)
result['assertion_status'] = 'PASS: counts, complete seed sets, finite values, bounds, duplicate consistency, moments/CI arithmetic, amplitude/count formulas, seven exact source hashes'
(ROOT / 'audit_results.json').write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding='utf-8')
for k, v in result.items():
    if k != 'duplicate_config_groups':
        print(k, json.dumps(v, ensure_ascii=False))

"""Create derived condition maps and an evidence-scoped seed ledger. No simulations."""
import csv
import hashlib
import io
import json
from collections import defaultdict
from pathlib import Path
import subprocess

OUT = Path(__file__).resolve().parent
REPO = OUT.parents[1]
B_SHA = '63d56b51cb74bca509b0fb7dc29d1d5997d3a1ae'
B_PATH = 'effective-social-input/run_B_E2/run_B_E2_results/B_E2_runs.csv'
B_SOURCE = REPO / 'reviews/sq1-b-e2-63d56b5/source' / B_PATH

def git(*args):
    return subprocess.check_output(['git', '-C', str(REPO), *args])

def write_csv(name, fields, rows):
    with (OUT / name).open('w', newline='', encoding='utf-8-sig') as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)

raw = B_SOURCE.read_bytes()
assert raw == git('show', B_SHA + ':' + B_PATH), 'B snapshot differs from pinned commit'
rows = list(csv.DictReader(io.StringIO(raw.decode('utf-8-sig'))))
metrics = [k for k in rows[0] if k not in ('experiment', 'strengthMode', 'selectionPolicy', 'k', 'seed')]

def named_key(r):
    return (r['strengthMode'], r['selectionPolicy'], int(r['k']))

def semantic_key(r):
    mode, policy, k = named_key(r)
    policy = 'balanced' if policy == 'angle-dispersed' else policy
    if k == 19:
        mode, policy = 'full-equivalent', 'all'
    return (mode, policy, k)

configs = {k: f'C{i:02d}' for i, k in enumerate(sorted({semantic_key(r) for r in rows}), 1)}
first = {}
members = defaultdict(list)
mapping = []
canonical = []
for line, r in enumerate(rows, 2):
    key = semantic_key(r)
    config_id = configs[key]
    run_key = (config_id, int(r['seed']))
    run_id = f'{config_id}-S{int(r["seed"]):03d}'
    members[config_id].append((line, r))
    reasons = []
    if r['selectionPolicy'] == 'angle-dispersed': reasons.append('angle-dispersed is code-identical to balanced')
    if int(r['k']) == 19: reasons.append('all 19 neighbors; original and fixed-total weights coincide')
    if run_key in first:
        prior_line, prior = first[run_key]
        assert all(r[m] == prior[m] for m in metrics), ('equivalence mismatch', line, prior_line)
        reasons.append('reuses canonical run; do not count as independent replication')
    else:
        first[run_key] = (line, r)
        canonical.append({'canonical_config_id': config_id, 'canonical_run_id': run_id,
                          'canonical_strength': key[0], 'canonical_policy': key[1],
                          'source_csv_line': line, **r})
    mapping.append({'source_commit': B_SHA, 'source_csv_line': line,
                    'experiment': r['experiment'], 'strengthMode': r['strengthMode'],
                    'selectionPolicy': r['selectionPolicy'], 'k': r['k'], 'seed': r['seed'],
                    'canonical_config_id': config_id, 'canonical_run_id': run_id,
                    'representative_csv_line': first[run_key][0],
                    'is_representative': line == first[run_key][0],
                    'reason': '; '.join(reasons) or 'first representative; original labels retained'})

write_csv('B_E2_CONDITION_MAP.csv', list(mapping[0]), mapping)
write_csv('B_E2_CANONICAL_RUNS.csv', list(canonical[0]), canonical)
catalog = []
for key, cid in sorted(configs.items(), key=lambda x: x[1]):
    mm = members[cid]
    mode, policy, k = key
    catalog.append({'canonical_config_id': cid, 'strength': mode, 'policy': policy, 'k': k,
                    'expected_total_weight': .012*k if mode == 'original' else .228,
                    'source_row_count': len(mm), 'unique_seed_count': len({r['seed'] for _, r in mm}),
                    'seeds': '1-20', 'experiment_aliases': ';'.join(sorted({r['experiment'] for _, r in mm})),
                    'policy_aliases': ';'.join(sorted({r['selectionPolicy'] for _, r in mm})),
                    'selection_information': 'distance-assisted' if policy in ('balanced', 'nearest') else 'no distance used for selection',
                    'selection_timing': 'all neighbors' if k == 19 else 'reselected every simulation step',
                    'evidence': 'exploratory/development; previously used seeds'})
write_csv('B_E2_CONFIG_CATALOG.csv', list(catalog[0]), catalog)

# Enumerate all currently fetched remote heads, pin hashes before reading.
prior_audit = OUT / 'ledger_audit.json'
if prior_audit.exists():
    refs = json.loads(prior_audit.read_text(encoding='utf-8'))['remote_ref_snapshot']
else:
    refs = git('for-each-ref', '--format=%(refname:short) %(objectname)', 'refs/remotes/origin').decode().splitlines()
    refs = dict(line.split(' ', 1) for line in refs if line.split(' ', 1)[0] != 'origin')
csv_blobs = {}
for ref, sha in refs.items():
    for line in git('ls-tree', '-r', sha).decode().splitlines():
        meta, path = line.split('\t', 1)
        if not path.endswith('.csv') or not path.startswith('effective-social-input/'):
            continue
        _, typ, blob = meta.split()
        if typ != 'blob': continue
        v = csv_blobs.setdefault((path, blob), {'path': path, 'blob': blob, 'refs': [], 'commits': []})
        v['refs'].append(ref)
        v['commits'].append(sha)

inventory = []
usage = []
ledger = defaultdict(list)
for index, ((path, blob), v) in enumerate(sorted(csv_blobs.items()), 1):
    dataset_id = f'D{index:03d}'
    data = git('cat-file', 'blob', blob)
    reader = csv.DictReader(io.StringIO(data.decode('utf-8-sig')))
    rr = list(reader)
    headers = reader.fieldnames or []
    seed_columns = [h for h in headers if h.lower() in ('seed', 'randomseed', 'random_seed')]
    observed = defaultdict(int)
    invalid = []
    for r in rr:
        for col in seed_columns:
            s = r.get(col, '')
            try:
                seed = int(s)
            except (TypeError, ValueError):
                invalid.append(s)
                continue
            observed[seed] += 1
    status = 'observed_result_seed_column' if seed_columns else 'no_explicit_seed_column_not_cleared'
    parent_summary = '/seed-results/' in path
    inventory.append({'dataset_id': dataset_id, 'path': path, 'git_blob': blob,
                      'sha256': hashlib.sha256(data).hexdigest(), 'bytes': len(data),
                      'refs': ';'.join(v['refs']), 'pinned_commits': ';'.join(v['commits']),
                      'row_count': len(rr), 'seed_columns': ';'.join(seed_columns),
                      'seed_count': len(observed), 'seeds': ';'.join(map(str, sorted(observed))),
                      'invalid_seed_cells': len(invalid), 'status': status,
                      'overlap_warning': 'per-seed export; may duplicate aggregate summary' if parent_summary else 'may overlap other exports/analyses; do not sum as trials'})
    for seed, n in sorted(observed.items()):
        event = {'seed': seed, 'dataset_id': dataset_id, 'path': path, 'git_blob': blob,
                 'observed_rows': n, 'status': 'already_observed_not_fresh',
                 'role': 'result or derived analysis; exact development/calibration role requires dataset documentation',
                 'replication_warning': 'row count is not independent trial count'}
        usage.append(event)
        ledger[seed].append(event)

write_csv('SOURCE_DATASET_INVENTORY.csv', list(inventory[0]), inventory)
write_csv('SQ1_SEED_USAGE.csv', list(usage[0]), usage)
seed_rows = [{'seed': seed, 'status': 'already_observed_not_fresh',
              'source_file_count': len(events), 'dataset_ids': ';'.join(e['dataset_id'] for e in events),
              'eligibility': 'do not designate as new held-out seed; reuse for development with disclosure',
              'independence_note': 'same numeric seed across different configurations is not automatically identical trajectory or independent replication'}
             for seed, events in sorted(ledger.items())]
write_csv('SQ1_SEED_LEDGER.csv', list(seed_rows[0]), seed_rows)

assert len(rows) == 480 and len(mapping) == 480 and len(canonical) == 240 and len(catalog) == 12
assert len({named_key(r) + (int(r['seed']),) for r in rows}) == 280
assert all(sum(r['canonical_config_id'] == cid for r in canonical) == 20 for cid in configs.values())
assert len({m['source_csv_line'] for m in mapping}) == 480
assert all(m['canonical_run_id'] in {r['canonical_run_id'] for r in canonical} for m in mapping)
assert not any(r['invalid_seed_cells'] for r in inventory)
assert all(s in ledger for s in range(1, 21))
def ranges(values):
    out = []
    for x in sorted(values):
        if out and x == out[-1][1] + 1: out[-1][1] = x
        else: out.append([x,x])
    return ', '.join(str(a) if a==b else f'{a}-{b}' for a,b in out)

audit = {'source_B_commit': B_SHA, 'source_B_sha256': hashlib.sha256(raw).hexdigest(),
         'remote_ref_snapshot': refs, 'original_rows': len(rows), 'named_run_keys': 280,
         'canonical_run_keys': len(canonical), 'canonical_configs': len(catalog),
         'csv_file_versions_inspected': len(inventory),
         'csv_file_versions_with_seed_column': sum(bool(r['seed_columns']) for r in inventory),
         'observed_seed_count': len(ledger), 'observed_seed_ranges': ranges(ledger),
         'coverage_limit': 'Committed remote CSV explicit seed columns only. Binary-only results, script defaults/plans, JSON-only provenance, unpushed/ignored files and unshared runs are not cleared. Absence is not proof of freshness. No new seeds reserved.',
         'checks': 'PASS: source identity, full row mapping, equivalent metrics identical, 12x20 canonical keys, explicit seed parsing, B seed coverage',
         'outputs': {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in OUT.glob('*.csv')}}
(OUT / 'ledger_audit.json').write_text(json.dumps(audit, ensure_ascii=False, indent=2), encoding='utf-8')
print(json.dumps(audit, ensure_ascii=False, indent=2))

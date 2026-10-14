"""Join the frozen canonical ledger to pinned MAT artifacts; do not alter inputs."""
import csv
import hashlib
import io
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
LEDGER = HERE.parent / 'sq1-c-condition-ledger-v1'
SOURCE = HERE.parent / 'sq1-whq-update-0036487'
COMMIT = '0036487818d891c23d357f7af49030f323f4c9af'

def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def read_csv(path):
    with path.open(encoding='utf-8-sig', newline='') as f:
        return list(csv.DictReader(f))

def write_new(path, data):
    if path.exists():
        assert path.read_bytes() == data, f'Refusing changed output: {path}'
    else:
        path.write_bytes(data)

def write_csv(name, rows):
    buf = io.StringIO(newline='')
    w = csv.DictWriter(buf, fieldnames=list(rows[0]))
    w.writeheader()
    w.writerows(rows)
    write_new(HERE / name, buf.getvalue().encode('utf-8-sig'))

def rel(path):
    return path.relative_to(REPO).as_posix()

def prepare():
    mapping = read_csv(LEDGER / 'B_E2_CONDITION_MAP.csv')
    canonical = read_csv(LEDGER / 'B_E2_CANONICAL_RUNS.csv')
    summary = read_csv(SOURCE / 'source/run_B_1_03_summary_results/B_E2_runs.csv')
    provenance = json.loads((SOURCE / 'source_audit.json').read_text(encoding='utf-8'))
    assert provenance['new_commit'] == COMMIT
    manifest = {r['file']: r for r in provenance['raw_manifest']}
    reps = {r['canonical_run_id']: r for r in canonical}
    assert len(mapping) == len(summary) == len(manifest) == 480
    assert len(canonical) == len(reps) == 240
    names, groups, rows = set(), defaultdict(list), []
    def filename(r):
        return f"{r['experiment']}_{r['strengthMode']}_{r['selectionPolicy']}_k{r['k']}_seed{r['seed']}.mat"
    for m in mapping:
        src = summary[int(m['source_csv_line']) - 2]
        for col in ('experiment', 'strengthMode', 'selectionPolicy', 'k', 'seed'):
            assert m[col] == src[col]
        rep = reps[m['canonical_run_id']]
        assert rep['source_csv_line'] == m['representative_csv_line']
        name = filename(m)
        path = SOURCE / 'raw' / name
        assert name not in names
        names.add(name)
        h = digest(path)
        assert h == manifest[name]['sha256']
        is_rep = m['source_csv_line'] == rep['source_csv_line']
        assert is_rep == (m['is_representative'] == 'True')
        reasons = []
        if m['experiment'] != rep['experiment']: reasons.append('experiment_label')
        if m['selectionPolicy'] != rep['selectionPolicy']: reasons.append('balanced_angle-dispersed_alias')
        if m['strengthMode'] != rep['strengthMode']: reasons.append('k19_strength_equivalence')
        row = dict(canonical_config_id=m['canonical_config_id'], canonical_run_id=m['canonical_run_id'],
                   canonical_strength=rep['canonical_strength'], canonical_policy=rep['canonical_policy'],
                   k=int(m['k']), seed=int(m['seed']), source_commit=COMMIT,
                   source_csv_line=int(m['source_csv_line']), experiment=m['experiment'],
                   strengthMode=m['strengthMode'], selectionPolicy=m['selectionPolicy'],
                   is_representative=int(is_rep), alias_reason=';'.join(reasons) or 'representative',
                   raw_path=rel(path), raw_sha256=h, raw_bytes=path.stat().st_size,
                   source_archive=manifest[name]['archive'], archive_member=manifest[name]['member'],
                   representative_raw_path=rel(SOURCE / 'raw' / filename(rep)))
        rows.append(row)
        groups[m['canonical_run_id']].append(row)
    assert names == set(manifest)
    for key, group in groups.items():
        assert sum(r['is_representative'] for r in group) == 1
    for config in {r['canonical_config_id'] for r in rows}:
        assert {r['seed'] for r in rows if r['canonical_config_id'] == config} == set(range(1, 21))
    write_csv('RAW_FILE_MAP.csv', rows)
    print(f'PREPARED: {len(rows)} verified files -> {len(groups)} canonical keys; hashes match.')
    return rows, groups

def finalize(rows, groups):
    verification = json.loads((HERE / 'alias_checks.json').read_text(encoding='utf-8'))
    checks = verification['checks']
    assert len(checks) == 480
    checked = {c['raw_path']: c for c in checks}
    assert len(checked) == 480 and set(checked) == {r['raw_path'] for r in rows}
    assert verification['allPassed'] and verification['nonRepresentativeComparisons'] == 240
    index = []
    for r in rows:
        c = checked[r['raw_path']]
        assert c['canonical_run_id'] == r['canonical_run_id']
        assert c['representative_raw_path'] == r['representative_raw_path']
        assert c['allPayloadExact'] and c['normalizedConfigExact']
        if not r['is_representative']: continue
        index.append(dict(r, source_label_count=len(groups[r['canonical_run_id']]),
                          alias_payload_verified=1, n_agents=c['n_agents'], n_neurons=c['n_neurons'],
                          n_steps=c['n_steps'], burn_in=c['burn_in'], dt=c['dt'],
                          time_series_length=c['time_series_length'], time_series_stride_steps=1,
                          time_series_dt=c['dt'], trajectory_frames=c['trajectory_frames'],
                          trajectory_stride_steps=c['trajectory_stride_steps'],
                          trajectory_dt=c['dt']*c['trajectory_stride_steps'],
                          time_series_fields=';'.join(c['time_series_fields']),
                          trajectory_fields=';'.join(c['trajectory_fields'])))
    index.sort(key=lambda r: (r['canonical_config_id'], r['seed']))
    assert len(index) == 240
    write_csv('CANONICAL_RAW_INDEX.csv', index)
    inputs = [LEDGER/'B_E2_CONDITION_MAP.csv', LEDGER/'B_E2_CANONICAL_RUNS.csv',
              SOURCE/'source_audit.json', SOURCE/'source/run_B_1_03_summary_results/B_E2_runs.csv']
    outputs = [HERE/'RAW_FILE_MAP.csv', HERE/'CANONICAL_RAW_INDEX.csv', HERE/'alias_checks.json']
    audit = dict(source_commit=COMMIT, evidence='software/data consistency audit; source experiments exploratory',
                 mapping_rows=len(rows), canonical_rows=len(index), alias_comparisons=240,
                 config_counts=dict(Counter(r['canonical_config_id'] for r in index)),
                 labels_per_canonical_run=dict(Counter(len(g) for g in groups.values())),
                 inputs={rel(p):digest(p) for p in inputs}, outputs={p.name:digest(p) for p in outputs},
                 scripts={p.name:digest(p) for p in [Path(__file__),HERE/'verify_aliases.m']},
                 payload_exact=True, normalized_config_exact=True)
    write_new(HERE/'audit.json', (json.dumps(audit, ensure_ascii=False, indent=2)+'\n').encode())
    print(json.dumps({k:v for k,v in audit.items() if k not in ('inputs','outputs','scripts')}, indent=2))

if __name__ == '__main__':
    rows, groups = prepare()
    if '--finalize' in sys.argv:
        finalize(rows, groups)

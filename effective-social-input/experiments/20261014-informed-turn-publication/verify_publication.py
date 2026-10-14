"""Read-only checks of saved development artifacts and their Git index bytes.

Writes a separate new publication report only when --out is provided.
Does not import the simulator or change any historical experiment output.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import subprocess
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
RESULTS = ROOT / 'effective-social-input/results/challenge2-informed-turn-v1-20261014'
SOURCE = ROOT / 'effective-social-input/challenge2-informed-turn-v1'
PROTOCOL = ROOT / 'effective-social-input/replication-protocols/20261014-informed-turn-n20-t18'


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run(check_index=False):
    audit = json.loads((RESULTS / 'audit.json').read_text(encoding='utf-8'))
    for rel, expected in audit['result_hashes'].items():
        assert sha256(RESULTS / rel.replace('\\', '/')) == expected, rel
    for rel, expected in audit['source_hashes'].items():
        assert sha256(SOURCE / rel.replace('\\', '/')) == expected, rel

    artifacts = [p for scope in [RESULTS, SOURCE, PROTOCOL]
                 for p in scope.rglob('*') if p.is_file() and '__pycache__' not in p.parts]
    assert len([p for p in RESULTS.rglob('*') if p.is_file()]) == 414
    assert max(p.stat().st_size for p in artifacts) < 100 * 1024**2

    # Build the directed graph using one-based human IDs, independently.
    graph = {}
    for i in range(1, 21):
        graph[i] = (({1, 2, 3} - {i}) | set(range(4, 17))) if i <= 3 else (
            {1, 2, 3} | {4 + ((i - 4 + r) % 17) for r in range(1, 13)})
    assert all(len(js) == 15 and i not in js for i, js in graph.items())
    selected = sum(map(len, graph.values()))
    hidden = sum(j <= 3 for i, js in graph.items() if i >= 4 for j in js)
    assert (selected, hidden, selected - hidden) == (300, 51, 249)

    stage = RESULTS / '06-development-comparison'
    raw_npz = list((stage / 'runs').glob('*.npz'))
    assert len(raw_npz) == 60
    samples = []
    for seed in [2026101401, 2026101402, 2026101403]:
        path = stage / 'runs' / f'seed-{seed}_turn_fresh_T-0.npz'
        with np.load(path) as z:
            assert z['xy'].shape == (801, 20, 2)
            assert z['velocity'].shape == (800, 20, 2)
            with path.with_suffix('.csv').open(encoding='utf-8-sig', newline='') as f:
                rows = list(csv.DictReader(f))
            assert len(rows) == 800
            n = 200
            vectors = []
            for vx, vy in z['velocity'][n]:
                speed = math.hypot(vx, vy)
                vectors.append((vx/speed, vy/speed) if speed > 1e-12 else (0., 0.))
            go = math.hypot(sum(x for x, _ in vectors), sum(y for _, y in vectors)) / 20
            error = abs(go - float(rows[n]['GO']))
            assert error < 1e-12
            samples.append(dict(seed=seed, zero_based_step=n, GO_abs_error=error))

    index_checks = 0
    if check_index:
        entries = subprocess.check_output(['git', 'ls-files', '-s', '-z'], cwd=ROOT)
        index = {}
        for entry in entries.split(b'\0'):
            if entry:
                info, path = entry.split(b'\t', 1)
                index[path.decode('utf-8')] = info.split()[1].decode('ascii')
        for p in artifacts:
            content = p.read_bytes()
            blob = hashlib.sha1(b'blob ' + str(len(content)).encode() + b'\0' + content).hexdigest()
            rel = p.relative_to(ROOT).as_posix()
            assert index.get(rel) == blob, ('Git index bytes differ or artifact missing', rel)
            index_checks += 1

    return dict(
        status='PASS', evidence='publication integrity and software checks; exploratory data unchanged',
        historical_result_hashes_verified=len(audit['result_hashes']),
        historical_source_hashes_verified=len(audit['source_hashes']),
        result_files=414, result_bytes=sum(p.stat().st_size for p in RESULTS.rglob('*') if p.is_file()),
        selected_edges=selected, masked_edges=hidden, fresh_edges=selected-hidden,
        final_stage_raw_runs=len(raw_npz), independent_GO_checks=samples,
        git_index_exact_byte_checks=index_checks,
        artifact_sha256={p.relative_to(ROOT).as_posix(): sha256(p) for p in artifacts},
        historical_status_note='Original audit commit/push/main flags describe experiment creation, not this publication.',
    )


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--check-index', action='store_true')
    parser.add_argument('--out', type=Path)
    args = parser.parse_args()
    result = run(args.check_index)
    if args.out:
        with args.out.open('x', encoding='utf-8') as f:
            json.dump(result, f, ensure_ascii=False, indent=2)
    print(json.dumps({k: v for k, v in result.items() if k != 'artifact_sha256'}, ensure_ascii=True))

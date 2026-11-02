"""Main reviewer: independent geometry, selection, noise and provenance audit."""
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
from scipy.io import loadmat


def wrap(x):
    return np.arctan2(np.sin(x), np.cos(x))


def audit(root):
    d = json.loads((root / 'design.json').read_text(encoding='utf-8'))
    manifest = json.loads((root / 'input_manifest.json').read_text(encoding='utf-8'))
    assert manifest['status'] == 'COMPLETED'
    assert len(manifest['inputs']) == 10 and len(d['cases']) == 15
    assert len(set(d['seeds'])) == 5
    arrays = {}
    max_geometry_error = 0.
    for entry in manifest['inputs']:
        path = root / entry['source']
        assert hashlib.file_digest(path.open('rb'), 'sha256').hexdigest() == entry['sha256']
        assert path.stat().st_size == entry['bytes']
        x = loadmat(path, squeeze_me=True)
        obs, truth = x['observed_bearing'], x['true_bearing']
        assert obs.shape == truth.shape == (2200, 140)
        assert np.isfinite(obs).all() and np.isfinite(truth).all()
        observer = x['observer_id'].astype(int) - 1
        target = x['target_id'].astype(int) - 1
        assert len(set(zip(observer, target))) == 140
        assert np.array_equal(observer, np.repeat(np.arange(20), 7))
        assert np.array_equal(target, (observer + np.tile(np.arange(1, 8), 20)) % 20)
        pos = x['positions']
        delta = pos[:, target, :] - pos[:, observer, :]
        delta -= 1000 * np.round(delta / 1000)
        independent = np.arctan2(delta[1], delta[0]).T
        error = np.max(np.abs(wrap(independent - truth)))
        assert error < 1e-12
        max_geometry_error = max(max_geometry_error, float(error))
        arrays[entry['source']] = x
    noise_max = 0.
    for seed in d['seeds']:
        seed = int(seed)
        base = arrays[f'raw/input/seed{seed}_base.mat']
        turn = arrays[f'raw/input/seed{seed}_turn.mat']
        assert np.array_equal(base['positions'][:, :, :1401], turn['positions'][:, :, :1401])
        assert np.array_equal(base['headings'][:, :1401], turn['headings'][:, :1401])
        assert np.array_equal(base['observed_bearing'][:1401], turn['observed_bearing'][:1401])
        bn = wrap(base['observed_bearing'] - base['true_bearing'])
        tn = wrap(turn['observed_bearing'] - turn['true_bearing'])
        max_diff = float(np.max(np.abs(bn - tn)))
        noise_max = max(noise_max, max_diff)
        assert max_diff < 3e-15
        rates = wrap(np.diff(base['true_bearing'], axis=0)) / d['dt']
        starts = np.arange(1201, 1682, 20)
        scores = np.array([np.sqrt(np.mean(rates[t-2:t+118] ** 2)) for t in starts])
        for label, index in [('S', np.argmin(scores)), ('V', np.argmax(scores))]:
            case = next(c for c in d['cases'] if c['seed'] == seed and c['stratum'] == label)
            assert case['t0'] == starts[index]
            assert abs(case['selection_score'] - scores[index]) < 1e-12
    audit = {
        'status': 'PASS', 'input_count': 10, 'case_count': 15,
        'independent_seeds': 5, 'directed_hidden_edges': 140,
        'geometry_max_wrapped_error_rad': max_geometry_error,
        'matched_noise_max_error_rad': noise_max,
        'truth_only_selection_independently_reproduced': True,
        'prefix_positions_headings_observations_exact_through_step': 1401,
    }
    dest = root / 'review'
    dest.mkdir(exist_ok=True)
    path = dest / 'input_audit.json'
    assert not path.exists()
    path.write_text(json.dumps(audit, indent=2), encoding='utf-8')
    print(json.dumps(audit))


if __name__ == '__main__':
    audit(Path(sys.argv[1]))

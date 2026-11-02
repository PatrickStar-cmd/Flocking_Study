"""Reviewer-run actual-observation parity against unchanged MATLAB functions."""
import json
from pathlib import Path
import sys

import numpy as np
from scipy.io import loadmat

from replay import replay


def verify(root):
    d = json.loads((root / 'design.json').read_text(encoding='utf-8'))
    refs = loadmat(root / 'raw/matlab_reference.mat', simplify_cells=True)['refs']
    maxima = {key: 0. for key in ('theta', 'p11', 'p12', 'p22', 'nis', 'age')}
    count = 0
    for ref in refs:
        case = next(c for c in d['cases'] if c['case_id'] == ref['case_id'])
        source = loadmat(root / case['source'], simplify_cells=True)
        obs = source['observed_bearing'][:, np.array(ref['edge_indices'], dtype=int)-1]
        missing = np.zeros(len(obs), dtype=bool)
        start = int(case['t0'])-1
        missing[start:start+int(ref['D'])] = True
        for method, q, tau in [('cv', d['q_cv'], None), ('gm', d['q_gm'], d['tau'])]:
            actual = replay(obs, missing, dt=d['dt'], q=q, R=d['R'], tau=tau)
            for key in maxima:
                expected = np.asarray(ref[method][key])
                result = np.asarray(actual[key])
                assert expected.shape == result.shape == obs.shape
                assert np.array_equal(np.isnan(expected), np.isnan(result)), (method, key)
                if key == 'theta':
                    difference = np.arctan2(np.sin(result-expected), np.cos(result-expected))
                    assert np.nanmax(np.abs(difference)) < 1e-10
                else:
                    difference = result-expected
                    assert np.allclose(result, expected, atol=1e-10, rtol=1e-10, equal_nan=True), (method, key)
                maxima[key] = max(maxima[key], float(np.nanmax(np.abs(difference))))
            count += 1
    report = {'status': 'PASS', 'method_case_count': count, 'edges_per_case': 3,
              'steps_per_edge': 2200, 'reference': 'stock MATLAB bearingTrackStep and bearingTrackStepGaussMarkov',
              'max_absolute_errors': maxima}
    path = root / 'review' / 'matlab_parity.json'
    assert not path.exists()
    path.write_text(json.dumps(report, indent=2), encoding='utf-8')
    print(json.dumps(report))


if __name__ == '__main__':
    verify(Path(sys.argv[1]))

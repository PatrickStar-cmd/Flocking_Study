"""Separate real bearing changes from minimum-image branch discontinuities."""
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from scipy.io import loadmat


def run(root):
    d = json.loads((root/'design.json').read_text(encoding='utf-8'))
    rows = []
    for case in d['cases']:
        x = loadmat(root/case['source'], simplify_cells=True)
        i = x['observer_id'].astype(int)-1
        j = x['target_id'].astype(int)-1
        delta = x['positions'][:, j, :]-x['positions'][:, i, :]
        delta -= 1000*np.round(delta/1000)
        # A nearest-periodic-image switch can reverse bearing at half-box distance.
        switches = np.any(np.abs(np.diff(delta, axis=2)) > 500, axis=0).T
        angles = np.diff(x['true_bearing'], axis=0)
        changes = np.arctan2(np.sin(angles), np.cos(angles))
        distances = np.linalg.norm(delta, axis=0).T
        for duration in d['durations']:
            start = int(case['t0'])-1
            sl = slice(start-1, start+int(duration)-1)
            switch = switches[sl]
            change = changes[sl]
            mask = ~switch
            rows.append({
                'case_id': case['case_id'], 'seed': int(case['seed']),
                'stratum': case['stratum'], 'D': duration,
                'edge_times': change.size, 'minimum_image_switch_count': int(switch.sum()),
                'near_nyquist_count': int((np.abs(change) >= .9*np.pi).sum()),
                'bearing_rate_rms': float(np.sqrt(np.mean(change**2))/d['dt']),
                'bearing_rate_rms_without_image_switches': float(np.sqrt(np.mean(change[mask]**2))/d['dt']),
                'minimum_range': float(distances[start:start+int(duration)].min()),
            })
    path = root/'review/geometry_events.csv'
    assert not path.exists()
    frame = pd.DataFrame(rows)
    frame.to_csv(path, index=False)
    print(frame.groupby(['stratum', 'D'])[['minimum_image_switch_count', 'near_nyquist_count']].sum().to_string())


if __name__ == '__main__':
    run(Path(sys.argv[1]))

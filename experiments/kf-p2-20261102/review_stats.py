"""Independent reviewer: recompute all reported numbers from raw predictions."""
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from scipy.io import loadmat


def close(actual, expected):
    assert np.allclose(actual, expected, rtol=1e-11, atol=1e-12, equal_nan=True), (actual, expected)


def run(root):
    design = json.loads((root/'design.json').read_text(encoding='utf-8'))
    manifest = json.loads((root/'raw/replay/manifest.json').read_text(encoding='utf-8'))
    summary = pd.read_csv(root/'analysis/summary.csv')
    contrasts = pd.read_csv(root/'analysis/contrasts.csv')
    groups = pd.read_csv(root/'analysis/grouped.csv')
    ages = pd.read_csv(root/'analysis/age_metrics.csv')
    assert len(summary) == 270 and len(contrasts) == 30 and len(groups) == 36
    assert len(ages) == 6300
    assert manifest['completed'] == 30 and manifest['errors'] == 0
    assert len(manifest['runs']) == 30
    counts = 0
    min_eigenvalue = np.inf
    geometry_rows = []
    for entry in manifest['runs']:
        assert entry['status'] == 'OK'
        path = root/'raw/replay'/entry['npz_file']
        assert hashlib.file_digest(path.open('rb'), 'sha256').hexdigest() == entry['npz_sha256']
        ca = next(c for c in design['cases'] if c['case_id'] == entry['case_id'])
        source_path = root/ca['source']
        assert hashlib.file_digest(source_path.open('rb'), 'sha256').hexdigest() == entry['source_sha256']
        source = loadmat(source_path, simplify_cells=True)
        raw = np.load(path, allow_pickle=False)
        duration = int(entry['D']); start = int(ca['t0'])-1
        steps = raw['window_steps'].astype(int)-1
        expected_steps = np.arange(start-200, start+duration+400)
        assert np.array_equal(steps, expected_steps)
        mask = raw['missing'].astype(bool)
        assert mask.sum() == duration and mask[200:200+duration].all()
        close(raw['truth'], source['true_bearing'][steps])
        observed = raw['observed_sanitized']
        assert np.isnan(observed[mask]).all()
        close(observed[~mask], source['observed_bearing'][steps][~mask])
        # Geometry stratification is supplementary, never a sample exclusion.
        i = source['observer_id'].astype(int)-1; j = source['target_id'].astype(int)-1
        delta = source['positions'][:, j, :]-source['positions'][:, i, :]
        delta -= 1000*np.round(delta/1000)
        jumps = np.any(np.abs(np.diff(delta, axis=2)) > 500, axis=0).T
        switch_edge = jumps[start-1:start+duration-1].any(axis=0)
        rmses = {}
        for method, key in [('Hold','Hold'), ('KF-CV','KF_CV'), ('KF-GM','KF_GM')]:
            err = raw['theta_'+key]-raw['truth']
            err = np.arctan2(np.sin(err), np.cos(err))
            close(err, raw['error_'+key])
            age = raw['age_'+key]
            assert np.all(np.isfinite(raw['theta_'+key]))
            if method == 'Hold':
                close(raw['theta_Hold'][200:200+duration], np.broadcast_to(source['observed_bearing'][start-1], (duration,140)))
                assert np.isnan(raw['p11_Hold']).all() and np.isnan(raw['nis_Hold']).all()
            else:
                cov = np.empty((*err.shape,2,2))
                cov[...,0,0] = raw['p11_'+key]; cov[...,0,1] = raw['p12_'+key]
                cov[...,1,0] = raw['p12_'+key]; cov[...,1,1] = raw['p22_'+key]
                eig = float(np.linalg.eigvalsh(cov).min())
                min_eigenvalue = min(min_eigenvalue, eig)
                assert eig >= -1e-12
                assert np.isnan(raw['nis_'+key][mask]).all()
            for window, sl in [('pre',slice(0,200)), ('occ',slice(200,200+duration)), ('recovery',slice(200+duration,None))]:
                row = summary[(summary.case_id == ca['case_id']) & (summary.D == duration) & (summary.method == method) & (summary.window == window)]
                assert len(row) == 1
                row = row.iloc[0]; e = err[sl]; a = age[sl]
                rmse = np.sqrt(np.mean(e*e))
                close(row.circular_rmse, rmse); close(row.circular_mae, np.mean(np.abs(e)))
                assert row.n_samples == e.size
                cutoff = a > design['max_age']
                close(row.cutoff_by_age_fraction, cutoff.mean())
                if method != 'Hold':
                    p = raw['p11_'+key][sl]
                    close(row.mean_P, p.mean())
                    close(row.mse_over_mean_P, np.mean(e*e)/p.mean())
                    close(row.approx95coverage, (np.abs(e) <= 1.96*np.sqrt(p)).mean())
                    close(row.wide_interval_fraction, (1.96*np.sqrt(p) >= np.pi).mean())
                    close(row.cutoff_by_P_fraction, (p > design['max_variance']).mean())
                    cutoff |= p > design['max_variance']
                    n = raw['nis_'+key][sl]; finite = n[np.isfinite(n)]
                    close(row.observed_fresh_meanNIS, finite.mean() if finite.size else np.nan)
                else:
                    assert np.isnan(row.approx95coverage) and np.isnan(row.mean_P)
                # On fresh steps the controller always uses a current observation.
                cutoff &= mask[sl,None]
                close(row.cutoff_union_fraction, cutoff.mean())
                active = ~cutoff
                assert np.array_equal(raw['active_'+key][sl], active)
                assert row.active_sample_count == active.sum()
                close(row.active_only_rmse, np.sqrt(np.mean(e[active]**2)))
                if window == 'occ':
                    rmses[method] = rmse
                    first = np.flatnonzero(cutoff.any(axis=1))
                    close(row.first_cutoff_age_step, first[0]+1 if len(first) else np.nan)
                counts += 1
            for index in range(duration):
                r = ages[(ages.case_id == ca['case_id']) & (ages.D == duration) & (ages.method == method) & (ages.age_step == index+1)]
                assert len(r) == 1
                close(r.iloc[0].circular_rmse, np.sqrt(np.mean(err[200+index]**2)))
            e = err[200:200+duration]
            geometry_rows.append({'case_id':ca['case_id'],'seed':int(ca['seed']),'stratum':ca['stratum'],'D':duration,'method':method,
                'edges_with_image_switch':int(switch_edge.sum()),
                'rmse_edges_without_image_switch':float(np.sqrt(np.mean(e[:,~switch_edge]**2))),
                'rmse_edges_with_image_switch':float(np.sqrt(np.mean(e[:,switch_edge]**2))) if switch_edge.any() else np.nan})
        row = contrasts[(contrasts.case_id == ca['case_id']) & (contrasts.D == duration)].iloc[0]
        for key, method in [('rmse_Hold','Hold'),('rmse_KF_CV','KF-CV'),('rmse_KF_GM','KF-GM')]:
            close(row[key], rmses[method])
        close(row.diff_KF_CV_minus_Hold, rmses['KF-CV']-rmses['Hold'])
        close(row.diff_KF_GM_minus_Hold, rmses['KF-GM']-rmses['Hold'])
        close(row.diff_KF_GM_minus_KF_CV, rmses['KF-GM']-rmses['KF-CV'])
    for _, row in groups.iterrows():
        subset = contrasts[(contrasts.stratum == row.stratum) & (contrasts.D == row.D)].sort_values('seed')
        assert len(subset) == subset.seed.nunique() == row.n_seeds == 5
        values = subset[row.metric].to_numpy()
        rng = np.random.default_rng(int(design['bootstrap_seed']))
        samples = values[rng.integers(0,5,size=(int(design['bootstrap_replicates']),5))].mean(axis=1)
        close(row['mean'], values.mean()); close(row['std'], values.std(ddof=1))
        close(row.ci_lower, np.percentile(samples,2.5)); close(row.ci_upper, np.percentile(samples,97.5))
    report = {'status':'PASS','raw_replays':30,'summary_rows_independently_recomputed':counts,
              'age_rmse_rows_independently_recomputed':len(ages),'paired_contrast_rows':30,
              'bootstrap_group_rows':36,'minimum_covariance_eigenvalue':min_eigenvalue,
              'all_sha256_verified':True,'independent_seeds_per_stratum':5}
    path = root/'review/statistical_audit.json'; assert not path.exists()
    path.write_text(json.dumps(report,indent=2),encoding='utf-8')
    geometry_path = root/'review/geometry_error_strata.csv'; assert not geometry_path.exists()
    pd.DataFrame(geometry_rows).to_csv(geometry_path,index=False)
    print(json.dumps(report))


if __name__ == '__main__':
    run(Path(sys.argv[1]))

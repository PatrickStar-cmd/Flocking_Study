#!/usr/bin/env python3
"""
analyze_kf_diagnostics.py

Streaming post-processing and statistical analysis for Flocking KF Diagnostic Experiment.
Complies strictly with Review Fixes 7, 8, 9, 11:
  - Primary response window is OCCLUSION (occ), secondary is RECOVERY (rec).
  - True bearing evaluated against EACH ALGORITHM'S OWN CLOSED-LOOP TRAJECTORY.
  - Distinguishes candidate vs actual injected estimator metrics.
  - summary.csv: Long-form schema partitioned by window (pre, occ, rec, common).
  - contrasts.csv: Paired seed-level contrasts (Hold vs KFmean primary, KFcov vs KFmean, vs Fresh, vs None).
  - estimator_diagnostics.csv: seed x D x mode x age breakdown (MAE, RMSE, P, coverage, fresh NIS).
  - Deterministic seed-level cluster bootstrap with fixed random seed.
  - Streaming memory architecture: large arrays (uRecord, positionsRecord) discarded immediately.
  - CLI: supports results_dir and --output <analysis_dir>, refuses overwriting existing analysis.
  - Validates full 45-run manifest (no duplicates, no omissions, no FAILED status).
"""

import os
import sys
import glob
import json
import argparse
import numpy as np
import pandas as pd
from scipy.io import loadmat
from scipy import stats

ANALYSIS_RANDOM_SEED = 20261101
T0 = 1401
DURATIONS = [20, 120]
EXPECTED_SEEDS = [2026110101, 2026110102, 2026110103, 2026110104, 2026110105]

def wrap_angle(angle):
    """Wrap angle to [-pi, pi]."""
    return (angle + np.pi) % (2 * np.pi) - np.pi

def circ_dist(a, b):
    """Shortest angular distance on circle."""
    return np.abs(wrap_angle(a - b))

def find_mat_files(results_dir):
    """Find MAT files in results_dir/raw/ or results_dir/."""
    raw_dir = os.path.join(results_dir, "raw")
    if os.path.isdir(raw_dir):
        files = glob.glob(os.path.join(raw_dir, "*.mat"))
        if files:
            return sorted(files)
    files = glob.glob(os.path.join(results_dir, "*.mat"))
    return sorted(files)

def inspect_manifest_and_files(results_dir, mat_files):
    """Check manifest status and ensure 45 complete runs without failure."""
    manifest_path = os.path.join(results_dir, "manifest.json")
    if not os.path.isfile(manifest_path):
        raise FileNotFoundError("A completed run manifest is required.")
    with open(manifest_path, "r", encoding="utf-8") as f:
        manifest = json.load(f)
    status = manifest.get("status", "")
    if status != "COMPLETED" or manifest.get("failed_count", 0):
        raise RuntimeError(f"Manifest is not successfully completed: {status}.")
    declared_seeds = np.atleast_1d(manifest.get("seeds", manifest.get("seed", []))).astype(int)
    entries = list(manifest.get("runs", {}).values())
    expected_count = 9 * len(declared_seeds)
    if len(set(declared_seeds)) != len(declared_seeds) or len(declared_seeds) not in (1, 5):
        raise ValueError("Expected one software-check seed or five development seeds.")
    if len(entries) != expected_count or len(mat_files) != expected_count:
        raise ValueError("Missing or extra condition files/manifest entries.")
    expected_paths = {os.path.normcase(os.path.abspath(os.path.join(results_dir, e['mat_file'])))
                      for e in entries if e.get('status') == 'COMPLETED'}
    actual_paths = {os.path.normcase(os.path.abspath(p)) for p in mat_files}
    if expected_paths != actual_paths or len(expected_paths) != expected_count:
        raise ValueError("Manifest entries must all be completed and match unique MAT files.")

    # Check for empty run list
    if not mat_files:
        raise FileNotFoundError(f"No MAT files found in {results_dir}")

def extract_run_summary(mat_path):
    """
    Stream a single MAT file, extract timeSeries and observation metrics,
    and return memory-light summary objects.
    """
    data = loadmat(mat_path, struct_as_record=False, squeeze_me=True)
    if 'result' in data:
        res = data['result']
    elif 'saveData' in data and hasattr(data['saveData'], 'result'):
        res = data['saveData'].result
    else:
        # Fallback to scanning fields
        found = None
        for k in ('result', 'res'):
            if k in data:
                found = data[k]
                break
        if found is None:
            raise KeyError(f"Could not find result struct in {mat_path}")
        res = found

    cfg = res.config
    seed = int(cfg.seed)
    dur = int(cfg.occlusionDurationSteps)
    mode_raw = str(cfg.predictionMode)

    # Normalize mode
    if dur == 0:
        mode = 'fresh'
    elif mode_raw in ('mean', 'kf-mean'):
        mode = 'kf-mean'
    elif mode_raw in ('kalman', 'kalman-unc'):
        mode = 'kalman'
    elif mode_raw in ('hold', 'hold-last'):
        mode = 'hold'
    elif mode_raw in ('none', 'drop', 'fresh-drop'):
        mode = 'none'
    else:
        mode = mode_raw

    ts = res.timeSeries
    N = int(cfg.NAgents)
    k_req = int(cfg.neighborCount)
    hidden_count = int(getattr(cfg, 'occlusionHiddenCount', 7))
    dt = float(cfg.dt)
    t0 = int(cfg.occlusionStartStep)

    # Extract time series arrays
    go = np.array(ts.globalOrder, dtype=np.float64)
    lo = np.array(ts.localOrder, dtype=np.float64)
    pd_norm = np.array(ts.meanPairDistanceNorm, dtype=np.float64)
    nd_norm = np.array(ts.meanNearestDistanceNorm, dtype=np.float64)
    an = np.array(ts.meanActiveNeighbors, dtype=np.float64)
    im = np.array(ts.meanInputMass, dtype=np.float64)
    bc = np.array(ts.meanBumpConcentration, dtype=np.float64)
    sp = np.array(ts.meanSpeed, dtype=np.float64)

    # Observation records for estimator diagnostics (streamed per-step)
    obs_diagnostics = []
    if hasattr(res, 'observationRecords') and dur > 0:
        obs_recs = res.observationRecords
        for k in range(1, dur + 1):
            step = t0 - 1 + k  # 1-based step
            rec = obs_recs[step - 1]

            true_b = np.array(rec.trueBearing)
            injected_b = np.array(rec.injectedBearing)
            fresh_mask = np.array(rec.fresh).astype(bool)
            active_mask = np.array(rec.active).astype(bool)
            selected_mask = np.isfinite(true_b)
            hidden_mask = selected_mask & (~fresh_mask)

            # Candidate bearing:
            # If KF mode, rawX contains tracker estimate regardless of truncation
            if mode in ('kalman', 'kf-mean') and hasattr(rec, 'rawX') and rec.rawX is not None:
                cand_b = np.array(rec.rawX)[:, :, 0]
            elif mode == 'hold':
                cand_b = injected_b
            else:
                cand_b = np.full_like(true_b, np.nan)

            # Raw P11 and cutoff reasons
            raw_p11 = np.array(rec.rawVariance) if hasattr(rec, 'rawVariance') else np.full_like(true_b, np.nan)
            cutoff = np.array(rec.cutoffReason) if hasattr(rec, 'cutoffReason') else np.full_like(true_b, -1)
            raw_nis = np.array(rec.rawNis) if hasattr(rec, 'rawNis') else np.full_like(true_b, np.nan)

            # Evaluate on hidden edges
            t_hid = true_b[hidden_mask]
            c_hid = cand_b[hidden_mask]
            inj_hid = injected_b[hidden_mask]
            p_hid = raw_p11[hidden_mask]
            cut_hid = cutoff[hidden_mask]
            act_hid = active_mask[hidden_mask]

            n_hid = len(t_hid)
            # Candidate errors
            valid_cand = np.isfinite(c_hid) & np.isfinite(t_hid)
            err_cand = circ_dist(c_hid[valid_cand], t_hid[valid_cand]) if np.any(valid_cand) else np.array([])
            mae_cand = np.mean(err_cand) if len(err_cand) > 0 else np.nan
            rmse_cand = np.sqrt(np.mean(err_cand**2)) if len(err_cand) > 0 else np.nan

            # Injected errors (only edges actually active and injected)
            valid_inj = act_hid & np.isfinite(inj_hid) & np.isfinite(t_hid)
            err_inj = circ_dist(inj_hid[valid_inj], t_hid[valid_inj]) if np.any(valid_inj) else np.array([])
            mae_inj = np.mean(err_inj) if len(err_inj) > 0 else np.nan
            rmse_inj = np.sqrt(np.mean(err_inj**2)) if len(err_inj) > 0 else np.nan

            # P and coverage (only applicable to KF modes)
            if mode in ('kalman', 'kf-mean'):
                valid_p = np.isfinite(p_hid)
                mean_p = np.mean(p_hid[valid_p]) if np.any(valid_p) else np.nan
                # Approximate 95% circular coverage: error <= 1.96 * sqrt(P)
                if len(err_cand) > 0 and np.any(valid_p):
                    sigma_est = np.sqrt(p_hid[valid_cand])
                    cov95 = np.mean(err_cand <= 1.96 * sigma_est)
                else:
                    cov95 = np.nan
            else:
                mean_p = np.nan
                cov95 = np.nan

            # Fresh NIS (only from visible edges)
            nis_fresh_vals = raw_nis[fresh_mask & np.isfinite(raw_nis)]
            mean_fresh_nis = (np.mean(nis_fresh_vals)
                              if mode in ('kalman', 'kf-mean') and len(nis_fresh_vals) > 0 else np.nan)

            # Cutoff counts
            age_cuts = int(np.sum(np.isin(cut_hid, [21, 23])))
            p_cuts = int(np.sum(np.isin(cut_hid, [22, 23])))

            obs_diagnostics.append({
                'seed': seed,
                'D': dur,
                'mode': mode,
                'step_in_occ': k,
                'age': k * dt,
                'n_hidden': n_hid,
                'mae_candidate': mae_cand,
                'rmse_candidate': rmse_cand,
                'mae_injected': mae_inj,
                'rmse_injected': rmse_inj,
                'mean_P11': mean_p,
                'coverage95_candidate': cov95,
                'fresh_nis_mean': mean_fresh_nis,
                'age_cutoff_count': age_cuts,
                'p_cutoff_count': p_cuts
                , 'age_trunc_fraction': age_cuts / n_hid if n_hid else np.nan
                , 'P_trunc_fraction': p_cuts / n_hid if n_hid else np.nan
                , 'joint_cutoff_count': int(np.sum(cut_hid == 23))
                , 'n_active_hidden': int(np.sum(valid_inj))
                , 'diagnostic_basis': 'own-closed-loop trajectory; not matched open-loop replay'
                , 'coverage95_note': 'Gaussian approximation; inspect circular errors when wide'
            })

    commit_value = str(data.get('commit', ''))
    run_id_value = str(data.get('run_id', os.path.basename(mat_path)))
    # Memory release: del large variables
    del res
    del data

    run_meta = {
        'commit': commit_value,
        'run_id': run_id_value,
        'seed': seed,
        'mode': mode,
        'D': dur,
        'N': N,
        'k': k_req,
        'hidden_count': hidden_count,
        'dt': dt,
        't0': t0,
        'result_path': os.path.relpath(mat_path, os.path.dirname(os.path.dirname(mat_path))).replace('\\', '/'),
        'go': go,
        'lo': lo,
        'pd_norm': pd_norm,
        'nd_norm': nd_norm,
        'an': an,
        'im': im,
        'bc': bc,
        'sp': sp,
        'obs_diagnostics': obs_diagnostics
    }
    return run_meta

def compute_window_means(arr, start_idx, end_idx):
    """Compute mean over [start_idx, end_idx] (0-based inclusive)."""
    return float(np.mean(arr[start_idx:end_idx + 1]))

def run_analysis(results_dir, output_dir=None):
    if output_dir is None:
        output_dir = os.path.join(results_dir, "analysis")

    if os.path.exists(output_dir) and os.listdir(output_dir):
        raise FileExistsError(f"Analysis output directory already exists and is not empty: {output_dir}. Refusing to overwrite.")
    os.makedirs(output_dir, exist_ok=True)

    mat_files = find_mat_files(results_dir)
    inspect_manifest_and_files(results_dir, mat_files)

    print(f"Loading and processing {len(mat_files)} MAT files in streaming mode...")
    runs = []
    runs_by_key = {}
    all_diagnostics = []

    for idx, mf in enumerate(mat_files):
        print(f"  [{idx+1}/{len(mat_files)}] Parsing {os.path.basename(mf)}...")
        meta = extract_run_summary(mf)
        runs.append(meta)
        key = (meta['seed'], meta['mode'], meta['D'])
        if key in runs_by_key:
            raise ValueError(f"Duplicate condition detected for key {key} in {mf}!")
        runs_by_key[key] = meta
        all_diagnostics.extend(meta['obs_diagnostics'])

    # Validation of complete A package
    seeds_found = sorted(list(set(r['seed'] for r in runs)))
    expected_keys = {(seed, 'fresh', 0) for seed in seeds_found}
    expected_keys |= {(seed, mode, duration) for seed in seeds_found
                      for mode in ('none', 'hold', 'kf-mean', 'kalman') for duration in DURATIONS}
    if set(runs_by_key) != expected_keys:
        raise ValueError('Each seed must have exactly Fresh plus four methods at both D values.')
    for run in runs:
        if run['t0'] != 1401 or run['dt'] != 0.3 or len(run['go']) != 2200:
            raise ValueError('Run configuration does not match the frozen P1 evaluation windows.')
        arrays = [run[k] for k in ('go','lo','pd_norm','nd_norm','an','im','bc','sp')]
        if any(len(a) != 2200 or not np.isfinite(a).all() for a in arrays):
            raise ValueError('Non-finite metrics or incomplete time series.')
    print(f"Unique seeds found: {seeds_found}")
    is_full_package = (len(seeds_found) == 5 and len(runs) == 45)
    if not is_full_package:
        print(f"[NOTE] Not a full 5-seed 45-run batch (n_seeds={len(seeds_found)}, total_runs={len(runs)}). Running exploratory mode.")

    # 1. Build summary.csv (Long table schema)
    summary_rows = []
    audit_path = os.path.join(results_dir, 'audit.json')
    with open(audit_path, encoding='utf-8') as stream:
        audit = json.load(stream)
    commit_id = audit['commit']
    if any(run['commit'] and run['commit'] != commit_id for run in runs):
        raise ValueError('Raw run commits do not match the run audit.')
    scenario = "allocentric_flocking_occlusion_kf"
    mask_rule = "index-fixed"

    for r in runs:
        seed = r['seed']
        mode = r['mode']
        dur = r['D']
        t0 = r['t0']
        dt = r['dt']

        # Determine target evaluation durations
        eval_durations = [dur] if dur > 0 else DURATIONS

        for eval_d in eval_durations:
            # Slices (0-based inclusive):
            # Pre: t0 - 200 : t0 - 1 -> [1200, 1399]
            pre_start, pre_end = t0 - 200 - 1, t0 - 1 - 1
            # Occ: t0 : t0 + eval_d - 1 -> [t0 - 1, t0 + eval_d - 2]
            occ_start, occ_end = t0 - 1, t0 + eval_d - 1 - 1
            # Rec: t0 + eval_d : t0 + eval_d + 399 -> [t0 + eval_d - 1, t0 + eval_d + 399 - 1]
            rec_start, rec_end = t0 + eval_d - 1, t0 + eval_d + 399 - 1
            # Common: 1401 : 1920 -> [1400, 1919] (520 steps)
            com_start, com_end = 1400, 1919

            windows = [
                ('pre', pre_start, pre_end),
                ('occ', occ_start, occ_end),
                ('rec', rec_start, rec_end),
                ('common', com_start, com_end)
            ]

            for w_name, w_s, w_e in windows:
                row = {
                    'run_id': r['run_id'],
                    'commit': commit_id,
                    'scenario': scenario,
                    'mask_rule': mask_rule,
                    'N': r['N'],
                    'k': r['k'],
                    'hidden_count': r['hidden_count'],
                    'effective_hidden_count': 0 if dur == 0 else r['hidden_count'],
                    'D_steps': dur,
                    'eval_D_steps': eval_d,
                    'dt': dt,
                    't0': t0,
                    'seed': seed,
                    'mode': mode,
                    'window': w_name,
                    'window_start_step': w_s + 1,
                    'window_end_step': w_e + 1,
                    'window_length_steps': w_e - w_s + 1,
                    'GO_mean': compute_window_means(r['go'], w_s, w_e),
                    'LO_mean': compute_window_means(r['lo'], w_s, w_e),
                    'pair_norm_mean': compute_window_means(r['pd_norm'], w_s, w_e),
                    'nearest_norm_mean': compute_window_means(r['nd_norm'], w_s, w_e),
                    'active_neighbors_mean': compute_window_means(r['an'], w_s, w_e),
                    'input_mass_mean': compute_window_means(r['im'], w_s, w_e),
                    'bump_concentration_mean': compute_window_means(r['bc'], w_s, w_e),
                    'speed_mean': compute_window_means(r['sp'], w_s, w_e),
                    'result_path': r['result_path']
                }
                # Common AUC definition: dt * sum(520 discrete step samples)
                if w_name == 'common':
                    row['common_auc_GO_discrete'] = dt * np.sum(r['go'][w_s:w_e + 1])
                    row['common_auc_pair_discrete'] = dt * np.sum(r['pd_norm'][w_s:w_e + 1])
                else:
                    row['common_auc_GO_discrete'] = np.nan
                    row['common_auc_pair_discrete'] = np.nan

                summary_rows.append(row)

    df_summary = pd.DataFrame(summary_rows)
    summary_path = os.path.join(output_dir, "summary.csv")
    df_summary.to_csv(summary_path, index=False)
    print(f"Saved: {summary_path}")

    # 2. Build contrasts.csv (Paired seed-level contrasts)
    # Primary window is OCCLUSION (occ), secondary is RECOVERY (rec)
    contrasts = []
    for seed in seeds_found:
        for dur in DURATIONS:
            r_fresh = runs_by_key.get((seed, 'fresh', 0))
            r_none = runs_by_key.get((seed, 'none', dur))
            r_hold = runs_by_key.get((seed, 'hold', dur))
            r_mean = runs_by_key.get((seed, 'kf-mean', dur))
            r_kalman = runs_by_key.get((seed, 'kalman', dur))

            if not (r_fresh and r_none and r_hold and r_mean and r_kalman):
                continue

            for w_name, w_s, w_e in [
                ('occ', T0 - 1, T0 + dur - 2),
                ('rec', T0 + dur - 1, T0 + dur + 399 - 1),
                ('common', 1400, 1919)
            ]:
                # Means in window
                go_fr = compute_window_means(r_fresh['go'], w_s, w_e)
                pd_fr = compute_window_means(r_fresh['pd_norm'], w_s, w_e)
                lo_fr = compute_window_means(r_fresh['lo'], w_s, w_e)

                go_no = compute_window_means(r_none['go'], w_s, w_e)
                pd_no = compute_window_means(r_none['pd_norm'], w_s, w_e)
                lo_no = compute_window_means(r_none['lo'], w_s, w_e)

                go_hd = compute_window_means(r_hold['go'], w_s, w_e)
                pd_hd = compute_window_means(r_hold['pd_norm'], w_s, w_e)
                lo_hd = compute_window_means(r_hold['lo'], w_s, w_e)

                go_km = compute_window_means(r_mean['go'], w_s, w_e)
                pd_km = compute_window_means(r_mean['pd_norm'], w_s, w_e)
                lo_km = compute_window_means(r_mean['lo'], w_s, w_e)

                go_kc = compute_window_means(r_kalman['go'], w_s, w_e)
                pd_kc = compute_window_means(r_kalman['pd_norm'], w_s, w_e)
                lo_kc = compute_window_means(r_kalman['lo'], w_s, w_e)

                # Baseline: Hold
                contrasts.append({
                    'seed': seed,
                    'D': dur,
                    'window': w_name,
                    'baseline': 'hold',
                    'method': 'kf-mean',
                    'is_primary_response': (w_name == 'occ'),
                    'benefitPair': pd_hd - pd_km,
                    'benefitGO': go_km - go_hd,
                    'benefitLO': lo_km - lo_hd
                })
                contrasts.append({
                    'seed': seed,
                    'D': dur,
                    'window': w_name,
                    'baseline': 'hold',
                    'method': 'kalman',
                    'is_primary_response': (w_name == 'occ'),
                    'benefitPair': pd_hd - pd_kc,
                    'benefitGO': go_kc - go_hd,
                    'benefitLO': lo_kc - lo_hd
                })

                # KFcov vs KFmean (Broadening comparison)
                contrasts.append({
                    'seed': seed,
                    'D': dur,
                    'window': w_name,
                    'baseline': 'kf-mean',
                    'method': 'kalman',
                    'is_primary_response': (w_name == 'occ'),
                    'benefitPair': pd_km - pd_kc,
                    'benefitGO': go_kc - go_km,
                    'benefitLO': lo_kc - lo_km
                })

                # Relative to Fresh (loss / degradation)
                contrasts.append({
                    'seed': seed,
                    'D': dur,
                    'window': w_name,
                    'baseline': 'fresh',
                    'method': 'none',
                    'is_primary_response': (w_name == 'occ'),
                    'benefitPair': pd_fr - pd_no,
                    'benefitGO': go_no - go_fr,
                    'benefitLO': lo_no - lo_fr
                })
                contrasts.append({
                    'seed': seed,
                    'D': dur,
                    'window': w_name,
                    'baseline': 'fresh',
                    'method': 'hold',
                    'is_primary_response': (w_name == 'occ'),
                    'benefitPair': pd_fr - pd_hd,
                    'benefitGO': go_hd - go_fr,
                    'benefitLO': lo_hd - lo_fr
                })
                contrasts.append({
                    'seed': seed,
                    'D': dur,
                    'window': w_name,
                    'baseline': 'fresh',
                    'method': 'kf-mean',
                    'is_primary_response': (w_name == 'occ'),
                    'benefitPair': pd_fr - pd_km,
                    'benefitGO': go_km - go_fr,
                    'benefitLO': lo_km - lo_fr
                })
                contrasts.append({
                    'seed': seed,
                    'D': dur,
                    'window': w_name,
                    'baseline': 'fresh',
                    'method': 'kalman',
                    'is_primary_response': (w_name == 'occ'),
                    'benefitPair': pd_fr - pd_kc,
                    'benefitGO': go_kc - go_fr,
                    'benefitLO': lo_kc - lo_fr
                })

    df_contrasts = pd.DataFrame(contrasts)
    df_contrasts['is_primary_response'] = ((df_contrasts['window'] == 'occ') &
        (df_contrasts['baseline'] == 'hold') & (df_contrasts['method'] == 'kf-mean'))
    contrasts_path = os.path.join(output_dir, "contrasts.csv")
    df_contrasts.to_csv(contrasts_path, index=False)
    print(f"Saved: {contrasts_path}")

    # 3. Seed-level grouped bootstrap and 95% CIs
    rng = np.random.default_rng(ANALYSIS_RANDOM_SEED)
    bootstrap_indices = rng.integers(0, len(seeds_found), size=(10000, len(seeds_found)))
    group_summary = []
    group_keys = ['D', 'window', 'baseline', 'method']

    for g_vals, g_df in df_contrasts.groupby(group_keys):
        g_df = g_df.sort_values('seed')
        if g_df['seed'].tolist() != seeds_found:
            raise ValueError('Every contrast must retain the full paired seed set.')
        dur, w_name, base_m, meth_m = g_vals
        n_seeds = len(g_df)

        for metric in ['benefitPair', 'benefitGO', 'benefitLO']:
            vals = g_df[metric].values
            mean_v = float(np.mean(vals))
            if n_seeds >= 2:
                std_v = float(np.std(vals, ddof=1))
                se_v = std_v / np.sqrt(n_seeds)
                t_crit = float(stats.t.ppf(0.975, df=n_seeds - 1))
                ci95_low = mean_v - t_crit * se_v
                ci95_high = mean_v + t_crit * se_v

                # Deterministic cluster bootstrap
                boot_means = vals[bootstrap_indices].mean(axis=1)
                boot_ci_low = float(np.percentile(boot_means, 2.5))
                boot_ci_high = float(np.percentile(boot_means, 97.5))
            else:
                std_v, se_v, ci95_low, ci95_high, boot_ci_low, boot_ci_high = (np.nan, np.nan, np.nan, np.nan, np.nan, np.nan)

            group_summary.append({
                'D': dur,
                'window': w_name,
                'baseline': base_m,
                'method': meth_m,
                'is_primary_response': (w_name == 'occ' and base_m == 'hold' and
                                        meth_m == 'kf-mean' and metric == 'benefitPair'),
                'metric': metric,
                'n_seeds': n_seeds,
                'mean': mean_v,
                'std': std_v,
                'se': se_v,
                'ci95_t_lower': ci95_low,
                'ci95_t_upper': ci95_high,
                'boot_ci95_lower': boot_ci_low,
                'boot_ci95_upper': boot_ci_high
            })

    df_group = pd.DataFrame(group_summary)
    group_path = os.path.join(output_dir, "contrasts_group_summary.csv")
    df_group.to_csv(group_path, index=False)
    print(f"Saved: {group_path}")

    # 4. Estimator Diagnostics CSV
    df_diag = pd.DataFrame(all_diagnostics)
    diag_path = os.path.join(output_dir, "estimator_diagnostics.csv")
    df_diag.to_csv(diag_path, index=False)
    with open(os.path.join(output_dir, 'analysis_audit.json'), 'x', encoding='utf-8') as stream:
        json.dump({'unique_trajectories': len(runs), 'seeds': seeds_found,
                   'primary': 'occ benefitPair: Hold pair minus KFmean pair',
                   'bootstrap_seed': ANALYSIS_RANDOM_SEED, 'bootstrap_replicates': 10000,
                   'resampling_unit': 'whole seed, identical index matrix across all contrasts and D',
                   'intervals': 'exploratory, pointwise 95%; no simultaneous or confirmatory claim',
                   'AUC': 'dt * sum of 520 post-update step samples',
                   'evidence': 'developmental exploratory'}, stream, indent=2)
    print(f"Saved: {diag_path}")
    print(f"Analysis completed successfully. All outputs written to: {output_dir}")

def main():
    parser = argparse.ArgumentParser(description="Streaming statistical analysis for Flocking KF Diagnostics.")
    parser.add_argument("results_dir", nargs="?", default=None, help="Path to results directory containing MAT files.")
    parser.add_argument("--output", dest="output_dir", default=None, help="Target directory for CSV outputs (refuses overwrite).")
    args = parser.parse_args()

    if args.results_dir is None:
        base_res = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "effective-social-input", "results"))
        candidates = sorted(glob.glob(os.path.join(base_res, "kf-diagnostic-v1-20261101-*")))
        if candidates:
            target_dir = candidates[-1]
        else:
            print("Error: No results directory specified and none found.")
            sys.exit(1)
    else:
        target_dir = os.path.abspath(args.results_dir)

    run_analysis(target_dir, args.output_dir)

if __name__ == "__main__":
    main()

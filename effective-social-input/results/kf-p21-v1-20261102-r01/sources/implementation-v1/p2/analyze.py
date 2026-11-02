"""experiments/kf-p2-20261102/analyze.py

Offline analysis and diagnostic visualization for Flocking Study P2 open-loop replay.
Generates comprehensive CSV summaries, paired estimator contrasts, grouped bootstrap
confidence intervals across seeds, and publication-quality diagnostic figures.

Output files are created under <RESULT_DIR>/analysis/ with strict non-overwrite guards.
"""

from __future__ import annotations

import argparse
import datetime
import hashlib
import json
import os
from pathlib import Path
import sys
from typing import Any, Dict, List, Optional, Tuple, Union

import matplotlib
matplotlib.use("Agg")  # Non-interactive headless backend
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.stats


def compute_sha256(filepath: Union[str, Path]) -> str:
    """Compute SHA256 hex digest of a file."""
    h = hashlib.sha256()
    with open(filepath, "rb") as f:
        while chunk := f.read(65536):
            h.update(chunk)
    return h.hexdigest()


def wrap_to_pi(angle: np.ndarray) -> np.ndarray:
    """Wrap angle differences to [-pi, pi)."""
    return np.arctan2(np.sin(angle), np.cos(angle))


def safe_save_csv(df: pd.DataFrame, path: Path) -> bool:
    """Save DataFrame to CSV, refusing overwrite if file already exists."""
    if path.exists():
        print(f"[SKIP CSV] Preserving existing file: {path.name}")
        return False
    df.to_csv(path, index=False)
    print(f"[WROTE CSV] {path.name} ({len(df)} rows)")
    return True


def safe_save_figure(fig: plt.Figure, path: Path) -> bool:
    """Save matplotlib figure to PNG, refusing overwrite if file already exists."""
    if path.exists():
        print(f"[SKIP FIG] Preserving existing figure: {path.name}")
        plt.close(fig)
        return False
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"[WROTE FIG] {path.name}")
    return True


def run_bootstrap_seeds(
    values: np.ndarray,
    n_replicates: int = 10000,
    rng_seed: int = 2026110200,
) -> Tuple[float, float, float, float]:
    """Compute bootstrap mean and 95% confidence interval by resampling seeds with replacement.

    Parameters
    ----------
    values : np.ndarray, shape [N_seeds]
        Summary metric for each seed (never treat edges/times as independent).
    n_replicates : int
        Number of bootstrap iterations (e.g. 10000).
    rng_seed : int
        Fixed RNG seed from design.json.

    Returns
    -------
    mean_val, std_val, ci_lower, ci_upper
    """
    arr = np.asarray(values, dtype=np.float64)
    n = len(arr)
    if n == 0 or np.all(np.isnan(arr)):
        return np.nan, np.nan, np.nan, np.nan

    mean_val = float(np.nanmean(arr))
    std_val = float(np.nanstd(arr, ddof=1)) if n > 1 else 0.0

    rng = np.random.default_rng(rng_seed)
    boot_indices = rng.integers(0, n, size=(n_replicates, n))
    boot_means = np.nanmean(arr[boot_indices], axis=1)

    ci_lower = float(np.percentile(boot_means, 2.5))
    ci_upper = float(np.percentile(boot_means, 97.5))
    return mean_val, std_val, ci_lower, ci_upper


def run_analysis_cli(result_dir_arg: str) -> None:
    """CLI execution for analyze.py."""
    result_dir = Path(result_dir_arg).resolve()
    if not result_dir.exists():
        sys.exit(f"Error: RESULT_DIR does not exist: {result_dir}")

    design_path = result_dir / "design.json"
    if not design_path.exists():
        sys.exit(f"Error: design.json not found in {result_dir}")

    with open(design_path, "r", encoding="utf-8") as f:
        design = json.load(f)

    dt = float(design["dt"])
    max_age = float(design["max_age"])
    max_variance = float(design["max_variance"])
    durations = [int(d) for d in design["durations"]]
    pre_steps = int(design.get("pre_steps", 200))
    recovery_steps = int(design.get("recovery_steps", 400))
    bootstrap_seed = int(design.get("bootstrap_seed", 2026110200))
    bootstrap_replicates = int(design.get("bootstrap_replicates", 10000))
    cases = design["cases"]

    replay_dir = result_dir / "raw" / "replay"
    if not replay_dir.exists():
        sys.exit(f"Error: Replay directory does not exist: {replay_dir}")

    analysis_dir = result_dir / "analysis"
    if analysis_dir.exists():
        raise FileExistsError("Analysis already exists; use a new experiment directory.")
    analysis_dir.mkdir(parents=True, exist_ok=True)
    figures_dir = analysis_dir / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)

    print(f"=== P2 Analyzer: Reading npz files from {replay_dir} ===")
    print(f"Output analysis directory: {analysis_dir}")

    # Discover and index npz files
    npz_data_store: Dict[Tuple[str, int], Dict[str, Any]] = {}
    missing_files = []
    with open(replay_dir/'manifest.json', encoding='utf-8') as f:
        manifest = json.load(f)
    if manifest['errors'] or manifest['completed'] != len(cases)*len(durations):
        raise ValueError("Cannot analyze an incomplete batch")
    expected_hashes = {r['npz_file']:r['npz_sha256'] for r in manifest['runs']}

    for case in cases:
        case_id = str(case["case_id"])
        for D in durations:
            npz_name = f"case_{case_id}_D_{D}.npz"
            npz_path = replay_dir / npz_name
            if not npz_path.exists():
                missing_files.append(npz_name)
                continue

            if compute_sha256(npz_path) != expected_hashes[npz_name]:
                raise ValueError(f"Raw SHA mismatch: {npz_name}")
            with np.load(npz_path, allow_pickle=False) as npz:
                loaded = {k: npz[k] for k in npz.files}
                loaded["_sha256"] = compute_sha256(npz_path)
                loaded["_path"] = str(npz_path)
                npz_data_store[(case_id, D)] = loaded

    if missing_files:
        raise FileNotFoundError(f"Missing replay files: {missing_files}")
        print(f"Warning: {len(missing_files)} expected npz files are missing:")
        for mf in missing_files[:5]:
            print(f"  - {mf}")
        if len(missing_files) > 5:
            print(f"  ... and {len(missing_files) - 5} more.")

    total_loaded = len(npz_data_store)
    print(f"Loaded {total_loaded} replay files for analysis.")
    if total_loaded == 0:
        sys.exit("Error: No replay files available to analyze.")

    methods = ["Hold", "KF-CV", "KF-GM"]
    method_keys = {"Hold": "Hold", "KF-CV": "KF_CV", "KF-GM": "KF_GM"}

    summary_rows: List[Dict[str, Any]] = []
    age_metric_rows: List[Dict[str, Any]] = []
    contrast_rows: List[Dict[str, Any]] = []

    # Calibration pooling for KF-CV and KF-GM
    calibration_pools: Dict[str, Dict[str, Dict[str, List[float]]]] = {
        m: {"occ": {"P": [], "sq_err": []}, "pre": {"P": [], "sq_err": []}}
        for m in ["KF-CV", "KF-GM"]
    }

    for (case_id, D), data in sorted(npz_data_store.items()):
        seed = int(data["seed"].item()) if data["seed"].ndim == 0 else int(data["seed"][0])
        stratum = str(data["stratum"].item()) if data["stratum"].ndim == 0 else str(data["stratum"][0])
        truth = data["truth"]  # [W, E]
        missing_win = data["missing"]  # [W] bool
        W, E = truth.shape

        occ_start = pre_steps
        occ_end = pre_steps + D
        rec_end = occ_end + recovery_steps

        window_slices = {
            "pre": slice(0, occ_start),
            "occ": slice(occ_start, occ_end),
            "recovery": slice(occ_end, rec_end),
        }

        # Angular speed of ground truth in this case
        # Wrap diff between consecutive steps: Delta truth = wrap(truth[t] - truth[t-1])
        diff_truth = wrap_to_pi(np.diff(truth, axis=0))
        speed_truth = diff_truth / dt  # [W-1, E]

        # Occ-level RMSE storage for paired contrast computation
        occ_rmse_by_method: Dict[str, float] = {}

        for method in methods:
            m_key = method_keys[method]
            err = data[f"error_{m_key}"]  # [W, E]
            p11 = data[f"p11_{m_key}"]  # [W, E]
            nis = data[f"nis_{m_key}"]  # [W, E]
            age = data[f"age_{m_key}"]  # [W, E]
            active = data[f"active_{m_key}"]  # [W, E] bool
            act_age = data[f"active_by_age_{m_key}"]  # [W, E] bool
            act_p = data[f"active_by_P_{m_key}"]  # [W, E] bool

            is_hold = (method == "Hold")

            # 1. Summary rows for each window: pre, occ, recovery
            for win_name, win_sl in window_slices.items():
                err_w = err[win_sl]  # [T_w, E]
                n_samples = err_w.size

                circ_mae = float(np.mean(np.abs(err_w)))
                circ_rmse = float(np.sqrt(np.mean(err_w**2)))

                if win_name == "occ":
                    occ_rmse_by_method[method] = circ_rmse

                if is_hold:
                    mean_p = np.nan
                    mse_over_mean_p = np.nan
                    approx95_cov = np.nan
                    wide_int_frac = np.nan
                    mean_nis = np.nan
                else:
                    p11_w = p11[win_sl]
                    mean_p = float(np.mean(p11_w))
                    mse_over_mean_p = float((circ_rmse**2) / mean_p) if mean_p > 0 else np.nan

                    # Approx 95% coverage: |err| <= 1.96 * sqrt(P)
                    cov_mask = np.abs(err_w) <= 1.96 * np.sqrt(np.maximum(0.0, p11_w))
                    approx95_cov = float(np.mean(cov_mask))

                    # Wide interval fraction: 1.96 * sqrt(P) >= pi
                    wide_mask = 1.96 * np.sqrt(np.maximum(0.0, p11_w)) >= np.pi
                    wide_int_frac = float(np.mean(wide_mask))

                    # Observed fresh mean NIS (only defined on fresh steps)
                    nis_w = nis[win_sl]
                    fin_nis = nis_w[np.isfinite(nis_w)]
                    mean_nis = float(np.mean(fin_nis)) if len(fin_nis) > 0 else np.nan

                    # Accumulate calibration data
                    if win_name in ["occ", "pre"]:
                        calibration_pools[method][win_name]["P"].extend(p11_w.ravel().tolist())
                        calibration_pools[method][win_name]["sq_err"].extend((err_w**2).ravel().tolist())

                # Gating & cutoff metrics over this window
                age_w = age[win_sl]
                cutoff_age_mask = age_w > max_age
                cutoff_by_age_frac = float(np.mean(cutoff_age_mask))

                if is_hold:
                    cutoff_by_p_frac = np.nan
                    cutoff_union_mask = cutoff_age_mask
                    cutoff_union_frac = cutoff_by_age_frac
                else:
                    p11_w = p11[win_sl]
                    cutoff_p_mask = p11_w > max_variance
                    cutoff_by_p_frac = float(np.mean(cutoff_p_mask))
                    cutoff_union_mask = cutoff_age_mask | cutoff_p_mask
                    cutoff_union_frac = float(np.mean(cutoff_union_mask))

                # First cutoff step during occlusion
                first_cutoff_step = np.nan
                if win_name == "occ":
                    # Check step-by-step 1..D
                    for step_idx in range(D):
                        if np.any(cutoff_union_mask[step_idx]):
                            first_cutoff_step = float(step_idx + 1)  # 1-based age step
                            break

                # Active-only RMSE
                act_w = active[win_sl]
                act_count = int(np.sum(act_w))
                if act_count > 0:
                    act_rmse = float(np.sqrt(np.mean(err_w[act_w] ** 2)))
                else:
                    act_rmse = np.nan

                # Truth angular speed metrics for this window
                # Account for diff slice (length W - 1)
                win_start_diff = max(0, win_sl.start - 1 if win_sl.start is not None else 0)
                win_end_diff = max(0, win_sl.stop - 1 if win_sl.stop is not None else W - 1)
                sp_w = speed_truth[win_start_diff:win_end_diff]
                diff_w = diff_truth[win_start_diff:win_end_diff]

                speed_rms = float(np.sqrt(np.mean(sp_w**2))) if sp_w.size > 0 else np.nan
                speed_max_abs = float(np.max(np.abs(sp_w))) if sp_w.size > 0 else np.nan
                near_nyquist_frac = float(np.mean(np.abs(diff_w) >= 0.9 * np.pi)) if diff_w.size > 0 else np.nan

                summary_rows.append({
                    "case_id": case_id,
                    "seed": seed,
                    "stratum": stratum,
                    "D": D,
                    "method": method,
                    "window": win_name,
                    "n_samples": n_samples,
                    "circular_mae": circ_mae,
                    "circular_rmse": circ_rmse,
                    "mean_P": mean_p,
                    "mse_over_mean_P": mse_over_mean_p,
                    "approx95coverage": approx95_cov,
                    "wide_interval_fraction": wide_int_frac,
                    "observed_fresh_meanNIS": mean_nis,
                    "cutoff_by_age_fraction": cutoff_by_age_frac,
                    "cutoff_by_P_fraction": cutoff_by_p_frac,
                    "cutoff_union_fraction": cutoff_union_frac,
                    "first_cutoff_age_step": first_cutoff_step,
                    "active_only_rmse": act_rmse,
                    "active_sample_count": act_count,
                    "truth_angular_speed_RMS": speed_rms,
                    "truth_angular_speed_max_abs": speed_max_abs,
                    "near_nyquist_fraction": near_nyquist_frac,
                })

            # 2. Per-age metrics (1..D) during occlusion
            occ_err = err[window_slices["occ"]]  # [D, E]
            occ_p11 = p11[window_slices["occ"]]  # [D, E]
            occ_age = age[window_slices["occ"]]  # [D, E]
            occ_act = active[window_slices["occ"]]  # [D, E]

            for d_step in range(1, D + 1):
                idx = d_step - 1
                e_step = occ_err[idx]
                d_rmse = float(np.sqrt(np.mean(e_step**2)))
                d_mae = float(np.mean(np.abs(e_step)))
                d_age_val = float(np.mean(occ_age[idx]))

                if is_hold:
                    d_p = np.nan
                    d_cov = np.nan
                    d_wide = np.nan
                    d_cut_p = np.nan
                else:
                    p_step = occ_p11[idx]
                    d_p = float(np.mean(p_step))
                    d_cov = float(np.mean(np.abs(e_step) <= 1.96 * np.sqrt(np.maximum(0.0, p_step))))
                    d_wide = float(np.mean(1.96 * np.sqrt(np.maximum(0.0, p_step)) >= np.pi))
                    d_cut_p = float(np.mean(p_step > max_variance))

                d_cut_age = float(np.mean(occ_age[idx] > max_age))
                d_cut_union = d_cut_age if is_hold else float(np.mean((occ_age[idx] > max_age) | (occ_p11[idx] > max_variance)))

                act_step = occ_act[idx]
                act_cnt = int(np.sum(act_step))
                act_r = float(np.sqrt(np.mean(e_step[act_step] ** 2))) if act_cnt > 0 else np.nan

                age_metric_rows.append({
                    "case_id": case_id,
                    "seed": seed,
                    "stratum": stratum,
                    "D": D,
                    "method": method,
                    "age_step": d_step,
                    "age_time": round(d_step * dt, 4),
                    "n_samples": E,
                    "circular_rmse": d_rmse,
                    "circular_mae": d_mae,
                    "mean_p11": d_p,
                    "approx95coverage": d_cov,
                    "wide_interval_fraction": d_wide,
                    "cutoff_by_age_fraction": d_cut_age,
                    "cutoff_by_P_fraction": d_cut_p,
                    "cutoff_union_fraction": d_cut_union,
                    "active_sample_count": act_cnt,
                    "active_only_rmse": act_r,
                })

        # 3. Paired contrasts for this case and D on occurrence window
        if all(m in occ_rmse_by_method for m in methods):
            rmse_h = occ_rmse_by_method["Hold"]
            rmse_cv = occ_rmse_by_method["KF-CV"]
            rmse_gm = occ_rmse_by_method["KF-GM"]

            contrast_rows.append({
                "case_id": case_id,
                "seed": seed,
                "stratum": stratum,
                "D": D,
                "rmse_Hold": rmse_h,
                "rmse_KF_CV": rmse_cv,
                "rmse_KF_GM": rmse_gm,
                "diff_KF_CV_minus_Hold": rmse_cv - rmse_h,
                "diff_KF_GM_minus_Hold": rmse_gm - rmse_h,
                "diff_KF_GM_minus_KF_CV": rmse_gm - rmse_cv,
            })

    # Convert to DataFrames
    df_summary = pd.DataFrame(summary_rows)
    df_age_metrics = pd.DataFrame(age_metric_rows)
    df_contrasts = pd.DataFrame(contrast_rows)

    # Save summary.csv, age_metrics.csv, contrasts.csv
    safe_save_csv(df_summary, analysis_dir / "summary.csv")
    safe_save_csv(df_age_metrics, analysis_dir / "age_metrics.csv")
    safe_save_csv(df_contrasts, analysis_dir / "contrasts.csv")

    # 4. Grouped analysis across seeds with bootstrap confidence intervals
    # "grouped per-stratum,D across5seeds means with bootstrap seed resampling intervals,
    # never treat edges/times as independent. 10000 bootstrap, fixed RNG seed from design"
    grouped_rows: List[Dict[str, Any]] = []

    metrics_to_group = [
        "rmse_Hold",
        "rmse_KF_CV",
        "rmse_KF_GM",
        "diff_KF_CV_minus_Hold",
        "diff_KF_GM_minus_Hold",
        "diff_KF_GM_minus_KF_CV",
    ]

    for stratum in sorted(df_contrasts["stratum"].unique()):
        for D in durations:
            sub = df_contrasts[(df_contrasts["stratum"] == stratum) & (df_contrasts["D"] == D)]
            if len(sub) == 0:
                continue

            for metric in metrics_to_group:
                seed_vals = sub[metric].to_numpy()
                mean_v, std_v, ci_lo, ci_hi = run_bootstrap_seeds(
                    seed_vals,
                    n_replicates=bootstrap_replicates,
                    rng_seed=bootstrap_seed,
                )
                grouped_rows.append({
                    "stratum": stratum,
                    "D": D,
                    "metric": metric,
                    "n_seeds": len(seed_vals),
                    "mean": mean_v,
                    "std": std_v,
                    "ci_lower": ci_lo,
                    "ci_upper": ci_hi,
                    "bootstrap_seed": bootstrap_seed,
                    "bootstrap_replicates": bootstrap_replicates,
                })

    df_grouped = pd.DataFrame(grouped_rows)
    safe_save_csv(df_grouped, analysis_dir / "grouped.csv")

    # 5. Calibration analysis CSV
    calib_rows: List[Dict[str, Any]] = []
    for method in ["KF-CV", "KF-GM"]:
        for win_name in ["occ", "pre"]:
            pool_p = np.array(calibration_pools[method][win_name]["P"])
            pool_sq = np.array(calibration_pools[method][win_name]["sq_err"])
            if len(pool_p) == 0:
                continue

            # Quantile bins (10 bins)
            try:
                quantiles = np.linspace(0.0, 1.0, 11)
                bins = np.percentile(pool_p, quantiles * 100.0)
                bins[0] -= 1e-12
                bins[-1] += 1e-12

                for b_idx in range(len(bins) - 1):
                    lo, hi = bins[b_idx], bins[b_idx + 1]
                    m_bin = (pool_p >= lo) & (pool_p < hi)
                    if np.sum(m_bin) == 0:
                        continue
                    b_p = pool_p[m_bin]
                    b_sq = pool_sq[m_bin]
                    cov_frac = np.mean(np.sqrt(b_sq) <= 1.96 * np.sqrt(np.maximum(0.0, b_p)))

                    calib_rows.append({
                        "method": method,
                        "window": win_name,
                        "bin_index": b_idx + 1,
                        "p11_min": float(np.min(b_p)),
                        "p11_max": float(np.max(b_p)),
                        "n_samples": int(np.sum(m_bin)),
                        "mean_p11": float(np.mean(b_p)),
                        "mean_squared_error": float(np.mean(b_sq)),
                        "approx95coverage": float(cov_frac),
                    })
            except Exception as e:
                print(f"Warning: Calibration binning failed for {method} {win_name}: {e}")

    df_calib = pd.DataFrame(calib_rows)
    safe_save_csv(df_calib, analysis_dir / "calibration.csv")

    # =========================================================================
    # Diagnostic Figures Generation (Strict English Labels)
    # =========================================================================

    # Figure 1: RMSE(age) per stratum / D with across-seed intervals
    try:
        fig1, axes1 = plt.subplots(3, 2, figsize=(12, 10), sharex="col", sharey="row")
        strata = ["S", "T", "V"]
        stratum_titles = {
            "S": "Stratum S (Slow Bearing)",
            "T": "Stratum T (Controlled Turn)",
            "V": "Stratum V (Fast Bearing)",
        }
        method_colors = {"Hold": "#7f7f7f", "KF-CV": "#1f77b4", "KF-GM": "#ff7f0e"}
        method_styles = {"Hold": "--", "KF-CV": "-", "KF-GM": "-."}

        for r_idx, strat in enumerate(strata):
            for c_idx, D in enumerate(durations):
                ax = axes1[r_idx, c_idx]
                sub_am = df_age_metrics[(df_age_metrics["stratum"] == strat) & (df_age_metrics["D"] == D)]

                for m in methods:
                    m_sub = sub_am[sub_am["method"] == m]
                    if len(m_sub) == 0:
                        continue

                    # Group by age_step across the seeds
                    agg = m_sub.groupby("age_step")["circular_rmse"].agg(["mean", "std", "count"]).reset_index()
                    steps = agg["age_step"].to_numpy()
                    means = agg["mean"].to_numpy()
                    stds = agg["std"].fillna(0.0).to_numpy()

                    ax.plot(
                        steps,
                        means,
                        label=m,
                        color=method_colors[m],
                        linestyle=method_styles[m],
                        linewidth=1.8,
                    )
                    ax.fill_between(
                        steps,
                        np.maximum(0.0, means - stds),
                        means + stds,
                        color=method_colors[m],
                        alpha=0.15,
                    )

                if r_idx == 0:
                    ax.set_title(f"D = {D} steps ({D * dt:.1f} model time)", fontsize=11, fontweight="bold")
                if c_idx == 0:
                    ax.set_ylabel(f"{stratum_titles[strat]}\nCircular RMSE (rad)", fontsize=10)
                if r_idx == 2:
                    ax.set_xlabel("Occlusion age (steps; dt = 0.3 model time)", fontsize=10)

                ax.grid(True, linestyle=":", alpha=0.6)
                if r_idx == 0 and c_idx == 0:
                    ax.legend(loc="upper left", frameon=True, fontsize=9)

        fig1.suptitle("RMSE vs Age (lines: 5-seed means; bands: +/-1 seed SD)", fontsize=13, fontweight="bold", y=0.99)
        safe_save_figure(fig1, figures_dir / "fig1_rmse_vs_age.png")
    except Exception as e:
        print(f"Error generating Figure 1: {e}")

    # Figure 2: P vs Squared Error (Calibration Curve)
    try:
        fig2, (ax2a, ax2b) = plt.subplots(1, 2, figsize=(12, 5))
        for ax, m, title in zip([ax2a, ax2b], ["KF-CV", "KF-GM"], ["KF-CV Calibration", "KF-GM Calibration"]):
            sub_c = df_calib[(df_calib["method"] == m) & (df_calib["window"] == "occ")]
            if len(sub_c) > 0:
                ax.plot(
                    sub_c["mean_p11"],
                    sub_c["mean_squared_error"],
                    marker="o",
                    color=method_colors[m],
                    label="Binned Empirical MSE",
                    linewidth=1.8,
                )
                min_v = min(sub_c["mean_p11"].min(), sub_c["mean_squared_error"].min())
                max_v = max(sub_c["mean_p11"].max(), sub_c["mean_squared_error"].max())
                ref_line = np.geomspace(max(1e-6, min_v), max_v, 50)
                ax.plot(ref_line, ref_line, "k--", alpha=0.7, label="Ideal (MSE = P)")
                ax.set_xscale("log")
                ax.set_yscale("log")

            ax.set_title(title, fontsize=11, fontweight="bold")
            ax.set_xlabel("Mean Predicted Variance P11 (rad^2)", fontsize=10)
            ax.set_ylabel("Empirical Mean Squared Error (rad^2)", fontsize=10)
            ax.grid(True, which="both", linestyle=":", alpha=0.5)
            ax.legend(loc="upper left", frameon=True, fontsize=9)

        fig2.suptitle("Descriptive pooled calibration (reused cases; no independent sample count)", fontsize=11, fontweight="bold")
        safe_save_figure(fig2, figures_dir / "fig2_p_vs_squared_error.png")
    except Exception as e:
        print(f"Error generating Figure 2: {e}")

    # Figure 3: Approx 95% Coverage vs Age
    try:
        fig3, (ax3a, ax3b) = plt.subplots(1, 2, figsize=(12, 4.5), sharey=True)
        for ax, D in zip([ax3a, ax3b], durations):
            sub_am = df_age_metrics[df_age_metrics["D"] == D]
            for m in ["KF-CV", "KF-GM"]:
                m_sub = sub_am[sub_am["method"] == m]
                if len(m_sub) == 0:
                    continue
                agg = m_sub.groupby("age_step")["approx95coverage"].mean().reset_index()
                ax.plot(
                    agg["age_step"],
                    agg["approx95coverage"] * 100.0,
                    label=m,
                    color=method_colors[m],
                    linewidth=1.8,
                )

            ax.axhline(95.0, color="r", linestyle="--", linewidth=1.2, label="Nominal 95% Target")
            ax.set_title(f"Occlusion Duration D = {D} steps", fontsize=11, fontweight="bold")
            ax.set_xlabel("Age (steps; dt = 0.3 model time)", fontsize=10)
            ax.set_ylabel("Empirical Coverage (%)" if ax == ax3a else "", fontsize=10)
            ax.set_ylim(-5.0, 105.0)
            ax.grid(True, linestyle=":", alpha=0.6)
            ax.legend(loc="lower left", frameon=True, fontsize=9)

        fig3.suptitle("Approx. Gaussian coverage (descriptive mean over S/T/V; 5 seeds)", fontsize=11, fontweight="bold")
        safe_save_figure(fig3, figures_dir / "fig3_coverage_vs_age.png")
    except Exception as e:
        print(f"Error generating Figure 3: {e}")

    # Figure 4: Fresh Normalized Innovation Squared (NIS) Distribution
    try:
        fig4, (ax4a, ax4b) = plt.subplots(1, 2, figsize=(12, 4.5), sharey=True)
        sub_fresh = df_summary[(df_summary["window"] == "pre") & (df_summary["method"] != "Hold") & (df_summary["D"] == durations[0])]

        for ax, m in zip([ax4a, ax4b], ["KF-CV", "KF-GM"]):
            m_data = sub_fresh[sub_fresh["method"] == m]["observed_fresh_meanNIS"].dropna()
            if len(m_data) > 0:
                ax.hist(m_data, bins=10, color=method_colors[m], alpha=0.7, edgecolor="k")
                ax.axvline(1.0, color="red", linestyle="--", linewidth=1.5, label="Expected Mean (NIS = 1)")
                ax.axvline(m_data.mean(), color="black", linestyle="-", linewidth=1.5, label=f"Mean: {m_data.mean():.2f}")

            ax.set_title(f"{m} Fresh NIS (Pre-Occlusion Window)", fontsize=11, fontweight="bold")
            ax.set_xlabel("Normalized Innovation Squared (NIS)", fontsize=10)
            ax.set_ylabel("Case summaries (3 strata x 5 reused seeds)" if ax == ax4a else "", fontsize=9)
            ax.grid(True, linestyle=":", alpha=0.5)
            ax.legend(loc="upper right", frameon=True, fontsize=9)

        fig4.suptitle("Innovation Consistency on Fresh Prefix Tracks", fontsize=13, fontweight="bold")
        safe_save_figure(fig4, figures_dir / "fig4_fresh_nis.png")
    except Exception as e:
        print(f"Error generating Figure 4: {e}")

    # Figure 5: Cutoff Fraction vs Age Step (D = 120)
    try:
        fig5, ax5 = plt.subplots(figsize=(8, 5))
        sub_d120 = df_age_metrics[df_age_metrics["D"] == 120]

        for m, ls in zip(["Hold", "KF-CV", "KF-GM"], ["--", "-", "-."]):
            m_sub = sub_d120[sub_d120["method"] == m]
            if len(m_sub) == 0:
                continue
            agg = m_sub.groupby("age_step")["cutoff_union_fraction"].mean().reset_index()
            ax5.plot(
                agg["age_step"],
                agg["cutoff_union_fraction"] * 100.0,
                label=f"{m} Union Cutoff",
                color=method_colors[m],
                linestyle=ls,
                linewidth=1.8,
            )

        ax5.set_title("Cutoff Fraction vs Occlusion Age (D = 120 steps)", fontsize=11, fontweight="bold")
        ax5.set_xlabel("Age (steps; max_age = 120 model time = 400 steps)", fontsize=10)
        ax5.set_ylabel("Cutoff Fraction (%)", fontsize=10)
        ax5.set_ylim(-5.0, 105.0)
        ax5.grid(True, linestyle=":", alpha=0.6)
        ax5.legend(loc="upper left", frameon=True, fontsize=9)
        safe_save_figure(fig5, figures_dir / "fig5_cutoff_fraction_vs_age.png")
    except Exception as e:
        print(f"Error generating Figure 5: {e}")

    # Figure 6: Paired RMSE Differences per Seed
    try:
        fig6, axes6 = plt.subplots(1, 2, figsize=(13, 5), sharey=True)
        contrast_keys = [
            ("diff_KF_CV_minus_Hold", "KF-CV - Hold", "#1f77b4"),
            ("diff_KF_GM_minus_Hold", "KF-GM - Hold", "#ff7f0e"),
            ("diff_KF_GM_minus_KF_CV", "KF-GM - KF-CV", "#2ca02c"),
        ]

        for ax, D in zip(axes6, durations):
            sub_c = df_contrasts[df_contrasts["D"] == D]
            x_ticks = []
            x_labels = []

            for i, (ck, label_c, col) in enumerate(contrast_keys):
                vals = sub_c[ck].to_numpy()
                x_pos = np.arange(len(vals)) + i * (len(vals) + 2)
                ax.scatter(x_pos, vals, color=col, label=label_c if ax == axes6[0] else "", s=40, zorder=3)
                mean_c = np.mean(vals)
                ax.hlines(mean_c, x_pos[0] - 0.4, x_pos[-1] + 0.4, colors=col, linestyles="-", linewidth=2.0)
                x_ticks.append(np.mean(x_pos))
                x_labels.append(label_c)

            ax.axhline(0.0, color="gray", linestyle="--", linewidth=1.2, zorder=1)
            ax.set_xticks(x_ticks)
            ax.set_xticklabels(x_labels, fontsize=9)
            ax.set_title(f"Duration D = {D} steps (Paired Seed Contrasts)", fontsize=11, fontweight="bold")
            ax.set_ylabel("RMSE Difference (rad)\n[Negative = Improvement]" if ax == axes6[0] else "", fontsize=10)
            ax.grid(True, linestyle=":", alpha=0.5)

        axes6[0].legend(loc="upper right", frameon=True, fontsize=9)
        fig6.suptitle("Descriptive S/T/V case contrasts (15 cases; 5 independent seeds)", fontsize=12, fontweight="bold")
        safe_save_figure(fig6, figures_dir / "fig6_paired_rmse_diff_per_seed.png")
    except Exception as e:
        print(f"Error generating Figure 6: {e}")

    # Figure 7: Fixed Representative Case Error Trajectories (First seed of each stratum, D=120)
    try:
        fig7, axes7 = plt.subplots(3, 1, figsize=(13, 9), sharex=True)
        for r_idx, strat in enumerate(strata):
            ax = axes7[r_idx]
            sub_cases = [c for c in cases if c["stratum"] == strat]
            if not sub_cases:
                continue
            rep_case = sub_cases[0]
            rep_id = str(rep_case["case_id"])

            rep_data = npz_data_store.get((rep_id, 120))
            if rep_data is None:
                continue

            # Plot mean error trajectory across 140 edges over the entire window
            steps_1based = rep_data["window_steps"]
            t_axis = (steps_1based - steps_1based[0]) * dt

            err_h = np.mean(np.abs(rep_data["error_Hold"]), axis=1)
            err_cv = np.mean(np.abs(rep_data["error_KF_CV"]), axis=1)
            err_gm = np.mean(np.abs(rep_data["error_KF_GM"]), axis=1)

            ax.plot(t_axis, err_h, label="Hold", color=method_colors["Hold"], linestyle="--", linewidth=1.4)
            ax.plot(t_axis, err_cv, label="KF-CV", color=method_colors["KF-CV"], linestyle="-", linewidth=1.5)
            ax.plot(t_axis, err_gm, label="KF-GM", color=method_colors["KF-GM"], linestyle="-.", linewidth=1.5)

            # Mark occlusion window
            t_occ_start = pre_steps * dt
            t_occ_end = (pre_steps + 120) * dt
            ax.axvspan(t_occ_start, t_occ_end, color="gray", alpha=0.15, label="Occlusion (D=120)")
            ax.axvline(t_occ_start, color="gray", linestyle=":", linewidth=1.0)
            ax.axvline(t_occ_end, color="gray", linestyle=":", linewidth=1.0)

            ax.set_title(
                f"{stratum_titles[strat]}: Representative Case {rep_id} (Seed {rep_case['seed']})",
                fontsize=10,
                fontweight="bold",
            )
            ax.set_ylabel("Mean |Error| (rad)", fontsize=9)
            ax.grid(True, linestyle=":", alpha=0.5)
            if r_idx == 0:
                ax.legend(loc="upper left", frameon=True, fontsize=9)

        axes7[-1].set_xlabel("Elapsed window time (model units; dt = 0.3)", fontsize=10)
        fig7.suptitle("Error Trajectories across Pre, Occlusion, and Recovery Windows", fontsize=13, fontweight="bold", y=0.99)
        safe_save_figure(fig7, figures_dir / "fig7_representative_case_trajectories.png")
    except Exception as e:
        print(f"Error generating Figure 7: {e}")

    # =========================================================================
    # Write Audit JSON
    # =========================================================================
    if len(list(figures_dir.glob('*.png'))) != 7:
        raise RuntimeError("Figure generation incomplete; logs and partial outputs preserved")
    audit_data = {
        "analysis_timestamp_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "result_dir": str(result_dir),
        "design_parameters": {
            "dt": dt,
            "max_age": max_age,
            "max_variance": max_variance,
            "durations": durations,
            "pre_steps": pre_steps,
            "recovery_steps": recovery_steps,
            "bootstrap_seed": bootstrap_seed,
            "bootstrap_replicates": bootstrap_replicates,
        },
        "python_environment": {
            "python_version": sys.version,
            "numpy_version": np.__version__,
            "pandas_version": pd.__version__,
            "matplotlib_version": matplotlib.__version__,
            "scipy_version": scipy.__version__,
        },
        "analyzed_runs_count": total_loaded,
        "summary_rows_count": len(df_summary),
        "age_metric_rows_count": len(df_age_metrics),
        "contrast_rows_count": len(df_contrasts),
        "grouped_rows_count": len(df_grouped),
        "analyzed_files": [
            {"case_id": c_id, "D": D_val, "path": d["_path"], "sha256": d["_sha256"]}
            for (c_id, D_val), d in sorted(npz_data_store.items())
        ],
    }

    audit_path = analysis_dir / "audit.json"
    with open(audit_path, "w", encoding="utf-8") as f:
        json.dump(audit_data, f, indent=2)
    print(f"[WROTE AUDIT] {audit_path.name}")

    print("=== Analysis Complete ===")
    print(f"Summary rows: {len(df_summary)}")
    print(f"Age metric rows: {len(df_age_metrics)}")
    print(f"Contrast rows: {len(df_contrasts)}")
    print(f"Grouped rows: {len(df_grouped)}")
    print(f"Calibration rows: {len(df_calib)}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="P2 Open-Loop Replay Diagnostic Analysis")
    parser.add_argument("result_dir", help="Path to experiment result directory containing design.json and raw/replay")
    args = parser.parse_args()
    run_analysis_cli(args.result_dir)

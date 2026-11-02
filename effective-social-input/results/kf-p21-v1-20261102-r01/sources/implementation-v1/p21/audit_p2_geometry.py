"""experiments/kf-p21-20261102/audit_p2_geometry.py

G0 audit script for P2 baseline reproduction and post-hoc geometric stratification.
Validates original P2 raw data, inputs, and summary rows to < 1e-10.
Recomputes geometry without filter access using geometry.py.
Outputs compressed geometry events, seed-level error tables, NIS tails, age coverage,
variance calibration bins, and diagnostic figures.

Usage:
  python audit_p2_geometry.py <OLD_P2_RESULT_DIR> <NEW_G0_DIR>
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
from pathlib import Path
import sys
from typing import Any, Dict, List, Optional, Tuple, Union

import matplotlib
matplotlib.use("Agg")  # Headless non-interactive plotting
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.io

# Import geometry from same directory
sys.path.insert(0, str(Path(__file__).resolve().parent))
import geometry


def compute_sha256(filepath: Union[str, Path]) -> str:
    """Compute SHA256 hex digest of a file in streaming chunks."""
    h = hashlib.sha256()
    with open(filepath, "rb") as f:
        while chunk := f.read(65536):
            h.update(chunk)
    return h.hexdigest()


def wrap_to_pi(angle: np.ndarray) -> np.ndarray:
    """Wrap angle differences to [-pi, pi)."""
    return np.arctan2(np.sin(angle), np.cos(angle))


def run_g0_audit(old_dir_arg: str, new_dir_arg: str) -> None:
    """Execute complete G0 baseline audit and geometric stratification."""
    old_p2_dir = Path(old_dir_arg).resolve()
    new_g0_dir = Path(new_dir_arg).resolve()

    if not old_p2_dir.exists():
        sys.exit(f"Error: OLD_P2_RESULT_DIR does not exist: {old_p2_dir}")

    # Guard: strictly refuse to overwrite existing output directory
    if new_g0_dir.exists():
        sys.exit(f"Error: Refusing to overwrite existing NEW_G0_DIR: {new_g0_dir}")

    new_g0_dir.mkdir(parents=True, exist_ok=False)
    figures_dir = new_g0_dir / "figures"
    figures_dir.mkdir(parents=True, exist_ok=False)

    print(f"=== P2.1 G0 Baseline Audit & Geometric Decoupling ===")
    print(f"Old P2 Directory: {old_p2_dir}")
    print(f"New G0 Directory: {new_g0_dir}")

    # 1. Load P2 design, manifests, and summary
    design_path = old_p2_dir / "design.json"
    if not design_path.exists():
        sys.exit(f"Error: design.json not found in {old_p2_dir}")
    with open(design_path, "r", encoding="utf-8") as f:
        design = json.load(f)

    dt = float(design["dt"])
    durations = [int(d) for d in design["durations"]]
    pre_steps = int(design.get("pre_steps", 200))
    recovery_steps = int(design.get("recovery_steps", 400))
    cases = design["cases"]

    input_manifest_path = old_p2_dir / "input_manifest.json"
    if not input_manifest_path.exists():
        sys.exit(f"Error: input_manifest.json not found in {old_p2_dir}")
    with open(input_manifest_path, "r", encoding="utf-8") as f:
        input_manifest = json.load(f)

    replay_manifest_path = old_p2_dir / "raw" / "replay" / "manifest.json"
    if not replay_manifest_path.exists():
        sys.exit(f"Error: replay manifest.json not found in {old_p2_dir}")
    with open(replay_manifest_path, "r", encoding="utf-8") as f:
        replay_manifest = json.load(f)

    summary_path = old_p2_dir / "analysis" / "summary.csv"
    if not summary_path.exists():
        sys.exit(f"Error: summary.csv not found in {old_p2_dir}")
    df_old_summary = pd.read_csv(summary_path)

    # 2. Verify SHA256 hashes of inputs and replays
    sha_audit: Dict[str, Any] = {"inputs": {}, "replays": {}}
    input_sha_ok = True
    for item in input_manifest.get("inputs", []):
        src_rel = item["source"]
        exp_sha = item["sha256"]
        fpath = old_p2_dir / src_rel
        if not fpath.exists():
            sha_audit["inputs"][src_rel] = {"status": "MISSING"}
            input_sha_ok = False
            continue
        actual_sha = compute_sha256(fpath)
        matched = (actual_sha == exp_sha)
        sha_audit["inputs"][src_rel] = {
            "expected": exp_sha,
            "actual": actual_sha,
            "matched": matched,
        }
        if not matched:
            input_sha_ok = False

    replay_sha_ok = True
    for item in replay_manifest.get("runs", []):
        r_file = item["npz_file"]
        exp_sha = item["npz_sha256"]
        fpath = old_p2_dir / "raw" / "replay" / r_file
        if not fpath.exists():
            sha_audit["replays"][r_file] = {"status": "MISSING"}
            replay_sha_ok = False
            continue
        actual_sha = compute_sha256(fpath)
        matched = (actual_sha == exp_sha)
        sha_audit["replays"][r_file] = {
            "expected": exp_sha,
            "actual": actual_sha,
            "matched": matched,
        }
        if not matched:
            replay_sha_ok = False

    # 3. Serial recovery of all 270 rows from summary.csv
    methods = ["Hold", "KF-CV", "KF-GM"]
    method_keys = {"Hold": "Hold", "KF-CV": "KF_CV", "KF-GM": "KF_GM"}
    window_slices = {
        "pre": lambda D: slice(0, pre_steps),
        "occ": lambda D: slice(pre_steps, pre_steps + D),
        "recovery": lambda D: slice(pre_steps + D, pre_steps + D + recovery_steps),
    }

    max_rmse_diff = 0.0
    max_nis_diff = 0.0
    rows_audited = 0
    all_rows_pass = True

    # Pre-index old summary for fast lookup
    old_summary_dict: Dict[Tuple[str, int, str, str], Dict[str, float]] = {}
    for _, row in df_old_summary.iterrows():
        key = (str(row["case_id"]), int(row["D"]), str(row["method"]), str(row["window"]))
        old_summary_dict[key] = {
            "circular_rmse": float(row["circular_rmse"]),
            "observed_fresh_meanNIS": float(row["observed_fresh_meanNIS"])
            if pd.notnull(row["observed_fresh_meanNIS"])
            else np.nan,
        }

    # Replay data cache during processing
    # To keep memory bounded, we process case-by-case serially
    geometry_events_path = new_g0_dir / "geometry_events.csv.gz"

    # Data structures for statistical aggregations
    counts_rows: List[Dict[str, Any]] = []
    error_seed_rows: List[Dict[str, Any]] = []
    nis_tail_rows: List[Dict[str, Any]] = []
    coverage_age_rows: List[Dict[str, Any]] = []
    p_bins_rows: List[Dict[str, Any]] = []

    # Clean rate storage for pilot tercile derivation
    pilot_g_clean_rates: List[float] = []

    # Calibration pooling for P bins
    p_pools: Dict[str, List[Tuple[float, float]]] = {
        "KF-CV": [],
        "KF-GM": [],
    }

    # Storage for representative edges for Figure G0-1
    fig_g0_1_data: Dict[str, Any] = {}

    # Storage for error and NIS samples for Figure G0-2
    fig_g0_2_error_samples: Dict[str, List[float]] = {
        label: [] for label in geometry.WINDOW_LABELS
    }
    fig_g0_2_nis_samples: Dict[str, List[float]] = {
        fl: [] for fl in geometry.FRESH_LABELS
    }

    # Open gzip output stream for geometry_events.csv.gz
    with gzip.open(geometry_events_path, "wt", encoding="utf-8") as gz_file:
        gz_header = [
            "case_id",
            "seed",
            "stratum",
            "D",
            "t0",
            "source",
            "step_rel",
            "step_global",
            "window_period",
            "edge",
            "obs_id",
            "tgt_id",
            "dx_min",
            "dy_min",
            "image_index_x",
            "image_index_y",
            "distance",
            "image_switch",
            "near_range",
            "near_half_box",
            "near_aliasing",
            "angular_rate_clean",
            "window_geom_label",
            "fresh_label",
        ]
        gz_file.write(",".join(gz_header) + "\n")

        # Serial loop over cases
        for case in cases:
            case_id = str(case["case_id"])
            seed = int(case["seed"])
            stratum = str(case["stratum"])
            source_rel = str(case["source"])
            t0 = int(case["t0"])

            source_path = old_p2_dir / source_rel
            if not source_path.exists():
                sys.exit(f"Error: Source MAT file missing: {source_path}")

            # Load full source MAT and recompute geometry using geometry.py
            mat_raw = scipy.io.loadmat(str(source_path), simplify_cells=True)
            g = geometry.geometry(mat_raw, arena_size=1000.0, dt=dt)
            fl_all = geometry.fresh_labels(g)

            T_full, E = g["distance"].shape
            obs_id_arr = np.asarray(mat_raw["observer_id"]).reshape(-1)
            tgt_id_arr = np.asarray(mat_raw["target_id"]).reshape(-1)

            for D in durations:
                npz_file = f"case_{case_id}_D_{D}.npz"
                npz_path = old_p2_dir / "raw" / "replay" / npz_file
                if not npz_path.exists():
                    sys.exit(f"Error: Replay NPZ missing: {npz_path}")

                with np.load(npz_path, allow_pickle=False) as npz:
                    truth = npz["truth"]  # [W, E]
                    missing = npz["missing"]  # [W] bool
                    W = truth.shape[0]

                    # Load method outputs
                    data_methods = {}
                    for m in methods:
                        m_k = method_keys[m]
                        data_methods[m] = {
                            "error": npz[f"error_{m_k}"],
                            "p11": npz[f"p11_{m_k}"],
                            "nis": npz[f"nis_{m_k}"],
                            "age": npz[f"age_{m_k}"],
                            "theta": npz[f"theta_{m_k}"],
                        }

                # 3a. Baseline audit verification for this (case_id, D)
                for m in methods:
                    for win_name in ["pre", "occ", "recovery"]:
                        sl = window_slices[win_name](D)
                        err_w = data_methods[m]["error"][sl]
                        rec_rmse = float(np.sqrt(np.mean(err_w**2)))

                        rec_nis = np.nan
                        if m != "Hold" and win_name in ["pre", "recovery"]:
                            nis_w = data_methods[m]["nis"][sl]
                            fin = nis_w[np.isfinite(nis_w)]
                            rec_nis = float(np.mean(fin)) if len(fin) > 0 else np.nan

                        old_val = old_summary_dict.get((case_id, D, m, win_name))
                        if old_val is None:
                            all_rows_pass = False
                            continue

                        diff_rmse = abs(rec_rmse - old_val["circular_rmse"])
                        max_rmse_diff = max(max_rmse_diff, diff_rmse)
                        if diff_rmse >= 1e-10:
                            all_rows_pass = False

                        if not np.isnan(rec_nis) and not np.isnan(old_val["observed_fresh_meanNIS"]):
                            diff_nis = abs(rec_nis - old_val["observed_fresh_meanNIS"])
                            max_nis_diff = max(max_nis_diff, diff_nis)
                            if diff_nis >= 1e-10:
                                all_rows_pass = False

                        rows_audited += 1

                # 3b. Window classification
                w_labels, w_rates = geometry.classify_window(g, t0=t0, D=D, pre_steps=pre_steps)
                # w_labels: [E], w_rates: [E]

                # Collect clean rates for pilot terciles
                for e in range(E):
                    if w_labels[e] == geometry.LABEL_CLEAN and np.isfinite(w_rates[e]):
                        pilot_g_clean_rates.append(float(w_rates[e]))

                # Counts per window label
                for lbl in geometry.WINDOW_LABELS:
                    count_lbl = int(np.sum(w_labels == lbl))
                    counts_rows.append({
                        "case_id": case_id,
                        "seed": seed,
                        "stratum": stratum,
                        "D": D,
                        "geometry_class": lbl,
                        "count_edges": count_lbl,
                        "total_edges": E,
                        "fraction": count_lbl / float(E),
                    })

                # Stream to geometry_events.csv.gz
                occ_start_global = (t0 - 1)
                win_start_global = occ_start_global - pre_steps

                d_min_win = g["d_min"][win_start_global : win_start_global + W]
                img_idx_win = g["image_index"][win_start_global : win_start_global + W]
                dist_win = g["distance"][win_start_global : win_start_global + W]
                switch_win = g["image_switch"][win_start_global : win_start_global + W]
                near_win = g["near_range"][win_start_global : win_start_global + W]
                half_box_win = g["near_half_box"][win_start_global : win_start_global + W]
                alias_win = g["near_aliasing"][win_start_global : win_start_global + W]
                rate_win = g["angular_rate_clean"][win_start_global : win_start_global + W]
                fresh_win = fl_all[win_start_global : win_start_global + W]

                for s_rel in range(W):
                    s_glob = win_start_global + s_rel
                    if s_rel < pre_steps:
                        period = "pre"
                    elif s_rel < pre_steps + D:
                        period = "occ"
                    else:
                        period = "recovery"

                    for e in range(E):
                        row_vals = [
                            case_id,
                            str(seed),
                            stratum,
                            str(D),
                            str(t0),
                            source_rel,
                            str(s_rel),
                            str(s_glob),
                            period,
                            str(e),
                            str(obs_id_arr[e]),
                            str(tgt_id_arr[e]),
                            f"{d_min_win[s_rel, e, 0]:.4f}",
                            f"{d_min_win[s_rel, e, 1]:.4f}",
                            f"{img_idx_win[s_rel, e, 0]:.1f}",
                            f"{img_idx_win[s_rel, e, 1]:.1f}",
                            f"{dist_win[s_rel, e]:.4f}",
                            "1" if switch_win[s_rel, e] else "0",
                            "1" if near_win[s_rel, e] else "0",
                            "1" if half_box_win[s_rel, e] else "0",
                            "1" if alias_win[s_rel, e] else "0",
                            f"{rate_win[s_rel, e]:.6f}" if np.isfinite(rate_win[s_rel, e]) else "nan",
                            str(w_labels[e]),
                            str(fresh_win[s_rel, e]),
                        ]
                        gz_file.write(",".join(row_vals) + "\n")

                # 3c. Error metrics stratified by window geometry in occlusion
                occ_sl = window_slices["occ"](D)
                for m in methods:
                    err_occ = data_methods[m]["error"][occ_sl]  # [D, E]

                    # All geometry
                    abs_err_all = np.abs(err_occ).ravel()
                    error_seed_rows.append({
                        "case_id": case_id,
                        "seed": seed,
                        "stratum": stratum,
                        "D": D,
                        "method": m,
                        "geometry_class": "all",
                        "n_samples": int(abs_err_all.size),
                        "circular_rmse": float(np.sqrt(np.mean(abs_err_all**2))),
                        "circular_mae": float(np.mean(abs_err_all)),
                        "p90": float(np.percentile(abs_err_all, 90)),
                        "p95": float(np.percentile(abs_err_all, 95)),
                        "p99": float(np.percentile(abs_err_all, 99)),
                        "max_abs_err": float(np.max(abs_err_all)),
                        "extreme_err_count": int(np.sum(abs_err_all > 1.0)),
                        "extreme_err_rate": float(np.mean(abs_err_all > 1.0)),
                    })

                    # Per geometry label
                    for lbl in geometry.WINDOW_LABELS:
                        mask_lbl = (w_labels == lbl)
                        if not np.any(mask_lbl):
                            continue
                        abs_err_lbl = np.abs(err_occ[:, mask_lbl]).ravel()
                        error_seed_rows.append({
                            "case_id": case_id,
                            "seed": seed,
                            "stratum": stratum,
                            "D": D,
                            "method": m,
                            "geometry_class": lbl,
                            "n_samples": int(abs_err_lbl.size),
                            "circular_rmse": float(np.sqrt(np.mean(abs_err_lbl**2))),
                            "circular_mae": float(np.mean(abs_err_lbl)),
                            "p90": float(np.percentile(abs_err_lbl, 90)),
                            "p95": float(np.percentile(abs_err_lbl, 95)),
                            "p99": float(np.percentile(abs_err_lbl, 99)),
                            "max_abs_err": float(np.max(abs_err_lbl)),
                            "extreme_err_count": int(np.sum(abs_err_lbl > 1.0)),
                            "extreme_err_rate": float(np.mean(abs_err_lbl > 1.0)),
                        })

                        # Collect samples for Figure G0-2 error tails (Hold and KF-CV)
                        if m == "KF-CV" and len(fig_g0_2_error_samples[lbl]) < 50000:
                            fig_g0_2_error_samples[lbl].extend(abs_err_lbl.tolist())

                # 3d. NIS stratification by fresh_labels in pre and recovery
                fresh_periods = np.zeros(W, dtype=bool)
                fresh_periods[:pre_steps] = True
                fresh_periods[pre_steps + D :] = True

                for m in ["KF-CV", "KF-GM"]:
                    nis_mat = data_methods[m]["nis"]  # [W, E]
                    fl_sub = fresh_win  # [W, E]

                    for fl_name in geometry.FRESH_LABELS:
                        mask_fl = (fl_sub == fl_name) & fresh_periods[:, None]
                        vals = nis_mat[mask_fl]
                        fin_vals = vals[np.isfinite(vals)]
                        if len(fin_vals) > 0:
                            nis_tail_rows.append({
                                "case_id": case_id,
                                "seed": seed,
                                "stratum": stratum,
                                "D": D,
                                "method": m,
                                "fresh_label": fl_name,
                                "n_samples": int(len(fin_vals)),
                                "mean_nis": float(np.mean(fin_vals)),
                                "median_nis": float(np.median(fin_vals)),
                                "p90_nis": float(np.percentile(fin_vals, 90)),
                                "p95_nis": float(np.percentile(fin_vals, 95)),
                                "p99_nis": float(np.percentile(fin_vals, 99)),
                            })
                            if m == "KF-CV" and len(fig_g0_2_nis_samples[fl_name]) < 50000:
                                fig_g0_2_nis_samples[fl_name].extend(fin_vals.tolist())

                # 3e. Age band calibration during occlusion
                # Bands: 1-20, 21-60, 61-100, 101-120
                age_bands = [
                    ("1-20", 1, 20),
                    ("21-60", 21, 60),
                    ("61-100", 61, 100),
                    ("101-120", 101, 120),
                ]

                for band_name, start_step, end_step in age_bands:
                    if start_step > D:
                        continue
                    act_end = min(D, end_step)
                    band_sl = slice(pre_steps + start_step - 1, pre_steps + act_end)

                    for m in ["KF-CV", "KF-GM"]:
                        e_b = data_methods[m]["error"][band_sl]  # [steps, E]
                        p_b = data_methods[m]["p11"][band_sl]  # [steps, E]

                        # Overall
                        cov = np.abs(e_b) <= 1.96 * np.sqrt(np.maximum(0.0, p_b))
                        wide = 1.96 * np.sqrt(np.maximum(0.0, p_b)) >= np.pi
                        sq_err = e_b**2

                        coverage_age_rows.append({
                            "case_id": case_id,
                            "seed": seed,
                            "stratum": stratum,
                            "D": D,
                            "method": m,
                            "age_band": band_name,
                            "geometry_class": "all",
                            "n_samples": int(e_b.size),
                            "mean_sq_err": float(np.mean(sq_err)),
                            "mean_P": float(np.mean(p_b)),
                            "ratio_sq_err_over_P": float(np.mean(sq_err) / np.mean(p_b))
                            if np.mean(p_b) > 0
                            else np.nan,
                            "approx95_coverage": float(np.mean(cov)),
                            "wide_interval_frac": float(np.mean(wide)),
                            "extreme_err_rate": float(np.mean(np.abs(e_b) > 1.0)),
                        })

                        # Collect samples for P pools
                        if m in p_pools:
                            p_flat = p_b.ravel()
                            sq_flat = sq_err.ravel()
                            for p_val, sq_val in zip(p_flat, sq_flat):
                                p_pools[m].append((float(p_val), float(sq_val)))

                # 3f. Capture Figure G0-1 data on first seed
                if seed == 2026110101 and D == 120 and "clean_edge" not in fig_g0_1_data:
                    clean_indices = np.where(w_labels == geometry.LABEL_CLEAN)[0]
                    switch_indices = np.where(w_labels == geometry.LABEL_SWITCH)[0]
                    if len(clean_indices) > 0 and len(switch_indices) > 0:
                        c_edge = int(clean_indices[0])
                        s_edge = int(switch_indices[0])
                        fig_g0_1_data["time_steps"] = (np.arange(W) - pre_steps) * dt
                        fig_g0_1_data["pre_steps"] = pre_steps
                        fig_g0_1_data["D"] = D
                        fig_g0_1_data["clean_edge"] = c_edge
                        fig_g0_1_data["switch_edge"] = s_edge
                        fig_g0_1_data["case_id"] = case_id
                        fig_g0_1_data["clean_truth"] = truth[:, c_edge]
                        fig_g0_1_data["switch_truth"] = truth[:, s_edge]
                        for m in methods:
                            fig_g0_1_data[f"clean_{m}"] = data_methods[m]["theta"][:, c_edge]
                            fig_g0_1_data[f"switch_{m}"] = data_methods[m]["theta"][:, s_edge]

    print(f"Audited {rows_audited} baseline rows. All match P2 summary: {all_rows_pass}")
    print(f"Max RMSE difference: {max_rmse_diff:.3e}")
    print(f"Max NIS difference:  {max_nis_diff:.3e}")

    # 4. Compute P bins calibration
    for m in ["KF-CV", "KF-GM"]:
        samples = p_pools[m]
        if not samples:
            continue
        p_arr = np.array([s[0] for s in samples])
        sq_arr = np.array([s[1] for s in samples])

        # 5 quantile bins
        q_edges = np.percentile(p_arr, [0, 20, 40, 60, 80, 100])
        for b_idx in range(5):
            lo, hi = q_edges[b_idx], q_edges[b_idx + 1]
            if b_idx == 4:
                mask = (p_arr >= lo) & (p_arr <= hi)
            else:
                mask = (p_arr >= lo) & (p_arr < hi)
            if not np.any(mask):
                continue
            sub_p = p_arr[mask]
            sub_sq = sq_arr[mask]
            ratio = sub_sq / np.maximum(1e-12, sub_p)
            p_bins_rows.append({
                "method": m,
                "bin_idx": b_idx + 1,
                "p_min": float(lo),
                "p_max": float(hi),
                "n_samples": int(len(sub_p)),
                "mean_P": float(np.mean(sub_p)),
                "mean_sq_err": float(np.mean(sub_sq)),
                "ratio_e2_over_P": float(np.mean(sub_sq) / np.mean(sub_p)),
                "ratio_median": float(np.median(ratio)),
                "ratio_p90": float(np.percentile(ratio, 90)),
                "ratio_p95": float(np.percentile(ratio, 95)),
                "ratio_p99": float(np.percentile(ratio, 99)),
            })

    # 5. Compute descriptive pilot clean rate terciles
    pilot_terciles = [np.nan, np.nan]
    if len(pilot_g_clean_rates) > 0:
        pilot_terciles = [
            float(np.percentile(pilot_g_clean_rates, 33.3333)),
            float(np.percentile(pilot_g_clean_rates, 66.6667)),
        ]
    print(f"G0 Pilot Clean Rate Terciles (descriptive only): q1={pilot_terciles[0]:.6f}, q2={pilot_terciles[1]:.6f}")

    # 6. Save audit.json and BASELINE_AUDIT.md
    audit_status = "PASS" if (all_rows_pass and input_sha_ok and replay_sha_ok and rows_audited == 270) else "FAILURE"
    audit_dict = {
        "status": audit_status,
        "rows_audited": rows_audited,
        "expected_rows": 270,
        "max_rmse_difference": max_rmse_diff,
        "max_nis_difference": max_nis_diff,
        "threshold": 1e-10,
        "input_sha_verified": input_sha_ok,
        "replay_sha_verified": replay_sha_ok,
        "g0_pilot_clean_rate_terciles": {
            "note": "Derived purely from old P2 5 pilot seeds for descriptive exploration. Strictly prohibited for development q calibration.",
            "q1_33pct": pilot_terciles[0],
            "q2_67pct": pilot_terciles[1],
            "n_clean_edges": len(pilot_g_clean_rates),
        },
    }
    with open(new_g0_dir / "audit.json", "w", encoding="utf-8") as f:
        json.dump(audit_dict, f, indent=2)

    baseline_md_content = f"""# P2.1 G0 基线审计与几何分层报告

**审计状态**：`{audit_status}`  
**核对时间**：2026-11-02  
**上游结果目录**：`{old_p2_dir}`  
**审计输出目录**：`{new_g0_dir}`  

---

## 一、P2 原始结果复算与 SHA 校验

- **复算行数**：`{rows_audited} / 270` 行（覆盖全部 15 cases × 2 durations × 3 methods × 3 windows）。
- **最大 RMSE 差异**：`{max_rmse_diff:.3e} rad`（严格 `< 1e-10`）。
- **最大 NIS 差异**：`{max_nis_diff:.3e}`（严格 `< 1e-10`）。
- **输入 MAT 文件 SHA 校验**：`{'全部通过 PASS' if input_sha_ok else '存在异常 FAIL'}`。
- **回放 NPZ 文件 SHA 校验**：`{'全部通过 PASS' if replay_sha_ok else '存在异常 FAIL'}`。

## 二、G0 探索性描述分位点（仅作先导展示，严禁用于开发集调参）

- **G-clean 有效样本边数**：`{len(pilot_g_clean_rates)}`
- **连续角速度变化率三分位**：
  - `q1 (33.3%)`：`{pilot_terciles[0]:.6f} rad/模型时间`
  - `q2 (66.7%)`：`{pilot_terciles[1]:.6f} rad/模型时间`
- **声明**：本分位阈值仅用于历史 5 个种子描述性审计，后续正式开发集实验必须独立从 12 个新开发种子真值分布中冻结，绝不交叉复用。

## 三、交付产物检查清单

1. `audit.json`：精确机器可读校验结果。
2. `geometry_events.csv.gz`：逐 case/edge/time 几何标签及最近镜像分量压缩明细。
3. `geometry_counts.csv`：各几何分层边数与样本占比统计。
4. `error_by_geometry_seed.csv`：按种子和几何分层方法误差（RMSE/MAE/P90/P95/P99/极端误差）。
5. `nis_tail_by_geometry.csv`：新鲜观测时刻按突跳邻域及洁净状态分层的厚尾 NIS。
6. `coverage_by_age.csv`：遮挡年龄分段经验覆盖率与方差失配比。
7. `p_bins.csv`：P 方差分箱校准分析。
8. `figures/fig_G0_1_trajectories.png`：典型几何清洁与镜像跳变方位轨迹。
9. `figures/fig_G0_2_error_nis_tails.png`：几何分层绝对误差尾部 ECDF 及 NIS 厚尾分布。
"""
    with open(new_g0_dir / "BASELINE_AUDIT.md", "w", encoding="utf-8") as f:
        f.write(baseline_md_content)

    # 7. Write DataFrames to CSV
    pd.DataFrame(counts_rows).to_csv(new_g0_dir / "geometry_counts.csv", index=False)
    pd.DataFrame(error_seed_rows).to_csv(new_g0_dir / "error_by_geometry_seed.csv", index=False)
    pd.DataFrame(nis_tail_rows).to_csv(new_g0_dir / "nis_tail_by_geometry.csv", index=False)
    pd.DataFrame(coverage_age_rows).to_csv(new_g0_dir / "coverage_by_age.csv", index=False)
    pd.DataFrame(p_bins_rows).to_csv(new_g0_dir / "p_bins.csv", index=False)

    # 8. Render figures
    # Figure G0-1: Trajectories
    if "clean_edge" in fig_g0_1_data:
        t_arr = fig_g0_1_data["time_steps"]
        D_val = fig_g0_1_data["D"]
        fig, axes = plt.subplots(2, 1, figsize=(10, 7), sharex=True)

        # Clean edge
        ax0 = axes[0]
        ax0.plot(t_arr, fig_g0_1_data["clean_truth"], "k-", lw=2, label="Truth")
        ax0.plot(t_arr, fig_g0_1_data["clean_Hold"], "b--", lw=1.5, label="Hold")
        ax0.plot(t_arr, fig_g0_1_data["clean_KF-CV"], "r-", lw=1.5, label="KF-CV")
        ax0.plot(t_arr, fig_g0_1_data["clean_KF-GM"], "g-.", lw=1.5, label="KF-GM")
        ax0.axvspan(0, D_val * dt, color="gray", alpha=0.2, label=f"Occlusion (D={D_val})")
        ax0.set_ylabel("Bearing (rad)")
        ax0.set_title(
            f"G0-1 Clean Edge (Edge {fig_g0_1_data['clean_edge']}, Case {fig_g0_1_data['case_id']})\n"
            f"[G0 Pilot 5 Reused Seeds | Posterior Geometry Labels | No Causal Claim]"
        )
        ax0.grid(True, ls=":", alpha=0.6)
        ax0.legend(loc="upper right", fontsize=8)

        # Switch edge
        ax1 = axes[1]
        ax1.plot(t_arr, fig_g0_1_data["switch_truth"], "k-", lw=2, label="Truth")
        ax1.plot(t_arr, fig_g0_1_data["switch_Hold"], "b--", lw=1.5, label="Hold")
        ax1.plot(t_arr, fig_g0_1_data["switch_KF-CV"], "r-", lw=1.5, label="KF-CV")
        ax1.plot(t_arr, fig_g0_1_data["switch_KF-GM"], "g-.", lw=1.5, label="KF-GM")
        ax1.axvspan(0, D_val * dt, color="gray", alpha=0.2, label=f"Occlusion (D={D_val})")
        ax1.set_xlabel("Time relative to occlusion start (model units)")
        ax1.set_ylabel("Bearing (rad)")
        ax1.set_title(
            f"G0-1 Switch Edge (Edge {fig_g0_1_data['switch_edge']}, Case {fig_g0_1_data['case_id']})\n"
            f"[G0 Pilot 5 Reused Seeds | Minimum-Image Branch Jump | No Causal Claim]"
        )
        ax1.grid(True, ls=":", alpha=0.6)
        ax1.legend(loc="upper right", fontsize=8)

        plt.tight_layout()
        fig.savefig(figures_dir / "fig_G0_1_trajectories.png", dpi=300)
        plt.close(fig)

    # Figure G0-2: Error Tails and NIS Tails
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # Error tail ECDF
    ax0 = axes[0]
    colors = {
        geometry.LABEL_CLEAN: "green",
        geometry.LABEL_SWITCH: "red",
        geometry.LABEL_CLOSE: "purple",
        geometry.LABEL_AMBIGUOUS: "gray",
    }
    for lbl in geometry.WINDOW_LABELS:
        vals = np.array(fig_g0_2_error_samples[lbl])
        if len(vals) > 0:
            s_vals = np.sort(vals)
            cdf = np.arange(1, len(s_vals) + 1) / float(len(s_vals))
            ax0.plot(s_vals, cdf, color=colors.get(lbl, "black"), lw=1.5, label=f"{lbl} (n={len(vals)})")
    ax0.set_xlabel("Absolute Occlusion Error |e| (rad)")
    ax0.set_ylabel("Empirical CDF")
    ax0.set_title(
        "G0-2 Error ECDF by Geometry Class (KF-CV)\n"
        "[G0 Pilot 5 Reused Seeds | Posterior Labels | No Causal Claim]"
    )
    ax0.grid(True, ls=":", alpha=0.6)
    ax0.legend(loc="lower right", fontsize=8)

    # NIS tail distribution
    ax1 = axes[1]
    nis_colors = {
        geometry.FRESH_CLEAN: "green",
        geometry.FRESH_SWITCH_EVENT: "red",
        geometry.FRESH_SWITCH_NEAR2: "orange",
        geometry.FRESH_CLOSE: "purple",
        geometry.FRESH_AMBIGUOUS: "gray",
    }
    for fl_name in [
        geometry.FRESH_CLEAN,
        geometry.FRESH_SWITCH_EVENT,
        geometry.FRESH_SWITCH_NEAR2,
        geometry.FRESH_CLOSE,
    ]:
        vals = np.array(fig_g0_2_nis_samples[fl_name])
        if len(vals) > 0:
            s_vals = np.sort(vals)
            cdf = np.arange(1, len(s_vals) + 1) / float(len(s_vals))
            ax1.plot(
                s_vals,
                1.0 - cdf,
                color=nis_colors.get(fl_name, "black"),
                lw=1.5,
                label=f"{fl_name} (n={len(vals)})",
            )
    ax1.set_xscale("log")
    ax1.set_yscale("log")
    ax1.set_xlabel("Normalized Innovation Squared (NIS)")
    ax1.set_ylabel("Tail Probability P(NIS > x)")
    ax1.set_title(
        "G0-2 Fresh NIS Tail by Geometry State (KF-CV)\n"
        "[G0 Pilot 5 Reused Seeds | Posterior Labels | No Causal Claim]"
    )
    ax1.grid(True, ls=":", alpha=0.6)
    ax1.legend(loc="upper right", fontsize=8)

    plt.tight_layout()
    fig.savefig(figures_dir / "fig_G0_2_error_nis_tails.png", dpi=300)
    plt.close(fig)

    print("=== G0 Audit Completed Successfully ===")


if __name__ == "__main__":
    if len(sys.argv) < 3:
        sys.exit("Usage: python audit_p2_geometry.py <OLD_P2_RESULT_DIR> <NEW_G0_DIR>")
    run_g0_audit(sys.argv[1], sys.argv[2])

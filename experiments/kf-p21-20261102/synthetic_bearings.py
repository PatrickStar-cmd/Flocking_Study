"""experiments/kf-p21-20261102/synthetic_bearings.py

Synthetic bearing benchmarks for controlled motion mechanism verification (G1).
Generates 8 canonical analytical configurations (A1, B1-B3, C1-C2, D1, E1) across
20 fixed random streams and 2 occlusion durations (D=20, 120), totaling 320 conditions.

Evaluates Hold, KF-CV (q=1e-4), and KF-GM (q=1e-5, tau=10) without future leakage.
Verifies minimum-image crossing geometry (E1) and near-range geometry (D1).
Standalone module with importable generator functions and CLI runner.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import sys
from typing import Any, Dict, List, Optional, Tuple, Union

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Import geometry
sys.path.insert(0, str(Path(__file__).resolve().parent))
import geometry

# Fixed experiment constants per PROTOCOL and RUN_RULES
SEEDS = [2026113000 + i for i in range(1, 21)]  # 2026113001 .. 2026113020
DURATIONS = [20, 120]
CONFIGS = ["A1", "B1", "B2", "B3", "C1", "C2", "D1", "E1"]
T_STEPS = 720
DT = 0.3
T0_1BASED = 201
PRE_STEPS = 200
R_VAR = 1e-4
THETA0 = 0.7
OMEGA0 = 0.01
ACCEL_VAL = 0.0005


def compute_sha256(filepath: Union[str, Path]) -> str:
    """Compute SHA256 hex digest of a file in streaming chunks."""
    h = hashlib.sha256()
    with open(filepath, "rb") as f:
        while chunk := f.read(65536):
            h.update(chunk)
    return h.hexdigest()


def wrap_to_pi(angle: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
    """Wrap angle differences to [-pi, pi)."""
    return np.arctan2(np.sin(angle), np.cos(angle))


def mod_2pi(angle: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
    """Wrap angle to [0, 2*pi)."""
    return np.mod(angle, 2.0 * np.pi)


def load_p2_replay_module() -> Any:
    """Import stock P2 replay module without modifying preexisting files."""
    p2_dir = Path(__file__).resolve().parent.parent / "kf-p2-20261102"
    replay_path = p2_dir / "replay.py"
    if not replay_path.exists():
        raise FileNotFoundError(f"P2 replay.py not found at {replay_path}")

    spec = importlib.util.spec_from_file_location("p2_replay_stock", replay_path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Could not load spec for {replay_path}")
    p2_mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(p2_mod)
    return p2_mod


def generate_synthetic_trajectory(
    config_id: str,
    seed: int,
    D: int,
    *,
    dt: float = DT,
    T: int = T_STEPS,
    t0_1based: int = T0_1BASED,
    R: float = R_VAR,
) -> Dict[str, np.ndarray]:
    """Generate ground truth, observations, and geometric flags for a synthetic condition.

    Parameters
    ----------
    config_id : str
        One of 'A1', 'B1', 'B2', 'B3', 'C1', 'C2', 'D1', 'E1'.
    seed : int
        RNG seed for observation noise.
    D : int
        Occlusion duration in steps (20 or 120).
    dt : float
        Sampling step in model time (default 0.3).
    T : int
        Total steps (default 720).
    t0_1based : int
        1-based step where occlusion begins (default 201).
    R : float
        Observation noise variance (default 1e-4).

    Returns
    -------
    dict with keys:
      'truth': [T] float64 in [0, 2*pi)
      'observations': [T] float64 in [0, 2*pi)
      'missing': [T] bool
      'noise': [T] float64
      'distance': [T] float64
      'image_switch': [T] bool
      'near_range': [T] bool
      'near_half_box': [T] bool
    """
    if config_id not in CONFIGS:
        raise ValueError(f"Unknown config_id: {config_id}. Expected one of {CONFIGS}")

    t0_0based = t0_1based - 1
    t_steps = np.arange(T, dtype=np.float64)
    t_model = (t_steps - t0_0based) * dt
    t_cross = (D * dt) / 2.0

    truth = np.zeros(T, dtype=np.float64)
    distance = np.full(T, 50.0, dtype=np.float64)  # Default distant
    image_switch = np.zeros(T, dtype=bool)
    near_range = np.zeros(T, dtype=bool)
    near_half_box = np.zeros(T, dtype=bool)

    if config_id == "A1":
        # Constant bearing
        truth = np.full(T, THETA0, dtype=np.float64)

    elif config_id == "B1":
        # Constant angular rate omega = 0.002
        truth = mod_2pi(THETA0 + 0.002 * t_model)

    elif config_id == "B2":
        # Constant angular rate omega = 0.01
        truth = mod_2pi(THETA0 + 0.01 * t_model)

    elif config_id == "B3":
        # Constant angular rate omega = 0.03
        truth = mod_2pi(THETA0 + 0.03 * t_model)

    elif config_id == "C1":
        # Angular acceleration: omega0=0.01, a=+0.0005
        truth = mod_2pi(THETA0 + OMEGA0 * t_model + 0.5 * ACCEL_VAL * (t_model**2))

    elif config_id == "C2":
        # Angular acceleration: omega0=0.01, a=-0.0005
        truth = mod_2pi(THETA0 + OMEGA0 * t_model + 0.5 * (-ACCEL_VAL) * (t_model**2))

    elif config_id == "D1":
        # Close passage without boundary jump
        # Target: x = 0.5 * (t_model - t_cross), y = 2.0. Observer: (0, 0).
        dx = 0.5 * (t_model - t_cross)
        dy = np.full(T, 2.0, dtype=np.float64)
        L = 1000.0

        idx_x = geometry.round_matlab(dx / L)
        idx_y = geometry.round_matlab(dy / L)
        dx_min = dx - L * idx_x
        dy_min = dy - L * idx_y

        truth = mod_2pi(np.arctan2(dy_min, dx_min))
        distance = np.sqrt(dx_min**2 + dy_min**2)
        near_range = distance < 5.0
        min_margin = np.minimum(
            np.abs(np.abs(dx_min) - (L / 2.0)),
            np.abs(np.abs(dy_min) - (L / 2.0)),
        )
        near_half_box = min_margin < 10.0
        diff_dmin = np.abs(np.diff(np.stack([dx_min, dy_min], axis=1), axis=0))
        image_switch[1:] = np.any(diff_dmin > (L / 2.0), axis=1)

    elif config_id == "E1":
        # Periodic minimum-image boundary switch
        # Target: x = 500.0 + 0.5 * (t_model - t_cross), y = 100.0. Observer: (0, 0).
        dx = 500.0 + 0.5 * (t_model - t_cross)
        dy = np.full(T, 100.0, dtype=np.float64)
        L = 1000.0

        idx_x = geometry.round_matlab(dx / L)
        idx_y = geometry.round_matlab(dy / L)
        dx_min = dx - L * idx_x
        dy_min = dy - L * idx_y

        truth = mod_2pi(np.arctan2(dy_min, dx_min))
        distance = np.sqrt(dx_min**2 + dy_min**2)
        near_range = distance < 5.0
        min_margin = np.minimum(
            np.abs(np.abs(dx_min) - (L / 2.0)),
            np.abs(np.abs(dy_min) - (L / 2.0)),
        )
        near_half_box = min_margin < 10.0
        diff_dmin = np.abs(np.diff(np.stack([dx_min, dy_min], axis=1), axis=0))
        image_switch[1:] = np.any(diff_dmin > (L / 2.0), axis=1)

    # Underlying noise: same 720 standard normal stream per seed and config
    rng = np.random.default_rng(seed)
    std_noise = rng.standard_normal(T)
    noise = std_noise * np.sqrt(R)
    observations = mod_2pi(truth + noise)

    missing = np.zeros(T, dtype=bool)
    missing[t0_0based : t0_0based + D] = True

    return {
        "truth": truth,
        "observations": observations,
        "missing": missing,
        "noise": noise,
        "distance": distance,
        "image_switch": image_switch,
        "near_range": near_range,
        "near_half_box": near_half_box,
    }


def run_synthetic_condition(
    config_id: str,
    seed: int,
    D: int,
    p2_mod: Any,
) -> Dict[str, Any]:
    """Run Hold, KF-CV, and KF-GM on a synthetic benchmark condition."""
    traj = generate_synthetic_trajectory(config_id, seed, D)

    obs_2d = traj["observations"][:, None]  # [T, 1]
    miss_2d = traj["missing"][:, None]  # [T, 1]

    # Run Hold
    res_hold = p2_mod.hold(obs_2d, miss_2d, dt=DT)

    # Run KF-CV default (q=1e-4)
    res_cv = p2_mod.replay(obs_2d, miss_2d, dt=DT, q=1e-4, R=R_VAR, tau=None)

    # Run KF-GM default (q=1e-5, tau=10)
    res_gm = p2_mod.replay(obs_2d, miss_2d, dt=DT, q=1e-5, R=R_VAR, tau=10.0)

    # Compute offline errors
    truth_1d = traj["truth"]
    err_hold = wrap_to_pi(res_hold["theta"][:, 0] - truth_1d)
    err_cv = wrap_to_pi(res_cv["theta"][:, 0] - truth_1d)
    err_gm = wrap_to_pi(res_gm["theta"][:, 0] - truth_1d)

    # Occlusion slice
    t0_0based = T0_1BASED - 1
    occ_sl = slice(t0_0based, t0_0based + D)

    # Compute occlusion summary metrics
    def calc_metrics(e: np.ndarray, p: np.ndarray) -> Dict[str, float]:
        sq = e**2
        rmse = float(np.sqrt(np.mean(sq)))
        mae = float(np.mean(np.abs(e)))
        p95 = float(np.percentile(np.abs(e), 95))
        max_e = float(np.max(np.abs(e)))
        if np.all(np.isnan(p)):
            cov = np.nan
            mean_p = np.nan
        else:
            cov = float(np.mean(np.abs(e) <= 1.96 * np.sqrt(np.maximum(0.0, p))))
            mean_p = float(np.mean(p))
        return {
            "rmse": rmse,
            "mae": mae,
            "p95": p95,
            "max_err": max_e,
            "mean_p": mean_p,
            "coverage95": cov,
        }

    met_hold = calc_metrics(err_hold[occ_sl], res_hold["p11"][occ_sl, 0])
    met_cv = calc_metrics(err_cv[occ_sl], res_cv["p11"][occ_sl, 0])
    met_gm = calc_metrics(err_gm[occ_sl], res_gm["p11"][occ_sl, 0])

    return {
        "traj": traj,
        "res_hold": res_hold,
        "res_cv": res_cv,
        "res_gm": res_gm,
        "error_hold": err_hold,
        "error_cv": err_cv,
        "error_gm": err_gm,
        "metrics_hold": met_hold,
        "metrics_cv": met_cv,
        "metrics_gm": met_gm,
    }


def run_all_synthetic(output_dir_arg: str) -> None:
    """Execute all 320 synthetic conditions and save NPZs, CSVs, and figures."""
    out_dir = Path(output_dir_arg).resolve()

    # Guard: strictly refuse overwrite
    if out_dir.exists():
        sys.exit(f"Error: Refusing to overwrite existing output directory: {out_dir}")

    out_dir.mkdir(parents=True, exist_ok=False)
    raw_dir = out_dir / "raw"
    raw_dir.mkdir(parents=True, exist_ok=False)
    fig_dir = out_dir / "figures"
    fig_dir.mkdir(parents=True, exist_ok=False)

    p2_mod = load_p2_replay_module()

    print(f"=== Synthetic Bearings Verification (G1) ===")
    print(f"Output directory: {out_dir}")
    print(f"Seeds: {len(SEEDS)} streams ({SEEDS[0]}..{SEEDS[-1]})")
    print(f"Configs: {len(CONFIGS)} ({', '.join(CONFIGS)})")
    print(f"Durations: {DURATIONS}")
    print(f"Total conditions: {len(SEEDS) * len(CONFIGS) * len(DURATIONS)}")

    summary_rows: List[Dict[str, Any]] = []
    age_metric_rows: List[Dict[str, Any]] = []
    geom_verify_rows: List[Dict[str, Any]] = []
    manifest_runs: List[Dict[str, Any]] = []

    # Store first stream D120 predictions for figure
    first_stream_data: Dict[str, Any] = {}

    for s_idx, seed in enumerate(SEEDS):
        for cfg in CONFIGS:
            for D in DURATIONS:
                res = run_synthetic_condition(cfg, seed, D, p2_mod)
                traj = res["traj"]

                # Save compressed NPZ per condition
                npz_name = f"case_{cfg}_seed{seed}_D{D}.npz"
                npz_path = raw_dir / npz_name

                np.savez_compressed(
                    npz_path,
                    truth=traj["truth"][:, None],
                    observations=traj["observations"][:, None],
                    missing=traj["missing"][:, None],
                    distance=traj["distance"][:, None],
                    image_switch=traj["image_switch"][:, None],
                    near_range=traj["near_range"][:, None],
                    near_half_box=traj["near_half_box"][:, None],
                    theta_Hold=res["res_hold"]["theta"],
                    p11_Hold=res["res_hold"]["p11"],
                    error_Hold=res["error_hold"][:, None],
                    theta_KF_CV=res["res_cv"]["theta"],
                    p11_KF_CV=res["res_cv"]["p11"],
                    nis_KF_CV=res["res_cv"]["nis"],
                    error_KF_CV=res["error_cv"][:, None],
                    theta_KF_GM=res["res_gm"]["theta"],
                    p11_KF_GM=res["res_gm"]["p11"],
                    nis_KF_GM=res["res_gm"]["nis"],
                    error_KF_GM=res["error_gm"][:, None],
                    config_id=cfg,
                    seed=seed,
                    D=D,
                    dt=DT,
                    t0=T0_1BASED,
                )

                file_sha = compute_sha256(npz_path)
                manifest_runs.append({
                    "config_id": cfg,
                    "seed": seed,
                    "D": D,
                    "npz_file": npz_name,
                    "npz_sha256": file_sha,
                })

                # Record summary metrics
                for m, met in [
                    ("Hold", res["metrics_hold"]),
                    ("KF-CV", res["metrics_cv"]),
                    ("KF-GM", res["metrics_gm"]),
                ]:
                    summary_rows.append({
                        "config_id": cfg,
                        "seed": seed,
                        "D": D,
                        "method": m,
                        "circular_rmse": met["rmse"],
                        "circular_mae": met["mae"],
                        "p95": met["p95"],
                        "max_err": met["max_err"],
                        "mean_P": met["mean_p"],
                        "approx95coverage": met["coverage95"],
                    })

                # Per-age metrics for error-vs-age plot
                t0_0 = T0_1BASED - 1
                for step_idx in range(1, D + 1):
                    t_cur = t0_0 + step_idx - 1
                    age_metric_rows.append({
                        "config_id": cfg,
                        "seed": seed,
                        "D": D,
                        "age_step": step_idx,
                        "age_time": round(step_idx * DT, 4),
                        "err_Hold": abs(float(res["error_hold"][t_cur])),
                        "err_KF_CV": abs(float(res["error_cv"][t_cur])),
                        "err_KF_GM": abs(float(res["error_gm"][t_cur])),
                    })

                # Verify geometry mechanism
                # D1 must verify near_range in occlusion; E1 must verify image_switch
                sw_count = int(np.sum(traj["image_switch"][t0_0 : t0_0 + D]))
                nr_count = int(np.sum(traj["near_range"][t0_0 : t0_0 + D]))
                geom_verify_rows.append({
                    "config_id": cfg,
                    "seed": seed,
                    "D": D,
                    "image_switch_count_occ": sw_count,
                    "near_range_count_occ": nr_count,
                    "min_distance_occ": float(np.min(traj["distance"][t0_0 : t0_0 + D])),
                })

                # Store first stream D120
                if seed == SEEDS[0] and D == 120:
                    first_stream_data[cfg] = {
                        "truth": traj["truth"],
                        "theta_Hold": res["res_hold"]["theta"][:, 0],
                        "theta_KF_CV": res["res_cv"]["theta"][:, 0],
                        "theta_KF_GM": res["res_gm"]["theta"][:, 0],
                        "missing": traj["missing"],
                    }

    # Save manifest
    manifest_data = {
        "status": "COMPLETED",
        "total_conditions": len(manifest_runs),
        "dt": DT,
        "T": T_STEPS,
        "t0": T0_1BASED,
        "R": R_VAR,
        "seeds": SEEDS,
        "configs": CONFIGS,
        "durations": DURATIONS,
        "runs": manifest_runs,
    }
    with open(out_dir / "manifest.json", "w", encoding="utf-8") as f:
        json.dump(manifest_data, f, indent=2)

    # Save CSV summaries
    df_summary = pd.DataFrame(summary_rows)
    df_summary.to_csv(out_dir / "synthetic_summary.csv", index=False)

    df_age = pd.DataFrame(age_metric_rows)
    df_age.to_csv(out_dir / "metrics_by_age.csv", index=False)

    df_geom = pd.DataFrame(geom_verify_rows)
    df_geom.to_csv(out_dir / "geometry_events_summary.csv", index=False)

    print("=== Geometry Verification Summary ===")
    d1_check = df_geom[df_geom["config_id"] == "D1"]
    e1_check = df_geom[df_geom["config_id"] == "E1"]
    print(f"D1 Near-range verified: all min_dist < 5: {np.all(d1_check['min_distance_occ'] < 5.0)}")
    print(f"D1 Image switches: {int(d1_check['image_switch_count_occ'].sum())} (expected 0)")
    print(f"E1 Image switches: all >= 1 switch in occ: {np.all(e1_check['image_switch_count_occ'] >= 1)}")
    print(f"E1 Near-range counts: {int(e1_check['near_range_count_occ'].sum())} (expected 0)")

    # Render Figures
    # 1. Error vs Age across Configs A-E for D=120
    df_age_120 = df_age[df_age["D"] == 120]
    grouped_age = (
        df_age_120.groupby(["config_id", "age_time"])[
            ["err_Hold", "err_KF_CV", "err_KF_GM"]
        ]
        .mean()
        .reset_index()
    )

    configs_to_plot = ["A1", "B2", "C1", "D1", "E1"]
    fig, axes = plt.subplots(1, 5, figsize=(18, 4), sharey=False)
    for idx, cfg in enumerate(configs_to_plot):
        ax = axes[idx]
        sub = grouped_age[grouped_age["config_id"] == cfg]
        ax.plot(sub["age_time"], sub["err_Hold"], "b--", lw=1.5, label="Hold")
        ax.plot(sub["age_time"], sub["err_KF_CV"], "r-", lw=1.5, label="KF-CV")
        ax.plot(sub["age_time"], sub["err_KF_GM"], "g-.", lw=1.5, label="KF-GM")
        ax.set_title(f"Config {cfg}", fontsize=11)
        ax.set_xlabel("Occlusion age (model time)")
        if idx == 0:
            ax.set_ylabel("Mean Absolute Error (rad)")
        ax.grid(True, ls=":", alpha=0.6)
        if idx == 0:
            ax.legend(loc="upper left", fontsize=8)

    plt.suptitle("Error vs Occlusion Age across Synthetic Configs (D=120, 20 Streams Mean)", y=1.02)
    plt.tight_layout()
    fig.savefig(fig_dir / "error_vs_age_A_E.png", dpi=300)
    plt.close(fig)

    # 2. Signed Prediction for First Stream D=120
    t_model_arr = (np.arange(T_STEPS) - (T0_1BASED - 1)) * DT
    fig, axes = plt.subplots(len(configs_to_plot), 1, figsize=(10, 10), sharex=True)
    for idx, cfg in enumerate(configs_to_plot):
        ax = axes[idx]
        data = first_stream_data[cfg]
        ax.plot(t_model_arr, data["truth"], "k-", lw=2, label="Truth")
        ax.plot(t_model_arr, data["theta_Hold"], "b--", lw=1.2, label="Hold")
        ax.plot(t_model_arr, data["theta_KF_CV"], "r-", lw=1.2, label="KF-CV")
        ax.plot(t_model_arr, data["theta_KF_GM"], "g-.", lw=1.2, label="KF-GM")
        ax.axvspan(0, 120 * DT, color="gray", alpha=0.2, label="Occlusion")
        ax.set_ylabel("Angle (rad)", fontsize=9)
        ax.set_title(f"Config {cfg} (Seed {SEEDS[0]}, D=120)", fontsize=10)
        ax.grid(True, ls=":", alpha=0.5)
        if idx == 0:
            ax.legend(loc="upper right", fontsize=7)

    axes[-1].set_xlabel("Time relative to occlusion (model units)")
    plt.tight_layout()
    fig.savefig(fig_dir / "signed_prediction_firststream_D120.png", dpi=300)
    plt.close(fig)

    print("=== Synthetic Bearings Completed Successfully ===")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit("Usage: python synthetic_bearings.py <NEW_OUTPUT_DIR>")
    run_all_synthetic(sys.argv[1])

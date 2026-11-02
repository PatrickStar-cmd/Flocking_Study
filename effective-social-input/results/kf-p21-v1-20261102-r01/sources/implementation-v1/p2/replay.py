"""experiments/kf-p2-20261102/replay.py

Open-loop estimator replay for Flocking Study P2 diagnostic experiment.
Implements exact Python ports of stock MATLAB bearingTrackStep (CV) and
bearingTrackStepGaussMarkov (GM), alongside the Hold-last-observation estimator.

Truth is strictly isolated from filter execution; angle errors are computed offline.
Outputs are saved to <RESULT_DIR>/raw/replay/*.npz with SHA256 manifest auditing.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import sys
from typing import Any, Dict, Optional, Tuple, Union

import numpy as np
import scipy.io


def wrap_to_pi(angle: np.ndarray) -> np.ndarray:
    """Wrap angle to [-pi, pi) using atan2(sin, cos)."""
    return np.arctan2(np.sin(angle), np.cos(angle))


def mod_2pi(angle: np.ndarray) -> np.ndarray:
    """Wrap angle to [0, 2*pi) matching MATLAB mod(angle, 2*pi)."""
    return np.mod(angle, 2.0 * np.pi)


def compute_sha256(filepath: Union[str, Path]) -> str:
    """Compute SHA256 hex digest of a file in streaming chunks."""
    h = hashlib.sha256()
    with open(filepath, "rb") as f:
        while chunk := f.read(65536):
            h.update(chunk)
    return h.hexdigest()


def replay(
    observations: np.ndarray,
    missing: np.ndarray,
    *,
    dt: float,
    q: float,
    R: float,
    tau: Optional[float] = None,
) -> Dict[str, np.ndarray]:
    """Open-loop Kalman Filter replay across time and edges.

    Exact port of stock MATLAB functions:
      - bearingTrackStep (tau is None or np.isinf(tau))
      - bearingTrackStepGaussMarkov (tau is finite positive float)

    Parameters
    ----------
    observations : np.ndarray, shape [T, E]
        Noisy bearing observations in radians. All elements must be finite.
    missing : np.ndarray, shape [T] or [T, E]
        Boolean array where True indicates hidden / occluded observation.
    dt : float
        Sampling time step, must be positive.
    q : float
        Process noise spectral density, must be non-negative.
    R : float
        Measurement noise variance, must be positive.
    tau : float, optional
        Velocity decorrelation time constant for Gauss-Markov model.
        If None or infinite, exact Constant Velocity (CV) model is used.

    Returns
    -------
    dict with keys:
      'theta': [T, E] float64 posterior/predicted bearing in [0, 2*pi)
      'p11':   [T, E] float64 bearing estimation variance
      'p12':   [T, E] float64 bearing-velocity covariance
      'p22':   [T, E] float64 angular velocity variance
      'nis':   [T, E] float64 normalized innovation squared (NaN on missing/init)
      'age':   [T, E] float64 cumulative age since last fresh observation
    """
    obs = np.asarray(observations, dtype=np.float64)
    if obs.ndim != 2:
        raise ValueError(f"observations must be 2D [T, E], got shape {obs.shape}")
    T, E = obs.shape
    if T == 0 or E == 0:
        raise ValueError("observations must be nonempty")
    if not np.all(np.isfinite(obs)):
        raise ValueError("All input observations must be finite numbers")

    miss = np.asarray(missing, dtype=bool)
    if miss.ndim == 1:
        if miss.shape[0] != T:
            raise ValueError(f"missing length {miss.shape[0]} does not match T={T}")
        miss = np.broadcast_to(miss[:, None], (T, E))
    elif miss.ndim == 2:
        if miss.shape != (T, E):
            raise ValueError(f"missing shape {miss.shape} does not match (T, E)={(T, E)}")
    else:
        raise ValueError(f"missing must be 1D [T] or 2D [T, E], got shape {miss.shape}")

    if dt <= 0.0 or not np.isfinite(dt):
        raise ValueError(f"dt must be finite positive float, got {dt}")
    if q < 0.0 or not np.isfinite(q):
        raise ValueError(f"q must be finite non-negative float, got {q}")
    if R <= 0.0 or not np.isfinite(R):
        raise ValueError(f"R must be finite positive float, got {R}")

    # Build transition matrix F and process noise covariance Q
    if tau is None or np.isinf(tau):
        # Constant Velocity (CV) model per stock bearingTrackStep.m
        # F = [1 dt; 0 1]
        # Q = q * [dt^3/3 dt^2/2; dt^2/2 dt]
        F = np.array([[1.0, dt], [0.0, 1.0]], dtype=np.float64)
        Q = q * np.array(
            [[dt**3 / 3.0, dt**2 / 2.0], [dt**2 / 2.0, dt]],
            dtype=np.float64,
        )
    else:
        if tau <= 0.0 or not np.isfinite(tau):
            raise ValueError(f"tau must be finite positive float, got {tau}")
        # Gauss-Markov (GM) model per stock bearingTrackStepGaussMarkov.m
        rho = np.exp(-dt / tau)
        F = np.array([[1.0, tau * (1.0 - rho)], [0.0, rho]], dtype=np.float64)
        q11 = tau**2 * dt - tau**3 * (2.0 * (1.0 - rho) - (1.0 - rho**2) / 2.0)
        q12 = tau**2 * (1.0 - rho) ** 2 / 2.0
        q22 = tau * (1.0 - rho**2) / 2.0
        Q = q * np.array([[q11, q12], [q12, q22]], dtype=np.float64)
        Q = 0.5 * (Q + Q.T)

    theta = np.zeros((T, E), dtype=np.float64)
    p11 = np.zeros((T, E), dtype=np.float64)
    p12 = np.zeros((T, E), dtype=np.float64)
    p22 = np.zeros((T, E), dtype=np.float64)
    nis = np.full((T, E), np.nan, dtype=np.float64)
    age = np.zeros((T, E), dtype=np.float64)

    # State vectors x: [E, 2], covariances P: [E, 2, 2], age_curr: [E]
    x = np.zeros((E, 2), dtype=np.float64)
    P = np.zeros((E, 2, 2), dtype=np.float64)
    age_curr = np.zeros(E, dtype=np.float64)
    initialized = np.zeros(E, dtype=bool)

    # Pre-allocate identity matrix for Joseph covariance update
    I2 = np.eye(2, dtype=np.float64)

    for t in range(T):
        y_t = obs[t]  # [E]
        m_t = miss[t]  # [E] bool

        # 1. Uninitialized edges receiving first finite observation
        first_init = (~initialized) & (~m_t)
        if np.any(first_init):
            # Stock init: x = [mod(observation, 2*pi); 0]; P = diag([R, 1])
            x[first_init, 0] = mod_2pi(y_t[first_init])
            x[first_init, 1] = 0.0
            P[first_init, 0, 0] = R
            P[first_init, 0, 1] = 0.0
            P[first_init, 1, 0] = 0.0
            P[first_init, 1, 1] = 1.0
            age_curr[first_init] = 0.0
            initialized[first_init] = True
            # For newly initialized tracks: age=0, nis=NaN, no prediction step

        # 2. Tracks that were already initialized prior to step t
        already_init = initialized & (~first_init)
        if np.any(already_init):
            # Predict step: x = F * x; P = F * P * F' + Q; age = age + dt
            # x[already_init] is [K, 2], F is [2, 2]: x @ F.T is [K, 2]
            x[already_init] = x[already_init] @ F.T
            # P[already_init] is [K, 2, 2]: F @ P @ F.T + Q
            P[already_init] = (
                np.matmul(F, np.matmul(P[already_init], F.T)) + Q[None, :, :]
            )
            age_curr[already_init] += dt

            # Update step for fresh observations
            fresh = already_init & (~m_t)
            if np.any(fresh):
                # Wrapped innovation: atan2(sin(y - x1), cos(y - x1))
                delta = y_t[fresh] - x[fresh, 0]
                innov = np.arctan2(np.sin(delta), np.cos(delta))

                # Kalman gain: K = P(:, 1) / (P(1, 1) + R)
                denom = P[fresh, 0, 0] + R
                K = P[fresh, :, 0] / denom[:, None]  # [K, 2]

                # NIS: innov^2 / (P(1, 1) + R)
                nis[t, fresh] = (innov**2) / denom
                age_curr[fresh] = 0.0

                # State update: x = x + K * innov
                x[fresh, 0] += K[:, 0] * innov
                x[fresh, 1] += K[:, 1] * innov

                # Joseph covariance update:
                # A = eye(2) - K * [1, 0]
                # P = A * P * A' + K * R * K'
                A = np.broadcast_to(I2[None, :, :], (np.sum(fresh), 2, 2)).copy()
                A[:, 0, 0] -= K[:, 0]
                A[:, 1, 0] -= K[:, 1]

                APA = np.matmul(A, np.matmul(P[fresh], A.swapaxes(1, 2)))
                KRK = R * np.matmul(K[:, :, None], K[:, None, :])
                P_up = APA + KRK
                P[fresh] = 0.5 * (P_up + P_up.swapaxes(1, 2))

            # Post-step angle wrapping and covariance symmetrizing
            x[already_init, 0] = mod_2pi(x[already_init, 0])
            P[already_init] = 0.5 * (
                P[already_init] + P[already_init].swapaxes(1, 2)
            )

        # Record step outputs
        theta[t] = np.where(initialized, x[:, 0], np.nan)
        p11[t] = np.where(initialized, P[:, 0, 0], np.nan)
        p12[t] = np.where(initialized, P[:, 0, 1], np.nan)
        p22[t] = np.where(initialized, P[:, 1, 1], np.nan)
        age[t] = np.where(initialized, age_curr, np.inf)

    return {
        "theta": theta,
        "p11": p11,
        "p12": p12,
        "p22": p22,
        "nis": nis,
        "age": age,
    }


def hold(
    observations: np.ndarray,
    missing: np.ndarray,
    *,
    dt: float,
) -> Dict[str, np.ndarray]:
    """Hold-last-observation estimator.

    Maintains the last actual observed bearing (in [0, 2*pi)).
    Age increments by dt on missing steps and resets to 0 on fresh steps.
    Covariance and NIS fields are defined as NaN (Hold has no calibrated covariance).

    Parameters
    ----------
    observations : np.ndarray, shape [T, E]
        Bearing observations in radians.
    missing : np.ndarray, shape [T] or [T, E]
        Boolean missing mask.
    dt : float
        Sampling step time in seconds.

    Returns
    -------
    dict with keys: 'theta', 'p11', 'p12', 'p22', 'nis', 'age'
    """
    obs = np.asarray(observations, dtype=np.float64)
    if obs.ndim != 2:
        raise ValueError(f"observations must be 2D [T, E], got shape {obs.shape}")
    T, E = obs.shape
    if T == 0 or E == 0:
        raise ValueError("observations must be nonempty")
    if not np.all(np.isfinite(obs)):
        raise ValueError("All input observations must be finite numbers")

    miss = np.asarray(missing, dtype=bool)
    if miss.ndim == 1:
        if miss.shape[0] != T:
            raise ValueError(f"missing length {miss.shape[0]} does not match T={T}")
        miss = np.broadcast_to(miss[:, None], (T, E))
    elif miss.ndim == 2:
        if miss.shape != (T, E):
            raise ValueError(f"missing shape {miss.shape} does not match (T, E)={(T, E)}")
    else:
        raise ValueError(f"missing must be 1D [T] or 2D [T, E], got shape {miss.shape}")

    if dt <= 0.0 or not np.isfinite(dt):
        raise ValueError(f"dt must be positive float, got {dt}")

    theta = np.zeros((T, E), dtype=np.float64)
    age = np.zeros((T, E), dtype=np.float64)
    p11 = np.full((T, E), np.nan, dtype=np.float64)
    p12 = np.full((T, E), np.nan, dtype=np.float64)
    p22 = np.full((T, E), np.nan, dtype=np.float64)
    nis = np.full((T, E), np.nan, dtype=np.float64)

    curr_theta = np.full(E, np.nan)
    curr_age = np.full(E, np.inf)
    for t in range(T):
        m_t = miss[t]
        curr_age = np.where(m_t, curr_age + dt, 0.0)
        curr_theta = np.where(m_t, curr_theta, mod_2pi(obs[t]))
        theta[t] = curr_theta
        age[t] = curr_age

    return {
        "theta": theta,
        "p11": p11,
        "p12": p12,
        "p22": p22,
        "nis": nis,
        "age": age,
    }


def compute_active_flags(
    age: np.ndarray,
    p11: np.ndarray,
    missing: np.ndarray,
    *,
    max_age: float,
    max_variance: float,
    is_hold: bool = False,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Compute active flags matching stock MATLAB gating criteria.

    Parameters
    ----------
    age : np.ndarray, shape [T, E]
        Cumulative estimator age in model time.
    p11 : np.ndarray, shape [T, E]
        Variance of angle estimate.
    missing : np.ndarray, shape [T, E]
        Missing observation indicator.
    max_age : float
        Maximum allowable age threshold.
    max_variance : float
        Maximum allowable variance threshold.
    is_hold : bool
        If True, variance threshold is not applicable.

    Returns
    -------
    active_by_age : [T, E] bool
    active_by_P : [T, E] bool
    active : [T, E] bool
    """
    # MATLAB uses exact '>' threshold check:
    # age > opt.maxAge -> cutoff; P(1,1) > opt.maxVariance -> cutoff
    active_by_age = age <= max_age

    if is_hold:
        active_by_P = np.full(age.shape, False, dtype=bool)
        # On fresh: always active. On hidden: gated by age only.
        active = np.where(missing, active_by_age, True)
    else:
        active_by_P = p11 <= max_variance
        # On fresh: always active. On hidden: gated by both age and P.
        active = np.where(missing, active_by_age & active_by_P, True)

    return active_by_age, active_by_P, active


def run_replay_cli(result_dir_arg: str) -> None:
    """CLI execution for replay.py."""
    result_dir = Path(result_dir_arg).resolve()
    if not result_dir.exists():
        sys.exit(f"Error: RESULT_DIR does not exist: {result_dir}")

    design_path = result_dir / "design.json"
    if not design_path.exists():
        sys.exit(f"Error: design.json not found in {result_dir}")

    with open(design_path, "r", encoding="utf-8") as f:
        design = json.load(f)

    dt = float(design["dt"])
    R = float(design["R"])
    q_cv = float(design["q_cv"])
    q_gm = float(design["q_gm"])
    tau = float(design["tau"])
    max_age = float(design["max_age"])
    max_variance = float(design["max_variance"])
    durations = [int(d) for d in design["durations"]]
    pre_steps = int(design.get("pre_steps", 200))
    recovery_steps = int(design.get("recovery_steps", 400))
    cases = design["cases"]

    out_dir = result_dir / "raw" / "replay"
    out_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = out_dir / "manifest.json"
    if manifest_path.exists() or any(out_dir.glob('*.npz')):
        raise FileExistsError("Replay outputs already exist; use a new experiment directory.")

    # Load existing manifest if present to preserve failure/skip records
    manifest_records = []
    if manifest_path.exists():
        try:
            with open(manifest_path, "r", encoding="utf-8") as f:
                manifest_records = json.load(f).get("runs", [])
        except Exception as e:
            print(f"Warning: Could not read existing manifest ({e}), starting fresh.")

    existing_manifest_map = {
        (r.get("case_id"), r.get("D")): r for r in manifest_records
    }

    total_tasks = len(cases) * len(durations)
    print(f"=== Replay: {len(cases)} cases x {len(durations)} durations = {total_tasks} runs ===")
    print(f"Result dir: {result_dir}")
    print(f"Parameters: dt={dt}, R={R}, q_cv={q_cv}, q_gm={q_gm}, tau={tau}")
    print(f"Thresholds: max_age={max_age} (model time), max_variance={max_variance}")

    completed_count = 0
    skipped_count = 0
    error_count = 0

    for c_idx, case in enumerate(cases, 1):
        case_id = case["case_id"]
        seed = int(case["seed"])
        stratum = case["stratum"]
        source_rel = case["source"]
        t0_1based = int(case["t0"])
        t0_0based = t0_1based - 1  # 0-based simulation observation step

        # Resolve source MAT path
        source_path = result_dir / source_rel
        if not source_path.exists():
            source_path = Path(source_rel)
        if not source_path.exists():
            source_path = Path.cwd() / source_rel
        if not source_path.exists():
            err_msg = f"Source MAT not found: {source_rel}"
            print(f"[ERROR] Case {case_id}: {err_msg}")
            error_count += len(durations)
            for D in durations:
                manifest_records.append({
                    "case_id": case_id,
                    "seed": seed,
                    "stratum": stratum,
                    "D": D,
                    "status": "ERROR_SOURCE_MISSING",
                    "error": err_msg,
                })
            continue

        source_sha256 = compute_sha256(source_path)
        with open(result_dir / 'input_manifest.json', encoding='utf-8') as f:
            source_manifest = json.load(f)
        expected = next(x for x in source_manifest['inputs'] if x['source'] == source_rel)
        if source_sha256 != expected['sha256']:
            raise ValueError(f"Input SHA mismatch: {source_rel}")

        # Load MAT source data
        try:
            mat = scipy.io.loadmat(str(source_path))
            true_bearing = np.asarray(mat["true_bearing"], dtype=np.float64)
            observed_bearing = np.asarray(mat["observed_bearing"], dtype=np.float64)
            observer_id = np.asarray(mat["observer_id"]).squeeze()
            target_id = np.asarray(mat["target_id"]).squeeze()
        except Exception as e:
            err_msg = f"Failed to load MAT file {source_path}: {e}"
            print(f"[ERROR] Case {case_id}: {err_msg}")
            error_count += len(durations)
            for D in durations:
                manifest_records.append({
                    "case_id": case_id,
                    "seed": seed,
                    "stratum": stratum,
                    "D": D,
                    "status": "ERROR_LOAD_FAILED",
                    "error": err_msg,
                })
            continue

        T_total, E = observed_bearing.shape

        for D in durations:
            out_filename = f"case_{case_id}_D_{D}.npz"
            out_path = out_dir / out_filename

            # Refuse overwrite: preserve existing outputs
            if out_path.exists():
                print(f"[{c_idx}/{len(cases)}] Case {case_id} D={D}: SKIPPED (file exists: {out_filename})")
                skipped_count += 1
                if (case_id, D) not in existing_manifest_map:
                    manifest_records.append({
                        "case_id": case_id,
                        "seed": seed,
                        "stratum": stratum,
                        "D": D,
                        "npz_file": out_filename,
                        "npz_sha256": compute_sha256(out_path),
                        "source_sha256": source_sha256,
                        "status": "EXISTS_SKIPPED",
                    })
                continue

            # Build missing mask for full trajectory [T_total]
            # Artificial open-loop mask hides ALL 140 edges on steps t0:t0+D-1 (1-based)
            occ_start = t0_0based
            occ_end = t0_0based + D
            if occ_start < 0 or occ_end > T_total:
                err_msg = f"Occlusion window [{occ_start}:{occ_end}] out of bounds for T={T_total}"
                print(f"[ERROR] Case {case_id} D={D}: {err_msg}")
                error_count += 1
                manifest_records.append({
                    "case_id": case_id,
                    "seed": seed,
                    "stratum": stratum,
                    "D": D,
                    "status": "ERROR_BOUNDS",
                    "error": err_msg,
                })
                continue

            missing_full = np.zeros(T_total, dtype=bool)
            missing_full[occ_start:occ_end] = True

            # Run estimators on full trajectory prefix from step 1
            # Truth is NEVER passed to the estimators!
            try:
                res_hold = hold(observed_bearing, missing_full, dt=dt)
                res_cv = replay(observed_bearing, missing_full, dt=dt, q=q_cv, R=R, tau=None)
                res_gm = replay(observed_bearing, missing_full, dt=dt, q=q_gm, R=R, tau=tau)
            except Exception as e:
                err_msg = f"Estimator execution error: {e}"
                print(f"[ERROR] Case {case_id} D={D}: {err_msg}")
                error_count += 1
                manifest_records.append({
                    "case_id": case_id,
                    "seed": seed,
                    "stratum": stratum,
                    "D": D,
                    "status": "ERROR_REPLAY_FAILED",
                    "error": err_msg,
                })
                continue

            # Slice window: pre_steps + occ D + recovery_steps
            win_start = occ_start - pre_steps
            win_end = occ_end + recovery_steps
            if win_start < 0 or win_end > T_total:
                err_msg = f"Window [{win_start}:{win_end}] exceeds trajectory bounds [0:{T_total}]"
                print(f"[ERROR] Case {case_id} D={D}: {err_msg}")
                error_count += 1
                manifest_records.append({
                    "case_id": case_id,
                    "seed": seed,
                    "stratum": stratum,
                    "D": D,
                    "status": "ERROR_WINDOW_BOUNDS",
                    "error": err_msg,
                })
                continue

            win_slice = slice(win_start, win_end)
            window_steps = np.arange(win_start + 1, win_end + 1, dtype=np.int32)  # 1-based steps

            missing_win = missing_full[win_slice]  # [W]
            missing_win_2d = np.broadcast_to(missing_win[:, None], (len(window_steps), E))

            truth_win = true_bearing[win_slice]  # [W, E]
            observed_sanitized_win = observed_bearing[win_slice].copy()
            observed_sanitized_win[missing_win, :] = np.nan

            # Active flags and angle errors for Hold
            act_age_hold, act_p_hold, act_hold = compute_active_flags(
                res_hold["age"][win_slice],
                res_hold["p11"][win_slice],
                missing_win_2d,
                max_age=max_age,
                max_variance=max_variance,
                is_hold=True,
            )
            err_hold = wrap_to_pi(res_hold["theta"][win_slice] - truth_win)

            # Active flags and angle errors for KF-CV
            act_age_cv, act_p_cv, act_cv = compute_active_flags(
                res_cv["age"][win_slice],
                res_cv["p11"][win_slice],
                missing_win_2d,
                max_age=max_age,
                max_variance=max_variance,
                is_hold=False,
            )
            err_cv = wrap_to_pi(res_cv["theta"][win_slice] - truth_win)

            # Active flags and angle errors for KF-GM
            act_age_gm, act_p_gm, act_gm = compute_active_flags(
                res_gm["age"][win_slice],
                res_gm["p11"][win_slice],
                missing_win_2d,
                max_age=max_age,
                max_variance=max_variance,
                is_hold=False,
            )
            err_gm = wrap_to_pi(res_gm["theta"][win_slice] - truth_win)

            # Save NPZ file
            save_dict = {
                # Metadata
                "case_id": np.array(case_id),
                "seed": np.array(seed, dtype=np.int64),
                "stratum": np.array(stratum),
                "D": np.array(D, dtype=np.int32),
                "t0": np.array(t0_1based, dtype=np.int32),
                "dt": np.array(dt, dtype=np.float64),
                "window_steps": window_steps,
                "pre_steps": np.array(pre_steps, dtype=np.int32),
                "recovery_steps": np.array(recovery_steps, dtype=np.int32),
                "max_age": np.array(max_age, dtype=np.float64),
                "max_variance": np.array(max_variance, dtype=np.float64),
                "observer_id": observer_id,
                "target_id": target_id,
                "methods": np.array(["Hold", "KF-CV", "KF-GM"]),
                # Shared observations and truth
                "missing": missing_win,
                "truth": truth_win,
                "observed_sanitized": observed_sanitized_win,
                # Hold
                "theta_Hold": res_hold["theta"][win_slice],
                "p11_Hold": res_hold["p11"][win_slice],
                "p12_Hold": res_hold["p12"][win_slice],
                "p22_Hold": res_hold["p22"][win_slice],
                "nis_Hold": res_hold["nis"][win_slice],
                "age_Hold": res_hold["age"][win_slice],
                "active_Hold": act_hold,
                "active_by_age_Hold": act_age_hold,
                "active_by_P_Hold": act_p_hold,
                "error_Hold": err_hold,
                # KF-CV
                "theta_KF_CV": res_cv["theta"][win_slice],
                "p11_KF_CV": res_cv["p11"][win_slice],
                "p12_KF_CV": res_cv["p12"][win_slice],
                "p22_KF_CV": res_cv["p22"][win_slice],
                "nis_KF_CV": res_cv["nis"][win_slice],
                "age_KF_CV": res_cv["age"][win_slice],
                "active_KF_CV": act_cv,
                "active_by_age_KF_CV": act_age_cv,
                "active_by_P_KF_CV": act_p_cv,
                "error_KF_CV": err_cv,
                # KF-GM
                "theta_KF_GM": res_gm["theta"][win_slice],
                "p11_KF_GM": res_gm["p11"][win_slice],
                "p12_KF_GM": res_gm["p12"][win_slice],
                "p22_KF_GM": res_gm["p22"][win_slice],
                "nis_KF_GM": res_gm["nis"][win_slice],
                "age_KF_GM": res_gm["age"][win_slice],
                "active_KF_GM": act_gm,
                "active_by_age_KF_GM": act_age_gm,
                "active_by_P_KF_GM": act_p_gm,
                "error_KF_GM": err_gm,
            }

            np.savez_compressed(out_path, **save_dict)
            npz_sha256 = compute_sha256(out_path)
            completed_count += 1
            print(
                f"[{c_idx}/{len(cases)}] Case {case_id} D={D} -> {out_filename} "
                f"(SHA256: {npz_sha256[:10]}...)"
            )

            manifest_records.append({
                "case_id": case_id,
                "seed": seed,
                "stratum": stratum,
                "D": D,
                "npz_file": out_filename,
                "npz_sha256": npz_sha256,
                "source_file": str(source_path.name),
                "source_sha256": source_sha256,
                "status": "OK",
            })

    # Save manifest
    with open(manifest_path, "w", encoding="utf-8") as f:
        json.dump(
            {
                "design_path": str(design_path),
                "total_runs": total_tasks,
                "completed": completed_count,
                "skipped": skipped_count,
                "errors": error_count,
                "runs": manifest_records,
            },
            f,
            indent=2,
        )

    print(f"=== Replay Finished: {completed_count} completed, {skipped_count} skipped, {error_count} errors ===")
    print(f"Manifest written to: {manifest_path}")
    if error_count or completed_count != total_tasks:
        raise RuntimeError("Incomplete replay; failure records preserved in manifest.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="P2 Open-Loop Estimator Replay")
    parser.add_argument("result_dir", help="Path to experiment result directory containing design.json")
    args = parser.parse_args()
    run_replay_cli(args.result_dir)

#!/usr/bin/env python3
"""
plot_kf_diagnostics.py

Generates the 8 official required diagnostic figures for Flocking KF Diagnostic Experiment (Package A),
plus optional estimator diagnostics.
Complies strictly with Review Fixes 10 & 11:
  - Exact file naming as mandated by project plan.
  - Generates P1_BumpConcentration.png (critical neural field bump concentration).
  - P1_PairDistance.png plots delta relative to same-seed Fresh control.
  - P1_GainScatter.png evaluates primary occlusion window (occ).
  - P1_FailureCases.png uses minimum seed (2026110101) with correct post-step index:
    state after step s corresponds to uRecord[:, :, s] (i.e. uRecord(:,:,s+1) in 1-based MATLAB).
  - Low-memory streaming: only loads necessary scalar timeSeries, only loads neural state for min seed.
  - Seed-level intervals use ddof=1 and t-distribution.
"""

import os
import sys
import glob
import argparse
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.io import loadmat
from scipy import stats

T0 = 1401
DURATIONS = [20, 120]
MODES = ['fresh', 'none', 'hold', 'kf-mean', 'kalman']
COLORS = {
    'fresh': '#2ca02c',   # Green
    'none': '#d62728',    # Red
    'hold': '#ff7f0e',    # Orange
    'kf-mean': '#1f77b4', # Blue
    'kalman': '#9467bd'   # Purple
}
LINESTYLES = {
    'fresh': '--',
    'none': ':',
    'hold': '-.',
    'kf-mean': '-',
    'kalman': '-'
}
LABELS = {
    'fresh': 'Fresh Control',
    'none': 'Drop (None)',
    'hold': 'Hold-Last',
    'kf-mean': 'KF-Mean',
    'kalman': 'KF-Cov (Broadening)'
}

def find_mat_files(results_dir):
    raw_dir = os.path.join(results_dir, "raw")
    if os.path.isdir(raw_dir):
        files = glob.glob(os.path.join(raw_dir, "*.mat"))
        if files:
            return sorted(files)
    return sorted(glob.glob(os.path.join(results_dir, "*.mat")))

def load_timeseries_streaming(results_dir):
    """
    Stream lightweight timeSeries from all MAT files without retaining heavy arrays.
    """
    mat_files = find_mat_files(results_dir)
    if not mat_files:
        raise FileNotFoundError(f"No MAT files found in {results_dir}")

    ts_data = {}
    seeds = set()

    for mf in mat_files:
        d = loadmat(mf, struct_as_record=False, squeeze_me=True)
        res = d.get('result', d.get('saveData', None))
        if res is not None and hasattr(res, 'result'):
            res = res.result
        if res is None:
            continue

        cfg = res.config
        seed = int(cfg.seed)
        dur = int(cfg.occlusionDurationSteps)
        mode_raw = str(cfg.predictionMode)

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

        seeds.add(seed)
        ts = res.timeSeries

        entry = {
            'seed': seed,
            'mode': mode,
            'D': dur,
            'go': np.array(ts.globalOrder, dtype=np.float64),
            'lo': np.array(ts.localOrder, dtype=np.float64),
            'pd_norm': np.array(ts.meanPairDistanceNorm, dtype=np.float64),
            'nd_norm': np.array(ts.meanNearestDistanceNorm, dtype=np.float64),
            'an': np.array(ts.meanActiveNeighbors, dtype=np.float64),
            'vn': np.array(ts.meanVisibleNeighbors, dtype=np.float64) if hasattr(ts, 'meanVisibleNeighbors') else np.full_like(ts.globalOrder, np.nan),
            'im': np.array(ts.meanInputMass, dtype=np.float64),
            'bc': np.array(ts.meanBumpConcentration, dtype=np.float64),
            'sp': np.array(ts.meanSpeed, dtype=np.float64),
            'mat_path': mf
        }
        key = (seed, mode, dur)
        ts_data[key] = entry
        del d
        del res

    sorted_seeds = sorted(list(seeds))
    return ts_data, sorted_seeds

def get_mean_and_ci(ts_data, seeds, mode, dur, field, t_slice):
    data_list = []
    for s in seeds:
        key = (s, 'fresh', 0) if mode == 'fresh' else (s, mode, dur)
        if key in ts_data:
            arr = ts_data[key][field][t_slice]
            data_list.append(arr)

    if not data_list:
        return None, None
    mat = np.array(data_list)
    mean_curve = np.mean(mat, axis=0)
    n = len(data_list)
    if n >= 2:
        se_curve = np.std(mat, axis=0, ddof=1) / np.sqrt(n)
        t_val = stats.t.ppf(0.975, df=n - 1)
        ci_curve = t_val * se_curve
    else:
        ci_curve = np.zeros_like(mean_curve)
    return mean_curve, ci_curve

# Figure 1: P1_Fresh_Drop_Hold_KFmean_KFcov_GO.png
def plot_p1_go(ts_data, seeds, out_dir):
    fig, axes = plt.subplots(2, 1, figsize=(11, 8), sharex=False)
    for idx, dur in enumerate(DURATIONS):
        ax = axes[idx]
        t_start, t_end = T0 - 200, T0 + dur + 400
        steps = np.arange(t_start, t_end + 1)
        t_slice = slice(t_start - 1, t_end)

        for mode in MODES:
            mean_c, ci_c = get_mean_and_ci(ts_data, seeds, mode, dur, 'go', t_slice)
            if mean_c is not None:
                lbl = LABELS[mode]
                ax.plot(steps, mean_c, label=lbl, color=COLORS[mode],
                        linestyle=LINESTYLES[mode], linewidth=2.0)
                if ci_c is not None and np.any(ci_c > 0):
                    ax.fill_between(steps, mean_c - ci_c, mean_c + ci_c,
                                    color=COLORS[mode], alpha=0.12)

        ax.axvline(T0, color='black', linestyle='--', linewidth=1.2, label='Occlusion Start (t0=1401)')
        ax.axvline(T0 + dur, color='black', linestyle=':', linewidth=1.2, label=f'Occlusion End (t0+{dur})')
        ax.axvspan(T0 - 200, T0 - 1, color='gray', alpha=0.08, label='Pre-window (200 steps)')
        ax.axvspan(T0, T0 + dur - 1, color='red', alpha=0.08, label=f'Occlusion ({dur} steps)')
        ax.axvspan(T0 + dur, T0 + dur + 399, color='blue', alpha=0.05, label='Recovery (400 steps)')

        ax.set_title(f"Global Order: D={dur} steps ({dur*0.3:.1f} model time units), n={len(seeds)} seeds", fontsize=12, fontweight='bold')
        ax.set_ylabel("Global Order (GO)", fontsize=11)
        ax.set_xlabel("Simulation step (dt = 0.3 model time units)", fontsize=11)
        ax.grid(True, alpha=0.3)
        if idx == 0:
            ax.legend(loc='lower left', fontsize=8, ncol=2)

    plt.tight_layout()
    p = os.path.join(out_dir, "P1_Fresh_Drop_Hold_KFmean_KFcov_GO.png")
    plt.savefig(p, dpi=300)
    plt.close()
    print(f"Saved: {p}")

# Figure 2: P1_PairDistance.png (Delta relative to same-seed Fresh)
def plot_p1_pair_distance(ts_data, seeds, out_dir):
    fig, axes = plt.subplots(2, 1, figsize=(11, 8), sharex=False)
    for idx, dur in enumerate(DURATIONS):
        ax = axes[idx]
        t_start, t_end = T0 - 200, T0 + dur + 400
        steps = np.arange(t_start, t_end + 1)
        t_slice = slice(t_start - 1, t_end)

        for mode in ['none', 'hold', 'kf-mean', 'kalman']:
            delta_list = []
            for s in seeds:
                key_fresh = (s, 'fresh', 0)
                key_mode = (s, mode, dur)
                if key_fresh in ts_data and key_mode in ts_data:
                    delta = ts_data[key_mode]['pd_norm'][t_slice] - ts_data[key_fresh]['pd_norm'][t_slice]
                    delta_list.append(delta)

            if delta_list:
                mat = np.array(delta_list)
                mean_d = np.mean(mat, axis=0)
                n = len(delta_list)
                se_d = np.std(mat, axis=0, ddof=1) / np.sqrt(n) if n >= 2 else np.zeros_like(mean_d)
                ci_d = stats.t.ppf(0.975, df=n - 1) * se_d if n >= 2 else np.zeros_like(mean_d)

                ax.plot(steps, mean_d, label=LABELS[mode], color=COLORS[mode],
                        linestyle=LINESTYLES[mode], linewidth=2.0)
                ax.fill_between(steps, mean_d - ci_d, mean_d + ci_d, color=COLORS[mode], alpha=0.12)

        ax.axhline(0, color='green', linestyle='--', linewidth=1.5, label='Fresh Baseline (Δ = 0)')
        ax.axvline(T0, color='black', linestyle='--', linewidth=1.2)
        ax.axvline(T0 + dur, color='black', linestyle=':', linewidth=1.2)
        ax.axvspan(T0, T0 + dur - 1, color='red', alpha=0.08)

        ax.set_title(f"Pair Distance Deviation Relative to Same-Seed Fresh: D = {dur} Steps", fontsize=12, fontweight='bold')
        ax.set_ylabel("Δ Pair Distance / L (+ means dispersed)", fontsize=11)
        ax.set_xlabel("Simulation Step", fontsize=11)
        ax.grid(True, alpha=0.3)
        if idx == 0:
            ax.legend(loc='upper left', fontsize=8, ncol=2)

    plt.tight_layout()
    p = os.path.join(out_dir, "P1_PairDistance.png")
    plt.savefig(p, dpi=300)
    plt.close()
    print(f"Saved: {p}")

# Figure 3: P1_LocalOrder.png
def plot_p1_local_order(ts_data, seeds, out_dir):
    fig, axes = plt.subplots(2, 1, figsize=(11, 8), sharex=False)
    for idx, dur in enumerate(DURATIONS):
        ax = axes[idx]
        t_start, t_end = T0 - 200, T0 + dur + 400
        steps = np.arange(t_start, t_end + 1)
        t_slice = slice(t_start - 1, t_end)

        for mode in MODES:
            mean_c, ci_c = get_mean_and_ci(ts_data, seeds, mode, dur, 'lo', t_slice)
            if mean_c is not None:
                ax.plot(steps, mean_c, label=LABELS[mode], color=COLORS[mode],
                        linestyle=LINESTYLES[mode], linewidth=2.0)
                if ci_c is not None and np.any(ci_c > 0):
                    ax.fill_between(steps, mean_c - ci_c, mean_c + ci_c, color=COLORS[mode], alpha=0.12)

        ax.axvline(T0, color='black', linestyle='--', linewidth=1.2)
        ax.axvline(T0 + dur, color='black', linestyle=':', linewidth=1.2)
        ax.axvspan(T0, T0 + dur - 1, color='red', alpha=0.08)
        ax.set_title(f"Local Order (LO, k=5): D = {dur} Steps", fontsize=12, fontweight='bold')
        ax.set_ylabel("Local Order (LO)", fontsize=11)
        ax.set_xlabel("Simulation Step", fontsize=11)
        ax.grid(True, alpha=0.3)
        if idx == 0:
            ax.legend(loc='lower left', fontsize=8, ncol=2)

    plt.tight_layout()
    p = os.path.join(out_dir, "P1_LocalOrder.png")
    plt.savefig(p, dpi=300)
    plt.close()
    print(f"Saved: {p}")

# Figure 4: P1_ActiveNeighbors.png (both fresh and active)
def plot_p1_active_neighbors(ts_data, seeds, out_dir):
    fig, axes = plt.subplots(2, 1, figsize=(11, 8), sharex=False)
    for idx, dur in enumerate(DURATIONS):
        ax = axes[idx]
        t_start, t_end = T0 - 50, T0 + dur + 100
        steps = np.arange(t_start, t_end + 1)
        t_slice = slice(t_start - 1, t_end)

        for mode in MODES:
            mean_an, _ = get_mean_and_ci(ts_data, seeds, mode, dur, 'an', t_slice)
            if mean_an is not None:
                ax.plot(steps, mean_an, label=f"{LABELS[mode]} (Active)", color=COLORS[mode],
                        linestyle='-', linewidth=2.0)

        # Also plot visible/fresh count for reference (drops to 8 during occlusion)
        mean_vn, _ = get_mean_and_ci(ts_data, seeds, 'kalman', dur, 'vn', t_slice)
        if mean_vn is not None and not np.all(np.isnan(mean_vn)):
            ax.plot(steps, mean_vn, label="Fresh Visible Neighbors (k_vis)", color='black',
                    linestyle=':', linewidth=1.8)

        ax.axvline(T0, color='black', linestyle='--', linewidth=1.2)
        ax.axvline(T0 + dur, color='black', linestyle=':', linewidth=1.2)
        ax.set_title(f"Active & Fresh Neighbor Count: D = {dur} Steps", fontsize=12, fontweight='bold')
        ax.set_ylabel("Neighbor Count (Out of k=15)", fontsize=11)
        ax.set_xlabel("Simulation Step", fontsize=11)
        ax.set_ylim(5, 17)
        ax.grid(True, alpha=0.3)
        if idx == 0:
            ax.legend(loc='lower left', fontsize=8, ncol=2)

    plt.tight_layout()
    p = os.path.join(out_dir, "P1_ActiveNeighbors.png")
    plt.savefig(p, dpi=300)
    plt.close()
    print(f"Saved: {p}")

# Figure 5: P1_InputMass.png
def plot_p1_input_mass(ts_data, seeds, out_dir):
    fig, axes = plt.subplots(2, 1, figsize=(11, 8), sharex=False)
    for idx, dur in enumerate(DURATIONS):
        ax = axes[idx]
        t_start, t_end = T0 - 50, T0 + dur + 100
        steps = np.arange(t_start, t_end + 1)
        t_slice = slice(t_start - 1, t_end)

        for mode in MODES:
            mean_c, ci_c = get_mean_and_ci(ts_data, seeds, mode, dur, 'im', t_slice)
            if mean_c is not None:
                ax.plot(steps, mean_c, label=LABELS[mode], color=COLORS[mode],
                        linestyle=LINESTYLES[mode], linewidth=2.0)
                if ci_c is not None and np.any(ci_c > 0):
                    ax.fill_between(steps, mean_c - ci_c, mean_c + ci_c, color=COLORS[mode], alpha=0.12)

        ax.axvline(T0, color='black', linestyle='--', linewidth=1.2)
        ax.axvline(T0 + dur, color='black', linestyle=':', linewidth=1.2)
        ax.set_title(f"Injected Social Input Mass (Integral of u_ext): D = {dur} Steps", fontsize=12, fontweight='bold')
        ax.set_ylabel("Integrated input (model units)", fontsize=11)
        ax.set_xlabel("Simulation Step", fontsize=11)
        ax.grid(True, alpha=0.3)
        if idx == 0:
            ax.legend(loc='lower left', fontsize=8, ncol=2)

    plt.tight_layout()
    p = os.path.join(out_dir, "P1_InputMass.png")
    plt.savefig(p, dpi=300)
    plt.close()
    print(f"Saved: {p}")

# Figure 6: P1_BumpConcentration.png (Mandated Package A Key Figure!)
def plot_p1_bump_concentration(ts_data, seeds, out_dir):
    fig, axes = plt.subplots(2, 1, figsize=(11, 8), sharex=False)
    for idx, dur in enumerate(DURATIONS):
        ax = axes[idx]
        t_start, t_end = T0 - 200, T0 + dur + 400
        steps = np.arange(t_start, t_end + 1)
        t_slice = slice(t_start - 1, t_end)

        for mode in MODES:
            mean_c, ci_c = get_mean_and_ci(ts_data, seeds, mode, dur, 'bc', t_slice)
            if mean_c is not None:
                ax.plot(steps, mean_c, label=LABELS[mode], color=COLORS[mode],
                        linestyle=LINESTYLES[mode], linewidth=2.0)
                if ci_c is not None and np.any(ci_c > 0):
                    ax.fill_between(steps, mean_c - ci_c, mean_c + ci_c, color=COLORS[mode], alpha=0.12)

        ax.axvline(T0, color='black', linestyle='--', linewidth=1.2)
        ax.axvline(T0 + dur, color='black', linestyle=':', linewidth=1.2)
        ax.axvspan(T0, T0 + dur - 1, color='red', alpha=0.08)
        ax.set_title(f"Neural Bump Concentration (Order Vector Length): D = {dur} Steps", fontsize=12, fontweight='bold')
        ax.set_ylabel("Bump Concentration", fontsize=11)
        ax.set_xlabel("Simulation Step", fontsize=11)
        ax.grid(True, alpha=0.3)
        if idx == 0:
            ax.legend(loc='lower left', fontsize=8, ncol=2)

    plt.tight_layout()
    p = os.path.join(out_dir, "P1_BumpConcentration.png")
    plt.savefig(p, dpi=300)
    plt.close()
    print(f"Saved: {p}")

# Figure 7: P1_GainScatter.png (Primary response window: OCCLUSION)
def plot_p1_gain_scatter(results_dir, out_dir):
    contrasts_file = os.path.join(out_dir, "contrasts.csv")
    if not os.path.exists(contrasts_file):
        contrasts_file = os.path.join(results_dir, "analysis", "contrasts.csv")
    if not os.path.exists(contrasts_file):
        contrasts_file = os.path.join(results_dir, "contrasts.csv")
    if not os.path.exists(contrasts_file):
        raise FileNotFoundError("contrasts.csv is required for the prespecified gain plot.")

    import pandas as pd
    df = pd.read_csv(contrasts_file)
    # Filter to primary response: OCCLUSION window
    occ_df = df[df['window'] == 'occ']
    if occ_df.empty:
        occ_df = df

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    ax_pair, ax_go = axes[0], axes[1]

    # Subplot 1: benefitPair in occ window (Hold vs KF-mean, Hold vs Kalman, KF-mean vs Kalman)
    # Subplot 2: benefitGO in occ window
    durations = [20, 120]
    bar_width = 0.25

    methods = ['kf-mean', 'kalman']
    for m_idx, m in enumerate(methods):
        sub_m = occ_df[(occ_df['baseline'] == 'hold') & (occ_df['method'] == m)]
        for d_idx, dur in enumerate(durations):
            sub_d = sub_m[sub_m['D'] == dur]
            if not sub_d.empty:
                x_idx = np.arange(len(sub_d)) + (m_idx * 2 + d_idx) * bar_width * 0.4
                ax_pair.scatter(x_idx, sub_d['benefitPair'], marker='o' if dur==20 else '^',
                                label=f"{m} (D={dur})", color=COLORS[m], alpha=0.7 + 0.2*d_idx)
                ax_go.scatter(x_idx, sub_d['benefitGO'], marker='o' if dur==20 else '^',
                              label=f"{m} (D={dur})", color=COLORS[m], alpha=0.7 + 0.2*d_idx)

    ax_pair.axhline(0, color='black', linestyle='-', linewidth=0.8)
    ax_pair.set_title("benefitPair in Occlusion Window [Hold - Method]\n(+ means Method maintains tighter cohesion)", fontsize=11, fontweight='bold')
    ax_pair.set_ylabel("Pair distance / L advantage (+ is benefit)", fontsize=10)
    ax_pair.set_xlabel("Seed Index", fontsize=10)
    ax_pair.grid(True, alpha=0.3)
    ax_pair.legend(fontsize=8)

    ax_go.axhline(0, color='black', linestyle='-', linewidth=0.8)
    ax_go.set_title("benefitGO in Occlusion Window [Method - Hold]\n(+ means Method maintains higher order)", fontsize=11, fontweight='bold')
    ax_go.set_ylabel("Δ Global Order Advantage (+ is benefit)", fontsize=10)
    ax_go.set_xlabel("Seed Index", fontsize=10)
    ax_go.grid(True, alpha=0.3)
    ax_go.legend(fontsize=8)
    seed_ids=sorted(occ_df['seed'].unique())
    for ax in axes:
        ax.set_xticks(np.arange(len(seed_ids)) + 0.15, [str(s) for s in seed_ids], rotation=45, ha='right')

    plt.tight_layout()
    p = os.path.join(out_dir, "P1_GainScatter.png")
    plt.savefig(p, dpi=300)
    plt.close()
    print(f"Saved: {p}")

# Figure 8: P1_FailureCases.png (Spatial positions, headings, and neural activity with correct post-step index)
def plot_p1_failure_cases(ts_data, seeds, out_dir):
    min_seed = min(seeds)
    # Find kalman D=120 or D=20 run for min_seed
    target_key = (min_seed, 'kalman', 120)
    if target_key not in ts_data:
        target_key = (min_seed, 'kalman', 20)
    if target_key not in ts_data:
        print("Target run for Fig 8 not found; skipping Fig 8.")
        return

    mat_path = ts_data[target_key]['mat_path']
    d = loadmat(mat_path, struct_as_record=False, squeeze_me=True)
    res = d.get('result', d.get('saveData', None))
    if hasattr(res, 'result'):
        res = res.result

    dur = int(res.config.occlusionDurationSteps)
    L = float(res.config.arenaSize)
    N = int(res.config.NAgents)
    nS = int(res.config.NNeurons)

    # 3 target steps:
    # 1. t0 - 1 (pre-occlusion end, step 1400)
    # 2. t0 + dur - 1 (occlusion end, step 1400+dur)
    # 3. t0 + dur + 399 (recovery end, step 1400+dur+400)
    s1 = T0 - 1
    s2 = T0 + dur - 1
    s3 = min(2200, T0 + dur + 399)

    target_steps = [s1, s2, s3]
    step_titles = [f"Pre-Occlusion (Step {s1})", f"End Occlusion (Step {s2})", f"End Recovery (Step {s3})"]

    # Trajectory positions and headings
    traj_x = np.array(res.trajectory.x)
    traj_y = np.array(res.trajectory.y)
    traj_h = np.array(res.trajectory.headings)

    # uRecord post-step index rule:
    # State after step s corresponds to uRecord[:, :, s] in Python 0-based indexing
    # (which was recorded at step s+1 before step s+1 push).
    # If s == 2200, use finalNeuralState.
    has_u_rec = hasattr(res, 'uRecord') and res.uRecord is not None
    final_u = np.array(res.finalNeuralState) if hasattr(res, 'finalNeuralState') else None

    fig, axes = plt.subplots(2, 3, figsize=(15, 8))

    for idx, s in enumerate(target_steps):
        ax_pos = axes[0, idx]
        ax_neu = axes[1, idx]

        # Top row: Spatial positions & headings
        rec_idx = s
        px = traj_x[:, rec_idx]
        py = traj_y[:, rec_idx]
        ph = traj_h[:, rec_idx]

        u_dir = np.cos(ph) * 25
        v_dir = np.sin(ph) * 25

        ax_pos.quiver(px, py, u_dir, v_dir, angles='xy', scale_units='xy', scale=1,
                      color='navy', alpha=0.8, width=0.005)
        ax_pos.scatter(px, py, color='crimson', s=35, zorder=5)
        ax_pos.set_xlim(0, L)
        ax_pos.set_ylim(0, L)
        ax_pos.set_aspect('equal')
        ax_pos.set_title(step_titles[idx], fontsize=11, fontweight='bold')
        ax_pos.set_xlabel("Arena X", fontsize=10)
        if idx == 0:
            ax_pos.set_ylabel("Arena Y", fontsize=10)
        ax_pos.grid(True, alpha=0.3)

        # Bottom row: Neural activity profile
        if s == 2200 and final_u is not None:
            u_step = final_u
        elif has_u_rec and s < res.uRecord.shape[2]:
            u_step = np.array(res.uRecord)[:, :, s]
        else:
            u_step = np.zeros((nS, N))

        # Plot positive neural output tanh(beta * u)+
        pos_out = np.maximum(0, np.tanh(1000.0 * u_step))
        im = ax_neu.imshow(pos_out, aspect='auto', origin='lower', cmap='viridis', vmin=0, vmax=1,
                           extent=(0.5, N+0.5, 0.5, nS+0.5))
        ax_neu.set_xticks([1,5,10,15,20])
        fig.colorbar(im, ax=ax_neu, fraction=0.046, pad=0.04, label='Positive activity')
        ax_neu.set_title(f"Neural Output (tanh(βu)+) at Step {s}", fontsize=10)
        ax_neu.set_xlabel("Agent Index (1..20)", fontsize=10)
        if idx == 0:
            ax_neu.set_ylabel("Neuron Angle Bin (1..100)", fontsize=10)

    plt.suptitle(f"Pre-specified Representative Case (Smallest Seed {min_seed}, KF-Cov, D={dur})",
                 fontsize=13, fontweight='bold', y=0.99)
    plt.tight_layout(rect=[0,0,1,0.95])
    p = os.path.join(out_dir, "P1_FailureCases.png")
    plt.savefig(p, dpi=300)
    plt.close()
    print(f"Saved: {p}")

def generate_all_plots(results_dir, output_dir=None):
    if output_dir is None:
        output_dir = os.path.join(results_dir, "analysis")
    os.makedirs(output_dir, exist_ok=True)
    expected = ['P1_Fresh_Drop_Hold_KFmean_KFcov_GO.png', 'P1_PairDistance.png',
                'P1_LocalOrder.png', 'P1_ActiveNeighbors.png', 'P1_InputMass.png',
                'P1_BumpConcentration.png', 'P1_GainScatter.png', 'P1_FailureCases.png']
    if any(os.path.exists(os.path.join(output_dir, name)) for name in expected):
        raise FileExistsError('Refuse to overwrite existing diagnostic figures.')
    from analyze_kf_diagnostics import inspect_manifest_and_files
    inspect_manifest_and_files(results_dir, find_mat_files(results_dir))
    if not os.path.isfile(os.path.join(output_dir, 'contrasts.csv')):
        raise FileNotFoundError('Run analysis into the same output directory before plotting.')

    print(f"Loading time series data from: {results_dir}")
    ts_data, seeds = load_timeseries_streaming(results_dir)
    print(f"Loaded {len(ts_data)} runs across seeds: {seeds}")

    # Generate the 8 official figures
    plot_p1_go(ts_data, seeds, output_dir)
    plot_p1_pair_distance(ts_data, seeds, output_dir)
    plot_p1_local_order(ts_data, seeds, output_dir)
    plot_p1_active_neighbors(ts_data, seeds, output_dir)
    plot_p1_input_mass(ts_data, seeds, output_dir)
    plot_p1_bump_concentration(ts_data, seeds, output_dir)
    plot_p1_gain_scatter(results_dir, output_dir)
    plot_p1_failure_cases(ts_data, seeds, output_dir)
    print("All 8 official figures generated successfully.")

def main():
    parser = argparse.ArgumentParser(description="Plot 8 official diagnostic figures for Flocking KF Diagnostics.")
    parser.add_argument("results_dir", nargs="?", default=None, help="Path to results directory containing MAT files.")
    parser.add_argument("--output", dest="output_dir", default=None, help="Directory to save figures.")
    args = parser.parse_args()

    if args.results_dir is None:
        base_res = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "effective-social-input", "results"))
        candidates = sorted(glob.glob(os.path.join(base_res, "kf-diagnostic-v1-20261101-*")))
        if candidates:
            target_dir = candidates[-1]
        else:
            print("Error: No results directory found.")
            sys.exit(1)
    else:
        target_dir = os.path.abspath(args.results_dir)

    generate_all_plots(target_dir, args.output_dir)

if __name__ == "__main__":
    main()

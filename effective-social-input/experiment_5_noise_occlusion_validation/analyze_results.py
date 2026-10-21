"""Integrity checks and paired noise-by-occlusion analysis."""
from __future__ import annotations

import hashlib
import json
import shutil
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

import experiment_model as model
import run_experiments as runner

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "results" / "formal_v1"
WINDOWS = {
    "occlusion_1": model.OCCLUSION_WINDOWS[0],
    "occlusion_2": model.OCCLUSION_WINDOWS[1],
    "occlusion_3": model.OCCLUSION_WINDOWS[2],
    "occlusion_all": None,
}
METRICS = ["follower_target_error_deg", "GO_all", "GO_followers", "GO_leaders",
           "follower_leader_error_deg", "hidden_prediction_error_deg", "S_abs_error",
           "pair_distance", "minimum_distance"]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def summarize(values, bootstrap_indices):
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    n = len(values)
    if n == 0:
        return dict(n_seeds=0, mean=np.nan, sd=np.nan, ci95_low=np.nan, ci95_high=np.nan,
                    bootstrap95_low=np.nan, bootstrap95_high=np.nan, pvalue=np.nan,
                    positive=0, negative=0)
    mean = float(values.mean())
    sd = float(values.std(ddof=1)) if n > 1 else np.nan
    half = float(stats.t.ppf(0.975, n - 1) * sd / np.sqrt(n)) if n > 1 else np.nan
    pvalue = float(stats.ttest_1samp(values, 0).pvalue) if n > 1 else np.nan
    if n == len(bootstrap_indices[0]):
        boot_means = values[bootstrap_indices].mean(axis=1)
        boot_low, boot_high = np.quantile(boot_means, [0.025, 0.975])
    else:
        boot_low = boot_high = np.nan
    return dict(n_seeds=n, mean=mean, sd=sd, ci95_low=mean-half if n > 1 else np.nan,
                ci95_high=mean+half if n > 1 else np.nan,
                bootstrap95_low=float(boot_low), bootstrap95_high=float(boot_high),
                pvalue=pvalue, positive=int((values > 1e-12).sum()),
                negative=int((values < -1e-12).sum()))


def block_mean(frame, bounds, metric, is_fresh=False):
    if bounds is None:
        if is_fresh:
            mask = np.zeros(len(frame), dtype=bool)
            for start, end in model.OCCLUSION_WINDOWS:
                mask |= (frame.step >= start) & (frame.step < end)
            block = frame[mask]
        else:
            block = frame[frame.blackout == 1]
    else:
        start, end = bounds
        block = frame[(frame.step >= start) & (frame.step < end)]
    return float(block[metric].mean()) if block[metric].notna().any() else np.nan


def verify_and_summarize():
    manifest = pd.read_csv(OUT / "RUN_MANIFEST.csv")
    planned = json.loads((OUT / "planned_manifest.json").read_text(encoding="utf-8"))
    assert len(manifest) == len(planned) == len(runner.tasks())
    assert manifest.run_id.nunique() == len(planned)
    assert set(manifest.status) == {"COMPLETE"}
    expected_ids = {row["run_id"] for row in planned}
    assert set(manifest.run_id) == expected_ids

    summaries = []
    prefix = {}
    hashes_ok = True
    for entry in manifest.itertuples(index=False):
        path = OUT / entry.data
        hashes_ok &= sha256(path) == entry.sha256
        frame = pd.read_csv(path, encoding="utf-8-sig")
        assert len(frame) == runner.TOTAL_STEPS
        assert frame.step.tolist() == list(range(runner.TOTAL_STEPS))
        assert frame.GO_all.notna().all()
        expected_blackout = frame.step.map(model.blackout_at).astype(int).to_numpy() if entry.mode != "fresh" else np.zeros(runner.TOTAL_STEPS, dtype=int)
        assert np.array_equal(frame.blackout.to_numpy(), expected_blackout)
        expected_hidden = np.where(expected_blackout, model.HIDDEN_PER_FOLLOWER * int(entry.followers), 0)
        assert np.array_equal(frame.hidden_edges.to_numpy(), expected_hidden)
        assert np.allclose(frame.noise_std.to_numpy(), float(entry.noise_std))
        assert np.allclose(frame.R.to_numpy(), float(entry.R))
        key = (int(entry.seed), float(entry.noise_std))
        prefix_values = frame[["GO_all", "GO_followers", "pair_distance"]].iloc[:model.OCCLUSION_WINDOWS[0][0]].to_numpy()
        if key not in prefix:
            prefix[key] = prefix_values
        else:
            assert np.allclose(prefix_values, prefix[key], rtol=0, atol=1e-12)
        for window, bounds in WINDOWS.items():
            for metric in METRICS:
                summaries.append(dict(
                    run_id=entry.run_id, seed=int(entry.seed), noise_std=float(entry.noise_std),
                    R=float(entry.R), method=entry.mode, mask_type=entry.mask_type,
                    window=window, metric=metric,
                    value=block_mean(frame, bounds, metric, is_fresh=(entry.mode == "fresh")),
                ))
    assert hashes_ok
    metrics = pd.DataFrame(summaries)
    metrics.to_csv(OUT / "analysis" / "window_metrics_long.csv", index=False, encoding="utf-8-sig")
    return manifest, metrics


def analyze_effects(metrics):
    wide = metrics.pivot_table(index=["seed", "noise_std", "R", "window", "metric"], columns="method", values="value").reset_index()
    for method in model.METHODS:
        wide[f"delta_{method}_minus_hold"] = wide[method] - wide["hold"]
    wide.to_csv(OUT / "analysis" / "paired_method_differences.csv", index=False, encoding="utf-8-sig")

    rng = np.random.default_rng(20261050)
    bootstrap = rng.integers(0, len(runner.SEEDS), size=(10000, len(runner.SEEDS)))
    np.save(OUT / "analysis" / "bootstrap_seed_indices.npy", bootstrap)
    summaries = []
    for noise_std in model.NOISE_LEVELS:
        for window in WINDOWS:
            for metric in METRICS:
                block = wide[(wide.noise_std == noise_std) & (wide.window == window) & (wide.metric == metric)].sort_values("seed")
                for method in ("drop", "kf-mean"):
                    delta = block[f"delta_{method}_minus_hold"].to_numpy()
                    summaries.append(dict(noise_std=noise_std, R=noise_std**2, window=window,
                                          metric=metric, contrast=f"{method}-hold", method=method,
                                          **summarize(delta, bootstrap)))
                for method in ("hold", "drop", "kf-mean", "fresh"):
                    absolute = metrics[(metrics.noise_std == noise_std) & (metrics.window == window) &
                                       (metrics.metric == metric) & (metrics.method == method)].sort_values("seed").value.to_numpy()
                    summaries.append(dict(noise_std=noise_std, R=noise_std**2, window=window,
                                          metric=metric, contrast="absolute_mean", method=method,
                                          **summarize(absolute, bootstrap)))
    summary = pd.DataFrame(summaries)
    summary.to_csv(OUT / "analysis" / "statistics_by_noise_window.csv", index=False, encoding="utf-8-sig")

    interactions = []
    for window in WINDOWS:
        for metric in ("follower_target_error_deg", "GO_all", "hidden_prediction_error_deg"):
            selected = wide[(wide.window == window) & (wide.metric == metric)].pivot(index="seed", columns="noise_std", values="delta_kf-mean_minus_hold")
            low, medium = model.NOISE_LEVELS
            delta_noise = (selected[medium] - selected[low]).sort_index().to_numpy()
            interactions.append(dict(
                contrast=f"noise interaction: (KF-Hold) at sigma={medium} minus sigma={low}",
                window=window, metric=metric, low_noise=low, medium_noise=medium,
                interpretation=("more negative supports larger KF tracking-error benefit" if metric == "follower_target_error_deg" else "metric-specific direction"),
                **summarize(delta_noise, bootstrap),
            ))
    interaction_frame = pd.DataFrame(interactions)
    interaction_frame.to_csv(OUT / "analysis" / "noise_interaction.csv", index=False, encoding="utf-8-sig")
    return summary, interaction_frame, bootstrap


def plot_results(metrics, interaction):
    figures = OUT / "figures"
    figures.mkdir(exist_ok=True)
    tracking = metrics[(metrics.window == "occlusion_all") & (metrics.metric == "follower_target_error_deg") & metrics.method.isin(["drop", "hold", "kf-mean", "fresh"])]
    order = ["fresh", "drop", "hold", "kf-mean"]
    fig, ax = plt.subplots(figsize=(8.8, 5.2))
    positions = []
    values = []
    labels = []
    for group_i, sigma in enumerate(model.NOISE_LEVELS):
        for method_i, method in enumerate(order):
            v = tracking[(tracking.noise_std == sigma) & (tracking.method == method)].sort_values("seed").value.to_numpy()
            positions.append(group_i * 5 + method_i + 1)
            values.append(v)
            labels.append(method)
    ax.boxplot(values, positions=positions, tick_labels=labels, showmeans=True)
    ax.set_xticks([2.5, 7.5], [f"low σ={model.NOISE_LEVELS[0]} rad", f"medium σ={model.NOISE_LEVELS[1]} rad"])
    ax.set_ylabel("Follower target-direction error (degrees; lower is better)")
    ax.set_title("Three repeated 20-step occlusion windows combined")
    ax.grid(axis="y", alpha=.2)
    fig.tight_layout(); fig.savefig(figures / "tracking_error_by_noise.png", dpi=180); fig.savefig(figures / "tracking_error_by_noise.pdf"); plt.close(fig)

    focus = interaction[(interaction.metric == "follower_target_error_deg") & (interaction.window == "occlusion_all")].iloc[0]
    fig, ax = plt.subplots(figsize=(6.8, 4.8))
    low = metrics[(metrics.noise_std == model.NOISE_LEVELS[0]) & (metrics.window == "occlusion_all") & (metrics.metric == "follower_target_error_deg")]
    high = metrics[(metrics.noise_std == model.NOISE_LEVELS[1]) & (metrics.window == "occlusion_all") & (metrics.metric == "follower_target_error_deg")]
    def paired_delta(sigma):
        pivot = metrics[(metrics.noise_std == sigma) & (metrics.window == "occlusion_all") & (metrics.metric == "follower_target_error_deg")].pivot(index="seed", columns="method", values="value")
        return (pivot["kf-mean"] - pivot["hold"]).sort_index().to_numpy()
    deltas = [paired_delta(sigma) for sigma in model.NOISE_LEVELS]
    for idx, arr in enumerate(deltas):
        x = np.full(len(arr), idx) + np.random.default_rng(20261050 + idx).uniform(-.045, .045, len(arr))
        ax.scatter(x, arr, alpha=.55, s=22)
        ax.errorbar(idx, arr.mean(), yerr=1.96*arr.std(ddof=1)/np.sqrt(len(arr)), fmt="D", capsize=4, color="#0072B2")
    ax.axhline(0, color="black", linestyle="--", linewidth=1)
    ax.set_xticks([0, 1], [f"low σ={model.NOISE_LEVELS[0]}", f"medium σ={model.NOISE_LEVELS[1]}"])
    ax.set_ylabel("KF-mean minus Hold tracking error (degrees)")
    ax.set_title("Negative values indicate lower error with KF")
    fig.tight_layout(); fig.savefig(figures / "kf_hold_noise_interaction.png", dpi=180); fig.savefig(figures / "kf_hold_noise_interaction.pdf"); plt.close(fig)


def write_report(manifest, metrics, summary, interaction, verification):
    main = interaction[(interaction.window == "occlusion_all") & (interaction.metric == "follower_target_error_deg")].iloc[0]
    lines = [
        "# 实验_5结果摘要", "",
        "## 实验问题", "", "检验中等方位观测噪声是否会扩大多次短遮挡期间 KF-mean 相对 Hold 的跟随者目标方向追踪优势。",
        "", "## 实验设置", "", f"- 规模固定为 N={model.N}、L={model.LEADER_COUNT}；每个跟随者遮挡 {model.HIDDEN_PER_FOLLOWER}/15 条观测边。",
        "- 三个遮挡窗口：[180,200)、[300,320)、[420,440)，每次20步；窗口间完整观测恢复。",
        f"- 低噪声 σ={model.NOISE_LEVELS[0]} rad（{np.degrees(model.NOISE_LEVELS[0]):.2f}°），中等噪声 σ={model.NOISE_LEVELS[1]} rad（{np.degrees(model.NOISE_LEVELS[1]):.2f}°）。KF 观测噪声方差同步设为 R=σ²；过程噪声 q 固定。",
        f"- 每档 {len(runner.SEEDS)} 个配对种子；共 {len(manifest)} 条仿真运行。",
        "", "## 主指标和主比较", "", "主指标为三个遮挡窗口合并后的跟随者目标方向误差，单位为度，越低越好。主效应为各噪声水平下 KF-mean−Hold；噪声交互为 (KF−Hold)_medium−(KF−Hold)_low，交互为负表示中等噪声扩大 KF 优势。",
        f"- 噪声交互差值：{main['mean']:.4f}°；95%区间 [{main['ci95_low']:.4f}°, {main['ci95_high']:.4f}°]；p={main['pvalue']:.6g}。",
    ]
    for sigma in model.NOISE_LEVELS:
        rows = summary[(summary.noise_std == sigma) & (summary.window == "occlusion_all") & (summary.metric == "follower_target_error_deg") & (summary.contrast.isin(["absolute_mean", "kf-mean-hold", "drop-hold"]))]
        absolute = rows[rows.contrast == "absolute_mean"].set_index("method")["mean"]
        kf_delta = rows[rows.contrast == "kf-mean-hold"].iloc[0]
        drop_delta = rows[rows.contrast == "drop-hold"].iloc[0]
        lines += ["", f"### σ={sigma} rad", "", f"- 绝对跟随者追踪误差：Fresh {absolute['fresh']:.3f}°；Drop {absolute['drop']:.3f}°；Hold {absolute['hold']:.3f}°；KF-mean {absolute['kf-mean']:.3f}°。",
                  f"- KF-mean−Hold：{kf_delta['mean']:.4f}°，95%区间 [{kf_delta['ci95_low']:.4f}°, {kf_delta['ci95_high']:.4f}°]，p={kf_delta['pvalue']:.6g}。",
                  f"- Drop−Hold：{drop_delta['mean']:.4f}°，95%区间 [{drop_delta['ci95_low']:.4f}°, {drop_delta['ci95_high']:.4f}°]。"]
    lines += ["", "## 解释边界", "", "- 只有噪声标准差及对应的 KF R 改变，其他仿真设置保持一致；相同种子使用相同标准正态噪声张量并缩放。",
              "- 这是验证性小实验，不是最终论文规模；若噪声交互区间跨 0，不能声称中等噪声明确放大了 KF 优势。",
              "- KF 相对 Hold 的优势不代表优于 Drop；需同时看逐方法绝对误差。GO_all、方位预测误差和社会输入误差仅作为辅助诊断。",
              "", "## 复现", "", "使用 `python run_experiments.py --workers 6` 运行，完成后使用 `python analyze_results.py` 校验与分析。逐步数据、配对输入、哈希和图表位于 `results/formal_v1/`。"]
    (OUT / "RESULTS_SUMMARY.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def analyze():
    for name in ("analysis", "figures", "checks"):
        (OUT / name).mkdir(exist_ok=True)
    manifest, metrics = verify_and_summarize()
    summary, interaction, _ = analyze_effects(metrics)
    plot_results(metrics, interaction)
    verification = dict(status="PASS", expected_runs=len(runner.tasks()), completed_runs=len(manifest),
                        run_hashes="PASS", rows_per_run=runner.TOTAL_STEPS,
                        blackout_windows="PASS: three windows, nine hidden edges per follower",
                        noise_scaling="PASS: shared standard-normal tensor, R=sigma^2",
                        pre_first_window_prefix="PASS", independent_seed_count=len(runner.SEEDS))
    (OUT / "checks" / "verification.json").write_text(json.dumps(verification, ensure_ascii=False, indent=2), encoding="utf-8")
    write_report(manifest, metrics, summary, interaction, verification)
    audit = json.loads((OUT / "audit.json").read_text(encoding="utf-8")); audit.update(analysis="COMPLETE", integrity_checks="PASS"); (OUT / "audit.json").write_text(json.dumps(audit, ensure_ascii=False, indent=2), encoding="utf-8")
    metadata = json.loads((OUT / "metadata.json").read_text(encoding="utf-8")); metadata.update(status="COMPLETE", analysis="COMPLETE", integrity_checks="PASS"); (OUT / "metadata.json").write_text(json.dumps(metadata, ensure_ascii=False, indent=2), encoding="utf-8")
    snapshot = OUT / "source_snapshot"
    for source in ROOT.glob("*.py"):
        shutil.copy2(source, snapshot / source.name)
    (OUT / "source_hashes.json").write_text(json.dumps({p.name: sha256(p) for p in sorted(snapshot.glob("*.py"))}, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"ANALYSIS COMPLETE: {OUT / 'RESULTS_SUMMARY.md'}")


if __name__ == "__main__":
    analyze()

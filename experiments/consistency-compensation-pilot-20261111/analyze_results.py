#!/usr/bin/env python3
"""Validate and analyze the paired intermittent-occlusion consistency pilot."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.io import loadmat


NUMERIC = (
    "S0",
    "C_occ_tail",
    "GO_occ_tail",
    "LO_occ_tail",
    "C_all_tail",
    "GO_all_tail",
    "LO_all_tail",
    "C_occ_all",
    "GO_occ_all",
    "LO_occ_all",
    "bearingMAE",
    "inputRMSE",
)
CONSISTENCY = (
    "C_occ_tail",
    "GO_occ_tail",
    "LO_occ_tail",
    "C_all_tail",
    "GO_all_tail",
    "LO_all_tail",
    "C_occ_all",
    "GO_occ_all",
    "LO_occ_all",
)
REQUIRED = ("seed", "split", "schedule", "remaining", "method", *NUMERIC)
PRIMARY = "C_occ_tail"
REPEATS = 4000


def write_csv(path: Path, rows: list[dict], columns: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)


def write_json(path: Path, value: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def comparable(value):
    if isinstance(value, np.ndarray):
        return [comparable(item) for item in value.flat]
    if isinstance(value, np.generic):
        return comparable(value.item())
    return value


def same_value(actual, expected) -> bool:
    if isinstance(actual, list) and isinstance(expected, list):
        return len(actual) == len(expected) and all(
            same_value(left, right) for left, right in zip(actual, expected)
        )
    if isinstance(actual, (int, float)) and isinstance(expected, (int, float)):
        return math.isclose(actual, expected, rel_tol=1e-12, abs_tol=1e-12)
    return actual == expected


def audit_forecast_params(root: Path, seeds: list[int]) -> dict:
    param_file = root / "results" / "test" / "forecast_params.json"
    expected = json.loads(param_file.read_text(encoding="utf-8"))
    before_hash = file_sha256(param_file)
    for seed in seeds:
        path = root / "results" / "test" / f"seed-{seed}" / "traces.mat"
        embedded = loadmat(path, variable_names=["params"], squeeze_me=True,
                           struct_as_record=False).get("params")
        if embedded is None or not hasattr(embedded, "_fieldnames"):
            raise ValueError(f"{path}: missing embedded calibration parameters")
        actual = {name: comparable(getattr(embedded, name)) for name in embedded._fieldnames}
        if set(actual) != set(expected) or not all(
            same_value(actual[name], expected[name]) for name in expected
        ):
            raise ValueError(f"{path}: embedded parameters differ from forecast_params.json")
    if file_sha256(param_file) != before_hash:
        raise ValueError("forecast_params.json changed during test audit")
    return {
        "verified_seed_count": len(seeds),
        "verified_seeds": seeds,
        "sha256": {
            "config.json": file_sha256(root / "config.json"),
            "runConsistencyCompensation2.m": file_sha256(root / "runConsistencyCompensation2.m"),
            "results/test/forecast_params.json": before_hash,
            "results/test/completion.json": file_sha256(root / "results" / "test" / "completion.json"),
        },
    }


def audit_raw_traces(root: Path, seeds: list[int], cells: dict, config: dict,
                     schedule: dict) -> int:
    period = np.r_[np.ones(schedule["occlusion_steps"], dtype=bool),
                   np.zeros(schedule["refresh_steps"], dtype=bool)]
    expected_mask = np.tile(period, schedule["cycles"])
    total_cells = 0
    for seed in seeds:
        path = root / "results" / "test" / f"seed-{seed}" / "traces.mat"
        traces = loadmat(path, variable_names=["traces"], squeeze_me=True,
                         struct_as_record=False).get("traces")
        if traces is None or np.size(traces) != len(config["remaining_counts"]) * len(config["methods"]):
            raise ValueError(f"{path}: missing or incomplete raw traces")
        observed = set()
        for item in np.ravel(traces):
            remaining, method = int(item.remaining), str(item.method)
            key = (seed, remaining, method)
            if key in observed or key not in cells:
                raise ValueError(f"{path}: duplicated or unplanned trace cell {key}")
            observed.add(key)
            series = np.asarray(item.series, dtype=float)
            if series.shape != (len(expected_mask), 10):
                raise ValueError(f"{path}: invalid series shape for {key}")
            if not np.array_equal(series[:, 9], expected_mask.astype(float)):
                raise ValueError(f"{path}: occlusion timing differs for {key}")
            if not np.allclose(series[:, 2], np.minimum(series[:, 0], series[:, 1]),
                               rtol=0, atol=1e-12):
                raise ValueError(f"{path}: C differs from min(GO, LO) for {key}")
            tail = np.zeros(len(expected_mask), dtype=bool)
            tail[-config["tail_steps"]:] = True
            occ_tail = expected_mask & tail
            raw = {
                "C_occ_tail": series[occ_tail, 2].mean(),
                "GO_occ_tail": series[occ_tail, 0].mean(),
                "LO_occ_tail": series[occ_tail, 1].mean(),
                "C_all_tail": series[tail, 2].mean(),
                "GO_all_tail": series[tail, 0].mean(),
                "LO_all_tail": series[tail, 1].mean(),
                "C_occ_all": series[expected_mask, 2].mean(),
                "GO_occ_all": series[expected_mask, 0].mean(),
                "LO_occ_all": series[expected_mask, 1].mean(),
                "inputRMSE": np.sqrt(np.mean(series[expected_mask, 8] ** 2)),
            }
            bearings = series[expected_mask, 7]
            raw["bearingMAE"] = bearings[np.isfinite(bearings)].mean() if np.isfinite(bearings).any() else math.nan
            summary = cells[key]
            for metric, value in raw.items():
                if not (math.isnan(value) and math.isnan(summary[metric])) and not same_value(
                    float(value), summary[metric]
                ):
                    raise ValueError(f"{path}: {metric} differs from CSV for {key}")
            total_cells += 1
    return total_cells


def parse_file(path: Path, seed: int, split: str, schedule: str, config: dict) -> list[dict]:
    if not path.is_file():
        raise FileNotFoundError(f"Missing required seed summary: {path}")
    if not (path.parent / "traces.mat").is_file():
        raise FileNotFoundError(f"Missing required raw traces beside {path}")
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        missing = set(REQUIRED) - set(reader.fieldnames or ())
        if missing:
            raise ValueError(f"{path}: missing columns {sorted(missing)}")
        rows = []
        for line, raw in enumerate(reader, 2):
            try:
                row = {
                    "seed": int(raw["seed"]),
                    "split": raw["split"].strip().lower(),
                    "schedule": raw["schedule"].strip(),
                    "remaining": int(raw["remaining"]),
                    "method": raw["method"].strip(),
                }
                row.update({
                    name: float(raw[name]) if raw[name].strip() else math.nan
                    for name in NUMERIC
                })
            except (TypeError, ValueError) as exc:
                raise ValueError(f"{path}:{line}: invalid summary value") from exc
            if row["seed"] != seed or row["split"] != split or row["schedule"] != schedule:
                raise ValueError(f"{path}:{line}: seed, split, or schedule does not match its file")
            if row["remaining"] not in config["remaining_counts"]:
                raise ValueError(f"{path}:{line}: unplanned remaining count")
            if row["method"] not in config["methods"]:
                raise ValueError(f"{path}:{line}: unplanned method")
            if not all(math.isfinite(row[name]) and -1e-9 <= row[name] <= 1 + 1e-9 for name in CONSISTENCY):
                raise ValueError(f"{path}:{line}: consistency outcome outside [0,1]")
            if not math.isfinite(row["S0"]) or row["S0"] < -1e-12:
                raise ValueError(f"{path}:{line}: invalid S0")
            strength_bound = row["remaining"] * config["model"]["single_edge_weight"]
            if row["S0"] > strength_bound + 1e-9:
                raise ValueError(f"{path}:{line}: S0 exceeds visible total-strength bound")
            if not math.isfinite(row["inputRMSE"]) or row["inputRMSE"] < 0:
                raise ValueError(f"{path}:{line}: invalid inputRMSE")
            if row["method"] in ("Hold", "Linear", "DampedHolt", "SG", "KF") and (
                not math.isfinite(row["bearingMAE"])
                or not 0 <= row["bearingMAE"] <= math.pi + 1e-9
            ):
                raise ValueError(f"{path}:{line}: invalid prediction bearingMAE")
            rows.append(row)
    return rows


def validate_run_folder(folder: Path, config: dict, expected_seeds: int, schedule: str) -> None:
    frozen = json.loads((folder / "config.json").read_text(encoding="utf-8"))
    if frozen != config:
        raise ValueError(f"{folder}: frozen result config differs from the analysis config")
    completion = json.loads((folder / "completion.json").read_text(encoding="utf-8"))
    if (completion.get("status") != "complete"
            or completion.get("seed_count") != expected_seeds
            or completion.get("schedule") != schedule):
        raise ValueError(f"{folder}: completion manifest is incomplete or inconsistent")


def validate_cells(rows: list[dict], seeds: list[int], remaining: list[int], methods: list[str], label: str) -> dict:
    keys = [(row["seed"], row["remaining"], row["method"]) for row in rows]
    if len(keys) != len(set(keys)):
        raise ValueError(f"{label}: duplicate seed/remaining/method cells")
    expected = {(seed, count, method) for seed in seeds for count in remaining for method in methods}
    observed = set(keys)
    if observed != expected:
        raise ValueError(f"{label}: missing={len(expected-observed)}, extra={len(observed-expected)}")
    by_cell = {(row["seed"], row["remaining"], row["method"]): row for row in rows}
    for seed in seeds:
        fresh = [by_cell[(seed, count, "Fresh")] for count in remaining]
        for metric in CONSISTENCY:
            if max(row[metric] for row in fresh) - min(row[metric] for row in fresh) > 1e-10:
                raise ValueError(f"{label}: Fresh {metric} varies across masks for seed {seed}")
        for count in remaining:
            s0 = [by_cell[(seed, count, method)]["S0"] for method in methods]
            if max(s0) - min(s0) > 1e-10:
                raise ValueError(f"{label}: first-event S0 varies across methods for seed {seed}, r={count}")
    return by_cell


def bootstrap(values: list[float], sample_indices: np.ndarray) -> tuple[float, float, float]:
    data = np.asarray(values, dtype=float)
    if not np.isfinite(data).all() or data.size != sample_indices.shape[1]:
        raise ValueError("Bootstrap requires one finite paired value per independent seed")
    draws = data[sample_indices].mean(axis=1)
    low, high = np.quantile(draws, [0.025, 0.975])
    return float(data.mean()), float(low), float(high)


def development_rows(root: Path, config: dict) -> tuple[dict, list[dict]]:
    first = config["seed_first"]
    seeds = list(range(first, first + config["development_seeds"]))
    all_rows = []
    maps = {}
    for schedule in config["schedule_candidates"]:
        name = schedule["name"]
        folder = root / "results" / "development" / name
        validate_run_folder(folder, config, len(seeds), name)
        rows = []
        for seed in seeds:
            path = folder / f"seed-{seed}" / "summary.csv"
            rows.extend(parse_file(path, seed, "development", name, config))
        maps[name] = validate_cells(rows, seeds, config["remaining_counts"], ["Fresh", "Drop"], name)
        all_rows.extend(rows)
    return maps, all_rows


def run_gate(root: Path, config: dict) -> dict:
    maps, rows = development_rows(root, config)
    seeds = list(range(config["seed_first"], config["seed_first"] + config["development_seeds"]))
    rng = np.random.default_rng(2026111101)
    indices = rng.integers(0, len(seeds), size=(REPEATS, len(seeds)))
    gate_rows = []
    per_seed = []
    passing = []
    for schedule in config["schedule_candidates"]:
        name = schedule["name"]
        cell = maps[name]
        for remaining in config["remaining_counts"]:
            delta = []
            for seed in seeds:
                fresh = cell[(seed, remaining, "Fresh")]
                drop = cell[(seed, remaining, "Drop")]
                value = fresh[PRIMARY] - drop[PRIMARY]
                delta.append(value)
                per_seed.append({
                    "seed": seed, "schedule": name, "remaining": remaining,
                    "S0": drop["S0"], "Fresh_C_occ_tail": fresh[PRIMARY],
                    "Drop_C_occ_tail": drop[PRIMARY], "Fresh_minus_Drop": value,
                })
            mean, low, high = bootstrap(delta, indices)
            positive = sum(value > 0 for value in delta)
            is_primary = remaining == config["primary_remaining"]
            passed = (
                is_primary
                and mean >= config["gate"]["min_mean_delta"]
                and positive >= config["gate"]["min_positive_seeds"]
            )
            gate_rows.append({
                "schedule": name, "remaining": remaining, "n_seeds": len(seeds),
                "mean_Fresh_minus_Drop": mean, "ci95_low": low, "ci95_high": high,
                "n_positive": positive, "passes_development_gate": passed,
            })
            if passed:
                passing.append(schedule)
    selected = min(passing, key=lambda item: item["occlusion_steps"])["name"] if passing else None
    gate = {
        "status": "pass" if selected else "fail",
        "selected_schedule": selected,
        "primary_remaining": config["primary_remaining"],
        "development_seed_count": len(seeds),
        "criteria": config["gate"],
        "metric": PRIMARY,
        "rows": len(rows),
        "schedule_results": gate_rows,
        "note": "Selection uses Fresh-Drop only; no compensation performance was inspected.",
    }
    out = root / "results" / "development" / "analysis"
    write_csv(out / "gate_summary.csv", gate_rows, list(gate_rows[0]))
    write_csv(out / "gate_by_seed.csv", per_seed, list(per_seed[0]))
    write_json(out / "gate.json", gate)
    lines = [
        "# 开发集：遮挡退化门槛",
        "",
        f"已检查 {len(seeds)} 个开发种子、{len(rows)} 条 Fresh/Drop 记录。只用 r={config['primary_remaining']} 选择时序；r=8 是描述性对照。",
        "C_occ_tail 是最后 1000 步中所有遮挡步的逐步 min(GO, LO) 均值，越高越好。",
        "开发集判据：Fresh−Drop 均值至少 0.05，且至少 5/6 种子为正；若两者都合格，取较短遮挡窗。",
        "",
        "| 时序 | r | Fresh−Drop 均值 [逐项 95% CI] | 正值种子 | 通过主门槛 |",
        "|---|---:|---:|---:|---|",
    ]
    for row in gate_rows:
        lines.append(
            f"| {row['schedule']} | {row['remaining']} | {row['mean_Fresh_minus_Drop']:+.4f} "
            f"[{row['ci95_low']:+.4f}, {row['ci95_high']:+.4f}] | "
            f"{row['n_positive']}/{row['n_seeds']} | {'是' if row['passes_development_gate'] else '否'} |"
        )
    lines += [
        "",
        f"锁定时序：**{selected or '无'}**。这只是开发集筛选；仍需 24 个独立测试种子确认 Fresh−Drop 退化。",
        "队友的 0.06–0.09 是永久遮挡设置下的候选敏感区域，不是本周期遮挡实验已验证的阈值。",
        "",
    ]
    (out / "README.md").write_text("\n".join(lines), encoding="utf-8")
    return gate


def test_rows(root: Path, config: dict, schedule: str) -> tuple[dict, list[dict]]:
    first = config["seed_first"] + config["development_seeds"]
    last = config["seed_first"] + config["seed_count"]
    seeds = list(range(first, last))
    folder = root / "results" / "test"
    validate_run_folder(folder, config, len(seeds), schedule)
    rows = []
    for seed in seeds:
        path = folder / f"seed-{seed}" / "summary.csv"
        rows.extend(parse_file(path, seed, "test", schedule, config))
    cell = validate_cells(rows, seeds, config["remaining_counts"], config["methods"], "test")
    return cell, rows


def plot_primary(path: Path, summaries: list[dict], remaining_counts: list[int], methods: list[str]) -> None:
    fig, axes = plt.subplots(1, len(remaining_counts), figsize=(max(8, 5 * len(remaining_counts)), 5), squeeze=False)
    color = {
        "Fresh": "#303030", "Hold": "#cb6d24", "Linear": "#007f87",
        "DampedHolt": "#786099", "SG": "#7b873d", "KF": "#2168ad",
        "StrengthOnly": "#ad558e",
    }
    for ax, remaining in zip(axes[0], remaining_counts):
        data = [row for row in summaries if row["remaining"] == remaining and row["metric"] == PRIMARY
                and row["reference"] == "Drop" and row["method"] in methods]
        data.sort(key=lambda row: methods.index(row["method"]))
        y = np.arange(len(data))
        means = np.asarray([row["mean_delta"] for row in data])
        lo = np.asarray([row["ci95_low"] for row in data])
        hi = np.asarray([row["ci95_high"] for row in data])
        for i, row in enumerate(data):
            ax.errorbar(row["mean_delta"], y[i],
                        xerr=np.array([[row["mean_delta"]-row["ci95_low"]],
                                       [row["ci95_high"]-row["mean_delta"]]]),
                        fmt="o", capsize=3, color=color.get(row["method"], "#444444"))
        ax.axvline(0, color="#555555", linewidth=1)
        ax.set(yticks=y, yticklabels=[row["method"] for row in data],
               xlabel="Paired C difference vs Drop (95% CI)", title=f"Visible edges r={remaining}")
        ax.invert_yaxis()
        ax.grid(axis="x", alpha=0.25)
    fig.suptitle("Final-window occlusion consistency: paired test seeds", fontsize=12)
    fig.tight_layout()
    fig.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def run_test(root: Path, config: dict) -> dict:
    gate_file = root / "results" / "development" / "analysis" / "gate.json"
    gate = json.loads(gate_file.read_text(encoding="utf-8"))
    if gate.get("status") != "pass" or not gate.get("selected_schedule"):
        raise ValueError("Development gate did not pass; test analysis is not authorized by the protocol")
    schedule = gate["selected_schedule"]
    if schedule not in {item["name"] for item in config["schedule_candidates"]}:
        raise ValueError("Gate selected an unknown schedule")
    cells, rows = test_rows(root, config, schedule)
    first = config["seed_first"] + config["development_seeds"]
    seeds = list(range(first, config["seed_first"] + config["seed_count"]))
    provenance = audit_forecast_params(root, seeds)
    selected = next(item for item in config["schedule_candidates"] if item["name"] == schedule)
    raw_cells_verified = audit_raw_traces(root, seeds, cells, config, selected)
    rng = np.random.default_rng(2026111102)
    indices = rng.integers(0, len(seeds), size=(REPEATS, len(seeds)))

    condition_rows = []
    for remaining in config["remaining_counts"]:
        for method in config["methods"]:
            records = [cells[(seed, remaining, method)] for seed in seeds]
            condition_rows.append({
                "remaining": remaining, "method": method, "n_seeds": len(seeds),
                "S0_mean": float(np.mean([record["S0"] for record in records])),
                **{name + "_mean": float(np.mean([record[name] for record in records])) for name in CONSISTENCY},
                "bearingMAE_mean": float(np.nanmean([record["bearingMAE"] for record in records]))
                if any(math.isfinite(record["bearingMAE"]) for record in records) else math.nan,
                "inputRMSE_mean": float(np.nanmean([record["inputRMSE"] for record in records]))
                if any(math.isfinite(record["inputRMSE"]) for record in records) else math.nan,
            })

    contrast_pairs = [(method, "Drop") for method in config["methods"] if method != "Drop"]
    contrast_pairs += [("KF", name) for name in ("Hold", "Linear", "DampedHolt", "SG", "StrengthOnly")]
    paired_seed = []
    summaries = []
    for remaining in config["remaining_counts"]:
        for method, reference in contrast_pairs:
            for metric in CONSISTENCY:
                deltas = [cells[(seed, remaining, method)][metric]
                          - cells[(seed, remaining, reference)][metric] for seed in seeds]
                mean, low, high = bootstrap(deltas, indices)
                summaries.append({
                    "remaining": remaining, "method": method, "reference": reference,
                    "metric": metric, "n_seeds": len(seeds),
                    "mean_delta": mean, "ci95_low": low, "ci95_high": high,
                    "n_positive": sum(value > 0 for value in deltas),
                })
                for seed, delta in zip(seeds, deltas):
                    paired_seed.append({
                        "seed": seed, "remaining": remaining, "S0": cells[(seed, remaining, "Drop")]["S0"],
                        "method": method, "reference": reference, "metric": metric, "delta": delta,
                    })

    main_damage = next(row for row in summaries if row["remaining"] == config["primary_remaining"]
                       and row["method"] == "Fresh" and row["reference"] == "Drop" and row["metric"] == PRIMARY)
    degradation_confirmed = (
        main_damage["mean_delta"] >= config["gate"]["min_mean_delta"]
        and main_damage["ci95_low"] > 0
    )
    out = root / "results" / "test" / "analysis"
    write_csv(out / "condition_means.csv", condition_rows, list(condition_rows[0]))
    write_csv(out / "paired_by_seed.csv", paired_seed, list(paired_seed[0]))
    write_csv(out / "paired_summary.csv", summaries, list(summaries[0]))
    plot_primary(out / "paired_primary.png", summaries, config["remaining_counts"],
                 [name for name in config["methods"] if name != "Drop"])
    audit = {
        "status": "COMPLETE",
        "selected_schedule": schedule,
        "development_gate": gate["status"],
        "expected_test_seeds": config["seed_count"] - config["development_seeds"],
        "test_seeds": len(seeds),
        "rows": len(rows),
        "bootstrap_repeats": REPEATS,
        "primary_remaining": config["primary_remaining"],
        "primary_metric": PRIMARY,
        "degradation_confirmed": degradation_confirmed,
        "primary_Fresh_minus_Drop": main_damage,
        "inference": "Exploratory pointwise paired bootstrap intervals; no multiple-comparison adjustment.",
        "provenance": provenance,
        "raw_trace_cells_verified": raw_cells_verified,
    }
    write_json(out / "audit.json", audit)

    lines = [
        "# 30 种子周期遮挡一致性实验：独立测试集",
        "",
        f"开发集锁定时序 {schedule}。本报告只分析 {len(seeds)} 个独立测试种子，另 6 个开发种子不计入效应区间。",
        "主指标是在最后 1000 步的遮挡时刻，对每步 min(GO, LO) 求均值，越高越好。",
        f"数据审计为 **COMPLETE**；主条件 r={config['primary_remaining']} 的退化确认：**{'是' if degradation_confirmed else '否'}**。",
        "",
        "## Fresh 与 Drop 的退化检查",
        "",
        "| 保留边数 | 实际 S0 均值 | Fresh−Drop 主指标 [95% CI] | 正值种子 |",
        "|---:|---:|---:|---:|",
    ]
    for remaining in config["remaining_counts"]:
        damage = next(row for row in summaries if row["remaining"] == remaining and row["method"] == "Fresh"
                      and row["reference"] == "Drop" and row["metric"] == PRIMARY)
        s0 = next(row["S0_mean"] for row in condition_rows if row["remaining"] == remaining and row["method"] == "Drop")
        lines.append(
            f"| {remaining} | {s0:.5f} | {damage['mean_delta']:+.4f} "
            f"[{damage['ci95_low']:+.4f}, {damage['ci95_high']:+.4f}] | "
            f"{damage['n_positive']}/{damage['n_seeds']} |"
        )
    lines += [
        "",
        "## 各方法相对 Drop 的主指标差值",
        "",
        "| 保留边数 | 方法 | 配对均值 [95% CI] | 正值种子 |",
        "|---:|---|---:|---:|",
    ]
    for remaining in config["remaining_counts"]:
        for method in config["methods"]:
            if method in ("Drop", "Fresh"):
                continue
            row = next(item for item in summaries if item["remaining"] == remaining and
                       item["method"] == method and item["reference"] == "Drop" and item["metric"] == PRIMARY)
            lines.append(
                f"| {remaining} | {method} | {row['mean_delta']:+.4f} "
                f"[{row['ci95_low']:+.4f}, {row['ci95_high']:+.4f}] | "
                f"{row['n_positive']}/{row['n_seeds']} |"
            )
    lines += [
        "",
        "## KF 与其他方法的直接比较",
        "",
        "| 保留边数 | 比较 | 配对均值 [95% CI] |",
        "|---:|---|---:|",
    ]
    for remaining in config["remaining_counts"]:
        for reference in ("Hold", "Linear", "DampedHolt", "SG", "StrengthOnly"):
            row = next(item for item in summaries if item["remaining"] == remaining and
                       item["method"] == "KF" and item["reference"] == reference and item["metric"] == PRIMARY)
            lines.append(
                f"| {remaining} | KF−{reference} | {row['mean_delta']:+.4f} "
                f"[{row['ci95_low']:+.4f}, {row['ci95_high']:+.4f}] |"
            )
    lines += [
        "",
        "## 解释边界",
        "",
        "- 主条件若未在测试集确认 Fresh−Drop 退化，不把补偿方法的正差值解释为恢复遮挡损失。",
        "- r=8 为次要对照，须查看实际 S0；不能用保留边数直接替代三因素社会输入经验量。",
        "- StrengthOnly 只恢复总权重，不预测隐藏方位，是输入幅值对照，不是 KF 的同类预测器。",
        "- 区间按种子配对重采样，均为逐项探索性区间；多方法比较未做多重检验校正。",
        "- 周期遮挡与队友的 6000 步永久遮挡不同；0.06–0.09 仍只是队友的候选敏感区，并非本实验已验证的阈值。",
        "- 详细 GO/LO、全时段 C 和输入误差在 condition_means.csv 与 paired_summary.csv；逐种子配对值在 paired_by_seed.csv。",
        "",
        "## 续跑与参数核查",
        "",
        f"`completion.json` 已确认 {len(seeds)} 个测试种子完整。测试种子：{', '.join(map(str, seeds))}。",
        f"逐一读取 {provenance['verified_seed_count']} 个 `traces.mat` 的内嵌参数，均与一次校准的 `forecast_params.json` 一致。",
        f"已从原始 MAT 轨迹复算 {raw_cells_verified} 个条件单元的 CSV 指标，并核对 6000 步遮挡时序。",
        "SHA256（最终文件快照）：",
        "",
    ]
    for name, digest in provenance["sha256"].items():
        lines.append(f"- `{name}`: `{digest}`")
    lines.append("")
    (out / "README.md").write_text("\n".join(lines), encoding="utf-8")
    return audit


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=("gate", "test", "all"))
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parent)
    args = parser.parse_args()
    root = args.root.resolve()
    config = json.loads((root / "config.json").read_text(encoding="utf-8"))
    if args.stage in ("gate", "all"):
        gate = run_gate(root, config)
        print(json.dumps({"stage": "gate", "status": gate["status"],
                          "selected_schedule": gate["selected_schedule"]}, ensure_ascii=False))
    if args.stage in ("test", "all"):
        if args.stage == "all" and gate["status"] != "pass":
            return
        audit = run_test(root, config)
        print(json.dumps({"stage": "test", "status": audit["status"],
                          "degradation_confirmed": audit["degradation_confirmed"]}, ensure_ascii=False))


if __name__ == "__main__":
    main()

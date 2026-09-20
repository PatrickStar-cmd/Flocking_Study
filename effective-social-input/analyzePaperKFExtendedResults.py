#!/usr/bin/env python3
"""Independent paired analysis for the extended KF occlusion simulation."""

from __future__ import annotations

import csv
import json
import math
import random
from collections import Counter, defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parent
RESULT_DIR = ROOT / "results" / "paper-kf-extended-v3" / "full"
SUMMARY = RESULT_DIR / "summary.csv"
BOOTSTRAP_REPLICATES = 5000
BOOTSTRAP_SEED = 20260915


def parse_summary():
    with SUMMARY.open(encoding="utf-8", newline="") as stream:
        rows = list(csv.DictReader(stream))
    for row in rows:
        row["seed"] = int(row["seed"])
        row["hiddenCount"] = int(row["hiddenCount"])
        row["T"] = float(row["T"])
        row["durationSteps"] = int(row["durationSteps"])
        for key, value in row.items():
            if key in {"mode"} or value == "":
                continue
            try:
                row[key] = float(value)
            except ValueError:
                pass
    return rows


def bootstrap_ci(values, seed_offset=0):
    values = [float(value) for value in values if math.isfinite(float(value))]
    if not values:
        return float("nan"), float("nan"), float("nan")
    mean = sum(values) / len(values)
    if len(values) == 1:
        return mean, mean, mean
    rng = random.Random(BOOTSTRAP_SEED + seed_offset)
    draws = []
    n = len(values)
    for _ in range(BOOTSTRAP_REPLICATES):
        draws.append(sum(values[rng.randrange(n)] for _ in range(n)) / n)
    draws.sort()
    return mean, draws[int(0.025 * BOOTSTRAP_REPLICATES)], draws[
        int(0.975 * BOOTSTRAP_REPLICATES) - 1
    ]


def bootstrap_relative_reduction(baseline, treated, seed_offset=0):
    pairs = [
        (float(base), float(value))
        for base, value in zip(baseline, treated)
        if math.isfinite(float(base)) and math.isfinite(float(value))
    ]
    if not pairs:
        return float("nan"), float("nan"), float("nan")
    base_mean = sum(base for base, _ in pairs) / len(pairs)
    treated_mean = sum(value for _, value in pairs) / len(pairs)
    if abs(base_mean) <= 1e-15:
        return float("nan"), float("nan"), float("nan")
    point = (base_mean - treated_mean) / abs(base_mean)
    if len(pairs) == 1:
        return point, point, point
    rng = random.Random(BOOTSTRAP_SEED + seed_offset)
    draws = []
    n = len(pairs)
    for _ in range(BOOTSTRAP_REPLICATES):
        sample = [pairs[rng.randrange(n)] for _ in range(n)]
        sample_base = sum(base for base, _ in sample) / n
        sample_treated = sum(value for _, value in sample) / n
        draws.append((sample_base - sample_treated) / max(1e-15, abs(sample_base)))
    draws.sort()
    return point, draws[int(0.025 * BOOTSTRAP_REPLICATES)], draws[
        int(0.975 * BOOTSTRAP_REPLICATES) - 1
    ]


def validate(rows):
    expected = 30 * 3 * 6 * 5
    if len(rows) != expected:
        raise ValueError(f"Expected {expected} rows, found {len(rows)}")
    counts = Counter((row["hiddenCount"], row["T"], row["mode"]) for row in rows)
    bad = {key: value for key, value in counts.items() if value != 30}
    if bad:
        raise ValueError(f"Unexpected group counts: {bad}")
    seeds = {row["seed"] for row in rows}
    if seeds != set(range(201, 231)):
        raise ValueError(f"Unexpected seed set: {sorted(seeds)}")
    print(f"validated rows={len(rows)} groups={len(counts)} seeds={len(seeds)}")


def write_group_summary(rows):
    metrics = [
        "eventGOLoss",
        "eventGODeviation",
        "postGOLoss",
        "postGODeviation",
        "eventPairDistanceIncrease",
        "eventPairDistanceDeviation",
        "postPairDistanceIncrease",
        "postPairDistanceDeviation",
        "meanPredictionError",
        "p90PredictionError",
        "recoveryTime",
    ]
    groups = defaultdict(list)
    for row in rows:
        groups[(row["hiddenCount"], row["T"], row["mode"])].append(row)
    output = []
    for (hidden, duration, mode), group in sorted(groups.items()):
        for metric in metrics:
            mean, lo, hi = bootstrap_ci([row[metric] for row in group])
            output.append(
                {
                    "hiddenCount": hidden,
                    "T": duration,
                    "mode": mode,
                    "metric": metric,
                    "n": len(group),
                    "mean": mean,
                    "ci95Low": lo,
                    "ci95High": hi,
                }
            )
    write_csv(RESULT_DIR / "bootstrap-group-summary.csv", output)


def paired_effects(rows):
    by_key = {
        (row["hiddenCount"], row["T"], row["seed"], row["mode"]): row for row in rows
    }
    metrics = [
        "eventGODeviation",
        "postPairDistanceDeviation",
        "eventGOLoss",
        "postPairDistanceIncrease",
    ]
    output = []
    for hidden in (7, 11, 13):
        for duration in (3.0, 6.0, 10.2, 18.0, 30.0, 60.0):
            for mode in ("hold", "mean", "kalman", "kalman-gm"):
                for metric in metrics:
                    baseline_values = []
                    treated_values = []
                    for seed in range(201, 231):
                        baseline = float(by_key[(hidden, duration, seed, "none")][metric])
                        treated = float(by_key[(hidden, duration, seed, mode)][metric])
                        baseline_values.append(baseline)
                        treated_values.append(treated)
                    values = [
                        treated - baseline
                        for baseline, treated in zip(baseline_values, treated_values)
                    ]
                    mean, lo, hi = bootstrap_ci(values, seed_offset=len(output))
                    rel_mean, rel_lo, rel_hi = bootstrap_relative_reduction(
                        baseline_values, treated_values, seed_offset=10000 + len(output)
                    )
                    output.append(
                        {
                            "hiddenCount": hidden,
                            "T": duration,
                            "mode": mode,
                            "metric": metric,
                            "nSeeds": len(values),
                            "meanDifference": mean,
                            "differenceCi95Low": lo,
                            "differenceCi95High": hi,
                            "meanRelativeReduction": rel_mean,
                            "relativeCi95Low": rel_lo,
                            "relativeCi95High": rel_hi,
                        }
                    )

                pred_values = []
                hold_values = []
                method_values = []
                for seed in range(201, 231):
                    hold = float(by_key[(hidden, duration, seed, "hold")]["meanPredictionError"])
                    method = by_key[(hidden, duration, seed, mode)]
                    if mode not in {"hold", "kalman", "kalman-gm"}:
                        continue
                    pred = float(method["meanPredictionError"])
                    pred_values.append(pred - hold)
                    hold_values.append(hold)
                    method_values.append(pred)
                if pred_values:
                    mean, lo, hi = bootstrap_ci(pred_values, seed_offset=20000 + len(output))
                    rel_mean, rel_lo, rel_hi = bootstrap_relative_reduction(
                        hold_values, method_values, seed_offset=30000 + len(output)
                    )
                    output.append(
                        {
                            "hiddenCount": hidden,
                            "T": duration,
                            "mode": mode,
                            "metric": "meanPredictionError_vs_hold",
                            "nSeeds": len(pred_values),
                            "meanDifference": mean,
                            "differenceCi95Low": lo,
                            "differenceCi95High": hi,
                            "meanRelativeReduction": rel_mean,
                            "relativeCi95Low": rel_lo,
                            "relativeCi95High": rel_hi,
                        }
                    )
    write_csv(RESULT_DIR / "paired-effects.csv", output)
    return output


def write_csv(path, rows):
    if not rows:
        raise ValueError(f"No rows available for {path}")
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def main():
    rows = parse_summary()
    validate(rows)
    write_group_summary(rows)
    effects = paired_effects(rows)

    selected = []
    for row in effects:
        if row["metric"] not in {
            "eventGODeviation",
            "postPairDistanceDeviation",
            "meanPredictionError_vs_hold",
        }:
            continue
        if row["hiddenCount"] not in {11, 13}:
            continue
        if row["T"] not in {10.2, 30.0, 60.0}:
            continue
        selected.append(row)
    audit = {
        "input": str(SUMMARY),
        "rows": len(rows),
        "seeds": [201, 230],
        "groups": 90,
        "bootstrapReplicates": BOOTSTRAP_REPLICATES,
        "bootstrapSeed": BOOTSTRAP_SEED,
        "selectedEffects": selected,
    }
    with (RESULT_DIR / "analysis-audit.json").open("w", encoding="utf-8") as stream:
        json.dump(audit, stream, ensure_ascii=False, indent=2)
    print(f"wrote grouped and paired analysis to {RESULT_DIR}")


if __name__ == "__main__":
    main()

# N20, 100-seed final statistical analysis

Source: ../effective-input-n20-100.csv. Reproduce with Python 3 standard library:

```text
python effective-social-input/finalStatistics.py
```

## Integrity

12,700 source rows contain 600 identical duplicate k=0 rows. Analysis retains
12,100 unique (N, seed, policy, strengthMode, k) observations: 121 groups with
seeds 1 through 100 each. No conflicting duplicates or missing seeds were found.
The source CSV is unchanged; its SHA256 is recorded in audit.json.
The six policy/mode variants of k=0 share a simulation and are not independent.
At k=19 all variants are full-input controls; do not pool them as new replicates.

## Calibration and interpretation

The existing MATLAB joint-score definition is reproduced, including linear
percentiles and the fifth ordered baseline score (100 baseline seeds).
Marginal limits: GO >= 0.3224604, LO >= 0.5292991, distance <= 0.2476612.
Joint-score cutoff: -0.1926198. Thus the actual joint acceptance limits are
GO >= 0.2550428, LO >= 0.4777576, distance <= 0.2653841.
This is baseline-relative acceptance, not a formal statistical equivalence test.

Baseline means: GO 0.654514, LO 0.786316, normalized distance 0.156806.
The calibrated joint rule accepts 96/100 baseline seeds, Wilson interval
[0.901628, 0.984337]. Direct conjunction of marginal limits accepts 88/100.
The calibration explicitly targets high in-sample baseline acceptance, so this
is not independent evidence of reliability. GO 0.255 is not high polarization.

## Thresholds

| Policy | Strength mode | First empirical k | Isotonic k | Recalibrated bootstrap 2.5-97.5 percentiles |
|---|---|---:|---:|---|
| random | fixed-total | 4 | 7 | 3-19 |
| nearest | fixed-total | 16 | 16 | 13-18 |
| balanced | fixed-total | 12 | 12 | 10-19 |
| random | original | 17 | 17 | 15-19 |
| nearest | original | 17 | 19 | 17-19 |
| balanced | original | 18 | 18 | 17-19 |

All pointwise Wilson-lower-bound thresholds are absent. Even 99/100 has a
two-sided 95% Wilson lower bound below 0.95. This does not prove failure; it
means a >=95% probability has not been established by that rule.

The bootstrap uses 1,000 resamples of entire seeds, preserving pairing across
conditions and recalibrating the baseline each time. No no-crossing replicates
occurred. Endpoints at 19 are influenced by the full-input baseline and its
in-sample calibration; do not interpret them as proof of a guaranteed threshold.
These are exploratory percentile intervals, not simultaneous confidence sets.
Isotonic fitting imposes monotonicity, which is not established physically.

For random/fixed-total, k=4 passes the calibrated rule 95/100 but the direct
conjunction only 76/100; k=7 gives 96/100 and 89/100 respectively.
Nearest/original gives 96/100 at k=17, 92/100 at k=18 and 96/100 at k=19,
explaining why its empirical and isotonic thresholds differ.

## Outputs and limits

group-success.csv contains original means, seed-paired differences from the
baseline, marginal success counts and Wilson intervals, direct-conjunction
counts, calibrated joint rates and isotonic rates. Single-metric passes use
the original marginal limits, not the relaxed effective joint limits.
thresholds.csv and audit.json record threshold estimates and audit metadata.
An independent PowerShell calculation using the equivalent three inequalities
matched all 121 joint counts.

The results support separate roles for amplitude and neighbor selection under
these parameters. Random selection redraws neighbors each step, so small k
does not imply a fixed small social network over time. The 100-seed baseline
is recalibrated and differs from the earlier 30-seed criterion; threshold shifts
cannot be attributed solely to improved Monte Carlo precision.

Next priority: freeze a scientifically justified success definition, evaluate
on independent seeds, and preserve time-resolved statistics for sustained
success and failure times. The source CSV contains time averages only; it
cannot establish collision safety or sustained collective motion. Neither
additional averaging nor the bootstrap can reconstruct that information.

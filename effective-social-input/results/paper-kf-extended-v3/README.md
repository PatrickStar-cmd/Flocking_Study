# Extended paired KF occlusion simulation

## Purpose

This experiment expands the original paper-aligned comparison while keeping the
old `results/paper-kf-simulation` files unchanged. It uses complete checkpoint
forking: control and all occlusion branches start from the same position,
heading, neural activity, preferred directions, selected directed edges, and
observation-noise stream.

## Design

- N=20, M=100, dt=0.3, L=1000, `v_s=0.05`, `sigma_s=0.4`, `beta=1000`.
- New seeds: 201--230. These seeds were not used in the earlier KF pilots.
- Fixed circular directed topology with k=15 selected inputs.
- Hidden directed inputs: 7, 11, or 13 of 15.
- Blackout horizons: T=3, 6, 10.2, 18, 30, and 60 model-time units.
- Methods: fresh-drop, hold-last, KF-mean, constant-velocity KF-uncertainty,
  and damped angular-velocity KF-GM-uncertainty.
- Measurement standard deviation: 0.01 rad. Observation variance R=1e-4.
- Constant-velocity process noise q=1e-4.
- Damped model: q=1e-5 and tau=10, fixed from the earlier offline development
  study rather than tuned on the new seeds.
- Three hundred recovery steps are simulated; the primary post window is the
  first 100 steps after blackout.

Each seed contributes 90 rows: 3 hidden counts x 6 horizons x 5 methods.
The complete table has 2,700 rows, 90 groups, and 30 independent seeds.

## Main results

The primary deviation measures are:

- `eventGODeviation`: mean absolute difference from the unoccluded control GO
  during the blackout.
- `postPairDistanceDeviation`: mean absolute difference from the unoccluded
  control pair distance in the first 100 post-blackout steps.

For the representative case in which 11 of 15 directed inputs are hidden:

| T | KF-mean GO reduction | KF-unc GO reduction | KF-GM GO reduction | KF-unc pair reduction |
|---:|---:|---:|---:|---:|
| 3 | 86.4% | 86.2% | 86.5% | 89.8% |
| 6 | 82.3% | 82.3% | 83.3% | 89.4% |
| 10.2 | 77.9% | 77.0% | 79.5% | 86.5% |
| 18 | 71.5% | 67.7% | 72.7% | 83.4% |
| 30 | 63.1% | 52.6% | 62.6% | 79.5% |
| 60 | 46.9% | 33.7% | 46.4% | 64.0% |

The reductions are paired against fresh-drop. Their nonparametric bootstrap
95% intervals do not cross zero in the displayed range. The same qualitative
pattern holds for 7/15 and 13/15 hidden inputs: compensation remains beneficial,
and the benefit decreases as the blackout becomes longer.

The KF variants have different strengths:

- KF-mean most consistently reduces GO deviation at long blackouts.
- KF-uncertainty most consistently protects post-blackout pair distance,
  especially for long blackouts and severe occlusion.
- KF-GM-uncertainty gives the best hidden-edge bearing prediction and combines
  much of the GO benefit of KF-mean with improved cohesion.

These statements compare the Kalman variants with fresh-drop. The table and
figures also include hold-last explicitly: at T=60, the GO-deviation reduction
of hold-last can be as large as or slightly larger than that of KF-mean or
KF-GM, whereas KF-uncertainty retains the clearest additional benefit in
post-blackout pair-distance recovery.

Relative to hold-last, KF-GM reduces hidden-edge mean bearing error by about
13%--19% for T=6--30 with bootstrap intervals excluding zero, and about 9%--10%
at T=60. At T=3 the mean reduction is 11%--13%, but the interval crosses zero
for two of the three occlusion loads, so the short-horizon result is directional
rather than significant. Constant-velocity KF is nearly neutral at short
horizons and improves relative to hold-last at longer horizons; its large T=60
reduction should be treated cautiously because angular error is already large.

## Files

- `summary.csv`: all 2,700 paired branch summaries.
- `paired-effects.csv`: paired differences and bootstrap reductions against
  fresh-drop, including hidden-edge prediction error versus hold-last.
- `bootstrap-group-summary.csv`: group means and bootstrap 95% intervals.
- `kf-performance-summary.png`: GO and pair-distance deviation versus T.
- `kf-prediction-error.png`: hidden-edge bearing error versus T.
- `kf-relative-improvement.png`: KF improvement relative to fresh-drop.
- `kf-versus-hold.png`: paired KF-minus-hold deviation difference with intervals.
- `curves.mat`: seed-averaged GO-loss and pair-distance curves.
- `metadata.json`: configuration snapshot.
- `analysis-audit.json`: validation and bootstrap metadata.

## Limits

This is a paired exploratory checkpoint experiment with new seeds, not a
threshold-based confirmatory success test. It does not model neighbor existence
probability, data-association errors, range sensing, physical collision radii,
or an explicit collision barrier. The prediction covariance has not been
recalibrated, and the damped model uses parameters selected in an earlier
offline development study. Do not describe these results as universal
stability or safety guarantees.

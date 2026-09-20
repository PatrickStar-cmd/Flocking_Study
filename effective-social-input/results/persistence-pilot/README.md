# Persistence pilot

This pilot is the first implementation of Stage 2 in `RESEARCH_ROADMAP.md`.
It reruns 40 simulations at N=20 with `recordStride=1`, preserving the original
dynamics and parameters. Conditions are full input, random/fixed-total k=7 and
k=15, and balanced/fixed-total k=15, with seeds 1-10.

The strict success thresholds are frozen from the independent split validation
(calibration seeds 1-50, evaluation seeds 51-100):

- GO >= 0.3413210143
- LO >= 0.5297546056
- normalized mean pair distance <= 0.2427761684

For each post-burn-in time step, success is the conjunction of all three
inequalities. The output reports the fraction of successful time steps, first
failure step, longest failure run, number of failure runs and recoveries.

The time-series means reproduce the existing CSV summaries exactly for every
matching seed and condition (maximum absolute difference 0). This confirms that
recording at every step did not change the simulator or its summary calculation.

Pilot means of sustained-success fraction:

| Condition | Mean fraction |
|---|---:|
| all/original k=19 | 0.8104 |
| random/fixed-total k=7 | 0.7348 |
| random/fixed-total k=15 | 0.6346 |
| balanced/fixed-total k=15 | 0.6566 |

These values are descriptive only: ten seeds are insufficient for inference, and
the strict thresholds were estimated from a related dataset. The surprising
ordering relative to burn-in means is evidence that sustained success is a
different endpoint, not evidence that k=7 outperforms k=15. The threshold and
time-window definition must be frozen before a larger persistence sweep.

`persistence-pilot.csv` and `.mat` are new outputs; no previous results were
overwritten.

# Challenge 2: informed turn under a controlled observation cut

This is a three-seed exploratory development experiment, not confirmation.
No historical MATLAB source or results are modified. The vectorized Python
field and CV-KF are validated against the repository MATLAB implementation.

Protocol frozen before the first run:

- Seeds: 2026101401--2026101403; NumPy PCG64, one serial worker.
- N=20, M=100, dt=.3, L=1000, beta=1000, speed scale=.05, sigma=.4,
  h0=.012, k=15. Each of the 17 followers observes all 3 leaders and the
  next 12 followers in a ring. Each leader observes 2 leaders and 13 followers.
- Construct an ordered initial neural bump, headings near zero, front leaders
  and followers on a perturbed grid. This tests response from an ordered
  constructed state, not spontaneous flock formation or sampled steady states.
- Warm prefix 100 steps. Leaders receive a common static cue during warmup,
  then a common .025 rad/time rotation through 45 degrees and a held endpoint.
  Initial candidate cue amplitude .08; turning begins 20 steps before masking.
  Turn command never enters follower control or the predictor.
- Mask all 51 follower-observes-leader edges, leaving the follower graph
  connected. Total 300 selected edges, 249 fresh during masking. This is a
  synthetic directed observation interruption, not obstacle geometry.
- T={6,10.2,18}; same command and total 500 steps across all conditions.
- Fresh, drop, drop-normalized, hold, KF-mean, KF-uncertainty. Straight-cue
  fresh/drop controls isolate the incremental effect of turning.
- Shared initialization, alpha, neural state, sensor noise, and complete
  pre-event history. Strict common-prefix and full T=0 equality checks.
- KF q=1e-4, R=1e-4, P0=diag(R,1), age limit120, variance limit1. Noise std=.01.
  All hidden edge-time points receive shadow KF diagnostics, including gated
  edges. Scored actual-input errors also report the active-edge count.
- Primary task diagnostics: leader/follower command error, follower turn,
  and paired increase in follower error during masking. Engineering gate:
  final-30-step fresh leader error <=10deg, follower error <=15deg, mean
  follower turn >=30deg, and drop event-error increase >=3deg. These declared
  development cutoffs are not physical or universal success thresholds.
- If the initial task gate fails, retain it. Adjust cue/geometry/lead time
  using fresh/drop task diagnostics only, not KF-versus-hold ranking. Give
  each changed configuration a new directory; retain all three seeds.
- GO and LO include all20. LO=self+5 current nearest by periodic distance.
  Fixed5 distance freezes the five nearest at the pre-event physical state;
  also report current nearest, all-pair, and leader--follower distances.
- S=sum w exp(cos(theta-heading)-1), A=sum w, C=S/A, over follower inputs.
  Fresh/drop/hold/mean use their injected centers. Uncertainty uses the
  posterior expectation of the empirical kernel, with a separately verified
  Fourier calculation. S_oracle uses current same-branch geometry only for
  evaluation. The task cue is excluded from social S/A/C. S is symmetric
  left/right and cannot replace the full direction field or task error.
- Raw trajectories, headings, velocities, final field/filter state, masks,
  noise and per-step CSV are saved. No seed selection, smoothing, or threshold
  fitting to KF outcomes. Statistical unit is the seed, not time or edges.

Run `reference-request`, the MATLAB `validate_reference(out)`, and
`reference-check` before `run`. Every stage output path must be new.
Analyze with `analyze_results.py --stage PATH`. Output metadata and audit
include hashes and environment versions. Reversion: stop using this directory;
if later committed, revert that isolated commit while preserving results.

## Development decisions (retained, not confirmation)

The initial Gaussian cue at amplitudes .08, .32, and 1.28 failed to make the
unoccluded task execute reliably. All three gate batches are retained. A
leader-only signed cosine cue, amplitude .32, then gave leader command errors
near3 degrees. This is an explicit exogenous steering intervention, distinct
from the positive Gaussian social inputs. Followers retain the original field.

The 500-step horizon was insufficient for follower catch-up; the final
comparison uses 800 steps with the same horizon for every method. At20 and40
pre-turn steps the masking effect was small. At100 pre-turn steps (30 time
units), followers still have 46--48 degree errors at masking onset and mean
leader-bearing change during the T=18 window is about18 degrees. The three
drop event-error increases are 2.70,3.41,2.42 degrees. Only1/3 passes the
original3-degree engineering gate; that gate is not relaxed. This setting was
frozen for KF/hold comparison without inspecting their ranking. The command
finishes about1.42 time units into this mask, while relative-bearing dynamics
and follower adjustment continue. Therefore this final batch primarily tests
the catch-up transient after an observed turn, rather than a long unseen leader
turn throughout the mask.

Final invocation:

```text
python run_experiment.py run --amplitude .32 --cue-shape cosine --lead-steps 100 --steps 800 --out NEW_DIRECTORY
python analyze_results.py --stage NEW_DIRECTORY
python verify_results.py --stage NEW_DIRECTORY
```

The program refuses existing stage directories. Figures include full-horizon
curves and event zooms. Earlier gate plots retain their original plotting
windows; their complete time series remain available in CSV. The final batch
also stores complete pre-mask physical/neural/filter checkpoints in each NPZ.

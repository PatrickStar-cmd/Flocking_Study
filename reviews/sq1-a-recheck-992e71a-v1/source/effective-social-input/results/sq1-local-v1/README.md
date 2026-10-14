# SQ1 local E1 replay v1

Formal model definition: `../../SQ1_LOCAL_INPUT_MODEL.md`.

Deterministic single-agent replay for input intensity and bearing structure. No group positions, random seeds, KF, range, or occlusion are used.

- Grid: NNeurons=100, dt=0.3, steps=120, sigma=0.4.
- Input: `I(alpha)=sum_j w_j exp(-d(alpha,theta_j)^2/(2 sigma^2))`.
- `inputMass` is the discrete angular integral; `z1/z2`, maximum gap and entropy describe bearing structure.
- The label-split case tests identical field and neural response after splitting one weighted bearing into duplicate labels.
- Fresh versus hold-last is a local time-pattern comparison, not a group performance or prediction claim.

Cases: 9. Results are exploratory development evidence and are not a population threshold or stability guarantee.

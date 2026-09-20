# Block metric pilot

This exploratory pilot reruns seeds 51-55 for four N=20 conditions and stores
every post-burn-in time block (100 or 300 simulation steps). The final partial
block is retained, so the 2000-step evaluation window produces 20 blocks at
size 100 and 7 blocks at size 300. Each row is one block, not one independent
replicate; inference must aggregate by seed and condition first.

The block endpoint is deliberately descriptive: it records the mean GO, mean
LO, and mean normalized pair distance within each fixed-time block, then applies
the frozen strict thresholds. It does not yet define a final success endpoint.

The pilot produced 540 rows: 5 seeds x 4 conditions x (20 + 7) blocks. The
`nBlocks` field is currently a per-row placeholder and must not be interpreted
as the number of blocks in the condition. The block length is the authoritative
field. This minor metadata defect does not affect the three block means or
`passBlocks` values.

Across all block rows (not an inferential statistic):

| Condition | 100-step passing blocks | 300-step passing blocks |
|---|---:|---:|
| all/original k=19 | 48/100 | 19/35 |
| random/fixed-total k=7 | 85/100 | 30/35 |
| random/fixed-total k=15 | 72/100 | 25/35 |
| balanced/fixed-total k=15 | 73/100 | 27/35 |

These aggregate block counts cannot be compared as independent trials and show
that a block-average endpoint can reverse the ranking of conditions under the
current threshold. That is a property to investigate, not a result to publish.

Before confirmation, implement an analysis that averages block metrics within
each seed, uses duration weighting for the final partial block, and defines the
endpoint from the scientific task or a frozen external criterion. Run synthetic
sequence tests for complete success, complete failure, one short failure and a
failure only at the end.

# r3 GitHub delivery

This directory contains the complete immutable r3 seed1 software acceptance
package. The effective result directory is
`results/sq1-fourbranch-protocol12-seed1-r3-results-v4/`.

The 96 files in `SHA256SUMS.csv` retain their original bytes. This delivery
note is an additional repository document, outside that original manifest.
Historical attempts and the duplicate ZIP are not included.

## Download the complete data

MAT files are stored using Git LFS because `T0_results.mat` exceeds GitHub's
100 MiB ordinary Git file limit. Install Git LFS before cloning, then run:

```powershell
git lfs install
git clone https://github.com/PatrickStar-cmd/Flocking_Study.git
cd Flocking_Study
git lfs pull
cd effective-social-input/results/sq1-fourbranch-protocol12-seed1-r3
./hash_package.ps1 -Verify
```

For the task branch before it is merged, clone with
`--branch feature/sq1-r3-software-package`. A Git LFS pointer is not a MAT
file: use `git lfs pull` when loading reports missing or invalid MAT files.

## Evidence and review

The original MATLAB R2023a run records 15 PASS, 0 FAIL, 0 NOT_RUN in
`tests.json`. These are software acceptance results for seed1, not claims
of scientific benefit. Formal 20-seed/560-trajectory execution remains
pending C's independent acceptance. Four seed1 T0 trajectories are complete;
the other 76 planned T0 trajectories are explicitly NOT_RUN.

Source file hashes, frozen protocol/A references, actual generation paths,
runtime, relative test logs, checkpoint comparisons and summary data are
included. Paths in historical audit records describe the generating machine;
the README and loading example support use from a new checkout.

Original code and historical repository results are unchanged. This package
is added on a task branch for a PR targeting main. Revert the package commit
if repository review requires withdrawing this delivery; retain local
historical experiment evidence.

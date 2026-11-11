# Inbreeding and disease-associated mortality — curated analysis code (2026)

This is a **reduced working candidate**, not a certified or executable end-to-end pipeline. It was constructed by choosing representative scripts from the Inbreeding and SimulationStudy snapshots. The original ZIP archives and Git history are authoritative recovery copies.

## Structure
- `R/`: candidate shared mortality distributions, NIMBLE versions, model-comparison and RJ helpers
- `scripts/`: selected legacy analysis/data preparation, survival model comparisons and results
- `tests/`: original demonstration/testing code plus a distribution validation starter
- `audit/`: exact original-file-to-disposition mapping

## Important limits
1. Some scripts still refer to their old absolute or relative folders and require dependency/path repairs.
2. Scripts need real data and R/NIMBLE to run; **no successful end-to-end rerun is claimed**.
3. The source NIMBLE Siler implementation was excluded due to a confirmed quantile upper-bound error. A tested replacement is required before fitting Siler in NIMBLE.
4. The historical fixed lifetime infection covariate **is not appropriate as a time-varying infection model**. This curated snapshot is only the starting point for a new workflow accepting full posterior infection trajectories.
5. Classification of omitted files as obsolete-for-active-worktree does **not** imply they are statistically incorrect.

## Before using on a server
Preserve the two original source folders (or ZIPs) and run the proposed installation on a new Git branch only. Never remove untracked data or large posterior outputs. Inspect `.gitignore` and compare manifests before committing.

## Scientific plan
- Map individual tattoo IDs and time points to posterior state trajectories and genetic metrics.
- Establish a new capture/censoring/survival data interface.
- Fit updated Gompertz–Makeham and Siler models with uncertainty propagation.
- Compare MLH effects, infection timing and mortality trajectories against published benchmarks.

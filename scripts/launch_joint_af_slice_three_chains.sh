#!/usr/bin/env bash
# Three independent AF_slice chains (6,7,8) of the SAME scientific model.
# Launches one tmux session per chain; never overwrites prior runs 1..5.
# Git revision must be pulled once BEFORE launching; no in-run git operations.
set -Eeuo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."
if [[ "$(git branch --show-current)" != "reanalysis-infection-2026" ]]; then
  echo "Wrong Git branch: expected reanalysis-infection-2026" >&2
  exit 2
fi
if ! command -v tmux >/dev/null 2>&1; then
  echo "tmux is required for independent sessions" >&2
  exit 2
fi
for chain in 6 7 8; do
  printf -v pad "%02d" "$chain"
  session="joint_sei_af_c${chain}"
  dir="outputs/joint_sei_survival/full_joint_20261117/chain_${pad}"
  if tmux has-session -t "$session" 2>/dev/null || \
     [[ -e "${dir}/chain_console.log" || \
        -e "${dir}/chain_started.rds" || \
        -e "${dir}/posterior_samples.rds" || \
        -e "${dir}/process_exit_status.txt" ]]; then
    echo "STOP: chain ${chain} has existing files or tmux session; no chains started." >&2
    exit 3
  fi
done
echo "Git revision: $(git rev-parse --short HEAD)"
echo "Starting 3 AF_slice chains, each with NITER=16000 NBURN=6000 THIN=5."
echo "All 1603 badgers; each chain saves its own posterior and log."
for chain in 6 7 8; do
  session="joint_sei_af_c${chain}"
  tmux new-session -d -s "$session" \
    "cd ~/Inbreeding && SAMPLER_SCHEME=af_slice NITER=16000 NBURN=6000 THIN=5 bash scripts/run_joint_full_chain.sh ${chain}"
  echo "Launched ${session}"
done
tmux ls

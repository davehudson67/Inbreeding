#!/usr/bin/env bash
# v4 joint model chain (v3 + removals censored).
# Usage: bash scripts/run_joint_v4_chain.sh 1
set -Eeuo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."
chain="${1:-1}"
if ! [[ "$chain" =~ ^[1-9][0-9]?$ ]]; then
  echo "CHAIN_ID must be 1..99" >&2
  exit 2
fi
printf -v label "%02d" "$chain"
export RUN_NAME="${RUN_NAME:-full_joint_v4_removals}"
run_dir="outputs/joint_sei_survival/${RUN_NAME}/chain_${label}"
mkdir -p "$run_dir"
log="${run_dir}/chain_console.log"
status_file="${run_dir}/process_exit_status.txt"
if [[ -e "$log" || -e "$status_file" || \
      -e "${run_dir}/chain_started.rds" || \
      -e "${run_dir}/posterior_samples.rds" ]]; then
  echo "Existing chain files found in ${run_dir}; refusing to overwrite." >&2
  exit 3
fi
export MODE=full CHAIN_ID="$chain" SAMPLER_SCHEME=default
export NITER="${NITER:-16000}"
export NBURN="${NBURN:-6000}"
export THIN="${THIN:-5}"
{
  printf "Joint model v4 chain %s\n" "$chain"
  printf "Git revision: %s\n" "$(git rev-parse HEAD)"
  printf "Started UTC: %s\n" "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  printf "NITER=%s NBURN=%s THIN=%s\n" "$NITER" "$NBURN" "$THIN"
  printf "Mortality results remain provisional pending death-cause audit.\n\n"
} > "$log"
set +e
Rscript scripts/joint_sei_survival_model_v4.R 2>&1 | tee -a "$log"
code=${PIPESTATUS[0]}
set -e
{
  printf "CHAIN_ID=%s\n" "$chain"
  printf "EXIT_STATUS=%s\n" "$code"
  printf "FINISHED_UTC=%s\n" "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
} > "$status_file"
printf "CHAIN %s EXIT STATUS: %s\n" "$chain" "$code" | tee -a "$log"
exit "$code"

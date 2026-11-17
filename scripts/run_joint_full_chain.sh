#!/usr/bin/env bash
# Full 1,603-badger chain of the SINGLE joint SEI + survival model.
# Usage: bash scripts/run_joint_full_chain.sh 1
set -Eeuo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."
chain="${1:-1}"
if ! [[ "$chain" =~ ^[1-9][0-9]?$ ]]; then
  echo "CHAIN_ID must be 1..99" >&2
  exit 2
fi
printf -v label "%02d" "$chain"
run_dir="outputs/joint_sei_survival/full_joint_20261117/chain_${label}"
mkdir -p "$run_dir"
log="${run_dir}/chain_console.log"
status_file="${run_dir}/process_exit_status.txt"
if [[ -e "$log" || -e "$status_file" || \
      -e "${run_dir}/chain_started.rds" || \
      -e "${run_dir}/posterior_samples.rds" ]]; then
  echo "Existing chain files found in ${run_dir}; refusing to overwrite." >&2
  exit 3
fi
export MODE=full CHAIN_ID="$chain"
export NITER="${NITER:-8000}"
export NBURN="${NBURN:-3000}"
export THIN="${THIN:-5}"
{
  printf "Single full joint model chain %s\n" "$chain"
  printf "Git revision: %s\n" "$(git rev-parse HEAD)"
  printf "Started UTC: %s\n" "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  printf "NITER=%s NBURN=%s THIN=%s\n" "$NITER" "$NBURN" "$THIN"
  printf "Mortality results remain provisional pending death-cause audit.\n\n"
} > "$log"
set +e
Rscript scripts/joint_sei_survival_model.R 2>&1 | tee -a "$log"
code=${PIPESTATUS[0]}
set -e
{
  printf "CHAIN_ID=%s\n" "$chain"
  printf "EXIT_STATUS=%s\n" "$code"
  printf "FINISHED_UTC=%s\n" "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
} > "$status_file"
printf "CHAIN %s EXIT STATUS: %s\n" "$chain" "$code" | tee -a "$log"
exit "$code"

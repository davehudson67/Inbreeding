#!/usr/bin/env bash
# ==========================================================================
# Stage 71: pilot-gated S->E acquisition vs E->I progression fitting.
# Does not touch existing GM/survival output. Invoke as:
#   bash scripts/71_launch_sei_transition_runs.sh
# ==========================================================================
set -euo pipefail
cd "$(dirname "$0")/.."
root="$PWD"
mkdir -p outputs/sei_transitions

if ! command -v Rscript >/dev/null 2>&1; then
  echo "Rscript not found on this server" >&2
  exit 1
fi
if ! command -v tmux >/dev/null 2>&1; then
  echo "tmux is required to launch persistent full runs" >&2
  exit 1
fi

echo "=== Checking R syntax (fit and aggregation) ==="
Rscript -e 'invisible(parse(file="scripts/69_fit_sei_transition_models.R")); invisible(parse(file="scripts/70_summarise_sei_transition_models.R")); cat("R SYNTAX PASS\n")'

echo "=== Running 25-draw pilot before longer jobs ==="
# Short pilot uses the same underlying data, Cox models and censoring.
# Only the short pilot output is overwritten if this is rerun.
MODE=pilot SHARD_ID=1 N_SHARDS=1 OVERWRITE=1 \
  Rscript scripts/69_fit_sei_transition_models.R \
  2>&1 | tee outputs/sei_transitions/pilot.log

echo "=== PILOT PASSED. Launching 3 disjoint full-draw shards ==="
for i in 1 2 3; do
  session="sei_transition_full_$i"
  filename="$(printf 'full_shard_%02d_of_03' "$i")"
  dest="outputs/sei_transitions/$filename.rds"
  log="outputs/sei_transitions/$filename.log"
  exitfile="outputs/sei_transitions/$filename.exit"

  if [[ -f "$dest" ]]; then
    echo "Shard $i already has saved results: $dest; leaving unchanged."
    continue
  fi
  if tmux has-session -t "=$session" 2>/dev/null; then
    echo "Shard $i tmux session is already running: $session"
    continue
  fi

  # $? is intentionally evaluated inside the tmux shell after Rscript.
  tmux new-session -d -s "$session" \
    "bash -lc 'cd \"$root\"; MODE=full SHARD_ID=$i N_SHARDS=3 Rscript scripts/69_fit_sei_transition_models.R > \"$log\" 2>&1; status=\$?; echo \$status > \"$exitfile\"; exit \$status'"
  echo "Launched shard $i in tmux session $session"
done

echo ""
echo "=== CURRENT SHARD STATUS ==="
for i in 1 2 3; do
  echo "----- SHARD $i -----"
  filename="$(printf 'full_shard_%02d_of_03' "$i")"
  log="outputs/sei_transitions/$filename.log"
  if [[ -f "$log" ]]; then
    tail -n 5 "$log"
  else
    echo "No log yet"
  fi
done
echo ""
echo "When ALL 3 shards exit with code 0 and save their RDS outputs:"
echo "  Rscript scripts/70_summarise_sei_transition_models.R"
echo "No Supabase access needed. No GM MCMC fitting is launched."

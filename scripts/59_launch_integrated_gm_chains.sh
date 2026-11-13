#!/usr/bin/env bash
# Launch the three full-cohort NIMBLE integrated-GM chains.
# Run explicitly on the user's server after benchmark approval.
set -euo pipefail
cd "$HOME/Inbreeding"
command -v tmux >/dev/null || {
  echo "tmux is required for resilient chain sessions." >&2
  exit 1
}
mkdir -p outputs/integrated_gm
test -f outputs/integrated_gm/full_cohort_benchmark.rds || {
  echo "Benchmark missing: run MODE=benchmark Stage 54 first." >&2
  exit 1
}
MAX_HOURS_PER_CHAIN="$(printenv MAX_HOURS_PER_CHAIN || echo 48)"
export MAX_HOURS_PER_CHAIN
Rscript - <<'RS'
x <- readRDS("outputs/integrated_gm/full_cohort_benchmark.rds")
f <- "data/model_inputs/integrated_gm_inputs_group_proxy_2026.rds"
if (!file.exists(f)) stop("Model input missing")
if (!identical(as.character(x$input_md5),
               as.character(unname(tools::md5sum(f)))))
  stop("Model input changed since benchmark; rerun MODE=benchmark")
hours <- x$projected_hours
cat("Benchmarked badgers:", x$n, "\n")
cat("Approximate projected hours per 12000-sweep chain:", hours, "\n")
maxhours <- as.numeric(Sys.getenv("MAX_HOURS_PER_CHAIN","48"))
if (!is.finite(hours) || !is.finite(maxhours) ||
    hours <= 0 || maxhours <= 0)
  stop("Invalid benchmark estimate or MAX_HOURS_PER_CHAIN")
if (hours > maxhours)
  stop("Estimated runtime exceeds MAX_HOURS_PER_CHAIN=",maxhours,
       "; examine speed before launching three chains")
if (x$n != 1603L)
  stop("Benchmark cohort size changed; review before running")
RS

for i in 1 2 3; do
  session="inb_gm_full_c$i"
  logfile="outputs/integrated_gm/full_chain_0$i.log"
  exitfile="outputs/integrated_gm/full_chain_0$i.exit"
  resultfile="outputs/integrated_gm/full_chain_0$i.rds"
  if [[ -e "$resultfile" ]]; then
    echo "Refusing to overwrite an existing completed chain: $resultfile" >&2
    exit 1
  fi
  if tmux has-session -t "$session" 2>/dev/null; then
    echo "Session $session is already running; not launching duplicate"
    continue
  fi
  tmux new-session -d -s "$session" \
    "bash -lc 'cd \"$PWD\"; MODE=full CHAIN_ID=$i NITER=12000 NBURN=4000 THIN=8 Rscript scripts/54_integrated_gm_survival_recapture.R > \"$logfile\" 2>&1; status=\$?; printf \"%s\n\" \"\$status\" > \"$exitfile\"; exit \$status'"
  echo "Launched $session -> $logfile"
done

echo
echo "Active integrated-GM sessions:"
tmux ls | grep -F inb_gm_full || true
echo "Log example: tail -n 30 -F outputs/integrated_gm/full_chain_01.log"
echo "After all chains: Rscript scripts/58_summarise_integrated_gm_chains.R"

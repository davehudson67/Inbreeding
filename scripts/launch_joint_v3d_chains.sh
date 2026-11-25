#!/usr/bin/env bash
# Joint model v3d (infected senescence): 20-badger compile gate, then 3 independent chains
# (tmux ${SESS}_c1..c3; NITER/NBURN/RUN_NAME/SESS overridable; THIN=5).
# Refuses to touch existing chain folders or tmux sessions.
set -Eeuo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."
repo="$(pwd)"
RUN_NAME="${RUN_NAME:-full_joint_v3d_infsen}"; SESS="${SESS:-joint_v3d}"
NITER="${NITER:-16000}"; NBURN="${NBURN:-6000}"
root="outputs/joint_sei_survival/$RUN_NAME"
# v3d lives as new files on joint-v3b-benton (v3 files unchanged there).
case "$(git branch --show-current)" in joint-v3d-infsen|joint-v3b-benton) ;;
  *) echo "STOP: expected branch joint-v3d-infsen or joint-v3b-benton" >&2; exit 2;; esac
command -v tmux >/dev/null || { echo "STOP: tmux required" >&2; exit 2; }
for c in 1 2 3; do
  if tmux has-session -t "${SESS}_c$c" 2>/dev/null ||
     [[ -e "$root/chain_0$c" ]]; then
    echo "STOP: ${SESS}_c$c session or $root/chain_0$c exists" >&2; exit 3
  fi
done
mkdir -p "$root"
gate_log="$root/compile_gate_$(date -u +%Y%m%dT%H%M%SZ).log"
echo "Git revision: $(git rev-parse --short HEAD)"
echo "Compile gate (20 badgers): $gate_log"
MODE=hierarchical_compile OVERWRITE=1 \
  Rscript scripts/joint_sei_survival_model_v3d.R > "$gate_log" 2>&1 || true
if ! grep -q "COMPILATION GATE: PASS" "$gate_log"; then
  tail -n 25 "$gate_log"; echo "STOP: compile gate did not pass" >&2; exit 4
fi
grep -E "posterior log density|COMPILATION GATE" "$gate_log"
for c in 1 2 3; do
  tmux new-session -d -s "${SESS}_c$c" "cd '$repo' && RUN_NAME=$RUN_NAME \
    NITER=$NITER NBURN=$NBURN THIN=5 bash scripts/run_joint_v3d_chain.sh $c"
  echo "Launched ${SESS}_c$c ($RUN_NAME, NITER=$NITER)"
done
tmux ls

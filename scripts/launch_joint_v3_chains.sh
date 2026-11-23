#!/usr/bin/env bash
# Joint model v3: 20-badger compile gate, then 3 independent chains
# (tmux joint_v3_c1..c3; NITER=16000 NBURN=6000 THIN=5; default samplers).
# Refuses to touch existing chain folders or tmux sessions.
set -Eeuo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."
repo="$(pwd)"
root="outputs/joint_sei_survival/full_joint_v3_detcentred"
[[ "$(git branch --show-current)" == "joint-v3-detcentred" ]] ||
  { echo "STOP: expected branch joint-v3-detcentred" >&2; exit 2; }
command -v tmux >/dev/null || { echo "STOP: tmux required" >&2; exit 2; }
for c in 1 2 3; do
  if tmux has-session -t "joint_v3_c$c" 2>/dev/null ||
     [[ -e "$root/chain_0$c" ]]; then
    echo "STOP: joint_v3_c$c session or $root/chain_0$c exists" >&2; exit 3
  fi
done
mkdir -p "$root"
gate_log="$root/compile_gate_$(date -u +%Y%m%dT%H%M%SZ).log"
echo "Git revision: $(git rev-parse --short HEAD)"
echo "Compile gate (20 badgers): $gate_log"
MODE=hierarchical_compile OVERWRITE=1 \
  Rscript scripts/joint_sei_survival_model_v3.R > "$gate_log" 2>&1 || true
if ! grep -q "COMPILATION GATE: PASS" "$gate_log"; then
  tail -n 25 "$gate_log"; echo "STOP: compile gate did not pass" >&2; exit 4
fi
grep -E "posterior log density|COMPILATION GATE" "$gate_log"
for c in 1 2 3; do
  tmux new-session -d -s "joint_v3_c$c" "cd '$repo' && \
    NITER=16000 NBURN=6000 THIN=5 bash scripts/run_joint_v3_chain.sh $c"
  echo "Launched joint_v3_c$c"
done
tmux ls

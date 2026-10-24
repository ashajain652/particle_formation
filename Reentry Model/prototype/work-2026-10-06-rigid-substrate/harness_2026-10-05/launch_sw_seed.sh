#!/bin/zsh
# launch_sw_seed.sh <run name> <plain|diag> <kn switch> <dt fine> <t max> <seed>: launch_sw.sh with `--seed <seed>`
# appended to the CLI arguments (fact 63: a non-default seed ends the run name in _seed-<n>), to measure the scatter
# between seeds at a fine step.
FF="/private/tmp/claude-501/-Users-ashajain-Documents-University-Documents--MIT-Graduate-Work-Research-Space-Sustainability-Particle-Wake-Evolution/930f32b6-e08a-42f0-a4ce-595da1b5cc42/scratchpad/fluxfix"
NAME="$1"; KIND="$2"; KN="$3"; DTF="$4"; TMAX="$5"; SEED="$6"
OUT="$FF/runs/$NAME"
if [ -e "$OUT" ]; then echo "refusing: $OUT exists"; exit 3; fi
case "$KIND" in plain) H="$FF/harness/dtswitch.py";; diag) H="$FF/harness/dtswitch_diag.py";; *) echo "kind?"; exit 2;; esac
mkdir -p "$OUT"
cd "$FF/meas"
echo "start $(date '+%F %T') $(pmset -g batt | head -1)" > "$FF/runs/$NAME.log"
caffeinate -i /Users/ashajain/miniforge3/envs/drama_env/bin/python "$H" "$FF/meas" "$OUT/record.json" "$KN" "$DTF" 10 -- \
  run --diameter 100 --altitude 77.500133 --velocity 7.5 --flight-path-angle -0.959331 --atmosphere us76 --thermal fem \
  --melt on --heating physics --removal girin --outdir "$OUT" --t-max "$TMAX" --frames-every 20 --seed "$SEED" >> "$FF/runs/$NAME.log" 2>&1
echo "exit $?" >> "$FF/runs/$NAME.log"
echo "end $(date '+%F %T') $(pmset -g batt | head -1)" >> "$FF/runs/$NAME.log"

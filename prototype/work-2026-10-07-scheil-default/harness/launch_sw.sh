#!/bin/zsh
# launch_sw.sh <run name> <plain|diag> <kn switch> <dt fine> <t max>: one instrumented 100 mm physics run of the frozen
# measurement copy fluxfix/meas (the runoff-flux fix) into fluxfix/runs/<run name>, under caffeinate, with the power
# source logged at the start and the end. The harness is the parent session's dtswitch.py / dtswitch_diag.py, copied
# unchanged into fluxfix/harness; the CLI arguments are those of dtswitch/launch.sh.
FF="/private/tmp/claude-501/-Users-ashajain-Documents-University-Documents--MIT-Graduate-Work-Research-Space-Sustainability-Particle-Wake-Evolution/930f32b6-e08a-42f0-a4ce-595da1b5cc42/scratchpad/fluxfix"
NAME="$1"; KIND="$2"; KN="$3"; DTF="$4"; TMAX="$5"
OUT="$FF/runs/$NAME"
if [ -e "$OUT" ]; then echo "refusing: $OUT exists"; exit 3; fi
case "$KIND" in plain) H="$FF/harness/dtswitch.py";; diag) H="$FF/harness/dtswitch_diag.py";; *) echo "kind?"; exit 2;; esac
mkdir -p "$OUT"
cd "$FF/meas"
echo "start $(date '+%F %T') $(pmset -g batt | head -1)" > "$FF/runs/$NAME.log"
caffeinate -i /Users/ashajain/miniforge3/envs/drama_env/bin/python "$H" "$FF/meas" "$OUT/record.json" "$KN" "$DTF" 10 -- \
  run --diameter 100 --altitude 77.500133 --velocity 7.5 --flight-path-angle -0.959331 --atmosphere us76 --thermal fem \
  --melt on --heating physics --removal girin --outdir "$OUT" --t-max "$TMAX" --frames-every 20 >> "$FF/runs/$NAME.log" 2>&1
echo "exit $?" >> "$FF/runs/$NAME.log"
echo "end $(date '+%F %T') $(pmset -g batt | head -1)" >> "$FF/runs/$NAME.log"

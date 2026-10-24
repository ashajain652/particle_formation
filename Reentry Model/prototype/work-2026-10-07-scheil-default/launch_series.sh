#!/bin/zsh
# launch_series.sh <run name> <kn switch> <dt fine> [extra CLI args...]: one instrumented 100 mm physics run to 120 s of
# the frozen measurement copy meas/ (Scheil default, rigid substrate on) into runs/<run name>, through the 2026-10-05
# harness (harness/dtswitch.py, unchanged), under caffeinate, the power source logged at the start and the end. The CLI
# arguments are those of the 2026-10-05 series (fluxfix/launch_sw.sh) with no --material, so the default applies.
W="${0:A:h}"
NAME="$1"; KN="$2"; DTF="$3"; shift 3
OUT="$W/runs/$NAME"
if [ -e "$OUT" ]; then echo "refusing: $OUT exists"; exit 3; fi
mkdir -p "$OUT"
cd "$W/meas"
echo "start $(date '+%F %T') $(pmset -g batt | head -1)" > "$W/runs/$NAME.log"
caffeinate -i /Users/ashajain/miniforge3/envs/drama_env/bin/python "$W/harness/dtswitch.py" "$W/meas" "$OUT/record.json" "$KN" "$DTF" 10 -- \
  run --diameter 100 --altitude 77.500133 --velocity 7.5 --flight-path-angle -0.959331 --atmosphere us76 --thermal fem \
  --melt on --heating physics --removal girin --outdir "$OUT" --t-max 120 --frames-every 20 "$@" >> "$W/runs/$NAME.log" 2>&1
echo "exit $?" >> "$W/runs/$NAME.log"
echo "end $(date '+%F %T') $(pmset -g batt | head -1)" >> "$W/runs/$NAME.log"

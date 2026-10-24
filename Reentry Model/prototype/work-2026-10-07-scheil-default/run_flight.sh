#!/bin/zsh
# run_flight.sh <name> <d100|d050> [extra args...]: one physics flight of the frozen copy meas/ (the 2026-10-06 command
# of work-2026-10-06-rigid-substrate/run_flight.sh, material named explicitly by the caller), under caffeinate, the
# power source logged at the start and the end. Used to check that the freeze-back fix leaves the 2026-10-06 flights alone.
PY=/Users/ashajain/miniforge3/envs/drama_env/bin/python
W="${0:A:h}"
name="$1"; case="$2"; shift 2
if [[ $case == d100 ]]; then geo=(--diameter 100 --altitude 77.500133 --t-max 120); else geo=(--diameter 50 --altitude 115); fi
mkdir -p "$W/recheck"
log="$W/recheck/$name.log"
{ echo "power at start: $(pmset -g batt | head -1)"; date; } > "$log"
cd "$W/meas" && caffeinate -i "$PY" -m reentry_model run $geo --velocity 7.5 --flight-path-angle -0.959331 --atmosphere us76 \
    --thermal fem --melt on --heating physics --outdir "$W/recheck" --name "$name" --quiet "$@" >> "$log" 2>&1
rc=$?
{ echo "exit $rc"; date; echo "power at end: $(pmset -g batt | head -1)"; } >> "$log"

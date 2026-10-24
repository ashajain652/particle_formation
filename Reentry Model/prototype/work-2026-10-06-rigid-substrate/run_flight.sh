#!/bin/zsh
# usage: run_flight.sh <code dir> <name> <d100|d050> [extra args...]
# Runs one physics flight under caffeinate, logging the power source at start and end.
PY=/Users/ashajain/miniforge3/envs/drama_env/bin/python
W="${0:A:h}"
code="$1"; name="$2"; case="$3"; shift 3
if [[ $case == d100 ]]; then geo=(--diameter 100 --altitude 77.500133 --t-max 120); else geo=(--diameter 50 --altitude 115); fi
log="$W/runs/$name.log"
{ echo "power at start: $(pmset -g batt | head -1)"; date; } > "$log"
cd "$W/$code" && caffeinate -i "$PY" -m reentry_model run $geo --velocity 7.5 --flight-path-angle -0.959331 --atmosphere us76 \
    --thermal fem --melt on --heating physics --material AA7075_range --outdir "$W/runs" --name "$name" --quiet "$@" >> "$log" 2>&1
rc=$?
{ echo "exit $rc"; date; echo "power at end: $(pmset -g batt | head -1)"; } >> "$log"

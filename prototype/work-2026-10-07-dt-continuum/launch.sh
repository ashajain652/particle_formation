#!/bin/zsh
# launch.sh <name> <d100|d050> [extra CLI args...]: one physics flight of the frozen copy meas/ with the model's own command
# line (no harness), under caffeinate, the power source logged at the start and the end.
W="${0:A:h}"
name="$1"; case="$2"; shift 2
if [[ $case == d100 ]]; then geo=(--diameter 100 --altitude 77.500133 --t-max 120); else geo=(--diameter 50 --altitude 115); fi
mkdir -p "$W/runs"
log="$W/runs/$name.log"
echo "start $(date '+%F %T') $(pmset -g batt | head -1)" > "$log"
cd "$W/meas" && caffeinate -i /Users/ashajain/miniforge3/envs/drama_env/bin/python -m reentry_model run $geo --velocity 7.5 \
  --flight-path-angle -0.959331 --atmosphere us76 --thermal fem --melt on --heating physics --removal girin \
  --outdir "$W/runs" --name "$name" --quiet "$@" >> "$log" 2>&1
echo "exit $?" >> "$log"
echo "end $(date '+%F %T') $(pmset -g batt | head -1)" >> "$log"

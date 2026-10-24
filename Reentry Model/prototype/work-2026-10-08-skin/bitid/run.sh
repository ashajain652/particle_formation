#!/bin/bash
P="/Users/ashajain/Documents/University Documents /MIT Graduate Work/Research/Space Sustainability/Particle Wake Evolution/prototype"
B="$P/work-2026-10-08-skin/bitid"
BASE="$P/work-2026-10-07-dt-continuum/code"; NEW="$P/work-2026-10-08-skin/code"
PY=/Users/ashajain/miniforge3/envs/drama_env/bin/python
COMMON=(run --velocity 7.5 --altitude 77.500133 --flight-path-angle -0.959331 --atmosphere us76 --thermal fem --melt on --heating physics)
echo "power at start: $(pmset -g batt | head -1)  $(date)"
run() { # label dir diameter extra...
  local label=$1 dir=$2 d=$3; shift 3
  (cd "$dir" && caffeinate -i $PY -m reentry_model "${COMMON[@]}" --diameter $d "$@" --outdir "$B/$label" > "$B/$label.log" 2>&1; echo "exit $?" >> "$B/$label.log")
}
run d100_orig "$BASE" 100 --t-max 55 & run d100_new "$NEW" 100 --t-max 55 & wait
run d50_orig "$BASE" 50 & run d50_new "$NEW" 50 & wait
for k in d100 d50; do echo "== $k"; $PY "$P/work-2026-10-07-skinplan/bitid/compare.py" "$B/${k}_orig" "$B/${k}_new"; done
echo "power at end: $(pmset -g batt | head -1)  $(date)"

#!/bin/bash
G="/Users/ashajain/Documents/University Documents /MIT Graduate Work/Research/Space Sustainability/Particle Wake Evolution/prototype/work-2026-10-08-gridgate"
PY=/Users/ashajain/miniforge3/envs/drama_env/bin/python
SEED=12345; OUT="$G/runs/hs14_seed$SEED"; mkdir -p "$OUT"
cd "$G/meas" && echo "power at start: $(pmset -g batt | head -1)" > "$OUT.log" && \
caffeinate -i "$PY" -m reentry_model run --diameter 100 --altitude 77.500133 --velocity 7.5 \
  --flight-path-angle -0.959331 --atmosphere us76 --thermal fem --melt on --heating physics --h-surface 1.4 \
  --seed $SEED --t-max 120 --outdir "$OUT" >> "$OUT.log" 2>&1; echo "exit $?" >> "$OUT.log"
echo "power at end: $(pmset -g batt | head -1)" >> "$OUT.log"

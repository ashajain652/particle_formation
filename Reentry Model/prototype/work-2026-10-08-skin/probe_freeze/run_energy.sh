#!/bin/bash
D="/Users/ashajain/Documents/University Documents /MIT Graduate Work/Research/Space Sustainability/Particle Wake Evolution/prototype/work-2026-10-08-skin/probe_freeze"
cd "/Users/ashajain/Documents/University Documents /MIT Graduate Work/Research/Space Sustainability/Particle Wake Evolution/prototype/work-2026-10-08-skin/code"
echo "power: $(pmset -g batt | head -1) $(date) head $(git rev-parse --short HEAD)"
caffeinate -i /Users/ashajain/miniforge3/envs/drama_env/bin/python -m reentry_model run --diameter 100 --velocity 7.5 --altitude 77.500133 \
  --flight-path-angle -0.959331 --atmosphere us76 --thermal fem --melt on --heating physics --t-max 75 \
  --surface-model skin --skin-preset dev --dt-continuum off --outdir "$D/run_energy_dt05"
echo "exit $?  $(date)"
/Users/ashajain/miniforge3/envs/drama_env/bin/python "$D/compare_freeze.py" "$D/run_dt05" "$D/run_energy_dt05"

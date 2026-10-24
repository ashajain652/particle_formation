#!/bin/sh
# The 100 mm Scheil frames flight on the frame-export amendment (2026-10-08): the MVP flight's flags at the 0.5 s step.
cd "$(dirname "$0")/meas" || exit 2
pmset -g batt | head -1
exec caffeinate -i /Users/ashajain/miniforge3/envs/drama_env/bin/python -m reentry_model run --diameter 100 --altitude 77.500133 \
  --velocity 7.5 --flight-path-angle -0.959331 --atmosphere us76 --thermal fem --melt on --heating physics \
  --material AA7075_scheil --removal girin --deep-runoff off --molten-cascade on --frames-every 1 --dt-continuum off \
  --outdir /Users/ashajain/MIT/particle_formation/reentry_model_output

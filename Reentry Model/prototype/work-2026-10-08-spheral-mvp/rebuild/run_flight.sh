#!/bin/bash
# run_flight.sh <tree> <diameter> <outdir> [extra args...]: one physics-mode melting flight of the prototype in <tree>
# (US76, winds off), as analysis/melt_verification.py's `physics` mode builds it, without a SESAM reference.
tree=$1; d=$2; out=$3; shift 3
case $d in 100) alt=77.500133;; 50) alt=115;; esac
cd "$tree" && exec /Users/ashajain/miniforge3/envs/drama_env/bin/python -m reentry_model run --diameter $d --altitude $alt \
  --velocity 7.5 --flight-path-angle -0.959331 --atmosphere us76 --thermal fem --melt on --heating physics \
  --material AA7075_range --removal girin --outdir "$out" "$@"

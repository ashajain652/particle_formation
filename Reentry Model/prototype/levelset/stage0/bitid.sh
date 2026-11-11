#!/bin/bash
# Stage 0 bit-identity runs: three flights whose history CSVs must be byte-identical before and after the
# shared-package refactor. Usage: bitid.sh <repo checkout> <output dir>
set -euo pipefail
REPO="$1"; OUT="$2"; PY=/Users/ashajain/miniforge3/envs/drama_env/bin/python
mkdir -p "$OUT"; cd "$REPO"
COMMON=(--diameter 100 --velocity 7.5 --altitude 77.500133 --flight-path-angle -0.959331 --outdir "$OUT" --quiet)
"$PY" -m reentry_model run "${COMMON[@]}" --name step1_msis                                     # Step 1, whole flight
"$PY" -m reentry_model run "${COMMON[@]}" --atmosphere us76 --thermal fem --heating physics \
    --t-max 20 --dt 0.5 --name step2_physics                                                  # Step 2, physics heating
"$PY" -m reentry_model run "${COMMON[@]}" --atmosphere us76 --thermal fem --heating sesam \
    --t-max 20 --dt 0.5 --name step2_sesam                                                    # Step 2, SESAM-equivalent heating

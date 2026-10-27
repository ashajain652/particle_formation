"""Wake column: the aluminium a demising droplet vaporises, followed from hot vapour to alumina nanoparticles.

Stage A slows the hot vapour to ambient speed (a slowing-down Monte Carlo that sets the plume's initial composition
and radius); Stage B is a 12-hour box model of a diffusing horizontal cylinder with Plane et al. (2021) gas chemistry,
sectional particles tallied by source molecule and route, uptake onto background meteoric smoke, coagulation and
settling. Design: docs/superpowers/specs/2026-10-27-wake-column-chemistry-design.md.

  species.py      the species table (data/wake_column/species.csv)
  reactions.py    the reaction table (data/wake_column/reactions.csv) and rate evaluation with validity clipping
  settings.py     the Settings, the Case and the run-name scheme
  ambient.py      the ambient-air columns (WACCM baseline, NRLMSIS + AFGL cross-check) and their interpolation
  analysis/       scripts that build the committed tables (WACCM extraction, MSIS + AFGL columns)
  tests/          the checks, run by pytest from Reentry Model/

All code works in SI units; literature units appear only in the committed tables and are converted on load.
"""
import os

from reentry_model import DATA_DIR, OUTPUT_DIR

DATA = os.path.join(DATA_DIR, "wake_column")
OUTPUT = os.path.join(OUTPUT_DIR, "wake_column_output")

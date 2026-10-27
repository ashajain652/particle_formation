"""Wake column: the aluminium a demising droplet vaporises, followed from hot vapour to alumina nanoparticles.

Stage A slows the hot vapour to ambient speed (a slowing-down Monte Carlo that sets the plume's initial composition
and radius); Stage B is a 12-hour box model of a diffusing horizontal cylinder with Plane et al. (2021) gas chemistry,
sectional particles tallied by source molecule and route, uptake onto background meteoric smoke, coagulation and
settling. Design: docs/superpowers/specs/2026-10-27-wake-column-chemistry-design.md.

  species.py      the species table (data/wake_column/species.csv)
  reactions.py    the reaction table (data/wake_column/reactions.csv) and rate evaluation with validity clipping
  settings.py     the Settings, the Case and the run-name scheme
  ambient.py      the ambient-air columns (WACCM baseline, NRLMSIS + AFGL cross-check) and their interpolation
  transport.py    molecular and eddy diffusion, the eddy switch, the cylinder volume and settling rate
  particles.py    size bins, particle size from composition, capture and formation rates, the smoke background
  trail.py        the droplet adapter: trail length, inventory, mean altitude and the emission samples
  slowing.py      Stage A, the slowing-down Monte Carlo (hot start: initial composition and radius)
  box.py          Stage B, the box model: gas chemistry, particles, smoke uptake, coagulation, settling, provenance
  run.py          one case from droplet to outputs (history.csv, size_distribution.csv, summary.json)
  analysis/       the sweep (sweep.py, 126 cases) and sensitivity set (sensitivity.py, 15 settings x 6 cases), their
                  plots (plot_sweep.py, plot_sensitivity.py; results.py holds what they share), and the scripts that
                  build the committed tables (extract_waccm_columns.py, build_msis_afgl_columns.py)
  tests/          the checks, run by pytest from Reentry Model/

Assumptions and their verification status: docs/model_assumptions.md section 9.

All code works in SI units; literature units appear only in the committed tables and are converted on load.
"""
import os

from reentry_model import DATA_DIR, OUTPUT_DIR

DATA = os.path.join(DATA_DIR, "wake_column")
OUTPUT = os.path.join(OUTPUT_DIR, "wake_column_output")

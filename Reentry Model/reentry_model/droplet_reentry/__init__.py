"""Demise of the molten AA7075 droplets shed by a demising fragment: trajectory, heating and vaporisation from release
to a steady temperature and terminal velocity, and optionally settling to the ground (promoted 2026-10-26 from
prototype/work-2026-10-21-droplet-demise, where REPORT.md, DATA_INDEX.md and VERTICAL_WINDS.md record its
development and results).

  model.py                 the droplet (DropletSimulator) and run_case, the entry point for one droplet
  singh_drag.py            Singh et al. (2021) general sphere drag, the default drag law
  particle_trajectory.py   air properties and the earlier composed drag law     } copied from the 2026-10-14
  droplet_heating_sweep.py Schaaf & Chambre free-molecule heat exchange          } settling prototype
  parent_trajectory.py     the parent's DRAMA trajectory (copied from DRAMA Sphere Experiments/sphere_sweep.py)
  analysis/                the sweeps and their plots, writing to OUTPUT (Reentry Model/output/droplet_reentry_output)
  tests/                   the checks, run by pytest from Reentry Model/
"""
import os

from reentry_model import OUTPUT_DIR

OUTPUT = os.path.join(OUTPUT_DIR, "droplet_reentry_output")

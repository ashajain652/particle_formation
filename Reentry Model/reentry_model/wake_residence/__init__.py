"""Wake residence: how long a follower (a solid sphere, or a molten droplet) stays inside its parent's reentry plume.

The parent leaves a plume behind it; the follower is inside it while its distance from the parent's track is less than
the plume radius emitted at the matching earlier time. The engine is general: it takes tracks and flow quantities, not
altitudes. Design: docs/superpowers/specs/2026-10-29-sphere-wake-residence-design.md.

  laws.py         the plume-radius laws (Taylor, Park, Feldman) and the regime choice from Kn, Re_d and Taylor's
                  pressure parameter
  tests/          the checks, run by pytest from Reentry Model/

All code works in SI units; literature units appear only inside the functions that convert them.
"""
import os

from reentry_model import OUTPUT_DIR

OUTPUT = os.path.join(OUTPUT_DIR, "wake_residence_output")

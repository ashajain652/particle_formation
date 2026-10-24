"""First-principles re-entry model of a solid sphere, Steps 1-2: trajectory and coupled 3D heat transfer
(specs under docs/superpowers/specs/)."""
import os

__version__ = "0.1.0"
# "Reentry Model/" holds this package, its data/, output/ and tests/; the repository holds the other projects
PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REPO_ROOT = os.path.dirname(PROJECT_ROOT)
DATA_DIR = os.path.join(PROJECT_ROOT, "data")
OUTPUT_DIR = os.path.join(PROJECT_ROOT, "output")
# the SESAM reference runs and DRAMA space-weather files belong to the DRAMA sphere experiments
DRAMA_DIR = os.path.join(REPO_ROOT, "DRAMA Sphere Experiments")
DRAMA_DATA_DIR = os.path.join(DRAMA_DIR, "data")

"""First-principles re-entry model of a solid sphere, Steps 1-3: trajectory, coupled 3D heat transfer, melting and
melt spraying (specs under docs/superpowers/specs/)."""
import os

__version__ = "0.1.0"
DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")

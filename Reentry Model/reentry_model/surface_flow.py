"""Moved to reentry_model.shared.surface_flow (stage 0 + S of the level-set design). This old path is kept as an alias:
importing it returns the shared module object itself, so attributes, private names, pickles and identity all match."""
import sys

from .shared import surface_flow as _shared

sys.modules[__name__] = _shared

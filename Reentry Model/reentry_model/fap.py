"""Moved to reentry_model.shared.fap (stage 0 + S of the level-set design). This old path is kept as an alias:
importing it returns the shared module object itself, so attributes, private names, pickles and identity all match."""
import sys

from .shared import fap as _shared

sys.modules[__name__] = _shared

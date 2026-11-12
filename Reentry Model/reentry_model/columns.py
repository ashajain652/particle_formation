"""Moved to reentry_model.skyscrapers.columns (skyscrapers subpackage, 2026-11-11). This old path is kept as an alias: importing it returns the moved module object itself, so attributes, private names, pickles and identity all match."""
import sys

from .skyscrapers import columns as _moved

sys.modules[__name__] = _moved

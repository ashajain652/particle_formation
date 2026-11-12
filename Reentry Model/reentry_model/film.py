"""Moved to reentry_model.skyscrapers.elements.film (skyscrapers subpackage, 2026-11-11), its lubrication formulas to reentry_model.shared.lubrication. This old path is kept as an alias: importing it returns the moved module object itself, so attributes, private names, pickles and identity all match."""
import sys

from .skyscrapers.elements import film as _moved

sys.modules[__name__] = _moved

"""FEniCSx (dolfinx) backend: Task 11 fills in the solver; until then the constructor only reports the missing library."""
from . import MissingBackend


class FenicsxThermalSolver:
    def __init__(self, **options):
        try:
            import dolfinx  # noqa: F401
        except ImportError as exc:
            raise MissingBackend("the fenicsx backend needs dolfinx, which is not importable here: create the separate "
                                 "conda environment fenicsx_env (spec section 12) and run with its interpreter") from exc
        raise NotImplementedError("FenicsxThermalSolver is implemented in Task 11")

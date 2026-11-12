"""Finite-element conduction solvers behind one protocol (spec section 8); `thermal_solver(name)` picks the backend.

Both backends solve rho c_p(T) dT/dt = div(k(T) grad T) with -k grad T . n = eps sigma (T^4 - T_amb^4) - q_conv on the
boundary, backward Euler in time, P1 tetrahedra in space, and report the same quantities: nodal temperatures,
stored enthalpy, radiated power. `dirichlet=(nodes, values)` in `step` is for the analytic tests only."""
from ..shared.thermal import SIGMA_SB, MissingBackend, StepResult, ThermalSolver  # noqa: F401

SOLVER_NAMES = ("skfem", "fenicsx")


def thermal_solver(name, **options):
    if name == "skfem":
        from .skfem_backend import SkfemThermalSolver
        return SkfemThermalSolver(**options)
    if name == "fenicsx":
        from .fenicsx_backend import FenicsxThermalSolver      # raises MissingBackend without dolfinx
        return FenicsxThermalSolver(**options)
    raise ValueError("thermal solver must be one of {}, got {!r}".format(SOLVER_NAMES, name))

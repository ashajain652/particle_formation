"""Finite-element conduction solvers behind one protocol (spec section 8); `thermal_solver(name)` picks the backend.

Both backends solve rho c_p(T) dT/dt = div(k(T) grad T) with -k grad T . n = eps sigma (T^4 - T_amb^4) - q_conv on the
boundary, backward Euler in time, P1 tetrahedra in space, and report the same quantities: nodal temperatures,
stored enthalpy, radiated power. `dirichlet=(nodes, values)` in `step` is for the analytic tests only."""
from dataclasses import dataclass
from typing import Protocol

import numpy as np

SIGMA_SB = 5.670374419e-8         # Stefan-Boltzmann [W/(m2 K4)]
SOLVER_NAMES = ("skfem", "fenicsx")


class MissingBackend(RuntimeError):
    """The selected backend's library is not importable in this interpreter."""


@dataclass
class StepResult:
    T: np.ndarray                 # nodal temperatures after the step [K]
    Q_conv: float                 # W, convective power applied over the step
    Q_rad: float                  # W, radiated power at the end of the step
    iterations: int               # Newton/Picard iterations taken
    Q_extra: float = 0.0          # W, nodal load applied over the step (Step 3 deferred melt energy)
    Q_dropped: float = 0.0        # W, load that fell on pinned (material-free) nodes and was not applied
    Q_interface: float = 0.0      # W, heat the skins passed in through their facets (Step 3 sub-plan 18)


class ThermalSolver(Protocol):
    def setup(self, mesh, material, emissivity) -> None: ...
    def set_temperature(self, T) -> None: ...                       # uniform float or nodal array
    def step(self, dt, q_conv, T_amb, dirichlet=None, nodal_load=None, interface=None) -> StepResult: ...
    def temperature(self) -> np.ndarray: ...
    def energy(self) -> float: ...                                  # stored enthalpy above material.T_REF [J]
    def element_energies(self, T=None) -> np.ndarray: ...           # phi_e rho V_e h(T_e) per element [J]
    def radiated_power(self, T_amb) -> float: ...
    def set_fractions(self, phi) -> None: ...                       # element material fractions (0 = dead), Step 3
    def set_film_mass(self, mass) -> None: ...                      # nodal mass of the melt film [kg], Step 3
    def nodal_capacity(self): ...                                   # J/K per node, material + film (Step 3)
    def film_weight(self): ...                                      # the film's share of each node's mass, 0-1 (Step 3)
    def state(self): ...                                            # what set_state needs to repeat a step (sub-plan 18)
    def set_state(self, state) -> None: ...


def thermal_solver(name, **options):
    if name == "skfem":
        from .skfem_backend import SkfemThermalSolver
        return SkfemThermalSolver(**options)
    if name == "fenicsx":
        from .fenicsx_backend import FenicsxThermalSolver      # raises MissingBackend without dolfinx
        return FenicsxThermalSolver(**options)
    raise ValueError("thermal solver must be one of {}, got {!r}".format(SOLVER_NAMES, name))

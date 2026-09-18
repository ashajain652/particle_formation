"""The body the trajectory carries, as a stepper (spec section 4). ConstantBody keeps mass and one temperature
(Step 1 behaviour, `--thermal none`); ThermalBody wraps a mesh, a material and a conduction solver and advances the
temperature field with each macro step's per-patch convective loads, keeping the heat bookkeeping SESAM reports."""
import math
from dataclasses import dataclass
from typing import Protocol

import numpy as np


class Body(Protocol):
    def mass(self, t) -> float: ...
    def advance(self, t, dt, loads) -> None: ...          # loads: heating.HeatingResult applied over [t - dt, t]
    def surface_temperature(self) -> np.ndarray: ...     # K per patch
    def mean_temperature(self) -> float: ...             # K, energy-equivalent (spec 6.4)
    def energy(self) -> float: ...                       # J stored above the material's reference temperature
    def field(self): ...                                 # nodal temperatures, or None


def sphere_mass(diameter_m, density_kgm3):
    """Mass of a solid sphere [kg]."""
    return density_kgm3 * 4.0 / 3.0 * math.pi * (diameter_m / 2.0) ** 3


@dataclass
class ConstantBody:
    mass_kg: float
    temperature_K: float = 300.0

    def mass(self, t):
        return self.mass_kg

    def advance(self, t, dt, loads):
        return None

    def surface_temperature(self):
        return np.array([self.temperature_K])

    def mean_temperature(self):
        return self.temperature_K

    def energy(self):
        return 0.0

    def field(self):
        return None


class ThermalBody:
    """Mesh + material + ThermalSolver. `v_hat` is the direction of motion in the body frame (fixed attitude: the
    stagnation patch is the one whose normal is along v_hat, spec 13.1)."""

    def __init__(self, mesh, material, solver, mass_kg, T0=300.0, emissivity=None, T_ambient=0.0, v_hat=(1.0, 0.0, 0.0)):
        self.mesh, self.material, self.solver, self.mass_kg = mesh, material, solver, mass_kg
        self.emissivity = material.emissivity if emissivity is None else float(emissivity)
        self.T_ambient = float(T_ambient)
        self.surface = mesh.surface()
        self.theta = self.surface.angles_to(v_hat)
        self.i_stag, self.i_back = self.surface.patch_toward(v_hat), self.surface.patch_toward(-np.asarray(v_hat, dtype=float))
        self.i_centre = mesh.centre_node()
        self.volume = mesh.volume()
        self.solver.setup(mesh, material, self.emissivity)
        self.solver.set_temperature(T0)
        self.energy0 = self.solver.energy()
        self.integrated_heat = 0.0       # J, integral of Q_conv dt (SESAM's integrated_heat_J)
        self.absorbed_heat = 0.0         # J, integral of (Q_conv - Q_rad) dt
        self.radiated_heat = 0.0         # J
        self.iterations = []             # Newton iterations per step
        self.last = None                 # thermal.StepResult of the last step

    def mass(self, t):
        return self.mass_kg

    def advance(self, t, dt, loads):
        res = self.solver.step(dt, loads.q_conv, self.T_ambient)
        self.integrated_heat += res.Q_conv * dt
        self.absorbed_heat += (res.Q_conv - res.Q_rad) * dt
        self.radiated_heat += res.Q_rad * dt
        self.iterations.append(res.iterations)
        self.last = res

    def surface_temperature(self):
        return self.surface.facet_mean(self.solver.temperature())

    def mean_temperature(self):
        return float(self.material.temperature_from_enthalpy(self.energy() / (self.material.rho * self.volume)))

    def energy(self):
        return self.solver.energy()

    def field(self):
        return self.solver.temperature()

    def radiated_power(self):
        return self.solver.radiated_power(self.T_ambient)

    def energy_balance_residual(self):
        """(E - E0 - absorbed heat) / absorbed heat: zero for an exact discrete balance."""
        return (self.energy() - self.energy0 - self.absorbed_heat) / self.absorbed_heat if self.absorbed_heat else 0.0

    def surface_stats(self):
        T, Tf = self.solver.temperature(), self.surface_temperature()
        return {"surface_T_max_K": float(Tf.max()), "surface_T_min_K": float(Tf.min()),
                "surface_T_mean_K": float((Tf * self.surface.areas).sum() / self.surface.area),
                "T_stagnation_K": float(Tf[self.i_stag]), "T_back_K": float(Tf[self.i_back]), "T_centre_K": float(T[self.i_centre])}

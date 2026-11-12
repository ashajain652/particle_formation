"""The body seam: the protocol the trajectory steps (Body), a solid sphere's mass and the constant-temperature body
of Step 1 (`--thermal none`)."""
import math
from dataclasses import dataclass
from typing import Protocol

import numpy as np


class Body(Protocol):
    def mass(self, t) -> float: ...
    def advance(self, t, dt, loads) -> None: ...          # loads: heating.HeatingResult applied over [t - dt, t]
    def reference_area(self): ...                        # m2 drag reference area, or None for the fixed pi D^2/4
    def reference_length(self): ...                      # m length of the body Knudsen number, or None for the fixed D
    def drag_shape_factor(self): ...                     # continuum C_D of the current shape / a sphere's (1 for a sphere)
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

    def reference_area(self):
        return None

    def reference_length(self):
        return None

    def drag_shape_factor(self):
        return 1.0

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

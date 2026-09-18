"""The body the trajectory carries. Step 1: constant mass and temperature. Step 2 replaces
ConstantBody with a thermal model whose on_step advances the temperature field."""
import math
from dataclasses import dataclass
from typing import Protocol


class Body(Protocol):
    def mass(self, t: float) -> float: ...
    def temperature(self, t: float) -> float: ...
    def on_step(self, t, state, freestream, aero) -> None: ...


def sphere_mass(diameter_m, density_kgm3):
    """Mass of a solid sphere [kg]."""
    return density_kgm3 * 4.0 / 3.0 * math.pi * (diameter_m / 2.0) ** 3


@dataclass
class ConstantBody:
    mass_kg: float
    temperature_K: float = 300.0

    def mass(self, t):
        return self.mass_kg

    def temperature(self, t):
        return self.temperature_K

    def on_step(self, t, state, freestream, aero):
        return None

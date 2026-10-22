"""Damping: light mass-proportional damping of the replay (Spheral M2 plan, Task 4 Step 2; D5; spec §9.1).

a = -2 zeta omega_1 v on every particle, added to the hydro acceleration in `evaluateDerivatives` (at the step's
midpoint state). Mass-proportional damping decays every mode, and rigid motion, as exp(-zeta omega_1 t): the rate is
zeta omega_1 whatever the mode's own frequency, so a mode of frequency omega has the damping ratio
zeta omega_1 / omega.

omega_1 = pi c_l / D (D5): c_l the longitudinal sound speed at the current scale factor, `material.thermal_sound_speed`
at the scaled free density (sqrt(K/(s rho_free) + 4 G / (3 s rho_free)), what the strength gives Spheral's step), by
default at 300 K (the cold body; `breathing_omega`), and D the body's extent along the flight axis
(`axial_extent`: the particles' axial span plus one spacing, each particle standing for a cell of width dx).
For the 100 mm sphere at 300 K and s = 1: c_l 6,146.9 m/s, omega_1 1.95e5 1/s; the lowest breathing mode of a
uniform elastic sphere (nu = 0.33) lies at 2.745 c_l / R, 1.75 times higher (Love; measured 1.3 % below that in
3D at 5 mm), so its damping ratio is zeta / 1.75. Measured decay rates give zeta within 2.3 % of the set 0.05 and
0.2 (3D and RZ, tests/test_spheral_frag_runner_packages.py).

Default zeta 0.05, bracket 0, 0.05, 0.2 (D5; `ZETA`, `ZETA_BRACKET`). The explicit treatment needs
2 zeta omega_1 dt << 1; `finalize` raises if a step exceeds `DT_FRACTION_MAX` (0.05): at Spheral's step
(dt c_l = 0.28 dx) 2 zeta omega_1 dt = 2 zeta pi 0.28 dx / D, 2.0e-3 at zeta 0.05 and 2.2 mm on the 100 mm body.
`dt()` votes no limit."""
from __future__ import annotations

import math

import numpy as np

from spheral_frag import material
from spheral_frag.runner.packages import LocalCount, add_acceleration, allmax, allmin, get_vector, state_velocity

ZETA = 0.05
ZETA_BRACKET = (0.0, 0.05, 0.2)
DT_FRACTION_MAX = 0.05          # largest 2 zeta omega_1 dt allowed ("<< 1"; explicit first-order damping error)

_CLASSES = {}


def longitudinal_speed(table, scale, T=material.T_ZERO):
    """c_l (m/s) at T and the scaled free density s rho_free(T): what PyThermalStrength gives Spheral's step."""
    T = np.asarray(T, dtype=float)
    rho = float(scale)*table.rho_free(T)
    cf = material.thermal_eos(table, rho, table.enthalpy(T), float(scale))["cs"]
    return material.thermal_sound_speed(table, T, rho, cf)


def axial_extent(body):
    """D (m): the global span of the particles' current axial positions plus one spacing (collective)."""
    x = get_vector(body.nodes.positions(), body.dim)[:, 0]
    lo = allmin(x.min() if len(x) else np.inf)
    hi = allmax(x.max() if len(x) else -np.inf)
    return hi - lo + body.dx


def breathing_omega(table, scale, extent, T=material.T_ZERO):
    """omega_1 = pi c_l / D (1/s), D5."""
    return math.pi*float(longitudinal_speed(table, scale, T))/float(extent)


def damping_class(S):
    if S in _CLASSES:
        return _CLASSES[S]

    class Damping(S.Physics):
        """a = -2 zeta omega_1 v (module docstring). Damping(body, zeta, omega1). Attribute `max_fraction`: the
        largest 2 zeta omega_1 dt seen on this rank."""

        def __init__(self, body, zeta, omega1):
            S.Physics.__init__(self)
            if zeta < 0.0 or not omega1 > 0.0:
                raise ValueError(f"need zeta >= 0 and omega1 > 0, got {zeta}, {omega1}")
            self.body = body
            self.dim = body.dim
            self.zeta, self.omega1 = float(zeta), float(omega1)
            self.c = 2.0*self.zeta*self.omega1
            self.count = LocalCount(body.nodes)
            self.max_fraction = 0.0

        def evaluateDerivatives(self, t, dt, db, state, derivs):
            if self.c == 0.0:
                return
            self.count.check("Damping")
            add_acceleration(S, derivs, self.dim, -self.c*state_velocity(S, state, self.dim))

        def finalize(self, t, dt, db, state, derivs):
            f = self.c*float(dt)
            self.max_fraction = max(self.max_fraction, f)
            if f > DT_FRACTION_MAX:
                raise RuntimeError(f"Damping: 2 zeta omega_1 dt = {f:.3e} exceeds {DT_FRACTION_MAX} (zeta "
                                   f"{self.zeta}, omega_1 {self.omega1:.4e} 1/s, dt {dt:.4e} s)")

        def dt(self, db, state, derivs, t):
            return (1.0e100, "Damping: no limit")

        def registerState(self, db, state):
            return

        def registerDerivatives(self, db, derivs):
            return

        def label(self):
            return "Damping"

    _CLASSES[S] = Damping
    return Damping


def make_damping(S, body, zeta, omega1):
    """A Damping package for the Spheral module S."""
    return damping_class(S)(body, zeta, omega1)

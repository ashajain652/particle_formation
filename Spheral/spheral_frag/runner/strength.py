"""PyThermalStrength: the strength of spec §8.3 without damage as a Python subclass of Spheral's StrengthModel (Spheral
M2 plan, Task 3; the mechanism of M0's PySteinbergGuinan, spheral_frag/m0/subclass_check.py).

  - shearModulus: `material.thermal_shear_modulus` at T, the table's G(T), zero at and above 50 % liquid.
  - yieldStrength: `material.flow_stress(table, T, max(plasticStrainRate, rate_min), law)` (D6; decision 4:
    rate_min = 1e-6 1/s, bracket 1e-7 and 1e-5). Laws: "flight" (Chen 2016 Eq. 20; provisional below 623 K,
    `material.FLOW_STRESS_PROVISIONAL`), "elastic", "constant:<Pa>". The elastic law's infinite yield is written as
    DBL_MAX (Spheral's PlasticStrainPolicy computes min(1, Y/sqrt(3 J2)), which is 1 either way; a finite value keeps
    inf out of the node fields and checkpoints). plasticStrainRate is the field Spheral hands the method (the
    previous update's (ps - ps0)/dt from PlasticStrainPolicy).
  - providesSoundSpeed is True and soundSpeed is the longitudinal speed `material.thermal_sound_speed`,
    sqrt(c_f^2 + 4 G(T)/(3 rho)) with c_f the EOS's sqrt(K/(s rho_free)) that SoundSpeedPolicy passes in, as Spheral's
    SteinbergGuinan provides it: GenericHydro's time step reads only the soundSpeed field, so without it the step
    would follow the bulk speed and the effective Courant number would be sqrt((K + 4G/3)/K) = 1.23 times the
    intended one (controller review of Task 3). Under scaling c goes as 1/sqrt(s). providesBulkModulus is False
    (the EOS's rho dP/drho is used).
  - T comes from the EOS's table path, `eos.temperature(u)`: the same exact inverse (and cache) as the EOS's.
  - Damage is not applied (M2 has none; M3's tearing law will).

Strength methods are called on internal values only (as the built-ins loop over numInternalElements).
`make_strength(S, eos, law, rate_min)` builds the class for the Spheral module S."""
from __future__ import annotations

import sys

import numpy as np

from spheral_frag import material
from spheral_frag.runner.eos import get_internal, put_internal

DBL_MAX = sys.float_info.max

_CLASSES = {}


def strength_class(S):
    """The PyThermalStrength class for the Spheral module S (built once per module)."""
    if S in _CLASSES:
        return _CLASSES[S]

    class PyThermalStrength(S.StrengthModel):
        """spec §8.3's flow stress and shear modulus (module docstring).

        PyThermalStrength(eos, law="flight", rate_min=1e-6), eos a PyThermalEOS (its table and T(u) are used)."""

        def __init__(self, eos, law="flight", rate_min=material.RATE_MIN):
            S.StrengthModel.__init__(self)
            material.parse_flow_law(law)               # ValueError on an unknown law
            if not float(rate_min) > 0.0:
                raise ValueError(f"rate_min must be positive, got {rate_min}")
            self.eos = eos
            self.table = eos.table
            self.law = law
            self.rate_min = float(rate_min)

        @property
        def provisional(self):
            return material.flow_stress_provisional(self.law)

        def providesSoundSpeed(self):
            return True

        def providesBulkModulus(self):
            return False

        def shear(self, u):
            return material.thermal_shear_modulus(self.table, self.eos.temperature(u))

        def yield_strength(self, u, plastic_strain_rate):
            Y = material.flow_stress(self.table, self.eos.temperature(u), plastic_strain_rate, self.law,
                                     self.rate_min)
            return np.minimum(Y, DBL_MAX)

        def shearModulus(self, shearModulus, density, specificThermalEnergy, pressure, damage):
            put_internal(shearModulus, self.shear(get_internal(specificThermalEnergy)))

        def yieldStrength(self, yieldStrength, density, specificThermalEnergy, pressure, plasticStrain,
                          plasticStrainRate, damage):
            put_internal(yieldStrength, self.yield_strength(get_internal(specificThermalEnergy),
                                                            get_internal(plasticStrainRate)))

        def soundSpeed(self, soundSpeed, density, specificThermalEnergy, pressure, fluidSoundSpeed, damage):
            # SoundSpeedPolicy passes the same field as soundSpeed and fluidSoundSpeed: read it before writing
            cf = np.array(get_internal(fluidSoundSpeed))
            put_internal(soundSpeed, material.thermal_sound_speed(
                self.table, self.eos.temperature(get_internal(specificThermalEnergy)), get_internal(density), cf))

    _CLASSES[S] = PyThermalStrength
    return PyThermalStrength


def make_strength(S, eos, law="flight", rate_min=material.RATE_MIN):
    """A PyThermalStrength for the Spheral module S on the PyThermalEOS `eos`."""
    return strength_class(S)(eos, law, rate_min)

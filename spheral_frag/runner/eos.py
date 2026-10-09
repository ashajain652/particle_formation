"""PyThermalEOS: the equation of state of spec §8.2 as a Python subclass of Spheral's SolidEquationOfState (Spheral M2
plan, Task 3; the mechanism of M0's PyMurnaghan, spheral_frag/m0/subclass_check.py).

The physics is `material.thermal_eos` (core, numpy only); this class only moves Spheral's fields to numpy and back.
With the material table of the prepared run and the mass scale factor s:

  - u, Spheral's specific thermal energy, is the table's enthalpy in J/kg from 300 K, latent heat included (spec
    §8.1; D2: imposed and unscaled). The table's zero is 300 K itself (`MaterialTable.enthalpy(300 K)` = 0 exactly;
    the finite element's 293 K zero is removed inside the table by the exact offset h_FE(300 K)), so the temperature
    is `table.temperature(u)`, the table's exact inverse, with no offset added here. Spheral's own temperature
    relations are never used (Review focus 2).
  - P = K(T) (rho/(s rho_free(T)) - 1); dP/drho = K/(s rho_free); the sound speed sqrt(dP/drho); the bulk modulus
    rho dP/drho = K(T) rho/(s rho_free), the same dP/drho (spec §2 traps, Review focus 3: the damage model reads the
    bulk modulus; it equals K(T) at the free density and differs from it by the strain eta - 1 elsewhere); dP/du the
    exact derivative of the law at fixed rho (P depends on u only through T: (dP/dT)/(dh/dT), with the table's
    piecewise-linear K, rho_free and piecewise-quadratic h, one-sided at their nodes; only porosity reads it in
    Spheral at 116c71f); gamma 1 + (s/rho) dP/du (a Grueneisen form with the physical density; Spheral's SolidSPH
    does not read it); c_p(T) the table's c_p without latent part; the entropy zero (undefined for this material; no
    package used here reads it).
  - Pressure limits, the external pressure and the pressure floor as Spheral's built-ins apply them
    (EquationOfStateInline.hh: P - P_ext, then below minimumPressure the floor or zero, above maximumPressure the
    maximum); unbounded by default.
  - Mass scaling (spec §9.1, Review focus 1): s is a constructor argument and the settable attribute `scale`. The
    reference density is s rho_free(300 K) and follows `scale`, so rho/referenceDensity (etamin, etamax, and the
    reference density SolidSPHRZ relaxes damaged material towards) is scale-invariant. The node list's own density
    limits rhoMin and rhoMax are the node list's (nodes.make_body sets them from the densities it is given, i.e. the
    scaled ones); whoever changes s on a live body (the switch, Task 11) must change those too.

Every method works on the whole field: the EOS methods on internal and ghost values (Spheral calls them on whole
fields), the T(u) inverse cached on the last energies seen (bitwise equality), the full state on the last (rho, u, s).

Spheral is imported by the caller: `make_eos(S, table, scale)` builds the class for the given Spheral module
(Spheral3d or SpheralRZ, one per process).

Measured 2026-10-09 (tests/test_spheral_frag_runner_material.py): every method equals the core bitwise on 4,000
states in 3D and RZ; with PyThermalStrength a step of the 2.2 mm 3D sphere costs about 6 % more than Murnaghan +
ConstantStrength on 18 processes (10 % on 8)."""
from __future__ import annotations

import struct
import sys

import numpy as np

from spheral_frag import material

DBL_MAX = sys.float_info.max
F64 = np.dtype("<f8")
ETA_RANGE = (0.2, 4.0)       # etamin, etamax: as nodes.RHO_RANGE_FACTORS (no method of this class reads them)

_CLASSES = {}


# ---------------------------------------------------------------------------------------------------------------
# Field <-> numpy (as M0's subclass_check: allValues/assign for whole fields, serialize/deserialize for internal values)
# ---------------------------------------------------------------------------------------------------------------

def get_all(field):
    """Internal and ghost values of a scalar Field."""
    return np.array(field.allValues(), dtype=F64)


def put_all(field, values):
    """Assign every value (internal and ghost) of a scalar Field."""
    field.assign(np.asarray(values, dtype=F64).tolist())


def get_internal(field):
    """The internal values of a scalar Field as a read-only numpy array (Field.serialize: name, count, values)."""
    b = field.serialize()
    n = field.numInternalElements
    off = 16 + len(field.name.encode())
    return np.frombuffer(b, dtype=F64, offset=off)[:n]


def put_internal(field, values):
    """Write the internal values of a scalar Field (ghosts untouched) through Field.deserialize."""
    name = field.name.encode()
    v = np.ascontiguousarray(values, dtype=F64)
    field.deserialize(struct.pack("<Q", len(name)) + name + struct.pack("<Q", v.size) + v.tobytes())


class _Cache:
    """The last result of f(*arrays, *scalars), reused while the arguments are bitwise the same."""

    def __init__(self, f):
        self.f, self.key, self.val = f, None, None

    def __call__(self, *args):
        key = self.key
        if key is not None and len(key) == len(args) and all(
                (a == b) if not isinstance(a, np.ndarray) else (a.shape == b.shape and np.array_equal(a, b))
                for a, b in zip(args, key)):
            return self.val
        self.val = self.f(*args)
        self.key = tuple(a.copy() if isinstance(a, np.ndarray) else a for a in args)
        return self.val

    def clear(self):
        self.key = self.val = None


def eos_class(S):
    """The PyThermalEOS class for the Spheral module S (built once per module)."""
    if S in _CLASSES:
        return _CLASSES[S]

    class PyThermalEOS(S.SolidEquationOfState):
        """spec §8.2's equation of state on the M1 material table (module docstring).

        PyThermalEOS(table, scale=1.0, constants=None, etamin=0.2, etamax=4.0, minimumPressure=-DBL_MAX,
        maximumPressure=DBL_MAX, externalPressure=0.0)"""

        def __init__(self, table, scale=1.0, constants=None, etamin=ETA_RANGE[0], etamax=ETA_RANGE[1],
                     minimumPressure=-DBL_MAX, maximumPressure=DBL_MAX, externalPressure=0.0):
            if constants is None:
                constants = S.MKS()
            scale = float(scale)
            if not scale > 0.0:
                raise ValueError(f"scale must be positive, got {scale}")
            self.rho_ref = float(table.rho_free(material.T_ZERO))      # rho_free(300 K), unscaled
            S.SolidEquationOfState.__init__(self, scale*self.rho_ref, etamin, etamax, constants, minimumPressure,
                                            maximumPressure, 0.0, S.MaterialPressureMinType.PressureFloor,
                                            externalPressure)
            self._constants = constants          # the C++ base keeps a reference to it
            self._floor = S.MaterialPressureMinType.PressureFloor
            self.table = table
            self._scale = scale
            self._T = _Cache(table.temperature)
            self._state = _Cache(lambda rho, u, s: material.thermal_eos(table, rho, u, s))

        # -- the scale factor (spec §9.1; the switch of Task 11 sets it)
        @property
        def scale(self):
            return self._scale

        @scale.setter
        def scale(self, s):
            s = float(s)
            if not s > 0.0:
                raise ValueError(f"scale must be positive, got {s}")
            self._scale = s
            self.referenceDensity = s*self.rho_ref

        # -- numpy
        def temperature(self, u):
            """K from u (J/kg from 300 K): the table's exact inverse."""
            return self._T(np.asarray(u, dtype=F64))

        def state(self, rho, u):
            """material.thermal_eos at the current scale."""
            return self._state(np.asarray(rho, dtype=F64), np.asarray(u, dtype=F64), self._scale)

        def limits(self, P):
            """EquationOfState::applyPressureLimits after the external pressure is taken off."""
            ext = self.externalPressure
            if ext != 0.0:
                P = P - ext
            pmin, pmax = self.minimumPressure, self.maximumPressure
            low = pmin if self.minimumPressureType == self._floor else 0.0
            return np.where(P < pmin, low, np.where(P > pmax, pmax, P))

        # -- EOSAbstractMethods
        def setPressure(self, Pressure, massDensity, specificThermalEnergy):
            st = self.state(get_all(massDensity), get_all(specificThermalEnergy))
            put_all(Pressure, self.limits(st["P"]))

        def setPressureAndDerivs(self, Pressure, dPdu, dPdrho, massDensity, specificThermalEnergy):
            st = self.state(get_all(massDensity), get_all(specificThermalEnergy))
            put_all(Pressure, self.limits(st["P"]))
            put_all(dPdu, st["dPdu"])
            put_all(dPdrho, st["dPdrho"])

        def setTemperature(self, temperature, massDensity, specificThermalEnergy):
            put_all(temperature, self.temperature(get_all(specificThermalEnergy)))

        def setSpecificThermalEnergy(self, specificThermalEnergy, massDensity, temperature):
            put_all(specificThermalEnergy, self.table.enthalpy(get_all(temperature)))

        def setSpecificHeat(self, specificHeat, massDensity, temperature):
            put_all(specificHeat, self.table.cp(get_all(temperature)))

        def setSoundSpeed(self, soundSpeed, massDensity, specificThermalEnergy):
            put_all(soundSpeed, self.state(get_all(massDensity), get_all(specificThermalEnergy))["cs"])

        def setGammaField(self, gamma, massDensity, specificThermalEnergy):
            put_all(gamma, self.state(get_all(massDensity), get_all(specificThermalEnergy))["gamma"])

        def setBulkModulus(self, bulkModulus, massDensity, specificThermalEnergy):
            put_all(bulkModulus, self.state(get_all(massDensity), get_all(specificThermalEnergy))["K"])

        def setEntropy(self, entropy, massDensity, specificThermalEnergy):
            put_all(entropy, np.zeros(entropy.numElements))

        def valid(self):
            return bool(self._scale > 0.0 and self.referenceDensity > 0.0 and self.etamin <= self.etamax)

    _CLASSES[S] = PyThermalEOS
    return PyThermalEOS


def make_eos(S, table, scale=1.0, **kw):
    """A PyThermalEOS for the Spheral module S on `table` (a material.MaterialTable) at scale factor s."""
    return eos_class(S)(table, scale, **kw)

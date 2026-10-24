"""Python subclasses of Spheral's equation-of-state, strength and damage classes against the built-ins (Spheral M0
plan, Task 7; spec section 12, runner check 2).

What is mimicked, and how Spheral allows it (all at 116c71f, read from /home/spheral/workspace in the image):

  - Equation of state: PyMurnaghan, a Python subclass of SolidEquationOfState copying MurnaghanEquationOfState
    (SolidMaterial/MurnaghanEquationOfState.cc), bench.py's EOS and route 2's (spec section 8.2), so that it is
    compared on the same baseline as the strength and damage classes; and PyLinearPolynomial, the plan's choice,
    copying LinearPolynomialEquationOfState, compared against its own built-in run. Every abstract method of
    EOSAbstractMethods (PYB11/Material/EOSAbstractMethods.py) is overridden: setPressure, setPressureAndDerivs,
    setTemperature, setSpecificThermalEnergy, setSpecificHeat, setSoundSpeed, setGammaField, setBulkModulus,
    setEntropy, valid. PYB11Generator writes a pybind11 trampoline for every bound class with virtuals
    (build_docker-gcc/build/src/PYB11/SolidMaterial/current_SpheralSolidMaterial/SpheralSolidMaterial.cc:13, class
    PYB11TrampolineSolidEquationOfState; trampolines exist for the built-in EOS and strength classes too).
  - Strength: PySteinbergGuinan, a Python subclass of StrengthModel copying SteinbergGuinanStrength
    (SolidMaterial/SteinbergGuinanStrength.cc) with bench.py's 6061-T6 parameters: shearModulus, yieldStrength,
    soundSpeed, meltSpecificEnergy, coldSpecificEnergy, providesSoundSpeed, providesBulkModulus overridden; the
    temperature goes through the EOS's virtual setTemperature as in the C++. The cold- and melt-energy fits are no
    obstacle: bench.py's are 0 and a constant, and NinthOrderPolynomialFit is a plain sum (PolynomialFitInline.hh:83).
  - Damage: PyDelegatingDamage, a Python subclass of DamageModel through its trampoline
    (current_SpheralDamage/SpheralDamage.cc:87, PYB11TrampolineDamageModel, which overrides every Physics and
    DamageModel virtual: evaluateDerivatives, computeScalarDDDt, dt, dtImplicit, registerState, registerDerivatives,
    initializeProblemStartup(Dependencies), preStepInitialize, initialize, finalizeDerivatives, postStateUpdate,
    finalize, applyGhostBoundaries, enforceBoundaries, requireConnectivity/Volumes/ReproducingKernels, maxResidual,
    extraEnergy/Momentum, registerAdditionalVisualizationState, label, dumpState, restoreState). It holds a C++
    ProbabilisticDamageModel that is not given to the integrator, forwards every virtual to it, and computes only
    the damage rate itself: computeScalarDDDt (Damage/DamageModel.cc:83-122), dD^(1/3)/dt = 0.4 c_l / (kernel extent
    h), as ProbabilisticDamageModel::evaluateDerivatives calls it (ProbabilisticDamageModel.cc:229-244). That is the
    method the tearing law of spec section 8.4 replaces: evaluateDerivatives is the per-step hook in which a damage
    model turns the state into a damage rate, and section 8.4's "where" (the brittle window), "when" (the critical
    strain) and "how" (growth at 0.4 c_l along the greatest tension) all become that rate. The Weibull flaw
    activation of ProbabilisticDamagePolicy::update (enrolled by registerState, ProbabilisticDamagePolicy.cc:117-227)
    and the strain policy will also have to give way to the tearing law's strain and threshold; registerState is
    overridable (it is forwarded here) and update policies have trampolines too (PYB11TrampolineUpdatePolicyBase in
    current_SpheralDataBase), but a Python policy is not tried here. The trampoline allows the delegation; one trap:
    Physics::appendBoundary is not virtual (Physics/Physics.hh:98), so the boundaries SpheralController gives the
    wrapper (in parallel the distributed boundary, SimulationControl/SpheralController.py:853-869) never reach the
    delegate, whose applyGhostBoundaries then fills no ghost damage. Without copying them over (_sync_boundaries) the
    2-process 3D run at 6 mm differed after 7 steps by 1e-11 of the largest velocity; with it, bitwise nothing.
    bench.py's built-in damage model is itself a Python subclass (ShadowProbabilisticDamageModel.py), so the
    baseline already pays a trampoline lookup on every virtual call.

Exact copies. Spheral is compiled by g++ 13 -O3 -std=c++20 (CMakeCache.txt), which on aarch64 contracts a*b + c into
fused multiply-adds: 290 FMADD/FMSUB/FNMSUB instructions in SteinbergGuinanStrengthInst.cc.o, 51 in
LinearPolynomialEquationOfStateInst.cc.o. A straightforward numpy port therefore differs from the built-in by an ULP
here and there. The classes here repeat the fusions read from the object code (fmadd/fmsub/fnmsub, formed exactly in
IEEE quad via numpy longdouble, which on aarch64 Linux has 113 bits; checked identical to glibc fma() on 200,000
triples including heavy cancellation), including those inside Spheral's closed-form eigenvalues of the damage
tensor (GeomSymmetricTensorInline.hh:2043-2090) and the 3D determinant of H in computeScalarDDDt. numpy's exp, pow,
sin, cos, atan2 and sqrt on this image are glibc's (checked identical on 200,000 arguments each). With --fma plain
the classes use unfused arithmetic instead, to measure what the emulation costs and how far ULP differences grow.

Measured (2026-10-07, run_subclass.sh, 18 processes, damage on, 10 warm-up + 100 timed steps, host 1-minute load
1.75-1.98 at every start):
  - Method by method on 2,000-4,000 synthetic points reaching every branch (both density clamps, melted and cold
    energies, damage eigenvalues outside [0, 1], excluded nodes): all 29 compared outputs bitwise identical to the
    built-ins in 3D and RZ.
  - End to end, built-in against Python for each class alone and for all three together, RZ at 1.1 mm (3,243
    particles, 166-181 per process) and 3D at 3.0 mm (19,381, 1,072-1,077 per process), on the same particles (a
    lattice id), the same decomposition and the same steps: position, velocity, deviatoric stress, pressure, damage,
    density, energy and plastic strain bitwise identical, and so the same dt sequence. With --fma plain the strength
    class (and all three) end 1.7e-11 (3D) and 3.0e-11 (RZ) of the largest velocity and deviatoric stress away.
  - Added cost per step (median against median; a repeat of the built-in run came out 2.3 % and 2.0 % slower, the
    noise floor): 3D, built-in 47.5 ms: EOS +5.0 %, strength +13.2 %, damage +3.5 %, all three +19.7 % (+12.1 % with
    plain arithmetic); RZ, built-in 3.45 ms: EOS +18.7 %, strength +41.1 %, damage +33.9 %, all three +79.6 %
    (+71.4 % plain). The Python cost is mostly per call (dozens of virtual calls per step, each converting Fields
    through serialize/assign), so it is a large fraction of the small RZ step at under 200 particles per process.

Run through the launcher from the repository root:

    spheral_frag/container/spheral -n 18 spheral_frag/m0/subclass_check.py run --geometry rz --dx-mm 1.1 \
        --variant strength                       # builtin, builtin_repeat, eos, strength, damage, all,
                                                 # linpoly_builtin, linpoly_python; --fma exact|plain
    spheral_frag/container/spheral spheral_frag/m0/subclass_check.py methods --geometry 3d --json FILE
    spheral_frag/container/spheral spheral_frag/m0/subclass_check.py compare     # -> OUT/summary.json

or the whole study with spheral_frag/m0/run_subclass.sh. A run writes OUT/<run name>/summary.json, steps.csv and
state.npz (the compared fields of every particle, ordered by lattice id); the run name encodes the whole
configuration and an existing summary.json skips the run. Timing as bench.py (the body of the controller's advance
loop between barriers, rank 0's view; the VM's load checked against --max-load unless --force; the host's passed in
with --host-load); the time spent inside the Python overrides is recorded per method (it includes the C++ work the
damage wrapper forwards to). Exit 0 ok, 1 a Spheral failure, 2 bad arguments or a machine that is not idle.
"""
import argparse
import hashlib
import json
import math
import os
import platform
import resource
import struct
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bench  # noqa: E402  -- Task 5's body, material constants and helpers

SPHERAL_COMMIT = bench.SPHERAL_COMMIT
VARIANTS = ("builtin", "builtin_repeat", "eos", "strength", "damage", "all", "linpoly_builtin", "linpoly_python")
DBL_MAX = sys.float_info.max
# Linear polynomial EOS of the linpoly variants and the method check: a1 = Murnaghan's K and the rest chosen non-zero
# (all of order K or a Gruneisen-like b0), so that every term of Spheral's polynomial does work. Not a material.
LINPOLY = dict(a0=0.0, a1=bench.MURNAGHAN["K"], a2=0.5*bench.MURNAGHAN["K"], a3=0.1*bench.MURNAGHAN["K"],
               b0=2.0, b1=0.5, b2=0.25, atomicWeight=bench.MURNAGHAN["atomicWeight"])
COLD_FIT = [0.0]*10
MELT_FIT = [bench.EMELT_PER_VOLUME] + [0.0]*9

# ------------------------------------------------------------------------------------------------- numerics
F64 = np.float64
LD = np.longdouble
if np.finfo(LD).nmant < 105:
    raise ImportError("subclass_check needs numpy longdouble with >= 106 mantissa bits (IEEE quad on aarch64 Linux) "
                      "to form exact fused multiply-adds")


def _ld(x):
    """x as IEEE quad. Every operand is converted: numpy 1.26 promotes a 0-d longdouble with a float64 array to
    float64 (value-based casting), which would silently round the product."""
    return np.asarray(x, dtype=LD)


FMA = {"mode": "exact"}     # "plain": unfused double arithmetic, as a straightforward numpy port would be


def fmadd(n, m, a):
    """a + n*m rounded once (AArch64 FMADD): the exact product (<= 106 bits) and the sum in IEEE quad (113 bits), then
    rounded to double. Checked against glibc fma() on 200,000 triples including heavy cancellation: identical."""
    if FMA["mode"] == "plain":
        return a + n*m
    return (_ld(n)*_ld(m) + _ld(a)).astype(F64)


def fmsub(n, m, a):
    """a - n*m rounded once (AArch64 FMSUB)."""
    if FMA["mode"] == "plain":
        return a - n*m
    return (_ld(a) - _ld(n)*_ld(m)).astype(F64)


def fnmsub(n, m, a):
    """n*m - a rounded once (AArch64 FNMSUB)."""
    if FMA["mode"] == "plain":
        return n*m - a
    return (_ld(n)*_ld(m) - _ld(a)).astype(F64)


def fuzzy_equal(a, b, fuzz=1.0e-15):
    """Spheral's fuzzyEqual (Utilities/SpheralFunctions.hh:30-46)."""
    return np.abs(a - b) <= fuzz*np.maximum(1.0, np.abs(a) + np.abs(b))


def sgn(x):
    return np.where(x < 0.0, -1.0, 1.0)


def cube_root(x):
    """FastMath::CubeRootHalley2 outside the Intel compiler: sgn(x) pow(|x|, 0.3333333333333333) (FastMath.hh:155)."""
    return sgn(x)*np.power(np.abs(x), 0.3333333333333333)


def poly9(c, x):
    """NinthOrderPolynomialFit::operator() (SolidMaterial/PolynomialFitInline.hh:83), in its source order. gcc contracts
    this sum into 27 fused multiply-adds per loop; that changes nothing here because every coefficient but c0 is zero
    in both fits used (a product with a zero coefficient is an exact zero), so this copy is exact only for such fits."""
    r = c[0] + c[1]*x
    p = x
    for k in range(2, 10):
        p = p*x if k > 2 else x*x
        r = r + c[k]*p
    return r


def eig_max_sym3(t):
    """GeomSymmetricTensor<3>::eigenValues().maxElement() (Geometry/GeomSymmetricTensorInline.hh:2063-2090, Eberly's
    closed form), with the fused multiply-adds gcc 13 -O3 emits for it on aarch64, transcribed from the object code of
    SteinbergGuinanStrengthInst.cc.o. t: (n, 6) as xx, xy, xz, yy, yz, zz."""
    xx, xy, xz, yy, yz, zz = (t[:, j] for j in range(6))
    tiny = 10.0*np.finfo(F64).eps
    fscale = np.maximum(tiny, np.max(np.abs(t), axis=1))
    fi = 1.0/fscale
    a00, a01, a02, a11, a12, a22 = xx*fi, xy*fi, xz*fi, yy*fi, yz*fi, zz*fi
    third = 1.0/3.0
    p0011 = a00*a11
    c1 = fmsub(a01, a01, p0011)
    c1 = fmadd(a00, a22, c1)
    c1 = fmsub(a02, a02, c1)
    c1 = fmadd(a11, a22, c1)
    c1 = fmsub(a12, a12, c1)
    c0 = ((a01 + a01)*a02)*a12
    c0 = fmadd(p0011, a22, c0)
    c0 = fmsub(a00*a12, a12, c0)
    c0 = fmsub(a02*a11, a02, c0)
    c0 = fmsub(a01*a22, a01, c0)
    c2 = (a00 + a11) + a22
    c2d3 = c2*third
    adiv3 = fmsub(c2, c2d3, c1)*third
    adiv3 = np.where(adiv3 < 0.0, adiv3, 0.0)
    k = fnmsub(c2d3, c2d3 + c2d3, c1)
    mbdiv2 = fmadd(c2d3, k, c0)*0.5
    q = fmadd(mbdiv2, mbdiv2, (adiv3*adiv3)*adiv3)
    q = np.where(q < 0.0, q, 0.0)
    mag = np.sqrt(-adiv3)
    angle = np.arctan2(np.sqrt(-q), mbdiv2)*third
    cs, sn = np.cos(angle), np.sin(angle)
    s3 = math.sqrt(3.0)
    e0 = fmadd(mag + mag, cs, c2d3)*fscale
    e1 = fmsub(fmadd(sn, s3, cs), mag, c2d3)*fscale
    e2 = fmsub(fmsub(sn, s3, cs), mag, c2d3)*fscale
    return np.maximum(e0, np.maximum(e1, e2))


def eig_max_sym2(t):
    """GeomSymmetricTensor<2>::eigenValues().maxElement() (GeomSymmetricTensorInline.hh:2043-2057), as gcc inlines it
    into SteinbergGuinanStrength<Dim<2>> (two fused operations). t: (n, 3) as xx, xy, yy."""
    xx, xy, yy = t[:, 0], t[:, 1], t[:, 2]
    b = xx + yy
    c = fnmsub(xx, yy, xy*xy)
    disc = fnmsub(b, b, c*4.0)
    s = np.where(disc > 0.0, sgn(b)*np.sqrt(np.where(disc > 0.0, disc, 0.0)), np.where(b < 0.0, -0.0, 0.0))
    q = (b + s)*0.5
    with np.errstate(divide="ignore", invalid="ignore"):
        e1 = c/q
    diag = np.abs(xy) < 1.0e-50
    return np.where(diag, np.maximum(xx, yy), np.maximum(q, e1))


def max_damage(D, ndim):
    """max(0, min(1, D.eigenValues().maxElement())) per point."""
    m = eig_max_sym3(D) if ndim == 3 else eig_max_sym2(D)
    return np.maximum(0.0, np.minimum(1.0, m))


# ------------------------------------------------------------------------------------------------- Field <-> numpy
NCOMP = {("Scalar", 2): 1, ("Scalar", 3): 1, ("Vector", 2): 2, ("Vector", 3): 3, ("SymTensor", 2): 3,
         ("SymTensor", 3): 6}


def get_internal(field, ncomp=1):
    """The internal values of a Field as a (read-only) numpy array, through Field.serialize (name, count, packed
    values): about 5 us for 1,300 doubles, where allValues() -> np.array takes 30 us. Packed elements are the raw
    C++ objects, so a 2D SymTensor is 4 doubles (xx, xy, yy, padding) and the stride is taken from the byte count."""
    b = field.serialize()
    n = field.numInternalElements
    off = 16 + len(field.name.encode())
    a = np.frombuffer(b, dtype=F64, offset=off)
    if ncomp == 1:
        return a[:n]
    return a.reshape(n, -1)[:, :ncomp]


def put_internal(field, values):
    """Write the internal values of a scalar Field (ghosts untouched) through Field.deserialize."""
    name = field.name.encode()
    v = np.ascontiguousarray(values, dtype=F64)
    field.deserialize(struct.pack("<Q", len(name)) + name + struct.pack("<Q", v.size) + v.tobytes())


def get_all(field):
    """Internal and ghost values of a scalar Field (allValues: a Python list)."""
    return np.array(field.allValues(), dtype=F64)


def put_all(field, values):
    """Assign every value of a scalar Field (Field.assign takes a list of exactly numNodes values)."""
    field.assign(np.asarray(values, dtype=F64).tolist())


class Clock:
    """Wall time spent inside the Python overrides: per method (inclusive of nested Python calls, e.g. the Python
    EOS's setTemperature inside the Python strength's shearModulus) and in total (outermost calls only)."""

    def __init__(self):
        self.t, self.n = {}, {}
        self.total = 0.0
        self.depth = 0

    def reset(self):
        self.t.clear()
        self.n.clear()
        self.total = 0.0


CLOCK = Clock()


def timed(fn):
    key = fn.__name__

    def wrapper(*args):
        CLOCK.depth += 1
        t0 = time.perf_counter()
        try:
            return fn(*args)
        finally:
            dt = time.perf_counter() - t0
            CLOCK.depth -= 1
            CLOCK.t[key] = CLOCK.t.get(key, 0.0) + dt
            CLOCK.n[key] = CLOCK.n.get(key, 0) + 1
            if CLOCK.depth == 0:
                CLOCK.total += dt
    wrapper.__name__ = key
    wrapper.__doc__ = fn.__doc__
    return wrapper


# ------------------------------------------------------------------------------------------------- the subclasses
def make_classes(S, ndim):
    """The four Python classes for one Spheral dimension module (Spheral3d or SpheralRZ, whose classes are Dim<2>)."""

    class _SolidEOSBase(S.SolidEquationOfState):
        def __init__(self, referenceDensity, etamin, etamax, constants, externalPressure=0.0,
                     minimumPressure=-DBL_MAX, maximumPressure=DBL_MAX, minimumPressureDamage=0.0,
                     minPressureType=None):
            if minPressureType is None:
                minPressureType = S.MaterialPressureMinType.PressureFloor
            S.SolidEquationOfState.__init__(self, referenceDensity, etamin, etamax, constants, minimumPressure,
                                            maximumPressure, minimumPressureDamage, minPressureType,
                                            externalPressure)
            self._constants = constants     # the C++ base keeps a reference to it
            self._R = constants.molarGasConstant
            self._floor = S.MaterialPressureMinType.PressureFloor

        def bounded_eta(self, rho):
            """SolidEquationOfState::boundedEta (SolidEquationOfStateInline.hh:123)."""
            return np.maximum(self.etamin, np.minimum(self.etamax, rho/self.referenceDensity))

        def limits(self, P):
            """EquationOfState::applyPressureLimits after the external pressure is taken off
            (Material/EquationOfStateInline.hh:98-106)."""
            pmin, pmax = self.minimumPressure, self.maximumPressure
            low = pmin if self.minimumPressureType == self._floor else 0.0
            return np.where(P < pmin, low, np.where(P > pmax, pmax, P))

        def valid(self):
            return bool(self.referenceDensity > 0.0 and self.etamin <= self.etamax)

    class PyMurnaghan(_SolidEOSBase):
        """MurnaghanEquationOfState (SolidMaterial/MurnaghanEquationOfState.cc at 116c71f) on whole arrays."""

        def __init__(self, referenceDensity, etamin, etamax, n, K, atomicWeight, constants, **kw):
            _SolidEOSBase.__init__(self, referenceDensity, etamin, etamax, constants, **kw)
            self.n, self.K, self.atomicWeight = n, K, atomicWeight
            self.cv = 3.0*self._R/atomicWeight
            self.nKi = K/n

        def _pressure(self, rho):
            eta = self.bounded_eta(rho)
            x = np.power(eta, self.n) - 1.0
            ext = self.externalPressure
            # gcc fuses mnKi*(eta^n - 1) - Pext into one FNMSUB (pressure(), MurnaghanEquationOfState.cc:197); with
            # Pext = 0 that is exactly the plain product.
            P = self.nKi*x if ext == 0.0 else fnmsub(self.nKi, x, ext)
            return np.where(fuzzy_equal(eta, self.etamin), 0.0, self.limits(P))

        def _dpdrho(self, rho):
            eta = self.bounded_eta(rho)
            r = np.maximum(0.0, np.power(eta, self.n - 1)*(self.K/self.referenceDensity))
            return np.where(fuzzy_equal(eta, self.etamin) | fuzzy_equal(eta, self.etamax), 0.0, r)

        def _gamma(self, rho):
            nden = (self.referenceDensity*self.bounded_eta(rho))/self.atomicWeight
            return 1.0 + self._R*nden/self.cv

        @timed
        def setPressure(self, Pressure, massDensity, specificThermalEnergy):
            put_all(Pressure, self._pressure(get_all(massDensity)))

        @timed
        def setPressureAndDerivs(self, Pressure, dPdu, dPdrho, massDensity, specificThermalEnergy):
            rho = get_all(massDensity)
            eta = self.bounded_eta(rho)
            put_all(Pressure, self._pressure(rho))
            put_all(dPdu, np.zeros_like(rho))
            put_all(dPdrho, self.n*self.nKi*np.power(eta, self.n - 1.0)/self.referenceDensity)

        @timed
        def setTemperature(self, temperature, massDensity, specificThermalEnergy):
            put_all(temperature, get_all(specificThermalEnergy)/self.cv + 300.0)

        @timed
        def setSpecificThermalEnergy(self, specificThermalEnergy, massDensity, temperature):
            put_all(specificThermalEnergy, (get_all(temperature) - 300.0)*self.cv)

        @timed
        def setSpecificHeat(self, specificHeat, massDensity, temperature):
            put_all(specificHeat, np.full(specificHeat.numElements, self.cv))

        @timed
        def setSoundSpeed(self, soundSpeed, massDensity, specificThermalEnergy):
            put_all(soundSpeed, np.sqrt(self._dpdrho(get_all(massDensity))))

        @timed
        def setGammaField(self, gamma, massDensity, specificThermalEnergy):
            put_all(gamma, self._gamma(get_all(massDensity)))

        @timed
        def setBulkModulus(self, bulkModulus, massDensity, specificThermalEnergy):
            rho = get_all(massDensity)
            put_all(bulkModulus, rho*self._dpdrho(rho))

        @timed
        def setEntropy(self, entropy, massDensity, specificThermalEnergy):
            rho = get_all(massDensity)
            with np.errstate(over="ignore"):
                x = np.power(rho, self._gamma(rho))
            put_all(entropy, self._pressure(rho)*(sgn(x)/np.maximum(1.0e-30, np.abs(x))))

    class PyLinearPolynomial(_SolidEOSBase):
        """LinearPolynomialEquationOfState (SolidMaterial/LinearPolynomialEquationOfState.cc at 116c71f), with the
        fused multiply-adds of its pressureAndDerivs (10) and computeDPDrho (7) transcribed from the object code."""

        def __init__(self, referenceDensity, etamin, etamax, a0, a1, a2, a3, b0, b1, b2, atomicWeight, constants,
                     **kw):
            _SolidEOSBase.__init__(self, referenceDensity, etamin, etamax, constants, **kw)
            self.a0, self.a1, self.a2, self.a3 = a0, a1, a2, a3
            self.b0, self.b1, self.b2 = b0, b1, b2
            self.atomicWeight = atomicWeight
            self.cv = 3.0*self._R/atomicWeight
            self.gamma0 = b0 + 1.0

        def _pad(self, rho, eps):
            """pressureAndDerivs: (P, dP/du, dP/drho), zero where eta is at etamin."""
            eta = self.bounded_eta(rho)
            mu = eta - 1.0
            a0, a1, a2, a3, b0, b1, b2 = self.a0, self.a1, self.a2, self.a3, self.b0, self.b1, self.b2
            p = fmadd(a1, mu, a0)
            p = fmadd(a2*mu, mu, p)
            p = fmadd((a3*mu)*mu, mu, p)
            dpdu = fmadd(b1, mu, b0)
            dpdu = fmadd(b2*mu, mu, dpdu)
            p = fmadd(dpdu, eps, p)
            P = self.limits(p - self.externalPressure)
            d = fmadd(a2 + a2, mu, a1)
            d = fmadd((a3*3.0)*mu, mu, d)
            d = fmadd(fmadd(b2 + b2, mu, b1), eps, d)
            dpdrho = d/self.referenceDensity
            at_min = fuzzy_equal(eta, self.etamin)
            return (np.where(at_min, 0.0, P), np.where(at_min, 0.0, dpdu), np.where(at_min, 0.0, dpdrho))

        def _dpdrho(self, rho, eps):
            """computeDPDrho."""
            rho0 = self.referenceDensity
            eta = self.bounded_eta(rho)
            mu = eta - 1.0
            a1, a2, a3, b0, b1, b2 = self.a1, self.a2, self.a3, self.b0, self.b1, self.b2
            d = fmadd(a2 + a2, mu, a1)
            d = fmadd((a3*3.0)*mu, mu, d)
            d = fmadd(fmadd(b2 + b2, mu, b1), eps, d)
            dpdrho_eps = np.abs(d)/rho0
            r = rho0*eta
            prho2 = self._pad(rho, eps)[0]/(r*r)
            dpdeps = fmadd(mu*b2, mu, fmadd(b1, mu, b0))
            return np.maximum(0.0, fmadd(prho2, dpdeps, dpdrho_eps))

        @timed
        def setPressure(self, Pressure, massDensity, specificThermalEnergy):
            put_all(Pressure, self._pad(get_all(massDensity), get_all(specificThermalEnergy))[0])

        @timed
        def setPressureAndDerivs(self, Pressure, dPdu, dPdrho, massDensity, specificThermalEnergy):
            P, du, dr = self._pad(get_all(massDensity), get_all(specificThermalEnergy))
            put_all(Pressure, P)
            put_all(dPdu, du)
            put_all(dPdrho, dr)

        @timed
        def setTemperature(self, temperature, massDensity, specificThermalEnergy):
            put_all(temperature, get_all(specificThermalEnergy)/self.cv + 300.0)

        @timed
        def setSpecificThermalEnergy(self, specificThermalEnergy, massDensity, temperature):
            put_all(specificThermalEnergy, (get_all(temperature) - 300.0)*self.cv)

        @timed
        def setSpecificHeat(self, specificHeat, massDensity, temperature):
            put_all(specificHeat, np.full(specificHeat.numElements, self.cv))

        @timed
        def setSoundSpeed(self, soundSpeed, massDensity, specificThermalEnergy):
            put_all(soundSpeed, np.sqrt(self._dpdrho(get_all(massDensity), get_all(specificThermalEnergy))))

        @timed
        def setGammaField(self, gamma, massDensity, specificThermalEnergy):
            put_all(gamma, np.full(gamma.numElements, self.gamma0))

        @timed
        def setBulkModulus(self, bulkModulus, massDensity, specificThermalEnergy):
            rho = get_all(massDensity)
            put_all(bulkModulus, rho*self._dpdrho(rho, get_all(specificThermalEnergy)))

        @timed
        def setEntropy(self, entropy, massDensity, specificThermalEnergy):
            rho = get_all(massDensity)
            x = np.power(rho, self.gamma0)
            P = self._pad(rho, get_all(specificThermalEnergy))[0]
            put_all(entropy, P*(sgn(x)/np.maximum(1.0e-30, np.abs(x))))

    class PySteinbergGuinan(S.StrengthModel):
        """SteinbergGuinanStrength (SolidMaterial/SteinbergGuinanStrength.cc at 116c71f) on whole arrays, with the
        fused multiply-adds gcc emits in its loops (two in shearModulus, one each in yieldStrength and soundSpeed,
        and those of the damage eigenvalue, see eig_max_sym3/2). The fits are passed as coefficient lists."""

        def __init__(self, eos, G0, A, B, Y0, Ymax, Yp, beta, gamma0, nhard, coldFit, meltFit, Gmax=1.0e100):
            S.StrengthModel.__init__(self)
            self.eos = eos
            self.G0, self.Gmax, self.A, self.B = G0, Gmax, A, B
            self.Y0, self.Ymax, self.Yp, self.beta, self.gamma0, self.nhard = Y0, Ymax, Yp, beta, gamma0, nhard
            self.cold, self.melt = list(coldFit), list(meltFit)
            self.constant = G0 > 0.0 and A == 0.0 and B == 0.0 and beta == 0.0 and nhard == 0.0
            self.ScalarField = S.ScalarField

        def providesSoundSpeed(self):
            return True

        def providesBulkModulus(self):
            return False

        def _eta(self, rho):
            e = self.eos
            return np.maximum(e.etamin, np.minimum(e.etamax, rho/e.referenceDensity))

        def _temperature(self, density, rho, eps):
            """computeTemperature: the EOS's temperature of eps - cold energy, less 300 K (internal values; the
            C++ builds the same two scratch Fields and calls the EOS through its virtual setTemperature)."""
            nl = density.nodeList()
            eps1 = self.ScalarField("new energy", nl)
            T = self.ScalarField("temperature", nl)
            put_internal(eps1, eps - poly9(self.cold, self._eta(rho) - 1.0)/self.eos.referenceDensity)
            self.eos.setTemperature(T, density, eps1)
            return get_internal(T) - 300.0

        def _fmelt(self, rho, eps):
            """meltAttenuation."""
            emelt = poly9(self.melt, self._eta(rho) - 1.0)/self.eos.referenceDensity
            melted = (eps > emelt) | fuzzy_equal(eps, emelt)
            with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
                f = np.exp((-self.Yp*eps)/(emelt - eps))
            return np.where(melted, 0.0, f)

        def _shear(self, density, specificThermalEnergy, pressure, damage):
            rho = get_internal(density)
            eps = get_internal(specificThermalEnergy)
            D = max_damage(get_internal(damage, NCOMP[("SymTensor", ndim)]), ndim)
            if self.constant:
                return (1.0 - D)*self.G0, D
            T = self._temperature(density, rho, eps)
            P = get_internal(pressure)
            fmelt = self._fmelt(rho, eps)
            x = fmsub(self.B, T, (self.A*P)/cube_root(self._eta(rho)) + 1.0)
            G = np.minimum(self.Gmax, self.G0*np.maximum(1.0e-10, fmelt*x))
            return fmadd(self.G0, fmelt*D, G*(1.0 - D)), D

        @timed
        def shearModulus(self, shearModulus, density, specificThermalEnergy, pressure, damage):
            put_internal(shearModulus, self._shear(density, specificThermalEnergy, pressure, damage)[0])

        @timed
        def yieldStrength(self, yieldStrength, density, specificThermalEnergy, pressure, plasticStrain,
                          plasticStrainRate, damage):
            G, D = self._shear(density, specificThermalEnergy, pressure, damage)
            if self.constant:
                put_internal(yieldStrength, (1.0 - D)*self.Y0)
                return
            ps = get_internal(plasticStrain)
            Yhard = np.minimum(self.Ymax, self.Y0*np.power(fmadd(self.beta, ps + self.gamma0, 1.0), self.nhard))
            put_internal(yieldStrength, ((Yhard*(1.0 - D))*G)/self.G0)

        @timed
        def soundSpeed(self, soundSpeed, density, specificThermalEnergy, pressure, fluidSoundSpeed, damage):
            mu = self._shear(density, specificThermalEnergy, pressure, damage)[0]
            cf = get_internal(fluidSoundSpeed)
            put_internal(soundSpeed, np.sqrt(fmadd(cf, cf, np.abs(mu*(4.0/3.0))/get_internal(density))))

        @timed
        def meltSpecificEnergy(self, meltSpecificEnergy, density, specificThermalEnergy):
            put_internal(meltSpecificEnergy,
                         poly9(self.melt, self._eta(get_internal(density)) - 1.0)/self.eos.referenceDensity)

        @timed
        def coldSpecificEnergy(self, coldSpecificEnergy, density, specificThermalEnergy):
            put_internal(coldSpecificEnergy,
                         poly9(self.cold, self._eta(get_internal(density)) - 1.0)/self.eos.referenceDensity)

    class PyDelegatingDamage(S.DamageModel):
        """A DamageModel subclass that holds a C++ ProbabilisticDamageModel (not given to the integrator) and forwards
        every Physics/DamageModel virtual to it, except computeScalarDDDt/evaluateDerivatives, which it computes
        itself in numpy: DamageModel::computeScalarDDDt (Damage/DamageModel.cc:83-122) as called by
        ProbabilisticDamageModel::evaluateDerivatives (ProbabilisticDamageModel.cc:229-244)."""

        def __init__(self, delegate):
            S.DamageModel.__init__(self, delegate.nodeList, delegate.kernel, delegate.crackGrowthMultiplier,
                                   delegate.damageCouplingAlgorithm)
            self.delegate = delegate

        # -- the method the tearing law replaces
        @timed
        def computeScalarDDDt(self, dataBase, state, time, dt, DDDt):
            d = self.delegate
            nodes = d.nodeList
            n = nodes.numInternalNodes
            if d.freezeDamage:
                put_internal(DDDt, np.zeros(n))
                return
            name = nodes.name
            cl = get_internal(state.scalarField(state.buildFieldKey(SolidFieldNames.longitudinalSoundSpeed, name)))
            H = get_internal(state.symTensorField(state.buildFieldKey(HydroFieldNames.H, name)),
                             NCOMP[("SymTensor", ndim)])
            A = d.crackGrowthMultiplier/d.kernel.kernelExtent
            if ndim == 3:
                xx, xy, xz, yy, yz, zz = (H[:, j] for j in range(6))
                det = (xy*yz)*xz
                det = fmadd(xx*yy, zz, det)
                det = fmadd(yz, xy*xz, det)
                det = fmsub(yz, xx*yz, det)
                det = fmsub(zz, xy*xy, det)
                det = fmsub(xz, yy*xz, det)
                hr = cube_root(det)
            else:
                det = fnmsub(H[:, 0], H[:, 2], H[:, 1]*H[:, 1])
                with np.errstate(invalid="ignore"):
                    hr = np.sqrt(det)
            r = (A*cl)*hr
            excl = d.excludeNodes
            if excl:
                r = r.copy()
                r[np.asarray(excl, dtype=np.int64)] = 0.0
            put_internal(DDDt, r)

        @timed
        def evaluateDerivatives(self, time, dt, dataBase, state, derivs):
            name = self.delegate.nodeList.name
            DDDt = derivs.scalarField(state.buildFieldKey("delta " + SolidFieldNames.scalarDamage, name))
            self.computeScalarDDDt(dataBase, state, time, dt, DDDt)

        # -- everything else is the delegate's (timed too, so the cost of forwarding is visible)
        def _sync_boundaries(self):
            """Physics::appendBoundary is not virtual (Physics/Physics.hh:98), so the boundaries SpheralController
            gives this package (in parallel at least the distributed boundary, SpheralController.py:853-869) never
            reach the delegate, whose applyGhostBoundaries would then fill no ghost damage. Copy them over."""
            mine = list(self.boundaryConditions)
            d = self.delegate
            if len(mine) != len(d.boundaryConditions) or not all(d.haveBoundary(bc) for bc in mine):
                d.clearBoundaries()
                for bc in mine:
                    d.appendBoundary(bc)

        @timed
        def label(self):
            return self.delegate.label()

        @timed
        def dt(self, dataBase, state, derivs, currentTime):
            return self.delegate.dt(dataBase, state, derivs, currentTime)

        @timed
        def dtImplicit(self, dataBase, state, derivs, currentTime):
            return self.delegate.dtImplicit(dataBase, state, derivs, currentTime)

        @timed
        def registerState(self, dataBase, state):
            self._sync_boundaries()
            self.delegate.registerState(dataBase, state)

        @timed
        def registerDerivatives(self, dataBase, derivs):
            self.delegate.registerDerivatives(dataBase, derivs)

        @timed
        def initializeProblemStartup(self, dataBase):
            self._sync_boundaries()
            self.delegate.initializeProblemStartup(dataBase)

        @timed
        def initializeProblemStartupDependencies(self, dataBase, state, derivs):
            self._sync_boundaries()
            self.delegate.initializeProblemStartupDependencies(dataBase, state, derivs)

        @timed
        def preStepInitialize(self, dataBase, state, derivs):
            self.delegate.preStepInitialize(dataBase, state, derivs)

        @timed
        def initialize(self, time, dt, dataBase, state, derivs):
            self._sync_boundaries()
            return self.delegate.initialize(time, dt, dataBase, state, derivs)

        @timed
        def finalizeDerivatives(self, time, dt, dataBase, state, derivs):
            self.delegate.finalizeDerivatives(time, dt, dataBase, state, derivs)

        @timed
        def postStateUpdate(self, time, dt, dataBase, state, derivs):
            return self.delegate.postStateUpdate(time, dt, dataBase, state, derivs)

        @timed
        def finalize(self, time, dt, dataBase, state, derivs):
            self._sync_boundaries()
            self.delegate.finalize(time, dt, dataBase, state, derivs)

        @timed
        def applyGhostBoundaries(self, state, derivs):
            self._sync_boundaries()
            self.delegate.applyGhostBoundaries(state, derivs)

        @timed
        def enforceBoundaries(self, state, derivs):
            self._sync_boundaries()
            self.delegate.enforceBoundaries(state, derivs)

        @timed
        def requireConnectivity(self):
            return self.delegate.requireConnectivity()

        @timed
        def requireVolumes(self):
            return self.delegate.requireVolumes()

        @timed
        def requireReproducingKernels(self):
            return self.delegate.requireReproducingKernels()

        @timed
        def maxResidual(self, dataBase, state1, state0, tol):
            return self.delegate.maxResidual(dataBase, state1, state0, tol)

        @timed
        def extraEnergy(self):
            return self.delegate.extraEnergy()

        @timed
        def extraMomentum(self):
            return self.delegate.extraMomentum()

        @timed
        def registerAdditionalVisualizationState(self, dataBase, state):
            self.delegate.registerAdditionalVisualizationState(dataBase, state)

        @timed
        def dumpState(self, file, pathName):
            self.delegate.dumpState(file, pathName)

        @timed
        def restoreState(self, file, pathName):
            self.delegate.restoreState(file, pathName)

    SolidFieldNames, HydroFieldNames = S.SolidFieldNames, S.HydroFieldNames
    return dict(PyMurnaghan=PyMurnaghan, PyLinearPolynomial=PyLinearPolynomial, PySteinbergGuinan=PySteinbergGuinan,
                PyDelegatingDamage=PyDelegatingDamage)


# ------------------------------------------------------------------------------------------------- material builders
def murnaghan(S, classes, python, units):
    m = bench.MURNAGHAN
    cls = classes["PyMurnaghan"] if python else S.MurnaghanEquationOfState
    return cls(bench.RHO0, bench.ETAMIN, bench.ETAMAX, m["n"], m["K"], m["atomicWeight"], units)


def linpoly(S, classes, python, units):
    p = LINPOLY
    cls = classes["PyLinearPolynomial"] if python else S.LinearPolynomialEquationOfState
    return cls(bench.RHO0, bench.ETAMIN, bench.ETAMAX, p["a0"], p["a1"], p["a2"], p["a3"], p["b0"], p["b1"], p["b2"],
               p["atomicWeight"], units)


def steinberg_guinan(S, classes, python, eos):
    g = bench.SG
    args = (eos, g["G0"], g["A"], g["B"], g["Y0"], g["Ymax"], g["Yp"], g["beta"], g["gamma0"], g["nhard"])
    if python:
        return classes["PySteinbergGuinan"](*args, COLD_FIT, MELT_FIT)
    return S.SteinbergGuinanStrength(*args, S.NinthOrderPolynomialFit(*COLD_FIT), S.NinthOrderPolynomialFit(*MELT_FIT))


def variant_classes(variant):
    """(python EOS, python strength, python damage, EOS kind) for a variant."""
    return {"builtin": (False, False, False, "murnaghan"), "builtin_repeat": (False, False, False, "murnaghan"),
            "eos": (True, False, False, "murnaghan"), "strength": (False, True, False, "murnaghan"),
            "damage": (False, False, True, "murnaghan"), "all": (True, True, True, "murnaghan"),
            "linpoly_builtin": (False, False, False, "linpoly"),
            "linpoly_python": (True, False, False, "linpoly")}[variant]


# ------------------------------------------------------------------------------------------------- run
def parse_run_args(argv):
    p = argparse.ArgumentParser(prog="subclass_check.py run", description="one run of the comparison")
    p.add_argument("--geometry", choices=("3d", "rz"), required=True)
    p.add_argument("--dx-mm", type=float, required=True)
    p.add_argument("--variant", choices=VARIANTS, required=True)
    p.add_argument("--steps", type=int, default=100)
    p.add_argument("--warmup", type=int, default=10)
    p.add_argument("--out", default="spheral_output/m0/subclass")
    p.add_argument("--radius-mm", type=float, default=50.0)
    p.add_argument("--max-load", type=float, default=2.0)
    p.add_argument("--force", action="store_true")
    p.add_argument("--host-load", default=None)
    p.add_argument("--fma", choices=("exact", "plain"), default="exact",
                   help="exact: the Python classes repeat gcc's fused multiply-adds (bitwise copies); plain: unfused "
                        "numpy arithmetic, as a straightforward port would be (cost of production-style Python)")
    a = p.parse_args(argv)
    if a.dx_mm <= 0 or a.radius_mm <= 0 or a.steps < 1 or a.warmup < 0 or a.dx_mm > a.radius_mm:
        p.error("need 0 < dx <= radius, steps >= 1, warmup >= 0")
    return a


def run_name(a, nproc):
    return (f"subclass_{a.geometry}_dx{a.dx_mm:.3f}mm_R{a.radius_mm:g}mm_n{nproc}_{a.variant}_fma{a.fma}_dmgon"
            f"_steps{a.steps}_wu{a.warmup}_{SPHERAL_COMMIT}")


def lattice_ids(pos, geometry, R, dx):
    """A global particle id from the lattice indices of the starting position (bench.lattice_generator's lattice):
    independent of the MPI distribution, so runs compare particle by particle."""
    k = int(math.floor(R/dx*(1.0 + 1e-12)))
    if geometry == "3d":
        idx = np.rint(pos/dx).astype(np.int64) + k
        n1 = 2*k + 1
        return (idx[:, 0]*n1 + idx[:, 1])*n1 + idx[:, 2]
    nr = int(math.ceil(R/dx - 0.5)) + 1
    iz = np.rint(pos[:, 0]/dx).astype(np.int64) + k
    jr = np.rint(pos[:, 1]/dx - 0.5).astype(np.int64)
    return iz*nr + jr


def run(argv):
    try:
        a = parse_run_args(argv)
    except SystemExit as e:
        return 2 if e.code else 0
    nproc = bench.mpi_size_from_env()
    rank = int(os.environ.get("OMPI_COMM_WORLD_RANK", os.environ.get("PMI_RANK", "0")))
    name = run_name(a, nproc)
    rundir = os.path.join(a.out, name)
    if os.path.exists(os.path.join(rundir, "summary.json")):
        if rank == 0:
            print(f"subclass_check: {name} exists, skipped")
        return 0

    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    load = comm.bcast(os.getloadavg() if comm.rank == 0 else None, root=0)
    if load[0] > a.max_load and not a.force:
        if comm.rank == 0:
            print(f"subclass_check: load average {load[0]:.2f} > {a.max_load} (the Linux VM's); not starting "
                  "(--force overrides)")
        return 2

    import mpi
    if a.geometry == "3d":
        import Spheral3d as S
        from PeanoHilbertDistributeNodes import distributeNodes3d as distributeNodes
        ndim = 3
    else:
        import SpheralRZ as S
        from PeanoHilbertDistributeNodes import distributeNodes2d as distributeNodes
        ndim = 2
    from SpheralController import SpheralController
    FMA["mode"] = a.fma
    classes = make_classes(S, ndim)
    py_eos, py_strength, py_damage, eos_kind = variant_classes(a.variant)

    R, dx = 1e-3*a.radius_mm, 1e-3*a.dx_mm
    units = S.MKS()
    eos = (murnaghan if eos_kind == "murnaghan" else linpoly)(S, classes, py_eos, units)
    strength = steinberg_guinan(S, classes, py_strength, eos)
    WT = S.TableKernel(S.WendlandC4Kernel(), 1000)
    big = 10.0*R
    box = (S.Vector(-big, -big, -big), S.Vector(big, big, big)) if ndim == 3 else (S.Vector(-big, -big),
                                                                                    S.Vector(big, big))
    nodes = S.makeSolidNodeList("AA7075", eos, strength, nPerh=bench.NPERH, kernelExtent=WT.kernelExtent,
                                hmin=1e-3*dx, hmax=100.0*dx, rhoMin=bench.ETAMIN*bench.RHO0,
                                rhoMax=bench.ETAMAX*bench.RHO0, xmin=box[0], xmax=box[1])
    gen, expected = bench.lattice_generator(a.geometry, R, dx)
    distributeNodes((nodes, gen))
    nlocal = nodes.numInternalNodes
    ntotal = mpi.allreduce(nlocal, mpi.SUM)
    nodes.specificThermalEnergy(S.ScalarField("eps0", nodes, 0.0))
    pos = nodes.positions()
    vel = nodes.velocity()
    for i in range(nlocal):
        vel[i] = pos[i]*bench.EDOT
    gid = lattice_ids(get_internal(pos, ndim).copy(), a.geometry, R, dx)
    all_gid = comm.gather(gid, root=0)
    fingerprint = None
    if comm.rank == 0:
        h = hashlib.sha256()
        for g in all_gid:
            h.update(np.ascontiguousarray(g).tobytes() + b"|")
        fingerprint = h.hexdigest()[:16]
        if len(np.unique(np.concatenate(all_gid))) != ntotal:
            print("subclass_check: lattice ids are not unique")
            return 1

    # As bench.main: SolidSPH, CheapSynchronousRK2, damage on.
    db = S.DataBase()
    db.appendNodeList(nodes)
    hydro = S.SPH(dataBase=db, W=WT, cfl=bench.CFL, compatibleEnergyEvolution=True, densityUpdate=S.IntegrateDensity,
                  XSPH=False, ASPH=False)
    integrator = S.CheapSynchronousRK2Integrator(db)
    integrator.appendPhysicsPackage(hydro)
    pdm = S.ProbabilisticDamageModel(materialName="aluminum", units=units, nodeList=nodes, kernel=WT,
                                     seed=bench.SEED)
    damage_pkg = classes["PyDelegatingDamage"](pdm) if py_damage else pdm
    integrator.appendPhysicsPackage(damage_pkg)
    cs0 = math.sqrt((bench.MURNAGHAN["K"] + (4.0/3.0)*bench.SG["G0"])/bench.RHO0)
    integrator.lastDt = 0.1*bench.CFL*dx/cs0
    integrator.dtMin = 1e-15
    integrator.dtMax = 1.0
    integrator.verbose = False
    controller = SpheralController(integrator, WT, statsStep=10**9, printStep=10**9, redistributeStep=None,
                                   restartStep=None, vizStep=None, vizTime=None, SPH=True)

    try:
        if a.warmup:
            controller.advance(1e40, a.warmup)
        mpi.barrier()
        CLOCK.reset()
        walls, dts = [], []
        for _ in range(a.steps):
            mpi.barrier()
            ts = time.perf_counter()
            integrator.step(1e40)
            controller.totalSteps += 1
            controller.doPeriodicWork()
            mpi.barrier()
            walls.append(time.perf_counter() - ts)
            dts.append(integrator.lastDt)
    except Exception as e:  # noqa: BLE001 -- any Spheral failure is a model failure (exit 1)
        print(f"subclass_check: Spheral failed on rank {mpi.rank}: {type(e).__name__}: {e}")
        return 1

    if nodes.numInternalNodes != nlocal:
        print("subclass_check: the particle count changed during the run")
        return 1
    fields = dict(
        position=get_internal(nodes.positions(), ndim), velocity=get_internal(nodes.velocity(), ndim),
        deviatoric_stress=get_internal(nodes.deviatoricStress(), NCOMP[("SymTensor", ndim)]),
        pressure=get_internal(hydro.pressure[0]), damage=get_internal(nodes.damage(), NCOMP[("SymTensor", ndim)]),
        density=get_internal(nodes.massDensity()), specific_thermal_energy=get_internal(nodes.specificThermalEnergy()),
        plastic_strain=get_internal(nodes.plasticStrain()))
    gathered = {k: comm.gather(np.array(v), root=0) for k, v in fields.items()}
    python_s = comm.gather(dict(CLOCK.t), root=0)
    python_n = comm.gather(dict(CLOCK.n), root=0)
    python_total = comm.gather(CLOCK.total, root=0)
    rss_total = mpi.allreduce(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss, mpi.SUM)*1024
    nmin = mpi.allreduce(nlocal, mpi.MIN)
    nmax = mpi.allreduce(nlocal, mpi.MAX)
    if comm.rank != 0:
        return 0

    order = np.argsort(np.concatenate(all_gid), kind="stable")
    arrays = {"gid": np.concatenate(all_gid)[order]}
    for k, parts in gathered.items():
        arrays[k] = np.concatenate(parts)[order]
    keys = sorted(set().union(*python_s))
    py_time = {k: dict(seconds_max_rank=max(d.get(k, 0.0) for d in python_s),
                       seconds_mean_rank=sum(d.get(k, 0.0) for d in python_s)/len(python_s),
                       calls_rank0=python_n[0].get(k, 0)) for k in keys}
    w = bench.stats(walls)
    summary = dict(
        run=name, spheral_commit=SPHERAL_COMMIT, geometry=a.geometry, dx_mm=a.dx_mm, radius_mm=a.radius_mm,
        variant=a.variant, fma=a.fma, python_eos=py_eos, python_strength=py_strength, python_damage=py_damage, eos=eos_kind,
        particles=ntotal, particles_expected=expected, particles_per_process_min=nmin,
        particles_per_process_max=nmax, processes=mpi.procs, omp_num_threads=os.environ.get("OMP_NUM_THREADS"),
        decomposition_fingerprint=fingerprint, steps=a.steps, warmup=a.warmup, step_wall_s=w,
        step_wall_total_s=sum(walls), dt_s=bench.stats(dts), time_reached_s=controller.time(),
        python_override_time_during_timed_steps=py_time,
        python_override_s_per_step_max_rank=max(python_total)/a.steps,
        python_override_s_per_step_mean_rank=sum(python_total)/len(python_total)/a.steps,
        peak_rss_bytes_total=rss_total, load_average_at_start=list(load),
        host_load_average_at_start=[float(x) for x in a.host_load.split()] if a.host_load else None,
        forced=bool(a.force and load[0] > a.max_load), machine=platform.machine(),
        material="timing only (bench.py's)", pin=json.load(open(bench.PIN)))
    os.makedirs(rundir, exist_ok=True)
    np.savez(os.path.join(rundir, "state.npz"), **arrays)
    with open(os.path.join(rundir, "steps.csv"), "w") as f:
        f.write("step,wall_s,dt_s\n")
        for i, (wt, d) in enumerate(zip(walls, dts)):
            f.write(f"{i + 1},{wt:.6e},{d:.6e}\n")
    with open(os.path.join(rundir, "summary.json"), "w") as f:
        json.dump(summary, f, indent=2)
        f.write("\n")
    print(f"subclass_check: {name}: {ntotal} particles on {mpi.procs} processes, median step {w['median']:.4g} s, "
          f"python overrides {summary['python_override_s_per_step_max_rank']*1e3:.3g} ms per step (slowest rank)")
    return 0


# ------------------------------------------------------------------------------------------------- comparison
def ordered_int(x):
    """float64 -> int64 whose difference counts the doubles between two values (the ULP distance)."""
    i = np.ascontiguousarray(x, dtype=F64).view(np.int64)
    return np.where(i < 0, np.int64(-0x8000000000000000) - i, i)


def field_diff(a, b):
    """Differences of b from a, point by point (rows) over every component: whether bitwise identical, how many points
    differ, the largest absolute difference, that over the field's largest magnitude, the largest per-point
    relative difference |b_i - a_i|/|a_i| (vector/tensor norms), and the largest distance in units in the last place
    over components of equal sign (a change of sign, or to or from zero, is counted separately: its ULP distance
    says nothing)."""
    a2 = np.asarray(a, F64).reshape(len(a), -1)
    b2 = np.asarray(b, F64).reshape(len(b), -1)
    d = np.abs(a2 - b2)
    scale = float(np.max(np.abs(a2))) if a2.size else 0.0
    na = np.sqrt(np.sum(a2*a2, axis=1))
    nd = np.sqrt(np.sum(d*d, axis=1))
    with np.errstate(divide="ignore", invalid="ignore"):
        rel = np.where(na > 0.0, nd/na, np.where(nd > 0.0, np.inf, 0.0))
    same_sign = (np.sign(a2) == np.sign(b2)) & (a2 != 0.0)
    ulp = np.abs(ordered_int(a2) - ordered_int(b2))[same_sign]
    return dict(identical=bool(np.array_equal(a2, b2)), points=int(len(a2)),
                points_differing=int(np.count_nonzero(np.any(d != 0, axis=1))),
                max_abs=float(np.max(d)) if d.size else 0.0, field_max_abs=scale,
                max_abs_over_field_max=float(np.max(d))/scale if scale > 0 else 0.0,
                max_pointwise_relative=float(np.max(rel)) if rel.size else 0.0,
                max_ulp_same_sign=int(np.max(ulp)) if ulp.size else 0,
                components_sign_or_zero_changed=int(np.count_nonzero(~same_sign & (a2 != b2))))


COMPARED = ("position", "velocity", "deviatoric_stress", "pressure", "damage", "density", "specific_thermal_energy",
            "plastic_strain")
BASELINE = {"builtin_repeat": "builtin", "eos": "builtin", "strength": "builtin", "damage": "builtin",
            "all": "builtin", "linpoly_python": "linpoly_builtin"}
ORDER = ("builtin_repeat", "eos", "strength", "damage", "all", "linpoly_python")


def run_label(s):
    """The variant, marked when its Python classes ran with unfused arithmetic."""
    python = s["python_eos"] or s["python_strength"] or s["python_damage"]
    return s["variant"] + ("+plainfma" if python and s.get("fma") == "plain" else "")


def compare(argv):
    p = argparse.ArgumentParser(prog="subclass_check.py compare")
    p.add_argument("--out", default="spheral_output/m0/subclass")
    p.add_argument("--summary", default=None, help="default OUT/summary.json")
    try:
        a = p.parse_args(argv)
    except SystemExit as e:
        return 2 if e.code else 0
    runs = {}
    for d in sorted(os.listdir(a.out)):
        sp = os.path.join(a.out, d, "summary.json")
        if os.path.isfile(sp) and d.startswith("subclass_"):
            s = json.load(open(sp))
            key = (s["geometry"], s["dx_mm"], s["processes"], s["steps"], s["warmup"], s["radius_mm"])
            runs.setdefault(key, {})[run_label(s)] = (s, os.path.join(a.out, d))
    out = dict(spheral_commit=SPHERAL_COMMIT, cases=[], methods={})
    for g in ("3d", "rz"):
        mp = os.path.join(a.out, f"methods_{g}.json")
        if os.path.isfile(mp):
            out["methods"][g] = json.load(open(mp))
    for key, vs in sorted(runs.items()):
        case = dict(geometry=key[0], dx_mm=key[1], processes=key[2], steps=key[3], warmup=key[4], radius_mm=key[5],
                    comparisons=[])
        if "builtin" in vs:
            b = vs["builtin"][0]
            case.update(particles=b["particles"], particles_per_process=[b["particles_per_process_min"],
                                                                         b["particles_per_process_max"]],
                        builtin_median_step_s=b["step_wall_s"]["median"])
        labels = sorted(vs, key=lambda x: (x.endswith("+plainfma"), ORDER.index(x.split("+")[0])
                                           if x.split("+")[0] in ORDER else 99))
        for v in labels:
            base = BASELINE.get(v.split("+")[0])
            if base is None or base not in vs:
                continue
            sv, dv = vs[v]
            sb, db_ = vs[base]
            A = np.load(os.path.join(db_, "state.npz"))
            B = np.load(os.path.join(dv, "state.npz"))
            same_ids = bool(np.array_equal(A["gid"], B["gid"]))
            mb, mv = sb["step_wall_s"]["median"], sv["step_wall_s"]["median"]
            case["comparisons"].append(dict(
                variant=v, baseline=base, same_particles=same_ids,
                same_decomposition=sv["decomposition_fingerprint"] == sb["decomposition_fingerprint"],
                same_dt_statistics=sv["dt_s"] == sb["dt_s"], time_reached_s=[sb["time_reached_s"], sv["time_reached_s"]],
                median_step_s=[mb, mv], added_cost_fraction=mv/mb - 1.0,
                python_override_s_per_step_max_rank=sv["python_override_s_per_step_max_rank"],
                python_override_fraction_of_builtin_step=sv["python_override_s_per_step_max_rank"]/mb,
                host_load_at_start=[sb.get("host_load_average_at_start"), sv.get("host_load_average_at_start")],
                fields={f: field_diff(A[f], B[f]) for f in COMPARED} if same_ids else None))
        out["cases"].append(case)
    path = a.summary or os.path.join(a.out, "summary.json")
    with open(path, "w") as f:
        json.dump(out, f, indent=2)
        f.write("\n")
    for g, m in out["methods"].items():
        print(f"methods {g}: largest difference of any method from its built-in {m['max_ulp_any_method']} ULP")
    for case in out["cases"]:
        print(f"{case['geometry']} dx {case['dx_mm']} mm, {case.get('particles')} particles, {case['processes']} "
              f"processes, built-in median step {case.get('builtin_median_step_s', float('nan')):.4g} s")
        for c in case["comparisons"]:
            fl = c["fields"] or {}
            worst = max((f["max_abs_over_field_max"] for f in fl.values()), default=float("nan"))
            same = all(f["identical"] for f in fl.values()) if fl else False
            print(f"  {c['variant']:>22} vs {c['baseline']:<15} cost {c['added_cost_fraction']:+7.1%} "
                  f"(python {c['python_override_fraction_of_builtin_step']:6.1%})  "
                  + ("bitwise identical" if same else f"largest difference / field max {worst:.3g}"))
    print(f"-> {path}")
    return 0


# ------------------------------------------------------------------------------------------------- method check
def methods(argv):
    """Each Python class against its built-in, method by method, on synthetic inputs that reach every branch."""
    p = argparse.ArgumentParser(prog="subclass_check.py methods")
    p.add_argument("--geometry", choices=("3d", "rz"), required=True)
    p.add_argument("--n", type=int, default=4000)
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--json", default=None, help="write the result here")
    try:
        a = p.parse_args(argv)
    except SystemExit as e:
        return 2 if e.code else 0
    if a.geometry == "3d":
        import Spheral3d as S
        ndim = 3
    else:
        import SpheralRZ as S
        ndim = 2
    classes = make_classes(S, ndim)
    units = S.MKS()
    rng = np.random.default_rng(a.seed)
    n = a.n
    res = {}

    keep = []

    def nodelist(eos, strength, label):
        box = (S.Vector(*([-1.0]*ndim)), S.Vector(*([1.0]*ndim)))
        nl = S.makeSolidNodeList(label, eos, strength, xmin=box[0], xmax=box[1])
        nl.numInternalNodes = n
        keep.extend([eos, strength, nl])
        return nl

    def scalar(nl, name, values):
        f = S.ScalarField(name, nl)
        put_internal(f, values)
        return f

    rho0 = bench.RHO0
    # Densities over the whole bounded range and beyond it (both clamps, the fuzzy etamin test), energies from cold to
    # well past the melt energy, plastic strain to 2, pressures of either sign.
    eta = np.concatenate([rng.uniform(0.1, 4.5, n - 6), [bench.ETAMIN, bench.ETAMAX, 1.0, 1.0 + 1e-16, 0.2, 4.0]])
    rho_v = eta*rho0
    emelt = bench.EMELT_PER_VOLUME/rho0
    eps_v = np.concatenate([rng.uniform(-0.1*emelt, 1.3*emelt, n - 3), [emelt, emelt*(1 + 1e-16), 0.0]])
    P_v = rng.uniform(-5e9, 5e9, n)
    ps_v = np.abs(rng.standard_normal(n))*np.where(rng.uniform(size=n) < 0.3, 0.0, 1.0)
    cf_v = rng.uniform(1e3, 8e3, n)
    T_v = rng.uniform(0.0, 2000.0, n)

    def random_damage():
        """Symmetric tensors with eigenvalues in [0, 1.2] (both clamps), a few exactly diagonal or zero."""
        m = ndim
        Q, _ = np.linalg.qr(rng.standard_normal((n, m, m)))
        lam = rng.uniform(0.0, 1.2, (n, m))
        lam[: n//10] = 1e-5
        T = np.einsum("nij,nj,nkj->nik", Q, lam, Q)
        T[n//10: n//10 + 20] = np.eye(m)*rng.uniform(0, 1, (20, 1, 1))
        T[-5:] = 0.0
        return T

    def set_damage(nl, T):
        D = nl.damage()
        for i in range(n):
            t = T[i]
            D[i] = S.SymTensor(*t.ravel())
        return D

    # -- equations of state
    for kind, builder in (("murnaghan", murnaghan), ("linpoly", linpoly)):
        eb = builder(S, classes, False, units)
        ep = builder(S, classes, True, units)
        nl = nodelist(eb, S.NullStrength(), "eos_" + kind)
        rho = scalar(nl, "rho", rho_v)
        eps = scalar(nl, "eps", eps_v)
        T = scalar(nl, "T", T_v)
        out = {}
        for m, args in (("setPressure", (rho, eps)), ("setTemperature", (rho, eps)),
                        ("setSpecificThermalEnergy", (rho, T)), ("setSpecificHeat", (rho, T)),
                        ("setSoundSpeed", (rho, eps)), ("setGammaField", (rho, eps)), ("setBulkModulus", (rho, eps)),
                        ("setEntropy", (rho, eps))):
            fb, fp = S.ScalarField("b", nl), S.ScalarField("p", nl)
            getattr(eb, m)(fb, *args)
            getattr(ep, m)(fp, *args)
            out[m] = field_diff(get_internal(fb), get_internal(fp))
        fbs = [S.ScalarField(f"b{i}", nl) for i in range(3)]
        fps = [S.ScalarField(f"p{i}", nl) for i in range(3)]
        eb.setPressureAndDerivs(*fbs, rho, eps)
        ep.setPressureAndDerivs(*fps, rho, eps)
        for i, lbl in enumerate(("P", "dPdu", "dPdrho")):
            out["setPressureAndDerivs." + lbl] = field_diff(get_internal(fbs[i]), get_internal(fps[i]))
        res["eos_" + kind] = out

    # -- strength, with the built-in Murnaghan EOS (as in the 'strength' variant)
    eos = murnaghan(S, classes, False, units)
    sb = steinberg_guinan(S, classes, False, eos)
    sp = steinberg_guinan(S, classes, True, eos)
    nl = nodelist(eos, sb, "strength")
    rho = scalar(nl, "rho", rho_v)
    eps = scalar(nl, "eps", eps_v)
    P = scalar(nl, "P", P_v)
    ps = scalar(nl, "ps", ps_v)
    psr = scalar(nl, "psr", np.zeros(n))
    cf = scalar(nl, "cf", cf_v)
    D = set_damage(nl, random_damage())
    out = {}
    for m, args in (("shearModulus", (rho, eps, P, D)), ("yieldStrength", (rho, eps, P, ps, psr, D)),
                    ("soundSpeed", (rho, eps, P, cf, D)), ("meltSpecificEnergy", (rho, eps)),
                    ("coldSpecificEnergy", (rho, eps))):
        fb, fp = S.ScalarField("b", nl), S.ScalarField("p", nl)
        getattr(sb, m)(fb, *args)
        getattr(sp, m)(fp, *args)
        out[m] = field_diff(get_internal(fb), get_internal(fp))
    res["strength_steinberg_guinan"] = out

    # -- damage rate: the delegate's computeScalarDDDt against the Python one, on a state with random H and c_l
    WT = S.TableKernel(S.WendlandC4Kernel(), 1000)
    nl2 = nodelist(eos, sb, "damage")
    db = S.DataBase()
    db.appendNodeList(nl2)
    pdm = S.ProbabilisticDamageModel(materialName="aluminum", units=units, nodeList=nl2, kernel=WT, seed=bench.SEED)
    pyd = classes["PyDelegatingDamage"](pdm)
    Hf = nl2.Hfield()
    Ht = random_damage()*rng.uniform(10.0, 2000.0, (n, 1, 1)) + np.eye(ndim)*5.0
    for i in range(n):
        Hf[i] = S.SymTensor(*Ht[i].ravel())
    state = S.State()
    state.enroll(Hf)
    cl = S.ScalarField(S.SolidFieldNames.longitudinalSoundSpeed, nl2)
    put_internal(cl, rng.uniform(1e3, 8e3, n))
    state.enroll(cl)
    out = {}
    for label, excl in (("computeScalarDDDt", []), ("computeScalarDDDt.excludeNodes", [0, 3, n - 1])):
        pdm.excludeNodes = S.vector_of_int(excl)
        fb, fp = S.ScalarField("b", nl2), S.ScalarField("p", nl2)
        pdm.computeScalarDDDt(db, state, 0.0, 1e-7, fb)
        pyd.computeScalarDDDt(db, state, 0.0, 1e-7, fp)
        out[label] = field_diff(get_internal(fb), get_internal(fp))
    res["damage_rate"] = out

    worst = max(max(d["max_ulp_same_sign"], 1 if d["components_sign_or_zero_changed"] else 0)
                for grp in res.values() for d in grp.values())
    summary = dict(geometry=a.geometry, n=n, seed=a.seed, max_ulp_any_method=worst, methods=res)
    if a.json:
        os.makedirs(os.path.dirname(os.path.abspath(a.json)), exist_ok=True)
        with open(a.json, "w") as f:
            json.dump(summary, f, indent=2)
            f.write("\n")
    for grp, d in res.items():
        for m, r in d.items():
            print(f"METHOD {a.geometry} {grp}.{m}: max ulp {r['max_ulp_same_sign']}, {r['points_differing']} of "
                  f"{r['points']} differ")
    print(f"METHODS {a.geometry} max ulp {worst}")
    return 0


def main(argv):
    if not argv or argv[0] not in ("run", "compare", "methods"):
        print("usage: subclass_check.py {run,compare,methods} ... (see the module docstring)")
        return 2
    return {"run": run, "compare": compare, "methods": methods}[argv[0]](argv[1:])


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))

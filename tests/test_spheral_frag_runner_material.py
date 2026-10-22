"""The material without damage (Spheral M2 plan, Task 3; spec §8.2-8.3): the flow-stress law in the core
(spheral_frag.material: chen_plastic, chen_semi_solid, flow_stress, thermal_eos, thermal_shear_modulus) and its
Spheral classes (spheral_frag.runner.eos.PyThermalEOS, spheral_frag.runner.strength.PyThermalStrength), on the
committed AA7075_scheil table fixture (identical, column by column, to the MVP flight's prepared table).

Core (drama_env):
  - Chen et al. (2016)'s own peak stresses at 600 C, 8.2 MPa at 1e-3 1/s and 22.2 MPa at 1 1/s (Step 4 plan Task 5),
    are his semi-solid branch at his measured liquid fraction there, f_l = 0.191 (both quotes give f_l = 0.1910 to
    four figures; his data reach f_l = 0.19, Step 4 plan §4). Not the plastic branch (19.4 and 49.5 MPa) and not this
    table's Scheil f_l at 600 C (0.233: 5.8 and 16.1 MPa).
  - The plastic branch extended to 300 K at 1e-3 1/s gives 449.7 MPa (plan fact 10: 450 MPa).
  - Continuity at 623 K (the same branch on both sides: the provisional placeholder is Chen's branch extended);
    monotone in rate everywhere and in T within each branch; the jump at the solidus (748 K, plastic -> semi-solid at
    f_l -> 0+) measured: +11.9 % at 1e-3 1/s, +23.1 % at the 1e-6 1/s floor (the semi-solid fit's f_l = 0 limit is
    stronger than the plastic fit; Chen fitted them separately).
  - The equation of state: P = 0 at the free density, dP/drho and dP/du against finite differences, scale invariance.

Under Spheral (marker `spheral`, tests/fixtures/spheral_frag/runner_material_check.py):
  - every method of both classes against the core on 4,000 synthetic states over solid, hot solid, mush, slurry and
    liquid at s = 1 and s = 1e7 (and an EOS built at 1 whose scale is then set to 1e7), 3D and RZ, to round-off
    (measured 2026-10-09: every EOS method and the shear modulus bitwise; the flight law's yield 2.2e-17 of its
    largest value, numpy 1.26's and 2.5's exp/arcsinh);
  - the strength provides the longitudinal sound speed sqrt(c_f^2 + 4G/(3 rho)) (controller review: Spheral's step
    reads only the sound-speed field, so the bulk speed alone would raise the effective Courant number 1.23 times):
    Spheral's step per spacing with the Python classes is the step of bench.py's Murnaghan + SteinbergGuinan (M0's
    measured 4.54e-5 s/m in 3D, 4.43e-5 in RZ) times the ratio of the longitudinal speeds (test docstring);
  - P, G and Y at (s rho, s) equal those at (rho, 1) to round-off (measured: G and Y bitwise, P 1.1e-14 of the
    largest |P|, i.e. 2e-16 of K, from rounding s rho / (s rho_free)); c^2 equals the finite-difference dP/drho to
    1e-8 (measured 1.8e-10, the difference's own round-off);
  - a uniform 5 mm lattice sphere at rest at 300 K and 700 K, built at its free density, stays stress-free for 100
    SolidSPH steps at s = 1 and 1e7: |P| and |S| below 1e-6 of K (measured: 3D exactly zero, velocities and energies
    unchanged bitwise; RZ at most 1.2e-15 of K, |v| 6e-12 m/s, |u - u0| 5e-23 J/kg).

Cost (measured 2026-10-09, before the strength provided the sound speed, with the fixture's `bench` mode on the body
stage's 2.2 mm 3D sphere, 49,173 particles, bench.py's expansion at 10 1/s, 10 warm-up and 40 timed steps, outputs
in spheral_output/m2/task3_bench/): the two Python classes against Murnaghan + ConstantStrength, medians per step,
on 18 processes 110.5, 115.0, 113.4 ms against 107.4, 104.1, 107.3 ms (+6 %, range +3 to +10 %); on 8 processes
186.2, 196.8 ms against 174.0, 173.1 ms (+10 %, range +7 to +14 %); M0's Python EOS + strength + damage cost +19.7 %
in 3D. The Linux VM's load rose from 0.8 to 7.9 during the series (the runs themselves). With the sound speed
provided (same body, 18 processes, two runs each): 112.8, 112.5 ms against Murnaghan + ConstantStrength 107.9, 105.8
ms (+5 %) and bench.py's Murnaghan + SteinbergGuinan 104.2, 103.4 ms (+8.5 %)."""
import math
import os
import shutil
import subprocess
import tempfile

import numpy as np
import pytest

from conftest import spheral_launcher
from spheral_frag.m0 import bench
from helpers import FIXTURES, REPO_ROOT

from spheral_frag import material

TABLE = os.path.join(FIXTURES, "spheral_frag", "material", "material_AA7075_scheil.npz")
SCRIPT = "tests/fixtures/spheral_frag/runner_material_check.py"
SCALE = 1.0e7
N_STATES = 4000
R, DX = 0.05, 0.005
T600C = 873.15
F_L_CHEN_600C = 0.191          # Chen et al. (2016)'s measured liquid fraction at 600 C, as his own stresses imply
# Round-off thresholds (a priori): a few ulps of the quantity's largest magnitude; the core runs under drama_env's
# numpy 2.5 and the classes under Spheral's numpy 1.26 (exp, log, arcsinh may round differently by an ulp).
ROUND = 1e-13


@pytest.fixture(scope="module")
def table():
    return material.MaterialTable.load(TABLE)


def branches(table, T):
    f_l = table.liquid_fraction(T)
    return dict(solid=T < 623.0, hot_solid=(T >= 623.0) & (f_l <= 0.0), mush=(f_l > 0.0) & (f_l < 0.5),
                slurry=(f_l >= 0.5) & (T < table.T_liquidus), liquid=T >= table.T_liquidus)


def synthetic_states(table, n=N_STATES, seed=5):
    """Temperatures spread over the five branches (solid, hot solid, mush, slurry, liquid) with their edges, densities
    within 2 % of the free density (and exactly at it), plastic strain rates zero or 1e-9 to 10 1/s."""
    rng = np.random.default_rng(seed)
    edges = [250.0, 623.0, table.T_solidus, table.T_half, table.T_liquidus, 1500.0]
    special = np.array([250.0, 300.0, 622.999999, 623.0, 623.000001, table.T_solidus, table.T_solidus + 1e-9, 760.0,
                        table.T_half - 1e-9, table.T_half, table.T_half + 1e-9, 900.0, table.T_liquidus, 1000.0,
                        1499.0, 1500.0])
    k = (n - len(special))//5
    T = np.concatenate([rng.uniform(edges[i], edges[i + 1], k) for i in range(5)])
    T = np.concatenate([T, rng.uniform(250.0, 1500.0, n - len(special) - len(T)), special])
    rho = table.rho_free(T)*(1.0 + rng.uniform(-0.02, 0.02, n))
    rho[::10] = table.rho_free(T[::10])
    psr = np.where(rng.uniform(size=n) < 0.3, 0.0, 10.0**rng.uniform(-9.0, 1.0, n))
    return dict(T=T, u=table.enthalpy(T), rho=rho, psr=psr)


# ---------------------------------------------------------------------------------------------------------------
# Core: the flow-stress law
# ---------------------------------------------------------------------------------------------------------------

def test_chen_peak_stresses_at_his_conditions(table):
    s3 = material.chen_semi_solid(T600C, F_L_CHEN_600C, 1e-3)
    s0 = material.chen_semi_solid(T600C, F_L_CHEN_600C, 1.0)
    assert abs(s3 - 8.2e6) < 0.05e6 and abs(s0 - 22.2e6) < 0.05e6      # his values, quoted to 0.1 MPa
    # the liquid fraction each quote implies, by bisection on the monotone branch: both 0.1910 (his data reach 0.19)
    for rate, target in ((1e-3, 8.2e6), (1.0, 22.2e6)):
        lo, hi = 0.0, 0.49
        for _ in range(60):
            mid = 0.5*(lo + hi)
            lo, hi = (mid, hi) if material.chen_semi_solid(T600C, mid, rate) > target else (lo, mid)
        assert abs(lo - F_L_CHEN_600C) < 5e-4, (rate, lo)
    # not the plastic branch, and not this table's Scheil liquid fraction at 600 C
    assert material.chen_plastic(T600C, 1e-3) == pytest.approx(19.44e6, rel=1e-3)
    assert material.chen_plastic(T600C, 1.0) == pytest.approx(49.51e6, rel=1e-3)
    f_l = float(table.liquid_fraction(T600C))
    assert f_l == pytest.approx(0.2334, abs=1e-4)
    assert material.flow_stress(table, T600C, 1e-3) == pytest.approx(5.76e6, rel=1e-3)
    assert material.flow_stress(table, T600C, 1.0) == pytest.approx(16.06e6, rel=1e-3)


def test_plastic_branch_at_300K_and_its_provisional_label():
    assert material.chen_plastic(300.0, 1e-3) == pytest.approx(449.7e6, abs=0.1e6)     # plan fact 10: 450 MPa
    assert material.flow_stress_provisional("flight") == (material.FLOW_STRESS_PROVISIONAL,)
    assert "623 K" in material.FLOW_STRESS_PROVISIONAL and material.flow_stress_provisional("elastic") == ()
    assert material.flow_stress_provisional("constant:1e6") == ()


def test_continuity_at_623K_and_the_solidus_jump(table):
    for rate in (1e-6, 1e-3, 1.0):
        lo, hi = material.flow_stress(table, [623.0 - 1e-9, 623.0 + 1e-9], rate)
        assert abs(hi/lo - 1.0) < 1e-10                  # continuous: the same branch on both sides
        assert material.flow_stress(table, 623.0, rate) == material.chen_plastic(623.0, rate)
    Ts = table.T_solidus
    jumps = {}
    for rate in (1e-6, 1e-3):
        below, above = material.flow_stress(table, [Ts, Ts + 1e-9], rate)
        jumps[rate] = above/below - 1.0
        assert below == material.chen_plastic(Ts, rate)
    # measured 2026-10-09: +23.1 % at the floor, +11.9 % at 1e-3 1/s (the semi-solid fit's f_l -> 0 limit)
    assert jumps[1e-6] == pytest.approx(0.231, abs=1e-3) and jumps[1e-3] == pytest.approx(0.119, abs=1e-3)


def test_monotone_in_rate_and_temperature(table):
    rates = np.logspace(-6, 1, 50)
    for T in (250.0, 300.0, 500.0, 623.0, 700.0, 748.0, 760.0, 850.0, 890.0, 895.0):
        Y = material.flow_stress(table, np.full(rates.shape, T), rates)
        assert np.all(np.diff(Y) > 0.0), T
    for rate in (1e-6, 1e-3, 1.0):
        for lo, hi in ((250.0, table.T_solidus), (table.T_solidus + 1e-6, table.T_half - 1e-6)):
            T = np.linspace(lo, hi, 2000)
            Y = material.flow_stress(table, T, rate)
            assert np.all(np.diff(Y) < 0.0), (rate, lo)
        assert np.all(material.flow_stress(table, np.linspace(table.T_half, 1500.0, 100), rate) == 0.0)


def test_rate_floor_laws_and_overflow(table):
    T = np.array([300.0, 700.0, 800.0])
    assert np.array_equal(material.flow_stress(table, T, 0.0), material.flow_stress(table, T, material.RATE_MIN))
    assert np.array_equal(material.flow_stress(table, T, 0.0, rate_min=1e-7),
                          material.flow_stress(table, T, 1e-7, rate_min=1e-7))
    assert material.RATE_MIN == 1e-6 and material.RATE_MIN_BRACKET == (1e-7, 1e-5)
    assert np.all(material.flow_stress(table, T, 1.0, law="elastic") == np.inf)
    assert np.all(material.flow_stress(table, T, 1.0, law="constant:2.5e7") == 2.5e7)
    for bad in ("plastic", "constant:", "constant:-1", "constant:abc", "constant:inf"):
        with pytest.raises(ValueError, match="law"):
            material.flow_stress(table, T, 1.0, law=bad)
    assert material.chen_plastic(300.0, 0.0) == 0.0 and material.chen_semi_solid(800.0, 0.1, 0.0) == 0.0
    assert material.chen_semi_solid(800.0, [0.5, 0.7], 1.0).tolist() == [0.0, 0.0]
    with np.errstate(all="raise"):                       # logs: exp(Q/RT) would overflow below about 26 K
        Y = material.chen_plastic(np.array([1.0, 10.0, 50.0]), 1e-3)
    assert np.all(np.isfinite(Y)) and np.all(np.diff(Y) < 0.0)
    # the large-z form z + ln 2 equals arcsinh(exp(z)) to round-off where it takes over (z = 20)
    assert abs(math.asinh(math.exp(20.0)) - (20.0 + math.log(2.0))) <= 4e-15


# ---------------------------------------------------------------------------------------------------------------
# Core: the equation of state
# ---------------------------------------------------------------------------------------------------------------

def test_thermal_eos_core(table):
    st = synthetic_states(table)
    T, u, rho = st["T"], st["u"], st["rho"]
    e = material.thermal_eos(table, rho, u)
    assert np.max(np.abs(e["T"] - T)) < 1e-9
    K, rf = table.bulk_modulus(e["T"]), table.rho_free(e["T"])
    assert np.array_equal(e["P"], K*(rho/rf - 1.0))
    free = material.thermal_eos(table, rf, u)
    assert np.max(np.abs(free["P"])) <= 1e-15*K.max()      # stress-free at the free density
    assert np.allclose(e["cs"]**2, e["dPdrho"], rtol=1e-15, atol=0)
    assert np.array_equal(e["K"], rho*e["dPdrho"])
    # dP/du against central differences in u, away from the table's nodes (one-sided there)
    hT = np.concatenate([table.h_T, table.T_grid, table.fl_T])
    away = np.min(np.abs(T[:, None] - hT[None, :]), axis=1) > 0.02
    du = 1e-5*table.dh_dT(T)                             # dT about 1e-5 K
    fd = (material.thermal_eos(table, rho, u + du)["P"] - material.thermal_eos(table, rho, u - du)["P"])/(2*du)
    scale = np.abs(e["dPdu"]) + 1e-3*np.abs(e["dPdu"]).max()
    assert np.max(np.abs(fd - e["dPdu"])[away]/scale[away]) < 1e-4
    # dh/dT against the enthalpy
    dT = 1e-4
    fdh = (table.enthalpy(T + dT) - table.enthalpy(T - dT))/(2*dT)
    away_h = np.min(np.abs(T[:, None] - table.h_T[None, :]), axis=1) > 2*dT
    assert np.max(np.abs(fdh - table.dh_dT(T))[away_h]/table.dh_dT(T)[away_h]) < 1e-6
    # scale invariance
    es = material.thermal_eos(table, SCALE*rho, u, SCALE)
    for k in ("P", "K", "gamma", "dPdu"):
        assert np.max(np.abs(es[k] - e[k])) <= ROUND*np.abs(e[k]).max(), k
    assert np.max(np.abs(es["dPdrho"]*SCALE/e["dPdrho"] - 1.0)) < 1e-15
    G = material.thermal_shear_modulus(table, T)
    assert np.all(G[table.liquid_fraction(T) >= 0.5] == 0.0) and np.all(G[T < table.T_solidus] > 2.6e10)
    # the longitudinal speed the strength provides: sqrt(c_f^2 + 4G/(3 rho)), = c_f where G = 0, 1/sqrt(s) scaling
    c = material.thermal_sound_speed(table, e["T"], rho, e["cs"])
    assert np.allclose(c**2, e["cs"]**2 + 4.0/3.0*material.thermal_shear_modulus(table, e["T"])/rho, rtol=1e-15)
    liquid = table.liquid_fraction(e["T"]) >= 0.5
    assert np.array_equal(c[liquid], e["cs"][liquid]) and np.all(c[~liquid] > e["cs"][~liquid])
    cs_ = material.thermal_sound_speed(table, e["T"], SCALE*rho, es["cs"])
    assert np.max(np.abs(cs_*math.sqrt(SCALE)/c - 1.0)) < 1e-14
    c300 = material.thermal_sound_speed(table, 300.0, table.rho_free(300.0),
                                        material.thermal_eos(table, table.rho_free(300.0), 0.0)["cs"])
    assert c300 == pytest.approx(6146.9, abs=0.1)        # K 70.29 GPa, G 26.95 GPa, rho 2,811.6 kg/m^3


# ---------------------------------------------------------------------------------------------------------------
# Under Spheral
# ---------------------------------------------------------------------------------------------------------------

def launch(*args, nproc=1, timeout=1800):
    cmd = [spheral_launcher()] + (["-n", str(nproc)] if nproc > 1 else []) + [SCRIPT, *args]
    r = subprocess.run(cmd, cwd=REPO_ROOT, capture_output=True, text=True, timeout=timeout)
    assert r.returncode == 0, r.stdout[-3000:] + r.stderr[-3000:]
    return r


def rel(p):
    return os.path.relpath(p, REPO_ROOT)


@pytest.fixture(scope="module")
def work():
    base = os.path.join(REPO_ROOT, "spheral_output", "tests")      # inside the repository: the container sees it
    os.makedirs(base, exist_ok=True)
    d = tempfile.mkdtemp(prefix="runner_material_", dir=base)
    yield d
    shutil.rmtree(d, ignore_errors=True)


@pytest.fixture(scope="module")
def methods(work, table):
    st = synthetic_states(table)
    inp = os.path.join(work, "states.npz")
    np.savez(inp, **st)
    out = {}
    for form in ("3d", "rz"):
        p = os.path.join(work, f"methods_{form}.npz")
        launch("methods", "--form", form, "--table", rel(TABLE), "--input", rel(inp), "--scale", repr(SCALE),
               "--out", rel(p))
        out[form] = dict(np.load(p))
    return st, out


def close(got, want, what):
    want = np.asarray(want, dtype=float)
    big = np.abs(want[np.isfinite(want)]).max() if np.isfinite(want).any() else 0.0
    assert got.shape == want.shape, what
    assert np.array_equal(np.isfinite(got), np.isfinite(want)), what
    fin = np.isfinite(want)
    err = np.max(np.abs(got[fin] - want[fin])) if fin.any() else 0.0
    assert err <= ROUND*big, f"{what}: {err:.3e} > {ROUND:.0e} of {big:.3e}"


@pytest.mark.spheral
@pytest.mark.parametrize("form", ("3d", "rz"))
def test_every_method_against_the_core(methods, table, form):
    st, out = methods
    z = out[form]
    T, u, rho, psr = st["T"], st["u"], st["rho"], st["psr"]
    br = branches(table, T)
    assert all(np.count_nonzero(m) >= 700 for m in br.values()), {k: int(m.sum()) for k, m in br.items()}
    rf0 = float(table.rho_free(material.T_ZERO))
    assert z["refdens"].tolist() == [rf0, SCALE*rf0, SCALE*rf0] and z["valid"].all()
    for tag, s in (("s1", 1.0), ("s", SCALE), ("set", SCALE)):
        e = material.thermal_eos(table, s*rho, u, s)
        close(z[f"{tag}.setPressure"], e["P"], f"{tag} P")
        close(z[f"{tag}.setPressureAndDerivs.P"], e["P"], f"{tag} P (derivs)")
        close(z[f"{tag}.setPressureAndDerivs.dPdu"], e["dPdu"], f"{tag} dPdu")
        close(z[f"{tag}.setPressureAndDerivs.dPdrho"], e["dPdrho"], f"{tag} dPdrho")
        close(z[f"{tag}.setTemperature"], table.temperature(u), f"{tag} T")
        close(z[f"{tag}.setSpecificThermalEnergy"], table.enthalpy(T), f"{tag} u")
        close(z[f"{tag}.setSpecificHeat"], table.cp(T), f"{tag} cp")
        close(z[f"{tag}.setSoundSpeed"], e["cs"], f"{tag} cs")
        close(z[f"{tag}.setGammaField"], e["gamma"], f"{tag} gamma")
        close(z[f"{tag}.setBulkModulus"], e["K"], f"{tag} K")
        assert not z[f"{tag}.setEntropy"].any()
        Tu = table.temperature(u)
        for law in ("flight", "elastic", "constant:2e7"):
            close(z[f"{tag}.shearModulus.{law}"], material.thermal_shear_modulus(table, Tu), f"{tag} G {law}")
            want = np.minimum(material.flow_stress(table, Tu, psr, law), np.finfo(float).max)
            close(z[f"{tag}.yieldStrength.{law}"], want, f"{tag} Y {law}")
            assert z[f"{tag}.provides.{law}"].tolist() == [True, False]
            c_l = material.thermal_sound_speed(table, Tu, s*rho, e["cs"])
            close(z[f"{tag}.soundSpeed.{law}"], c_l, f"{tag} c_l {law}")
            close(z[f"{tag}.soundSpeed_aliased.{law}"], c_l, f"{tag} c_l aliased {law}")
        # c^2 against the finite-difference dP/drho (P is linear in rho: exact up to round-off of the difference)
        fd = (z[f"{tag}.P_plus"] - z[f"{tag}.P_minus"])/(2e-6*s*rho)
        assert np.max(np.abs(fd/z[f"{tag}.setSoundSpeed"]**2 - 1.0)) < 1e-8
    # the yield strength reaches every branch: Chen plastic, semi-solid, zero
    Y = z["s1.yieldStrength.flight"]
    assert np.all(Y[br["solid"] | br["hot_solid"]] > 1e7) and np.all(Y[br["slurry"] | br["liquid"]] == 0.0)
    assert np.all(Y[br["mush"]] >= 0.0) and np.count_nonzero(Y[br["mush"]] > 0.0) > 600


@pytest.mark.spheral
@pytest.mark.parametrize("form", ("3d", "rz"))
def test_scale_invariance(methods, form):
    _, out = methods
    z = out[form]
    for q in ("setPressure", "setBulkModulus", "setGammaField", "setPressureAndDerivs.dPdu", "setTemperature",
              "shearModulus.flight", "yieldStrength.flight", "yieldStrength.constant:2e7"):
        for tag in ("s", "set"):
            ref = z[f"s1.{q}"]
            assert np.max(np.abs(z[f"{tag}.{q}"] - ref)) <= ROUND*np.abs(ref).max(), (q, tag)
    for tag in ("s", "set"):
        assert np.max(np.abs(z[f"{tag}.setSoundSpeed"]*math.sqrt(SCALE)/z["s1.setSoundSpeed"] - 1.0)) < 1e-14
        c1 = z["s1.soundSpeed.flight"]
        assert np.max(np.abs(z[f"{tag}.soundSpeed.flight"]*math.sqrt(SCALE)/c1 - 1.0)) < 1e-14


def lattice_sphere(form):
    k = int(math.floor(R/DX*(1 + 1e-12)))
    if form == "3d":
        g = np.arange(-k, k + 1)*DX
        X = np.stack(np.meshgrid(g, g, g, indexing="ij"), -1).reshape(-1, 3)
    else:
        X = np.stack(np.meshgrid(np.arange(-k, k + 1)*DX, (np.arange(0, k + 1) + 0.5)*DX, indexing="ij"),
                     -1).reshape(-1, 2)
    X = X[np.einsum("ij,ij->i", X, X) <= (R*(1 + 1e-9))**2]
    V = np.full(len(X), DX**3) if form == "3d" else 2*math.pi*X[:, 1]*DX**2
    return dict(x_ref=X, V=V, pid=np.arange(len(X), dtype=np.int64) + 1, dx=np.array(DX))


@pytest.mark.spheral
@pytest.mark.parametrize("scale", (1.0, SCALE))
@pytest.mark.parametrize("T", (300.0, 700.0))
@pytest.mark.parametrize("form", ("3d", "rz"))
def test_uniform_body_at_rest_stays_stress_free(work, table, form, T, scale):
    inp = os.path.join(work, f"lattice_{form}.npz")
    if not os.path.exists(inp):
        np.savez(inp, **lattice_sphere(form))
    out = os.path.join(work, f"rest_{form}_{T:.0f}_{scale:.0e}.npz")
    launch("rest", "--form", form, "--table", rel(TABLE), "--input", rel(inp), "--scale", repr(scale),
           "--T", repr(T), "--steps", "100", "--out", rel(out), nproc=4 if form == "3d" else 1)
    z = np.load(out)
    K = float(z["K"])
    assert int(z["steps"]) == 100 and float(z["time"]) > 0.0
    assert float(z["rho0"]) == scale*float(table.rho_free(T))
    assert np.max(np.abs(z["P"])) <= 1e-6*K
    assert np.max(np.abs(z["S"])) <= 1e-6*K
    # nothing heats a body at rest: 3D bitwise; RZ by round-off of its hoop terms (measured up to 5e-23 J/kg at
    # u = 0 and s = 1); threshold 1e-9 J/kg, i.e. about 1e-12 K
    assert np.max(np.abs(z["u"] - float(z["u0"]))) <= 1e-9
    if form == "3d":
        assert np.all(z["u"] == float(z["u0"]))
    print(f"{form} T {T:.0f} K s {scale:.0e}: max|P|/K {np.abs(z['P']).max()/K:.2e}, max|S|/K "
          f"{np.abs(z['S']).max()/K:.2e}, max|v| {np.abs(z['v']).max():.2e} m/s, max|u - u0| "
          f"{np.abs(z['u'] - float(z['u0'])).max():.2e} J/kg, t {float(z['time']):.3e} s")


# M0's measured step (spec §2, data/spheral/m0_benchmark.csv): dt = 4.54e-5 s/m x dx in 3D, 4.43e-5 in RZ, with
# bench.py's Murnaghan + SteinbergGuinan (K 70.3 GPa, G 27.6 GPa, rho 2,813: c_l = 6,170 m/s).
M0_DT_PER_DX = {"3d": 4.54e-5, "rz": 4.43e-5}


@pytest.mark.spheral
@pytest.mark.parametrize("form", ("3d", "rz"))
def test_time_step_follows_the_longitudinal_speed(work, table, form):
    """Spheral's own step with the Python classes against bench.py's Murnaghan + SteinbergGuinan (M0's measured step)
    on the 5 mm sphere at rest at 300 K, 100 steps. The ratio of the steps is the ratio of the longitudinal speeds:
    the table's K 70.29 GPa, G 26.95 GPa at rho_free 2,811.6 kg/m^3 give c_l = 6,146.9 m/s; bench.py's 70.3 and 27.6
    GPa at that density 6,170.9 m/s (its Murnaghan's c_f is sqrt(K/rho0)), a ratio of 1.0039. Measured 2026-10-09:
    3D 4.5575e-5 against 4.5378e-5 s/m (1.0043; also on the 2.2 mm 3D body at 18 processes: 4.5575e-5 against
    4.5384e-5, dt 0.1003 against 0.0998 us), RZ 4.4471e-5 against 4.4276e-5 (1.0044); with the bulk speed alone
    (bench.py's Murnaghan + ConstantStrength, no sound speed from the strength) 5.60e-5 and 5.47e-5, 1.234 times
    longer.
    Thresholds: the ratio within 1e-3 of the moduli's (the 4e-4 left is SteinbergGuinan's pressure term and the
    kernel's sampling of a non-uniform h; measured); the step within 1 % of M0's measured value scaled by the ratio."""
    inp = os.path.join(work, f"lattice_{form}.npz")
    if not os.path.exists(inp):
        np.savez(inp, **lattice_sphere(form))
    dt = {}
    for mat in ("thermal", "sg"):
        out = os.path.join(work, f"dt_{form}_{mat}.npz")
        launch("rest", "--form", form, "--table", rel(TABLE), "--input", rel(inp), "--material", mat, "--T", "300",
               "--steps", "100", "--out", rel(out), nproc=4)
        z = np.load(out)
        dt[mat] = float(np.median(z["dts"][20:]))/float(z["dx"])
    T = material.T_ZERO
    rho = float(table.rho_free(T))
    c_ours = math.sqrt((float(table.bulk_modulus(T)) + 4.0/3.0*float(table.shear_modulus(T)))/rho)
    c_m0 = math.sqrt(bench.MURNAGHAN["K"]/bench.RHO0 + 4.0/3.0*bench.SG["G0"]/rho)
    ratio = c_m0/c_ours
    assert ratio == pytest.approx(1.0039, abs=1e-4)
    assert dt["thermal"]/dt["sg"] == pytest.approx(ratio, abs=1e-3)
    assert dt["thermal"] == pytest.approx(M0_DT_PER_DX[form]*ratio, rel=1e-2)
    assert dt["sg"] == pytest.approx(M0_DT_PER_DX[form], rel=1e-2)

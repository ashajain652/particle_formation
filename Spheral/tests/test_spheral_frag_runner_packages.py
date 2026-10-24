"""The replay's physics packages (Spheral M2 plan, Task 4; spec §7, §9.1; D1, D2, D5): ImposedTemperature, Damping,
SurfaceLoads and Deceleration (spheral_frag/runner/packages/), under Spheral through
tests/fixtures/spheral_frag/runner_packages_check.py, on lattice spheres of radius 50 mm (the body stage's lattice:
3D with a node at the centre, RZ rows at r = (j + 1/2) dx) with the committed AA7075_scheil table.

Measured 2026-10-09 (the numbers the thresholds come from):

  - ImposedTemperature, 5 mm sphere, three synthetic frames about 20 steps apart: the imposed field equals the
    table's linear interpolation at every step to round-off, the frame columns bitwise where Spheral's time lands on
    the frame time, u = h(T) and T(u) = T to round-off. What the equation of state sees at the step's midpoint is
    the temperature at the step's start plus Spheral's own half-step heating: max |T(u_mid) - T_imp(t_mid)| 1.2508
    K (3D) and 1.2167 K (RZ) against the largest imposed half-step change 1.2475 and 1.2131 K (the 3e-3 K left is
    the work of the thermal expansion these large synthetic jumps drive). On the real flight (RZ 1.1 mm replay body,
    interval k 20 -> 21, s = 9.43e6 for 10 crossings, 400 steps, dt 1.48e-4 s) the lag is 1.44e-3 K, exactly the
    imposed half-step change (the interval's largest change 9.7 K); temperature.py explains why no package hook can
    do better without recomputing the pressure policies.
  - Damping, 5 mm sphere started with v = x (1/s), 600 steps (7 periods): the fitted decay rate of the breathing
    mode, less the undamped run's (the artificial viscosity's: 2.85e3 1/s in 3D, -2.6e2 in RZ), divided by
    omega_1 (1.839e5 1/s in 3D, D = 105 mm; 2.033e5 in RZ, D = 95 mm) gives zeta = 0.0506 (set 0.05) and 0.2046
    (0.2) in 3D, 0.0503 and 0.2027 in RZ: within 2.3 %. The breathing frequency 3.36e5 1/s in 3D (Love's 2.745 c_l
    / R = 3.375e5). 2 zeta omega_1 dt at most 0.018 (zeta 0.2, 5 mm).
  - Surface detection at THRESHOLD 0.5 and normals at 2 h, 2.2 mm: 3D 6,194 surface particles, none deeper than
    1.07 dx, all within 0.5 dx marked, normals within 2.39 deg of radial; RZ 63 rings, depth at most 0.86 dx,
    normals within 4.82 deg. Area shares sum to 4 pi R^2 to 2.2e-16. A uniform 10 kPa: net force 1.6e-15 (3D) and
    2.0e-17 (RZ) of p pi R^2, momentum after 20 steps 6.6e-16 of p pi R^2 t (3D; RZ axial 1e-20 kg m/s).
  - Frame 182's tables (lee-extended at the defaults) on the 2.2 mm sphere, s = 1e7: the particles' drag 30.676 N
    (3D) and 29.829 N (RZ) against the continuum integral of the same table over the sphere, 2 pi R^2
    int (p cos + tau sin) sin dtheta = 30.717 N: -0.13 % and -2.9 %. The plan's literal exposure rule
    ("most_upstream") gives 23.565 N (-23 %) and 26.109 N (-15 %): it puts the base pressure (3 % of p_stag =
    231 Pa, against the table's 19 Pa) on the whole lee and drops the shoulder's windward loads.
    Deceleration: 3D momentum after 50 steps 1.6e-15 of sum F t; the packages' own balance sum m (a - a_frame) at
    round-off in both forms. RZ's momentum is not at round-off: Spheral's SolidSPHRZ itself does not conserve axial
    momentum (a sheared 2.2 mm sphere without any package: 4.3e-7 after 1 step, 6.0e-4 after 10, 3.5e-2 after 100;
    3D 4e-15), so under the frame's loads the RZ body gains 1.4e-2 of sum F t in 50 steps; the damping bounds the
    drift it causes.

Package cost per step (fixture `cost` mode, the body stage's replay bodies k 0-47 at their frame-0 temperatures,
s for 10 crossings per 0.5 s, flight law, the frame-0/1 tables of the fixture; 30 timed steps after 5 warm-up;
spheral_output/m2/task4_cost/): the time inside each package's hooks per step, the largest over the ranks, is the
reliable measure; the step medians moved by +-10 ms between configurations and with the order they ran in (the
bare step measured 82.7, 122.3, 99.9 and 111.7 ms in four runs, against 109-112 ms with all four packages).
3D 2.2 mm, 49,173 particles, 18 processes: ImposedTemperature 0.21-0.25 ms, Damping 0.34-0.40, SurfaceLoads
0.17-0.32, Deceleration 0.19-0.21; all four 0.96 ms, 0.9 % of the ~110 ms step; SurfaceLoads' construction (surface
detection, exposure, shares, tables) 0.19 s once per interval. RZ 1.1 mm, 3,243 particles: on 4 processes 0.05,
0.07, 0.04-0.05, 0.04-0.05 ms, all four 0.19 ms of a 5.5 ms step (3.5 %), construction 0.015 s; on 1 process all
four 0.41 ms of 12.2 ms (3.4 %), construction 0.027 s.
"""
import math
import os
import shutil
import subprocess
import tempfile

import numpy as np
import pytest

from conftest import spheral_launcher
from helpers import FIXTURES, REPO_ROOT

from spheral_frag import loads as core_loads, material

TABLE = os.path.join(FIXTURES, "spheral_frag", "material", "material_AA7075_scheil.npz")
LOADS = os.path.join(FIXTURES, "spheral_frag", "loads", "loads_mvp.npz")
SCRIPT = "tests/fixtures/spheral_frag/runner_packages_check.py"
R = 0.05
P_UNIFORM = 1.0e4
FORMS = ("3d", "rz")
PROCS = {"3d": 4, "rz": 1}


def launch(*args, nproc=1, timeout=1800):
    cmd = [spheral_launcher()] + (["-n", str(nproc)] if nproc > 1 else []) + [SCRIPT, *args]
    r = subprocess.run(cmd, cwd=REPO_ROOT, capture_output=True, text=True, timeout=timeout)
    assert r.returncode == 0, r.stdout[-3000:] + r.stderr[-3000:]
    return r


def rel(p):
    return os.path.relpath(p, REPO_ROOT)


def lattice_sphere(form, dx):
    k = int(math.floor(R/dx*(1 + 1e-12)))
    if form == "3d":
        g = np.arange(-k, k + 1)*dx
        X = np.stack(np.meshgrid(g, g, g, indexing="ij"), -1).reshape(-1, 3)
    else:
        X = np.stack(np.meshgrid(np.arange(-k, k + 1)*dx, (np.arange(0, k + 1) + 0.5)*dx, indexing="ij"),
                     -1).reshape(-1, 2)
    X = X[np.einsum("ij,ij->i", X, X) <= (R*(1 + 1e-9))**2]
    V = np.full(len(X), dx**3) if form == "3d" else 2*math.pi*X[:, 1]*dx**2
    return dict(x_ref=X, V=V, pid=np.arange(len(X), dtype=np.int64) + 1, dx=np.array(dx))


@pytest.fixture(scope="module")
def table():
    return material.MaterialTable.load(TABLE)


@pytest.fixture(scope="module")
def work():
    base = os.path.join(REPO_ROOT, "spheral_output", "tests")      # inside the repository: the container sees it
    os.makedirs(base, exist_ok=True)
    d = tempfile.mkdtemp(prefix="runner_packages_", dir=base)
    yield d
    shutil.rmtree(d, ignore_errors=True)


def lattice_file(work, form, dx):
    p = os.path.join(work, f"lattice_{form}_{dx*1e3:.1f}mm.npz")
    if not os.path.exists(p):
        np.savez(p, **lattice_sphere(form, dx))
    return p


# ---------------------------------------------------------------------------------------------------------------
# ImposedTemperature
# ---------------------------------------------------------------------------------------------------------------

@pytest.mark.spheral
@pytest.mark.parametrize("form", FORMS)
def test_imposed_temperature_reproduces_the_table(work, table, form):
    d = lattice_sphere(form, 0.005)
    r = np.sqrt(np.einsum("ij,ij->i", d["x_ref"], d["x_ref"]))
    T0 = 300.0 + 100.0*(r/R)**2
    T1 = T0 + 30.0 + 20.0*np.sin(7.0*r/R)
    T2 = T1 + 50.0*(r/R)**3
    d["T_table"] = np.asfortranarray(np.stack([T0, T1, T2], 1))
    tau = 20*0.28*0.005/6147.0                     # about 20 of Spheral's steps per interval at 5 mm, s = 1
    d["times"] = 1.0 + np.array([0.0, tau, 2*tau])  # the replay's times are flight times, not zero-based
    inp = os.path.join(work, f"temperature_{form}.npz")
    np.savez(inp, **d)
    out = os.path.join(work, f"temperature_{form}_out.npz")
    launch("temperature", "--form", form, "--table", rel(TABLE), "--input", rel(inp), "--out", rel(out),
           nproc=PROCS[form])
    z = np.load(out)
    assert np.array_equal(z["pid"], d["pid"])
    times = d["times"]
    hit = 0
    for t, T, u in zip(z["t"], z["T"], z["u"]):
        j = min(max(int(np.searchsorted(times, t, side="right")) - 1, 0), len(times) - 2)
        w = min(max((t - times[j])/(times[j + 1] - times[j]), 0.0), 1.0)
        want = (1.0 - w)*d["T_table"][:, j] + w*d["T_table"][:, j + 1]
        assert np.max(np.abs(T - want)) <= 1e-10, t           # round-off of w (t ~ 1 s, steps ~ 1e-7 s)
        assert np.max(np.abs(u - table.enthalpy(T))) <= 1e-9*np.abs(u).max()
        assert np.max(np.abs(table.temperature(u) - T)) <= 1e-9
        for k, tk in enumerate(times):
            if t == tk:                                       # exactly at a frame time: the column bitwise
                assert np.array_equal(T, d["T_table"][:, k])
                hit += 1
    assert hit >= 1 and abs(z["t"][-1] - times[-1]) <= 1e-12
    n_eval, seen, half = z["monitor"]
    # the midpoint sees the step-start temperature plus Spheral's own half-step heating (temperature.py): within
    # 1 % of the largest imposed half-step change (measured 0.3 %)
    assert n_eval == z["steps"] + 1 and half > 0.0 and abs(seen - half) <= 0.01*half
    print(f"{form}: {int(z['steps'])} steps, max |T(u_mid) - T_imp(t_mid)| {seen:.4f} K, largest imposed half-step "
          f"change {half:.4f} K")


# ---------------------------------------------------------------------------------------------------------------
# Damping
# ---------------------------------------------------------------------------------------------------------------

def decay_rate(rows):
    from scipy.optimize import curve_fit
    t, q = rows[:, 0], rows[:, 1]

    def model(t, A, lam, w, ph, c):
        return A*np.exp(-lam*t)*np.cos(w*t + ph) + c

    p, _ = curve_fit(model, t, q, p0=[q[0], 1e3, 3.3e5, 0.0, 0.0], maxfev=20000)
    return p[1], p[2]


@pytest.mark.spheral
@pytest.mark.parametrize("form", FORMS)
def test_damping_decays_the_breathing_mode_at_zeta(work, form):
    inp = lattice_file(work, form, 0.005)
    lam, z = {}, None
    for zeta in (0.0, 0.05, 0.2):
        out = os.path.join(work, f"breathing_{form}_{zeta}.npz")
        launch("breathing", "--form", form, "--table", rel(TABLE), "--input", rel(inp), "--zeta", repr(zeta),
               "--steps", "600", "--out", rel(out), nproc=PROCS[form])
        z = np.load(out)
        lam[zeta], w = decay_rate(z["rows"])
        assert float(z["max_fraction"]) <= 0.05
    omega1 = float(z["omega1"])
    # omega_1 = pi c_l / D with D the axial span + dx (D5): 3D 105 mm, RZ 95 mm at 5 mm
    D = {"3d": 0.105, "rz": 0.095}[form]
    assert omega1 == pytest.approx(math.pi*float(z["c_l"])/D, rel=1e-12)
    # Love's lowest breathing mode of the uniform sphere, 2.745 c_l / R: measured -1.3 % (3D), -4.1 % (RZ) at 5 mm
    assert w/omega1 == pytest.approx(2.745*D/(math.pi*R), rel=0.06)
    for zeta in (0.05, 0.2):
        got = (lam[zeta] - lam[0.0])/omega1                  # the artificial viscosity's own decay taken off
        assert got == pytest.approx(zeta, rel=0.10), (zeta, got)
        print(f"{form}: zeta {zeta}: measured {got:.4f} ({got/zeta - 1:+.2%}); undamped decay {lam[0.0]:.3e} 1/s")


# ---------------------------------------------------------------------------------------------------------------
# SurfaceLoads and Deceleration
# ---------------------------------------------------------------------------------------------------------------

@pytest.fixture(scope="module")
def uniform(work):
    out = {}
    for form in FORMS:
        p = os.path.join(work, f"uniform_{form}.npz")
        launch("loads", "--form", form, "--table", rel(TABLE), "--input", rel(lattice_file(work, form, 0.0022)),
               "--tables", "uniform", "--pressure", repr(P_UNIFORM), "--steps", "20", "--out", rel(p),
               nproc={"3d": 8, "rz": 1}[form])
        out[form] = dict(np.load(p))
    return out


@pytest.mark.spheral
@pytest.mark.parametrize("form", FORMS)
def test_surface_detection_on_the_lattice_sphere(uniform, form):
    z = uniform[form]
    dx = 0.0022
    x, s, n = z["st.x"], z["st.surface"], z["st.normal"]
    r = np.sqrt(np.einsum("ij,ij->i", x, x))
    depth = (R - r)/dx
    assert np.all(s[depth < 0.5])                       # the outer shell is marked
    assert depth[s].max() <= 1.5                        # and nothing deeper than 1.5 dx (measured 1.07, 0.86)
    radial = x[s]/r[s][:, None]
    ang = np.degrees(np.arccos(np.clip(np.einsum("ij,ij->i", n[s], radial), -1.0, 1.0)))
    assert ang.max() <= 10.0                            # measured 2.39 deg (3D), 4.82 deg (RZ)
    assert int(z["n_surface"]) == int(s.sum()) == {"3d": 6194, "rz": 63}[form]
    assert int(z["n_exposed"]) == int(z["n_surface"])   # a convex body: every surface particle on the first body
    th = z["st.theta_deg"][s]
    assert np.allclose(th, np.degrees(np.arccos(np.clip(n[s][:, 0], -1, 1))))
    print(f"{form}: {int(s.sum())} surface particles, deepest {depth[s].max():.2f} dx, normals within "
          f"{ang.max():.2f} deg of radial")


@pytest.mark.spheral
@pytest.mark.parametrize("form", FORMS)
def test_area_shares_and_a_uniform_pressure(uniform, form):
    z = uniform[form]
    A = 4.0*math.pi*R**2
    assert abs(float(z["area_sum"])/A - 1.0) <= 1e-14
    assert abs(np.sum(z["st.area"])/A - 1.0) <= 1e-14
    assert np.all(z["st.area"][~z["st.surface"]] == 0.0)
    assert np.all(z["st.p_k"][z["st.surface"]] == P_UNIFORM)
    ref = P_UNIFORM*math.pi*R**2
    assert np.max(np.abs(z["force_k"])) <= 1e-3*ref       # measured 1.6e-15 (3D), 2e-17 (RZ)
    # no motion of the body as a whole: momentum along the axis (RZ) or every component (3D)
    p1 = z["p1"] if form == "3d" else z["p1"][:1]
    assert np.max(np.abs(p1)) <= 1e-12*ref*float(z["t"])
    # the pressure's forces: -p A n, summed per particle
    F = z["st.force_k"][z["st.surface"]]
    want = -P_UNIFORM*z["st.area"][z["st.surface"]][:, None]*z["st.normal"][z["st.surface"]]
    assert np.max(np.abs(F - want)) <= 1e-12*np.abs(want).max()


def continuum_drag(table):
    th = np.linspace(0.0, math.pi, 200001)
    p, tau = core_loads.evaluate(table, np.degrees(th))
    f = (p*np.cos(th) + tau*np.sin(th))*np.sin(th)
    return 2.0*math.pi*R**2*float(np.sum(0.5*(f[1:] + f[:-1])*np.diff(th)))


@pytest.fixture(scope="module")
def frame182(work):
    out = {}
    for form in FORMS:
        for exposure in ("first_body", "most_upstream"):
            p = os.path.join(work, f"frame182_{form}_{exposure}.npz")
            launch("loads", "--form", form, "--table", rel(TABLE), "--input",
                   rel(lattice_file(work, form, 0.0022)), "--tables", "fixture", "--fixture", rel(LOADS),
                   "--frames", "182,183", "--decel", "--exposure", exposure, "--scale", "1e7", "--steps", "50",
                   "--out", rel(p), nproc={"3d": 8, "rz": 1}[form])
            out[form, exposure] = dict(np.load(p))
    return out


def fixture_tables():
    z = np.load(LOADS)
    tabs = core_loads.unstack({k: z[k] for k in ("theta_deg", "edges_deg", "p", "tau", "area", "lee",
                                                  "theta_last_deg", "p_stag")})
    return {int(k): core_loads.with_lee(t) for k, t in zip(z["k"], tabs)}


@pytest.mark.spheral
@pytest.mark.parametrize("form", FORMS)
def test_frame_loads_drag_and_the_deceleration(frame182, form):
    tabs = fixture_tables()
    D_cont = [continuum_drag(tabs[182]), continuum_drag(tabs[183])]
    assert D_cont[0] == pytest.approx(30.717, abs=1e-3)
    z = frame182[form, "first_body"]
    lit = frame182[form, "most_upstream"]
    drag = -z["force_k"][:, 0]
    # the particles' drag against the continuum integral of the same table over the sphere (measured -0.13 % in 3D,
    # -2.9 % in RZ at 2.2 mm); the plan's literal exposure rule misses by -23 % and -15 %
    tol = {"3d": 0.005, "rz": 0.04}[form]
    for i in range(2):
        assert abs(drag[i]/D_cont[i] - 1.0) <= tol
        assert -lit["force_k"][i, 0]/D_cont[i] - 1.0 < -0.1
    # the deceleration: the packages' own sum m (a_loads - a_frame) vanishes to round-off, and sum m a_loads = sum F
    F = np.abs(z["F_t"][0])
    assert abs(z["bal_decel"][0]) <= 1e-12*F
    assert abs(z["sum_ma"][0] - z["F_t"][0]) <= 1e-12*F
    h = z["decel_history"]
    dim = z["force_k"].shape[1]
    assert len(h) == int(z["steps"]) == 50
    assert np.allclose(h[:, 1 + dim:], h[:, 1:1 + dim]/float(z["mass"]), rtol=1e-15, atol=0)
    # each row is the forces at the step's midpoint, linear between the frames
    w = h[:, 0]/0.5
    want = (1 - w)[:, None]*z["force_k"][0] + w[:, None]*z["force_k"][1]
    assert np.max(np.abs(h[:, 1:1 + dim] - want)) <= 1e-12*F
    pt = abs(drag[0])*float(z["t"])
    if form == "3d":
        assert np.max(np.abs(z["p1"])) <= 1e-12*pt          # measured 1.6e-15: at rest under the frame's loads
    else:
        # SolidSPHRZ does not conserve axial momentum itself (module docstring): measured 1.4e-2 after 50 steps
        assert abs(z["p1"][0]) <= 0.05*pt
    assert float(z["max_displacement_dx"]) < 0.01
    print(f"{form}: drag {drag[0]:.4f} N vs continuum {D_cont[0]:.4f} N ({drag[0]/D_cont[0] - 1:+.2%}); most "
          f"upstream {-lit['force_k'][0, 0]:.4f} N; momentum / (F t) {abs(z['p1'][0])/pt:.2e}")


@pytest.mark.spheral
def test_exposure_shadows_and_wakes(work):
    """A 20 mm lattice sphere at the origin with a 6 mm one 40 mm upstream (+x) and a 6 mm one 40 mm downstream:
    under "first_body" the upstream sphere is exposed whole, the main sphere is exposed except on the lines the
    upstream one covers (its nose and its lee there), and the downstream sphere is hidden on the lines the main
    sphere covers, i.e. entirely; under "most_upstream" only one surface particle per line is exposed."""
    dx = 0.002

    def ball(c, rad):
        k = int(math.ceil(rad/dx)) + 1
        g = np.arange(-k, k + 1)*dx
        X = np.stack(np.meshgrid(g, g, g, indexing="ij"), -1).reshape(-1, 3)
        X = X[np.einsum("ij,ij->i", X, X) <= rad**2*(1 + 1e-9)]
        return X + np.array([c, 0.0, 0.0])

    parts = [ball(0.0, 0.02), ball(0.04, 0.006), ball(-0.04, 0.006)]
    x = np.concatenate(parts)
    which = np.concatenate([np.full(len(p), i) for i, p in enumerate(parts)])
    inp = os.path.join(work, "exposure.npz")
    np.savez(inp, x=x, V=np.full(len(x), dx**3), dx=np.array(dx))
    out = os.path.join(work, "exposure_out.npz")
    launch("exposure", "--form", "3d", "--table", rel(TABLE), "--input", rel(inp), "--out", rel(out))
    z = np.load(out)
    s, e = z["surface"], z["exposed.first_body"]
    lines = np.rint(x[:, 1:]/dx).astype(int)
    shadow = {tuple(v) for v in lines[(which == 1)]}
    in_shadow = np.array([tuple(v) in shadow for v in lines])
    assert np.all(e[(which == 1) & s])
    assert np.all(e[(which == 0) & s & ~in_shadow]) and not np.any(e[(which == 0) & in_shadow])
    assert not np.any(e[which == 2])
    lit = z["exposed.most_upstream"]
    assert np.count_nonzero(lit) == len({tuple(v) for v in lines[s]})
    print(f"surface {int(s.sum())}, first_body exposed {int(e.sum())}, most_upstream {int(lit.sum())}")

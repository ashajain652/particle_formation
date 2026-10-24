"""analysis/spheral_m2_elastic_fe.py: the finite-element reference of check 6 (Spheral M2 plan, Task 8 Step 1), on
the coarse 10 mm ball only (about 2,300 tetrahedra, 11,700 dofs; a few seconds)."""
import importlib.util
import os

import numpy as np
import pytest

from helpers import REPO_ROOT

LOADS_FIXTURE = os.path.join(REPO_ROOT, "tests", "fixtures", "spheral_frag", "loads", "loads_mvp.npz")
E, NU, RHO = 71.7e9, 0.33, 2811.618          # aa7075_mechanical.json at 300 K; rho_free(300 K) of the MVP table
H = 0.010                                    # m


def _module():
    spec = importlib.util.spec_from_file_location("spheral_m2_elastic_fe",
                                                  os.path.join(REPO_ROOT, "analysis", "spheral_m2_elastic_fe.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def fe_mod():
    return _module()


@pytest.fixture(scope="module")
def prob(fe_mod):
    return fe_mod.Problem(fe_mod.ball_mesh(fe_mod.RADIUS, H), E, NU, RHO)


@pytest.fixture(scope="module")
def frame182(fe_mod):
    with np.load(LOADS_FIXTURE) as z:
        arrays = {k: z[k] for k in z.files}
    i = int(np.flatnonzero(arrays["k"] == 182)[0])
    return fe_mod.table_from_arrays(arrays, i, lee=True)


def _interior_points(n=200, seed=1):
    rng = np.random.default_rng(seed)
    d = rng.normal(size=(n, 3))
    d /= np.linalg.norm(d, axis=1)[:, None]
    return d*(0.049*rng.random(n)**(1.0/3.0))[:, None]


def test_traction_is_the_particle_packages_force_per_area(fe_mod, frame182):
    """The same definition as SurfaceLoads: traction = surface_forces(n, area 1, p(theta), tau(theta))."""
    from spheral_frag import loads as core_loads
    from spheral_frag.runner.packages import loads as pkg
    th = np.radians(np.linspace(0.0, 180.0, 37))
    ph = np.radians(np.linspace(0.0, 330.0, 12))
    T, P = np.meshgrid(th, ph)
    n = np.stack([np.cos(T), np.sin(T)*np.cos(P), np.sin(T)*np.sin(P)], axis=-1).reshape(-1, 3)
    p, tau = core_loads.evaluate(frame182, pkg.inclination_deg(n))
    expect = pkg.surface_forces(n, np.ones(len(n)), p, tau)
    got = fe_mod.traction_function(frame182)(n.T).T
    assert np.abs(got - expect).max() <= 1e-12*np.abs(expect).max()
    # the drag direction: at 45 deg from the nose both pressure and shear push toward the tail (-x)
    k = np.argmin(np.abs(np.degrees(np.arccos(n[:, 0])) - 45.0))
    assert got[k, 0] < 0.0


def test_uniform_pressure_gives_minus_p_identity(fe_mod, prob):
    """A sphere under uniform pressure p: sigma = -p I everywhere. On the closed faceted ball the discrete load has
    zero net force and moment, and P2 contains the constant-strain solution, so this holds to the solver tolerance
    (CG rtol 1e-12; measured 8.4e-13 p at worst over 200 interior points, threshold 1e-8 p)."""
    p0 = 1.0e4
    f, I = fe_mod.traction_load(prob, lambda n: -p0*np.asarray(n))
    assert np.abs(I["force"]).max() <= 1e-12*p0*I["area"]
    fb, a, alpha, before, after = fe_mod.balance(prob, f)
    assert np.abs(a).max() <= 1e-12*p0*I["area"]/prob.mass_total
    pts = _interior_points()
    for solver in (fe_mod.solve_weak, fe_mod.solve_pinned):
        u = solver(prob, fb)[0]
        s = fe_mod.stress_at(prob.basis, u, pts, prob.lam, prob.mu)
        expect = np.array([-p0, -p0, -p0, 0.0, 0.0, 0.0])
        assert np.abs(s - expect).max() <= 1e-8*p0, solver.__name__


def test_frame182_load_is_balanced_and_rigid_constraint_free(fe_mod, prob, frame182):
    """Frame 182's tractions plus -rho (a + alpha x X): net force and moment zero to round-off (measured 1.9e-15 of
    the drag and 1.3e-16 of drag R; threshold 1e-12); a_x = -drag / M; the stresses of the weak (Lagrange) and the
    pinned (3-2-1) rigid-mode constraints agree (measured 9e-13 of the largest von Mises; threshold 1e-9); the
    surface integral of the drag is within 1.5 % of the table's continuum integral on the exact sphere (measured
    -0.85 % at 10 mm, -0.20 % at 5 mm: O(h^2))."""
    f, I = fe_mod.traction_load(prob, fe_mod.traction_function(frame182))
    fb, a, alpha, before, after = fe_mod.balance(prob, f)
    D = I["drag"]
    assert D > 0.0
    assert np.linalg.norm(after[:3]) <= 1e-12*D
    assert np.linalg.norm(after[3:]) <= 1e-12*D*fe_mod.RADIUS
    assert a[0] == pytest.approx(-D/prob.mass_total, rel=1e-9)
    assert abs(D/fe_mod.continuum_drag(frame182) - 1.0) <= 1.5e-2
    uA = fe_mod.solve_weak(prob, fb)[0]
    uB, _, reaction = fe_mod.solve_pinned(prob, fb)
    assert np.linalg.norm(reaction[:3]) <= 1e-9*D                     # the supports carry nothing
    pts = _interior_points()
    sA = fe_mod.stress_at(prob.basis, uA, pts, prob.lam, prob.mu)
    sB = fe_mod.stress_at(prob.basis, uB, pts, prob.lam, prob.mu)
    assert fe_mod.von_mises(sA - sB).max() <= 1e-9*fe_mod.von_mises(sA).max()
    # the weak constraint holds: no mean displacement or rotation
    assert np.abs(prob.MR.T @ uA).max() <= 1e-9*np.abs(prob.MR.T @ uB).max()

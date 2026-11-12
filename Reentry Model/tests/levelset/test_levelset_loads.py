import math
from datetime import datetime

import numpy as np
import pytest

from reentry_model import atmosphere, trajectory as tj
from reentry_model.levelset import loads
from reentry_model.levelset.boxmesh import BoxBase, build
from reentry_model.levelset.cutcell import patches
from reentry_model.levelset.levelset import redistance
from reentry_model.levelset.markers import Markers
from reentry_model.shared import heating, nose, surface_flow

V_HAT = np.array([1.0, 0.0, 0.0])


def _sphere_markers(R=25e-3, centre=(0.0, 0.0, 0.0), n=400):
    rng = np.random.default_rng(3)
    u = rng.normal(size=(n, 3))
    u /= np.linalg.norm(u, axis=1)[:, None]
    c = np.array(centre)
    return Markers(c + R * u, u, np.full(n, 4 * math.pi * R * R / n), np.full((n, 2), 1 / R), np.arange(n),
                   np.arange(n), 1e-3)


def _state(rho=2.727e-5, T=203.6, V=7500.0, kn=0.0298, ma=26.2):
    fs = atmosphere.Freestream(rho, T, rho / 28.9644e-3 * 8.314462618 * T, 28.9644 * 1.66053906660e-27, np.zeros(3))
    return tj.AeroState(77.5e3, 0.0, 0.0, fs, np.array([V, 0.0, 0.0]), V, kn, ma, 0.9, np.array([-V * 1e-3, 0, 0]),
                        0.5 * rho * V * V)


def test_a_sphere_fits_its_own_radius():
    m = _sphere_markers(25e-3, (0.003, -0.002, 0.001))
    fit = loads.fit_nose(m, np.array([0.003, -0.002, 0.001]), V_HAT, "current", fallback=0.05)
    assert fit.fitted == pytest.approx(25e-3, abs=1e-6) and fit.used == fit.fitted
    assert fit.n_points > 100
    assert loads.fit_nose(m, np.zeros(3), V_HAT, "initial", fallback=0.05).used == 0.05


def test_a_flat_face_is_capped_at_1_67_transverse_radii():
    rng = np.random.default_rng(1)
    r, ang = 10e-3 * np.sqrt(rng.random(300)), 2 * math.pi * rng.random(300)
    x = np.column_stack([np.zeros(300), r * np.cos(ang), r * np.sin(ang)])
    n = np.tile(V_HAT, (300, 1))
    m = Markers(x, n, np.full(300, 1e-6), np.zeros((300, 2)), np.arange(300), np.arange(300), 1e-3)
    fit = loads.fit_nose(m, np.array([-0.02, 0.0, 0.0]), V_HAT, "current", fallback=0.05)
    assert fit.bound == pytest.approx(nose.NOSE_CAP_FACTOR * r.max(), rel=1e-12)
    assert fit.used == fit.bound and not math.isfinite(fit.fitted)


def test_a_size_feedback_other_than_the_named_ones_is_refused():
    m = _sphere_markers()
    with pytest.raises(ValueError, match="size_feedback"):
        loads.fit_nose(m, np.zeros(3), V_HAT, "bogus", fallback=0.05)


def test_loads_are_the_shared_laws():
    m = _sphere_markers()
    theta = loads.marker_angles(m, V_HAT)
    assert np.allclose(np.cos(theta), m.normal @ V_HAT)
    lee = loads.leeward(None, None, None, m, V_HAT, test="centroid", centroid=np.zeros(3))
    ang = nose.load_angles(theta, ~lee)
    fit = loads.fit_nose(m, np.zeros(3), V_HAT, "current", fallback=0.05)
    T_top = np.linspace(500.0, 900.0, m.n)
    st = _state()
    rho_l = 2400.0
    for h in (heating.PhysicsHeating(), heating.SesamEquivalentHeating()):
        flow = surface_flow.SurfaceFlow()
        gl = loads.GasLoads(h, flow, rho_l)
        q, res = gl.evaluate(st, m, T_top, ang, fit)
        ref_q = h.evaluate(st, ang, T_top, fit.used).q_conv
        ref_res = flow.evaluate(st, ang, fit.used, rho_l, T_top)
        assert np.array_equal(q, ref_q)
        assert np.array_equal(res.tau, ref_res.tau) and np.array_equal(res.dp_ds, ref_res.dp_ds)
        if isinstance(h, heating.PhysicsHeating):      # (the SESAM-equivalent mode is uniform by design)
            assert q[ang < 0.4].min() > 0.0 and (q[ang > 0.5 * math.pi] == 0.0).all()


def test_leeward_centroid_matches_shared_rule():
    m = _sphere_markers(centre=(0.0, 0.0, 0.0))
    c = np.array([0.004, 0.0, 0.0])
    lee = loads.leeward(None, None, None, m, V_HAT, test="centroid", centroid=c)
    assert np.array_equal(lee, ~nose.ahead_of(m.x, c, V_HAT))
    with pytest.raises(ValueError):
        loads.leeward(None, None, None, m, V_HAT, test="nonsense", centroid=c)


def test_on_a_sphere_shadow_and_facing_agree():
    C = np.array([1.5e-3] * 3)
    mesh = build(BoxBase.around((0, 0, 0), (3e-3, 3e-3, 3e-3), 0.125e-3))
    phi = redistance(mesh, np.linalg.norm(mesh.points - C, axis=1) - 1.0e-3)
    P = patches(mesh, phi)
    m = Markers.create(mesh, phi)
    m.assign(mesh, P)
    lee = loads.leeward(mesh, phi, P, m, V_HAT, test="shadow")
    theta = loads.marker_angles(m, V_HAT)
    facing = m.normal[:, 0] > 0.2
    away = m.normal[:, 0] < -0.2
    assert lee[away].all() and not lee[facing].any()
    assert (np.abs(theta[~lee]) <= 0.5 * math.pi + 0.2).all()

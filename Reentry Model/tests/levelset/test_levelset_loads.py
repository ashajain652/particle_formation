import math
from datetime import datetime

import numpy as np
import pytest

from reentry_model.levelset import loads
from reentry_model.levelset.boxmesh import BoxBase, build
from reentry_model.levelset.cutcell import patches
from reentry_model.levelset.levelset import redistance
from reentry_model.levelset.markers import Markers
from reentry_model.shared import heating, nose, surface_flow

from cases import gas_state

V_HAT = np.array([1.0, 0.0, 0.0])


def _sphere_markers(R=25e-3, centre=(0.0, 0.0, 0.0), n=400):
    rng = np.random.default_rng(3)
    u = rng.normal(size=(n, 3))
    u /= np.linalg.norm(u, axis=1)[:, None]
    c = np.array(centre)
    return Markers(c + R * u, u, np.full(n, 4 * math.pi * R * R / n), np.full((n, 2), 1 / R), np.arange(n),
                   np.arange(n), 1e-3)


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
    st = gas_state()
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


# -- the loads by the area each marker owns (Task 11's fix, item 2, and its round 1, C1; cited findings F10, F2) ------
BEAM = np.array([0.0, 0.0, 1.0])     # a beam arriving along -z: each patch's projection max(n . BEAM, 0)
COVERAGE_EXACT = 1e-9                # the heat the markers take against the flux over phi's patches, relative (C1,
                                     # cited finding F2)
TILT = np.array([0.3, -0.2, 0.15])   # rad: the tilted box's rotation vector


def _box_surface(half=(5e-3, 5e-3, 4e-3), edge=0.8e-3, fine=0.7e-3, tilt=None):
    """phi of a box with its edges rounded to `edge` [m] (the body's phi over a box shrunk by its stub), on a 2.8 mm base
    refined to `fine` near it, its centre off the grid, turned by the rotation vector `tilt` [rad] if given, with
    markers seeded one per 2 mm voxel and tessellated."""
    from scipy.spatial.transform import Rotation
    from cases import refined, seeded
    off = np.array([0.13e-3, 0.21e-3, 0.17e-3])
    half = np.asarray(half)
    rot = np.eye(3) if tilt is None else Rotation.from_rotvec(tilt).as_matrix()

    def level(x):
        q = np.abs((x - off) @ rot) - (half - edge)
        return np.linalg.norm(np.maximum(q, 0.0), axis=1) + np.minimum(q.max(axis=1), 0.0) - edge
    lo, hi = -half - 2.8e-3, half + 2.8e-3
    mesh = refined(BoxBase.around(lo, hi, 2.8e-3), level, fine)
    phi = redistance(mesh, level(mesh.points), band=4e-3)
    mk = seeded(mesh, phi, 2e-3, lo)
    mk.tessellate(mesh, phi)
    return mesh, phi, mk, off, half


def test_a_uniform_beam_on_a_flat_face_gives_every_face_marker_the_same_flux():
    """A beam along -z on a box's flat top face: every marker whose part of phi lies wholly on the face takes exactly the
    beam's flux by its pieces, and so does the shared heating's flux at the pieces of a face facing the flow."""
    mesh, phi, mk, off, half = _box_surface()
    P = patches(mesh, phi)
    on_face = loads.owned(mk, P) @ (P.normal[:, 2] < 1.0 - 1e-9).astype(float) == 0.0
    q = loads.patch_mean(mk, P, np.maximum(P.normal @ BEAM, 0.0))
    face = on_face & (mk.normal[:, 2] > 0.99)
    assert face.sum() >= 4
    np.testing.assert_allclose(q[face], 1.0, rtol=0.0, atol=1e-14)
    gl = loads.GasLoads(heating.PhysicsHeating(), surface_flow.SurfaceFlow(), 2400.0)
    fit = loads.NoseFit(float("nan"), 0.05, 0.05, 0, False)
    theta = np.arccos(np.clip(P.normal @ BEAM, -1.0, 1.0))
    qg, _ = loads.by_patches(gl, gas_state(), mk, P, np.full(mk.n, 700.0), theta, fit, BEAM)
    assert np.ptp(qg[face]) <= 1e-12 * qg[face].max() and qg[face].min() > 0.0


@pytest.mark.parametrize("tilt", [None, TILT])
def test_the_box_s_markers_take_its_projected_heat_exactly(tilt):
    """The q6p box's shape (8 x 10 x 10 mm, edges rounded to the stub's 0.8 mm) under a beam along -z, upright and
    turned so that no face is aligned with the beam or the grid: the markers, each its mean over the area of phi it owns
    on that area (Task 11's fix round 1, C1, cited finding F2: one partition of phi), take the beam over phi's patches,
    each patch's flux times its area, to COVERAGE_EXACT (the review's M3, cited finding F5: by the tangent-plane cells'
    areas the upright box's markers took 1.013 of it, within a band of 0.95 to 1.053 that could not see the error)."""
    mesh, phi, mk, off, half = _box_surface(half=(5e-3, 5e-3, 4e-3), tilt=tilt)
    P = patches(mesh, phi)
    proj = np.maximum(P.normal @ BEAM, 0.0)
    want = (P.area * proj).sum()
    own = (mk.area * np.maximum(mk.normal @ BEAM, 0.0)).sum() / want
    by_patch = (mk.area * loads.patch_mean(mk, P, proj)).sum() / want
    print("box ({}): the markers take {:.3f} of the projected heat by their own normals, {:.15f} by what they "
          "own".format("tilted" if tilt is not None else "upright", own, by_patch))
    assert abs(by_patch - 1.0) <= COVERAGE_EXACT
    assert mk.area.sum() == pytest.approx(P.area.sum(), rel=1e-12)


SPHERE_MEDIAN, SPHERE_TOTAL = 0.01, 0.005   # the patch loads against the own-normal loads on a smooth sphere: the median
                                            # marker's and the windward total's relative change at most these


def test_on_a_smooth_sphere_the_patch_loads_are_the_own_normal_loads_within_a_percent():
    """On the 50 mm sphere (radius 25 mm, 1.4 mm cells at its surface, markers 2 mm apart) at a gas state, the shared
    heating by the patches each marker owns differs from that at the marker's own normal by under SPHERE_MEDIAN for the
    median windward marker (within 60 degrees of the stagnation point) and by under SPHERE_TOTAL in their total.
    Measured once: median 0.8 per cent, 95th percentile 3.2 per cent, largest 5.0 per cent, total -0.04 per cent; a
    marker seeded on its voxel's largest patch sits off its cell's centre (its normal 0.76 degrees from its cell's mean
    at the median, 2.3 at most), and near 60 degrees the heating changes by tan(theta) per radian of angle."""
    from cases import refined, seeded
    R = 25e-3
    level = lambda x: np.linalg.norm(x, axis=1) - R
    mesh = refined(BoxBase.around((-R - 2.8e-3,) * 3, (R + 2.8e-3,) * 3, 2.8e-3), level, 1.4e-3)
    phi = redistance(mesh, level(mesh.points), band=4e-3)
    mk = seeded(mesh, phi, 2e-3, (-R - 2.8e-3,) * 3)
    mk.tessellate(mesh, phi)
    P = patches(mesh, phi)
    gl = loads.GasLoads(heating.PhysicsHeating(), surface_flow.SurfaceFlow(), 2400.0)
    fit = loads.NoseFit(R, R, R, mk.n, False)
    st, T = gas_state(), np.full(mk.n, 700.0)
    q_own, _ = gl.evaluate(st, mk, T, loads.marker_angles(mk, V_HAT), fit)
    q_p, flow = loads.by_patches(gl, st, mk, P, T, np.arccos(np.clip(P.normal @ V_HAT, -1.0, 1.0)), fit, V_HAT)
    wind = loads.marker_angles(mk, V_HAT) < math.radians(60.0)
    rel = np.abs(q_p[wind] / q_own[wind] - 1.0)
    tot = (mk.area[wind] * q_p[wind]).sum() / (mk.area[wind] * q_own[wind]).sum() - 1.0
    print("sphere: {} windward markers, the patch loads differ by at most {:.2e} (median {:.2e}), the total by "
          "{:.2e}".format(int(wind.sum()), rel.max(), np.median(rel), tot))
    assert np.median(rel) < SPHERE_MEDIAN and abs(tot) < SPHERE_TOTAL
    assert np.all(np.isfinite(flow.tau)) and len(flow.tau) == mk.n

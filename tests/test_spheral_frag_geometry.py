"""Point location, centre crossing and the thickness map (Spheral M1 plan, Task 5; check "Thickness map").

On the committed synthetic runs of tests/fixtures/spheral_frag/ (their README lists the answers), read back with
pyvista through `spheral_frag_synthetic.read_frame_pyvista` (the tests' independent reader; Task 4's
`frames.read_frame` is not used, so this file tests geometry alone). Thresholds:

- locator: P1 reproduces a linear field exactly, so 1e-12 (round-off on values of order 3); the surface slack is
  the module's SURFACE_TOL_M = 1e-12 m (plan Task 5, Step 1).
- thickness, sphere: the plan expected epsilon_sphere below 0.5 % (the facet-centroid offset h^2/8R with h = 8 mm).
  Measured on the fixture: 0.816 % (mean 0.369 %), because gmsh's facets are larger than h (centroids 0.25-0.70 % of R
  inside the sphere, the nearest facet plane at r_in = 49.546 mm). The check is therefore made against the faceted
  geometry, exactly: for an inscribed convex polyhedron the ray from a facet at plane distance d whose line passes b
  from the centre exits between the ball of radius r_in and the sphere, d + sqrt(r_in^2 - b^2) <= t <=
  d + sqrt(R^2 - b^2) (to 1e-12 m) -- and epsilon_sphere is bounded a priori by 1 - r_in/R + b^2/(4 R r_in) (0.907 %
  here), both computed from the fixture's own facets, not tuned.
- thickness, dumbbell: the neck's 567 lateral patches 2 r_neck within one surface cell (h = 1.5 mm, plan); the
  sphere patches within the same faceted bounds about their own sphere, and within one far cell (h_far = 5 mm) of 2R
  (measured worst 4.35 mm short, under three 13 mm facets at the left sphere's pole).
- slab: an exact box, so every patch's thickness is the box's extent along its normal to 1e-12 m.
- brute force and vtk agree to 1e-12 m (plan; measured 0.0, bitwise: both evaluate the same Moller-Trumbore test)."""
import json
import os
import time

import numpy as np
import pytest

import spheral_frag_synthetic as syn
from helpers import REPO_ROOT
from spheral_frag import geometry
from spheral_frag.geometry import TetLocator, centre_crossed, thickness_map, thickness_summary, thin_patches

pytest.importorskip("pyvista")

FIXTURE_DIR = os.path.join(REPO_ROOT, "tests", "fixtures", "spheral_frag")
RUNS = {"sphere_run": syn.SPHERE_RUN_NAME, "dumbbell_frame": syn.DUMBBELL_RUN_NAME, "slab_frame": syn.SLAB_RUN_NAME}
EXACT = 1e-12            # m, or the value scale of order 1-10 of the linear test field: round-off


def run_dir(fixture):
    return os.path.join(FIXTURE_DIR, fixture, RUNS[fixture])


@pytest.fixture(scope="module")
def frames():
    """Every committed frame, read once: {(fixture, k): Frame}."""
    out = {}
    for fixture, n in (("sphere_run", 3), ("dumbbell_frame", 1), ("slab_frame", 1)):
        for k in range(n):
            out[(fixture, k)] = syn.read_frame_pyvista(run_dir(fixture), k)[0]
    return out


@pytest.fixture(scope="module")
def thickness(frames):
    """Brute-force thickness and patch geometry of the four distinct surfaces, computed once."""
    out = {}
    for key in (("sphere_run", 0), ("sphere_run", 2), ("dumbbell_frame", 0), ("slab_frame", 0)):
        f = frames[key]
        area, normal, centroid, theta = syn.patch_geometry(f.points, f.faces)
        t0 = time.perf_counter()
        t = thickness_map(f.points, f.faces, normal, centroid, method="bruteforce")
        out[key] = dict(t=t, normal=normal, centroid=centroid, seconds=time.perf_counter() - t0)
    return out


def sphere_dead_elements():
    with open(os.path.join(FIXTURE_DIR, "sphere_run", syn.SPHERE_RUN_NAME + ".json")) as fh:
        return json.load(fh)["synthetic"]["dead_elements_frame2"]


def faceted_bounds(centroid, normal, centre, R, r_in):
    """(lower, upper) bounds of the inward chord from each facet of a convex polyhedron inscribed in the sphere
    (centre, R) and containing the ball (centre, r_in)."""
    c = centroid - centre
    d = np.einsum("ij,ij->i", c, normal)                          # facet plane distance from the centre
    b2 = np.einsum("ij,ij->i", c, c) - d ** 2                     # squared distance of the ray's line from the centre
    return d + np.sqrt(np.maximum(r_in ** 2 - b2, 0.0)), d + np.sqrt(R ** 2 - b2), d, b2


# ------------------------------------------------------------------------------------------------------- locator
def test_locator_finds_every_tet_centroid_in_its_own_tet(frames):
    for key in (("sphere_run", 0), ("sphere_run", 2), ("dumbbell_frame", 0), ("slab_frame", 0)):
        f = frames[key]
        loc = TetLocator.from_frame(f)
        idx, w = loc.locate(f.points[f.tets].mean(axis=1))
        np.testing.assert_array_equal(idx, np.arange(len(f.tets)))
        np.testing.assert_allclose(w, 0.25, rtol=0, atol=EXACT)


def test_locator_interpolates_a_linear_field_exactly_inside_and_nan_outside(frames):
    f = frames[("sphere_run", 0)]
    loc = TetLocator(f.points, f.tets)
    _, normal, centroid, _ = syn.patch_geometry(f.points, f.faces)
    r_in = np.einsum("ij,ij->i", centroid, normal).min()            # the faceted body contains this ball
    rng = np.random.default_rng(5)
    u = rng.normal(size=(20000, 3))
    x = (u / np.linalg.norm(u, axis=1, keepdims=True)) * (0.999 * r_in * rng.random(20000) ** (1 / 3))[:, None]
    linear = lambda p: 3.0 + 20.0 * p[:, 0] - 10.0 * p[:, 1] + 5.0 * p[:, 2]
    idx, w = loc.locate(x)
    assert (idx >= 0).all()
    np.testing.assert_allclose(w.sum(axis=1), 1.0, rtol=0, atol=EXACT)
    assert (w >= -1e-9).all()
    np.testing.assert_allclose(loc.interpolate(x, linear(f.points)), linear(x), rtol=0, atol=EXACT)
    two = loc.interpolate(x, np.column_stack([linear(f.points), 2.0 * f.points[:, 0]]))     # several fields at once
    np.testing.assert_allclose(two[:, 1], 2.0 * x[:, 0], rtol=0, atol=EXACT)
    # outside: beyond the circumscribed sphere, and far outside the grid
    y = (u[:2000] / np.linalg.norm(u[:2000], axis=1, keepdims=True)) * 0.0501
    y = np.vstack([y, [[1.0, 0.0, 0.0], [-5.0, 2.0, 1.0]]])
    idx, w = loc.locate(y)
    assert (idx == -1).all() and np.isnan(w).all()
    assert np.isnan(loc.interpolate(y, linear(f.points))).all()
    assert loc.contains(x).all() and not loc.contains(y).any()


def test_locator_ties_go_to_the_lowest_tet_and_the_surface_slack_is_1e_12_m(frames):
    f = frames[("slab_frame", 0)]                                  # an exact box: faces lie on x, y, z planes
    loc = TetLocator(f.points, f.tets)
    used = np.unique(f.tets)
    lowest = np.full(len(f.points), np.iinfo(np.int64).max)
    np.minimum.at(lowest, f.tets.ravel(), np.repeat(np.arange(len(f.tets)), 4))
    idx, _ = loc.locate(f.points[used])                           # a node lies in every tet that uses it
    np.testing.assert_array_equal(idx, lowest[used])
    # the centroid of an interior face lies in both tets sharing it
    from reentry_model.mesh import _FACE_OF_VERTEX
    tf = np.sort(f.tets[:, _FACE_OF_VERTEX].reshape(-1, 3), axis=1)
    owner = np.repeat(np.arange(len(f.tets)), 4)
    order = np.lexsort(tf.T[::-1])
    tf, owner = tf[order], owner[order]
    shared = np.nonzero((tf[1:] == tf[:-1]).all(axis=1))[0]
    idx, _ = loc.locate(f.points[tf[shared]].mean(axis=1))
    np.testing.assert_array_equal(idx, np.minimum(owner[shared], owner[shared + 1]))
    # surface slack on the top face z = H
    H = f.points[:, 2].max()
    p = np.array([[0.0071, 0.0013, H], [0.0071, 0.0013, H + 0.5e-12], [0.0071, 0.0013, H + 1e-11],
                  [0.0071, 0.0013, H + 1e-6]])
    idx, _ = loc.locate(p)
    assert (idx[:2] >= 0).all() and (idx[2:] == -1).all()


def test_locator_answer_does_not_depend_on_the_bucket_size(frames):
    f = frames[("dumbbell_frame", 0)]
    rng = np.random.default_rng(11)
    lo, hi = f.points.min(axis=0), f.points.max(axis=0)
    x = lo + (hi - lo) * rng.random((20000, 3))
    ref = TetLocator(f.points, f.tets).locate(x)
    assert 0 < (ref[0] >= 0).sum() < len(x)
    for cell in (1.5e-3, 6e-3, 1.5e-2):
        idx, w = TetLocator(f.points, f.tets, cell_size=cell).locate(x)
        np.testing.assert_array_equal(idx, ref[0])
        np.testing.assert_allclose(w, ref[1], rtol=0, atol=EXACT)


def test_centre_crossing_flags_exactly_the_dead_elements(frames):
    f1, f2 = frames[("sphere_run", 1)], frames[("sphere_run", 2)]
    dead = sphere_dead_elements()
    assert len(f1.tets) - len(f2.tets) == len(dead) == 10
    loc1, loc2 = TetLocator.from_frame(f1), TetLocator.from_frame(f2)
    crossed = centre_crossed(loc2, f1.points[f1.tets].mean(axis=1))          # one centre per frame-1 element
    np.testing.assert_array_equal(np.nonzero(crossed)[0], np.sort(dead))
    rng = np.random.default_rng(3)                                           # and random centres in frame 1's body
    x = rng.uniform(-0.05, 0.05, size=(200000, 3))
    owner = loc1.locate(x)[0]
    x = x[owner >= 0]
    in_dead = np.isin(owner[owner >= 0], dead)
    assert in_dead.sum() > 50
    np.testing.assert_array_equal(centre_crossed(loc2, x), in_dead)


# ------------------------------------------------------------------------------------------------------ thickness
def test_sphere_thickness_lies_within_the_faceted_bounds_and_records_epsilon_sphere(thickness):
    R = syn.SPHERE_R
    for k in (0, 2):
        th = thickness[("sphere_run", k)]
        t = th["t"]
        assert np.isfinite(t).all()                                      # closed, staircase included
        if k == 2:
            continue      # the staircase is not an inscribed convex polyhedron: closed and finite is the check
        d_all = np.einsum("ij,ij->i", th["centroid"], th["normal"])
        r_in = d_all.min()
        lower, upper, d, b2 = faceted_bounds(th["centroid"], th["normal"], np.zeros(3), R, r_in)
        assert (t >= lower - EXACT).all() and (t <= upper + EXACT).all()
        eps_sphere = np.abs(t / (2 * R) - 1).max()
        bound = 1 - r_in / R + b2.max() / (4 * R * r_in)
        # measured 2026-10-07 on the committed fixture: eps_sphere 0.816 %, bound 0.907 % (r_in 49.546 mm)
        assert eps_sphere <= bound
        assert eps_sphere == pytest.approx(0.008159, abs=5e-6)
        assert r_in == pytest.approx(0.049546, abs=1e-6)


def test_dumbbell_neck_is_2_r_neck_and_the_spheres_2R(frames, thickness):
    f = frames[("dumbbell_frame", 0)]
    d = syn.DUMBBELL
    t = thickness[("dumbbell_frame", 0)]["t"]
    normal, centroid = thickness[("dumbbell_frame", 0)]["normal"], thickness[("dumbbell_frame", 0)]["centroid"]
    neck = syn.dumbbell_neck_patches(f, d["R"], d["r_neck"], d["L_neck"])
    spheres = syn.dumbbell_sphere_patches(f, d["R"], d["r_neck"], d["L_neck"])
    assert neck.sum() == 567 and spheres.sum() == 1471
    assert np.isfinite(t).all()
    np.testing.assert_allclose(t[neck], 2 * d["r_neck"], rtol=0, atol=d["h"])        # within one surface cell
    assert t[neck].max() <= 2 * d["r_neck"] + EXACT                                   # inscribed: never thicker
    a = 0.5 * d["L_neck"] + d["R"]
    for side in (-1.0, 1.0):
        centre = np.array([side * a, 0.0, 0.0])
        on = (np.abs(np.linalg.norm(f.points[f.faces] - centre, axis=2) - d["R"]) < 1e-9).all(axis=1)
        sel = spheres & (np.sign(centroid[:, 0]) == side)
        assert (sel <= on).all()
        r_in = np.einsum("ij,ij->i", centroid[on] - centre, normal[on]).min()
        lower, upper, _, _ = faceted_bounds(centroid[sel], normal[sel], centre, d["R"], r_in)
        assert (t[sel] >= lower - EXACT).all() and (t[sel] <= upper + EXACT).all()
    np.testing.assert_allclose(t[spheres], 2 * d["R"], rtol=0, atol=d["h_far"])      # within one far cell
    # measured 2026-10-07: neck 7.879-7.992 mm; spheres 35.65-39.98 mm
    assert t[neck].min() == pytest.approx(7.879e-3, abs=1e-6) and t[spheres].min() == pytest.approx(35.649e-3, abs=1e-6)


def test_thin_patches_flags_exactly_the_neck_at_2_2_mm_and_nothing_at_1_5_mm(frames, thickness):
    f = frames[("dumbbell_frame", 0)]
    t = thickness[("dumbbell_frame", 0)]["t"]
    neck = syn.dumbbell_neck_patches(f)
    np.testing.assert_array_equal(thin_patches(t, 2.2e-3), neck)        # 8 mm < 4 x 2.2 = 8.8 mm
    assert not thin_patches(t, 1.5e-3).any()                           # 4 x 1.5 = 6 mm < 7.88 mm
    assert not thin_patches(np.array([np.nan]), 1.0).any()


def test_slab_thickness_is_the_box_extent_along_each_normal(frames, thickness):
    f = frames[("slab_frame", 0)]
    th = thickness[("slab_frame", 0)]
    extent = f.points.max(axis=0) - f.points.min(axis=0)               # 20, 4, 12 mm
    axis = np.argmax(np.abs(th["normal"]), axis=1)
    np.testing.assert_allclose(np.abs(th["normal"])[np.arange(len(axis)), axis], 1.0, rtol=0, atol=EXACT)
    np.testing.assert_allclose(th["t"], extent[axis], rtol=0, atol=EXACT)


def test_bruteforce_and_vtk_agree(frames, thickness):
    pytest.importorskip("vtkmodules")
    for key, th in thickness.items():
        f = frames[key]
        t_vtk = thickness_map(f.points, f.faces, th["normal"], th["centroid"], method="vtk")
        np.testing.assert_allclose(t_vtk, th["t"], rtol=0, atol=EXACT)


def test_thickness_does_not_depend_on_the_winding_of_the_triangles_hit(frames, thickness):
    f = frames[("sphere_run", 2)]
    th = thickness[("sphere_run", 2)]
    faces = f.faces.copy()
    flip = np.random.default_rng(1).random(len(faces)) < 0.5            # wind half the triangles inward
    faces[flip] = faces[flip][:, [0, 2, 1]]
    t = thickness_map(f.points, faces, th["normal"], th["centroid"], method="bruteforce")   # outward normals given
    np.testing.assert_allclose(t, th["t"], rtol=0, atol=EXACT)        # same hits; vertex order changes round-off


def test_escaped_rays_are_nan_and_counted(frames, thickness):
    f = frames[("sphere_run", 0)]
    th = thickness[("sphere_run", 0)]
    keep = th["centroid"][:, 0] < 0.04                                  # an open surface: the +x cap removed
    t = thickness_map(f.points, f.faces[keep], th["normal"][keep], th["centroid"][keep], method="bruteforce")
    lost = np.isnan(t)
    assert lost.any() and (th["centroid"][keep][lost, 0] < -0.03).all()  # only rays aimed at the hole escape
    s = thickness_summary(t)
    assert s["n_nan"] == int(lost.sum()) and s["min_m"] == pytest.approx(np.nanmin(t))
    # a normal passed inward-wrong sends its ray out of the body
    n = th["normal"].copy()
    n[0] *= -1.0
    t = thickness_map(f.points, f.faces, n, th["centroid"], method="bruteforce")
    assert np.isnan(t[0]) and np.isfinite(t[1:]).all()
    assert thickness_summary(th["t"]) == {"min_m": float(th["t"].min()), "n_nan": 0}


def test_auto_method_and_cost_per_frame(frames, thickness):
    f = frames[("dumbbell_frame", 0)]
    th = thickness[("dumbbell_frame", 0)]
    assert len(f.faces) <= geometry.BRUTEFORCE_MAX_FACES
    np.testing.assert_array_equal(thickness_map(f.points, f.faces, th["normal"], th["centroid"]), th["t"])
    with pytest.raises(ValueError):
        thickness_map(f.points, f.faces, th["normal"], th["centroid"], method="obb")
    # cost per frame on the fixtures (measured 2026-10-07: 0.05 s sphere, 0.12 s dumbbell, 0.07 s slab); a loose
    # ceiling so that a quadratic blow-up of the brute force is caught, not a slow machine
    assert all(v["seconds"] < 5.0 for v in thickness.values())

"""Point location, centre crossing, the thickness map (Spheral M1 plan, Task 5; check "Thickness map"), and the
slurry depth and the three zones (Task 6; check "Zones").

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
- brute force and vtk agree to 1e-12 m (plan; measured 0.0, bitwise: both evaluate the same Moller-Trumbore test).
- depths, slab: f_l is linear in depth below every exact top patch, so the crossing is exact up to round-off (the
  plan's "within ds/10" is the looser statement of the same; measured 5e-18 m), and 1e-12 m is used. Below a top
  patch whose cube spans two columns the P1 value is 0.5 + (sum lambda_i D_i - d)/G, so its crossing lies between
  the two columns' depths (a priori).
- depths, sphere: the plan's "within one element size" (h = 8 mm, the fixture's surface size); measured on the
  committed fixture: slurry depth 3.662-9.029 mm against 3.051 mm on frame 1 (error +0.61 to +5.98 mm) and on frame
  2's 1,184 patches on the sphere 4.211-9.635 mm against 5.236 mm (-1.02 to +4.40 mm); liquid depth 0 against 2.010 mm
  (P1 cannot hold a liquid layer thinner than the first element: the nodes below the surface are under the liquidus).
  A priori, independent of the fixture: where P1 f_l = 0.5 the containing tetrahedron has nodes on both sides of the
  f_l = 0.5 sphere, so |r(x) - r_half| <= its diameter (+ ds for the sampling)."""
import json
import os
import time

import numpy as np
import pytest

import spheral_frag_synthetic as syn
from helpers import REPO_ROOT
from spheral_frag import frames as fr
from spheral_frag import geometry
from spheral_frag.geometry import (ZONE_BULK, ZONE_INVALID, ZONE_NONE, ZONE_RUNOFF, ZONE_SKIN, TetLocator,
                                   bulk_slurry_area, centre_crossed, classify_zones, layer_depths, thickness_map,
                                   thickness_summary, thin_patches, zone_summary)

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


# ------------------------------------------------------------------------------------------- depths and zones
SLAB_COLS = range(len(syn.SLAB_SLURRY_DEPTH))


@pytest.fixture(scope="module")
def oriented(frames):
    """The frames as the API requires them: wound outward (frames.orient_outward) with their derived patch arrays,
    a locator and the thickness map. The fixtures are written outward, so nothing is flipped here; on the real
    export about 18 % of the triangles are (Task 5)."""
    out = {}
    for key in (("sphere_run", 0), ("sphere_run", 1), ("sphere_run", 2), ("slab_frame", 0)):
        f, flipped = fr.orient_outward(frames[key])
        assert len(flipped) == 0
        d = fr.derived_patch_arrays(f)
        th = thickness_map(f.points, f.faces, d["normal"], d["centroid"], method="bruteforce")
        out[key] = dict(frame=f, d=d, loc=TetLocator.from_frame(f), thickness=th)
    return out


@pytest.fixture(scope="module")
def slab_depths(oriented):
    o = oriented[("slab_frame", 0)]
    liq, slu, info = layer_depths(o["frame"], o["loc"], o["d"]["normal"], o["d"]["centroid"],
                                  max_depth=o["thickness"], return_info=True)
    top, col, exact = syn.slab_exact_patches(o["frame"])
    return dict(liquid=liq, slurry=slu, info=info, top=top, col=col, exact=exact)


def slab_zones(o, sd, **kw):
    p = o["frame"].patch
    return classify_zones(p["film_thickness"], sd["liquid"], sd["slurry"], p["delta_m"],
                          deep_thickness=p["deep_thickness"], **kw)


def per_column(values, sd):
    """The value on each column's exact top patches, checked to be one value per column (to round-off)."""
    out = []
    for c in SLAB_COLS:
        v = np.asarray(values)[sd["exact"] & (sd["col"] == c)]
        assert len(v) == (32 if c == 4 else 24)
        if v.dtype.kind == "f" and np.isnan(v).any():
            assert np.isnan(v).all()
            out.append(np.nan)
            continue
        np.testing.assert_allclose(v, v[0], rtol=0, atol=EXACT)
        out.append(v[0])
    return np.array(out)


def test_slab_slurry_depth_is_exact_below_every_exact_top_patch(oriented, slab_depths):
    sd = slab_depths
    np.testing.assert_allclose(per_column(sd["slurry"], sd), syn.SLAB_SLURRY_DEPTH, rtol=0, atol=EXACT)
    assert (sd["liquid"][sd["top"]] == 0.0).all()                    # the surface f_l stays below 1 everywhere
    # a top patch over two columns: its crossing lies between the two columns' depths (P1, a priori)
    D = np.asarray(syn.SLAB_SLURRY_DEPTH)
    o = oriented[("slab_frame", 0)]
    cube = np.floor(o["d"]["centroid"][:, 0] / syn.SLAB["h"]).astype(int)
    for i in np.flatnonzero(sd["top"] & ~sd["exact"]):
        a, b = D[min(cube[i] // 4, 4)], D[min((cube[i] + 1) // 4, 4)]
        assert min(a, b) - EXACT <= sd["slurry"][i] <= max(a, b) + EXACT
    # only patches whose surface value exceeds 0.5 are marched (cost scales with the molten area)
    f = o["frame"]
    assert sd["info"]["n_marched"] == np.count_nonzero(f.node["f_l"][f.faces].mean(axis=1) > 0.5)
    assert sd["info"]["n_wrong_way"] == 0
    assert (sd["slurry"][f.node["f_l"][f.faces].mean(axis=1) <= 0.5] == 0.0).all()


def test_slab_zones_at_film_limits_2_and_3_mm_with_both_zone_3_readings(oriented, slab_depths):
    o, sd = oriented[("slab_frame", 0)], slab_depths
    z = slab_zones(o, sd, film_limit=2e-3, dx=1.1e-3)
    np.testing.assert_allclose(per_column(z.layer, sd), [0.25e-3, 1.1e-3, 2.6e-3, 4.1e-3, 1.1e-3], rtol=0, atol=EXACT)
    np.testing.assert_array_equal(per_column(z.zone, sd), [ZONE_SKIN, ZONE_RUNOFF, ZONE_BULK, ZONE_BULK, ZONE_RUNOFF])
    np.testing.assert_allclose(per_column(z.skin, sd), [1e-4, 1e-4, 1e-4, 1e-4, np.nan], rtol=0, atol=EXACT)
    np.testing.assert_allclose(per_column(z.runoff, sd), [0.15e-3, 1.0e-3, 1.9e-3, 1.9e-3, 1.1e-3], rtol=0, atol=EXACT)
    np.testing.assert_allclose(per_column(z.bulk, sd), [0.0, 0.0, 0.6e-3, 2.1e-3, 0.0], rtol=0, atol=EXACT)
    # decision 6: the whole layer (zone 3, `layer`) or only the part below the film limit (`bulk`)
    np.testing.assert_array_equal(per_column(z.under_resolved, sd), [False, False, True, True, False])   # < 4.4 mm
    np.testing.assert_array_equal(per_column(z.under_resolved_below_limit, sd), [False, False, True, True, False])
    z05 = slab_zones(o, sd, film_limit=2e-3, dx=0.5e-3)                                                 # 2 mm
    np.testing.assert_array_equal(per_column(z05.under_resolved, sd), [False] * 5)
    np.testing.assert_array_equal(per_column(z05.under_resolved_below_limit, sd), [False, False, True, False, False])
    assert not slab_zones(o, sd).under_resolved.any()                                                   # dx None
    z3 = slab_zones(o, sd, film_limit=3e-3)
    np.testing.assert_array_equal(per_column(z3.zone, sd), [ZONE_SKIN, ZONE_RUNOFF, ZONE_RUNOFF, ZONE_BULK,
                                                             ZONE_RUNOFF])
    np.testing.assert_allclose(per_column(z3.bulk, sd), [0.0, 0.0, 0.0, 1.1e-3, 0.0], rtol=0, atol=EXACT)
    # the three parts add up to the layer (skin <= delta_m < film limit here); patches without film or melt: none
    for zz in (z, z3):
        total = np.where(np.isnan(zz.skin), 0.0, zz.skin) + zz.runoff + zz.bulk
        np.testing.assert_allclose(total[sd["top"]], zz.layer[sd["top"]], rtol=0, atol=EXACT)
    side = ~sd["top"] & (sd["slurry"] == 0.0)
    assert (z.zone[side] == ZONE_NONE).all() and side.any()
    # route-2 input: the zone-3 area, identical under both readings (zone 3 exactly where bulk > 0)
    area = o["d"]["area"]
    assert ((z.zone == ZONE_BULK) == (z.bulk > 0)).all()
    assert bulk_slurry_area(z, area) == pytest.approx(float(area[z.zone == ZONE_BULK].sum()), rel=1e-15)
    s = zone_summary(z, area)
    assert s["bulk_area_m2"] == bulk_slurry_area(z, area) and s["film_limit_m"] == 2e-3
    assert s["n_patches"]["bulk"] == int(np.count_nonzero(z.zone == ZONE_BULK))
    assert s["bulk_volume_below_limit_m3"] < s["bulk_volume_whole_layer_m3"]


def test_include_deep_adds_deep_thickness_to_the_layer_and_nothing_else(oriented, slab_depths):
    o, sd = oriented[("slab_frame", 0)], slab_depths
    z0 = slab_zones(o, sd)
    z1 = slab_zones(o, sd, include_deep=True)
    assert not z0.include_deep and z1.include_deep
    np.testing.assert_allclose(per_column(z1.layer, sd), [0.45e-3, 1.3e-3, 2.8e-3, 4.3e-3, 1.3e-3], rtol=0, atol=EXACT)
    np.testing.assert_array_equal(z1.layer, z0.layer + o["frame"].patch["deep_thickness"])
    np.testing.assert_array_equal(z1.skin, z0.skin)
    # the rule applied to the thicker layer: column 0 (0.45 mm > delta_m 0.3 mm) leaves the skin zone
    np.testing.assert_array_equal(per_column(z1.zone, sd), [ZONE_RUNOFF, ZONE_RUNOFF, ZONE_BULK, ZONE_BULK,
                                                             ZONE_RUNOFF])
    np.testing.assert_allclose(per_column(z1.bulk, sd), [0.0, 0.0, 0.8e-3, 2.3e-3, 0.0], rtol=0, atol=EXACT)
    with pytest.raises(ValueError):
        classify_zones(0.0, 0.0, 0.0, 0.0, include_deep=True)


def test_classify_zones_rule_on_scalars():
    # film, liquid, slurry, delta_m -> zone at the 2 mm limit
    cases = [((0.0, 0.0, 0.0, 3e-4), ZONE_NONE), ((1e-4, 0.0, 0.0, 3e-4), ZONE_SKIN),
             ((1e-4, 0.0, 1.5e-4, 3e-4), ZONE_SKIN), ((0.0, 0.0, 3e-4, 3e-4), ZONE_SKIN),
             ((1e-4, 0.0, 2e-4 + 1e-9, 3e-4), ZONE_RUNOFF),
             ((1e-4, 0.0, 0.0, np.nan), ZONE_RUNOFF), ((0.0, 0.0, 2e-3, 3e-4), ZONE_RUNOFF),
             ((1e-4, 0.0, 1.9e-3 + 1e-9, 3e-4), ZONE_BULK), ((1e-4, 0.0, 5e-3, np.nan), ZONE_BULK),
             ((1e-4, 0.0, np.nan, 3e-4), ZONE_INVALID)]
    for (film, liq, slu, dm), zone in cases:
        z = classify_zones(np.array([film]), np.array([liq]), np.array([slu]), np.array([dm]))
        assert z.zone[0] == zone, (film, liq, slu, dm)
    z = classify_zones(np.array([1e-4]), np.array([5e-4]), np.array([8e-4]), np.array([3e-4]))
    assert z.skin[0] == 3e-4 and z.runoff[0] == pytest.approx(6e-4, abs=1e-15)   # skin capped at delta_m
    assert z.zone.dtype == np.int8


def test_layer_depths_thresholds_f_l_source_cap_and_wrong_way_normals(oriented, slab_depths):
    o, sd = oriented[("slab_frame", 0)], slab_depths
    f, d, loc = o["frame"], o["d"], o["loc"]
    # the material's f_l as a callable of the nodal temperature (Task 3's MaterialTable.liquid_fraction) gives the
    # same answer as the frame's own nodal f_l when it is the same law
    liq, slu = layer_depths(f, loc, d["normal"], d["centroid"], max_depth=o["thickness"],
                            f_l=syn.synthetic_liquid_fraction)
    np.testing.assert_allclose(slu, sd["slurry"], rtol=0, atol=EXACT)
    # another threshold: f_l > 0.6 at D - 0.1 G = D - 2.4 mm (columns 2 and 3), none where the surface is below 0.6
    _, slu6 = layer_depths(f, loc, d["normal"], d["centroid"], slurry_fl=0.6)
    np.testing.assert_allclose(per_column(slu6, sd), [0.0, 0.0, 0.1e-3, 1.6e-3, 0.0], rtol=0, atol=EXACT)
    # a liquid threshold the slab reaches: f_l >= 0.6 is "liquid" with slurry_fl 0.55 -> D - 2.4 mm again
    liq6, slu55 = layer_depths(f, loc, d["normal"], d["centroid"], slurry_fl=0.55, liquid_fl=0.6)
    np.testing.assert_allclose(per_column(liq6, sd), [0.0, 0.0, 0.1e-3, 1.6e-3, 0.0], rtol=0, atol=EXACT)
    np.testing.assert_allclose(per_column(slu55, sd), [0.0, 0.0, 1.3e-3, 2.8e-3, 0.0], rtol=0, atol=EXACT)
    assert (liq6 <= slu55).all()
    with pytest.raises(ValueError):
        layer_depths(f, loc, d["normal"], d["centroid"], slurry_fl=0.6, liquid_fl=0.5)
    # the cap: both depths at most max_depth, and the rays that hit it are flagged
    liq_c, slu_c, info = layer_depths(f, loc, d["normal"], d["centroid"], max_depth=2e-3, return_info=True)
    np.testing.assert_allclose(per_column(slu_c, sd), np.minimum(syn.SLAB_SLURRY_DEPTH, 2e-3), rtol=0, atol=EXACT)
    assert info["left_body"][sd["exact"] & np.isin(sd["col"], [3])].all()
    assert not info["left_body"][sd["exact"] & np.isin(sd["col"], [0, 1, 4])].any()
    # without a cap the march still ends where the ray leaves the body (at the last sample inside, within ds)
    _, slu_n = layer_depths(f, loc, d["normal"], d["centroid"])
    np.testing.assert_allclose(slu_n[sd["top"]], sd["slurry"][sd["top"]], rtol=0, atol=EXACT)
    capped = np.isclose(sd["slurry"], o["thickness"], rtol=0, atol=EXACT) & (sd["slurry"] > 0)
    assert capped.any() and (np.abs(slu_n[capped] - o["thickness"][capped]) <= geometry.DEPTH_DS_M).all()
    # a normal pointing out of the body (the file's winding, not orient_outward's) is detected: NaN, counted
    n = d["normal"].copy()
    wrong = np.flatnonzero(sd["exact"] & (sd["col"] == 2))[:5]
    n[wrong] *= -1.0
    liq_w, slu_w, info = layer_depths(f, loc, n, d["centroid"], return_info=True)
    assert info["n_wrong_way"] == len(wrong)
    assert np.isnan(slu_w[wrong]).all() and np.isnan(liq_w[wrong]).all()
    assert np.isfinite(np.delete(slu_w, wrong)).all()
    z = classify_zones(f.patch["film_thickness"], liq_w, slu_w, f.patch["delta_m"])
    assert (z.zone[wrong] == ZONE_INVALID).all()
    with pytest.raises(ValueError):
        layer_depths(f, loc, d["normal"], d["centroid"], f_l=np.zeros(3))


def test_sphere_slurry_depth_within_one_element_of_the_radial_profile(oriented):
    """Frames 1 and 2 of the sphere run: T = 300 + (T_wall - 300)(r/R)^2 with the linear f_l (README)."""
    R = syn.SPHERE_R
    answers = syn.sphere_answers()
    o0 = oriented[("sphere_run", 0)]
    _, slu0, info0 = layer_depths(o0["frame"], o0["loc"], o0["d"]["normal"], o0["d"]["centroid"], return_info=True)
    assert info0["n_marched"] == 0 and (slu0 == 0.0).all()             # 300 K: nothing marched
    measured = {1: (0.61e-3, 5.98e-3), 2: (-1.02e-3, 4.40e-3)}         # 2026-10-07, committed fixture
    for k in (1, 2):
        o = oriented[("sphere_run", k)]
        f, d, loc = o["frame"], o["d"], o["loc"]
        liq, slu, info = layer_depths(f, loc, d["normal"], d["centroid"], max_depth=o["thickness"], return_info=True)
        on = (np.abs(np.linalg.norm(f.points[f.faces], axis=2) - R) < 1e-9).all(axis=1)   # not the staircase
        assert on.sum() == (1194 if k == 1 else 1184)
        err = slu[on] - answers["slurry_depth_m"][k]
        assert np.abs(err).max() <= syn.SPHERE_H                       # within one element size (plan)
        assert err.min() == pytest.approx(measured[k][0], abs=1e-5) and err.max() == pytest.approx(measured[k][1],
                                                                                                    abs=1e-5)
        # a priori: the crossing's tetrahedron straddles the f_l = 0.5 sphere
        r_half = R * np.sqrt((0.5 * (syn.T_SOLIDUS + syn.T_LIQUIDUS) - 300.0) / (syn.SPHERE_T_WALL[k] - 300.0))
        ok = ~info["left_body"]
        x = d["centroid"][ok] - slu[ok][:, None] * d["normal"][ok]
        idx, _ = loc.locate(x)
        assert (idx >= 0).all()
        xt = f.points[f.tets[idx]]
        diam = np.max([np.linalg.norm(xt[:, i] - xt[:, j], axis=1) for i in range(4) for j in range(i + 1, 4)], axis=0)
        assert (np.abs(np.linalg.norm(x, axis=1) - r_half) <= diam + geometry.DEPTH_DS_M).all()
        # liquid: P1 holds none below the surface (the first interior nodes are under the liquidus), against the
        # continuum's 2.010 mm on frame 2
        assert (liq[on] <= 1e-9).all() and (liq <= slu).all()
    assert answers["liquid_depth_m"][2] == pytest.approx(2.010e-3, abs=1e-6)

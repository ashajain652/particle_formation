"""Tests of the level set's markers (design sections 9.1, 9.2, 9.6 and 21.5, stage 3 Task 5): markers on phi carrying one
skin each, snapped to the zero contour, given the area of phi each owns (`Markers.tessellate`: ownership with the
same-side and same-piece rules; plan 03b Task 6 retired the tangent-plane Voronoi faces), split, merged by the hand-over
rule, trapped and released, with the skins' books exact."""
import numpy as np
import pytest

from reentry_model.levelset.boxmesh import BoxBase, build
from reentry_model.levelset.cutcell import patches
from reentry_model.levelset.levelset import pieces, redistance
from reentry_model.levelset import markers as markers_mod
from reentry_model.levelset.markers import Markers
from reentry_model.levelset.skin import Skins, SkinSettings

from cases import DESIGN_FINE, DESIGN_SPACING, box_sdf, refined, seeded
from materials import ConstMaterial, LinearMeltMaterial

R = 25e-3                  # m, the sphere
H0 = 2e-3                  # m, base cells of the coarse test meshes
SETTINGS = SkinSettings(cell=20e-6)
FLAT_DEPTH = 0.8e-3        # m, the flat columns' depth in these tests


def _mesh(half=31e-3, h0=H0):
    return build(BoxBase.around((-half,) * 3, (half,) * 3, h0))


def _sphere(mesh, c=(0.0, 0.0, 0.0), r=R):
    return np.linalg.norm(mesh.points - np.asarray(c), axis=1) - r


def _filled(markers, mat, rng=None):
    """Markers' skins filled with metal: flat columns of FLAT_DEPTH of solid at 600 K, or random profiles."""
    sk = markers.skins
    if rng is None:
        sk.fill(mat.rho * sk.area * FLAT_DEPTH, 600.0)
        return sk
    n, K = sk.n, 30
    count = rng.integers(3, K + 1, n)
    m = rng.uniform(0.2, 2.0, (n, K)) * mat.rho * sk.area[:, None] * 20e-6
    T = rng.uniform(600.0, 1000.0, (n, K))
    sk.install(m, m * mat.enthalpy(T), count, rng.normal(0.0, 1.0, (n, K, 3)))
    return sk


def _arrays(sk):
    out = [a.copy() for a in (sk.m, sk.H, sk.U, sk.count, sk.area, sk.kappa)]
    return out + ([sk.profile.copy(), sk.profile_levels.copy()] if sk.profile is not None else [])


def test_markers_are_created_one_per_patch_on_the_zero_contour():
    mesh = _mesh()
    phi = _sphere(mesh)
    P = patches(mesh, phi)
    mk = Markers.create(mesh, phi, None, material=ConstMaterial(), settings=SETTINGS)
    assert mk.n == len(P.cell)
    assert np.array_equal(mk.patch, np.arange(mk.n)) and np.array_equal(mk.skin_index, np.arange(mk.n))
    assert np.allclose(mk.area, P.area)
    assert np.abs(np.linalg.norm(mk.x, axis=1) - R).max() < 0.05 * H0   # the P1 contour lies within h^2/(8R) of the sphere
    assert np.all(np.einsum("ij,ij->i", mk.normal, mk.x / R) > 0.99)
    assert np.allclose(mk.kappa, 1.0 / R, rtol=0.1)
    assert mk.spacing == pytest.approx(np.sqrt(P.area.mean()))
    assert mk.skins.n == mk.n and np.array_equal(mk.skins.area, mk.area)
    # a selection by mask makes a marker on each selected patch only
    sel = P.centroid[:, 2] > 0.0
    up = Markers.create(mesh, phi, sel)
    assert up.n == sel.sum() and up.skins is None
    # snap carries a marker displaced along the normal back onto the contour
    mk.x = mk.x * (1.0 + 0.3 * H0 / R)
    mk.snap(mesh, phi)
    assert np.abs(np.linalg.norm(mk.x, axis=1) - R).max() < 0.05 * H0   # the P1 contour lies within h^2/(8R) of the sphere


def test_areas_sum_to_the_surface_area():
    """Every marker on the sphere owns a part of phi, the areas sum to phi's, the ownership `tessellate` returns holds
    them as its rows, and each marker's longest direction (`axis`, which a split follows) is a unit vector in its tangent
    plane. (Plan 03b Task 6: the faces, their closure and the cap went with the tangent-plane cells; the division's
    boundaries replace them, tested in test_levelset_division.py.)"""
    mesh = _mesh()
    phi = _sphere(mesh)
    P = patches(mesh, phi)
    mk = Markers.create(mesh, phi, None)
    W = mk.tessellate(mesh, phi)
    assert mk.area.sum() == pytest.approx(P.area.sum(), rel=1e-12)
    assert np.all(mk.area > 0.0)
    assert W.shape == (mk.n, len(P.area))
    np.testing.assert_allclose(np.asarray(W.sum(axis=1)).ravel(), mk.area, rtol=1e-14, atol=0.0)
    np.testing.assert_allclose(np.linalg.norm(mk.axis, axis=1), 1.0, rtol=1e-12)
    assert np.abs(np.einsum("ij,ij->i", mk.axis, mk.normal)).max() < 1e-12
    for gone in ("faces", "cell_area", "capped"):
        assert not hasattr(mk, gone)


CLOSE_GAP = 0.7e-3         # m: the gas between the two cubes of the same-piece test (test_levelset_division's pair)


def test_a_marker_owns_area_only_on_its_own_piece():
    """Plan 03b Task 6 (the controller's ruling T1-i carried to the markers' areas): two 6 mm cubes 0.7 mm apart across
    gas, markers seeded on both as the body seeds them. Facing faces lie side by side across the gap, so the nearest
    marker to a point of one cube's face may stand on the other (the plain ownership crosses the gap); the markers'
    areas (`tessellate`) do not: the areas of each cube's markers sum to that cube's part of phi to round-off, and every
    piece of phi a marker owns lies on its own cube."""
    side = np.full(3, 6e-3)
    off = 0.13e-3
    a_lo = np.array([off - 6e-3 - 0.5 * CLOSE_GAP, off - 3e-3, off - 3e-3])
    b_lo = np.array([off + 0.5 * CLOSE_GAP, off - 3e-3, off - 3e-3])
    fa, fb = box_sdf(a_lo, a_lo + side), box_sdf(b_lo, b_lo + side)
    sdf = lambda p: np.minimum(fa(p), fb(p))
    base = BoxBase.around((-9.8e-3, -5.8e-3, -5.8e-3), (9.8e-3, 5.8e-3, 5.8e-3), 2.8e-3)
    mesh = refined(base, sdf, 0.35e-3)
    phi = sdf(mesh.points)
    mk = seeded(mesh, phi, DESIGN_SPACING, base.lo)
    P = patches(mesh, phi)
    on_b = mk.x[:, 0] >= off
    tri_b = P.centroid[:, 0] >= off
    plain = markers_mod.owned(mk.x, P, markers_mod.OWN_SHARE * DESIGN_SPACING).tocoo()
    assert np.any(on_b[plain.row] != tri_b[plain.col])                  # the nearest marker alone crosses the gap
    W = mk.tessellate(mesh, phi).tocoo()
    label, n_pieces = pieces(mesh, phi)
    print("two cubes {:.1f} mm apart: {} pieces, {} markers; areas {:.3e} and {:.3e} m2 against their cubes' {:.3e} and "
          "{:.3e}".format(CLOSE_GAP * 1e3, n_pieces, mk.n, mk.area[~on_b].sum(), mk.area[on_b].sum(),
                          P.area[~tri_b].sum(), P.area[tri_b].sum()))
    assert n_pieces == 2
    assert not np.any(on_b[W.row] != tri_b[W.col])
    assert mk.area[on_b].sum() == pytest.approx(P.area[tri_b].sum(), rel=1e-12)
    assert mk.area[~on_b].sum() == pytest.approx(P.area[~tri_b].sum(), rel=1e-12)


def test_rigid_motion_leaves_skins_unchanged():
    """Markers carried with a rigid motion of phi keep their skins bit for bit -- the cells, the areas, the curvatures
    and the measured profile (plan 03b Task 5: every skin's cross-section is a measured one; here each column narrowing
    by 30 % to H, so that the rows are not all alike) -- and their areas follow phi. (The division's own invariance is
    Task 1's to test: dividing this 25 mm sphere's thousands of markers took half an hour of the unit tier.)"""
    mat = LinearMeltMaterial()
    mesh = _mesh(half=33e-3)
    phi = _sphere(mesh)
    mk = Markers.create(mesh, phi, None, material=mat, settings=SETTINGS)
    _filled(mk, mat, np.random.default_rng(3))
    mk.tessellate(mesh, phi)
    levels = np.linspace(0.0, SETTINGS.depth, 25)
    mk.skins.set_profile(levels, mk.skins.area[:, None] * (1.0 - 0.3 * levels / SETTINGS.depth)[None, :])
    A0, skins0 = mk.area.copy(), _arrays(mk.skins)
    # translation by whole base cells: the mesh-relative geometry is identical
    shift = np.array([2 * H0, -H0, 3 * H0])
    mk.carry(lambda x: np.broadcast_to(shift, x.shape), 1.0)
    phi_t = _sphere(mesh, c=shift)
    mk.snap(mesh, phi_t)
    mk.tessellate(mesh, phi_t)
    for a, b in zip(skins0, _arrays(mk.skins)):
        assert np.array_equal(a, b)
    assert np.abs(mk.area / A0 - 1.0).max() < 1e-10
    assert (mk.splits, mk.merges, mk.trapped) == (0, 0, 0)
    # a rotation about the centre (phi unchanged): a Kuhn mesh is not rotation-symmetric, so 1 % on the areas (the areas
    # of phi each marker owns, Task 11's fix round 1, C1 (cited finding F2): 0.50 % at most measured once, 0.31 % for
    # the tangent-plane cells)
    c = shift
    from scipy.spatial.transform import Rotation
    rot = Rotation.from_rotvec([0.1, 0.2, 0.25]).as_matrix()
    mk.carry(lambda x: (x - c) @ rot.T + c - x, 1.0)      # the exact rotation as one step's displacement
    mk.snap(mesh, phi_t)
    mk.tessellate(mesh, phi_t)
    for a, b in zip(skins0, _arrays(mk.skins)):
        assert np.array_equal(a, b)
    assert np.abs(mk.area / A0 - 1.0).max() < 0.01


OWNED_EXACT = 1e-12        # each marker's area against its Voronoi cell's on a plane, relative
OWNED_ROTATION = 5e-3      # the largest change of a marker's area under a rigid rotation at the body's spacing


def test_each_marker_s_area_is_the_area_of_phi_nearest_it():
    """Task 11's fix round 1, C1 (the controller's ruling; cited finding F2): one partition of phi, each marker's area
    the area of phi nearer it than any other marker. On a plane tilted across a 0.25 mm mesh, 150 markers scattered on
    it: the area of every marker whose Voronoi cell in the plane (scipy's) lies 1 mm or more inside the mesh is that
    cell's area to OWNED_EXACT, and the areas sum to phi's to round-off."""
    from scipy.spatial import Voronoi
    nrm = np.array([0.3, -0.2, 1.0])
    nrm /= np.linalg.norm(nrm)
    e1 = np.cross(nrm, [1.0, 0.0, 0.0])
    e1 /= np.linalg.norm(e1)
    e2 = np.cross(nrm, e1)
    c = np.array([0.11e-3, -0.07e-3, 0.05e-3])
    mesh = build(BoxBase.around((-5e-3,) * 3, (5e-3,) * 3, 0.25e-3))
    phi = (mesh.points - c) @ nrm
    rng = np.random.default_rng(5)
    uv = rng.uniform(-4.5e-3, 4.5e-3, (150, 2))
    x = c + uv[:, :1] * e1 + uv[:, 1:] * e2
    mk = Markers(x, np.tile(nrm, (len(x), 1)), np.zeros(len(x)), np.zeros((len(x), 2)), np.arange(len(x)),
                 np.zeros(len(x), np.int64), 0.8e-3)
    mk.tessellate(mesh, phi)
    P = patches(mesh, phi)
    vor = Voronoi(uv)
    inner = []
    for i, r in enumerate(vor.point_region):
        reg = vor.regions[r]
        if not reg or -1 in reg:
            continue
        poly = vor.vertices[reg]
        X = c + poly[:, :1] * e1 + poly[:, 1:] * e2
        if np.all(np.abs(X) <= 4e-3):
            ctr = poly.mean(axis=0)
            ang = np.argsort(np.arctan2(*(poly - ctr).T[::-1]))
            q = poly[ang]
            inner.append((i, 0.5 * abs(np.dot(q[:, 0], np.roll(q[:, 1], -1)) - np.dot(q[:, 1], np.roll(q[:, 0], -1)))))
    idx = np.array([i for i, _ in inner])
    want = np.array([v for _, v in inner])
    err = np.abs(mk.area[idx] / want - 1.0)
    print("owned areas on a plane: {} of {} markers checked, within {:.1e} of their Voronoi cells".format(
        len(idx), mk.n, err.max()))
    assert len(idx) >= 50
    assert err.max() <= OWNED_EXACT
    assert mk.area.sum() == pytest.approx(P.area.sum(), rel=1e-12)


def _ownership_at_ffae441(x, tris, tri_patch, n_patch, size):
    """`markers.ownership` as it stood at ffae441, before the same-side rule and the shared boundaries (plan 03b Task 1),
    copied unchanged: the reference its default call must reproduce bit for bit."""
    from scipy.sparse import coo_matrix, csr_matrix
    from scipy.spatial import cKDTree
    x = np.asarray(x, dtype=float).reshape(-1, 3)
    n = len(x)
    t = np.asarray(tris, dtype=float).reshape(-1, 3, 3)
    tp = np.asarray(tri_patch, dtype=np.int64)
    if not n or not len(t):
        return csr_matrix((n, n_patch))
    tree = cKDTree(x)
    m = min(markers_mod.OWN_CANDIDATES, n)
    a = 0.5 * np.linalg.norm(np.cross(t[:, 1] - t[:, 0], t[:, 2] - t[:, 0]), axis=1)
    _, o = tree.query(t.reshape(-1, 3))
    o = np.asarray(o, dtype=np.int64).reshape(-1, 3)
    rows, cols, vals = [], [], []
    while len(t):
        same = (o[:, 0] == o[:, 1]) & (o[:, 0] == o[:, 2])
        rows.append(o[same, 0]), cols.append(tp[same]), vals.append(a[same])
        t, tp, o, a = t[~same], tp[~same], o[~same], a[~same]
        if not len(t):
            break
        c = t.mean(axis=1)
        r = np.linalg.norm(t - c[:, None, :], axis=2).max(axis=1)
        d, k = tree.query(c, k=min(m + 1, n))
        d, k = d.reshape(len(t), -1), np.asarray(k, dtype=np.int64).reshape(len(t), -1)
        reach = d[:, 0] + 2.0 * r
        live = (d[:, :m] <= reach[:, None]).sum(axis=1)
        last = np.sqrt(2.0 * a) <= size
        done = last | (d.shape[1] == m) | (d[:, -1] > reach)
        if done.any():
            share = markers_mod._divide(t[done], a[done], x, k[done, :m], np.where(last[done], m, live[done]))
            keep = share > 0.0
            rows.append(k[done, :m][keep]), cols.append(np.broadcast_to(tp[done][:, None], share.shape)[keep])
            vals.append(share[keep])
        t, tp, o, a = t[~done], tp[~done], o[~done], a[~done]
        if not len(t):
            break
        A, B, C = t[:, 0], t[:, 1], t[:, 2]
        AB, BC, CA = 0.5 * (A + B), 0.5 * (B + C), 0.5 * (C + A)
        _, om = tree.query(np.concatenate([AB, BC, CA]))
        oAB, oBC, oCA = np.asarray(om, dtype=np.int64).reshape(3, -1)
        oA, oB, oC = o[:, 0], o[:, 1], o[:, 2]
        t = np.concatenate([np.stack(v, axis=1) for v in ((A, AB, CA), (AB, B, BC), (CA, BC, C), (AB, BC, CA))])
        o = np.concatenate([np.stack(v, axis=1) for v in ((oA, oAB, oCA), (oAB, oB, oBC), (oCA, oBC, oC),
                                                          (oAB, oBC, oCA))])
        tp, a = np.tile(tp, 4), np.tile(0.25 * a, 4)
    W = coo_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(n, n_patch)).tocsr()
    W.sum_duplicates()
    return W


def test_ownership_with_no_normals_is_unchanged():
    """Plan 03b Task 1: `ownership` without normals and without boundaries is today's (ffae441's), bit for bit, on the
    test sphere with markers seeded 2 mm apart; asking for the boundaries leaves the areas as they were."""
    mesh = _mesh()
    phi = _sphere(mesh)
    P = patches(mesh, phi)
    mk = seeded(mesh, phi, DESIGN_SPACING, mesh.base.lo)
    size = markers_mod.OWN_SHARE * DESIGN_SPACING
    old = _ownership_at_ffae441(mk.x, P.tris, P.tri_patch, len(P.area), size)
    new = markers_mod.ownership(mk.x, P.tris, P.tri_patch, len(P.area), size)
    both, (i, j, length) = markers_mod.ownership(mk.x, P.tris, P.tri_patch, len(P.area), size, boundaries=True)
    assert old.nnz > len(P.area)                                           # some patches are divided
    for W in (new, both):
        assert W.shape == old.shape
        for name in ("indptr", "indices", "data"):
            assert np.array_equal(getattr(W, name), getattr(old, name))
    assert len(i) > mk.n and np.all(i < j) and np.all(length > 0.0)


TIE_OFF = 0.13e-3          # m: the flat face of the tie test off the mesh's node planes (the stage 3 box's offset)


def test_a_boundary_along_a_mesh_edge_is_counted_once():
    """Plan 03b Task 6: markers a and b mirror-symmetric about a plane of the mesh's nodes (x = -1 and +1 mm about x = 0),
    so that the plane halfway between them runs along the edges of a flat face's triangles and the corners on it are
    exactly as near one as the other, and a third, c, on that plane 2 mm along it (z = 2 mm), so that the triangles
    beside the edge near the three's meeting point are divided among them. On the face 12 mm across the boundaries are
    the plane's: a's with b from z = -6 mm to the meeting point at z = 0.75 mm, 6.75 mm, and a's and b's with c, each
    3 sqrt(5) mm to the face's edge; the division of the triangles gives each to 1e-12, with and without the same-side
    rule (before, an edge counted in both triangles beside it or in neither, by the candidates' order and the
    nearest-marker query's choice at its corners: on a 10 mm sphere centred on the mesh's planes one column's boundary
    fluxes missed closing by 22 %); and the areas are the ones asking for no boundaries gives, bit for bit."""
    base = BoxBase((-6e-3, -2.5e-3, -6e-3), 0.5e-3, (24, 6, 24))
    mesh = build(base)
    phi = mesh.points[:, 1] - TIE_OFF
    P = patches(mesh, phi)
    size = markers_mod.OWN_SHARE * DESIGN_SPACING
    x = np.array([[-1e-3, TIE_OFF, 0.0], [1e-3, TIE_OFF, 0.0], [0.0, TIE_OFF, 2e-3]])
    nrm = np.tile([0.0, 1.0, 0.0], (3, 1))
    assert np.any(P.tris[:, :, 0] == 0.0)                                     # corners on the halfway plane
    want = {(0, 1): 6.75e-3, (0, 2): 3.0 * np.sqrt(5.0) * 1e-3, (1, 2): 3.0 * np.sqrt(5.0) * 1e-3}
    for kw in ({}, dict(normals=nrm, tri_normal=P.normal[P.tri_patch])):
        W0 = markers_mod.ownership(x, P.tris, P.tri_patch, len(P.area), size, **kw)
        W, (i, j, length) = markers_mod.ownership(x, P.tris, P.tri_patch, len(P.area), size, boundaries=True, **kw)
        got = {(int(a), int(b)): float(v) for a, b, v in zip(i, j, length)}
        print("three markers{}: boundaries {} mm".format(" (same side)" if kw else "",
                                                       {k: round(v * 1e3, 9) for k, v in got.items()}))
        assert set(got) == set(want)
        for k, v in want.items():
            assert got[k] == pytest.approx(v, rel=1e-12), k
        for name in ("indptr", "indices", "data"):
            assert np.array_equal(getattr(W, name), getattr(W0, name))


def test_a_turned_sphere_s_markers_keep_their_areas_at_the_body_s_spacing():
    """Task 11's fix round 1, C1 (cited finding F2): a 10 mm sphere on the body's design cells (0.7 mm near the
    surface), markers seeded as the body seeds them, 2 mm apart, turned rigidly about its centre (phi unchanged, the
    markers carried and snapped): no marker's area, the area of phi it owns, changes by more than OWNED_ROTATION. The
    change is the faceted surface's own, the cut-cell patches meeting the turned markers differently. (Measured once:
    median 0.023 per cent, 95th percentile 0.071, largest 0.21; Task 6 measured 0.31 per cent at most for the
    tangent-plane cells on a 25 mm sphere with one marker per patch, where the owned areas now change by 0.50 per cent
    at most, `test_rigid_motion_leaves_skins_ unchanged`.)"""
    from scipy.spatial.transform import Rotation
    Rs = 10e-3
    base = BoxBase.around((-Rs - 5.6e-3,) * 3, (Rs + 5.6e-3,) * 3, 2.8e-3)
    mesh = refined(base, lambda x: np.linalg.norm(x, axis=1) - Rs, DESIGN_FINE)
    phi = redistance(mesh, np.linalg.norm(mesh.points, axis=1) - Rs, band=4e-3)
    mk = seeded(mesh, phi, DESIGN_SPACING, base.lo)
    mk.tessellate(mesh, phi)
    A0 = mk.area.copy()
    rot = Rotation.from_rotvec([0.1, 0.2, 0.25]).as_matrix()
    mk.carry(lambda x: x @ rot.T - x, 1.0)
    mk.snap(mesh, phi)
    mk.tessellate(mesh, phi)
    turned = np.abs(mk.area / A0 - 1.0)
    print("turned: {} markers, area change median {:.2e}, 95th {:.2e}, largest {:.2e}".format(
        mk.n, np.median(turned), np.percentile(turned, 95), turned.max()))
    assert turned.max() <= OWNED_ROTATION


def test_uniform_stretching_thins_skins_and_fills_with_the_implied_metal():
    """Markers carried radially from a 25 mm sphere to a 26 mm one keep their metal spread thinner over their larger
    areas; refilled at their bases to the flat column they now stand for (`Skins.exchange_base`, the fill of plan 03b
    Task 5, which replaced the draw to a minimum mass), each takes exactly the metal the larger column implies."""
    mat = ConstMaterial()
    mesh = _mesh()
    phi = _sphere(mesh)
    mk = Markers.create(mesh, phi, None, material=mat, settings=SETTINGS)
    mk.tessellate(mesh, phi)
    sk = mk.skins
    sk.area = mk.area[np.argsort(mk.skin_index)].copy()
    _filled(mk, mat)
    A_old, M_old, thick_old = sk.area.copy(), sk.m.sum(axis=1), sk.thickness()
    mk.carry(lambda x: 1e-3 * x / np.linalg.norm(x, axis=1)[:, None], 1.0)      # 25 to 26 mm, radially
    phi26 = _sphere(mesh, r=26e-3)
    mk.snap(mesh, phi26)
    mk.tessellate(mesh, phi26)
    sk.area[mk.skin_index] = mk.area
    ratio = sk.area / A_old
    assert np.median(ratio) == pytest.approx((26 / 25) ** 2, rel=5e-4)
    assert np.abs(ratio / (26 / 25) ** 2 - 1.0).max() < 0.01
    assert np.allclose(sk.thickness(), thick_old / ratio, rtol=1e-12)        # the same metal, spread thinner
    implied = sk.area * FLAT_DEPTH - A_old * FLAT_DEPTH                     # the new flat column's volume less the old
    ex = sk.exchange_base(np.arange(sk.n), implied, 600.0)
    drawn, H = ex.mass, ex.H
    assert np.allclose(drawn, mat.rho * implied, rtol=1e-12, atol=0.0)
    assert np.allclose(sk.m.sum(axis=1), mat.rho * sk.area * FLAT_DEPTH, rtol=1e-12, atol=0.0)
    assert np.allclose(sk.m.sum(axis=1), M_old + drawn, rtol=1e-14, atol=0.0)
    assert np.allclose(H, drawn * mat.enthalpy(600.0), rtol=1e-12)


def _books(sk):
    M, E, P = sk.totals()
    return M.sum(), E.sum(), P.sum(axis=0)


def test_split_and_merge_conserve_mass_and_enthalpy():
    mat = LinearMeltMaterial()
    mesh = _mesh()
    phi = _sphere(mesh)
    mk = Markers.create(mesh, phi, None, material=mat, settings=SETTINGS)
    _filled(mk, mat, np.random.default_rng(11))
    mk.tessellate(mesh, phi)
    mk.skins.area[mk.skin_index] = mk.area                 # the caller gives the skins their markers' areas
    M0, E0, P0 = _books(mk.skins)
    n0, A0 = mk.n, mk.area.sum()
    target = 0.6 * np.median(mk.area)
    big = np.flatnonzero(mk.area > 2 * target)
    assert 0 < len(big) < n0
    x_big, m_big = mk.x[big].copy(), mk.skins.m[mk.skin_index[big]].copy()
    mk.split(target)
    assert mk.splits == len(big) and mk.n == n0 + len(big) and mk.skins.n == mk.n
    assert len(np.unique(mk.skin_index)) == mk.n
    M1, E1, P1 = _books(mk.skins)
    assert M1 == pytest.approx(M0, rel=1e-14) and E1 == pytest.approx(E0, rel=1e-14)
    assert np.allclose(P1, P0, rtol=0, atol=1e-14 * np.abs(mk.skins.m.sum()))
    assert mk.area.sum() == pytest.approx(A0, rel=1e-14)
    new = np.arange(n0, mk.n)
    gap = np.linalg.norm(mk.x[big] - mk.x[new], axis=1)
    assert np.allclose(gap, 0.5 * np.sqrt(target), rtol=1e-12)
    assert np.allclose(0.5 * (mk.x[big] + mk.x[new]), x_big, rtol=0, atol=1e-15)
    # each half keeps the profile with half the mass of every cell
    assert np.array_equal(mk.skins.m[mk.skin_index[big]], 0.5 * m_big)
    assert np.array_equal(mk.skins.m[mk.skin_index[new]][:, :m_big.shape[1]], 0.5 * m_big)
    # the split pairs are half a spacing apart: a merge at twice that spacing joins them again
    n1 = mk.n
    mk.merge(2.0 * np.sqrt(target))
    assert mk.merges >= len(big) and mk.n == n1 - mk.merges and mk.skins.n == mk.n
    assert np.array_equal(np.sort(mk.skin_index), np.arange(mk.n))
    M2, E2, P2 = _books(mk.skins)
    assert M2 == pytest.approx(M0, rel=1e-14) and E2 == pytest.approx(E0, rel=1e-14)
    assert np.allclose(P2, P0, rtol=0, atol=1e-14 * np.abs(mk.skins.m.sum()))
    assert mk.area.sum() == pytest.approx(A0, rel=1e-14)
    assert np.allclose(mk.skins.area[mk.skin_index], mk.area, rtol=1e-14)
    assert mk.trapped == 0


MERGE_MUSH = (885.0, 905.0)   # K: the merge test's columns, within LinearMeltMaterial's mushy range (specific volume
                              # linear in enthalpy, so mixing keeps volume; density varying with depth and column)


def test_a_merge_keeps_each_depth_s_metal_at_its_depth():
    """Plan 03b Task 6 (the controller's ruling T6-c; design 21.5: metal keeps its depth when its owner changes): two
    markers 0.2 mm apart, each with a full column of a flat measured profile (cross-sections 3 and 1 mm2, the skin depth
    deep) of mushy metal whose temperature, and so density, varies with depth differently in each (MERGE_MUSH), moving
    differently, merge (`Markers.merge`, closer than a third of the spacing): the merged column's mass, enthalpy and
    momentum are the pair's to 1e-14, every cell keeps its depth span (1e-12 of the depth) and holds the pair's metal from
    that span (1e-12 of its mass). (By the hand-over rule's relative mass depth, the rule before, the denser column's
    metal shifts in depth against the lighter's.)"""
    mat = LinearMeltMaterial()
    depth = SETTINGS.depth
    levels = np.linspace(0.0, depth, 25)
    area = np.array([3e-6, 1e-6])
    sk = Skins(mat, area, np.zeros((2, 2)), SETTINGS)
    sk.set_profile(levels, area[:, None] * np.ones(len(levels))[None, :])
    lo, hi = MERGE_MUSH
    sk.fill(sk.density(np.full(2, hi)) * sk.column_volume(), hi)
    sk.regrid()
    dz = sk.dz()
    z = np.concatenate([np.zeros((2, 1)), np.cumsum(dz, axis=1)], axis=1)
    mid = 0.5 * (z[:, 1:] + z[:, :-1])
    T = np.where(sk.valid, hi - (hi - lo) * (mid / depth) ** np.array([[0.5], [2.0]]), lo)
    span = np.diff(sk._volume_at(z), axis=1)                     # each cell's span's volume, kept full at its T
    m = np.where(sk.valid, sk.density(T) * span, 0.0)
    U = np.zeros(m.shape + (3,))
    U[0, :, 0], U[1, :, 1] = 0.3, -0.7
    sk.install(m, m * sk.props.enthalpy(T), sk.count, U)
    assert np.array_equal(sk.count[0], sk.count[1]) and np.abs(z[0] - z[1]).max() <= 1e-15    # the same cells in depth
    K = int(sk.count[0])
    held, spans = sk.m[:, :K].sum(axis=0), np.cumsum(sk.dz(), axis=1)[0, :K]
    M0, E0, P0 = _books(sk)
    x = np.array([[0.0, 0.0, 0.0], [0.2e-3, 0.0, 0.0]])
    mk = Markers(x, np.tile([0.0, 0.0, 1.0], (2, 1)), area, np.zeros((2, 2)), np.arange(2), np.full(2, -1), 1e-3,
                 skins=sk)
    assert mk.merge(1e-3) == 1 and mk.n == 1 and sk.n == 1
    M1, E1, P1 = _books(sk)
    depth_miss = float(np.abs(np.cumsum(sk.dz(), axis=1)[0, :K] - spans).max() / depth)
    level_miss = float(np.abs(sk.m[0, :K] - held).max() / sk.m[0].sum())
    print("merge by depth: {} cells; books mass {:.1e} enthalpy {:.1e}; cell depths {:.1e} of H, cell metal {:.1e} of the "
          "column".format(sk.count[0], M1 / M0 - 1, E1 / E0 - 1, depth_miss, level_miss))
    assert sk.count[0] == K
    assert M1 == pytest.approx(M0, rel=1e-14) and E1 == pytest.approx(E0, rel=1e-14)
    assert np.abs(P1 - P0).max() <= 1e-14 * np.abs(P0).max()
    assert depth_miss <= 1e-12 and level_miss <= 1e-12


def test_the_hand_over_rule_mixes_layers_at_the_same_relative_depth():
    mat = LinearMeltMaterial()
    sk = Skins(mat, np.array([4e-6, 2e-6]), np.zeros((2, 2)), SETTINGS)
    m = np.array([[1.0, 1.0, 2.0, 0.0], [1.0, 3.0, 0.0, 0.0]]) * 1e-6
    T = np.array([[1000.0, 800.0, 600.0, 0.0], [900.0, 500.0, 0.0, 0.0]])
    U = np.zeros((2, 4, 3))
    U[0, :, 0], U[1, :, 1] = 1.0, -2.0
    sk.install(m, m * mat.enthalpy(T), [3, 2], U)
    M0, E0, P0 = (np.sum(t, axis=0) for t in sk.totals())
    sk.combine([0], [1])
    assert sk.count[1] == 0 and sk.area[0] == pytest.approx(6e-6) and sk.area[1] == 0.0
    M1, E1, P1 = (np.sum(t, axis=0) for t in sk.totals())
    assert M1 == pytest.approx(M0, rel=1e-14) and E1 == pytest.approx(E0, rel=1e-14)
    assert np.allclose(P1, P0, rtol=1e-14)
    # boundaries at the mass fractions 0.25, 0.5 (column 0) and 0.25 (column 1): cells [0, .25], [.25, .5], [.5, 1]
    assert sk.count[0] == 3
    assert np.allclose(sk.m[0, :3], [2e-6, 2e-6, 4e-6], rtol=1e-14)
    top = sk.H[0, 0]
    assert top == pytest.approx(1e-6 * mat.enthalpy(1000.0) + 1e-6 * mat.enthalpy(900.0), rel=1e-14)
    second = 1e-6 * mat.enthalpy(800.0) + 1e-6 * mat.enthalpy(500.0)
    assert sk.H[0, 1] == pytest.approx(second, rel=1e-14)
    gone = sk.remove([1])
    assert sk.n == 1 and all(np.all(g == 0.0) for g in gone)


def test_a_trapped_marker_returns_its_skin_exactly():
    """Plan 03b, planning ruling 6: the markers trapped in a fold (a second sphere overlapping the first) are released,
    each giving its skin to the nearest marker that stays, by depth (`Skins.combine(..., by="depth")`), not to the volume:
    the skins' mass, enthalpy and momentum are unchanged to round-off, each keeper holds its own skin and exactly what
    the released markers nearest it gave, and every other marker keeps its own skin."""
    mat = LinearMeltMaterial()
    r, h = 10e-3, 1e-3
    mesh = build(BoxBase.around((-12e-3, -12e-3, -12e-3), (24e-3, 12e-3, 12e-3), h))
    phi_a = _sphere(mesh, r=r)
    mk = Markers.create(mesh, phi_a, None, material=mat, settings=SETTINGS)
    _filled(mk, mat, np.random.default_rng(5))
    M0, E0, P0 = _books(mk.skins)
    rows0 = [a.copy() for a in mk.skins.totals()]
    x0 = mk.x.copy()
    c_b = np.array([1.2 * r, 0.0, 0.0])
    phi = np.minimum(phi_a, _sphere(mesh, c=c_b, r=r))      # a second sphere overlapping the first
    idx = mk.trapped_markers(mesh, phi, h)
    depth_b = np.linalg.norm(x0 - c_b, axis=1) - r
    assert np.all(depth_b[idx] < 0.0)                       # only markers inside the other sphere are trapped
    assert np.all(np.isin(np.flatnonzero(depth_b < -h), idx))
    assert np.all(~np.isin(np.flatnonzero(depth_b > 0.5 * h), idx))
    assert len(idx) > 10
    rows = mk.skin_index[idx]
    keep = np.setdiff1d(np.arange(len(x0)), idx)
    x, keeper, (mass, H, P) = mk.release(idx)
    assert mk.trapped == len(idx) and mk.n == len(x0) - len(idx) and mk.skins.n == mk.n
    assert np.array_equal(x, x0[idx])
    assert np.array_equal(mk.x, x0[keep])
    assert np.allclose(mass, rows0[0][rows], rtol=1e-14, atol=0.0)
    assert np.allclose(H, rows0[1][rows], rtol=1e-14, atol=0.0)
    assert np.allclose(P, rows0[2][rows], rtol=1e-14, atol=0.0)
    dist = np.linalg.norm(x0[keep][None, :, :] - x[:, None, :], axis=2)
    np.testing.assert_allclose(dist[np.arange(len(x)), keeper], dist.min(axis=1), rtol=1e-12, atol=0.0)   # nearest
    M1, E1, P1 = _books(mk.skins)
    assert M1 == pytest.approx(M0, rel=1e-14) and E1 == pytest.approx(E0, rel=1e-14)
    assert np.allclose(P1, P0, rtol=0, atol=1e-14 * np.abs(mk.skins.m.sum()))
    # each keeper: its own skin and what was given it; every other survivor its own skin alone
    got = np.bincount(keeper, mass, minlength=mk.n)
    totals = mk.skins.totals()[0][mk.skin_index]
    np.testing.assert_allclose(totals, rows0[0][keep] + got, rtol=1e-14, atol=0.0)
    untouched = np.flatnonzero(got == 0.0)              # their skins as they were (the padding grew: sums to round-off)
    np.testing.assert_allclose(totals[untouched], rows0[0][keep][untouched], rtol=1e-15, atol=0.0)


def test_assign_finds_the_patch_whose_cut_cell_holds_each_marker():
    mesh = _mesh()
    phi = _sphere(mesh)
    P = patches(mesh, phi)
    mk = Markers.create(mesh, phi, None)
    perm = np.random.default_rng(2).permutation(mk.n)
    mk.x = mk.x[perm]
    mk.assign(mesh, P)
    assert np.array_equal(mk.patch, perm)


def test_snap_converges_at_a_kink_of_the_p1_contour():
    """A P1 phi with a kink at a layer of nodes (slope 3 in the layer holding the zero contour, 1 in the layer above,
    the root 0.2 mm below the kink): the recovered-normal snap alone overshoots across the kink and back and cycles
    for ever 1.2 mm off phi; the cell-gradient Newton steps of the body's last snap round land on the contour (Task 9, fix round 4)."""
    from reentry_model.levelset.markers import SNAP_PASSES, _field
    from reentry_model.levelset.transfer import interpolate
    h = 1e-3
    mesh = build(BoxBase.around((-3e-3, -3e-3, -2e-3), (3e-3, 3e-3, 4e-3), h))
    z = mesh.points[:, 2]
    r = 0.8e-3                                                       # the zero contour, z = r, in the layer [0, h]
    phi = np.where(z <= h + 1e-12, 3.0 * (z - r), 3.0 * (h - r) + (z - h))
    x0 = np.array([[0.3e-3, -0.2e-3, 0.4e-3], [1.1e-3, 0.6e-3, 0.4e-3]])
    x = x0.copy()                                                    # the old snap: recovered-normal steps only
    for _ in range(20 * SNAP_PASSES):
        val, n = _field(mesh, phi, x)
        x = x - val[:, None] * n
    assert np.allclose(np.abs(interpolate(mesh, phi, x)), 1.2e-3, rtol=1e-9)   # still cycling, 1.2 mm off
    mk = Markers(x0, np.tile([0.0, 0.0, 1.0], (2, 1)), np.ones(2), np.zeros((2, 2)), [0, 1], [-1, -1], h)
    mk.snap(mesh, phi)                                               # the default reach: too far off for Newton
    assert np.allclose(mk.x, x0, atol=1e-15)                         # (back where it started, an even pass count)
    mk.snap(mesh, phi, reach=np.inf)                                 # the body's last snap round
    assert np.abs(interpolate(mesh, phi, mk.x)).max() <= 1e-15
    assert np.allclose(mk.x[:, 2], r, atol=1e-15) and np.allclose(mk.x[:, :2], x0[:, :2], atol=1e-15)
    assert np.allclose(mk.normal, [0.0, 0.0, 1.0])


def test_snap_converges_where_phi_is_far_from_a_distance():
    """A P1 phi rising only 0.2 m per metre along its normal (as in two cells of the full-size cube's second step, where
    |grad phi| was 0.20 and 0.88 and the recovered normal 30 and 73 degrees off the cell's gradient): each
    recovered-normal step removes only a fifth of the distance, so 40 of them from 1.6 mm leave the marker some 2e-7 m
    off, far above the body's 1e-10 m; one snap whose cell-gradient Newton steps reach it lands on the contour (fix round 4)."""
    from reentry_model.levelset.markers import SNAP_PASSES, _field
    from reentry_model.levelset.transfer import interpolate
    h = 1e-3
    mesh = build(BoxBase.around((-3e-3, -3e-3, -3e-3), (3e-3, 3e-3, 3e-3), h))
    phi = 0.2 * (mesh.points[:, 2] - 0.3e-3)
    x0 = np.array([[0.3e-3, -0.2e-3, 1.9e-3]])
    x = x0.copy()
    for _ in range(20 * SNAP_PASSES):
        val, n = _field(mesh, phi, x)
        x = x - val[:, None] * n
    assert abs(x[0, 2] - 0.3e-3) == pytest.approx(1.6e-3 * 0.8 ** 40, rel=1e-6)    # about 2.1e-7 m off
    mk = Markers(x0, [[0.0, 0.0, 1.0]], [1.0], [[0.0, 0.0]], [0], [-1], h)
    mk.snap(mesh, phi, reach=np.inf)
    assert abs(mk.x[0, 2] - 0.3e-3) <= 1e-15 and np.allclose(mk.x[0, :2], x0[0, :2], atol=1e-15)


def test_markers_pickled_before_task_9b_load_with_their_positions_as_the_last_tessellation():
    """Task 9b (fix round 1): markers pickled before `x_tess` and `normal_tess` existed (a state saved between steps,
    when they equal the positions and normals) load with them filled from `x` and `normal`; markers that have them keep
    their own."""
    import pickle
    mesh = _mesh()
    phi = _sphere(mesh)
    mk = Markers.create(mesh, phi, None)
    mk.tessellate(mesh, phi)
    del mk.x_tess, mk.normal_tess                         # as a pickle from before Task 9b holds them
    old = pickle.loads(pickle.dumps(mk))
    assert np.array_equal(old.x_tess, old.x) and np.array_equal(old.normal_tess, old.normal)
    assert old.x_tess is not old.x
    mk.x_tess, mk.normal_tess = mk.x + 1e-6, -mk.normal
    new = pickle.loads(pickle.dumps(mk))
    assert np.array_equal(new.x_tess, mk.x + 1e-6) and np.array_equal(new.normal_tess, -mk.normal)
    # markers pickled with the retired modes (stage 3 Task 10 to plan 03b Task 5) load without them, and so do markers
    # pickled with the retired tangent-plane faces, cells and cap flags (to plan 03b Task 6)
    mk.top, mk.hold = np.zeros(mk.n, bool), np.zeros(mk.n, np.int64)
    mk.faces, mk.cell_area, mk.capped = (np.zeros(0), np.zeros(0), np.zeros(0), np.zeros((0, 3))), mk.area.copy(), \
        np.zeros(mk.n, bool)
    retired = pickle.loads(pickle.dumps(mk))
    for gone in ("top", "hold", "faces", "cell_area", "capped"):
        assert not hasattr(retired, gone)


def test_a_released_marker_gives_its_skin_to_the_nearest_marker_on_its_own_side():
    """Plan 03b Task 6 (optional, the review's note on `Markers.release`, cited finding F22): a released marker's keeper
    is the nearest marker that may own its part of phi by the same-side rule (normals not opposed,
    `markers.OPPOSED_NORMALS`), as ownership elsewhere chooses. Two rows of markers 1 mm apart facing each other across
    0.4 mm of gas (a slot's walls, normals +z and -z): a marker of the lower row released gives its skin to its own
    row's nearest marker, 1 mm away, not to the upper row's, 0.4 mm away; the books are exact."""
    mat = LinearMeltMaterial()
    n = 5
    row = np.arange(n) * 1e-3
    x = np.concatenate([np.column_stack([row, np.zeros(n), np.zeros(n)]),
                        np.column_stack([row, np.zeros(n), np.full(n, 0.4e-3)])])
    normal = np.concatenate([np.tile([0.0, 0.0, 1.0], (n, 1)), np.tile([0.0, 0.0, -1.0], (n, 1))])
    sk = Skins(mat, np.full(2 * n, 1e-6), np.zeros((2 * n, 2)), SETTINGS)
    sk.fill(mat.rho * 1e-6 * FLAT_DEPTH * (1.0 + np.arange(2 * n)), 600.0)
    mk = Markers(x, normal, np.full(2 * n, 1e-6), np.zeros((2 * n, 2)), np.arange(2 * n), np.full(2 * n, -1), 1e-3,
                 skins=sk)
    M0 = sk.totals()[0].copy()
    _, keeper, (mass, _, _) = mk.release(np.array([2]))
    survivors = np.r_[np.arange(2), np.arange(3, 2 * n)]
    print("released marker 2: its keeper is marker {} at {} mm".format(survivors[keeper[0]],
                                                                         np.round(x[survivors[keeper[0]]] * 1e3, 2)))
    assert survivors[keeper[0]] in (1, 3)                              # its own row's neighbour, not marker 7 above it
    assert mass[0] == pytest.approx(M0[2], rel=1e-14)
    assert sk.m.sum() == pytest.approx(M0.sum(), rel=1e-14)

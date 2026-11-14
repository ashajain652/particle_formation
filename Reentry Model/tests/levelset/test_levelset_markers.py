"""Tests of the level set's markers (design sections 9.1, 9.2 and 9.6, stage 3 Task 5): markers on phi carrying one skin
each, snapped to the zero contour, tessellated by a Voronoi construction in each marker's tangent plane (Ruling 2),
split, merged by the hand-over rule and trapped, with the skins' books exact."""
import numpy as np
import pytest

from reentry_model.levelset.boxmesh import BoxBase, build
from reentry_model.levelset.cutcell import patches
from reentry_model.levelset.levelset import redistance
from reentry_model.levelset import markers as markers_mod
from reentry_model.levelset.markers import Markers
from reentry_model.levelset.skin import Skins, SkinSettings

from cases import DESIGN_FINE, DESIGN_SPACING, CAPPED_SHARE, OWNERSHIP_FACTOR, nearest_ratio, refined, seeded
from materials import ConstMaterial, LinearMeltMaterial

R = 25e-3                  # m, the sphere
H0 = 2e-3                  # m, base cells of the coarse test meshes
SETTINGS = SkinSettings(cell=20e-6)


def _mesh(half=31e-3, h0=H0):
    return build(BoxBase.around((-half,) * 3, (half,) * 3, h0))


def _sphere(mesh, c=(0.0, 0.0, 0.0), r=R):
    return np.linalg.norm(mesh.points - np.asarray(c), axis=1) - r


def _filled(markers, mat, rng=None):
    """Markers' skins filled with metal: flat columns of 0.8 mm of solid at 600 K, or random profiles."""
    sk = markers.skins
    if rng is None:
        sk.fill(mat.rho * sk.area * SETTINGS.stub, 600.0)
        return sk
    n, K = sk.n, 30
    count = rng.integers(3, K + 1, n)
    m = rng.uniform(0.2, 2.0, (n, K)) * mat.rho * sk.area[:, None] * 20e-6
    T = rng.uniform(600.0, 1000.0, (n, K))
    sk.install(m, m * mat.enthalpy(T), count, rng.normal(0.0, 1.0, (n, K, 3)))
    return sk


def _arrays(sk):
    return [a.copy() for a in (sk.m, sk.H, sk.U, sk.count, sk.area, sk.kappa)]


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
    mesh = _mesh()
    phi = _sphere(mesh)
    P = patches(mesh, phi)
    mk = Markers.create(mesh, phi, None)
    faces = mk.tessellate(mesh, phi)
    assert mk.area.sum() == pytest.approx(P.area.sum(), rel=1e-12)
    assert np.all(mk.area > 0.0)
    i, j, length, normal = faces
    assert np.all(i < j) and np.all(length > 0.0)
    assert np.allclose(np.linalg.norm(normal, axis=1), 1.0, rtol=1e-12)
    # a face's normal lies in the surface and points from i toward j
    assert np.all(np.einsum("ij,ij->i", normal, mk.x[j] - mk.x[i]) > 0.0)
    assert np.abs(np.einsum("ij,ij->i", normal, mk.normal[i] + mk.normal[j])).max() < 0.2
    # every marker is a neighbour of some other, and the cap binds almost nowhere on a closed surface
    assert len(np.setdiff1d(np.arange(mk.n), np.concatenate([i, j]))) == 0
    assert mk.capped.mean() < 0.02
    # each cell closes: the outward face vectors of an uncapped marker nearly cancel
    vec = np.zeros((mk.n, 3))
    np.add.at(vec, i, length[:, None] * normal)
    np.add.at(vec, j, -length[:, None] * normal)
    per = np.bincount(i, length, mk.n) + np.bincount(j, length, mk.n)
    ok = ~mk.capped
    assert np.median(np.linalg.norm(vec[ok], axis=1) / per[ok]) < 0.05


def test_rigid_motion_leaves_skins_unchanged():
    mat = LinearMeltMaterial()
    mesh = _mesh(half=33e-3)
    phi = _sphere(mesh)
    mk = Markers.create(mesh, phi, None, material=mat, settings=SETTINGS)
    _filled(mk, mat, np.random.default_rng(3))
    mk.tessellate(mesh, phi)
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
    # of phi each marker owns, Task 11's fix round 1, C1: 0.50 % at most measured once, 0.31 % for the tangent-plane cells)
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
    """Task 11's fix round 1, C1 (the controller's ruling): one partition of phi, each marker's area the area of phi
    nearer it than any other marker. On a plane tilted across a 0.25 mm mesh, 150 markers scattered on it: the area of
    every marker whose Voronoi cell in the plane (scipy's) lies 1 mm or more inside the mesh is that cell's area to
    OWNED_EXACT, and the areas sum to phi's to round-off."""
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


def test_a_turned_sphere_s_markers_keep_their_areas_at_the_body_s_spacing():
    """Task 11's fix round 1, C1: a 10 mm sphere on the body's design cells (0.7 mm near the surface), markers seeded as
    the body seeds them, 2 mm apart, turned rigidly about its centre (phi unchanged, the markers carried and snapped): no
    marker's area, the area of phi it owns, changes by more than OWNED_ROTATION. The change is the faceted surface's own,
    the cut-cell patches meeting the turned markers differently. (Measured once: median 0.023 per cent, 95th percentile
    0.071, largest 0.21; Task 6 measured 0.31 per cent at most for the tangent-plane cells on a 25 mm sphere with one
    marker per patch, where the owned areas now change by 0.50 per cent at most, `test_rigid_motion_leaves_skins_
    unchanged`.)"""
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


def test_uniform_stretching_thins_skins_and_draws_the_implied_metal():
    mat = ConstMaterial()
    mesh = _mesh()
    phi = _sphere(mesh)
    mk = Markers.create(mesh, phi, None, material=mat, settings=SETTINGS)
    mk.tessellate(mesh, phi)
    sk = mk.skins
    sk.area = mk.area[np.argsort(mk.skin_index)].copy()
    sk.kappa[:] = 0.0                                     # flat columns: the stub's volume is A0 times the stub
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
    m_min = sk.minimum_mass()
    drawn, H = sk.draw(600.0, 0.0, m_min)
    implied = sk.area * SETTINGS.stub - A_old * SETTINGS.stub                 # the new stub volume minus the old
    assert np.allclose(drawn, mat.rho * implied, rtol=1e-12, atol=0.0)
    assert np.allclose(sk.m.sum(axis=1), m_min, rtol=1e-12, atol=0.0)
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
    x, mass, H, P = mk.release(idx)
    assert mk.trapped == len(idx) and mk.n == len(x0) - len(idx) and mk.skins.n == mk.n
    assert np.array_equal(x, x0[idx])
    assert np.allclose(mass, rows0[0][rows], rtol=1e-14, atol=0.0)
    assert np.allclose(H, rows0[1][rows], rtol=1e-14, atol=0.0)
    assert np.allclose(P, rows0[2][rows], rtol=1e-14, atol=0.0)
    M1, E1, P1 = _books(mk.skins)
    assert M1 + mass.sum() == pytest.approx(M0, rel=1e-14)
    assert E1 + H.sum() == pytest.approx(E0, rel=1e-14)
    assert np.allclose(P1 + P.sum(axis=0), P0, rtol=0, atol=1e-14 * M0)
    # the survivors keep their own skins
    keep = np.setdiff1d(np.arange(len(x0)), idx)
    assert np.array_equal(mk.x, x0[keep])
    assert np.array_equal(mk.skins.totals()[0][mk.skin_index], rows0[0][keep])


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


# -- the neighbours on meshes whose patches are much smaller than the spacing (Task 9b) -------------------------------
def test_the_tessellation_keeps_its_neighbours_on_a_sphere_at_design_size_cells():
    """Task 9b (the ruling on the neighbour test): a 10 mm sphere on the body's design cells (0.7 mm near the surface)
    with markers seeded as the body seeds them, one per 2 mm voxel. The patches are pieces of cut cells, about a sixth of
    the spacing across; the neighbour test bounds the path along the patch graph, not its hops, so the true neighbours
    are found: under 1 % of the cells are bounded by the cap, and every area is within a factor of 1.5 of the area
    nearest the marker. (Two hops of the patch graph, the rule before, capped most cells on the test cube and left areas
    0.38 to 2.67 times their nearest-marker areas.)"""
    R = 10e-3
    base = BoxBase.around((-R - 5.6e-3,) * 3, (R + 5.6e-3,) * 3, 2.8e-3)
    mesh = refined(base, lambda x: np.linalg.norm(x, axis=1) - R, DESIGN_FINE)
    phi = np.linalg.norm(mesh.points, axis=1) - R
    mk = seeded(mesh, phi, DESIGN_SPACING, base.lo)
    mk.tessellate(mesh, phi)
    r = nearest_ratio(mk, mesh, phi)
    print("design-size sphere: {} markers, {} capped; areas {:.3f} to {:.3f} times the nearest-marker areas".format(
        mk.n, int(mk.capped.sum()), r.min(), r.max()))
    assert mk.capped.mean() < CAPPED_SHARE
    assert r.min() >= 1.0 / OWNERSHIP_FACTOR and r.max() <= OWNERSHIP_FACTOR


def slot_body():
    """A block 7 x 7 x 5 mm with a slot 1.5 mm wide cut 5 mm deep into it from its +y face through its whole height, on
    0.25 mm cells near the surface; markers seeded at a 1 mm spacing (a reach of 3 mm, wider than the slot)."""
    base = BoxBase.around((-5e-3,) * 3, (5e-3,) * 3, 1e-3)

    def level(x):
        block = np.max(np.abs(x) - np.array([3.5e-3, 3.5e-3, 2.5e-3]), axis=1)
        cut = np.maximum(np.abs(x[:, 0]) - 0.75e-3, -1.5e-3 - x[:, 1])          # negative inside the slot
        return np.maximum(block, -cut)
    mesh = refined(base, level, 0.25e-3)
    phi = redistance(mesh, level(mesh.points))
    return mesh, phi, seeded(mesh, phi, 1e-3, base.lo)


SLOT_FLOOR = -1.5e-3       # m: the slot's floor (y); its walls rise from it to the block's +y face at 3.5 mm


def test_two_sheets_across_a_thin_gap_are_not_neighbours():
    """Task 9b (fix round 1): the slot's two walls face each other 1.5 mm apart, half the 3 mm reach, and meet only round
    its floor, so that low on the walls the path along the surface between them is within the path bound. Over the walls'
    whole height no marker on one wall is a neighbour of one on the other (their normals oppose and each lies in front
    of the other), though dozens of such pairs lie within the reach; every wall marker keeps a neighbour on its own wall;
    and the wall markers more than 1 mm above the floor have areas within a factor of 1.5 of the area nearest them. The
    lower ones, by the concave corner at the floor, are reported and not bounded. (With the path bound alone, 23 pairs
    crossed the gap, with wall markers 0.42 to 2.67 mm above the floor, and cut the lowest ones' areas to 0.137 of the
    nearest-marker area.)"""
    from reentry_model.levelset.markers import NEIGHBOUR_REACH
    mesh, phi, mk = slot_body()
    faces = mk.tessellate(mesh, phi)
    x, n = mk.x, mk.normal
    wa = (np.abs(x[:, 0] + 0.75e-3) < 0.1e-3) & (n[:, 0] > 0.9)
    wb = (np.abs(x[:, 0] - 0.75e-3) < 0.1e-3) & (n[:, 0] < -0.9)
    across = (wa[faces.i] & wb[faces.j]) | (wb[faces.i] & wa[faces.j])
    within = (np.linalg.norm(x[wa][:, None] - x[wb][None], axis=2) < NEIGHBOUR_REACH * mk.spacing).sum()
    own = np.zeros(mk.n, bool)
    for side in (wa, wb):
        same = side[faces.i] & side[faces.j]
        own[faces.i[same]] = own[faces.j[same]] = True
    r = nearest_ratio(mk, mesh, phi)
    height = x[:, 1] - SLOT_FLOOR
    wall = wa | wb
    high, low = wall & (height > 1e-3), wall & (height <= 1e-3)
    print("slot: walls of {} and {} markers {:.2f} to {:.2f} mm above the floor, {} pairs across within the reach, {} of "
          "them neighbours; areas over nearest-marker areas {:.3f} to {:.3f} more than 1 mm above the floor, {:.3f} to "
          "{:.3f} by the floor's corner ({} markers)".format(
              int(wa.sum()), int(wb.sum()), height[wall].min() * 1e3, height[wall].max() * 1e3, int(within),
              int(across.sum()), r[high].min(), r[high].max(), r[low].min(), r[low].max(), int(low.sum())))
    assert wa.sum() >= 10 and wb.sum() >= 10 and within >= 30 and height[wall].min() < 0.5e-3
    assert not across.any()
    assert own[wall].all()
    assert r[high].min() >= 1.0 / OWNERSHIP_FACTOR and r[high].max() <= OWNERSHIP_FACTOR


def test_a_concave_right_angle_corner_keeps_its_neighbours():
    """Task 9b (fix round 1): a step, a floor (normal +z) meeting a wall (normal +x) at a concave right angle, markers
    1 mm apart. Pairs across the corner lie each in front of the other, as the slot's walls do, but their normals are at
    about a right angle (measured 0.34 to 0.47 for the dot product: the recovered normals blend within a cell of the
    corner), far above OPPOSED_NORMALS (-0.5): the neighbour test keeps them, faces join the floor's and the wall's
    markers by the corner, and none of their cells is capped."""
    from reentry_model.levelset.markers import OPPOSED_NORMALS
    base = BoxBase.around((-4e-3, -3e-3, -3e-3), (4e-3, 3e-3, 4e-3), 1e-3)

    def level(x):
        floor = np.maximum(x[:, 2], -x[:, 2] - 2.5e-3)                   # the slab 0 > z > -2.5 mm
        step = np.maximum(np.maximum(x[:, 0], x[:, 2] - 3e-3), -x[:, 2] - 2.5e-3)   # and the block x < 0, z < 3 mm
        box = np.maximum(np.abs(x[:, 0]) - 3.5e-3, np.abs(x[:, 1]) - 2.5e-3)
        return np.maximum(np.minimum(floor, step), box)
    mesh = refined(base, lambda x: level(x + 0.13e-3), 0.25e-3)
    phi = redistance(mesh, level(mesh.points + 0.13e-3))
    mk = seeded(mesh, phi, 1e-3, base.lo)
    faces = mk.tessellate(mesh, phi)
    x, n = mk.x + 0.13e-3, mk.normal
    fl = (np.abs(x[:, 2]) < 0.1e-3) & (n[:, 2] > 0.9) & (x[:, 0] > 0.0) & (x[:, 0] < 1.2e-3)
    wl = (np.abs(x[:, 0]) < 0.1e-3) & (n[:, 0] > 0.9) & (x[:, 2] > 0.0) & (x[:, 2] < 1.2e-3)
    across = np.flatnonzero((fl[faces.i] & wl[faces.j]) | (wl[faces.i] & fl[faces.j]))
    i, j = faces.i[across], faces.j[across]
    d = mk.x[j] - mk.x[i]
    dot = np.einsum("ij,ij->i", n[i], n[j])
    front = (np.einsum("ij,ij->i", n[i], d) > 0.0) & (np.einsum("ij,ij->i", n[j], d) < 0.0)
    print("concave corner: {} floor and {} wall markers by it, {} faces across it (normals' dot products {}), {} "
          "capped".format(int(fl.sum()), int(wl.sum()), len(across), np.round(dot, 2), int(mk.capped[fl | wl].sum())))
    assert fl.sum() >= 3 and wl.sum() >= 3 and len(across) >= 3
    assert front.all() and dot.min() > OPPOSED_NORMALS + 0.5
    assert not mk.capped[fl | wl].any()


def test_patch_paths_match_an_unrestricted_dijkstra(monkeypatch):
    """Task 9b (fix round 1, minor 2): `patch_paths`, Dijkstra from its sources in spatial chunks each on the patches
    within the bound of its sources, gives every pair's path along the patch graph exactly as Dijkstra over the whole
    graph does, and np.inf beyond the bound: on a 10 mm sphere on 0.7 mm cells, 2000 random pairs of patches up to the
    sphere's diameter apart, a 6 mm bound, and chunks of 7 sources (so that most sources share a chunk with others)."""
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import dijkstra
    from reentry_model.levelset import markers as M
    from reentry_model.levelset.thin import patch_adjacency
    R = 10e-3
    base = BoxBase.around((-R - 5.6e-3,) * 3, (R + 5.6e-3,) * 3, 2.8e-3)
    mesh = refined(base, lambda x: np.linalg.norm(x, axis=1) - R, DESIGN_FINE)
    phi = np.linalg.norm(mesh.points, axis=1) - R
    P = patches(mesh, phi)
    rng = np.random.default_rng(12)
    a, b = rng.integers(0, len(P.cell), 2000), rng.integers(0, len(P.cell), 2000)
    limit = 6e-3
    monkeypatch.setattr(M, "PATH_CHUNK", 7)
    got = M.patch_paths(mesh, P, a, b, limit)
    A = patch_adjacency(mesh, P).tocoo()
    w = np.maximum(np.linalg.norm(P.centroid[A.row] - P.centroid[A.col], axis=1), 1e-300)
    src = np.unique(a)
    D = dijkstra(coo_matrix((w, (A.row, A.col)), shape=A.shape).tocsr(), indices=src, limit=limit)
    want = D[np.searchsorted(src, a), b]
    within = np.isfinite(want)
    print("patch paths: {} pairs, {} within the bound, largest difference {:.1e} m".format(
        len(a), int(within.sum()), np.abs(got[within] - want[within]).max()))
    assert 100 < within.sum() < len(a)
    assert np.array_equal(np.isfinite(got), within)
    np.testing.assert_allclose(got[within], want[within], rtol=1e-12, atol=0.0)


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


def test_a_merge_across_modes_keeps_the_top_zone_marker_s_mode_and_hold():
    """Stage 3 Task 10, fix round 1 (the controller's ruling (4)): two markers closer than a third of the spacing merge,
    their skins mixed by the hand-over rule; where one is top-zone and the other full-depth, the merged marker is
    top-zone with the top-zone marker's hold, whichever of the two the merge keeps (the step end then gives the volume
    what the mixed column holds below its fine zone); where both share a mode, the kept marker's mode and hold stay.
    Three pairs 0.1 mm apart on a flat row 1 mm apart: (full, top), (top, full) and (full, full)."""
    mat = ConstMaterial()
    x = np.array([[0, 0, 0], [1e-4, 0, 0], [3e-3, 0, 0], [3.1e-3, 0, 0], [6e-3, 0, 0], [6.1e-3, 0, 0]], float)
    n = len(x)
    sk = Skins(mat, np.full(n, 1e-6), np.zeros((n, 2)), SETTINGS)
    sk.fill(mat.rho * 1e-6 * 0.8e-3, 600.0)
    mk = Markers(x, np.tile([0.0, 0.0, 1.0], (n, 1)), np.full(n, 1e-6), np.zeros((n, 2)), np.arange(n),
                 np.full(n, -1), 1e-3, skins=sk)
    mk.top = np.array([False, True, True, False, False, False])
    mk.hold = np.array([2, 1, 0, 2, 1, 2])
    assert mk.merge(1e-3) == 3 and mk.n == 3
    order = np.argsort(mk.x[:, 0])
    np.testing.assert_array_equal(mk.top[order], [True, True, False])
    np.testing.assert_array_equal(mk.hold[order], [1, 0, 1])

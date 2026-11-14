import numpy as np
import pytest

from reentry_model.levelset.boxmesh import BoxBase, build
from reentry_model.levelset.cutcell import (Patches, cell_gradients, inside_subtets, inside_volumes, patches,
                                           subtet_volumes, subtets, volumes_below)

H0 = 1.0e-3


def test_inside_and_outside_parts_add_to_each_cell_for_random_fields():
    m = build(BoxBase.around((0, 0, 0), (3e-3, 3e-3, 3e-3), H0))
    phi = np.random.default_rng(3).normal(size=len(m.points))
    total = inside_volumes(m, phi) + inside_volumes(m, -phi)
    assert np.allclose(total, m.volumes, rtol=1e-12)


def test_a_plane_gives_the_exact_volume_area_and_normal():
    base = BoxBase.around((0, 0, 0), (4e-3, 3e-3, 2e-3), H0)
    m = build(base, lambda x, s: np.where(np.abs(x[:, 0] - 1.37e-3) < 0.6e-3, 0.25e-3, np.inf))
    n = np.array([0.3, 0.5, 0.81]); n /= np.linalg.norm(n)
    d = n @ np.array([2e-3, 1.5e-3, 1e-3])
    phi = m.points @ n - d
    vol = inside_volumes(m, phi).sum()
    P = patches(m, phi)
    assert np.allclose(P.normal, n, atol=1e-12)
    # exact answers from a fine plane slicing of the same box, independently by Monte-Carlo-free formula:
    flat = build(base)
    assert vol == pytest.approx(inside_volumes(flat, phi).sum(), rel=1e-12)
    assert P.area.sum() == pytest.approx(patches(flat, phi).area.sum(), rel=1e-12)


def test_a_sphere_volume_and_area_converge_at_second_order():
    R, c = 1.1e-3, np.array([1.5e-3, 1.5e-3, 1.5e-3])
    errs = []
    for h in (0.25e-3, 0.125e-3):
        m = build(BoxBase.around((0, 0, 0), (3e-3, 3e-3, 3e-3), h))
        phi = np.linalg.norm(m.points - c, axis=1) - R
        V = inside_volumes(m, phi).sum()
        A = patches(m, phi).area.sum()
        errs.append((abs(V / (4 / 3 * np.pi * R ** 3) - 1), abs(A / (4 * np.pi * R ** 2) - 1)))
    (v1, a1), (v2, a2) = errs
    assert v2 < 0.01 and a2 < 0.005          # measured 0.65 % and 0.33 % at R/h = 8.8
    assert v1 / v2 > 3.5 and a1 / a2 > 3.5          # measured 4.0 and 4.0 per halving


def test_patch_normals_point_out_of_the_metal_and_centroids_lie_on_the_surface():
    R, c = 1.1e-3, np.array([1.5e-3, 1.5e-3, 1.5e-3])
    m = build(BoxBase.around((0, 0, 0), (3e-3, 3e-3, 3e-3), 0.25e-3))
    phi = np.linalg.norm(m.points - c, axis=1) - R
    P = patches(m, phi)
    radial = (P.centroid - c) / np.linalg.norm(P.centroid - c, axis=1, keepdims=True)
    assert (np.einsum("ij,ij->i", P.normal, radial) > 0.9).all()
    assert np.abs(np.linalg.norm(P.centroid - c, axis=1) - R).max() < 0.25e-3 ** 2 / R


def test_a_contour_through_nodes_leaves_no_zero_area_patches():
    m = build(BoxBase.around((0, 0, 0), (2e-3, 2e-3, 2e-3), 0.25e-3))
    phi = m.points[:, 0] - 1.0e-3                    # the plane x = 1 mm passes through a layer of nodes
    P = patches(m, phi)
    assert (P.area > 0).all()
    assert P.area.sum() == pytest.approx(4e-6, rel=1e-12)


def test_subtets_reproduce_the_cut_volumes():
    L = 14e-3
    m = build(BoxBase.around((-L,) * 3, (L,) * 3, 2e-3))
    phi = np.linalg.norm(m.points, axis=1) - 10e-3
    nneg = (phi[m.tets] < 0).sum(axis=1)
    cut = np.flatnonzero((nneg > 0) & (nneg < 4))
    assert set(np.unique(nneg[cut])) == {1, 2, 3}
    S, owner = inside_subtets(m, phi, cut)
    v = np.bincount(owner, subtet_volumes(S), minlength=len(cut))
    assert np.allclose(v, inside_volumes(m, phi, cut), rtol=1e-12, atol=1e-24)


# The cut-cell geometry at 770259e, copied verbatim (plan 03b Task 2, fix round 1): `inside_volumes`, `inside_subtets`
# and `patches` now call the coordinate-generic `volumes_below`, `subtets` and `contour`, which the level set's move by
# volumes also uses, and must give exactly what they gave before.

def _old_tvol(a, b, c, d):
    return np.abs(np.einsum("ij,ij->i", b - a, np.cross(c - a, d - a))) / 6.0


def _old_tri_area(a, b, c):
    return 0.5 * np.linalg.norm(np.cross(b - a, c - a), axis=1)


def _old_inside_volumes(mesh, phi, cells=None):
    """Volume of {phi < 0} in each cell (all cells, or the given ones)."""
    cells = np.arange(len(mesh.tets)) if cells is None else np.asarray(cells)
    T = mesh.tets[cells]
    f = phi[T]
    P = mesh.points[T]
    neg = f < 0
    nneg = neg.sum(axis=1)
    V = _old_tvol(P[:, 0], P[:, 1], P[:, 2], P[:, 3])
    out = np.where(nneg == 4, V, 0.0)
    for want, flip in ((1, False), (3, True)):
        m = nneg == want
        if not m.any():
            continue
        ff, PP = f[m], P[m]
        lone = np.argmax(~neg[m] if flip else neg[m], axis=1)       # the single vertex on its own side
        r = np.arange(len(ff))
        fa = ff[r, lone]
        frac = np.ones(len(ff))
        for j in range(4):
            other = j != lone
            with np.errstate(divide="ignore", invalid="ignore"):
                t = fa / (fa - ff[:, j])
            frac = np.where(other, frac * t, frac)
        corner = V[m] * frac
        out[m] = (V[m] - corner) if flip else corner
    m = nneg == 2
    if m.any():
        ff, PP = f[m], P[m]
        order = np.argsort(~neg[m], axis=1, kind="stable")          # two negatives first
        r = np.arange(len(ff))[:, None]
        fo, Po = ff[r, order], PP[r, order]
        a, b, c, d = (Po[:, i] for i in range(4))
        fa, fb, fc, fd = (fo[:, i] for i in range(4))
        cut = lambda x, y, fx, fy: x + (fx / (fx - fy))[:, None] * (y - x)
        pac, pad, pbc, pbd = cut(a, c, fa, fc), cut(a, d, fa, fd), cut(b, c, fb, fc), cut(b, d, fb, fd)
        # prism bottom (a, pac, pad), top (b, pbc, pbd): three tetrahedra
        out[m] = _old_tvol(a, pac, pad, b) + _old_tvol(pac, pad, b, pbc) + _old_tvol(pad, b, pbc, pbd)
    return out


def _old_inside_subtets(mesh, phi, cells):
    """The metal part {phi < 0} of each given cut cell as sub-tetrahedra.

    Returns (S, owner): S is (n, 4, 3) vertices, owner (n,) indexes into `cells`.  One tetrahedron for a lone inside
    vertex, three for the 2-2 prism (the split used by `inside_volumes`) and three for the wedge left by a lone
    outside vertex.  Nodes with phi exactly zero count as outside, as in `inside_volumes`."""
    cells = np.asarray(cells)
    T = mesh.tets[cells]
    f, P = phi[T], mesh.points[T]
    neg = f < 0
    nneg = neg.sum(axis=1)
    cut = lambda x, y, fx, fy: x + (fx / (fx - fy))[:, None] * (y - x)
    rest = np.array([[j for j in range(4) if j != lone] for lone in range(4)])
    out, owner = [], []
    m = np.flatnonzero(nneg == 1)
    if len(m):
        lone = np.argmax(neg[m], axis=1)
        others = rest[lone]
        r = np.arange(len(m))
        a, fa = P[m][r, lone], f[m][r, lone]
        pts = [cut(a, P[m][r, others[:, i]], fa, f[m][r, others[:, i]]) for i in range(3)]
        out.append(np.stack([a] + pts, axis=1)); owner.append(m)
    m = np.flatnonzero(nneg == 3)
    if len(m):
        lone = np.argmax(~neg[m], axis=1)                           # the single outside vertex
        others = rest[lone]
        r = np.arange(len(m))
        d, fd = P[m][r, lone], f[m][r, lone]
        a, b, c = (P[m][r, others[:, i]] for i in range(3))
        fa, fb, fc = (f[m][r, others[:, i]] for i in range(3))
        pa, pb, pc = cut(a, d, fa, fd), cut(b, d, fb, fd), cut(c, d, fc, fd)
        for tet in ((a, b, c, pa), (b, c, pa, pb), (c, pa, pb, pc)):
            out.append(np.stack(tet, axis=1)); owner.append(m)
    m = np.flatnonzero(nneg == 2)
    if len(m):
        order = np.argsort(~neg[m], axis=1, kind="stable")          # two negatives first
        rr = np.arange(len(m))[:, None]
        Po, fo = P[m][rr, order], f[m][rr, order]
        a, b, c, d = (Po[:, i] for i in range(4))
        fa, fb, fc, fd = (fo[:, i] for i in range(4))
        pac, pad, pbc, pbd = cut(a, c, fa, fc), cut(a, d, fa, fd), cut(b, c, fb, fc), cut(b, d, fb, fd)
        for tet in ((a, pac, pad, b), (pac, pad, b, pbc), (pad, b, pbc, pbd)):
            out.append(np.stack(tet, axis=1)); owner.append(m)
    if not out:
        return np.zeros((0, 4, 3)), np.zeros(0, int)
    return np.concatenate(out), np.concatenate(owner)


def _old_patches(mesh, phi, cells=None):
    """The patches of every cut cell (or of the cut cells among `cells`)."""
    cells = np.arange(len(mesh.tets)) if cells is None else np.asarray(cells)
    nneg_all = (phi[mesh.tets[cells]] < 0).sum(axis=1)
    cut_cells = cells[(nneg_all > 0) & (nneg_all < 4)]
    T = mesh.tets[cut_cells]
    f, P = phi[T], mesh.points[T]
    neg = f < 0
    nneg = neg.sum(axis=1)
    g = cell_gradients(mesh, phi, cut_cells)
    normal = g / np.linalg.norm(g, axis=1, keepdims=True)
    tris, tri_patch = [], []
    cut = lambda x, y, fx, fy: x + (fx / (fx - fy))[:, None] * (y - x)
    for want in (1, 3):
        m = np.flatnonzero(nneg == want)
        if not len(m):
            continue
        lone = np.argmax(neg[m] if want == 1 else ~neg[m], axis=1)
        others = np.array([[j for j in range(4) if j != l] for l in range(4)])[lone]
        rr = np.arange(len(m))
        pa, fa = P[m][rr, lone], f[m][rr, lone]
        pts = [cut(pa, P[m][rr, others[:, i]], fa, f[m][rr, others[:, i]]) for i in range(3)]
        tris.append(np.stack(pts, axis=1))
        tri_patch.append(m)
    m = np.flatnonzero(nneg == 2)
    if len(m):
        order = np.argsort(~neg[m], axis=1, kind="stable")
        rr = np.arange(len(m))[:, None]
        Po, fo = P[m][rr, order], f[m][rr, order]
        a, b, c, d = (Po[:, i] for i in range(4))
        fa, fb, fc, fd = (fo[:, i] for i in range(4))
        pac, pbc, pbd, pad = cut(a, c, fa, fc), cut(b, c, fb, fc), cut(b, d, fb, fd), cut(a, d, fa, fd)
        tris.append(np.stack([pac, pbc, pbd], axis=1)); tri_patch.append(m)
        tris.append(np.stack([pac, pbd, pad], axis=1)); tri_patch.append(m)
    tris = np.concatenate(tris) if tris else np.zeros((0, 3, 3))
    tri_patch = np.concatenate(tri_patch) if tri_patch else np.zeros(0, int)
    ta = _old_tri_area(tris[:, 0], tris[:, 1], tris[:, 2])
    area = np.bincount(tri_patch, ta, minlength=len(cut_cells))
    cen = np.stack([np.bincount(tri_patch, ta * tris[:, :, k].mean(axis=1), minlength=len(cut_cells))
                    for k in range(3)], axis=1) / np.maximum(area, 1e-300)[:, None]
    # a contour through a node or along an edge leaves zero-area pieces; their neighbours carry that surface
    keep = area > 0
    renum = np.cumsum(keep) - 1
    tk = keep[tri_patch]
    return Patches(cut_cells[keep], area[keep], cen[keep], normal[keep], tris[tk], renum[tri_patch[tk]])



def _fields(m):
    """Fields on the mesh m that cut its cells every way: a sphere's distance, and a random field with a tenth of its
    nodes exactly zero (and some at -0.0), so that corners at exactly zero count as outside."""
    rng = np.random.default_rng(11)
    sphere = np.linalg.norm(m.points - np.array([1.4e-3, 1.6e-3, 1.5e-3]), axis=1) - 1.1e-3
    noisy = rng.normal(size=len(m.points))
    noisy[rng.random(len(m.points)) < 0.1] = 0.0
    noisy[rng.random(len(m.points)) < 0.02] = -0.0
    return sphere, noisy


def test_the_generic_cut_cell_geometry_gives_exactly_what_it_gave_before():
    """inside_volumes, inside_subtets and patches, now built on volumes_below, subtets and contour, are bit for bit
    the functions they were, on every cell and on a subset, for fields with corners at exactly zero; and subtets with
    `whole` keeps the wholly negative tetrahedra too, its parts' volumes summing to volumes_below."""
    m = build(BoxBase.around((0, 0, 0), (3e-3, 3e-3, 3e-3), 0.5e-3),
              lambda x, s: np.where(np.abs(x[:, 0] - 1.37e-3) < 0.6e-3, 0.125e-3, np.inf))
    some = np.random.default_rng(5).choice(len(m.tets), size=len(m.tets) // 3, replace=False)
    sphere, noisy = _fields(m)
    assert (noisy == 0.0).sum() > len(noisy) // 20 and np.signbit(noisy[noisy == 0.0]).any()
    for phi in (sphere, noisy):
        assert {1, 2, 3} <= set(np.unique((phi[m.tets] < 0).sum(axis=1)).tolist())
        for cells in (None, some):
            assert np.array_equal(inside_volumes(m, phi, cells), _old_inside_volumes(m, phi, cells))
            new, old = patches(m, phi, cells), _old_patches(m, phi, cells)
            for field in ("cell", "area", "centroid", "normal", "tris", "tri_patch"):
                assert np.array_equal(getattr(new, field), getattr(old, field)), field
        nneg = (phi[m.tets] < 0).sum(axis=1)
        cut = np.flatnonzero((nneg > 0) & (nneg < 4))
        for cells in (cut, cut[::3], cut[:0]):
            S, owner = inside_subtets(m, phi, cells)
            S0, owner0 = _old_inside_subtets(m, phi, cells)
            assert np.array_equal(S, S0) and np.array_equal(owner, owner0) and owner.dtype == owner0.dtype
        T = m.tets
        S, owner = subtets(m.points[T], phi[T], whole=True)
        v = np.bincount(owner, subtet_volumes(S), minlength=len(T))
        assert np.allclose(v, volumes_below(m.points[T], phi[T]), rtol=1e-12, atol=1e-24)

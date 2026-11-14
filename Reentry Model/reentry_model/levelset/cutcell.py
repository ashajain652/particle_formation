"""Cut cells of a P1 level set on tetrahedra (design §7.2): the volume of {phi < 0} in every cell and, in every cell
the zero contour crosses, its flat piece (the patch): area, centroid and outward normal grad(phi)/|grad(phi)|.

A node with phi exactly zero counts as outside (phi treated as +0), so every cell is cut consistently.

The constructions are written once, for tetrahedra given by their corners and a linear function's values there
(`volumes_below`, `subtets`, `contour`), and the mesh's cells use them (`inside_volumes`, `inside_subtets`,
`patches`). The level set's move by volumes (`levelset.move_by_volumes`) applies the same ones to the parts of cells it
divides among the markers, so that its groups' volumes are the books' to round-off (plan 03b Task 2, fix round 1)."""
from dataclasses import dataclass

import numpy as np


@dataclass
class Patches:
    cell: np.ndarray       # (P,) index of the cut cell
    area: np.ndarray       # (P,) m2
    centroid: np.ndarray   # (P, 3)
    normal: np.ndarray     # (P, 3) unit, pointing out of the metal (toward phi > 0)
    tris: np.ndarray       # (Q, 3, 3) the patch triangles (a quad is two), for distances and output
    tri_patch: np.ndarray  # (Q,) patch index of each triangle


def _tvol(a, b, c, d):
    return np.abs(np.einsum("ij,ij->i", b - a, np.cross(c - a, d - a))) / 6.0


def _tri_area(a, b, c):
    return 0.5 * np.linalg.norm(np.cross(b - a, c - a), axis=1)


_REST = np.array([[j for j in range(4) if j != lone] for lone in range(4)])   # a tetrahedron's corners but each one


def _cut(x, y, fx, fy):
    """The points between x and y (k, D) where the linear function with the values fx and fy (k,) there is zero."""
    return x + (fx / (fx - fy))[:, None] * (y - x)


def cell_gradients(mesh, phi, cells=None):
    """Gradient of the P1 field in each cell (all, or the given ones), from the cached inverse edge matrices."""
    from .transfer import _geom
    _, Jinv, _ = _geom(mesh)
    T = mesh.tets if cells is None else mesh.tets[cells]
    Ji = Jinv if cells is None else Jinv[cells]
    dphi = phi[T[:, 1:]] - phi[T[:, :1]]
    return np.einsum("nji,nj->ni", Ji, dphi)


def volumes_below(P, f):
    """The volume of {f < 0} in each tetrahedron with corners P (m, 4, 3), f linear on it with the values f (m, 4) at the
    corners (a corner at exactly zero outside): the whole, the corner cut off at a lone corner (or the whole less it),
    and the 2-2 prism as three tetrahedra."""
    neg = f < 0
    nneg = neg.sum(axis=1)
    V = _tvol(P[:, 0], P[:, 1], P[:, 2], P[:, 3])
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
        pac, pad, pbc, pbd = _cut(a, c, fa, fc), _cut(a, d, fa, fd), _cut(b, c, fb, fc), _cut(b, d, fb, fd)
        # prism bottom (a, pac, pad), top (b, pbc, pbd): three tetrahedra
        out[m] = _tvol(a, pac, pad, b) + _tvol(pac, pad, b, pbc) + _tvol(pad, b, pbc, pbd)
    return out


def inside_volumes(mesh, phi, cells=None):
    """Volume of {phi < 0} in each cell (all cells, or the given ones)."""
    cells = np.arange(len(mesh.tets)) if cells is None else np.asarray(cells)
    T = mesh.tets[cells]
    return volumes_below(mesh.points[T], phi[T])


def subtet_volumes(S):
    """Volumes of sub-tetrahedra given as an (n, 4, 3) array of vertices."""
    return _tvol(S[:, 0], S[:, 1], S[:, 2], S[:, 3])


def subtets(S, f, whole=False):
    """The part {f < 0} of each tetrahedron with corners S (m, 4, D) in any coordinates (physical, or barycentric in a
    cell), f linear on it with the values f (m, 4) at the corners, as tetrahedra: (corners (k, 4, D), the index (k,) of
    the tetrahedron each came from). One for a lone inside corner, three for the 2-2 prism (the split `volumes_below`
    uses) and three for the wedge a lone outside corner leaves; a corner at exactly zero counts as outside. The wholly
    negative tetrahedra are left out, or with `whole` returned whole after the rest."""
    neg = f < 0
    nneg = neg.sum(axis=1)
    out, owner = [], []
    m = np.flatnonzero(nneg == 1)
    if len(m):
        lone = np.argmax(neg[m], axis=1)
        others = _REST[lone]
        r = np.arange(len(m))
        a, fa = S[m][r, lone], f[m][r, lone]
        pts = [_cut(a, S[m][r, others[:, i]], fa, f[m][r, others[:, i]]) for i in range(3)]
        out.append(np.stack([a] + pts, axis=1)); owner.append(m)
    m = np.flatnonzero(nneg == 3)
    if len(m):
        lone = np.argmax(~neg[m], axis=1)                           # the single outside vertex
        others = _REST[lone]
        r = np.arange(len(m))
        d, fd = S[m][r, lone], f[m][r, lone]
        a, b, c = (S[m][r, others[:, i]] for i in range(3))
        fa, fb, fc = (f[m][r, others[:, i]] for i in range(3))
        pa, pb, pc = _cut(a, d, fa, fd), _cut(b, d, fb, fd), _cut(c, d, fc, fd)
        for tet in ((a, b, c, pa), (b, c, pa, pb), (c, pa, pb, pc)):
            out.append(np.stack(tet, axis=1)); owner.append(m)
    m = np.flatnonzero(nneg == 2)
    if len(m):
        order = np.argsort(~neg[m], axis=1, kind="stable")          # two negatives first
        rr = np.arange(len(m))[:, None]
        Po, fo = S[m][rr, order], f[m][rr, order]
        a, b, c, d = (Po[:, i] for i in range(4))
        fa, fb, fc, fd = (fo[:, i] for i in range(4))
        pac, pad, pbc, pbd = _cut(a, c, fa, fc), _cut(a, d, fa, fd), _cut(b, c, fb, fc), _cut(b, d, fb, fd)
        for tet in ((a, pac, pad, b), (pac, pad, b, pbc), (pad, b, pbc, pbd)):
            out.append(np.stack(tet, axis=1)); owner.append(m)
    if whole:
        m = np.flatnonzero(nneg == 4)
        out.append(S[m]); owner.append(m)
    if not out:
        return np.zeros((0, 4, S.shape[2])), np.zeros(0, int)
    return np.concatenate(out), np.concatenate(owner)


def inside_subtets(mesh, phi, cells):
    """The metal part {phi < 0} of each given cut cell as sub-tetrahedra (`subtets`).

    Returns (S, owner): S is (n, 4, 3) vertices, owner (n,) indexes into `cells`.  One tetrahedron for a lone inside
    vertex, three for the 2-2 prism (the split used by `inside_volumes`) and three for the wedge left by a lone
    outside vertex.  Nodes with phi exactly zero count as outside, as in `inside_volumes`."""
    cells = np.asarray(cells)
    T = mesh.tets[cells]
    return subtets(mesh.points[T], phi[T])


def contour(S, f):
    """The zero contour of the linear function with the values f (m, 4) at the corners S (m, 4, D) of tetrahedra, in
    any coordinates, as triangles: (corners (Q, 3, D), the index (Q,) of the tetrahedron each lies in). One for a lone
    corner on its side, two for the 2-2 quad; a corner at exactly zero counts as outside."""
    neg = f < 0
    nneg = neg.sum(axis=1)
    tris, owner = [], []
    for want in (1, 3):
        m = np.flatnonzero(nneg == want)
        if not len(m):
            continue
        lone = np.argmax(neg[m] if want == 1 else ~neg[m], axis=1)
        others = _REST[lone]
        rr = np.arange(len(m))
        pa, fa = S[m][rr, lone], f[m][rr, lone]
        pts = [_cut(pa, S[m][rr, others[:, i]], fa, f[m][rr, others[:, i]]) for i in range(3)]
        tris.append(np.stack(pts, axis=1))
        owner.append(m)
    m = np.flatnonzero(nneg == 2)
    if len(m):
        order = np.argsort(~neg[m], axis=1, kind="stable")
        rr = np.arange(len(m))[:, None]
        Po, fo = S[m][rr, order], f[m][rr, order]
        a, b, c, d = (Po[:, i] for i in range(4))
        fa, fb, fc, fd = (fo[:, i] for i in range(4))
        pac, pbc, pbd, pad = _cut(a, c, fa, fc), _cut(b, c, fb, fc), _cut(b, d, fb, fd), _cut(a, d, fa, fd)
        tris.append(np.stack([pac, pbc, pbd], axis=1)); owner.append(m)
        tris.append(np.stack([pac, pbd, pad], axis=1)); owner.append(m)
    tris = np.concatenate(tris) if tris else np.zeros((0, 3, S.shape[2]))
    owner = np.concatenate(owner) if owner else np.zeros(0, int)
    return tris, owner


def patches(mesh, phi, cells=None):
    """The patches of every cut cell (or of the cut cells among `cells`)."""
    cells = np.arange(len(mesh.tets)) if cells is None else np.asarray(cells)
    nneg_all = (phi[mesh.tets[cells]] < 0).sum(axis=1)
    cut_cells = cells[(nneg_all > 0) & (nneg_all < 4)]
    T = mesh.tets[cut_cells]
    g = cell_gradients(mesh, phi, cut_cells)
    normal = g / np.linalg.norm(g, axis=1, keepdims=True)
    tris, tri_patch = contour(mesh.points[T], phi[T])
    ta = _tri_area(tris[:, 0], tris[:, 1], tris[:, 2])
    area = np.bincount(tri_patch, ta, minlength=len(cut_cells))
    cen = np.stack([np.bincount(tri_patch, ta * tris[:, :, k].mean(axis=1), minlength=len(cut_cells))
                    for k in range(3)], axis=1) / np.maximum(area, 1e-300)[:, None]
    # a contour through a node or along an edge leaves zero-area pieces; their neighbours carry that surface
    keep = area > 0
    renum = np.cumsum(keep) - 1
    tk = keep[tri_patch]
    return Patches(cut_cells[keep], area[keep], cen[keep], normal[keep], tris[tk], renum[tri_patch[tk]])

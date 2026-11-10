"""Cut cells of a P1 level set on tetrahedra (design §7.2): the volume of {phi < 0} in every cell and, in every cell
the zero contour crosses, its flat piece (the patch): area, centroid and outward normal grad(phi)/|grad(phi)|.

A node with phi exactly zero counts as outside (phi treated as +0), so every cell is cut consistently."""
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


def cell_gradients(mesh, phi, cells=None):
    """Gradient of the P1 field in each cell (all, or the given ones), from the cached inverse edge matrices."""
    from .transfer import _geom
    _, Jinv, _ = _geom(mesh)
    T = mesh.tets if cells is None else mesh.tets[cells]
    Ji = Jinv if cells is None else Jinv[cells]
    dphi = phi[T[:, 1:]] - phi[T[:, :1]]
    return np.einsum("nji,nj->ni", Ji, dphi)


def inside_volumes(mesh, phi, cells=None):
    """Volume of {phi < 0} in each cell (all cells, or the given ones)."""
    cells = np.arange(len(mesh.tets)) if cells is None else np.asarray(cells)
    T = mesh.tets[cells]
    f = phi[T]
    P = mesh.points[T]
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
        cut = lambda x, y, fx, fy: x + (fx / (fx - fy))[:, None] * (y - x)
        pac, pad, pbc, pbd = cut(a, c, fa, fc), cut(a, d, fa, fd), cut(b, c, fb, fc), cut(b, d, fb, fd)
        # prism bottom (a, pac, pad), top (b, pbc, pbd): three tetrahedra
        out[m] = _tvol(a, pac, pad, b) + _tvol(pac, pad, b, pbc) + _tvol(pad, b, pbc, pbd)
    return out


def patches(mesh, phi, cells=None):
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
    ta = _tri_area(tris[:, 0], tris[:, 1], tris[:, 2])
    area = np.bincount(tri_patch, ta, minlength=len(cut_cells))
    cen = np.stack([np.bincount(tri_patch, ta * tris[:, :, k].mean(axis=1), minlength=len(cut_cells))
                    for k in range(3)], axis=1) / np.maximum(area, 1e-300)[:, None]
    # a contour through a node or along an edge leaves zero-area pieces; their neighbours carry that surface
    keep = area > 0
    renum = np.cumsum(keep) - 1
    tk = keep[tri_patch]
    return Patches(cut_cells[keep], area[keep], cen[keep], normal[keep], tris[tk], renum[tri_patch[tk]])

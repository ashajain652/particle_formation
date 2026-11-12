"""Adaptive tetrahedral box mesh for the level-set architecture (design §6.1 as amended in §20.5).

The box is a regular grid of cubes of edge h0, each split into six Kuhn (Freudenthal) tetrahedra along its main
diagonal. Refinement is Maubach's newest-vertex bisection with tags, which keeps the mesh conforming (no hanging
nodes) and the tetrahedra within three similarity classes. Coarsening is done by rebuilding from the base grid with
the current targets, so every mesh is a deterministic function of (base, targets).

A leaf at bisection level L has the cube-edge size h0 / 2**(L // 3); three bisections halve it. A target size s is
met when that size is at most s, so the achieved size is the largest h0 / 2**k not above s.
"""
from dataclasses import dataclass
from itertools import permutations

import numpy as np


@dataclass(frozen=True)
class BoxBase:
    lo: tuple          # lower corner (m)
    h0: float          # cube edge (m)
    n: tuple           # cubes along x, y, z

    @classmethod
    def around(cls, lo, hi, h0):
        lo, hi = np.asarray(lo, float), np.asarray(hi, float)
        n = np.maximum(1, np.ceil((hi - lo) / h0 - 1e-9)).astype(int)
        return cls(tuple(lo), float(h0), tuple(int(v) for v in n))

    @property
    def hi(self):
        return tuple(np.asarray(self.lo) + self.h0 * np.asarray(self.n))


class AdaptiveMesh:
    """Leaves of the bisection forest. tets rows are Maubach-ordered (x0, x1, x2, x3); tag in 1..3; level from 0."""

    def __init__(self, base, points, tets, tag, level):
        self.base, self.points, self.tets, self.tag, self.level = base, points, tets, tag, level
        self._cache = {}

    # --- geometry -------------------------------------------------------------------------------------------
    @property
    def size(self):
        return self.base.h0 * 2.0 ** (-(self.level // 3))

    @property
    def centroids(self):
        return self.points[self.tets].mean(axis=1)

    @property
    def volumes(self):
        if "vol" not in self._cache:
            p = self.points[self.tets]
            self._cache["vol"] = np.abs(np.einsum("ij,ij->i", p[:, 1] - p[:, 0],
                                                  np.cross(p[:, 2] - p[:, 0], p[:, 3] - p[:, 0]))) / 6.0
        return self._cache["vol"]

    def node_cells(self):
        """CSR matrix nodes x cells (1 where the node is a vertex of the cell), cached."""
        if "nc" not in self._cache:
            from scipy.sparse import coo_matrix
            m = len(self.tets)
            self._cache["nc"] = coo_matrix((np.ones(4 * m, np.int8), (self.tets.ravel(), np.repeat(np.arange(m), 4))),
                                           shape=(len(self.points), m)).tocsr()
        return self._cache["nc"]

    def ring(self, cells, rings):
        """The cells within `rings` rings (through shared nodes) of the given cells, sorted."""
        cells = np.unique(np.asarray(cells))
        nc = self.node_cells()
        for _ in range(rings):
            nodes = np.unique(self.tets[cells])
            cells = np.unique(nc[nodes].indices)
        return cells

    def to_skfem(self):
        from skfem import MeshTet
        return MeshTet(self.points.T.copy(), self.tets.T.copy())


def _edge_key(a, b, n_max):
    lo, hi = np.minimum(a, b), np.maximum(a, b)
    return lo.astype(np.int64) * np.int64(n_max) + hi.astype(np.int64)


def base_mesh(base):
    nx, ny, nz = base.n
    ii, jj, kk = np.meshgrid(np.arange(nx + 1), np.arange(ny + 1), np.arange(nz + 1), indexing="ij")
    pts = np.column_stack([ii.ravel(), jj.ravel(), kk.ravel()]).astype(float) * base.h0 + np.asarray(base.lo)
    nid = lambda i, j, k: (i * (ny + 1) + j) * (nz + 1) + k
    ci, cj, ck = [a.ravel() for a in np.meshgrid(np.arange(nx), np.arange(ny), np.arange(nz), indexing="ij")]
    tets = []
    for perm in permutations(range(3)):
        o = np.zeros((len(ci), 3), int)
        verts = [nid(ci, cj, ck)]
        for ax in perm:
            o[:, ax] += 1
            verts.append(nid(ci + o[:, 0], cj + o[:, 1], ck + o[:, 2]))
        tets.append(np.column_stack(verts))
    tets = np.vstack(tets).astype(np.int64)
    m = len(tets)
    return AdaptiveMesh(base, pts, tets, np.full(m, 3, np.int8), np.zeros(m, np.int16))


class _Builder:
    def __init__(self, mesh):
        self.pts = [mesh.points]
        self.npts = len(mesh.points)
        self.tets, self.tag, self.level = mesh.tets, mesh.tag, mesh.level
        self.fresh = np.ones(len(mesh.tets), bool)     # leaves whose target has not been evaluated yet
        self.nmax = 1 << 31
        self.mid_keys = np.empty(0, np.int64)
        self.mid_ids = np.empty(0, np.int64)

    def _points(self):
        if len(self.pts) > 1:
            self.pts = [np.vstack(self.pts)]
        return self.pts[0]

    def _midpoints(self, a, b):
        keys = _edge_key(a, b, self.nmax)
        out = np.empty(len(keys), np.int64)
        pos = np.searchsorted(self.mid_keys, keys)
        pos_c = np.minimum(pos, max(len(self.mid_keys) - 1, 0))
        known = (len(self.mid_keys) > 0) & (self.mid_keys[pos_c] == keys) if len(self.mid_keys) else np.zeros(len(keys), bool)
        out[known] = self.mid_ids[pos_c[known]]
        new_keys, inv = np.unique(keys[~known], return_inverse=True)
        if len(new_keys):
            pts = self._points()
            lo = new_keys // self.nmax
            hi = new_keys % self.nmax
            new_ids = self.npts + np.arange(len(new_keys))
            self.pts.append(0.5 * (pts[lo] + pts[hi]))
            self.npts += len(new_keys)
            out[~known] = new_ids[inv]
            keys_all = np.concatenate([self.mid_keys, new_keys])
            ids_all = np.concatenate([self.mid_ids, new_ids])
            order = np.argsort(keys_all, kind="stable")
            self.mid_keys, self.mid_ids = keys_all[order], ids_all[order]
        return out

    def bisect(self, sel):
        T, k = self.tets[sel], self.tag[sel].astype(int)
        r = np.arange(len(T))
        z = self._midpoints(T[:, 0], T[r, k])
        c1, c2 = np.empty_like(T), np.empty_like(T)
        for kk in (1, 2, 3):
            m = k == kk
            if not m.any():
                continue
            t = T[m]
            zz = z[m]
            # child1 = (x0..x_{k-1}, z, x_{k+1}..x3); child2 = (x1..x_k, z, x_{k+1}..x3)
            a = [t[:, j] for j in range(kk)] + [zz] + [t[:, j] for j in range(kk + 1, 4)]
            b = [t[:, j] for j in range(1, kk + 1)] + [zz] + [t[:, j] for j in range(kk + 1, 4)]
            c1[m], c2[m] = np.column_stack(a), np.column_stack(b)
        newtag = np.where(k > 1, k - 1, 3).astype(np.int8)
        newlev = (self.level[sel] + 1).astype(np.int16)
        keep = np.ones(len(self.tets), bool)
        keep[sel] = False
        self.fresh = np.concatenate([self.fresh[keep], np.ones(2 * len(sel), bool)])
        self.tets = np.vstack([self.tets[keep], c1, c2])
        self.tag = np.concatenate([self.tag[keep], newtag, newtag])
        self.level = np.concatenate([self.level[keep], newlev, newlev])

    def hanging(self):
        if len(self.mid_keys) == 0:
            return np.zeros(len(self.tets), bool)
        T = self.tets
        bad = np.zeros(len(T), bool)
        for i, j in ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)):
            keys = _edge_key(T[:, i], T[:, j], self.nmax)
            pos = np.minimum(np.searchsorted(self.mid_keys, keys), len(self.mid_keys) - 1)
            bad |= self.mid_keys[pos] == keys
        return bad


def build(base, target_size=None, max_level=30, max_passes=400):
    """Refine the base grid until every leaf meets target_size(centroids, sizes) -> target sizes (m).

    target_size is called on the leaves created by the previous pass only (a leaf that met its target keeps meeting
    it); it returns, per leaf, the largest allowed cube-edge size (np.inf for none). Closure bisections keep the mesh
    conforming. Returns an AdaptiveMesh."""
    b = _Builder(base_mesh(base))
    for _ in range(max_passes):
        pts = b._points()
        mark = np.zeros(len(b.tets), bool)
        if target_size is not None and b.fresh.any():
            f = np.flatnonzero(b.fresh)
            size = base.h0 * 2.0 ** (-(b.level[f] // 3))
            want = np.asarray(target_size(pts[b.tets[f]].mean(axis=1), size), float)
            with np.errstate(divide="ignore"):
                need = np.where(want < base.h0, 3 * np.ceil(np.log2(base.h0 / np.maximum(want, 1e-300)) - 1e-12), 0)
            mark[f] = b.level[f] < np.minimum(need, max_level)
            b.fresh[:] = False
        mark |= b.hanging()
        if not mark.any():
            break
        b.bisect(np.flatnonzero(mark))
    else:
        raise RuntimeError("adaptive build did not converge in {} passes".format(max_passes))
    pts = b._points()
    return AdaptiveMesh(base, pts.copy(), b.tets.copy(), b.tag.copy(), b.level.copy())

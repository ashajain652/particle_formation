"""Geometry of a finite-element frame: point location, the thickness map and centre crossing (M1 plan, Task 5).

Spec §6.3 (fields on particles), §9.1 (removal by centre crossing) and §9.2 trigger 3 (a patch thinner than a few
particle spacings), on the arrays of `contract.Frame` (finite-element numbering, or the compact numbering of a
prepared frame; both work, nothing here depends on unreferenced nodes).

- `TetLocator`: which active tetrahedron contains a point, and its barycentric weights; a uniform grid of buckets
  over the tetrahedra's bounding boxes. The body is the union of the active tetrahedra (an element with phi < 1 is
  located like any other: the frame's geometry does not shrink with phi, plan fact 6). Ties on a shared face go to
  the lowest tetrahedron index; a point within `SURFACE_TOL_M` of the surface counts as inside.
- `centre_crossed`: a particle centre no longer inside the next frame's body (spec §9.1).
- `thickness_map`: from each patch centroid, a ray along the inward normal to its first intersection with any other
  surface triangle. The triangle test (Moller-Trumbore, no back-face culling) does not depend on the winding of the
  triangles it hits, and a hit on an edge or vertex shared by several triangles (non-manifold surfaces after element
  death, plan fact 7) gives the same first distance from any of them. **The inward direction is `-normals`**: the
  caller passes outward normals (Task 4's `derived_patch_arrays`, after its orientation check). A patch whose normal
  points outward-wrong sends its ray out of the body; on a closed body that ray usually escapes, giving NaN, which
  `thickness_summary` counts.
- `thin_patches`: the trigger: thickness below `n_spacings` particle spacings.

The module imports numpy and the standard library only. `thickness_map(method="vtk")` imports vtk inside the
function (prepare side only); the brute force is numpy and is the default up to `BRUTEFORCE_MAX_FACES` faces."""
from __future__ import annotations

import math

import numpy as np

SURFACE_TOL_M = 1e-12          # m: points this close outside the body count as inside (plan Task 5, Step 1)
BRUTEFORCE_MAX_FACES = 5000    # thickness_map(method="auto") uses the brute force up to this many faces
_RAY_UV_TOL = 1e-12            # barycentric slack of the ray-triangle test: a ray through an edge hits both sides
_VTK_TOL = 1e-9               # vtk cell-locator candidate tolerance, relative to the bounding-box diagonal
_PAIRS_PER_CHUNK = 2_000_000   # point-tetrahedron pairs (locator) or ray-triangle pairs (brute force) per chunk


# ------------------------------------------------------------------------------------------------------- locator
class TetLocator:
    """Point location in the active tetrahedra `tets` (m, 4) on `points` (n, 3).

    Buckets: a uniform grid of cubic cells of side `cell_size` (default: the bounding box volume per tetrahedron,
    cube-rooted, so about one tetrahedron per cell on average); each tetrahedron is listed in every cell its
    bounding box (widened by `SURFACE_TOL_M`) touches, in increasing tetrahedron index, so the first containing
    candidate of a point is the lowest-index one."""

    def __init__(self, points, tets, cell_size=None):
        self.points = np.ascontiguousarray(points, dtype=np.float64)
        self.tets = np.ascontiguousarray(tets, dtype=np.int64)
        if self.tets.ndim != 2 or self.tets.shape[1] != 4 or len(self.tets) == 0:
            raise ValueError("tets must be a non-empty (m, 4) array")
        x = self.points[self.tets]                                   # (m, 4, 3)
        self._x0 = x[:, 0].copy()
        J = np.stack([x[:, 1] - x[:, 0], x[:, 2] - x[:, 0], x[:, 3] - x[:, 0]], axis=2)   # columns = edges
        self._Jinv = np.linalg.inv(J)                                # rows: grad of lambda_1..3
        grad = np.concatenate([-self._Jinv.sum(axis=1, keepdims=True), self._Jinv], axis=1)   # (m, 4, 3)
        # lambda_i = (signed distance from face i) * |grad lambda_i|: a distance slack of SURFACE_TOL_M
        self._tol = SURFACE_TOL_M * np.linalg.norm(grad, axis=2)    # (m, 4)

        lo = x.min(axis=1) - 2.0 * SURFACE_TOL_M
        hi = x.max(axis=1) + 2.0 * SURFACE_TOL_M
        self._lo = lo.min(axis=0)
        extent = hi.max(axis=0) - self._lo
        if cell_size is None:
            cell_size = (float(np.prod(np.maximum(extent, 1e-300))) / len(self.tets)) ** (1.0 / 3.0)
        self.cell_size = float(cell_size)
        self._dims = np.maximum(np.floor(extent / self.cell_size).astype(np.int64) + 1, 1)
        i0 = self._index(lo)
        i1 = self._index(hi)
        span = i1 - i0 + 1                                           # (m, 3) cells per axis
        count = span.prod(axis=1)
        tet_of_pair = np.repeat(np.arange(len(self.tets), dtype=np.int64), count)
        local = np.arange(int(count.sum()), dtype=np.int64) - np.repeat(np.cumsum(count) - count, count)
        sp = span[tet_of_pair]
        ix = i0[tet_of_pair, 0] + local % sp[:, 0]
        iy = i0[tet_of_pair, 1] + (local // sp[:, 0]) % sp[:, 1]
        iz = i0[tet_of_pair, 2] + local // (sp[:, 0] * sp[:, 1])
        cell = (iz * self._dims[1] + iy) * self._dims[0] + ix
        order = np.lexsort((tet_of_pair, cell))                      # by cell, then by tetrahedron index
        self._items = tet_of_pair[order]
        n_cells = int(self._dims.prod())
        self._start = np.zeros(n_cells + 1, dtype=np.int64)
        np.cumsum(np.bincount(cell, minlength=n_cells), out=self._start[1:])

    @classmethod
    def from_frame(cls, frame, cell_size=None) -> "TetLocator":
        return cls(frame.points, frame.tets, cell_size)

    def _index(self, x):
        return np.floor((x - self._lo) / self.cell_size).astype(np.int64)

    def _cells(self, x):
        """Linear cell index of each point, -1 outside the grid (then outside every tetrahedron)."""
        i = self._index(x)
        ok = ((i >= 0) & (i < self._dims)).all(axis=1)
        i = np.where(ok[:, None], i, 0)
        return np.where(ok, (i[:, 2] * self._dims[1] + i[:, 1]) * self._dims[0] + i[:, 0], -1)

    def _barycentric(self, x, t):
        l123 = np.einsum("pij,pj->pi", self._Jinv[t], x - self._x0[t])
        return np.concatenate([1.0 - l123.sum(axis=1, keepdims=True), l123], axis=1)

    def locate(self, x) -> tuple[np.ndarray, np.ndarray]:
        """(tetrahedron index or -1, barycentric weights (n, 4); NaN where outside) for points x (n, 3)."""
        x = np.ascontiguousarray(np.asarray(x, dtype=np.float64).reshape(-1, 3))
        n = len(x)
        found = np.full(n, -1, dtype=np.int64)
        weights = np.full((n, 4), np.nan)
        cell = self._cells(x)
        pts = np.nonzero(cell >= 0)[0]
        start = self._start[cell[pts]]
        count = self._start[cell[pts] + 1] - start
        keep = count > 0
        pts, start, count = pts[keep], start[keep], count[keep]
        csum = np.cumsum(count)
        lo = 0
        while lo < len(pts):                                         # chunks of at most ~_PAIRS_PER_CHUNK pairs
            base = csum[lo - 1] if lo else 0
            hi = max(int(np.searchsorted(csum, base + _PAIRS_PER_CHUNK, side="right")), lo + 1)
            p, s, c = pts[lo:hi], start[lo:hi], count[lo:hi]
            pair_pt = np.repeat(p, c)
            offs = np.arange(int(c.sum()), dtype=np.int64) - np.repeat(np.cumsum(c) - c, c)
            pair_tet = self._items[np.repeat(s, c) + offs]
            lam = self._barycentric(x[pair_pt], pair_tet)
            inside = (lam >= -self._tol[pair_tet]).all(axis=1)
            # candidates of a point are contiguous and in increasing tetrahedron index: the first inside wins
            hit_pt, first = np.unique(pair_pt[inside], return_index=True)
            sel = np.nonzero(inside)[0][first]
            found[hit_pt] = pair_tet[sel]
            weights[hit_pt] = lam[sel]
            lo = hi
        return found, weights

    def contains(self, x) -> np.ndarray:
        return self.locate(x)[0] >= 0

    def interpolate(self, x, nodal) -> np.ndarray:
        """P1 interpolation of the nodal field(s) `nodal` (n_nodes,) or (n_nodes, k) at x; NaN outside."""
        idx, w = self.locate(x)
        nodal = np.asarray(nodal, dtype=np.float64)
        out = np.full((len(idx),) + nodal.shape[1:], np.nan)
        ok = idx >= 0
        vals = nodal[self.tets[idx[ok]]]                             # (p, 4, ...)
        out[ok] = np.einsum("pi,pi...->p...", w[ok], vals)
        return out


def centre_crossed(locator_next, x) -> np.ndarray:
    """True where a particle centre x (n, 3) is no longer inside the next frame's body (spec §9.1)."""
    return locator_next.locate(x)[0] < 0


# ------------------------------------------------------------------------------------------------------ thickness
def _ray_hits(origins, dirs, v0, e1, e2, uv_tol=_RAY_UV_TOL):
    """Moller-Trumbore for every (ray, triangle) pair: the ray parameter t (rays (r, 3), triangles (f, 3) -> (r, f)),
    NaN where it misses or the ray lies in the triangle's plane. No back-face culling: winding does not matter."""
    P = np.cross(dirs[:, None, :], e2[None, :, :])                   # (r, f, 3)
    det = np.einsum("fk,rfk->rf", e1, P)
    scale = np.linalg.norm(e1, axis=1) * np.linalg.norm(e2, axis=1)
    ok = np.abs(det) > 1e-12 * scale[None, :]
    inv = np.where(ok, 1.0 / np.where(ok, det, 1.0), 0.0)
    Tv = origins[:, None, :] - v0[None, :, :]
    u = np.einsum("rfk,rfk->rf", Tv, P) * inv
    Q = np.cross(Tv, e1[None, :, :])
    v = np.einsum("rk,rfk->rf", dirs, Q) * inv
    t = np.einsum("fk,rfk->rf", e2, Q) * inv
    hit = ok & (u >= -uv_tol) & (v >= -uv_tol) & (u + v <= 1.0 + uv_tol)
    return np.where(hit, t, np.nan)


def _rays(points, faces, normals, centroids, eps):
    points = np.asarray(points, dtype=np.float64)
    faces = np.asarray(faces, dtype=np.int64)
    n = np.asarray(normals, dtype=np.float64)
    n = n / np.linalg.norm(n, axis=1, keepdims=True)
    c = np.asarray(centroids, dtype=np.float64)
    v0 = points[faces[:, 0]]
    return c - eps * n, -n, v0, points[faces[:, 1]] - v0, points[faces[:, 2]] - v0


def _thickness_bruteforce(points, faces, normals, centroids, eps=1e-7) -> np.ndarray:
    """Every ray against every triangle (numpy; O(f^2), for meshes under about 5,000 faces and the tests)."""
    origins, dirs, v0, e1, e2 = _rays(points, faces, normals, centroids, eps)
    f = len(origins)
    out = np.full(f, np.nan)
    block = max(1, _PAIRS_PER_CHUNK // max(f, 1) // 4)
    for lo in range(0, f, block):
        hi = min(f, lo + block)
        t = _ray_hits(origins[lo:hi], dirs[lo:hi], v0, e1, e2)
        t[np.arange(hi - lo), np.arange(lo, hi)] = np.nan            # not the patch's own triangle
        t[~(t > 0.0)] = np.nan                                        # in front of the start point only
        good = ~np.isnan(t).all(axis=1)
        out[lo:hi][good] = np.nanmin(t[good], axis=1) + eps           # distance from the centroid
    return out


def _thickness_vtk(points, faces, normals, centroids, eps=1e-7) -> np.ndarray:
    """vtk's `vtkStaticCellLocator` lists the candidate triangles along each ray (a segment twice the bounding-box
    diagonal long); the distance is then recomputed with the brute force's own Moller-Trumbore test, so the two agree
    to round-off and the locator's tolerance only widens the candidate list. Prepare side only (imports vtk).

    Not `vtkOBBTree` as the plan first said: in drama_env's VTK 9.7.1 its `IntersectWithLine(p0, p1, points, ids)`
    reports no intersection at all, and the overload with a tolerance returns cell lists that omit the cell actually
    hit and vary between calls (measured 2026-10-07 on the sphere fixture); the static cell locator reports exactly
    the cell hit."""
    from vtkmodules.util.numpy_support import numpy_to_vtk, numpy_to_vtkIdTypeArray
    from vtkmodules.vtkCommonCore import vtkIdList, vtkPoints
    from vtkmodules.vtkCommonDataModel import vtkCellArray, vtkPolyData, vtkStaticCellLocator

    points = np.ascontiguousarray(points, dtype=np.float64)
    faces = np.ascontiguousarray(faces, dtype=np.int64)
    origins, dirs, v0, e1, e2 = _rays(points, faces, normals, centroids, eps)
    vpts = vtkPoints()
    vpts.SetData(numpy_to_vtk(points, deep=True))
    cells = vtkCellArray()
    offsets = np.arange(0, 3 * len(faces) + 1, 3, dtype=np.int64)
    cells.SetData(numpy_to_vtkIdTypeArray(offsets, deep=True),
                  numpy_to_vtkIdTypeArray(faces.ravel().copy(), deep=True))
    poly = vtkPolyData()
    poly.SetPoints(vpts)
    poly.SetPolys(cells)
    locator = vtkStaticCellLocator()
    locator.SetDataSet(poly)
    locator.BuildLocator()
    diag = float(np.linalg.norm(points.max(axis=0) - points.min(axis=0)))
    length = 2.0 * diag
    tol = _VTK_TOL * diag
    hit_pts, ids = vtkPoints(), vtkIdList()
    out = np.full(len(faces), np.nan)
    for i in range(len(faces)):
        o = origins[i]
        locator.IntersectWithLine(o.tolist(), (o + length * dirs[i]).tolist(), tol, hit_pts, ids)
        cand = np.array([ids.GetId(j) for j in range(ids.GetNumberOfIds())], dtype=np.int64)
        cand = cand[cand != i]
        if len(cand) == 0:
            continue
        t = _ray_hits(o[None, :], dirs[i][None, :], v0[cand], e1[cand], e2[cand])[0]
        t = t[t > 0.0]
        if len(t):
            out[i] = t.min() + eps
    return out


def thickness_map(points, faces, normals, centroids, eps=1e-7, method="auto") -> np.ndarray:
    """Thickness per patch, m (spec §9.2 trigger 3): the distance from the patch centroid, along the inward normal
    `-normals`, to the first other surface triangle (the ray starts `eps` inside, the distance is from the centroid).
    NaN where the ray escapes (an open surface, or a patch whose normal points inward). `normals` must be outward;
    the hit test itself is independent of the winding of the triangles hit.

    method: "bruteforce" (numpy, O(f^2)), "vtk" (vtkStaticCellLocator candidates, prepare side only), or "auto"
    (the brute force up to BRUTEFORCE_MAX_FACES faces, vtk above). Both give the same distances to round-off."""
    if method == "auto":
        method = "bruteforce" if len(faces) <= BRUTEFORCE_MAX_FACES else "vtk"
    if method == "bruteforce":
        return _thickness_bruteforce(points, faces, normals, centroids, eps)
    if method == "vtk":
        return _thickness_vtk(points, faces, normals, centroids, eps)
    raise ValueError("method must be 'auto', 'bruteforce' or 'vtk', got {!r}".format(method))


def thin_patches(thickness, dx, n_spacings=4.0) -> np.ndarray:
    """True where the patch is thinner than `n_spacings` particle spacings `dx` (spec §9.2 trigger 3). A NaN
    thickness (escaped ray) is not flagged; `thickness_summary` counts it."""
    with np.errstate(invalid="ignore"):
        return np.asarray(thickness, dtype=np.float64) < n_spacings * float(dx)


def thickness_summary(thickness) -> dict:
    """The prepare JSON's per-frame thickness entry: thinnest finite patch and the number of escaped rays."""
    t = np.asarray(thickness, dtype=np.float64)
    finite = t[np.isfinite(t)]
    return {"min_m": float(finite.min()) if len(finite) else math.nan,
            "n_nan": int(np.count_nonzero(~np.isfinite(t)))}

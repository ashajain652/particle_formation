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
from dataclasses import dataclass

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


# ---------------------------------------------------------------------------------------------- depths and zones
# Task 6 (spec §10, §9.2 trigger 2; check "Zones").
#
# Thresholds on the liquid fraction. Spec §10 defines the sprayable skin by the liquidus ("liquid, above 908 K") and
# the slurry depth by "the depth where the material becomes more than half solid", i.e. f_l = 0.5 (895.1 K on the
# Scheil curve; spec §2). Both are thresholds on f_l, so they hold for any material: the finite element's own nodal
# f_l (the frame's `liquid_fraction`) is marched by default, and Task 3's `MaterialTable.liquid_fraction` can be
# passed instead as a callable of the nodal temperature (`layer_depths(..., f_l=table.liquid_fraction)`).
SLURRY_FL = 0.5                # f_l above this: more than half liquid (spec §10 "more than half solid" below it)
LIQUID_FL = 1.0 - 1e-9         # f_l at or above this: liquid (spec §10 "above 908 K"; plan Task 6, Step 1)
DEPTH_DS_M = 5e-5              # m, the march's sample spacing: a fifth of the 0.25 mm outer cell (plan Task 6)
_PROBE_M = 1e-9                # m: the orientation probe's depth below the centroid
_SAMPLES_PER_BATCH = 64        # samples evaluated per marched ray and pass

ZONE_INVALID, ZONE_NONE, ZONE_SKIN, ZONE_RUNOFF, ZONE_BULK = -1, 0, 1, 2, 3
FILM_LIMIT_M = 2e-3            # m, spec §10's default film limit (bracket 3 mm; naming.BRACKET_DEFAULTS film_limit_mm)


def _crossing(s0, s1, f0, f1, thr):
    """Depth where the linear interpolant between samples (s0, f0) and (s1, f1) reaches thr (f0 above, f1 not)."""
    with np.errstate(invalid="ignore", divide="ignore"):
        w = np.where(f0 != f1, (f0 - thr) / (f0 - f1), 1.0)
    return s0 + np.clip(w, 0.0, 1.0) * (s1 - s0)


def layer_depths(frame, locator, normals, centroids, ds=DEPTH_DS_M, max_depth=None, f_l=None,
                 slurry_fl=SLURRY_FL, liquid_fl=LIQUID_FL, return_info=False):
    """(liquid_depth, slurry_depth) in m per patch of `frame`, measured in the mesh along each patch's inward normal
    (spec §10: "casting rays inward to the depth where the material becomes more than half solid").

    **The frame must be wound outward and `normals` must be its outward unit normals** -- that is,
    `frames.orient_outward(frame)` then `frames.derived_patch_arrays(...)["normal"]`. The real Step 3 export winds
    some triangles inward (about 18 % on the prototype's frames, Task 5), so the file's own winding must never be
    used here. A patch whose normal points out of the body is detected (a probe `_PROBE_M` below the centroid lies
    outside every tetrahedron) and gets NaN in both depths, counted in `info["n_wrong_way"]`.

    liquid_depth  contiguous f_l >= liquid_fl from the surface;
    slurry_depth  contiguous f_l > slurry_fl from the surface.
    The surface value is the P1 value at the centroid (the mean of the face's three nodes). A patch whose surface
    value fails a threshold has zero depth at once; only patches above `slurry_fl` are marched, so the cost scales
    with the molten area. The ray is sampled every `ds` from the centroid; the field at each sample is the P1
    interpolant of the nodal f_l in the active tetrahedra (`locator`, a `TetLocator` of the same frame), and the
    crossing is placed by linear interpolation between the last sample that passes and the first that fails --
    exact wherever f_l is linear along the ray between them, otherwise within ds. A sample outside the body ends
    the layer at the last inside sample (the ray left the body: the whole remaining thickness is layer).

    max_depth: None, a scalar or a per-patch array (the thickness map, Task 5; NaN = no cap): the march stops there
    and both depths are capped at it. Without it the march stops when the ray leaves the body or after the bounding
    box's diagonal.
    f_l: None (the frame's nodal `f_l`), a nodal array, or a callable of the nodal temperature `frame.node["T"]`
    (Task 3's `MaterialTable.liquid_fraction`). The thresholds are arguments so that a material's own definition
    of "liquid" or "half solid" can be passed.

    With return_info=True, also a dict: n_marched, n_wrong_way, n_samples, and left_body (bool per patch) /
    n_left_body: rays whose slurry layer reached the body's far side or max_depth instead of a crossing.
    Numpy only (uses `TetLocator`); prepare side by role."""
    if not liquid_fl > slurry_fl:
        raise ValueError("liquid_fl must exceed slurry_fl")
    faces = np.asarray(frame.faces, dtype=np.int64)
    nf = len(faces)
    if f_l is None:
        nodal = np.asarray(frame.node["f_l"], dtype=np.float64)
    elif callable(f_l):
        nodal = np.asarray(f_l(np.asarray(frame.node["T"], dtype=np.float64)), dtype=np.float64)
    else:
        nodal = np.asarray(f_l, dtype=np.float64)
    if nodal.shape != (len(frame.points),):
        raise ValueError("nodal f_l must have one value per node of the frame")
    n = np.asarray(normals, dtype=np.float64)
    n = n / np.linalg.norm(n, axis=1, keepdims=True)
    c = np.asarray(centroids, dtype=np.float64)
    if max_depth is None:
        cap = np.full(nf, np.inf)
    else:
        cap = np.broadcast_to(np.asarray(max_depth, dtype=np.float64), (nf,)).copy()
        cap[~np.isfinite(cap)] = np.inf
    span = float(np.linalg.norm(frame.points.max(axis=0) - frame.points.min(axis=0)))
    cap = np.minimum(cap, span)

    f_surf = nodal[faces].mean(axis=1)
    liquid = np.zeros(nf)
    slurry = np.zeros(nf)
    info = {"n_marched": 0, "n_wrong_way": 0, "n_left_body": 0, "n_samples": 0, "left_body": np.zeros(nf, bool)}
    march = np.flatnonzero(f_surf > slurry_fl)
    if len(march):
        probe = locator.locate(c[march] - _PROBE_M * n[march])[0]
        wrong = march[probe < 0]
        liquid[wrong] = np.nan
        slurry[wrong] = np.nan
        info["n_wrong_way"] = int(len(wrong))
        march = march[probe >= 0]
    info["n_marched"] = int(len(march))

    # state per marched ray: the last sample's depth and value, whether each layer is still open
    s_prev = np.zeros(len(march))
    f_prev = f_surf[march].copy()
    liq_open = f_prev >= liquid_fl
    slu_open = np.ones(len(march), bool)
    liq_d = np.zeros(len(march))
    slu_d = np.zeros(len(march))
    j0 = 1
    active = np.arange(len(march))
    while len(active):
        rays = march[active]
        steps = np.arange(j0, j0 + _SAMPLES_PER_BATCH, dtype=np.float64) * ds          # (b,)
        s = np.minimum(steps[None, :], cap[rays][:, None])                              # (a, b): last at the cap
        x = c[rays][:, None, :] - s[:, :, None] * n[rays][:, None, :]
        vals = locator.interpolate(x.reshape(-1, 3), nodal).reshape(s.shape)            # NaN outside the body
        info["n_samples"] += int(vals.size)
        sp, fp = s_prev[active].copy(), f_prev[active].copy()
        lo_, so_ = liq_open[active].copy(), slu_open[active].copy()
        ld, sd = liq_d[active].copy(), slu_d[active].copy()
        left = np.zeros(len(active), bool)
        for b in range(s.shape[1]):
            sb, fb = s[:, b], vals[:, b]
            outside = np.isnan(fb)
            # liquid layer
            close = lo_ & ~outside & ~(fb >= liquid_fl)
            ld[close] = _crossing(sp[close], sb[close], fp[close], fb[close], liquid_fl)
            stop = lo_ & outside
            ld[stop] = sp[stop]
            lo_ &= ~(close | stop)
            # slurry layer
            close = so_ & ~outside & ~(fb > slurry_fl)
            sd[close] = _crossing(sp[close], sb[close], fp[close], fb[close], slurry_fl)
            stop = so_ & outside
            sd[stop] = sp[stop]
            left |= stop
            so_ &= ~(close | stop)
            # the cap: a layer still open at the cap ends there
            at_cap = so_ & (sb >= cap[rays])
            ld[at_cap & lo_] = sb[at_cap & lo_]
            lo_ &= ~at_cap
            sd[at_cap] = sb[at_cap]
            left |= at_cap
            so_ &= ~at_cap
            upd = ~outside
            sp[upd], fp[upd] = sb[upd], fb[upd]
        s_prev[active], f_prev[active] = sp, fp
        liq_open[active], slu_open[active] = lo_, so_
        liq_d[active], slu_d[active] = ld, sd
        info["left_body"][rays[left]] = True
        active = active[so_]                       # the liquid layer is inside the slurry layer: it is closed too
        j0 += _SAMPLES_PER_BATCH
    info["n_left_body"] = int(np.count_nonzero(info["left_body"]))
    liquid[march] = np.minimum(liq_d, cap[march])
    slurry[march] = np.minimum(slu_d, cap[march])
    return (liquid, slurry, info) if return_info else (liquid, slurry)


@dataclass
class ZoneMap:
    """Spec §10's three zones per patch (all lengths in m). Decision 6 of the M1 review: for a layer thicker than the
    film limit both readings are reported -- the whole layer is Spheral's (`zone == ZONE_BULK`, thickness `layer`,
    resolution `under_resolved`), or only the part below the film limit (`bulk`, resolution
    `under_resolved_below_limit`); the choice is made before M4."""
    layer: np.ndarray          # film_thickness + slurry_depth (+ deep_thickness when include_deep)
    skin: np.ndarray           # min(film + liquid_depth, delta_m) where delta_m is finite, else NaN
    runoff: np.ndarray         # max(min(layer, film_limit) - skin, 0), skin taken as 0 where undefined
    bulk: np.ndarray           # max(layer - film_limit, 0)
    zone: np.ndarray           # int8: the deepest zone the layer reaches (ZONE_NONE where layer = 0, ZONE_INVALID NaN)
    under_resolved: np.ndarray             # bool: zone 3 and layer < n_spacings dx (whole-layer reading)
    under_resolved_below_limit: np.ndarray  # bool: zone 3 and bulk < n_spacings dx (below-the-limit reading)
    film_limit: float
    include_deep: bool


def classify_zones(film_thickness, liquid_depth, slurry_depth, delta_m, film_limit=FILM_LIMIT_M, dx=None,
                   n_spacings=4.0, deep_thickness=None, include_deep=False) -> ZoneMap:
    """Spec §10's rule as the M1 plan (Task 6, Step 2) states it, per patch:

    zone 1 (sprayable skin)  layer <= delta_m;
    zone 2 (thin runoff)     delta_m < layer <= film_limit (where delta_m is NaN: 0 < layer <= film_limit);
    zone 3 (bulk, Spheral)   layer > film_limit.
    The layer is the film plus the contiguous more-than-half-liquid depth (spec §9.2 trigger 2's "region more than
    half liquid"). Note that zone 1 is decided on that layer, while `skin` is the liquid part (film + liquid depth)
    within delta_m (spec §10: "liquid (above 908 K) within the conjugate depth"): a slurry layer thinner than delta_m
    is zone 1 with a skin of only its film (the inconsistency the plan flags; both quantities are reported).

    `deep_thickness` (a mass per area, m_d / (rho_l A)) is excluded by default: at the 0.5 s step it is a time-step
    artefact (spec §2, §6.1; decision 6). `include_deep=True` adds it to the layer and changes nothing else.
    `dx` (particle spacing, m): flags zone-3 patches thinner than n_spacings dx under both readings; None flags
    nothing. NaN in the film, depths or (with include_deep) deep_thickness gives ZONE_INVALID and NaN lengths."""
    film = np.asarray(film_thickness, dtype=np.float64)
    liq = np.asarray(liquid_depth, dtype=np.float64)
    slu = np.asarray(slurry_depth, dtype=np.float64)
    dm = np.asarray(delta_m, dtype=np.float64)
    layer = film + slu
    if include_deep:
        if deep_thickness is None:
            raise ValueError("include_deep=True needs deep_thickness")
        layer = layer + np.asarray(deep_thickness, dtype=np.float64)
    limit = float(film_limit)
    girin = np.isfinite(dm)
    valid = np.isfinite(layer) & np.isfinite(liq)
    skin = np.where(girin & valid, np.minimum(film + liq, dm), np.nan)
    with np.errstate(invalid="ignore"):
        runoff = np.maximum(np.minimum(layer, limit) - np.where(girin, skin, 0.0), 0.0)
        bulk = np.maximum(layer - limit, 0.0)
        zone = np.full(layer.shape, ZONE_NONE, dtype=np.int8)
        zone[layer > 0.0] = ZONE_RUNOFF
        zone[girin & (layer > 0.0) & (layer <= dm)] = ZONE_SKIN
        zone[layer > limit] = ZONE_BULK
    zone[~valid] = ZONE_INVALID
    runoff[~valid] = np.nan
    bulk[~valid] = np.nan
    if dx is None:
        under = np.zeros(layer.shape, bool)
        under_b = np.zeros(layer.shape, bool)
    else:
        thin = n_spacings * float(dx)
        is3 = zone == ZONE_BULK
        under = is3 & (layer < thin)
        under_b = is3 & (bulk < thin)
    return ZoneMap(layer=layer, skin=skin, runoff=runoff, bulk=bulk, zone=zone, under_resolved=under,
                   under_resolved_below_limit=under_b, film_limit=limit, include_deep=bool(include_deep))


def bulk_slurry_area(zones, area) -> float:
    """m^2 of zone-3 patches (spec §9.2 trigger 2's input, per frame). The same under both of decision 6's readings:
    a patch is zone 3 exactly where `bulk` > 0."""
    a = np.asarray(area, dtype=np.float64)
    return float(math.fsum(a[zones.zone == ZONE_BULK]))


def zone_summary(zones, area) -> dict:
    """Per-frame numbers for the prepare JSON: zone-3 area, the largest layer, the patch count per zone, and the
    zone-3 volume under both readings (whole layer: sum A layer; below the limit: sum A bulk)."""
    a = np.asarray(area, dtype=np.float64)
    is3 = zones.zone == ZONE_BULK
    finite = zones.layer[np.isfinite(zones.layer)]
    return {"film_limit_m": zones.film_limit, "include_deep": zones.include_deep,
            "bulk_area_m2": bulk_slurry_area(zones, a),
            "max_layer_m": float(finite.max()) if len(finite) else math.nan,
            "n_patches": {name: int(np.count_nonzero(zones.zone == z)) for name, z in
                          (("invalid", ZONE_INVALID), ("none", ZONE_NONE), ("skin", ZONE_SKIN),
                           ("runoff", ZONE_RUNOFF), ("bulk", ZONE_BULK))},
            "bulk_volume_whole_layer_m3": float(math.fsum(a[is3] * zones.layer[is3])),
            "bulk_volume_below_limit_m3": float(math.fsum(a[is3] * zones.bulk[is3])),
            "n_under_resolved": int(np.count_nonzero(zones.under_resolved)),
            "n_under_resolved_below_limit": int(np.count_nonzero(zones.under_resolved_below_limit))}

"""The skins' upstream-ray shelter test: which top-surface corners are sheltered from the oncoming gas.

From each corner of the top surface (the skins' gas boundary) a ray starts at the corner moved 0.05 sqrt(area) [m] out
along its surface normal (a twentieth of the corner's cell width, so the ray clears the corner's own facets) and runs
along `+v_hat`, up the flow (`v_hat` is the direction of motion, so the gas arrives from `+v_hat`). A corner is sheltered
when that ray crosses any triangle of the top surface at a positive distance along `v_hat`. Edges and vertices of a
triangle count as crossings, and the corner's own triangles are not skipped, so a corner on a face that leans over itself
is sheltered by it. This is the design's rule 2 (`docs/superpowers/specs/2026-11-13-fitted-sphere-loads-design.md`) and
mirrors the level-set shadow test in `reentry_model/levelset/shadow.py`.

`upstream_hits` finds the triangles a ray could cross with a uniform grid on the plane across `v_hat`, cells about one
median projected edge length [m] wide, each triangle listed in every cell its projected bounding box overlaps.
`upstream_hits_plain` tests every ray against every triangle's bounding box and is kept as the reference; both use the
same pair arithmetic, so they agree exactly.
"""
import numpy as np

OFFSET_FRACTION = 0.05


def _basis(v):
    """Unit v and an orthonormal pair (e1, e2) spanning the plane across it."""
    v = np.asarray(v, float)
    v = v / np.linalg.norm(v)
    e1 = np.cross(v, [0, 0, 1.0]) if abs(v[2]) < 0.9 else np.cross(v, [0, 1.0, 0])
    e1 /= np.linalg.norm(e1)
    return v, e1, np.cross(v, e1)


def _setup(points, faces, normals, areas, v_hat, candidates):
    v, e1, e2 = _basis(v_hat)
    n = len(points)
    cand = np.ones(n, bool) if candidates is None else np.asarray(candidates, bool)
    O = points + normals * (OFFSET_FRACTION * np.sqrt(areas))[:, None]
    o2 = np.stack([O @ e1, O @ e2], 1)
    ox = O @ v
    T = points[faces]
    t2 = np.stack([T @ e1, T @ e2], -1)
    tx = T @ v
    return cand, o2, ox, t2, tx


def _cross(q, mi, fi, lo, hi, t2, tx, ox):
    """For ray/triangle pairs (ray index mi into q and ox, triangle index fi): True where the point is inside the
    triangle's bounding box, inside the triangle (inclusive) and the hit is further along v_hat than the ray's start."""
    p = q[mi]
    keep = ((p >= lo[fi]) & (p <= hi[fi])).all(1)
    mi, fi, p = mi[keep], fi[keep], p[keep]
    a, b, c = t2[fi, 0], t2[fi, 1], t2[fi, 2]
    den = (b[:, 1] - c[:, 1]) * (a[:, 0] - c[:, 0]) + (c[:, 0] - b[:, 0]) * (a[:, 1] - c[:, 1])
    ok = np.abs(den) > 1e-30
    d = np.where(ok, den, 1)
    l1 = np.where(ok, ((b[:, 1] - c[:, 1]) * (p[:, 0] - c[:, 0]) + (c[:, 0] - b[:, 0]) * (p[:, 1] - c[:, 1])) / d, -1)
    l2 = np.where(ok, ((c[:, 1] - a[:, 1]) * (p[:, 0] - c[:, 0]) + (a[:, 0] - c[:, 0]) * (p[:, 1] - c[:, 1])) / d, -1)
    l3 = 1 - l1 - l2
    inside = ok & (l1 >= 0) & (l2 >= 0) & (l3 >= 0)
    xhit = l1 * tx[fi, 0] + l2 * tx[fi, 1] + l3 * tx[fi, 2]
    return mi[inside & (xhit > ox[mi])]


def upstream_hits_plain(points, faces, normals, areas, v_hat, candidates=None):
    """bool (n,): True where the corner's upstream ray crosses a top-surface triangle. Tests every candidate against every
    triangle's bounding box in chunks; the reference for `upstream_hits`."""
    cand, o2, ox, t2, tx = _setup(points, faces, normals, areas, v_hat, candidates)
    lo, hi = t2.min(1), t2.max(1)
    idx = np.nonzero(cand)[0]
    hit = np.zeros(len(points), bool)
    nf = len(faces)
    for s in range(0, len(idx), 128):
        sel = idx[s:s + 128]
        q = o2[sel]
        inbox = ((q[:, None, :] >= lo[None]) & (q[:, None, :] <= hi[None])).all(-1)
        mi, fi = np.nonzero(inbox)
        hit[sel[np.unique(_cross(q, mi, fi, lo, hi, t2, tx, ox[sel]))]] = True
    return hit


def upstream_hits(points, faces, normals, areas, v_hat, candidates=None):
    """bool (n,): as `upstream_hits_plain`, found with a grid on the plane across `v_hat`."""
    cand, o2, ox, t2, tx = _setup(points, faces, normals, areas, v_hat, candidates)
    hit = np.zeros(len(points), bool)
    idx = np.nonzero(cand)[0]
    if len(idx) == 0:
        return hit
    lo, hi = t2.min(1), t2.max(1)
    edges = np.concatenate([np.linalg.norm(t2[:, 1] - t2[:, 0], axis=1), np.linalg.norm(t2[:, 2] - t2[:, 1], axis=1),
                            np.linalg.norm(t2[:, 0] - t2[:, 2], axis=1)])
    h = float(np.median(edges[edges > 0])) if (edges > 0).any() else 1.0
    org = lo.min(0)
    nc = np.maximum(np.floor((hi.max(0) - org) / h).astype(int) + 1, 1)
    c_lo = np.floor((lo - org) / h).astype(int)
    c_hi = np.floor((hi - org) / h).astype(int)
    # list each triangle in every cell its projected bounding box overlaps (CSR: cell -> triangles)
    span = c_hi - c_lo + 1
    tri = np.repeat(np.arange(len(faces)), span[:, 0] * span[:, 1])
    k = np.arange(len(tri)) - np.repeat(np.cumsum(span[:, 0] * span[:, 1]) - span[:, 0] * span[:, 1], span[:, 0] * span[:, 1])
    cx = c_lo[tri, 0] + k % span[tri, 0]
    cy = c_lo[tri, 1] + k // span[tri, 0]
    cell = cx * nc[1] + cy
    order = np.argsort(cell, kind="stable")
    tri_sorted = tri[order]
    start = np.searchsorted(cell[order], np.arange(nc[0] * nc[1] + 1))
    q = o2[idx]
    qc = np.clip(np.floor((q - org) / h).astype(int), 0, nc - 1)
    qcell = qc[:, 0] * nc[1] + qc[:, 1]
    cnt = start[qcell + 1] - start[qcell]
    mi = np.repeat(np.arange(len(idx)), cnt)
    off = np.arange(len(mi)) - np.repeat(np.cumsum(cnt) - cnt, cnt)
    fi = tri_sorted[start[qcell][mi] + off]
    hit[idx[np.unique(_cross(q, mi, fi, lo, hi, t2, tx, ox[idx]))]] = True
    return hit

"""Operations on a P1 level set (negative inside the metal) on an AdaptiveMesh: reset to a signed distance, motion
under a prescribed velocity, connected pieces, and the exact volume correction per piece (design §6.2-§6.4; §20.3:
the evolved level set is the skin base beta)."""
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree

from .cutcell import inside_volumes, patches
from .transfer import interpolate


def _segment_distance(p, a, b):
    ab = b - a
    t = np.einsum("ij,ij->i", p - a, ab) / np.maximum(np.einsum("ij,ij->i", ab, ab), 1e-300)
    return np.linalg.norm(p - (a + np.clip(t, 0.0, 1.0)[:, None] * ab), axis=1)


def point_triangle_distance(p, a, b, c):
    """Distance from points p (n,3) to triangles (a, b, c) (n,3 each). The nearest of the three edges, or the
    perpendicular distance where the foot falls inside the triangle; degenerate (zero-area) triangles reduce to
    their edges, so a contour passing through a node gives no NaN."""
    d = np.minimum(np.minimum(_segment_distance(p, a, b), _segment_distance(p, b, c)), _segment_distance(p, c, a))
    n = np.cross(b - a, c - a)
    nn = np.einsum("ij,ij->i", n, n)
    ok = nn > 1e-30 * np.maximum(np.einsum("ij,ij->i", b - a, b - a), 1e-300) ** 2
    with np.errstate(divide="ignore", invalid="ignore"):
        h = np.einsum("ij,ij->i", p - a, n) / nn
        foot = p - h[:, None] * n
        u = np.einsum("ij,ij->i", np.cross(b - foot, c - foot), n) / nn
        v = np.einsum("ij,ij->i", np.cross(c - foot, a - foot), n) / nn
    inside = ok & (u >= 0) & (v >= 0) & (u + v <= 1)
    d[inside] = np.minimum(d[inside], np.abs(h[inside]) * np.sqrt(nn[inside]))
    return d


def redistance(mesh, phi, k=8, band=None):
    """Signed distance to the zero contour of phi, keeping phi's sign. Nodal values are exact distances to the
    piecewise-flat contour, searched among the k triangles with the nearest centroids. With `band` (m), only nodes
    within band of a triangle centroid are computed; the rest get +-band, never more than their true distance, so
    sphere tracing on the result stays safe."""
    P = patches(mesh, phi)
    if len(P.tris) == 0:
        return np.where(phi < 0, -np.inf, np.inf)
    tri = P.tris
    kd = cKDTree(tri.mean(axis=1))
    kk = min(k, len(tri))
    d = np.full(len(mesh.points), np.inf if band is None else float(band))
    if band is None:
        sel = np.arange(len(mesh.points))
        _, idx = kd.query(mesh.points, k=kk)
    else:
        near, _ = kd.query(mesh.points, k=1, distance_upper_bound=band)
        sel = np.flatnonzero(np.isfinite(near))
        _, idx = kd.query(mesh.points[sel], k=kk)
    idx = idx.reshape(len(sel), kk)
    ds = np.full(len(sel), np.inf)
    for j in range(kk):
        t = tri[idx[:, j]]
        ds = np.minimum(ds, point_triangle_distance(mesh.points[sel], t[:, 0], t[:, 1], t[:, 2]))
    d[sel] = ds if band is None else np.minimum(ds, band)
    return np.where(phi < 0, -d, d)


def advect(mesh, phi, velocity, dt):
    """Semi-Lagrangian step: phi_new(x) = phi(x - u(x) dt), then reset to a distance."""
    x = mesh.points
    return redistance(mesh, interpolate(mesh, phi, x - velocity(x) * dt))


def _face_pairs(mesh, cells):
    T = mesh.tets[cells]
    f = np.vstack([T[:, [0, 1, 2]], T[:, [0, 1, 3]], T[:, [0, 2, 3]], T[:, [1, 2, 3]]])
    f.sort(axis=1)
    owner = np.tile(np.arange(len(cells)), 4)
    _, inv = np.unique(f, axis=0, return_inverse=True)
    inv = inv.ravel()
    order = np.argsort(inv, kind="stable")
    s_inv, s_own = inv[order], owner[order]
    same = s_inv[1:] == s_inv[:-1]
    return s_own[:-1][same], s_own[1:][same]


def pieces(mesh, phi):
    """Connected pieces of metal (cells holding any metal, joined through shared faces).
    Returns (label per cell, -1 where no metal; number of pieces), pieces ordered by decreasing volume."""
    vol = inside_volumes(mesh, phi)
    cells = np.flatnonzero(vol > 0)
    label = np.full(len(mesh.tets), -1)
    if len(cells) == 0:
        return label, 0
    a, b = _face_pairs(mesh, cells)
    g = coo_matrix((np.ones(len(a)), (a, b)), shape=(len(cells), len(cells)))
    n, lab = connected_components(g, directed=False)
    pv = np.bincount(lab, vol[cells], minlength=n)
    rank = np.empty(n, int)
    rank[np.argsort(-pv, kind="stable")] = np.arange(n)
    label[cells] = rank[lab]
    return label, n


def piece_volumes(mesh, phi, label, n):
    vol = inside_volumes(mesh, phi)
    m = label >= 0
    return np.bincount(label[m], vol[m], minlength=n)


def correct_volumes(mesh, phi, label, targets, rings=2, rtol=1e-14, max_iter=60):
    """Shift phi uniformly over each piece (its cells' nodes, dilated by `rings` rings of cells) so that the
    piece's volume equals targets[p] to round-off. Returns (phi, largest relative correction)."""
    phi = phi.copy()
    worst = 0.0
    targets = np.asarray(targets, float)
    if (targets <= 0).any():
        raise ValueError("volume target of piece {} is not positive; remove the piece instead".format(
            int(np.flatnonzero(targets <= 0)[0])))
    for p, target in enumerate(targets):
        others = np.zeros(len(mesh.points), bool)
        others[mesh.tets[(label >= 0) & (label != p)].ravel()] = True       # nodes of every other piece
        nodes = np.unique(mesh.tets[label == p])
        for _ in range(rings):
            touch = np.isin(mesh.tets, nodes).any(axis=1)
            nodes = np.unique(mesh.tets[touch])
        nodes = nodes[~others[nodes]]
        cells = np.flatnonzero(np.isin(mesh.tets, nodes).any(axis=1) & ((label < 0) | (label == p)))
        base = phi[nodes].copy()
        lo, hi, c = None, None, 0.0
        V0 = inside_volumes(mesh, phi, cells).sum()
        for _ in range(max_iter):
            phi[nodes] = base - c
            V = inside_volumes(mesh, phi, cells).sum()
            err = target - V
            if abs(err) <= rtol * target:
                break
            if err > 0:
                lo = c
            else:
                hi = c
            A = patches_area(mesh, phi, cells)
            c_new = c + err / max(A, 1e-300)
            if (lo is not None and c_new <= lo) or (hi is not None and c_new >= hi):
                c_new = 0.5 * (lo + hi) if (lo is not None and hi is not None) else c + 2 * (c_new - c)
            c = c_new
        else:
            raise RuntimeError("volume correction of piece {} did not converge".format(p))
        worst = max(worst, abs(target - V0) / target)
    return phi, worst


def patches_area(mesh, phi, cells):
    return float(patches(mesh, phi, cells).area.sum())

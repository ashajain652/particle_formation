"""Operations on a P1 level set (negative inside the metal) on an AdaptiveMesh: reset to a signed distance, motion
under a prescribed velocity, connected pieces, and the exact volume correction per piece (design §6.2-§6.4; §20.3:
the evolved level set is the skin base beta)."""
import numpy as np
from scipy.optimize import brentq
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
    within band of a triangle are computed (exactly, by the search below); the rest get +-band, never more than their
    true distance."""
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
        r_max = np.linalg.norm(tri - tri.mean(axis=1, keepdims=True), axis=2).max()   # centroid to farthest vertex
        near, _ = kd.query(mesh.points, k=1, distance_upper_bound=band + r_max)
        sel = np.flatnonzero(np.isfinite(near))
        _, idx = kd.query(mesh.points[sel], k=kk)
    idx = idx.reshape(len(sel), kk)
    ds = np.full(len(sel), np.inf)
    for j in range(kk):
        t = tri[idx[:, j]]
        ds = np.minimum(ds, point_triangle_distance(mesh.points[sel], t[:, 0], t[:, 1], t[:, 2]))
    d[sel] = ds if band is None else np.minimum(ds, band)
    return np.where(phi < 0, -d, d)


def advect(mesh, phi, velocity, dt, band=None):
    """Semi-Lagrangian step: phi_new(x) = phi(x - u(x) dt), then reset to a distance."""
    x = mesh.points
    return redistance(mesh, interpolate(mesh, phi, x - velocity(x) * dt), band=band)


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


SWEEPS_MAX = 20     # sweeps over the pieces at most in `correct_volumes`
BRENT_XTOL = 1e-18  # m: Brent's method's tolerance on a piece's shift (`correct_volumes`), far below any cell
BRENT_ITER = 200    # its iterations at most


def correct_volumes(mesh, phi, label, targets, rings=2, rtol=1e-14, max_iter=60, misses=None):
    """Shift phi uniformly over each piece (its cells' nodes, dilated by `rings` rings of cells) so that the
    piece's volume equals targets[p] to round-off. Returns (phi, largest relative correction); with `misses` a list,
    each piece's miss before the correction, as a share of its target, is appended to it in label order (stage 3 Task
    11's fix round 2: the body's, piece 0's, beside the worst).

    Pieces whose dilations meet (Task 9, fix round 4): each node is shifted with the first piece (in label order) that
    reaches it, and each cell holding no metal before is counted with the first piece that reaches it, so no metal is
    counted for two pieces; the pieces are corrected in turn, sweep after sweep (at most SWEEPS_MAX), until every piece
    holds its target at once, since a later piece's shift moves the contour in a cell an earlier piece counts. (Counted
    for both, the metal a shift grew in a cell between two pieces 0.25 to 0.4 mm apart on 0.125 mm cells was booked
    twice and the total fell short by up to 3e-5; on the full-size cube's second step, with a 2e-9 m3 piece beside the
    body, a short 1e-14 m3 a move made the step's mass residual -9.4e-9.) One piece is corrected as before. A piece
    that does not converge raises with its target, the volume it held, its last shift and bracket, the iterations, its
    claimed nodes and cells, and the miss before the correction (stage 3 Task 11's review, M6).

    Each piece's shift is found by Newton steps on the contour's area, kept within the bracket of shifts that held too
    little and too much: bisection of the bracket where a Newton step leaves it, and Brent's method on it where a Newton
    step fails to halve the step before the last (rtsafe's test; stage 3 Task 11's fix round 1). A speck a base move makes beside the body,
    asked to grow some 1e5 times, first steps far beyond the nodes it may shift; there its volume no longer responds
    while the cut cells at its claim's edge keep a patch area, and the bare Newton steps crawled for all 60 iterations
    (the q6p box's step to 1.5 s and the shoulder's to 2.0 s stopped so)."""
    phi = phi.copy()
    targets = np.asarray(targets, float)
    if (targets <= 0).any():
        raise ValueError("volume target of piece {} is not positive; remove the piece instead".format(
            int(np.flatnonzero(targets <= 0)[0])))
    claimed_nodes = np.zeros(len(mesh.points), bool)
    claimed_cells = np.zeros(len(mesh.tets), bool)
    sets = []
    for p in range(len(targets)):
        others = np.zeros(len(mesh.points), bool)
        others[mesh.tets[(label >= 0) & (label != p)].ravel()] = True       # nodes of every other piece
        nodes = np.unique(mesh.tets[label == p])
        for _ in range(rings):
            touch = np.isin(mesh.tets, nodes).any(axis=1)
            nodes = np.unique(mesh.tets[touch])
        nodes = nodes[~others[nodes] & ~claimed_nodes[nodes]]
        claimed_nodes[nodes] = True
        cells = np.flatnonzero(np.isin(mesh.tets, nodes).any(axis=1) & ((label < 0) | (label == p)) & ~claimed_cells)
        claimed_cells[cells] = True
        sets.append((nodes, cells, phi[nodes].copy()))
    shift = np.zeros(len(targets))
    worst = 0.0
    for sweep in range(SWEEPS_MAX):
        for p, target in enumerate(targets):
            nodes, cells, base = sets[p]
            lo, hi, c = None, None, shift[p]
            step_before = step = np.inf                    # the last two steps (the bisection safeguard below)
            for _ in range(max_iter):
                phi[nodes] = base - c
                V, c_held = inside_volumes(mesh, phi, cells).sum(), c
                err = target - V
                if sweep == 0 and c == 0.0:
                    worst = max(worst, abs(err) / target)
                    if misses is not None:
                        misses.append(abs(err) / target)
                if abs(err) <= rtol * target:
                    break
                if err > 0:
                    lo = c
                else:
                    hi = c
                A = patches_area(mesh, phi, cells)
                c_new = c + err / max(A, 1e-300)
                bracketed = lo is not None and hi is not None
                if (lo is not None and c_new <= lo) or (hi is not None and c_new >= hi):
                    c_new = 0.5 * (lo + hi) if bracketed else c + 2 * (c_new - c)
                elif bracketed and abs(c_new - c) > 0.5 * abs(step_before):
                    # Newton not halving the step before the last: the contour's area is not the volume's slope here
                    # (cut cells at the claim's edge keep an area that does not respond), so Brent's method on the bracket
                    def held(x, nodes=nodes, base=base, cells=cells, target=target):
                        phi[nodes] = base - x
                        return target - inside_volumes(mesh, phi, cells).sum()
                    c_new = brentq(held, lo, hi, xtol=BRENT_XTOL, rtol=4.0 * np.finfo(float).eps, maxiter=BRENT_ITER)
                step_before, step = step, c_new - c
                c = c_new
            else:
                raise RuntimeError("volume correction of piece {} did not converge: target {:.6e} m3, held {:.6e} m3 "
                                   "at the shift {:.6e} m, bracket [{}, {}] m, {} iterations (sweep {}), {} nodes "
                                   "and {} cells claimed, the largest miss before the correction {:.3e}".format(
                                       p, target, V, c_held, "none" if lo is None else "{:.6e}".format(lo),
                                       "none" if hi is None else "{:.6e}".format(hi), max_iter, sweep, len(nodes),
                                       len(cells), worst))
            shift[p] = c
        if len(targets) == 1 or all(abs(t - inside_volumes(mesh, phi, sets[p][1]).sum()) <= rtol * t
                                    for p, t in enumerate(targets)):
            return phi, worst
    raise RuntimeError("the volume correction did not hold every piece's target at once in {} sweeps".format(SWEEPS_MAX))


def patches_area(mesh, phi, cells):
    return float(patches(mesh, phi, cells).area.sum())

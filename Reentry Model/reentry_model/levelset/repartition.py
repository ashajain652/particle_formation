"""The re-partition of the skins at the step end (level-set stage 3 Task 9b, the controller's rulings of 2026-11-12): the
area a tessellation hands out afresh between the markers, and the skins' metal following it.

The transfer. `transfer_matrix` follows the area a tessellation hands out afresh, for the body's step end
(`stepend.StepEnd._repartition`). The surface before and after is sampled (`surface_samples`). Each sample of the old
surface belongs to the marker nearest it at its old position (`Markers.x_tess`); each sample of the new surface to its
old owner, the nearest marker at its old position carried along its old normal by the normal part of its move (so that
a receding face keeps its owners), and to its new owner, the nearest at its new position. A marker's stretch is the new
surface's area its old owners hold over the old surface's, each summed over the marker and its neighbours (the new
faces) before the ratio, then scaled, on each piece of the faces' graph, so that the old areas times the stretch sum to
the new areas exactly (a marker whose neighbourhood holds no sample takes the piece's ratio of new to old area). The
transfer starts from the new-surface samples whose old and new owners are one marker or neighbours, each row scaled to
its old area times its stretch (a row with none keeps that area itself); its columns are then brought to the new areas
by a potential flow across the faces (`_fit_columns`), which moves area only within rows and turns no entry negative.
Rows and columns hold to round-off, and when nothing moves the transfer is the old areas on its diagonal. Measured once
each (Task 9b): a cold step end of the test cube, markers moved 13 nm, passed 3.8e-7 of the area between markers; a
10 mm sphere's markers carried onto 10.4 mm found stretches of 1.0749 to 1.0879 against (10.4/10)^2 = 1.0816, with
0.35 % of the area passed between markers.

The skins follow it (`apply_transfer`): each old column is scaled by its stretch (its top area becomes its old area
times its stretch, its metal unchanged; design 9.2's draw or return at its base follows in the next part, by the draw to
the minimum mass and the hand-down), cut into the fractions of its row of the transfer (`Skins.cut`), and each new column
is the hand-over mix (`Skins.combine`) of the pieces it receives.

Numbers here are numerical settings; no physics."""
import math

import numpy as np
from scipy.sparse import coo_matrix, diags, identity
from scipy.sparse.csgraph import connected_components
from scipy.sparse.linalg import spsolve
from scipy.spatial import cKDTree

from .cutcell import patches as cut_patches

SAMPLE_SHARE = 1.0 / 16.0   # of the markers' target spacing: the surface samples' size for a re-partition
SAMPLE_LEVEL_MAX = 5        # midpoint subdivisions of a patch triangle at most for the surface samples
                            # (`surface_samples`)
FIT_PASSES = 4              # flows at most that bring the re-partition's columns to the new areas (`_fit_columns`)
FIT_RTOL = 1e-14            # of a piece's mean area: the fit stops when no column is farther from its area than this
FIT_FAIL = 1e-12            # of a piece's mean area: a fit still farther after FIT_PASSES raises (fix round 1: of the
                            # mean, not of the column's own area, which may be zero)


def _subdivision(level):
    """Barycentric centroids (4**level, 3) of the triangles a triangle is cut into by its edge midpoints `level` times."""
    tri = np.eye(3)[None, :, :]
    for _ in range(level):
        a, b, c = tri[:, 0], tri[:, 1], tri[:, 2]
        ab, bc, ca = 0.5 * (a + b), 0.5 * (b + c), 0.5 * (c + a)
        tri = np.concatenate([np.stack(t, axis=1) for t in ((a, ab, ca), (ab, b, bc), (ca, bc, c), (ab, bc, ca))])
    return tri.mean(axis=1)


def surface_samples(mesh, phi, size):
    """Points sampling phi's zero contour, with the area each stands for: every patch triangle cut into 4**L equal
    triangles by its edge midpoints, L the least level (at most SAMPLE_LEVEL_MAX) at which they are at most `size` [m]
    across (the square root of twice their area), each at its centroid. Returns (points (k, 3), areas (k,) [m2]); the
    areas sum to the patches' area to round-off."""
    t = cut_patches(mesh, phi).tris
    a = 0.5 * np.linalg.norm(np.cross(t[:, 1] - t[:, 0], t[:, 2] - t[:, 0]), axis=1)
    level = np.clip(np.ceil(np.log2(np.maximum(np.sqrt(2.0 * a), 1e-300) / size)), 0, SAMPLE_LEVEL_MAX).astype(int)
    y, w = [np.zeros((0, 3))], [np.zeros(0)]
    for L in np.unique(level):
        sel = level == L
        bary = _subdivision(int(L))
        y.append(np.einsum("sk,qkd->qsd", bary, t[sel]).reshape(-1, 3))
        w.append(np.repeat(a[sel] / len(bary), len(bary)))
    return np.concatenate(y), np.concatenate(w)


def transfer_matrix(old_x, old_normal, old_area, old_samples, x, area, faces, samples):
    """The step end's re-partition (module docstring): n markers, the same before and after, each
    owning before the area `old_area` [m2] (its skin's) about `old_x` (n, 3) with the normal `old_normal` on the
    surface sampled by `old_samples`, and after the tessellated `area` about `x` on the surface sampled by `samples`
    (each a (points, areas) pair from `surface_samples`), with the tessellation's `faces`. Returns (T, stretch): T
    (n, n) csr [m2], T[j, i] the area of marker j's old column that marker i's new column receives, its rows summing to
    old_area * stretch and its columns to `area`, each to round-off; `stretch` (n,) each old column's stretch."""
    old_x, old_normal, x = (np.asarray(v, dtype=float).reshape(-1, 3) for v in (old_x, old_normal, x))
    old_area, area = np.asarray(old_area, dtype=float), np.asarray(area, dtype=float)
    n = len(x)
    (y0, w0), (y1, w1) = old_samples, samples
    carried = old_x + np.einsum("ij,ij->i", x - old_x, old_normal)[:, None] * old_normal
    _, own0 = cKDTree(old_x).query(y0)                      # the old surface's samples by their old owners
    _, own1 = cKDTree(carried).query(y1)                    # the new surface's by theirs, carried along the old normal
    _, new1 = cKDTree(x).query(y1)                          # and by their new owners
    fi, fj = (np.zeros(0, np.int64),) * 2 if faces is None else (faces.i, faces.j)
    F = coo_matrix((np.ones(len(fi)), (fi, fj)), shape=(n, n))
    W = ((F + F.T + identity(n)) > 0).astype(float).tocsr()
    ncomp, comp = connected_components(W, directed=False)
    with np.errstate(divide="ignore", invalid="ignore"):   # each marker's stretch, smoothed over its neighbours
        s = (W @ np.bincount(own1, w1, minlength=n)) / (W @ np.bincount(own0, w0, minlength=n))
    R = np.zeros(n)
    for c in range(ncomp):                                  # normalised to the covered area's change, piece by piece
        m = np.flatnonzero(comp == c)
        A_new, A_old = math.fsum(area[m]), math.fsum(old_area[m])
        sc = np.where(np.isfinite(s[m]) & (s[m] > 0.0), s[m], A_new / A_old if A_old > 0.0 else 1.0)
        R[m] = old_area[m] * sc
        total = math.fsum(R[m])
        if total > 0.0:
            R[m] *= A_new / total
            big = m[np.argmax(R[m])]                       # the rows' total to the new areas' to round-off: what the
            R[big] += A_new - math.fsum(R[m])              # scaling's rounding left, on the largest row
        else:
            R[m] = 0.0
    stretch = np.divide(R, old_area, out=np.ones(n), where=old_area > 0.0)
    K = coo_matrix((w1, (own1, new1)), shape=(n, n)).tocsr().multiply(W).tocsr()   # between neighbours or itself
    rs = np.asarray(K.sum(axis=1)).ravel()
    T = (diags(np.divide(R, rs, out=np.zeros(n), where=rs > 0.0)) @ K
         + diags(np.where(rs > 0.0, 0.0, R))).tocsr()          # rows to their stretched old areas
    length = np.zeros(0) if faces is None else faces.length
    return _fit_columns(T, R, area, fi, fj, length, x, comp), stretch


def _fit_columns(T, rows, area, fi, fj, length, x, comp):
    """T (csr) with each column's sum brought to `area` and each row's to `rows`: the columns' shortfalls, r = area less the
    column sums, are met by a flow of area between the columns of neighbours across their faces, the potential flow of
    the faces' graph (conductance the face length over the markers' distance; one potential held per piece `comp`)
    with L p = -r, so that column i passes g = c (p_i - p_k) to column k; the columns are visited from the highest
    potential down, and what one passes on it takes from its own column's metal first (row i's entry T[i, i], moved to
    T[i, k]) and the rest from what it received, in proportion (row j's T[j, i] moved to T[j, k]): every move stays
    within its row and nothing turns negative, since a column passes on less than it holds (it keeps its area).
    Repeated on the round-off left, FIT_PASSES times at most, until no column but the one whose potential is held is
    further from its area than FIT_RTOL of its piece's mean area (a column whose new area is zero, too, is then met); a
    fit with any column still FIT_FAIL of that short raises. Before each pass every row is scaled back to `rows`; the
    sums are exactly rounded (`math.fsum`), so the tolerances hold for the true sums. The held column carries what no
    flow can move, the piece's total less its new areas' total: the rounding of every move and scaling, summed (fix round
    1: on 443 markers turned until 94 % of the area changed hands, 6e-14 of a mean area, against n eps = 1e-13)."""
    n = len(area)
    cols = [dict() for _ in range(n)]
    coo = T.tocoo()
    for j, i, v in zip(coo.row.tolist(), coo.col.tolist(), coo.data.tolist()):
        if v > 0.0:
            cols[i][j] = cols[i].get(j, 0.0) + v

    def shortfall():                                       # each column's sum exactly rounded (math.fsum)
        return area - np.array([math.fsum(c.values()) for c in cols])

    def rescale_rows():                                    # each row back to its target, its entries in proportion
        held = [[] for _ in range(n)]
        for c in cols:
            for j, v in c.items():
                held[j].append(v)
        f = [rows[j] / math.fsum(h) if h else 1.0 for j, h in enumerate(held)]
        for c in cols:
            for j in c:
                c[j] *= f[j]
    scale = (np.bincount(comp, area) / np.bincount(comp))[comp]   # each column's piece's mean area
    if len(fi):
        g = length / np.linalg.norm(x[fi] - x[fj], axis=1)
        C = coo_matrix((g, (fi, fj)), shape=(n, n))
        C = (C + C.T).tocsr()
        L = (diags(np.asarray(C.sum(axis=1)).ravel()) - C).tocsr()
        free = np.ones(n, bool)
        free[np.unique(comp, return_index=True)[1]] = False   # one potential held per piece
        Lf = L[free][:, free].tocsc()
    else:
        free = np.zeros(n, bool)
    for _ in range(FIT_PASSES):
        rescale_rows()
        r = shortfall()
        if np.all(np.abs(r[free]) <= FIT_RTOL * scale[free]):
            break
        p = np.zeros(n)
        if free.any():
            p[free] = spsolve(Lf, -r[free])
        out = [[] for _ in range(n)]
        for a, b, q in zip(fi.tolist(), fj.tolist(), (g * (p[fi] - p[fj])).tolist()):
            if q > 0.0:
                out[a].append((b, q))
            elif q < 0.0:
                out[b].append((a, -q))
        for i in np.argsort(-p, kind="stable").tolist():
            if not out[i]:
                continue
            col = cols[i]
            O = sum(q for _, q in out[i])
            own = col.get(i, 0.0)
            take = min(own, O)
            for k, q in out[i]:
                cols[k][i] = cols[k].get(i, 0.0) + q * (take / O)
            col[i] = own - take if take < own else 0.0
            others = [(h, v) for h, v in col.items() if h != i and v > 0.0]
            if O > take and others:
                frac = min((O - take) / sum(v for _, v in others), 1.0)
                for h, v in others:
                    moved = v * frac
                    col[h] = v - moved
                    for k, q in out[i]:
                        cols[k][h] = cols[k].get(h, 0.0) + moved * (q / O)
    r = shortfall()
    if np.any(np.abs(r) > FIT_FAIL * scale):
        raise RuntimeError("the re-partition's columns did not meet the new areas: {:.1e} of their mean area left".format(
            float(np.max(np.abs(r) / scale))))
    rows = [j for c in cols for j in c]
    at = [i for i, c in enumerate(cols) for _ in c]
    vals = [v for c in cols for v in c.values()]
    fitted = coo_matrix((vals, (rows, at)), shape=(n, n)).tocsr()
    fitted.eliminate_zeros()
    return fitted


def apply_transfer(skins, rows, old_area, T, stretch):
    """The skins follow a re-partition: each old column (skin row rows[j] of marker j, `skins` a `Skins`) scaled by its
    stretch -- its top area old_area[j] * stretch[j], its metal unchanged -- then cut into the fractions T[j, i] / (row
    j's sum) of its profile (`Skins.cut`), and each new column the hand-over mix (`Skins.combine`) of what it receives,
    one piece per receiver at a time. The pieces are appended after every marker's row and removed once emptied. Mass,
    enthalpy and momentum add to round-off; the caller gives the columns their new areas and curvatures."""
    skins.area[rows] = old_area * stretch
    R = np.asarray(T.sum(axis=1)).ravel()
    coo = T.tocoo()
    give = (coo.row != coo.col) & (coo.data > 0.0)
    donor, receiver = coo.row[give], coo.col[give]
    if len(donor):
        pieces = skins.cut(rows[donor], coo.data[give] / R[donor])
        order = np.argsort(receiver, kind="stable")
        donor, receiver, pieces = donor[order], receiver[order], pieces[order]
        first = np.r_[0, np.flatnonzero(np.diff(receiver)) + 1]
        rank = np.arange(len(receiver)) - np.repeat(first, np.diff(np.r_[first, len(receiver)]))
        for k in range(int(rank.max()) + 1):            # one piece per receiver at a time: distinct pairs
            at = rank == k
            skins.combine(rows[receiver[at]], pieces[at])
        skins.remove(pieces)                            # emptied by combine; appended after every marker's row

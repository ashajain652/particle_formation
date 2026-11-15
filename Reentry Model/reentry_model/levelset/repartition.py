"""The re-partition of the skins at the step end (level-set stage 3 Task 9b, the controller's rulings of 2026-11-12; plan
03b Task 6, the controller's ruling P4): the layer a tessellation and the division that follows it hand out afresh
between the markers' columns, and the skins' metal following it by depth.

The transfer. `transfer_matrix` follows the columns' volumes a step end's settle hands out afresh, for the body's step
end (`stepend.StepEnd._repartition`). The surface before and after is sampled (`surface_samples`). Each sample of the
old surface belongs to the marker nearest it at its old position (`Markers.x_tess`); each sample of the new surface to
its old owner, the nearest marker at its old position carried along its old normal by the normal part of its move (so
that a receding face keeps its owners), and to its new owner, the nearest at its new position. A marker's stretch is
the new surface's area its old owners hold over the old surface's, each summed over the marker and its neighbours (the
columns sharing a boundary with its column in the new division) before the ratio, then scaled, on each piece of the
boundaries' graph, so that the old columns' volumes times the stretch sum to the new columns' volumes exactly (a marker
whose neighbourhood holds no sample takes the piece's ratio of new to old volume). The transfer starts from the
new-surface samples whose old and new owners are one marker or neighbours, each row scaled to its old volume times its
stretch (a row with none keeps that volume itself); its columns are then brought to the new volumes by a potential flow
across the division's boundaries (`_fit_columns`, each boundary's conductance its area over the layer, the integral of
its length over depth, over the markers' distance), which moves volume only within rows and turns no entry negative.
Rows and columns hold to round-off, and when nothing moves the transfer is the old volumes on its diagonal.

Why volumes (plan 03b Task 6): the division's per-piece closure makes its columns' volumes sum to the layer's volume
exactly, before and after a settle on the same surface, so the stretch of an unchanged surface is one to round-off; and
where the columns are flat each new column receives the volume its new division gives it, so the fill that follows moves
nothing through the base. (By the markers' areas of phi, with the pieces mixed by mass fraction, the fill after a
re-partition of a flat box's top face moved 2.3e-3 of its inner columns' volume through the base, plan 03b Task 5,
measured once.) Measured before Task 6, on areas and the tangent-plane faces, once each (Task 9b): a cold step end of the
test cube, markers moved 13 nm, passed 3.8e-7 of the area between markers; a 10 mm sphere's markers carried onto 10.4 mm
found stretches of 1.0749 to 1.0879 against (10.4/10)^2 = 1.0816, with 0.35 % of the area passed between markers.

The skins follow it (`apply_transfer`): each old column is scaled by its stretch -- its top area and its measured
cross-section, its metal unchanged (the stretch carries the profile: a stretched surface spreads the same metal over
more layer; the next fill draws the rest at its base) -- cut into the fractions of its row of the transfer
(`Skins.cut`, each piece that fraction of the column's metal and of its profile), and each new column takes the pieces
it receives by depth (`Skins.combine(..., by="depth")`: each piece's metal joins the new column's cells at the depths it
held, design 21.5's rule that metal stays where it is and only its owner changes).

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
FIT_RTOL = 1e-14            # of a piece's mean column (its volume): the fit stops when no column is farther from its
                            # own than this
FIT_FAIL = 1e-12            # of a piece's mean column: a fit still farther after FIT_PASSES raises (fix round 1: of the
                            # mean, not of the column's own, which may be zero)


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


def _boundaries(division):
    """The division's shared boundaries with an area over the layer: (i, j) markers (E,), each boundary's area [m2] (the
    trapezoid integral of its lengths over the levels)."""
    lv = np.asarray(division.levels, dtype=float)
    L = np.asarray(division.length, dtype=float)
    area = (0.5 * (L[:, :-1] + L[:, 1:]) * np.diff(lv)).sum(axis=1) if len(L) else np.zeros(0)
    keep = area > 0.0
    return np.asarray(division.i)[keep], np.asarray(division.j)[keep], area[keep]


def transfer_matrix(old_x, old_normal, old_volume, old_samples, x, volume, division, samples):
    """The step end's re-partition (module docstring): n markers, the same before and after, each one's old column of
    the volume `old_volume` [m3] (its skin's profile's) owning the surface about `old_x` (n, 3) with the normal
    `old_normal` on the surface sampled by `old_samples`, and after the new column of the volume `volume` about `x` on the
    surface sampled by `samples` (each a (points, areas) pair from `surface_samples`), of the new `division` of the layer
    among the markers (`division.Division`, its shared boundaries the columns' neighbours). Returns (T, stretch): T (n,
    n) csr [m3], T[j, i] the volume of marker j's old column that marker i's new column receives, its rows summing to
    old_volume * stretch and its columns to `volume`, each to round-off; `stretch` (n,) each old column's stretch."""
    old_x, old_normal, x = (np.asarray(v, dtype=float).reshape(-1, 3) for v in (old_x, old_normal, x))
    old_volume, volume = np.asarray(old_volume, dtype=float), np.asarray(volume, dtype=float)
    n = len(x)
    (y0, w0), (y1, w1) = old_samples, samples
    carried = old_x + np.einsum("ij,ij->i", x - old_x, old_normal)[:, None] * old_normal
    _, own0 = cKDTree(old_x).query(y0)                      # the old surface's samples by their old owners
    _, own1 = cKDTree(carried).query(y1)                    # the new surface's by theirs, carried along the old normal
    _, new1 = cKDTree(x).query(y1)                          # and by their new owners
    fi, fj, face = _boundaries(division)
    F = coo_matrix((np.ones(len(fi)), (fi, fj)), shape=(n, n))
    W = ((F + F.T + identity(n)) > 0).astype(float).tocsr()
    ncomp, comp = connected_components(W, directed=False)
    with np.errstate(divide="ignore", invalid="ignore"):   # each marker's stretch, smoothed over its neighbours
        s = (W @ np.bincount(own1, w1, minlength=n)) / (W @ np.bincount(own0, w0, minlength=n))
    R = np.zeros(n)
    for c in range(ncomp):                                  # normalised to the new columns' volume, piece by piece
        m = np.flatnonzero(comp == c)
        V_new, V_old = math.fsum(volume[m]), math.fsum(old_volume[m])
        sc = np.where(np.isfinite(s[m]) & (s[m] > 0.0), s[m], V_new / V_old if V_old > 0.0 else 1.0)
        R[m] = old_volume[m] * sc
        total = math.fsum(R[m])
        if total > 0.0:
            R[m] *= V_new / total
            big = m[np.argmax(R[m])]                       # the rows' total to the new volumes' to round-off: what the
            R[big] += V_new - math.fsum(R[m])              # scaling's rounding left, on the largest row
        else:
            R[m] = 0.0
    stretch = np.divide(R, old_volume, out=np.ones(n), where=old_volume > 0.0)
    K = coo_matrix((w1, (own1, new1)), shape=(n, n)).tocsr().multiply(W).tocsr()   # between neighbours or itself
    rs = np.asarray(K.sum(axis=1)).ravel()
    T = (diags(np.divide(R, rs, out=np.zeros(n), where=rs > 0.0)) @ K
         + diags(np.where(rs > 0.0, 0.0, R))).tocsr()          # rows to their stretched old volumes
    return _fit_columns(T, R, volume, fi, fj, face, x, comp), stretch


def _fit_columns(T, rows, area, fi, fj, length, x, comp):
    """T (csr) with each column's sum brought to `area` and each row's to `rows` (the columns' measure: the new columns'
    volumes for `transfer_matrix`): the columns' shortfalls, r = area less the column sums, are met by a flow between the
    columns of neighbours across their faces, the potential flow of the faces' graph (conductance the face's `length` --
    for `transfer_matrix` the boundary's area over the layer -- over the markers' distance; one potential held per piece
    `comp`) with L p = -r, so that column i passes g = c (p_i - p_k) to column k; the columns are visited from the
    highest potential down, and what one passes on it takes from its own column's metal first (row i's entry T[i, i],
    moved to T[i, k]) and the rest from what it received, in proportion (row j's T[j, i] moved to T[j, k]): every move
    stays within its row and nothing turns negative, since a column passes on less than it holds (it keeps its measure).
    Repeated on the round-off left, FIT_PASSES times at most, until no column but the one whose potential is held is
    further from its measure than FIT_RTOL of its piece's mean (a column whose new measure is zero, too, is then met); a
    fit with any column still FIT_FAIL of that short raises. Before each pass every row is scaled back to `rows`; the
    sums are exactly rounded (`math.fsum`), so the tolerances hold for the true sums. The held column carries what no
    flow can move, the piece's total less its new measures' total: the rounding of every move and scaling, summed (fix
    round 1: on 443 markers turned until 94 % of the area changed hands, 6e-14 of a mean area, against n eps = 1e-13)."""
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
    scale = (np.bincount(comp, area) / np.bincount(comp))[comp]   # each column's piece's mean measure
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
        raise RuntimeError("the re-partition's columns did not meet their new measures: {:.1e} of their mean left".format(
            float(np.max(np.abs(r) / scale))))
    rows = [j for c in cols for j in c]
    at = [i for i, c in enumerate(cols) for _ in c]
    vals = [v for c in cols for v in c.values()]
    fitted = coo_matrix((vals, (rows, at)), shape=(n, n)).tocsr()
    fitted.eliminate_zeros()
    return fitted


def apply_transfer(skins, rows, T, stretch):
    """The skins follow a re-partition (module docstring): each old column (skin row rows[j] of marker j, `skins` a
    `Skins`) scaled by its stretch -- its top area and its measured cross-section (the stretch carries the profile), its
    metal unchanged -- then cut into the fractions T[j, i] / (row j's sum) of itself (`Skins.cut`), and each new column
    takes what it receives by depth (`Skins.combine(..., by="depth")`), one piece per receiver at a time. The pieces are
    appended after every marker's row and removed once emptied. Mass, enthalpy and momentum add to round-off; the caller
    gives the columns their new areas, curvatures and profiles."""
    stretch = np.asarray(stretch, dtype=float)
    skins.area[rows] = skins.area[rows] * stretch
    if skins.profile is not None:
        skins.profile[rows] = skins.profile[rows] * stretch[:, None]
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
            skins.combine(rows[receiver[at]], pieces[at], by="depth")
        skins.remove(pieces)                            # emptied by combine; appended after every marker's row

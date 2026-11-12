"""The thin test of design §20.5: per surface patch, the molten depth h (phi to psi), the largest principal curvature
of phi at the patch and of psi at its foot, and the slope of the depth along the surface; with hysteresis (a patch
changes mode only after `hold` consecutive opposite answers, together with a connected group of `min_group`)."""
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

from .criteria import BIG, finite_values, hessian_field, nodal_average, principal_curvatures
from .cutcell import cell_gradients
from .transfer import locate


def thin_candidate(h, kmax, slope, c_kappa=0.1, c_slope=0.2):
    """Thin when h <= c_kappa / kmax (h at most c_kappa times the smallest radius of curvature) and slope <= c_slope.
    kmax is the largest |principal curvature| of either surface: a neck has small mean but large principal
    curvature."""
    return (np.asarray(h) * np.asarray(kmax) <= c_kappa) & (np.asarray(slope) <= c_slope)


def patch_measures(mesh, phi, psi, P):
    """(h, kmax, slope) at the patches P of phi. h is the nearest distance from the patch to the zero contour of psi
    (never more than the depth along the normal that the design names; an approximation revisited in stage 3).
    Where psi is infinite at the patch (no rigid metal: fully molten) h is infinite (thick), the slope is 0 and
    kmax comes from phi alone."""
    x, n = P.centroid, P.normal
    cell, b = locate(mesh, x)
    T = mesh.tets[cell]
    h = np.maximum(np.einsum("ni,ni->n", b, finite_values(psi)[T]), 0.0)
    h = np.where(h > 0.5 * BIG, np.inf, h)
    molten = np.isinf(h)
    psi_c = np.clip(psi, -1e3, 1e3)                         # finite copy for the derivative recoveries
    near = mesh.ring(cell, 1)
    gpsi = nodal_average(mesh, cell_gradients(mesh, psi_c, near), near)
    g = np.einsum("ni,nik->nk", b, gpsi[T])
    slope = np.linalg.norm(g - np.einsum("nk,nk->n", g, n)[:, None] * n, axis=1)
    slope = np.where(molten, 0.0, slope)
    k1, k2 = principal_curvatures(mesh, phi, x)
    kmax = np.max(np.abs(np.stack([k1, k2])), axis=0)
    use = ~molten & (h > 0)
    if use.any():
        foot = x[use] - h[use][:, None] * n[use]
        q1, q2 = principal_curvatures(mesh, psi_c, foot)
        kmax[use] = np.maximum(kmax[use], np.maximum(np.abs(q1), np.abs(q2)))
    return h, kmax, slope


def patch_adjacency(mesh, P):
    """Patches are neighbours when their cells share a node."""
    T = mesh.tets[P.cell]
    n = len(P.cell)
    inc = coo_matrix((np.ones(4 * n), (np.repeat(np.arange(n), 4), T.ravel())), shape=(n, len(mesh.points))).tocsr()
    A = (inc @ inc.T).tocoo()
    keep = A.row != A.col
    return coo_matrix((np.ones(keep.sum()), (A.row[keep], A.col[keep])), shape=(n, n)).tocsr()


class ModeState:
    """Thin (True) or thick (False) per patch, with hysteresis. Stage 3 keys it by marker instead of patch."""

    def __init__(self, thin, hold=3, min_group=7):
        self.thin = np.asarray(thin, bool).copy()
        self.count = np.zeros(len(self.thin), int)
        self.hold, self.min_group = hold, min_group

    def update(self, candidate, adjacency):
        candidate = np.asarray(candidate, bool)
        differ = candidate != self.thin
        self.count = np.where(differ, self.count + 1, 0)
        ready = np.flatnonzero(differ & (self.count >= self.hold))
        if len(ready) == 0:
            return np.zeros(len(self.thin), bool)
        sub = adjacency[ready][:, ready]
        _, lab = connected_components(sub, directed=False)
        size = np.bincount(lab)
        switch = ready[size[lab] >= self.min_group]
        self.thin[switch] = ~self.thin[switch]
        self.count[switch] = 0
        changed = np.zeros(len(self.thin), bool)
        changed[switch] = True
        return changed

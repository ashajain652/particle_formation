"""The melt film on the surface patches: lubrication velocity and flux, and the explicit upwind runoff transport on
the patch graph (spec Step 3 section 8).

State: m_f per patch [kg], thickness b = m_f / (rho_l A). Lubrication solution with the patch's shear tau, driving
gradient G, thickness b, melt boundary-layer thickness delta_m (from Girin's Eq. 2, spray.melt_layer) and the
liquid viscosity mu_l:
    thin (b <= delta_m):  V_s = tau b / mu_l + G b^2 / (2 mu_l),   q = tau b^2 / (2 mu_l) + G b^3 / (3 mu_l)
    thick (b > delta_m):  V_s = tau delta_m / mu_l,               q = V_s delta_m / 2 + G b^3 / (3 mu_l)
(q per unit width [m^2/s], signed along the surface direction t away from the stagnation point; the velocity
gradient for the instability is V_s / b or V_s / delta_m). Where delta_m is undefined (free-molecular patches) the
film is thin-branch (Couette) by definition.

Runoff: for every edge shared by patches i and j, flux = q_donor * l_edge * (t_donor . n_edge)^+ with n_edge the
in-plane edge normal pointing out of the donor; both directions are evaluated (a negative q reverses the flow toward
the nose). The transport over a macro step is the linearly implicit upwind scheme (I + dt C) m_new = m_old with the
edge coefficients c = q l (t . n)^+ / (A b) [1/s] taken at the start of the step: unconditionally stable, positive,
and mass-conserving to round-off (every edge flux leaves one patch and enters another); its steady state is the
exact nonlinear one. The spec's explicit sub-stepped
scheme was replaced on 2026-09-20: films driven by the free-molecular shear near the equator move at ~10 m/s and
cross the hemisphere many times per 0.5 s step, so an explicit CFL needed 1e4-1e5 sub-steps per macro step."""
from dataclasses import dataclass

import numpy as np


def lubrication(tau, G, b, delta_m, mu_l, b_layer=None):
    """(V_s, q, shear_rate, thick) per patch for the thin/thick lubrication branches (module docstring); `thick` is the
    branch mask, true where the sheared layer does not reach the bottom of the liquid.

    `b` is the mobile film -- the mass that is actually available to move -- while `b_layer`, if given, is the depth of
    liquid beneath the wall (film plus the contiguous molten material under it, `MeltingBody.liquid_layer_depth`). The
    branch is decided on the layer, because whether the gas shear penetrates the whole liquid or only its top delta_m
    is a property of the liquid's depth, not of how much of it the melt bookkeeping has mobilised; the fluxes stay on
    the film. Deciding the branch on the film instead put 79-99.8 % of the patches on the thin branch where the layer
    says 94-100 % are thick (measured on the 100 mm physics flight, 2026-09-22).
    """
    tau, G, b = (np.asarray(x, dtype=float) for x in (tau, G, b))
    delta_m = np.asarray(delta_m, dtype=float)
    layer = b if b_layer is None else np.maximum(np.asarray(b_layer, dtype=float), b)
    thick = np.isfinite(delta_m) & (layer > delta_m)
    V_thin = tau * b / mu_l + G * b * b / (2.0 * mu_l)
    q_thin = tau * b * b / (2.0 * mu_l) + G * b ** 3 / (3.0 * mu_l)
    d = np.where(thick, delta_m, 0.0)
    V_thick = tau * d / mu_l
    q_thick = V_thick * d / 2.0 + G * b ** 3 / (3.0 * mu_l)
    V = np.where(thick, V_thick, V_thin)
    q = np.where(thick, q_thick, q_thin)
    with np.errstate(divide="ignore", invalid="ignore"):
        rate = np.where(thick, V / np.where(thick, d, 1.0), np.where(b > 0.0, V / np.where(b > 0.0, b, 1.0), 0.0))
    return V, q, rate, thick


class Runoff:
    """Explicit upwind transport of the film mass on a SurfaceMesh (rebuilt whenever the surface changes)."""

    def __init__(self, surface, points, windward=None):
        e, i, j = surface.edges()
        if windward is not None:                               # the flow stops at the equator: no edge leads into a leeward patch
            keep = windward[i] & windward[j]
            e, i, j = e[keep], i[keep], j[keep]
        x1, x2 = points[e[:, 0]], points[e[:, 1]]
        self.length = np.linalg.norm(x2 - x1, axis=1)
        mid = 0.5 * (x1 + x2)
        self.i, self.j = i, j
        self.n_i = self._edge_normal(x2 - x1, surface.normals[i], mid - surface.centroids[i])
        self.n_j = self._edge_normal(x2 - x1, surface.normals[j], mid - surface.centroids[j])
        self.n_patches = surface.n_patches

    @staticmethod
    def _edge_normal(edge, normal, outward):
        n = np.cross(edge, normal)
        n = n / np.maximum(np.linalg.norm(n, axis=1), 1e-300)[:, None]
        sign = np.sign(np.einsum("ij,ij->i", n, outward))
        return n * np.where(sign == 0.0, 1.0, sign)[:, None]

    def edge_coefficients(self, q, b, t_hat, areas):
        """Emptying-rate coefficients c_ij, c_ji [1/s] of every edge: the flux i -> j is m_i c_ij with
        c_ij = q_i l (t_i . n_ij)^+ / (A_i b_i) (zero where the patch is dry); q, b per patch."""
        with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
            rate = np.where(b > 0.0, q / (areas * np.where(b > 0.0, b, 1.0)), 0.0)          # 1/m per unit edge length
        rate = np.where(np.isfinite(rate), rate, 0.0)                                       # a degenerate (sliver) patch moves nothing
        c_ij = np.maximum(0.0, rate[self.i] * self.length * np.einsum("ij,ij->i", t_hat[self.i], self.n_i))
        c_ji = np.maximum(0.0, rate[self.j] * self.length * np.einsum("ij,ij->i", t_hat[self.j], self.n_j))
        return c_ij, c_ji

    def transport(self, m_f, q_of_thickness, t_hat, rho_l, areas, dt, substeps=4):
        """Advance m_f over dt in `substeps` linearly implicit upwind steps (I + dt_s C) m_new = m_old, the edge
        coefficients C taken at the start of each sub-step; each M-matrix system is solved directly, so the scheme is
        unconditionally stable, positive and conservative to round-off (the columns of C sum to zero), and its
        steady state C(m) m = 0 is the exact nonlinear one. Dry patches have no coefficient, so a wetting front
        advances one patch per sub-step (a documented limit; the film that matters is stripped where it forms). A
        Picard iteration on the fully implicit form does not contract when dt x c >> 1, which is the case for micron
        films at 30 Pa on millimetre patches (measured 2026-09-20). `q_of_thickness(b)` returns the signed flux per
        patch for thickness b. Returns (m_f, n_solves, mass that arrived on another patch [kg])."""
        import scipy.sparse as sp
        import scipy.sparse.linalg as spla
        m = np.array(m_f, dtype=float)
        if not (m > 0.0).any() or len(self.i) == 0:
            return m, 0, 0.0
        moved, dt_s, n = 0.0, dt / substeps, 0
        for n in range(1, substeps + 1):
            b = m / (rho_l * areas)
            c_ij, c_ji = self.edge_coefficients(q_of_thickness(b), b, t_hat, areas)
            out = np.bincount(self.i, c_ij, self.n_patches) + np.bincount(self.j, c_ji, self.n_patches)
            rows = np.concatenate([np.arange(self.n_patches), self.j, self.i])
            cols = np.concatenate([np.arange(self.n_patches), self.i, self.j])
            vals = np.concatenate([1.0 + dt_s * out, -dt_s * c_ij, -dt_s * c_ji])
            A = sp.csc_matrix((vals, (rows, cols)), shape=(self.n_patches, self.n_patches))
            m_new = np.maximum(spla.spsolve(A, m), 0.0)
            moved += float(np.maximum(m_new - m, 0.0).sum())
            m = m_new
        return m, n, moved

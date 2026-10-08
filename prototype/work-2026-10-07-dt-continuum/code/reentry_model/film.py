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

A thick patch whose film is thinner than delta_m (amendment of 2026-10-05): the branch is decided on the liquid layer
beneath the wall (plan fact 28), so a patch is thick wherever film, molten material and deep liquid together are
deeper than delta_m, however little of it is film. The shear-driven velocity is linear within the conjugate layer,
V_s at the surface and zero at depth delta_m, and the film is the top b of that layer, so the film carries
    q = V_s b (1 - b / (2 delta_m)) + G b^3 / (3 mu_l)     for b < delta_m,
which is V_s delta_m / 2 at b = delta_m (continuous) and V_s b as b -> 0, so the runoff's emptying rate q / (A b) stays
below V_s l / A. Handing the whole layer's V_s delta_m / 2 to a vanishing film made that rate grow like 1/b, emptied
the film by a factor c dt in one implicit sub-step, squared its depth from sub-step to sub-step and broke the direct
solve (NaN film on every patch at a 0.0125 s macro step); it also moved such a film as fast as a whole conjugate layer.

The liquid below the film (amendment of 2026-10-02): where contiguous liquid lies beneath the film, the column's
runoff flux needs two depth limits in one expression -- the shear-driven part over min(delta_m, h), the pressure-
gradient and deceleration-driven part over the whole liquid depth h (plan fact 29):
    q_column = tau d^2 / (2 mu_l) + G h^3 / (3 mu_l),    d = min(delta_m, h).
The film keeps its own share above; `deep_flux` is the rest, G (h^3 - b^3) / (3 mu_l), carried by the liquid
beneath it, which the shear does not reach.

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
    says 94-100 % are thick (measured on the 100 mm physics flight, 2026-09-22). On a thick patch the film is the top
    b of the conjugate layer and carries only its own part of that layer's shear-driven flux, V_s b (1 - b/(2 delta_m))
    while b < delta_m (module docstring, amendment of 2026-10-05); V_s, the shear rate and the branch are unchanged.
    """
    tau, G, b = (np.asarray(x, dtype=float) for x in (tau, G, b))
    delta_m = np.asarray(delta_m, dtype=float)
    layer = b if b_layer is None else np.maximum(np.asarray(b_layer, dtype=float), b)
    thick = np.isfinite(delta_m) & (layer > delta_m)
    V_thin = tau * b / mu_l + G * b * b / (2.0 * mu_l)
    q_thin = tau * b * b / (2.0 * mu_l) + G * b ** 3 / (3.0 * mu_l)
    d = np.where(thick, delta_m, 0.0)
    V_thick = tau * d / mu_l
    # the film is the top b of the conjugate layer, whose velocity falls linearly from V_thick at the surface to zero
    # at depth delta_m, so it carries V b (1 - b / (2 delta_m)) while it is thinner than delta_m and the whole layer's
    # V delta_m / 2 once it fills it (amendment of 2026-10-05: the whole layer's flux handed to a vanishing film made
    # the runoff's emptying rate grow like 1/b and the transport break down)
    with np.errstate(divide="ignore", invalid="ignore"):
        q_shear = np.where(b < d, V_thick * b * (1.0 - b / (2.0 * np.where(thick, d, 1.0))), V_thick * d / 2.0)
    q_thick = q_shear + G * b ** 3 / (3.0 * mu_l)
    V = np.where(thick, V_thick, V_thin)
    q = np.where(thick, q_thick, q_thin)
    with np.errstate(divide="ignore", invalid="ignore"):
        rate = np.where(thick, V / np.where(thick, d, 1.0), np.where(b > 0.0, V / np.where(b > 0.0, b, 1.0), 0.0))
    return V, q, rate, thick


def deep_flux(G, b, h, mu_l):
    """Flux per unit width [m^2/s] of the liquid beneath the film: the pressure-gradient and deceleration-driven part
    of a column of contiguous liquid of depth h whose top b is the film,
        q_deep = G (h^3 - b^3) / (3 mu_l)            (zero where h <= b; signed like G, negative toward the nose).
    G h^3 / (3 mu_l) is the Poiseuille flux of a half-channel of depth h -- no slip at the solid, no stress from this
    part at the free surface -- and `lubrication` already gives the film G b^3 / (3 mu_l) of it, so with the film's
    shear-driven share over the conjugate depth the column carries exactly plan fact 29's
    tau d^2 / (2 mu_l) + G h^3 / (3 mu_l). There is no shear term here: the gas shear reaches only Girin's conjugate
    depth, and that layer is the film's (amendment of 2026-10-02). Since the amendment of 2026-10-05 the film carries
    only its own part of the conjugate layer's shear flux, so the column total is fact 29's where the film fills the
    conjugate layer (b >= delta_m) and short by tau (delta_m - b)^2 / (2 mu_l) where it does not: that is the shear
    flux of the part of the conjugate layer held in the elements, which no account moves. The emptying rate this flux
    gives the deep liquid, q_deep / (A h_D) <= G (b + h_D)^2 / (mu_l A) per unit edge length, stays bounded as h_D -> 0."""
    G, b, h = (np.asarray(x, dtype=float) for x in (G, b, h))
    return G * (np.maximum(h, b) ** 3 - b ** 3) / (3.0 * mu_l)


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

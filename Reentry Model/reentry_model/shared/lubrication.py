"""The melt film's lubrication velocity and flux (spec Step 3 section 8), shared by both surface models; the
runoff transport on the patch graph is `reentry_model.skyscrapers.elements.film.Runoff`.

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
beneath it, which the shear does not reach."""
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

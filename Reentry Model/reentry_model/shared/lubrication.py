"""The lubrication solution of a sheared melt film (moved from `reentry_model/film.py` for level-set stage 3, Task 1:
the physics both architectures read; the skins' `film.lubrication` is this function). Per patch, with the gas shear
tau [Pa], the driving gradient G [Pa/m], the film depth b [m], the melt boundary-layer depth delta_m [m] (Girin's
Eq. 2, `spray.melt_layer`; nan off Girin's closure) and the liquid viscosity mu_l [Pa s]:
    thin (layer <= delta_m):  V_s = tau b / mu_l + G b^2 / (2 mu_l),   q = tau b^2 / (2 mu_l) + G b^3 / (3 mu_l)
    thick (layer > delta_m):  V_s = tau delta_m / mu_l,
                              q = V_s b (1 - b / (2 delta_m)) + G b^3 / (3 mu_l) for b < delta_m, V_s delta_m / 2 + ... beyond
(q per unit width [m^2/s], signed along the surface direction away from the stagnation point). The branch is decided
on the liquid layer beneath the wall (`b_layer`), the fluxes stay on the film; the amendments of 2026-10-02 and
2026-10-05 are told in `reentry_model/film.py`'s docstring, which keeps the runoff that uses this solution."""
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

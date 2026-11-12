"""Gas loads per marker (level-set stage 3, Task 7): the angle of each marker's normal to the flow, which markers are
leeward, the fitted nose radius, and the convective heating and surface-flow result, all from the shared laws
(`shared.heating`, `shared.surface_flow`, `shared.nose`); this module holds the geometry that feeds them.

The nose sphere (design ledger row 1, Asha 2026-11-11; decision D1 of the plan): one sphere fitted, unweighted, to the
markers (which lie on phi, the gas surface) ahead of the plane through the body's mass centroid at right angles to the
flight direction, as the skyscrapers fit one sphere to their windward cap. The bound is `NOSE_CAP_FACTOR` times the
largest transverse distance (from the flight axis through the mass centroid) of those points; the radius the loads read
comes from `shared.nose.nose_radius`. A fit whose centred points are coplanar to rounding is a flat face, and a fit whose
centre lies on the flow side of the points' mean is concave (as `SkinBody._fit_nose` decides them).

Leeward: `shadow` (default), the stage 1 ray cast on each marker's patch of phi, or `centroid`, the skyscrapers'
compatibility rule (behind the centroid's plane)."""
from collections import namedtuple

import numpy as np

from . import shadow as shadow_mod
from ..shared import nose

# numerical settings
SIZE_FEEDBACK_NAMES = ("current", "initial")   # as reentry_model/body.py (not imported: the level set never reads body.py)
LEEWARD_NAMES = ("shadow", "centroid")
FLAT_FACE_SV_RATIO = 1e-9      # smallest/largest singular value of the centred points below which they are coplanar
                               # (SkinBody._fit_nose)

NoseFit = namedtuple("NoseFit", "fitted bound used n_points concave")


def marker_angles(markers, v_hat):
    """Angle [rad] between each marker's outward normal and the direction of motion `v_hat` (0 at the stagnation point)."""
    v = np.asarray(v_hat, float)
    v = v / np.linalg.norm(v)
    return np.arccos(np.clip(markers.normal @ v, -1.0, 1.0))


def leeward(mesh, phi, P, markers, v_hat, test="shadow", centroid=None):
    """Which markers are leeward, bool (n,). `shadow`: the ray test `shadow.shadowed(mesh, phi, P, v_hat)` once over phi's
    patches `P`, read at each marker's patch (a marker without a patch, -1, is not shadowed). `centroid`: not ahead of
    the plane through `centroid` (`shared.nose.ahead_of`)."""
    if test not in LEEWARD_NAMES:
        raise ValueError("leeward test must be one of {}, got {!r}".format(LEEWARD_NAMES, test))
    v = np.asarray(v_hat, float)
    v = v / np.linalg.norm(v)
    if test == "centroid":
        if centroid is None:
            raise ValueError("the centroid test needs the mass centroid")
        return ~nose.ahead_of(markers.x, np.asarray(centroid, float), v)
    hit = shadow_mod.shadowed(mesh, phi, P, v)
    patch = markers.patch
    return np.where(patch >= 0, hit[np.maximum(patch, 0)], False)


def fit_nose(markers, centroid, v_hat, size_feedback="current", fallback=float("nan")):
    """The nose fit (`NoseFit`: fitted radius or nan, bound, radius used, points fitted, concave) from the markers ahead of
    the plane through `centroid`; `fallback` [m] is the radius read without a fit (the body's initial radius)."""
    if size_feedback not in SIZE_FEEDBACK_NAMES:
        raise ValueError("size_feedback must be one of {}, got {!r}".format(SIZE_FEEDBACK_NAMES, size_feedback))
    v = np.asarray(v_hat, float)
    v = v / np.linalg.norm(v)
    c = np.asarray(centroid, float)
    pts = markers.x[nose.ahead_of(markers.x, c, v)]
    n, fitted, concave, bound = len(pts), float("nan"), False, float("nan")
    if n:
        X = pts - c
        r_t = np.sqrt(np.maximum(np.einsum("ij,ij->i", X, X) - (X @ v) ** 2, 0.0))
        bound = float(nose.NOSE_CAP_FACTOR * r_t.max())
    if n >= 4:
        sv = np.linalg.svd(pts - pts.mean(axis=0), compute_uv=False)
        centre, r = nose.fit_sphere(pts)
        if sv[0] > 0.0 and sv[-1] < FLAT_FACE_SV_RATIO * sv[0]:
            concave = True
        elif np.isfinite(r) and r > 0.0:
            fitted = r
            concave = bool((centre - pts.mean(axis=0)) @ v > 0.0)
    used = nose.nose_radius(fitted, bound, concave, n, fallback, size_feedback)
    return NoseFit(fitted, bound, float(used), n, concave)


class GasLoads:
    """The shared heating law `heating` and surface flow `surface_flow` applied to the markers; `rho_liquid` [kg/m3] is the
    film density the surface flow's inertial term reads (the shared Material's liquid density)."""

    def __init__(self, heating, surface_flow, rho_liquid):
        self.heating, self.surface_flow, self.rho_liquid = heating, surface_flow, float(rho_liquid)

    def evaluate(self, state, markers, T_top, theta, nose_fit):
        """(q_conv [W/m2 per marker], SurfaceFlowResult) at the aero `state`, the angles `theta` [rad] the windward loads
        read (after `shared.nose.load_angles`), the top-cell temperatures `T_top` [K] and the `NoseFit`."""
        T_top = np.asarray(T_top, float)
        q = self.heating.evaluate(state, theta, T_top, nose_fit.used).q_conv
        flow = self.surface_flow.evaluate(state, theta, nose_fit.used, self.rho_liquid, T_top)
        return q, flow

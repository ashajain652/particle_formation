"""The windward nose sphere and the leeward rules that every windward load reads (moved and extracted for level-set
stage 3, Task 1, so that both architectures read one copy): the algebraic sphere fit and its cap constants (from
`reentry_model/body.py`), the rule that picks the radius the loads read (from `SkinBody.nose_fit`), the test of which
surface points lie ahead of the plane through the mass centroid (from `SkinBody._leeward`) and the angles the windward
consumers read under that test (from `SkinBody.load_angles`). The skins call these functions; the arithmetic is the
skins' own, unchanged (sub-plan 22 Tasks 5g and 5j; Asha, 2026-11-10)."""
import math

import numpy as np

NOSE_CAP_ANGLE = 30.0            # deg: the windward cap fitted for the nose radius (depth (1 - cos 30 deg) R_t behind the front)
NOSE_CAP_FACTOR = 1.67           # a flat face of radius R_t heats like a sphere of 1.67 R_t (its stagnation velocity gradient is
                                 # ~0.6 x a sphere's of the same radius, Boison & Curtiss 1959): the cap on the fitted radius

SIZE_FEEDBACK_NAMES = ("current", "initial")   # the settings of `size_feedback` (moved here from `body.py`, which re-exports it)
FLAT_FACE_SV_RATIO = 1e-9        # smallest/largest singular value of the centred cap points below which they are coplanar: a flat face

# The leeward test (sub-plan 22 Task 5j; Asha, 2026-11-10). `facing`: a surface point is windward if its normal is within 90
# degrees of the direction of motion. `centroid`: it is windward only if it also lies ahead of the plane through the body's
# mass centroid at right angles to the direction of motion (`ahead_of`).
LEEWARD_TESTS = ("facing", "centroid")
LEEWARD_ANGLE = 0.5 * np.pi + 1e-6   # rad: the angle given a point that faces the flow but lies behind the centroid's plane,
                                     # just above 90 degrees, which the surface flow and the heating treat as leeward


def fit_sphere(points):
    """Algebraic least-squares sphere through `points` (n, 3): (centre, radius)."""
    x = np.asarray(points, dtype=float)
    A = np.column_stack([2.0 * x, np.ones(len(x))])
    p = np.linalg.lstsq(A, (x * x).sum(axis=1), rcond=None)[0]
    c = p[:3]
    return c, float(math.sqrt(max(p[3] + c @ c, 0.0)))


def fit_cap(points, v_hat):
    """The fit of a windward cap `points` (n, 3) for the flight direction `v_hat`: (fitted radius [m] or nan, concave). Fewer
    than 4 points is no fit. Points coplanar to rounding (smallest singular value of the centred points below
    FLAT_FACE_SV_RATIO of the largest) are a flat face: concave, no fitted radius. Otherwise the sphere of `fit_sphere`,
    used only if its radius is finite and positive, is concave where its centre lies on the flow side of the points' mean
    point. The decision as `SkinBody._fit_nose` made it inline until level-set stage 3 Task 7."""
    n, fitted, concave = len(points), float("nan"), False
    if n >= 4:
        sv = np.linalg.svd(points - points.mean(axis=0), compute_uv=False)
        centre, r = fit_sphere(points)
        if sv[0] > 0.0 and sv[-1] < FLAT_FACE_SV_RATIO * sv[0]:
            concave = True                               # coplanar to rounding: a flat face, which no sphere describes
        elif np.isfinite(r) and r > 0.0:
            fitted = r
            concave = bool((centre - points.mean(axis=0)) @ v_hat > 0.0)
    return fitted, concave


def nose_radius(fitted, bound, concave, n_points, fallback, size_feedback):
    """The nose radius [m] the windward loads read, from the windward cap's fit: `fitted` [m] the fitted sphere's radius
    (nan where there is no fit), `bound` [m] NOSE_CAP_FACTOR times the windward transverse radius (nan where nothing is
    windward), `concave` whether the cap is concave or coplanar, `n_points` the points fitted and `fallback` [m] the
    radius read without a fit (the body's initial radius). Under `size_feedback` "current": the bound where the cap is
    concave or coplanar (a flat face heats like a sphere of 1.67 times its radius, Boison and Curtiss 1959), the smaller
    of the fit and the bound where it is convex (Asha, 2026-11-10, option B; the fit alone where there is no bound), and
    `fallback` with fewer than 4 points or no usable value; under any other setting ("initial") `fallback` throughout.
    The rule as `SkinBody.nose_fit` applied it inline until level-set stage 3 Task 1."""
    used = fallback
    if size_feedback == "current" and n_points >= 4:     # fewer than 4 points never has a fit (SkinBody._fit_nose)
        if concave and np.isfinite(bound) and bound > 0.0:
            used = bound
        elif np.isfinite(fitted) and not concave:
            used = min(fitted, bound) if np.isfinite(bound) and bound > 0.0 else fitted
    return used


def ahead_of(points, centre, v_hat):
    """Which `points` (n, 3) lie ahead of the plane through `centre` at right angles to the flight direction `v_hat`, bool
    (n,): strictly ahead, a point on the plane is behind it (`SkinBody._leeward`'s test)."""
    return (points - centre) @ v_hat > 0.0


def load_angles(theta, windward):
    """The angles to the flow [rad] the windward consumers read under the `centroid` leeward test: `theta` where `windward`,
    and the larger of `theta` and LEEWARD_ANGLE (90 degrees plus a microradian) where not, so a point behind the
    centroid's plane gets no shear, no spray and no Rayleigh-Taylor mode whichever way it faces (Asha, 2026-11-10)."""
    return np.where(windward, theta, np.maximum(theta, LEEWARD_ANGLE))

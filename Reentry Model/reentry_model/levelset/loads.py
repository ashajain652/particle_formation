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
compatibility rule (behind the centroid's plane).

By the patches each marker stands for (stage 3 Task 11's fix, item 2; `by_patches`): the loads are evaluated at phi's
patches -- each at its own angle to the flow and leeward test, with the top temperature of the marker that owns it
(`owners`, the nearest marker) -- and each marker takes the area-weighted mean over the patches it owns, on its own area.
A marker's own normal stands poorly for a cell that wraps an edge sharper than the marker spacing: on the q6p box the
markers' own normals caught 74 per cent of the heat its projected area receives. The surface flow's fields along the
meridian (MERIDIAN_FIELDS) are averaged as vectors, each patch's along its own meridian, and the marker takes the mean
vector's part along its own meridian (in its tangent plane); categorical fields take the value covering the most of the
marker's area; the body-scale scalars are the evaluation's own. A marker that owns no patch takes the patch nearest it."""
import dataclasses
from collections import namedtuple
from types import SimpleNamespace

import numpy as np
from scipy.spatial import cKDTree

from . import shadow as shadow_mod
from ..shared import nose

MERIDIAN_FIELDS = ("tau", "tau_continuum", "tau_fm", "G", "dp_ds", "u_e", "u_eff")   # the surface flow's per-point fields
                  # that are magnitudes along the meridian away from the stagnation point (averaged as vectors)

# numerical settings
LEEWARD_NAMES = ("shadow", "centroid")

NoseFit = namedtuple("NoseFit", "fitted bound used n_points concave")


def owners(markers, points):
    """The marker that owns each point on phi (k, 3) [m]: the nearest marker. The one ownership of phi by the markers that
    the gas loads (each marker the mean over the patches it owns, `by_patches`) and the volume's groups beneath them
    (`skinbody.LevelSetSkinBody._beta_owners`) read, as the re-partition's surface samples do (stage 3 Task 11's fix,
    items 1 and 2): the skin above a patch and the group below it are one marker by construction. Returns marker
    indices (k,)."""
    _, m = cKDTree(np.asarray(markers.x, dtype=float)).query(np.asarray(points, dtype=float))
    return np.asarray(m, dtype=np.int64)


class Points:
    """Points on phi the loads are evaluated at, as markers carry them: positions `x` (k, 3) [m] and unit outward
    normals `normal` (k, 3) (phi's patches for `by_patches`)."""

    def __init__(self, x, normal):
        self.x, self.normal = np.asarray(x, dtype=float), np.asarray(normal, dtype=float)

    @property
    def n(self):
        return len(self.x)


def _meridian(normals, v_hat):
    """Unit vectors (k, 3) in each tangent plane along the meridian away from the stagnation point (the tangent part of
    -v_hat, normalised; zero at the stagnation and rear points), as the runoff's pushes read it."""
    v = np.asarray(v_hat, dtype=float) / np.linalg.norm(v_hat)
    away = -(v[None, :] - (normals @ v)[:, None] * normals)
    size = np.linalg.norm(away, axis=1)
    return np.divide(away, size[:, None], out=np.zeros_like(away), where=size[:, None] > 1e-12)


def patch_mean(markers, P, values, owner=None):
    """Each marker's area-weighted mean (n,) of per-patch `values` (P,) over the patches of phi `P` it owns (`owners`
    unless `owner` (P,) is given), finite values only; a marker that owns no patch, or none with a finite value, takes
    the value of the patch nearest it."""
    owner = owners(markers, P.centroid) if owner is None else np.asarray(owner, dtype=np.int64)
    values = np.asarray(values, dtype=float)
    ok = np.isfinite(values)
    w = np.bincount(owner, np.where(ok, P.area, 0.0), minlength=markers.n)
    s = np.bincount(owner, np.where(ok, P.area * np.where(ok, values, 0.0), 0.0), minlength=markers.n)
    out = np.divide(s, w, out=np.full(markers.n, np.nan), where=w > 0.0)
    none = w <= 0.0
    if none.any():
        _, near = cKDTree(P.centroid).query(markers.x[none])
        out[none] = values[near]
    return out


def _patch_mode(markers, P, values, owner):
    """Each marker's value (n,) of a categorical per-patch field: the one covering the most of the area it owns."""
    values = np.asarray(values)
    cats = np.unique(values)
    cover = np.stack([np.bincount(owner, np.where(values == c, P.area, 0.0), minlength=markers.n) for c in cats])
    out = cats[np.argmax(cover, axis=0)]
    none = cover.sum(axis=0) <= 0.0
    if none.any():
        _, near = cKDTree(P.centroid).query(markers.x[none])
        out[none] = values[near]
    return out


def _patch_meridian_mean(markers, P, values, owner, v_hat):
    """Each marker's magnitude along its meridian (n,) of a per-patch field along the meridian: the area-weighted mean,
    over the patches it owns, of each patch's value times its own meridian direction, its part along the marker's."""
    vec = np.asarray(values, dtype=float)[:, None] * _meridian(P.normal, v_hat)
    mean = np.stack([patch_mean(markers, P, vec[:, k], owner) for k in range(3)], axis=1)
    return np.einsum("ij,ij->i", mean, _meridian(markers.normal, v_hat))


def by_patches(loads, state, markers, P, T_top, theta, nose_fit, v_hat, owner=None):
    """The gas loads per marker by the patches of phi each owns (module docstring): `loads.evaluate` at the patches `P`
    (as `Points`) with each patch's angle `theta` (P,) [rad] (after its leeward test, `nose.load_angles`) and the top
    temperature [K] of its owner (`T_top` per marker), then each marker's area-weighted mean on its own area. Returns
    (q_conv [W/m2] per marker, the surface flow's result per marker: its per-point fields averaged, its body-scale fields
    the evaluation's own)."""
    owner = owners(markers, P.centroid) if owner is None else np.asarray(owner, dtype=np.int64)
    q_p, flow_p = loads.evaluate(state, Points(P.centroid, P.normal), np.asarray(T_top, dtype=float)[owner],
                                 np.asarray(theta, dtype=float), nose_fit)
    q = patch_mean(markers, P, q_p, owner)
    k = len(P.area)
    out = {}
    items = ([(f.name, getattr(flow_p, f.name)) for f in dataclasses.fields(flow_p)]
             if dataclasses.is_dataclass(flow_p) else list(vars(flow_p).items()))
    for name, value in items:
        if not isinstance(value, np.ndarray) or value.ndim != 1 or len(value) != k:
            continue                                     # body-scale fields: the evaluation's own
        if name in MERIDIAN_FIELDS:
            out[name] = _patch_meridian_mean(markers, P, value, owner, v_hat)
        elif value.dtype.kind in "biu":
            out[name] = _patch_mode(markers, P, value, owner)
        else:
            out[name] = patch_mean(markers, P, value, owner)
    flow = dataclasses.replace(flow_p, **out) if dataclasses.is_dataclass(flow_p) else SimpleNamespace(
        **{**vars(flow_p), **out})
    return q, flow


def leeward_patches(mesh, phi, P, v_hat, test="shadow", centroid=None):
    """Which of phi's patches `P` are leeward, bool (P,): the shadow test's ray cast, or not ahead of the plane through
    `centroid` (the `centroid` test), as `leeward` decides for markers."""
    if test not in LEEWARD_NAMES:
        raise ValueError("leeward test must be one of {}, got {!r}".format(LEEWARD_NAMES, test))
    v = np.asarray(v_hat, float)
    v = v / np.linalg.norm(v)
    if test == "centroid":
        if centroid is None:
            raise ValueError("the centroid test needs the mass centroid")
        return ~nose.ahead_of(P.centroid, np.asarray(centroid, float), v)
    return shadow_mod.shadowed(mesh, phi, P, v)


def marker_angles(markers, v_hat):
    """Angle [rad] between each marker's outward normal and the direction of motion `v_hat` (0 at the stagnation point)."""
    v = np.asarray(v_hat, float)
    v = v / np.linalg.norm(v)
    return np.arccos(np.clip(markers.normal @ v, -1.0, 1.0))


def leeward(mesh, phi, P, markers, v_hat, test="shadow", centroid=None, patch_lee=None):
    """Which markers are leeward, bool (n,). `shadow`: the ray test `shadow.shadowed(mesh, phi, P, v_hat)` once over phi's
    patches `P` (or `patch_lee`, `leeward_patches`' result, if given), read at each marker's patch (a marker without a
    patch, -1, is not shadowed). `centroid`: not ahead of the plane through `centroid` (`shared.nose.ahead_of`)."""
    if test not in LEEWARD_NAMES:
        raise ValueError("leeward test must be one of {}, got {!r}".format(LEEWARD_NAMES, test))
    v = np.asarray(v_hat, float)
    v = v / np.linalg.norm(v)
    if test == "centroid":
        if centroid is None:
            raise ValueError("the centroid test needs the mass centroid")
        return ~nose.ahead_of(markers.x, np.asarray(centroid, float), v)
    hit = shadow_mod.shadowed(mesh, phi, P, v) if patch_lee is None else np.asarray(patch_lee, bool)
    patch = markers.patch
    return np.where(patch >= 0, hit[np.maximum(patch, 0)], False)


def fit_nose(markers, centroid, v_hat, size_feedback="current", fallback=float("nan")):
    """The nose fit (`NoseFit`: fitted radius or nan, bound, radius used, points fitted, concave) from the markers ahead of
    the plane through `centroid`; `fallback` [m] is the radius read without a fit (the body's initial radius)."""
    if size_feedback not in nose.SIZE_FEEDBACK_NAMES:
        raise ValueError("size_feedback must be one of {}, got {!r}".format(nose.SIZE_FEEDBACK_NAMES, size_feedback))
    v = np.asarray(v_hat, float)
    v = v / np.linalg.norm(v)
    c = np.asarray(centroid, float)
    pts = markers.x[nose.ahead_of(markers.x, c, v)]
    n, bound = len(pts), float("nan")
    if n:
        X = pts - c
        r_t = np.sqrt(np.maximum(np.einsum("ij,ij->i", X, X) - (X @ v) ** 2, 0.0))
        bound = float(nose.NOSE_CAP_FACTOR * r_t.max())
    fitted, concave = nose.fit_cap(pts, v)
    used = nose.nose_radius(fitted, bound, concave, n, fallback, size_feedback)
    return NoseFit(fitted, bound, float(used), n, concave)


class GasLoads:
    """The shared heating law `heating` and surface flow `surface_flow` applied to the markers; `rho_liquid` [kg/m3] is the
    film density the surface flow's inertial term reads (the shared Material's liquid density)."""

    def __init__(self, heating, surface_flow, rho_liquid):
        self.heating, self.surface_flow, self.rho_liquid = heating, surface_flow, float(rho_liquid)

    def evaluate(self, state, markers, T_top, theta, nose_fit):
        """(q_conv [W/m2 per point], SurfaceFlowResult) at the aero `state` at the points `markers` (markers, or phi's
        patches as `Points` for `by_patches`; the shared laws read only the angles), the angles `theta` [rad] the windward
        loads read (after `shared.nose.load_angles`), the top-cell temperatures `T_top` [K] and the `NoseFit`."""
        T_top = np.asarray(T_top, float)
        q = self.heating.evaluate(state, theta, T_top, nose_fit.used).q_conv
        flow = self.surface_flow.evaluate(state, theta, nose_fit.used, self.rho_liquid, T_top)
        return q, flow

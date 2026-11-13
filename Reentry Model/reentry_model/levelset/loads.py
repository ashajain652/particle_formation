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

By the area of phi each marker owns (stage 3 Task 11's fix, item 2, and its round 1, C1; `by_patches`): the loads are
evaluated at the pieces of phi's patches each marker owns (`owned`, the one partition of phi: each point belongs to its
nearest marker, and a marker's area is what it owns) -- each piece at its patch's angle to the flow and leeward test,
with the top temperature of its owner -- and each marker takes the area-weighted mean over its pieces, on its own area,
so that the heat the markers take is the flux over phi's patches exactly. A marker's own normal stands poorly for a cell
that wraps an edge sharper than the marker spacing: on the q6p box the markers' own normals caught 74 per cent of the
heat its projected area receives. The surface flow's fields along the meridian (MERIDIAN_FIELDS) are averaged as
vectors, each piece's along its patch's meridian, and the marker takes the mean vector's part along its own meridian
(in its tangent plane); categorical fields take the value covering the most of the marker's area; the body-scale scalars
are the evaluation's own. A marker that owns no piece takes the patch nearest it."""
import dataclasses
from collections import namedtuple
from types import SimpleNamespace

import numpy as np
from scipy.spatial import cKDTree

from . import markers as markers_mod
from . import shadow as shadow_mod
from ..shared import nose

MERIDIAN_FIELDS = ("tau", "tau_continuum", "tau_fm", "G", "dp_ds", "u_e", "u_eff")   # the surface flow's per-point fields
                  # that are magnitudes along the meridian away from the stagnation point (averaged as vectors)

# numerical settings
LEEWARD_NAMES = ("shadow", "centroid")

NoseFit = namedtuple("NoseFit", "fitted bound used n_points concave")


def owned(markers, P):
    """The area of each of phi's patches `P` each marker owns, W (n, P) csr [m2] (`markers.ownership`, the one partition of
    phi: each point of phi belongs to its nearest marker; stage 3 Task 11's fix round 1, C1). The markers' areas are its
    row sums (`Markers.tessellate`) and the volume's groups beneath them read it (`skinbody.LevelSetSkinBody._beta_owners`),
    so that the skin above a piece of phi and the group below it are one marker by construction."""
    return markers_mod.owned(markers.x, P, markers_mod.OWN_SHARE * markers.spacing)


def _pieces(W):
    """The pieces of phi a W (`owned`) holds: (owner marker (k,), patch (k,), area (k,) [m2])."""
    C = W.tocoo()
    return np.asarray(C.row, dtype=np.int64), np.asarray(C.col, dtype=np.int64), np.asarray(C.data, dtype=float)


def _mean(markers, P, owner, patch, area, values):
    """Each marker's area-weighted mean (n,) of `values` at the pieces (`_pieces`), finite values only; a marker that owns
    no piece with a finite value takes the finite value of the piece whose patch is nearest it."""
    n = markers.n
    values = np.asarray(values, dtype=float)
    ok = np.isfinite(values)
    w = np.bincount(owner, np.where(ok, area, 0.0), minlength=n)
    s = np.bincount(owner, np.where(ok, area * np.where(ok, values, 0.0), 0.0), minlength=n)
    out = np.divide(s, w, out=np.full(n, np.nan), where=w > 0.0)
    none = w <= 0.0
    if none.any() and ok.any():
        _, near = cKDTree(P.centroid[patch[ok]]).query(markers.x[none])
        out[none] = values[ok][near]
    return out


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


def patch_mean(markers, P, values, W=None):
    """Each marker's area-weighted mean (n,) of per-patch `values` (P,) over the area of phi's patches `P` it owns (W,
    `owned`, computed if not given), finite values only; a marker that owns none with a finite value takes the finite
    value of the patch nearest it."""
    owner, patch, area = _pieces(owned(markers, P) if W is None else W)
    return _mean(markers, P, owner, patch, area, np.asarray(values, dtype=float)[patch])


def _mode(markers, P, owner, patch, area, values):
    """Each marker's value (n,) of a categorical field at the pieces: the one covering the most of the area it owns."""
    values = np.asarray(values)
    cats = np.unique(values)
    cover = np.stack([np.bincount(owner, np.where(values == c, area, 0.0), minlength=markers.n) for c in cats])
    out = cats[np.argmax(cover, axis=0)]
    none = cover.sum(axis=0) <= 0.0
    if none.any():
        _, near = cKDTree(P.centroid[patch]).query(markers.x[none])
        out[none] = values[near]
    return out


def _meridian_mean(markers, P, owner, patch, area, values, v_hat):
    """Each marker's magnitude along its meridian (n,) of a field along the meridian at the pieces: the area-weighted mean
    of each piece's value times its patch's own meridian direction, its part along the marker's."""
    vec = np.asarray(values, dtype=float)[:, None] * _meridian(P.normal[patch], v_hat)
    mean = np.stack([_mean(markers, P, owner, patch, area, vec[:, k]) for k in range(3)], axis=1)
    return np.einsum("ij,ij->i", mean, _meridian(markers.normal, v_hat))


def by_patches(loads, state, markers, P, T_top, theta, nose_fit, v_hat, W=None):
    """The gas loads per marker by the area of phi each owns (module docstring): `loads.evaluate` at each piece of phi a
    marker owns (W, `owned`, computed if not given: a patch's centroid and normal with its angle `theta` (P,) [rad],
    after its leeward test, `nose.load_angles`, and the top temperature [K] of the piece's owner, `T_top` per marker),
    then each marker's area-weighted mean over its pieces, on its own area. Returns (q_conv [W/m2] per marker, the surface
    flow's result per marker: its per-point fields averaged, its body-scale fields the evaluation's own)."""
    owner, patch, area = _pieces(owned(markers, P) if W is None else W)
    q_k, flow_k = loads.evaluate(state, Points(P.centroid[patch], P.normal[patch]),
                                 np.asarray(T_top, dtype=float)[owner], np.asarray(theta, dtype=float)[patch], nose_fit)
    q = _mean(markers, P, owner, patch, area, q_k)
    k = len(patch)
    out = {}
    items = ([(f.name, getattr(flow_k, f.name)) for f in dataclasses.fields(flow_k)]
             if dataclasses.is_dataclass(flow_k) else list(vars(flow_k).items()))
    for name, value in items:
        if not isinstance(value, np.ndarray) or value.ndim != 1 or len(value) != k:
            continue                                     # body-scale fields: the evaluation's own
        if name in MERIDIAN_FIELDS:
            out[name] = _meridian_mean(markers, P, owner, patch, area, value, v_hat)
        elif value.dtype.kind in "biu":
            out[name] = _mode(markers, P, owner, patch, area, value)
        else:
            out[name] = _mean(markers, P, owner, patch, area, value)
    flow = dataclasses.replace(flow_k, **out) if dataclasses.is_dataclass(flow_k) else SimpleNamespace(
        **{**vars(flow_k), **out})
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

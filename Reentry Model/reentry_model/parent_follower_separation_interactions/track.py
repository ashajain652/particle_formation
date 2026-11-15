"""Track: a time-ordered, interpolated path with its air-relative speed, path length and flow state.

Both bodies of the overlap engine are Tracks. The parent's Track carries the flow quantities that decide the plume
(density, temperature, pressure, viscosity, mean free path, parent diameter, C_D A, equal-area radius) so that
`emission_state(tau)` can hand `laws.py` the state at the moment the plume was emitted; the follower's Track needs
only positions. Nothing here takes or knows an altitude, and nothing is imported from the trajectory, droplet, body,
aero or atmosphere modules: a Track is built from arrays or from a callable (design spec section 4.0).

What a Track promises
  * Construction validates every array (times, positions, velocities, speeds, flow). Times must be finite and
    strictly increasing; everything else must be finite (speeds also non-negative). The ValueError names the array
    and the first bad index, so a NaN or a duplicate time never becomes a silently wrong spline, and `laws.classify`
    (which does not test for NaN) never receives one from a Track.
  * Positions are a cubic Hermite spline when velocities are given (positions and velocities both honoured at the
    nodes) and a cubic spline with not-a-knot ends otherwise. Velocity is the spline's derivative.
  * Queries outside [t0, t1] raise ValueError instead of extrapolating a cubic (round-off of 1e-9 of the span is
    clamped). Scalar and array times are both accepted; a scalar returns a plain float or a (3,) array.
  * Speed is air-relative. By default it is the magnitude of the spline velocity (ground speed); an explicit `speed`
    array (for example the parent's speed relative to a moving wake gas) is interpolated linearly between nodes.
    Flow arrays are also interpolated linearly, so a density or temperature never overshoots its neighbours.
  * `path_length(t)` is the integral of speed from t0, by 5-point Gauss-Legendre on each interval (exact for the
    linear speed of an explicit array and accurate to well below a millimetre on a smooth orbit arc), so it is
    non-decreasing. Interval totals are stored at construction; a call evaluates one interval's partial integral.

Units: SI throughout (seconds, metres, metres per second).
"""
import math

import numpy as np
from scipy.interpolate import CubicHermiteSpline, CubicSpline

from reentry_model.parent_follower_separation_interactions.laws import EmissionState

FLOW_KEYS = ("rho", "T", "p", "mu", "mfp", "d", "cda", "r0")   # the EmissionState fields other than U
_GL_X, _GL_W = np.polynomial.legendre.leggauss(5)
_RANGE_TOL = 1e-9          # of the track's time span: round-off allowed outside [t0, t1], clamped
_MAX_SAMPLES = 50_000_000  # guard against a spacing small enough to exhaust memory (3 doubles per sample)


def _first_non_finite(name, a):
    """Raise ValueError naming `name` and the first index along axis 0 where `a` holds a NaN or infinity."""
    bad = ~np.isfinite(a)
    if bad.any():
        idx = np.argwhere(bad)[0]
        where = f"index {idx[0]}" + (f", component {idx[1]}" if len(idx) > 1 else "")
        raise ValueError(f"{name} has a non-finite value ({a[tuple(idx)]}) at {where}")


def _as_array(name, a, shape):
    """Float array of the required shape, or ValueError naming the array."""
    try:
        out = np.array(a, dtype=float)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} cannot be converted to a float array: {exc}") from None
    if out.shape != shape:
        raise ValueError(f"{name} has shape {out.shape}, expected {shape}")
    return out


class Track:
    """An interpolated path r(t), with air-relative speed, path length and (optionally) per-sample flow state.

    Build one with `Track.from_arrays` or `Track.from_callable`; the constructor is not part of the interface.
    """

    def __init__(self, t, r, v, speed, flow):
        self._t = t
        self._r = r
        self._v = v
        self._speed_nodes = speed          # None means "magnitude of the spline velocity"
        self._flow = flow                  # None or dict key -> (N,) array
        self.t0 = float(t[0])
        self.t1 = float(t[-1])
        if v is not None:
            self._pos = CubicHermiteSpline(t, r, v, axis=0)
        else:
            self._pos = CubicSpline(t, r, axis=0)
        self._vel = self._pos.derivative()
        self._pos_fn = None                # from_callable: exact position / velocity callables
        self._vel_fn = None
        # Path length at every node: Gauss-Legendre on each interval, once, vectorised.
        a, b = t[:-1], t[1:]
        self._s_nodes = np.concatenate([[0.0], np.cumsum(self._gl(a, b - a))])

    # ---- construction ----

    @classmethod
    def from_arrays(cls, t, r, v=None, speed=None, flow=None):
        """Track through samples: t (N,), r (N, 3), optionally v (N, 3), air-relative speed (N,) and a flow dict.

        `flow` holds (N,) arrays for exactly the keys rho, T, p, mu, mfp, d, cda and r0. Raises ValueError naming
        the array and the first bad index for non-increasing times or any non-finite value.
        """
        t = np.array(t, dtype=float)
        if t.ndim != 1:
            raise ValueError(f"t has shape {t.shape}, expected (N,)")
        n = len(t)
        if n < 2:
            raise ValueError(f"a Track needs at least 2 samples, got {n}")
        _first_non_finite("t", t)
        dt = np.diff(t)
        if (dt <= 0).any():
            i = int(np.argmax(dt <= 0)) + 1
            raise ValueError(f"t is not strictly increasing at index {i} (t[{i - 1}] = {t[i - 1]!r}, t[{i}] = {t[i]!r})")
        r = _as_array("r", r, (n, 3))
        _first_non_finite("r", r)
        if v is not None:
            v = _as_array("v", v, (n, 3))
            _first_non_finite("v", v)
        if speed is not None:
            speed = _as_array("speed", speed, (n,))
            _first_non_finite("speed", speed)
            if (speed < 0).any():
                i = int(np.argmax(speed < 0))
                raise ValueError(f"speed is negative ({speed[i]!r}) at index {i}")
        if flow is not None:
            if set(flow) != set(FLOW_KEYS):
                missing, extra = sorted(set(FLOW_KEYS) - set(flow)), sorted(set(flow) - set(FLOW_KEYS))
                raise ValueError(f"flow must have exactly the keys {FLOW_KEYS}; missing {missing}, unexpected {extra}")
            flow = {k: _as_array(f"flow[{k!r}]", flow[k], (n,)) for k in FLOW_KEYS}
            for k in FLOW_KEYS:
                _first_non_finite(k, flow[k])
        return cls(t, r, v, speed, flow)

    @classmethod
    def from_callable(cls, position, t0, t1, velocity=None, speed=None, flow=None, n_samples=2001):
        """Track from `position(t)` -> (3,) on [t0, t1], sampled at n_samples evenly spaced times.

        The callable itself answers `position` queries (and `velocity` if one is given); the samples feed the
        interpolant that `samples`, `path_length` and `interpolation_error` use. `speed` is a callable returning a
        scalar air-relative speed; `flow` maps the keys rho ... r0 to callables of t (or to (n_samples,) arrays).
        """
        t0, t1 = float(t0), float(t1)
        if not (math.isfinite(t0) and math.isfinite(t1) and t1 > t0):
            raise ValueError(f"need finite t0 < t1, got t0 = {t0!r}, t1 = {t1!r}")
        if int(n_samples) < 2:
            raise ValueError(f"n_samples must be at least 2, got {n_samples}")
        t = np.linspace(t0, t1, int(n_samples))

        def sample(fn, shape_tail):
            return np.array([np.asarray(fn(x), dtype=float).reshape(shape_tail) for x in t])

        def sample_or_array(x):
            return sample(x, ()) if callable(x) else x

        r = sample(position, (3,))
        v = sample(velocity, (3,)) if velocity is not None else None
        spd = sample(speed, ()) if speed is not None else None
        fl = {k: sample_or_array(f) for k, f in flow.items()} if flow is not None else None
        tr = cls.from_arrays(t, r, v=v, speed=spd, flow=fl)
        tr._pos_fn, tr._vel_fn = position, velocity
        return tr

    # ---- internals ----

    def _prep(self, t):
        """(times as a flat float array clamped to [t0, t1], was_scalar, original shape); ValueError outside."""
        a = np.asarray(t, dtype=float)
        flat = a.reshape(-1)
        tol = _RANGE_TOL * (self.t1 - self.t0)
        if not np.all((flat >= self.t0 - tol) & (flat <= self.t1 + tol)):       # False for NaN too
            bad = flat[~((flat >= self.t0 - tol) & (flat <= self.t1 + tol))][0]
            raise ValueError(f"time {bad!r} is outside the track's range [{self.t0!r}, {self.t1!r}] (or not finite)")
        return np.clip(flat, self.t0, self.t1), a.ndim == 0, a.shape

    def _speed_flat(self, tt):
        """Air-relative speed at clamped times `tt` (any shape): explicit nodes interpolated, else |spline velocity|."""
        if self._speed_nodes is not None:
            return np.interp(tt, self._t, self._speed_nodes)
        return np.linalg.norm(self._vel(tt), axis=-1)

    def _gl(self, a, h):
        """Integral of speed over [a, a + h] (elementwise), 5-point Gauss-Legendre; a and h are same-shape arrays."""
        pts = a[:, None] + h[:, None] * (0.5 * (_GL_X + 1.0))
        sp = self._speed_flat(pts.reshape(-1)).reshape(pts.shape)
        return 0.5 * h * (sp @ _GL_W)

    # ---- queries ----

    def position(self, t):
        """Position (3,) at scalar t, or (M, 3) for an array of times."""
        tt, scalar, shape = self._prep(t)
        if self._pos_fn is not None:
            out = np.array([np.asarray(self._pos_fn(x), dtype=float).reshape(3) for x in tt])
        else:
            out = np.atleast_2d(self._pos(tt))
        return out[0] if scalar else out.reshape(shape + (3,))

    def velocity(self, t):
        """Velocity (3,) at scalar t, or (M, 3) for an array of times (the spline derivative unless one was given)."""
        tt, scalar, shape = self._prep(t)
        if self._vel_fn is not None:
            out = np.array([np.asarray(self._vel_fn(x), dtype=float).reshape(3) for x in tt])
        else:
            out = np.atleast_2d(self._vel(tt))
        return out[0] if scalar else out.reshape(shape + (3,))

    def speed(self, t):
        """Air-relative speed in m/s: a float at scalar t, an array for an array of times."""
        tt, scalar, shape = self._prep(t)
        out = self._speed_flat(tt)
        return float(out[0]) if scalar else out.reshape(shape)

    def path_length(self, t):
        """Path length in metres from t0 (integral of speed dt); non-decreasing in t."""
        tt, scalar, shape = self._prep(t)
        i = np.clip(np.searchsorted(self._t, tt, side="right") - 1, 0, len(self._t) - 2)
        a = self._t[i]
        out = self._s_nodes[i] + self._gl(a, tt - a)
        out = np.maximum(out, 0.0)
        return float(out[0]) if scalar else out.reshape(shape)

    def emission_state(self, tau):
        """`laws.EmissionState` at emission time tau: flow linearly interpolated, U = speed(tau).

        A scalar tau gives float fields; an array of times gives array fields. Raises ValueError if the track
        was built without flow.
        """
        if self._flow is None:
            raise ValueError("this Track has no flow arrays (flow=None), so it cannot give an emission state")
        tt, scalar, shape = self._prep(tau)
        vals = {k: np.interp(tt, self._t, self._flow[k]) for k in FLOW_KEYS}
        vals["U"] = self._speed_flat(tt)
        if scalar:
            return EmissionState(**{k: float(x[0]) for k, x in vals.items()})
        return EmissionState(**{k: x.reshape(shape) for k, x in vals.items()})

    def samples(self, spacing_m):
        """(tau, xyz): times and interpolated positions at roughly `spacing_m` metres of path between samples.

        Both ends of the track are included. Fully vectorised (a 1000 km track at 0.7 m spacing is 1.4 million
        points). Positions come from the interpolant even for a Track built from a callable.
        """
        sp = float(spacing_m)
        if not (math.isfinite(sp) and sp > 0):
            raise ValueError(f"spacing must be a positive finite number of metres, got {spacing_m!r}")
        total = float(self._s_nodes[-1])
        n = max(int(math.ceil(total / sp)) + 1, 2)
        if n > _MAX_SAMPLES:
            raise ValueError(f"spacing {sp} m over a {total:.0f} m track needs {n} samples (limit {_MAX_SAMPLES})")
        if total == 0.0:
            tau = np.array([self.t0, self.t1])
        else:
            target = np.linspace(0.0, total, n)
            tau = np.interp(target, self._s_nodes, self._t)          # exact if speed were constant per interval
            # Refine by Newton steps on a Hermite model of path length (S' = speed), cheap to evaluate in bulk.
            s_model = CubicHermiteSpline(self._t, self._s_nodes, self._speed_flat(self._t))
            ds = s_model.derivative()
            for _ in range(3):
                slope = np.maximum(ds(tau), 1e-12 * max(1.0, total / (self.t1 - self.t0)))
                tau = np.clip(tau - (s_model(tau) - target) / slope, self.t0, self.t1)
            tau[0], tau[-1] = self.t0, self.t1
            tau = np.maximum.accumulate(tau)
        return tau, np.atleast_2d(self._pos(tau))

    def interpolation_error(self):
        """Largest distance in metres between this interpolant and one rebuilt from every other sample.

        The even-numbered samples (and the last) are kept; the interpolants are compared at the dropped ones. It
        estimates how far the sampling is from resolving the path, not the error against the true path; a coarser
        track gives a larger value. Zero for a track of fewer than 3 samples (nothing to drop).
        """
        n = len(self._t)
        if n < 3:
            return 0.0
        keep = np.unique(np.concatenate([np.arange(0, n, 2), [n - 1]]))
        drop = np.setdiff1d(np.arange(n), keep)
        if len(drop) == 0:
            return 0.0
        if self._v is not None:
            reduced = CubicHermiteSpline(self._t[keep], self._r[keep], self._v[keep], axis=0)
        else:
            reduced = CubicSpline(self._t[keep], self._r[keep], axis=0)
        return float(np.max(np.linalg.norm(reduced(self._t[drop]) - self._pos(self._t[drop]), axis=1)))

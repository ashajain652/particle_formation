"""Overlap engine: is a follower inside the parent's plume, by how much, and over which stretches of time.

The plume is the union of balls of radius r_b centred on every point P(tau) of the parent's track, each ball growing
as the parent moves on (design spec section 4.4). A follower at position x at time t is inside if

    there is an emission time tau in [max(t0, t_r - tau_pre), t] with  |x - P(tau)| <= r_b(s(tau, t)),

where s(tau, t) is the parent's (air-relative) path length from tau to t and r_b comes from `laws.plume_radius`, with
the law chosen per emission point from the parent's own flow state (spec section 4.2). The **margin** is the largest
value over admissible tau of r_b - |x - P(tau)|: positive inside, zero on the edge, negative outside. Nothing here
takes or knows an altitude, and nothing is imported from the trajectory, droplet, body, aero or atmosphere modules
(spec section 4.0); a test inspects the imports.

The window is anchored at the follower's release time t_r (spec section 4.4): every point of track emitted from
tau_pre (60 s) before release up to now counts, however long the follower lingers. `Plume.margin(x, t, t_release)`
takes t_r, and `residence` passes its t_start as t_r on every call. Without t_release the window slides with t,
[max(t0, t - tau_pre), t]: the plume of the last tau_pre seconds, for single queries with no release in view. A
follower that stays near its release point longer than tau_pre would see a false exit with the sliding window.

Margin search (`Plume.margin`). The parent track is sampled once per Plume at a spacing of min(r0/2, 50 m) of path
(the smallest r0 on the track), widened if needed so that there are at most 200,000 samples (a 50 mm parent over
100 km would otherwise need 8 million), and the samples go into a k-d tree. The spacing only chooses candidates;
the refinement between neighbouring samples is exact whatever it is. For a query, the 16 nearest samples whose times are
admissible are taken (k doubles until at least one is admissible, or every sample has been tried). Each candidate's
neighbouring sample intervals are refined with a bounded scalar search on -(r_b - |x - P(tau)|), using the exact
`parent.position` (for a Track built from a callable this is the callable, not the spline the samples came from), and
tau = t (the parent's own position, where the ball has radius r0) is always evaluated. Before refining an interval,
an upper bound on the margin inside it is formed from the distance to the chord between its ends, the largest possible
bow of the track away from that chord (from the speed at the ends) and the larger end radius; intervals whose bound is
below the best value found at the sample points are skipped, and each run of adjacent surviving intervals gets one
bounded search (relative tolerance 1e-7 of the run). That bound assumes the radius between two neighbouring samples (a
fraction of a metre apart) is no larger than at their ends, which holds for every law at a fixed flow state.
`Plume.margin_bruteforce` is the reference: a dense grid over the whole window, its best point refined in the same
way. A test requires the two to agree to 1e-6 m.

Two limitations of the search, stated rather than hidden: it finds the best emission point among those near the
follower, so a plume that grows faster than the parent moves (dr_b/ds above 1) could have a better point farther
away; and when every nearby candidate has no law (rarefied without a fallback, or a parent slower than
`u_min_parent`) the margin is NaN even if a distant part of the track has a law. The regime then reported is that of
the nearest candidate.

Residence search (`residence`). The margin is sampled from t_start to t_end and every change of inside/outside is
bracketed and refined with `brentq` to `tol` (1 ms by default; bisection on inside/outside where a bracket end has
no law). The finest step is `scan_dt`: by default the time for the follower-parent relative speed at the start (at
least 1 m/s) to cover r0/4, clipped to [1e-4, 0.05] s. Where the margin is far from zero the step grows, to at most
`max_step` (0.05 s), to about half the time the margin needs to reach zero at the rate it can change: the smaller of
the follower-parent relative speed and the follower's own speed plus the plume's growth rate at the best emission
point. A step is accepted only if the margin changed by at most half its value, otherwise it is cut by four, so every
sign change is bracketed by a step of `scan_dt`. A fixed grid at `scan_dt` (`max_step=scan_dt`) gives the same
crossings; the adaptive step exists because a stationary follower behind a 6 km/s parent would otherwise need 2e5
margin evaluations for 20 s.

Results. `ResidenceResult` lists every inside interval (an interval still open at t_end has end None and the case is
censored), the first-interval residence time measured from t_start (0 if the follower starts outside, NaN if the start
point has no law), the total time inside, the state at each exit, and per interval a `ValidityRecord` of the regimes,
laws, Knudsen and Reynolds numbers, pressure parameter and range flags met there. Flags are recorded, never used to
refuse a computation. A follower that moves ahead of the parent from release is inside for the time its margin takes
to go negative (spec section 4.4), never "never inside".

The front plane (controller ruling "Task 6b"). The plume exists only behind the plane through the parent's current
position P(t) normal to its velocity: with proj = (x - P(t)) . v_hat(t) (positive ahead of the parent), the margin is
min(union-of-balls margin, -proj), and `MarginResult.front_limited` says when -proj is the smaller term. Without the
plane, the ball emitted at tau = t reaches r0 ahead of the parent and carries full wake gas, so a follower in a
forward-moving wake (wake_flow "both") could overtake the parent inside its own plume; the parent is not an obstacle
in this model. A follower that catches the parent would in reality hit it or be pushed aside by its shock, which is
not modelled; instead the exit through the plane is labelled (`ExitState.through == "front"`). A first interval that
ends that way has reason "reached_parent" when the follower fell behind the parent and later caught it up (the time
it reaches the parent), and "ahead_of_parent" when the exit is within 2 tol of t_start: the follower never fell behind
(a body that decelerates less than the parent, or one released on or ahead of the plane). The wake gas of
`Plume.wake_state` is ramped down to ambient over the last r0 before the plane (see there). A NaN union margin (no
law) stays NaN: the plane never makes a no-law point inside. At release the follower sits on the plane, so the start
of `residence` counts as inside when the union margin is non-negative and the follower is less than
`START_PLANE_TOL_M` (1e-6 m) ahead of the plane: round-off of a release built through geodetic coordinates (3e-9 m
for the twin), and small even for millimetre bodies.

Units: SI throughout (metres, seconds).
"""
from dataclasses import dataclass, field, fields, replace
import math

import numpy as np
from scipy.optimize import brentq, minimize_scalar
from scipy.spatial import cKDTree

from reentry_model.parent_follower_separation_interactions.laws import (EmissionState, LawSettings, RadiusResult, knudsen, plume_radius,
                                               pressure_parameter, reynolds)
from reentry_model.parent_follower_separation_interactions.track import Track
from reentry_model.parent_follower_separation_interactions.wake_gas import Ambient, GasState, ambient_state, check_mode, gas_state

K_NEAREST = 16                     # admissible samples refined per margin query (brief, Task 3)
MAX_SPACING_M = 50.0               # tree sample spacing is min(r0 / 2, this) ...
MAX_TREE_SAMPLES = 200_000         # ... widened to keep at most this many samples
SCAN_DT_RANGE = (1e-4, 0.05)       # s: clip of the default finest scan step
MAX_STEP_S = 0.05                  # s: default largest scan step
SCAN_MIN_SPEED = 1.0               # m/s: floor on the relative speed used for steps
REASONS = ("ok", "censored", "starts_outside", "rarefied_no_law", "too_slow_no_law", "reached_parent",
           "ahead_of_parent")
START_PLANE_TOL_M = 1e-6           # m: a follower this close ahead of the front plane at t_start counts as on it

_PENALTY = 1e30                    # objective value where no law applies (bounded search minimises -margin)
_REL_XTOL = 1e-7                   # bounded-search tolerance, as a fraction of the interval length
_SPEED_SLACK = 1e-4                # relative allowance on the speed bound used for the track's bow
_STATE_FIELDS = tuple(f.name for f in fields(EmissionState))


@dataclass(frozen=True)
class MarginResult:
    """The margin at one point and time, and the emission point that gives it.

    margin    r_b - distance in m, maximised over admissible tau; NaN when no admissible candidate has a radius
    tau_star  emission time of the best point, s (of the nearest candidate when the margin is NaN)
    distance  |x - P(tau_star)|, m
    s_behind  parent path length from tau_star to t, m (the distance behind the emission point the law is fed)
    r_b       plume radius there, m (NaN when no law applies)
    regime, laws, flags  from `laws.plume_radius` at tau_star
    emission  the `EmissionState` at tau_star (with the Plume's r0 override applied)
    front_limited  the front-plane term -proj, not the union of balls, set the margin
    behind_m  -proj = how far behind the parent's front plane the point is, m (negative ahead; NaN for a stopped parent)

    `margin` is min(union-of-balls margin, behind_m) (module docstring, "The front plane"); the other fields describe
    the union-of-balls best point.
    """
    margin: float
    tau_star: float
    distance: float
    s_behind: float
    r_b: float
    regime: str
    laws: tuple
    flags: dict
    emission: EmissionState
    front_limited: bool = False
    behind_m: float = math.nan


@dataclass(frozen=True)
class ExitState:
    """The follower at one exit: time [s], position (3,) [m], its air-relative speed [m/s], its lag behind the parent
    along the parent's velocity [m], its distance from the best emission point [m], the plume radius there [m] and that
    point's distance behind the parent along the track [m]; `through` is "front" when the follower left through the
    parent's front plane (it reached the parent) and "edge" when it left through the side of the plume."""
    t: float
    position: np.ndarray
    follower_speed: float
    lag_m: float
    distance_m: float
    r_b: float
    s_behind: float
    through: str = "edge"


@dataclass(frozen=True)
class ValidityRecord:
    """What held over one inside interval, from the margin evaluations in it.

    regimes and laws are sets; kn, re_d and p_param are (min, max) over the emission points (empty if none had a law);
    flags maps each `laws` flag name to whether it was raised anywhere in the interval.
    """
    interval: tuple
    regimes: set
    laws: set
    kn: tuple
    re_d: tuple
    p_param: tuple
    flags: dict


@dataclass(frozen=True)
class ResidenceResult:
    """Every inside interval of a follower and what the spec asks to be reported about them.

    intervals       [(entry, exit)], exit None for an interval still open at t_end
    residence_time  first interval's length from t_start: 0 if the follower starts outside, NaN if the start point
                    has no law, t_end - t_start if it never leaves (then a lower bound)
    total_inside    summed length of all intervals, s (open intervals counted to t_end)
    censored        the follower is still inside at t_end
    exits           an `ExitState` per closed interval
    validity        a `ValidityRecord` per interval
    interp_error_m  the follower Track's `interpolation_error()`, m; the parent's is reported separately by callers
    reason          "ok" (started inside and left through the plume's edge), "reached_parent" (started inside, fell
                    behind, and the first interval ended at the parent's front plane later than 2 tol after t_start),
                    "ahead_of_parent" (the first interval ended at the front plane within 2 tol of t_start: the
                    follower never fell behind the parent), "censored" (started inside and never left),
                    "starts_outside", "rarefied_no_law" or "too_slow_no_law" (no law at the start point)
    """
    intervals: list
    residence_time: float
    total_inside: float
    censored: bool
    exits: list
    validity: list
    interp_error_m: float
    reason: str


@dataclass
class Plume:
    """The parent's plume: its track (with flow arrays), the law settings, which bound of a law range to use, how far
    back emission counts (tau_pre, s) and an optional r0 override (metres, or a callable of the emission time tau)."""
    parent: Track
    settings: LawSettings = field(default_factory=LawSettings)
    bound: str = "low"
    tau_pre: float = 60.0
    r0: object = None
    _cache: dict = field(default_factory=dict, init=False, repr=False, compare=False)

    def __post_init__(self):
        if self.bound not in ("low", "high"):
            raise ValueError(f"bound must be 'low' or 'high', got {self.bound!r}")
        if not (math.isfinite(self.tau_pre) and self.tau_pre >= 0.0):
            raise ValueError(f"tau_pre must be a finite, non-negative number of seconds, got {self.tau_pre!r}")
        if self.r0 is not None and not callable(self.r0):
            self.r0 = _check_r0(self.r0, "r0")
        self.emission(self.parent.t0)                  # fails now, clearly, if the parent carries no flow arrays

    # ---- emission state and radius ----

    def _r0_at(self, tau):
        return _check_r0(self.r0(tau) if callable(self.r0) else self.r0, f"r0 at tau = {tau!r}")

    def emission(self, tau: float) -> EmissionState:
        """The parent's flow state at emission time tau, with the r0 override applied."""
        e = self.parent.emission_state(float(tau))
        return e if self.r0 is None else replace(e, r0=self._r0_at(float(tau)))

    def _emissions(self, taus):
        """`EmissionState` objects at an array of emission times (one vectorised flow lookup)."""
        es = self.parent.emission_state(taus)
        cols = {k: np.asarray(getattr(es, k), dtype=float).reshape(-1) for k in _STATE_FIELDS}
        if self.r0 is not None:
            cols["r0"] = np.array([self._r0_at(float(x)) for x in np.reshape(taus, -1)])
        return [EmissionState(**{k: float(cols[k][i]) for k in _STATE_FIELDS}) for i in range(len(cols["r0"]))]

    def radius(self, tau: float, t: float) -> RadiusResult:
        """Radius at time t of the ball emitted at tau, at x = parent.path_length(t) - parent.path_length(tau)."""
        tau, t = float(tau), float(t)
        if tau > t:
            raise ValueError(f"emission time tau = {tau!r} is after the query time t = {t!r}")
        s = max(self.parent.path_length(t) - self.parent.path_length(tau), 0.0)   # path length may dip by an ulp
        return plume_radius(s, self.emission(tau), self.settings, self.bound)

    # ---- the sample tree ----

    @property
    def spacing(self) -> float:
        """Tree sample spacing in metres: min(smallest r0 on the track / 2, 50 m), but at least the track's path
        length / 200,000."""
        if "spacing" not in self._cache:
            tt = np.linspace(self.parent.t0, self.parent.t1, 2001)
            r0 = (np.asarray(self.parent.emission_state(tt).r0) if self.r0 is None
                  else np.array([self._r0_at(float(x)) for x in tt]))
            total = self.parent.path_length(self.parent.t1)
            self._cache["spacing"] = float(max(min(np.min(r0) / 2.0, MAX_SPACING_M), total / MAX_TREE_SAMPLES))
        return self._cache["spacing"]

    def _index(self):
        """(sample times, k-d tree of sample positions), built once on first use."""
        if "tree" not in self._cache:
            tau_s, xyz = self.parent.samples(self.spacing)
            self._cache["tree"] = (tau_s, cKDTree(xyz))
        return self._cache["tree"]

    # ---- margin ----

    def _window(self, x, t, t_release):
        x = np.asarray(x, dtype=float).reshape(3)
        if not np.all(np.isfinite(x)):
            raise ValueError(f"follower position {x!r} is not finite")
        S_t = self.parent.path_length(float(t))          # raises ValueError outside the parent's time range
        t = min(max(float(t), self.parent.t0), self.parent.t1)
        if t_release is None:
            return x, t, S_t, max(self.parent.t0, t - self.tau_pre)
        t_release = float(t_release)
        if not (math.isfinite(t_release) and t_release <= t + 1e-9 * (self.parent.t1 - self.parent.t0)):
            raise ValueError(f"t_release = {t_release!r} must be finite and not after the query time t = {t!r}")
        return x, t, S_t, max(self.parent.t0, min(t_release, t) - self.tau_pre)

    def _values(self, x, S_t, taus):
        """At an array of emission times: positions (M, 3), distances, radii and the RadiusResults."""
        P = np.atleast_2d(self.parent.position(taus))
        s = np.maximum(S_t - np.asarray(self.parent.path_length(taus)), 0.0)
        rr = [plume_radius(float(si), e, self.settings, self.bound) for si, e in zip(s, self._emissions(taus))]
        dist = np.linalg.norm(x - P, axis=1)
        return P, dist, np.array([r.r for r in rr])

    def _objective(self, x, S_t, tau):
        """-(r_b - |x - P(tau)|) at one emission time, or a large penalty where no law applies."""
        s = max(S_t - self.parent.path_length(tau), 0.0)
        r = plume_radius(s, self.emission(tau), self.settings, self.bound).r
        if math.isnan(r):
            return _PENALTY
        return -(r - float(np.linalg.norm(x - self.parent.position(tau))))

    def _refine(self, x, S_t, a, b):
        """(tau, margin) of the best point in [a, b] by a bounded scalar search (in a local variable, so the
        tolerance is relative to the interval, not to the size of t). Margin NaN if no law was met."""
        L = b - a
        res = minimize_scalar(lambda u: self._objective(x, S_t, a + u), bounds=(0.0, L), method="bounded",
                              options={"xatol": L * _REL_XTOL})
        return a + float(res.x), (-float(res.fun) if res.fun < 0.5 * _PENALTY else math.nan)

    def _result(self, x, t, S_t, tau, margin_nan=False):
        e = self.emission(tau)
        s = max(S_t - self.parent.path_length(tau), 0.0)
        rr = plume_radius(s, e, self.settings, self.bound)
        d = float(np.linalg.norm(x - self.parent.position(tau)))
        m = math.nan if (margin_nan or math.isnan(rr.r)) else rr.r - d
        return MarginResult(m, float(tau), d, float(s), rr.r, rr.regime, rr.laws, rr.flags, e)

    def _best(self, x, t, S_t, nodes, g, dist):
        """MarginResult at the best of the evaluated nodes and refined points, or at the nearest node if none had a
        law."""
        if np.all(np.isnan(g)):
            return self._front(x, t, self._result(x, t, S_t, float(nodes[int(np.argmin(dist))]), margin_nan=True))
        return self._front(x, t, self._result(x, t, S_t, float(nodes[int(np.nanargmax(g))])))

    def _front(self, x, t, m):
        """Apply the front plane: margin = min(union margin, -proj); a NaN (no-law) margin stays NaN."""
        v = np.asarray(self.parent.velocity(t), dtype=float).reshape(3)
        speed = float(np.linalg.norm(v))
        if not speed > 0.0:                              # a stopped parent has no front
            return m
        behind = -float(np.dot(x - np.asarray(self.parent.position(t), dtype=float).reshape(3), v / speed))
        if not math.isnan(m.margin) and behind < m.margin:
            return replace(m, margin=behind, front_limited=True, behind_m=behind)
        return replace(m, behind_m=behind)

    def margin(self, x, t, t_release=None) -> MarginResult:
        """Margin of the point x (3,) at time t: the k-d tree search described in the module docstring.

        Admissible emission times are [max(t0, t_release - tau_pre), t] for a follower released at t_release (spec
        section 4.4), or [max(t0, t - tau_pre), t] (a window sliding with t) when t_release is None.
        """
        x, t, S_t, lo = self._window(x, t, t_release)
        tau_s, tree = self._index()
        n = len(tau_s)
        # k starts at half the k the last query needed (a follower far ahead of the parent needs thousands of
        # samples before one is admissible, and keeps needing them); the 16 nearest admissible are the same either way
        k = min(max(K_NEAREST, self._cache.get("k_last", 0) // 2), n)
        while True:
            _, idx = tree.query(x, k=k)
            idx = np.atleast_1d(idx)
            ok = idx[(tau_s[idx] >= lo) & (tau_s[idx] <= t)][:K_NEAREST]   # nearest first
            if ok.size or k >= n:
                break
            k = min(2 * k, n)
        self._cache["k_last"] = k
        if ok.size:
            j = np.unique(np.concatenate([ok - 1, ok]))
            j = j[(j >= 0) & (j < n - 1)]
            a, b = np.clip(tau_s[j], lo, t), np.clip(tau_s[j + 1], lo, t)
        else:                                            # no sample inside the window: refine the whole window
            a, b = np.array([lo]), np.array([t])
        keep = b > a
        a, b = a[keep], b[keep]
        nodes = np.unique(np.concatenate([a, b, [t]]))
        P, dist, r = self._values(x, S_t, nodes)
        g = r - dist
        best = np.nanmax(g) if not np.all(np.isnan(g)) else -math.inf
        cand_tau, cand_g = list(nodes), list(g)
        if a.size:
            # Upper bound on the margin inside each interval: larger end radius, minus the distance to the chord,
            # plus the farthest a curve of the interval's length (from the end speeds) can bow away from the chord.
            ia, ib = np.searchsorted(nodes, a), np.searchsorted(nodes, b)
            speed = np.linalg.norm(np.atleast_2d(self.parent.velocity(nodes)), axis=1)
            pa, chord = P[ia], P[ib] - P[ia]
            c2 = np.einsum("ij,ij->i", chord, chord)
            u = np.clip(np.einsum("ij,ij->i", x - pa, chord) / np.where(c2 > 0.0, c2, 1.0), 0.0, 1.0)
            d_seg = np.linalg.norm(x - (pa + u[:, None] * chord), axis=1)
            arc = np.maximum(speed[ia], speed[ib]) * (b - a) * (1.0 + _SPEED_SLACK)
            bow = 0.5 * np.sqrt(np.maximum(arc * arc - c2, 0.0))
            rmax = np.fmax(r[ia], r[ib])                 # NaN only when neither end has a law: never refined
            live = rmax * (1.0 + 1e-6) + 1e-9 - d_seg + bow >= best
            # Refine each run of adjacent surviving intervals with one bounded search.
            m, nm = 0, len(a)
            while m < nm:
                if not live[m]:
                    m += 1
                    continue
                e = m
                while e + 1 < nm and live[e + 1] and b[e] == a[e + 1]:
                    e += 1
                tau_r, g_r = self._refine(x, S_t, float(a[m]), float(b[e]))
                cand_tau.append(tau_r)
                cand_g.append(g_r)
                m = e + 1
        cand_d = np.concatenate([dist, np.full(len(cand_tau) - len(nodes), math.inf)])   # nearest: nodes only
        return self._best(x, t, S_t, np.array(cand_tau), np.array(cand_g), cand_d)

    def wake_state(self, x, t, ambient: Ambient, mode: str, k_E: float = 1.0,
                   t_release=None) -> tuple[GasState, MarginResult]:
        """The wake gas the follower at x sees at time t, and the margin result it came from (spec section 4.3a).

        The gas state is evaluated at the best emission point of `margin` (same window, including `t_release`): its
        distance from the follower is r_axis, its plume radius r_b, and the gas moves along the parent's unit velocity
        there. `ambient` is the undisturbed air at the follower (numbers passed in; no atmosphere module is used). The
        state is plain ambient air, with no gas velocity, when `mode` is "off", when no emission point has a law (the
        margin is NaN: a rarefied or too-slow parent has no wake gas state), or beyond the profiles' cutoff.

        Front plane (ruling "Task 6b"): the wake excess -- the gas velocity and q = 1 - T_inf/T, the fraction by which
        the gas is hotter and thinner than ambient -- is multiplied by w = smoothstep(clip(behind_m / r0, 0, 1)),
        smoothstep(z) = 3 z^2 - 2 z^3, with r0 from the emission state: exactly ambient at and ahead of the parent's
        front plane, the unramped `gas_state` from one r0 behind it, and continuous with a continuous slope between.
        """
        check_mode(mode)
        m = self.margin(x, t, t_release=t_release)
        if mode == "off" or math.isnan(m.margin) or not math.isfinite(m.r_b):
            return ambient_state(ambient), m
        v = np.asarray(self.parent.velocity(m.tau_star), dtype=float).reshape(3)
        speed = float(np.linalg.norm(v))
        if not speed > 0.0:                              # a stopped parent has no direction to push the gas along
            return ambient_state(ambient), m
        gas = gas_state(m.distance, m.r_b, m.emission, v / speed, ambient, mode, k_E=k_E)
        return _ramp(gas, ambient, m.behind_m / m.emission.r0), m

    def margin_bruteforce(self, x, t, n: int = 20001, t_release=None) -> MarginResult:
        """Reference margin (for tests): n evenly spaced emission times over the whole window, the best one then
        refined over its two neighbouring grid intervals by the same bounded search. Window as in `margin`."""
        x, t, S_t, lo = self._window(x, t, t_release)
        taus = np.linspace(lo, t, int(n)) if t > lo else np.array([t])
        _, dist, r = self._values(x, S_t, taus)
        g = r - dist
        cand_tau, cand_g = list(taus), list(g)
        if not np.all(np.isnan(g)) and len(taus) > 1:
            kb = int(np.nanargmax(g))
            tau_r, g_r = self._refine(x, S_t, float(taus[max(kb - 1, 0)]), float(taus[min(kb + 1, len(taus) - 1)]))
            cand_tau.append(tau_r)
            cand_g.append(g_r)
        cand_d = np.concatenate([dist, np.full(len(cand_tau) - len(taus), math.inf)])
        return self._best(x, t, S_t, np.array(cand_tau), np.array(cand_g), cand_d)


def _ramp(gas: GasState, ambient: Ambient, z: float) -> GasState:
    """Scale the wake excess of `gas` by smoothstep(clip(z, 0, 1)): its velocity, and q = 1 - T_inf/T (temperature
    and density follow from q at constant pressure, as in `gas_state`). w = 1 returns `gas` itself and w = 0 exactly
    ambient; a NaN z (a stopped parent, which has no front) leaves `gas` unchanged."""
    if math.isnan(z):
        return gas
    zc = min(max(z, 0.0), 1.0)
    w = zc * zc * (3.0 - 2.0 * zc)
    if w >= 1.0:
        return gas
    if w <= 0.0:
        return ambient_state(ambient)
    q = w * (1.0 - ambient.T / gas.T)
    return GasState(ambient.rho * (1.0 - q), ambient.T / (1.0 - q), gas.p, w * np.asarray(gas.u, dtype=float))


def _check_r0(value, what):
    v = float(value)
    if not (math.isfinite(v) and v > 0.0):
        raise ValueError(f"{what} must be a finite positive radius in metres, got {value!r}")
    return v


# ------------------------------------------------------------------------------------------------- residence

class _NoLaw(Exception):
    """Raised inside a root search when the margin has no law at a trial time."""


def _inside(m: MarginResult) -> bool:
    return (not math.isnan(m.margin)) and m.margin >= 0.0


def _span(values):
    v = [x for x in values if math.isfinite(x)]
    return (min(v), max(v)) if v else ()


def residence(plume: Plume, follower: Track, t_start=None, t_end=None, tol: float = 1e-3, scan_dt=None,
              max_step=None) -> ResidenceResult:
    """Every interval during which `follower` is inside `plume`, over [t_start, t_end] (default: the follower's range).

    t_start is the follower's release time: every margin is taken over emission times [max(t0, t_start - tau_pre), t]
    (spec section 4.4). Raises ValueError, naming both ranges, when the window is not inside the parent's time range
    or the follower's.
    `scan_dt` is the finest scan step and `max_step` the largest (default 0.05 s, never below scan_dt); see the module
    docstring for the search.
    """
    par = plume.parent
    t_start = follower.t0 if t_start is None else float(t_start)
    t_end = follower.t1 if t_end is None else float(t_end)
    if not (math.isfinite(t_start) and math.isfinite(t_end) and t_end > t_start):
        raise ValueError(f"need finite t_start < t_end, got t_start = {t_start!r}, t_end = {t_end!r}")
    for name, tr in (("parent", par), ("follower", follower)):
        slack = 1e-9 * (tr.t1 - tr.t0)
        if t_start < tr.t0 - slack or t_end > tr.t1 + slack:
            raise ValueError(f"the window [{t_start!r}, {t_end!r}] s is not inside the {name} track's time range "
                             f"[{tr.t0!r}, {tr.t1!r}] s")
    if not (math.isfinite(tol) and tol > 0.0):
        raise ValueError(f"tol must be a positive number of seconds, got {tol!r}")
    if scan_dt is None:
        v_rel = max(float(np.linalg.norm(follower.velocity(t_start) - par.velocity(t_start))), SCAN_MIN_SPEED)
        scan_dt = min(max(plume.emission(t_start).r0 / 4.0 / v_rel, SCAN_DT_RANGE[0]), SCAN_DT_RANGE[1])
    scan_dt = float(scan_dt)
    if not (math.isfinite(scan_dt) and scan_dt > 0.0):
        raise ValueError(f"scan_dt must be a positive number of seconds, got {scan_dt!r}")
    max_step = max(MAX_STEP_S if max_step is None else float(max_step), scan_dt)

    evals = {}                                           # t -> MarginResult, every evaluation made

    def M(tt):
        tt = float(tt)
        if tt not in evals:
            evals[tt] = plume.margin(follower.position(tt), tt, t_release=t_start)
        return evals[tt]

    def rate(tt, m):
        """How fast the margin can change, m/s (see the module docstring)."""
        vf = follower.velocity(tt)
        v_rel = float(np.linalg.norm(vf - par.velocity(tt)))
        growth = 0.0
        if math.isfinite(m.r_b):
            ds = max(1.0, 1e-3 * m.s_behind)
            r2 = plume_radius(m.s_behind + ds, m.emission, plume.settings, plume.bound).r
            if math.isfinite(r2):
                growth = max(r2 - m.r_b, 0.0) / ds * par.speed(tt)
        return max(min(v_rel, float(np.linalg.norm(vf)) + growth), SCAN_MIN_SPEED)

    def inside(tt, m):
        """Inside, with the release on the front plane: at t_start a follower whose union-of-balls margin is
        non-negative and which is less than START_PLANE_TOL_M ahead of the plane counts as inside."""
        if _inside(m):
            return True
        return (tt == t_start and m.front_limited and m.r_b - m.distance >= 0.0
                and m.behind_m >= -START_PLANE_TOL_M)

    def crossing(a, ma, b, mb):
        """Time of the change of inside/outside in [a, b], to tol."""
        ins_a = inside(a, ma)
        if not (math.isnan(ma.margin) or math.isnan(mb.margin)) and (ma.margin >= 0.0) != (mb.margin >= 0.0):
            def f(tt):
                v = M(tt).margin
                if math.isnan(v):
                    raise _NoLaw
                return v
            try:
                root = brentq(f, a, b, xtol=tol)
                M(root)
                return root
            except _NoLaw:
                pass
        while b - a > tol:                               # bisection on inside/outside
            mid = 0.5 * (a + b)
            if _inside(M(mid)) == ins_a:
                a = mid
            else:
                b = mid
        root = 0.5 * (a + b)
        M(root)
        return root

    # ---- scan ----
    t, m = t_start, M(t_start)
    crossings = []                                       # (time, entering)
    while t < t_end:
        h = max_step if math.isnan(m.margin) else min(max(0.5 * abs(m.margin) / rate(t, m), scan_dt), max_step)
        while True:
            t1 = min(t + h, t_end)
            m1 = M(t1)
            if (t1 - t) <= scan_dt * (1.0 + 1e-9) or math.isnan(m.margin):
                break
            if not math.isnan(m1.margin) and abs(m1.margin - m.margin) <= 0.5 * abs(m.margin):
                break
            h = max(h / 4.0, scan_dt)
        if inside(t, m) != _inside(m1):
            crossings.append((crossing(t, m, t1, m1), _inside(m1)))
        t, m = t1, m1

    # ---- intervals ----
    m0 = M(t_start)
    intervals, open_at = [], (t_start if inside(t_start, m0) else None)
    for tc, entering in crossings:
        if entering:
            open_at = tc
        elif open_at is not None:
            intervals.append((open_at, tc))
            open_at = None
    censored = open_at is not None
    if censored:
        intervals.append((open_at, None))

    if math.isnan(m0.margin):
        residence_time = math.nan
        reason = "too_slow_no_law" if m0.regime == "too_slow" else "rarefied_no_law"
    elif not inside(t_start, m0):
        residence_time, reason = 0.0, "starts_outside"
    else:
        first_end = intervals[0][1]
        residence_time = (t_end if first_end is None else first_end) - t_start
        if first_end is None:
            reason = "censored"
        elif M(first_end).front_limited:
            reason = "ahead_of_parent" if first_end - t_start <= 2.0 * tol else "reached_parent"
        else:
            reason = "ok"
    total_inside = sum((t_end if e is None else e) - s for s, e in intervals)

    exits = []
    for s, e in intervals:
        if e is None:
            continue
        me, x = M(e), follower.position(e)
        vp = par.velocity(e)
        vn = float(np.linalg.norm(vp))
        lag = float(np.dot(par.position(e) - x, vp / vn)) if vn > 0.0 else math.nan
        exits.append(ExitState(e, x, follower.speed(e), lag, me.distance, me.r_b, me.s_behind,
                               "front" if me.front_limited else "edge"))

    validity = []
    for s, e in intervals:
        e_eff = t_end if e is None else e
        ms = [mr for tt, mr in evals.items() if s <= tt <= e_eff and math.isfinite(mr.r_b)]
        flags = {}
        for mr in ms:
            for name, raised in mr.flags.items():
                flags[name] = flags.get(name, False) or bool(raised)
        validity.append(ValidityRecord(
            (s, e), {mr.regime for mr in ms}, {law for mr in ms for law in mr.laws},
            _span([knudsen(mr.emission) for mr in ms]), _span([reynolds(mr.emission) for mr in ms]),
            _span([pressure_parameter(mr.emission) for mr in ms]), flags))

    return ResidenceResult(intervals, residence_time, total_inside, censored, exits, validity,
                           follower.interpolation_error(), reason)

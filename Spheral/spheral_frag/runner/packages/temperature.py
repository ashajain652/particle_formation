"""ImposedTemperature: the replay's temperatures, imposed (Spheral M2 plan, Task 4 Step 1; D1, D2; spec §9.1).

Each particle's temperature comes from the body stage's table (`temperature.npy`, one row per particle, one column per
frame, sampled at the particle's reference position, D1), linear in time between the two frame columns that bracket
t; its specific thermal energy is set to u = table.enthalpy(T), the M1 table's enthalpy from 300 K (D2: imposed and
unscaled; the table's zero is 300 K, so this is D2's h(T) - h(300 K) with no offset). Spheral's own energy
evolution (compatible energy: the pair work, artificial-viscosity and plastic heating) is overwritten.

**Hooks.** u is set once when the package is built (at the interval's start time) and in `finalize` after every
step, at the step's end time t + dt. Within a CheapSynchronousRK2 step the forces are evaluated once, at the
midpoint, from the state advanced by dt/2 with the previous midpoint's derivatives; the equation of state and the
strength read u there through their state-update policies (PressurePolicy, SoundSpeedPolicy, the strength's
policies), which run inside `state.update`, before any package hook. So no hook short of recomputing those policies
can make the midpoint see h(T(t + dt/2)): `initialize` and `postStateUpdate` come after the pressure is computed.
What the midpoint sees is h(T(t)) + (dt/2) du/dt|mech, i.e. the imposed temperature half a step late plus half a
step of Spheral's own heating. `monitor=True` measures it: in `evaluateDerivatives` (at the midpoint time) the
largest |table.temperature(u) - T_imposed(t_mid)| over the local particles, kept per step in `self.monitor_rows`
(t_mid, max difference, max |T_imposed(t_mid) - T_imposed(t_start)|). Measured (2026-10-09; test module
docstring): the midpoint lags by the imposed half-step change and nothing else to within Spheral's own heating; on
the real flight (RZ 1.1 mm, interval 20 -> 21, s = 9.43e6) 1.44e-3 K, against a 9.7 K change over the interval.

**Rows.** The table's rows follow the body's `pid` array (body.npz); the local particles' rows are looked up by their
pid field at construction. Columns follow `frames_time_s`. The two bracketing columns' local rows are cached per
interval of the table, so a step reads no file. A particle with a NaN temperature at t (one the body stage removed
before t) is an error: the driver deletes particles at the boundary of their `remove_k`.

At a frame time the result is the frame's column bitwise (weight exactly 0 or 1); between frames
(1 - w) T_j + w T_{j+1}. Spheral's time is the flight time (packages/__init__.py)."""
from __future__ import annotations

import os

import numpy as np

from spheral_frag.runner.eos import get_internal, put_internal
from spheral_frag.runner.packages import LocalCount, allmax

TIME_TOL = 1e-9          # s: how far t may lie outside the table's time range (round-off of Spheral's time sum)

_CLASSES = {}


class TemperatureHistory:
    """T(pid, t) from a (n, n_frames) table, linear in time between columns. `pid` (n,) is the table's row order
    (body.npz's pid), `times` (n_frames,) strictly increasing seconds."""

    def __init__(self, temperature, times, pid):
        self.T = temperature
        self.times = np.asarray(times, dtype=np.float64)
        self.pid = np.asarray(pid, dtype=np.int64)
        if self.T.ndim != 2 or self.T.shape != (len(self.pid), len(self.times)):
            raise ValueError(f"temperature table {self.T.shape} does not match {len(self.pid)} pids and "
                             f"{len(self.times)} frame times")
        if len(self.times) < 1 or np.any(np.diff(self.times) <= 0.0):
            raise ValueError("frame times must be strictly increasing")
        self._order = np.argsort(self.pid, kind="stable")
        self._cache = None

    def rows(self, pids):
        """Table rows of the given pids (KeyError-like ValueError for an unknown pid)."""
        pids = np.asarray(pids, dtype=np.int64)
        i = np.searchsorted(self.pid, pids, sorter=self._order)
        i = np.clip(i, 0, len(self.pid) - 1)
        r = self._order[i]
        if len(pids) and not np.array_equal(self.pid[r], pids):
            raise ValueError("a particle's pid is not in the temperature table")
        return r

    def bracket(self, t):
        """(j, w): columns j, j + 1 and the weight of j + 1 at time t; w exactly 0 or 1 at a frame time."""
        times = self.times
        if t < times[0] - TIME_TOL or t > times[-1] + TIME_TOL:
            raise ValueError(f"t = {t!r} s lies outside the table's frames {times[0]!r} .. {times[-1]!r} s")
        if len(times) == 1:
            return 0, 0.0
        j = int(np.searchsorted(times, t, side="right")) - 1
        j = min(max(j, 0), len(times) - 2)
        if t <= times[j]:
            return j, 0.0
        if t >= times[j + 1]:
            return j, 1.0
        return j, (t - times[j])/(times[j + 1] - times[j])

    def columns(self, j, rows):
        """Columns j and j + 1 (j alone if there is one frame) at the rows, cached per (j, rows)."""
        c = self._cache
        if c is not None and c[0] == j and c[1] is rows:
            return c[2], c[3]
        a = np.asarray(self.T[:, j], dtype=np.float64)[rows]
        b = np.asarray(self.T[:, j + 1], dtype=np.float64)[rows] if j + 1 < self.T.shape[1] else a
        self._cache = (j, rows, a, b)
        return a, b

    def at(self, t, rows):
        """T (K) at time t for the rows."""
        j, w = self.bracket(t)
        a, b = self.columns(j, rows)
        if w == 0.0:
            return a.copy()
        if w == 1.0:
            return b.copy()
        return (1.0 - w)*a + w*b


def load_history(body_dir, mmap=True):
    """The body stage's temperature table of `body_dir` (temperature.npy memory-mapped, body.npz's pid and
    frames_time_s) as a TemperatureHistory, and the body.npz arrays."""
    with np.load(os.path.join(body_dir, "body.npz"), allow_pickle=False) as z:
        b = {k: z[k] for k in z.files}
    T = np.load(os.path.join(body_dir, "temperature.npy"), mmap_mode="r" if mmap else None)
    return TemperatureHistory(T, b["frames_time_s"], b["pid"]), b


def temperature_class(S):
    """The ImposedTemperature class for the Spheral module S."""
    if S in _CLASSES:
        return _CLASSES[S]

    class ImposedTemperature(S.Physics):
        """u = table.enthalpy(T_imposed(t)) at construction (time t0) and after every step (module docstring).

        ImposedTemperature(body, history, table, t0, monitor=False). Attributes: `temperature` (a ScalarField
        "temperature" on the node list, the imposed T for output), `rows` (the local particles' table rows),
        `monitor_rows` (list of (t_mid, max |T(u) - T_imp(t_mid)|, max |T_imp(t_mid) - T_imp(t_last)|) per
        evaluation when monitor is on, local maxima reduced over ranks by `monitor_summary`)."""

        def __init__(self, body, history, table, t0, monitor=False):
            S.Physics.__init__(self)
            self.body, self.history, self.table = body, history, table
            self.nodes = body.nodes
            self.count = LocalCount(body.nodes)
            self.rows = history.rows(body.local_pids())
            self.temperature = S.ScalarField("temperature", body.nodes)
            self.monitor = bool(monitor)
            self.monitor_rows = []
            self.t_last = float(t0)
            self.apply(float(t0))

        # -- the imposed state
        def imposed(self, t):
            T = self.history.at(t, self.rows)
            if not np.all(np.isfinite(T)):
                bad = int(np.count_nonzero(~np.isfinite(T)))
                raise RuntimeError(f"ImposedTemperature: {bad} local particles have no temperature at t = {t!r} s "
                                   "(removed by the body stage before t; delete them at their remove_k)")
            return T

        def apply(self, t):
            """Set u = h(T(t)) and the temperature field on the local particles."""
            self.count.check("ImposedTemperature")
            T = self.imposed(t)
            put_internal(self.nodes.specificThermalEnergy(), self.table.enthalpy(T))
            put_internal(self.temperature, T)
            self.t_last = float(t)
            return T

        def monitor_summary(self):
            """(n evaluations, the largest |T(u) - T_imp(t_mid)| over every rank and evaluation, the largest
            imposed half-step change) -- collective."""
            d = max((r[1] for r in self.monitor_rows), default=0.0)
            c = max((r[2] for r in self.monitor_rows), default=0.0)
            return len(self.monitor_rows), allmax(d), allmax(c)

        # -- Physics
        def evaluateDerivatives(self, t, dt, db, state, derivs):
            if not self.monitor:
                return
            self.count.check("ImposedTemperature")
            u = get_internal(self.nodes.specificThermalEnergy())
            seen = self.table.temperature(np.array(u))
            T_mid = self.imposed(float(t))
            T_start = self.history.at(self.t_last, self.rows)
            d = float(np.max(np.abs(seen - T_mid))) if len(seen) else 0.0
            c = float(np.max(np.abs(T_mid - T_start))) if len(seen) else 0.0
            self.monitor_rows.append((float(t), d, c))

        def finalize(self, t, dt, db, state, derivs):
            self.apply(float(t))

        def dt(self, db, state, derivs, t):
            return (1.0e100, "ImposedTemperature: no limit")

        def registerState(self, db, state):
            return

        def registerDerivatives(self, db, derivs):
            return

        def label(self):
            return "ImposedTemperature"

    _CLASSES[S] = ImposedTemperature
    return ImposedTemperature


def make_imposed_temperature(S, body, history, table, t0, monitor=False):
    """An ImposedTemperature for the Spheral module S (module docstring)."""
    return temperature_class(S)(body, history, table, t0, monitor)

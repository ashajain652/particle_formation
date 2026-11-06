"""The segmented trail: the droplet's trail as a row of equal-length columns, each holding the aluminium really
released there, added back together at a common time since release.

One column of length L holding all the aluminium (the single-column model, Settings.n_segments = 0) gives every point
of the trail the same age and the same concentration. In fact the droplet lays the trail over about a second, releases
most of its aluminium just after its peak temperature, and leaves a long thin tail. Here:

  * `cut_trail` cuts the path flown out to `LOSS_FRACTION` (0.999) of the aluminium loss into n equal lengths. Each
    segment holds the aluminium lost over its own length (from a monotone curve through the droplet's cumulative
    loss), and the last one also holds the remaining 0.1 %, so no aluminium is dropped. A stretch with no release makes
    no segment. Each segment carries its own emission samples (equal-mass slices of its own release), so it gets its
    own Stage A hot start, altitude and ambient air; `t_laid_s` is the aluminium-weighted mean release time.
  * Each segment runs as its own Box. It appears whole at t_laid_s, and its clock starts there; actual laying takes
    about 0.03 s near the peak of a 200 um release, so time series are not to be read in detail before the droplet
    finishes. Its volume is its length's share of the whole trail's volume (Box.volume_share), so the end correction
    erf(L / 2r) belongs to the whole trail, not to the segment, whose neighbours fill the space along the trail.
  * `TrailSolution` adds the segments at the same time since release, each read at its own age (interpolated in log
    age). It is a Solution, so the history, size-distribution and summary writers of run.py use it unchanged. Shares
    are of the droplet's whole aluminium; `unreleased(i)` is the share still in segments not yet laid.

Exchange between segments is not solved. `along_trail` adds it afterwards by spreading each segment's atoms along the
trail with the kernel that goes with the model's radius, exp(-(x - x')^2 / r^2) / (sqrt(pi) r), in closed form (a top
hat convolved with a Gaussian). That is right when the shares in each segment do not depend on concentration once the
segments overlap, and when spreading is equal along and across the trail. A test with four times less and four times
more aluminium at 70 km (2026-11-05, one case) gave identical gas ratios at 12 hours; only plume particles, nearly
finished long before the segments overlap, scale with concentration. The fully coupled alternative is kept as a
back-up in docs/milestones/2026-11-05-coupled-trail-segments.md.
"""
from dataclasses import dataclass
from types import SimpleNamespace

import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.special import erf

from reentry_model.wake_column.box import Box, Solution
from reentry_model.wake_column.slowing import HotStartResult, slow_down
from reentry_model.wake_column.species import AL_SPECIES
from reentry_model.wake_column.trail import EmissionSample, M_AL_ATOM, Trail, _at_loss
from reentry_model.wake_column.transport import volume

LOSS_FRACTION = 0.999                     # the path to this share of the aluminium loss is cut into segments
DENSE_POINTS = 200001                     # points of the cumulative-loss curve along the path


@dataclass(frozen=True)
class Segment:
    """One piece of the trail: its position [m along the path], when it appears [s after release], its Trail (L_m is
    the WHOLE trail's length; N_Al, z_mean_km and samples are its own) and its share of the whole trail's volume."""
    index: int
    x0_m: float
    x1_m: float
    t_laid_s: float
    trail: Trail
    volume_share: float


def cut_trail(rows, d0_um, h0_km, n_segments, n_samples=20, loss_fraction=LOSS_FRACTION):
    """The Segments of a droplet run from its history `rows` (module docstring), in order along the path."""
    if n_segments < 1:
        raise ValueError(f"cut_trail needs at least one segment, got {n_segments}")
    t = np.array([r["t_s"] for r in rows], dtype=float)
    h = np.array([r["h_km"] for r in rows], dtype=float)
    V = np.array([r["V_ms"] for r in rows], dtype=float)
    Td = np.array([r["T_K"] for r in rows], dtype=float)
    mAl = np.array([r["m_Al_kg"] for r in rows], dtype=float)
    lost = np.maximum.accumulate(mAl[0] - mAl)
    total = lost[-1]
    path = np.concatenate(([0.0], np.cumsum(0.5 * (V[1:] + V[:-1]) * np.diff(t))))

    x = np.linspace(0.0, path[-1], DENSE_POINTS)
    cum = np.maximum.accumulate(np.clip(PchipInterpolator(path, lost / total)(x), 0.0, 1.0))
    whole = float(_at_loss(cum, x, loss_fraction)[0])
    edges = np.linspace(0.0, whole, n_segments + 1)
    at_edge = np.interp(edges, x, cum)
    at_edge[0], at_edge[-1] = 0.0, 1.0            # the last segment holds the tail past loss_fraction

    segments = []
    for i in range(n_segments):
        share = at_edge[i + 1] - at_edge[i]
        if not share > 0.0:
            continue
        points = _at_loss(cum, x, at_edge[i] + (np.arange(n_samples) + 0.5) / n_samples * share)
        ts, hs, Vs, Ts = (np.interp(points, path, a) for a in (t, h, V, Td))
        samples = tuple(EmissionSample(float(ts[k]), float(hs[k]), float(Vs[k]), float(Ts[k]), 1.0 / n_samples)
                        for k in range(n_samples))
        mass = share * total
        trail = Trail(d0_um, h0_km, whole, float(mass / M_AL_ATOM), float(mass), 0.0, 0.0, float(hs.mean()), samples)
        segments.append(Segment(len(segments), float(edges[i]), float(edges[i + 1]), float(ts.mean()), trail,
                                float((edges[i + 1] - edges[i]) / whole)))
    return tuple(segments)


def solve_segments(segments, column, settings, species, channels, time_of_day):
    """Stage A and the box of every segment: a list of (Segment, Solution, HotStartResult)."""
    out = []
    for seg in segments:
        hot = slow_down(seg.trail.samples, column, settings, species, channels)
        amb = column.at(seg.trail.z_mean_km)
        sol = Box(seg.trail, hot, amb, time_of_day, settings, species=species,
                  volume_share=seg.volume_share).solve()
        out.append((seg, sol, hot))
    return out


def combined_hot(solved):
    """The aluminium-weighted hot start of a segmented trail, for the summary."""
    w = np.array([seg.trail.N_Al for seg, _, _ in solved], dtype=float)
    w /= w.sum()
    hots = [h for _, _, h in solved]
    mean = lambda field: float(sum(a * getattr(h, field) for a, h in zip(w, hots)))
    fractions = {k: float(sum(a * h.fractions[k] for a, h in zip(w, hots))) for k in hots[0].fractions}
    return HotStartResult(fractions=fractions, r0_m=float(np.sqrt(sum(a * h.r0_m ** 2 for a, h in zip(w, hots)))),
                          mean_collisions=mean("mean_collisions"), mean_stop_time_s=mean("mean_stop_time_s"),
                          share_above_ionisation=mean("share_above_ionisation"),
                          share_above_AlO_bond=mean("share_above_AlO_bond"), per_sample=())


def state_at_age(sol, age):
    """The state rows of a Solution at ages `age` [s], linear in log age between outputs (linear from t = 0)."""
    age = np.asarray(age, dtype=float)
    t = sol.t
    k = np.clip(np.searchsorted(t, age, side="right") - 1, 0, len(t) - 2)
    a, b = t[k], t[k + 1]
    with np.errstate(divide="ignore", invalid="ignore"):
        f = np.where(a > 0.0, np.log(np.maximum(age, a) / a) / np.log(b / a), age / b)
    f = np.clip(f, 0.0, 1.0)[:, None]
    return (1.0 - f) * sol.Y[k] + f * sol.Y[k + 1]


class TrailSolution(Solution):
    """The segments added at the output times `t` [s since release]: a Solution whose state is the sum of the
    segments' extensive state, each at its own age, whose r2 is the aluminium-weighted mean over the segments laid."""

    def __init__(self, parts, t):
        parts = [(p[0], p[1]) for p in parts]
        t = np.asarray(t, dtype=float)
        first = parts[0][1].box
        lay = first.layout
        Y = np.zeros((len(t), lay.size))
        Vs, r2, held = np.zeros(len(t)), np.zeros(len(t)), np.zeros(len(t))
        for seg, sol in parts:
            on = t >= seg.t_laid_s
            if not on.any():
                continue
            rows = state_at_age(sol, t[on] - seg.t_laid_s)
            box = sol.box
            Y[on] += rows
            Vs[on] += [box.volume_share * volume(float(np.sqrt(q)), box.L) for q in rows[:, lay.r2.start]]
            r2[on] += box.inventory * rows[:, lay.r2.start]
            held[on] += box.inventory
        Y[:, lay.r2.start] = np.divide(r2, held, out=np.zeros_like(r2), where=held > 0.0)
        total = sum(sol.box.inventory for _, sol in parts)
        warnings = sorted({str(w) for _, sol in parts for w in sol.box.warnings})
        T_K = sum(sol.box.inventory * sol.box.amb.T_K for _, sol in parts) / total
        facade = SimpleNamespace(layout=lay, bins=first.bins, inventory=total, L=first.L, volume_share=1.0,
                                 amb=SimpleNamespace(T_K=T_K), warnings=warnings)
        super().__init__(t, Y, facade)
        self._V, self._held = Vs, held

    def V(self, i):
        """Plume volume (m^3) at output i: the sum over the segments laid."""
        return float(self._V[i])

    def concentration(self, name, i):
        """Number density (m^-3) of gas species `name` at output i over the segments laid (0 before the first)."""
        return self.count(name, i) / self._V[i] if self._V[i] > 0.0 else 0.0

    def unreleased(self, i):
        """Share of the droplet's aluminium in segments not yet laid at output i."""
        return 1.0 - float(self._held[i]) / self.box.inventory


def _line_contents(sol, row):
    """{name: atoms} of one segment's state row: each aluminium gas species, plume particles and smoke uptake. Settled
    matter has left the plume and is not spread."""
    u = sol.box.layout.unpack(row)
    out = {name: float(u["gas"][j]) for j, name in enumerate(AL_SPECIES)}
    out["particles"] = float(u["tally"].sum())
    out["smoke"] = float(u["smoke"].sum())
    return out


def along_trail(parts, t_s, x_m):
    """The trail at time t_s [s since release] at positions x_m [m along the path] with spreading along the trail
    (module docstring): {"line": {name: atoms per metre}, "axis": {name: concentration on the trail's axis, m^-3}}."""
    x = np.asarray(x_m, dtype=float)
    names = list(AL_SPECIES) + ["particles", "smoke"]
    line = {n: np.zeros_like(x) for n in names}
    axis = {n: np.zeros_like(x) for n in names}
    for p in parts:
        seg, sol = p[0], p[1]
        if t_s < seg.t_laid_s:
            continue
        row = state_at_age(sol, np.array([t_s - seg.t_laid_s]))[0]
        r = float(np.sqrt(row[sol.box.layout.r2.start]))
        kernel = 0.5 * (erf((x - seg.x0_m) / r) - erf((x - seg.x1_m) / r)) / (seg.x1_m - seg.x0_m)
        for name, atoms in _line_contents(sol, row).items():
            line[name] += atoms * kernel
            axis[name] += atoms * kernel / (np.pi * r * r)
    return {"line": line, "axis": axis}

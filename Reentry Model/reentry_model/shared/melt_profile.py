"""The melt layer's speed profile as straight segments (kinked-profile dispersion design,
docs/superpowers/specs/2026-11-11-kinked-profile-dispersion-design.md, section 6.4 and stage 3), on plain arrays.

An architecture hands over its melt as `LayerColumns`: n columns of stacked cells, top first, each with its centre's
depth, thickness, speed along the top cell's direction, temperature, mass, density and viscosity, the contiguous
flowing and fully liquid runs from the top, and the depth of the rigid front (the skins' `skinbody.layer_columns`; the
level set can build the same arrays). `fit_profiles` turns each column whose top cell flows and moves into a continuous
piecewise-linear speed profile, a `dispersion.KinkedProfiles` with breakpoints at cell centres:

- the surface speed is extrapolated linearly from the top two cell centres, or with one flowing cell from the line
  through its centre and zero at the front;
- the fit passes through the surface speed, each breakpoint cell's own speed and zero at the front, and every flowing
  cell above the front lies within FIT_TOL times the top cell's speed of it, with at most N_MAX breakpoints. Breakpoints
  are inserted greedily at the cell farthest from the fit, then pruned where dropping one keeps the tolerance. [plan
  call: greedy insertion plus pruning stands in for the spec's "fewest breakpoints".]
- each segment's density is the mass-weighted density of the cells whose centres lie in its closed depth interval.

The floor test (`dispersion.place_wall`) reads `nu_below`, the kinematic viscosity of the cell just below each
breakpoint; `couette_viscosity` is the layer's effective viscosity above a depth. The physics is the spec's; the
greedy fit is the builder's (spec section 8)."""
from dataclasses import dataclass, fields

import numpy as np

from .dispersion import KinkedProfiles

FIT_TOL = 0.03      # every fitted cell within this times |the top cell's speed| of the fit; assumed until stage 6
N_MAX = 6           # breakpoints per column at most; assumed until stage 6

_FLOATS = ("z", "dz", "u", "T", "m", "rho", "mu", "front")
_BOOLS = ("flowing", "liquid")


@dataclass
class LayerColumns:
    """n columns of layered melt, top first, each array (n, K) but `front` (n,). Beyond each column's cells: dz and m
    0, mu inf, flowing and liquid False, and z, u, T and rho NaN (`skinbody.layer_columns` writes them so). The
    constructor takes the arrays as keywords, so `LayerColumns(**arrays)` builds one from a saved file's arrays."""
    z: np.ndarray           # cell-centre depth below the column's top, m
    dz: np.ndarray          # cell thickness, m
    u: np.ndarray           # each cell's speed along the top cell's direction, m/s
    T: np.ndarray           # K
    m: np.ndarray           # kg
    rho: np.ndarray         # kg/m3
    mu: np.ndarray          # Pa s at the cell's shear rate; inf where the cell does not flow
    flowing: np.ndarray     # bool: the contiguous flowing cells from the top
    liquid: np.ndarray      # bool: the contiguous fully liquid cells from the top
    front: np.ndarray       # (n,) the depth of the rigid front, m

    def __post_init__(self):
        for name in _FLOATS:
            setattr(self, name, np.asarray(getattr(self, name), dtype=float))
        for name in _BOOLS:
            setattr(self, name, np.asarray(getattr(self, name), dtype=bool))

    def rows(self, idx):
        """The columns idx (an index, an index array or a boolean mask), as a new LayerColumns of copies."""
        if np.ndim(idx) == 0:
            idx = [idx]
        return LayerColumns(**{f.name: getattr(self, f.name)[idx] for f in fields(self)})


@dataclass
class Fit:
    """The straight-segment fits of n columns. `profiles` has B = N_MAX (the fit's n_max) breakpoints at most, its wall
    the rigid front at speed 0; `cell` (n, B) is the cell at each breakpoint, -1 beyond each column's count. Columns
    that are not `reducible` have count 0, NaN depths, speeds, wall and densities, and no cells."""
    profiles: KinkedProfiles
    cell: np.ndarray


def reducible(cols):
    """(n,) bool: the columns with a profile to fit -- the top cell flows and moves, and the front lies below its
    centre. The others take today's theoretical path (spec 6.4, decision 2)."""
    return cols.flowing[:, 0] & (np.abs(cols.u[:, 0]) > 0.0) & (cols.front > cols.z[:, 0])


def _fitted(cols):
    """The cells the fit holds to its tolerance, (n, K): the flowing run's cells whose centres lie above the front."""
    z = np.where(cols.flowing, cols.z, np.inf)
    return cols.flowing & (z < cols.front[:, None])


def _surface_speed(z, u, fitted, front):
    """The surface speed (m/s), extrapolated linearly from the top two cell centres; with one fitted cell, from the line
    through its centre and zero at the front."""
    two = fitted[:, 1] if z.shape[1] > 1 else np.zeros(len(z), dtype=bool)
    z0, u0 = z[:, 0], u[:, 0]
    z1, u1 = (z[:, 1], u[:, 1]) if z.shape[1] > 1 else (z0, u0)
    gap = np.where(two, z1 - z0, front - z0)
    slope = np.where(two, u1 - u0, -u0) / gap
    return u0 - z0 * slope


def _evaluate(z, u, node, u_top, front):
    """The broken line through (0, u_top), the node cells' (z, u) and (front, 0), at every cell's centre, (n, K); z and
    u must be finite and z increasing along each row up to the last node."""
    n, K = z.shape
    j = np.arange(K)[None, :]
    rr = np.arange(n)[:, None]
    above = np.maximum.accumulate(np.where(node, j, -1), axis=1)                    # the node at or above each cell
    below = np.minimum.accumulate(np.where(node, j, K)[:, ::-1], axis=1)[:, ::-1]   # the node at or below it
    a, b = np.maximum(above, 0), np.minimum(below, K - 1)
    z_a = np.where(above >= 0, z[rr, a], 0.0)
    u_a = np.where(above >= 0, u[rr, a], u_top[:, None])
    z_b = np.where(below < K, z[rr, b], front[:, None])
    u_b = np.where(below < K, u[rr, b], 0.0)
    span = z_b - z_a
    w = np.divide(z - z_a, span, out=np.zeros_like(z), where=span > 0.0)
    return np.where(node, u, u_a + w * (u_b - u_a))


def _misfit(z, u, fitted, node, u_top, front):
    return np.where(fitted, np.abs(u - _evaluate(z, u, node, u_top, front)), 0.0)


def _greedy_then_prune(z, u, fitted, u_top, front, limit, n_max):
    """The breakpoint cells (n, K) bool: inserted one at a time at the fitted cell farthest from the fit until every
    fitted cell lies within `limit` (n,) or n_max are placed; then, while some column can drop a breakpoint and keep
    every fitted cell within its limit, it drops the one whose removal leaves the smallest misfit."""
    n, K = z.shape
    rows = np.arange(n)
    node = np.zeros((n, K), dtype=bool)
    for _ in range(n_max):
        e = _misfit(z, u, fitted, node, u_top, front)
        worst = np.argmax(e, axis=1)
        add = e[rows, worst] > limit
        if not add.any():
            break
        node[rows[add], worst[add]] = True
    while node.any():
        rank = np.cumsum(node, axis=1)
        best, drop = np.full(n, np.inf), np.full(n, -1)
        for s in range(1, int(rank[:, -1].max()) + 1):
            here = node & (rank == s)
            has = here.any(axis=1)
            cell = np.argmax(here, axis=1)
            trial = node.copy()
            trial[rows[has], cell[has]] = False
            e = _misfit(z, u, fitted, trial, u_top, front).max(axis=1)
            take = has & (e <= limit) & (e < best)
            best, drop = np.where(take, e, best), np.where(take, cell, drop)
        go = drop >= 0
        if not go.any():
            break
        node[rows[go], drop[go]] = False
    return node


def _segment_densities(cols, depth, count, front):
    """Each segment's mass-weighted density (kg/m3), (n, B + 1): the cells whose centres lie in its closed depth
    interval, from the surface (segment 0) through the breakpoints to the front (segment count); NaN beyond."""
    n, B = depth.shape
    z = np.where(np.isfinite(cols.z) & (cols.dz > 0.0), cols.z, np.nan)
    mass = np.where(np.isfinite(z), cols.m, 0.0)
    mrho = np.where(np.isfinite(z), cols.m * np.where(np.isfinite(cols.rho), cols.rho, 0.0), 0.0)
    tops = np.concatenate([np.zeros((n, 1)), depth], axis=1)
    bots = np.concatenate([depth, np.full((n, 1), np.nan)], axis=1)
    j = np.arange(B + 1)[None, :]
    bots = np.where(j == count[:, None], front[:, None], bots)
    out = np.full((n, B + 1), np.nan)
    for s in range(B + 1):
        live = s <= count
        inside = (z >= tops[:, s, None]) & (z <= bots[:, s, None]) & live[:, None]
        m_sum = np.where(inside, mass, 0.0).sum(axis=1)
        mr_sum = np.where(inside, mrho, 0.0).sum(axis=1)
        out[:, s] = np.divide(mr_sum, m_sum, out=np.full(n, np.nan), where=m_sum > 0.0)
    return out


def fit_profiles(cols, tol=FIT_TOL, n_max=N_MAX):
    """The straight-segment fit of every column (module docstring), as a `Fit` whose profiles have B = n_max: on the
    `reducible` columns, breakpoints at cell centres with their cells' own speeds, the surface speed extrapolated, the
    wall the front at speed 0 and each segment's mass-weighted density; the other columns get count 0 and NaN. The
    tolerance is tol times |the top cell's speed|, over every flowing cell whose centre lies above the front; a column
    that cannot reach it stops at n_max breakpoints."""
    n, K = cols.z.shape
    count = np.zeros(n, dtype=int)
    depth = np.full((n, n_max), np.nan)
    u_break = np.full((n, n_max), np.nan)
    u_top = np.full(n, np.nan)
    wall = np.full(n, np.nan)
    u_wall = np.full(n, np.nan)
    cell = np.full((n, n_max), -1)
    r = np.flatnonzero(reducible(cols))
    if r.size:
        sub = cols.rows(r)
        fitted = _fitted(sub)
        front = sub.front
        z = np.where(fitted, sub.z, 0.0)
        u = np.where(fitted, sub.u, 0.0)
        top = _surface_speed(z, u, fitted, front)
        limit = tol * np.abs(sub.u[:, 0])
        node = _greedy_then_prune(z, u, fitted, top, front, limit, n_max)
        c = node.sum(axis=1)
        order = np.argsort(~node, axis=1, kind="stable")[:, :n_max]      # the breakpoint cells first, top down
        if order.shape[1] < n_max:
            order = np.concatenate([order, np.zeros((len(r), n_max - order.shape[1]), dtype=int)], axis=1)
        used = np.arange(n_max)[None, :] < c[:, None]
        rr = np.arange(len(r))[:, None]
        count[r] = c
        cell[r] = np.where(used, order, -1)
        depth[r] = np.where(used, z[rr, order], np.nan)
        u_break[r] = np.where(used, u[rr, order], np.nan)
        u_top[r] = top
        wall[r] = front
        u_wall[r] = 0.0
    rho = np.full((n, n_max + 1), np.nan)
    if r.size:
        rho[r] = _segment_densities(cols.rows(r), depth[r], count[r], wall[r])
    profiles = KinkedProfiles(count=count, depth=depth, u_top=u_top, u_break=u_break, wall=wall, u_wall=u_wall, rho=rho)
    return Fit(profiles=profiles, cell=cell)


def nu_below(cols, cell):
    """(n, B) m2/s: the kinematic viscosity mu / rho of the cell just below each breakpoint's cell (`Fit.cell`), the
    least viscous part of the would-be floor (spec 6.3); inf where that cell does not flow or lies beyond the column,
    NaN where cell is -1 (padding)."""
    cell = np.asarray(cell, dtype=int)
    n, K = cols.mu.shape
    below = cell + 1
    inside = (cell >= 0) & (below < K)
    rr = np.arange(n)[:, None]
    b = np.clip(below, 0, K - 1)
    mu, rho = cols.mu[rr, b], cols.rho[rr, b]
    flows = inside & np.isfinite(mu) & np.isfinite(rho) & (rho > 0.0)
    nu = np.divide(mu, rho, out=np.full(cell.shape, np.inf), where=flows)
    return np.where(cell >= 0, nu, np.nan)


def couette_viscosity(cols, depth):
    """(nu_eff (n,) m2/s, rho (n,) kg/m3) of the cells whose centres lie above `depth` (m, (n,)): the series (Couette)
    viscosity mu_eff = sum dz / sum (dz / mu), their mass-weighted density, and nu_eff = mu_eff / rho. NaN where no
    cell lies above the depth; inf where every one of them is rigid."""
    depth = np.asarray(depth, dtype=float)
    n = cols.z.shape[0]
    above = (cols.dz > 0.0) & (cols.z < depth[:, None])
    flows = above & np.isfinite(cols.mu)
    dz = np.where(above, cols.dz, 0.0)
    thick = dz.sum(axis=1)
    compliance = (dz / np.where(flows, cols.mu, np.inf)).sum(axis=1)          # a rigid cell adds no compliance
    m = np.where(above, cols.m, 0.0)
    m_sum = m.sum(axis=1)
    rho = np.divide((m * np.where(above, cols.rho, 0.0)).sum(axis=1), m_sum, out=np.full(n, np.nan), where=m_sum > 0.0)
    mu = np.divide(thick, compliance, out=np.where(thick > 0.0, np.inf, np.nan), where=compliance > 0.0)
    nu = np.divide(mu, rho, out=np.full(n, np.nan), where=rho > 0.0)
    return nu, rho

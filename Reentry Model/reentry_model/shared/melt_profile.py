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
breakpoint; `couette_viscosity` is the layer's effective viscosity above a depth. `threshold_run` is the release run's
rule (spec 6.6, decision 5): the flowing cells from the top at least as warm as a threshold temperature, which the
skins' cascade takes at the mode's rate (`skin.mode_run`). The physics is the spec's; the greedy fit is the builder's
(spec section 8).

The mode inputs and the summary (spec 6.4-6.8): `profile_modes` puts the fit, the floor test, the broken-line solver and
the layer's properties above the wall together into what the shared spray model's profile path reads
(`spray.ProfileModes`) and what the records and the history report (`ProfileRecord`); `select_run` picks each column's
threshold temperature by the mode that took it, `run_properties` gives the droplets' properties from the run's cells,
`profile_rows` the records' appended columns (`spray.PROFILE_COLUMNS`) and `summary` the history's (`DISPERSION_COLUMNS`)."""
from dataclasses import dataclass, fields

import numpy as np

from . import spray
from .dispersion import KinkedProfiles, above_wall, fastest_waves, kinks, place_wall

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
    rho = np.full((n, n_max + 1), np.nan)
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
        rho[r] = _segment_densities(sub, depth[r], c, front)
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


def threshold_run(T, flowing, liquid, T_run):
    """Each column's release run (spec 6.6, decisions 3 and 5), (n, K) bool, top first: where its threshold temperature
    T_run [K, (n,)] is a number, the contiguous flowing cells from the top at least that warm -- liquid and slurry
    together, whatever their temperature against the liquidus -- so that the first cell colder than it or not flowing
    ends the run and a top that fails gives an empty one (a plug's top colder than the cell at its release depth: nothing
    leaves); where T_run is NaN, the run of `liquid` cells from the top, today's. T [K], flowing and liquid are (n, K);
    T may be NaN beyond a column's cells, where flowing and liquid must be False. The release is heat-limited: a cell
    joins the run by warming to the threshold, not by the number of steps (the slurry cascade's T_slurry rule)."""
    T = np.asarray(T, dtype=float)
    T_run = np.broadcast_to(np.asarray(T_run, dtype=float), T.shape[:1])
    warm = np.asarray(flowing, dtype=bool) & (T >= T_run[:, None])
    take = np.where(np.isnan(T_run)[:, None], np.asarray(liquid, dtype=bool), warm)
    return np.cumprod(take, axis=1).astype(bool)


# -- the mode inputs and the summary (spec 6.4-6.8) ---------------------------------------------------------------------
# the history's columns on the profile path (spec 6.8), per macro step, over the resolved (wet windward) columns: the
# share that are floor-like, the count that are fluid-like, the medians of the depth to the wall [mm], of the floor
# Reynolds number and of the recorded Weber number, and the mass [kg] sprayed by the gradient mode and by the gas-driven
# mode
DISPERSION_COLUMNS = ("floor_like_share", "fluid_like_columns", "wall_depth_median_mm", "floor_re_median",
                      "we_s_profile_median", "sprayed_gradient_kg", "sprayed_gas_driven_kg")


@dataclass
class ProfileRecord:
    """What the profile path reports per column (spec 6.6 and 6.8), all (n,) arrays, NaN on the columns it does not
    resolve: `wall_depth` [m], the depth to the wall; `wall_floor`, 1 where the wall is a floor-like kink and 0 where it
    is the rigid front; `floor_re`, the floor Reynolds number at the wall, or at the first kink of a fluid-like column
    (`dispersion.Wall.re_f`; NaN with no kink); `kink_depth` [m] and `kink_T` [K], the first kink's depth and its cell's
    temperature (NaN with none); `release_depth` [m], the gradient mode's (NaN where it is stable or the layer above the
    wall is one straight segment); `growth` [1/s], its growth rate (0 there); and the threshold temperatures [K] of the
    two modes' release runs: `T_gradient`, that of the cell the release depth falls in, or the gas-driven mode's where
    it reaches the wall (NaN without a release depth), and `T_driven`, that of the deepest cell whose centre lies
    strictly above the wall (NaN with none)."""
    wall_depth: np.ndarray
    wall_floor: np.ndarray
    floor_re: np.ndarray
    kink_depth: np.ndarray
    kink_T: np.ndarray
    release_depth: np.ndarray
    growth: np.ndarray
    T_gradient: np.ndarray
    T_driven: np.ndarray


def _temperature_at(cols, depth):
    """(n,) K: the temperature of the cell each depth [m, (n,)] falls in, the cells stacked from the top (cell k from the
    sum of the thicknesses above it to that sum plus its own); NaN for a NaN depth or one below the column."""
    bottom = np.cumsum(cols.dz, axis=1)
    k = (bottom <= np.asarray(depth, dtype=float)[:, None]).sum(axis=1)
    inside = np.isfinite(depth) & (k < (cols.dz > 0.0).sum(axis=1))
    return np.where(inside, cols.T[np.arange(len(k)), np.minimum(k, cols.T.shape[1] - 1)], np.nan)


def _temperature_above(cols, depth):
    """(n,) K: the temperature of the deepest cell whose centre lies strictly above each depth [m, (n,)]; NaN with none."""
    z = np.where(np.isfinite(cols.z) & (cols.dz > 0.0), cols.z, np.inf)
    count = (z < np.asarray(depth, dtype=float)[:, None]).sum(axis=1)
    return np.where(count > 0, cols.T[np.arange(len(count)), np.maximum(count - 1, 0)], np.nan)


def _run_depth(cols, T_run):
    """(n,) m: the depth at the step's start of each column's release run at the threshold T_run (`threshold_run`)."""
    return np.where(threshold_run(cols.T, cols.flowing, cols.liquid, T_run), cols.dz, 0.0).sum(axis=1)


def profile_modes(cols, windward, sigma):
    """The profile path's inputs to the shared spray model and its record (spec 6.3-6.7): (`spray.ProfileModes`,
    `ProfileRecord`) for the columns `cols` (`LayerColumns`), on the windward ones (`windward`, (n,) bool) whose profile
    is `reducible` -- the resolved columns; every other column takes today's theoretical path (spec 6.4, decision 2),
    with NaN in its modes and record but a growth of 0. At the surface tension `sigma` [N/m, a scalar or (n,)]:

    - the straight-segment fit (`fit_profiles`), the floor test with each breakpoint's `nu_below` and the profile cut
      at its wall (`dispersion.place_wall`, `dispersion.above_wall`);
    - the fastest wave on the layer above the wall (`dispersion.fastest_waves`) only where that layer has a breakpoint:
      one straight segment is stable, growth 0, so a straight film on a floor never sprays by the gradient mode (the
      spec's anti-goal);
    - the layer's Couette viscosity above the wall (`couette_viscosity`), and the recorded Weber number (spec 6.7): the
      top segment's density times the square of its velocity jump (the surface's speed less the first breakpoint's, or
      the wall's where the layer is straight) times its depth, over sigma;
    - the two modes' threshold temperatures (`ProfileRecord`) and, as the crest cap's depths, the depths at the step's
      start of the runs they give (`threshold_run`; NaN gives the liquid run's)."""
    n = cols.front.size
    sigma = np.broadcast_to(np.asarray(sigma, dtype=float), (n,))
    resolved = reducible(cols) & np.asarray(windward, dtype=bool)
    growth = np.zeros(n)
    k, wall_depth, nu_eff, we_s, cap_gradient, cap_driven = (np.full(n, np.nan) for _ in range(6))
    rec = ProfileRecord(*(np.full(n, np.nan) for _ in fields(ProfileRecord)))
    r = np.flatnonzero(resolved)
    if r.size:
        sub, m = cols.rows(r), r.size
        rows = np.arange(m)
        fit = fit_profiles(sub)
        p = fit.profiles
        wall = place_wall(p, nu_below(sub, fit.cell))
        cut = above_wall(p, wall)
        g, kk, release = np.zeros(m), np.full(m, np.nan), np.full(m, np.nan)
        curved = np.flatnonzero(cut.count > 0)
        if curved.size:
            waves = fastest_waves(cut.rows(curved), sigma[r][curved])
            g[curved], kk[curved], release[curved] = waves.growth, waves.k, waves.release_depth
        straight = cut.count == 0
        d_top = np.where(straight, cut.wall, cut.depth[:, 0])
        u_low = np.where(straight, cut.u_wall, cut.u_break[:, 0])
        is_kink = kinks(p)
        has = is_kink.any(axis=1)
        first = np.argmax(is_kink, axis=1)
        T_driven = _temperature_above(sub, wall.depth)
        T_gradient = np.where(release >= wall.depth, T_driven, _temperature_at(sub, release))
        growth[r], k[r], wall_depth[r] = g, kk, wall.depth
        nu_eff[r] = couette_viscosity(sub, wall.depth)[0]
        we_s[r] = cut.rho[:, 0] * (cut.u_top - u_low) ** 2 * d_top / sigma[r]
        cap_gradient[r], cap_driven[r] = _run_depth(sub, T_gradient), _run_depth(sub, T_driven)
        for name, value in (("wall_depth", wall.depth), ("wall_floor", np.where(wall.floor, 1.0, 0.0)),
                            ("floor_re", wall.re_f), ("kink_depth", np.where(has, p.depth[rows, first], np.nan)),
                            ("kink_T", np.where(has, sub.T[rows, np.maximum(fit.cell[rows, first], 0)], np.nan)),
                            ("release_depth", release), ("growth", g), ("T_gradient", T_gradient),
                            ("T_driven", T_driven)):
            getattr(rec, name)[r] = value
    modes = spray.ProfileModes(resolved=resolved, growth=growth, k=k, wall_depth=wall_depth, nu_eff=nu_eff, we_s=we_s,
                               cap_gradient=cap_gradient, cap_driven=cap_driven)
    return modes, rec


def select_run(branch, modes, record):
    """Each column's threshold temperature T_run [K, (n,)] for the macro step's release (spec 6.6): the gradient mode's
    (`record.T_gradient`) where a resolved column took branch 5, the gas-driven mode's (`record.T_driven`) where it took
    the thin or the rarefied branch, and NaN -- today's liquid run -- elsewhere: on the columns that are not resolved, and
    where a resolved column took no mode or the Rayleigh-Taylor mode, a liquid pool's. `branch` is the mode per column;
    a caller passes -1 where a mode releases nothing, so that only the mode that sprays fixes a run."""
    branch = np.asarray(branch)
    on = np.asarray(modes.resolved, dtype=bool)
    driven = on & np.isin(branch, (spray.BRANCH_THIN, spray.BRANCH_RAREFIED))
    return np.where(on & (branch == spray.BRANCH_PROFILE), record.T_gradient, np.where(driven, record.T_driven, np.nan))


def run_properties(cols, run):
    """The droplets' properties from the run's cells (spec 6.6), (rho [kg/m3, (n,)], mu [Pa s, (n,)]): the density and
    the viscosity (at the cells' shear rate, `cols.mu`) of the cells `run` [(n, K) bool] marks, mass-weighted; NaN for
    an empty run. [plan call: the run at the step's start stands in for the released cells.]"""
    run = np.asarray(run, dtype=bool)
    m = np.where(run, cols.m, 0.0)
    total = m.sum(axis=1)

    def weighted(x):
        return np.divide((m * np.where(run, x, 0.0)).sum(axis=1), total, out=np.full(total.size, np.nan),
                         where=total > 0.0)
    return weighted(cols.rho), weighted(cols.mu)


def profile_rows(record, delta_m, v_s):
    """(n, len(spray.PROFILE_COLUMNS)) the records' appended columns (spec 6.8), in their order: the record's wall depth
    [m], wall kind, floor Reynolds number, first kink's depth [m] and temperature [K], release depth [m] and growth rate
    [1/s], and Girin's theoretical delta_m [m] and V_s [m/s] (`delta_m`, `v_s`, (n,))."""
    n = np.shape(record.wall_depth)[0]
    return np.column_stack([record.wall_depth, record.wall_floor, record.floor_re, record.kink_depth, record.kink_T,
                            record.release_depth, record.growth, np.broadcast_to(np.asarray(delta_m, dtype=float), (n,)),
                            np.broadcast_to(np.asarray(v_s, dtype=float), (n,))]).astype(float)


def _median(x):
    x = x[np.isfinite(x)]
    return float(np.median(x)) if x.size else float("nan")


def summary(modes, record, branch, mass):
    """The history's DISPERSION_COLUMNS for a macro step (spec 6.8) over its resolved columns: the share that are
    floor-like (NaN with none), the count that are fluid-like, the medians of the wall's depth [mm], of the floor
    Reynolds number (over the columns that have one) and of the recorded Weber number (NaN with none), and the mass
    [kg] each mode sprayed: `mass` [kg, (n,)], the liquid's mode's per column, summed where a resolved column took
    branch 5 (the gradient mode) and where it took the thin or rarefied branch (the gas-driven mode). `modes` None --
    no step evaluated on the profile path yet -- is a step with none."""
    if modes is None:
        floor = wall = re = we = np.zeros(0)
        gradient = driven = 0.0
    else:
        on = np.asarray(modes.resolved, dtype=bool)
        floor, wall, re = (np.asarray(x, dtype=float)[on] for x in (record.wall_floor, record.wall_depth, record.floor_re))
        we = np.asarray(modes.we_s, dtype=float)[on]
        branch, mass = np.asarray(branch), np.asarray(mass, dtype=float)
        gradient = float(mass[on & (branch == spray.BRANCH_PROFILE)].sum())
        driven = float(mass[on & np.isin(branch, (spray.BRANCH_THIN, spray.BRANCH_RAREFIED))].sum())
    share = float(np.mean(floor == 1.0)) if floor.size else float("nan")
    return dict(zip(DISPERSION_COLUMNS, (share, int((floor == 0.0).sum()), _median(wall) * 1e3, _median(re), _median(we),
                                         gradient, driven)))

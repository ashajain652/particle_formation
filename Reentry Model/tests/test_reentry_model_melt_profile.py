"""melt_profile.py: the straight-segment fit of each melt column's speed profile (kinked-profile dispersion design,
spec 6.4 and stage 3), on plain arrays. The synthetic columns are K = 40 cells of 10 um, centres at (i + 1/2) dz, the
first 30 flowing and the rigid front at 0.3 mm unless a test says otherwise."""
import dataclasses
import os
import warnings

import numpy as np
import pytest

from reentry_model import PROJECT_ROOT
from reentry_model.shared import dispersion as dp
from reentry_model.shared import melt_profile as mp
from reentry_model.shared import rheology

K = 40
DZ = 10e-6
N_FLOW = 30
FRONT = 0.3e-3
Z = (np.arange(K) + 0.5) * DZ
FIXTURE = os.path.join(PROJECT_ROOT, "data", "reference_values", "d050_layer_columns.npz")


def columns(u, n_flow=N_FLOW, front=FRONT, rho=2400.0, mu=1.3e-3):
    """LayerColumns of synthetic columns, one per row of u (m/s, (n, K) or (K,)): the first n_flow cells (per column or
    one for all) flow at viscosity mu and the rest are rigid and at rest; rho and mu per cell or one for all; each
    cell's mass is its density times its thickness (a column of 1 m2)."""
    u = np.atleast_2d(np.asarray(u, dtype=float))
    n = len(u)
    flowing = np.arange(K)[None, :] < np.broadcast_to(np.asarray(n_flow), (n,))[:, None]
    rho = np.array(np.broadcast_to(np.asarray(rho, dtype=float), (n, K)))
    mu = np.where(flowing, np.broadcast_to(np.asarray(mu, dtype=float), (n, K)), np.inf)
    dz = np.full((n, K), DZ)
    return mp.LayerColumns(z=np.tile(Z, (n, 1)), dz=dz, u=np.where(flowing, u, 0.0), T=np.where(flowing, 950.0, 880.0),
                           m=rho * dz, rho=rho, mu=mu, flowing=flowing, liquid=flowing.copy(),
                           front=np.array(np.broadcast_to(np.asarray(front, dtype=float), (n,))))


def fit_speeds(fit, cols, i):
    """The fit of column i evaluated at its cells' centres: the surface, its breakpoints and the zero at the front."""
    p = fit.profiles
    c = p.count[i]
    x = np.concatenate([[0.0], p.depth[i, :c], [p.wall[i]]])
    y = np.concatenate([[p.u_top[i]], p.u_break[i, :c], [p.u_wall[i]]])
    return np.interp(cols.z[i], x, y)


def fitted_cells(cols):
    return cols.flowing & (np.nan_to_num(cols.z, nan=np.inf) < cols.front[:, None])


def test_the_fit_constants_are_the_specs():
    assert mp.FIT_TOL == 0.03 and mp.N_MAX == 6


def test_girin_linear_profile_gives_one_breakpoint_at_its_depth():
    """Girin's profile, u = 3 (1 - z/D) above D = 0.2 mm and at rest below it: the line from the surface speed to the
    zero at the front misses the cell at 195 um by 0.975 m/s; with a breakpoint there every cell lies within 3 % of the
    top cell's speed (the cell at 205 um is off by 0.068 m/s, against 0.088)."""
    D = 0.2e-3
    cols = columns(3.0 * np.clip(1.0 - Z / D, 0.0, None))
    fit = mp.fit_profiles(cols)
    p = fit.profiles
    assert p.depth.shape == (1, mp.N_MAX) and p.rho.shape == (1, mp.N_MAX + 1) and fit.cell.shape == (1, mp.N_MAX)
    assert p.count.tolist() == [1]
    assert abs(p.depth[0, 0] - D) <= DZ
    assert fit.cell[0, 0] == 19 and (fit.cell[0, 1:] == -1).all()
    assert p.u_top[0] == pytest.approx(3.0, rel=1e-9)
    assert p.u_break[0, 0] == cols.u[0, 19]
    assert p.wall[0] == FRONT and p.u_wall[0] == 0.0
    assert np.isnan(p.depth[0, 1:]).all() and np.isnan(p.u_break[0, 1:]).all() and np.isnan(p.rho[0, 2:]).all()


def test_a_layer_moving_as_one_piece_has_no_kink():
    """A uniform layer has one breakpoint, at the last flowing cell's centre, forced by the zero at the front; the
    segment above it does not shear, so it is no kink [plan call: the spec's "gives none" read as "no kink"]."""
    cols = columns(np.full(K, 0.5))
    fit = mp.fit_profiles(cols)
    p = fit.profiles
    assert p.count.tolist() == [1] and fit.cell[0, 0] == N_FLOW - 1
    assert p.depth[0, 0] == cols.z[0, N_FLOW - 1]
    assert p.u_top[0] == 0.5 and p.u_break[0, 0] == 0.5
    assert not dp.kinks(p).any()


def test_reducible_needs_a_moving_flowing_top_above_the_front():
    """Only a column whose top cell flows and moves, with the front below its centre, has a profile to fit (spec 6.4's
    fallback: the others take today's theoretical path)."""
    u = np.tile(np.where(Z < FRONT, 1.0, 0.0), (5, 1))
    u[1] = 0.0                                                     # top at rest
    cols = columns(u, n_flow=[N_FLOW, N_FLOW, 0, N_FLOW, N_FLOW], front=[FRONT, FRONT, 0.0, 4e-6, FRONT])
    cols.z[4] = np.nan                                             # a column with no cells
    cols.dz[4], cols.m[4] = 0.0, 0.0
    cols.flowing[4] = False
    assert mp.reducible(cols).tolist() == [True, False, False, False, False]


def test_a_plug_column_and_a_one_cell_column_fit_without_warnings():
    """Review Focus 1 and 2, under numpy's errors and warnings raised: a plug column, whose speed rises from 0.41 m/s at
    the top cell to 0.5 m/s 50 um down and falls to rest only below 100 um, fits within its tolerance with a top segment
    sloping the other way (the speed rising with depth); a one-cell column takes its surface speed from the line through
    that cell's centre and the zero at the front, and needs no breakpoint; a layer whose last flowing cell carries a kink
    has a rigid cell beneath it, so its kink has an infinite viscosity below, Re_f = 0 and a floor-like wall there; a top
    at rest and a top that does not flow give unresolved rows (count 0, NaN)."""
    plug = np.where(Z < 50e-6, 0.5 - 0.1 * (1.0 - Z / 50e-6), np.where(Z < 100e-6, 0.5, 0.5 * (FRONT - Z) / 200e-6))
    one = np.where(Z < DZ, 0.2, 0.0)
    kinked = 1.0 - 0.95 * Z / Z[N_FLOW - 1]
    u = np.vstack([plug, one, kinked, np.zeros(K), np.full(K, 0.3)])
    cols = columns(u, n_flow=[N_FLOW, 1, N_FLOW, N_FLOW, 0], front=[FRONT, 12e-6, 0.5e-3, FRONT, 0.0])
    with np.errstate(all="raise"), warnings.catch_warnings():
        warnings.simplefilter("error")
        red = mp.reducible(cols)
        fit = mp.fit_profiles(cols)
        p = fit.profiles
        nu = mp.nu_below(cols, fit.cell)
        kinks = dp.kinks(p)
        re = dp.floor_reynolds(p, nu)
        wall = dp.place_wall(p, nu)
        cut = dp.above_wall(p, wall)
        nu_eff, rho_eff = mp.couette_viscosity(cols, wall.depth)
        slopes = p.slopes()
    assert red.tolist() == [True, True, True, False, False]
    # the plug: within its tolerance, the top segment's slope negative (y up, the speed rising with depth)
    assert 1 <= p.count[0] <= mp.N_MAX
    err = np.abs(cols.u[0] - fit_speeds(fit, cols, 0))[fitted_cells(cols)[0]]
    assert err.max() <= mp.FIT_TOL * abs(cols.u[0, 0])
    assert p.u_top[0] == pytest.approx(0.40, rel=1e-12) and slopes[0, 0] < 0.0
    # the one-cell column: the surface speed on the line through (5 um, 0.2 m/s) and (12 um, 0)
    assert p.count[1] == 0 and (fit.cell[1] == -1).all()
    assert p.u_top[1] == pytest.approx(0.2 * 12.0 / 7.0, rel=1e-12)
    assert p.wall[1] == 12e-6 and p.u_wall[1] == 0.0 and np.isnan(p.rho[1, 1:]).all()
    assert p.rho[1, 0] == pytest.approx(2400.0, rel=1e-14)
    assert np.isfinite(slopes[1, 0]) and not kinks[1].any() and not wall.floor[1] and wall.depth[1] == 12e-6
    # the kink at the last flowing cell: shear 3.2e3 1/s above it and 244 1/s below, a rigid cell under it
    assert p.count[2] == 1 and fit.cell[2, 0] == N_FLOW - 1 and kinks[2, 0]
    assert nu[2, 0] == np.inf and re[2, 0] == 0.0
    assert wall.floor[2] and wall.index[2] == 0 and wall.depth[2] == Z[N_FLOW - 1] and cut.count[2] == 0
    assert nu_eff[2] == pytest.approx(1.3e-3 / 2400.0, rel=1e-13) and rho_eff[2] == pytest.approx(2400.0, rel=1e-14)
    # unresolved rows: count 0 and NaN, no breakpoint cells
    for i in (3, 4):
        assert p.count[i] == 0 and (fit.cell[i] == -1).all()
        for a in (p.depth[i], p.u_break[i], p.rho[i], [p.u_top[i], p.wall[i], p.u_wall[i]]):
            assert np.isnan(a).all()
        assert np.isnan(nu[i]).all() and np.isnan(nu_eff[i]) and np.isnan(rho_eff[i])
    for i in range(3):
        assert np.isfinite(p.u_top[i]) and np.isfinite(p.wall[i]) and np.isfinite(wall.depth[i])


def two_layers():
    """One column of two layers: ten cells of 2300 kg/m3 at 1.3 mPa s over twenty of slurry at 2500 kg/m3 and 0.5 Pa s,
    over rigid cells at 2600 kg/m3; the speed falls linearly from 3 m/s at the surface to 1 m/s at the tenth cell's
    centre (95 um), then to rest at the front."""
    top = np.arange(K) < 10
    rho = np.where(top, 2300.0, np.where(np.arange(K) < N_FLOW, 2500.0, 2600.0))
    mu = np.where(top, 1.3e-3, 0.5)
    u = np.where(Z <= Z[9], 3.0 - 2.0 * Z / Z[9], (FRONT - Z) / (FRONT - Z[9]))
    return columns(u, rho=rho, mu=mu)


def test_segment_densities_are_mass_weighted():
    """A segment's density is the mass-weighted density of the cells whose centres lie in its closed depth interval: the
    breakpoint's cell belongs to both segments it joins."""
    cols = two_layers()
    fit = mp.fit_profiles(cols)
    p = fit.profiles
    assert p.count.tolist() == [1] and fit.cell[0, 0] == 9 and p.u_top[0] == pytest.approx(3.0, rel=1e-12)
    m, rho = cols.m[0], cols.rho[0]
    lower = slice(9, N_FLOW)
    assert p.rho[0, 0] == pytest.approx(2300.0, rel=1e-14)
    expected = (m[lower] * rho[lower]).sum() / m[lower].sum()       # (2300^2 + 20 2500^2) / (2300 + 20 2500)
    assert p.rho[0, 1] == pytest.approx(expected, rel=1e-13)
    assert p.rho[0, 1] != pytest.approx(rho[lower].mean(), rel=1e-6)    # not the plain (volume) mean
    assert np.isnan(p.rho[0, 2:]).all()


def test_couette_viscosity_of_two_layers():
    """Over the cells whose centres lie above the depth: mu_eff = sum dz / sum (dz / mu), the mass-weighted density, and
    nu_eff = mu_eff / rho."""
    cols = two_layers()
    nu, rho = mp.couette_viscosity(cols, np.array([FRONT]))
    mu_eff = 300e-6 / (100e-6 / 1.3e-3 + 200e-6 / 0.5)
    rho_eff = (10 * 2300.0 ** 2 + 20 * 2500.0 ** 2) / (10 * 2300.0 + 20 * 2500.0)
    assert rho[0] == pytest.approx(rho_eff, rel=1e-13)
    assert nu[0] == pytest.approx(mu_eff / rho_eff, rel=1e-13)
    nu, rho = mp.couette_viscosity(cols, np.array([Z[10]]))           # strictly above: the top layer alone
    assert nu[0] == pytest.approx(1.3e-3 / 2300.0, rel=1e-13) and rho[0] == pytest.approx(2300.0, rel=1e-14)
    with np.errstate(all="raise"):
        nu, rho = mp.couette_viscosity(cols.rows([0, 0]), np.array([np.nan, 0.0]))
    assert np.isnan(nu).all() and np.isnan(rho).all()


def test_nu_below_reads_the_cell_under_each_breakpoint():
    """mu / rho of the cell just below each breakpoint's cell: inf where that cell is rigid or beyond the column, NaN in
    the padding."""
    cols = two_layers()
    cell = np.array([[3, 9, N_FLOW - 1, K - 1, -1, -1]])
    with np.errstate(all="raise"):
        nu = mp.nu_below(cols, cell)
    assert nu.shape == (1, 6)
    assert nu[0, 0] == 1.3e-3 / 2300.0 and nu[0, 1] == 0.5 / 2500.0
    assert nu[0, 2] == np.inf and nu[0, 3] == np.inf
    assert np.isnan(nu[0, 4:]).all()


def test_layer_columns_take_their_fields_as_keywords_and_select_rows():
    cols = two_layers().rows([0, 0, 0])
    assert cols.z.shape == (3, K) and cols.front.shape == (3,)
    again = mp.LayerColumns(**{f.name: getattr(cols, f.name) for f in dataclasses.fields(mp.LayerColumns)})
    sub = again.rows(1)
    assert sub.u.shape == (1, K) and sub.flowing.dtype == bool and np.array_equal(sub.mu, cols.mu[:1])


# -- the d050 fixture (spec stage 3) -------------------------------------------------------------------------------------
def load_fixture():
    """The fixture's LayerColumns and its other arrays (theta, tau, group), ruling R3: the npz keeps them beside the
    LayerColumns fields, and loading filters to the dataclass's own fields. `group` is 0 for a thick column, 1 for a
    plug column and 2 for another wet windward column."""
    names = {f.name for f in dataclasses.fields(mp.LayerColumns)}
    with np.load(FIXTURE) as data:
        cols = mp.LayerColumns(**{k: data[k] for k in data.files if k in names})
        extra = {k: data[k] for k in data.files if k not in names}
    return cols, extra


def isotherm_depth(cols, T_level):
    """The depth (m) where each column's contiguous cells above T_level end, interpolated between cell centres, as
    `Skins.depth_above` reads it: 0 if the top cell is not above it, the whole column if every cell is."""
    n, K_ = cols.z.shape
    valid = cols.dz > 0.0
    T = np.where(valid, cols.T, -np.inf)
    k = np.cumprod(valid & (T > T_level), axis=1).sum(axis=1)
    rows = np.arange(n)
    Ka, Kb = np.clip(k - 1, 0, K_ - 1), np.clip(k, 0, K_ - 1)
    Ta, Tb = T[rows, Ka], T[rows, Kb]
    frac = np.clip(np.divide(Ta - T_level, Ta - Tb, out=np.zeros(n), where=(Ta > Tb) & np.isfinite(Tb)), 0.0, 1.0)
    za, zb = cols.z[rows, Ka], np.where(valid[rows, Kb], cols.z[rows, Kb], cols.z[rows, Ka])
    count = valid.sum(axis=1)
    return np.where(k == 0, 0.0, np.where(k >= count, cols.dz.sum(axis=1), za + frac * (zb - za)))


def test_every_fit_holds_its_tolerance_or_its_cap_on_the_d050_fixture():
    """On the 36 columns of the d050 state at 178.5 s (`analysis/kinked_dispersion/extract_fixture.py`): every column is
    reducible, every fit holds every fitted cell within FIT_TOL of its top cell's speed or has N_MAX breakpoints, and the
    breakpoint depths increase from the surface to above the front. The plug rows (group 1) are genuine plug columns:
    their fitted top segment slopes the other way, the speed rising with depth (Review Focus 1). Printed, as a
    measurement (spec stage 3), not a pass/fail: the medians of the first kink's depth and of the wall's depth
    (`dispersion.place_wall` with `nu_below`), each over the depth of the isotherm at the bridge's base (906.4 K). The
    builder check measured 0.94 for the first kink on the 622 thick columns."""
    cols, extra = load_fixture()
    n = cols.front.size
    assert n == 36 and set(extra) == {"theta", "tau", "group"}
    assert all(a.shape == (n,) for a in extra.values())
    group = extra["group"]
    assert np.bincount(group, minlength=3).tolist() == [24, 6, 6]
    with np.errstate(all="raise"), warnings.catch_warnings():
        warnings.simplefilter("error")
        assert mp.reducible(cols).all()
        fit = mp.fit_profiles(cols)
        nu = mp.nu_below(cols, fit.cell)
        wall = dp.place_wall(fit.profiles, nu)
        kinks = dp.kinks(fit.profiles)
    p = fit.profiles
    fitted = fitted_cells(cols)
    within = np.zeros(n, dtype=bool)
    for i in range(n):
        c = p.count[i]
        err = np.abs(cols.u[i] - fit_speeds(fit, cols, i))[fitted[i]].max()
        within[i] = err <= mp.FIT_TOL * abs(cols.u[i, 0]) * (1.0 + 1e-12)
        assert within[i] or c == mp.N_MAX, i
        d = p.depth[i, :c]
        assert (np.diff(np.concatenate([[0.0], d, [p.wall[i]]])) > 0.0).all(), i
        assert np.array_equal(d, cols.z[i, fit.cell[i, :c]]) and p.wall[i] == cols.front[i] and p.u_wall[i] == 0.0
    assert (p.slopes()[group == 1, 0] < 0.0).all()
    d906 = isotherm_depth(cols, rheology.BRIDGE_LOW)
    has = kinks.any(axis=1)
    first = np.where(has, p.depth[np.arange(n), np.argmax(kinks, axis=1)], np.nan)
    ok = d906 > 0.0
    kink_ratio = first[ok & has] / d906[ok & has]
    wall_ratio = wall.depth[ok] / d906[ok]
    print("\nd050 fixture, {} columns ({} with the 906.4 K isotherm below the top cell's centre): columns by breakpoint "
          "count 0 to {} {}; {} within tolerance, the rest at the cap; {} with a kink".format(
              n, int(ok.sum()), mp.N_MAX, np.bincount(p.count, minlength=mp.N_MAX + 1).tolist(), int(within.sum()),
              int(has.sum())))
    print("  median first kink depth / 906.4 K depth: {:.3f} (p10 {:.3g}, p90 {:.3g}, {} columns)".format(
        np.median(kink_ratio), *np.percentile(kink_ratio, [10, 90]), kink_ratio.size))
    print("  median wall depth / 906.4 K depth: {:.3f} (p10 {:.3g}, p90 {:.3g}, {} columns, {} of them floor-like)".format(
        np.median(wall_ratio), *np.percentile(wall_ratio, [10, 90]), wall_ratio.size, int((wall.floor & ok).sum())))
    assert np.isfinite(kink_ratio).all() and np.isfinite(wall_ratio).all() and kink_ratio.size > 0

"""Tests of the level set's spraying per marker (design sections 9.5 and 20.6; level-set stage 3 Task 8): Girin's plan on
the markers' skins through the shared `SprayModel`, the release cascade from the skins' top cells, the droplet records
through the shared writer, and Girin and Kopyt's 1994 front-surface case through the level-set path; and the kinked-profile
spraying (plan 03b Task 10): the table path byte for byte as before, and the profile path against the skyscrapers' on the
same cells. The skyscrapers' `skin.Skins`, `cascade_mass`, `release_top` and `skinbody.girin_plan` and `layer_columns`
are imported here only to compare against (decision D4): the package itself never imports them."""
import csv
import json
import os
from types import SimpleNamespace

import numpy as np
import pytest

from reentry_model.levelset import spraying
from reentry_model.levelset.markers import Markers
from reentry_model.levelset.skin import Skins, SkinSettings
from reentry_model.shared import droplets, spray, surface_flow as sf
from reentry_model.shared.lubrication import lubrication
from reentry_model.shared.material import LiquidProperties, Material

from helpers import PROJECT_ROOT
from materials import LinearMeltMaterial

LIQUID, SLURRY, RIGID = 950.0, 902.0, 880.0    # K: all liquid; 30 % solid, flowing; rigid (AA7075's T_rigid 895.1 K)
CELL = 1.0e-5                                  # m: each stacked cell holds 10 um of solid metal's mass
SETTINGS = SkinSettings(cell=10e-6, fine_zone=0.4e-3)
H = 1.0e-3                                     # m: the markers' spacing on the flat grid
SOLID_TOL = 1.0e-9                             # the D4 test's tolerance on the solid released (its docstring)
TAU_RENEWAL = 20.0                             # Pa: the gas shear of the renewal test (supercritical: We_s above 3.08
                                               # on the thick branch)


def aa7075():
    return Material.from_drama_json("AA7075_scheil")


def row_markers(n, h=H, order=None):
    """n markers on a straight row h apart in z = 0, normals +z, square cells of area h^2; `order` the skin row of each
    marker (default the identity). Their division's boundaries are `row_division`'s."""
    x = np.column_stack([np.arange(n) * h, np.zeros(n), np.zeros(n)])
    skin_index = np.arange(n) if order is None else np.asarray(order)
    return Markers(x, np.tile([0.0, 0.0, 1.0], (n, 1)), np.full(n, h * h), np.zeros((n, 2)), skin_index,
                   np.full(n, -1), h)


def row_division(n, h=H):
    """What `spray_plan` reads of the division of `row_markers`' layer (plan 03b Task 6: the wave's domain graph is the
    division's boundaries at phi): a boundary of length h at every level between neighbours along the row, as markers'
    indices i < j."""
    return SimpleNamespace(i=np.arange(n - 1), j=np.arange(1, n), length=np.full((max(n - 1, 0), 2), h))


def stacked(mk, temps, mat=None, rng=None, settings=SETTINGS):
    """Skins on the markers: marker i's skin holds a cell of 10 um of solid metal's mass at each temperature of temps[i],
    top first; every cell moving with a random velocity when `rng` is given."""
    mat = aa7075() if mat is None else mat
    order = np.argsort(mk.skin_index)                    # the marker of each skin row
    sk = Skins(mat, mk.area[order], mk.kappa[order], settings)
    temps = [temps[m] for m in order]
    count = np.array([len(t) for t in temps])
    K = int(count.max())
    T = np.array([list(t) + [300.0] * (K - len(t)) for t in temps])
    valid = np.arange(K)[None, :] < count[:, None]
    m = np.where(valid, mat.rho * sk.area[:, None] * CELL, 0.0)
    U = None if rng is None else np.where(valid[..., None], rng.normal(0.0, 0.3, (sk.n, K, 3)), 0.0)
    sk.install(m, m * sk.props.enthalpy(T), count, U)
    return sk


def gas(n, closure=sf.CLOSURE_GIRIN, tau=60.0, G=0.0, u=2000.0, rho=7e-4, mu=1.6e-4, delta_a=7e-3, mach=1.0,
        decel=30.0):
    """The fields of a SurfaceFlowResult the spraying reads, one value for every marker, on the shock-layer branch."""
    full = lambda v: np.full(n, v, dtype=float)
    return SimpleNamespace(closure=np.full(n, closure), branch=sf.BRANCH_SHOCK_LAYER, u_eff=full(u), u_e=full(u),
                           rho_e=full(rho), mu_e=full(mu), delta_a=full(delta_a), mach_e=full(mach), tau=full(tau),
                           G=full(G), kn_local=full(1e-3), deceleration=decel)


STATE = SimpleNamespace(freestream=SimpleNamespace(rho=8e-5), V=7200.0, ma=24.0, h=71e3)


def model_of(mat, **kw):
    kw.setdefault("spray_gates", "2026-11-10")
    return spray.SprayModel(LiquidProperties(mat.liquid.rho, mat.liquid.mu, mat.liquid.sigma), **kw)


MIXED = [[LIQUID] * 30 + [RIGID] * 10,                       # deep liquid on rigid metal
         [LIQUID] * 4 + [RIGID] * 10,                        # a thin liquid film on rigid metal
         [LIQUID] * 6 + [SLURRY] * 20 + [RIGID] * 10,         # liquid on slurry
         [SLURRY] * 3 + [LIQUID] * 5 + [RIGID] * 10,          # slurry exposed over liquid
         [RIGID] * 3 + [LIQUID] * 5 + [RIGID] * 10,           # a floating crust
         [RIGID] * 12]                                       # rigid metal alone


def test_marker_release_equals_the_shared_model():
    """For six markers with given depths at a fixed gas state, the plan's mdot, r and branch are the shared
    SprayModel.evaluate's called directly with the same arrays, bit for bit -- the markers' skins in reversed rows."""
    n = len(MIXED)
    mk = row_markers(n, order=np.arange(n)[::-1])
    mat = aa7075()
    sk = stacked(mk, MIXED, mat)
    flow, model = gas(n), model_of(mat)
    theta = np.linspace(0.1, 0.6, n)
    windward, w_n = theta <= 0.5 * np.pi, flow.deceleration * np.cos(theta)
    plan = spraying.spray_plan(sk, mk, flow, STATE, model, windward, w_n, 0.5, 0.05, 0.05, theta=theta,
                               division=row_division(n))
    # the same arrays by hand, in the markers' order
    r = mk.skin_index
    F, L = sk.depth_above(mat.T_rigid)[r], sk.depth_above(mat.T_liquidus)[r]
    T = sk.T()
    flowing, liquid = sk.flowing(T), sk.liquid(T)
    n_flow, n_liq = np.cumprod(flowing, axis=1).sum(axis=1)[r], np.cumprod(liquid, axis=1).sum(axis=1)[r]
    delta_m = spray.spray_layer(flow, model.liquid, model.slip_band, theta, 0.05)
    v_s = lubrication(flow.tau, flow.G, F, delta_m, model.liquid.mu, b_layer=F)[0]
    A = sk.area[r]
    extent = spray.region_extent(np.column_stack([row_division(n).i, row_division(n).j]), L > 0.0, A)
    m_f = (sk.m * np.cumprod(liquid, axis=1))[r].sum(axis=1)
    want = model.evaluate(flow, STATE, F, delta_m, v_s, windward, 0.5, A, m_f, 0.05, b_layer=F, extent=extent,
                          deceleration_n=w_n, regime_layer=F, on_slurry=n_flow > n_liq, rt_layer=L, cap_depth=L)
    assert np.array_equal(plan.result.mdot, want.mdot) and np.array_equal(plan.rate, want.mdot)
    assert np.array_equal(plan.result.r, want.r, equal_nan=True)
    assert np.array_equal(plan.result.branch, want.branch)
    assert (want.mdot > 0.0).sum() >= 2 and len(set(want.branch.tolist())) >= 2    # the case is not trivial


def test_the_wave_s_region_is_the_division_s_boundaries_at_phi(monkeypatch):
    """Plan 03b Task 6: the wave's domain graph (`spray.region_extent`'s edges, for the liquid region and the slurry's
    flowing one) is the division's shared boundaries that have a length at phi, as markers' indices; a boundary the
    division finds only below phi (two columns meeting at depth beyond an edge) joins no region. Six markers in a row,
    the boundary between the third and fourth of no length at phi: the plan's edges are the other four."""
    mat = aa7075()
    n = len(MIXED)
    mk = row_markers(n)
    sk = stacked(mk, MIXED, mat)
    d = row_division(n)
    d.length[2, 0] = 0.0                                  # the third and fourth columns meet only below phi
    seen = []
    extent = spraying.spray_mod.region_extent

    def spy(edges, wet, area):
        seen.append(np.array(edges))
        return extent(edges, wet, area)
    monkeypatch.setattr(spraying.spray_mod, "region_extent", spy)
    theta = np.linspace(0.1, 0.6, n)
    spraying.spray_plan(sk, mk, gas(n), STATE, model_of(mat), np.ones(n, bool), np.full(n, 30.0), 0.5, 0.05, 0.05,
                        slurry=True, theta=theta, division=d)
    want = np.array([[0, 1], [1, 2], [3, 4], [4, 5]])
    assert len(seen) >= 1 and all(np.array_equal(e, want) for e in seen)
    with pytest.raises(ValueError, match="division"):
        spraying.spray_plan(sk, mk, gas(n), STATE, model_of(mat), np.ones(n, bool), np.full(n, 30.0), 0.5, 0.05, 0.05,
                            theta=theta)


def test_the_slip_band_needs_the_markers_angles():
    """With the model's slip band on, the band's depth reads the markers' angles: a plan without them is refused rather
    than silently given Girin's delta_m."""
    mat = aa7075()
    n = len(MIXED)
    mk = row_markers(n)
    sk = stacked(mk, MIXED, mat)
    model = model_of(mat, slip_band=True)
    with pytest.raises(ValueError, match="theta"):
        spraying.spray_plan(sk, mk, gas(n), STATE, model, np.ones(n, bool), np.full(n, 30.0), 0.5, 0.05, 0.05,
                            division=row_division(n))


SLURRY_COLUMNS = MIXED[:5] + [[LIQUID] * 2 + [907.8] * 300 + [RIGID] * 5, [907.5] * 200 + [RIGID] * 5,
                              [LIQUID] * 3 + [SLURRY] * 40 + [RIGID] * 5]


def test_the_slurry_mode_matches_the_skyscrapers_on_the_same_columns():
    """The exposed slurry's own mode, read from the level set's skins (`spraying.spray_plan` with the slurry) and from the
    skyscrapers' skins holding the same cells (`skinbody.slurry_modes`, as their `girin_plan` calls it), at one gas state
    on Girin's closure: the same rate, droplet radius, viscosity, shear rate, delta_m, layer depth and branch. Both call
    the shared `spray.slurry_modes`; what differs is each architecture's reading of its cells (each skin's own property
    table for the temperatures and depths), so the numbers agree to 1e-9; the column of deep warm slurry under liquid
    sprays on the thick branch, and the others reach a branch the viscous film or the critical Weber number holds."""
    from reentry_model import skin as sky, skinbody
    mat = aa7075()
    n = len(SLURRY_COLUMNS)
    mk = row_markers(n)
    ours = stacked(mk, SLURRY_COLUMNS, mat, rng=np.random.default_rng(2))
    theirs = sky.Skins(mat, ours.area.copy(), sky.SkinSettings(rheology=ours.settings.rheology))
    theirs._install(ours.m.copy(), ours.H.copy(), ours.count.copy(), U=ours.U.copy())
    flow, model = gas(n, tau=200.0, G=-3e4), model_of(mat)
    theta = np.linspace(0.05, 0.9, n)
    windward, w_n = theta <= 0.5 * np.pi, flow.deceleration * np.cos(theta)
    plan = spraying.spray_plan(ours, mk, flow, STATE, model, windward, w_n, 0.5, 0.05, 0.04, slurry=True,
                               theta=theta, division=row_division(n))
    F = theirs.depth_above(mat.T_rigid)
    flowing = theirs.valid & np.isfinite(sky.cell_viscosity(theirs, 0.0))
    liquid = flowing & (theirs.T >= mat.T_liquidus)
    n_flow, n_liq = np.cumprod(flowing, axis=1).sum(axis=1), np.cumprod(liquid, axis=1).sum(axis=1)
    edges = np.column_stack([row_division(n).i, row_division(n).j])
    rate, modes = skinbody.slurry_modes(theirs, flow, STATE, model, windward, w_n, theirs.area, edges, 0.05, 0.5, F, n_liq,
                                        n_flow, theta=theta, nose_radius=0.04)
    assert (rate > 0.0).sum() >= 1 and (modes["branch"] >= 0).sum() >= 3 and np.all(modes["branch"][rate > 0.0] == spray.BRANCH_THICK)
    np.testing.assert_allclose(plan.rate_slurry, rate, rtol=1e-9, atol=0.0)
    np.testing.assert_array_equal(plan.slurry["branch"], modes["branch"])
    for key in ("r", "mu", "gamma", "delta_m", "depth", "rho", "T"):
        np.testing.assert_allclose(plan.slurry[key], modes[key], rtol=1e-9, atol=0.0, equal_nan=True, err_msg=key)


def test_the_exposed_slurry_mode_matches_the_skyscrapers_on_one_column():
    """The controller's ruling (c) of Task 10: one column -- two liquid cells over 3 mm of warm slurry (907.8 K) over rigid
    metal, on one wide marker with no neighbours -- read by the level set (`spray_plan` with the slurry) and by the skyscrapers'
    `skinbody.slurry_modes` on their skins holding the same cells: the same slurry rate, radius, viscosity, shear rate,
    delta_m, layer depth, density, temperature and branch (to 1e-9, each skin finding its temperatures from its own
    property table), on the thick branch. Both call the shared `spray.slurry_modes`: this pins the wiring of one column,
    the region's extent of a lone marker included."""
    from reentry_model import skin as sky, skinbody
    mat = aa7075()
    mk = row_markers(1, h=0.05)                                      # a wide marker: the slurry's wave fits on it alone
    ours = stacked(mk, [[LIQUID] * 2 + [907.8] * 300 + [RIGID] * 5], mat)
    theirs = sky.Skins(mat, ours.area.copy(), sky.SkinSettings(rheology=ours.settings.rheology))
    theirs._install(ours.m.copy(), ours.H.copy(), ours.count.copy(), U=ours.U.copy())
    flow, model = gas(1, tau=200.0, G=-3e4), model_of(mat)
    theta = np.array([0.2])
    windward, w_n = np.array([True]), flow.deceleration * np.cos(theta)
    plan = spraying.spray_plan(ours, mk, flow, STATE, model, windward, w_n, 0.5, 0.05, 0.04, slurry=True, theta=theta,
                               division=row_division(1))
    F = theirs.depth_above(mat.T_rigid)
    flowing = theirs.valid & np.isfinite(sky.cell_viscosity(theirs, 0.0))
    liquid = flowing & (theirs.T >= mat.T_liquidus)
    n_flow, n_liq = np.cumprod(flowing, axis=1).sum(axis=1), np.cumprod(liquid, axis=1).sum(axis=1)
    rate, modes = skinbody.slurry_modes(theirs, flow, STATE, model, windward, w_n, theirs.area, np.zeros((0, 2), int),
                                        0.05, 0.5, F, n_liq, n_flow, theta=theta, nose_radius=0.04)
    assert rate[0] > 0.0 and modes["branch"][0] == spray.BRANCH_THICK
    np.testing.assert_allclose(plan.rate_slurry, rate, rtol=1e-9, atol=0.0)
    np.testing.assert_array_equal(plan.slurry["branch"], modes["branch"])
    for key in ("r", "mu", "gamma", "delta_m", "depth", "rho", "T"):
        np.testing.assert_allclose(plan.slurry[key], modes[key], rtol=1e-9, atol=0.0, err_msg=key)


TABLE_FIXTURE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "fixtures", "spraying_table_path.npz")


def table_path_case():
    """The table path's records and released masses (ruling P5 of plan 03b Task 10, the "before" of
    `test_the_table_path_is_unchanged`): on SLURRY_COLUMNS' skins, every cell moving, the markers in reversed rows, two
    plans for a 0.5 s macro step on Girin's theory (a spray model whose `dispersion` is the default, "table") -- the
    liquid's mode alone at the default gas, and with the slurry at the slurry tests' gas -- each released over four flow
    sub-steps of 0.125 s and recorded at 0.5 s. Returns {name: array}: per case the records (`SprayPlan.records`), the
    mass each sub-step released per skin row, the plan's rates, books, droplets and size histograms, and the step
    summary's names and values (`step_summary`)."""
    mat = aa7075()
    n = len(SLURRY_COLUMNS)
    theta = np.linspace(0.05, 0.9, n)
    out = {}
    for case, slurry, flow in (("liquid", False, gas(n)), ("slurry", True, gas(n, tau=200.0, G=-3e4))):
        mk = row_markers(n, order=np.arange(n)[::-1])
        sk = stacked(mk, SLURRY_COLUMNS, mat, rng=np.random.default_rng(17))
        plan = spraying.spray_plan(sk, mk, flow, STATE, model_of(mat), theta <= 0.5 * np.pi,
                                   flow.deceleration * np.cos(theta), 0.5, 0.05, 0.04, slurry=slurry, theta=theta,
                                   centre=np.array([0.0, 0.0, -0.05]), v_hat=np.array([0.0, 0.0, -1.0]),
                                   division=row_division(n))
        out[case + "_released"] = np.array([plan.release(sk, 0.125).mass for _ in range(4)])
        out[case + "_records"] = plan.records(0.5)
        for name in ("rate", "rate_slurry", "mass", "mass_slurry", "enthalpy", "enthalpy_slurry", "momentum",
                     "momentum_slurry", "number", "number_slurry", "hist_n", "hist_m"):
            out[case + "_" + name] = np.asarray(getattr(plan, name), dtype=float)
        summary = plan.step_summary(sk, n_before=1.0)
        out[case + "_summary_names"] = np.array(sorted(summary))
        out[case + "_summary"] = np.array([summary[k] for k in sorted(summary)], dtype=float)
    return out


def write_table_fixture(path=TABLE_FIXTURE):
    """Write `table_path_case` to `path` (not a test). The committed file was written on the package as it stood at plan
    03b Task 10's BASE, 964b415, before spraying.py was edited (the commit that added the file changed no code), from
    "Reentry Model/": `PYTHONPATH=tests:tests/levelset python -c "import test_levelset_spraying as t;
    t.write_table_fixture()"`."""
    np.savez(path, **table_path_case())


def test_the_table_path_is_unchanged():
    """Plan 03b Task 10 (ruling P5): on Girin's theory, the spray model's default `dispersion` "table", the plans'
    records, the masses each flow sub-step released, the plans' books, droplets and histograms and the step summary are
    byte for byte those recorded at the task's BASE, before the profile path was added (`table_path_case`,
    TABLE_FIXTURE): same names, dtypes, shapes and bytes. Both cases release from both modes, so the records hold rows
    of each."""
    want = np.load(TABLE_FIXTURE)
    got = table_path_case()
    assert sorted(got) == sorted(want.files)
    for key in want.files:
        a, b = got[key], want[key]
        assert a.dtype == b.dtype and a.shape == b.shape and a.tobytes() == b.tobytes(), key
    assert got["liquid_mass"].sum() > 0.0 and got["slurry_mass"].sum() > 0.0 and got["slurry_mass_slurry"].sum() > 0.0
    assert got["liquid_records"].shape[1] == len(spray.SOURCE_COLUMNS)


def test_the_cascade_takes_liquid_first_and_never_rigid_metal():
    """Asked for everything over a long time, the release takes the run of liquid cells from the top; with the slurry
    every flowing cell down to the first rigid one; never rigid metal, which is as it was in every column."""
    mat = aa7075()
    for slurry, cells in ((False, [30, 4, 6, 0, 0, 0]), (True, [30, 4, 26, 8, 0, 0])):
        mk = row_markers(len(MIXED))
        sk = stacked(mk, MIXED, mat, rng=np.random.default_rng(3))
        m0 = sk.m[:, 0].copy()
        T = sk.T()
        rigid0 = np.where(sk.valid & ~sk.flowing(T), sk.m, 0.0).sum(axis=1)
        big = np.full(sk.n, 1e6)
        dm, _ = spraying.cascade_mass(sk, big, big if slurry else np.zeros(sk.n), sk.area, 1.0, slurry)
        book = spraying.release_top(sk, dm, slurry)
        np.testing.assert_allclose(book.mass, np.array(cells) * m0, rtol=1e-13, atol=0.0)
        T = sk.T()
        rigid = np.where(sk.valid & ~sk.flowing(T), sk.m, 0.0).sum(axis=1)
        np.testing.assert_allclose(rigid, rigid0, rtol=1e-15, atol=0.0)
        assert np.all(book.by_phase.mass[1] == 0.0) if not slurry else book.by_phase.mass[1][2] > 0.0
    # in time order: the liquid's rate takes two liquid cells in two thirds of the step, the slurry's half a cell after
    mk = row_markers(1)
    sk = stacked(mk, [[LIQUID] * 2 + [SLURRY] * 3 + [RIGID] * 3], mat)
    m = sk.m[0, 0]
    dt = 0.01
    dm, from_slurry = spraying.cascade_mass(sk, np.array([3 * m / (sk.area[0] * dt)]),
                                            np.array([1.5 * m / (sk.area[0] * dt)]), sk.area, dt, True)
    np.testing.assert_allclose([dm[0], from_slurry[0]], [2.5 * m, 0.5 * m], rtol=1e-12)


def test_the_cascade_matches_the_skyscrapers_release_on_one_column():
    """Decision D4: the same columns of cells -- liquid over slurry over rigid metal, slurry over liquid, a crust, a warm
    slurry cell -- released by both architectures' cascades over one macro step of eight flow sub-steps, at a liquid
    rate and the slurry's own with a slurry floor: the same masses, enthalpies, momenta and solid, in total and by the cells
    they came from (liquid and the others), and the same mass taken at the slurry's rate, within 1e-12. The solid is
    within SOLID_TOL: it is the cells' liquid fractions at their temperatures, and each architecture's skins find a
    temperature from the enthalpy their own way (the level set on its 0.01 K property table, the skyscrapers through the
    material), which moves the solid of a partly solid cell by 6e-10 of itself (measured 2026-11-12)."""
    from reentry_model import skin as sky
    mat = aa7075()
    temps = [[LIQUID] * 2 + [SLURRY] * 3 + [RIGID] * 3, [SLURRY] * 2 + [LIQUID] * 2 + [RIGID] * 3,
             [LIQUID] * 2 + [RIGID] * 3, [RIGID] * 4, [LIQUID] * 2 + [907.0] + [SLURRY] * 2 + [RIGID] * 2,
             [RIGID] * 2 + [LIQUID] * 3 + [RIGID] * 2]
    mk = row_markers(len(temps))
    for slurry in (False, True):
        ours = stacked(mk, temps, mat, rng=np.random.default_rng(11))
        theirs = sky.Skins(mat, ours.area.copy(), sky.SkinSettings(rheology=ours.settings.rheology))
        theirs._install(ours.m.copy(), ours.H.copy(), ours.count.copy(), U=ours.U.copy())
        m = ours.m[:, 0]
        dt, steps = 0.5, 8
        h = dt / steps
        rate, rate_s = 0.7 * m / (ours.area * h), 0.4 * m / (ours.area * h)
        floor = np.full(ours.n, 905.0)
        fields = ("mass", "enthalpy", "momentum", "solid")
        pieces = lambda bk, ds: [getattr(bk, f) for f in fields] + [getattr(bk.by_phase, f) for f in fields] + [ds]
        got = ref = None
        for _ in range(steps):
            dm, ds = spraying.cascade_mass(ours, rate, rate_s, ours.area, h, slurry, floor)
            b = pieces(spraying.release_top(ours, dm, slurry), ds)
            dm2, ds2 = sky.cascade_mass(theirs, rate, rate_s, theirs.area, h, slurry, floor)
            b2 = pieces(sky.release_top(theirs, dm2, slurry), ds2)
            got = b if got is None else [a + x for a, x in zip(got, b)]
            ref = b2 if ref is None else [a + x for a, x in zip(ref, b2)]
        assert got[0].sum() > 0.0 and (got[-1].sum() > 0.0) == slurry
        assert (got[3].sum() > 0.0) == slurry and (got[4][1].sum() > 0.0) == slurry   # the solid, the slurry's share
        for name, a, b in zip(fields + tuple("phase " + f for f in fields) + ("released_slurry",), got, ref):
            scale = np.abs(b).max()
            tol = SOLID_TOL if name.endswith("solid") else 1e-12
            np.testing.assert_allclose(a, b, rtol=tol, atol=tol * scale if scale > 0.0 else 0.0, err_msg=name)
        np.testing.assert_allclose(ours.m.sum(axis=1), theirs.m.sum(axis=1), rtol=1e-12)


def test_release_is_not_capped_by_a_renewal_time():
    """The ruling of 2026-11-11: no renewal limit. A 2 mm liquid column on rigid metal, its plan evaluated for a 0.5 s
    macro step (Girin's thick branch, delta_m about 0.19 mm), is released over four flow sub-steps of 2.5 ms: the
    shared rate asks far more than design 9.5's clause would allow over that time (rho_l delta_m A per renewal time
    delta_m^2/nu_l, about 6.9 kg/(m2 s), against about 300), and the release takes the shared rate's mass, mdot A t,
    within 1e-12 -- about half the column, so neither the column nor any cap binds."""
    mat = aa7075()
    mk = row_markers(1, h=0.05)                                      # a wide marker: the wave fits
    sk = stacked(mk, [[LIQUID] * 200 + [RIGID] * 20], mat)
    flow, model = gas(1, tau=TAU_RENEWAL), model_of(mat)
    theta = np.array([0.3])
    plan = spraying.spray_plan(sk, mk, flow, STATE, model, np.array([True]), flow.deceleration * np.cos(theta), 0.5,
                               0.05, 0.05, theta=theta, division=row_division(1))
    assert plan.result.branch[0] == spray.BRANCH_THICK
    liq = model.liquid
    d = plan.result.delta_m[0]
    A = sk.area[0]
    h, steps = 2.5e-3, 4
    renewal = liq.rho * d * A * (steps * h) / (d * d / (liq.mu / liq.rho))   # kg the clause would allow
    asked = plan.rate[0] * A * steps * h
    M0 = sk.totals()[0][0]
    column = sk.m[0, :200].sum()
    assert asked > 10.0 * renewal and 0.3 * column < asked < 0.9 * column
    total = sum(plan.release(sk, h).mass.sum() for _ in range(steps))
    assert total == pytest.approx(asked, rel=1e-12)
    assert plan.mass[0] == pytest.approx(asked, rel=1e-12) and M0 - sk.totals()[0][0] == pytest.approx(asked, rel=1e-12)


def test_release_books_close():
    """On mixed skins moving in every cell, a release with and without the slurry: each skin's mass, enthalpy and momentum
    fall by exactly what the book says left (1e-14 of the largest), and the plan's record is the books' sum."""
    mat = aa7075()
    n = len(MIXED)
    for slurry in (False, True):
        mk = row_markers(n, order=np.arange(n)[::-1])
        sk = stacked(mk, MIXED, mat, rng=np.random.default_rng(5))
        flow, model = gas(n), model_of(mat)
        theta = np.linspace(0.1, 0.6, n)
        plan = spraying.spray_plan(sk, mk, flow, STATE, model, theta <= 0.5 * np.pi, 30.0 * np.cos(theta), 0.5, 0.05,
                                   0.05, slurry=slurry, theta=theta, division=row_division(n))
        plan.rate[:] = np.maximum(plan.rate, 50.0)               # every skin releases something it may
        if slurry:
            plan.rate_slurry[:] = 20.0
            plan.slurry["T"][:] = SLURRY - 1.0
        M0, E0, P0 = sk.totals()
        books = [plan.release(sk, 0.01) for _ in range(5)]
        M, E, P = sk.totals()
        mass = sum(b.mass for b in books)
        assert mass.sum() > 0.0
        np.testing.assert_allclose(M, M0 - mass, rtol=0.0, atol=1e-14 * M0.max())
        np.testing.assert_allclose(E, E0 - sum(b.enthalpy for b in books), rtol=0.0, atol=1e-14 * np.abs(E0).max())
        np.testing.assert_allclose(P, P0 - sum(b.momentum for b in books), rtol=0.0, atol=1e-14 * np.abs(P0).max())
        r = mk.skin_index
        np.testing.assert_allclose(plan.mass + plan.mass_slurry, mass[r], rtol=1e-14, atol=0.0)
        assert (plan.mass_slurry.sum() > 0.0) == slurry


class History:
    """A run history as the shared writer reads it: `columns` by name and its length in rows."""

    def __init__(self, columns):
        self.columns = columns

    def __len__(self):
        return len(self.columns["time_s"])


def test_droplet_records_go_through_the_shared_writer(tmp_path):
    """A fake history of three macro steps, each a plan released over its step: the shared `write_particles` writes
    particles.npz, particles_summary.csv and size_distribution.csv from the plans' rows and summaries, and the released
    mass in each equals the books'."""
    mat = aa7075()
    n = len(MIXED)
    mk = row_markers(n)
    sk = stacked(mk, MIXED, mat, rng=np.random.default_rng(9))
    model = model_of(mat)
    theta = np.linspace(0.1, 0.6, n)
    rows, history, booked, number = [], {k: [] for _, k in droplets.SKIN_SUMMARY_COLUMNS}, 0.0, 0.0
    for step in range(3):
        plan = spraying.spray_plan(sk, mk, gas(n), STATE, model, theta <= 0.5 * np.pi, 30.0 * np.cos(theta), 0.5,
                                   0.05, 0.05, slurry=True, theta=theta, centre=np.array([0.0, 0.0, -0.05]),
                                   v_hat=np.array([0.0, 0.0, -1.0]), division=row_division(n))
        for _ in range(4):
            booked += plan.release(sk, 0.125).mass.sum()
        block = plan.records(0.5 * (step + 1))
        rows.append(block)
        summary = plan.step_summary(sk, n_before=number)
        number = summary["n_released"]
        assert set(history) <= set(summary)
        for k in history:
            history[k].append(summary[k])
    assert booked > 0.0
    hist = History({k: np.array(v) for k, v in history.items()})
    paths = droplets.write_particles(str(tmp_path), SimpleNamespace(source_rows=rows), hist, window=1.0)
    assert set(paths) == {"particles", "particles_summary", "size_distribution"}
    cols = np.load(paths["particles"])
    assert cols["dm_kg"].sum() == pytest.approx(booked, rel=1e-12)
    with open(paths["particles_summary"]) as fh:
        table = list(csv.DictReader(fh))
    assert len(table) == 3 and sum(float(r["released_mass_kg"]) for r in table) == pytest.approx(booked, rel=1e-8)
    with open(paths["size_distribution"]) as fh:
        dist = [r for r in csv.DictReader(fh)]
    total = [r for r in dist if float(r["window_start_s"]) == 0.0 and float(r["window_end_s"]) == 1.5]
    assert sum(float(r["dm_kg"]) for r in total) == pytest.approx(booked, rel=1e-8)


class Iron:
    """Girin and Kopyt's (1994) liquid: rho_1 7.8 g/cm3, sigma_t 1200 dyn/cm (their Tables 1-2); mu an iron melt's."""
    rho, mu, sigma = 7800.0, 5.8e-3, 1.2


@pytest.mark.parametrize("dispersion", spray.DISPERSION_NAMES)
def test_girin_cases_through_the_level_set_path(dispersion):
    """Girin and Kopyt (1994) Table 2, the front-surface Rayleigh-Taylor mode, fed through `spray_plan` on one flat
    marker at the stagnation point holding a pool of their liquid one wavelength deep (so tanh(k h) is 1 to 1e-5), off
    Girin's closure (no shear mode): the plan's wavelength is their lambda* (Eq. 14) within 2 % and its growth time
    their tau* within 1 %, at each of the three decelerations. On either dispersion (plan 03b Task 10: Girin's published
    cases are unchanged on the profile path, whose pool at rest resolves no profile and keeps the table path, the
    release and growth time the table path's to the bit)."""
    with open(os.path.join(PROJECT_ROOT, "data", "reference_values", "girin1994_tables.json")) as fh:
        t = json.load(fh)["table2"]
    mat = LinearMeltMaterial(rho=Iron.rho, rho_l=Iron.rho, mu=Iron.mu, sigma=Iron.sigma)
    model = spray.SprayModel(LiquidProperties(Iron.rho, Iron.mu, Iron.sigma), spray_gates="2026-11-10",
                             dispersion=dispersion)
    table = spray.SprayModel(LiquidProperties(Iron.rho, Iron.mu, Iron.sigma), spray_gates="2026-11-10")
    for W, lam_ref, tau_ref in zip(t["W_ms2"], t["lambda_star_m_eq14"], t["tau_star_s"]):
        plans = []
        for m in (model, table):
            mk = row_markers(1, h=1.0)                                 # 1 m2: the pool admits the fastest wave
            sk = Skins(mat, mk.area, mk.kappa, SkinSettings(cell=20e-6, fine_zone=0.4e-3))
            sk.fill(Iron.rho * lam_ref, 1000.0)                        # one wavelength deep, liquid
            flow = gas(1, closure=sf.CLOSURE_COUETTE, tau=0.0, decel=W)
            theta = np.zeros(1)
            plans.append(spraying.spray_plan(sk, mk, flow, STATE, m, np.array([True]), np.array([W]), 0.5, 0.05, 0.05,
                                             theta=theta, division=row_division(1)))
        plan, base = plans
        res = plan.result
        assert res.branch[0] == spray.BRANCH_RT and plan.rate[0] > 0.0
        assert res.wavelength[0] == pytest.approx(lam_ref, rel=2e-2)
        assert res.growth[0] == pytest.approx(tau_ref, rel=1e-2)
        assert np.array_equal(plan.rate, base.rate) and np.array_equal(res.growth, base.result.growth)
        assert (plan.resolved is None) if dispersion == "table" else not plan.resolved.any()


# -- the profile path (plan 03b Task 10) ----------------------------------------------------------------------------------
def profile_column(n_fast, n_slow, n_slurry, n_rigid, u_top, u_kink, u_slurry, T_slurry=SLURRY):
    """One skin's temperatures and speeds [m/s] along x, top first: n_fast liquid cells sheared from u_top to u_kink,
    n_slow liquid cells from u_kink to u_slurry, n_slurry slurry cells at T_slurry creeping from u_slurry to rest at the
    rigid front, and n_rigid rigid cells at rest -- a profile with a kink where the liquid's shear drops."""
    temps = [LIQUID] * (n_fast + n_slow) + [T_slurry] * n_slurry + [RIGID] * n_rigid
    u = np.concatenate([np.linspace(u_top, u_kink, n_fast, endpoint=False),
                        np.linspace(u_kink, u_slurry, n_slow, endpoint=False),
                        np.linspace(u_slurry, 0.0, n_slurry + 1)[:-1], np.zeros(n_rigid)])
    return temps, u


def moving(mk, temps, u, mat, settings=SETTINGS):
    """The level set's skins on `mk` holding `temps` (`stacked`) with the speeds `u` [m/s] along x (one list per marker,
    top first)."""
    ours = stacked(mk, temps, mat, settings=settings)
    U = np.zeros(ours.m.shape + (3,))
    for i, row in enumerate(mk.skin_index):
        U[row, :len(u[i]), 0] = u[i]
    ours.U = np.where(ours.valid[..., None], U, 0.0)
    return ours


def both_skins(mk, temps, u, mat, settings=SETTINGS):
    """`moving`'s skins, and the skyscrapers' skins holding the same cells."""
    from reentry_model import skin as sky
    ours = moving(mk, temps, u, mat, settings)
    theirs = sky.Skins(mat, ours.area.copy(), sky.SkinSettings(rheology=ours.settings.rheology))
    theirs._install(ours.m.copy(), ours.H.copy(), ours.count.copy(), U=ours.U.copy())
    return ours, theirs


def test_layer_columns_match_the_skyscrapers():
    """`spraying.layer_columns` reads the skins' cells as the skyscrapers' `skinbody.layer_columns` reads theirs, on the
    same cells (MIXED, every cell moving at about a millimetre a second, a floating crust among them) under both
    rheology laws: the same depths, thicknesses, speeds, temperatures, masses, densities and fronts, the same flowing and
    liquid runs, and the same viscosities at the cells' shear rates. Under the bridge law the slurry's viscosity depends
    on its rate, here some 30 to 150 1/s, below the law's cap (rheology.GAMMA_MAX, 367 1/s), so this pins the shear rates
    the momentum step starts from (`_shear_rates`, the skyscrapers' `skin.cell_shear_rates`): those viscosities move
    by more than 10 % from their values at rest. Each skin finds its temperatures from its own property table, so the
    floats agree to 1e-9 (the slurry tests' tolerance)."""
    from dataclasses import replace

    from reentry_model import skinbody as sky_body
    mat = aa7075()
    n = len(MIXED)
    mk = row_markers(n)
    rng = np.random.default_rng(23)
    u = [rng.normal(0.0, 1e-3, len(t)) for t in MIXED]
    tau = np.linspace(5.0, 200.0, n)
    for law in ("composite", "bridge"):
        ours, theirs = both_skins(mk, MIXED, u, mat, replace(SETTINGS, rheology=law))
        got, want = spraying.layer_columns(ours, tau), sky_body.layer_columns(theirs, tau)
        for name in ("flowing", "liquid"):
            assert np.array_equal(getattr(got, name), getattr(want, name)), (law, name)
        for name in ("z", "dz", "u", "T", "m", "rho", "mu", "front"):
            np.testing.assert_allclose(getattr(got, name), getattr(want, name), rtol=1e-9, atol=0.0, equal_nan=True,
                                       err_msg=law + " " + name)
        assert got.flowing.sum() >= 30 and np.isfinite(got.mu[got.flowing]).all()
    slurry = got.flowing & ~got.liquid                                         # the bridge law's, from the last pass
    rest = ours.viscosity(0.0)[slurry]
    assert slurry.sum() >= 20 and np.abs(got.mu[slurry] / rest - 1.0).max() > 0.1


def assert_plans_agree(plan, ref):
    """A level-set plan and the skyscrapers' on the same cells (both on the profile path): the same branches, resolved
    columns and slurry branches; the thresholds T_run within 1e-12 (NaN where no profile mode sprays); the rates, the
    droplets' density and viscosity, the slurry's own modes, the records' appended columns and the profile modes within
    1e-9 (each skin finds its temperatures from its own property table)."""
    assert np.array_equal(plan.result.branch, ref.result.branch) and np.array_equal(plan.resolved, ref.resolved)
    assert np.array_equal(plan.slurry["branch"], ref.slurry["branch"])
    np.testing.assert_allclose(plan.T_run, ref.T_run, rtol=1e-12, atol=0.0, equal_nan=True)
    pairs = [("rate", plan.rate, ref.rate), ("rate_slurry", plan.rate_slurry, ref.rate_slurry),
             ("rho", plan.props.rho, ref.props.rho), ("mu", plan.props.mu, ref.props.mu),
             ("profile", plan.profile, ref.profile)]
    pairs += [("slurry " + k, plan.slurry[k], ref.slurry[k]) for k in ("r", "mu", "gamma", "delta_m", "depth", "rho", "T")]
    pairs += [("modes " + k, getattr(plan.modes, k), getattr(ref.modes, k)) for k in
              ("growth", "k", "wall_depth", "nu_eff", "we_s", "cap_gradient", "cap_driven")]
    for name, a, b in pairs:
        np.testing.assert_allclose(a, b, rtol=1e-9, atol=0.0, equal_nan=True, err_msg=name)


# one column each: (the column's profile, the slurry's own mode on, the flow sub-step [s]); the first is taken by the
# gas-driven mode on the rigid front, its run at the slurry's temperature (the slurry beneath the liquid released at the
# mode's rate); the second by the gradient mode, its run the liquid, with the exposed slurry's own mode beneath it
PROFILE_CASES = ((profile_column(4, 16, 30, 10, 4.0, 1.0, 0.05), False, 4e-4, spray.BRANCH_THIN, SLURRY),
                 (profile_column(4, 16, 60, 5, 8.0, 1.0, 0.05, 907.8), True, 1e-4, spray.BRANCH_PROFILE, LIQUID))
PROFILE_STEPS = 6     # flow sub-steps released: the cascade reaches into the slurry and stops short of the column's end
PROFILE_SOLID_TOL = 1e-8   # the solid released, relative: each skin's own property table gives the liquid fraction of the
                           # second case's 907.8 K slurry, 1.3 % solid, 3.4e-9 of its solid apart (measured 2026-11-15;
                           # 8e-11 at 902 K), past the D4 test's SOLID_TOL


def test_the_profile_release_matches_the_skyscrapers_on_one_column():
    """Plan 03b Task 10, decision D4 on the profile path: one wide marker's skin and the skyscrapers' column holding the
    same cells, a kinked profile (PROFILE_CASES), planned by `spraying.spray_plan` and by their `skinbody.girin_plan` on
    a spray model with `dispersion` "profile", then released over PROFILE_STEPS flow sub-steps -- the level set by its
    plan's `release`, the skyscrapers by their `skin.cascade_mass` and `skin.release_top` with their own plan's
    thresholds and slurry floor, as their flow does with its `Cascade`, at the level set's rates. The plans agree: the
    same branch, threshold T_run and resolved column, the rates and the droplets' density and viscosity within 1e-9. The
    releases agree within 1e-12: mass, enthalpy and momentum in total, of each mode, and taken at the slurry's rate (the
    solid within PROFILE_SOLID_TOL, by each skin's own liquid fractions, as in D4). The records agree: 36 columns, the same rows and branches, and every value within
    1e-9 (each skin finds its temperatures from its own property table: about 2e-10 at most, measured) but the solid
    fraction, within PROFILE_SOLID_TOL. In the first case
    the gas-driven mode's run reaches into the slurry (heat-limited at its temperature); in the second the gradient mode
    takes the liquid and the exposed slurry's own mode the slurry beneath."""
    from reentry_model import skin as sky, skinbody as sky_body
    mat = aa7075()
    model = model_of(mat, dispersion="profile")
    centre, v_hat = np.array([0.0, 0.0, -0.05]), np.array([0.0, 0.0, -1.0])
    for (temps, u), slurry, h, branch, T_run in PROFILE_CASES:
        mk = row_markers(1, h=0.05)                                  # a wide marker: the waves fit on it alone
        ours, theirs = both_skins(mk, [temps], [u], mat)
        flow = gas(1, tau=200.0, G=-3e4)
        theta = np.array([0.2])
        windward, w_n = np.array([True]), flow.deceleration * np.cos(theta)
        plan = spraying.spray_plan(ours, mk, flow, STATE, model, windward, w_n, 0.5, 0.05, 0.04, slurry=slurry,
                                   theta=theta, centre=centre, v_hat=v_hat, division=row_division(1))
        ref = sky_body.girin_plan(theirs, flow, STATE, model, theta, theirs.area, np.zeros((0, 2), int), 0.05, 0.5,
                                  slurry=slurry, points=mk.x.copy(), centre=centre, v_hat=v_hat, nose_radius=0.04)
        assert plan.result.branch[0] == ref.result.branch[0] == branch and plan.resolved[0] and ref.resolved[0]
        assert plan.T_run[0] == pytest.approx(T_run, abs=1e-9) and ref.T_run[0] == pytest.approx(T_run, abs=1e-9)
        assert (plan.rate_slurry[0] > 0.0) == slurry
        assert_plans_agree(plan, ref)
        m0 = ours.m.sum()
        floor = np.where(np.isfinite(ref.slurry["T"]), ref.slurry["T"], np.inf)
        got = {k: 0.0 for k in ("mass", "enthalpy", "momentum", "solid", "slurry")}
        want = dict(got)
        for _ in range(PROFILE_STEPS):
            book = plan.release(ours, h)
            dm, from_slurry = sky.cascade_mass(theirs, plan.rate, plan.rate_slurry, theirs.area, h, slurry, floor,
                                               ref.T_run)
            rb = sky.release_top(theirs, dm, slurry, ref.T_run)
            ph = rb.by_phase
            for key, a in (("mass", ph.mass), ("enthalpy", ph.enthalpy), ("momentum", ph.momentum),
                           ("solid", ph.solid)):
                setattr(ref, key, getattr(ref, key) + a[0])
                setattr(ref, key + "_slurry", getattr(ref, key + "_slurry") + a[1])
                got[key] = got[key] + getattr(book, key)
                want[key] = want[key] + getattr(rb, key)
            got["slurry"] = got["slurry"] + book.released_slurry
            want["slurry"] = want["slurry"] + from_slurry
        released = plan.mass.sum() + plan.mass_slurry.sum()
        assert 0.0 < released < 0.9 * m0 and (plan.mass_slurry.sum() > 0.0) == slurry
        assert plan.solid.sum() > 0.0 if not slurry else plan.solid_slurry.sum() > 0.0   # slurry left in either mode
        for key in ("mass", "enthalpy", "momentum", "solid", "slurry"):
            tol = PROFILE_SOLID_TOL if key == "solid" else 1e-12
            np.testing.assert_allclose(got[key], want[key], rtol=tol, atol=0.0, err_msg=key)
        for key in ("mass", "enthalpy", "momentum", "solid"):
            for name in (key, key + "_slurry"):
                tol = PROFILE_SOLID_TOL if key == "solid" else 1e-12
                np.testing.assert_allclose(getattr(plan, name), getattr(ref, name), rtol=tol, atol=0.0, err_msg=name)
        block = plan.records(0.5)
        ref.count(ref.props.rho)
        block_ref = ref.source_block(0.5)
        names = spray.record_columns(36)
        assert block.shape == block_ref.shape == (2 if slurry else 1, 36)
        assert np.array_equal(np.isnan(block), np.isnan(block_ref))
        for name in ("closure", "branch", "breakup", "floor_like_wall"):
            j = names.index(name)
            assert np.array_equal(block[:, j], block_ref[:, j]), name
        solid = names.index("solid_fraction")                       # by each skin's own liquid fractions
        np.testing.assert_allclose(block[:, solid], block_ref[:, solid], rtol=PROFILE_SOLID_TOL, atol=0.0)
        rest = [j for j in range(len(names)) if j != solid]
        np.testing.assert_allclose(block[:, rest], block_ref[:, rest], rtol=1e-9, atol=0.0, equal_nan=True)


def test_the_profile_plan_matches_the_skyscrapers_on_a_row():
    """The profile path's wiring across markers, against the skyscrapers' `girin_plan` on columns holding the same cells
    on the same graph (`assert_plans_agree`): a row of 1 mm markers holding PROFILE_CASES' two kinked skins, a moving
    slurry with no liquid, and a moving slurry over liquid, with the exposed slurry's mode on. The slurry over liquid
    has no liquid at its top, so its wave fits only in the region of the layers above the walls (the liquid's region
    would give it none), and it sprays; the gas-driven skin's run takes all its flowing metal, so no slurry is left for
    the slurry's own mode beneath it. The same markers with their skins in reversed rows (the final review's M12, cited
    finding F36; the table path's tests do so) give the same plan by marker (`assert_plans_agree`), and a release of a
    flow sub-step takes each marker's metal from its own skin row, the same masses, enthalpies and momenta to 1e-12."""
    from reentry_model import skinbody as sky_body
    mat = aa7075()
    slurry_top = ([SLURRY] * 60 + [RIGID] * 10, np.concatenate([np.linspace(0.8, 0.0, 60, endpoint=False),
                                                                 np.zeros(10)]))
    over_liquid = ([SLURRY] * 3 + [LIQUID] * 30 + [SLURRY] * 30 + [RIGID] * 10,
                   np.concatenate([np.linspace(2.0, 1.0, 3, endpoint=False), np.linspace(1.0, 0.5, 30, endpoint=False),
                                   np.linspace(0.5, 0.0, 30, endpoint=False), np.zeros(10)]))
    columns = [PROFILE_CASES[1][0], slurry_top, PROFILE_CASES[0][0], over_liquid]
    n, h = len(columns), 1e-3
    mk = row_markers(n, h=h)
    ours, theirs = both_skins(mk, [c[0] for c in columns], [c[1] for c in columns], mat)
    flow = gas(n, tau=200.0, G=-3e4)
    theta = np.linspace(0.1, 0.5, n)
    windward, w_n = theta <= 0.5 * np.pi, flow.deceleration * np.cos(theta)
    model = model_of(mat, dispersion="profile")
    d = row_division(n, h=h)
    plan = spraying.spray_plan(ours, mk, flow, STATE, model, windward, w_n, 0.5, 0.05, 0.04, slurry=True, theta=theta,
                               division=d)
    ref = sky_body.girin_plan(theirs, flow, STATE, model, theta, theirs.area, np.column_stack([d.i, d.j]), 0.05, 0.5,
                              slurry=True, points=mk.x.copy(), nose_radius=0.04)
    assert plan.resolved.all() and np.array_equal(plan.liquid > 0.0, [True, False, True, False])
    assert plan.rate[3] > 0.0 and plan.rate_slurry[0] > 0.0 and plan.slurry["branch"][2] == -1
    assert_plans_agree(plan, ref)
    mk_r = row_markers(n, h=h, order=np.arange(n)[::-1])                       # the skins in reversed rows
    ours_r, _ = both_skins(mk_r, [c[0] for c in columns], [c[1] for c in columns], mat)
    plan_r = spraying.spray_plan(ours_r, mk_r, flow, STATE, model, windward, w_n, 0.5, 0.05, 0.04, slurry=True,
                                 theta=theta, division=d)
    assert_plans_agree(plan_r, plan)
    r = mk_r.skin_index
    got, want = plan_r.release(ours_r, 1e-4), plan.release(ours, 1e-4)
    assert want.mass.sum() > 0.0
    for key in ("mass", "enthalpy", "momentum"):
        a, b = np.asarray(getattr(got, key))[r], np.asarray(getattr(want, key))
        np.testing.assert_allclose(a, b, rtol=1e-12, atol=1e-12 * np.abs(b).max(), err_msg=key)


def test_a_skin_flowing_to_its_base_takes_its_base_as_the_wall():
    """Design 21.7: where the flowing layer reaches the skin's base the skin's base acts as the wall. A skin that flows
    to its base -- liquid throughout, its speed falling linearly from the top to rest at the base, so its profile has no
    kink -- has its front at the column's own depth (`layer_columns`: `depth_above(T_rigid)`, the whole thickness), and
    the profile path resolves it with the wall at that depth, the rigid front's kind (not a floor-like kink) and no kink."""
    mat = aa7075()
    mk = row_markers(1, h=0.05)
    n_cells = 60
    sk = stacked(mk, [[LIQUID] * n_cells], mat)
    depth = sk.thickness()[0]
    z = (np.arange(n_cells) + 0.5) * depth / n_cells                  # equal liquid cells: the centres
    sk.U = np.zeros(sk.m.shape + (3,))
    sk.U[0, :, 0] = 3.0 * (1.0 - z / depth)
    cols = spraying.layer_columns(sk, np.array([60.0]))
    assert cols.front[0] == depth and cols.flowing[0].all()
    flow = gas(1)
    plan = spraying.spray_plan(sk, mk, flow, STATE, model_of(mat, dispersion="profile"), np.array([True]),
                               flow.deceleration * np.cos([0.2]), 0.5, 0.05, 0.04, theta=np.array([0.2]),
                               division=row_division(1))
    assert plan.resolved[0]
    assert plan.modes.wall_depth[0] == pytest.approx(depth, rel=1e-12)
    assert plan.record.floor_like_wall[0] == 0.0 and np.isnan(plan.record.kink_depth[0])


def test_profile_records_are_36_wide(tmp_path):
    """On the profile path the records carry spray.PROFILE_COLUMNS after the 27 SOURCE_COLUMNS: each liquid's row its
    marker's `profile` (melt_profile.profile_rows), each slurry's row the same with no release depth and no growth rate
    (ruling T10-a, the skyscrapers' rule); the body names its records `spray.record_columns(36)` on the profile path and
    27 on the table path or without a spray model (`LevelSetSkinBody.record_names`), and the shared writer writes the 36.
    Three wide markers: PROFILE_CASES' two kinked skins (the gas-driven mode's, the gradient mode's with the exposed
    slurry's own mode beneath) and a deep liquid pool at rest, which no profile resolves (its profile columns NaN but
    Girin's delta_m and V_s)."""
    from reentry_model.levelset.skinbody import LevelSetSkinBody
    mat = aa7075()
    (a, u_a), (b, u_b) = PROFILE_CASES[0][0], PROFILE_CASES[1][0]
    pool = [LIQUID] * 30 + [RIGID] * 10
    mk = row_markers(3, h=0.05)
    sk = moving(mk, [a, b, pool], [u_a, u_b, np.zeros(len(pool))], mat)
    flow = gas(3, tau=200.0, G=-3e4)
    theta = np.array([0.2, 0.3, 0.4])
    plan = spraying.spray_plan(sk, mk, flow, STATE, model_of(mat, dispersion="profile"), theta <= 0.5 * np.pi,
                               flow.deceleration * np.cos(theta), 0.5, 0.05, 0.04, slurry=True, theta=theta,
                               division=row_division(3, h=0.05))
    assert np.array_equal(plan.resolved, [True, True, False])
    for _ in range(PROFILE_STEPS):
        plan.release(sk, 1e-4)
    block = plan.records(0.5)
    width = len(spray.SOURCE_COLUMNS)
    assert block.shape == (4, width + len(spray.PROFILE_COLUMNS)) and block.shape[1] == 36
    liquid, slurry = np.flatnonzero(plan.mass > 0.0), np.flatnonzero(plan.mass_slurry > 0.0)
    assert np.array_equal(liquid, [0, 1, 2]) and np.array_equal(slurry, [1])
    np.testing.assert_array_equal(block[:3, width:], plan.profile[liquid])
    assert np.isnan(block[2, width:width + 7]).all() and np.isfinite(block[:3, width + 7:]).all()
    blank = [spray.PROFILE_COLUMNS.index("release_depth_m"), spray.PROFILE_COLUMNS.index("growth_rate_per_s")]
    kept = [j for j in range(len(spray.PROFILE_COLUMNS)) if j not in blank]
    assert np.isfinite(block[1, width + np.array(blank)]).all() and np.isnan(block[3, width + np.array(blank)]).all()
    np.testing.assert_array_equal(block[3, width + np.array(kept)], plan.profile[1, kept])
    body = LevelSetSkinBody.__new__(LevelSetSkinBody)
    for model, names in ((model_of(mat, dispersion="profile"), spray.record_columns(36)),
                         (model_of(mat), spray.record_columns(27)), (None, spray.record_columns(27))):
        body.spray_model = model
        assert body.record_names() == names
    body.spray_model, body.source_rows = model_of(mat, dispersion="profile"), [block]
    history = History({"time_s": np.array([0.5])})
    cols = np.load(droplets.write_particles(str(tmp_path), body, history, window=1.0)["particles"])
    assert list(cols.files) == spray.record_columns(36)


def test_the_step_summary_gains_the_dispersion_columns():
    """On the profile path the step summary ends with the shared melt_profile.DISPERSION_COLUMNS, the summary of the
    plan's modes and record (`melt_profile.summary`, the skyscrapers' `dispersion_stats`); a step that made no plan gives
    an empty step's (the body's `_summary`), so every step of a profile run has the same columns. The table path's and a
    body without a spray model have none of them."""
    from reentry_model.levelset.skinbody import LevelSetSkinBody
    from reentry_model.shared import melt_profile
    mat = aa7075()
    n = len(SLURRY_COLUMNS)
    theta = np.linspace(0.05, 0.9, n)
    summaries = {}
    for dispersion in spray.DISPERSION_NAMES:
        mk = row_markers(n)
        sk = stacked(mk, SLURRY_COLUMNS, mat, rng=np.random.default_rng(17))
        flow = gas(n, tau=200.0, G=-3e4)
        plan = spraying.spray_plan(sk, mk, flow, STATE, model_of(mat, dispersion=dispersion), theta <= 0.5 * np.pi,
                                   flow.deceleration * np.cos(theta), 0.5, 0.05, 0.04, theta=theta,
                                   division=row_division(n))
        plan.release(sk, 0.1)
        plan.records(0.5)
        summaries[dispersion] = (plan, plan.step_summary(sk))
    plan, out = summaries["profile"]
    want = melt_profile.summary(plan.modes, plan.record, plan.result.branch, plan.mass)
    assert list(out)[-len(melt_profile.DISPERSION_COLUMNS):] == list(melt_profile.DISPERSION_COLUMNS)
    assert {k: out[k] for k in want} == pytest.approx(want, nan_ok=True)
    assert out["sprayed_gradient_kg"] + out["sprayed_gas_driven_kg"] > 0.0
    table = summaries["table"][1]
    assert not set(melt_profile.DISPERSION_COLUMNS) & set(table)
    assert set(out) - set(melt_profile.DISPERSION_COLUMNS) == set(table)
    body = LevelSetSkinBody.__new__(LevelSetSkinBody)
    body.skins, body.material, body.time, body.n_released = sk, mat, 0.5, 0.0
    for model, has in ((model_of(mat, dispersion="profile"), True), (model_of(mat), False), (None, False)):
        body.spray_model = model
        empty = body._summary(None, None)
        assert set(melt_profile.DISPERSION_COLUMNS) <= set(empty) if has else \
            not set(melt_profile.DISPERSION_COLUMNS) & set(empty)
        if has:
            assert {k: empty[k] for k in melt_profile.DISPERSION_COLUMNS} == pytest.approx(
                melt_profile.summary(None, None, None, None), nan_ok=True)

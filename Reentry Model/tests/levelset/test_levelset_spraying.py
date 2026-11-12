"""Tests of the level set's spraying per marker (design sections 9.5 and 20.6; level-set stage 3 Task 8): Girin's plan on
the markers' skins through the shared `SprayModel`, the release cascade from the skins' top cells, the droplet records
through the shared writer, and Girin and Kopyt's 1994 front-surface case through the level-set path. The skyscrapers'
`skin.Skins`, `cascade_mass` and `release_top` are imported here only to compare against (decision D4): the package
itself never imports them."""
import csv
import json
import os
from types import SimpleNamespace

import numpy as np
import pytest

from reentry_model.levelset import spraying
from reentry_model.levelset.markers import Faces, Markers
from reentry_model.levelset.skin import Skins, SkinSettings
from reentry_model.shared import droplets, spray, surface_flow as sf
from reentry_model.shared.lubrication import lubrication
from reentry_model.shared.material import LiquidProperties, Material

from helpers import PROJECT_ROOT
from materials import LinearMeltMaterial

LIQUID, SLURRY, RIGID = 950.0, 902.0, 880.0    # K: all liquid; 30 % solid, flowing; rigid (AA7075's T_rigid 895.1 K)
CELL = 1.0e-5                                  # m: each stacked cell holds 10 um of solid metal's mass
SETTINGS = SkinSettings(cell=10e-6, fine_zone=0.4e-3, stub=0.4e-3)
H = 1.0e-3                                     # m: the markers' spacing on the flat grid
SOLID_TOL = 1.0e-9                             # the D4 test's tolerance on the solid released (its docstring)
TAU_RENEWAL = 20.0                             # Pa: the gas shear of the renewal test (supercritical: We_s above 3.08
                                               # on the thick branch)


def aa7075():
    return Material.from_drama_json("AA7075_scheil")


def row_markers(n, h=H, order=None):
    """n markers on a straight row h apart in z = 0, normals +z, square cells of area h^2, a face of length h between
    neighbours; `order` the skin row of each marker (default the identity)."""
    x = np.column_stack([np.arange(n) * h, np.zeros(n), np.zeros(n)])
    skin_index = np.arange(n) if order is None else np.asarray(order)
    mk = Markers(x, np.tile([0.0, 0.0, 1.0], (n, 1)), np.full(n, h * h), np.zeros((n, 2)), skin_index,
                 np.full(n, -1), h)
    i, j = np.arange(n - 1), np.arange(1, n)
    mk.faces = Faces(i, j, np.full(n - 1, h), np.tile([1.0, 0.0, 0.0], (n - 1, 1)))
    return mk


def stacked(mk, temps, mat=None, rng=None):
    """Skins on the markers: marker i's skin holds a cell of 10 um of solid metal's mass at each temperature of temps[i],
    top first; every cell moving with a random velocity when `rng` is given."""
    mat = aa7075() if mat is None else mat
    order = np.argsort(mk.skin_index)                    # the marker of each skin row
    sk = Skins(mat, mk.area[order], mk.kappa[order], SETTINGS)
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
    plan = spraying.spray_plan(sk, mk, flow, STATE, model, windward, w_n, 0.5, 0.05, 0.05, theta=theta)
    # the same arrays by hand, in the markers' order
    r = mk.skin_index
    F, L = sk.depth_above(mat.T_rigid)[r], sk.depth_above(mat.T_liquidus)[r]
    T = sk.T()
    flowing, liquid = sk.flowing(T), sk.liquid(T)
    n_flow, n_liq = np.cumprod(flowing, axis=1).sum(axis=1)[r], np.cumprod(liquid, axis=1).sum(axis=1)[r]
    delta_m = spray.spray_layer(flow, model.liquid, model.slip_band, theta, 0.05)
    v_s = lubrication(flow.tau, flow.G, F, delta_m, model.liquid.mu, b_layer=F)[0]
    A = sk.area[r]
    extent = spray.region_extent(np.column_stack([mk.faces.i, mk.faces.j]), L > 0.0, A)
    m_f = (sk.m * np.cumprod(liquid, axis=1))[r].sum(axis=1)
    want = model.evaluate(flow, STATE, F, delta_m, v_s, windward, 0.5, A, m_f, 0.05, b_layer=F, extent=extent,
                          deceleration_n=w_n, regime_layer=F, on_slurry=n_flow > n_liq, rt_layer=L, cap_depth=L)
    assert np.array_equal(plan.result.mdot, want.mdot) and np.array_equal(plan.rate, want.mdot)
    assert np.array_equal(plan.result.r, want.r, equal_nan=True)
    assert np.array_equal(plan.result.branch, want.branch)
    assert (want.mdot > 0.0).sum() >= 2 and len(set(want.branch.tolist())) >= 2    # the case is not trivial


def test_the_slip_band_needs_the_markers_angles():
    """With the model's slip band on, the band's depth reads the markers' angles: a plan without them is refused rather
    than silently given Girin's delta_m."""
    mat = aa7075()
    n = len(MIXED)
    mk = row_markers(n)
    sk = stacked(mk, MIXED, mat)
    model = model_of(mat, slip_band=True)
    with pytest.raises(ValueError, match="theta"):
        spraying.spray_plan(sk, mk, gas(n), STATE, model, np.ones(n, bool), np.full(n, 30.0), 0.5, 0.05, 0.05)


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
    theirs = sky.Skins(mat, ours.area.copy())
    theirs._install(ours.m.copy(), ours.H.copy(), ours.count.copy(), U=ours.U.copy())
    flow, model = gas(n, tau=200.0, G=-3e4), model_of(mat)
    theta = np.linspace(0.05, 0.9, n)
    windward, w_n = theta <= 0.5 * np.pi, flow.deceleration * np.cos(theta)
    plan = spraying.spray_plan(ours, mk, flow, STATE, model, windward, w_n, 0.5, 0.05, 0.04, slurry=True,
                               theta=theta)
    F = theirs.depth_above(mat.T_rigid)
    flowing = theirs.valid & np.isfinite(sky.cell_viscosity(theirs, 0.0))
    liquid = flowing & (theirs.T >= mat.T_liquidus)
    n_flow, n_liq = np.cumprod(flowing, axis=1).sum(axis=1), np.cumprod(liquid, axis=1).sum(axis=1)
    edges = np.column_stack([mk.faces.i, mk.faces.j])
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
    theirs = sky.Skins(mat, ours.area.copy())
    theirs._install(ours.m.copy(), ours.H.copy(), ours.count.copy(), U=ours.U.copy())
    flow, model = gas(1, tau=200.0, G=-3e4), model_of(mat)
    theta = np.array([0.2])
    windward, w_n = np.array([True]), flow.deceleration * np.cos(theta)
    plan = spraying.spray_plan(ours, mk, flow, STATE, model, windward, w_n, 0.5, 0.05, 0.04, slurry=True, theta=theta)
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


TOP_ZONE = [[LIQUID] * 40, [LIQUID] * 10 + [SLURRY] * 30]      # top-zone skins: the fine zone over thick melt


def test_a_top_zone_skin_sprays_by_its_depth_from_phi_to_psi():
    """Design 20.6 and decision D2 (the controller's ruling (b) of Task 10): on a top-zone skin -- the top 0.4 mm over
    thick melt the volume holds -- the plan reads the depth from phi to psi (b_mode "phi_psi") for the melt depth b (the
    regime test and the lubrication's layer with it), the cap depth, the Rayleigh-Taylor pool and the slurry's layer, its
    cap and the extent its wave must fit in (b at least B_MIN); a full-depth skin beside it (nan) reads its own. The
    liquid's mode is the shared model's called with those depths and the skin's own film mass; the release takes only
    the skin's own top cells; and the skin's base draw then refills it from the melt beneath, at the melt's temperature,
    to its fine zone's mass."""
    from reentry_model.levelset import modes
    mat = aa7075()
    cols = MIXED[:2] + TOP_ZONE
    n = len(cols)
    top = np.array([False, False, True, True])
    d = np.array([np.nan, np.nan, 3.0e-3, 2.0e-3])                        # phi to psi at the top-zone markers [m]
                                                                          # (read only where top)
    mk = row_markers(n, h=0.05)
    sk = stacked(mk, cols, mat)
    flow, model = gas(n, tau=TAU_RENEWAL), model_of(mat)
    theta = np.full(n, 0.3)
    windward, w_n = np.ones(n, bool), flow.deceleration * np.cos(theta)
    plan = spraying.spray_plan(sk, mk, flow, STATE, model, windward, w_n, 0.5, 0.05, 0.05, slurry=True,
                               b_mode="phi_psi", depth_phi_psi=d, top=top, theta=theta)
    # the mode is given, not read off the depths: without it, or with no depth at a top-zone marker, the plan refuses
    with pytest.raises(ValueError, match="top"):
        spraying.spray_plan(sk, mk, flow, STATE, model, windward, w_n, 0.5, 0.05, 0.05, b_mode="phi_psi",
                            depth_phi_psi=d, theta=theta)
    with pytest.raises(ValueError, match="depth"):
        spraying.spray_plan(sk, mk, flow, STATE, model, windward, w_n, 0.5, 0.05, 0.05, b_mode="phi_psi",
                            depth_phi_psi=np.where(top, np.nan, 1e-3), top=top, theta=theta)
    F, L = sk.depth_above(mat.T_rigid), sk.depth_above(mat.T_liquidus)
    b = np.where(top, d, F)
    np.testing.assert_array_equal(plan.flowing, b)
    np.testing.assert_array_equal(plan.cap_depth, b)                      # with the slurry: the flowing depth, or d
    T = sk.T()
    flowing, liquid = sk.flowing(T), sk.liquid(T)
    run = np.cumprod(flowing, axis=1).astype(bool)
    delta_m = spray.spray_layer(flow, model.liquid, model.slip_band, theta, 0.05)
    v_s = lubrication(flow.tau, flow.G, b, delta_m, model.liquid.mu, b_layer=b)[0]
    extent = spray.region_extent(np.column_stack([mk.faces.i, mk.faces.j]), L > 0.0, sk.area)
    n_flow, n_liq = np.cumprod(flowing, axis=1).sum(axis=1), np.cumprod(liquid, axis=1).sum(axis=1)
    want = model.evaluate(flow, STATE, b, delta_m, v_s, windward, 0.5, sk.area, (sk.m * run).sum(axis=1), 0.05,
                          b_layer=b, extent=extent, deceleration_n=w_n, regime_layer=b, on_slurry=n_flow > n_liq,
                          rt_layer=np.where(top, d, L), cap_depth=b)
    np.testing.assert_array_equal(plan.result.mdot, want.mdot)
    np.testing.assert_array_equal(plan.result.branch, want.branch)
    assert plan.rate[2] > 0.0 and plan.result.branch[2] == spray.BRANCH_THICK      # 3 mm of melt: the thick branch
    # the slurry exposed on the second top-zone skin: its layer runs from the slurry cell's top down to psi
    dz = sk.dz(T)
    assert plan.on_slurry[3]
    assert plan.slurry["depth"][3] == pytest.approx(d[3] - dz[3, :10].sum(), rel=1e-12)
    # the release takes the skin's own top cells only, and the base draw refills it from the melt beneath
    M0 = sk.m.sum(axis=1)
    released = sum(plan.release(sk, 0.0125).mass for _ in range(8))
    assert released[2] > 0.0 and np.all(released <= M0 * (1 + 1e-14))
    keep = modes.minimum_mass(sk, top)
    assert sk.m[2].sum() < keep[2]
    drawn, dH = sk.draw(LIQUID, 0.0, keep, cap=modes.draw_cap(sk, top))
    assert drawn[2] > 0.0 and sk.m[2].sum() == pytest.approx(keep[2], rel=1e-14)
    assert dH[2] / drawn[2] == pytest.approx(float(sk.props.enthalpy(np.array([LIQUID]))[0]), rel=1e-12)


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
        theirs = sky.Skins(mat, ours.area.copy())
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
                               0.05, 0.05, theta=theta)
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
                                   0.05, slurry=slurry, theta=theta)
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
                                   v_hat=np.array([0.0, 0.0, -1.0]))
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


def test_girin_cases_through_the_level_set_path():
    """Girin and Kopyt (1994) Table 2, the front-surface Rayleigh-Taylor mode, fed through `spray_plan` on one flat
    marker at the stagnation point holding a pool of their liquid one wavelength deep (so tanh(k h) is 1 to 1e-5), off
    Girin's closure (no shear mode): the plan's wavelength is their lambda* (Eq. 14) within 2 % and its growth time
    their tau* within 1 %, at each of the three decelerations."""
    with open(os.path.join(PROJECT_ROOT, "data", "reference_values", "girin1994_tables.json")) as fh:
        t = json.load(fh)["table2"]
    mat = LinearMeltMaterial(rho=Iron.rho, rho_l=Iron.rho, mu=Iron.mu, sigma=Iron.sigma)
    model = spray.SprayModel(LiquidProperties(Iron.rho, Iron.mu, Iron.sigma), spray_gates="2026-11-10")
    for W, lam_ref, tau_ref in zip(t["W_ms2"], t["lambda_star_m_eq14"], t["tau_star_s"]):
        mk = row_markers(1, h=1.0)                                     # 1 m2: the pool admits the fastest wave
        sk = Skins(mat, mk.area, mk.kappa, SkinSettings(cell=20e-6, fine_zone=0.4e-3, stub=0.4e-3))
        sk.fill(Iron.rho * lam_ref, 1000.0)                            # one wavelength deep, liquid
        flow = gas(1, closure=sf.CLOSURE_COUETTE, tau=0.0, decel=W)
        theta = np.zeros(1)
        plan = spraying.spray_plan(sk, mk, flow, STATE, model, np.array([True]), np.array([W]), 0.5, 0.05, 0.05,
                                   theta=theta)
        res = plan.result
        assert res.branch[0] == spray.BRANCH_RT and plan.rate[0] > 0.0
        assert res.wavelength[0] == pytest.approx(lam_ref, rel=2e-2)
        assert res.growth[0] == pytest.approx(tau_ref, rel=1e-2)

"""spray.py: delta_m, the thin-film and Rayleigh-Taylor modes against Girin & Kopyt (1994), the thick branch, release
bookkeeping, size caps and histograms."""
import json
import os

import numpy as np
import pytest

from reentry_model import dispersion, material, spray, surface_flow as sf
from helpers import REPO_ROOT

LIQ = material.LiquidProperties(2400.0, 1.3e-3, 0.86)


class Iron:
    rho, mu, sigma = 7800.0, 5.8e-3, 1.2


def test_thin_film_mode_reproduces_the_1994_table_1_scaling():
    """r_d = lambda*/4 and tau_d = tau* for the six (rho_2, V0) cases: one shock-layer factor on rho_2 V0^2 (fitted on
    the first entry, 7.8) reproduces every r_d within 0.5 % and tau_d within 1 %; the printed mass rates are
    rho_1 r_d / (2 tau_d) within 1 %."""
    t = json.load(open(os.path.join(REPO_ROOT, "data", "reference_values", "girin1994_tables.json")))["table1"]
    rho2, V0 = np.array(t["rho2_kgm3"]), np.array(t["V0_ms"])
    lam, tau = spray.thin_film_mode(np.full((3, 2), t["mach"]), rho2[:, None] * V0[None, :] ** 2, Iron)
    r_tab, tau_tab, m_tab = (np.array(t[k]) for k in ("r_d_m", "tau_d_s", "m_kgm2s"))
    f_s = (lam / 4.0)[0, 0] / r_tab[0, 0]
    assert 7.5 < f_s < 8.1
    assert np.abs(lam / 4.0 / f_s / r_tab - 1.0).max() < 5e-3 and np.abs(tau / f_s ** 1.5 / tau_tab - 1.0).max() < 1e-2
    assert np.abs(Iron.rho * (lam / 4.0 / f_s) / (2.0 * tau / f_s ** 1.5) / m_tab - 1.0).max() < 1e-2
    assert spray.CAPILLARY_TAU == pytest.approx(0.798, rel=1e-3)                                  # tau* = 2 x 2 pi / omega_cap
    assert np.isinf(spray.thin_film_mode(np.array([3.0]), np.array([0.0]), Iron)[0][0])


def test_rayleigh_taylor_mode_reproduces_the_1994_table_2():
    t = json.load(open(os.path.join(REPO_ROOT, "data", "reference_values", "girin1994_tables.json")))["table2"]
    for W, lam_ref, tau_ref in zip(t["W_ms2"], t["lambda_star_m_eq14"], t["tau_star_s"]):
        active, lam, tau = spray.rayleigh_taylor(W, np.array([1e-3, 1e-5]), Iron)
        assert lam[0] == pytest.approx(lam_ref, rel=2e-2) and tau[0] == pytest.approx(tau_ref, rel=1e-2)
        assert active[0] == (W * 1e-6 * Iron.rho > 3.0 * Iron.sigma) and not active[1]
    assert not spray.rayleigh_taylor(30.0, np.array([1e-4]), LIQ)[0][0]                          # our deceleration: inactive


class FakeFlow:
    """The fields spray.py reads: a shock-layer branch with the Girin closure unless told otherwise."""

    def __init__(self, n, closure=sf.CLOSURE_GIRIN, branch=sf.BRANCH_SHOCK_LAYER, u=2000.0, rho=7e-4, mu=1.6e-4,
                 delta_a=7e-3, mach=1.0, decel=30.0):
        self.closure = np.full(n, closure)
        self.branch = branch
        self.u_eff = self.u_e = np.full(n, u)
        self.rho_e, self.mu_e, self.delta_a, self.mach_e = np.full(n, rho), np.full(n, mu), np.full(n, delta_a), np.full(n, mach)
        self.deceleration = decel


class FakeState:
    class freestream:
        rho = 8e-5
    V = 7200.0
    ma = 24.0
    h = 71e3


def test_melt_layer_and_branches():
    flow = FakeFlow(3)
    delta_m, factor = spray.melt_layer(flow, LIQ)
    alpha, mu = 7e-4 / 2400.0, 1.6e-4 / 1.3e-3
    assert np.allclose(delta_m, (alpha / mu ** 2) ** (1 / 3) * 7e-3) and np.allclose(factor, (alpha * mu) ** (1 / 3) / (1 + (alpha * mu) ** (1 / 3)))
    assert 5e-5 < delta_m[0] < 2e-4                                                                 # ~0.1 mm (spec estimate)
    model = spray.SprayModel(LIQ)
    b = np.array([1e-3, 2e-5, 0.0])                                                                 # thick, thin, dry
    v_s = np.array([5.0, 20.0, 0.0])          # both films supercritical: Girin's criterion now gates every branch
    areas, m_f = np.full(3, 1e-5), b * LIQ.rho * 1e-5
    res = model.evaluate(flow, FakeState(), b, delta_m, v_s, np.array([True, True, True]), 0.5, areas, m_f)
    assert res.branch.tolist() == [spray.BRANCH_THICK, spray.BRANCH_THIN, -1]
    we_s = LIQ.rho * 25.0 * delta_m[0] / LIQ.sigma
    assert res.we_s[0] == pytest.approx(we_s) and res.unstable[0] == (we_s > 4.62)
    d_f, im_f, _ = dispersion.DispersionTable()(we_s)
    lam = 2 * np.pi * delta_m[0] / d_f
    assert res.r[0] == pytest.approx(min(0.17 * lam, (3 * m_f[0] / (4 * np.pi * LIQ.rho)) ** (1 / 3), 0.0125))
    assert res.mdot[0] == pytest.approx(LIQ.rho * np.pi * (0.17 * lam) ** 2 / (lam * 1.1 * delta_m[0] / (5.0 * im_f)))
    lam_t, tau_t = spray.thin_film_mode(np.array([1.0]), np.array([7e-4 * 2000.0 ** 2]), LIQ)
    assert res.r[1] == pytest.approx(min(lam_t[0] / 4.0, (3 * m_f[1] / (4 * np.pi * LIQ.rho)) ** (1 / 3), 0.0125))
    assert res.mdot[1] == pytest.approx(LIQ.rho * min(2e-5, lam_t[0] / 8.0) / tau_t[0])
    assert res.dm[1] == pytest.approx(min(res.mdot[1] * 1e-5 * 0.5, m_f[1])) and res.dm[2] == 0.0
    assert res.dn[1] == pytest.approx(res.dm[1] / (4 / 3 * np.pi * LIQ.rho * res.r[1] ** 3))
    assert not res.rt_active
    # Couette closure (the gate declined to certify the patch): no delta_m, so the thin mode with the edge state
    flow_c = FakeFlow(2, closure=sf.CLOSURE_COUETTE)
    dm_c, factor_c = spray.melt_layer(flow_c, LIQ)
    assert np.all(np.isnan(dm_c)) and np.all(factor_c == 0.0)
    res = model.evaluate(flow_c, FakeState(), np.array([2e-5, 1e-3]), dm_c, np.array([20.0, 20.0]), np.array([True, True]), 0.5, areas[:2], m_f[:2] + 1e-9)
    assert res.branch.tolist() == [spray.BRANCH_THIN, spray.BRANCH_THIN]          # a thick film cannot take Girin's branch here
    # free-molecular branch: no edge state exists, the freestream drives the mode
    flow_fm = FakeFlow(2, closure=sf.CLOSURE_COUETTE, branch=sf.BRANCH_FREE_MOLECULAR)
    res = model.evaluate(flow_fm, FakeState(), np.array([2e-5, 1e-3]), dm_c, np.array([20.0, 20.0]), np.array([True, True]), 0.5, areas[:2], m_f[:2] + 1e-9)
    assert res.branch.tolist() == [spray.BRANCH_RAREFIED, spray.BRANCH_RAREFIED]
    lam_r, _ = spray.thin_film_mode(np.array([24.0]), np.array([8e-5 * 7200.0 ** 2]), LIQ)
    assert res.r[0] == pytest.approx(min(lam_r[0] / 4.0, (3 * (m_f[0] + 1e-9) / (4 * np.pi * LIQ.rho)) ** (1 / 3), 0.0125))
    # leeward or dry: nothing
    res = model.evaluate(flow, FakeState(), b, delta_m, v_s, np.array([False, False, False]), 0.5, areas, m_f)
    assert not res.unstable.any() and res.dm.sum() == 0.0
    with pytest.raises(ValueError):
        spray.SprayModel(LIQ, k_r=0.0)


def test_source_rows_and_histogram():
    flow = FakeFlow(2)
    delta_m, _ = spray.melt_layer(flow, LIQ)
    b = np.array([1e-3, 1e-3])
    m_f = b * LIQ.rho * 1e-5
    res = spray.SprayModel(LIQ).evaluate(flow, FakeState(), b, delta_m, np.array([5.0, 5.0]), np.array([True, True]), 0.5, np.full(2, 1e-5), m_f)
    rows = spray.source_rows(10.0, 71e3, 7200.0, np.radians([30.0, 60.0]), np.zeros((2, 3)), np.tile([0.0, 1.0, 0.0], (2, 1)), flow, res, b, LIQ, FakeState())
    assert len(rows) == 2 and len(rows[0]) == len(spray.SOURCE_COLUMNS)
    row = dict(zip(spray.SOURCE_COLUMNS, rows[0]))
    assert row["time_s"] == 10.0 and row["altitude_km"] == 71.0 and row["theta_deg"] == pytest.approx(30.0) and row["branch"] == spray.BRANCH_THICK
    assert row["closure"] == sf.CLOSURE_GIRIN
    assert row["we_d"] == pytest.approx(8e-5 * 7200.0 ** 2 * 2 * row["r_m"] / LIQ.sigma) and row["oh"] == pytest.approx(LIQ.mu / np.sqrt(LIQ.rho * LIQ.sigma * 2 * row["r_m"]))
    assert row["breakup"] == float(row["we_d"] > 12.0) and row["dm_kg"] == res.dm[0] and row["v_s_ms"] == 5.0
    n_hist, m_hist = spray.histogram(np.array([2e-6, 5e-5, 5e-5, np.nan, 2e-2]), np.array([1.0, 2.0, 3.0, 4.0, 5.0]), np.array([1.0, 1.0, 1.0, 1.0, 1.0]))
    assert n_hist.sum() == 6.0 and m_hist.sum() == 3.0 and n_hist.shape == (spray.N_BINS,)          # nan and out-of-range dropped
    assert spray.BIN_EDGES[0] == 1e-6 and spray.BIN_EDGES[-1] == 1e-2


def test_girins_critical_weber_number_gates_the_thin_branch_too():
    """Girin & Kopyt's inviscid side mode is unstable at every wavelength and carries no threshold of its own, so the
    thin branch used to strip film wherever a wavelength existed -- at any angle, however small the shear. Girin's
    (2017) criterion We_s > We_cr now gates it as well, measured over the liquid the shear can actually reach,
    min(delta_m, layer): on a film shallower than the conjugate layer that is the film's own depth."""
    flow = FakeFlow(2, closure=sf.CLOSURE_COUETTE)                # no delta_m: the thin branch
    model = spray.SprayModel(LIQ)
    delta_m, _ = spray.melt_layer(flow, LIQ)
    b = np.array([2e-5, 2e-5])
    areas, m_f = np.full(2, 1e-5), b * LIQ.rho * 1e-5
    v_s = np.array([0.5, 20.0])                                   # subcritical, supercritical
    res = model.evaluate(flow, FakeState(), b, delta_m, v_s, np.array([True, True]), 0.5, areas, m_f)
    we = LIQ.rho * v_s ** 2 * b / LIQ.sigma                       # layer = b here, and delta_m is nan
    assert res.we_s == pytest.approx(we) and we[0] < 4.62 < we[1]
    assert res.branch.tolist() == [spray.BRANCH_THIN, spray.BRANCH_THIN]          # classified either way
    assert not res.unstable[0] and res.dm[0] == 0.0 and np.isnan(res.r[0])        # but stable: nothing leaves
    assert res.unstable[1] and res.dm[1] > 0.0                                    # and the supercritical one sprays
    lam, tau = spray.thin_film_mode(np.array([1.0]), np.array([7e-4 * 2000.0 ** 2]), LIQ)
    assert np.isfinite(lam[0]) and np.isfinite(tau[0])            # a wavelength existed for both -- only We decided


def test_a_wave_longer_than_its_melt_region_cannot_form():
    """A mode cannot develop in a puddle narrower than itself, and the domain is the contiguous molten region, not the
    mesh facet. The thin-branch wavelength diverges as the edge dynamic pressure falls, so near the stagnation point it
    is the fit that has to decide (2026-09-24)."""
    flow = FakeFlow(2, closure=sf.CLOSURE_COUETTE)
    model = spray.SprayModel(LIQ)
    delta_m, _ = spray.melt_layer(flow, LIQ)
    b, v_s = np.full(2, 2e-5), np.full(2, 20.0)                   # both supercritical
    areas, m_f = np.full(2, 1e-5), b * LIQ.rho * 1e-5
    lam, _ = spray.thin_film_mode(np.array([1.0]), np.array([7e-4 * 2000.0 ** 2]), LIQ)
    tight, loose = 0.5 * lam[0], 2.0 * lam[0]
    res = model.evaluate(flow, FakeState(), b, delta_m, v_s, np.array([True, True]), 0.5, areas, m_f,
                         extent=np.array([tight, loose]))
    assert not res.unstable[0] and res.dm[0] == 0.0               # the wave is twice as wide as its melt
    assert res.unstable[1] and res.dm[1] > 0.0                    # here it fits
    assert spray.wave_fits(lam, None).all()                       # no extent given disables the test
    assert not spray.wave_fits(lam, np.zeros(1))[0]               # a region of no extent admits nothing


def test_the_rayleigh_taylor_driving_is_the_deceleration_normal_to_the_film():
    """Girin & Kopyt's dispersion relation carries the deceleration as W sin(Theta), Theta being the angle between the
    relative velocity of the media and the deceleration vector; on this body that is W cos(phi). It is the whole
    deceleration at the stagnation point and zero at the equator, where W lies in the surface and cannot lift the
    interface. Passing the whole deceleration instead over-drove every patch off the nose by 1/cos(phi)."""
    W, h = 95.0, np.full(4, 6e-3)          # deep enough that the flank still qualifies at cos 60 = 1/2:
    #                                       W cos(phi) h^2 rho_l > 3 Sigma needs h > 4.76 mm there, 3.37 mm at the nose
    phi = np.radians([0.0, 60.0, 90.0, 150.0])                    # nose, flank, equator, leeward
    w_n = W * np.maximum(np.cos(phi), 0.0)
    on, lam, tau = spray.rayleigh_taylor(w_n, h, LIQ)
    assert on[0] and on[1] and not on[2] and not on[3]            # driven at the nose, dead at and past the equator
    assert lam[0] == pytest.approx(2 * np.pi * np.sqrt(3 * LIQ.sigma / (W * LIQ.rho)))
    assert lam[1] == pytest.approx(lam[0] * np.sqrt(2.0))          # cos 60 = 1/2, and lambda* goes as W^-1/2
    assert lam[2] > 1e3                                           # cos 90 is round-off: kilometres of wavelength, no mode
    assert not np.isfinite(lam[3]) and not np.isfinite(tau[3])     # leeward, clamped to zero: no mode at all
    assert spray.rayleigh_taylor(W, h, LIQ)[0].tolist() == [True] * 4          # a scalar still broadcasts
    # and the bounded form takes the same normal component
    onb, lamb, taub = spray.rayleigh_taylor_bounded(w_n, h, np.full(4, 0.2), LIQ)
    assert onb[0] and not onb[2] and not onb[3] and taub[0] < taub[1] < np.inf


def test_the_front_surface_rayleigh_taylor_mode_is_applied_where_it_is_the_faster_mode():
    """Girin & Kopyt assign the front surface to their aperiodic Rayleigh-Taylor mode, not to the side mode, and with
    Girin's critical angle gating the shear branches the stagnation cap otherwise has no mechanism at all. The mode is
    applied where it passes *both* criteria -- the depth threshold W cos(phi) h^2 rho_l > 3 Sigma, which is what confines
    it to the nose, and an admissible wave in the molten region -- and where its growth time beats the shear mode's,
    since two modes of one interface cannot both break it up."""
    flow = FakeFlow(3, closure=sf.CLOSURE_COUETTE)            # no delta_m: the shear mode here is the weak thin one
    model = spray.SprayModel(LIQ)
    delta_m, _ = spray.melt_layer(flow, LIQ)
    b = np.full(3, 1e-4)
    layer = np.array([5e-3, 1e-3, 5e-3])                      # deep, too shallow for the threshold, deep
    v_s = np.array([0.1, 0.1, 40.0])                          # stable shear, stable shear, strong shear
    w_n = np.full(3, 95.0)                                    # at the stagnation point cos(phi) = 1
    areas, m_f = np.full(3, 1e-5), b * LIQ.rho * 1e-5
    deep, _, _ = spray.rayleigh_taylor(w_n, layer, LIQ)
    assert deep.tolist() == [True, False, True]               # h > lambda*/2 pi only for the deep pair
    res = model.evaluate(flow, FakeState(), b, delta_m, v_s, np.full(3, True), 0.5, areas, m_f,
                         b_layer=layer, deceleration_n=w_n)
    assert res.branch.tolist() == [spray.BRANCH_RT, spray.BRANCH_THIN, spray.BRANCH_THIN]
    assert res.dm[0] > 0.0 and res.dm[1] == 0.0 and res.dm[2] > 0.0
    _, lam_rt, tau_rt = spray.rayleigh_taylor_bounded(w_n, layer, np.full(3, np.inf), LIQ)
    assert res.growth[0] == pytest.approx(tau_rt[0])          # the winning mode's growth time is the RT one
    assert res.growth[2] < res.growth[0]                      # and the shear mode is the faster one on patch 2
    assert res.mdot[0] == pytest.approx(LIQ.rho * min(b[0], lam_rt[0] / 2.0) / tau_rt[0])
    cap = (3 * m_f[0] / (4 * np.pi * LIQ.rho)) ** (1 / 3)
    assert lam_rt[0] > 1e-2 and res.r[0] == pytest.approx(cap)   # lambda* is centimetres; the film-mass cap binds
    # with the mode reported but not applied -- every run before 2026-09-24 -- the deep stable patch releases nothing
    off = spray.SprayModel(LIQ, rt_spray=False).evaluate(
        flow, FakeState(), b, delta_m, v_s, np.full(3, True), 0.5, areas, m_f, b_layer=layer, deceleration_n=w_n)
    assert off.branch.tolist() == [spray.BRANCH_THIN] * 3 and off.dm[0] == 0.0 and off.dm[2] > 0.0
    assert off.dm[2] == pytest.approx(res.dm[2])              # the shear-dominated patch is untouched either way


def test_the_applied_rayleigh_taylor_mode_needs_the_wave_to_fit_its_melt_region():
    """The lateral bound applies to this mode as it does to the shear ones: a pool narrower than half the marginal
    wavelength admits nothing at all, however deep it is."""
    flow = FakeFlow(2, closure=sf.CLOSURE_COUETTE)
    model = spray.SprayModel(LIQ)
    delta_m, _ = spray.melt_layer(flow, LIQ)
    b, layer, v_s = np.full(2, 1e-4), np.full(2, 5e-3), np.full(2, 0.1)
    w_n = np.full(2, 95.0)
    areas, m_f = np.full(2, 1e-5), b * LIQ.rho * 1e-5
    lam_c = 2 * np.pi * np.sqrt(LIQ.sigma / (95.0 * LIQ.rho))          # the marginal wavelength
    res = model.evaluate(flow, FakeState(), b, delta_m, v_s, np.full(2, True), 0.5, areas, m_f,
                         b_layer=layer, extent=np.array([0.4 * lam_c, 20.0 * lam_c]), deceleration_n=w_n)
    assert res.branch.tolist() == [spray.BRANCH_THIN, spray.BRANCH_RT]   # too narrow, then wide enough
    assert res.dm[0] == 0.0 and res.dm[1] > 0.0


def test_regime_2_gives_a_deep_film_with_a_viscous_melt_the_step_profile_mode():
    """Regime 2: where delta_m has formed and the melt is the more viscous medium, the patch takes Girin's
    Kelvin-Helmholtz cell -- the 1994 step-profile mode -- and not regime 3's gradient instability.

    delta_m is supplied directly rather than through `melt_layer`, and it has to be: Eq. (2) puts the predicted melt
    boundary layer at metres for a melt this viscous, which is Girin's own observation that for a high-viscosity melt the
    layer "becomes comparable with the radius of the meteoroid remnant". On a body this size such a melt therefore never
    passes stage one and lands in regime 1 instead, which the next test measures. Here the branch itself is exercised on
    the state stage two would hand it."""
    nu_gas = 1.6e-4 / 7e-4                                     # the edge state FakeFlow carries, 0.229 m2/s
    MELT_VISCOUS = material.LiquidProperties(2400.0, 10.0 * nu_gas * 2400.0, 0.86)   # nu_melt = 10 x nu_gas
    assert MELT_VISCOUS.mu / MELT_VISCOUS.rho > nu_gas
    flow = FakeFlow(2)
    model = spray.SprayModel(MELT_VISCOUS)
    delta_m = np.full(2, 1e-4)                                 # a conjugate depth the film can exceed
    b, v_s = np.full(2, 1e-3), np.full(2, 10.0)                # deep film; v_s large enough to clear We_cr
    areas, m_f = np.full(2, 1e-5), np.full(2, 1e-3 * 2400.0 * 1e-5)
    res = model.evaluate(flow, FakeState(), b, delta_m, v_s, np.full(2, True), 0.5, areas, m_f)
    assert res.branch.tolist() == [spray.BRANCH_REGIME2] * 2    # not BRANCH_THICK, and the RT mode does not take it
    assert res.we_s[0] > model.we_critical                     # Girin's gate still applies, as on every branch
    # the wavelength and growth time are Girin & Kopyt's Eqs. (11) and (12) -- the thin branch's pair, same mode
    lam = 1.5 * flow.mach_e[0] * MELT_VISCOUS.sigma / (flow.rho_e[0] * flow.u_eff[0] ** 2)
    tau = spray.CAPILLARY_TAU * lam ** 1.5 * np.sqrt(MELT_VISCOUS.rho / MELT_VISCOUS.sigma)
    assert res.growth[0] == pytest.approx(tau, rel=1e-12)
    assert res.r[0] == pytest.approx(spray.K_R * lam, rel=1e-12)            # r = k_r lambda*, below both caps here
    # the rate is Girin 2017's torus shedding on that wavelength, with no min(b, lambda*/8) factor: one torus of
    # cross-section r = k_r lambda* per wavelength per growth time, mdot = rho_l pi r^2 / (lambda* tau*)
    mdot = MELT_VISCOUS.rho * np.pi * (spray.K_R * lam) ** 2 / (lam * tau)
    assert res.mdot[0] == pytest.approx(mdot, rel=1e-12)
    # "not depth-limited" means the rate does not reference the film depth at all, as on the thick branch: doubling the
    # film leaves it unchanged, where the thin branch's min(b, lambda*/8) factor follows a film shallower than the wave
    deeper = model.evaluate(flow, FakeState(), 2.0 * b, delta_m, v_s, np.full(2, True), 0.5, areas, 2.0 * m_f)
    assert deeper.branch.tolist() == [spray.BRANCH_REGIME2] * 2
    assert deeper.mdot[0] == pytest.approx(res.mdot[0], rel=1e-12)
    # mdot is the instability's demand; the film is what limits the release, exactly as on every other branch
    assert mdot * areas[0] * 0.5 > m_f[0] and res.dm[0] == pytest.approx(m_f[0])


def test_a_kelvin_helmholtz_melt_fails_girins_first_stage_on_this_body():
    """A melt viscous enough for regime 2 cannot satisfy stage one here, so it routes to regime 1 -- correctly.

    delta_m/delta_a = (alpha/mu^2)^(1/3) = ((nu_melt/nu_gas)^2 rho_l/rho_e)^(1/3), so at the viscosity threshold itself
    the predicted melt boundary layer is already (rho_l/rho_e)^(1/3) ~ 151 air boundary-layer thicknesses, of order a
    metre at delta_a = 7 mm, and no film a 100 mm sphere can hold accommodates it. Stage one then says the layer never
    formed, which is Girin's dominant-ablation case and his Girin & Kopyt (1994) mode: the thin branch. That is what his
    own ordering produces and not a gap in it -- recorded so the regime-2 branch is not later mistaken for dead code that
    ought to have fired."""
    flow = FakeFlow(1)
    nu_gas = flow.mu_e[0] / flow.rho_e[0]
    for ratio in (1.0, 10.0, 100.0):
        liq = material.LiquidProperties(2400.0, ratio * nu_gas * 2400.0, 0.86)
        delta_m, _ = spray.melt_layer(flow, liq)
        expected = (ratio ** 2 * liq.rho / flow.rho_e[0]) ** (1.0 / 3.0)            # delta_m / delta_a
        assert delta_m[0] / flow.delta_a[0] == pytest.approx(expected, rel=1e-12)
        assert delta_m[0] > 1.0                                  # metres of liquid would be needed to out-deepen it
    # at the threshold itself the factor is (rho_l/rho_e)^(1/3), and no film the body can hold reaches it
    liq = material.LiquidProperties(2400.0, nu_gas * 2400.0, 0.86)
    delta_m, _ = spray.melt_layer(flow, liq)
    assert delta_m[0] / flow.delta_a[0] == pytest.approx((liq.rho / flow.rho_e[0]) ** (1.0 / 3.0), rel=1e-12)
    b = np.full(1, 1e-2)                                         # 10 mm of melt: a tenth of the sphere's radius
    res = spray.SprayModel(liq).evaluate(flow, FakeState(), b, delta_m, np.full(1, 10.0), np.full(1, True), 0.5,
                                        np.full(1, 1e-5), b * liq.rho * 1e-5)
    assert res.branch.tolist() == [spray.BRANCH_THIN]            # regime 1, the stabilised-core case, as Girin has it


def test_regime_2_never_fires_for_liquid_aluminium():
    """The regression guard: for the real melt the branch is inert, so every earlier run is reproduced exactly.

    Liquid aluminium has nu_melt = 5.42e-7 m2/s against a post-shock edge nu_gas measured at 0.018-3.75 m2/s over the
    two physics flights, a ratio of 3.3e4 to 6.9e6, so the melt is far the less viscous medium everywhere and a formed
    layer is always regime 3 -- Girin's gradient instability, the thick branch -- exactly as before 2026-09-25 (0 of
    533 960 windward wet patch-steps on the 100 mm flight and 0 of 42 268 on the 50 mm one were assigned regime 2)."""
    flow = FakeFlow(2)
    assert LIQ.mu / LIQ.rho < flow.mu_e[0] / flow.rho_e[0]       # aluminium is far the less viscous medium
    delta_m, _ = spray.melt_layer(flow, LIQ)
    b, v_s = np.full(2, 1e-3), np.full(2, 10.0)
    assert (b > delta_m).all()                                   # a genuinely deep film by Girin's own conjugate depth
    areas, m_f = np.full(2, 1e-5), b * LIQ.rho * 1e-5
    res = spray.SprayModel(LIQ).evaluate(flow, FakeState(), b, delta_m, v_s, np.full(2, True), 0.5, areas, m_f)
    assert res.branch.tolist() == [spray.BRANCH_THICK] * 2


# ---------------------------------------------------------------------------------------------------------------
# Step 3 amendment of 2026-10-06: the thin branch needs a rigid substrate

def test_a_thin_film_on_slurry_deeper_than_the_conjugate_depth_takes_the_thick_branch():
    """Girin's dominant ablation -- the thin branch, Girin & Kopyt's (1994) mode -- is the case where "the rigid core
    still stabilises the disturbances". Slurry (more than half liquid) is no rigid core, so the regime test reads the
    non-rigid layer, the film plus everything above 50 % liquid beneath it (`regime_layer`), against delta_m. A film on
    slurry deeper than delta_m is deep melt and takes Girin's (2017) thick branch with his Weber number on delta_m; one
    whose slurry ends within delta_m keeps the thin branch, its Weber number unchanged. Under Girin's closure the
    `on_slurry` flag holds nothing back: the depth decides."""
    flow = FakeFlow(2)
    model = spray.SprayModel(LIQ)
    delta_m, _ = spray.melt_layer(flow, LIQ)
    b = np.full(2, 2e-5)                                    # a film well inside the conjugate depth, nothing molten beneath
    regime = np.array([0.5, 3.0]) * delta_m                 # slurry ending within delta_m / slurry deeper than delta_m
    v_s = np.full(2, 20.0)
    areas, m_f = np.full(2, 1e-5), b * LIQ.rho * 1e-5
    res = model.evaluate(flow, FakeState(), b, delta_m, v_s, np.ones(2, bool), 0.5, areas, m_f, b_layer=b,
                         regime_layer=regime, on_slurry=np.array([True, True]))
    assert res.branch.tolist() == [spray.BRANCH_THIN, spray.BRANCH_THICK]
    assert res.we_s[1] == pytest.approx(LIQ.rho * 400.0 * delta_m[1] / LIQ.sigma)       # Girin's, on delta_m
    assert res.we_s[0] == pytest.approx(LIQ.rho * 400.0 * b[0] / LIQ.sigma)             # the thin branch's, on the film
    assert res.dm[0] > 0.0                                   # rigid material within the shear's reach: the thin mode acts
    old = model.evaluate(flow, FakeState(), b, delta_m, v_s, np.ones(2, bool), 0.5, areas, m_f, b_layer=b)
    assert old.branch.tolist() == [spray.BRANCH_THIN, spray.BRANCH_THIN]                 # the liquid layer alone: both thin


def test_the_rayleigh_taylor_criteria_keep_the_liquid_layer():
    """The slurry counts for the regime test only: the Rayleigh-Taylor depth criterion W h^2 rho_l > 3 Sigma describes a
    liquid pool, so it still reads the liquid layer (`b_layer`), not the non-rigid one."""
    flow = FakeFlow(1, decel=1e5)                            # deep enough a regime layer would pass the depth criterion
    model = spray.SprayModel(LIQ)
    delta_m, _ = spray.melt_layer(flow, LIQ)
    b, regime = np.array([2e-5]), 3.0 * delta_m
    assert spray.rayleigh_taylor(1e5, regime, LIQ)[0][0] and not spray.rayleigh_taylor(1e5, b, LIQ)[0][0]
    res = model.evaluate(flow, FakeState(), b, delta_m, np.array([20.0]), np.ones(1, bool), 0.5, np.full(1, 1e-5),
                         b * LIQ.rho * 1e-5, b_layer=b, regime_layer=regime)
    assert not res.rt_active and res.branch[0] == spray.BRANCH_THICK


def test_without_a_conjugate_depth_a_film_on_slurry_takes_no_shear_mode():
    """Off Girin's closure there is no conjugate depth and so no deep-melt mode, and the thin mode needs a rigid wall:
    a film whose base is slurry (`on_slurry`) is not sprayed by either, under the Couette closure and in the
    free-molecular branch alike. A film on a rigid wall keeps the thin (or rarefied) mode."""
    model = spray.SprayModel(LIQ)
    for branch, code in ((sf.BRANCH_SHOCK_LAYER, spray.BRANCH_THIN), (sf.BRANCH_FREE_MOLECULAR, spray.BRANCH_RAREFIED)):
        flow = FakeFlow(2, closure=sf.CLOSURE_COUETTE, branch=branch)
        delta_m, _ = spray.melt_layer(flow, LIQ)
        b = np.full(2, 2e-5)
        areas, m_f = np.full(2, 1e-5), b * LIQ.rho * 1e-5
        res = model.evaluate(flow, FakeState(), b, delta_m, np.full(2, 20.0), np.ones(2, bool), 0.5, areas, m_f,
                             on_slurry=np.array([True, False]))
        assert res.branch.tolist() == [-1, code]
        assert res.dm[0] == 0.0 and not res.unstable[0] and res.dm[1] > 0.0

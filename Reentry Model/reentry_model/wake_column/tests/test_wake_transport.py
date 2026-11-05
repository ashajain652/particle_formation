import math

import numpy as np
import pytest

from reentry_model.wake_column import transport as tr, species as sp
from reentry_model.wake_column.ambient import AmbientColumn
from reentry_model.wake_column.settings import Settings


def test_eddy_profile_is_megner_reference():
    e = tr.EddyProfile.load()
    Ks = [e.K(z) for z in np.arange(50, 100, 2.5)]
    assert 5.0 <= max(Ks) <= 15.0 and min(Ks) >= 0.05


def test_eddy_profile_interpolates_log_linearly_and_refuses_extrapolation():
    e = tr.EddyProfile(np.array([50.0, 60.0]), np.array([1.0, 100.0]))
    assert e.K(55.0) == pytest.approx(10.0, rel=1e-12)           # geometric mean at the midpoint
    assert e.K(50.0) == pytest.approx(1.0) and e.K(60.0) == pytest.approx(100.0)
    for z in (49.9, 60.1, float("nan")):
        with pytest.raises(ValueError, match="outside"):
            e.K(z)
    with pytest.raises(ValueError, match="outside"):
        tr.EddyProfile.load().K(5.0)


def test_molecular_diffusion_matches_spec_magnitudes():      # spec 4.1: ~0.01 m2/s at 50 km, ~3 at 90 km
    s, col = sp.load_species(), AmbientColumn.load("waccm", "night")
    assert 0.003 < tr.D_molecular(s["AlO"], col.at(50.0)) < 0.03
    assert 1.0 < tr.D_molecular(s["AlO"], col.at(90.0)) < 10.0


def test_molecular_diffusion_formula_on_a_hand_state():
    # D = 3/(8 n sigma12^2) sqrt(kT / (2 pi mu)), air 3.711 A / 28.96 amu, AlO 3.5 A / 42.981 amu
    s = sp.load_species()["AlO"]
    amb = AmbientColumn.load("waccm", "night").at(70.0)
    k, amu = 1.380649e-23, 1.66053906660e-27
    sig = 0.5 * (3.711 + s.sigma_A) * 1e-10
    mu = 28.96 * s.mass_amu / (28.96 + s.mass_amu) * amu
    want = 3.0 / (8.0 * amb.n_total_m3 * sig ** 2) * math.sqrt(k * amb.T_K / (2 * math.pi * mu))
    assert tr.D_molecular(s, amb) == pytest.approx(want, rel=1e-6)


def test_dissipation_table_covers_the_ambient_columns():
    e = tr.DissipationProfile.load()
    for z, want in ((50.0, 1e-3), (60.0, 2e-3), (70.0, 1e-3), (80.0, 1e-2), (90.0, 3e-2)):
        assert e.epsilon(z) == pytest.approx(want)                       # the review's central values, W/kg
    assert e.epsilon(75.0) == pytest.approx(math.sqrt(1e-3 * 1e-2))      # log-linear between rows
    assert e.epsilon(45.0) > 0 and e.epsilon(100.0) > 0
    with pytest.raises(ValueError, match="44.0 km"):
        e.epsilon(44.0)


def test_richardson_K_reproduces_the_t_cubed_law():
    """r^2 = g epsilon t^3 / 3 (a cloud's cross-section, from <l^2> = g epsilon t^3) solves d(r^2)/dt = 4 K_R(r)."""
    g, eps = 0.52, 1e-3
    for t in (10.0, 100.0, 1000.0):
        r = math.sqrt(g * eps * t ** 3 / 3.0)
        assert 4.0 * tr.richardson_K(r, eps, g) == pytest.approx(g * eps * t ** 2, rel=1e-12)


def test_total_diffusion_grows_with_size_up_to_the_capped_eddy_part():
    s, amb = sp.load_species()["AlO"], AmbientColumn.load("waccm", "night").at(70.0)
    K, eps = 5.0, 1e-3
    Dm = tr.D_molecular(s, amb)
    st = Settings()
    kr = lambda r, e=eps: tr.richardson_K(r, e, st.richardson_g)
    assert tr.D_total(1.0, s, amb, K, eps, st) == pytest.approx(Dm + kr(1.0))             # small plume: Richardson
    assert tr.D_total(1000.0, s, amb, K, eps, st) == pytest.approx(Dm + K)                # wide plume: capped
    assert tr.D_total(1.0, s, amb, K, eps, Settings(epsilon_scale=10.0)) == pytest.approx(Dm + kr(1.0, 10 * eps))
    assert tr.D_total(1000.0, s, amb, K, eps, Settings(eddy_scale=0.1)) == pytest.approx(Dm + 0.5)
    assert tr.D_total(1000.0, s, amb, K, eps, Settings(eddy_scale=0.0)) == pytest.approx(Dm)  # molecular only
    r = np.logspace(-3, 4, 400)
    D = np.array([tr.D_total(x, s, amb, K, eps, st) for x in r])
    assert np.all(np.diff(D) >= 0)                                                        # no step anywhere
    assert np.max(np.diff(np.log(D))) < 0.1


def test_volume_limits():
    assert tr.volume(1.0, 1e3) == pytest.approx(math.pi * 1e3, rel=1e-9)            # cylinder
    assert tr.volume(1e4, 1.0) == pytest.approx(math.pi ** 1.5 * 1e12, rel=1e-3)    # sphere
    assert tr.dlnV_dr2(1.0, 1e6) == pytest.approx(1.0, rel=1e-9)                     # d ln(pi r^2 L)/d r^2 = 1/r^2
    assert tr.dlnV_dr2(1.0, 1e-6) == pytest.approx(1.5, rel=1e-6)                    # d ln(r^3)/d r^2 = 1.5/r^2


def test_volume_finite_and_monotone_over_extreme_ratios():  # Review Focus 3
    r = np.logspace(-6, 6, 241)
    V = np.array([tr.volume(x, 1.0) for x in r]); g = np.array([tr.dlnV_dr2(x, 1.0) for x in r])
    assert np.all(np.isfinite(V)) and np.all(np.diff(V) > 0)
    assert np.all(np.isfinite(g)) and np.all((g * r**2 > 0.99) & (g * r**2 < 1.51))


def test_dlnV_dr2_matches_a_finite_difference_of_ln_volume():
    L = 7.0
    for r in (0.5, 2.0, 3.5, 6.0, 20.0):                    # x = L/(2r) from 8.75 down to 0.175
        h = 1e-5 * r * r
        fd = (math.log(tr.volume(math.sqrt(r * r + h), L)) - math.log(tr.volume(math.sqrt(r * r - h), L))) / (2 * h)
        assert tr.dlnV_dr2(r, L) == pytest.approx(fd, rel=1e-6)


def test_dlnV_dr2_continuous_across_the_series_switch():
    L = 1.0
    r_sw = L / (2 * 1e-4)                                   # L/(2r) = 1e-4
    lo, hi = tr.dlnV_dr2(r_sw * (1 + 1e-9), L), tr.dlnV_dr2(r_sw * (1 - 1e-9), L)
    assert lo == pytest.approx(hi, rel=1e-6)


def test_settling_rate():
    assert tr.settling_rate(0.03, 2.0) == pytest.approx(2 * 0.03 / (math.pi * 2.0))


def _smoke_rows():
    import csv
    import os
    from reentry_model.wake_column import DATA
    with open(os.path.join(DATA, "megner2006_smoke.csv"), newline="") as f:
        return [(float(r["z_km"]), float(r["r_nm"]), float(r["dN_dlnr_m3"]))
                for r in csv.DictReader(l for l in f if not l.startswith("#"))]


def test_smoke_table_shape():
    """Altitudes 45-100 km; at each, the visible Fig. 2b radii (0.3 nm and up) decrease steeply with radius, and
    the first row is the assumed deficit bin at Megner's smallest radius, 0.2 nm."""
    rows = _smoke_rows()
    zs = sorted({z for z, _, _ in rows})
    assert zs[0] == 45.0 and zs[-1] == 100.0
    for z in zs:
        r, d = np.array([(a, b) for zz, a, b in rows if zz == z]).T
        assert r[0] == pytest.approx(0.2) and np.all(np.diff(r) > 0) and np.all(d > 0)
        assert np.all(np.diff(d[1:]) < 0)


@pytest.mark.parametrize("z, n_fig2a_cm3", [(60.0, 1742.0), (75.0, 20287.0), (90.0, 47784.0)])
def test_smoke_table_integrates_to_fig2a_number_density(z, n_fig2a_cm3):
    """Sum of dN/dln r times the bin width in ln r (Megner's 35 geometric bins, ln(200)/34 each) over the table's
    radii equals the Fig. 2a total number density read at that altitude (cm^-3 recorded in the CSV header)."""
    r, d = np.array([(a, b) for zz, a, b in _smoke_rows() if zz == z]).T
    n_cm3 = d.sum() * (math.log(200.0) / 34) * 1e-6
    assert n_cm3 == pytest.approx(n_fig2a_cm3, rel=0.01)

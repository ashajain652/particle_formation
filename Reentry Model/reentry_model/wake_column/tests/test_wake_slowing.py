"""Stage A, the slowing-down Monte Carlo (spec section 4.4, tests 9-12, Review Focus 4)."""
import numpy as np
import pytest
from reentry_model.wake_column import slowing as sl, species as sp, reactions as rx
from reentry_model.wake_column.ambient import AmbientColumn
from reentry_model.wake_column.settings import Settings
from reentry_model.wake_column.trail import EmissionSample

S = sp.load_species()
CH = sl.calibrate(rx.load_reactions(stage=("hot", "both")), S)
NIGHT = AmbientColumn.load("waccm", "night")


def one(V, T, h=80.0):
    return (EmissionSample(0.0, h, V, T, 1.0),)


def test_calibrated_probabilities():
    P = {c.rid: c.P0 for c in CH}
    assert P["R1"] == pytest.approx(0.88, abs=0.03) and P["G42"] == pytest.approx(0.082, abs=0.01)
    assert P["D1"] == 1.0


def test_mean_energy_loss_per_collision_hard_sphere():        # spec test 9
    st = Settings(n_test_particles=20000)
    r = sl.slow_down(one(7600.0, 300.0), NIGHT, st, S, CH, partners_at_rest=True, reactions_on=False)
    assert r.per_sample[0]["mean_fractional_loss_vs_N2"] == pytest.approx(0.4998, abs=0.01)


def test_stopping_distance_scales_inversely_with_density():   # spec test 9
    st = Settings(n_test_particles=5000)
    a = sl.slow_down(one(7600.0, 1500.0, 70.0), NIGHT, st, S, CH, reactions_on=False)
    b = sl.slow_down(one(7600.0, 1500.0, 80.0), NIGHT, st, S, CH, reactions_on=False)
    ratio = NIGHT.at(70.0).n_total_m3 / NIGHT.at(80.0).n_total_m3
    assert b.r0_m / a.r0_m == pytest.approx(ratio, rel=0.15)


@pytest.mark.parametrize("rid", ["R1", "G42", "R3", "G33"])
def test_line_of_centers_reproduces_thermal_rate(rid):         # spec test 10
    ch = next(c for c in CH if c.rid == rid)
    T = next(x for x in rx.load_reactions(stage=("hot", "both")) if x.rid == rid).T_max
    rng = np.random.default_rng(1)
    assert sl.sample_rate(ch, T, S, 400000, rng) == pytest.approx(sl.loc_rate(ch, T, S), rel=0.05)


def test_cold_start_has_no_threshold_reactions():             # spec test 11
    T_amb = NIGHT.at(80.0).T_K
    r = sl.slow_down(one(0.0, T_amb), NIGHT, Settings(n_test_particles=5000), S, CH)
    assert r.fractions["OAlO"] == 0.0 and r.share_above_AlO_bond == 0.0
    lam = 1.0 / (NIGHT.at(80.0).n_total_m3 * np.pi * (3.711e-10) ** 2)
    assert 0.1 * lam < r.r0_m < 10 * lam


def test_packet_conservation_and_hot_products():               # spec test 12
    r = sl.slow_down(one(7600.0, 1800.0), NIGHT, Settings(n_test_particles=5000), S, CH)
    assert sum(r.fractions.values()) == pytest.approx(1.0, abs=1e-12)
    assert set(r.fractions) == {"Al", "AlO", "OAlO"}
    # The plan's bound here was < 1e-3, which ignored the isotropic thermal emission velocity of decision 15: at
    # T_d = 1800 K about 15 % of atoms leave faster than 8.6 km/s in the air frame, and about 1.7 % of Al collisions
    # exceed the 578 kJ/mol ionisation energy (an independent sampling of first collisions agrees). The droplet speed
    # alone (403 kJ/mol for Al + air at 7.6 km/s) does not reach it, which the cold-droplet case below checks.
    assert r.fractions["OAlO"] > 0.0 and r.share_above_ionisation < 0.05
    cold = sl.slow_down(one(7600.0, 300.0), NIGHT, Settings(n_test_particles=5000), S, CH)
    assert cold.share_above_ionisation < 1e-3


def test_partner_mix_follows_density_times_cross_section():    # pins the null-collision acceptance
    # With partners at rest and reactions off, partner j collides at n_j pi sigma_j^2 |v|, so its share of the
    # collisions is n_j sigma_j^2 / sum_k n_k sigma_k^2 (N2 standing for the remainder of the air). Accepting with the
    # whole majorant instead of the candidate's own term would weight each partner by its share twice (O2 about 5 %
    # instead of 19 %).
    r = sl.slow_down(one(7600.0, 300.0), NIGHT, Settings(n_test_particles=20000), S, CH, partners_at_rest=True,
                     reactions_on=False)
    a = NIGHT.at(80.0)
    n = {p: a.n[p] for p in ("O2", "O", "H2O")}
    n["N2"] = a.n_total_m3 - sum(n.values())
    w = {p: n[p] * (S["Al"].sigma_A + S[p].sigma_A) ** 2 for p in n}
    counts = r.per_sample[0]["collisions_by_partner"]
    total = sum(counts.values())
    for p in ("N2", "O2"):
        assert counts[p] / total == pytest.approx(w[p] / sum(w.values()), rel=0.03)


def test_degenerate_inputs_stay_finite():                      # Review Focus 4
    import dataclasses
    a = NIGHT.at(60.0)
    no_O = dataclasses.replace(a, n={**a.n, "O": 0.0})
    col = type("C", (), {"at": lambda self, z: no_O})()
    r = sl.slow_down(one(0.0, 908.0, 60.0), col, Settings(n_test_particles=2000), S, CH)
    assert np.isfinite(r.r0_m) and all(np.isfinite(v) for v in r.fractions.values())

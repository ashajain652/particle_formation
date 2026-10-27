"""The box model's particle terms (Task 8): formation by all 28 pairs, growth, oxygen uptake, uptake onto background
smoke, free-molecule coagulation, Epstein settling and the provenance tallies.

Spec section 7 tests 1 (aluminium conservation), 5 (self-preserving coagulation), 6 (settling decay) and 7
(provenance closure), plus the reporting of every pathway and route and a check of the analytic Jacobian.
"""
import numpy as np
import pytest

from reentry_model.wake_column.box import oxidation_g
from reentry_model.wake_column.tests.box_helpers import make_box


def test_aluminium_conserved_in_full_run():                     # spec test 1
    for tod in ("day", "night"):
        sol = make_box(time_of_day=tod, N_Al=1e14, t_end_s=3600.0).solve()
        for i in range(len(sol.t)):
            assert sum(sol.al_budget(i).values()) == pytest.approx(1.0, abs=1e-6)


def test_coagulation_reaches_self_preserving_width():          # spec test 5 (decision 8 fallback below)
    b = make_box(N_Al=0.0, chemistry=False, settling=False, smoke=False, freeze_radius=True,
                 initial_bins={3: 1e15}, r0=1.0, L=1.0, t_end_s=1e6)
    sol = b.solve(); bins = sol.bins(-1)
    assert bins["N"].sum() < 1e-3 * 1e15
    lnd = np.log(bins["d_m"]); w = bins["N"] / bins["N"].sum(); mu = (w * lnd).sum()
    gsd = np.exp(np.sqrt((w * (lnd - mu) ** 2).sum()))
    print(f"self-preserving geometric standard deviation {gsd:.4f} at bin ratio {b.bins.bin_ratio:.6g}")
    assert gsd == pytest.approx(1.45, abs=0.03)


def test_settling_only_decays_at_two_vs_over_pi_r():            # spec test 6
    b = make_box(N_Al=0.0, chemistry=False, coagulation=False, smoke=False, freeze_radius=True,
                 initial_bins={10: 1e6}, r0=2.0, t_end_s=1e3)
    sol = b.solve()
    rate = b.settling_rate_of_bin(10)
    for i in range(len(sol.t)):
        assert sol.bins(i)["N"][10] == pytest.approx(1e6 * np.exp(-rate * sol.t[i]), rel=1e-4)


def test_provenance_closes_in_every_bin():                      # spec test 7
    sol = make_box(time_of_day="day", N_Al=1e14, t_end_s=600.0).solve()
    for i in range(len(sol.t)):
        b = sol.bins(i)
        tallied = sum(b["tally"][s][r] for s in b["tally"] for r in b["tally"][s])
        assert np.allclose(tallied, b["N"] * b["n_atoms"], rtol=1e-6, atol=1e-6 * 1e14)
        assert sum(sol.provenance(i).values()) == pytest.approx(
            sum(v for k, v in sol.al_budget(i).items() if k in ("particles", "smoke", "settled")), abs=1e-9)


def test_every_pathway_and_the_smoke_route_are_reported():
    sol = make_box(time_of_day="night", N_Al=1e14, t_end_s=600.0).solve()
    prov = sol.provenance(-1)                                     # every (source, route) key is always present
    assert len(prov) == 7 * 3
    assert sum(v for (_, r), v in prov.items() if r == "formation") > 0
    assert sum(v for (_, r), v in prov.items() if r == "smoke") > 0


def test_seeded_bins_carry_one_oxygen_per_aluminium():
    b = make_box(N_Al=0.0, initial_bins={3: 1e10, 5: 2e9})
    u = b.layout.unpack(b.y0() * b.scale)
    n = b.bins.n_atoms
    assert u["O"][3] == pytest.approx(1e10 * n[3]) and u["O"][5] == pytest.approx(2e9 * n[5])
    assert not u["H"].any()


def test_analytic_jacobian_matches_finite_differences():
    """The analytic Jacobian (kernels frozen at the current composition) against a central difference of the
    right-hand side, at a state where the composition dependence is locally flat (H/Al above 3, so every particle is
    hydroxide-sized) and every process is active. The ambient-gas rows are left out: they carry no particle term
    (oxygen uptake does not consume the gas, decision 9), and their difference quotients are round-off against the
    large entrainment term."""
    b = make_box(time_of_day="day", z=76.0, N_Al=3e13, L=500.0, r0=0.5, start="AlO")
    lay, K = b.layout, b.bins.K
    rng = np.random.default_rng(1)
    x = b.y0() * b.scale
    x[lay.al] = 3e13 * rng.uniform(0.05, 0.3, lay.n_al)
    N = 1e8 * rng.uniform(0.5, 1.5, K) / b.bins.n_atoms
    x[lay.N] = N
    tally = rng.uniform(0.0, 1.0, (K, 7, 2))
    tally *= (N * b.bins.n_atoms / tally.sum(axis=(1, 2)))[:, None, None]
    x[lay.tally] = tally.ravel()
    al = N * b.bins.n_atoms
    x[lay.O] = al * np.where(np.arange(K) % 2 == 0, 1.2, 1.49)     # inside the flat part and the ramp of g
    x[lay.H] = 4.0 * al
    x[lay.r2] = 0.7 ** 2
    y = x / b.scale
    J = b.jac(0.0, y)
    J = J.toarray() if hasattr(J, "toarray") else np.asarray(J)
    Jfd = np.zeros_like(J)
    for j in range(lay.size):
        h = 1e-6 * max(abs(y[j]), 1e-6)
        yp, ym = y.copy(), y.copy()
        yp[j] += h; ym[j] -= h
        Jfd[:, j] = (b.rhs(0.0, yp) - b.rhs(0.0, ym)) / (2 * h)
    rows = np.setdiff1d(np.arange(lay.size), np.arange(lay.size)[lay.ambient])
    J, Jfd = J[rows], Jfd[rows]
    tol = 1e-4 * np.abs(Jfd).max(axis=1, keepdims=True) + 1e-12
    assert np.all(np.abs(J - Jfd) <= tol + 1e-4 * np.abs(Jfd))


def test_oxygen_uptake_factor_runs_from_one_to_zero_with_a_continuous_slope():
    x = np.array([1.0, 1.485, 1.4925, 1.5, 1.6])
    g, dg = oxidation_g(x)
    assert g.tolist() == pytest.approx([1.0, 1.0, 0.5, 0.0, 0.0])
    assert dg[[0, 1, 3, 4]].tolist() == [0.0, 0.0, 0.0, 0.0]
    eps = 1e-12                                                 # the slope tends to 0 at both ends of the ramp
    assert abs(oxidation_g(1.485 + eps)[1]) < 1e-5 and abs(oxidation_g(1.5 - eps)[1]) < 1e-5


def test_seeded_particles_take_up_oxygen_to_one_and_a_half_per_aluminium():
    b = make_box(N_Al=0.0, chemistry=False, coagulation=False, settling=False, smoke=False, freeze_radius=True,
                 initial_bins={5: 1e6}, t_end_s=1.0)
    sol = b.solve()
    assert sol.bins(0)["O_per_Al"][5] == pytest.approx(1.0)
    assert sol.bins(-1)["O_per_Al"][5] == pytest.approx(1.5, abs=1e-3)

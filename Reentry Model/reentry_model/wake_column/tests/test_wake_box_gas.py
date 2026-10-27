"""The box model's gas phase (Task 7): geometry, entrainment and the neutral network at ambient temperature.

Spec section 7 tests 2 (inert tracer), 3 (ambient stability) and 4 (single-reaction lifetimes; decision 12 puts the
total air density in R6's third body), and conservation of gas-phase aluminium with the full network by day and night.
"""
import numpy as np
import pytest

from reentry_model.wake_column import transport
from reentry_model.wake_column.tests.box_helpers import make_box, S


def test_inert_tracer_follows_r0_over_r_squared():            # spec test 2
    b = make_box(chemistry=False, particles=False, eddy_scale=0.0, t_end_s=1e3)
    sol = b.solve()
    D = transport.D_molecular(S["AlO"], b.amb)
    for i, t in enumerate(sol.t):
        r2 = 0.05**2 + 4 * D * t
        assert sol.r(i) ** 2 == pytest.approx(r2, rel=1e-6)
        assert sol.concentration("Al", i) * sol.V(i) == pytest.approx(1e12, rel=1e-6)


def test_ambient_stays_ambient_without_aluminium():             # spec test 3
    b = make_box(N_Al=0.0, t_end_s=1e3); sol = b.solve()
    for name in ("O", "O2", "O3", "H", "H2O", "CO2"):
        assert sol.concentration(name, -1) == pytest.approx(b.amb.n[name], rel=1e-9)


@pytest.mark.parametrize("z", [70.0, 80.0, 90.0])
def test_AlO_lifetime_against_R6(z):                            # spec test 4 (decision 12)
    b = make_box(z=z, start="AlO", reactions=("R6",), freeze_radius=True, particles=False, t_end_s=1e4)
    tau = 1.0 / (b.reaction("R6").k(b.amb.T_K) * b.amb.n_total_m3 * b.amb.n["O2"])
    print(f"AlO lifetime against R6 at {z:g} km (night WACCM, T = {b.amb.T_K:.1f} K): {tau:.4g} s")
    sol = b.solve()
    i = np.searchsorted(sol.t, tau)
    assert sol.al_budget(i)["AlO"] == pytest.approx(np.exp(-sol.t[i] / tau), rel=1e-3)


def test_Al_lifetime_against_R1():
    b = make_box(reactions=("R1",), freeze_radius=True, particles=False, t_end_s=1.0)
    tau = 1.0 / (b.reaction("R1").k(b.amb.T_K) * b.amb.n["O2"])
    sol = b.solve(); i = np.searchsorted(sol.t, tau)
    assert sol.al_budget(i)["Al"] == pytest.approx(np.exp(-sol.t[i] / tau), rel=1e-3)


def test_gas_aluminium_is_conserved_with_full_chemistry():
    for tod in ("day", "night"):
        sol = make_box(time_of_day=tod, particles=False, t_end_s=3600.0).solve()
        for i in range(len(sol.t)):
            assert sum(sol.al_budget(i).values()) == pytest.approx(1.0, abs=1e-6)

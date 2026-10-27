"""Checks of the droplet adapter (trail.py): the trail length, inventory, mean altitude and emission samples that
Stage A and the box model take from one droplet run."""
import pytest

from reentry_model.wake_column.trail import EmissionSample, Trail, trail_for, trail_from_rows

M_AL_ATOM = 26.9815385e-3 / 6.02214076e23


def _rows(ts, hs, vs, Ts, mAl, mZn=None, mMg=None):
    n = len(ts)
    mZn = mZn or [0.0] * n
    mMg = mMg or [0.0] * n
    return [dict(t_s=ts[i], h_km=hs[i], V_ms=vs[i], T_K=Ts[i], m_Al_kg=mAl[i], m_Zn_kg=mZn[i], m_Mg_kg=mMg[i])
            for i in range(n)]


@pytest.mark.parametrize("d, h, L", [(1.0, 80.0, 51.0), (10.0, 60.0, 40.0)])   # spec 4.2 table (2026-10-27 runs)
def test_trail_length_matches_design_measurements(parent, d, h, L):
    t = trail_for(d, h, parent=parent)
    assert t.L_m == pytest.approx(L, rel=0.05)
    assert abs(t.z_mean_km - h) < 0.05
    assert sum(s.weight for s in t.samples) == pytest.approx(1.0)


def test_inventory_is_the_whole_aluminium_loss(parent):
    from reentry_model.droplet_reentry.model import run_case
    fs, rows = run_case(10.0, 80.0, parent=parent)
    t = trail_from_rows(rows, 10.0, 80.0)
    lost = rows[0]["m_Al_kg"] - rows[-1]["m_Al_kg"]
    assert t.m_Al_kg == pytest.approx(lost, rel=1e-9)
    assert t.N_Al == pytest.approx(lost / (26.9815385e-3 / 6.02214076e23), rel=1e-9)
    assert all(100.0 < s.V_ms < 8000.0 and 900.0 < s.T_d_K < 2000.0 for s in t.samples)


def test_no_rows_gives_empty_trail():
    t = trail_from_rows(None, 1.0, 50.0)
    assert t.N_Al == 0 and t.samples == ()


def test_no_rows_keeps_release_altitude_and_zero_masses():
    t = trail_from_rows(None, 1.0, 50.0)
    assert (t.L_m, t.m_Al_kg, t.m_Zn_kg, t.m_Mg_kg, t.z_mean_km) == (0.0, 0.0, 0.0, 0.0, 50.0)


def test_rows_without_aluminium_loss_give_empty_trail_but_keep_zinc_and_magnesium():
    rows = _rows([0, 1, 2], [70, 69, 68], [5000, 4000, 3000], [900, 1000, 1100], [1e-15] * 3,
                 mZn=[3e-16, 2e-16, 1e-16], mMg=[2e-16, 2e-16, 0.5e-16])
    t = trail_from_rows(rows, 1.0, 70.0)
    assert t.N_Al == 0 and t.m_Al_kg == 0.0 and t.L_m == 0.0 and t.samples == ()
    assert t.z_mean_km == 70.0
    assert t.m_Zn_kg == pytest.approx(2e-16) and t.m_Mg_kg == pytest.approx(1.5e-16)


def test_synthetic_linear_loss_gives_exact_length_altitude_and_samples():
    # Constant speed 1000 m/s, aluminium lost linearly from 1e-15 to 0 over 10 s, altitude falling 80 -> 70 km.
    ts = [0.0, 2.0, 4.0, 6.0, 8.0, 10.0]
    rows = _rows(ts, [80 - t for t in ts], [1000.0] * 6, [1000.0 + 100 * t for t in ts],
                 [1e-15 * (1 - t / 10) for t in ts])
    t = trail_from_rows(rows, 1.0, 80.0, n_samples=5)
    assert t.L_m == pytest.approx(9000.0)                    # 90 % of the loss is reached at t = 9 s
    assert t.m_Al_kg == pytest.approx(1e-15)
    assert t.z_mean_km == pytest.approx(75.0)
    assert [s.t_s for s in t.samples] == pytest.approx([1.0, 3.0, 5.0, 7.0, 9.0])
    assert [s.T_d_K for s in t.samples] == pytest.approx([1100.0, 1300.0, 1500.0, 1700.0, 1900.0])
    assert all(s.weight == pytest.approx(0.2) for s in t.samples)
    assert isinstance(t.samples[0], EmissionSample) and isinstance(t, Trail)


def test_flat_stretches_of_the_loss_history_do_not_break_interpolation():
    # The aluminium is lost in the first step only; later rows (a long coast) carry no loss.
    rows = _rows([0, 1, 2, 3], [80, 79, 78, 77], [4000, 3000, 2000, 1000], [900, 1200, 1200, 1200],
                 [1e-15, 0.0, 0.0, 0.0])
    t = trail_from_rows(rows, 1.0, 80.0, n_samples=4)
    assert t.L_m == pytest.approx(0.9 * 3500.0)
    assert t.z_mean_km == pytest.approx(79.5)
    assert [s.V_ms for s in t.samples] == pytest.approx([3875.0, 3625.0, 3375.0, 3125.0])


def test_trail_for_returns_empty_trail_for_a_weber_skipped_case(parent):
    t = trail_for(5000.0, 60.0, parent=parent)   # release Weber number 31 > critical
    assert t.N_Al == 0 and t.samples == ()

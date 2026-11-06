"""The segmented trail: cutting a droplet's trail into equal lengths holding their real aluminium, adding the segments
at a common time since release, and the along-trail spreading applied after the runs."""
import csv
import math

import numpy as np
import pytest

from reentry_model.wake_column.tests.box_helpers import make_box
from reentry_model.wake_column import run
from reentry_model.wake_column.segments import Segment, TrailSolution, along_trail, cut_trail, state_at_age
from reentry_model.wake_column.settings import Case, Settings
from reentry_model.wake_column.transport import volume
from reentry_model.wake_column.trail import Trail


def rows_with_bump(n=2001, t_end=10.0, V=1000.0, m0=1e-9):
    """Droplet rows at constant speed whose aluminium leaves in a Gaussian bump around t = 3 s (x = 3 km)."""
    t = np.linspace(0.0, t_end, n)
    cum = 0.5 * (1.0 + np.vectorize(math.erf)((t - 3.0) / 0.5))
    cum = (cum - cum[0]) / (cum[-1] - cum[0])
    return [{"t_s": float(a), "h_km": 70.0 - 0.001 * float(a), "V_ms": V, "T_K": 1500.0 + float(a),
             "m_Al_kg": m0 * (1.0 - 0.5 * float(c))} for a, c in zip(t, cum)]


def test_cut_conserves_the_aluminium():
    segs = cut_trail(rows_with_bump(), 200.0, 70.0, 24, n_samples=4)
    assert sum(s.trail.m_Al_kg for s in segs) == pytest.approx(0.5e-9, rel=1e-12)


def test_cut_gives_equal_lengths_out_to_the_whole_trail():
    segs = cut_trail(rows_with_bump(), 200.0, 70.0, 10, n_samples=4)
    whole = segs[0].trail.L_m
    assert whole == pytest.approx(3000.0 + 500.0 * 2.185, rel=0.01)       # erf(z) = 0.998 at z = 2.185
    lengths = [s.x1_m - s.x0_m for s in segs]
    assert lengths == pytest.approx([whole / 10] * len(segs), rel=1e-9)
    assert sum(s.volume_share for s in segs) == pytest.approx(len(segs) / 10, rel=1e-9)


def test_mass_follows_the_release_not_equal_shares():
    segs = cut_trail(rows_with_bump(), 200.0, 70.0, 20, n_samples=4)
    shares = np.array([s.trail.m_Al_kg for s in segs]) / 0.5e-9
    peak = segs[int(np.argmax(shares))]
    assert peak.x0_m <= 3000.0 <= peak.x1_m
    assert shares.max() > 3.0 / len(segs)


def test_each_segment_carries_its_own_flight_state():
    segs = cut_trail(rows_with_bump(), 200.0, 70.0, 10, n_samples=4)
    laid = [s.t_laid_s for s in segs]
    assert laid == sorted(laid)
    for s in segs:
        assert len(s.trail.samples) == 4
        assert sum(q.weight for q in s.trail.samples) == pytest.approx(1.0)
        assert s.t_laid_s == pytest.approx(np.mean([q.t_s for q in s.trail.samples]))
        assert s.trail.z_mean_km == pytest.approx(70.0 - 0.001 * s.t_laid_s, abs=1e-6)
    late = [s for s in segs if s.x0_m >= 3000.0]
    assert late[0].trail.samples[0].T_d_K > 1500.0 + 2.9


def test_a_stretch_with_no_release_makes_no_segment():
    rows = rows_with_bump()
    segs = cut_trail(rows, 200.0, 70.0, 200, n_samples=2)
    assert all(s.trail.N_Al > 0.0 for s in segs)
    assert len(segs) < 200


def test_box_volume_share_scales_the_volume():
    a = make_box(L=100.0, r0=0.05, t_end_s=1.0, freeze_radius=True, reactions=(), particles=False)
    b = make_box(L=100.0, r0=0.05, t_end_s=1.0, freeze_radius=True, reactions=(), particles=False, volume_share=0.25)
    assert b.V0 == pytest.approx(0.25 * a.V0)
    assert b.solve().V(0) == pytest.approx(0.25 * volume(0.05, 100.0))


def test_a_quarter_of_the_atoms_in_a_quarter_of_the_volume_evolves_the_same():
    """Same concentrations, so the same shares, including the second-order formation that depends on them."""
    kw = dict(z=70.0, L=1000.0, r0=0.01, start="AlO", t_end_s=1.0, time_of_day="night")
    a = make_box(N_Al=1e16, **kw).solve()
    b = make_box(N_Al=2.5e15, volume_share=0.25, **kw).solve()
    for key, v in a.al_budget(len(a.t) - 1).items():
        assert b.al_budget(len(b.t) - 1)[key] == pytest.approx(v, rel=1e-3, abs=1e-9)
    assert a.al_budget(len(a.t) - 1)["particles"] > 1e-4


def two_part_solution(t=None):
    kw = dict(z=70.0, L=1000.0, r0=0.01, start="AlO", t_end_s=100.0, time_of_day="night")
    parts = []
    for i, (N, laid) in enumerate(((3e15, 0.1), (1e15, 0.3))):
        tr = Trail(1.0, 70.0, 1000.0, N, 0.0, 0.0, 0.0, 70.0, ())
        seg = Segment(index=i, x0_m=500.0 * i, x1_m=500.0 * (i + 1), t_laid_s=laid, trail=tr, volume_share=0.5)
        parts.append((seg, make_box(N_Al=N, volume_share=0.5, **kw).solve()))
    grid = parts[0][1].t if t is None else t
    return parts, TrailSolution(parts, grid)


def test_aluminium_is_conserved_at_every_time_since_release():
    _, sol = two_part_solution()
    for i in range(len(sol.t)):
        b = sol.al_budget(i)
        total = sum(b.values()) + sol.unreleased(i)
        assert total == pytest.approx(1.0, abs=1e-5)


def test_a_segment_counts_only_after_it_is_laid():
    _, sol = two_part_solution(np.array([0.05, 0.2, 0.5]))
    assert sol.unreleased(0) == pytest.approx(1.0)
    assert sol.unreleased(1) == pytest.approx(0.25)
    assert sol.unreleased(2) == pytest.approx(0.0)


def test_the_sum_reads_each_segment_at_its_own_age():
    parts, sol = two_part_solution(np.array([50.1]))
    seg, part = parts[0]
    k = int(np.argmin(np.abs(part.t - 50.0)))
    _, near = two_part_solution(np.array([part.t[k] + 0.1]))
    other = parts[1][1]
    expect = part.count("AlO", k) + float(np.interp(np.log(part.t[k] - 0.2), np.log(other.t[1:]),
                                                    [other.count("AlO", i) for i in range(1, len(other.t))]))
    assert near.count("AlO", 0) == pytest.approx(expect, rel=1e-3)


def alo_at(part, age):
    """AlO molecules of one segment's Solution at `age`, read as the sum over segments reads them."""
    return float(state_at_age(part, np.array([age]))[0][part.box.layout.index("AlO")])


def test_along_trail_keeps_every_atom():
    parts, _ = two_part_solution()
    x = np.linspace(-2e3, 3e3, 2001)
    p = along_trail(parts, 50.0, x)
    assert np.trapezoid(p["line"]["AlO"], x) == pytest.approx(
        sum(alo_at(q, 50.0 - s.t_laid_s) for s, q in parts), rel=1e-6)


def test_along_trail_is_a_top_hat_when_the_plume_is_narrow():
    parts, _ = two_part_solution()
    x = np.array([250.0, 750.0])
    p = along_trail(parts, 0.5, x)
    a, b = parts[0], parts[1]
    assert p["line"]["AlO"][0] / p["line"]["AlO"][1] == pytest.approx(alo_at(a[1], 0.4) / alo_at(b[1], 0.2), rel=1e-6)


def test_axis_concentration_of_one_segment_matches_its_box():
    kw = dict(z=70.0, L=500.0, r0=0.01, start="AlO", t_end_s=100.0, time_of_day="night")
    sol = make_box(N_Al=1e15, **kw).solve()
    tr = Trail(1.0, 70.0, 500.0, 1e15, 0.0, 0.0, 0.0, 70.0, ())
    seg = Segment(index=0, x0_m=0.0, x1_m=500.0, t_laid_s=0.0, trail=tr, volume_share=1.0)
    i = len(sol.t) - 1
    p = along_trail([(seg, sol)], sol.t[i], np.array([250.0]))
    assert p["axis"]["AlO"][0] == pytest.approx(sol.concentration("AlO", i), rel=1e-9)


def test_spreading_matches_a_top_hat_convolved_with_the_radius():
    """No chemistry and a fixed 50 m radius, so the exact answer is the closed form for a 100 m top hat."""
    tr = Trail(1.0, 70.0, 100.0, 1e15, 0.0, 0.0, 0.0, 70.0, ())
    sol = make_box(N_Al=1e15, z=70.0, L=100.0, r0=50.0, start="AlO", t_end_s=10.0, time_of_day="night",
                   freeze_radius=True, chemistry=False, particles=False).solve()
    seg = Segment(index=0, x0_m=0.0, x1_m=100.0, t_laid_s=0.0, trail=tr, volume_share=1.0)
    p = along_trail([(seg, sol)], 5.0, np.array([50.0, 150.0]))
    assert p["line"]["AlO"][0] == pytest.approx(1e13 * math.erf(1.0), rel=1e-6)
    assert p["line"]["AlO"][1] == pytest.approx(1e13 * 0.5 * (math.erf(3.0) - math.erf(1.0)), rel=1e-6)


def test_run_writes_the_segmented_outputs(tmp_path, parent):
    st = Settings(t_end_s=60.0, n_test_particles=400, n_emission_samples=2, n_segments=3)
    s = run.run(Case(1.0, 80.0, "night"), st, out_dir=tmp_path)
    d = tmp_path / s["run_name"]
    for name in ("history.csv", "size_distribution.csv", "segments.csv", "along_trail.csv", "summary.json"):
        assert (d / name).exists()
    assert s["n_segments"] == 3 and s["trail_length_m"] > s["L_m"]
    with open(d / "segments.csv", newline="") as f:
        table = list(csv.DictReader(f))
    assert 1 <= len(table) <= 3
    assert sum(float(r["share_Al"]) for r in table) == pytest.approx(1.0, rel=1e-9)
    assert sum(float(r["share_Al"]) * float(r["condensed"]) for r in table) == pytest.approx(
        s["condensed_fraction_end"], rel=1e-3)
    with open(d / "history.csv", newline="") as f:
        hist = list(csv.DictReader(f))
    gas = ["frac_" + n for n in ("Al", "AlO", "OAlO", "OAlO2", "AlCO3", "AlOH", "AlOH2")]
    for row in hist:
        total = sum(float(row[c]) for c in gas) + float(row["frac_particles"]) + float(row["frac_smoke"]) \
            + float(row["frac_settled"]) + float(row["frac_unreleased"])
        assert total == pytest.approx(1.0, abs=1e-4)
    assert float(hist[0]["frac_unreleased"]) == pytest.approx(1.0)
    assert float(hist[-1]["frac_unreleased"]) == 0.0
    with open(d / "along_trail.csv", newline="") as f:
        prof = list(csv.DictReader(f))
    for c in ("t_s", "x_m", "line_AlOH_per_m", "c_AlOH_m3", "line_particles_per_m"):
        assert c in prof[0]


def test_settings_name_a_non_default_segment_count():
    assert Settings().n_segments == 24
    assert "seg48" in Settings(n_segments=48).tag()
    assert "seg0" in Settings(n_segments=0).tag()
    with pytest.raises(ValueError):
        Settings(n_segments=-1)

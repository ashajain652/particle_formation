"""Tests for the overlap engine (overlap.py): the union-of-balls margin against a brute-force reference, every inside
interval, censoring, the ahead-of-parent rule, the window check and the per-interval validity record.

The helpers build turbulent-flow tracks (Re_d about 1e7, pressure parameter about 53 cm Hg-cm), so the radius is
Taylor's unless a test says otherwise. Taylor's radius behind this parent is about 10.6 t^(1/3) m at t seconds behind
a 6 km/s parent, which sets the lateral speeds the hand-computed cases use.
"""
import ast, inspect, math, numpy as np, pytest
from reentry_model.wake_residence import overlap
from reentry_model.wake_residence.overlap import Plume, residence
from reentry_model.wake_residence.track import Track
from reentry_model.wake_residence.laws import LawSettings, taylor_radius

U = 6000.0
TURB = dict(rho=1e-2, T=250.0, p=300.0, mu=1.6e-5, mfp=1e-5, d=2.813, cda=5.5, r0=1.4065)

def line_parent(T=20.0, n=2001, flow=TURB):
    t = np.linspace(0, T, n); r = np.c_[U*t, 0*t, 0*t]
    return Track.from_arrays(t, r, v=np.tile([U, 0, 0], (n, 1)), flow={k: np.full(n, v) for k, v in flow.items()})

def follower(t, r):
    return Track.from_arrays(np.asarray(t), np.asarray(r))

@pytest.mark.parametrize("modname", ["overlap", "laws", "track", "wake_gas"])
def test_engine_imports_nothing_physical(modname):
    # wake_gas.py arrives in Task 4; until then its case is skipped rather than failed
    mod = pytest.importorskip("reentry_model.wake_residence." + modname)
    tree = ast.parse(inspect.getsource(mod))
    names = {a.name for n in ast.walk(tree) if isinstance(n, (ast.Import, ast.ImportFrom))
             for a in n.names} | {n.module or "" for n in ast.walk(tree) if isinstance(n, ast.ImportFrom)}
    for bad in ("trajectory", "droplet_reentry", "parent_trajectory", "body", "aero", "atmosphere"):
        assert not any(bad in x.split(".") for x in names), (modname, bad)

def test_regime_functions_take_no_altitude():
    from reentry_model.wake_residence import laws
    for fn in (laws.classify, laws.plume_radius, overlap.Plume.margin):
        assert not any("alt" in p or p in ("h", "h_km") for p in inspect.signature(fn).parameters)

def test_fast_margin_matches_bruteforce():
    p = Plume(line_parent()); rng = np.random.default_rng(1)
    for _ in range(50):
        t = rng.uniform(2, 20); x = np.array([U*t - rng.uniform(0, 5000), rng.uniform(-30, 30), rng.uniform(-30, 30)])
        assert p.margin(x, t).margin == pytest.approx(p.margin_bruteforce(x, t).margin, abs=1e-6)

def test_straight_divergence_matches_hand_computation():
    # follower stays at the release point (x = 0) while the parent moves on; lateral offset grows at 5 m/s
    # (at 1 m/s Taylor's radius, about 10.6 t^(1/3) m, stays larger than the offset until about 34 s, past the track)
    t = np.linspace(0, 20, 401); f = follower(t, np.c_[0*t, 5.0*t, 0*t])
    res = residence(Plume(line_parent()), f, t_start=0.0)
    from scipy.optimize import brentq
    # the best emission point is the release point itself (tau = 0): radius at s = U t, distance 5 t
    expected = brentq(lambda tt: taylor_radius(U*tt, 5.5, 1.4065, 0.66, 1/3) - 5.0*tt, 0.01, 19.9)
    assert res.residence_time == pytest.approx(expected, abs=1e-3)

def test_curved_parent_matches_bruteforce_exit():
    R = 6.4e6; w = U / R; t = np.linspace(0, 20, 4001)
    r = np.c_[R*np.cos(w*t), R*np.sin(w*t), 0*t]; v = np.c_[-U*np.sin(w*t), U*np.cos(w*t), 0*t]
    par = Track.from_arrays(t, r, v=v, flow={k: np.full(t.size, x) for k, x in TURB.items()})
    f = follower(t, np.c_[R*np.ones_like(t), 0*t, 2.0*t])           # stays at the start point, rises 2 m/s
    p = Plume(par); res = residence(p, f, t_start=0.0)
    te = res.intervals[0][1]
    assert p.margin_bruteforce(f.position(te), te).margin == pytest.approx(0.0, abs=1e-3)

def test_follower_starting_outside_enters_later_and_reenters():
    t = np.linspace(0, 20, 2001)
    y = 22 + 10*np.cos(0.8*t)        # 32 m at t=0 (outside), 12 m near t=4 (inside), 32 m near t=8 (outside), 12 m near t=12
    res = residence(Plume(line_parent()), follower(t, np.c_[0*t, y, 0*t]), t_start=0.0)
    assert res.residence_time == 0.0 and res.reason == "starts_outside"
    assert len(res.intervals) >= 2

def test_r0_changing_with_emission_time():
    # lateral offset grows at 5 m/s so that both cases exit before the track ends (at 1 m/s both are censored at 20 s)
    t = np.linspace(0, 20, 401); f = follower(t, np.c_[0*t, 5.0*t, 0*t])
    park = LawSettings(turbulent_law="park")                         # Park's radius scales with r0^(2/3)
    fixed = residence(Plume(line_parent(), park, r0=1.4065), f, t_start=0.0).residence_time
    shrinking = residence(Plume(line_parent(), park, r0=lambda tau: 1.4065 * 0.5), f, t_start=0.0).residence_time
    assert shrinking < fixed

def test_subsampled_follower_gives_the_same_exit_and_a_larger_error_estimate():
    # the 0.2 sin(3t) term keeps the path from being a quadratic, which a cubic spline reproduces exactly (both error
    # estimates would then be round-off, about 1e-15 m, and their order a coin toss)
    t = np.linspace(0, 20, 801); r = np.c_[0*t, 0.05*t**2 + 0.5*t + 0.2*np.sin(3*t), 0*t]
    fine = residence(Plume(line_parent()), follower(t, r), t_start=0.0)
    coarse = residence(Plume(line_parent()), follower(t[::2], r[::2]), t_start=0.0)
    assert abs(fine.residence_time - coarse.residence_time) < 1e-3
    assert coarse.interp_error_m >= fine.interp_error_m

def test_follower_identical_to_parent_is_censored():
    par = line_parent(); res = residence(Plume(par), par, t_start=1.0, t_end=19.0)
    assert res.censored and res.intervals[-1][1] is None and res.reason == "censored"

def test_follower_ahead_of_parent_has_small_positive_residence():           # Review Focus 4
    t = np.linspace(0, 20, 2001); f = follower(t, np.c_[1.05*U*t, 0*t, 0*t])
    res = residence(Plume(line_parent()), f, t_start=0.0)
    assert 0.0 < res.residence_time < 0.05 and res.reason == "ok"

def test_window_outside_parent_raises():                                    # Review Focus 2
    t = np.linspace(-5, 20, 101); f = follower(t, np.c_[0*t, 0*t, 0*t])
    with pytest.raises(ValueError, match="parent"):
        residence(Plume(line_parent()), f, t_start=-5.0)

def test_law_switching_along_one_track_tags_intervals():
    n = 3001; t = np.linspace(0, 30, n)
    flow = {k: np.full(n, v) for k, v in TURB.items()}
    flow["mfp"] = np.where(t < 10, 0.1, 1e-5)                         # rarefied, then continuum
    flow["rho"] = np.where(t < 20, 1e-5, 1e-2); flow["p"] = np.where(t < 20, 0.5, 300.0)   # laminar, then turbulent
    par = Track.from_arrays(t, np.c_[U*t, 0*t, 0*t], v=np.tile([U, 0, 0], (n, 1)), flow=flow)
    p = Plume(par)
    regimes = [p.margin(par.position(tt) + [0, 0.5, 0], tt).regime for tt in (5.0, 15.0, 25.0)]
    assert regimes == ["rarefied", "laminar", "turbulent"]

def test_rarefied_start_reports_no_law():
    flow = dict(TURB, mfp=0.1)
    t = np.linspace(0, 20, 401)
    res = residence(Plume(line_parent(flow=flow)), follower(t, np.c_[0*t, 0.1*t, 0*t]), t_start=0.0)
    assert math.isnan(res.residence_time) and res.reason == "rarefied_no_law"

def test_validity_record_is_populated():
    t = np.linspace(0, 20, 401)
    res = residence(Plume(line_parent()), follower(t, np.c_[0*t, 1.0*t, 0*t]), t_start=0.0)
    v = res.validity[0]
    assert v.regimes == {"turbulent"} and v.kn and v.re_d and v.p_param and "taylor" in v.laws


# ---- checks beyond the brief: the adaptive scan, both ends of the window, and the too-slow parent ----

def test_adaptive_scan_finds_the_same_crossings_as_a_uniform_grid():
    t = np.linspace(0, 20, 2001); f = follower(t, np.c_[0*t, 22 + 10*np.cos(0.8*t), 0*t])
    p = Plume(line_parent())
    adaptive = residence(p, f, t_start=0.0, t_end=11.0)
    uniform = residence(p, f, t_start=0.0, t_end=11.0, scan_dt=0.01, max_step=0.01)   # the brief's fixed grid
    assert len(adaptive.intervals) == len(uniform.intervals) == 2                # in, out, in again (censored)
    for (a0, a1), (u0, u1) in zip(adaptive.intervals, uniform.intervals):
        assert a0 == pytest.approx(u0, abs=1e-3) and (a1 is None) == (u1 is None)
        if a1 is not None:
            assert a1 == pytest.approx(u1, abs=1e-3)

def test_window_ending_after_parent_raises_and_a_later_start_is_normal():
    t = np.linspace(0, 25, 101); f = follower(t, np.c_[0*t, 0*t, 0*t])
    with pytest.raises(ValueError, match=r"\[0\.0, 25\.0\].*parent.*\[0\.0, 20\.0\]"):
        residence(Plume(line_parent()), f)
    t = np.linspace(5, 20, 101); late = residence(Plume(line_parent()), follower(t, np.c_[U*t, 0*t, 0*t]))
    assert late.reason == "censored" and late.intervals[0][0] == 5.0

def test_too_slow_parent_reports_no_law():
    n = 201; t = np.linspace(0, 20, n)
    par = Track.from_arrays(t, np.c_[0.5*t, 0*t, 0*t], v=np.tile([0.5, 0, 0], (n, 1)),
                            flow={k: np.full(n, v) for k, v in TURB.items()})       # 0.5 m/s, below u_min_parent
    res = residence(Plume(par), follower(t, np.c_[0.5*t, 0*t, 0*t]), t_start=0.0)
    assert math.isnan(res.residence_time) and res.reason == "too_slow_no_law" and res.intervals == []

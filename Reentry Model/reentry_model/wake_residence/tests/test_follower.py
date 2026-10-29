"""Followers in the wake (Task 6, spec 4.3 and 4.3a): the WakeAtmosphere wrapper, sphere followers and sphere parents."""
import time, numpy as np, pytest
from reentry_model.atmosphere import US76TableAtmosphere
from reentry_model.wake_residence.adapters import (WakeAtmosphere, sphere_simulator, dense_solve, release_state)
from reentry_model.wake_residence.overlap import Plume
from reentry_model.wake_residence.laws import LawSettings

def release(twin, h_km=60.0):
    t_r = twin.t_at_altitude(h_km); return t_r, release_state(twin.track, t_r, twin.epoch)

def test_off_returns_the_wrapped_freestream(twin):
    base = US76TableAtmosphere(); w = WakeAtmosphere(base, Plume(twin.track), 0.0, "off")
    a, b = w.state(1.0, 70e3, 0.5, 0.2), base.state(1.0, 70e3, 0.5, 0.2)
    assert (a.rho, a.T, a.p, a.m_bar) == (b.rho, b.T, b.p, b.m_bar) and np.array_equal(a.wind_enu, b.wind_enu)

def test_off_run_reproduces_plain_simulator(twin):
    t_r, st = release(twin)
    plain = dense_solve(sphere_simulator(st, 0.05, US76TableAtmosphere()), 5.0)
    wrapped = dense_solve(sphere_simulator(st, 0.05, WakeAtmosphere(US76TableAtmosphere(), Plume(twin.track), t_r, "off")), 5.0)
    assert np.allclose(plain.sol.sol(5.0), wrapped.sol.sol(5.0), rtol=0, atol=1e-9)

def test_tailwind_and_thin_gas_reduce_deceleration(twin):
    t_r, st = release(twin)
    speeds = {}
    for mode in ("off", "both"):
        run = dense_solve(sphere_simulator(st, 0.05, WakeAtmosphere(US76TableAtmosphere(), Plume(twin.track), t_r, mode)), 1.0)
        speeds[mode] = np.linalg.norm(run.sol.sol(1.0)[3:])
    assert speeds["both"] > speeds["off"]

def test_timing_of_one_wake_run_is_recorded(twin, capsys):
    t_r, st = release(twin)
    w = WakeAtmosphere(US76TableAtmosphere(), Plume(twin.track), t_r, "both")
    t0 = time.perf_counter(); dense_solve(sphere_simulator(st, 0.05, w), 2.0)
    print("wake run: {:.1f} s for 2 s of flight, {} atmosphere calls".format(time.perf_counter() - t0, w.n_calls))
    assert w.n_calls > 0


# ---- beyond the brief's four tests: the conversions and the tracks the later tasks rely on ----

def test_wake_gas_velocity_is_added_to_the_wind_in_east_north_up(twin):
    from reentry_model import earth
    from reentry_model.atmosphere import StaticProfileWind
    from reentry_model.wake_residence.wake_gas import Ambient
    t_r = twin.t_at_altitude(70.0)
    base = US76TableAtmosphere(wind=StaticProfileWind())
    plume = Plume(twin.track)
    w = WakeAtmosphere(base, plume, t_r, "both")
    dt = 0.05                                            # 50 ms after release, 300 m behind, 0.5 m to the side
    side = np.cross(twin.track.velocity(t_r + dt), twin.track.position(t_r + dt))
    x = twin.track.position(t_r + dt) - 300.0 * twin.track.velocity(t_r + dt) / twin.track.speed(t_r + dt) \
        + 0.5 * side / np.linalg.norm(side)
    h, lat, lon = earth.ecef_to_geodetic(x)
    fs, b = w.state(dt, h, lat, lon), base.state(dt, h, lat, lon)
    gas, m = plume.wake_state(earth.geodetic_to_ecef(h, lat, lon), t_r + dt, Ambient(b.p, b.T, b.rho), "both",
                              t_release=t_r)
    assert m.margin > 0.0 and np.linalg.norm(gas.u) > 1.0 and np.any(b.wind_enu != 0.0)
    e, n, u = earth.enu_basis(lat, lon)
    added = fs.wind_enu - b.wind_enu
    assert np.allclose(added[0] * e + added[1] * n + added[2] * u, gas.u, rtol=1e-12, atol=1e-9)
    assert (fs.rho, fs.T, fs.p, fs.m_bar) == (gas.rho, gas.T, gas.p, b.m_bar) and fs.rho < b.rho


def test_far_from_the_plume_the_wrapped_freestream_is_returned(twin):
    t_r = twin.t_at_altitude(70.0)
    base = US76TableAtmosphere()
    w = WakeAtmosphere(base, Plume(twin.track), t_r, "both")
    a = w.state(0.0, 70e3, 0.0, 0.0)                     # a point thousands of km from the track
    b = base.state(0.0, 70e3, 0.0, 0.0)
    assert (a.rho, a.T, a.p, a.m_bar) == (b.rho, b.T, b.p, b.m_bar) and np.array_equal(a.wind_enu, b.wind_enu)
    assert w.last_margin.margin < 0.0 and w.n_calls == 1


def test_bad_mode_and_time_beyond_the_parent_raise(twin):
    with pytest.raises(ValueError, match="wake_flow"):
        WakeAtmosphere(US76TableAtmosphere(), Plume(twin.track), 0.0, "on")
    t_r = twin.t_at_altitude(60.0)
    w = WakeAtmosphere(US76TableAtmosphere(), Plume(twin.track), t_r, "both")
    assert w.follower_t_max == pytest.approx(twin.track.t1 - t_r)
    with pytest.raises(ValueError, match="outside the parent track's range"):
        w.state(w.follower_t_max + 1.0, 50e3, 0.0, 0.0)


def test_solution_track_is_the_run_on_the_parent_clock(twin):
    from reentry_model.wake_residence.adapters import solution_track
    t_r, st = release(twin)
    run = dense_solve(sphere_simulator(st, 0.05, US76TableAtmosphere()), 1.0)
    tr = solution_track(run, t_r)
    assert tr.t0 == pytest.approx(t_r) and tr.t1 == pytest.approx(t_r + run.t_end)
    for tt in (0.0, 0.3337, 1.0):
        assert np.allclose(tr.position(t_r + tt), run.sol.sol(tt)[:3], rtol=0, atol=1e-6)
    assert np.allclose(tr.position(t_r), twin.track.position(t_r), rtol=0, atol=1e-6)   # released on the parent
    assert len(tr._t) == 501                             # every 2 ms over 1 s, both ends included


def test_sphere_parent_track_carries_the_spheres_own_drag(twin):
    from reentry_model import aero, earth
    from reentry_model.wake_residence.adapters import sphere_parent_track
    t_r, st = release(twin)
    atm = US76TableAtmosphere()
    sim = sphere_simulator(st, 0.05, atm)
    run = dense_solve(sim, 1.0)
    tr = sphere_parent_track(run, sim, t_r, atm, cadence=0.01)
    e = tr.emission_state(t_r + 0.5)
    a = sim.aero_state(0.5, run.sol.sol(0.5)[:3], run.sol.sol(0.5)[3:])
    assert e.d == pytest.approx(0.05) and e.r0 == pytest.approx(0.025)
    assert e.cda == pytest.approx(a.cd * sim.area, rel=1e-6) and e.U == pytest.approx(a.V, rel=1e-6)
    assert e.rho == pytest.approx(a.freestream.rho, rel=1e-6)
    plume = Plume(tr)                                     # a sphere parent's plume can be built and queried
    assert plume.margin(tr.position(t_r + 0.5), t_r + 0.5).margin == pytest.approx(0.025, rel=1e-6)

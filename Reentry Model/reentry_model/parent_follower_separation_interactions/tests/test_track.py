import numpy as np, pytest
from reentry_model.parent_follower_separation_interactions.track import Track

FLOW_KEYS = ("rho", "T", "p", "mu", "mfp", "d", "cda", "r0")


def circle(n):
    t = np.linspace(0, 10, n); R, w = 6.4e6, 7500 / 6.4e6
    r = np.c_[R*np.cos(w*t), R*np.sin(w*t), 0*t]; v = np.c_[-R*w*np.sin(w*t), R*w*np.cos(w*t), 0*t]
    return t, r, v


def test_straight_line_is_exact_and_path_length_is_speed_times_time():
    t = np.linspace(0, 2, 5); r = np.c_[7000*t, 0*t, 0*t]
    tr = Track.from_arrays(t, r, v=np.tile([7000, 0, 0], (5, 1)))
    assert np.allclose(tr.position(1.234), [8638.0, 0, 0], atol=1e-9)
    assert tr.path_length(1.5) == pytest.approx(10500.0, rel=1e-9)


def test_hermite_on_orbit_arc_is_millimetre_accurate():
    t, r, v = circle(101)
    tr = Track.from_arrays(t, r, v=v)
    tt = 3.21; R, w = 6.4e6, 7500 / 6.4e6
    assert np.linalg.norm(tr.position(tt) - [R*np.cos(w*tt), R*np.sin(w*tt), 0]) < 1e-3


def test_interpolation_error_flags_coarse_sampling():
    t, r, v = circle(6)
    assert Track.from_arrays(t, r).interpolation_error() > Track.from_arrays(*circle(201)[:2]).interpolation_error()


@pytest.mark.parametrize("bad", ["dup", "nan"])
def test_bad_input_names_index(bad):
    t = np.array([0.0, 1.0, 2.0, 3.0]); r = np.zeros((4, 3)); r[:, 0] = t
    if bad == "dup": t[2] = 1.0
    else: r[3, 1] = np.nan
    with pytest.raises(ValueError, match="index"):
        Track.from_arrays(t, r)


def test_emission_state_interpolates_flow():
    t = np.array([0.0, 1.0]); r = np.c_[[0.0, 1.0], [0, 0], [0, 0]]
    flow = {k: np.array([1.0, 3.0]) for k in ("rho", "T", "p", "mu", "mfp", "d", "cda", "r0")}
    s = Track.from_arrays(t, r, flow=flow).emission_state(0.5)
    assert s.rho == pytest.approx(2.0) and s.U == pytest.approx(1.0)


# ---- validation: every array, the first bad index and the array name in the message ----

def _good(n=5):
    t = np.linspace(0.0, 4.0, n)
    r = np.c_[100.0*t, 0*t, 0*t]
    return t, r, np.tile([100.0, 0, 0], (n, 1)), np.full(n, 100.0), {k: np.ones(n) for k in FLOW_KEYS}


def test_decreasing_time_names_array_and_first_index():
    t, r, *_ = _good(); t[3] = 0.5
    with pytest.raises(ValueError, match=r"t .*index 3"):
        Track.from_arrays(t, r)


@pytest.mark.parametrize("which, where", [("t", (2,)), ("r", (3, 1)), ("v", (1, 2)), ("speed", (4,)), ("flow:mu", (2,))])
@pytest.mark.parametrize("value", [np.nan, np.inf])
def test_non_finite_anywhere_names_array_and_index(which, where, value):
    t, r, v, speed, flow = _good()
    if which == "t": t[where] = value
    elif which == "r": r[where] = value
    elif which == "v": v[where] = value
    elif which == "speed": speed[where] = value
    else: flow["mu"][where] = value
    name = "mu" if which == "flow:mu" else which
    with pytest.raises(ValueError, match=rf"{name}.*index {where[0]}"):
        Track.from_arrays(t, r, v=v, speed=speed, flow=flow)


def test_first_bad_index_is_reported_when_several_are_bad():
    t, r, *_ = _good(); r[2, 0] = np.nan; r[4, 0] = np.nan
    with pytest.raises(ValueError, match="index 2"):
        Track.from_arrays(t, r)


def test_negative_speed_and_shape_and_flow_key_errors():
    t, r, v, speed, flow = _good()
    s = speed.copy(); s[1] = -1.0
    with pytest.raises(ValueError, match="index 1"):
        Track.from_arrays(t, r, speed=s)
    with pytest.raises(ValueError, match="r"):
        Track.from_arrays(t, r[:, :2])
    with pytest.raises(ValueError, match="flow"):
        Track.from_arrays(t, r, flow={k: flow[k] for k in FLOW_KEYS[:-1]})        # missing r0
    with pytest.raises(ValueError, match="flow"):
        Track.from_arrays(t, r, flow={**flow, "extra": np.ones(5)})
    with pytest.raises(ValueError, match="flow"):
        Track.from_arrays(t, r, flow={**flow, "rho": np.ones(4)})
    with pytest.raises(ValueError, match="at least 2"):
        Track.from_arrays(t[:1], r[:1])


# ---- queries ----

def test_queries_outside_the_track_raise_instead_of_extrapolating():
    t, r, v, *_ = _good()
    tr = Track.from_arrays(t, r, v=v)
    for f in (tr.position, tr.velocity, tr.speed, tr.path_length):
        with pytest.raises(ValueError, match=r"\[0.0, 4.0\]"):
            f(4.5)
        with pytest.raises(ValueError):
            f(-0.1)
        with pytest.raises(ValueError):
            f(float("nan"))
    assert np.allclose(tr.position(4.0 + 1e-13), [400.0, 0, 0])      # round-off at the end is tolerated


def test_vectorised_calls_match_scalar_calls():
    t, r, v = circle(101)
    tr = Track.from_arrays(t, r, v=v)
    tt = np.array([0.0, 1.234, 5.5, 9.99, 10.0])
    assert tr.position(tt).shape == (5, 3) and tr.velocity(tt).shape == (5, 3)
    assert tr.speed(tt).shape == (5,) and tr.path_length(tt).shape == (5,)
    for i, x in enumerate(tt):
        assert np.allclose(tr.position(x), tr.position(tt)[i], atol=1e-9)
        assert tr.path_length(x) == pytest.approx(tr.path_length(tt)[i], rel=1e-12)
        assert tr.speed(x) == pytest.approx(tr.speed(tt)[i], rel=1e-12)
    assert tr.position(1.0).shape == (3,) and isinstance(tr.speed(1.0), float)


def test_path_length_on_orbit_is_arc_length_and_monotone():
    t, r, v = circle(101)
    tr = Track.from_arrays(t, r, v=v)
    tt = np.linspace(0, 10, 777)
    s = tr.path_length(tt)
    assert np.all(np.diff(s) >= 0)
    assert s[0] == 0.0 and s[-1] == pytest.approx(75000.0, rel=1e-9)
    assert tr.path_length(3.21) == pytest.approx(7500 * 3.21, rel=1e-9)


def test_default_speed_without_velocity_is_the_spline_derivative():
    t, r, v = circle(201)
    tr = Track.from_arrays(t, r)
    assert tr.speed(4.321) == pytest.approx(7500.0, rel=1e-6)
    assert np.allclose(tr.velocity(4.321), v[0] * 0 + [-7500*np.sin(7500/6.4e6*4.321), 7500*np.cos(7500/6.4e6*4.321), 0], rtol=1e-5, atol=1e-2)


def test_explicit_speed_is_air_relative_and_drives_path_length_and_U():
    t = np.linspace(0, 2, 5); r = np.c_[7000*t, 0*t, 0*t]
    speed = np.full(5, 3000.0)                          # a tailwind: air-relative speed below ground speed
    flow = {k: np.ones(5) for k in FLOW_KEYS}
    tr = Track.from_arrays(t, r, speed=speed, flow=flow)
    assert tr.speed(0.7) == pytest.approx(3000.0)
    assert tr.path_length(2.0) == pytest.approx(6000.0)
    assert tr.emission_state(0.7).U == pytest.approx(3000.0)


def test_emission_state_without_flow_raises_and_vectorised_state_works():
    t, r, v, speed, flow = _good()
    with pytest.raises(ValueError, match="flow"):
        Track.from_arrays(t, r, v=v).emission_state(1.0)
    s = Track.from_arrays(t, r, v=v, flow=flow).emission_state(np.array([0.5, 1.5]))
    assert s.rho.shape == (2,) and np.allclose(s.U, 100.0)


# ---- samples ----

def test_samples_have_roughly_the_requested_spacing_and_are_ordered():
    t, r, v = circle(101)
    tr = Track.from_arrays(t, r, v=v)
    tau, xyz = tr.samples(250.0)
    assert tau[0] == 0.0 and tau[-1] == 10.0 and np.all(np.diff(tau) > 0)
    d = np.linalg.norm(np.diff(xyz, axis=0), axis=1)
    assert np.all(d == pytest.approx(250.0, rel=1e-3))
    assert len(tau) == 301
    assert np.allclose(xyz[17], tr.position(tau[17]))


def test_samples_follow_path_length_when_speed_varies():
    t = np.linspace(0, 4, 9); r = np.c_[t**2, 0*t, 0*t]          # speed 2t, so equal path steps mean unequal time steps
    tr = Track.from_arrays(t, r, v=np.c_[2*t, 0*t, 0*t])
    tau, xyz = tr.samples(1.0)
    assert np.allclose(tr.path_length(tau), np.linspace(0, 16, len(tau)), atol=1e-6)
    assert np.allclose(xyz[:, 0], tau**2, atol=1e-9)


def test_samples_scale_to_a_long_track_and_reject_bad_spacing():
    t = np.linspace(0, 150.0, 31); r = np.c_[7000*t, 0*t, 0*t]   # about 1000 km
    tr = Track.from_arrays(t, r, v=np.tile([7000.0, 0, 0], (31, 1)))
    tau, xyz = tr.samples(0.7)
    assert len(tau) == int(np.ceil(1.05e6 / 0.7)) + 1 and xyz.shape == (len(tau), 3)
    for bad in (0.0, -1.0, float("nan")):
        with pytest.raises(ValueError, match="spacing"):
            tr.samples(bad)


def test_samples_of_a_stationary_track_are_the_two_ends():
    t = np.array([0.0, 1.0, 2.0]); r = np.zeros((3, 3))
    tau, xyz = Track.from_arrays(t, r, v=np.zeros((3, 3))).samples(1.0)
    assert np.allclose(tau, [0.0, 2.0]) and xyz.shape == (2, 3)


# ---- interpolation error ----

def test_interpolation_error_is_tiny_for_a_straight_line_and_zero_for_two_nodes():
    t = np.linspace(0, 2, 9); r = np.c_[7000*t, 0*t, 0*t]
    assert Track.from_arrays(t, r).interpolation_error() < 1e-6
    assert Track.from_arrays(t[:2], r[:2]).interpolation_error() == 0.0


def test_interpolation_error_with_velocity_uses_hermite_and_is_smaller_than_without():
    t, r, v = circle(21)
    assert Track.from_arrays(t, r, v=v).interpolation_error() < Track.from_arrays(t, r).interpolation_error()


# ---- from_callable ----

def test_from_callable_uses_the_callable_for_positions_and_samples_feed_path_length():
    R, w = 6.4e6, 7500 / 6.4e6
    pos = lambda t: np.array([R*np.cos(w*t), R*np.sin(w*t), 0.0])
    tr = Track.from_callable(pos, 0.0, 10.0, n_samples=201)
    assert tr.t0 == 0.0 and tr.t1 == 10.0
    assert np.array_equal(tr.position(3.21), pos(3.21))               # the callable itself, not an interpolant
    assert tr.path_length(10.0) == pytest.approx(75000.0, rel=1e-6)
    tau, xyz = tr.samples(500.0)
    assert xyz.shape == (len(tau), 3) and abs(len(tau) - 151) <= 1   # 150 intervals of 500 m; the ceil may add one
    assert tr.position(np.array([1.0, 2.0])).shape == (2, 3)


def test_from_callable_with_velocity_speed_and_flow_callables():
    pos = lambda t: np.array([7000.0*t, 0.0, 0.0])
    vel = lambda t: np.array([7000.0, 0.0, 0.0])
    flow = {k: (lambda t, k=k: 1.0 + t) for k in FLOW_KEYS}
    tr = Track.from_callable(pos, 0.0, 2.0, velocity=vel, speed=lambda t: 3000.0, flow=flow, n_samples=11)
    assert np.allclose(tr.velocity(1.0), [7000.0, 0, 0])
    assert tr.speed(0.3) == pytest.approx(3000.0)
    assert tr.path_length(2.0) == pytest.approx(6000.0)
    s = tr.emission_state(0.5)
    assert s.rho == pytest.approx(1.5) and s.U == pytest.approx(3000.0)


def test_from_callable_rejects_bad_window_and_non_finite_output():
    pos = lambda t: np.array([t, 0.0, 0.0])
    with pytest.raises(ValueError, match="t1"):
        Track.from_callable(pos, 1.0, 1.0)
    with pytest.raises(ValueError, match="r.*index"):
        Track.from_callable(lambda t: np.array([t, np.nan if t > 0.5 else 0.0, 0.0]), 0.0, 1.0, n_samples=11)


def test_track_module_imports_none_of_the_forbidden_modules():
    import reentry_model.parent_follower_separation_interactions.track as m, re
    src = open(m.__file__).read()
    for name in ("trajectory", "droplet_reentry", "parent_trajectory", "body", "aero", "atmosphere"):
        assert not re.search(rf"^\s*(from|import)\s+\S*\b{name}\b", src, re.M)

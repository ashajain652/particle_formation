"""dispersion.py: Girin's Eq. (1) solved numerically -- onset, asymptotes, the cached table (spec Step 3 section 9)."""
import os

import numpy as np
import pytest

from reentry_model import dispersion as dp


def test_onset_and_asymptotes():
    assert np.isnan(dp.fastest_mode(3.0)[0]) and dp.fastest_mode(3.0)[1] == 0.0        # stable at We_s = 3.00
    d, im, re = dp.fastest_mode(3.08)
    assert 0.0 < im and 0.0 < d < 0.1                                                      # unstable at 3.08, long waves
    d, im, re = dp.fastest_mode(1e4)
    assert d == pytest.approx(1.225, rel=1e-2) and im == pytest.approx(0.24, rel=4e-2)    # Girin's Fig. 5 asymptotes (measured 1.226 / 0.247)
    for we in (20.0, 100.0, 1e4):
        d, im, re = dp.fastest_mode(we)
        assert 1.5 <= re / im <= 1.7                                                       # Re Omega_f = 1.5-1.7 Im Omega_f (Girin: 1.5)
    d, im, _ = dp.fastest_mode(4.62)
    assert d == pytest.approx(0.58, rel=2e-2) and im == pytest.approx(0.049, rel=3e-2)
    d, im, _ = dp.fastest_mode(10.0)
    assert d == pytest.approx(1.34, rel=2e-2) and im == pytest.approx(0.18, rel=3e-2)


def test_cubic_coefficients_reproduce_the_relation():
    """A root of the cubic satisfies Eq. (1) written out."""
    delta, we = 1.2, 20.0
    a2, a1, a0 = dp.cubic_coefficients(delta, we)
    roots = np.roots([1.0, a2, a1, a0])
    K = (1.0 - np.exp(-2.0 * delta)) / (2.0 * delta)
    for w in roots:
        lhs = (w - delta) * ((w - delta) * (w + K * delta) + (1.0 - K) * delta)
        rhs = delta ** 3 / we * (w - K * delta)
        assert abs(lhs - rhs) < 1e-9


def test_table_is_cached_and_interpolates(tmp_path):
    path = str(tmp_path / "disp.json")
    t = dp.DispersionTable(path)
    assert os.path.isfile(path) and 3.0 < t.we_onset < 3.2
    d, im, re = t(np.array([2.0, 4.62, 100.0, 1e5]))
    assert np.isnan(d[0]) and im[0] == 0.0 and d[1] == pytest.approx(0.58, rel=3e-2) and d[3] == pytest.approx(1.226, rel=1e-2)
    assert im[2] == pytest.approx(0.241, rel=2e-2)
    packaged = dp.DispersionTable()
    assert np.allclose(packaged.we, t.we) and np.allclose(packaged.im_omega_f, t.im_omega_f, atol=1e-6)     # the committed table is current


# The broken-line solver (spec 6.2, stage 1). The helpers below work in each source's own units: lengths in layer
# depths, speeds in surface speeds and densities in the top layer's, so a wavenumber k is Girin's Delta and a growth
# rate k Im c is his Im Omega.

def girin(we):            # Girin's units: delta_m = V_s = rho = 1, sigma = 1 / We_s
    return dp.KinkedProfiles.one([1.0], [1.0, 0.0, 0.0]), 1.0 / we


def slurry_over_wall(beta, eps, r=1.0, we=1e4):   # notes fulltext_girin_papers.md, lines 79-106
    S1 = 1.0 / (1.0 + eps * (beta - 1.0)); u_h = eps * (beta - 1.0) * S1
    return dp.KinkedProfiles.one([1.0], [1.0, u_h, 0.0], wall=beta, rho=[1.0, r]), 1.0 / we, S1, u_h


def ghosh(K2, d=0.4):     # Ghosh 2017 in units of V and H: kink at depth d, top slope K2 times the bottom's, sigma 0
    S1 = 1.0 / (d + (1.0 - d) / K2)
    return dp.KinkedProfiles.one([d], [1.0, (1.0 - d) * S1 / K2, 0.0], wall=1.0), 0.0


def stack(*profiles):
    """One batch from single profiles of any breakpoint count, padded with NaN to the largest count."""
    B = max(p.depth.shape[1] for p in profiles)

    def pad(a, width):
        out = np.full((a.shape[0], width), np.nan)
        out[:, :a.shape[1]] = a
        return out
    return dp.KinkedProfiles(count=np.concatenate([p.count for p in profiles]),
                             depth=np.concatenate([pad(p.depth, B) for p in profiles]),
                             u_top=np.concatenate([p.u_top for p in profiles]),
                             u_break=np.concatenate([pad(p.u_break, B) for p in profiles]),
                             wall=np.concatenate([p.wall for p in profiles]),
                             u_wall=np.concatenate([p.u_wall for p in profiles]),
                             rho=np.concatenate([pad(p.rho, B + 1) for p in profiles]))


def assert_same_roots(got, want, atol):
    """Each wanted root has its own computed root within atol (matched one to one, so a double root counts twice)."""
    got, want = list(got), list(want)
    assert len(got) == len(want), (got, want)
    for w in want:
        j = int(np.argmin([abs(g - w) for g in got]))
        assert abs(got[j] - w) <= atol, (w, got)
        got.pop(j)


def test_kinked_profiles_build_validate_and_give_slopes():
    assert dp.DISPERSION_NAMES == ("table", "profile") and dp.DEFAULT_DISPERSION == "profile"
    p = dp.KinkedProfiles.one([0.5, 1.0], [2.0, 1.0, 0.5, 0.0], wall=2.0, rho=[1.0, 2.0, 3.0])
    assert p.count.tolist() == [2] and p.depth.shape == (1, 2) and p.u_break.shape == (1, 2) and p.rho.shape == (1, 3)
    assert np.allclose(p.slopes(), [[2.0, 1.0, 0.5]])                       # dU/dy with y up: speed falling with depth
    assert np.allclose(dp.KinkedProfiles.one([0.5], [1.0, 2.0, 0.0], wall=1.0).slopes(), [[-2.0, 4.0]])   # a plug
    unbounded = dp.KinkedProfiles.one([1.0], [1.0, 0.0, 0.0])
    assert np.isinf(unbounded.wall[0]) and np.allclose(unbounded.slopes(), [[1.0, 0.0]])
    assert np.allclose(dp.KinkedProfiles.one([1.0], [1.0, 0.5, 0.5], rho=[1.0, 2.0]).rho, [[1.0, 2.0]])
    for depths, speeds, wall in (([1.0, 0.5], [1.0, 0.5, 0.2, 0.0], 2.0),    # the depths do not increase
                                 ([1.0], [1.0, 0.5, 0.0], 1.0),              # the wall is not below the last breakpoint
                                 ([1.0], [1.0, 0.5, 0.2], np.inf),           # an unbounded bottom that is not uniform
                                 ([], [1.0, 1.0], np.inf)):                  # no breakpoint and no wall
        with pytest.raises(ValueError):
            dp.KinkedProfiles.one(depths, speeds, wall=wall)
    # The dataclass itself accepts an unresolved row (count 0 and NaN everywhere), which a later fit builds.
    nan = np.array([np.nan])
    unresolved = dp.KinkedProfiles(count=np.array([0]), depth=np.full((1, 2), np.nan), u_top=nan, u_break=np.full((1, 2), np.nan),
                                   wall=nan, u_wall=nan, rho=np.full((1, 3), np.nan))
    batch = stack(p, unresolved)
    assert batch.count.tolist() == [2, 0] and batch.rows([1]).count.tolist() == [0]
    assert np.allclose(batch.rows(np.array([True, False])).slopes(), p.slopes())
    with np.errstate(all="raise"):                                           # callers filter them, but nothing faults
        assert np.isnan(dp.broken_line_eigenvalues(batch, [1.0, 3.0], 0.01)[1]).all()
        w = dp.fastest_waves(batch, 0.01)
    assert w.growth[1] == 0.0 and np.isnan(w.k[1]) and np.isnan(w.release_depth[1])


def test_girin_profile_gives_the_roots_of_his_cubic():
    for delta, we in ((0.5, 5.0), (1.2, 20.0), (1.225, 1e4), (3.0, 50.0)):
        prof, sigma = girin(we)
        c = dp.broken_line_eigenvalues(prof, [delta], sigma)[0, 0]
        assert c.shape == (3,) and np.isfinite(c).sum() == 3
        assert_same_roots(c * delta, np.roots([1.0, *dp.cubic_coefficients(delta, we)]), atol=1e-10)


def test_girin_profile_fastest_wave_matches_the_cached_table():
    prof, sigma = girin(1e4)
    w = dp.fastest_waves(prof, sigma)
    delta_table, im_table = dp.fastest_mode(1e4)[:2]
    assert abs(w.k[0] - delta_table) <= 0.002                                  # the table's grid spacing
    assert im_table - 1e-9 <= w.growth[0] <= im_table + 1e-4
    assert w.growth[0] == pytest.approx(w.k[0] * w.c[0].imag, rel=1e-12)


def test_girin_threshold_is_three():
    prof, sigma = girin(3.00)
    assert np.all(dp.broken_line_c(prof, np.logspace(-4, 0.7, 400), sigma).imag <= 1e-10)
    # Measured on the cubic, 2026-11-11: at We_s 3.01 the unstable band is Delta below 0.0044, with Im Omega at most
    # 3.2e-6, below the fastest-wave grid, so the eigenvalue function is called directly.
    prof, sigma = girin(3.01)
    assert dp.broken_line_c(prof, 0.0033, sigma)[0, 0].imag > 1e-4


def test_slurry_over_a_wall_gives_the_closed_form_roots():
    """Every finite eigenvalue c makes the closed form, multiplied by (1 - c)(u_h - c), vanish to 1e-10 relative to its
    largest term: (1 - c)(Delta tanh Delta + q) - [S1 + Delta^2 / (We (1 - c))](1 + (q / Delta) tanh Delta) = 0 with
    q = r Delta coth(Delta (beta - 1)) + S1 (1 - r eps) / (u_h - c)."""
    we = 1e4
    for beta, eps, r in ((1.5, 0.0, 1.0), (2.0, 1e-3, 1.0), (5.0, 1e-3, 1.1), (20.0, 1e-3, 0.95)):
        prof, sigma, S1, u_h = slurry_over_wall(beta, eps, r, we)
        for delta in (0.5, 1.2, 3.0):
            cs = dp.broken_line_eigenvalues(prof, [delta], sigma)[0, 0]
            assert np.isfinite(cs).sum() == 3
            T, C, a = np.tanh(delta), 1.0 / np.tanh(delta * (beta - 1.0)), S1 * (1.0 - r * eps)
            for c in cs:
                w, v = 1.0 - c, u_h - c
                terms = [w * w * v * delta * T, w * w * v * r * delta * C, w * w * a, -S1 * w * v, -S1 * w * r * C * T * v,
                         -S1 * w * a * T / delta, -delta ** 2 / we * v, -delta ** 2 / we * r * C * T * v,
                         -delta ** 2 / we * a * T / delta]
                assert abs(sum(terms)) <= 1e-10 * max(abs(t) for t in terms), (beta, eps, r, delta, c)


def test_slurry_over_a_wall_fastest_waves():
    """The side session's values, at equal densities and no slurry shear (its notes give beta = 20 with a shear ratio of
    1e-3 as 0.242, outside this tolerance)."""
    for beta, k, growth in ((1.5, 1.21, 0.196), (2.0, 1.24, 0.234), (20.0, 1.225, 0.247)):
        prof, sigma, _, _ = slurry_over_wall(beta, 0.0)
        w = dp.fastest_waves(prof, sigma)
        assert abs(w.k[0] - k) <= 0.01 and abs(w.growth[0] - growth) <= 0.002, (beta, w.k[0], w.growth[0])


def test_ghosh_limits():
    w = dp.fastest_waves(*ghosh(1.5))
    assert abs(w.growth[0] - 0.229) <= 0.002 and abs(w.k[0] - 2.6) <= 0.15
    for K2, growth in ((1.1, 0.109), (2.5, 0.349), (3.0, 0.383)):
        assert abs(dp.fastest_waves(*ghosh(K2)).growth[0] - growth) <= 0.003, K2
    # More viscous fluid on top (a smaller top slope), and a single straight film on a wall, are stable (Ghosh, Eq. 20).
    k = np.logspace(-2, 2, 200)
    for prof in (ghosh(0.5)[0], dp.KinkedProfiles.one([], [1.0, 0.0], wall=1.0)):
        for sigma in (0.0, 0.01, 1.0):
            assert np.all(dp.broken_line_c(prof, k, sigma).imag <= 1e-10)
            assert dp.fastest_waves(prof, sigma).growth[0] == 0.0


def test_release_depth_on_girin_profile():
    prof, sigma = girin(1e4)
    w = dp.fastest_waves(prof, sigma)
    y = np.linspace(0.0, 4.0, 4001)
    phi, disp = dp.wave_structure(prof, w.k, w.c, y)
    assert phi[0, 0] == pytest.approx(1.0, abs=1e-12) and disp[0, 0] == pytest.approx(1.0, abs=1e-12)
    i = int(np.argmax(disp[0]))
    assert disp[0, i] == pytest.approx(1.6, abs=0.05) and y[i] == pytest.approx(0.67, abs=0.05)
    assert w.release_depth[0] == pytest.approx(1.82, abs=0.02)
    a = np.abs(phi[0])
    j = int(np.argmax(a <= np.exp(-1.0)))
    y_phi = y[j - 1] + (a[j - 1] - np.exp(-1.0)) * (y[j] - y[j - 1]) / (a[j - 1] - a[j])
    assert y_phi == pytest.approx(1.29, abs=0.02)


def test_a_wall_closing_on_the_kink_kills_the_growth():
    near = dp.fastest_waves(*slurry_over_wall(1.001, 1e-3)[:2])
    far = dp.fastest_waves(*slurry_over_wall(20.0, 1e-3)[:2])
    assert far.growth[0] > 0.2 and near.growth[0] < 1e-2 * far.growth[0]


def test_deep_walls_and_short_waves_stay_accurate():
    """A wall 100/k below the kink (ruling R1: at depth 1 + 100/k, since a wall at 100/k would meet the kink at k = 100)
    changes nothing; segments 1e-3 and 10 thick at kH = 1000 give finite eigenvalues without a floating-point fault."""
    with np.errstate(all="raise"):
        for we in (20.0, 1e4):
            open_bottom, sigma = girin(we)
            for k in (1.0, 10.0, 100.0):
                walled = dp.KinkedProfiles.one([1.0], [1.0, 0.0, 0.0], wall=1.0 + 100.0 / k)
                want = dp.broken_line_eigenvalues(open_bottom, [k], sigma)[0, 0]
                got = dp.broken_line_eigenvalues(walled, [k], sigma)[0, 0]
                assert_same_roots(got, want, atol=1e-8 * np.max(np.abs(want)))
        thin_thick = dp.KinkedProfiles.one([1e-3], [1.0, 0.5, 0.0], wall=10.0 + 1e-3)
        c = dp.broken_line_eigenvalues(thin_thick, [1000.0 / 10.001], 1e-6)
        assert np.isfinite(c).all()
        w = dp.fastest_waves(thin_thick, 1e-6)          # its grid reaches kH = 5e4; the top layer's We_s is 250
        assert np.isfinite([w.k[0], w.growth[0], w.release_depth[0]]).all() and 0.0 < w.release_depth[0] <= 10.001


def test_degenerate_pencils_give_finite_eigenvalues_only():
    """Review Focus 5: zero surface tension, a wall 1e-12 below a breakpoint and breakpoints 1e-9 apart (smooth across
    the gap, and with a jump across it) give finite eigenvalues; NaN appears only as padding."""
    profiles = stack(dp.KinkedProfiles.one([0.5], [1.0, 0.4, 0.0], wall=0.5 + 1e-12),
                     dp.KinkedProfiles.one([0.5, 0.5 + 1e-9], [1.0, 0.5, 0.5 - 1e-9, 0.0], wall=1.0),
                     dp.KinkedProfiles.one([0.5, 0.5 + 1e-9], [1.0, 0.5, 0.2, 0.0], wall=1.0),
                     dp.KinkedProfiles.one([], [1.0, 0.0], wall=1.0))
    with np.errstate(all="raise"):
        c = dp.broken_line_eigenvalues(profiles, np.logspace(-2, 3, 60), 0.0)
        w = dp.fastest_waves(profiles, 0.0)
    assert c.shape == (4, 60, 4) and not np.isinf(c).any()
    for i, n in enumerate(profiles.count):
        assert np.isfinite(c[i, :, :n + 2]).all() and np.isnan(c[i, :, n + 2:]).all()
    assert np.isfinite(w.growth).all()


def test_the_generalised_route_agrees_with_the_standard_one(monkeypatch):
    """The per-column scipy.linalg.eig(A, B) route, used where B is ill conditioned, gives the batched route's values."""
    profiles = stack(girin(20.0)[0], slurry_over_wall(2.0, 1e-3, 1.1)[0], ghosh(1.5)[0],
                     dp.KinkedProfiles.one([0.3, 0.5], [1.0, 0.6, 0.1, 0.0], wall=1.5, rho=[1.0, 1.2, 1.5]))
    k = np.logspace(-1.5, 1.5, 7)
    sigma = np.array([0.05, 1e-4, 0.0, 0.02])
    standard = dp.broken_line_eigenvalues(profiles, k, sigma)
    monkeypatch.setattr(dp, "PENCIL_COND_MAX", 0.0)
    general = dp.broken_line_eigenvalues(profiles, k, sigma)
    assert np.array_equal(np.isnan(standard), np.isnan(general))
    for i, n in enumerate(profiles.count):
        for j in range(len(k)):
            want = standard[i, j, :n + 2]
            assert_same_roots(general[i, j, :n + 2], want, atol=1e-10 * np.max(np.abs(want)))

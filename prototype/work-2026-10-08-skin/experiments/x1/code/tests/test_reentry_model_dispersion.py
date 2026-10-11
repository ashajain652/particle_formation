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

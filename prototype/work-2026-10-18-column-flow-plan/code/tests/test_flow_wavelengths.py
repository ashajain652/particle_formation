"""The share of the columns' flow at short wavelengths (sub-plan 21 Task 8, `analysis/flow_wavelengths.py`): the
measure of the depth-uniform pressure's limitation (layered-skin design section 14, assumption 1)."""
import os
import sys

import numpy as np
import pytest

from reentry_model import columns

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "analysis"))
import flow_wavelengths as fw  # noqa: E402


def cylinder(N=40, M=4, dx=1e-3):
    """An open cylinder of N x M square facets of side dx, each split in two: flat to itself and periodic around, so
    the cotangent Laplacian is the five-point one (weights 1 on the squares' sides, 0 on their diagonals, half on the
    open ends' edges, whose columns have half the area) and a sampled cosine around it is one of its modes exactly."""
    R = dx / (2.0 * np.sin(np.pi / N))
    a = 2.0 * np.pi * np.arange(N) / N
    points = np.array([(R * np.cos(t), R * np.sin(t), r * dx) for r in range(M + 1) for t in a])
    at = lambda r, j: r * N + j % N
    faces = [f for r in range(M) for j in range(N)
             for f in ((at(r, j), at(r, j + 1), at(r + 1, j + 1)), (at(r, j), at(r + 1, j + 1), at(r + 1, j)))]
    return columns.Columns(points, np.array(faces))


def test_the_overstatement_factors_are_the_inertial_and_viscous_layers_rate_ratios():
    """sqrt(kh / tanh kh) (spec 14.1) and the lubrication rate over the Stokes rate of a viscous layer on a rigid base,
    by hand at kh = 1: sqrt(1 / tanh 1) = 1.14588 and (2/3)(cosh 2 + 3) / (sinh 2 - 2) = 2.77112; both 1 at kh = 0, the
    viscous one's series 1 + 9/5 (kh)^2 meeting the formula at kh = 0.05 to 3e-7; 1.05 at kh = 0.560 (KH_SHORT) and
    0.167."""
    assert fw.inertial_factor(1.0) == pytest.approx(np.sqrt(1.0 / np.tanh(1.0)), rel=1e-15)
    assert fw.viscous_factor(1.0) == pytest.approx((2.0 / 3.0) * (np.cosh(2.0) + 3.0) / (np.sinh(2.0) - 2.0), rel=1e-14)
    assert fw.inertial_factor(0.0) == 1.0 and fw.viscous_factor(0.0) == 1.0
    assert fw.viscous_factor(0.05) == pytest.approx(fw.viscous_factor(0.0500001), abs=1e-6)
    assert fw.inertial_factor(fw.KH_SHORT) == pytest.approx(1.05, abs=1e-3)
    assert fw.kh_at(fw.viscous_factor, 1.05) == pytest.approx(0.167, abs=1e-3)


@pytest.mark.parametrize("m", [1, 4, 8])
def test_a_single_sinusoid_on_a_flat_layer_has_the_share_and_excess_of_its_own_wavenumber(m):
    """cos(2 pi m j / N) around the cylinder of 40 facets of 1 mm on a layer 1.2 mm deep: all its power is in its own
    mode, of the five-point Laplacian's wavenumber k = (2 / dx) sin(pi m / N) -- 156.9, 618.0 and 1 175.6 1/m, kh 0.19,
    0.74 and 1.41 -- so its short-wave share is 0 below KH_SHORT and 1 above, and its mean excesses are the factors'
    at that kh less one, to 1e-15 and 3e-15 (measured 2026-10-24); its power is its area-weighted square, half the
    cylinder's 160 mm2."""
    c = cylinder()
    j = np.arange(c.n) % 40
    res = fw.wavelength_share(c, np.cos(2.0 * np.pi * m * j / 40), np.full(c.n, 1.2e-3))
    kh = 2.0 / 1e-3 * np.sin(np.pi * m / 40) * 1.2e-3
    assert res["share"] == pytest.approx(1.0 if kh > fw.KH_SHORT else 0.0, abs=1e-15)    # 2.5e-30 leaks at m = 1
    assert res["inertial"] == pytest.approx(float(fw.inertial_factor(kh)) - 1.0, rel=1e-12, abs=1e-15)
    assert res["viscous"] == pytest.approx(float(fw.viscous_factor(kh)) - 1.0, rel=1e-12, abs=1e-15)
    assert res["power"].sum() == pytest.approx(80e-6, rel=1e-12)
    assert res["k"][np.argmax(res["power"])] == pytest.approx(kh / 1.2e-3, rel=1e-12)


def test_the_share_of_two_waves_is_the_short_one_s_power_over_both_and_a_vector_field_counts_every_component():
    """a cos of m = 1 (kh 0.19) plus b cos of m = 8 (kh 1.41): the share is b^2 / (a^2 + b^2) and the inertial excess
    the factors' power-weighted mean, by hand; a field of three Cartesian components sums their powers, so the same
    waves along x and z give the same share; and over a region the modes are the region's own, its columns' depths
    weighted by each mode's square."""
    c = cylinder()
    j = np.arange(c.n) % 40
    lo, hi = np.cos(2.0 * np.pi * j / 40), np.cos(2.0 * np.pi * 8 * j / 40)
    a, b = 2.0, 0.5
    res = fw.wavelength_share(c, a * lo + b * hi, np.full(c.n, 1.2e-3))
    kl, kh_ = (2e3 * np.sin(np.pi * m / 40) * 1.2e-3 for m in (1, 8))
    assert res["share"] == pytest.approx(b ** 2 / (a ** 2 + b ** 2), rel=1e-12)
    mean = (a ** 2 * (fw.inertial_factor(kl) - 1.0) + b ** 2 * (fw.inertial_factor(kh_) - 1.0)) / (a ** 2 + b ** 2)
    assert res["inertial"] == pytest.approx(float(mean), rel=1e-12)
    vec = np.stack([a * lo, np.zeros(c.n), b * hi], axis=1)
    assert fw.wavelength_share(c, vec, np.full(c.n, 1.2e-3))["share"] == pytest.approx(res["share"], rel=1e-12)
    ring = np.arange(c.n) < 120                                       # the first three rings: a cylinder 2 mm long
    depth = np.where(np.arange(c.n) < 80, 1e-3, 2e-3)
    sub = fw.wavelength_share(c, hi, depth, region=ring)
    k, phi, idx = fw.modes(c, ring)
    assert len(idx) == 120 and sub["h"] == pytest.approx((c.area[idx][:, None] * phi ** 2 * depth[idx][:, None]).sum(0))

"""Limits of singh_drag (from the prototype's check_singh_drag.py) and its agreement with the earlier composed law of
particle_trajectory in fast free-molecule flow."""
import math

import pytest

from reentry_model.droplet_reentry import particle_trajectory as pt
from reentry_model.droplet_reentry import singh_drag as sg


def test_stokes_limit():
    cd, _ = sg.drag_coefficient(1e-8, 1e-4, 1.0)
    assert cd / 2.4e5 == pytest.approx(1.0, abs=0.01)


def test_free_molecule_limit_is_schaaf_chambre():
    cd, g = sg.drag_coefficient(10.0, 1e-3, 1.0)
    assert cd / g["cd_fm"] == pytest.approx(1.0, abs=1e-3)


def test_hypersonic_continuum():
    cd, _ = sg.drag_coefficient(20.0, 1e7, 1.0)
    assert cd == pytest.approx(0.9, abs=0.05)


@pytest.mark.parametrize("re", [1.0, 10.0, 30.0])
def test_low_mach_against_schiller_naumann(re):
    cd, _ = sg.drag_coefficient(0.01, re, 1.0)
    assert cd / (24 / re * (1 + 0.15 * re ** 0.687)) == pytest.approx(1.0, abs=0.1)


def test_continuous_at_mach_1():
    lo, _ = sg.drag_coefficient(0.999999, 50.0, 1.0)
    hi, _ = sg.drag_coefficient(1.000001, 50.0, 1.0)
    assert hi / lo == pytest.approx(1.0, abs=1e-4)


@pytest.mark.parametrize("tw", [220.0, 1800.0])
def test_matches_composed_law_in_fast_free_molecule_flow(tw):
    """Air at 220 K and 2 Pa (about 77 km), Knudsen number 100, speed ratio 25: the two laws are the same physics."""
    T, p = 220.0, 2.0
    rho, mu = p / (pt.R_AIR * T), pt.viscosity(T)
    lam = mu / rho * math.sqrt(math.pi / (2 * pt.R_AIR * T))
    d, V = lam / 100.0, 25.0 * math.sqrt(2 * pt.R_AIR * T)
    c_s, _ = sg.drag_coefficient(V / math.sqrt(pt.GAMMA * pt.R_AIR * T), rho * V * d / mu, tw / T)
    c_c, _ = pt.drag_coefficient(V, rho, T, p, d, tw / T, "blend")
    assert c_s / c_c == pytest.approx(1.0, abs=0.01)

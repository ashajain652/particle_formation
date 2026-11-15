"""The twin of the DRAMA parent (design spec section 4.3): the ballistic-coefficient table read from the DRAMA file,
the 1 percent speed check against DRAMA above 60 km, a release that starts exactly on the twin's track, and a
release altitude outside the twin's range (Review Focus 5)."""
import numpy as np
import pytest

from reentry_model.wake_residence.adapters import drama_beta_table, release_state, twin_vs_drama
from reentry_model.wake_residence.overlap import Plume


def test_beta_table_matches_builder_check(parent):
    b = drama_beta_table(parent)
    assert 120 < b.beta_at(85e3) < 160 and 20 < b.beta_at(60e3) < 40
    assert b.mass_at(90e3) == pytest.approx(740.0, rel=1e-3) and b.mass_at(66e3) < 300


def test_twin_reproduces_drama_speed_within_one_percent_above_60km(twin, parent):
    rows = [r for r in twin_vs_drama(twin, parent) if r["drama_altitude_km"] > 60.0]
    worst = max(abs(r["speed_rel_diff"]) for r in rows)
    print("twin vs DRAMA worst speed difference above 60 km:", worst)
    assert worst <= 0.01


def test_release_starts_exactly_on_the_track(twin):
    t_r = twin.t_at_altitude(80.0)
    st = release_state(twin.track, t_r, twin.epoch)
    from reentry_model.shared import earth
    x = earth.geodetic_to_ecef(st.altitude, st.lat, st.lon)
    m = Plume(twin.track).margin(x, t_r)
    # on the track at the parent's position: the union of balls gives r_b, and the front plane (Task 6b) caps it at 0
    assert m.r_b - m.distance == pytest.approx(m.r_b, abs=1e-3) and m.distance < 1e-3
    assert m.margin == pytest.approx(0.0, abs=1e-3) and m.front_limited


def test_altitude_outside_track_raises(twin):
    with pytest.raises(ValueError, match="km"):
        twin.t_at_altitude(130.0)

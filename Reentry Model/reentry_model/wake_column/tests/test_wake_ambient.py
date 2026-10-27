import pytest

from reentry_model.wake_column.ambient import AmbientColumn


@pytest.mark.parametrize("src", ["waccm", "msis_afgl"])
def test_columns_cover_45_to_100_km_and_close(src):
    for tod in ("day", "night"):
        col = AmbientColumn.load(src, tod)
        for z in (45.0, 60.0, 76.0, 90.0, 100.0):
            a = col.at(z)
            assert 130 < a.T_K < 320
            assert sum(a.n.values()) == pytest.approx(a.n_total_m3, rel=0.02)
            assert all(v >= 0 for v in a.n.values())


def test_day_night_atomic_oxygen_contrast_below_80_km():
    day, night = AmbientColumn.load("waccm", "day").at(70.0), AmbientColumn.load("waccm", "night").at(70.0)
    assert day.n["O"] / night.n["O"] > 100
    assert day.n["H2O"] == pytest.approx(night.n["H2O"], rel=0.3)


def test_waccm_and_msis_agree_on_total_density_within_factor_2():
    for z in (60.0, 76.0, 90.0):
        a, b = AmbientColumn.load("waccm", "night").at(z), AmbientColumn.load("msis_afgl", "night").at(z)
        assert 0.5 < a.n_total_m3 / b.n_total_m3 < 2.0


def test_outside_table_raises_with_altitude_in_message():   # Review Focus 2
    with pytest.raises(ValueError, match="44.0 km"):
        AmbientColumn.load("waccm", "day").at(44.0)

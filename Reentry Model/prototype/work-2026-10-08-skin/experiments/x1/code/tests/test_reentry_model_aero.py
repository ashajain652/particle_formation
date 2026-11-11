"""aero.py: Knudsen/Mach numbers and SESAM's sphere drag (tables + Knudsen bridging)."""
import csv
import math
import os

import numpy as np
import pytest

from reentry_model import aero
from reentry_model import constants as c

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REFERENCE_DIR = os.path.join(REPO_ROOT, "data", "reference_runs")
REFERENCE_CSVS = sorted(f for f in os.listdir(REFERENCE_DIR) if f.endswith(".csv"))


def test_knudsen_reproduces_sesam_at_the_break_off_altitude():
    # SESAM: rho = 2.727e-5 kg/m3 at 77.5 km (US76), D = 100 mm -> knudsen column 0.0298
    kn = aero.knudsen(2.727e-5, c.M_BAR_AIR, 0.100)
    assert kn == pytest.approx(0.0298, rel=0.02)
    assert aero.mean_free_path(2.727e-5, c.M_BAR_AIR) == pytest.approx(2.98e-3, rel=0.02)


def test_speed_of_sound_and_mach_at_sea_level():
    assert aero.speed_of_sound(288.15, c.M_BAR_AIR) == pytest.approx(340.3, abs=0.3)
    assert aero.mach(7500.0, 288.15, c.M_BAR_AIR) == pytest.approx(7500.0 / 340.3, rel=1e-3)


class TestTables:
    def test_loads_the_bundled_atdb(self):
        t = aero.SphereDragTables.from_json()
        assert list(t.mach) == [5.0, 10.0, 15.0, 20.0, 25.0, 30.0]
        assert t.cd_free_molecular(5.0) == 2.360635 and t.cd_continuum(30.0) == 0.913599

    def test_interpolates_and_clamps_above_mach_30(self):
        t = aero.SphereDragTables.from_json()
        assert t.cd_continuum(27.5) == pytest.approx(0.5 * (0.913411 + 0.913599))
        assert t.cd_free_molecular(21.07) == pytest.approx(2.062249 + (21.07 - 20.0) / 5.0 * (2.046842 - 2.062249))
        assert t.cd_continuum(35.0) == 0.913599 and t.cd_free_molecular(40.0) == 2.036933

    def test_below_mach_5_matches_sesam_clamping(self):
        t = aero.SphereDragTables.from_json()
        assert t.cd_continuum(3.0) == 0.898818 and t.cd_continuum(1.02) == 0.898818
        assert t.cd_continuum(1.0) == pytest.approx(0.75 * 0.898818) and t.cd_continuum(0.98) == 0.5 * 0.898818   # the Ma-1 step smoothed over +-2 % (Step 3)
        assert t.cd_continuum(0.5) == pytest.approx(0.449409, abs=1e-6)     # exactly half the Ma-5 value
        assert t.cd_free_molecular(2.0) == 2.360635


class TestBridging:
    def test_sesam_table_limits_midpoint_and_monotonicity(self):
        f = aero.SesamTable()
        assert f(1e-3) == 0.0 and f(10.0) == 1.0 and f(41.0) == 1.0
        assert f(0.143) == pytest.approx(0.50, abs=0.03)
        kns = np.logspace(-3, 2, 200)
        values = [f(k) for k in kns]
        assert all(b >= a for a, b in zip(values, values[1:]))

    def test_sesam_erf_values(self):
        f = aero.SesamErf()
        assert f(0.01) == pytest.approx(0.0026, abs=5e-4)
        assert f(0.143) == pytest.approx(0.50, abs=0.01)
        assert f(1.0) == pytest.approx(0.98, abs=0.01)

    def test_sin2_textbook_form(self):
        f = aero.Sin2()
        assert f(0.001) == 0.0 and f(0.01) == 0.0 and f(1.0) == 1.0 and f(5.0) == 1.0
        assert f(0.1) == pytest.approx(0.5)

    def test_matting_is_not_available_yet(self):
        with pytest.raises(NotImplementedError):
            aero.Matting()(0.1)

    def test_factory(self):
        assert isinstance(aero.bridging_by_name("sesam-table"), aero.SesamTable)
        assert isinstance(aero.bridging_by_name("textbook"), aero.Sin2)
        with pytest.raises(ValueError):
            aero.bridging_by_name("legge")


class TestDragCoefficient:
    def test_limits(self):
        t = aero.SphereDragTables.from_json()
        assert aero.drag_coefficient(1e-4, 0.5, t, aero.SesamTable()) == pytest.approx(0.449409, abs=1e-6)
        assert aero.drag_coefficient(41.0, 21.07, t, aero.SesamTable()) == pytest.approx(2.059, abs=2e-3)   # R50 first row
        assert aero.drag_coefficient(0.0596, 26.23, t, aero.SesamTable()) == pytest.approx(1.125, abs=0.02) # 50 mm sweep, t = 0
        assert aero.drag_coefficient(0.0, 26.0, t, aero.SesamTable()) == t.cd_continuum(26.0)

    @pytest.mark.parametrize("name", REFERENCE_CSVS)
    def test_reproduces_the_reference_drag_column(self, name):
        t = aero.SphereDragTables.from_json()
        f = aero.SesamTable()
        rows = list(csv.DictReader(open(os.path.join(REFERENCE_DIR, name))))
        hyper, low = [], []
        for r in rows:
            kn, ma, cd = float(r["knudsen"]), float(r["mach"]), float(r["drag"])
            model = aero.drag_coefficient(kn, ma, t, f)
            if ma >= 5.0:
                hyper.append(model - cd)
            elif abs(ma - 1.0) > 0.02:                 # SESAM's subsonic switch sits at Ma = 1 (3-decimal columns)
                low.append(model - cd)
        hyper = np.array(hyper)
        assert math.sqrt(np.mean(hyper ** 2)) <= 0.010, name
        assert np.abs(hyper).max() <= 0.030, name
        if low:                                        # the melting references (Step 3) demise while still hypersonic
            assert np.abs(low).max() <= 1.5e-3, name   # the reference prints C_D with 3 decimals

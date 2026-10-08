"""atmosphere.py: US76 table, NRLMSISE-00 (pymsis), SESAM replay, winds."""
import math
import os
from datetime import date

import numpy as np
import pytest

from reentry_model import aero, atmosphere, fap, sesam_io
from reentry_model import constants as c

R100 = "sphere_d100.00mm_T0300.0K_v07.50000kms_h077.500km_mAA7075_nomelt_msis_nowind"
R50 = "sphere_d050.00mm_T0300.0K_v07.50000kms_h115.000km_mAA7075_nomelt_msis_nowind"
LAT, LON = math.radians(29.546067), math.radians(-82.134333)


def ref(name):
    return sesam_io.load_reference(os.path.join(sesam_io.REFERENCE_DIR, name + ".csv"))


class TestUS76:
    def test_sea_level_and_70km(self):
        atm = atmosphere.US76TableAtmosphere()
        s0 = atm.state(0.0, 0.0, LAT, LON)
        assert s0.rho == pytest.approx(1.225, abs=1e-3) and s0.T == pytest.approx(288.15, abs=0.01)
        assert s0.p == pytest.approx(101325.0, abs=1.0) and s0.m_bar == c.M_BAR_AIR
        assert np.all(s0.wind_enu == 0.0)
        s70 = atm.state(0.0, 70000.0, LAT, LON)
        assert s70.rho == pytest.approx(8.283e-5, rel=2e-3) and s70.T == pytest.approx(219.6, abs=0.1)

    def test_matches_sesam_static_density_at_break_off(self):
        atm = atmosphere.US76TableAtmosphere()
        assert atm.state(0.0, 77500.133, LAT, LON).rho == pytest.approx(2.727e-5, rel=0.01)

    def test_log_linear_between_rows_and_top_of_table(self):
        atm = atmosphere.US76TableAtmosphere()
        a, m, b = (atm.state(0.0, h, LAT, LON).rho for h in (70000.0, 70050.0, 70100.0))
        assert m == pytest.approx(math.sqrt(a * b), rel=1e-6)
        with pytest.raises(ValueError):
            atm.state(0.0, 150001.0, LAT, LON)

    def test_static_wind_profile(self):
        wind = atmosphere.StaticProfileWind()
        w = wind.wind_enu(60000.0)
        assert w[0] == pytest.approx(10.90, abs=0.01) and w[1] == pytest.approx(-1.30, abs=0.01) and w[2] == 0.0
        atm = atmosphere.US76TableAtmosphere(wind=wind)
        assert np.allclose(atm.state(0.0, 60000.0, LAT, LON).wind_enu, w)


class TestNRLMSISE00:
    def test_reproduces_sesam_density_at_the_reference_starts(self):
        # Informative band, not a tight tolerance: measured ratios s.rho/r.density[0] are 1.057 (R100,
        # 77.5 km) and 0.858 (R50, 115 km) (2026-09-17), and are insensitive to the solar-index convention
        # (tried F10.7 170-246, F10.7a 170-234, Ap 8-56, storm-mode 3-hourly ap, pymsis historical indices -
        # all give 1.05-1.07 and 0.84-0.86) - the gap is SESAM's NRLMSISE-00 implementation/switches vs.
        # pymsis, not a fap.py convention bug; its effect on the trajectory is measured in Task 10.
        solar = fap.solar_indices(fap.load_fap_day(fap.DEFAULT_FAP_DAY), date(2024, 8, 1))
        for name in (R100, R50):
            r = ref(name)
            atm = atmosphere.NRLMSISE00Atmosphere(r.initial.epoch, solar)
            s = atm.state(0.0, r.initial.altitude, r.initial.lat, r.initial.lon)
            ratio = s.rho / r.density[0]
            assert 0.80 <= ratio <= 1.25, (name, s.rho, r.density[0], ratio)
            assert s.T > 150.0 and s.p > 0.0

    def test_mean_molecular_mass_falls_off_above_100km(self):
        solar = fap.SolarIndices(f107=234.0, f107a=194.0, ap=19.0)
        atm = atmosphere.NRLMSISE00Atmosphere(ref(R50).initial.epoch, solar)
        low = atm.state(0.0, 77500.0, LAT, LON).m_bar / c.ATOMIC_MASS_UNIT
        high = atm.state(0.0, 115000.0, LAT, LON).m_bar / c.ATOMIC_MASS_UNIT
        assert low == pytest.approx(28.96, abs=0.2)
        assert 25.0 < high < 28.0                       # SESAM's Kn column implies ~26.6 u at 115 km


class TestReplay:
    @pytest.mark.parametrize("name", [R100, R50])
    def test_reproduces_density_knudsen_and_mach_columns(self, name):
        r = ref(name)
        atm = atmosphere.ReplayAtmosphere(r)
        for i in range(0, len(r.time), 25):
            s = atm.state(r.time[i], r.altitude[i], r.lat[i], r.lon[i])
            assert s.rho == pytest.approx(r.density[i], rel=1e-6)
            if r.knudsen[i] >= 1e-3:                     # SESAM prints knudsen with 5 decimals: below this it's <3 sig figs
                assert aero.knudsen(s.rho, s.m_bar, r.diameter) == pytest.approx(r.knudsen[i], rel=0.02)
            if r.mach[i] > 0.3:
                assert aero.mach(r.velocity[i], s.T, s.m_bar) == pytest.approx(r.mach[i], rel=0.01)


def test_vacuum():
    s = atmosphere.VacuumAtmosphere().state(0.0, 1e5, LAT, LON)
    assert s.rho == 0.0 and s.p == 0.0 and s.T == 200.0 and np.all(s.wind_enu == 0.0)

"""fap.py: DRAMA fap_day.dat -> NRLMSISE-00 solar/geomagnetic inputs."""
from datetime import date

import pytest

from reentry_model import fap


def test_parses_the_repo_file_and_the_reference_epoch():
    records = fap.load_fap_day(fap.DEFAULT_FAP_DAY)
    assert records[date(2024, 8, 1)].f10 == 234.0
    assert records[date(2024, 8, 1)].f3m == 194.0
    assert records[date(2024, 8, 1)].ap == 19.0
    assert records[date(2024, 7, 31)].f10 == 234.0 and records[date(2024, 8, 2)].f10 == 246.0


def test_solar_indices_use_previous_day_f107_and_same_day_average_and_ap():
    records = fap.load_fap_day(fap.DEFAULT_FAP_DAY)
    s = fap.solar_indices(records, date(2024, 8, 2))
    assert s.f107 == 234.0          # F10.7 of 2024-08-01 (previous day)
    assert s.f107a == 194.0         # 81-day mean column of 2024-08-02
    assert s.ap == 8.0              # daily Ap of 2024-08-02


def test_missing_day_raises(tmp_path):
    p = tmp_path / "fap_day.dat"
    p.write_text("#d/mm/yyyy F10 F3M SSN Ap  3-hr Kp-Indices\n01/08/2024 234 194 259 019 5+4+3+2+3o2+2+2o\n")
    records = fap.load_fap_day(str(p))
    assert list(records) == [date(2024, 8, 1)]
    with pytest.raises(KeyError):
        fap.solar_indices(records, date(2024, 8, 1))      # needs 2024-07-31 for f107

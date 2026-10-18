"""DRAMA space-weather file (fap_day.dat) -> the inputs NRLMSISE-00 wants.

File format (one line per day): `dd/mm/yyyy F10 F3M SSN Ap kp-string`, `#` comment lines.
F10 = daily F10.7, F3M = its 81-day (3-month) mean, Ap = daily Ap.
"""
import os
from dataclasses import dataclass
from datetime import date, datetime, timedelta

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DEFAULT_FAP_DAY = os.path.join(REPO_ROOT, "data", "fap_day.dat")


@dataclass(frozen=True)
class FapRecord:
    f10: float
    f3m: float
    ssn: int
    ap: float


@dataclass(frozen=True)
class SolarIndices:
    f107: float      # F10.7 of the previous day (NRLMSISE-00 convention)
    f107a: float     # 81-day running mean
    ap: float        # daily Ap


def load_fap_day(path):
    """{date: FapRecord} for every data line of a fap_day.dat file."""
    records = {}
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.split()
            day = datetime.strptime(parts[0], "%d/%m/%Y").date()
            records[day] = FapRecord(float(parts[1]), float(parts[2]), int(parts[3]), float(parts[4]))
    return records


def solar_indices(records, day):
    """NRLMSISE-00 inputs for `day`; KeyError if the day or the previous day is missing."""
    previous = records[day - timedelta(days=1)]
    today = records[day]
    return SolarIndices(f107=previous.f10, f107a=today.f3m, ap=today.ap)

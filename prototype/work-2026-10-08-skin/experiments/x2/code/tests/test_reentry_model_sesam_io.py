"""sesam_io.py: a wrapper run (CSV + JSON) as a Reference."""
import math
import os
from datetime import datetime

import numpy as np
import pytest

from reentry_model import sesam_io

R50 = "sphere_d050.00mm_T0300.0K_v07.50000kms_h115.000km_mAA7075_nomelt_msis_nowind"
R100 = "sphere_d100.00mm_T0300.0K_v07.50000kms_h077.500km_mAA7075_nomelt_msis"


def test_loads_r50_with_si_units_and_initial_state():
    ref = sesam_io.load_reference(os.path.join(sesam_io.REFERENCE_DIR, R50 + ".csv"))
    assert ref.name == R50 and ref.use_wind is False and ref.atmosphere == "nrlmsise"
    assert len(ref.time) > 500 and ref.time[0] == 0.0 and np.all(np.diff(ref.time) > 0)
    assert ref.altitude[0] == pytest.approx(115000.0, abs=1.0) and ref.altitude[-1] == 0.0
    assert ref.velocity[0] == pytest.approx(7500.0) and ref.diameter == 0.05
    assert ref.material_density == 2813.0
    assert ref.initial.velocity == 7500.0 and ref.initial.altitude == 115000.0
    assert ref.initial.flight_path == pytest.approx(math.radians(-0.959331))
    assert ref.initial.heading == pytest.approx(math.radians(347.168296))
    assert ref.initial.lat == pytest.approx(math.radians(29.546067))
    assert ref.initial.epoch == datetime(2024, 8, 1, 12, 53, 7)
    assert ref.knudsen[0] == pytest.approx(40.95647) and ref.drag[0] == pytest.approx(2.059)
    assert len(ref.sha256) == 64


def test_json_path_defaults_to_the_csv_stem_and_must_exist(tmp_path):
    ref = sesam_io.load_reference(os.path.join(sesam_io.REFERENCE_DIR, R100 + ".csv"))
    assert ref.json_path.endswith(R100 + ".json") and ref.use_wind is True
    orphan = tmp_path / "x.csv"
    orphan.write_text("time_s,altitude_km\n0,1\n")
    with pytest.raises(FileNotFoundError):
        sesam_io.load_reference(str(orphan))

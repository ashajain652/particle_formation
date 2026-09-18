"""compare.py: model-vs-SESAM metrics and plots on a synthetic history with known offsets."""
import math
import os

import numpy as np
import pytest

from reentry_model import compare, sesam_io
from reentry_model import trajectory as tj

R100 = "sphere_d100.00mm_T0300.0K_v07.50000kms_h077.500km_mAA7075_nomelt_msis_nowind"


def synthetic_history(ref, dv=5.0, dh=-20.0, truncate=None):
    n = len(ref.time) if truncate is None else truncate
    cols = {k: np.zeros(n) for k in tj.CSV_COLUMNS}
    cols["time_s"] = ref.time[:n].copy()
    cols["velocity_kms"] = (ref.velocity[:n] + dv) / 1e3
    cols["altitude_km"] = (ref.altitude[:n] + dh) / 1e3
    cols["flight_path_deg"] = np.degrees(ref.flight_path[:n])
    cols["heading_deg"] = np.degrees(ref.heading[:n])
    cols["lat_deg"] = np.degrees(ref.lat[:n])
    cols["lon_deg"] = np.degrees(ref.lon[:n])
    cols["knudsen"] = ref.knudsen[:n]
    cols["drag"] = ref.drag[:n]
    cols["mach"] = ref.mach[:n]
    cols["density_kgm3"] = ref.density[:n]
    return tj.History(cols, np.zeros((n, 6)), "ground", {})


@pytest.fixture(scope="module")
def ref():
    return sesam_io.load_reference(os.path.join(sesam_io.REFERENCE_DIR, R100 + ".csv"))


def test_metrics_recover_known_offsets(ref):
    m = compare.metrics(synthetic_history(ref), ref)
    assert m["n_points"] == len(ref.time)
    assert m["all"]["dV_max_ms"] == pytest.approx(5.0, abs=1e-6) and m["all"]["dV_rms_ms"] == pytest.approx(5.0, abs=1e-6)
    assert m["all"]["dh_max_m"] == pytest.approx(20.0, abs=1e-6)
    assert m["hypersonic"]["n_points"] < m["all"]["n_points"]
    assert m["hypersonic"]["dV_rel_max"] == pytest.approx(5.0 / 1000.0, abs=3e-4)   # worst at the last V_ref > 1 km/s
    assert m["d_end_time_s"] == 0.0 and m["final_velocity_model_ms"] == pytest.approx(ref.velocity[-1] + 5.0)
    assert m["d_final_velocity_ms"] == pytest.approx(5.0)


def test_alignment_uses_only_reference_times_inside_the_model_flight(ref):
    short = synthetic_history(ref, truncate=100)
    a = compare.align(short, ref)
    assert a["t"][-1] <= short.columns["time_s"][-1] and len(a["t"]) == 100
    m = compare.metrics(short, ref)
    assert m["d_end_time_s"] < 0.0


def test_phase_metrics_are_nan_when_the_phase_is_empty(ref):
    a = compare.align(synthetic_history(ref), ref)
    phase = compare._phase(a, np.zeros(len(a["t"]), dtype=bool))
    assert phase["n_points"] == 0
    assert math.isnan(phase["dV_max_ms"])
    assert math.isnan(phase["dV_rms_ms"])
    assert math.isnan(phase["dV_rel_max"])
    assert math.isnan(phase["dV_rel_rms"])
    assert math.isnan(phase["dh_max_m"])
    assert math.isnan(phase["dh_rms_m"])


def test_plots_are_written(ref, tmp_path):
    paths = compare.plot_all(synthetic_history(ref), ref, str(tmp_path), "synthetic")
    assert [os.path.basename(p) for p in paths] == list(compare.PLOT_NAMES)
    for p in paths:
        assert os.path.getsize(p) > 5000

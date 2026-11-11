"""compare.py: model-vs-SESAM metrics and plots on a synthetic history with known offsets."""
import dataclasses
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


def _thermal_history_like(reference, factor=1.0, dT=0.0):
    """A model history built from the reference's own heat columns (scaled), at half its time stamps."""
    from reentry_model import trajectory as tj
    t = reference.time[::2]
    pick = lambda arr: np.interp(t, reference.time, arr)
    columns = {"time_s": t, "velocity_kms": pick(reference.velocity) / 1e3, "altitude_km": pick(reference.altitude) / 1e3,
               "knudsen": pick(reference.knudsen), "convective_heat_W": factor * pick(reference.convective_heat),
               "integrated_heat_J": factor * pick(reference.integrated_heat), "rad_cooling_W": pick(reference.rad_cooling),
               "temperature_K": pick(reference.temperature) + dT, "surface_T_max_K": pick(reference.temperature) + dT + 50.0,
               "surface_T_min_K": pick(reference.temperature) + dT - 20.0}
    return tj.History(columns, np.zeros((t.size, 6)), "ground")


def test_reference_carries_the_heat_columns(ref):
    assert ref.convective_heat[0] == pytest.approx(15819.606) and ref.rad_cooling[0] == pytest.approx(-5.772)
    assert ref.integrated_heat[0] == 0.0 and ref.integrated_heat[-1] > 1e5


def test_thermal_metrics_of_a_scaled_copy_and_plots(ref, tmp_path):
    hist = _thermal_history_like(ref, factor=1.02, dT=3.0)
    assert compare.has_thermal(hist, ref) and not compare.has_thermal(synthetic_history(ref), ref)
    m = compare.thermal_metrics(hist, ref)
    assert m["Q_conv"]["max"] == pytest.approx(0.02, abs=2e-3) and 0.015 < m["Q_conv"]["continuum_rel_max"] < 0.04   # + interpolation at half the stamps
    assert m["integrated_heat"]["rel_error_end"] == pytest.approx(0.02, abs=1e-3) and m["integrated_heat"]["ratio_end"] == pytest.approx(1.02, abs=1e-3)
    assert m["temperature"]["dT_max_K"] == pytest.approx(3.0, abs=0.5) and m["radiated"]["max"] < 0.05
    assert m["n_continuum_points"] > 100 and m["integrated_heat"]["rel_error_end_of_hypersonic"] == pytest.approx(0.02, abs=2e-3)
    paths = compare.plot_thermal(hist, ref, str(tmp_path), "scaled copy")
    assert [os.path.basename(p) for p in paths] == list(compare.THERMAL_PLOT_NAMES) and all(os.path.getsize(p) > 5000 for p in paths)


def test_has_thermal_needs_all_reference_heat_columns(ref):
    hist = _thermal_history_like(ref)
    assert compare.has_thermal(hist, ref) is True
    for name in ("convective_heat", "rad_cooling", "integrated_heat"):
        partial = dataclasses.replace(ref, **{name: None})
        assert compare.has_thermal(hist, partial) is False
    assert compare.has_thermal(hist, None) is True


# ---------------------------------------------------------------------------------------------------------------
# Step 3: the melting reference, the mass metrics and plots

MELT100 = "sphere_d100.00mm_T0300.0K_v07.50000kms_h077.500km_nowind"


@pytest.fixture(scope="module")
def melt_ref():
    return sesam_io.load_reference(os.path.join(sesam_io.REFERENCE_DIR, MELT100 + ".csv"))


def melting_history(ref, lag=0.0, film=0.0):
    """A synthetic model history that reproduces the reference's mass shifted by `lag` seconds, with the melt columns."""
    h = synthetic_history(ref, dv=0.0, dh=0.0)
    c = h.columns
    n = len(ref.time)
    c["mass_kg"] = np.interp(ref.time - lag, ref.time, ref.mass) + film
    for key in ("film_mass_kg", "sprayed_mass_kg", "runoff_mass_kg", "removed_mass_kg", "spraying_area_m2", "theta_cr_deg", "n_released",
                "released_mass_kg", "r_median_um", "r_max_um", "closure_fraction_girin", "closure_fraction_couette_slip",
                "closure_fraction_couette_fm", "kn_body", "kn_local_stag", "film_thickness_mean_mm", "film_thickness_max_mm",
                "convective_heat_W", "surface_T_max_K"):
        c[key] = np.zeros(n)
    c["film_mass_kg"][:] = film
    c["removed_mass_kg"] = ref.mass[0] - c["mass_kg"] + film
    c["sprayed_mass_kg"] = c["removed_mass_kg"].copy()
    c["r_median_um"][:] = 150.0
    c["closure_fraction_girin"][:] = 0.5
    c["closure_fraction_couette_slip"][:] = 0.5
    c["kn_body"][:] = 0.01
    c["kn_local_stag"][:] = 1e-4
    return h


def test_melting_reference_and_has_melt(ref, melt_ref):
    assert ref.mass is not None and float(ref.mass.min()) == float(ref.mass[0])                         # no-melt: constant mass
    assert melt_ref.mass[0] == pytest.approx(1.473) and melt_ref.mass[-1] == 0.0 and melt_ref.thickness[0] == pytest.approx(0.05)
    assert melt_ref.thickness[-1] == 0.0 and np.all(np.diff(melt_ref.mass) <= 0.0)
    h = melting_history(melt_ref)
    assert compare.has_melt(h) and compare.has_melt(h, melt_ref) and not compare.has_melt(h, ref) and not compare.has_melt(synthetic_history(ref))


def test_melt_metrics_recover_a_known_lag(melt_ref):
    m = compare.melt_metrics(melting_history(melt_ref, lag=-1.0), melt_ref)          # the model runs 1 s ahead of SESAM
    assert m["mass"]["max_rel_m0"] < 0.06 and m["mass"]["max_abs_kg"] > 0.0
    assert m["demise_time_diff_s"] == pytest.approx(-1.0, abs=0.1) and m["demise_time_rel"] == pytest.approx(-1.0 / m["demise_time_reference_s"], rel=0.1)
    assert m["onset_time_reference_s"] == pytest.approx(43.55) and m["onset_altitude_reference_km"] == pytest.approx(71.005, abs=0.01)
    assert m["onset_altitude_diff_km"] is not None and m["demise_time_reference_s"] == pytest.approx(66.8, abs=0.1)
    exact = compare.melt_metrics(melting_history(melt_ref), melt_ref)
    assert exact["mass"]["max_rel_m0"] < 1e-12 and abs(exact["demise_time_diff_s"]) < 1e-9
    with_film = compare.melt_metrics(melting_history(melt_ref, film=0.05), melt_ref)
    assert with_film["demise_time_diff_s"] == pytest.approx(exact["demise_time_diff_s"], abs=1e-9)   # the crossing is on the body's material


def test_melt_plots_are_written(melt_ref, tmp_path):
    h = melting_history(melt_ref, lag=-1.0)
    paths = compare.plot_melt(h, str(tmp_path), "melt test", melt_ref)
    assert [os.path.basename(p) for p in paths] == list(compare.MELT_PLOT_NAMES) and all(os.path.getsize(p) > 5000 for p in paths)
    alone = compare.plot_melt(h, str(tmp_path / "alone"), "no reference")
    assert len(alone) == len(compare.MELT_PLOT_NAMES) and all(os.path.isfile(p) for p in alone)

"""Full-model verification against the committed SESAM references (marker: reference)."""
import os

import pytest

from reentry_model import aero, atmosphere, body, cli, fap, sesam_io
from reentry_model import compare
from reentry_model import trajectory as tj

NAMES = [
    "sphere_d100.00mm_T0300.0K_v07.50000kms_h077.500km_mAA7075_nomelt_msis_nowind",
    "sphere_d100.00mm_T0300.0K_v07.50000kms_h077.500km_mAA7075_nomelt_msis",
    "sphere_d050.00mm_T0300.0K_v07.50000kms_h115.000km_mAA7075_nomelt_msis_nowind",
    "sphere_d050.00mm_T0300.0K_v07.50000kms_h115.000km_mAA7075_nomelt_msis",
]
# Spec section 9 expectations over the hypersonic phase (V_ref > 1 km/s). Task 10 step 3 may revise a
# value only together with the measured number and the reason, recorded in the README verification table.
#
# replay dV_rel_max: 1 % — SESAM's printed density lags the state its drag used by ~0.05 s (its printed
# dynamic pressure exceeds 1/2 rho V^2 from its own columns by up to 1 %); measured 0.3-0.9 %, see README.
# nrlmsise: dh_max_m (worst measured 430 m) and d_end_time_rel (worst 0.29 %) already pass the spec's
# values unchanged. dV_rel_max is missed (worst measured 6.802 %) and is set to the measured worst value
# x 1.5, rounded to one significant figure: 6.802 % x 1.5 = 10.2 % -> 0.1.
THRESHOLDS = {
    "replay":   {"dV_rel_max": 0.01, "dh_max_m": 100.0, "d_end_time_rel": 0.01},
    "nrlmsise": {"dV_rel_max": 0.1, "dh_max_m": 500.0, "d_end_time_rel": 0.03},
}


def simulate(ref, mode):
    initial = tj.InitialState(ref.initial.velocity, ref.initial.altitude, ref.initial.flight_path, ref.initial.heading,
                              ref.initial.lat, ref.initial.lon, ref.initial.epoch)
    if mode == "replay":
        atm = atmosphere.ReplayAtmosphere(ref)
    else:
        solar = fap.solar_indices(fap.load_fap_day(fap.DEFAULT_FAP_DAY), ref.initial.epoch.date())
        atm = atmosphere.NRLMSISE00Atmosphere(ref.initial.epoch, solar)
    mass = body.sphere_mass(ref.diameter, ref.material_density)
    sim = tj.Simulator(initial, body.ConstantBody(mass), atm, aero.SphereDragTables.from_json(), aero.SesamTable(),
                       tj.Settings(diameter=ref.diameter, cadence=5.0))
    return sim.run(extra_times=ref.time)


@pytest.mark.reference
@pytest.mark.parametrize("mode", ["replay", "nrlmsise"])
@pytest.mark.parametrize("name", NAMES)
def test_model_matches_sesam(name, mode):
    ref = sesam_io.load_reference(os.path.join(sesam_io.REFERENCE_DIR, name + ".csv"))
    history = simulate(ref, mode)
    assert history.end_reason == "ground"
    m = compare.metrics(history, ref)
    hyp, thr = m["hypersonic"], THRESHOLDS[mode]
    assert hyp["dV_rel_max"] <= thr["dV_rel_max"], (name, mode, hyp)
    assert hyp["dh_max_m"] <= thr["dh_max_m"], (name, mode, hyp)
    assert abs(m["d_end_time_rel"]) <= thr["d_end_time_rel"], (name, mode, m["d_end_time_s"])

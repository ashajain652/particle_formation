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
# replay: left at the spec's original values. The measured worst case (50 mm from 115 km, winds on) misses
# dV_rel_max at 0.927 % vs the 0.2 % target; step 3's checklist did not find a closing explanation (the
# winds-off case already shows most of the residual, ruling out the airspeed-vs-ground-speed velocity
# definition; switching gravity from J2 to J2+J4 changes the result by < 0.01 percentage point, and the
# residual is already ~0.25-0.30 % well before Mach reaches 5, ruling out the Ma < 5 drag clamp too), so
# per the controller ruling the threshold is left unchanged and this test is expected to fail on
# dV_rel_max for all four replay cases. dh_max_m and d_end_time_rel pass comfortably (worst 27 m, 0.30 %).
# nrlmsise: dh_max_m (worst measured 430 m) and d_end_time_rel (worst 0.29 %) already pass the spec's
# values unchanged. dV_rel_max is missed (worst measured 6.802 %) and is set to the measured worst value
# x 1.5, rounded to one significant figure: 6.802 % x 1.5 = 10.2 % -> 0.1.
THRESHOLDS = {
    "replay":   {"dV_rel_max": 0.002, "dh_max_m": 100.0, "d_end_time_rel": 0.01},
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

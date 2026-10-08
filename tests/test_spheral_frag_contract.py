"""spheral_frag.contract: the frame contract is one consistent table (M1 plan, Task 1; review focus 1)."""
import pytest

from spheral_frag import contract
from spheral_frag.contract import FE_FIELDS


def test_keys_and_names_unique():
    keys = [f.key for f in FE_FIELDS]
    assert len(keys) == len(set(keys))
    names = [(f.location, f.name) for f in FE_FIELDS]
    assert len(names) == len(set(names))
    # also unique across locations, so a name alone identifies the item in messages and the prepare JSON
    assert len({f.name for f in FE_FIELDS}) == len(FE_FIELDS)


def test_locations_and_flags_are_well_formed():
    for f in FE_FIELDS:
        assert f.location in contract.LOCATIONS, f
        assert f.units and f.consumer and f.name, f
        assert isinstance(f.required, bool) and isinstance(f.nan_allowed, bool) and isinstance(f.confirmed, bool)


def test_every_required_patch_field_has_a_consumer():
    req = contract.required("patch")
    assert {f.key for f in req} >= {"faces", "q_conv", "film_thickness", "film_T", "p_w", "tau", "release_rate",
                                     "delta_m"}
    for f in req:
        assert f.consumer.strip() and f.consumer != "informational", f


def test_fe_name_round_trip():
    for f in FE_FIELDS:
        assert contract.fe_name(f.key) == f.name
        assert contract.key_of(f.location, f.name) == f.key
        assert contract.field_spec(f.key) is f
    with pytest.raises(KeyError):
        contract.fe_name("no_such_item")
    with pytest.raises(ValueError):
        contract.required("everywhere")


def test_units_of_the_per_step_fields():
    """Review focus 3: release_rate per macro step, deep_thickness a mass per area, p_w = 0 means not evaluated."""
    assert "per macro step" in contract.field_spec("release_rate").units
    assert "mass per area" in contract.field_spec("deep_thickness").units
    assert "0 = not evaluated" in contract.field_spec("p_w").units
    assert {f.key for f in FE_FIELDS if f.nan_allowed} == {"delta_m", "r_droplet", "kn_local"}


def test_history_and_run_items():
    hist = [f.name for f in contract.required("history")]
    assert hist[:10] == ["time_s", "altitude_km", "velocity_kms", "mass_kg", "flight_path_deg", "lat_deg", "lon_deg",
                         "density_kgm3", "dynamic_pressure_Pa", "load_factor_g"]
    assert set(hist[10:]) == {"p_w_stag_Pa", "film_mass_kg", "deep_mass_kg", "sprayed_mass_kg"}
    run = {f.name for f in contract.required("run")}
    assert {"run_name", "inputs.diameter_mm", "inputs.mass_kg", "settings.atmosphere", "settings.wind",
            "settings.macro_step_s", "settings.frames_every", "settings.seed", "v_hat"} <= run
    assert contract.V_HAT == (1.0, 0.0, 0.0)


def test_nothing_confirmed_before_task_10():
    """`confirmed` is set only by checking the real flight's frames (Task 10)."""
    assert contract.unconfirmed() == [f.key for f in FE_FIELDS]


def test_prepared_frame_schema_covers_the_patch_fields():
    arrays = contract.PREPARED_FRAME_ARRAYS
    for f in contract.data_fields("patch"):
        assert "patch_" + f.key in arrays
    for name, (dtype, shape, units) in arrays.items():
        assert dtype in ("f8", "i8", "i4"), name
        assert isinstance(shape, tuple) and all(s in ("n", "m", "f") or isinstance(s, int) for s in shape), name
        assert units, name
    assert arrays["points"][1] == ("n", 3) and arrays["tets"] == ("i4", ("m", 4), arrays["tets"][2])
    assert contract.PREPARED_FRAME_OPTIONAL <= set(arrays)
    assert "patch_p_w" not in contract.PREPARED_FRAME_OPTIONAL and "patch_q_rad" in contract.PREPARED_FRAME_OPTIONAL

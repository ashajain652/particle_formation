"""The package's bundled data and the committed SESAM reference runs."""
import csv
import json
import os

import pytest

import reentry_model
from reentry_model import constants as c

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REFERENCE_DIR = os.path.join(REPO_ROOT, "data", "reference_runs")

REFERENCE_NAMES = [
    "sphere_d100.00mm_T0300.0K_v07.50000kms_h077.500km_mAA7075_nomelt_msis",
    "sphere_d100.00mm_T0300.0K_v07.50000kms_h077.500km_mAA7075_nomelt_msis_nowind",
    "sphere_d050.00mm_T0300.0K_v07.50000kms_h115.000km_mAA7075_nomelt_msis",
    "sphere_d050.00mm_T0300.0K_v07.50000kms_h115.000km_mAA7075_nomelt_msis_nowind",
]


def test_atdb_tables_have_the_six_mach_points():
    with open(os.path.join(reentry_model.DATA_DIR, "atdb_sphere.json")) as fh:
        atdb = json.load(fh)
    assert atdb["mach"] == [5.0, 10.0, 15.0, 20.0, 25.0, 30.0]
    assert atdb["cd_free_molecular"][0] == 2.360635 and atdb["cd_continuum"][-1] == 0.913599
    assert len(atdb["heat_flux_factor_continuum"]) == 6


def test_us76_table_is_the_standard_atmosphere():
    path = os.path.join(reentry_model.DATA_DIR, "us76_static_environment.csv")
    rows = [r for r in csv.reader(open(path)) if r and not r[0].startswith("#")]
    assert len(rows) == 1501                      # 0..150 km every 100 m
    first = [float(v) for v in rows[0]]
    assert first[0] == 0.0 and abs(first[1] - 101325.0) < 1.0 and abs(first[2] - 1.225) < 1e-3
    at70 = [float(v) for v in rows[700]]
    assert at70[0] == 70000.0 and abs(at70[2] - 8.283e-5) / 8.283e-5 < 2e-3   # US76: 8.283e-5 kg/m3


def test_reference_runs_are_committed_and_consistent():
    for name in REFERENCE_NAMES:
        csv_path = os.path.join(REFERENCE_DIR, name + ".csv")
        json_path = os.path.join(REFERENCE_DIR, name + ".json")
        assert os.path.isfile(csv_path) and os.path.isfile(json_path), name
        doc = json.load(open(json_path))
        assert doc["run_name"] == name and doc["status"] == "ok"
        assert doc["inputs"]["material"] == "AA7075_nomelt" and doc["inputs"]["atmosphere"] == "nrlmsise"
        rows = list(csv.DictReader(open(csv_path)))
        assert len(rows) > 100 and float(rows[-1]["altitude_km"]) == 0.0
        assert abs(float(rows[-1]["mass_kg"]) - doc["inputs"]["initial_mass_kg"]) < 1e-3   # no mass loss


def test_constants():
    assert c.WGS84_A == 6378137.0 and abs(c.WGS84_F - 1 / 298.257223563) < 1e-15
    assert c.MU_EARTH == 3.986004418e14 and c.OMEGA_EARTH == 7.2921159e-5
    assert abs(c.M_BAR_AIR / c.ATOMIC_MASS_UNIT - 28.9644) < 1e-6
    assert c.HARD_SPHERE_DIAMETER == 3.65e-10 and c.GAMMA_AIR == 1.4


MELTING_REFERENCE_NAMES = ["sphere_d100.00mm_T0300.0K_v07.50000kms_h077.500km_nowind",
                           "sphere_d050.00mm_T0300.0K_v07.50000kms_h115.000km_nowind"]


def test_melting_reference_runs_are_committed_and_consistent():
    """The two Step 3 melting references: DRAMA's drama-AA7075 (850 K, 400 kJ/kg) on the US76 table, winds off; the
    sphere melts completely (mass 0 at the end, 'demised')."""
    for name, onset_km, end_km in zip(MELTING_REFERENCE_NAMES, (71.005, 77.104), (66.476, 73.094)):
        csv_path = os.path.join(REFERENCE_DIR, name + ".csv")
        doc = json.load(open(os.path.join(REFERENCE_DIR, name + ".json")))
        assert doc["run_name"] == name and doc["status"] == "ok" and doc["inputs"]["material"] == "drama-AA7075"
        assert doc["inputs"]["atmosphere"] == "static" and doc["inputs"]["use_wind"] is False
        assert doc["results"]["melting_temperature_K"] == 850.0 and abs(doc["results"]["altitude_first_melt_km"] - onset_km) < 1e-3
        assert abs(doc["results"]["altitude_last_melt_km"] - end_km) < 1e-3
        rows = list(csv.DictReader(open(csv_path)))
        assert float(rows[-1]["mass_kg"]) == 0.0 and float(rows[-1]["thick_mm"]) == 0.0 and float(rows[-1]["altitude_km"]) > 60.0
        assert abs(float(rows[0]["mass_kg"]) - doc["inputs"]["initial_mass_kg"]) < 1e-3


def test_disc_aerothermal_table():
    """data/atdb_disc.json: ATDB_CYLINDER at zero angle of attack and the thinnest aspect ratio -- the flat-face limit
    of the sphere-to-disc family (Step 3). HTG's continuum database is modified Newtonian, so the sphere entry is
    0.49897 of the disc entry at every Mach number, and the free-molecular entries differ by only 3 %."""
    sphere = json.load(open(os.path.join(REPO_ROOT, "reentry_model", "data", "atdb_sphere.json")))
    disc = json.load(open(os.path.join(REPO_ROOT, "reentry_model", "data", "atdb_disc.json")))
    assert disc["mach"] == sphere["mach"] == [5.0, 10.0, 15.0, 20.0, 25.0, 30.0]
    assert "ATDB_CYLINDER" in disc["_provenance"] and "angle of attack 0" in disc["_provenance"]
    cd_c, cd_fm = disc["cd_continuum"], disc["cd_free_molecular"]
    assert cd_c[1] == pytest.approx(1.824148) and cd_fm[1] == pytest.approx(2.202072)
    assert all(abs(s / d - 0.49897) < 1e-5 for s, d in zip(sphere["cd_continuum"], cd_c))
    assert all(0.99 < d / s < 1.08 for s, d in zip(sphere["cd_free_molecular"], cd_fm))
    assert all(v == pytest.approx(0.329561) for v in disc["heat_flux_factor_continuum"])       # Mach-independent, like the sphere's

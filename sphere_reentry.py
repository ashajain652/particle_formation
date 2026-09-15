#!/usr/bin/env python3
"""
sphere_reentry.py -- one solid AA7075 sphere through ESA DRAMA / SARA (SESAM)

Runs SESAM through the pyDRAMA package for a single sphere fragment starting
from a geodetic state (altitude, latitude, longitude, velocity, flight-path
angle, heading) and writes

    <outdir>/<run_name>.csv   full-resolution history (trajectory + aerothermal, joined on time)
    <outdir>/<run_name>.json  inputs, statistics, warnings, provenance

Usage
-----
    python sphere_reentry.py --velocity 7.5 --altitude 77.5 --temperature 300 --diameter 50
    python sphere_reentry.py ... --flight-path-angle -0.96 --heading 347.2 --lat 29.5 --lon -82.1 \
                                 --epoch 2024-08-01T12:53:07
    python sphere_reentry.py ... --dry-run          # print the SESAM config and run name, run nothing

Units: velocity km/s, altitude km, temperature K, diameter mm, angles deg.
Requires DRAMA 4.1.4 and pyDRAMA (conda env drama_env); see README.md.
Design: docs/superpowers/specs/2026-09-13-sphere-reentry-sweep-design.md
"""
from __future__ import annotations

import argparse
import csv
import glob
import json
import math
import os
import re
import shutil
import socket
import sys
import time
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from datetime import datetime, timezone

SCRIPT_VERSION = "1.0.0"
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DRAMA_INSTALL_PATH = "/Applications/DRAMA-4.1.4"

# --- material: drama-AA7075 exactly as shipped in DRAMA 4.1.4 (the default) ----
MATERIAL_NAME = "drama-AA7075"
RHO_AA7075 = 2813.0            # [kg/m3]
T_MELT_AA7075 = 850.0          # [K]
MELT_TOLERANCE_K = 0.5         # |T - T_melt| <= tol counts as "at melting temperature"
# DRAMA's own material database; --material NAME is looked up here (metals only)
MATERIAL_DB_PATH = os.path.join(DRAMA_INSTALL_PATH, "TOOLS", "material_database.xml")

# --- run control ------------------------------------------------------------------
ENERGY_THRESHOLD_J = 1e-9      # SESAM drops fragments below this kinetic energy; parent run used 15 J
CSV_MASS_RESOLUTION_KG = 1e-3  # SESAM prints mass with 3 decimals
COARSE_MASS_LIMIT_KG = 0.05    # warn when the initial mass is below this (CSV mass column too coarse)
DEMISE_MASS_FRACTION = 0.05    # final/initial mass below this -> outcome "demised"
OBJECT_NAME = "sphere"
OBJECT_UUID = "5b3e0000-0000-4000-8000-000000000001"
PARENT_EPOCH = datetime(2024, 8, 1, 12, 0, 0)   # epoch of the parent satellite run
DEFAULT_TIMEOUT_S = 600
DEFAULT_OUTDIR = os.path.join(SCRIPT_DIR, "sphere_sweep_output", "runs")
DEFAULT_RAW_DIR = os.path.join(SCRIPT_DIR, "sphere_sweep_output", "raw")
DEFAULT_FAP_DAY = os.path.join(SCRIPT_DIR, "data", "fap_day.dat")
DEFAULT_FAP_MON = os.path.join(SCRIPT_DIR, "data", "fap_mon.dat")


# =============================================================================
# 1. Inputs
# =============================================================================

class MaterialError(ValueError):
    """A material name or definition that cannot be used."""


@dataclass(frozen=True)
class Material:
    """The sphere's material: what SESAM is told, plus the two numbers this script needs itself.

    density_kgm3 sets the sphere's mass and final radius; melting_temperature_K sets the
    "time at melting temperature" statistic. For a built-in DRAMA material `definition`
    is None (SESAM resolves the name from its own materials.xml); for a custom material
    it holds the full DRAMA material dict that is injected as materialList.
    """
    name: str
    density_kgm3: float
    melting_temperature_K: float
    source: str = "builtin"          # "builtin" | "custom_file"
    definition: dict = None


DEFAULT_MATERIAL = Material(MATERIAL_NAME, RHO_AA7075, T_MELT_AA7075)


def list_builtin_materials(db_path: str = MATERIAL_DB_PATH) -> dict:
    """{name: Material} for every <metalMaterial> in DRAMA's material database."""
    if not os.path.isfile(db_path):
        raise MaterialError("DRAMA material database not found: {}".format(db_path))
    try:
        root = ET.parse(db_path).getroot()
    except ET.ParseError as exc:
        raise MaterialError("cannot parse DRAMA material database {}: {}".format(db_path, exc))
    materials = {}
    for el in root.findall("metalMaterial"):
        name = (el.findtext("name") or "").strip()
        density = el.findtext("density")
        melting = el.findtext("meltingTemperature")
        if not name or density is None or melting is None:
            continue
        materials[name] = Material(name, float(density), float(melting))
    return materials


def load_builtin_material(name: str, db_path: str = MATERIAL_DB_PATH) -> Material:
    """The built-in DRAMA metal called `name`; MaterialError listing the valid names otherwise."""
    materials = list_builtin_materials(db_path)
    if name not in materials:
        raise MaterialError("unknown built-in material {!r}; DRAMA's metals are: {}".format(
            name, ", ".join(sorted(materials))))
    return materials[name]


# the required keys of a metal in pyDRAMA's sara_materials.schema.json
CUSTOM_MATERIAL_REQUIRED = ("name", "density", "specificHeatCapacity", "meltingHeat",
                            "meltingTemperature", "emissivity", "heatConductivity",
                            "oxideActivationTemperature", "oxideEmissivity",
                            "oxideHeatOfFormation", "oxideReactionProbability")


def load_material_file(path: str) -> Material:
    """A custom metal from a JSON file in DRAMA's material format (see examples/)."""
    try:
        with open(path) as fh:
            data = json.load(fh)
    except OSError as exc:
        raise MaterialError("cannot read material file {}: {}".format(path, exc))
    except ValueError as exc:
        raise MaterialError("material file {} is not valid JSON: {}".format(path, exc))
    if not isinstance(data, dict):
        raise MaterialError("material file {} must contain one JSON object".format(path))
    missing = [k for k in CUSTOM_MATERIAL_REQUIRED if k not in data]
    if missing:
        raise MaterialError("material file {} is missing required fields: {}".format(
            path, ", ".join(missing)))
    definition = dict(data)
    definition.setdefault("materialType", "metal")
    if definition["materialType"] != "metal":
        raise MaterialError("material file {}: only materialType 'metal' is supported for a solid "
                            "sphere, got {!r}".format(path, definition["materialType"]))
    try:
        density = float(definition["density"])
        melting = float(definition["meltingTemperature"])
    except (TypeError, ValueError) as exc:
        raise MaterialError("material file {}: density/meltingTemperature must be numbers ({})".format(path, exc))
    if density <= 0.0:
        raise MaterialError("material file {}: density must be > 0, got {}".format(path, density))
    if melting <= 0.0:
        raise MaterialError("material file {}: meltingTemperature must be > 0, got {}".format(path, melting))
    return Material(str(definition["name"]), density, melting, source="custom_file", definition=definition)


def sphere_mass_kg(diameter_mm: float, density_kgm3: float = RHO_AA7075) -> float:
    """Mass of a solid sphere: rho * 4/3 * pi * r^3 (AA7075 density unless given)."""
    r = diameter_mm / 2000.0
    return density_kgm3 * 4.0 / 3.0 * math.pi * r ** 3


def sphere_cross_section_m2(diameter_mm: float) -> float:
    r = diameter_mm / 2000.0
    return math.pi * r * r


@dataclass(frozen=True)
class SphereRun:
    """All inputs of one sphere run (CLI units)."""
    velocity_kms: float
    altitude_km: float
    temperature_K: float
    diameter_mm: float
    flight_path_deg: float = 0.0
    heading_deg: float = 0.0
    lat_deg: float = 0.0
    lon_deg: float = 0.0
    epoch: datetime = PARENT_EPOCH
    material: Material = DEFAULT_MATERIAL

    @property
    def radius_m(self) -> float:
        return self.diameter_mm / 2000.0

    @property
    def mass_kg(self) -> float:
        return sphere_mass_kg(self.diameter_mm, self.material.density_kgm3)

    @property
    def cross_section_m2(self) -> float:
        return sphere_cross_section_m2(self.diameter_mm)


def material_slug(name: str) -> str:
    """Filename-safe form of a material name (letters, digits, '-' and '.' kept; the rest -> '_')."""
    return "".join(c if c.isalnum() or c in "-." else "_" for c in name)


def run_name(run: SphereRun) -> str:
    """Fixed-width name so that lexicographic order equals parameter order (spec 4.3).

    The default material keeps the original name (so existing sweeps stay resumable);
    any other material appends "_m<slug>" so two materials can never share output files.
    """
    name = "sphere_d{:06.2f}mm_T{:06.1f}K_v{:08.5f}kms_h{:07.3f}km".format(
        run.diameter_mm, run.temperature_K, run.velocity_kms, run.altitude_km)
    if run.material.name != DEFAULT_MATERIAL.name:
        name += "_m" + material_slug(run.material.name)
    return name


def read_lines(path: str) -> list:
    """Non-blank lines of a text file without trailing newlines (DRAMA fap_day/fap_mon content)."""
    with open(path) as fh:
        return [ln.rstrip("\n") for ln in fh if ln.strip()]


# =============================================================================
# 2. SESAM configuration (spec 4.4) -- mirrors the parent satellite run except energyThreshold
# =============================================================================

def build_config(run: SphereRun) -> dict:
    """The complete pyDRAMA SARA configuration for one sphere (pass it as config=[cfg])."""
    cfg = _build_config(run)
    if run.material.definition is not None:
        cfg["materialList"] = [json.loads(json.dumps(run.material.definition))]
    return cfg


def _build_config(run: SphereRun) -> dict:
    obj = {
        "name": OBJECT_NAME,
        "uniqueID": OBJECT_UUID,
        "primitive": {"sphere": {"radius": run.radius_m}},
        "mass": run.mass_kg,
        "material": run.material.name,
        "solid": True,
        "relativePosition": {"cartX": 0.0, "cartY": 0.0, "cartZ": 0.0,
                             "yaw": 0.0, "pitch": 0.0, "roll": 0.0},
        "scalingFactors": {"drag": 1.0, "lift": 1.0, "sideForce": 1.0,
                           "averageHeatFlux": 1.0, "averageHeatFluxCT": 1.0,
                           "averageHeatFluxTR": 1.0, "averageHeatFluxFM": 1.0},
        "attitude": "tumbling",
        "quantity": 1,
        "temperature": run.temperature_K,
    }
    return {
        # ---- general -------------------------------------------------------
        "runID": "SPHERE",
        "comment1": run_name(run),
        "comment2": "solid AA7075 sphere fragment",
        "beginDate": run.epoch,
        "runMode": "reentry-only",
        "monteCarlo": False,
        # ---- initial state: geodetic = altitude, lat, lon, velocity, flight path, heading
        "coordinateSystem": "geodetic",
        "initialDate": run.epoch,
        "element1": run.altitude_km,
        "element2": run.lat_deg,
        "element3": run.lon_deg,
        "element4": run.velocity_kms,
        "element5": run.flight_path_deg,
        "element6": run.heading_deg % 360.0,
        # ---- spacecraft-level aerodynamics / attitude ----------------------
        "assumedCrossSection": run.cross_section_m2,
        "dragCoefficient": 2.2,
        "reflectivityCoefficient": 1.3,
        "attitude": "tumbling",
        "fragmentsAttitudeAfterBreakup": "inherited",
        "globalSpacecraftTemperature": run.temperature_K,
        # ---- environment (same keys/values as the parent run) ---------------
        "densityScalingFactor": 1.0,
        "dynamicEnvironment": True,
        "useWind": True,
        "solarActivityFromFile": True,
        "useEnvironmentCSV": False,
        "ap": 8,
        "f107a": 170,
        # ---- numerics / output ----------------------------------------------
        "voxelatorMode": 1,
        "energyThreshold": ENERGY_THRESHOLD_J,   # parent used 15 J; see spec section 2
        "plotVisibilityMaps": False,
        "plotObjectTrajectories": False,
        "propagationWithOscar": False,
        # ---- the model: a built-in material is resolved from DRAMA's own materials.xml;
        #      a custom one is injected as materialList (which REPLACES that table for the run)
        "objects": [obj],
    }


# =============================================================================
# 3. SESAM output parsing (spec 4.6)
# =============================================================================

AERO_COLUMNS = ["time", "altitude", "temp", "mass", "thick", "convectiveHeat",
                "radiativeHeat", "oxidationHeat", "radCooling", "integratedHeat",
                "visibilityFactor"]
TRAJ_COLUMNS = ["time", "altitude", "lat", "lon", "velocity", "downrange", "drag",
                "lift", "side", "knudsen", "mach", "path", "heading", "density",
                "dynamicPressure", "loadFactor"]


def read_sara_table(path: str, columns: list) -> list:
    """Rows of a SESAM whitespace table; '#' lines and rows with the wrong width are skipped."""
    rows = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) != len(columns):
                continue
            try:
                rows.append({c: float(p) for c, p in zip(columns, parts)})
            except ValueError:
                continue
    return rows


_SESAM_VERSION_RE = re.compile(r"SESAM\s+([0-9][0-9.]*[0-9])")


def parse_sesam_version(path):
    """SESAM version from the '#  ---- DRAMA ( SESAM 2.3.0 ) ----' header line, or None."""
    if not path or not os.path.isfile(path):
        return None
    with open(path, errors="replace") as fh:
        for line in fh:
            if not line.startswith("#"):
                break
            m = _SESAM_VERSION_RE.search(line)
            if m:
                return m.group(1)
    return None


def parse_impacting_fragments(path):
    """First <fragment> of PySara.ImpactingFragments.xml (full-precision impact state), or None."""
    if not path or not os.path.isfile(path):
        return None
    try:
        root = ET.parse(path).getroot()
    except ET.ParseError:
        return None
    frag = root.find("fragment")
    if frag is None:
        return None

    def number(tag):
        el = frag.find(tag)
        return float(el.text) if el is not None and el.text else None

    velocity_ms = number("velocity")
    epoch = frag.find("epoch")
    return {
        "mass_kg": number("mass"),
        "velocity_kms": velocity_ms / 1000.0 if velocity_ms is not None else None,
        "lat_deg": number("latitude"),
        "lon_deg": number("longitude"),
        "epoch": epoch.text if epoch is not None else None,
    }


_EOL_RE = re.compile(r"reached end of life because (?:of )?(.+?)\.")
_EVENT_END_RE = re.compile(r"EVENT end \S+ ([0-9.eE+-]+)")
_EVENT_START_RE = re.compile(r"EVENT start \S+ ([0-9.eE+-]+)")


def parse_sesam_log(path):
    """End-of-life reason and the 6-decimal EVENT start/end masses from sesam.log."""
    info = {"end_of_life_reason": "unknown", "event_end_mass_kg": None, "event_start_mass_kg": None}
    if not path or not os.path.isfile(path):
        return info
    with open(path, errors="replace") as fh:
        text = fh.read()
    m = _EOL_RE.search(text)
    if m:
        info["end_of_life_reason"] = m.group(1).strip()
    m = _EVENT_END_RE.search(text)
    if m:
        info["event_end_mass_kg"] = float(m.group(1))
    m = _EVENT_START_RE.search(text)
    if m:
        info["event_start_mass_kg"] = float(m.group(1))
    return info


def find_output_files(raw_dir):
    """Locate the four SESAM files anywhere under raw_dir (pyDRAMA writes run_0/reentry/)."""
    def first(pattern):
        hits = sorted(glob.glob(os.path.join(raw_dir, "**", pattern), recursive=True))
        return hits[0] if hits else None
    return {
        "aero": first("*_AeroThermalHistory.txt"),
        "traj": first("*_Trajectory.txt"),
        "fragments": first("*ImpactingFragments.xml"),
        "log": first("sesam.log"),
    }


# =============================================================================
# 4. Merge + CSV (spec 4.6 / 4.7)
# =============================================================================

CSV_COLUMNS = [
    "time_s", "altitude_km", "velocity_kms", "temperature_K", "mass_kg", "thick_mm",
    "lat_deg", "lon_deg", "downrange_km", "flight_path_deg", "heading_deg",
    "drag", "lift", "side", "knudsen", "mach", "density_kgm3", "dynamic_pressure_Pa",
    "load_factor_g", "convective_heat_W", "radiative_heat_W", "oxidation_heat_W",
    "rad_cooling_W", "integrated_heat_J", "visibility_factor",
]
# CSV column -> trajectory-file column ("time_s" / "altitude_km" handled explicitly)
_TRAJ_MAP = {"velocity_kms": "velocity", "lat_deg": "lat", "lon_deg": "lon",
             "downrange_km": "downrange", "flight_path_deg": "path", "heading_deg": "heading",
             "drag": "drag", "lift": "lift", "side": "side", "knudsen": "knudsen", "mach": "mach",
             "density_kgm3": "density", "dynamic_pressure_Pa": "dynamicPressure",
             "load_factor_g": "loadFactor"}
# CSV column -> aerothermal-file column
_AERO_MAP = {"temperature_K": "temp", "mass_kg": "mass", "thick_mm": "thick",
             "convective_heat_W": "convectiveHeat", "radiative_heat_W": "radiativeHeat",
             "oxidation_heat_W": "oxidationHeat", "rad_cooling_W": "radCooling",
             "integrated_heat_J": "integratedHeat", "visibility_factor": "visibilityFactor"}


def _merged_row(t, aero, traj):
    row = {c: None for c in CSV_COLUMNS}
    row["time_s"] = t
    if traj is not None:
        row["altitude_km"] = traj["altitude"]
        for col, src in _TRAJ_MAP.items():
            row[col] = traj[src]
    if aero is not None:
        if row["altitude_km"] is None:
            row["altitude_km"] = aero["altitude"]
        for col, src in _AERO_MAP.items():
            row[col] = aero[src]
    return row


def merge_histories(aero_rows, traj_rows):
    """Join the two SESAM histories on time. Identical grids zip; otherwise outer-join + warning."""
    warnings = []
    a_times = [r["time"] for r in aero_rows]
    t_times = [r["time"] for r in traj_rows]
    if a_times == t_times:
        rows = [_merged_row(t, a, tr) for t, a, tr in zip(a_times, aero_rows, traj_rows)]
    else:
        warnings.append("time grids differ (aero {} rows, traj {} rows)".format(len(aero_rows), len(traj_rows)))
        by_aero = {r["time"]: r for r in aero_rows}
        by_traj = {r["time"]: r for r in traj_rows}
        rows = [_merged_row(t, by_aero.get(t), by_traj.get(t)) for t in sorted(set(by_aero) | set(by_traj))]
    return rows, warnings


def write_csv(path, rows):
    """One row per time step; floats via str() (shortest round-trip repr), None as empty."""
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=CSV_COLUMNS)
        writer.writeheader()
        for row in rows:
            writer.writerow({c: ("" if row.get(c) is None else row[c]) for c in CSV_COLUMNS})


# =============================================================================
# 5. Statistics (spec 4.8)
# =============================================================================

def radius_mm_from_mass(mass_kg, material: Material = DEFAULT_MATERIAL):
    """Radius of a solid sphere of the given mass and material; 0 for zero mass."""
    if mass_kg is None or mass_kg <= 0.0:
        return 0.0
    return 1000.0 * (3.0 * mass_kg / (4.0 * math.pi * material.density_kgm3)) ** (1.0 / 3.0)


def classify_outcome(reason, final_mass_kg, initial_mass_kg):
    if reason == "ground impact":
        return "survived"
    if reason == "ballooning":
        return "demised"
    if (final_mass_kg is not None and initial_mass_kg > 0.0
            and final_mass_kg < DEMISE_MASS_FRACTION * initial_mass_kg):
        return "demised"
    return "other"


def _at_melt(temp, melting_temperature_K=T_MELT_AA7075):
    return temp is not None and abs(temp - melting_temperature_K) <= MELT_TOLERANCE_K


def _last_with(rows, key):
    for row in reversed(rows):
        if row.get(key) is not None:
            return row
    return None


def compute_stats(run, rows, fragments, log_info):
    """The statistics block of the JSON plus a list of warning strings."""
    results = {}
    warnings = []

    # -- maximum temperature (first occurrence) ------------------------------------
    with_temp = [r for r in rows if r["temperature_K"] is not None]
    if with_temp:
        hottest = max(with_temp, key=lambda r: r["temperature_K"])   # max() keeps the first maximum
        results["max_temperature_K"] = hottest["temperature_K"]
        results["time_of_max_temperature_s"] = hottest["time_s"]
        results["altitude_of_max_temperature_km"] = hottest["altitude_km"]
    else:
        results["max_temperature_K"] = None
        results["time_of_max_temperature_s"] = None
        results["altitude_of_max_temperature_km"] = None

    # -- time at melting temperature -------------------------------------------------
    t_melt = run.material.melting_temperature_K
    melt_rows = [r for r in rows if _at_melt(r["temperature_K"], t_melt)]
    duration = 0.0
    for prev, cur in zip(rows, rows[1:]):
        if _at_melt(prev["temperature_K"], t_melt) and _at_melt(cur["temperature_K"], t_melt):
            duration += cur["time_s"] - prev["time_s"]
    results["melting_temperature_K"] = t_melt
    results["melt_tolerance_K"] = MELT_TOLERANCE_K
    results["time_at_melting_temperature_s"] = duration
    results["n_rows_at_melt"] = len(melt_rows)
    results["first_time_at_melt_s"] = melt_rows[0]["time_s"] if melt_rows else None
    results["last_time_at_melt_s"] = melt_rows[-1]["time_s"] if melt_rows else None
    results["altitude_first_melt_km"] = melt_rows[0]["altitude_km"] if melt_rows else None
    results["altitude_last_melt_km"] = melt_rows[-1]["altitude_km"] if melt_rows else None

    # -- final mass: XML -> log -> history file ---------------------------------------
    initial_mass = run.mass_kg
    last_mass_row = _last_with(rows, "mass_kg")
    if fragments is not None and fragments.get("mass_kg") is not None:
        final_mass, source = fragments["mass_kg"], "impacting_fragments_xml"
    elif log_info.get("event_end_mass_kg") is not None:
        final_mass, source = log_info["event_end_mass_kg"], "sesam_log_event_end"
    elif last_mass_row is not None:
        final_mass, source = last_mass_row["mass_kg"], "history_file"
    else:
        final_mass, source = None, None
    results["initial_mass_kg"] = initial_mass
    results["final_mass_kg"] = final_mass
    results["final_mass_source"] = source
    results["mass_loss_fraction"] = (1.0 - final_mass / initial_mass) if final_mass is not None else None
    results["final_radius_mm"] = radius_mm_from_mass(final_mass, run.material) if final_mass is not None else None
    last_thick_row = _last_with(rows, "thick_mm")
    results["final_thick_mm"] = last_thick_row["thick_mm"] if last_thick_row else None

    # -- final velocity and trajectory end -------------------------------------------
    last_traj = _last_with(rows, "velocity_kms")
    if fragments is not None and fragments.get("velocity_kms") is not None:
        results["final_velocity_kms"] = fragments["velocity_kms"]
        results["final_velocity_source"] = "impacting_fragments_xml"
    elif last_traj is not None:
        results["final_velocity_kms"] = last_traj["velocity_kms"]
        results["final_velocity_source"] = "trajectory_file"
    else:
        results["final_velocity_kms"] = None
        results["final_velocity_source"] = None
    last = rows[-1] if rows else None
    results["final_time_s"] = last["time_s"] if last else None
    results["final_altitude_km"] = last["altitude_km"] if last else None
    results["final_latitude_deg"] = last_traj["lat_deg"] if last_traj else None
    results["final_longitude_deg"] = last_traj["lon_deg"] if last_traj else None
    results["downrange_km"] = last_traj["downrange_km"] if last_traj else None

    # -- end of life --------------------------------------------------------------------
    reason = log_info.get("end_of_life_reason", "unknown")
    results["end_of_life_reason"] = reason
    results["outcome"] = classify_outcome(reason, final_mass, initial_mass)
    results["n_rows"] = len(rows)

    # -- warnings -----------------------------------------------------------------------
    if initial_mass < COARSE_MASS_LIMIT_KG:
        warnings.append("CSV mass column ({:g} kg resolution) is coarse for initial mass {:.6g} kg".format(
            CSV_MASS_RESOLUTION_KG, initial_mass))
    if len(rows) == 1:
        warnings.append("history has a single row")
    if reason == "ground impact" and fragments is None:
        warnings.append("ground impact reported but no ImpactingFragments.xml entry")
    if reason == "unknown":
        warnings.append("end-of-life reason not found in sesam.log")
    return results, warnings


# =============================================================================
# 6. Running SESAM through pyDRAMA (spec 4.5 / 4.9 / 4.10)
# =============================================================================

class DramaNotAvailable(RuntimeError):
    """pyDRAMA (package 'drama') cannot be imported."""


def _import_sara():
    os.environ.setdefault("DRAMA_INSTALL_PATH", DRAMA_INSTALL_PATH)
    try:
        import drama
        from drama import sara
    except ImportError as exc:
        raise DramaNotAvailable(
            "pyDRAMA not importable ({}). Install it into this interpreter with:\n"
            "    pip install '{}/TOOLS/drama_python_package'".format(exc, DRAMA_INSTALL_PATH))
    return sara, getattr(drama, "__version__", None)


def _log_tail(text, n=40):
    if text is None:
        return ""
    if isinstance(text, bytes):
        text = text.decode(errors="replace")
    return "\n".join(str(text).splitlines()[-n:])


def _base_document(run, cfg, name, csv_path, run_raw, pydrama_version):
    settings = {k: v for k, v in cfg.items() if k != "objects"}
    return {
        "schema_version": 1,
        "run_name": name,
        "status": None,
        "error": None,
        "inputs": {
            "diameter_mm": run.diameter_mm,
            "radius_m": run.radius_m,
            "initial_mass_kg": run.mass_kg,
            "cross_section_m2": run.cross_section_m2,
            "initial_temperature_K": run.temperature_K,
            "initial_velocity_kms": run.velocity_kms,
            "initial_altitude_km": run.altitude_km,
            "flight_path_angle_deg": run.flight_path_deg,
            "heading_deg": run.heading_deg % 360.0,
            "latitude_deg": run.lat_deg,
            "longitude_deg": run.lon_deg,
            "epoch_utc": run.epoch.strftime("%Y-%m-%dT%H:%M:%S"),
            "material": MATERIAL_NAME,
            "material_density_kgm3": RHO_AA7075,
            "melting_temperature_K": T_MELT_AA7075,
            "melt_tolerance_K": MELT_TOLERANCE_K,
            "energy_threshold_J": ENERGY_THRESHOLD_J,
            "sesam_settings": json.loads(json.dumps(settings, default=str)),
        },
        "results": None,
        "warnings": [],
        "files": {"csv": os.path.abspath(csv_path), "raw_dir": run_raw},
        "provenance": {
            "drama_install_path": os.environ.get("DRAMA_INSTALL_PATH", DRAMA_INSTALL_PATH),
            "sesam_version": None,
            "pydrama_version": pydrama_version,
            "script_version": SCRIPT_VERSION,
            "wall_time_s": None,
            "created_utc": None,
            "hostname": socket.gethostname(),
        },
    }


def _finish(doc, json_path, status, t0, error=None):
    doc["status"] = status
    doc["error"] = error
    doc["provenance"]["wall_time_s"] = round(time.time() - t0, 3)
    doc["provenance"]["created_utc"] = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    with open(json_path, "w") as fh:
        json.dump(doc, fh, indent=2, default=str)
    return doc


def run_sphere(run, outdir, raw_dir, timeout, keep_raw, fap_day_lines, fap_mon_lines):
    """Run SESAM for `run`; write <outdir>/<run_name>.csv and .json; return the JSON document."""
    sara, pydrama_version = _import_sara()
    name = run_name(run)
    os.makedirs(outdir, exist_ok=True)
    os.makedirs(raw_dir, exist_ok=True)
    csv_path = os.path.join(outdir, name + ".csv")
    json_path = os.path.join(outdir, name + ".json")
    run_raw = os.path.join(raw_dir, name)
    if os.path.isdir(run_raw):
        shutil.rmtree(run_raw)          # pyDRAMA's copytree refuses an existing destination
    if os.path.isfile(csv_path):
        os.remove(csv_path)
    cfg = build_config(run)
    doc = _base_document(run, cfg, name, csv_path, run_raw, pydrama_version)
    t0 = time.time()
    try:
        results = sara.run(config=[cfg], save_output_dirs=run_raw, keep_output_files="all",
                           fap_day_content=fap_day_lines, fap_mon_content=fap_mon_lines,
                           parallel=False, timeout=timeout, log_level="ERROR", spell_check=False)
    except Exception as exc:  # noqa: BLE001 -- anything pyDRAMA raises is a failed run
        return _finish(doc, json_path, "error", t0, "pyDRAMA raised {!r}".format(exc))

    errors = results.get("errors") or []
    if errors or not results.get("results"):
        status, messages = "error", []
        for err in errors:
            text = str(err.get("status"))
            if "timeout" in text.lower():
                status = "timeout"
            messages.append(text + "\n" + _log_tail(err.get("reentry_logfile") or err.get("logfile")))
        if not errors:
            messages.append("pyDRAMA returned no result")
        return _finish(doc, json_path, status, t0, "\n".join(messages))

    files = find_output_files(run_raw)
    if not files["aero"] or not files["traj"]:
        return _finish(doc, json_path, "error", t0, "history files not found under {}\n{}".format(
            run_raw, _log_tail(results["results"][0].get("reentry_logfile"))))

    aero = read_sara_table(files["aero"], AERO_COLUMNS)
    traj = read_sara_table(files["traj"], TRAJ_COLUMNS)
    rows, warnings = merge_histories(aero, traj)
    if not rows:
        return _finish(doc, json_path, "error", t0, "history files contain no data rows")
    fragments = parse_impacting_fragments(files["fragments"])
    log_info = parse_sesam_log(files["log"])
    stats, more_warnings = compute_stats(run, rows, fragments, log_info)
    write_csv(csv_path, rows)
    doc["results"] = stats
    doc["warnings"] = warnings + more_warnings
    doc["provenance"]["sesam_version"] = parse_sesam_version(files["aero"])
    if keep_raw:
        doc["files"]["raw_dir"] = run_raw
    else:
        shutil.rmtree(run_raw, ignore_errors=True)
        doc["files"]["raw_dir"] = None
    return _finish(doc, json_path, "ok", t0)

# =============================================================================
# 7. Command line (spec 4.1)
# =============================================================================

def parse_epoch(text):
    """ISO-8601 UTC timestamp (a trailing 'Z' is accepted) -> naive datetime."""
    text = text.strip()
    if text.endswith("Z"):
        text = text[:-1]
    try:
        return datetime.fromisoformat(text)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("invalid epoch {!r}: {}".format(text, exc))


def build_parser():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--velocity", type=float, required=True, help="initial velocity [km/s]")
    p.add_argument("--altitude", type=float, required=True, help="initial geodetic altitude [km]")
    p.add_argument("--temperature", type=float, required=True, help="initial bulk temperature [K]")
    p.add_argument("--diameter", type=float, required=True, help="sphere diameter [mm]")
    p.add_argument("--flight-path-angle", type=float, default=0.0, help="[deg], negative = descending (default 0)")
    p.add_argument("--heading", type=float, default=0.0, help="[deg] (default 0)")
    p.add_argument("--lat", type=float, default=0.0, help="geodetic latitude [deg] (default 0)")
    p.add_argument("--lon", type=float, default=0.0, help="longitude [deg] (default 0)")
    p.add_argument("--epoch", type=parse_epoch, default=PARENT_EPOCH,
                   help="ISO-8601 UTC (default %(default)s, the parent run's epoch)")
    p.add_argument("--outdir", default=DEFAULT_OUTDIR, help="CSV/JSON directory (default %(default)s)")
    p.add_argument("--raw-dir", default=DEFAULT_RAW_DIR, help="raw DRAMA trees root (default %(default)s)")
    p.add_argument("--timeout", type=int, default=DEFAULT_TIMEOUT_S, help="SESAM timeout [s] (default %(default)s)")
    p.add_argument("--keep-raw", action="store_true", help="keep the raw DRAMA tree on success")
    p.add_argument("--dry-run", action="store_true", help="print the config and run name; run nothing")
    p.add_argument("--quiet", action="store_true", help="no console output except errors")
    p.add_argument("--fap-day", default=DEFAULT_FAP_DAY, help="DRAMA fap_day.dat (default %(default)s)")
    p.add_argument("--fap-mon", default=DEFAULT_FAP_MON, help="DRAMA fap_mon.dat (default %(default)s)")
    return p


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.velocity <= 0.0:
        parser.error("--velocity must be > 0")
    if args.altitude < 0.0:
        parser.error("--altitude must be >= 0")
    if args.temperature <= 0.0:
        parser.error("--temperature must be > 0")
    if args.diameter <= 0.0:
        parser.error("--diameter must be > 0")
    run = SphereRun(velocity_kms=args.velocity, altitude_km=args.altitude,
                    temperature_K=args.temperature, diameter_mm=args.diameter,
                    flight_path_deg=args.flight_path_angle, heading_deg=args.heading,
                    lat_deg=args.lat, lon_deg=args.lon, epoch=args.epoch)
    if args.dry_run:
        print(json.dumps(build_config(run), indent=2, default=str))
        print("run name: " + run_name(run))
        return 0
    for label, path in (("--fap-day", args.fap_day), ("--fap-mon", args.fap_mon)):
        if not os.path.isfile(path):
            parser.error("{} file not found: {}".format(label, path))
    try:
        doc = run_sphere(run, args.outdir, args.raw_dir, args.timeout, args.keep_raw,
                         read_lines(args.fap_day), read_lines(args.fap_mon))
    except DramaNotAvailable as exc:
        print("ERROR: {}".format(exc), file=sys.stderr)
        return 2
    if doc["status"] != "ok":
        print("ERROR ({}): {}".format(doc["status"], doc["error"]), file=sys.stderr)
        return 1
    if not args.quiet:
        r = doc["results"]
        print("{}: {} after {:.1f} s ({}); Tmax {:.1f} K, {:.1f} s at melt, final mass {:.6g} kg "
              "(r = {:.3f} mm), final velocity {:.4f} km/s".format(
                  doc["run_name"], r["outcome"], r["final_time_s"], r["end_of_life_reason"],
                  r["max_temperature_K"], r["time_at_melting_temperature_s"], r["final_mass_kg"],
                  r["final_radius_mm"], r["final_velocity_kms"]))
        print("  csv  -> {}".format(doc["files"]["csv"]))
        print("  json -> {}".format(os.path.abspath(os.path.join(args.outdir, doc["run_name"] + ".json"))))
    return 0


if __name__ == "__main__":
    sys.exit(main())

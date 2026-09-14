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

# --- material: drama-AA7075 exactly as shipped in DRAMA 4.1.4 -----------------
MATERIAL_NAME = "drama-AA7075"
RHO_AA7075 = 2813.0            # [kg/m3]
T_MELT_AA7075 = 850.0          # [K]
MELT_TOLERANCE_K = 0.5         # |T - T_melt| <= tol counts as "at melting temperature"

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

def sphere_mass_kg(diameter_mm: float) -> float:
    """Mass of a solid AA7075 sphere: rho * 4/3 * pi * r^3."""
    r = diameter_mm / 2000.0
    return RHO_AA7075 * 4.0 / 3.0 * math.pi * r ** 3


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

    @property
    def radius_m(self) -> float:
        return self.diameter_mm / 2000.0

    @property
    def mass_kg(self) -> float:
        return sphere_mass_kg(self.diameter_mm)

    @property
    def cross_section_m2(self) -> float:
        return sphere_cross_section_m2(self.diameter_mm)


def run_name(run: SphereRun) -> str:
    """Fixed-width name so that lexicographic order equals parameter order (spec 4.3)."""
    return "sphere_d{:06.2f}mm_T{:06.1f}K_v{:08.5f}kms_h{:07.3f}km".format(
        run.diameter_mm, run.temperature_K, run.velocity_kms, run.altitude_km)


def read_lines(path: str) -> list:
    """Non-blank lines of a text file without trailing newlines (DRAMA fap_day/fap_mon content)."""
    with open(path) as fh:
        return [ln.rstrip("\n") for ln in fh if ln.strip()]


# =============================================================================
# 2. SESAM configuration (spec 4.4) -- mirrors the parent satellite run except energyThreshold
# =============================================================================

def build_config(run: SphereRun) -> dict:
    """The complete pyDRAMA SARA configuration for one sphere (pass it as config=[cfg])."""
    obj = {
        "name": OBJECT_NAME,
        "uniqueID": OBJECT_UUID,
        "primitive": {"sphere": {"radius": run.radius_m}},
        "mass": run.mass_kg,
        "material": MATERIAL_NAME,
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
        # ---- the model (built-in materials.xml is used: no materialList) ----
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

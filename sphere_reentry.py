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

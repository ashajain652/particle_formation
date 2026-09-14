#!/usr/bin/env python3
"""
sphere_sweep.py -- run sphere_reentry.py over a grid of diameters, temperatures and velocities

The sphere's remaining initial-state values (altitude, latitude, longitude,
flight-path angle, heading, epoch) are inherited from a parent-satellite SESAM
run (a DRAMA GUI dmf_output.json) at the point of its final descent where the
parent reached the sphere's velocity. The velocity grid runs from 7.5 km/s down
to the parent's minimum velocity (its impact velocity), which is read from the
parent file and logged.

Usage
-----
    python sphere_sweep.py --dry-run                       # show the matrix, run nothing
    python sphere_sweep.py                                 # everything pending in one confirmed batch
    python sphere_sweep.py --batch-size 5000               # confirm + choose cores before every batch
    python sphere_sweep.py --diameters 50 --temperatures 300 --velocities 7.5,0.5 --yes --cores 2

Outputs (under --outdir, default sphere_sweep_output/):
    runs/<run_name>.csv|.json   per run (written by sphere_reentry.py)
    raw/<run_name>/             raw DRAMA tree, failed runs only
    sweep_manifest.json         every matrix point with its state and status
    sweep_summary.csv           one row per point with the statistics
    sweep.log                   log of prompts, answers, failures

Requires conda env drama_env (pyDRAMA, numpy, tqdm); see README.md.
Design: docs/superpowers/specs/2026-09-13-sphere-reentry-sweep-design.md
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import logging
import math
import os
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass, asdict
from datetime import datetime, timedelta, timezone

import numpy as np
from tqdm import tqdm

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
if SCRIPT_DIR not in sys.path:
    sys.path.insert(0, SCRIPT_DIR)
import sphere_reentry as sr  # noqa: E402

SPHERE_REENTRY = os.path.join(SCRIPT_DIR, "sphere_reentry.py")
DEFAULT_PARENT = os.path.join(SCRIPT_DIR, "Generic_Satellite Reentry", "output", "dmf_output.json")
DEFAULT_PARENT_OBJECT = "Main Body"
DEFAULT_OUTDIR = os.path.join(SCRIPT_DIR, "sphere_sweep_output")
DEFAULT_TIMEOUT_S = 600
V_TOP_KMS = 7.5                    # fixed upper bound of the velocity grid
N_VELOCITIES = 100                 # linspace(V_TOP_KMS, v_min_parent, N_VELOCITIES)
DIAMETERS_MM = [float(d) for d in range(5, 101, 5)]          # 5 .. 100 mm step 5
TEMPERATURES_K = [float(t) for t in range(300, 751, 10)]     # 300 .. 750 K step 10
SECONDS_PER_RUN_ESTIMATE = 0.3     # per run per core, before any run has been timed
STDERR_TAIL_LINES = 20
EPOCH_ARG_FORMAT = "%Y-%m-%dT%H:%M:%S"
PARENT_EPOCH_FORMAT = "%Y-%m-%dT%H:%M:%S.%fZ"

log = logging.getLogger("sweep")


# =============================================================================
# 1. Parent trajectory (spec 5.2)
# =============================================================================

def sha256_of_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def descending_branch_start(velocities):
    """Index one past the last velocity increase (0 when velocity never increases)."""
    start = 0
    for i in range(len(velocities) - 1):
        if velocities[i + 1] > velocities[i]:
            start = i + 1
    return start


@dataclass
class ParentTrajectory:
    path: str
    sha256: str
    object_name: str
    epoch: datetime
    rows: list            # full trajectory rows (dicts straight from dmf_output.json)
    branch_start: int     # first row of the final descending branch

    @property
    def branch(self):
        return self.rows[self.branch_start:]

    @property
    def v_max(self):
        return self.branch[0]["velocity"]

    @property
    def v_min(self):
        return self.branch[-1]["velocity"]

    @property
    def v_min_row(self):
        return len(self.rows) - 1

    def describe_v_min(self):
        last = self.branch[-1]
        return "v_min = {:.6f} km/s at t = {:.3f} s, altitude {:.3f} km (parent row {})".format(
            last["velocity"], last["time"], last["altitude"], self.v_min_row)


def load_parent(path, object_name):
    """Read a DRAMA GUI dmf_output.json and select one object's trajectory."""
    with open(path) as fh:
        doc = json.load(fh)
    satellites = doc["singleModuleOutputs"]["satellites"]
    if len(satellites) != 1:
        raise ValueError("expected exactly one satellite in {}, found {}: {}".format(
            path, len(satellites), [s.get("name") for s in satellites]))
    phases = satellites[0]["missionPhases"]
    if len(phases) != 1:
        raise ValueError("expected exactly one mission phase, found {}: {}".format(
            len(phases), [p.get("name") for p in phases]))
    epochs = phases[0]["epochs"]
    if len(epochs) != 1:
        raise ValueError("expected exactly one epoch, found {}: {}".format(
            len(epochs), [e.get("epoch") for e in epochs]))
    epoch = datetime.strptime(epochs[0]["epoch"], PARENT_EPOCH_FORMAT)
    results = epochs[0]["analysisModules"][0]["results"]
    wanted = object_name.replace(" ", "_")
    trajectory_keys = [k for k in results if k.endswith("_Trajectory.txt")]
    keys = [k for k in trajectory_keys if len(k.split(".")) >= 2 and k.split(".")[1] == wanted]
    if len(keys) != 1:
        raise ValueError("expected one trajectory for object {!r}, found {}; available: {}".format(
            object_name, keys, trajectory_keys))
    rows = results[keys[0]]
    if len(rows) < 2:
        raise ValueError("parent trajectory {} has fewer than 2 rows".format(keys[0]))
    start = descending_branch_start([r["velocity"] for r in rows])
    if len(rows) - start < 2:
        raise ValueError("descending branch of the parent trajectory has fewer than 2 rows")
    return ParentTrajectory(path=path, sha256=sha256_of_file(path), object_name=object_name,
                            epoch=epoch, rows=rows, branch_start=start)


def lerp(a, b, f):
    return a + f * (b - a)


def angle_lerp(a_deg, b_deg, f):
    """Interpolate two angles through their sines/cosines; result in (-180, 180]."""
    ra, rb = math.radians(a_deg), math.radians(b_deg)
    s = lerp(math.sin(ra), math.sin(rb), f)
    c = lerp(math.cos(ra), math.cos(rb), f)
    return math.degrees(math.atan2(s, c))


def wrap_longitude(lon):
    return ((lon + 180.0) % 360.0) - 180.0


def interpolate_state(parent, v_kms):
    """Parent state at velocity v on the descending branch (linear in velocity), or None outside it."""
    if v_kms > parent.v_max or v_kms < parent.v_min:
        return None
    branch = parent.branch
    for a, b in zip(branch, branch[1:]):
        if a["velocity"] >= v_kms >= b["velocity"]:
            dv = a["velocity"] - b["velocity"]
            f = (a["velocity"] - v_kms) / dv if dv else 0.0
            t = lerp(a["time"], b["time"], f)
            return {
                "altitude_km": lerp(a["altitude"], b["altitude"], f),
                "lat_deg": lerp(a["lat"], b["lat"], f),
                "lon_deg": wrap_longitude(angle_lerp(a["lon"], b["lon"], f)),
                "flight_path_deg": lerp(a["path"], b["path"], f),
                "heading_deg": angle_lerp(a["heading"], b["heading"], f) % 360.0,
                "time_s": t,
                "epoch": parent.epoch + timedelta(seconds=t),
            }
    return None


# =============================================================================
# 2. Grid and matrix (spec 5.3)
# =============================================================================

def velocity_grid(v_min_parent):
    """100 velocities from 7.5 km/s down to the parent's minimum velocity, both included."""
    return [float(v) for v in np.linspace(V_TOP_KMS, v_min_parent, N_VELOCITIES)]


def parse_float_list(text):
    values = [float(x) for x in text.split(",") if x.strip()]
    if not values:
        raise argparse.ArgumentTypeError("expected a comma-separated list of numbers")
    return values


def format_arg(x):
    """How every numeric argument is passed to sphere_reentry.py."""
    return "{:.6f}".format(x)


@dataclass
class Point:
    diameter_mm: float
    temperature_K: float
    velocity_kms: float
    altitude_km: float = None
    flight_path_deg: float = None
    heading_deg: float = None
    lat_deg: float = None
    lon_deg: float = None
    epoch_utc: str = None
    run_name: str = None
    status: str = "pending"        # pending | skipped | done | failed
    skip_reason: str = None
    returncode: int = None
    wall_time_s: float = None
    stderr_tail: str = None


def sphere_run_for(point):
    """The SphereRun script 1 will construct from the formatted CLI strings (so run names agree)."""
    return sr.SphereRun(
        velocity_kms=float(format_arg(point.velocity_kms)),
        altitude_km=float(format_arg(point.altitude_km)),
        temperature_K=float(format_arg(point.temperature_K)),
        diameter_mm=float(format_arg(point.diameter_mm)),
        flight_path_deg=float(format_arg(point.flight_path_deg)),
        heading_deg=float(format_arg(point.heading_deg)),
        lat_deg=float(format_arg(point.lat_deg)),
        lon_deg=float(format_arg(point.lon_deg)),
        epoch=datetime.strptime(point.epoch_utc, EPOCH_ARG_FORMAT),
    )


def build_points(parent, diameters, temperatures, velocities, limit=None):
    """All matrix points ordered diameter -> temperature -> velocity (descending), with parent states."""
    points = []
    for d in sorted(diameters):
        for T in sorted(temperatures):
            for v in sorted(velocities, reverse=True):
                point = Point(diameter_mm=d, temperature_K=T, velocity_kms=v)
                state = interpolate_state(parent, v)
                if state is None:
                    point.status = "skipped"
                    point.skip_reason = "parent never reaches velocity"
                else:
                    point.altitude_km = state["altitude_km"]
                    point.flight_path_deg = state["flight_path_deg"]
                    point.heading_deg = state["heading_deg"]
                    point.lat_deg = state["lat_deg"]
                    point.lon_deg = state["lon_deg"]
                    point.epoch_utc = state["epoch"].strftime(EPOCH_ARG_FORMAT)
                    point.run_name = sr.run_name(sphere_run_for(point))
                points.append(point)
    if limit is not None:
        points = points[:limit]
    return points

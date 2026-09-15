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
    python sphere_sweep.py --material drama-TiAl6v4          # every sphere in another DRAMA metal
    python sphere_sweep.py --material-file examples/material_al_li_2195.json   # or a custom metal

The material is one setting for the whole sweep (it is not a grid axis). A non-default
material adds a "_m<name>" suffix to every run name, so sweeps of different materials can
share an output directory without overwriting each other; it is recorded in the manifest.

Outputs (under --outdir, default sphere_sweep_output/):
    runs/<run_name>.csv|.json   per run (written by sphere_reentry.py)
    raw/<run_name>/             raw DRAMA tree, failed runs only
    sweep_manifest_<material>.json   every matrix point with its state and status
    sweep_summary_<material>.csv     one row per point with the statistics
    sweep.log                        log of prompts, answers, failures (all invocations)

<material> is the sweep's material with DRAMA's "drama-" prefix dropped (sweep_summary_AA7075.csv
by default), so sweeps of different materials never overwrite each other's manifest or summary.

Requires conda env drama_env (pyDRAMA, numpy, tqdm); see README.md.
Design: docs/superpowers/specs/2026-09-13-sphere-reentry-sweep-design.md
"""
from __future__ import annotations

import argparse
import csv
import functools
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


def sphere_run_for(point, material=sr.DEFAULT_MATERIAL):
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
        material=material,
    )


def build_points(parent, diameters, temperatures, velocities, limit=None, material=sr.DEFAULT_MATERIAL):
    """All matrix points ordered diameter -> temperature -> velocity (descending), with parent states.

    `material` is the sweep's single material; it only affects the run names (a non-default
    material gets a suffix so it can never share output files with another material).
    """
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
                    point.run_name = sr.run_name(sphere_run_for(point, material))
                points.append(point)
    # Two points that round to the same run_name would run concurrently under the thread-pool
    # executor and clobber each other's output files (run A's JSON paired with run B's CSV),
    # not just harmlessly overwrite each other sequentially. Unreachable on the default grid;
    # reachable via close-together explicit --velocities overrides after format_arg's rounding.
    seen = {}
    for p in points:
        if p.run_name is not None:
            if p.run_name in seen:
                other = seen[p.run_name]
                raise ValueError(
                    "duplicate run_name {!r}: points ({}, {}, {}) and ({}, {}, {}) would run "
                    "concurrently and clobber each other's output files".format(
                        p.run_name, other.diameter_mm, other.temperature_K, other.velocity_kms,
                        p.diameter_mm, p.temperature_K, p.velocity_kms))
            seen[p.run_name] = p
    if limit is not None:
        points = points[:limit]
    return points


# =============================================================================
# 3. Manifest and resume (spec 5.4)
# =============================================================================

def utc_now_str():
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def manifest_document(parent, grid, settings, points, created_utc=None):
    now = utc_now_str()
    last = parent.branch[-1]
    return {
        "created_utc": created_utc or now,
        "updated_utc": now,
        "parent": {
            "path": os.path.abspath(parent.path),
            "sha256": parent.sha256,
            "object": parent.object_name,
            "epoch_utc": parent.epoch.strftime(EPOCH_ARG_FORMAT),
            "branch_first_row": parent.branch_start,
            "branch_last_row": len(parent.rows) - 1,
            "v_min_kms": parent.v_min,
            "v_min_time_s": last["time"],
            "v_min_altitude_km": last["altitude"],
            "v_max_kms": parent.v_max,
        },
        "grid": grid,
        "settings": settings,
        "points": [asdict(p) for p in points],
    }


def write_manifest(path, doc):
    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(doc, fh, indent=2)
    os.replace(tmp, path)


def load_manifest(path):
    with open(path) as fh:
        return json.load(fh)


def _json_status(path):
    """status field of a run JSON, or None when the file is missing/unreadable."""
    if not os.path.isfile(path):
        return None
    try:
        with open(path) as fh:
            return json.load(fh).get("status") or "unreadable"
    except (OSError, ValueError):
        return "unreadable"


def apply_resume(points, runs_dir, force=False, retry_failed=False):
    """Mark points done/failed/pending from the run JSONs already in runs_dir."""
    for point in points:
        if point.status == "skipped":
            continue
        status = _json_status(os.path.join(runs_dir, point.run_name + ".json"))
        if status is None:
            point.status = "pending"
        elif status == "ok":
            point.status = "pending" if force else "done"
        else:
            point.status = "pending" if (force or retry_failed) else "failed"


# =============================================================================
# 4. Aggregate summary (spec 5.6)
# =============================================================================

SUMMARY_COLUMNS = [
    "diameter_mm", "initial_temperature_K", "initial_velocity_kms", "initial_altitude_km",
    "flight_path_angle_deg", "heading_deg", "latitude_deg", "longitude_deg", "epoch_utc",
    "status", "skip_reason", "max_temperature_K", "final_mass_kg", "final_mass_source",
    "mass_loss_fraction", "time_at_melting_temperature_s", "final_velocity_kms",
    "final_radius_mm", "final_altitude_km", "end_of_life_reason", "outcome", "wall_time_s", "run_name",
]
_RESULT_COLUMNS = ["max_temperature_K", "final_mass_kg", "final_mass_source", "mass_loss_fraction",
                   "time_at_melting_temperature_s", "final_velocity_kms", "final_radius_mm",
                   "final_altitude_km", "end_of_life_reason", "outcome"]


def _blank(value):
    return "" if value is None else value


def summary_row(point, runs_dir):
    row = {c: "" for c in SUMMARY_COLUMNS}
    row.update({
        "diameter_mm": point.diameter_mm,
        "initial_temperature_K": point.temperature_K,
        "initial_velocity_kms": point.velocity_kms,
        "initial_altitude_km": _blank(point.altitude_km),
        "flight_path_angle_deg": _blank(point.flight_path_deg),
        "heading_deg": _blank(point.heading_deg),
        "latitude_deg": _blank(point.lat_deg),
        "longitude_deg": _blank(point.lon_deg),
        "epoch_utc": _blank(point.epoch_utc),
        "status": point.status,
        "skip_reason": _blank(point.skip_reason),
        "wall_time_s": _blank(point.wall_time_s),
        "run_name": _blank(point.run_name),
    })
    if point.status == "done":
        try:
            with open(os.path.join(runs_dir, point.run_name + ".json")) as fh:
                results = json.load(fh).get("results") or {}
        except (OSError, ValueError):
            results = {}
        for column in _RESULT_COLUMNS:
            row[column] = _blank(results.get(column))
    return row


def write_summary(path, points, runs_dir):
    tmp = path + ".tmp"
    with open(tmp, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=SUMMARY_COLUMNS)
        writer.writeheader()
        for point in points:
            writer.writerow(summary_row(point, runs_dir))
    os.replace(tmp, path)


# =============================================================================
# 5. Execution: one subprocess per point, batches, prompts, progress (spec 5.5)
# =============================================================================

def material_cli_args(material_name, material_file):
    """The sphere_reentry.py arguments that select the sweep's material ([] for the default)."""
    if material_file:
        return ["--material-file", os.path.abspath(material_file)]
    if material_name != sr.MATERIAL_NAME:
        return ["--material", material_name]
    return []


def build_command(point, runs_dir, raw_dir, timeout, python=None, material_args=()):
    return [python or sys.executable, SPHERE_REENTRY,
            "--velocity", format_arg(point.velocity_kms),
            "--altitude", format_arg(point.altitude_km),
            "--temperature", format_arg(point.temperature_K),
            "--diameter", format_arg(point.diameter_mm),
            "--flight-path-angle", format_arg(point.flight_path_deg),
            "--heading", format_arg(point.heading_deg),
            "--lat", format_arg(point.lat_deg),
            "--lon", format_arg(point.lon_deg),
            "--epoch", point.epoch_utc,
            "--outdir", runs_dir,
            "--raw-dir", raw_dir,
            "--timeout", str(int(timeout)),
            *material_args,
            "--quiet"]


def run_point(point, runs_dir, raw_dir, timeout, material_args=()):
    """Run sphere_reentry.py for one point; mark it done/failed from the exit code and its JSON."""
    cmd = build_command(point, runs_dir, raw_dir, timeout, material_args=material_args)
    t0 = time.time()
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout + 60)
        point.returncode = proc.returncode
        stderr = proc.stderr or ""
    except subprocess.TimeoutExpired:
        point.returncode = -1
        stderr = "sweep: subprocess timed out after {} s".format(timeout + 60)
    point.wall_time_s = round(time.time() - t0, 3)
    tail = "\n".join(stderr.splitlines()[-STDERR_TAIL_LINES:])
    point.stderr_tail = tail or None
    json_ok = _json_status(os.path.join(runs_dir, point.run_name + ".json")) == "ok"
    point.status = "done" if (point.returncode == 0 and json_ok) else "failed"
    return point


def split_batches(points, batch_size):
    points = list(points)
    if not points:
        return []
    if not batch_size:
        return [points]
    return [points[i:i + batch_size] for i in range(0, len(points), batch_size)]


def default_cores():
    return max(1, (os.cpu_count() or 2) - 1)


def format_duration(seconds):
    seconds = int(round(seconds))
    hours, rest = divmod(seconds, 3600)
    minutes, secs = divmod(rest, 60)
    return "{:d}:{:02d}:{:02d}".format(hours, minutes, secs)


def confirm_batch(k, n_batches, batch, seconds_per_run, cores_for_estimate, assume_yes, ask=None):
    """Describe the batch, then ask 'Run this batch? [Y/n]' unless assume_yes."""
    ask = ask or input        # resolved at call time (tests monkeypatch builtins.input)
    estimate = len(batch) * seconds_per_run / max(1, cores_for_estimate)
    log.info("batch %d/%d: %d runs, %s .. %s, estimated %s at %d cores",
             k, n_batches, len(batch), batch[0].run_name, batch[-1].run_name,
             format_duration(estimate), cores_for_estimate)
    if assume_yes:
        return True
    while True:
        answer = ask("Run this batch? [Y/n] ").strip().lower()
        if answer in ("", "y", "yes"):
            log.info("user confirmed batch %d", k)
            return True
        if answer in ("n", "no", "q", "quit"):
            log.info("user declined batch %d", k)
            return False


def ask_cores(default, preset=None, ask=None):
    """Cores for the next batch: --cores preset, or a prompt with `default` (1..cpu_count)."""
    ask = ask or input        # resolved at call time (tests monkeypatch builtins.input)
    max_cores = os.cpu_count() or 1
    if preset is not None:
        if 1 <= preset <= max_cores:
            return preset
        raise ValueError("--cores must be between 1 and {}".format(max_cores))
    while True:
        answer = ask("Cores for this batch [{}]: ".format(default)).strip()
        if answer == "":
            return default
        try:
            cores = int(answer)
        except ValueError:
            cores = 0
        if 1 <= cores <= max_cores:
            return cores
        print("please enter an integer between 1 and {}".format(max_cores))


class TqdmHandler(logging.Handler):
    """Console log handler that writes through tqdm so the progress bar is not garbled."""

    def emit(self, record):
        try:
            tqdm.write(self.format(record))
        except Exception:  # noqa: BLE001
            self.handleError(record)


def make_progress(total, initial):
    return tqdm(total=total, initial=initial, desc="sweep", unit="run", dynamic_ncols=True)


def _run_and_mark(point, runs_dir, raw_dir, timeout, progress, runner):
    if progress is not None:
        progress.set_postfix_str("d={:g}mm T={:g}K v={:.5f}km/s".format(
            point.diameter_mm, point.temperature_K, point.velocity_kms), refresh=False)
    return runner(point, runs_dir, raw_dir, timeout)


def run_batch(batch, cores, runs_dir, raw_dir, timeout, progress=None, runner=None):
    """Run one batch on `cores` worker threads. Returns (n_ok, n_failed, wall_seconds).

    KeyboardInterrupt cancels queued points (they stay pending), lets the running
    subprocesses finish, then propagates.
    """
    runner = runner or run_point
    t0 = time.time()
    n_ok = n_failed = 0
    executor = ThreadPoolExecutor(max_workers=cores)
    try:
        futures = [executor.submit(_run_and_mark, p, runs_dir, raw_dir, timeout, progress, runner)
                   for p in batch]
        for future in as_completed(futures):
            point = future.result()
            if point.status == "done":
                n_ok += 1
            else:
                n_failed += 1
                log.warning("FAILED %s (rc=%s): %s", point.run_name, point.returncode, point.stderr_tail)
            if progress is not None:
                progress.update(1)
    except KeyboardInterrupt:
        executor.shutdown(wait=True, cancel_futures=True)
        raise
    executor.shutdown(wait=True)
    return n_ok, n_failed, time.time() - t0


# =============================================================================
# 6. Command line (spec 5.1, 5.5, 5.7)
# =============================================================================

def build_parser():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--parent", default=DEFAULT_PARENT, help="parent dmf_output.json (default %(default)s)")
    p.add_argument("--parent-object", default=DEFAULT_PARENT_OBJECT, help="object whose trajectory is used")
    p.add_argument("--outdir", default=DEFAULT_OUTDIR, help="sweep output directory (default %(default)s)")
    p.add_argument("--batch-size", type=int, default=None, help="runs per batch (default: all pending runs in one batch)")
    p.add_argument("--cores", type=int, default=None, help="cores for every batch (default: ask before each batch)")
    p.add_argument("--yes", action="store_true", help="do not ask for confirmation before a batch")
    p.add_argument("--diameters", type=parse_float_list, default=None, help="comma-separated diameters [mm]")
    p.add_argument("--temperatures", type=parse_float_list, default=None, help="comma-separated temperatures [K]")
    p.add_argument("--velocities", type=parse_float_list, default=None, help="comma-separated velocities [km/s]")
    p.add_argument("--limit", type=int, default=None, help="only the first N points of the matrix")
    p.add_argument("--dry-run", action="store_true", help="show the matrix and batch plan; run nothing")
    p.add_argument("--force", action="store_true", help="re-run every point, even completed ones")
    p.add_argument("--retry-failed", action="store_true", help="re-run points that failed before")
    p.add_argument("--timeout", type=int, default=DEFAULT_TIMEOUT_S, help="per-run SESAM timeout [s]")
    material = p.add_mutually_exclusive_group()
    material.add_argument("--material", default=sr.MATERIAL_NAME,
                          help="a metal from DRAMA's material database for every sphere (default %(default)s)")
    material.add_argument("--material-file", default=None,
                          help="JSON file defining a custom metal (DRAMA material format) for every sphere")
    return p


def summary_slug(material_name):
    """Filename part for a material: 'drama-AA7075' -> 'AA7075', 'user-moltenAA7075' unchanged."""
    name = material_name[len("drama-"):] if material_name.startswith("drama-") else material_name
    return sr.material_slug(name)


def manifest_path(outdir, material_name):
    return os.path.join(outdir, "sweep_manifest_{}.json".format(summary_slug(material_name)))


def summary_path(outdir, material_name):
    return os.path.join(outdir, "sweep_summary_{}.csv".format(summary_slug(material_name)))


def material_settings(material, material_file):
    """The manifest's record of the sweep's material (name, source, numbers, file + hash)."""
    return {
        "name": material.name,
        "source": material.source,
        "density_kgm3": material.density_kgm3,
        "melting_temperature_K": material.melting_temperature_K,
        "file": os.path.abspath(material_file) if material_file else None,
        "file_sha256": sha256_of_file(material_file) if material_file else None,
    }


def setup_logging(outdir, dry_run):
    log.setLevel(logging.INFO)
    for handler in list(log.handlers):
        handler.close()
    log.handlers.clear()
    log.propagate = True
    fmt = logging.Formatter("%(asctime)s %(levelname)s %(message)s")
    console = TqdmHandler()
    console.setFormatter(fmt)
    log.addHandler(console)
    if not dry_run:
        os.makedirs(outdir, exist_ok=True)
        file_handler = logging.FileHandler(os.path.join(outdir, "sweep.log"))
        file_handler.setFormatter(fmt)
        log.addHandler(file_handler)


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    max_cores = os.cpu_count() or 1
    if args.batch_size is not None and args.batch_size < 1:
        parser.error("--batch-size must be >= 1")
    if args.limit is not None and args.limit < 1:
        parser.error("--limit must be >= 1")
    if args.cores is not None and not 1 <= args.cores <= max_cores:
        parser.error("--cores must be between 1 and {}".format(max_cores))

    runs_dir = os.path.join(args.outdir, "runs")
    raw_dir = os.path.join(args.outdir, "raw")
    setup_logging(args.outdir, args.dry_run)
    log.info("sphere_sweep start: parent=%s object=%r outdir=%s batch_size=%s cores=%s yes=%s force=%s retry_failed=%s",
             args.parent, args.parent_object, args.outdir, args.batch_size, args.cores, args.yes,
             args.force, args.retry_failed)

    # -- parent trajectory and the velocity grid's lower bound ------------------------
    try:
        parent = load_parent(args.parent, args.parent_object)
    except (OSError, ValueError, KeyError) as exc:
        log.error("cannot load parent trajectory from %s: %s", args.parent, exc)
        return 2
    log.info("parent %s (sha256 %s) object %r epoch %s; descending branch rows %d..%d, v_max %.6f km/s",
             parent.path, parent.sha256[:12], parent.object_name,
             parent.epoch.strftime(EPOCH_ARG_FORMAT), parent.branch_start, len(parent.rows) - 1, parent.v_max)
    log.info("%s", parent.describe_v_min())

    try:
        material = sr.resolve_material(args.material, args.material_file)
    except sr.MaterialError as exc:
        log.error("cannot use material: %s", exc)
        return 2
    material_args = material_cli_args(args.material, args.material_file)
    log.info("material %s (%s): density %.1f kg/m3, melting temperature %.1f K%s",
             material.name, material.source, material.density_kgm3, material.melting_temperature_K,
             " from {}".format(os.path.abspath(args.material_file)) if args.material_file else "")

    diameters = args.diameters or DIAMETERS_MM
    temperatures = args.temperatures or TEMPERATURES_K
    velocities = args.velocities or velocity_grid(parent.v_min)
    if args.velocities is None:
        log.info("velocity grid: linspace(%.6f, %.6f, %d) km/s (lower bound = parent v_min)",
                 V_TOP_KMS, parent.v_min, N_VELOCITIES)
    grid = {"diameters_mm": diameters, "temperatures_K": temperatures, "velocities_kms": velocities}

    # -- matrix and resume ------------------------------------------------------------
    points = build_points(parent, diameters, temperatures, velocities, limit=args.limit, material=material)
    if not args.dry_run:
        os.makedirs(runs_dir, exist_ok=True)
        os.makedirs(raw_dir, exist_ok=True)
    apply_resume(points, runs_dir, force=args.force, retry_failed=args.retry_failed)
    skipped = [p for p in points if p.status == "skipped"]
    pending = [p for p in points if p.status == "pending"]
    n_done = sum(p.status == "done" for p in points)
    n_failed_before = sum(p.status == "failed" for p in points)
    batches = split_batches(pending, args.batch_size)
    log.info("matrix: %d points = %d diameters x %d temperatures x %d velocities%s; "
             "%d skipped, %d done, %d failed (not retried), %d pending in %d batch(es)",
             len(points), len(diameters), len(temperatures), len(velocities),
             " (limited to {})".format(args.limit) if args.limit else "",
             len(skipped), n_done, n_failed_before, len(pending), len(batches))
    for velocity in sorted({p.velocity_kms for p in skipped}, reverse=True):
        count = sum(p.velocity_kms == velocity for p in skipped)
        log.info("skipped v=%.6f km/s for %d point(s): parent never reaches velocity (range %.6f..%.6f)",
                 velocity, count, parent.v_min, parent.v_max)

    if args.dry_run:
        for k, batch in enumerate(batches, 1):
            log.info("batch %d/%d: %d runs, %s .. %s", k, len(batches), len(batch),
                     batch[0].run_name, batch[-1].run_name)
        log.info("dry run: nothing executed, nothing written")
        return 0

    manifest_file = manifest_path(args.outdir, material.name)
    summary_file = summary_path(args.outdir, material.name)
    settings = {"timeout_s": args.timeout, "batch_size": args.batch_size,
                "sphere_reentry_version": sr.SCRIPT_VERSION,
                "material": material_settings(material, args.material_file)}
    created = utc_now_str()

    def save():
        write_manifest(manifest_file, manifest_document(parent, grid, settings, points, created))
        write_summary(summary_file, points, runs_dir)

    save()
    if not pending:
        log.info("nothing to do: %d done, %d failed, %d skipped (use --retry-failed / --force to re-run)",
                 n_done, n_failed_before, len(skipped))
        return 0

    # -- batches ------------------------------------------------------------------------
    completed = 0
    total_ok = total_failed = 0
    total_wall = 0.0
    seconds_per_run = SECONDS_PER_RUN_ESTIMATE      # per run per core
    exit_code = 0
    try:
        for k, batch in enumerate(batches, 1):
            if not confirm_batch(k, len(batches), batch, seconds_per_run,
                                 args.cores or default_cores(), args.yes):
                log.info("stopped before batch %d by user", k)
                break
            cores = ask_cores(default_cores(), preset=args.cores)
            log.info("batch %d/%d: %d runs on %d cores", k, len(batches), len(batch), cores)
            progress = make_progress(total=len(pending), initial=completed)
            try:
                n_ok, n_failed, wall = run_batch(batch, cores, runs_dir, raw_dir, args.timeout, progress,
                                                 runner=functools.partial(run_point, material_args=material_args))
            finally:
                progress.close()
            completed += n_ok + n_failed
            total_ok += n_ok
            total_failed += n_failed
            total_wall += wall
            if n_ok + n_failed:
                seconds_per_run = wall * cores / (n_ok + n_failed)
            log.info("batch %d/%d done: %d ok, %d failed, %s (%.2f s/run/core)",
                     k, len(batches), n_ok, n_failed, format_duration(wall), seconds_per_run)
            save()
    except KeyboardInterrupt:
        log.warning("interrupted: saving manifest and summary; re-run to resume")
        exit_code = 130
    save()
    still_pending = sum(p.status == "pending" for p in points)
    log.info("sweep finished: %d ok, %d failed, %d pending, total run time %s",
             total_ok, total_failed, still_pending, format_duration(total_wall))
    return exit_code


if __name__ == "__main__":
    sys.exit(main())

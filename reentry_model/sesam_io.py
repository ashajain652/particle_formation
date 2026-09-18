"""A SESAM run written by sphere_reentry.py (history CSV + run JSON) as a Reference in SI units."""
import csv
import hashlib
import json
import math
import os
from dataclasses import dataclass
from datetime import datetime

import numpy as np

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REFERENCE_DIR = os.path.join(REPO_ROOT, "data", "reference_runs")


@dataclass(frozen=True)
class InitialStateSpec:
    velocity: float      # m/s, relative to the rotating atmosphere
    altitude: float      # m, geodetic
    flight_path: float   # rad, negative = descending
    heading: float       # rad, clockwise from north
    lat: float           # rad, geodetic
    lon: float           # rad
    epoch: datetime


@dataclass
class Reference:
    name: str
    csv_path: str
    json_path: str
    sha256: str
    diameter: float
    material_density: float
    atmosphere: str
    use_wind: bool
    initial: InitialStateSpec
    time: np.ndarray
    altitude: np.ndarray
    velocity: np.ndarray
    lat: np.ndarray
    lon: np.ndarray
    flight_path: np.ndarray
    heading: np.ndarray
    downrange: np.ndarray
    drag: np.ndarray
    knudsen: np.ndarray
    mach: np.ndarray
    density: np.ndarray
    dynamic_pressure: np.ndarray
    temperature: np.ndarray
    convective_heat: np.ndarray = None      # W (SESAM's convective_heat_W), None if the CSV lacks the column
    rad_cooling: np.ndarray = None          # W, negative (SESAM's rad_cooling_W)
    integrated_heat: np.ndarray = None      # J, integral of the convective heat (SESAM's integrated_heat_J)


def _sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def load_reference(csv_path, json_path=None):
    """Read <run>.csv and <run>.json (default: same stem) into a Reference."""
    if json_path is None:
        json_path = os.path.splitext(csv_path)[0] + ".json"
    if not os.path.isfile(json_path):
        raise FileNotFoundError("reference JSON not found next to {}: {}".format(csv_path, json_path))
    with open(json_path) as fh:
        doc = json.load(fh)
    inp = doc["inputs"]
    initial = InitialStateSpec(
        velocity=float(inp["initial_velocity_kms"]) * 1e3,
        altitude=float(inp["initial_altitude_km"]) * 1e3,
        flight_path=math.radians(float(inp["flight_path_angle_deg"])),
        heading=math.radians(float(inp["heading_deg"])),
        lat=math.radians(float(inp["latitude_deg"])),
        lon=math.radians(float(inp["longitude_deg"])),
        epoch=datetime.strptime(inp["epoch_utc"], "%Y-%m-%dT%H:%M:%S"),
    )
    with open(csv_path) as fh:
        rows = list(csv.DictReader(fh))
    col = lambda key, scale=1.0: np.array([float(r[key]) for r in rows]) * scale
    optional = lambda key: col(key) if key in rows[0] else None
    return Reference(
        name=doc["run_name"], csv_path=os.path.abspath(csv_path), json_path=os.path.abspath(json_path),
        sha256=_sha256(csv_path),
        diameter=float(inp["diameter_mm"]) * 1e-3, material_density=float(inp["material_density_kgm3"]),
        atmosphere=inp.get("atmosphere", "static"), use_wind=bool(inp.get("use_wind", False)),
        initial=initial,
        time=col("time_s"), altitude=col("altitude_km", 1e3), velocity=col("velocity_kms", 1e3),
        lat=np.radians(col("lat_deg")), lon=np.radians(col("lon_deg")),
        flight_path=np.radians(col("flight_path_deg")), heading=np.radians(col("heading_deg")),
        downrange=col("downrange_km", 1e3), drag=col("drag"), knudsen=col("knudsen"), mach=col("mach"),
        density=col("density_kgm3"), dynamic_pressure=col("dynamic_pressure_Pa"), temperature=col("temperature_K"),
        convective_heat=optional("convective_heat_W"), rad_cooling=optional("rad_cooling_W"), integrated_heat=optional("integrated_heat_J"),
    )

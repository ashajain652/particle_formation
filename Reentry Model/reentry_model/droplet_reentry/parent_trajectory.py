"""The parent satellite's DRAMA trajectory, from which the droplets are released.

Copied 2026-10-26 from DRAMA Sphere Experiments/sphere_sweep.py (section 1, "Parent trajectory"): load_parent, the
ParentTrajectory it returns, and the interpolation helpers, unchanged. The rest of sphere_sweep.py drives SESAM
runs and is not needed here.

DEFAULT_PARENT is the committed snapshot of the parent's DRAMA GUI output (Generic_Satellite Reentry/output/
dmf_output.json, run 2026-09-13), kept in Reentry Model/data/ because the DRAMA GUI rewrites the original folder on
every run.
"""
import hashlib
import json
import math
import os
from dataclasses import dataclass
from datetime import datetime, timedelta  # noqa: F401  (timedelta: used by the callers through this module)

from reentry_model import DATA_DIR

DEFAULT_PARENT = os.path.join(DATA_DIR, "generic_satellite_reentry", "dmf_output.json")
DEFAULT_PARENT_OBJECT = "Main Body"
PARENT_EPOCH_FORMAT = "%Y-%m-%dT%H:%M:%S.%fZ"


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


def load_parent(path=DEFAULT_PARENT, object_name=DEFAULT_PARENT_OBJECT):
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

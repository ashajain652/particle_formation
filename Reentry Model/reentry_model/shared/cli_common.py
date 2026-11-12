"""The run command's architecture-independent parts, shared by both architectures' run drivers: the flight's
defaults and arguments (diameter through name), the atmosphere factory, the git commit and provenance record, and the
--architecture pre-parse that routes a run to its architecture. Moved unchanged from cli.py (combined refactor, Task 6)."""
import argparse
import importlib.metadata
import subprocess
import sys
from datetime import datetime

import numpy as np

from . import aero, atmosphere, fap, sesam_io
from .earth import GRAVITY_MODELS
from .. import OUTPUT_DIR, REPO_ROOT, __version__

DEFAULT_OUTDIR = OUTPUT_DIR
# the reference cases' break-off state (parent satellite, see README "Reference runs")
DEFAULT_HEADING_DEG = 347.168296
DEFAULT_LAT_DEG = 29.546067
DEFAULT_LON_DEG = -82.134333
DEFAULT_EPOCH = "2024-08-01T12:53:07"
DEFAULT_MATERIAL_DENSITY = 2813.0          # drama-AA7075
ATMOSPHERES = ("nrlmsise", "us76")         # plus replay:<path>
WINDS = ("none", "static")


def parse_epoch(text):
    try:
        return datetime.strptime(text, "%Y-%m-%dT%H:%M:%S")
    except ValueError:
        raise argparse.ArgumentTypeError("epoch must be YYYY-MM-DDTHH:MM:SS, got {!r}".format(text))



def make_atmosphere(spec, epoch, wind_name):
    """(atmosphere object, short name, provenance dict) for an --atmosphere value."""
    wind = atmosphere.NoWind() if wind_name == "none" else atmosphere.StaticProfileWind()
    if spec == "us76":
        return atmosphere.US76TableAtmosphere(wind=wind), "us76", {}
    if spec == "nrlmsise":
        solar = fap.solar_indices(fap.load_fap_day(fap.DEFAULT_FAP_DAY), epoch.date())
        return atmosphere.NRLMSISE00Atmosphere(epoch, solar, wind), "nrlmsise", \
            {"f107": solar.f107, "f107a": solar.f107a, "ap": solar.ap, "fap_day": fap.DEFAULT_FAP_DAY}
    if spec.startswith("replay:"):
        ref = sesam_io.load_reference(spec[len("replay:"):])
        return atmosphere.ReplayAtmosphere(ref, wind), "replay", {"replay_reference": ref.name, "replay_sha256": ref.sha256}
    raise ValueError("--atmosphere must be one of {} or replay:<sesam.csv>, got {!r}".format(ATMOSPHERES, spec))


def git_commit():
    try:
        return subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=REPO_ROOT, stderr=subprocess.DEVNULL).decode().strip()
    except Exception:
        return "unknown"


def provenance():
    import scipy
    try:
        import pymsis
        pymsis_version = pymsis.__version__
    except ImportError:
        pymsis_version = None
    package_names = {"skfem": "scikit-fem", "gmsh": "gmsh", "pyamg": "pyamg", "pyvista": "pyvista",
                      "cantera": "cantera", "dolfinx": "fenics-dolfinx"}
    versions = {}
    for key, pkg in package_names.items():
        try:
            versions[key] = importlib.metadata.version(pkg)
        except importlib.metadata.PackageNotFoundError:
            versions[key] = None
    return {"package_version": __version__, "git_commit": git_commit(), "numpy": np.__version__,
            "scipy": scipy.__version__, "pymsis": pymsis_version, "python": sys.version.split()[0], **versions}


def add_flight_arguments(parser):
    """The flight, atmosphere and output arguments every architecture's `run` takes, from --diameter through --name."""
    parser.add_argument("--diameter", type=float, required=True, help="sphere diameter [mm]")
    parser.add_argument("--velocity", type=float, required=True, help="initial velocity [km/s], relative to the rotating atmosphere")
    parser.add_argument("--altitude", type=float, required=True, help="initial geodetic altitude [km]")
    parser.add_argument("--flight-path-angle", type=float, default=0.0, help="[deg], negative = descending")
    parser.add_argument("--heading", type=float, default=DEFAULT_HEADING_DEG, help="[deg] clockwise from north (default %(default)s)")
    parser.add_argument("--lat", type=float, default=DEFAULT_LAT_DEG, help="geodetic latitude [deg] (default %(default)s)")
    parser.add_argument("--lon", type=float, default=DEFAULT_LON_DEG, help="longitude [deg] (default %(default)s)")
    parser.add_argument("--epoch", type=parse_epoch, default=parse_epoch(DEFAULT_EPOCH), help="UTC YYYY-MM-DDTHH:MM:SS (default {})".format(DEFAULT_EPOCH))
    parser.add_argument("--material-density", type=float, default=DEFAULT_MATERIAL_DENSITY, help="[kg/m3] (default %(default)s)")
    parser.add_argument("--temperature", type=float, default=300.0,
                        help="initial temperature [K]: recorded only with --thermal none, the field's initial condition with --thermal fem (default %(default)s)")
    parser.add_argument("--atmosphere", default="nrlmsise", help="nrlmsise (default) | us76 | replay:<sesam.csv>")
    parser.add_argument("--wind", choices=WINDS, default="none")
    parser.add_argument("--bridging", choices=aero.BRIDGING_NAMES, default="sesam-table")
    parser.add_argument("--gravity", choices=GRAVITY_MODELS, default="j2")
    parser.add_argument("--rtol", type=float, default=1e-9)
    parser.add_argument("--cadence", type=float, default=1.0, help="history sample spacing [s] (default %(default)s)")
    parser.add_argument("--t-max", type=float, default=3600.0, help="[s] (default %(default)s)")
    parser.add_argument("--reference", default=None,
                        help="SESAM run CSV to compare against (with --thermal none the model is also sampled at its times; "
                             "with --thermal fem the macro-step history is interpolated)")
    parser.add_argument("--outdir", default=DEFAULT_OUTDIR)
    parser.add_argument("--name", default=None, help="run name (default: model_d..mm_v..kms_h..km_<atmosphere>_<bridging>_<wind>[_fem-<heating>])")


ARCHITECTURES = ("skyscrapers", "levelset")
DEFAULT_ARCHITECTURE = "skyscrapers"


def architecture_of(run_argv):
    """The --architecture value in a `run` command's arguments (the arguments after `run`): `--architecture X` or
    `--architecture=X`, default DEFAULT_ARCHITECTURE. An unknown or missing value ends the program with exit code 2 and the
    `run` usage, as the full parser would."""
    pre = argparse.ArgumentParser(prog="reentry_model run", add_help=False)
    pre.add_argument("--architecture", choices=ARCHITECTURES, default=DEFAULT_ARCHITECTURE)
    known, _ = pre.parse_known_args(run_argv)
    return known.architecture

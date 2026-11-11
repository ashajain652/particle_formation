"""Command line of the re-entry model.

    python -m reentry_model run --diameter 100 --velocity 7.5 --altitude 77.500133 --flight-path-angle -0.959331 \
        [--atmosphere nrlmsise|us76|replay:<sesam.csv>] [--reference <sesam.csv>] [--outdir ...]
        [--thermal fem --heating sesam|physics ... --animate]          (Step 2: coupled 3D conduction)
        [--melt on --material AA7075_scheil --removal girin|instant ...] (Step 3: melting, film, spraying)
        [--surface-model skin|elements --skin-flow on|off --slurry-spray on|off ...]
                                     (the layered skins, the default with --thermal fem, their flow and spraying the
                                      default with --melt on; the element model with --surface-model elements)
    python -m reentry_model compare --model <model.csv> --reference <sesam.csv> [--outdir ...]

Exit codes: 0 ok, 1 the flight escaped / integration failed, 2 bad input or a missing optional library
(cantera for --heating physics, dolfinx for --thermal-solver fenicsx: create the fenicsx_env environment).
"""
import argparse
import importlib.metadata
import math
import os
import subprocess
import sys
from datetime import datetime

import numpy as np

from . import __version__, aero, atmosphere, body, compare, coupled, dispersion, fap, heating, material, mesh, movemesh, sesam_io, skin, skinbody, spray, surface_flow, thermal, viz
from . import trajectory as tj
from .earth import GRAVITY_MODELS
from .thermal import MissingBackend

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DEFAULT_OUTDIR = os.path.join(REPO_ROOT, "reentry_model_output")
# the reference cases' break-off state (parent satellite, see README "Reference runs")
DEFAULT_HEADING_DEG = 347.168296
DEFAULT_LAT_DEG = 29.546067
DEFAULT_LON_DEG = -82.134333
DEFAULT_EPOCH = "2024-08-01T12:53:07"
DEFAULT_MATERIAL_DENSITY = 2813.0          # drama-AA7075
ATMOSPHERES = ("nrlmsise", "us76")         # plus replay:<path>
WINDS = ("none", "static")
THERMAL_MODES = ("none", "fem")
MELT_MODES = ("off", "on")
# numpy's global generator is seeded with this at the start of every run (amendment of 2026-10-05): pyamg draws from it
# the starting vector of the spectral-radius estimate that weights the prolongation smoother of every hierarchy it builds,
# so an unseeded run differed from the next from the first solve on. 12345 is the seed the measurement harness of the
# amendments of 2026-10-02 and 2026-10-03 set, so the model's default reproduces their runs.
DEFAULT_SEED = 12345
# The macro step once the flight enters the continuum regime, for the model proper (amendment of 2026-10-07): from the
# first step whose body Knudsen number is below --kn-body-shock, where Girin's closure begins, the step shortens to this,
# latched. The switched-step series converges from it -- the masses and the branch split within the scatter between
# seeds, the droplet count within 2.7 % of the 0.00625 s step's (plan facts 92-94) -- while the default 0.5 s step makes
# six times fewer droplets. Before the switch the default step serves.
DEFAULT_DT_CONTINUUM = 0.0125
# The surface models a run can choose (sub-plan 20): the layered skins on the corners of a moving mesh
# (skinbody.SkinBody), the default of a thermal run since sub-plan 22 Task 3 gave them Girin's spraying (spec 10.1), or
# the element model, kept unchanged for comparison (spec 9 and 11.1 level 3), whose runs end in `_surface-elements`.
# Sub-plan 18's one-dimensional skins are retired (sub-plan 22 Task 3).
SURFACE_MODELS = ("skin", "elements")
DEFAULT_SURFACE_MODEL = "skin"
# The fine band under a skin run's surface (Asha, 2026-10-14, option 1): the 3D model's elements are finer to this depth
# while the surface keeps its triangulation, so that the 3D model follows the temperature just below the rigid front
# (on 2 mm elements throughout the plates melted 2-6 % late), and every rebuild makes the band again
# (movemesh.refined_mesh). Under a fixed 2 mm surface gmsh fills the band coarser than asked: 0.7 mm asked matches 1 mm
# elements throughout on both plates' onset, 1 mm asked gets about half way (measured 2026-10-15); on the 100 mm sphere
# it takes the 3D model from 177 000 to 1.21 million elements and a 3D step from 0.16 to 1.17 s.
FINE_BAND_MM, FINE_SIZE_MM = 4.0, 0.7


def parse_epoch(text):
    try:
        return datetime.strptime(text, "%Y-%m-%dT%H:%M:%S")
    except ValueError:
        raise argparse.ArgumentTypeError("epoch must be YYYY-MM-DDTHH:MM:SS, got {!r}".format(text))


def parse_seed(text):
    """--seed: an integer in [0, 2**32 - 1], the range numpy's np.random.seed accepts."""
    try:
        seed = int(text)
    except ValueError:
        seed = -1
    if not 0 <= seed < 2 ** 32:
        raise argparse.ArgumentTypeError("seed must be an integer in [0, 2**32 - 1], got {!r}".format(text))
    return seed


def parse_dt_continuum(text):
    """--dt-continuum: `off`, or a macro step [s] > 0."""
    if text == "off":
        return text
    try:
        value = float(text)
    except ValueError:
        value = -1.0
    if not value > 0.0:
        raise argparse.ArgumentTypeError("dt-continuum must be off or a step > 0 s, got {!r}".format(text))
    return value


def continuum_step(args):
    """The macro step after the continuum switch [s], or None for no switch: --dt-continuum if given (`off`: none), else
    DEFAULT_DT_CONTINUUM for the model proper -- melting with Girin removal -- and none for the instant-removal device
    (no film; its SESAM thresholds were measured at the default step) or a run that does not melt."""
    if args.dt_continuum is not None:
        return None if args.dt_continuum == "off" else args.dt_continuum
    return DEFAULT_DT_CONTINUUM if (args.thermal == "fem" and args.melt == "on" and removal_of(args) == "girin") else None


def surface_model_of(args):
    """The run's surface model: --surface-model if given, else the layered skins for a thermal run (spec 10.1; sub-plan
    22 Task 3); None without a thermal model, which has no surface to model."""
    if args.surface_model is not None:
        return args.surface_model
    return DEFAULT_SURFACE_MODEL if args.thermal == "fem" else None


def skin_flow_of(args):
    """--skin-flow if given, else the run's own (sub-plan 22 Task 3b; Asha, 2026-10-29): on for the skins' melting runs,
    which flow and spray by Girin's model, but off where --removal instant is given, the heat-only skins' removal; off
    for a run that does not melt or does not use the skins, where the flow does not apply."""
    if args.skin_flow is not None:
        return args.skin_flow
    melting_skins = surface_model_of(args) == "skin" and args.thermal == "fem" and args.melt == "on"
    return "on" if melting_skins and args.removal != "instant" else "off"


def removal_of(args):
    """--removal if given, else the run's own (sub-plan 22 Task 3): Girin's, for the element model and for the skins with
    the flow on, whose liquid only spraying removes (spec 10.1: --skin-flow on needs --removal girin); instant removal for
    the heat-only skins, whose liquid leaves as it forms (stage 2's stand-in for flow and spray, spec 11.2)."""
    if args.removal is not None:
        return args.removal
    return "instant" if surface_model_of(args) == "skin" and skin_flow_of(args) == "off" else "girin"


def model_run_name(diameter_m, velocity_ms, altitude_m, atmosphere_name, bridging_name, wind_name, heating_name=None, melt=None,
                   deep_runoff=True, molten_cascade=True, seed=DEFAULT_SEED, rigid_substrate=True, dt_continuum="default",
                   surface_model=DEFAULT_SURFACE_MODEL, skin_suffix=""):
    """The run name encodes the configuration, so runs of different settings can share an output directory. A melting
    run with the deep runoff off (amendment of 2026-10-02) ends in `_deeprunoff-off`, one with the molten cascade off
    (amendment of 2026-10-03) in `_moltencascade-off`, one with the rigid substrate off (amendment of 2026-10-06) in
    `_rigidsubstrate-off`, one whose continuum step is not its removal mode's default (amendment of 2026-10-07) in
    `_dtcontinuum-off` or `_dtcontinuum-<s>`; a thermal run with the element model in `_surface-elements`, the layered
    skins being the default (spec 10.1; sub-plan 22 Task 3), and one with the skins in `skin_suffix`, their non-default
    settings (`skin_name_suffix`, empty for the defaults); and any run with a seed other than DEFAULT_SEED (amendment of
    2026-10-05) in `_seed-<n>`. Without a thermal model there is no surface model and no suffix for it."""
    dtc_default = DEFAULT_DT_CONTINUUM if melt == "girin" else None
    dtc = dtc_default if dt_continuum == "default" else dt_continuum
    name = "model_d{:06.2f}mm_v{:08.5f}kms_h{:07.3f}km_{}_{}_{}".format(
        diameter_m * 1e3, velocity_ms / 1e3, altitude_m / 1e3, atmosphere_name, bridging_name, wind_name)
    return (name + ("_fem-" + heating_name if heating_name else "") + ("_melt-" + melt if melt else "")
            + ("_deeprunoff-off" if melt and not deep_runoff else "") + ("_moltencascade-off" if melt and not molten_cascade else "")
            + ("_rigidsubstrate-off" if melt and not rigid_substrate else "")
            + ("" if not melt or dtc == dtc_default else "_dtcontinuum-off" if dtc is None else "_dtcontinuum-{:g}".format(dtc))
            + (("_surface-elements" if surface_model == "elements" else skin_suffix) if heating_name else "")
            + ("_seed-{}".format(seed) if seed != DEFAULT_SEED else ""))


def skin_name_suffix(min_thickness_mm, cell_um, growth, substep_ms, preset, deform_mm, swell=True,
                     fine_band_mm=FINE_BAND_MM, fine_size_mm=FINE_SIZE_MM, fine_zone_mm=0.4, flow=True, flow_cycles=None,
                     slurry_spray=False):
    """The run-name suffix of a skin run's non-default settings (layered-skin design section 10.1): `_skin-dev` for the
    dev preset, `_skinmin<mm>mm` for another minimum thickness, `_skinfine<mm>mm` for another fine zone (both changed by
    option B', 2026-10-16: 0.8 and 0.4 mm), `_cell<um>um` for an explicit cell other than the
    preset's, `_growth<ratio>`, `_sub<ms>ms` and `_deform<mm>mm` for another growth, sub-step or deforming depth,
    `_noswell` for skins held at the solid's density, `_fine<mm>x<mm>mm` for another fine band (`_fine0` for none);
    `_skinflow-off` for the heat-only skins of a melting run (`flow` False: the flow, sub-plan 21 Task 6c, is on by default
    with melting since sub-plan 22 Task 3b; None for a run that does not melt, which the flow does not concern), and with
    the flow on `_flowcycles<n>` where the flow sub-steps in each heat sub-step are fixed at n rather than chosen each heat
    sub-step (Task 6d) and `_slurryspray-on` where the exposed slurry sprays too (spec 10.1; sub-plan 22 Task 3); empty
    for the defaults."""
    out = "_skin-dev" if preset == "dev" else ""
    if min_thickness_mm != 0.8:
        out += "_skinmin{:g}mm".format(min_thickness_mm)
    if fine_zone_mm != 0.4:
        out += "_skinfine{:g}mm".format(fine_zone_mm)
    if cell_um is not None and abs(cell_um * 1e-6 - skin.SKIN_PRESETS[preset]) > 1e-12:
        out += "_cell{:g}um".format(cell_um)
    if growth != 1.10:
        out += "_growth{:g}".format(growth)
    if substep_ms != 10.0:
        out += "_sub{:g}ms".format(substep_ms)
    if deform_mm != 15.0:
        out += "_deform{:g}mm".format(deform_mm)
    if not swell:
        out += "_noswell"
    if fine_band_mm <= 0.0:
        out += "_fine0"
    elif (fine_band_mm, fine_size_mm) != (FINE_BAND_MM, FINE_SIZE_MM):
        out += "_fine{:g}x{:g}mm".format(fine_band_mm, fine_size_mm)
    if flow is False:
        out += "_skinflow-off"
    if flow_cycles is not None:
        out += "_flowcycles{}".format(flow_cycles)
    if slurry_spray:
        out += "_slurryspray-on"
    return out


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


def build_parser():
    p = argparse.ArgumentParser(prog="reentry_model", description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="command", required=True)

    r = sub.add_parser("run", help="integrate one sphere trajectory")
    r.add_argument("--diameter", type=float, required=True, help="sphere diameter [mm]")
    r.add_argument("--velocity", type=float, required=True, help="initial velocity [km/s], relative to the rotating atmosphere")
    r.add_argument("--altitude", type=float, required=True, help="initial geodetic altitude [km]")
    r.add_argument("--flight-path-angle", type=float, default=0.0, help="[deg], negative = descending")
    r.add_argument("--heading", type=float, default=DEFAULT_HEADING_DEG, help="[deg] clockwise from north (default %(default)s)")
    r.add_argument("--lat", type=float, default=DEFAULT_LAT_DEG, help="geodetic latitude [deg] (default %(default)s)")
    r.add_argument("--lon", type=float, default=DEFAULT_LON_DEG, help="longitude [deg] (default %(default)s)")
    r.add_argument("--epoch", type=parse_epoch, default=parse_epoch(DEFAULT_EPOCH), help="UTC YYYY-MM-DDTHH:MM:SS (default {})".format(DEFAULT_EPOCH))
    r.add_argument("--material-density", type=float, default=DEFAULT_MATERIAL_DENSITY, help="[kg/m3] (default %(default)s)")
    r.add_argument("--temperature", type=float, default=300.0,
                   help="initial temperature [K]: recorded only with --thermal none, the field's initial condition with --thermal fem (default %(default)s)")
    r.add_argument("--atmosphere", default="nrlmsise", help="nrlmsise (default) | us76 | replay:<sesam.csv>")
    r.add_argument("--wind", choices=WINDS, default="none")
    r.add_argument("--bridging", choices=aero.BRIDGING_NAMES, default="sesam-table")
    r.add_argument("--gravity", choices=GRAVITY_MODELS, default="j2")
    r.add_argument("--rtol", type=float, default=1e-9)
    r.add_argument("--cadence", type=float, default=1.0, help="history sample spacing [s] (default %(default)s)")
    r.add_argument("--t-max", type=float, default=3600.0, help="[s] (default %(default)s)")
    r.add_argument("--reference", default=None,
                   help="SESAM run CSV to compare against (with --thermal none the model is also sampled at its times; "
                        "with --thermal fem the macro-step history is interpolated)")
    r.add_argument("--outdir", default=DEFAULT_OUTDIR)
    r.add_argument("--name", default=None, help="run name (default: model_d..mm_v..kms_h..km_<atmosphere>_<bridging>_<wind>[_fem-<heating>])")
    r.add_argument("--seed", type=parse_seed, default=DEFAULT_SEED,
                   help="seed of numpy's global random generator, set at the start of the run (default %(default)s): pyamg draws "
                        "the starting vectors of its spectral-radius estimates from it, so a fixed seed makes a run reproducible "
                        "bit for bit and another seed changes it only through that round-off, which the melting model amplifies "
                        "into the run-to-run scatter; an integer in [0, 2**32 - 1], recorded in the run JSON, and a seed other "
                        "than the default ends the run name in _seed-<n>")
    r.add_argument("--quiet", action="store_true")
    th = r.add_argument_group("thermal model (Step 2)")
    th.add_argument("--thermal", choices=THERMAL_MODES, default="none", help="none: Step 1 trajectory only (default); fem: coupled 3D conduction")
    th.add_argument("--heating", choices=heating.HEATING_NAMES, default="physics",
                    help="physics (default) or sesam: SESAM's uniform 0.27471 x q_stag on every patch -- a verification device, not physical")
    th.add_argument("--stagnation", choices=heating.STAGNATION_NAMES, default="fay-riddell")
    th.add_argument("--bridging-heat", choices=heating.BRIDGING_HEAT_NAMES, default="matting")
    th.add_argument("--matting-n", type=float, default=1.0, help="Matting exponent n (default %(default)s)")
    th.add_argument("--accommodation", type=float, default=0.8, help="free-molecular energy accommodation A_cq (default %(default)s)")
    th.add_argument("--catalycity", type=float, default=1.0, help="wall catalycity 0..1 in Fay-Riddell (default %(default)s)")
    th.add_argument("--material", default=None,
                    help="AA7075_nomelt | AA7075 | AA7075_range | AA7075_scheil | AA7075-empiricaldata | <DRAMA material JSON> "
                         "(default: AA7075_nomelt, or AA7075_scheil with --melt on)")
    th.add_argument("--emissivity", type=float, default=None, help="override the material's emissivity")
    th.add_argument("--k-scale", type=float, default=1.0,
                    help="multiplies the material's conductivity (default 1): a verification device for the near-isothermal (lumped) limit, not physical")
    th.add_argument("--t-ambient", type=float, default=0.0, help="radiation background [K] (default %(default)s, SESAM's)")
    th.add_argument("--mesh-size", type=float, default=1.0, help="multiplies --h-surface and --h-core (default %(default)s)")
    th.add_argument("--h-surface", type=float, default=mesh.DEFAULT_H_SURFACE * 1e3, help="surface element size [mm] (default %(default)s)")
    th.add_argument("--h-core", type=float, default=mesh.DEFAULT_H_CORE * 1e3, help="core element size [mm] (default %(default)s)")
    th.add_argument("--thermal-solver", choices=thermal.SOLVER_NAMES, default="skfem")
    th.add_argument("--linear-solver", choices=("direct", "amg"), default="amg")
    th.add_argument("--consistent-mass", action="store_true", help="consistent capacity matrix instead of the lumped nodal-enthalpy one (skfem only, analytic checks)")
    th.add_argument("--dt", type=float, default=0.5, help="macro step [s] (default %(default)s)")
    th.add_argument("--frames-every", type=int, default=0, help="VTK frame every n macro steps (default 0: none; 10 with --animate/--stills)")
    th.add_argument("--animate", action="store_true", help="MP4/GIF of the surface temperature and of the meridional cross-section, plus stills")
    th.add_argument("--stills", action="store_true", help="only the stills (start, peak heating, peak surface T, end; surface and section)")
    me = r.add_argument_group("melting and spraying (Step 3)")
    me.add_argument("--melt", choices=MELT_MODES, default="off", help="off: Step 2 behaviour (default); on: melting, melt film, spraying")
    me.add_argument("--removal", choices=body.REMOVAL_NAMES, default=None,
                    help="girin: Girin's spraying -- the element model's film, runoff and spraying, or the skins' flow and "
                         "spraying from their top cells (default, except for the heat-only skins); instant: liquid removed "
                         "as it forms -- the lumped-melting verification device, not physical (the default of the "
                         "heat-only skins, --skin-flow off)")
    me.add_argument("--runoff", choices=("on", "off"), default="on", help="film runoff transport along the surface (default on)")
    me.add_argument("--deep-runoff", choices=("on", "off"), default="on",
                    help="move the contiguous liquid below Girin's conjugate depth along the surface by the pressure-gradient and "
                         "deceleration-driven lubrication flux over the whole liquid depth, never spraying it; it becomes film only "
                         "from the top, as the film thins below the conjugate depth (default on; needs --runoff on). off: that "
                         "liquid stays in its elements until the surface recedes to it, as in every run before 2026-10-02")
    me.add_argument("--molten-cascade", choices=("on", "off"), default="on",
                    help="when an element death exposes an element that is fully molten (its mean temperature at or above the top "
                         "of the feed ramp, the test the molten layer is measured by), feed it whole within the same macro step, so "
                         "that it dies in turn and the surface recedes through molten material by as many elements as are molten "
                         "(default on). off: an exposed element waits for the next step's feed, so the surface recedes through "
                         "molten material by one element per macro step, as in every run before the molten-cascade amendment of "
                         "2026-10-03")
    me.add_argument("--rigid-substrate", choices=("on", "off"), default="on",
                    help="the thin spraying branch needs a rigid substrate: the regime test compares Girin's conjugate depth "
                         "with the film plus everything more than half liquid beneath it (slurry included), and where there "
                         "is no conjugate depth a film on slurry does not spray (default on). off: the regime test reads the "
                         "film plus the fully molten material, as in every run before the amendment of 2026-10-06")
    me.add_argument("--surface-model", choices=SURFACE_MODELS, default=None,
                    help="skin: a layered skin on every corner of the surface holds the outer metal in fine cells above a "
                         "moving 3D mesh (the layered-skin design; default with --thermal fem since sub-plan 22). elements: "
                         "the surface elements melt and feed the film, as in every run before it, kept for comparison; "
                         "its runs end in _surface-elements")
    me.add_argument("--skin-min-thickness", type=float, default=0.8,
                    help="the skins' minimum thickness [mm of solid]: where the non-rigid layer is thinner, a skin's bottom "
                         "is rigid metal (default %(default)s)")
    me.add_argument("--skin-fine-zone", type=float, default=0.4,
                    help="the depth below a skin's top held in fine cells [mm]; below it each cell is --skin-growth times "
                         "thicker than the one above (default %(default)s)")
    me.add_argument("--skin-cell", type=float, default=None, help="skin fine cell [um] (default: the preset's)")
    me.add_argument("--skin-growth", type=float, default=1.10,
                    help="each skin cell below the fine zone this much thicker than the one above (default %(default)s)")
    me.add_argument("--skin-substep", type=float, default=10.0, help="skin sub-step [ms] (default %(default)s)")
    me.add_argument("--fine-band", type=float, default=FINE_BAND_MM,
                    help="with --surface-model skin: the 3D model's elements are --fine-size to this depth [mm] below the "
                         "surface, which keeps its triangulation, and every rebuild makes the band again; 0: none "
                         "(default %(default)s)")
    me.add_argument("--fine-size", type=float, default=FINE_SIZE_MM,
                    help="the elements' size asked for in the fine band [mm] (default %(default)s; capped at the surface's)")
    me.add_argument("--skin-swell", choices=("on", "off"), default="on",
                    help="on: each skin cell at its own phase's density, so melting swells the skin (default). off: every "
                         "cell at the solid's density, the verification setting of levels 2 and 2b, whose references hold "
                         "one density")
    me.add_argument("--skin-flow", choices=("on", "off"), default=None,
                    help="with --surface-model skin and --melt on: on, the skins' liquid and slurry flow over the surface "
                         "(sub-plan 21) and Girin's model sprays them from the top cells (sub-plan 22), with --removal "
                         "girin, whose continuum step (--dt-continuum) follows (default, unless --removal instant is "
                         "given); off: the heat-only skins, liquid removed as it forms with --removal instant, named "
                         "_skinflow-off")
    me.add_argument("--slurry-spray", choices=("on", "off"), default="off",
                    help="with --skin-flow on: on, the exposed slurry sprays too, once no liquid covers it, by its own mode "
                         "(the cascade of the layered-skin design, section 8); it sprays only within about a kelvin of the "
                         "liquidus, on a viscosity law with no data there (default off)")
    me.add_argument("--skin-flow-cycles", type=int, default=None,
                    help="with --skin-flow on: flow sub-steps in each skin sub-step, fixed (default: chosen each skin "
                         "sub-step from the flow's Courant number and growth, {} to {})".format(skin.FLOW_CYCLES_MIN,
                                                                                            skin.FLOW_CYCLES_MAX))
    me.add_argument("--skin-preset", choices=tuple(skin.SKIN_PRESETS), default="production",
                    help="skin cell-size preset: production 10 um, dev 20 um (default %(default)s)")
    me.add_argument("--deform-depth", type=float, default=15.0,
                    help="depth of the moving mesh's deforming zone [mm] (default %(default)s, the dense band's)")
    me.add_argument("--rt-spray", choices=("on", "off"), default="on",
                    help="apply Girin & Kopyt's front-surface Rayleigh-Taylor mode as a release mechanism, on the patches where it "
                         "grows faster than the shear mode (default on). off: evaluate and report it but release nothing through it, "
                         "which is how every run before 2026-09-24 behaved")
    me.add_argument("--rarefied-shear", choices=surface_flow.RAREFIED_SHEAR_NAMES, default="slip",
                    help="within the shock-layer branch: slip (default) uses Maxwell slip below Kn_local 0.1 and free-molecular shear above; "
                         "bridged blends them with SESAM's f(Kn). The merged branch always bridges")
    me.add_argument("--gamma-pm", type=float, default=surface_flow.GAMMA_PM,
                    help="effective ratio of specific heats of the Prandtl-Meyer expansion that sets the wall pressure beyond the sonic "
                         "point (default %(default)s; 1.4 frozen air, ~1.15 dissociated -- a factor ~2 on p_w at 90 deg)")
    me.add_argument("--dt-continuum", type=parse_dt_continuum, default=None,
                    help="macro step [s] from the first step whose body Knudsen number is below --kn-body-shock, or, "
                         "with the skins, from the first step at whose end they hold liquid if that comes first, latched "
                         "(default {} with --removal girin, where the droplet population converges; off with instant "
                         "removal). off: the --dt step throughout, as every run before 2026-10-07".format(DEFAULT_DT_CONTINUUM))
    me.add_argument("--kn-body-shock", type=float, default=surface_flow.KN_BODY_SHOCK,
                    help="body Knudsen number below which a distinct bow shock is assumed and the shock-layer construction is used "
                         "(default %(default)s); above it the flow is treated as merged and every melt closure is flagged")
    me.add_argument("--we-critical", type=float, default=dispersion.WE_CRITICAL_THEORY,
                    help="critical surface Weber number of Girin's thick branch (default %(default)s); Girin & Kopyt's thin, "
                         "rarefied and regime-2 branches are gated by the wave fit and the viscous film instead")
    me.add_argument("--kr", type=float, default=spray.K_R, help="droplet radius / wavelength (default %(default)s)")
    me.add_argument("--kt", type=float, default=spray.K_T, help="release period / growth time (default %(default)s)")
    me.add_argument("--prism-layers", type=int, default=None, help="prism layers under the surface (default 4 with --melt on, 0 otherwise)")
    me.add_argument("--layer-thickness", type=float, default=mesh.DEFAULT_LAYER_THICKNESS * 1e3, help="outermost layer thickness [mm] (default %(default)s, growth 2)")
    me.add_argument("--demise-fraction", type=float, default=0.01, help="the run ends when the body mass falls below this fraction of the initial (default %(default)s)")
    me.add_argument("--size-feedback", choices=body.SIZE_FEEDBACK_NAMES, default=None,
                    help="current: Kn on the equivalent diameter of the remaining mass and the stagnation radius fitted to the windward cap "
                         "(default with --removal girin); initial: D0 and R0 throughout, SESAM's convention (default with --removal instant)")
    me.add_argument("--particles", dest="particles", action="store_true", default=True,
                    help="write the particle source table and the size distribution of a spraying run (default)")
    me.add_argument("--no-particles", dest="particles", action="store_false")

    c = sub.add_parser("compare", help="metrics and plots for an existing model history")
    c.add_argument("--model", required=True, help="model history CSV")
    c.add_argument("--reference", required=True, help="SESAM run CSV")
    c.add_argument("--outdir", default=DEFAULT_OUTDIR)
    c.add_argument("--title", default=None)
    c.add_argument("--quiet", action="store_true")
    return p


def build_thermal(args, settings, mass):
    """(ThermalBody or MeltingBody, HeatingModel, settings-provenance dict) for --thermal fem [--melt on]."""
    radius = settings.diameter / 2.0
    melting = args.melt == "on"
    skins = args.surface_model == "skin"
    h_surface, h_core = args.h_surface * 1e-3 * args.mesh_size, args.h_core * 1e-3 * args.mesh_size
    layers = args.prism_layers if args.prism_layers is not None else (mesh.DEFAULT_LAYERS if melting and not skins else 0)
    band = mesh.DEFAULT_BAND if melting or skins else 0.0    # non-melting runs keep the Step 2 field: the committed Step 2
                                                    # verification numbers were measured on it (fact 12); the skins'
                                                    # moving mesh is rebuilt with the dense band's size field
    the_mesh = mesh.sphere_mesh(radius, h_surface, h_core, layers=layers, layer_thickness=args.layer_thickness * 1e-3,
                               band=band)
    if skins and args.fine_band > 0.0:                  # the fine band (FINE_BAND_MM), made again at every rebuild
        the_mesh = movemesh.refined_mesh(the_mesh, args.fine_band * 1e-3, min(args.fine_size * 1e-3, h_surface))
    # the melting default is Scheil's curve (decided 2026-09-27, made the default on 2026-10-07); the linear range stays
    # selectable by name
    material_name = args.material or ("AA7075_scheil" if melting else "AA7075_nomelt")
    mat = material.Material.from_drama_json(material_name)
    if args.k_scale != 1.0:
        mat.k_table = mat.k_table * args.k_scale                  # verification device (near-isothermal body), not physical
    solver = thermal.thermal_solver(args.thermal_solver, linear_solver=args.linear_solver, lumped_mass=not args.consistent_mass)
    if skins:
        # with the flow on the liquid stays, flows and is sprayed from the top cells by Girin's model (sub-plan 22), the
        # exposed slurry too with --slurry-spray on; instant removal is instead of both
        flow = args.skin_flow == "on"
        flow_model = surface_flow.SurfaceFlow(rarefied_shear=args.rarefied_shear, gamma_pm=args.gamma_pm,
                                              kn_body_shock=args.kn_body_shock) if flow else None
        the_body = skinbody.SkinBody(the_mesh, mat, solver, mass, skin_settings(args), T0=args.temperature,
                                     emissivity=args.emissivity, T_ambient=args.t_ambient, deform_depth=args.deform_depth * 1e-3,
                                     remove=melting and not flow, demise_fraction=args.demise_fraction, flow=flow,
                                     flow_cycles=args.skin_flow_cycles, surface_flow=flow_model,
                                     spray=spray_model_of(args, mat) if flow else None,
                                     slurry_spray=args.slurry_spray == "on", particles=args.particles)
    elif melting:
        flow = surface_flow.SurfaceFlow(rarefied_shear=args.rarefied_shear, gamma_pm=args.gamma_pm, kn_body_shock=args.kn_body_shock)
        spray_model = spray_model_of(args, mat)
        size_feedback = args.size_feedback or ("current" if args.removal == "girin" else "initial")
        melt_settings = body.MeltSettings(removal=args.removal, runoff=args.runoff == "on", demise_fraction=args.demise_fraction,
                                          particles=args.particles, size_feedback=size_feedback,
                                          deep_runoff=args.deep_runoff == "on", molten_cascade=args.molten_cascade == "on",
                                          rigid_substrate=args.rigid_substrate == "on")
        the_body = body.MeltingBody(the_mesh, mat, solver, mass, flow, spray_model, melt_settings, T0=args.temperature,
                                    emissivity=args.emissivity, T_ambient=args.t_ambient)
    else:
        the_body = body.ThermalBody(the_mesh, mat, solver, mass, T0=args.temperature, emissivity=args.emissivity, T_ambient=args.t_ambient)
    if args.heating == "sesam":
        heating_model = heating.SesamEquivalentHeating()
    else:
        heating_model = heating.PhysicsHeating(stagnation=args.stagnation, bridging=args.bridging_heat, matting_n=args.matting_n,
                                               accommodation=args.accommodation, catalycity=args.catalycity)
    info = {"heating": args.heating, "material": mat.name, "material_file": os.path.abspath(material.MATERIAL_NAMES.get(material_name, material_name)),
            "emissivity": the_body.emissivity, "t_ambient_K": args.t_ambient, "mesh_file": the_mesh.params["path"],
            "h_surface_mm": h_surface * 1e3, "h_core_mm": h_core * 1e3, "prism_layers": layers, "layer_thickness_mm": args.layer_thickness,
            "n_nodes": the_mesh.n_nodes, "n_elements": the_mesh.n_elements, "n_patches": the_body.surface.n_patches,
            "thermal_solver": args.thermal_solver, "linear_solver": args.linear_solver, "lumped_mass": not args.consistent_mass,
            "k_scale": args.k_scale, "melt": args.melt, "surface_model": args.surface_model}
    if skins:
        info["skin"] = {"min_thickness_mm": args.skin_min_thickness, "fine_zone_mm": args.skin_fine_zone,
                        "cell_um": skin_settings(args).cell * 1e6,
                        "growth": args.skin_growth, "substep_ms": args.skin_substep, "preset": args.skin_preset,
                        "deform_depth_mm": args.deform_depth, "n_columns": the_body.columns.n,
                        "fine_band_mm": args.fine_band, "fine_size_mm": min(args.fine_size, h_surface * 1e3),
                        "swell": args.skin_swell == "on", "flow": args.skin_flow,
                        "flow_cycles": "adaptive" if the_body.flow_cycles is None else the_body.flow_cycles,
                        "slurry_spray": args.slurry_spray}
    if args.heating == "physics":
        info.update({"stagnation": args.stagnation, "bridging_heat": args.bridging_heat, "matting_n": args.matting_n,
                     "accommodation": args.accommodation, "catalycity": args.catalycity})
    if melting and skins:
        info.update({"removal": args.removal, "T_solidus_K": mat.T_solidus, "T_liquidus_K": mat.T_liquidus,
                     "latent_heat_Jkg": mat.latent_heat})
        if the_body.spray is not None:                  # Girin's spraying from the top cells (sub-plan 22)
            info.update({"rt_spray": args.rt_spray, "rarefied_shear": args.rarefied_shear, "we_critical": args.we_critical,
                         "k_r": args.kr, "k_t": args.kt, "gamma_pm": args.gamma_pm, "kn_body_shock": args.kn_body_shock,
                         "particles": args.particles,
                         "liquid": {"rho": mat.liquid.rho, "mu": mat.liquid.mu, "sigma": mat.liquid.sigma}})
    elif melting:
        info.update({"removal": args.removal, "runoff": args.runoff, "deep_runoff": args.deep_runoff,
                     "molten_cascade": args.molten_cascade, "rigid_substrate": args.rigid_substrate, "rt_spray": args.rt_spray,
                     "rarefied_shear": args.rarefied_shear, "we_critical": args.we_critical,
                     "k_r": args.kr, "k_t": args.kt, "demise_fraction": args.demise_fraction, "particles": args.particles,
                     "size_feedback": size_feedback, "gamma_pm": args.gamma_pm, "kn_body_shock": args.kn_body_shock,
                     "liquid": {"rho": mat.liquid.rho, "mu": mat.liquid.mu, "sigma": mat.liquid.sigma},
                     "T_solidus_K": mat.T_solidus, "T_liquidus_K": mat.T_liquidus, "latent_heat_Jkg": mat.latent_heat})
    return the_body, heating_model, info


def spray_model_of(args, mat):
    """The shared spray physics with the command line's constants, for either surface model."""
    return spray.SprayModel(mat.liquid, k_r=args.kr, k_t=args.kt, we_critical=args.we_critical, rt_spray=args.rt_spray == "on")


def skin_settings(args):
    """The layered skins' settings from the command line (sub-plan 20)."""
    cell = args.skin_cell * 1e-6 if args.skin_cell is not None else skin.SKIN_PRESETS[args.skin_preset]
    return skin.SkinSettings(min_thickness=args.skin_min_thickness * 1e-3, cell=cell, growth=args.skin_growth,
                             substep=args.skin_substep * 1e-3, swell=args.skin_swell == "on",
                             fine_zone=args.skin_fine_zone * 1e-3)


def run_name_of(args, diameter_m, velocity_ms, altitude_m, atmosphere_name, dt_continuum):
    """The run's name (`model_run_name`) from its command line, its surface model and removal resolved (`surface_model_of`,
    `removal_of`) and its continuum step (`continuum_step`); --name overrides it."""
    return args.name or model_run_name(diameter_m, velocity_ms, altitude_m, atmosphere_name, args.bridging, args.wind,
                                       args.heating if args.thermal == "fem" else None,
                                       removal_of(args) if args.melt == "on" else None,
                                       deep_runoff=args.deep_runoff == "on", molten_cascade=args.molten_cascade == "on",
                                       seed=args.seed, rigid_substrate=args.rigid_substrate == "on",
                                       dt_continuum=dt_continuum, surface_model=surface_model_of(args),
                                       skin_suffix=skin_name_suffix(args.skin_min_thickness, args.skin_cell, args.skin_growth,
                                                                    args.skin_substep, args.skin_preset, args.deform_depth,
                                                                    args.skin_swell == "on", args.fine_band, args.fine_size,
                                                                    args.skin_fine_zone,
                                                                    skin_flow_of(args) == "on" if args.melt == "on" else None,
                                                                    args.skin_flow_cycles, args.slurry_spray == "on"))


def cmd_run(args, parser):
    # first, before anything that could draw from it (amendment of 2026-10-05; plan facts 52 and 62): the only draws in a
    # run are pyamg's, at every hierarchy the skfem backend builds, and with the generator seeded here two runs of one
    # build agree bit for bit
    np.random.seed(args.seed)
    if args.surface_model == "skin" and args.thermal != "fem":
        parser.error("--surface-model skin needs --thermal fem")
    # the surface model, the flow and the removal not given take the run's own (sub-plan 22 Tasks 3 and 3b): the skins for
    # a thermal run, their flow on when they melt unless instant removal is given, and Girin's removal except for the
    # heat-only skins
    args.surface_model, args.skin_flow, args.removal = surface_model_of(args), skin_flow_of(args), removal_of(args)
    for label, value in (("--diameter", args.diameter), ("--velocity", args.velocity), ("--material-density", args.material_density),
                         ("--cadence", args.cadence), ("--t-max", args.t_max)):
        if value <= 0.0:
            parser.error("{} must be > 0".format(label))
    if args.altitude < 0.0:
        parser.error("--altitude must be >= 0")
    if args.atmosphere not in ATMOSPHERES and not args.atmosphere.startswith("replay:"):
        parser.error("--atmosphere must be one of {} or replay:<sesam.csv>".format(ATMOSPHERES))
    for label, value in (("--dt", args.dt), ("--mesh-size", args.mesh_size), ("--h-surface", args.h_surface), ("--h-core", args.h_core), ("--k-scale", args.k_scale)):
        if value <= 0.0:
            parser.error("{} must be > 0".format(label))
    if not 0.0 <= args.catalycity <= 1.0:
        parser.error("--catalycity must be within [0, 1]")
    if args.thermal == "none" and (args.animate or args.stills or args.frames_every):
        parser.error("--animate/--stills/--frames-every need --thermal fem")
    if args.melt == "on" and args.thermal != "fem":
        parser.error("--melt on needs --thermal fem")
    if not 1.0 < args.gamma_pm < 2.0:
        parser.error("--gamma-pm must be within (1, 2)")
    for label, value in (("--we-critical", args.we_critical), ("--kr", args.kr), ("--kt", args.kt), ("--layer-thickness", args.layer_thickness),
                         ("--kn-body-shock", args.kn_body_shock)):
        if value <= 0.0:
            parser.error("{} must be > 0".format(label))
    if args.prism_layers is not None and args.prism_layers < 0:
        parser.error("--prism-layers must be >= 0")
    if not 0.0 < args.demise_fraction < 1.0:
        parser.error("--demise-fraction must be within (0, 1)")
    if args.dt_continuum is not None and args.melt != "on":
        parser.error("--dt-continuum needs --melt on")
    dt_continuum = continuum_step(args)
    if dt_continuum is not None and dt_continuum > args.dt:
        parser.error("--dt-continuum must be no longer than --dt")
    for label, value in (("--skin-min-thickness", args.skin_min_thickness), ("--skin-fine-zone", args.skin_fine_zone),
                         ("--skin-substep", args.skin_substep),
                         ("--skin-cell", 1.0 if args.skin_cell is None else args.skin_cell),
                         ("--deform-depth", args.deform_depth)):
        if value <= 0.0:
            parser.error("{} must be > 0".format(label))
    if args.skin_growth < 1.0:
        parser.error("--skin-growth must be at least 1")
    if args.fine_band < 0.0 or args.fine_size <= 0.0:
        parser.error("--fine-band must be >= 0 and --fine-size > 0")
    cell_um = args.skin_cell if args.skin_cell is not None else skin.SKIN_PRESETS[args.skin_preset] * 1e6
    if 2.0 * cell_um * 1e-3 > min(args.skin_min_thickness, args.skin_fine_zone) * (1.0 + 1e-12):
        parser.error("--skin-cell must be at most half of --skin-min-thickness and of --skin-fine-zone")
    if args.skin_flow == "on" and not (args.surface_model == "skin" and args.melt == "on"):
        parser.error("--skin-flow on needs --surface-model skin and --melt on: it moves the skins' molten metal")
    if args.skin_flow_cycles is not None and args.skin_flow != "on":
        parser.error("--skin-flow-cycles needs --skin-flow on")
    if args.skin_flow_cycles is not None and args.skin_flow_cycles < 1:
        parser.error("--skin-flow-cycles must be at least 1")
    if args.surface_model == "skin" and args.melt == "on" and args.skin_flow == "on" and args.removal != "girin":
        parser.error("--skin-flow on keeps the liquid, which only spraying removes: it needs --removal girin, not instant "
                     "removal, which stands in for flow and spraying in the heat-only skins (spec 11.2)")
    if args.surface_model == "skin" and args.melt == "on" and args.skin_flow == "off" and args.removal != "instant":
        parser.error("Girin's spraying from the skins needs their flow: the heat-only skins (--skin-flow off) remove their "
                     "liquid as it forms, --removal instant")
    if args.slurry_spray == "on" and args.skin_flow != "on":
        parser.error("--slurry-spray on needs --skin-flow on (with --surface-model skin and --melt on): the exposed slurry "
                     "sprays by the skins' spraying from their top cells")
    if args.surface_model == "skin":
        for flag, given in (("--runoff off", args.runoff == "off"), ("--deep-runoff off", args.deep_runoff == "off"),
                            ("--molten-cascade off", args.molten_cascade == "off"),
                            ("--rigid-substrate off", args.rigid_substrate == "off"),
                            ("--size-feedback", args.size_feedback is not None)):
            if given:
                parser.error("{} is the element model's setting: give --surface-model elements (the layered skins, the "
                             "default since sub-plan 22, have no film, runoff, deep account or element deaths, and their "
                             "regime test always reads their own depths)".format(flag))
    if args.surface_model == "skin" and (args.animate or args.stills):
        parser.error("--animate/--stills read the element model's frames; with --surface-model skin use --frames-every")
    if args.surface_model == "skin" and args.prism_layers is not None and args.prism_layers > 0:
        parser.error("--prism-layers {} is refused with --surface-model skin: the moving mesh's rebuild cannot reproduce prism "
                     "layers (a skin run uses 0, the default)".format(args.prism_layers))

    initial = tj.InitialState(velocity=args.velocity * 1e3, altitude=args.altitude * 1e3,
                              flight_path=math.radians(args.flight_path_angle), heading=math.radians(args.heading),
                              lat=math.radians(args.lat), lon=math.radians(args.lon), epoch=args.epoch)
    settings = tj.Settings(diameter=args.diameter * 1e-3, gravity=args.gravity, rtol=args.rtol,
                           cadence=args.cadence, t_max=args.t_max)
    mass = body.sphere_mass(settings.diameter, args.material_density)
    atm, atm_name, atm_info = make_atmosphere(args.atmosphere, args.epoch, args.wind)
    reference = sesam_io.load_reference(args.reference) if args.reference else None
    name = run_name_of(args, settings.diameter, initial.velocity, initial.altitude, atm_name, dt_continuum)
    run_dir = os.path.join(args.outdir, name)
    os.makedirs(args.outdir, exist_ok=True)
    thermal_info, the_body = {}, body.ConstantBody(mass, args.temperature)
    if args.thermal == "fem":
        the_body, heating_model, thermal_info = build_thermal(args, settings, mass)
    sim = tj.Simulator(initial, the_body, atm, aero.SphereDragTables.from_json(), aero.bridging_by_name(args.bridging), settings)
    if args.thermal == "fem":
        frames_every = args.frames_every or (10 if (args.animate or args.stills) else 0)
        run = coupled.CoupledRun(sim, the_body, heating_model,
                                 coupled.CoupledSettings(dt=args.dt, frames_every=frames_every, output_dir=os.path.join(run_dir, "vtk"),
                                                         dt_continuum=dt_continuum, kn_switch=args.kn_body_shock,
                                                         switch_on_melt=args.surface_model == "skin"))
        history = run.run()
        thermal_info["macro_step_s"], thermal_info["frames_every"] = args.dt, frames_every
        thermal_info["dt_continuum_s"] = dt_continuum
    else:
        history = sim.run(extra_times=reference.time if reference is not None else None)
    csv_path = os.path.join(args.outdir, name + ".csv")
    json_path = os.path.join(args.outdir, name + ".json")
    tj.write_history_csv(history, csv_path)
    prov = provenance()
    # spec section 8: the reference file(s)' SHA-256 also live under provenance, alongside the copies
    # already recorded under "comparison" (--reference) and "settings" (replay_reference/replay_sha256).
    prov["reference_sha256"] = reference.sha256 if reference is not None else None
    prov["replay_sha256"] = atm_info.get("replay_sha256")
    doc = {
        "schema_version": 1,
        "run_name": name,
        "inputs": {"diameter_mm": args.diameter, "initial_velocity_kms": args.velocity, "initial_altitude_km": args.altitude,
                   "flight_path_angle_deg": args.flight_path_angle, "heading_deg": args.heading, "latitude_deg": args.lat,
                   "longitude_deg": args.lon, "epoch_utc": args.epoch.strftime("%Y-%m-%dT%H:%M:%S"),
                   "material_density_kgm3": args.material_density, "mass_kg": mass, "initial_temperature_K": args.temperature},
        "settings": {"atmosphere": atm_name, "wind": args.wind, "bridging": args.bridging, "gravity": args.gravity,
                     "rtol": args.rtol, "atol_position_m": settings.atol_position, "atol_velocity_ms": settings.atol_velocity,
                     "cadence_s": args.cadence, "t_max_s": args.t_max, "seed": args.seed, **atm_info, "thermal": args.thermal,
                     **thermal_info},
        "results": history.results,
        "comparison": None,
        "provenance": prov,
        "files": {"csv": os.path.abspath(csv_path)},
    }
    melting = args.thermal == "fem" and args.melt == "on" and args.surface_model != "skin"
    # the skins' spraying flight (sub-plan 22 Task 3): the element model's droplet files, from the same writer; the
    # per-step summary waits for Task 4, which gives the skin history its droplet columns, and is recorded as not written
    spraying_skins = args.thermal == "fem" and args.melt == "on" and args.surface_model == "skin" and args.skin_flow == "on"
    if melting:
        doc["files"].update({k: os.path.abspath(v) for k, v in coupled.write_particles(run_dir, the_body, history).items()} if args.particles
                            else {})
        doc["files"]["melt_plots"] = [os.path.abspath(p) for p in compare.plot_melt(history, run_dir, name, reference if compare.has_melt(history, reference) else None)]
    elif spraying_skins and args.particles:
        doc["files"].update({k: os.path.abspath(v) for k, v in coupled.write_particles(run_dir, the_body, history).items()})
        doc["files"]["particles_summary"] = None
    if reference is not None:
        plots = compare.plot_all(history, reference, run_dir, name)
        doc["comparison"] = {"reference": reference.name, "reference_csv": reference.csv_path,
                             "reference_sha256": reference.sha256, "metrics": compare.metrics(history, reference),
                             "plots": [os.path.abspath(p) for p in plots]}
        if compare.has_thermal(history, reference):
            doc["comparison"]["thermal_metrics"] = compare.thermal_metrics(history, reference)
            doc["comparison"]["plots"] += [os.path.abspath(p) for p in compare.plot_thermal(history, reference, run_dir, name)]
        if compare.has_melt(history, reference):
            doc["comparison"]["melt_metrics"] = compare.melt_metrics(history, reference)
    if args.thermal == "fem" and (args.animate or args.stills):
        vtk_dir = os.path.join(run_dir, "vtk")
        out = viz.animate(vtk_dir, history, settings.diameter / 2.0, animation=args.animate, melting=melting)
        section = viz.animate_section(vtk_dir, history, settings.diameter / 2.0, animation=args.animate, melting=melting)
        doc["files"]["animation"], doc["files"]["section"] = out.get("animation"), section.get("animation")
        doc["files"]["stills"] = out["stills"] + section["stills"]
        if melting:
            film = viz.animate_film(vtk_dir, history, settings.diameter / 2.0, animation=args.animate)
            doc["files"]["film"] = film.get("animation")
            doc["files"]["stills"] += film["stills"]
    if args.thermal == "fem":
        doc["files"]["vtk_dir"] = os.path.abspath(os.path.join(run_dir, "vtk")) if history.results.get("n_frames") else None
    tj.write_run_json(json_path, doc)
    if not args.quiet:
        res = history.results
        print("{} (seed {}): {} at t = {:.1f} s, final V {:.4f} km/s, Kn {:.3g} -> {:.3g}, {} RHS evaluations in {:.1f} s".format(
            name, args.seed, res["end_reason"], res["final_time_s"], res["final_velocity_kms"], res["knudsen_start"],
            history.columns["knudsen"][-1], res["rhs_evaluations"], res["runtime_s"]))
        if args.thermal == "fem":
            print("  thermal: peak surface T {:.0f} K at t = {:.0f} s, peak mean T {:.0f} K, integrated heat {:.3g} J, "
                  "energy balance residual {:.1e}, {} macro steps, {:.1f} Newton iterations/step".format(
                      res["peak_surface_T_K"], res["time_of_peak_surface_T_s"], res["peak_mean_T_K"], res["integrated_heat_J"],
                      res["energy_balance_residual"], res["n_macro_steps"], res["mean_newton_iterations"]))
        if args.thermal == "fem" and args.surface_model == "skin":
            print("  skins: {} columns, removed {:.4g} kg, drawn {:.4g} kg, {} rebuilds, {} interface repeats, {} base faults; "
                  "wall time {:.0f} s in the skins, {:.0f} s in the 3D solve".format(
                      res["skin_columns"], res["removed_mass_kg"], res["drawn_mass_kg"], res["rebuilds"],
                      res["interface_repeats"], res["jump_faults"], res["wall_skins_s"], res["wall_3d_s"]))
        if spraying_skins:
            print("  spraying: released {:.4g} kg as {:.3g} droplets in {} source rows{}; no particles_summary.csv until "
                  "sub-plan 22 Task 4 gives the skin history its droplet columns".format(
                      the_body.released_mass, the_body.droplet_number, sum(len(b) for b in the_body.source_rows),
                      " (particles.npz, size_distribution.csv)" if args.particles else ""))
        if melting:
            print("  melt: onset {} km, spraying onset {} km, demise {} km at t = {} s; sprayed {:.4f} kg of {:.4f}, film left {:.4f} kg, "
                  "{:.3g} droplets (median r {} um), {} source rows, melt balance {:.1e}".format(
                      _fmt(res["melt_onset_altitude_km"]), _fmt(res["spraying_onset_altitude_km"]), _fmt(res["demise_altitude_km"]),
                      _fmt(res["demise_time_s"]), res["sprayed_mass_kg"], res["initial_mass_kg"], res["film_mass_kg"], res["n_released"],
                      _fmt(res["r_median_um"]), res["n_source_rows"], res["melt_energy_balance_residual"]))
            if reference is not None and "melt_metrics" in doc["comparison"]:
                mm = doc["comparison"]["melt_metrics"]
                print("  mass vs SESAM: max |dm| {:.2%} of m0, melt onset {:+.2f} km, 1 %-mass time {:+.1f} s ({:+.1%})".format(
                    mm["mass"]["max_rel_m0"], mm["onset_altitude_diff_km"], mm["demise_time_diff_s"], mm["demise_time_rel"]))
        if reference is not None:
            hyp = doc["comparison"]["metrics"]["hypersonic"]
            print("  vs {}: hypersonic max |dV| {:.1f} m/s ({:.3%}), max |dh| {:.0f} m; end time {:+.1f} s".format(
                reference.name, hyp["dV_max_ms"], hyp["dV_rel_max"], hyp["dh_max_m"], doc["comparison"]["metrics"]["d_end_time_s"]))
            if "thermal_metrics" in doc["comparison"]:
                tm = doc["comparison"]["thermal_metrics"]
                print("  heat vs SESAM: Q_conv max {:.2%} of peak ({:.2%} point-wise, continuum), integrated heat {:+.2%} (hypersonic) "
                      "{:+.2%} (end), |dT_eq| max {:.1f} K ({:.2%}), radiated max {:.2%} of peak".format(
                          tm["Q_conv"]["max"], tm["Q_conv"]["continuum_rel_max"], tm["integrated_heat"]["rel_error_end_of_hypersonic"],
                          tm["integrated_heat"]["rel_error_end"], tm["temperature"]["dT_max_K"], tm["temperature"]["dT_rel_max"], tm["radiated"]["max"]))
        print("  csv  -> {}\n  json -> {}".format(os.path.abspath(csv_path), os.path.abspath(json_path)))
    # spec section 7: exit 1 when the flight escaped rather than reaching the ground or t_max;
    # the CSV/JSON are already written above so the escaped trajectory is still available.
    if history.end_reason == "escape":
        return 1
    return 0


def _fmt(x):
    return "n/a" if x is None else "{:.1f}".format(x)


def cmd_compare(args, parser):
    history = tj.read_history_csv(args.model)
    reference = sesam_io.load_reference(args.reference)
    stem = os.path.splitext(os.path.basename(args.model))[0]
    title = args.title or "{} vs {}".format(stem, reference.name)
    os.makedirs(args.outdir, exist_ok=True)
    plots = compare.plot_all(history, reference, args.outdir, title)
    doc = {"model": os.path.abspath(args.model), "reference": reference.name, "reference_sha256": reference.sha256,
           "metrics": compare.metrics(history, reference), "plots": [os.path.abspath(p) for p in plots]}
    if compare.has_thermal(history, reference):
        doc["thermal_metrics"] = compare.thermal_metrics(history, reference)
        doc["plots"] += [os.path.abspath(p) for p in compare.plot_thermal(history, reference, args.outdir, title)]
    if compare.has_melt(history, reference):
        doc["melt_metrics"] = compare.melt_metrics(history, reference)
        doc["plots"] += [os.path.abspath(p) for p in compare.plot_melt(history, args.outdir, title, reference)]
    out = os.path.join(args.outdir, stem + "_vs_reference.json")
    tj.write_run_json(out, doc)
    if not args.quiet:
        hyp = doc["metrics"]["hypersonic"]
        print("{}: hypersonic max |dV| {:.1f} m/s ({:.3%}), max |dh| {:.0f} m -> {}".format(title, hyp["dV_max_ms"], hyp["dV_rel_max"], hyp["dh_max_m"], out))
    return 0


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        if args.command == "run":
            return cmd_run(args, parser)
        return cmd_compare(args, parser)
    except KeyError as exc:
        # e.g. fap.solar_indices(): the run epoch (or the day before it) has no record in
        # data/fap_day.dat; could also be a reference CSV missing an expected column.
        key = exc.args[0] if exc.args else exc
        print("ERROR: missing key {!r} (epoch outside the fap file, or a reference CSV without that column)".format(key),
              file=sys.stderr)
        return 2
    except (ValueError, FileNotFoundError) as exc:
        print("ERROR: {}".format(exc), file=sys.stderr)
        return 2
    except MissingBackend as exc:
        print("ERROR: {}".format(exc), file=sys.stderr)
        return 2
    except ModuleNotFoundError as exc:
        print("ERROR: missing optional library {!r}: install requirements-step2.txt into drama_env (dolfinx: the separate "
              "fenicsx_env environment)".format(exc.name), file=sys.stderr)
        return 2
    except RuntimeError as exc:
        print("ERROR: {}".format(exc), file=sys.stderr)
        return 1

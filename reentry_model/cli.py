"""Command line of the re-entry model.

    python -m reentry_model run --diameter 100 --velocity 7.5 --altitude 77.500133 --flight-path-angle -0.959331 \
        [--atmosphere nrlmsise|us76|replay:<sesam.csv>] [--reference <sesam.csv>] [--outdir ...]
        [--thermal fem --heating sesam|physics ... --animate]          (Step 2: coupled 3D conduction)
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

from . import __version__, aero, atmosphere, body, compare, coupled, fap, heating, material, mesh, sesam_io, thermal, viz
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


def parse_epoch(text):
    try:
        return datetime.strptime(text, "%Y-%m-%dT%H:%M:%S")
    except ValueError:
        raise argparse.ArgumentTypeError("epoch must be YYYY-MM-DDTHH:MM:SS, got {!r}".format(text))


def model_run_name(diameter_m, velocity_ms, altitude_m, atmosphere_name, bridging_name, wind_name, heating_name=None):
    name = "model_d{:06.2f}mm_v{:08.5f}kms_h{:07.3f}km_{}_{}_{}".format(
        diameter_m * 1e3, velocity_ms / 1e3, altitude_m / 1e3, atmosphere_name, bridging_name, wind_name)
    return name + ("_fem-" + heating_name if heating_name else "")


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
    th.add_argument("--material", default=material.DEFAULT_MATERIAL, help="DRAMA material JSON (default: packaged AA7075_nomelt)")
    th.add_argument("--emissivity", type=float, default=None, help="override the material's emissivity")
    th.add_argument("--t-ambient", type=float, default=0.0, help="radiation background [K] (default %(default)s, SESAM's)")
    th.add_argument("--mesh-size", type=float, default=1.0, help="multiplies --h-surface and --h-core (default %(default)s)")
    th.add_argument("--h-surface", type=float, default=mesh.DEFAULT_H_SURFACE * 1e3, help="surface element size [mm] (default %(default)s)")
    th.add_argument("--h-core", type=float, default=mesh.DEFAULT_H_CORE * 1e3, help="core element size [mm] (default %(default)s)")
    th.add_argument("--thermal-solver", choices=thermal.SOLVER_NAMES, default="skfem")
    th.add_argument("--linear-solver", choices=("direct", "amg"), default="amg")
    th.add_argument("--lumped-mass", action="store_true")
    th.add_argument("--dt", type=float, default=0.5, help="macro step [s] (default %(default)s)")
    th.add_argument("--frames-every", type=int, default=0, help="VTK frame every n macro steps (default 0: none; 10 with --animate/--stills)")
    th.add_argument("--animate", action="store_true", help="MP4/GIF of the surface temperature and of the meridional cross-section, plus stills")
    th.add_argument("--stills", action="store_true", help="only the stills (start, peak heating, peak surface T, end; surface and section)")

    c = sub.add_parser("compare", help="metrics and plots for an existing model history")
    c.add_argument("--model", required=True, help="model history CSV")
    c.add_argument("--reference", required=True, help="SESAM run CSV")
    c.add_argument("--outdir", default=DEFAULT_OUTDIR)
    c.add_argument("--title", default=None)
    c.add_argument("--quiet", action="store_true")
    return p


def build_thermal(args, settings, mass):
    """(ThermalBody, HeatingModel, settings-provenance dict) for --thermal fem."""
    radius = settings.diameter / 2.0
    h_surface, h_core = args.h_surface * 1e-3 * args.mesh_size, args.h_core * 1e-3 * args.mesh_size
    # band = 0 is the Step 2 field (graded wall to centre over the whole radius). Every run in this package is
    # non-melting, and the committed Step 2 verification numbers were measured on it: the 15 mm default band
    # gives 177 363 tets against 87 632 on the 100 mm reference (measured 2026-09-27). Step 3's Task 13 adds
    # the --band-thickness flag and turns it on for melting runs.
    the_mesh = mesh.sphere_mesh(radius, h_surface, h_core, band=0.0)
    mat = material.Material.from_drama_json(args.material)
    solver = thermal.thermal_solver(args.thermal_solver, linear_solver=args.linear_solver, lumped_mass=args.lumped_mass)
    the_body = body.ThermalBody(the_mesh, mat, solver, mass, T0=args.temperature, emissivity=args.emissivity, T_ambient=args.t_ambient)
    if args.heating == "sesam":
        heating_model = heating.SesamEquivalentHeating()
    else:
        heating_model = heating.PhysicsHeating(stagnation=args.stagnation, bridging=args.bridging_heat, matting_n=args.matting_n,
                                               accommodation=args.accommodation, catalycity=args.catalycity)
    info = {"heating": args.heating, "material": mat.name, "material_file": os.path.abspath(args.material), "emissivity": the_body.emissivity,
            "t_ambient_K": args.t_ambient, "mesh_file": the_mesh.params["path"], "h_surface_mm": h_surface * 1e3, "h_core_mm": h_core * 1e3,
            "n_nodes": the_mesh.n_nodes, "n_elements": the_mesh.n_elements, "n_patches": the_body.surface.n_patches,
            "thermal_solver": args.thermal_solver, "linear_solver": args.linear_solver, "lumped_mass": args.lumped_mass}
    if args.heating == "physics":
        info.update({"stagnation": args.stagnation, "bridging_heat": args.bridging_heat, "matting_n": args.matting_n,
                     "accommodation": args.accommodation, "catalycity": args.catalycity})
    return the_body, heating_model, info


def cmd_run(args, parser):
    for label, value in (("--diameter", args.diameter), ("--velocity", args.velocity), ("--material-density", args.material_density),
                         ("--cadence", args.cadence), ("--t-max", args.t_max)):
        if value <= 0.0:
            parser.error("{} must be > 0".format(label))
    if args.altitude < 0.0:
        parser.error("--altitude must be >= 0")
    if args.atmosphere not in ATMOSPHERES and not args.atmosphere.startswith("replay:"):
        parser.error("--atmosphere must be one of {} or replay:<sesam.csv>".format(ATMOSPHERES))
    for label, value in (("--dt", args.dt), ("--mesh-size", args.mesh_size), ("--h-surface", args.h_surface), ("--h-core", args.h_core)):
        if value <= 0.0:
            parser.error("{} must be > 0".format(label))
    if not 0.0 <= args.catalycity <= 1.0:
        parser.error("--catalycity must be within [0, 1]")
    if args.thermal == "none" and (args.animate or args.stills or args.frames_every):
        parser.error("--animate/--stills/--frames-every need --thermal fem")

    initial = tj.InitialState(velocity=args.velocity * 1e3, altitude=args.altitude * 1e3,
                              flight_path=math.radians(args.flight_path_angle), heading=math.radians(args.heading),
                              lat=math.radians(args.lat), lon=math.radians(args.lon), epoch=args.epoch)
    settings = tj.Settings(diameter=args.diameter * 1e-3, gravity=args.gravity, rtol=args.rtol,
                           cadence=args.cadence, t_max=args.t_max)
    mass = body.sphere_mass(settings.diameter, args.material_density)
    atm, atm_name, atm_info = make_atmosphere(args.atmosphere, args.epoch, args.wind)
    reference = sesam_io.load_reference(args.reference) if args.reference else None
    name = args.name or model_run_name(settings.diameter, initial.velocity, initial.altitude, atm_name, args.bridging, args.wind,
                                       args.heating if args.thermal == "fem" else None)
    run_dir = os.path.join(args.outdir, name)
    os.makedirs(args.outdir, exist_ok=True)
    thermal_info, the_body = {}, body.ConstantBody(mass, args.temperature)
    if args.thermal == "fem":
        the_body, heating_model, thermal_info = build_thermal(args, settings, mass)
    sim = tj.Simulator(initial, the_body, atm, aero.SphereDragTables.from_json(), aero.bridging_by_name(args.bridging), settings)
    if args.thermal == "fem":
        frames_every = args.frames_every or (10 if (args.animate or args.stills) else 0)
        run = coupled.CoupledRun(sim, the_body, heating_model,
                                 coupled.CoupledSettings(dt=args.dt, frames_every=frames_every, output_dir=os.path.join(run_dir, "vtk")))
        history = run.run()
        thermal_info["macro_step_s"], thermal_info["frames_every"] = args.dt, frames_every
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
                     "cadence_s": args.cadence, "t_max_s": args.t_max, **atm_info, "thermal": args.thermal, **thermal_info},
        "results": history.results,
        "comparison": None,
        "provenance": prov,
        "files": {"csv": os.path.abspath(csv_path)},
    }
    if reference is not None:
        plots = compare.plot_all(history, reference, run_dir, name)
        doc["comparison"] = {"reference": reference.name, "reference_csv": reference.csv_path,
                             "reference_sha256": reference.sha256, "metrics": compare.metrics(history, reference),
                             "plots": [os.path.abspath(p) for p in plots]}
        if compare.has_thermal(history, reference):
            doc["comparison"]["thermal_metrics"] = compare.thermal_metrics(history, reference)
            doc["comparison"]["plots"] += [os.path.abspath(p) for p in compare.plot_thermal(history, reference, run_dir, name)]
    if args.thermal == "fem" and (args.animate or args.stills):
        vtk_dir = os.path.join(run_dir, "vtk")
        out = viz.animate(vtk_dir, history, settings.diameter / 2.0, animation=args.animate)
        section = viz.animate_section(vtk_dir, history, settings.diameter / 2.0, animation=args.animate)
        doc["files"]["animation"], doc["files"]["section"] = out.get("animation"), section.get("animation")
        doc["files"]["stills"] = out["stills"] + section["stills"]
    if args.thermal == "fem":
        doc["files"]["vtk_dir"] = os.path.abspath(os.path.join(run_dir, "vtk")) if history.results.get("n_frames") else None
    tj.write_run_json(json_path, doc)
    if not args.quiet:
        res = history.results
        print("{}: {} at t = {:.1f} s, final V {:.4f} km/s, Kn {:.3g} -> {:.3g}, {} RHS evaluations in {:.1f} s".format(
            name, res["end_reason"], res["final_time_s"], res["final_velocity_kms"], res["knudsen_start"],
            history.columns["knudsen"][-1], res["rhs_evaluations"], res["runtime_s"]))
        if args.thermal == "fem":
            print("  thermal: peak surface T {:.0f} K at t = {:.0f} s, peak mean T {:.0f} K, integrated heat {:.3g} J, "
                  "energy balance residual {:.1e}, {} macro steps, {:.1f} Newton iterations/step".format(
                      res["peak_surface_T_K"], res["time_of_peak_surface_T_s"], res["peak_mean_T_K"], res["integrated_heat_J"],
                      res["energy_balance_residual"], res["n_macro_steps"], res["mean_newton_iterations"]))
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

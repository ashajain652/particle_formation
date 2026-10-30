"""The sensitivity set (spec 4.6): each variation of one parameter at a time against its baseline.

    "$PY" -m reentry_model.parent_follower_separation_interactions.analysis.sensitivity [--diameters mm ...]
        [--altitudes km ...] [--theta deg] [--bounds low high] [--parameters name ...] [--no-rarefied]
        [--workers N] [--out DIR] [--force]

Reference set: spheres of 10, 50 and 200 mm at 80 and 60 km released at `--theta` (default 90 degrees, the far wake),
each at both bounds. The variations (`variations`) are Taylor's K by the factors 0.5/0.7 and 0.9/0.7 (Park's original
0.5/0.7/0.9 range applied to Taylor's K), Taylor's exponent 0.30 and 0.40, the turbulent law Park instead of Taylor, the
Feldman hot-axis temperature 2000 and 5000 K and core-speed fraction 0.75, the parent frontal area 3.1 and 12.4 m2
(the spec's 6.2 m2 is the baseline 6.215 m2), the ballistic-coefficient table scaled by 0.9 and 1.1, the sphere density
2700 and 2900 kg/m3, and `wake_flow` `velocity`, `density` and `both` against the baseline `off` (so each effect's share of the change can be read). The rarefied set (`rarefied_variations`) is 10 mm at
90 km at theta 180 and 120 degrees with `shock_standoff="billig"`, `pressure_gradient` off, `follower_wall_T_K` 1000 K
and `shock_thickness_mfp` 1 and 10; where 10 mm is not a rarefied follower at that theta the largest default diameter
that is takes its place (`rarefied_diameter`; the choice is printed and is in the `case` name).

Files under `--out` (default OUTPUT): `sensitivity/runs/` (the cases, resumable as in the sweep), `sensitivity.csv`
(case, parameter, value, bound, residence_time, change_vs_baseline in seconds, shock_layer_exit, t_shock_layer_exit,
then the relative change, release kind, follower regime, plume regime at release and the reason) and
`sensitivity_ranking.csv` (per parameter, the largest absolute change in residence time over all cases).
"""
import argparse
import concurrent.futures
import dataclasses
import math
import os
import sys
from dataclasses import dataclass
from typing import Callable

from reentry_model.parent_follower_separation_interactions import OUTPUT
from reentry_model.parent_follower_separation_interactions.analysis.sweep import (
    DEFAULT_BOUNDS, DEFAULT_DIAMETERS, get_twin, run_or_load, write_csv)
from reentry_model.parent_follower_separation_interactions.laws import LawSettings
from reentry_model.parent_follower_separation_interactions.scenarios import ScenarioSettings

REFERENCE_DIAMETERS = (10.0, 50.0, 200.0)       # mm
REFERENCE_ALTITUDES = (80.0, 60.0)              # km
RAREFIED_ALTITUDE = 90.0                        # km
RAREFIED_THETAS = (180.0, 120.0)                # degrees
RAREFIED_DIAMETER = 10.0                        # mm, used where it is a rarefied follower
TAYLOR_K = LawSettings().taylor_K               # the baseline coefficient the factors multiply
COLUMNS = ("case", "parameter", "value", "bound", "residence_time", "change_vs_baseline", "shock_layer_exit",
           "t_shock_layer_exit", "relative_change", "set", "diameter_mm", "altitude_km", "theta_deg", "release_kind",
           "follower_regime", "regime_at_release", "reason", "censored")


@dataclass(frozen=True)
class Variation:
    """One parameter set to one value: `apply(settings)` returns the settings with it changed."""
    parameter: str
    value: object
    apply: Callable


def _law(**kw):
    return lambda s: dataclasses.replace(s, law=dataclasses.replace(s.law, **kw))


def _shock(**kw):
    return lambda s: dataclasses.replace(s, shock=dataclasses.replace(s.shock, **kw))


def _top(**kw):
    return lambda s: dataclasses.replace(s, **kw)


def variations():
    """The one-at-a-time variations of the reference set (without the baseline), in the order of spec 4.6."""
    out = [Variation("taylor_K", TAYLOR_K * f, _law(taylor_K=TAYLOR_K * f)) for f in (0.5 / 0.7, 0.9 / 0.7)]
    out += [Variation("taylor_n", n, _law(taylor_n=n)) for n in (0.30, 0.40)]
    out += [Variation("turbulent_law", "park", _law(turbulent_law="park"))]
    out += [Variation("T_axis_hot", T, _law(T_axis_hot=T)) for T in (2000.0, 5000.0)]
    out += [Variation("feldman_f", 0.75, _law(feldman_f=0.75))]
    out += [Variation("frontal_area", a, _top(frontal_area=a)) for a in (3.1, 12.4)]
    out += [Variation("beta_scale", b, _top(beta_scale=b)) for b in (0.9, 1.1)]
    out += [Variation("density", rho, _top(density=rho)) for rho in (2700.0, 2900.0)]
    out += [Variation("wake_flow", m, _top(wake_flow=m)) for m in ("velocity", "density", "both")]
    return out


def rarefied_variations():
    """The variations of the rarefied follower's shock-layer model (without the baseline)."""
    return [Variation("shock_standoff", "billig", _shock(shock_standoff="billig")),
            Variation("pressure_gradient", False, _shock(pressure_gradient=False)),
            Variation("follower_wall_T_K", 1000.0, _shock(follower_wall_T_K=1000.0)),
            Variation("shock_thickness_mfp", 1.0, _shock(shock_thickness_mfp=1.0)),
            Variation("shock_thickness_mfp", 10.0, _shock(shock_thickness_mfp=10.0))]


def rarefied_diameter(h_km, theta_deg, bound="low", out_dir=None, candidates=DEFAULT_DIAMETERS, force=False):
    """(diameter, why) for the rarefied set at this altitude and theta. The diameter [mm] is RAREFIED_DIAMETER if its
    follower is rarefied (transitional or free-molecular) there, else the largest of `candidates` that is, else None.
    `why` says what RAREFIED_DIAMETER was instead (None when it is used): "a continuum follower" or "a case with no plume
    law (<reason>)". Each test runs the baseline case (in `out_dir`, so it is kept); an exception from one propagates."""
    out_dir = out_dir if out_dir is not None else os.path.join(OUTPUT, "sensitivity")
    base = ScenarioSettings(bound=bound)
    why = None
    for d in [RAREFIED_DIAMETER] + sorted(candidates, reverse=True):
        s, _ = run_or_load(out_dir, d, h_km, theta_deg, base, force)
        regime = s.get("follower_regime")
        if regime in ("transitional", "free_molecular"):
            return d, why
        if d == RAREFIED_DIAMETER:
            why = ("a continuum follower" if regime == "continuum"
                   else f"a case with no plume law ({s.get('reason')})")
    return None, why


def _case_name(d_mm, h_km, theta):
    return f"sphere{d_mm:g}mm_{h_km:g}km_th{theta:g}"


def _run_job(job):
    """(job, summary or None, error message): one case in this process."""
    out_dir, d_mm, h_km, theta, settings = job[:5]
    try:
        summary, _ = run_or_load(out_dir, d_mm, h_km, theta, settings, job[8])
        return job, summary, ""
    except Exception as exc:
        return job, None, f"{type(exc).__name__}: {exc}"


def _jobs(out_dir, diameters, altitudes, theta, bounds, variation_list, rarefied, force=False):
    """Job tuples (out_dir, d, h, theta, settings, set, parameter, value, force) for the baseline and each variation of every
    reference case, and of every rarefied case."""
    jobs = []
    def add(d, h, th, vs, name):
        for b in bounds:
            base = ScenarioSettings(bound=b)
            jobs.append((out_dir, d, h, th, base, name, "baseline", "", force))
            for v in vs:
                jobs.append((out_dir, d, h, th, v.apply(base), name, v.parameter, v.value, force))
    for d in diameters:
        for h in altitudes:
            add(d, h, theta, variation_list, "reference")
    for h, th, d in rarefied:
        add(d, h, th, rarefied_variations(), "rarefied")
    return jobs


def _rows(results):
    """Rows of sensitivity.csv in job order, with each change taken against the baseline of the same case and bound."""
    base = {}
    for job, summary, _ in results:
        if job[6] == "baseline" and summary is not None:
            base[(job[1:4], job[4].bound)] = summary["residence_time"]
    rows = []
    for job, summary, error in results:
        _, d, h, th, settings, name, parameter, value, _ = job
        if summary is None:
            continue
        b = base.get(((d, h, th), settings.bound))
        res = summary["residence_time"]
        change = res - b if b is not None and res is not None else math.nan
        rel = change / b if b not in (None, 0.0) and not math.isnan(change) else math.nan
        rows.append({"case": _case_name(d, h, th), "parameter": parameter, "value": value, "bound": settings.bound,
                     "residence_time": res, "change_vs_baseline": change, "relative_change": rel, "set": name,
                     "diameter_mm": d, "altitude_km": h, "theta_deg": th,
                     **{k: summary.get(k) for k in ("shock_layer_exit", "t_shock_layer_exit", "release_kind",
                                                    "follower_regime", "regime_at_release", "reason", "censored")}})
    return rows


def ranking(rows):
    """Per parameter, the largest absolute change of residence time over its rows, largest first, with the case."""
    best = {}
    for r in rows:
        c = r["change_vs_baseline"]
        if r["parameter"] == "baseline" or c is None or (isinstance(c, float) and math.isnan(c)):
            continue
        if r["parameter"] not in best or abs(c) > abs(best[r["parameter"]]["change_vs_baseline"]):
            best[r["parameter"]] = r
    return [{"parameter": p, "max_abs_change_s": abs(r["change_vs_baseline"]), "change_s": r["change_vs_baseline"],
             "case": r["case"], "bound": r["bound"], "value": r["value"]}
            for p, r in sorted(best.items(), key=lambda kv: -abs(kv[1]["change_vs_baseline"]))]


def _parser():
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0], prog="sensitivity")
    p.add_argument("--diameters", nargs="+", type=float, default=list(REFERENCE_DIAMETERS), help="sphere diameters, mm")
    p.add_argument("--altitudes", nargs="+", type=float, default=list(REFERENCE_ALTITUDES), help="altitudes, km")
    p.add_argument("--theta", type=float, default=90.0, help="release angle of the reference set, degrees")
    p.add_argument("--bounds", nargs="+", choices=DEFAULT_BOUNDS, default=list(DEFAULT_BOUNDS))
    p.add_argument("--parameters", nargs="+", default=None, help="only these parameters (the baseline always runs)")
    p.add_argument("--no-rarefied", action="store_true", help="skip the rarefied set")
    p.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 2) - 1))
    p.add_argument("--out", default=OUTPUT, help="output folder")
    p.add_argument("--force", action="store_true", help="rerun every case, even if its JSON is current")
    return p


def main(argv=None):
    """Run the sensitivity set. Returns the number of cases that raised."""
    args = _parser().parse_args(argv)
    out = os.fspath(args.out)
    case_dir = os.path.join(out, "sensitivity")
    os.makedirs(case_dir, exist_ok=True)
    keep = (lambda p: True) if args.parameters is None else (lambda p: p in args.parameters)
    var = [v for v in variations() if keep(v.parameter)]
    rare = []
    if not args.no_rarefied:
        for th in RAREFIED_THETAS:
            try:
                d, why = rarefied_diameter(RAREFIED_ALTITUDE, th, args.bounds[0], case_dir, force=args.force)
            except Exception as exc:                       # the rarefied set is skipped, the reference set still runs
                print(f"rarefied set at {RAREFIED_ALTITUDE:g} km, theta {th:g} skipped: {type(exc).__name__}: {exc}")
                continue
            if d is None:
                print(f"no default diameter is a rarefied follower at {RAREFIED_ALTITUDE:g} km, theta {th:g}")
            else:
                if d != RAREFIED_DIAMETER:
                    print(f"{RAREFIED_DIAMETER:g} mm is {why} at {RAREFIED_ALTITUDE:g} km, theta {th:g}; "
                          f"the rarefied set uses {d:g} mm there")
                rare.append((RAREFIED_ALTITUDE, th, d))
    jobs = _jobs(case_dir, args.diameters, args.altitudes, args.theta, args.bounds, var, rare, args.force)
    if args.parameters is not None:                                  # keep the rarefied variations out unless asked for
        jobs = [j for j in jobs if j[5] == "reference" or j[6] == "baseline" or keep(j[6])]
    get_twin()
    print(f"{len(jobs)} cases, {args.workers} worker(s)")
    if args.workers <= 1 or len(jobs) <= 1:
        results = [_run_job(j) for j in jobs]
    else:
        results = [None] * len(jobs)
        with concurrent.futures.ProcessPoolExecutor(max_workers=min(args.workers, len(jobs))) as pool:
            futures = {pool.submit(_run_job, j): i for i, j in enumerate(jobs)}
            for fut in concurrent.futures.as_completed(futures):
                i = futures[fut]
                try:
                    results[i] = fut.result()
                except Exception as exc:                    # a worker that died takes its case with it
                    results[i] = (jobs[i], None, f"{type(exc).__name__}: {exc}")
    rows = _rows(results)
    write_csv(os.path.join(out, "sensitivity.csv"), rows, COLUMNS)
    rank = ranking(rows)
    write_csv(os.path.join(out, "sensitivity_ranking.csv"), rank,
              ["parameter", "max_abs_change_s", "change_s", "case", "bound", "value"])
    failures = [(j, e) for j, s, e in results if s is None]
    print(f"{len(rows)} rows in sensitivity.csv; parameters by largest change in residence time:")
    for r in rank:
        print(f"  {r['parameter']:>20}  {r['max_abs_change_s']:.4g} s  ({r['case']}, {r['bound']}, value {r['value']})")
    for j, e in failures:
        print(f"FAILED {_case_name(*j[1:4])} {j[6]} {j[7]}: {e}")
    return len(failures)


if __name__ == "__main__":
    sys.exit(1 if main() else 0)

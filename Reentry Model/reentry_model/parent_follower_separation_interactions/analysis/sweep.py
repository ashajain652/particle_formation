"""The sweep (spec 4.5): every (diameter, release altitude, release angle, bound, wake_flow) of a sphere released from the
twin of the Generic Spacecraft, classified (AFFIx zone, the follower's own regime) and, where spec 4.3c-d integrates it,
run by `scenarios.sphere_in_spacecraft_plume`.

    "$PY" -m reentry_model.parent_follower_separation_interactions.analysis.sweep [--diameters mm ...] [--altitudes km ...]
        [--thetas deg ...] [--bounds low high] [--wake-flow off both] [--workers N] [--out DIR] [--force]

Work units. The far-wake integration is shared by every far-wake theta of one (diameter, altitude, bound, wake_flow)
through an in-process memo in `scenarios`, so one worker takes one such combination with all its thetas, one after
another; the units go to a process pool (default: CPU count less one workers; one worker runs in this process).
Each worker builds its own twin (`get_twin`, once per process).

Files under `--out` (default `parent_follower_separation_interactions.OUTPUT`):
  runs/<run name>.json, .csv   one case each (`scenarios.write_case`). A case whose JSON exists is not run again, unless
                               that JSON was made by another commit or from a dirty working tree (a stale case: it is
                               rerun and a line is printed) or `--force` is given.
  summary.csv                  one row per (diameter, altitude, theta, bound, wake_flow), sorted in that order (bound
                               low before high, wake_flow in the order off, velocity, density, both), read back from
                               the JSON files. It holds no wall-clock value or git state, so two sweeps of the same
                               settings give byte-identical files.
  crossover.csv                per altitude, wake_flow and bound, from the far-wake runs: the largest diameter whose
                               follower is still behind the parent 1 s after release (positive lag) and the smallest
                               that has passed it (negative lag).
"""
import argparse
import concurrent.futures
import csv
import functools
import json
import math
import os
import sys
from collections import Counter

from reentry_model.parent_follower_separation_interactions import OUTPUT
from reentry_model.parent_follower_separation_interactions.adapters import PARENT_FRONTAL_AREA, twin_parent
from reentry_model.parent_follower_separation_interactions import scenarios
from reentry_model.parent_follower_separation_interactions.scenarios import (
    RELEASES, NO_THETA_RELEASES, ScenarioSettings, run_name, sphere_in_spacecraft_plume, write_case)
from reentry_model.parent_follower_separation_interactions.wake_gas import WAKE_FLOW_MODES

DEFAULT_DIAMETERS = (1.0, 5.0, 10.0, 25.0, 50.0, 100.0, 200.0, 400.0)      # mm
DEFAULT_ALTITUDES = (100.0, 90.0, 80.0, 70.0, 60.0, 50.0)                    # km
DEFAULT_THETAS = (90.0, 120.0, 135.0, 145.0, 160.0, 180.0)                   # degrees
DEFAULT_BOUNDS = ("low", "high")
DEFAULT_WAKE_FLOWS = ("off", "both")
LAG_TIME_S = 1.0                      # s after release at which the lag behind the parent is read for the crossover

# Columns of summary.csv that are always present (blank where a case does not have them), in this order; any other
# scalar value in the case files follows, sorted by name.
SUMMARY_COLUMNS = (
    "name", "diameter_mm", "altitude_km", "theta_deg", "bound", "wake_flow",
    "affix_zone", "release_kind", "follower_regime", "affix_status", "transitional_follower",
    "residence_time", "censored", "reason", "residence_reason", "total_inside", "t_first_entry", "first_interval",
    "n_intervals", "residence_exceeds_t_fz", "t_release", "follower_t_end", "follower_end_reason",
    "exit_altitude_m", "exit_speed", "exit_lag_m", "exit_distance_m", "exit_r_b", "exit_s", "exit_through",
    "regime_at_exit", "laws_at_exit", "regime_at_release", "laws_at_release",
    "plume_distance_out_of_range", "plume_speed_out_of_range", "plume_feldman_extrapolation",
    "kn_f_free_stream", "kn_f_release", "kn_f_min_in_layer", "kn_f_exit", "kn_f_end", "follower_continuum_in_layer",
    "shock_layer_exit", "t_shock_layer_exit", "flux_mismatch", "flux_mismatch_above_tol",
    "surface_match_phi_deg", "shoulder_pressure_ratio",
    "t_hat", "t_fz", "theta_inn", "theta_bif", "affix_r")
SUMMARY_EXCLUDED = frozenset({"kind", "config", "limitation", "wall_time_s", "field_build_time_s",
                              "git_commit", "git_dirty", "parent_sha256"})     # wall-clock values and git state


@functools.lru_cache(maxsize=None)
def get_twin(frontal_area=PARENT_FRONTAL_AREA, beta_scale=1.0):
    """The twin for these settings, built once per process (it takes a second or two)."""
    return twin_parent(frontal_area=frontal_area, beta_scale=beta_scale)


def csv_value(v):
    """Text of a value in a CSV file: blank for None and NaN, `+`-joined lists, otherwise Python's own text."""
    if v is None or (isinstance(v, float) and math.isnan(v)):
        return ""
    if isinstance(v, (list, tuple)):
        return "+".join(str(x) for x in v)
    return v


def write_csv(path, rows, fieldnames, comment=None):
    """Write dict rows to a CSV file (blank for missing and NaN values), after `comment` lines starting with '#'."""
    with open(path, "w", newline="") as fh:
        for line in (comment or []):
            fh.write("# " + line + "\n")
        w = csv.DictWriter(fh, fieldnames=list(fieldnames), extrasaction="ignore", lineterminator="\n")
        w.writeheader()
        for row in rows:
            w.writerow({k: csv_value(row.get(k)) for k in fieldnames})


def case_files(out_dir, d_mm, h_km, theta_deg, settings):
    """{release kind: JSON path} of the files that exist for this case. The release kind is part of a case's name and is
    only known once the case is classified, so a case may have left files under several kinds (an earlier run whose
    classification differed); every kind that can have a theta in its name is looked for."""
    found = {}
    for release in RELEASES:
        if release in NO_THETA_RELEASES:
            continue
        path = os.path.join(out_dir, "runs", run_name("sphere", d_mm, h_km, settings, release=release,
                                                      theta_deg=theta_deg) + ".json")
        if os.path.exists(path):
            found[release] = path
    return found


def remove_other_kinds(out_dir, keep_path, d_mm, h_km, theta_deg, settings):
    """Delete the JSON and history CSV of this case under every release kind except the one whose JSON is `keep_path`
    (orphans of an earlier classification, which would otherwise be rerun, and summarised, forever)."""
    for path in case_files(out_dir, d_mm, h_km, theta_deg, settings).values():
        if path != keep_path:
            os.remove(path)
            if os.path.exists(path[:-5] + ".csv"):
                os.remove(path[:-5] + ".csv")
            print(f"removed orphaned case {os.path.basename(path)[:-5]} (its release kind changed)", flush=True)


def git_state():
    """(HEAD commit, working tree dirty) as `scenarios` records them in a case's JSON. A sweep reads it once and passes
    it to every worker, so every case of a run is compared with, and stamped with, the same state."""
    return scenarios._git_commit(), scenarios._git_dirty()


def stale_reason(summary, git):
    """Why a case's JSON cannot be trusted for the current code (`git` = (commit, dirty)), or None. It cannot if the
    current tree has uncommitted changes (or git cannot say), or the JSON was made from such a tree or by another commit."""
    commit, dirty = git
    if dirty is not False or commit is None:
        return f"the current working tree is dirty or git cannot say (commit {str(commit)[:10]}, dirty {dirty!r})"
    if summary.get("git_dirty") is not False:
        return f"made from a dirty or unknown working tree (git_dirty = {summary.get('git_dirty')!r})"
    if summary.get("git_commit") != commit:
        return f"made by commit {str(summary.get('git_commit'))[:10]}, not {str(commit)[:10]}"
    return None


def _current_case(path, git, resume=True):
    """(summary, None) if the JSON at `path` is current, else (None, reason). Its name must be its file name. For a
    resume (`resume` True) it must also not be stale (`stale_reason`); for a summary it must have been recorded in
    exactly the state `git` of this run, which is what every case of the run was run or resumed in, so a run from a
    dirty tree still summarises its own cases."""
    with open(path) as fh:
        summary = json.load(fh)
    if summary.get("name") != os.path.basename(path)[:-5]:
        return None, f"its name {summary.get('name')!r} is not its file name"
    if resume:
        why = stale_reason(summary, git)
    else:
        why = None if (summary.get("git_commit"), summary.get("git_dirty")) == tuple(git) else "recorded in another state"
    return (summary, None) if why is None else (None, why)


def run_or_load(out_dir, d_mm, h_km, theta_deg, settings, force=False, git=None):
    """(summary dict, ran): the case's summary from its JSON file if one exists under any release kind and is current,
    else run it (in this process, with this process's twin for the settings), write it, and read it back, so that a
    fresh and a resumed case look the same. A stale JSON (`stale_reason`) is rerun with a printed line; `force` reruns
    every case, as does a case with files under more than one release kind. After a run, the files of the case under other
    release kinds are deleted (`remove_other_kinds`). `git` is the
    (commit, dirty) of the run (default: read now); it decides what is stale and is what the new JSON records."""
    git = git if git is not None else git_state()
    files = case_files(out_dir, d_mm, h_km, theta_deg, settings)
    if not force and len(files) > 1:                  # one case, files under several kinds: which is right is unknowable
        print(f"case has files under several release kinds ({', '.join(files)}); rerun", flush=True)
    elif not force:
        for path in files.values():
            summary, why = _current_case(path, git)
            if summary is not None:
                return summary, False
            print(f"stale case {os.path.basename(path)[:-5]}: {why}; rerun", flush=True)
    twin = get_twin(settings.frontal_area, settings.beta_scale)
    case = sphere_in_spacecraft_plume(d_mm, h_km, theta_deg, settings, twin=twin)
    case.summary["git_commit"], case.summary["git_dirty"] = git
    _, path = write_case(case, out_dir)
    remove_other_kinds(out_dir, path, d_mm, h_km, theta_deg, settings)
    with open(path) as fh:
        return json.load(fh), True


def _run_unit(unit):
    """One work unit: all thetas of one (diameter, altitude, bound, wake_flow), run consecutively in this process so the
    far-wake integration is shared. Returns [(theta, status, message)], status "ran", "resumed" or "failed"."""
    d_mm, h_km, bound, wake_flow, thetas, out_dir, force, git = unit
    settings = ScenarioSettings(bound=bound, wake_flow=wake_flow)
    result = []
    for theta in thetas:
        try:
            _, ran = run_or_load(out_dir, d_mm, h_km, theta, settings, force, git)
            result.append((theta, "ran" if ran else "resumed", ""))
        except Exception as exc:                                   # one bad case must not stop the sweep
            result.append((theta, "failed", f"{type(exc).__name__}: {exc}"))
    return result


# ------------------------------------------------------------------------------------------------- summary files

def _sort_key(row):
    return (row["diameter_mm"], row["altitude_km"], row["theta_deg"], DEFAULT_BOUNDS.index(row["bound"]),
            WAKE_FLOW_MODES.index(row["wake_flow"]))


def summary_rows(out_dir, cases, git=None):
    """Rows of summary.csv for the (diameter, altitude, theta, bound, wake_flow) tuples in `cases`, sorted by those five
    fields in that order. Only a current JSON is summarised (`_current_case`: its name is its file name and it was recorded in
    the state `git` of this run), so an orphan or stale file never reaches the table (a case with files under several release kinds is left out)."""
    git = git if git is not None else git_state()
    rows = []
    for d_mm, h_km, theta, bound, wake_flow in cases:
        files = case_files(out_dir, d_mm, h_km, theta, ScenarioSettings(bound=bound, wake_flow=wake_flow))
        if len(files) != 1:                          # none, or several release kinds (ambiguous): not summarised
            continue
        summary, _ = _current_case(next(iter(files.values())), git, resume=False)
        if summary is not None:
            rows.append({k: v for k, v in summary.items() if k not in SUMMARY_EXCLUDED})
    rows.sort(key=_sort_key)
    return rows


def summary_columns(rows):
    """SUMMARY_COLUMNS, then every other scalar or list value found in the rows, sorted by name."""
    extra = sorted({k for r in rows for k, v in r.items() if k not in SUMMARY_COLUMNS and not isinstance(v, dict)})
    return list(SUMMARY_COLUMNS) + extra


def crossover_bracket(lag_by_diameter):
    """(largest diameter with positive lag, smallest with negative lag, monotone) from {diameter: lag}. A positive lag
    is a follower still behind the parent. A bracket end is None when no diameter qualifies. `monotone` is False when
    the sign of the lag changes more than once as the diameter grows (the bracket is then not a crossover)."""
    ds = sorted(lag_by_diameter)
    pos = [d for d in ds if lag_by_diameter[d] > 0.0]
    neg = [d for d in ds if lag_by_diameter[d] < 0.0]
    signs = [lag_by_diameter[d] > 0.0 for d in ds if lag_by_diameter[d] != 0.0]
    changes = sum(1 for a, b in zip(signs, signs[1:]) if a != b)
    return (max(pos) if pos else None), (min(neg) if neg else None), changes <= 1


def _lag_at(csv_path, t_since_release):
    """Lag behind the parent [m] at the history row nearest `t_since_release`, or None if the history ends before it."""
    with open(csv_path) as fh:
        rows = [(float(r["t_since_release"]), float(r["lag_m"])) for r in csv.DictReader(fh)]
    if not rows or rows[-1][0] < t_since_release:
        return None
    return min(rows, key=lambda r: abs(r[0] - t_since_release))[1]


def crossover_rows(out_dir, rows):
    """One row per (altitude, wake_flow, bound) with at least one far-wake run: the bracket of `crossover_bracket`
    from the lag at LAG_TIME_S after release of each far-wake diameter (all far-wake thetas share one integration, so
    the first theta found is used)."""
    lags = {}
    for r in rows:
        if r.get("release_kind") != "far_wake":
            continue
        key = (r["altitude_km"], r["wake_flow"], r["bound"])
        per_d = lags.setdefault(key, {})
        if r["diameter_mm"] in per_d:
            continue
        csv_path = os.path.join(out_dir, "runs", run_name(
            "sphere", r["diameter_mm"], r["altitude_km"], ScenarioSettings(bound=r["bound"], wake_flow=r["wake_flow"]),
            release="far_wake", theta_deg=r["theta_deg"]) + ".csv")
        lag = _lag_at(csv_path, LAG_TIME_S) if os.path.exists(csv_path) else None
        if lag is not None:
            per_d[r["diameter_mm"]] = lag
    out = []
    for (h, flow, bound), per_d in sorted(lags.items(), key=lambda kv: (-kv[0][0], WAKE_FLOW_MODES.index(kv[0][1]),
                                                                       DEFAULT_BOUNDS.index(kv[0][2]))):
        low, high, monotone = crossover_bracket(per_d)
        out.append({"altitude_km": h, "wake_flow": flow, "bound": bound, "lag_time_s": LAG_TIME_S,
                    "largest_diameter_positive_lag_mm": low, "lag_there_m": per_d.get(low),
                    "smallest_diameter_negative_lag_mm": high, "lag_there_negative_m": per_d.get(high),
                    "n_far_wake_diameters": len(per_d), "monotone": monotone})
    return out


# ------------------------------------------------------------------------------------------------------ command

def _parser():
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0], prog="sweep")
    p.add_argument("--diameters", nargs="+", type=float, default=list(DEFAULT_DIAMETERS), help="sphere diameters, mm")
    p.add_argument("--altitudes", nargs="+", type=float, default=list(DEFAULT_ALTITUDES), help="release altitudes, km")
    p.add_argument("--thetas", nargs="+", type=float, default=list(DEFAULT_THETAS), help="release angles, degrees")
    p.add_argument("--bounds", nargs="+", choices=DEFAULT_BOUNDS, default=list(DEFAULT_BOUNDS))
    p.add_argument("--wake-flow", nargs="+", choices=WAKE_FLOW_MODES, default=list(DEFAULT_WAKE_FLOWS), dest="wake_flow")
    p.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 2) - 1), help="worker processes")
    p.add_argument("--out", default=OUTPUT, help="output folder")
    p.add_argument("--force", action="store_true", help="rerun every case, even if its JSON is current")
    return p


def kept_altitudes(altitudes, twin):
    """The altitudes inside the twin's range (the others are dropped with a printed message), in the order given."""
    kept = []
    for h in altitudes:
        try:
            twin.t_at_altitude(h)
            kept.append(h)
        except ValueError as exc:
            print(f"dropping altitude {h:g} km: {exc}")
    return kept


def regime_counts(rows):
    """{altitude_km: Counter of follower regimes} over the distinct (diameter, theta) of the rows; a case with no plume
    law has the regime "no_plume_law"."""
    seen, counts = set(), {}
    for r in rows:
        k = (r["altitude_km"], r["diameter_mm"], r["theta_deg"])
        if k in seen:
            continue
        seen.add(k)
        counts.setdefault(r["altitude_km"], Counter())[r.get("follower_regime") or "no_plume_law"] += 1
    return counts


def main(argv=None):
    """Run the sweep. Returns the number of cases that raised (0 when all went well)."""
    args = _parser().parse_args(argv)
    out = os.fspath(args.out)
    os.makedirs(os.path.join(out, "runs"), exist_ok=True)
    altitudes = kept_altitudes(args.altitudes, get_twin())
    git = git_state()
    if git[1] is not False:
        print(f"WARNING: the working tree is dirty or git cannot say (dirty = {git[1]!r}): every existing case counts as "
              f"stale and is rerun, and the new cases are recorded as made from an unclean tree", flush=True)
    thetas = sorted(args.thetas)
    units = [(d, h, b, f, tuple(thetas), out, args.force, git)
             for d in sorted(args.diameters) for h in altitudes for b in args.bounds for f in args.wake_flow]
    print(f"{len(units)} work units, {len(units) * len(thetas)} cases, {args.workers} worker(s), output {out}")
    failures, n_ran, n_resumed = [], 0, 0

    def collect(unit, result):
        nonlocal n_ran, n_resumed
        for theta, status, message in result:
            n_ran += status == "ran"
            n_resumed += status == "resumed"
            if status == "failed":
                failures.append((unit[:4] + (theta,), message))

    if args.workers <= 1 or len(units) <= 1:
        for i, unit in enumerate(units, 1):
            collect(unit, _run_unit(unit))
            print(f"[{i}/{len(units)}] d {unit[0]:g} mm, h {unit[1]:g} km, {unit[2]}, {unit[3]} done", flush=True)
    else:
        with concurrent.futures.ProcessPoolExecutor(max_workers=min(args.workers, len(units))) as pool:
            futures = {pool.submit(_run_unit, u): u for u in units}
            for i, fut in enumerate(concurrent.futures.as_completed(futures), 1):
                unit = futures[fut]
                try:
                    collect(unit, fut.result())
                except Exception as exc:                            # a worker that died takes its whole unit with it
                    failures.extend((unit[:4] + (t,), f"{type(exc).__name__}: {exc}") for t in unit[4])
                print(f"[{i}/{len(units)}] d {unit[0]:g} mm, h {unit[1]:g} km, {unit[2]}, {unit[3]} done", flush=True)

    failed = {c for c, _ in failures}                  # a failed case may leave a stale JSON: it is not summarised
    cases = [(d, h, t, b, f) for d in args.diameters for h in altitudes for t in thetas for b in args.bounds
             for f in args.wake_flow if (d, h, b, f, t) not in failed]
    rows = summary_rows(out, cases, git)
    write_csv(os.path.join(out, "summary.csv"), rows, summary_columns(rows))
    write_csv(os.path.join(out, "crossover.csv"), crossover_rows(out, rows),
              ["altitude_km", "wake_flow", "bound", "lag_time_s", "largest_diameter_positive_lag_mm", "lag_there_m",
               "smallest_diameter_negative_lag_mm", "lag_there_negative_m", "n_far_wake_diameters", "monotone"])
    print(f"{n_ran} cases run, {n_resumed} resumed, {len(failures)} failed; {len(rows)} rows in summary.csv")
    for h, c in sorted(regime_counts(rows).items(), reverse=True):
        print(f"  {h:g} km follower regimes (distinct diameter and theta): " +
              ", ".join(f"{k} {c[k]}" for k in ("continuum", "transitional", "free_molecular", "no_plume_law")))
    print("  AFFIx zones: " + ", ".join(f"{k} {v}" for k, v in sorted(Counter(r.get("affix_zone") or "none" for r in rows).items())))
    print("  release kinds: " + ", ".join(f"{k} {v}" for k, v in sorted(Counter(r.get("release_kind") or "none" for r in rows).items())))
    exits = Counter(r["shock_layer_exit"] for r in rows if r.get("shock_layer_exit"))
    print("  shock-layer outcomes of rarefied followers: " + (", ".join(f"{k} {v}" for k, v in sorted(exits.items())) or "none"))
    for case, message in failures:
        print(f"FAILED d {case[0]:g} mm, h {case[1]:g} km, {case[2]}, {case[3]}, theta {case[4]:g}: {message}")
    return len(failures)


if __name__ == "__main__":
    sys.exit(1 if main() else 0)

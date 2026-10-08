"""`python -m spheral_frag analyse`: a Spheral run's outputs -> the fragment record, the debris log and the mass
accounts (spec §4.2, §11; M1 plan, Task 9 -- the skeleton M3 fills).

    "$PY" -m spheral_frag analyse --run spheral_output/runs/<run name> [--prepared <dir>] [--min-particles N]
          [--outdir spheral_output/analyse] [--force] [--quiet]

Reads the run directory through Task 8's runner -> analyse contract (`record`): `meta.json`, every
`checks/check_<n>.npz` and `removed.npz`; and from the prepared inputs the run was built from (`--prepared`,
default `spheral_output/prepare/<meta.prepared>`) the material table and the flight table. Then, check by check:

- classifies the particles (`record.classify_groups`: main body, fragments with attached dust, debris) and books
  the accounts (`record.accounts_at`: main + fragments + debris + film account = the starting mass, summed exactly);
- follows every fragment across checks by its permanent IDs (`record.match_groups`), from the first check at which
  it was a separate group; once it is clear of the main body (`record.is_clear`) it is released, and its row of the
  record is built from the particles of that first separate check (spec §11.1-11.2: release is the first check at
  which the group was separate, decided once it is clear);
- logs each sub-floor group and each isolated dust particle once, at the first check it appears as debris.

Writes `<outdir>/<analysis name>/`: `fragments.csv` (record.FRAGMENT_COLUMNS), `debris.csv` (record.DEBRIS_COLUMNS)
and last `analyse.json` (accounts at every check, the largest residual, counts, fragments followed but never
clear, the provenance of the run and of the prepared inputs). The analysis name is the run's name
(`naming.analyse_name`); `--min-particles` other than the run's bracket renames it with that bracket's suffix, since
it changes the result.

Route and mechanism of a row (spec §11.2 Origin): `meta.json`'s `route` and `mechanism` when the runner writes them,
otherwise from the run's mode (replay and separation: route "separation", mechanism "tearing"; slurry:
"slurry_breakup"; detachment: "neck_failure"; synthetic: "synthetic", "ring_release") -- an M1 default until M3's
runner records each group's own origin.

Exit codes: 0 ok; 1 accounts not closing within ACCOUNTS_RTOL of the starting mass, or unreadable run outputs;
2 bad arguments or a missing run or prepared directory. Resumable: an existing `analyse.json` whose run files match
(sha256 of meta.json and removed.npz, the check list) is skipped; `--force` rebuilds.

Numpy and the standard library at module level."""
from __future__ import annotations

import argparse
import datetime
import json
import math
import os
import sys
import time

import numpy as np

from . import __version__, material, naming, record
from .prepare import REPO_ROOT, git_commit, jsonable, sha256_file, write_json

DEFAULT_OUTDIR = os.path.join(REPO_ROOT, "spheral_output", "analyse")
DEFAULT_PREPARED = os.path.join(REPO_ROOT, "spheral_output", "prepare")
ANALYSE_SCHEMA = "spheral_frag.analyse"
ANALYSE_SCHEMA_VERSION = 1
ACCOUNTS_RTOL = 1e-12            # plan Task 9: accounts must close within 1e-12 m0 (Task 8 measured 2e-18)

ORIGIN_BY_MODE = {"replay": ("separation", "tearing"), "separation": ("separation", "tearing"),
                  "slurry": ("slurry", "slurry_breakup"), "detachment": ("detachment", "neck_failure"),
                  "synthetic": ("synthetic", "ring_release")}


def _say(quiet, *a):
    if not quiet:
        print(*a, flush=True)


def _resolve_prepared(meta, arg):
    if arg:
        return os.path.abspath(arg)
    for cand in (meta.get("prepared_dir"), os.path.join(DEFAULT_PREPARED, meta["prepared"])):
        if cand and os.path.isdir(cand):
            return os.path.abspath(cand)
    return os.path.join(DEFAULT_PREPARED, meta["prepared"])


def _source_signature(run_dir, checks):
    removed = os.path.join(run_dir, "removed.npz")
    return {"meta_sha256": sha256_file(os.path.join(run_dir, "meta.json")),
            "removed_sha256": sha256_file(removed) if os.path.isfile(removed) else None,
            "checks": [n for n, _ in checks]}


def follow(run_dir, checks, removed, meta, min_particles, link_h, quiet=True):
    """The check-by-check pass (module docstring). Returns (accounts per check; released fragments as (track, first
    separate check, label at that check, check at which it was clear); fragments followed but never clear, with the
    check at which they stopped being separate fragments or None if still separate at the last check; debris items as
    (first check, kind, permanent IDs))."""
    m0 = meta["initial_mass_kg"]
    accounts, released, debris, ended = [], [], [], []
    tracks = {}                     # track id -> {"first": n, "label": label at n, "ids": ids at the last check}
    next_track = 0
    logged = set()                  # permanent IDs already in the debris log
    for n, path in checks:
        p = record.read_check(path)
        if p.check != n:
            raise ValueError("check file {} holds check {}".format(path, p.check))
        groups = record.classify_groups(p, min_particles=min_particles, link_h=link_h)
        acc = record.accounts_at(m0, p, groups, removed, n)
        acc["n_particles"] = p.n
        accounts.append(acc)
        # fragments: follow the open tracks by their IDs, open new ones, release the clear ones
        curr = record.group_ids(p, groups["fragments"])
        prev = {tid: tr["ids"] for tid, tr in tracks.items()}
        mapping = record.match_groups(prev, curr)
        taken = {}                                    # current label -> track
        for tid in sorted(mapping):
            taken.setdefault(mapping[tid], tid)       # two tracks merged: the older one keeps the group
        for tid in [t for t in tracks if t not in taken.values()]:
            tr = tracks.pop(tid)                      # no longer a separate fragment (merged back, below the floor)
            if not tr.get("released"):
                ended.append({"first_check": tr["first"], "label_at_first": tr["label"], "n_ids": len(tr["ids"]),
                              "ended_check": n})
        for lab in sorted(curr):
            if lab not in taken:
                tracks[next_track] = {"first": n, "label": lab}
                taken[lab] = next_track
                next_track += 1
        for lab, tid in sorted(taken.items()):
            tr = tracks[tid]
            tr["ids"] = curr[lab]                     # a released fragment is still followed, so it is not re-opened
            if not tr.get("released") and record.is_clear(p, groups["fragments"][lab], groups["main"]):
                released.append((tid, tr["first"], tr["label"], n))
                tr["released"] = True
        # debris: sub-floor groups and isolated dust, each logged once
        for lab in sorted(groups["debris_groups"]):
            sel = groups["debris_groups"][lab]
            ids = p.id[sel]
            if not set(ids.tolist()) <= logged:
                debris.append((n, "small_group", np.sort(ids)))
                logged.update(ids.tolist())
        for i in np.flatnonzero(groups["debris_dust"]):
            pid = int(p.id[i])
            if pid not in logged:
                debris.append((n, "dust", np.array([pid])))
                logged.add(pid)
        _say(quiet, "  check {:5d}  t {:8.3f} s  particles {:7d}  fragments {:3d}  residual {:.3e}".format(
            n, p.t, p.n, len(groups["fragments"]), acc["rel_residual"]))
    unreleased = ended + [{"first_check": tr["first"], "label_at_first": tr["label"], "n_ids": len(tr["ids"]),
                           "ended_check": None} for tr in tracks.values() if not tr.get("released")]
    return accounts, released, unreleased, debris


def analyse(args) -> int:
    t_start = time.perf_counter()
    quiet = getattr(args, "quiet", False)
    run_dir = os.path.abspath(args.run)
    meta_path = os.path.join(run_dir, "meta.json")
    if not os.path.isfile(meta_path):
        print("no Spheral run at {} (meta.json missing)".format(run_dir), file=sys.stderr)
        return 2
    try:
        meta = record.read_meta(meta_path)
    except (ValueError, KeyError) as e:
        print("meta.json: {}".format(e), file=sys.stderr)
        return 1
    parsed = naming.parse_run_name(meta["run_name"])
    brackets = dict(naming.BRACKET_DEFAULTS, **meta["brackets"])
    if args.min_particles is not None:
        if args.min_particles < 1:
            print("--min-particles must be a positive integer", file=sys.stderr)
            return 2
        brackets["min_particles"] = int(args.min_particles)
    name = naming.analyse_name(naming.run_name(parsed["prepared"], parsed["mode"], parsed["frames"], parsed["form"],
                                               parsed["dx_mm"], brackets, parsed["seed"]))
    prepared = _resolve_prepared(meta, args.prepared)
    for need in ("prepare.json", "material_table.npz", "flight.npz"):
        if not os.path.isfile(os.path.join(prepared, need)):
            print("prepared inputs {}: missing {}".format(prepared, need), file=sys.stderr)
            return 2
    checks = record.list_checks(run_dir)
    if not checks:
        print("run {} has no checks/check_<n>.npz".format(run_dir), file=sys.stderr)
        return 1
    out = os.path.join(os.path.abspath(args.outdir), name)
    json_path = os.path.join(out, "analyse.json")
    source = _source_signature(run_dir, checks)
    if os.path.isfile(json_path) and not args.force:
        try:
            with open(json_path, encoding="utf-8") as fh:
                old = json.load(fh)
        except (OSError, ValueError):
            old = None
        if old and old.get("source") == source:
            print("analysed already ({}): {}".format("passed" if old.get("passed") else "FAILED", out))
            return 0 if old.get("passed") else 1

    with open(os.path.join(prepared, "prepare.json"), encoding="utf-8") as fh:
        prep = json.load(fh)
    table = material.MaterialTable.load(os.path.join(prepared, "material_table.npz"))
    flight = record.FlightTable.load(os.path.join(prepared, "flight.npz"))
    removed_path = os.path.join(run_dir, "removed.npz")
    removed = record.read_removed(removed_path) if os.path.isfile(removed_path) else record.Removal.empty()
    route, mechanism = meta.get("route"), meta.get("mechanism")
    if route is None or mechanism is None:
        r, m = ORIGIN_BY_MODE[meta["mode"]]
        route, mechanism = route or r, mechanism or m
    _say(quiet, "analyse {} ({} checks) -> {}".format(meta["run_name"], len(checks), out))
    try:
        accounts, released, unreleased, debris = follow(run_dir, checks, removed, meta, brackets["min_particles"],
                                                        brackets["link_h"], quiet)
    except ValueError as e:
        print("run outputs: {}".format(e), file=sys.stderr)
        return 1

    form, dx = meta["form"], float(meta["dx_mm"]) * 1e-3
    window = "k{:05d}-{:05d}".format(*parsed["frames"])
    provenance = {"run": meta["run_name"], "window": window, "dx_mm": meta["dx_mm"], "seed": meta["seed"],
                  "brackets": brackets}
    paths = dict(checks)
    rows = []
    cache = {}

    def particles(n):
        if n not in cache:
            cache.clear()
            cache[n] = record.read_check(paths[n])
        return cache[n]

    for number, (tid, n_first, label, n_clear) in enumerate(sorted(released, key=lambda r: (r[1], r[0]))):
        p = particles(n_first)
        groups = record.classify_groups(p, min_particles=brackets["min_particles"], link_h=brackets["link_h"])
        sel = groups["fragments"][label]
        row = record.fragment_row(p, sel, groups["fragment_dust"][label], table, flight, route, mechanism,
                                  provenance, form=form, number=number)
        rows.append(row)
    debris_rows = []
    for number, (n, kind, ids) in enumerate(debris):
        p = particles(n)
        sel = np.isin(p.id, ids)
        debris_rows.append(record.debris_row(p, sel, flight, dx, table=table, form=form, number=number, kind=kind))

    os.makedirs(out, exist_ok=True)
    record.write_rows(os.path.join(out, "fragments.csv"), rows, record.FRAGMENT_COLUMNS)
    record.write_rows(os.path.join(out, "debris.csv"), debris_rows, record.DEBRIS_COLUMNS)
    worst = max(abs(a["rel_residual"]) for a in accounts)
    passed = worst <= ACCOUNTS_RTOL
    m0 = meta["initial_mass_kg"]
    doc = {
        "schema": ANALYSE_SCHEMA, "schema_version": ANALYSE_SCHEMA_VERSION, "name": name,
        "created_utc": datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "spheral_frag_version": __version__, "git": git_commit(), "command": " ".join(
            ["python -m spheral_frag"] + list(getattr(args, "argv", None) or ["analyse"])),
        "run": {"dir": run_dir, "meta": meta}, "source": source,
        "prepared": {"dir": prepared, "name": prep.get("name"), "created_utc": prep.get("created_utc"),
                     "prepare_json_sha256": sha256_file(os.path.join(prepared, "prepare.json")),
                     "fe_run": {k: prep.get("fe_run", {}).get(k) for k in ("name", "json_sha256", "csv_sha256",
                                                                            "package")},
                     "material_table_sha256": sha256_file(os.path.join(prepared, "material_table.npz"))},
        "settings": {"min_particles": brackets["min_particles"], "link_h": brackets["link_h"],
                     "dust_damage": record.DUST_DAMAGE, "clear_h": record.CLEAR_H, "route": route,
                     "mechanism": mechanism, "accounts_rtol": ACCOUNTS_RTOL},
        "counts": {"n_checks": len(checks), "n_fragments": len(rows), "n_unreleased": len(unreleased),
                   "n_debris_rows": len(debris_rows), "n_removals": len(removed),
                   "n_removals_by_reason": {name_: int(np.count_nonzero(removed.reason == code))
                                            for code, name_ in record.REASONS.items()},
                   "fragment_mass_kg": math.fsum(r["mass_kg"] for r in rows),
                   "debris_mass_kg": math.fsum(r["mass_kg"] for r in debris_rows)},
        "released": [{"fragment": i, "release_check": r[1], "clear_check": r[3]}
                     for i, r in enumerate(sorted(released, key=lambda r: (r[1], r[0])))],
        "unreleased": unreleased,
        "accounts": accounts,
        "max_abs_rel_residual": worst, "max_abs_residual_kg": worst * m0,
        "passed": passed, "files": {"fragments": "fragments.csv", "debris": "debris.csv"},
        "timing_s": time.perf_counter() - t_start,
    }
    write_json(json_path, doc)
    _say(quiet, "  {} fragments, {} debris rows, largest account residual {:.3e} m0: {}".format(
        len(rows), len(debris_rows), worst, "closed" if passed else "NOT CLOSED"))
    if not passed:
        print("accounts do not close: largest |residual| {:.3e} m0 > {:g}".format(worst, ACCOUNTS_RTOL),
              file=sys.stderr)
        return 1
    return 0


def add_arguments(p: argparse.ArgumentParser) -> None:
    p.add_argument("--run", required=True, help="Spheral run directory (spheral_output/runs/<run name>)")
    p.add_argument("--prepared", default=None, help="prepared inputs (default: spheral_output/prepare/<prepared>)")
    p.add_argument("--min-particles", type=int, default=None, help="resolution floor (default: the run's bracket, "
                   "{})".format(naming.BRACKET_DEFAULTS["min_particles"]))
    p.add_argument("--outdir", default=DEFAULT_OUTDIR, help="parent of the analysis directory (default "
                   "spheral_output/analyse)")
    p.add_argument("--force", action="store_true", help="rebuild even if analysed already")
    p.add_argument("--quiet", action="store_true", help="print only the summary and errors")

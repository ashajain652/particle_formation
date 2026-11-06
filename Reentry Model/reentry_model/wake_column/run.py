"""One case from droplet to outputs: droplet run, trail, Stage A slowing-down, the 12-hour box, and the files.

`run` chains the pieces for one `Case` and writes `<out_dir>/<run_name>/` (spec section 6):

  history.csv            one row per output time: plume radius and volume, ambient temperature, the share of the
                         aluminium inventory in each gas species, in plume particles, background smoke and settled
                         matter, in every (source molecule, route) slice, and the ambient O, O2, O3, H and H2O
                         concentrations;
  size_distribution.csv  the particle bins at 1e-3, 1, 60, 600, 3600 and 43200 s (the nearest output times that do
                         not pass the end of the run): number, diameter, composition and the aluminium atoms of
                         every (source, route) slice. The composition columns of a bin holding under
                         MASK_FRACTION of the inventory in aluminium are left empty (its O/Al and H/Al are set by a
                         handful of atoms and are noise); its diameter is then the alumina diameter;
  summary.json           the case, settings, trail, hot start and outcome; written last, so a run that stopped
                         part-way leaves no summary and a sweep repeats it.

With Settings.n_segments > 0 (the default) the trail is a row of segments (segments.py). The files above then hold the
whole trail at each time since release, history.csv gains frac_unreleased (the share in segments not yet laid), and two
files are added:

  segments.csv           one row per segment: its position, release time, altitude, aluminium share, hot start and its
                         own shares at the end of the run;
  along_trail.csv        at 1, 60, 600, 3600 and 43200 s (those within the run) and PROFILE_POINTS positions: atoms per
                         metre of trail (line_*) and concentration on the axis (c_*) of each aluminium gas species,
                         plume particles and smoke uptake, with spreading along the trail applied after the runs.

"Condensed" (decision 17) is everything that left the gas phase: plume particles + background-smoke uptake + settled
matter. A finished run (its summary exists) is returned as it stands unless `overwrite`.
"""
import dataclasses
import json
import math
import os
import subprocess
import time

import numpy as np

from reentry_model.wake_column import OUTPUT
from reentry_model.wake_column.ambient import AmbientColumn
from reentry_model.wake_column.box import PROVENANCE_ROUTES, ROUTES, SOURCES, Box
from reentry_model.wake_column.particles import diameter
from reentry_model.wake_column.reactions import load_reactions
from reentry_model.wake_column.segments import TrailSolution, along_trail, combined_hot, cut_trail, solve_segments
from reentry_model.wake_column.settings import Case, Settings, run_name
from reentry_model.wake_column.slowing import calibrate, slow_down
from reentry_model.wake_column.species import AL_SPECIES, load_species
from reentry_model.wake_column.trail import rows_for, trail_for, trail_from_rows

MASK_FRACTION = 1e-9                                 # bins holding less aluminium than this share of the inventory
SIZE_TIMES_S = (1e-3, 1.0, 60.0, 600.0, 3600.0, 43200.0)
HISTORY_AMBIENT = ("O", "O2", "O3", "H", "H2O")
ALUMINA_LIKE_BELOW = 0.5                             # H/Al below this counts as alumina-like
STEADY_TOLERANCE = 0.01                              # decision 18: change per decade, as a share of the inventory
SKIPPED_NO_AL = "no aluminium vaporised"
PROFILE_TIMES_S = (1.0, 60.0, 600.0, 3600.0, 43200.0)
PROFILE_POINTS = 201
PROFILE_REACH = 3.0                                  # along_trail.csv extends this many radii past the trail's ends


# ---- steady state (decision 18) -------------------------------------------------------------------------------------

def steady_state_time(t, reservoirs):
    """The earliest output time t_ss such that, for every output time t_j >= t_ss with 10 t_j <= t_end, each
    reservoir (a mapping name -> array of inventory fractions at the times `t`) changes by less than 0.01 between
    t_j and 10 t_j; None when no output time qualifies (including a run shorter than one decade past its first
    point, since a decade of evidence is required). The value at 10 t_j is interpolated linearly in log t."""
    t = np.asarray(t, dtype=float)
    pos = t > 0.0
    tp = t[pos]
    if tp.size == 0:
        return None
    t_end = t[-1]
    cand = np.nonzero(10.0 * tp <= t_end * (1.0 + 1e-9))[0]
    if cand.size == 0:
        return None
    logt = np.log(tp)
    quiet = np.ones(cand.size, dtype=bool)
    for f in reservoirs.values():
        f = np.asarray(f, dtype=float)[pos]
        later = np.interp(logt[cand] + math.log(10.0), logt, f)
        quiet &= np.abs(later - f[cand]) < STEADY_TOLERANCE
    if not quiet[-1]:
        return None
    bad = np.nonzero(~quiet)[0]
    first = 0 if bad.size == 0 else int(bad[-1]) + 1
    return float(tp[cand[first]])


# ---- summary --------------------------------------------------------------------------------------------------------

def _clean(x):
    """JSON-ready copy of x: numpy values as Python ones, NaN and infinity as None, tuples as lists, keys as text."""
    if isinstance(x, dict):
        return {str(k): _clean(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [_clean(v) for v in x]
    if isinstance(x, np.ndarray):
        return _clean(x.tolist())
    if isinstance(x, np.generic):
        return _clean(x.item())
    if isinstance(x, float):
        return x if math.isfinite(x) else None
    return x


def _git_commit():
    """The repository's HEAD commit (git rev-parse HEAD run beside this file), or None when git cannot say."""
    try:
        out = subprocess.run(["git", "rev-parse", "HEAD"], cwd=os.path.dirname(os.path.abspath(__file__)),
                             capture_output=True, text=True, timeout=30)
        sha = out.stdout.strip()
        return sha if out.returncode == 0 and sha else None
    except (OSError, subprocess.SubprocessError):
        return None


def _git_dirty():
    """True when the Reentry Model project folder has uncommitted changes (git status --porcelain over it, untracked
    files included, so a run made from edited code is flagged), False when it is clean, None when git cannot say."""
    project = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    try:
        out = subprocess.run(["git", "status", "--porcelain", "--", project], cwd=project,
                             capture_output=True, text=True, timeout=30)
        return bool(out.stdout.strip()) if out.returncode == 0 else None
    except (OSError, subprocess.SubprocessError):
        return None


def _case_settings(case, settings):
    return {"case": dataclasses.asdict(case), "settings": dataclasses.asdict(settings),
            "run_name": run_name(case, settings)}


def _trail_block(trail):
    return {"L_m": trail.L_m, "z_mean_km": trail.z_mean_km, "N_Al": trail.N_Al, "m_Al_kg": trail.m_Al_kg,
            "m_Zn_kg": trail.m_Zn_kg, "m_Mg_kg": trail.m_Mg_kg}


def _condensed(budget):
    """Everything that left the gas phase (decision 17) from an al_budget dictionary: particles + smoke + settled."""
    return budget["particles"] + budget["smoke"] + budget["settled"]


def _first_time_reaching(t, values, level):
    hit = np.nonzero(np.asarray(values) >= level)[0]
    return float(t[hit[0]]) if hit.size else None


def summarise(sol, trail, hot, case, settings):
    """The summary dictionary of a solved case (everything of summary.json except `runtime_s`, which `run` adds)."""
    n_t = len(sol.t)
    last = n_t - 1
    budgets = [sol.al_budget(i) for i in range(n_t)]
    condensed = np.array([_condensed(bud) for bud in budgets])
    final = budgets[last]

    # mean diameters and the alumina-like share over the plume particles at the final time
    b = sol.bins(last)
    N, d = b["N"], b["d_m"]
    atoms = sum(np.asarray(b["tally"][s][r]) for s in SOURCES for r in ROUTES)
    d_num = d_vol = share = None
    if N.sum() > 0.0:
        d_num = float((N * d).sum() / N.sum()) * 1e9
        d_vol = float((N * d ** 4).sum() / (N * d ** 3).sum()) * 1e9      # volume-weighted mean diameter
    if atoms.sum() > 0.0:
        share = float(atoms[b["H_per_Al"] < ALUMINA_LIKE_BELOW].sum() / atoms.sum())

    reservoirs = {name: np.array([bud[name] for bud in budgets]) for name in AL_SPECIES}
    for name in ("particles", "smoke", "settled"):
        reservoirs[name] = np.array([bud[name] for bud in budgets])

    prov = sol.provenance(last)
    out = _case_settings(case, settings)
    out.update(_trail_block(trail))
    out.update({
        "hot_fractions": dict(hot.fractions), "r0_m": hot.r0_m,
        "share_above_ionisation": hot.share_above_ionisation, "share_above_AlO_bond": hot.share_above_AlO_bond,
        "condensed_fraction_end": float(condensed[last]),
        "t50_s": _first_time_reaching(sol.t, condensed, 0.5),
        "t90_s": _first_time_reaching(sol.t, condensed, 0.9),
        "d_number_mean_nm": d_num, "d_volume_mean_nm": d_vol, "alumina_like_share": share,
        "provenance_end": {f"{s}/{r}": v for (s, r), v in prov.items()},
        "settled_fraction_end": final["settled"], "smoke_fraction_end": final["smoke"],
        "particle_fraction_end": final["particles"],
        "steady_state_time_s": steady_state_time(sol.t, reservoirs),
        "rate_warnings": list(sol.box.warnings),
        "git_commit": _git_commit(), "git_dirty": _git_dirty(), "runtime_s": None})
    return out


# ---- files ----------------------------------------------------------------------------------------------------------

def _cell(v):
    """Text of one CSV cell: empty for None, integers as integers, floats at full precision."""
    if v is None:
        return ""
    if isinstance(v, str):
        return v
    if isinstance(v, (int, np.integer)):
        return str(int(v))
    return repr(float(v))


def _write_csv(path, header, rows):
    with open(path, "w", newline="") as f:
        f.write(",".join(header) + "\n")
        for row in rows:
            f.write(",".join(_cell(v) for v in row) + "\n")


def write_history(path, sol):
    box = sol.box
    prov_keys = [(s, r) for s in SOURCES for r in PROVENANCE_ROUTES]
    laying = getattr(sol, "unreleased", None)
    header = (["t_s", "r_m", "V_m3", "T_K"] + [f"frac_{s}" for s in AL_SPECIES]
              + ["frac_particles", "frac_smoke", "frac_settled", "frac_condensed"]
              + [f"frac_{s}_{r}" for s, r in prov_keys] + [f"c_{s}_m3" for s in HISTORY_AMBIENT]
              + (["frac_unreleased"] if laying else []))
    rows = []
    for i, t in enumerate(sol.t):
        bud = sol.al_budget(i)
        prov = sol.provenance(i)
        rows.append([t, sol.r(i), sol.V(i), box.amb.T_K] + [bud[s] for s in AL_SPECIES]
                    + [bud["particles"], bud["smoke"], bud["settled"], _condensed(bud)]
                    + [prov[k] for k in prov_keys] + [sol.concentration(s, i) for s in HISTORY_AMBIENT]
                    + ([laying(i)] if laying else []))
    _write_csv(path, header, rows)


def size_indices(t, t_end):
    """Indices of the output times nearest to SIZE_TIMES_S (those not past the end of the run), without repeats."""
    t = np.asarray(t, dtype=float)
    idx = []
    for target in SIZE_TIMES_S:
        if target > t_end * (1.0 + 1e-9):
            continue
        k = int(np.argmin(np.abs(np.log(np.maximum(t, 1e-300)) - math.log(target))))
        if k not in idx:
            idx.append(k)
    return idx


def write_size_distribution(path, sol, settings):
    box = sol.box
    pairs = [(s, r) for s in SOURCES for r in ROUTES]
    header = (["t_s", "bin", "n_atoms", "d_nm", "N", "O_per_Al", "H_per_Al", "alumina_like"]
              + [f"Al_{s}_{r}" for s, r in pairs])
    rows = []
    for i in size_indices(sol.t, settings.t_end_s):
        b = sol.bins(i)
        atoms = sum(np.asarray(b["tally"][s][r]) for s, r in pairs)
        for k in range(len(b["N"])):
            if atoms[k] < MASK_FRACTION * box.inventory:
                comp = [None, None, None]
                d_nm = diameter(b["n_atoms"][k], 0.0) * 1e9
            else:
                comp = [b["O_per_Al"][k], b["H_per_Al"][k], float(b["H_per_Al"][k] < ALUMINA_LIKE_BELOW)]
                d_nm = b["d_m"][k] * 1e9
            rows.append([sol.t[i], int(k), b["n_atoms"][k], d_nm, b["N"][k]] + comp
                        + [b["tally"][s][r][k] for s, r in pairs])
    _write_csv(path, header, rows)


def write_segments(path, solved):
    """One row per segment: position, release, hot start and the segment's own shares at the end of its run."""
    total = sum(seg.trail.N_Al for seg, _, _ in solved)
    header = (["segment", "x0_m", "x1_m", "t_laid_s", "z_km", "share_Al", "N_Al", "r0_m", "V_ms", "T_droplet_K"]
              + [f"hot_{s}" for s in sorted(solved[0][2].fractions)] + [f"frac_{s}" for s in AL_SPECIES]
              + ["particles", "smoke", "settled", "condensed"])
    rows = []
    for seg, sol, hot in solved:
        bud = sol.al_budget(len(sol.t) - 1)
        s = seg.trail.samples
        rows.append([seg.index, seg.x0_m, seg.x1_m, seg.t_laid_s, seg.trail.z_mean_km, seg.trail.N_Al / total,
                     seg.trail.N_Al, hot.r0_m, sum(q.V_ms for q in s) / len(s), sum(q.T_d_K for q in s) / len(s)]
                    + [hot.fractions[k] for k in sorted(hot.fractions)] + [bud[k] for k in AL_SPECIES]
                    + [bud["particles"], bud["smoke"], bud["settled"], _condensed(bud)])
    _write_csv(path, header, rows)


def write_along_trail(path, solved, sol, settings):
    """The trail with spreading along it at PROFILE_TIMES_S (those within the run), on PROFILE_POINTS positions."""
    names = list(AL_SPECIES) + ["particles", "smoke"]
    header = ["t_s", "x_m"] + [f"line_{n}_per_m" for n in names] + [f"c_{n}_m3" for n in names]
    whole = max(seg.x1_m for seg, _, _ in solved)
    rows = []
    for t in PROFILE_TIMES_S:
        if t > settings.t_end_s * (1.0 + 1e-9):
            continue
        i = int(np.argmin(np.abs(sol.t - t)))
        reach = PROFILE_REACH * sol.r(i)
        x = np.linspace(-reach, whole + reach, PROFILE_POINTS)
        p = along_trail(solved, t, x)
        for k in range(len(x)):
            rows.append([t, float(x[k])] + [float(p["line"][n][k]) for n in names]
                        + [float(p["axis"][n][k]) for n in names])
    _write_csv(path, header, rows)


def _write_summary(path, summary):
    tmp = path + ".tmp"
    with open(tmp, "w") as f:
        json.dump(_clean(summary), f, indent=2, allow_nan=False)
        f.write("\n")
    os.replace(tmp, path)


# ---- the case -------------------------------------------------------------------------------------------------------

def run(case: Case, settings=None, out_dir=OUTPUT, trail=None, overwrite=False):
    """Run one case and write its outputs under out_dir/<run_name>/; returns the summary dictionary (the finished
    run's own when a valid summary.json exists and `overwrite` is false; an unreadable or truncated one is replaced
    by a fresh run). `trail=` is for tests: it is not encoded in the run name."""
    settings = settings if settings is not None else Settings()
    name = run_name(case, settings)
    d = os.path.join(out_dir, name)
    summary_path = os.path.join(d, "summary.json")
    if os.path.exists(summary_path) and not overwrite:
        try:
            with open(summary_path) as f:
                done = json.load(f)
            if isinstance(done, dict) and "run_name" in done:
                return done
        except (ValueError, OSError):
            pass                                     # truncated or hand-edited: fall through and run again

    t0 = time.perf_counter()
    rows = None
    if trail is None and settings.n_segments > 0:
        rows = rows_for(case.d0_um, case.h0_km)
        trail = trail_from_rows(rows, case.d0_um, case.h0_km, n_samples=settings.n_emission_samples)
    elif trail is None:
        trail = trail_for(case.d0_um, case.h0_km, settings.n_emission_samples)
    os.makedirs(d, exist_ok=True)
    if not trail.N_Al > 0.0:
        summary = _case_settings(case, settings)
        summary.update(_trail_block(trail))
        summary.update({"skipped": SKIPPED_NO_AL, "git_commit": _git_commit(),
                        "git_dirty": _git_dirty(),
                        "runtime_s": time.perf_counter() - t0})
        _write_summary(summary_path, summary)
        return _clean(summary)

    species = load_species()
    column = AmbientColumn.load(settings.ambient, case.time_of_day)
    amb = column.at(trail.z_mean_km)
    channels = calibrate(load_reactions(stage=("hot", "both")), species)
    extra = {}
    if rows is not None:
        segments = cut_trail(rows, case.d0_um, case.h0_km, settings.n_segments, settings.n_emission_samples)
        solved = solve_segments(segments, column, settings, species, channels, case.time_of_day)
        hot = combined_hot(solved)
        sol = TrailSolution(solved, solved[0][1].box.time_grid())
        write_segments(os.path.join(d, "segments.csv"), solved)
        write_along_trail(os.path.join(d, "along_trail.csv"), solved, sol, settings)
        extra = {"n_segments": len(segments), "trail_length_m": segments[0].trail.L_m}
    else:
        hot = slow_down(trail.samples, column, settings, species, channels)
        sol = Box(trail, hot, amb, case.time_of_day, settings, species=species).solve()

    write_history(os.path.join(d, "history.csv"), sol)
    write_size_distribution(os.path.join(d, "size_distribution.csv"), sol, settings)
    summary = summarise(sol, trail, hot, case, settings)
    summary.update(extra)
    summary["runtime_s"] = time.perf_counter() - t0
    _write_summary(summary_path, summary)            # last: a crash above leaves no summary, so a sweep re-runs it
    return _clean(summary)

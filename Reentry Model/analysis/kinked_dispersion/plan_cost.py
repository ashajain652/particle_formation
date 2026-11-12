"""The cost check of the kinked-profile dispersion (Task 8, spec section 5 and stage 4): how much longer a macro step's
spray plan takes on the `profile` path than on `table`, against the wall time of a macro step of the flight the state
came from. The budget is 5 % of a macro step's wall time.

    "$PY" -m analysis.kinked_dispersion.plan_cost --state <flight>.state.pkl      # from Reentry Model/

The plan (`skinbody.girin_plan`) is evaluated once per macro step, at its start. It is called REPEATS times with the
saved body's spray model on `table` and REPEATS times on `profile`, alternately so that a change in the machine's load
falls on both, on the body's own inputs: its skins, its last surface-flow evaluation (`body.last_flow`), the aero state,
angles and top areas its last plan was evaluated with (`body.spray_plan`), its columns' graph, its transverse radius, its
nose radius, its exposed-slurry setting and the run's macro step. The added time is the median time on `profile` less the
median on `table` [s].

The macro step's wall time is read from the flight's rows file (`--rows`, by default the state's name with `.rows.csv`
for `.state.pkl`): the increase, between consecutive rows from --start to --end s of flight, of the wall time spent in
the skins, the 3D solve, the mesh's moves and its rebuilds (wall_skins_s + wall_3d_s + wall_move_s + wall_rebuild_s,
which accumulate over the flight), one increase per macro step; its median [s]. The flight's other wall time (the
trajectory, the writer) is left out, which makes the macro step shorter and the check stricter.

The spray model's setting is changed on the loaded copy only; the state on disk is read, never written."""
import argparse
import csv
import os
import time

import numpy as np

from reentry_model.skyscrapers import restart, skinbody

REPEATS = 5                       # calls of the plan on each path
BUDGET = 0.05                     # the added time's budget, as a fraction of a macro step's wall time (spec section 5)
WALL_COLUMNS = ("wall_skins_s", "wall_3d_s", "wall_move_s", "wall_rebuild_s")
WINDOW_S = (177.5, 179.5)         # s of flight: the d050 flight's 0.5 s macro steps around its saved state (178.5 s)


def set_dispersion(model, setting):
    """Put the spray model on `setting` as the command line does: `table` is the class default, stored on no instance."""
    if setting == "table":
        vars(model).pop("dispersion", None)
    else:
        model.dispersion = setting


def plan_once(body, dt):
    """One evaluation of the macro step's spray plan on the body's own inputs: (seconds, the plan)."""
    p = body.spray_plan
    start = time.perf_counter()
    plan = skinbody.girin_plan(body.skins, body.last_flow, p.state, body.spray, p.theta, p.area, body.columns.edges,
                               body.transverse_radius(), dt, body.slurry_spray, nose_radius=body.nose_radius())
    return time.perf_counter() - start, plan


def plan_times(body, dt, repeats=REPEATS):
    """{setting: [seconds of each call]} and the last plan on each setting, the settings called alternately."""
    times, plans = {"table": [], "profile": []}, {}
    for _ in range(repeats):
        for setting in ("table", "profile"):
            set_dispersion(body.spray, setting)
            seconds, plans[setting] = plan_once(body, dt)
            times[setting].append(seconds)
    set_dispersion(body.spray, "table")
    return times, plans


def step_wall_times(rows_path, start=WINDOW_S[0], end=WINDOW_S[1]):
    """[(t0, t1, seconds)]: the increase of the summed WALL_COLUMNS between consecutive rows whose times lie within
    [start, end] s, one per macro step."""
    with open(rows_path, newline="") as fh:
        rows = [r for r in csv.DictReader(fh) if start - 1e-9 <= float(r["time_s"]) <= end + 1e-9]
    total = [(float(r["time_s"]), sum(float(r[k]) for k in WALL_COLUMNS)) for r in rows]
    return [(t0, t1, w1 - w0) for (t0, w0), (t1, w1) in zip(total, total[1:])]


def report(state_path, rows_path, start=WINDOW_S[0], end=WINDOW_S[1], repeats=REPEATS):
    state = restart.load(state_path)
    run = state["run"]
    body, dt = run.body, run.settings.dt
    if body.spray is None or body.spray_plan is None or body.last_flow is None:
        raise SystemExit("{} has no spray plan to repeat (a flight of spraying skins saved after a macro step)".format(
            state_path))
    print("state {} at {:.6g} s (step {}): {} columns, macro step {:g} s, spray model saved on {}".format(
        state_path, state["time_s"], state["step"], body.columns.n, dt, body.spray.dispersion))
    times, plans = plan_times(body, dt, repeats)
    med = {k: float(np.median(v)) for k, v in times.items()}
    added = med["profile"] - med["table"]
    for setting in ("table", "profile"):
        print("girin_plan on {}: {} s; median {:.3f} s".format(
            setting, ", ".join("{:.3f}".format(t) for t in times[setting]), med[setting]))
    resolved = plans["profile"].resolved
    print("profile path: {} columns resolved, {} of them floor-like".format(
        int(resolved.sum()), int(np.nansum(plans["profile"].record.wall_floor[resolved]))))
    print("added time of the plan on profile over table (median less median): {:.3f} s".format(added))

    steps = step_wall_times(rows_path, start, end)
    if not steps:
        raise SystemExit("{} has fewer than two rows from {:g} to {:g} s".format(rows_path, start, end))
    for t0, t1, w in steps:
        print("macro step {:g} to {:g} s: {:.1f} s of wall time ({})".format(t0, t1, w, " + ".join(WALL_COLUMNS)))
    wall = float(np.median([w for _, _, w in steps]))
    share = added / wall
    print("macro step's wall time, median of {} steps: {:.1f} s".format(len(steps), wall))
    print("added time {:.3f} s is {:.3f} % of it, against a budget of {:g} %: {}".format(
        added, 100.0 * share, 100.0 * BUDGET, "PASS" if share < BUDGET else "FAIL"))
    return added, wall, share


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--state", required=True, help="a spraying skin flight's saved state (restart.save); read only")
    p.add_argument("--rows", default=None, help="the flight's rows file (default: the state's name with .rows.csv)")
    p.add_argument("--start", type=float, default=WINDOW_S[0], help="[s] of flight, first row of the wall-time window")
    p.add_argument("--end", type=float, default=WINDOW_S[1], help="[s] of flight, last row of the wall-time window")
    p.add_argument("--repeats", type=int, default=REPEATS, help="calls of the plan on each path")
    a = p.parse_args(argv)
    rows = a.rows
    if rows is None:
        if not a.state.endswith(".state.pkl"):
            p.error("--rows is needed for a state whose name does not end in .state.pkl")
        rows = a.state[:-len(".state.pkl")] + ".rows.csv"
    if not os.path.isfile(rows):
        p.error("no rows file at {}".format(rows))
    report(a.state, rows, a.start, a.end, a.repeats)


if __name__ == "__main__":
    main()

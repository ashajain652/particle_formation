"""Spheral M2 Task 6b: the replay's crossings-error monitor (plan decision 7) validated against check 4 and measured
on the real flight.

    "$PY" analysis/spheral_m2_monitor.py [--only run|analyse]

**The estimator** (`run.crossings_error`): on a monitored interval the replay reruns the accepted interval at twice
the crossings and takes err = 2 max_i vM(sigma_N,i - sigma_2N,i) / max_i vM(sigma_N,i), with the RMS and the 99th
percentile over particles in place of the maximum for context. The factor 2 assumes that the quasi-static lag falls
as 1/N (sigma_N - sigma_inf = 2 (sigma_N - sigma_2N)).

**Validation against check 4** (`validate`). Check 4's heated sphere (the ramp in over interval 0 -> 1, then held;
`analysis/spheral_m2_check4.py`) has the exact solution and runs at 5, 10 and 20 crossings for four forms and
spacings, so at k = 1 (the end of the ramp, the monitor's first interval) per case:

  est_10      the estimator at N = 10 (sigma_10 against sigma_20; what the monitor reports on interval 0 -> 1)
  est_5       the estimator at N = 5 (sigma_5 against sigma_10)
  ratio_5_10  est_5 / est_10: 2 if the lag falls as 1/N in the same spatial pattern
  pressure, pressure_5   the hydrostatic part the von Mises estimate leaves out, 2 max |P_N - P_2N| / max vM(sigma_N)
              (`run.pressure_error`), at N = 10 and 5
  exact_N     the same per-particle measure against the exact solution, max vM(sigma_N - sigma_exact) / max
              vM(sigma_N), N = 5, 10, 20: the whole error, dominated by the surface particles and SPH's particle
              scatter (0.60-0.63 everywhere; 0.16-0.18 inside R - 2 dx), so per particle the exact solution cannot
              resolve the lag; check 4's binned interior errors (below) can
  limit, limit_basis   max(est, pressure) and which decided: what the monitor's limit keys on (controller, after
              Task 6b), computed here from the stored estimates
  binned_lag  2 x the binned interior difference |sigma_10 - sigma_20| over sigma_0 (check 4's binned measure,
              `spheral_m2_check4.pair_diff`), and check 4's interior errors at 10 and 20 and their difference
              (the "lag" check 4 measured: 0.6 % at RZ 1.1 mm)

and at k = 3 (held for 1 s) the same where the runs exist. The monitor itself is run on check 4's RZ 2.2 mm case
(`--monitor-every 1`, 10 crossings, damping 0.05, the check's other settings): its interval 0 -> 1 estimate must equal
est_10 computed here from check 4's 10- and 20-crossing runs (the same two states, bitwise), and its later intervals
(held) show what it reports once the ramp is over. The two monitor runs were made while the limit keyed on the von
Mises estimate alone; the `limit` column applies the present rule (the larger of the two) to their stored estimates,
which needs no rerun. Under it check 4's RZ 2.2 mm ramp interval (limit 6.1 %, pressure) would have switched to 20
crossings, so that run's later rows describe the run at 10; the flight (largest 1.45 %, von Mises) is unaffected.

**The real flight** (`flight`): the MVP flight's RZ 2.2 mm replay body, frames 0-20, `--monitor-every 2` (the default
is 10; 2 for the measurement), every other setting at its default. Per monitored interval the estimates, the
interval's and the rerun's wall time; the overhead = monitor wall / accepted wall over the run.

Outputs: `data/spheral/m2_monitor.csv` (committed; one row per validation case and per monitored flight interval),
the runs in `spheral_output/m2/monitor/runs/`. Run under `caffeinate -i -s` (the Mac's sleep suspends Docker's VM)."""
from __future__ import annotations

import argparse
import csv
import json
import math
import os
import subprocess
import sys

import numpy as np

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, REPO)
sys.path.insert(0, os.path.join(REPO, "analysis"))

import spheral_m2_check4 as c4  # noqa: E402
from spheral_frag import analytic, material, naming  # noqa: E402
from spheral_frag.runner import checkpoint, run  # noqa: E402

OUT = os.path.join("spheral_output", "m2", "monitor")
RUNS_OUT = os.path.join(OUT, "runs")
CSV = os.path.join("data", "spheral", "m2_monitor.csv")
LAUNCHER = os.path.join("spheral_frag", "container", "spheral")
FLIGHT_BODY = os.path.join("spheral_output", "body", naming.body_name(os.path.basename(c4.PREPARED), 0, 47, "rz", 2.2))
FLIGHT_FRAMES = (0, 20)
FLIGHT_EVERY = 2
CHECK4_EVERY = 1

COLUMNS = ["case", "form", "dx_mm", "k0", "k1", "t1_s", "crossings", "monitor_crossings", "est", "rms", "p99",
           "pressure", "limit", "limit_basis", "est_5", "pressure_5", "ratio_5_10", "exact_5", "exact_10", "exact_20", "binned_lag", "interior_err_10",
           "interior_err_20", "interior_diff_10_20", "monitor_run_est", "wall_s", "monitor_wall_s",
           "ke_internal_ratio", "monitor_ke_internal_ratio", "run"]


def check4_monitor_name():
    b = naming.canonical_brackets({"strength": "elastic", "sampler": "radial", "loads": "off",
                                   "monitor_every": CHECK4_EVERY})
    return naming.run_name(os.path.basename(c4.PREPARED), "replay", c4.FRAMES, "rz", 2.2, b, 1)


def flight_name():
    b = naming.canonical_brackets({"monitor_every": FLIGHT_EVERY})
    return naming.run_name(os.path.basename(c4.PREPARED), "replay", FLIGHT_FRAMES, "rz", 2.2, b, 1)


def launch(args, name):
    d = os.path.join(REPO, RUNS_OUT, name)
    if os.path.isfile(os.path.join(d, "summary.json")):
        return
    cmd = [LAUNCHER, "spheral_frag/runner/run.py", "replay", *args, "--out", RUNS_OUT]
    if os.path.isdir(d):
        cmd.append("--resume")
    print(" ".join(cmd), flush=True)
    r = subprocess.run(cmd, cwd=REPO, capture_output=True, text=True)
    os.makedirs(os.path.join(REPO, OUT, "logs"), exist_ok=True)
    with open(os.path.join(REPO, OUT, "logs", name[-60:] + ".log"), "w") as fh:
        fh.write(" ".join(cmd) + "\n" + r.stdout + "\n--- stderr\n" + r.stderr)
    if r.returncode != 0:
        raise SystemExit("run failed ({}): {}".format(r.returncode, r.stderr[-2000:]))


def run_all():
    launch(["--body", c4.body_dir("rz", 2.2), "--prepared", c4.PREPARED, "--frames", "{}-{}".format(*c4.FRAMES),
            "--strength", "elastic", "--loads", "off", "--ke-fraction", "1e300",
            "--monitor-every", str(CHECK4_EVERY)], check4_monitor_name())
    launch(["--body", FLIGHT_BODY, "--prepared", c4.PREPARED, "--frames", "{}-{}".format(*FLIGHT_FRAMES),
            "--monitor-every", str(FLIGHT_EVERY)], flight_name())


# ------------------------------------------------------------------------------------------------------ measuring
def exact_tensor(st, E, nu, alpha):
    """The exact Timoshenko-Goodier tensor at each particle's reference position, in stress_state's frame."""
    x = np.asarray(st["x_ref"], dtype=float)
    p = np.zeros((len(x), 3))
    p[:, :x.shape[1]] = x
    r = np.sqrt(np.einsum("ij,ij->i", p, p))
    n = np.zeros_like(p)
    ok = r > 0.0
    n[ok] = p[ok]/r[ok, None]
    n[~ok, 0] = 1.0
    srr, stt = analytic.sphere_quadratic_stress(r, c4.R, E, nu, alpha, c4.DELTA_K)
    return stt[:, None, None]*np.eye(3)[None] + (srr - stt)[:, None, None]*np.einsum("ni,nj->nij", n, n)


def state(d, k, table):
    return run.stress_state(checkpoint.read(os.path.join(d, "checkpoints", checkpoint.filename(k))), table)


def validate(table):
    E, nu, alpha = c4.constants(table)
    s0 = analytic.quadratic_peak(E, nu, alpha, c4.DELTA_K)
    c4rows = {}
    with open(os.path.join(REPO, c4.CSV), newline="") as fh:
        for r in csv.DictReader(fh):
            if r["scaled"] == "1" and r["damping"] == "0.05":
                c4rows[(r["form"], float(r["dx_mm"]), int(r["crossings"]), int(r["k"]))] = r
    mon_dir = os.path.join(REPO, RUNS_OUT, check4_monitor_name())
    mon_rows = {int(r["k0"]): r for r in run.read_intervals(os.path.join(mon_dir, "intervals.csv"))}
    rows = []
    for form, dx, _ in c4.FORMS:
        for k in c4.COMPARE_K:
            st = {n: state(c4.run_dir(form, dx, n, 0.05, True), k, table) for n in (5, 10, 20)}
            est10 = run.crossings_error(st[10], st[20])
            est5 = run.crossings_error(st[5], st[10])
            pr10, pr5 = run.pressure_error(st[10], st[20]), run.pressure_error(st[5], st[10])
            ex = exact_tensor(st[10], E, nu, alpha)
            exact = {}
            for n in (5, 10, 20):
                o = np.argsort(st[n]["pid"])
                assert np.array_equal(st[n]["pid"][o], st[10]["pid"])
                exact[n] = float(run.von_mises(st[n]["sigma"][o] - ex).max()/run.von_mises(st[n]["sigma"]).max())
            meas = {n: c4.measure(c4.run_dir(form, dx, n, 0.05, True), k, table, E, nu, alpha, dx*1e-3)
                    for n in (10, 20)}
            lag = 2.0*c4.pair_diff(meas[10], meas[20], dx*1e-3)[0]
            i10 = float(c4rows[(form, dx, 10, k)]["interior_err"])
            i20 = float(c4rows[(form, dx, 20, k)]["interior_err"])
            mr = mon_rows.get(k - 1) if (form, dx) == ("rz", 2.2) else None
            rows.append(dict(case="check4", form=form, dx_mm=dx, k0=k - 1, k1=k, t1_s=0.5*k, crossings=10,
                             monitor_crossings=20, est=est10[0], rms=est10[1], p99=est10[2], pressure=pr10,
                             est_5=est5[0], pressure_5=pr5,
                             ratio_5_10=est5[0]/est10[0], exact_5=exact[5], exact_10=exact[10], exact_20=exact[20],
                             binned_lag=lag, interior_err_10=i10, interior_err_20=i20,
                             interior_diff_10_20=i10 - i20,
                             monitor_run_est=float(mr["crossings_err_est"]) if mr else math.nan,
                             wall_s=float(mr["wall_s"]) if mr else math.nan,
                             monitor_wall_s=float(mr["monitor_wall_s"]) if mr else math.nan,
                             ke_internal_ratio=float(c4rows[(form, dx, 10, k)]["ke_internal_ratio"]),
                             monitor_ke_internal_ratio=float(c4rows[(form, dx, 20, k)]["ke_internal_ratio"]),
                             run=os.path.basename(c4.run_dir(form, dx, 10, 0.05, True))))
    # the monitor run on check 4's RZ 2.2 mm case: every interval as the monitor reported it
    for k0, r in sorted(mon_rows.items()):
        rows.append(dict(case="check4_monitor", form="rz", dx_mm=2.2, k0=k0, k1=int(r["k1"]), t1_s=float(r["t1_s"]),
                         crossings=int(r["crossings"]), monitor_crossings=int(r["monitor_crossings"]),
                         est=float(r["crossings_err_est"]), rms=float(r["crossings_err_rms"]),
                         p99=float(r["crossings_err_p99"]), pressure=float(r["crossings_err_pressure"]),
                         wall_s=float(r["wall_s"]),
                         monitor_wall_s=float(r["monitor_wall_s"]), ke_internal_ratio=float(r["ke_internal_ratio"]),
                         monitor_ke_internal_ratio=float(r["monitor_ke_internal_ratio"]), run=check4_monitor_name()))
    return rows, s0


def flight():
    d = os.path.join(REPO, RUNS_OUT, flight_name())
    rows = []
    with open(os.path.join(d, "summary.json")) as fh:
        summ = json.load(fh)
    for r in run.read_intervals(os.path.join(d, "intervals.csv")):
        rows.append(dict(case="flight", form="rz", dx_mm=2.2, k0=int(r["k0"]), k1=int(r["k1"]), t1_s=float(r["t1_s"]),
                         crossings=int(r["crossings"]), monitor_crossings=int(r["monitor_crossings"]),
                         est=float(r["crossings_err_est"]), rms=float(r["crossings_err_rms"]),
                         p99=float(r["crossings_err_p99"]), pressure=float(r["crossings_err_pressure"]),
                         wall_s=float(r["wall_s"]),
                         monitor_wall_s=float(r["monitor_wall_s"]), ke_internal_ratio=float(r["ke_internal_ratio"]),
                         monitor_ke_internal_ratio=float(r["monitor_ke_internal_ratio"]), run=flight_name()))
    return rows, summ


def analyse():
    table = material.MaterialTable.load(os.path.join(REPO, c4.PREPARED))
    vrows, s0 = validate(table)
    frows, summ = flight()
    rows = vrows + frows
    for r in rows:      # the limit keys on the larger of the two (controller, after Task 6b), from the stored values
        r["limit"], r["limit_basis"] = run.limit_error(r["est"], r["pressure"])
    path = os.path.join(REPO, CSV)
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(COLUMNS)
        for r in rows:
            w.writerow([repr(float(r[c])) if isinstance(r.get(c), (float, np.floating)) else r.get(c, "nan")
                        for c in COLUMNS])
    print("wrote {} ({} rows); sigma_0 {:.4g} MPa".format(CSV, len(rows), s0/1e6))
    print("\ncheck 4 validation (N = 10 against 20; fractions)")
    hdr = ("frm", "dx", "k", "est10", "rms", "p99", "pres10", "limit10", "est5", "pres5", "5/10", "exact10", "binlag",
           "int10", "int20", "mon.run")
    print("{:3s} {:>4s} {:>2s}".format(*hdr[:3]) + "".join(" {:>8s}".format(h) for h in hdr[3:]))
    for r in vrows:
        if r["case"] != "check4":
            continue
        vals = [r[c] for c in ("est", "rms", "p99", "pressure", "limit", "est_5", "pressure_5", "ratio_5_10",
                               "exact_10", "binned_lag", "interior_err_10", "interior_err_20", "monitor_run_est")]
        print("{:3s} {:4.1f} {:2d}".format(r["form"], r["dx_mm"], r["k1"]) + "".join(" {:8.4f}".format(v)
                                                                                for v in vals))
    print("\ncheck 4 RZ 2.2 mm with --monitor-every 1")
    for r in vrows:
        if r["case"] == "check4_monitor":
            print("  k {}->{}  est {:.4f}  rms {:.4f}  p99 {:.4f}  pressure {:.4f}  limit {:.4f} ({})  wall {:.1f} s  "
                  "monitor {:.1f} s".format(r["k0"], r["k1"], r["est"], r["rms"], r["p99"], r["pressure"], r["limit"],
                                            r["limit_basis"], r["wall_s"],
                          r["monitor_wall_s"]))
    print("\nflight RZ 2.2 mm, frames {}-{}, --monitor-every {}".format(*FLIGHT_FRAMES, FLIGHT_EVERY))
    for r in frows:
        if r["monitor_crossings"]:
            print("  k {:2d}->{:2d} t {:5.2f} s  N {}  est {:.4f}  rms {:.4f}  p99 {:.4f}  pressure {:.4f}  limit "
                  "{:.4f} ({})  wall {:.1f} s  monitor {:.1f} s".format(
                      r["k0"], r["k1"], r["t1_s"], r["crossings"], r["est"], r["rms"], r["p99"], r["pressure"],
                      r["limit"], r["limit_basis"], r["wall_s"], r["monitor_wall_s"]))
    wall = math.fsum(r["wall_s"] for r in frows)
    mon = math.fsum(r["monitor_wall_s"] for r in frows)
    print("  accepted {:.1f} s, monitor {:.1f} s: overhead {:.1%} at every {} ({} monitored of {}); at the default 10 "
          "about {:.1%}".format(wall, mon, mon/wall, FLIGHT_EVERY, sum(1 for r in frows if r["monitor_crossings"]),
                                len(frows), (mon/max(1, sum(1 for r in frows if r["monitor_crossings"])))*sum(
                                    run.monitor_due(j, 10) for j in range(len(frows)))/wall))
    print("  max limit estimate {:.4f} ({}); the run's summary.json (von Mises rule): switches {}".format(
        max(r["limit"] for r in frows if r["monitor_crossings"]),
        max((r for r in frows if r["monitor_crossings"]), key=lambda r: r["limit"])["limit_basis"],
        summ["crossings_monitor"]["switches"]))
    return rows


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--only", choices=("run", "analyse"))
    a = p.parse_args(argv)
    os.chdir(REPO)
    if a.only != "analyse":
        run_all()
    if a.only != "run":
        analyse()
    return 0


if __name__ == "__main__":
    sys.exit(main())

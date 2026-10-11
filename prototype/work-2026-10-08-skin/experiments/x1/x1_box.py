"""Experiment X1 (Asha's decision 2026-10-09, option a): two candidate fixes for the skins' stall in the box, measured.

Runs the box harness of analysis/box/box_ablation.py (imported, not edited: its run_specs, run_one, classify and
collect_variant) on the experiment copy experiments/x1/code instead of the frozen meas/ copy, with the two MeltSettings
switches of the copy:
  draw_through  -- a skin whose owner element is exhausted within a step goes on drawing from the element behind it
  cover_exposed -- faces a death exposes get no skin by the solidus test (only by the hand-over, or when liquid arrives)

Variant q6p (6 MW/m2 heated by projection, 2 s, Scheil, 700 K), dev preset (20 um cells), sub-step 10 ms.

    python x1_box.py --list
    python x1_box.py base__dt0.05 both__h1mm ...      # run the named runs, one JSON each under runs/
    python x1_box.py --collect                         # runs/*.json -> results.json, summary.md, x1_box.png

Each run records, beyond box_ablation's series, per step: the skins on the lateral interior facing the heat
(n_x < -0.1) and not facing it (n_x >= -0.1), the lateral liquid fed, the skins created, the metal drawn through, and
the faces the deaths exposed.
"""
import glob
import json
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
CODE = os.path.join(HERE, "code")
BOX = os.path.normpath(os.path.join(HERE, "..", "..", "analysis", "box"))
RUNS = os.path.join(HERE, "runs")
sys.dont_write_bytecode = True
sys.path.insert(0, CODE)
import reentry_model  # noqa: E402  the experiment copy, imported before box_ablation puts meas/ on the path

sys.path.insert(0, BOX)
import box_ablation as bx  # noqa: E402

import numpy as np  # noqa: E402

from reentry_model import body  # noqa: E402

for _m in (reentry_model, bx.body, bx.skin, bx.mesh, bx.thermal, bx.material, bx.heating, body):
    assert os.path.realpath(_m.__file__).startswith(os.path.realpath(CODE)), _m.__file__

TAG = "q6p"
V = bx.VARIANTS[TAG]
CONFIGS = {"base": (False, False), "draw_through": (True, False), "cover_exposed": (False, True), "both": (True, True)}
MACRO = bx.MACRO
FEED_ZERO = 0.0          # kg: a step whose lateral feed is at or below this fed nothing (exactly zero)
SLOW_RATE = 1.1e-3       # m/s of lateral recession: a quarter of the fine box's mean rate between 1 and 2 s (4.40 mm/s); a
                         # step feeding slower than this is "slow" (a second, softer measure of the stall)
ORIG_RUNS = os.path.join(BOX, "runs")
COLOURS = {"base": "tab:blue", "draw_through": "tab:orange", "cover_exposed": "tab:green", "both": "tab:red"}


def specs():
    base = bx.run_specs(TAG, V["q"], V["rule"])
    out = {}
    for cfg, (dthr, cov) in CONFIGS.items():
        for dt in MACRO:
            s = dict(base["skin_dev_dt{:g}".format(dt)])
            s.update(draw_through=dthr, cover_exposed=cov, config=cfg, orig="q6p__skin_dev_dt{:g}".format(dt))
            out["{}__dt{:g}".format(cfg, dt)] = s
        if cfg in ("base", "both"):
            for hh in (1.0, 0.5):
                s = dict(base["diag_skin_dev_dt0.05_h{:g}mm".format(hh)])
                s.update(draw_through=dthr, cover_exposed=cov, config=cfg, orig="q6p__diag_skin_dev_dt0.05_h{:g}mm".format(hh))
                out["{}__h{:g}mm".format(cfg, hh)] = s
    return out


_orig_make_body = bx.make_body
_last = {}


def make_body(spec):
    """box_ablation's body, the switches set on its settings (nothing in the construction reads them), and its advance
    wrapped to record the X1 diagnostics after every step (before the harness removes the liquid)."""
    b, mat = _orig_make_body(spec)
    b.settings.draw_through = bool(spec.get("draw_through", False))
    b.settings.cover_exposed = bool(spec.get("cover_exposed", False))
    lx, ly, lz, h = spec["box"]
    rec = {k: [] for k in ("t", "skins_lat_away", "skins_lat_away_exposed", "skins_lat_away_back", "skins_lat_heated", "skins_total", "skins_created", "fed_lat_kg",
                           "drawn_kg", "drawn_through_kg", "handover_kg", "exposed_faces", "skins_on_exposed", "n_active")}
    inner = b.advance

    def advance(t, dt, loads, state=None):
        inner(t, dt, loads, state)
        _, _, _, lateral = bx.classify(b, spec["box"])
        s = b.surface
        rec["t"].append(t)
        if b.skins is not None and b.skins.n:
            p = b.patch_of_skin
            lat = lateral[p]
            nx = s.normals[p, 0]
            on_exp = np.isin(b.skins.face_id, np.fromiter(b.exposed_faces, dtype=np.int64, count=len(b.exposed_faces)))
            away = lat & (nx >= -0.1)
            rec["skins_lat_away"].append(int(away.sum()))
            rec["skins_lat_away_exposed"].append(int((away & on_exp).sum()))       # on faces a death exposed
            rec["skins_lat_away_back"].append(int((away & (s.centroids[p, 0] > lx * (1.0 - 1e-9))).sum()))   # the x = lx face
            rec["skins_lat_heated"].append(int((lat & (nx < -0.1)).sum()))
            rec["skins_on_exposed"].append(int(on_exp.sum()))
        else:
            rec["skins_lat_away"].append(0)
            rec["skins_lat_away_exposed"].append(0)
            rec["skins_lat_away_back"].append(0)
            rec["skins_lat_heated"].append(0)
            rec["skins_on_exposed"].append(0)
        rec["skins_total"].append(int(b.skins.n) if b.skins is not None else 0)
        rec["skins_created"].append(int(b.skins_created))
        rec["fed_lat_kg"].append(float((b.m_f + b.m_d)[lateral].sum()))
        rec["drawn_kg"].append(float(b.skin_drawn_mass))
        rec["drawn_through_kg"].append(float(b.skin_drawn_through_mass))
        rec["handover_kg"].append(float(b.skin_handover_mass))
        rec["exposed_faces"].append(len(b.exposed_faces))
        rec["n_active"].append(int(b.mesh.n_active))

    b.advance = advance
    _last["body"], _last["rec"] = b, rec
    return b, mat


bx.make_body = make_body     # run_one looks make_body up in its module at call time


def run(name, spec):
    print("== X1 {} ({}) ==".format(name, spec["config"]), flush=True)
    try:
        r = bx.run_one(name, spec, log=lambda s: print(s, flush=True))
        b, x = _last["body"], _last["rec"]
        r["x1"] = x
        r["mass_balance_rel"] = abs(r["mass_balance"]) / r["mass0_kg"]
        r["skin_drawn_mass_kg"] = float(b.skin_drawn_mass)
        r["skin_drawn_through_mass_kg"] = float(b.skin_drawn_through_mass)
        r["skin_handover_mass_kg"] = float(b.skin_handover_mass)
        r["exposed_faces_total"] = len(b.exposed_faces)
        r["code"] = os.path.dirname(reentry_model.__file__)
    except Exception as exc:
        import traceback
        traceback.print_exc()
        r = {"name": name, "tag": TAG, "status": "failed", "error": "{}: {}".format(type(exc).__name__, exc)}
    r["config"] = spec["config"]
    r["x1_spec"] = {k: (list(v) if isinstance(v, tuple) else v) for k, v in spec.items()}
    os.makedirs(RUNS, exist_ok=True)
    with open(os.path.join(RUNS, name + ".json"), "w") as f:
        json.dump(r, f, indent=1)


# -- collection -----------------------------------------------------------------------------------------------------
def stall(t, fed, dt):
    """The longest interval with no lateral feed after the recession starts: (duration [s], from [s], to [s]). A step
    feeds nothing when its lateral liquid is at or below FEED_ZERO; the interval runs from the end of the last step that
    fed to the end of the last that did not."""
    fed = np.asarray(fed)
    on = np.flatnonzero(fed > FEED_ZERO)
    if not on.size:
        return None, None, None
    best, start = (0, None), None
    for k in range(on[0] + 1, len(fed)):
        if fed[k] <= FEED_ZERO:
            start = k if start is None else start
            n = k - start + 1
            if n > best[0]:
                best = (n, start)
        else:
            start = None
    if best[0] == 0:
        return 0.0, None, None
    n, k0 = best
    return n * dt, t[k0] - dt, t[k0 + n - 1]


def slow(t, fed, dt, A_lat):
    """The longest interval after the recession starts in which every step's lateral recession rate is below SLOW_RATE:
    (duration [s], from [s], to [s])."""
    rate = np.asarray(fed) / (A_lat * bx.RHO_METAL) / dt
    on = np.flatnonzero(np.asarray(fed) > FEED_ZERO)
    if not on.size:
        return None, None, None
    best, start = (0, None), None
    for k in range(on[0] + 1, len(rate)):
        if rate[k] < SLOW_RATE:
            start = k if start is None else start
            if k - start + 1 > best[0]:
                best = (k - start + 1, start)
        else:
            start = None
    if best[0] == 0:
        return 0.0, None, None
    n, k0 = best
    return n * dt, t[k0] - dt, t[k0 + n - 1]


def at_time(t, v, x):
    i = int(np.argmin(np.abs(np.asarray(t) - x)))
    return v[i] if abs(t[i] - x) < 1e-9 else None


def collect():
    with open(os.path.join(ORIG_RUNS, "q6p__fine.json")) as f:
        fine = json.load(f)
    recs = {}
    for path in sorted(glob.glob(os.path.join(RUNS, "*.json"))):
        with open(path) as f:
            r = json.load(f)
        recs[r["name"]] = r
    res = bx.collect_variant(TAG, dict(recs, fine=fine))
    sp = specs()
    out = {"fine_recession_mm_at_times": res["fine_recession_mm_at_times"], "compare_times_s": res["compare_times_s"],
           "criterion": "recession within 5 % of the fine box's at 0.5, 1.0, 1.5 and 2.0 s (where the fine box has >= 0.05 mm)",
           "runs": {}}
    for name, r in recs.items():
        o = res["runs"][name]
        if r.get("status") != "ok":
            out["runs"][name] = {"status": "failed", "error": r.get("error"), "config": r.get("config")}
            continue
        S, x = r["series"], r["x1"]
        dt = r["spec"]["dt"]
        d, a, b_ = stall(S["t"], x["fed_lat_kg"], dt)
        lx_, ly_, lz_, hh_ = r["spec"]["box"]
        ds, as_, bs_ = slow(S["t"], x["fed_lat_kg"], dt, (ly_ - 2 * hh_) * (lz_ - 2 * hh_))
        rel = o["r_lat_mm_rel_diff_at_times"]
        e = {"status": "ok", "config": r["config"], "dt": dt, "box": r["spec"]["box"],
             "r_lat_mm_at_times": o["r_lat_mm_at_times"], "rel_diff_at_times": rel,
             "max_abs_rel_diff_4times": o["r_lat_mm_max_abs_rel_diff"],
             "max_abs_rel_diff_curve": o["r_lat_mm_max_abs_rel_diff_curve"],
             "max_abs_diff_curve_mm": o["r_lat_mm_max_abs_diff_curve_mm"],
             "meets_5pct": all(abs(v) <= 0.05 for v in rel if v is not None),
             "stall_s": d, "stall_from_s": a, "stall_to_s": b_,
             "slow_s": ds, "slow_from_s": as_, "slow_to_s": bs_,
             "onset_s": S["t"][int(np.flatnonzero(np.asarray(x["fed_lat_kg"]) > FEED_ZERO)[0])] if np.any(np.asarray(x["fed_lat_kg"]) > FEED_ZERO) else None,
             "worst_energy_residual": r["worst_energy_residual"], "mass_balance_rel": r["mass_balance_rel"],
             "repeats_per_step": r["interface_repeats_total"] / r["n_steps"], "max_shortfall_kg": r["max_shortfall_kg"],
             "skins_lat_away_at_1s": at_time(x["t"], x["skins_lat_away"], 1.0),
             "skins_lat_away_at_2s": at_time(x["t"], x["skins_lat_away"], 2.0),
             "skins_lat_away_exposed_at_1s": at_time(x["t"], x["skins_lat_away_exposed"], 1.0),
             "skins_lat_away_exposed_at_2s": at_time(x["t"], x["skins_lat_away_exposed"], 2.0),
             "skins_lat_away_back_at_1s": at_time(x["t"], x["skins_lat_away_back"], 1.0),
             "skins_lat_away_back_at_2s": at_time(x["t"], x["skins_lat_away_back"], 2.0),
             "skins_lat_heated_at_1s": at_time(x["t"], x["skins_lat_heated"], 1.0),
             "skins_lat_heated_at_2s": at_time(x["t"], x["skins_lat_heated"], 2.0),
             "skins_created": r["skins_created"], "n_active_final": S["n_active"][-1], "n_elements": r["n_elements"],
             "run_time_s": r["run_time_s"], "skin_wall_s": r["skin_wall_s"],
             "skin_drawn_mass_kg": r["skin_drawn_mass_kg"], "skin_drawn_through_mass_kg": r["skin_drawn_through_mass_kg"],
             "skin_handover_mass_kg": r["skin_handover_mass_kg"], "exposed_faces_total": r["exposed_faces_total"],
             "T_top_lat_at_times": o["T_top_at_times"], "T_plateau_std_K": o.get("T_plateau_std_K"),
             "T_top_max_K": o.get("T_top_max_K")}
        # the base runs against box_ablation's own runs of the same code and settings
        orig = os.path.join(ORIG_RUNS, sp[name]["orig"] + ".json") if name in sp else None
        if r["config"] == "base" and orig and os.path.exists(orig):
            with open(orig) as f:
                q = json.load(f)["series"]
            e["reproduces_box_ablation"] = {
                "file": os.path.basename(orig),
                "max_abs_diff_r_lat_mm": float(np.max(np.abs(np.array(S["r_lat_mm"]) - np.array(q["r_lat_mm"])))),
                "max_abs_diff_r_total_mm": float(np.max(np.abs(np.array(S["r_total_mm"]) - np.array(q["r_total_mm"])))),
                "max_abs_diff_T_top_lat_K": float(np.nanmax(np.abs(np.array(S["T_top_lat"]) - np.array(q["T_top_lat"])))),
                "bit_identical_r_lat": S["r_lat_mm"] == q["r_lat_mm"], "bit_identical_n_skins": S["n_skins"] == q["n_skins"]}
        out["runs"][name] = e
    with open(os.path.join(HERE, "results.json"), "w") as f:
        json.dump(out, f, indent=1)
    return out, recs, fine


def pct(x):
    return "n/a" if x is None else "{:+.1f} %".format(100 * x)


def verdict(e):
    rel = [(t, v) for t, v in zip(bx.COMPARE_TIMES, e["rel_diff_at_times"]) if v is not None]
    worst_t, worst = max(rel, key=lambda tv: abs(tv[1]))
    good = [t for t, v in rel if abs(v) <= 0.05]
    if len(good) == len(rel):
        return "Meets the 5 % criterion: its largest difference is {} at {:g} s.".format(pct(worst), worst_t)
    return ("Fails the 5 % criterion: its largest difference is {} at {:g} s, and it is within 5 % only at {}.".format(
        pct(worst), worst_t, ", ".join("{:g} s".format(t) for t in good) if good else "none of the compared times"))


def write_summary(out):
    L = ["# Experiment X1: the two switches in the box (variant q6p, dev preset, 10 ms sub-step)", "",
         "Generated by `x1_box.py --collect` from `runs/`. Lateral recession (`r_lat_mm`, the liquid removed from the "
         "lateral interior over its projected area and 2813 kg/m3) against the fine box's `r_total_mm` "
         "(analysis/box/runs/q6p__fine.json), which is {} mm at 0.5 / 1.0 / 1.5 / 2.0 s. Relative differences are (coarse - fine) / fine "
         "where the fine box has at least 0.05 mm, so not at 0.5 s.".format(" / ".join("{:.4f}".format(v) for v in out["fine_recession_mm_at_times"])), ""]
    L += ["| run | recession at 0.5 / 1.0 / 1.5 / 2.0 s (mm) | rel. diff at 1.0 / 1.5 / 2.0 s | largest abs rel. diff (4 times) | "
          "largest abs diff over the curve (mm) | stall: no lateral feed (s, from-to) | slow: below 1.1 mm/s (s, from-to) | worst energy residual | mass books (rel.) | repeats per step | "
          "max shortfall (kg) | skins lateral, not facing heat (n_x >= -0.1), at 1 / 2 s [on exposed faces; on the back face x = 8 mm] | skins created | active elements at end | run time (s) |",
          "|" + "---|" * 15]
    order = [n for c in CONFIGS for n in ["{}__dt{:g}".format(c, dt) for dt in MACRO] + ["{}__h1mm".format(c), "{}__h0.5mm".format(c)]
             if n in out["runs"]]
    for n in order:
        e = out["runs"][n]
        if e["status"] != "ok":
            L.append("| {} | FAILED: {} |".format(n, e.get("error")) + " |" * 13)
            continue
        st = "{:.2f} ({:.2f}-{:.2f})".format(e["stall_s"], e["stall_from_s"], e["stall_to_s"]) if e["stall_from_s"] is not None else (
            "none" if e["stall_s"] == 0.0 else "no feed")
        sl = "{:.2f} ({:.2f}-{:.2f})".format(e["slow_s"], e["slow_from_s"], e["slow_to_s"]) if e["slow_from_s"] is not None else (
            "none" if e["slow_s"] == 0.0 else "no feed")
        L.append("| {} | {} | {} | {} | {:.4f} | {} | {} | {:.1e} | {:.1e} | {:.2f} | {:.1e} | {} / {} [{} / {}; {} / {}] | {} | {} / {} | {:.0f} |".format(
            n, " / ".join("{:.4f}".format(v) for v in e["r_lat_mm_at_times"]),
            " / ".join(pct(v) for v in e["rel_diff_at_times"][1:]),
            "{:.1f} %".format(100 * e["max_abs_rel_diff_4times"]), e["max_abs_diff_curve_mm"], st, sl, e["worst_energy_residual"],
            e["mass_balance_rel"], e["repeats_per_step"], e["max_shortfall_kg"], e["skins_lat_away_at_1s"], e["skins_lat_away_at_2s"],
            e["skins_lat_away_exposed_at_1s"], e["skins_lat_away_exposed_at_2s"], e["skins_lat_away_back_at_1s"], e["skins_lat_away_back_at_2s"],
            e["skins_created"], e["n_active_final"], e["n_elements"], e["run_time_s"]))
    L += ["", "## Verdict per configuration against the 5 % criterion (the four macro steps, 20 um cells, 2 mm mesh)", ""]
    for c in CONFIGS:
        es = [(dt, out["runs"].get("{}__dt{:g}".format(c, dt))) for dt in MACRO]
        es = [(dt, e) for dt, e in es if e is not None and e["status"] == "ok"]
        if not es:
            continue
        met = [dt for dt, e in es if e["meets_5pct"]]
        L.append("- `{}`: {} of the four macro steps meet the criterion{}; its largest differences at the compared times are {}.".format(
            c, "none" if not met else len(met), "" if not met else " ({})".format(", ".join("{:g} s".format(d) for d in met)),
            ", ".join("{:.1f} % at {:g} s".format(100 * e["max_abs_rel_diff_4times"], dt) for dt, e in es)))
    L += ["", "## Verdicts per run against the 5 % criterion", ""]
    for n in order:
        e = out["runs"][n]
        L.append("- `{}`: {}".format(n, verdict(e) if e["status"] == "ok" else "The run failed: {}.".format(e.get("error"))))
    L += ["", "## Do the base runs reproduce box_ablation's?", ""]
    for n in order:
        e = out["runs"][n]
        rep = e.get("reproduces_box_ablation")
        if rep:
            L.append("- `{}` against `{}`: the largest difference in the lateral recession series is {:.3g} mm and in the total "
                     "recession {:.3g} mm; the lateral recession series is {}bit for bit the same, and the skin count series {}.".format(
                         n, rep["file"], rep["max_abs_diff_r_lat_mm"], rep["max_abs_diff_r_total_mm"],
                         "" if rep["bit_identical_r_lat"] else "not ", "is the same" if rep["bit_identical_n_skins"] else "differs"))
    with open(os.path.join(HERE, "summary.md"), "w") as f:
        f.write("\n".join(L) + "\n")


def plot(recs, fine):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception:
        return
    panels = ["dt{:g}".format(dt) for dt in MACRO] + ["h1mm", "h0.5mm"]
    fig, ax = plt.subplots(2, 3, figsize=(15, 8.5), squeeze=False)
    ft = [0.0] + fine["series"]["t"]
    fr = [0.0] + fine["series"]["r_total_mm"]
    for i, pnl in enumerate(panels):
        a = ax[i // 3][i % 3]
        a.plot(ft, fr, "k-", lw=2.5, label="fine box")
        for c in CONFIGS:
            r = recs.get("{}__{}".format(c, pnl))
            if r is None or r.get("status") != "ok":
                continue
            a.plot([0.0] + r["series"]["t"], [0.0] + r["series"]["r_lat_mm"], lw=1.3, label=c, color=COLOURS[c])
        a.set_title("q6p, dev, " + pnl.replace("dt", "macro step ").replace("h", "h = ", 1) + (" s" if pnl.startswith("dt") else ", 0.05 s"))
        a.set_xlabel("time (s)")
        a.set_ylabel("lateral recession (mm)")
        a.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "x1_box.png"), dpi=120)


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("names", nargs="*")
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--collect", action="store_true")
    a = ap.parse_args()
    sp = specs()
    if a.list:
        print("\n".join(sp))
        return
    if a.collect:
        out, recs, fine = collect()
        write_summary(out)
        plot(recs, fine)
        return
    for name in a.names:
        if name not in sp:
            raise SystemExit("unknown run " + name)
        run(name, sp[name])


if __name__ == "__main__":
    main()

"""Windward surface loads on a sphere: Holman & Boyd (2011) DSMC against Daniels & Schwartzentruber (2026) and the
re-entry model's current equations, all at Holman & Boyd's conditions (0.3048 m sphere, Mach 25 air, 7108 m/s,
200 K freestream, 1500 K wall; freestream Knudsen numbers 0.002, 0.01, 0.05, 0.25).

Data (research_notes/Film separation at a convex shoulder/data/):
  hb_dsmc.csv      DSMC curves digitized from Holman & Boyd Figs 5-6, 8-9, 11-12, 14-15
  ds_model.csv     Daniels & Schwartzentruber Appendix A, perfect gas, fully diffuse wall (no heating model)
  model_loads.csv  reentry_model.shared SurfaceFlow + PhysicsHeating with a default flight's options
Coefficients: Cp = (p - p_inf)/(0.5 rho U^2), C_tau = tau/(0.5 rho U^2), Cq = q/(0.5 rho U^3).
"""
import csv
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
DATA = HERE.parent / "reference"
CASES = (0.002, 0.01, 0.05, 0.25)
SOURCES = (  # (file, label, colour, line style, width)
    ("hb_dsmc.csv", "Holman & Boyd 2011, DSMC (digitized)", "#2a78d6", "-", 2.2),
    ("ds_model.csv", "Daniels & Schwartzentruber 2026 model", "#eb6834", "--", 1.8),
    ("model_loads.csv", "Re-entry model, current equations", "#1baf7a", "-.", 1.8),
)
# Candidate fix, shear row only: local-similarity wall shear (Lees 1956, cold wall, variable rho*mu) on the model's own
# edge state, passed through the model's slip and bridging (data/similarity_shear.py). Same colour as the re-entry
# model (it is that model with tau_c replaced); solid line with markers to tell the variants apart.
SIMILARITY = (("similarity_shear.csv", "Re-entry model with similarity-based shear", "#1baf7a", "-", 2.0, "o"),
              ("similarity_hb_pressure.csv", "Similarity shear on Holman & Boyd's pressure, equilibrium gas",
               "#1baf7a", ":", 2.0, "s"),
              ("similarity_hb_pressure_frozen.csv", "Similarity shear on Holman & Boyd's pressure, unreacted gas",
               "#1baf7a", "--", 2.0, "^"))
ROWS = (("Cp", "Pressure coefficient", "log", (1e-3, 3)),
        ("Ctau", "Shear coefficient", "linear", None),
        ("Cq", "Heat-flux coefficient", "log", (3e-4, 1.5)))


def load(name):
    out = defaultdict(lambda: defaultdict(list))
    with open(DATA / name) as f:
        for r in csv.DictReader(f):
            th = float(r["theta_deg"])
            if th > 90:
                continue
            for key in ("Cp", "Ctau", "Cq"):
                if r[key] not in ("", None):
                    out[float(r["case_kn"])][key].append((th, float(r[key])))
    return out


def main():
    data = {name: load(name) for name, *_ in SOURCES}
    sim = {name: defaultdict(list) for name, *_ in SIMILARITY}
    for name, *_ in SIMILARITY:
        with open(DATA / name) as f:
            for r in csv.DictReader(f):
                sim[name][float(r["case_kn"])].append((float(r["theta_deg"]), float(r["Ctau"])))
    plt.rcParams.update({"font.size": 9, "axes.edgecolor": "#8a8984", "axes.labelcolor": "#0b0b0b",
                         "xtick.color": "#52514e", "ytick.color": "#52514e"})
    fig, axes = plt.subplots(3, 4, figsize=(13.5, 9.2), sharex=True)
    fig.patch.set_facecolor("#fcfcfb")
    for i, (key, ylabel, scale, ylim) in enumerate(ROWS):
        for j, kn in enumerate(CASES):
            ax = axes[i, j]
            ax.set_facecolor("#fcfcfb")
            for name, label, colour, ls, lw in SOURCES:
                pts = data[name][kn].get(key, [])
                if scale == "log":
                    pts = [(t, v) for t, v in pts if v > 0]
                if pts:
                    ax.plot(*zip(*pts), color=colour, ls=ls, lw=lw, label=label)
            if key == "Ctau":
                for k, (name, label, colour, ls, lw, mk) in enumerate(SIMILARITY):
                    ax.plot(*zip(*sim[name][kn]), color=colour, ls=ls, lw=lw, marker=mk, ms=4,
                            markevery=((3 * k) % 10, 10), label=label)
            ax.set_yscale(scale)
            if ylim:
                ax.set_ylim(*ylim)
            else:
                ax.set_ylim(bottom=0)
            ax.set_xlim(0, 90)
            ax.set_xticks([0, 15, 30, 45, 60, 75, 90])
            ax.grid(True, which="major", color="#e4e3df", lw=0.6)
            for s in ("top", "right"):
                ax.spines[s].set_visible(False)
            if i == 0:
                ax.set_title(f"Freestream Knudsen number {kn:g}", fontsize=10, color="#0b0b0b")
            if j == 0:
                ax.set_ylabel(ylabel)
            if i == 2:
                ax.set_xlabel("Angle from stagnation point (deg)")
    axes[2, 0].text(3, 4.5e-4, "Daniels & Schwartzentruber\ngive no heating model", fontsize=8, color="#52514e")
    handles, labels = axes[1, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=2, frameon=False, bbox_to_anchor=(0.5, 0.97))
    fig.suptitle("Windward surface loads, 0.305 m sphere in Mach 25 air at 7.1 km/s with a 1500 K wall "
                 "(Holman & Boyd 2011 conditions)", fontsize=11.5, y=0.99)
    fig.tight_layout(rect=(0, 0, 1, 0.89))
    out = HERE.parents[2] / "output" / "gas_loads" / "windward_loads_comparison.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=150, facecolor=fig.get_facecolor())
    print("saved", out)
    # values quoted in the write-up
    for kn in CASES:
        for key in ("Cp", "Ctau", "Cq"):
            row = []
            for name, label, *_ in SOURCES:
                d = dict(data[name][kn].get(key, []))
                row.append(" ".join(f"{d[t]:.3g}" if t in d else "-" for t in (0, 30, 45, 60, 75, 90)))
            print(f"Kn {kn:<5} {key:<4} HB[{row[0]}]  DS[{row[1]}]  MODEL[{row[2]}]")


if __name__ == "__main__":
    main()

"""Leeward surface loads on a sphere (90-180 deg): Holman & Boyd (2011) DSMC against Daniels & Schwartzentruber (2026),
at Holman & Boyd's conditions (0.3048 m sphere, Mach 25 air at 7108 m/s, 200 K freestream, 1500 K wall; freestream
Knudsen numbers 0.002, 0.01, 0.05, 0.25). Data in research_notes/Film separation at a convex shoulder/data/.

Pressure is plotted as p / p_inf so both sources and the model's current lee value (p = p_inf) share one scale.
Holman & Boyd's plotted pressure coefficient is read as p / (0.5 rho U^2): their rear values fall to 3e-5 (Kn 0.05)
and 3e-7 (Kn 0.25), far below p_inf / (0.5 rho U^2) = 0.0023, which only p/q can do (an almost empty wake at speed
ratio 21). Daniels & Schwartzentruber's Cp is (p - p_inf)/(0.5 rho U^2) as their paper defines it.
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
Q_OVER_PINF = 0.5 * 1.4 * (7108.0 / (1.4 * 287.05 * 200.0) ** 0.5) ** 2      # 0.5 gamma M^2 = 440
SOURCES = (
    ("hb_dsmc.csv", "Holman & Boyd 2011, DSMC (digitized)", "#2a78d6", "-", 2.2, lambda cp: cp * Q_OVER_PINF),
    ("ds_model.csv", "Daniels & Schwartzentruber 2026 model", "#eb6834", "--", 1.8, lambda cp: 1.0 + cp * Q_OVER_PINF),
)
ROWS = (("p", "Pressure / freestream pressure", "log", (1e-4, 100)),
        ("Ctau", "Shear coefficient", "linear", None),
        ("Cq", "Heat-flux coefficient", "log", (1e-4, 0.1)))


def load(name, p_of_cp):
    out = defaultdict(lambda: defaultdict(list))
    with open(DATA / name) as f:
        for r in csv.DictReader(f):
            th = float(r["theta_deg"])
            if th < 90:
                continue
            kn = float(r["case_kn"])
            if r["Cp"]:
                out[kn]["p"].append((th, p_of_cp(float(r["Cp"]))))
            for key in ("Ctau", "Cq"):
                if r[key]:
                    out[kn][key].append((th, float(r[key])))
    return out


def main():
    data = {name: load(name, f) for name, _, _, _, _, f in SOURCES}
    plt.rcParams.update({"font.size": 9, "axes.edgecolor": "#8a8984", "xtick.color": "#52514e", "ytick.color": "#52514e"})
    fig, axes = plt.subplots(3, 4, figsize=(13.5, 9.2), sharex=True)
    fig.patch.set_facecolor("#fcfcfb")
    for i, (key, ylabel, scale, ylim) in enumerate(ROWS):
        for j, kn in enumerate(CASES):
            ax = axes[i, j]
            ax.set_facecolor("#fcfcfb")
            if key == "p":
                ax.axhline(1.0, color="#8a8984", lw=1.2, ls=(0, (1, 2)),
                           label="Freestream pressure (the model's current lee value)")
            for name, label, colour, ls, lw, _ in SOURCES:
                pts = data[name][kn].get(key, [])
                if scale == "log":
                    pts = [(t, v) for t, v in pts if v > 0]
                if pts:
                    ax.plot(*zip(*pts), color=colour, ls=ls, lw=lw, label=label)
            ax.set_yscale(scale)
            if ylim:
                ax.set_ylim(*ylim)
            else:
                ax.set_ylim(bottom=0)
            ax.set_xlim(90, 180)
            ax.set_xticks([90, 105, 120, 135, 150, 165, 180])
            ax.grid(True, which="major", color="#e4e3df", lw=0.6)
            for s in ("top", "right"):
                ax.spines[s].set_visible(False)
            if i == 0:
                ax.set_title(f"Freestream Knudsen number {kn:g}", fontsize=10)
            if j == 0:
                ax.set_ylabel(ylabel)
            if i == 2:
                ax.set_xlabel("Angle from stagnation point (deg)")
    axes[2, 0].text(92, 1.6e-4, "Daniels & Schwartzentruber give no heating model;\nthe re-entry model sets lee "
                    "shear and heating to zero", fontsize=8, color="#52514e")
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=3, frameon=False, bbox_to_anchor=(0.5, 0.955))
    fig.suptitle("Leeward surface loads, 0.305 m sphere in Mach 25 air at 7.1 km/s with a 1500 K wall "
                 "(Holman & Boyd 2011 conditions)", fontsize=11.5, y=0.99)
    fig.tight_layout(rect=(0, 0, 1, 0.925))
    out = HERE.parents[2] / "output" / "gas_loads" / "leeward_loads_comparison.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=150, facecolor=fig.get_facecolor())
    print("saved", out)
    for kn in CASES:
        for key in ("p", "Ctau", "Cq"):
            row = []
            for name, *_ in SOURCES:
                d = dict(data[name][kn].get(key, []))
                row.append(" ".join(f"{d[t]:.3g}" if t in d else "-" for t in (90, 100, 110, 120, 135, 150, 165, 175)))
            print(f"Kn {kn:<5} {key:<4} HB[{row[0]}]  DS[{row[1]}]")


if __name__ == "__main__":
    main()

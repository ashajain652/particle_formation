"""Heating round the whole sphere (0-180 deg) at Holman & Boyd's conditions: their DSMC, the re-entry model's current
heating (Fay-Riddell equilibrium, fully catalytic, Lees shape, Matting bridge; zero leeward), and the similarity-based
heating on the design's decided pressure with unreacted post-shock gas and the model's Matting bridge (computed until
the boundary layer separates on the rising lee pressure). Top row absolute heat-flux coefficient; bottom row the same
divided by its stagnation value, to separate the shape from the magnitude. Data: research_notes/Film separation at a
convex shoulder/data/ (hb_dsmc.csv, model_loads.csv, similarity_heat.csv from similarity_heat.py).
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


def read(name, key):
    out = defaultdict(list)
    with open(DATA / name) as f:
        for r in csv.DictReader(f):
            if r[key] and float(r[key]) > 0:
                out[float(r["case_kn"])].append((float(r["theta_deg"]), float(r[key])))
    return out


def main():
    sources = (
        (read("hb_dsmc.csv", "Cq"), "Holman & Boyd 2011, DSMC (digitized)", dict(color="#2a78d6", ls="-", lw=2.2)),
        (read("model_loads.csv", "Cq"), "Re-entry model, current heating (zero leeward)",
         dict(color="#1baf7a", ls="-.", lw=1.8)),
        (read("similarity_heat.csv", "Cq"), "Similarity-based heating, unreacted gas, model's Matting bridge (to separation)",
         dict(color="#1baf7a", ls="-", lw=2.0, marker="o", ms=4, markevery=15)),
    )
    plt.rcParams.update({"font.size": 9, "axes.edgecolor": "#8a8984", "xtick.color": "#52514e", "ytick.color": "#52514e"})
    fig, axes = plt.subplots(2, 4, figsize=(13.5, 7.0), sharex=True)
    fig.patch.set_facecolor("#fcfcfb")
    for j, kn in enumerate(CASES):
        for i in range(2):
            ax = axes[i, j]
            ax.set_facecolor("#fcfcfb")
            ax.axvline(90, color="#c3c2b7", lw=1.0)
            for data, label, style in sources:
                pts = data[kn]
                if not pts:
                    continue
                if i == 1:
                    q0 = dict(pts).get(0.0) or pts[0][1]
                    pts = [(t, v / q0) for t, v in pts]
                ax.plot(*zip(*pts), label=label, **style)
            ax.set_yscale("log")
            ax.set_ylim(*((1e-5, 2) if i == 0 else (1e-5, 2)))
            ax.set_xlim(0, 180); ax.set_xticks([0, 30, 60, 90, 120, 150, 180])
            ax.grid(True, color="#e4e3df", lw=0.6)
            for s in ("top", "right"):
                ax.spines[s].set_visible(False)
            if i == 0:
                ax.set_title(f"Freestream Knudsen number {kn:g}", fontsize=10)
            if j == 0:
                ax.set_ylabel("Heat-flux coefficient" if i == 0 else "Heat flux / stagnation value")
            if i == 1:
                ax.set_xlabel("Angle from stagnation point (deg)")
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=2, frameon=False, bbox_to_anchor=(0.5, 0.955))
    fig.suptitle("Heating round the sphere, 0.305 m sphere in Mach 25 air at 7.1 km/s with a 1500 K wall "
                 "(Holman & Boyd 2011 conditions)", fontsize=11.5, y=0.99)
    fig.tight_layout(rect=(0, 0, 1, 0.88))
    out = HERE.parents[2] / "output" / "gas_loads" / "heating_comparison.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=150, facecolor=fig.get_facecolor())
    print("saved", out)
    for kn in CASES:
        row = []
        for data, *_ in sources:
            d = dict(data[kn])
            row.append(" ".join(f"{d[t]:.3g}" if t in d else "-" for t in (0.0, 30.0, 60.0, 90.0, 100.0, 120.0, 140.0, 160.0)))
        print(f"Kn {kn:<5} HB[{row[0]}]  MODEL[{row[1]}]  SIM[{row[2]}]")


if __name__ == "__main__":
    main()

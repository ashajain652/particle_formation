"""Leeward shear options (90-180 deg) at Holman & Boyd's conditions: their DSMC, Daniels & Schwartzentruber, and the
similarity-based shear run round the body on the design's decided pressure (the model's Prandtl-Meyer expansion
continued past 90 deg until it meets Daniels & Schwartzentruber's lee pressure; unreacted post-shock gas; separated,
zero shear, once the rising lee pressure stalls the boundary layer). Data: research_notes/Film separation at a convex
shoulder/data/ (hb_dsmc.csv, ds_model.csv, similarity_lee.csv from similarity_lee.py).
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
Q_OVER_PINF = 0.5 * 1.4 * (7108.0 / (1.4 * 287.05 * 200.0) ** 0.5) ** 2      # 440


def read(name, pkey, pfun, tkey):
    out = defaultdict(lambda: defaultdict(list))
    with open(DATA / name) as f:
        for r in csv.DictReader(f):
            th = float(r["theta_deg"])
            if th < 90:
                continue
            kn = float(r["case_kn"])
            if r[pkey]:
                out[kn]["p"].append((th, pfun(float(r[pkey]))))
            if r[tkey]:
                out[kn]["tau"].append((th, float(r[tkey])))
    return out


def main():
    sources = (
        (read("hb_dsmc.csv", "Cp", lambda c: c * Q_OVER_PINF, "Ctau"), "Holman & Boyd 2011, DSMC (digitized)",
         dict(color="#2a78d6", ls="-", lw=2.2)),
        (read("ds_model.csv", "Cp", lambda c: 1 + c * Q_OVER_PINF, "Ctau"), "Daniels & Schwartzentruber 2026 model",
         dict(color="#eb6834", ls="--", lw=1.8)),
        (read("similarity_lee.csv", "p_over_pinf", lambda v: v, "Ctau"),
         "Similarity shear on the decided lee pressure (Prandtl-Meyer joined to Daniels), unreacted gas",
         dict(color="#1baf7a", ls="-", lw=2.0, marker="o", ms=4, markevery=10)),
    )
    plt.rcParams.update({"font.size": 9, "axes.edgecolor": "#8a8984", "xtick.color": "#52514e", "ytick.color": "#52514e"})
    fig, axes = plt.subplots(2, 4, figsize=(13.5, 6.8), sharex=True)
    fig.patch.set_facecolor("#fcfcfb")
    for j, kn in enumerate(CASES):
        for i, key in enumerate(("p", "tau")):
            ax = axes[i, j]
            ax.set_facecolor("#fcfcfb")
            if key == "p":
                ax.axhline(1.0, color="#8a8984", lw=1.2, ls=(0, (1, 2)), label="Freestream pressure")
            for data, label, style in sources:
                pts = [(t, v) for t, v in data[kn][key] if (v > 0 or key == "tau")]
                if pts:
                    ax.plot(*zip(*pts), label=label, **style)
            if key == "p":
                ax.set_yscale("log"); ax.set_ylim(1e-4, 100)
            else:
                ax.set_ylim(bottom=0)
            ax.set_xlim(90, 180); ax.set_xticks([90, 105, 120, 135, 150, 165, 180])
            ax.grid(True, color="#e4e3df", lw=0.6)
            for s in ("top", "right"):
                ax.spines[s].set_visible(False)
            if i == 0:
                ax.set_title(f"Freestream Knudsen number {kn:g}" + ("\n" if kn < 0.05 else "\n(design uses Daniels here)"),
                             fontsize=10)
            if j == 0:
                ax.set_ylabel("Pressure / freestream pressure" if key == "p" else "Shear coefficient")
            if i == 1:
                ax.set_xlabel("Angle from stagnation point (deg)")
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=2, frameon=False, bbox_to_anchor=(0.5, 0.955))
    fig.suptitle("Leeward shear options, 0.305 m sphere in Mach 25 air at 7.1 km/s with a 1500 K wall "
                 "(Holman & Boyd 2011 conditions)", fontsize=11.5, y=0.99)
    fig.tight_layout(rect=(0, 0, 1, 0.87))
    out = HERE.parents[2] / "output" / "gas_loads" / "leeward_shear_options.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=150, facecolor=fig.get_facecolor())
    print("saved", out)
    for kn in CASES:
        row = []
        for data, label, _ in sources:
            d = dict(data[kn]["tau"])
            row.append(" ".join(f"{d[t]:.4f}" if t in d else "-" for t in (90, 100, 110, 120, 135, 150)))
        print(f"Kn {kn:<5} shear  HB[{row[0]}]  DS[{row[1]}]  SIM[{row[2]}]")


if __name__ == "__main__":
    main()

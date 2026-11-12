#!/usr/bin/env python3
"""The composite near-liquidus viscosity law of AA7075 against the present (bridge) law: the figure and the table of the
near-liquidus viscosity design (docs/superpowers/specs/2026-11-12-near-liquidus-viscosity-design.md, sections 3.2-3.7).

    "$PY" analysis/near_liquidus_viscosity.py [--outdir output/near_liquidus_viscosity]

Writes `law.png` (viscosity over the liquid's against solid fraction, 0 to 0.5: the law for internal solid fractions
0.5, 1 and 0.25; Li et al.'s five per-rate fits at their measured solid fractions; the present law at 0 and 367 1/s) and
`law_tables.txt` (spec 3.6's table). Every curve comes from `reentry_model.shared.rheology`; nothing of the law is copied
here, only Li et al.'s Table 3 data for reference. The law's assumptions and uncertainties are those of spec 3.7 and of
the module docstring of `reentry_model/shared/rheology.py`: nothing below 10 % solid has been measured in aluminium, the
internal solid fraction is assumed, and the jump at 10 % solid is unresolved."""
import argparse
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, REPO_ROOT)
from reentry_model import OUTPUT_DIR  # noqa: E402
from reentry_model.shared import material as M, rheology as R  # noqa: E402

MU_L = 1.3e-3                                   # Pa s, AA7075's liquid
LI_RATES = np.array([61.24, 122.47, 244.94, 367.41, 489.88])        # 1/s, Li et al. 2014 Table 3
LI_A = np.array([0.66711, 0.59623, 0.33326, 0.14323, 0.04346])      # Pa s
LI_B = np.array([4.78, 3.52, 3.45, 3.66, 3.23])
LI_F = np.array([0.1, 0.2, 0.3, 0.4, 0.5])                          # their measured solid fractions
INK, INK2, GRID, SURFACE, REF = "#0b0b0b", "#52514e", "#e6e5e1", "#fcfcfb", "#9a9893"
BLUES = ["#86b6ef", "#5598e7", "#2a78d6", "#1c5cab", "#104281"]
TABLE_F = [0.01, 0.025, 0.038, 0.057, 0.0975, 0.1 - 1e-9, 0.1, 0.2, 0.3, 0.4, 0.45, 0.4999]
TABLE_FI = [(0.5, "f_i = 0.5 (the law)"), (1.0, "f_i = 1 (bracket)"), (0.25, "f_i = 0.25 (bracket)")]


def tables():
    head = "viscosity / liquid's viscosity (1.3 mPa s), from rheology.composite_viscosity (spec 3.6)\n"
    head += "{:<22s}".format("solid fraction") + "".join("{:>9s}".format("{:.4g}".format(f) if f < 0.0999999 or f >= 0.1 else "<0.1")
                                                        for f in TABLE_F) + "\n"
    rows = [head]
    for fi, name in TABLE_FI:
        rows.append("{:<22s}".format(name) + "".join("{:9.3f}".format(R.composite_viscosity(f, MU_L, fi) / MU_L)
                                                     for f in TABLE_F) + "\n")
    rows.append("\nThe jump at 10 % solid (unresolved, spec 3.5): ")
    rows.append(", ".join("{} to {:.1f} times the liquid (a factor of {:.1f})".format(
        "{:.2f}".format(R.composite_viscosity(0.1 - 1e-12, MU_L, fi) / MU_L), R.composite_viscosity(0.1, MU_L, fi) / MU_L,
        R.composite_viscosity(0.1, MU_L, fi) / R.composite_viscosity(0.1 - 1e-12, MU_L, fi)) for fi, _ in TABLE_FI) + "\n")
    rows.append("\nAssumptions and uncertainties (spec 3.7): nothing below 10 % solid has been measured in aluminium; the "
                "internal solid fraction is assumed (0.5), bracketed by 1 and 0.25; the structure in flight is assumed to be "
                "Li et al.'s most broken-down one and the viscosity independent of rate; the jump at 10 % solid is "
                "unresolved. See the module docstring of reentry_model/shared/rheology.py.\n")
    return "".join(rows)


def figure(path):
    mat = M.Material.from_drama_json(M.MATERIAL_NAMES["AA7075_scheil"])
    f_b = R.bridge_solid_fraction(mat)
    fs = np.linspace(0.0, 0.499, 1500)
    T = np.linspace(895.0, 910.0, 6001)
    f_T = 1.0 - np.asarray(mat.liquid_fraction(T), dtype=float)
    fig, ax = plt.subplots(figsize=(9.0, 6.0), facecolor=SURFACE)
    ax.set_facecolor(SURFACE)
    ax.axvspan(0.0, 0.1, color=GRID, alpha=0.7, lw=0)
    ax.text(0.05, 6e3, "assumed:\nno aluminium data\n(spec 3.7)", ha="center", va="center", fontsize=8.5, color=INK2)
    for fi, c, ls, lw in ((0.5, BLUES[3], "-", 2.2), (1.0, BLUES[1], (0, (4, 2)), 1.5), (0.25, BLUES[4], (0, (1, 1.5)), 1.5)):
        eta = R.composite_viscosity(fs, MU_L, fi) / MU_L
        lo = fs < R.LI_SOLID_MIN
        ax.plot(fs[lo], eta[lo], color=c, ls=ls, lw=lw, label="composite law, f_i = {:g}{}".format(fi, " (default)" if fi == 0.5 else ""))
        ax.plot(fs[~lo], eta[~lo], color=c, ls=ls, lw=lw)
    ax.plot([0.1, 0.1], [R.composite_viscosity(0.1 - 1e-12, MU_L) / MU_L, R.composite_viscosity(0.1, MU_L) / MU_L],
            color=INK, lw=1.0, ls=":")
    ax.annotate("jump unresolved", xy=(0.1, 9.0), xytext=(0.115, 3.0), fontsize=8.5, color=INK, arrowprops=dict(arrowstyle="-", color=INK2))
    for k, (g, a, b) in enumerate(zip(LI_RATES, LI_A, LI_B)):
        ax.plot(LI_F, a * np.exp(b * LI_F) / MU_L, "o", ms=5, color=BLUES[k], mec=SURFACE, mew=1.2,
                label="Li et al. 2014 fit at {:.0f} 1/s, at measured points".format(g))
    for g, c, ls in ((0.0, "#c4501f", "-"), (367.0, "#8f3a15", (0, (4, 2)))):
        eta = R.law_viscosity("bridge", T, f_T, g, MU_L, f_b) / MU_L
        ok = np.isfinite(eta) & (f_T < 0.5)
        o = np.argsort(f_T[ok])
        ax.plot(f_T[ok][o], eta[ok][o], color=c, ls=ls, lw=1.4, label="present (bridge) law at {:.0f} 1/s".format(g))
    ax.set_yscale("log")
    ax.set_xlim(0.0, 0.5)
    ax.set_ylim(0.9, 1.5e4)
    ax.grid(True, which="major", color=GRID, lw=0.8)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.set_xlabel("solid fraction", color=INK, fontsize=10)
    ax.set_ylabel("viscosity / liquid's viscosity (1.3 mPa s)", color=INK, fontsize=10)
    ax.set_title("AA7075 near-liquidus viscosity: composite law against the present law", loc="left", color=INK, fontsize=11)
    ax.legend(loc="lower right", fontsize=8, frameon=False)
    fig.tight_layout()
    fig.savefig(path, dpi=160, facecolor=SURFACE)
    plt.close(fig)


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--outdir", default=os.path.join(OUTPUT_DIR, "near_liquidus_viscosity"))
    a = p.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    text = tables()
    with open(os.path.join(a.outdir, "law_tables.txt"), "w") as fh:
        fh.write(text)
    figure(os.path.join(a.outdir, "law.png"))
    print(text)
    print("wrote", a.outdir)


if __name__ == "__main__":
    main()

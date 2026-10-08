"""Girin's published spraying cases against the model's transcription of his theory (spec Step 3 sections 13.2-13.3).

    "$PY" analysis/girin_reference.py [--outdir reentry_model_output/verification_melt/girin] [--re-density ambient|shock|both]

Runs the three variants of Girin (2017) Table 1 in Girin-as-published mode (reentry_model.girin_case) and checks the
Girin & Kopyt (1994) Tables 1-2 with the model's thin-film and Rayleigh-Taylor modes (reentry_model.spray). Writes
girin2017.json / girin1994.json (every number next to its published value and ratio), the size distributions dn(r),
dM(r) per variant with the table values marked, the mass-loss law, and log-log plots of the 1994 tables with
residuals, plus girin_summary.md. Thresholds (README, Step 3): exact tier GI and phi_cr within 2 %; integrated tier
t_f, N and r_med within 30 % with the ambient-density Reynolds number; t_s.d. reported."""
import argparse
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, REPO_ROOT)
from reentry_model import girin_case, spray  # noqa: E402
from reentry_model.compare import INK, MODEL_COLOR, MUTED, SECOND, apply_rcparams, strip_top_right_spines  # noqa: E402

DEFAULT_OUTDIR = os.path.join(REPO_ROOT, "reentry_model_output", "verification_melt", "girin")
TABLE_2017 = os.path.join(REPO_ROOT, "data", "reference_values", "girin2017_table1.json")
TABLE_1994 = os.path.join(REPO_ROOT, "data", "reference_values", "girin1994_tables.json")


def ratio(model, ref):
    return None if (model is None or ref in (None, 0)) else float(model / ref)


def run_2017(outdir, densities):
    tab = json.load(open(TABLE_2017))
    results = {}
    for name, d in tab["variants"].items():
        v = girin_case.GirinVariant(name, d["V_ms"], d["rho_m"], d["mu_m"], d["sigma"])
        results[name] = {}
        for rd in densities:
            o = girin_case.run_variant(v, re_density=rd)
            entry = {"GI": (o["GI"], d["GI"]), "phi_cr_deg": (o["phi_cr_deg_eq3"], d["phi_cr_deg"]), "t_f_us": (o["t_f_min_us_tau0"], d["t_f_us"]),
                     "N": (o["N"], d["N"]), "r_med_um": (o["r_med_um"], d["r_med_um"]), "r_min_um": (o["r_min_um"], d["r_min_um"]),
                     "r_max_um": (o["r_max_um"], d["r_max_um"]), "t_sd_ms": (o["t_sd_ms"], d["t_sd_ms"]), "tau_end": (o["tau_end"], d["tau_end"]),
                     "t_ch_ms": (o["t_ch_s"] * 1e3, d["t_ch_ms"]), "Re_inf": (o["Re_inf"], d["Re_inf"])}
            results[name][rd] = {k: {"model": m, "published": p, "ratio": ratio(m, p)} for k, (m, p) in entry.items()}
            results[name][rd]["r_um_90deg_tau0"] = o["r_um_90deg_tau0"]
            # plots: size distributions and the mass-loss law
            apply_rcparams(plt)
            fig, (ax, bx, cx) = plt.subplots(1, 3, figsize=(14, 4.2))
            if o["radii"].size:
                n_hist, m_hist = spray.histogram(o["radii"], o["counts"], o["counts"] * 4.0 / 3.0 * np.pi * v.rho_m * o["radii"] ** 3)
                mid = np.sqrt(spray.BIN_EDGES[:-1] * spray.BIN_EDGES[1:]) * 1e6
                ax.loglog(mid[n_hist > 0], n_hist[n_hist > 0], color=MODEL_COLOR, lw=1.4, label="model")
                bx.loglog(mid[m_hist > 0], m_hist[m_hist > 0] * 1e3, color=MODEL_COLOR, lw=1.4, label="model")
                for value, label in ((d["r_med_um"], "published r_med"), (d["r_min_um"], "published range"), (d["r_max_um"], None)):
                    if value:
                        ax.axvline(value, color=INK, lw=0.8, ls=":", label=label)
                        bx.axvline(value, color=INK, lw=0.8, ls=":")
                mh = o["mass_history"]
                cx.plot(mh[:, 0] * 1e3, mh[:, 1] / mh[0, 1], color=MODEL_COLOR, lw=1.4, label="model")
                t_sd = mh[-1, 0]
                cx.plot(mh[:, 0] * 1e3, (1.0 - np.minimum(mh[:, 0] / t_sd, 1.0)) ** 3, color=MUTED, lw=1.0, ls="--", label="(1 - t/t_s.d.)^3")
                cx.axvline(d["t_sd_ms"], color=INK, lw=0.8, ls=":", label="published t_s.d.")
            ax.set_xlabel("r [um]"); ax.set_ylabel("dn per bin"); ax.legend(frameon=False)
            bx.set_xlabel("r [um]"); bx.set_ylabel("dM per bin [g]")
            cx.set_xlabel("t [ms]"); cx.set_ylabel("m / m0"); cx.legend(frameon=False)
            ax.set_title("Girin 2017 variant {} ({}), Re on {} density".format(name, d["description"], rd), color=SECOND, fontsize=10)
            for a in (ax, bx, cx):
                strip_top_right_spines(a)
            fig.tight_layout(); fig.savefig(os.path.join(outdir, "girin2017_variant_{}_{}.png".format(name, rd)), dpi=150); plt.close(fig)
    return results


def run_1994(outdir):
    tab = json.load(open(TABLE_1994))
    t1, t2 = tab["table1"], tab["table2"]

    class Liquid:
        rho, sigma = t1["rho1_kgm3"], t1["sigma_Nm"]
        mu = 0.0
    rho2, V0 = np.array(t1["rho2_kgm3"]), np.array(t1["V0_ms"])
    r_tab, tau_tab, m_tab = (np.array(t1[k]) for k in ("r_d_m", "tau_d_s", "m_kgm2s"))
    lam, tau = spray.thin_film_mode(np.full((3, 2), t1["mach"]), rho2[:, None] * V0[None, :] ** 2, Liquid)
    f_s = float((lam / 4.0)[0, 0] / r_tab[0, 0])                     # their shock-layer factor, fitted on the first entry
    r_model, tau_model = lam / 4.0 / f_s, tau / f_s ** 1.5
    m_model = Liquid.rho * r_model / (2.0 * tau_model)
    W = np.array(t2["W_ms2"])
    lam_rt, tau_rt = zip(*[spray.rayleigh_taylor(w, np.array([1.0]), Liquid)[1:] for w in W])
    lam_rt, tau_rt = np.array([x[0] for x in lam_rt]), np.array([x[0] for x in tau_rt])
    res = {"table1": {"shock_layer_factor": f_s, "r_d_ratio": (r_model / r_tab).tolist(), "tau_d_ratio": (tau_model / tau_tab).tolist(),
                      "m_ratio": (m_model / m_tab).tolist(), "r_d_model_m": r_model.tolist(), "tau_d_model_s": tau_model.tolist()},
           "table2": {"lambda_star_model_m": lam_rt.tolist(), "lambda_star_ratio_to_eq14_column": (lam_rt / np.array(t2["lambda_star_m_eq14"])).tolist(),
                      "lambda_star_ratio_to_printed": (lam_rt / np.array(t2["lambda_star_m_as_printed"])).tolist(),
                      "tau_star_model_s": tau_rt.tolist(), "tau_star_ratio": (tau_rt / np.array(t2["tau_star_s"])).tolist()}}
    apply_rcparams(plt)
    fig, axes = plt.subplots(2, 3, figsize=(14, 7), gridspec_kw={"height_ratios": [3, 1]})
    for j, (name, model, published, unit) in enumerate((("r_d", r_model, r_tab, "m"), ("tau_d", tau_model, tau_tab, "s"), ("m", m_model, m_tab, "kg/m2/s"))):
        ax, rx = axes[0, j], axes[1, j]
        for k, V in enumerate(V0):
            ax.loglog(rho2, published[:, k], "o", color=INK, label="published, V0 = {:.0f} km/s".format(V / 1e3))
            ax.loglog(rho2, model[:, k], "--", color=MODEL_COLOR, label="model x f_s" if k == 0 else None)
            rx.semilogx(rho2, 100.0 * (model[:, k] / published[:, k] - 1.0), "s-", color=MODEL_COLOR if k == 0 else MUTED, lw=0.8)
        ax.set_ylabel("{} [{}]".format(name, unit)); ax.legend(frameon=False, fontsize=8); ax.set_title("Girin & Kopyt 1994 Table 1", color=SECOND, fontsize=10)
        rx.axhline(0.0, color=MUTED, lw=0.6); rx.set_xlabel("rho_2 [kg/m3]"); rx.set_ylabel("model/published - 1 [%]")
        strip_top_right_spines(ax); strip_top_right_spines(rx)
    fig.tight_layout(); fig.savefig(os.path.join(outdir, "girin1994_table1.png"), dpi=150); plt.close(fig)
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(10, 4.2))
    ax.loglog(W, t2["lambda_star_m_eq14"], "o", color=INK, label="Eq. (14) (printed x 10)")
    ax.loglog(W, t2["lambda_star_m_as_printed"], "x", color=MUTED, label="Table 2 as printed")
    ax.loglog(W, lam_rt, "--", color=MODEL_COLOR, label="model")
    ax.set_xlabel("W [m/s2]"); ax.set_ylabel("lambda* [m]"); ax.legend(frameon=False, fontsize=8); ax.set_title("Girin & Kopyt 1994 Table 2", color=SECOND, fontsize=10)
    bx.loglog(W, t2["tau_star_s"], "o", color=INK, label="published"); bx.loglog(W, tau_rt, "--", color=MODEL_COLOR, label="model")
    bx.set_xlabel("W [m/s2]"); bx.set_ylabel("tau* [s]"); bx.legend(frameon=False, fontsize=8)
    strip_top_right_spines(ax); strip_top_right_spines(bx)
    fig.tight_layout(); fig.savefig(os.path.join(outdir, "girin1994_table2.png"), dpi=150); plt.close(fig)
    return res


def summary_markdown(r2017, r1994, densities):
    lines = ["| variant | Re density | GI | phi_cr [deg] | t_f [us] | N | r_med [um] | range [um] | t_s.d. [ms] |", "|---|---|---|---|---|---|---|---|---|"]
    fmt = lambda e, f="{:.3g}": "{} ({}{})".format(f.format(e["model"]) if e["model"] is not None else "n/a", f.format(e["published"]) if e["published"] is not None else "n/a",
                                                  ", x{:.2f}".format(e["ratio"]) if e["ratio"] else "")
    for name, per in r2017.items():
        for rd in densities:
            e = per[rd]
            rng = "{}-{}".format("{:.1f}".format(e["r_min_um"]["model"]) if e["r_min_um"]["model"] else "n/a", "{:.1f}".format(e["r_max_um"]["model"]) if e["r_max_um"]["model"] else "n/a")
            rng += " ({}-{})".format(e["r_min_um"]["published"] or "n/a", e["r_max_um"]["published"] or "n/a")
            lines.append("| {} | {} | {} | {} | {} | {} | {} | {} | {} |".format(name, rd, fmt(e["GI"], "{:.2f}"), fmt(e["phi_cr_deg"], "{:.1f}"), fmt(e["t_f_us"], "{:.1f}"),
                                                                       fmt(e["N"], "{:.2e}"), fmt(e["r_med_um"], "{:.1f}"), rng, fmt(e["t_sd_ms"], "{:.2f}")))
    lines += ["", "Values in parentheses: published (Girin 2017 Table 1) and model/published.", "",
              "Girin & Kopyt 1994 Table 1: shock-layer factor {:.2f}; r_d ratios {}; tau_d ratios {}; m ratios {}".format(
                  r1994["table1"]["shock_layer_factor"], np.round(r1994["table1"]["r_d_ratio"], 3).tolist(), np.round(r1994["table1"]["tau_d_ratio"], 3).tolist(),
                  np.round(r1994["table1"]["m_ratio"], 3).tolist()),
              "Girin & Kopyt 1994 Table 2: lambda* / Eq. (14) {}; lambda* / printed {}; tau* ratios {}".format(
                  np.round(r1994["table2"]["lambda_star_ratio_to_eq14_column"], 3).tolist(), np.round(r1994["table2"]["lambda_star_ratio_to_printed"], 3).tolist(),
                  np.round(r1994["table2"]["tau_star_ratio"], 3).tolist())]
    return "\n".join(lines) + "\n"


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--outdir", default=DEFAULT_OUTDIR)
    p.add_argument("--re-density", choices=("ambient", "shock", "both"), default="both")
    args = p.parse_args(argv)
    os.makedirs(args.outdir, exist_ok=True)
    densities = ("ambient", "shock") if args.re_density == "both" else (args.re_density,)
    r2017 = run_2017(args.outdir, densities)
    r1994 = run_1994(args.outdir)
    json.dump(r2017, open(os.path.join(args.outdir, "girin2017.json"), "w"), indent=2)
    json.dump(r1994, open(os.path.join(args.outdir, "girin1994.json"), "w"), indent=2)
    md = summary_markdown(r2017, r1994, densities)
    open(os.path.join(args.outdir, "girin_summary.md"), "w").write(md)
    print(md)
    return 0


if __name__ == "__main__":
    sys.exit(main())

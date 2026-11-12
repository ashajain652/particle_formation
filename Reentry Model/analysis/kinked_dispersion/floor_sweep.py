"""The viscous check's sweeps (spec stage 6, Task 10, Step 3): the Orr-Sommerfeld solver (`orr_sommerfeld`) on sampled
flight columns and on synthetic two-layer columns across a floor Reynolds number of 1, against the inviscid broken-line
solver and the floor test of the profile path (`reentry_model.shared.dispersion`). Analysis only; no model code changes.

    "$PY" -m analysis.kinked_dispersion.floor_sweep extract     # once: sample columns from the saved states
    "$PY" -m analysis.kinked_dispersion.floor_sweep run [--workers 4]
    "$PY" -m analysis.kinked_dispersion.floor_sweep plot          # from Reentry Model/

`extract` takes the d050 fixture's 36 columns (`data/reference_values/d050_layer_columns.npz`, d050 hand-down off at
178.5 s) and SAMPLE columns from each further saved state in STATES (each loaded once, read only), chosen as the fixture's
were (`extract_fixture`): the spec's thick set, plug columns and other wet windward columns, evenly along angle or plug
strength. It writes their `melt_profile.LayerColumns` with the state, group and angle of each to FLIGHT (small, committed),
so the sweep runs without the states.

`run` evaluates every flight column at its own speeds, and every FAST_EVERY-th one also with its speeds scaled until the
Weber number of the layer above its first kink, rho (u_top - u_kink)^2 d / sigma, is WE_TARGET (a probe of the faster
shear the spec expects under hard deceleration, where the inviscid gradient mode exists; the cells' viscosities are held,
though the slurry thins with shear, so the scaled column is a what-if); and every synthetic column of the grid below.
Each case records:

- the straight-segment fit and the floor test at the first kink: the floor Reynolds number Re_f = |u_top - u_kink| d /
  (1.225 nu_below), the slope ratio (shear below the kink over shear above it), and whether the column is floor-like
  under the defaults FLOOR_RE = 1.0 and KINK_DROP = 0.2;
- the inviscid broken-line growth of the fastest wave on the whole fitted profile, the slurry moving with the wave
  (`g_full`), and on the profile above the wall the floor test places (`g_cut`), with the profile path's viscous gate
  (growth above 2 nu_eff k^2, nu_eff the layer's Couette viscosity above the wall) applied to it (`g_model`);
- for flight columns, also the floor Reynolds number the column would have at the threshold Weber number 3
  (`re_f_star` = Re_f (3 / We_kink)^(1/2), speeds scaled together): below about 4 the floor can matter for that column
  once the mode exists; above it, the slurry is already fluid-like when surface tension first lets the mode grow;
- the viscous growth of the fastest wave in two bands of wavenumber k, with d the depth of the first kink (the top
  segment's depth where there is none): the gradient mode's band, k d from GRAD_BAND[0] to K_SPAN[1] times d over the
  thinnest fitted segment (`g_visc`, `k_visc`; `edge` where the fastest lies at an end of the band), and the long waves
  below it, from K_SPAN[0] over the column's depth to k d = GRAD_BAND[0] (`g_long`, `k_long`), where a sheared film with
  no gravity has a slow inertial mode of its own (found on a uniform column with no kink at all); and the viscous growth
  at the inviscid fastest wavenumber (`g_visc_at_kfull`). The scaled flight case, whose solves cost up to minutes each
  (its fast top cells need many more points), has instead the largest viscous growth at FAST_K times the inviscid
  fastest wavenumber (at k d = 1.225 where the inviscid profile is stable).

The synthetic columns are two-layer Couette columns under a free surface, in the top layer's units (depth d = 1, its
velocity jump 1, density 1): a top layer of viscosity 1 / Re_1 over a slurry of viscosity M / Re_1 down to a wall at
depth b; the shear stress is continuous at the kink, so the slope ratio is 1 / M, the slurry's speed at the kink
(b - 1) / M, and Re_f = Re_1 / (1.225 M); the surface tension is 1 / We. `synthetic_jobs` lists the grid: for each Re_1
of RE1, M from M_LO(Re_1) over M_DECADES decades at M_PER_DECADE points per decade, so that Re_f runs from about 80 down
to about 0.08 (or M from 1 where Re_1 is small); b = 15 and We = 1e3, with b = 3 and We = 30 as sensitivities at
Re_1 = 1e3. A slurry at high Reynolds number (M near 1 with Re_1 large) is left out: it is deep in the fluid-like regime,
where the Girin check already shows the viscous growth equal to the inviscid one, and it costs a minute a wavenumber.

Results go to OUT_DIR (`flight.csv`, `synthetic.csv`, `resolution.csv`, `timing.json`); `plot` draws the figures there."""
import os

os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "2")      # several workers share the machine
os.environ.setdefault("OMP_NUM_THREADS", "2")

import argparse  # noqa: E402
import csv  # noqa: E402
import gc  # noqa: E402
import json  # noqa: E402
import time  # noqa: E402
from dataclasses import fields  # noqa: E402
from multiprocessing import Pool  # noqa: E402

import numpy as np  # noqa: E402

from analysis.kinked_dispersion import orr_sommerfeld as osm  # noqa: E402
from reentry_model import PROJECT_ROOT  # noqa: E402
from reentry_model.shared import dispersion, melt_profile  # noqa: E402

OUT_DIR = os.path.join(PROJECT_ROOT, "reports", "kinked_dispersion_orr_sommerfeld")
FIXTURE = os.path.join(PROJECT_ROOT, "data", "reference_values", "d050_layer_columns.npz")
FLIGHT = os.path.join(OUT_DIR, "flight_columns.npz")
_SUBPLAN20 = os.path.join(os.path.dirname(os.path.dirname(PROJECT_ROOT)), "subplan20", "prototype",
                          "work-2026-10-28-skin-spraying-plan", "code", "reentry_model_output")
STATES = {                                           # label: saved state (read only); those without flow add nothing
    "d050_snapshot": os.path.join(_SUBPLAN20, "task5_snapshots", "d050__physics__skin.state.pkl"),
    "d100_snapshot": os.path.join(_SUBPLAN20, "task5_snapshots", "d100__physics__skin.state.pkl"),
    "d100_handdown_off": os.path.join(_SUBPLAN20, "task5_compare", "d100__handdown-off.state.pkl"),
    "d050_parts_melting": os.path.join(_SUBPLAN20, "task5_compare", "d050__parts-melting.state.pkl"),
    "d100_parts_melting": os.path.join(_SUBPLAN20, "task5_compare", "d100__parts-melting.state.pkl"),
    "d050_spread60": os.path.join(_SUBPLAN20, "task5_spreading", "d050__spread60.state.pkl"),
    "d100_spread60": os.path.join(_SUBPLAN20, "task5_spreading", "d100__spread60.state.pkl"),
}
FIXTURE_LABEL = "d050_fixture_178.5s"
SAMPLE = (4, 2, 2)                  # thick, plug and other columns from each further state
SIGMA = 0.80                        # N/m, the liquid's and SLURRY_SIGMA (spec 6.6)
WE_TARGET = 30.0                    # the flight columns' second case: speeds scaled until the kink's Weber number is this
WE_THRESHOLD = 3.0                  # the broken-line solver's threshold on Girin's profile (3.00)
FAST_EVERY = 4                      # every FAST_EVERY-th flight column is also run at the scaled speeds
FAST_K = (0.5, 0.71, 1.0, 1.41)     # the scaled case's wavenumbers, times the inviscid fastest
OWN_GRAD, OWN_LONG = (6, 4), (4, 2)  # the own-speed case's search grids (points, golden-section steps)
K_SPAN = dispersion.K_SPAN          # the wavenumber range of the broken-line solver's search, (0.02, 5)
GRAD_BAND = 0.3                     # k d at which the gradient mode's band starts and the long waves end
GRAD_GRID = (9, 7)                  # grid points and golden-section steps of the viscous search in the gradient band
LONG_GRID = (5, 4)                  # the same in the long-wave band
RE1 = (3.0, 10.0, 30.0, 1e2, 1e3, 1e4)   # the top layer's Reynolds number, velocity jump times depth over nu
M_DECADES, M_PER_DECADE = 3.0, 3    # the span of M = mu_slurry / mu_top per Re_1, and its points per decade
DEPTH, WEBER = 15.0, 1e3            # b / d (15 as the d050 thick columns, 830 / 55 micrometres) and rho dU^2 d / sigma
SENSITIVITY = ((1e3, 3.0, 1e3), (1e3, 15.0, 30.0))   # (Re_1, b, We): a shallow slurry; a Weber number of 30
RES_CHECK = 6                       # flight cases re-solved at doubled resolution
SYN_RES_CHECK = 4                   # synthetic cases re-solved at doubled resolution


# ---- extraction -------------------------------------------------------------------------------------------------------

def _layer_fields():
    return [f.name for f in fields(melt_profile.LayerColumns)]


def load_fixture():
    d = np.load(FIXTURE)
    cols = melt_profile.LayerColumns(**{f: d[f] for f in _layer_fields()})
    return cols, d["group"], d["theta"]


def _sample_state(path, label):
    """The sampled columns of one saved state, with their groups and angles (as `extract_fixture`)."""
    from analysis.kinked_dispersion import extract_fixture as ef
    from analysis.kinked_dispersion.shear_check import forcing, select_columns
    from reentry_model.skyscrapers import restart, skin, skinbody
    state = restart.load(path)
    body = state["run"].body
    sk = body.skins
    shear, _, _, _ = forcing(body)
    tau = np.linalg.norm(shear, axis=1)
    theta = np.asarray(body.load_angles(), dtype=float)
    cols = skinbody.layer_columns(sk, tau)
    thick = select_columns(body)
    wet = skin.flowing_tops(sk) & (np.degrees(theta) < ef.WINDWARD_DEG)
    plug, peak = ef.plug_columns(np.linalg.norm(sk.U, axis=2), cols.flowing, wet)
    is_plug = np.zeros(sk.n, dtype=bool)
    is_plug[plug] = True
    t_idx = np.flatnonzero(thick & wet & ~is_plug)
    other = np.flatnonzero(wet & ~thick & ~is_plug)
    n_thick, n_plug, n_other = SAMPLE
    chosen = [ef.evenly(t_idx, theta[t_idx], n_thick), ef.evenly(plug, peak, n_plug),
              ef.evenly(other, theta[other], n_other + max(0, n_thick - len(t_idx)))]
    idx = np.concatenate(chosen).astype(int)
    group = np.repeat([ef.GROUP_THICK, ef.GROUP_PLUG, ef.GROUP_OTHER], [len(c) for c in chosen])
    print("{}: t = {} s; wet windward {}, thick {}, plug {}, other {}; sampled {}".format(
        label, state.get("time_s"), int(wet.sum()), len(t_idx), len(plug), len(other), [len(c) for c in chosen]),
        flush=True)
    sub = cols.rows(idx)
    t = float(state.get("time_s", np.nan))
    del state, body, sk, cols
    gc.collect()
    return sub, group, theta[idx], t


def extract():
    """Write FLIGHT: the fixture's columns and SAMPLE more from each state of STATES that has wet windward columns, cut to
    their deepest cell."""
    cols, group, theta = load_fixture()
    parts = [(FIXTURE_LABEL, cols, group, theta, 178.5)]
    for label, path in STATES.items():
        sub, g, th, t = _sample_state(path, label)
        if len(g):
            parts.append(("{}_{:.1f}s".format(label, t), sub, g, th, t))
    K = max(int((p[1].dz > 0.0).sum(axis=1).max()) for p in parts)
    out = {f: [] for f in _layer_fields()}
    labels, groups, thetas = [], [], []
    for label, c, g, th, _ in parts:
        for f in _layer_fields():
            a = getattr(c, f)
            if a.ndim == 2:
                pad = {"dz": 0.0, "m": 0.0, "mu": np.inf, "flowing": False, "liquid": False}.get(f, np.nan)
                a = np.concatenate([a[:, :K], np.full((a.shape[0], max(0, K - a.shape[1])), pad, dtype=a.dtype)],
                                   axis=1)
            out[f].append(a)
        labels += [label] * len(g)
        groups.append(g)
        thetas.append(th)
    arrays = {f: np.concatenate(v) for f, v in out.items()}
    os.makedirs(OUT_DIR, exist_ok=True)
    np.savez_compressed(FLIGHT, state=np.array(labels), group=np.concatenate(groups), theta=np.concatenate(thetas),
                        **arrays)
    print("wrote {}: {} columns, {:.3f} MB".format(FLIGHT, len(labels), os.path.getsize(FLIGHT) / 1e6))


def load_flight():
    d = np.load(FLIGHT)
    cols = melt_profile.LayerColumns(**{f: d[f] for f in _layer_fields()})
    return cols, d["state"], d["group"], d["theta"]


# ---- one flight column ------------------------------------------------------------------------------------------------

def _first_kink(p, nu):
    """(index, Re_f, slope ratio) at each profile's first kink (NaN / -1 without one)."""
    kk = dispersion.kinks(p)
    re = dispersion.floor_reynolds(p, nu)
    s = np.abs(p.slopes())
    rr = np.arange(p.count.size)
    j = np.where(kk.any(axis=1), np.argmax(kk, axis=1), -1)
    jj = np.maximum(j, 0)
    ratio = np.where(j >= 0, s[rr, jj + 1] / s[rr, jj], np.nan)
    return j, np.where(j >= 0, re[rr, jj], np.nan), ratio


def _scaled(cols, s):
    c = cols.rows(np.arange(cols.front.size))
    c.u = cols.u * s
    return c


def inviscid(cols, sigma=SIGMA, floor_re=dispersion.FLOOR_RE, kink_drop=dispersion.KINK_DROP):
    """The profile path's view of LayerColumns `cols` (all reducible): a dict of (n,) arrays -- the fit's first kink
    (depth, Re_f, slope ratio), the wall, the inviscid growth on the whole profile and above the wall, and the gate."""
    fit = melt_profile.fit_profiles(cols)
    p = fit.profiles
    nu = melt_profile.nu_below(cols, fit.cell)
    j, re_f, ratio = _first_kink(p, nu)
    wall = dispersion.place_wall(p, nu, floor_re=floor_re, kink_drop=kink_drop)
    cut = dispersion.above_wall(p, wall)
    n = p.count.size
    full = dispersion.fastest_waves(p, sigma)
    g_cut, k_cut = np.zeros(n), np.full(n, np.nan)
    curved = np.flatnonzero(cut.count > 0)
    if curved.size:
        w = dispersion.fastest_waves(cut.rows(curved), sigma)
        g_cut[curved], k_cut[curved] = w.growth, w.k
    nu_eff = melt_profile.couette_viscosity(cols, wall.depth)[0]
    gate = 2.0 * nu_eff * np.where(np.isfinite(k_cut), k_cut, 0.0) ** 2
    top, bot, _, _ = p.segments()
    seg = bot - top
    s_min = np.nanmin(np.where(np.isfinite(seg) & (seg > 0.0), seg, np.nan), axis=1)
    rr = np.arange(n)
    d_top = np.where(p.count > 0, p.depth[:, 0], p.wall)
    u_low = np.where(p.count > 0, p.u_break[:, 0], p.u_wall)
    we_top = p.rho[:, 0] * (p.u_top - u_low) ** 2 * d_top / sigma
    jj = np.maximum(j, 0)
    d_k = np.where(j >= 0, p.depth[rr, jj], d_top)
    u_k = np.where(j >= 0, p.u_break[rr, jj], u_low)
    we_kink = p.rho[:, 0] * (p.u_top - u_k) ** 2 * d_k / sigma
    return dict(count=p.count, u_top=p.u_top, d_top=d_top, we_top=we_top, we_kink=we_kink,
                re_f_star=re_f * np.sqrt(WE_THRESHOLD / we_kink),
                kink_depth=np.where(j >= 0, p.depth[rr, np.maximum(j, 0)], np.nan),
                re_f=re_f, slope_ratio=ratio, floor=wall.floor, wall_depth=wall.depth, g_full=full.growth,
                k_full=full.k, g_cut=g_cut, k_cut=k_cut, gate=gate, g_model=np.where(g_cut > gate, g_cut, 0.0),
                s_min=s_min, front=p.wall)


def _two_bands(y, U, mu, rho, sigma, d, k_lo, k_hi, u_faces=None):
    """The viscous fastest wave in the gradient mode's band (k d from GRAD_BAND up to k_hi) and in the long waves
    (k_lo up to k d = GRAD_BAND), each an `orr_sommerfeld.Fastest` (the second None if the band is empty)."""
    k_mid = GRAD_BAND / d
    grad = osm.fastest(y, U, mu, rho, sigma, (max(k_mid, k_lo), k_hi), u_faces=u_faces, n_grid=GRAD_GRID[0],
                       n_refine=GRAD_GRID[1])
    long = None
    if k_lo < 0.99 * k_mid:
        long = osm.fastest(y, U, mu, rho, sigma, (k_lo, k_mid), u_faces=u_faces, n_grid=LONG_GRID[0],
                           n_refine=LONG_GRID[1])
    return grad, long


def _visc_fields(grad, long):
    return dict(g_visc=grad.growth, k_visc=grad.k, c_visc_re=grad.c.real, c_visc_im=grad.c.imag, edge=int(grad.edge),
                g_long=long.growth if long else np.nan, k_long=long.k if long else np.nan)


def flight_case(args):
    """One flight column, at its own speeds ("own") or scaled to the kink's Weber number WE_TARGET ("fast"): the
    inviscid view and the viscous growth (a dict, one CSV row)."""
    i, mode, label, group, theta = args
    cols, _, _, _ = load_flight()
    c = cols.rows([i])
    s = 1.0
    if mode == "fast":
        we = inviscid(c)["we_kink"][0]
        s = float(np.sqrt(WE_TARGET / we))
        c = _scaled(c, s)
    inv = {k: (v[0].item() if hasattr(v[0], "item") else v[0]) for k, v in inviscid(c).items()}
    y, U, mu, rho = osm.column_cells(c, 0)
    depth = y[-1]
    s_min = inv["s_min"] if np.isfinite(inv["s_min"]) else float(np.diff(y).min())
    d = inv["kink_depth"] if np.isfinite(inv["kink_depth"]) else inv["d_top"]
    t0 = time.perf_counter()
    at_full = np.nan
    if mode == "own":
        k_mid = GRAD_BAND / d
        grad = osm.fastest(y, U, mu, rho, SIGMA, (k_mid, K_SPAN[1] / s_min), n_grid=OWN_GRAD[0], n_refine=OWN_GRAD[1])
        long = osm.fastest(y, U, mu, rho, SIGMA, (K_SPAN[0] / depth, k_mid), n_grid=OWN_LONG[0],
                           n_refine=OWN_LONG[1]) if K_SPAN[0] / depth < 0.99 * k_mid else None
        visc = _visc_fields(grad, long)
    else:
        ks = inv["k_full"] * np.array(FAST_K) if np.isfinite(inv["k_full"]) else np.array([dispersion.FASTEST_KD / d])
        cs = [osm.growth_rate(k, y, U, mu, rho, SIGMA) for k in ks]
        g = np.array([k * cc.imag for k, cc in zip(ks, cs)])
        b = int(np.argmax(g))
        visc = dict(g_visc=float(g[b]), k_visc=float(ks[b]), c_visc_re=cs[b].real, c_visc_im=cs[b].imag,
                    edge=int(b in (0, len(ks) - 1)), g_long=np.nan, k_long=np.nan)
        if np.isfinite(inv["k_full"]):
            at_full = float(g[list(FAST_K).index(1.0)])
    if mode == "own" and np.isfinite(inv["k_full"]):
        at_full = inv["k_full"] * osm.growth_rate(inv["k_full"], y, U, mu, rho, SIGMA).imag
    row = dict(state=label, column=i, group=int(group), theta_deg=float(np.degrees(theta)), mode=mode, scale=s,
               cells=U.size, depth_m=depth, d_ref=d, **inv, **visc, g_visc_at_kfull=at_full,
               seconds=time.perf_counter() - t0)
    print("flight {} col {} {} (x{:.3g}): Re_f {:.3g} (at We 3: {:.3g}), ratio {:.3g}, We_kink {:.3g}, g_full {:.4g}, "
          "g_model {:.4g}, g_visc {:.4g} at k d {:.3g} ({:.0f} s)".format(
              label, i, mode, s, inv["re_f"], inv["re_f_star"], inv["slope_ratio"], inv["we_kink"], inv["g_full"],
              inv["g_model"], visc["g_visc"], visc["k_visc"] * d, row["seconds"]), flush=True)
    return row


def resolution_case(args):
    """A flight case's viscous eigenvalue at its fastest wavenumber at resolutions 1 and 2."""
    i, s, label, k = args
    cols, _, _, _ = load_flight()
    c = _scaled(cols.rows([i]), s)
    y, U, mu, rho = osm.column_cells(c, 0)
    c1 = osm.growth_rate(k, y, U, mu, rho, SIGMA)
    c2 = osm.growth_rate(k, y, U, mu, rho, SIGMA, resolution=2)
    return dict(kind="flight", case="{} col {} x{:g}".format(label, i, s), k=k, c1_re=c1.real, c1_im=c1.imag,
                c2_re=c2.real, c2_im=c2.imag, growth1=k * c1.imag, growth2=k * c2.imag, speed=float(abs(U[0])),
                diff_over_speed=abs(c2 - c1) / abs(U[0]))


# ---- one synthetic column ---------------------------------------------------------------------------------------------

def synthetic_column(re1, M, b, n_top=6, n_slurry=8):
    """(y, U, mu, rho, u_faces) of a two-layer Couette column in the top layer's units (module docstring); the slurry's
    cells thin geometrically toward the kink."""
    uk = (b - 1.0) / M
    w = 1.6 ** np.arange(n_slurry)
    edges = np.concatenate([np.linspace(0.0, 1.0, n_top + 1), 1.0 + np.cumsum(w) / w.sum() * (b - 1.0)])
    uf = np.where(edges <= 1.0, uk + 1.0 - edges, uk * (b - edges) / (b - 1.0))
    U = 0.5 * (uf[:-1] + uf[1:])
    mu = np.where(edges[:-1] < 1.0, 1.0 / re1, M / re1)
    return edges, U, mu, np.ones(U.size), uf


def synthetic_jobs():
    """(Re_1, M, b, We) of every synthetic case (module docstring)."""
    jobs = []
    for re1, b, we in [(r, DEPTH, WEBER) for r in RE1] + list(SENSITIVITY):
        m_lo = max(1.0, re1 / (dispersion.FASTEST_KD * 80.0))
        for j in range(int(M_DECADES * M_PER_DECADE) + 1):
            jobs.append((re1, float(np.round(m_lo * 10.0 ** (j / M_PER_DECADE), 6)), b, we))
    return jobs


def synthetic_case(args):
    re1, M, b, we = args
    y, U, mu, rho, uf = synthetic_column(re1, M, b)
    sigma = 1.0 / we
    uk = (b - 1.0) / M
    prof = dispersion.KinkedProfiles.one([1.0], [uk + 1.0, uk, 0.0], wall=b, rho=1.0)
    full = dispersion.fastest_waves(prof, sigma)
    nu_eff = b / (1.0 * re1 + (b - 1.0) * re1 / M)          # the whole column's Couette viscosity (rho = 1)
    re_f = re1 / (dispersion.FASTEST_KD * M)
    t0 = time.perf_counter()
    grad, long = _two_bands(y, U, mu, rho, sigma, 1.0, K_SPAN[0] / b, K_SPAN[1], u_faces=uf)
    at_full = np.nan
    if np.isfinite(full.k[0]):
        at_full = full.k[0] * osm.growth_rate(full.k[0], y, U, mu, rho, sigma, u_faces=uf).imag
    row = dict(re1=re1, M=M, b=b, we=we, re_f=re_f, slope_ratio=1.0 / M, g_full=float(full.growth[0]),
               k_full=float(full.k[0]), gate_fluid=2.0 * nu_eff * (full.k[0] if np.isfinite(full.k[0]) else 0.0) ** 2,
               **_visc_fields(grad, long), g_visc_at_kfull=at_full, seconds=time.perf_counter() - t0)
    print("synthetic Re1 {:g} M {:g} b {:g} We {:g}: Re_f {:.3g}, g_full {:.4g}, g_visc {:.4g} at k {:.3g}{}, g_long "
          "{:.3g} (at k_full {:.4g}; {:.0f} s)".format(re1, M, b, we, re_f, row["g_full"], grad.growth, grad.k,
                                                        " (edge)" if grad.edge else "", row["g_long"], at_full,
                                                        row["seconds"]), flush=True)
    return row


def synthetic_resolution(args):
    re1, M, b, we, k = args
    y, U, mu, rho, uf = synthetic_column(re1, M, b)
    c1 = osm.growth_rate(k, y, U, mu, rho, 1.0 / we, u_faces=uf)
    c2 = osm.growth_rate(k, y, U, mu, rho, 1.0 / we, u_faces=uf, resolution=2)
    return dict(kind="synthetic", case="Re1 {:g} M {:g} b {:g} We {:g}".format(re1, M, b, we), k=k, c1_re=c1.real,
                c1_im=c1.imag, c2_re=c2.real, c2_im=c2.imag, growth1=k * c1.imag, growth2=k * c2.imag,
                speed=float(abs(uf[0])), diff_over_speed=abs(c2 - c1) / abs(uf[0]))


def _cell(v):
    """A CSV cell: floats in full (repr), numpy scalars as plain numbers, booleans as 0 / 1."""
    if isinstance(v, (bool, np.bool_)):
        return int(v)
    if isinstance(v, (float, np.floating)):
        return repr(float(v))
    if isinstance(v, np.integer):
        return int(v)
    return v


def _write(rows, name):
    path = os.path.join(OUT_DIR, name)
    keys = list(rows[0])
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        for r in rows:
            w.writerow({k: _cell(v) for k, v in r.items()})
    return path


def resolution_picks(fl):
    """The flight cases re-solved at doubled resolution: RES_CHECK own-speed cases evenly through the list, and the
    scaled case with the fewest cells times speed factor (a large scaled column at doubled resolution needs a dense solve
    of over ten thousand unknowns: one ran for over 25 minutes without finishing, 2026-11-12)."""
    if isinstance(fl, dict):
        rows = [{k: fl[k][j] for k in fl} for j in range(len(fl["mode"]))]
    else:
        rows = fl
    own = sorted([r for r in rows if r["mode"] == "own"], key=lambda r: r["column"])
    pick = [own[j] for j in np.unique(np.linspace(0, len(own) - 1, RES_CHECK).astype(int))]
    fast = [r for r in rows if r["mode"] == "fast"]
    if fast:
        pick.append(min(fast, key=lambda r: r["cells"] * r["scale"]))
    return [(int(r["column"]), float(r["scale"]), str(r["state"]), float(r["k_visc"])) for r in pick]


def resolution_only(workers=4):
    """Re-solve `resolution_picks` of an existing flight.csv at doubled resolution (resolution_flight.csv)."""
    t0 = time.time()
    with Pool(workers) as pool:
        res = pool.map(resolution_case, resolution_picks(_read("flight.csv")), chunksize=1)
    _write(res, "resolution_flight.csv")
    print("flight resolution checks: {:.0f} s".format(time.time() - t0))


def run(workers=4, only=None):
    os.makedirs(OUT_DIR, exist_ok=True)
    cols, states, groups, thetas = load_flight()
    reducible = melt_profile.reducible(cols)
    timing = {}
    res = []
    with Pool(workers) as pool:
        if only in (None, "synthetic"):
            t0 = time.time()
            jobs = synthetic_jobs()
            jobs.sort(key=lambda j: j[1] / j[0])             # the costliest (a slurry at high Reynolds number) first
            syn = pool.map(synthetic_case, jobs, chunksize=1)
            syn.sort(key=lambda r: (r["b"], -r["we"], r["re1"], r["M"]))
            timing["synthetic_s"] = time.time() - t0
            timing["synthetic_cases"] = len(jobs)
            _write(syn, "synthetic.csv")
            pick = [syn[j] for j in np.unique(np.linspace(0, len(syn) - 1, SYN_RES_CHECK + 2).astype(int))[1:-1]]
            t0 = time.time()
            res += pool.map(synthetic_resolution, [(r["re1"], r["M"], r["b"], r["we"], r["k_visc"]) for r in pick],
                            chunksize=1)
            timing["synthetic_resolution_s"] = time.time() - t0
        if only in (None, "flight"):
            t0 = time.time()
            red = np.flatnonzero(reducible)
            jobs = [(i, "own", states[i], groups[i], thetas[i]) for i in red]
            jobs += [(i, "fast", states[i], groups[i], thetas[i]) for i in red[::FAST_EVERY]]
            jobs.sort(key=lambda j: (j[1] != "fast", -int(cols.flowing[j[0]].sum())))    # the costliest first
            fl = pool.map(flight_case, jobs, chunksize=1)
            timing["flight_s"] = time.time() - t0
            timing["flight_cases"] = len(jobs)
            fl.sort(key=lambda r: (r["mode"] != "own", r["column"]))
            _write(fl, "flight.csv")
            t0 = time.time()
            res += pool.map(resolution_case, resolution_picks(fl), chunksize=1)
            timing["flight_resolution_s"] = time.time() - t0
    if res:
        _write(res, "resolution.csv" if only is None else "resolution_{}.csv".format(only))
    timing["workers"] = workers
    with open(os.path.join(OUT_DIR, "timing.json" if only is None else "timing_{}.json".format(only)), "w") as fh:
        json.dump(timing, fh, indent=1)
    print(timing)


# ---- analysis and figures ---------------------------------------------------------------------------------------------

INK, INK2, MUTED, GRID, SURFACE = "#0b0b0b", "#52514e", "#8a8984", "#e4e3df", "#fcfcfb"
BLUES = ("#86b6ef", "#5598e7", "#2a78d6", "#1c5cab", "#104281", "#0d366b")    # sequential blue, light to dark
SERIES = ("#2a78d6", "#eb6834", "#1baf7a")                                    # categorical slots 1-3
CANDIDATE_FLOOR_RE = (0.1, 0.3, 1.0, 3.0, 10.0, 30.0)
CANDIDATE_KINK_DROP = (0.05, 0.1, 0.2, 0.5, 1.0)
SURVIVES = 0.5          # the viscous gradient mode "survives" where it keeps at least this share of the inviscid growth


def _read(name):
    path = os.path.join(OUT_DIR, name)
    if not os.path.isfile(path):
        return None
    with open(path, newline="") as fh:
        rows = list(csv.DictReader(fh))
    out = {}
    for k in rows[0]:
        vals = [r[k] for r in rows]
        vals = [v[len("np.float64("):-1] if v.startswith("np.float64(") else v for v in vals]
        try:
            out[k] = np.array([float({"True": 1, "False": 0}.get(v, v)) if v not in ("", "None") else np.nan
                               for v in vals])
        except ValueError:
            out[k] = np.array(vals)
    return out


def crossing(x, y, level):
    """The x (log-interpolated) where y, read along falling x, first drops below level; NaN if it never does."""
    order = np.argsort(-x)
    x, y = x[order], y[order]
    for j in range(1, len(x)):
        if y[j - 1] >= level > y[j]:
            t = (y[j - 1] - level) / (y[j - 1] - y[j])
            return float(np.exp(np.log(x[j - 1]) + t * (np.log(x[j]) - np.log(x[j - 1]))))
    return float("nan")


def model_growth(syn, floor_re, kink_drop, gate=True):
    """The profile path's growth on each synthetic column under (floor_re, kink_drop): 0 where the kink is a floor
    (the straight top layer on a floor is stable), otherwise the inviscid growth on the whole profile where it beats the
    viscous gate 2 nu_eff k^2 (with `gate`), else 0."""
    floor = (syn["slope_ratio"] <= kink_drop) & (syn["re_f"] < floor_re)
    fluid = np.where(syn["g_full"] > syn["gate_fluid"], syn["g_full"], 0.0) if gate else syn["g_full"]
    return np.where(floor, 0.0, fluid)


def analyse():
    """summary.json: per Re_1, the floor Reynolds numbers at which the viscous growth falls to SURVIVES, 0.1 and 0.01
    of the inviscid; and, for each candidate (FLOOR_RE, KINK_DROP), how often the profile path's verdict (sprays by the
    gradient mode or not) agrees with the viscous one on the main grid, and its mean growth error over the inviscid
    growth."""
    syn = _read("synthetic.csv")
    out = {"survives": SURVIVES, "series": {}, "candidates": []}
    main = (syn["b"] == DEPTH) & (syn["we"] == WEBER)
    groups = [("Re1 {:g}".format(r), main & (syn["re1"] == r)) for r in RE1]
    groups += [("Re1 {:g} b {:g} We {:g}".format(*c), (syn["re1"] == c[0]) & (syn["b"] == c[1]) & (syn["we"] == c[2]))
               for c in SENSITIVITY]
    for name, m in groups:
        ok = m & (syn["g_full"] > 0.0)
        ratio = syn["g_visc"][ok] / syn["g_full"][ok]
        out["series"][name] = dict(
            re_f_half=crossing(syn["re_f"][ok], ratio, SURVIVES), re_f_tenth=crossing(syn["re_f"][ok], ratio, 0.1),
            re_f_hundredth=crossing(syn["re_f"][ok], ratio, 0.01), ratio_max=float(np.max(ratio)) if ratio.size else None,
            ratio_min=float(np.min(ratio)) if ratio.size else None)
    # the one-parameter curve Re_f / (Re_f + c) fitted to the ratio on the main grid where the top layer is close to
    # inviscid (Re_1 >= 30) and the kink a real one (M >= 2, a slope ratio of 0.5 or less; at M = 1 there is no kink and
    # both growths are small), over Re_f from 0.3 to 80 (below it the residual film mode of the high-Re_1 series lifts
    # the ratio off the curve)
    fitm = main & (syn["g_full"] > 0.0) & (syn["re1"] >= 30.0) & (syn["re_f"] >= 0.3) & (syn["M"] >= 2.0)
    rf, rt = syn["re_f"][fitm], syn["g_visc"][fitm] / syn["g_full"][fitm]
    cs = np.logspace(-0.5, 1.5, 401)
    err = [np.sqrt(np.mean((rf / (rf + c) - rt) ** 2)) for c in cs]
    jbest = int(np.argmin(err))
    out["smooth_fit"] = dict(c=float(cs[jbest]), rms=float(err[jbest]), cases=int(fitm.sum()),
                             max_abs=float(np.max(np.abs(rf / (rf + cs[jbest]) - rt))))
    m = main & (syn["g_full"] > 0.0)
    truth = syn["g_visc"][m] >= SURVIVES * syn["g_full"][m]
    for fr in CANDIDATE_FLOOR_RE:
        for kd in CANDIDATE_KINK_DROP:
            g = model_growth(syn, fr, kd)[m]
            verdict = g > 0.0
            out["candidates"].append(dict(floor_re=fr, kink_drop=kd, agree=int((verdict == truth).sum()),
                                          cases=int(m.sum()),
                                          false_spray=int((verdict & ~truth).sum()),
                                          missed_spray=int((~verdict & truth).sum()),
                                          mean_abs_error=float(np.mean(np.abs(g - syn["g_visc"][m]) / syn["g_full"][m]))))
    # the smooth alternative on the same cases: the inviscid growth times Re_f / (Re_f + c), with and without the gate
    c = out["smooth_fit"]["c"]
    factor = syn["re_f"] / (syn["re_f"] + c)
    for name, g in (("smooth_gated", np.where(syn["g_full"] > syn["gate_fluid"], syn["g_full"], 0.0) * factor),
                    ("smooth_ungated", syn["g_full"] * factor)):
        out[name] = dict(mean_abs_error=float(np.mean(np.abs(g[m] - syn["g_visc"][m]) / syn["g_full"][m])),
                         mean_abs_error_re1_ge_30=float(np.mean((np.abs(g - syn["g_visc"]) / syn["g_full"])[
                             m & (syn["re1"] >= 30.0)])) if np.any(m & (syn["re1"] >= 30.0)) else None)
    for cand in out["candidates"]:
        g = model_growth(syn, cand["floor_re"], cand["kink_drop"])
        cand["mean_abs_error_re1_ge_30"] = float(np.mean((np.abs(g - syn["g_visc"]) / syn["g_full"])[
            m & (syn["re1"] >= 30.0)]))
        g = model_growth(syn, cand["floor_re"], cand["kink_drop"], gate=False)
        cand["ungated_agree"] = int(((g[m] > 0.0) == truth).sum())
        cand["ungated_mean_abs_error"] = float(np.mean(np.abs(g[m] - syn["g_visc"][m]) / syn["g_full"][m]))
    fl = _read("flight.csv")
    if fl is not None:
        out["flight"] = {}
        for name in ("own", "fast"):
            f = fl["mode"] == name
            gated = f & (fl["g_cut"] > 0.0) & (fl["g_model"] == 0.0)
            out["flight"][name + "_gate"] = dict(
                killed_by_gate=int(gated.sum()),
                viscous_over_inviscid_there=[float(x) for x in np.round((fl["g_visc"] / fl["g_cut"])[gated], 3)],
                gate_over_growth_there=[float(x) for x in np.round((fl["gate"] / fl["g_cut"])[gated], 2)])
            out["flight"][name] = dict(
                cases=int(f.sum()), floor_like=int(np.sum(fl["floor"][f] == 1.0)),
                inviscid_unstable=int(np.sum(fl["g_full"][f] > 0.0)), model_sprays=int(np.sum(fl["g_model"][f] > 0.0)),
                viscous_positive=int(np.sum(fl["g_visc"][f] > 0.0)),
                g_visc_median=float(np.nanmedian(fl["g_visc"][f])), g_visc_max=float(np.nanmax(fl["g_visc"][f])),
                re_f_median=float(np.nanmedian(fl["re_f"][f])), we_kink_median=float(np.nanmedian(fl["we_kink"][f])),
                re_f_star_quartiles=[float(q) for q in np.nanpercentile(fl["re_f_star"][f], [10, 25, 50, 75, 90])])
    if fl is not None:
        out["flight_candidates"] = flight_candidates(fl)
    with open(os.path.join(OUT_DIR, "summary.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    return out


def flight_candidates(fl):
    """On the scaled flight cases, the profile path's inviscid growth under each candidate (FLOOR_RE, KINK_DROP), with
    and without the viscous gate, against the viscous growth (which does not depend on the constants): how often the
    verdict -- the gradient mode keeps at least SURVIVES of the whole profile's inviscid growth, or not -- agrees, and
    the mean of |model - viscous| over the whole profile's inviscid growth. Also the walls' median depth."""
    cols, _, _, _ = load_flight()
    f = np.flatnonzero((fl["mode"] == "fast") & (fl["g_full"] > 0.0))
    if f.size == 0:
        return []
    idx = fl["column"][f].astype(int)
    c = cols.rows(idx)
    c.u = cols.u[idx] * fl["scale"][f][:, None]
    g_full, g_visc = fl["g_full"][f], fl["g_visc"][f]
    truth = g_visc >= SURVIVES * g_full
    out = []
    for fr in CANDIDATE_FLOOR_RE + (3.5,):
        for kd in CANDIDATE_KINK_DROP:
            inv = inviscid(c, floor_re=fr, kink_drop=kd)
            for gate in (True, False):
                g = inv["g_model"] if gate else inv["g_cut"]
                out.append(dict(floor_re=fr, kink_drop=kd, gate=gate, cases=int(f.size),
                                agree=int(((g >= SURVIVES * g_full) == truth).sum()),
                                mean_abs_error=float(np.mean(np.abs(g - g_visc) / g_full)),
                                floor_like=int(inv["floor"].sum()),
                                wall_median_um=float(np.median(inv["wall_depth"]) * 1e6)))
    return out


def _style(ax):
    ax.set_facecolor(SURFACE)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(MUTED)
    ax.tick_params(colors=INK2, labelsize=9)
    ax.grid(True, which="major", color=GRID, lw=0.6)
    ax.set_axisbelow(True)


def plot():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    summary = analyse()
    syn = _read("synthetic.csv")
    plt.rcParams.update({"font.size": 9, "axes.labelcolor": INK, "text.color": INK, "figure.facecolor": SURFACE})

    # Figure 1: the synthetic two-layer columns, viscous over inviscid growth against Re_f and the slope ratio
    fig, axes = plt.subplots(1, 2, figsize=(10.0, 4.2), sharey=True)
    main = (syn["b"] == DEPTH) & (syn["we"] == WEBER)
    for ax, key, xlabel, ref, reflabel in (
            (axes[0], "re_f", "floor Reynolds number at the kink, Re_f", dispersion.FLOOR_RE, "FLOOR_RE = 1.0"),
            (axes[1], "slope_ratio", "slope ratio at the kink (shear below over shear above)", dispersion.KINK_DROP,
             "KINK_DROP = 0.2")):
        _style(ax)
        for colour, re1 in zip(BLUES, RE1):
            m = main & (syn["re1"] == re1) & (syn["g_full"] > 0.0) & (syn["M"] >= 2.0)
            o = np.argsort(syn[key][m])
            ax.plot(syn[key][m][o], (syn["g_visc"] / syn["g_full"])[m][o], "-o", color=colour, lw=2, ms=4.5,
                    mec=SURFACE, mew=1.0, label="Re_1 = {:g}".format(re1))
        for (re1, b, we), ls in zip(SENSITIVITY, ("--", ":")):
            m = (syn["re1"] == re1) & (syn["b"] == b) & (syn["we"] == we) & (syn["g_full"] > 0.0) & (syn["M"] >= 2.0)
            o = np.argsort(syn[key][m])
            ax.plot(syn[key][m][o], (syn["g_visc"] / syn["g_full"])[m][o], ls, color=INK2, lw=1.4,
                    label="Re_1 = {:g}, b = {:g} d, We = {:g}".format(re1, b, we))
        if key == "re_f":
            xs = np.logspace(-2.7, 2.0, 200)
            c = summary["smooth_fit"]["c"]
            ax.plot(xs, xs / (xs + c), color=SERIES[1], lw=1.2, ls="-", alpha=0.9,
                    label="Re_f / (Re_f + {:.1f}), fitted".format(c))
        ax.axvline(ref, color=INK2, lw=1.0, ls=(0, (4, 3)))
        ax.text(ref, 1.06, reflabel + " (assumed)", color=INK2, fontsize=8, ha="center", va="bottom")
        ax.axhline(SURVIVES, color=MUTED, lw=0.8, ls=":")
        ax.set_xscale("log")
        ax.set_xlabel(xlabel)
        ax.set_ylim(-0.05, 1.15)
    axes[0].set_ylabel("viscous growth / inviscid growth (whole profile)")
    axes[0].legend(frameon=False, fontsize=7.5, loc="upper left")
    axes[1].text(0.02, 0.97, "columns with a viscosity contrast M of 2 or more\n(M below 2 is no real kink)",
                 transform=axes[1].transAxes, ha="left", va="top", color=INK2, fontsize=8)
    fig.suptitle("Two-layer Couette columns under a free surface: where the slurry stops moving with the wave",
                 color=INK, fontsize=10)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT_DIR, "synthetic_growth_ratio.png"), dpi=130)
    plt.close(fig)

    # Figure 2: the flight columns -- at their own speeds, and scaled to a kink Weber number of WE_TARGET
    fl = _read("flight.csv")
    if fl is not None:
        fig, axes = plt.subplots(1, 2, figsize=(10.0, 4.2))
        ax = axes[0]
        _style(ax)
        own = (fl["mode"] == "own") & np.isfinite(fl["re_f"])
        floor = fl["floor"] == 1.0
        g = fl["g_visc"]
        for sel, filled, label in ((own & floor & (g > 0), True, "floor-like under the defaults"),
                                   (own & ~floor & (g > 0), False, "fluid-like under the defaults")):
            ax.scatter(fl["re_f"][sel], g[sel], s=28, facecolor=SERIES[0] if filled else SURFACE, edgecolor=SERIES[0],
                       lw=1.4, zorder=3, label=label)
        low = own & (g <= 0)
        if low.any():
            ax.scatter(fl["re_f"][low], np.full(low.sum(), 1e-3), s=22, marker="v", color=INK2, zorder=3,
                       label="no growing wave (drawn at 1e-3)")
        ax.axvline(dispersion.FLOOR_RE, color=INK2, lw=1.0, ls=(0, (4, 3)))
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("floor Reynolds number at the first kink, Re_f")
        ax.set_ylabel("viscous growth of the fastest wave, 1/s")
        ax.set_title("at their own speeds (gradient band)", color=INK, fontsize=9)
        ax.legend(frameon=False, fontsize=7.5, loc="upper left")
        ax = axes[1]
        _style(ax)
        main = (syn["b"] == DEPTH) & (syn["we"] == WEBER) & (syn["g_full"] > 0) & (syn["M"] >= 2.0)
        ax.scatter(syn["re_f"][main], (syn["g_visc"] / syn["g_full"])[main], s=12, color=GRID, edgecolor=MUTED,
                   lw=0.5, zorder=2, label="synthetic two-layer columns")
        xs = np.logspace(-2.0, 2.3, 200)
        c = summary["smooth_fit"]["c"]
        ax.plot(xs, xs / (xs + c), color=SERIES[1], lw=1.2, zorder=2, label="Re_f / (Re_f + {:.1f})".format(c))
        fast = (fl["mode"] == "fast") & (fl["g_full"] > 0) & np.isfinite(fl["re_f"])
        ax.scatter(fl["re_f"][fast], (fl["g_visc"] / fl["g_full"])[fast], s=34, facecolor=SERIES[0],
                   edgecolor=SURFACE, lw=1.0, zorder=3,
                   label="flight columns at kink We = {:g}".format(WE_TARGET))
        ax.axvline(dispersion.FLOOR_RE, color=INK2, lw=1.0, ls=(0, (4, 3)))
        ax.set_xscale("log")
        ax.set_ylim(-0.05, 1.15)
        ax.set_xlabel("floor Reynolds number at the first kink, Re_f")
        ax.set_ylabel("viscous growth / inviscid growth (whole profile)")
        ax.set_title("speeds scaled to a kink Weber number of {:g}".format(WE_TARGET), color=INK, fontsize=9)
        ax.legend(frameon=False, fontsize=7.5, loc="upper left")
        fig.suptitle("Sampled flight columns (d050 at 178.5 s, d100 at 32.5 s)", color=INK, fontsize=10)
        fig.tight_layout()
        fig.savefig(os.path.join(OUT_DIR, "flight_growth.png"), dpi=130)
        plt.close(fig)
    print(json.dumps(summary["series"], indent=1))


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("action", choices=("extract", "run", "resolution", "plot"))
    p.add_argument("--workers", type=int, default=4)
    p.add_argument("--only", choices=("synthetic", "flight"), default=None)
    a = p.parse_args(argv)
    if a.action == "extract":
        extract()
    elif a.action == "run":
        run(a.workers, a.only)
    elif a.action == "resolution":
        resolution_only(a.workers)
    else:
        plot()


if __name__ == "__main__":
    main()

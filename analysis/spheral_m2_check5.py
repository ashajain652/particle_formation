"""Spheral M2 check 5: a heated bar held at both ends (M2 plan, Task 7; spec §12 check 5).

    "$PY" analysis/spheral_m2_check5.py [--only fe|run|analyse] [--max-load 2.0] [--force]

**The case.** A bar along the flight axis x, 100 mm long, 20 mm square (3D) or a 10 mm-radius rod on the axis (RZ),
at 2.2 and 1.1 mm (`body --geometry bar`: particles at cell centres, so the bar is n dx long and wide -- 99 x 19.8
mm and RZ radius 11 mm at 2.2 mm, 100.1 x 19.8 mm and radius 9.9 mm at 1.1 mm); the particles within two spacings
of each end held at rest (the runner's Hold package; run.py "Held particles": exact, the held particles' largest
|v| and |x - x_ref| are 0.0 in every run). Uniform heating from 300 K by 300 K over 10 intervals (the body's
`--synthetic uniform`, 30 K per 0.5 s interval, the MVP flight's frame times 0-10), `--strength constant:100e6`,
`--loads off`, the mass-scaled replay at the defaults (10 crossings, damping 0.05, the crossings monitor every 10:
intervals 0->1 and 9->10). `--ke-fraction 1e300` as in check 4: the study measures the default crossings as given
and records the KE ratios instead of letting the slowness test repeat intervals.

**The material over 300-600 K** (checked by `constants`): the table's E 71.7 GPa, nu 0.33, K, G are constant, the
yield 100 MPa by construction, and the free density gives the free linear strain exactly alpha_eff (T - 300 K),
alpha_eff = alpha / (1 + alpha (300 K - 293 K)) = 2.33962e-5 1/K (rho_free = rho_s / (1 + alpha (T - 293))^3), to
4e-16. So the closed forms hold with constant properties: E alpha_eff dT reaches Y at dT_y = 59.61 K, inside the
second interval.

**The expectation** (plan Task 7 Step 2), in the middle third (|x| < n_x dx / 6): sigma_xx = -E alpha_eff dT until
-Y, then -Y; the lateral stresses zero; the equivalent plastic strain alpha_eff dT - Y/E once yielded. Two
qualifications, both measured here:

1. *The elastic phase is not uniaxial over the whole bar.* Near the held ends the lateral expansion is prevented,
   and a laterally constrained element pushes axially 1/(1 - 2 nu) = 2.9 times harder than a free one; since the
   axial force is the same along the bar, the middle carries more than E alpha_eff dT. The finite-element reference
   (`bar_fe`: linear thermoelasticity, P2 on gmsh tetrahedra, the bar between the held layers -- n_x dx - 4 dx long,
   its end faces clamped -- with the particle bar's section, sampled at the middle third's particle positions)
   gives the ratio sigma_xx / (-E alpha_eff dT) of the middle (`fe_ratio`): the elastic expectation is
   fe_ratio (-E alpha_eff dT).
2. *Perfect plasticity does not fix where the bar flows.* Once the whole bar is at -Y, equilibrium and the yield
   condition hold for any distribution of plastic shortening along it; only its integral is fixed (the held length
   does not change). The bar yields first near the ends, where the constrained zones raise the stress, and that
   early flow is never undone, so the middle third's plastic strain lags alpha_eff dT - Y/E by a near-constant
   offset while its *increment* per interval, once all of the middle yields, is the uniaxial alpha_eff 30 K =
   7.02e-4. Measured: the increment per interval (`ps_increment`) against 7.02e-4, and the definition check below.

**Spheral's plastic strain** (src/Strength/PlasticStrainPolicy.cc at the pinned commit): after each stage's
update, with J2 = S:S/2 (RZ: the 2D double dot plus (S_xx + S_rr)^2 for the hoop element, which equals S_zz^2
because SolidSPHRZ keeps S traceless in 3D -- the checkpoints' trace is at most 6e-17 of |S|), f = min(1, Y /
sqrt(3 J2)), S <- f S, and eps_p += (1 - f) sqrt(3 J2) / (3 G): the radial return's equivalent plastic strain
increment (the trial von Mises stress's excess over Y is 3 G d eps_p). In uniaxial flow (plastic strain
diag(-a, a/2, a/2), incompressible) the equivalent strain sqrt(2/3 eps_p:eps_p) is a, the axial plastic shortening.
**The definition check** (`ps_kinematic`): in the middle third the axial strain eps_xx (the least-squares slope of
the axial displacement against the reference position) is the elastic, thermal and plastic strains together,
eps_xx = sigma_xx / E + alpha_eff dT - eps_p, so eps_p = alpha_eff dT + mean sigma_xx / E - eps_xx, against the
particles' mean eps_p (`ps_mean`).

**Measures** per run and frame k = 1 ... 10 (dT = 30 k K), over the middle third's particles (volume-weighted means
and spreads; RZ: rings): sigma_xx's mean, standard deviation, min and max, its error against the expectation
(-min(fe_ratio E alpha_eff dT, Y) -- the plateau is -Y whatever fe_ratio) over Y; the lateral stresses (3D yy and zz;
RZ the radial and hoop) mean and largest |.| over Y; eps_p's mean, sd, min, max, the expectation alpha_eff dT - Y/E,
the increment from the previous frame, the kinematic eps_p; over the held length (every non-held particle) eps_p's
mean (`ps_free_mean`); the internal KE ratio, the crossings monitor's estimates (intervals.csv), the held
particles' largest speed and displacement, the cost.

**Measured** (2026-10-10; the CSV; thresholds in tests/test_spheral_frag_runner_check5.py). The holds are exact
(0.0 in every row). fe_ratio 1.0540 (RZ 1.1 mm), 1.0660 and 1.0622 (3D 2.2, 1.1 mm), converged to 7e-4 between
the 1.5 and 1.0 mm meshes. Stable runs (RZ 1.1 mm, 3D 2.2 and 1.1 mm): the elastic middle (k 1) within -0.45 %,
-1.46 % and -0.39 % of Y of fe_ratio (-E alpha_eff dT) (the uniaxial value alone is 6-10 % short); the plateau
(k >= 3) at -Y within 1.4e-3 Y (spread 0.45 %, 3.1 %, 0.6 % of Y); the lateral stresses' mean within 1.7e-3 Y,
their largest |.| 1.8 %, 10.8 % (the 9-particle section's corners), 3.1 % of Y; at 1.1 mm the plastic-strain
increment per interval 0.86-1.04 times 7.02e-4 (k >= 5), the middle's eps_p at k 10 4.92e-3 (RZ) and 5.10e-3 (3D)
against 5.62e-3 (lag 7.1e-4 and 5.3e-4, the ends' early flow; 3D 2.2 mm 2.90e-3, lag 2.7e-3: the coarse bar's flow
concentrates 25-40 mm from the centre), and Spheral's eps_p equals the kinematic one to 4.7 %, 7.3 % (11 % at
2.2 mm) of the expectation; internal KE/E at most 1.2e-4; the monitor's estimates 1.2-3.8 % (no switch). Cost per
interval: RZ 1.1 mm 10 s (4 processes), 3D 2.2 mm 63 s and 3D 1.1 mm 282 s (18 processes; 29,484 particles, 84
ms/step, 3,251 steps), with the monitor's two reruns at 20 crossings on top. The RZ 2.2 mm rod -- five particle rows across its radius -- is unstable in
SolidSPHRZ under the axial compression: the axis rows slide against the outer ones from the second interval on
(KE/E 0.15-0.2 with 10 crossings, worse at 20 and 40; the same unscaled, `--mass-scaling off`, and the same with
`--strength elastic`, where the particles end interpenetrated at -11 GPa), so it is recorded and flagged
(`collapsed`), not held to the thresholds. RZ 1.1 mm (nine rows) and 3D at both spacings are stable.

**Running it.** Under `caffeinate -i -s` (the host's sleep suspends Docker's VM). Before each run on more than 4
processes the driver waits for the VM's load to fall to half of --max-load (as check 4).

Outputs: `data/spheral/m2_check5.csv` (one row per run and frame; committed), the finite-element reference
`spheral_output/m2/check5/fe.json`, the runs in `spheral_output/m2/check5/runs/`, plots in `spheral_output/m2/check5/`.
Resumable: bodies, the FE reference and runs whose outputs exist are skipped.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import os
import subprocess
import sys
import time

import numpy as np

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, REPO)
sys.path.insert(0, os.path.join(REPO, "analysis"))

from spheral_frag import body as body_stage, material, naming  # noqa: E402
from spheral_frag.runner import checkpoint, run  # noqa: E402

PREPARED = os.path.join("spheral_output", "prepare", "prep_model_d100.00mm_v07.50000kms_h077.500km_us76_sesam-table_"
                        "none_fem-physics_melt-girin_material-AA7075_scheil_deeprunoff-off_k00000-01254")
FRAMES = (0, 10)
OUT = os.path.join("spheral_output", "m2", "check5")
RUNS_OUT = os.path.join(OUT, "runs")
FE_JSON = os.path.join(OUT, "fe.json")
CSV = os.path.join("data", "spheral", "m2_check5.csv")
LAUNCHER = os.path.join("spheral_frag", "container", "spheral")
STRENGTH = "constant:100e6"
Y = 100e6                                   # Pa, the constant law's yield
DT_PER_FRAME = body_stage.UNIFORM_DELTA_K/body_stage.UNIFORM_RAMP_FRAMES     # 30 K per frame
HOLD = body_stage.BAR_HOLD_SPACINGS
FE_H_MM = (1.5, 1.0)                        # the finite-element reference's two mesh sizes (convergence)

FORMS = [("rz", 2.2, 1), ("rz", 1.1, 4), ("3d", 2.2, 18), ("3d", 1.1, 18)]     # form, dx mm, processes
COLLAPSED = {("rz", 2.2)}                  # measured unstable (module docstring): recorded, flagged

COLUMNS = ["form", "dx_mm", "k", "dT_K", "collapsed", "n_particles", "n_mid", "n_held", "processes", "scale",
           "fe_ratio", "sxx_expected_Pa", "sxx_mean_Pa", "sxx_sd_Pa", "sxx_min_Pa", "sxx_max_Pa", "sxx_err_over_Y",
           "lat_mean_over_Y", "lat_max_over_Y", "ps_expected", "ps_mean", "ps_sd", "ps_min", "ps_max",
           "ps_increment", "ps_increment_expected", "eps_xx_mid", "ps_kinematic", "ps_free_mean",
           "ke_internal_ratio", "monitored", "crossings_err_est", "crossings_err_pressure", "crossings_err_limit",
           "run_crossings", "held_max_speed_ms", "held_max_disp_m", "steps", "wall_s", "wall_per_step_ms", "run"]


def body_dir(form, dx):
    prep = os.path.basename(PREPARED)
    return naming.body_dir(os.path.join("spheral_output", "body"),
                           naming.body_name(prep, *FRAMES, form, dx, synthetic="uniform", geometry="bar"))


def run_name(form, dx):
    b = naming.canonical_brackets({"strength": STRENGTH, "sampler": "uniform", "loads": "off", "geometry": "bar"})
    return naming.run_name(os.path.basename(PREPARED), "replay", FRAMES, form, dx, b, 1)


def run_dir(form, dx):
    return naming.run_dir(os.path.join(REPO, RUNS_OUT), run_name(form, dx))


# ------------------------------------------------------------------------------------------------------- material
def constants(table):
    """(E, nu, alpha_eff): constant over 300-600 K in this table (checked), alpha_eff the free density's linear
    expansion referred to 300 K (module docstring)."""
    T = np.linspace(300.0, 600.0, 301)
    for f in ("young", "poisson", "bulk_modulus", "shear_modulus"):
        if np.ptp(getattr(table, f)(T)) != 0.0:
            raise RuntimeError(f"{f} is not constant between 300 and 600 K in this table")
    alpha = float(table.header["mechanical"]["entries"]["alpha"]["value"])
    a_eff = alpha/(1.0 + alpha*(300.0 - material.T_REF_EXPANSION))
    lin = (table.rho_free(300.0)/table.rho_free(T))**(1.0/3.0) - 1.0
    if float(np.max(np.abs(lin - a_eff*(T - 300.0)))) > 1e-12:
        raise RuntimeError("the table's free density is not the expected expansion law")
    if float(np.ptp(material.flow_stress(table, np.array([300.0, 600.0]), np.array([1e-3, 1e-3]), STRENGTH))) != 0.0:
        raise RuntimeError("the constant law is not constant")
    return float(table.young(300.0)), float(table.poisson(300.0)), a_eff


# ---------------------------------------------------------------------------------------- finite-element reference
def bar_mesh(form, length, width, h):
    """gmsh tetrahedra of the clamped bar: a box length x width x width (3d) or a cylinder of diameter `width`
    (rz), centred at the origin, axis x."""
    import gmsh
    import skfem as fe
    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.model.add("bar")
        if form == "3d":
            vol = gmsh.model.occ.addBox(-length/2, -width/2, -width/2, length, width, width)
        else:
            vol = gmsh.model.occ.addCylinder(-length/2, 0.0, 0.0, length, 0.0, 0.0, width/2)
        gmsh.model.occ.synchronize()
        gmsh.option.setNumber("Mesh.MeshSizeMin", h)
        gmsh.option.setNumber("Mesh.MeshSizeMax", h)
        for option in ("Mesh.MeshSizeFromPoints", "Mesh.MeshSizeFromCurvature", "Mesh.MeshSizeExtendFromBoundary"):
            gmsh.option.setNumber(option, 0)
        gmsh.option.setNumber("Mesh.Algorithm3D", 10)
        gmsh.model.mesh.generate(3)
        tags, coords, _ = gmsh.model.mesh.getNodes()
        _, _, enodes = gmsh.model.mesh.getElements(3, vol)
    finally:
        gmsh.finalize()
    tags = np.asarray(tags, dtype=np.int64)
    p = np.asarray(coords, dtype=float).reshape(-1, 3)
    t_tags = np.asarray(enodes[0], dtype=np.int64).reshape(-1, 4)
    used, inv = np.unique(t_tags, return_inverse=True)
    lookup = np.full(tags.max() + 1, -1, dtype=np.int64)
    lookup[tags] = np.arange(len(tags))
    return fe.MeshTet(p[lookup[used]].T.copy(), inv.reshape(-1, 4).T.copy())


def bar_fe(form, length, width, points, E, nu, h):
    """Linear thermoelasticity of the bar with both end faces clamped and a unit free strain (e_T = 1): sigma (n, 6)
    at `points` (xx, yy, zz, xy, yz, xz), per unit E e_T. P2, CG with pyamg."""
    import scipy.sparse.linalg as spla
    import skfem as fe
    from skfem.helpers import ddot, sym_grad
    from skfem.models.elasticity import lame_parameters, linear_elasticity
    import spheral_m2_elastic_fe as fe6
    mesh = bar_mesh(form, length, width, h)
    basis = fe.Basis(mesh, fe.ElementVector(fe.ElementTetP2()))
    lam, mu = lame_parameters(1.0, nu)                       # E = 1: stresses per unit E
    K = fe.asm(linear_elasticity(lam, mu), basis)
    three_k = 3.0*lam + 2.0*mu                               # 3 K: the thermal stress 3 K e_T I

    @fe.LinearForm
    def thermal(v, w):
        return three_k*(sym_grad(v)[0, 0] + sym_grad(v)[1, 1] + sym_grad(v)[2, 2])

    f = fe.asm(thermal, basis)
    clamp = basis.get_dofs(lambda x: np.abs(np.abs(x[0]) - length/2) < 1e-9*length).all()
    I = np.setdiff1d(np.arange(basis.N), clamp)
    KI = K[I][:, I].tocsr()
    import pyamg
    ml = pyamg.smoothed_aggregation_solver(KI, symmetry="hermitian", max_coarse=500)
    u = np.zeros(basis.N)
    uI, info = spla.cg(KI, f[I], M=ml.aspreconditioner(cycle="V"), rtol=1e-12, atol=0.0, maxiter=5000)
    if info != 0:
        raise RuntimeError("CG did not converge")
    u[I] = uI
    s = fe6.stress_at(basis, u, points, lam, mu)
    s[:, :3] -= three_k                                      # sigma = C:(eps - e_T I)
    return s, basis.N


def middle_points(form, dx):
    """The body's middle-third reference points as 3D points (RZ at azimuth 0: (x, r, 0)), and the volumes."""
    x_ref, x_sample, info = body_stage.bar_points(form, dx, body_stage.BAR_LENGTH_MM*1e-3,
                                                  body_stage.BAR_WIDTH_MM*1e-3)
    V = body_stage.reference_volumes(form, x_ref, dx, 1.0)
    mid = np.abs(x_ref[:, 0]) < info["length_m"]/6.0
    return x_sample[mid], V[mid], info


def fe_reference(E, nu, force=False):
    """fe_ratio per form and spacing: the middle third's volume-weighted mean sigma_xx / (-E e_T) on the finer mesh,
    with the coarser mesh's value and the lateral stresses (fe.json)."""
    path = os.path.join(REPO, FE_JSON)
    if os.path.isfile(path) and not force:
        with open(path) as fh:
            return json.load(fh)
    out = {}
    for form, dx_mm, _ in FORMS:
        dx = dx_mm*1e-3
        pts, V, info = middle_points(form, dx)
        L = info["length_m"] - 2*HOLD*dx                    # between the held layers
        W = info["width_m"] if form == "3d" else 2.0*info["radius_m"]
        # the clamped faces are at |x| = L/2, half a spacing outside the outermost free particles; every point
        # (RZ rows at r = (j + 1/2) dx < n_r dx) lies inside the mesh
        res = {}
        for h in FE_H_MM:
            t = time.perf_counter()
            s, ndof = bar_fe(form, L, W, pts, E, nu, h*1e-3)
            ratio = -float(np.sum(V*s[:, 0])/np.sum(V))
            lat = float(np.max(np.abs(s[:, 1:3])))
            res[h] = dict(ratio=ratio, lateral_max=lat, ndof=int(ndof), seconds=time.perf_counter() - t)
            print("  FE {} {} mm h {} mm: ratio {:.5f}, lateral max {:.2e} E e_T, {} dofs, {:.1f} s".format(
                form, dx_mm, h, ratio, lat, ndof, res[h]["seconds"]), flush=True)
        out["{}_{}".format(form, dx_mm)] = dict(form=form, dx_mm=dx_mm, length_m=L, width_m=W,
                                                 ratio=res[FE_H_MM[-1]]["ratio"],
                                                 ratio_coarse=res[FE_H_MM[0]]["ratio"],
                                                 lateral_max=res[FE_H_MM[-1]]["lateral_max"],
                                                 meshes={str(h): r for h, r in res.items()})
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path + ".tmp", "w") as fh:
        json.dump(out, fh, indent=2)
    os.replace(path + ".tmp", path)
    return out


# ------------------------------------------------------------------------------------------------------- running
def build_bodies():
    for form, dx, _ in FORMS:
        if os.path.isfile(os.path.join(REPO, body_dir(form, dx), "body.json")):
            continue
        cmd = [sys.executable, "-m", "spheral_frag", "body", "--prepared", PREPARED, "--frames", "{}-{}".format(*FRAMES),
               "--form", form, "--dx-mm", str(dx), "--synthetic", "uniform", "--geometry", "bar", "--quiet"]
        subprocess.run(cmd, cwd=REPO, check=True)


def launch(form, dx, procs, max_load, force):
    import spheral_m2_check4 as c4
    d = run_dir(form, dx)
    if os.path.isfile(os.path.join(d, "summary.json")):
        return "done"
    if procs > 4 and not force:
        c4.wait_for_load(0.5*max_load)
    cmd = [LAUNCHER] + (["-n", str(procs)] if procs > 1 else []) + [
        "spheral_frag/runner/run.py", "replay", "--body", body_dir(form, dx), "--prepared", PREPARED,
        "--frames", "{}-{}".format(*FRAMES), "--strength", STRENGTH, "--loads", "off", "--ke-fraction", "1e300",
        "--out", RUNS_OUT, "--max-load", repr(max_load)]
    if os.path.isdir(d):
        cmd.append("--resume")
    if force:
        cmd.append("--force")
    t = time.perf_counter()
    r = subprocess.run(cmd, cwd=REPO, capture_output=True, text=True)
    os.makedirs(os.path.join(REPO, OUT, "logs"), exist_ok=True)
    with open(os.path.join(REPO, OUT, "logs", "{}_dx{}.log".format(form, dx)), "w") as fh:
        fh.write(" ".join(cmd) + "\n" + r.stdout + "\n--- stderr\n" + r.stderr)
    if r.returncode != 0:
        print("  FAILED ({}): {}".format(r.returncode, r.stderr.strip().splitlines()[-1:] if r.stderr else ""),
              flush=True)
        return "failed"
    return "ran in {:.0f} s".format(time.perf_counter() - t)


# ----------------------------------------------------------------------------------------------------- measuring
def wmean(a, w):
    return float(np.sum(w*a)/np.sum(w))


def wsd(a, w):
    m = wmean(a, w)
    return math.sqrt(float(np.sum(w*(a - m)**2)/np.sum(w)))


def measure(d, k, table, E, alpha, dx, form):
    """The middle third's measures at frame k of run directory d (module docstring)."""
    ck = checkpoint.read(os.path.join(d, "checkpoints", checkpoint.filename(k)))
    st = run.stress_state(ck, table)
    x_ref, x = ck["x_ref"], ck["x"]
    half = float(np.max(x_ref[:, 0])) + 0.5*dx
    mid = np.abs(x_ref[:, 0]) < half/3.0
    free = np.abs(x_ref[:, 0]) < half - HOLD*dx
    V = st["V"]
    w = V[mid]
    s = st["sigma"][mid]
    sxx = s[:, 0, 0]
    lat = np.concatenate([s[:, 1, 1], s[:, 2, 2]])
    ps = ck["plastic_strain"]
    xm = x_ref[mid, 0]
    ux = x[mid, 0] - xm
    A = np.stack([xm, np.ones_like(xm)], axis=1)*np.sqrt(w)[:, None]
    slope = float(np.linalg.lstsq(A, ux*np.sqrt(w), rcond=None)[0][0])
    dT = DT_PER_FRAME*k
    return dict(n_particles=len(ck["pid"]), n_mid=int(mid.sum()), sxx_mean=wmean(sxx, w), sxx_sd=wsd(sxx, w),
                sxx_min=float(sxx.min()), sxx_max=float(sxx.max()),
                lat_mean=float(np.mean([wmean(s[:, 1, 1], w), wmean(s[:, 2, 2], w)])),
                lat_max=float(np.abs(lat).max()), ps_mean=wmean(ps[mid], w), ps_sd=wsd(ps[mid], w),
                ps_min=float(ps[mid].min()), ps_max=float(ps[mid].max()), eps_xx=slope,
                ps_kin=alpha*dT + wmean(sxx, w)/E - slope, ps_free=wmean(ps[free], V[free]))


def analyse():
    table = material.MaterialTable.load(os.path.join(REPO, PREPARED))
    E, nu, alpha = constants(table)
    fe_ref = fe_reference(E, nu)
    print("E {:.4g} Pa, nu {}, alpha_eff {:.6e} 1/K, dT_y {:.3f} K".format(E, nu, alpha, Y/(E*alpha)))
    rows = []
    for form, dx_mm, procs in FORMS:
        d = run_dir(form, dx_mm)
        if not os.path.isfile(os.path.join(d, "summary.json")):
            print("  missing: {}".format(naming.run_name_from_dir(d)))
            continue
        with open(os.path.join(d, "summary.json")) as fh:
            summ = json.load(fh)
        iv = {int(r["k1"]): r for r in run.read_intervals(os.path.join(d, "intervals.csv"))}
        ratio = fe_ref["{}_{}".format(form, dx_mm)]["ratio"]
        prev = 0.0
        for k in range(FRAMES[0] + 1, FRAMES[1] + 1):
            m = measure(d, k, table, E, alpha, dx_mm*1e-3, form)
            r = iv[k]
            dT = DT_PER_FRAME*k
            sxx_exp = -min(ratio*E*alpha*dT, Y)
            ps_exp = max(alpha*dT - Y/E, 0.0)

            def num(c):
                v = r.get(c, "")
                return float(v) if v not in ("", None) else math.nan
            rows.append(dict(form=form, dx_mm=dx_mm, k=k, dT_K=dT, collapsed=int((form, dx_mm) in COLLAPSED),
                             n_particles=m["n_particles"], n_mid=m["n_mid"], n_held=int(num("n_held")),
                             processes=summ["processes"], scale=num("scale"), fe_ratio=ratio,
                             sxx_expected_Pa=sxx_exp, sxx_mean_Pa=m["sxx_mean"], sxx_sd_Pa=m["sxx_sd"],
                             sxx_min_Pa=m["sxx_min"], sxx_max_Pa=m["sxx_max"],
                             sxx_err_over_Y=(m["sxx_mean"] - sxx_exp)/Y, lat_mean_over_Y=m["lat_mean"]/Y,
                             lat_max_over_Y=m["lat_max"]/Y, ps_expected=ps_exp, ps_mean=m["ps_mean"],
                             ps_sd=m["ps_sd"], ps_min=m["ps_min"], ps_max=m["ps_max"],
                             ps_increment=m["ps_mean"] - prev, ps_increment_expected=(
                                 max(alpha*dT - Y/E, 0.0) - max(alpha*(dT - DT_PER_FRAME) - Y/E, 0.0)),
                             eps_xx_mid=m["eps_xx"], ps_kinematic=m["ps_kin"], ps_free_mean=m["ps_free"],
                             ke_internal_ratio=num("ke_internal_ratio"), monitored=int(num("monitored")),
                             crossings_err_est=num("crossings_err_est"),
                             crossings_err_pressure=num("crossings_err_pressure"),
                             crossings_err_limit=num("crossings_err_limit"), run_crossings=int(num("run_crossings")),
                             held_max_speed_ms=num("held_max_speed_ms"), held_max_disp_m=num("held_max_disp_m"),
                             steps=int(num("steps")), wall_s=num("wall_s"),
                             wall_per_step_ms=1e3*num("wall_per_step_s"), run=naming.run_name_from_dir(d)))
            prev = m["ps_mean"]
    path = os.path.join(REPO, CSV)
    with open(path, "w", newline="") as fh:
        wr = csv.writer(fh)
        wr.writerow(COLUMNS)
        for row in rows:
            wr.writerow([repr(float(row[c])) if isinstance(row[c], (float, np.floating)) else row[c]
                         for c in COLUMNS])
    print("wrote {} ({} rows)".format(CSV, len(rows)))
    report(rows)
    plots(rows, E, alpha)
    return rows


def report(rows):
    print("\n{:3s} {:>4s} {:>2s} {:>4s} {:>8s} {:>8s} {:>7s} {:>8s} {:>9s} {:>9s} {:>9s} {:>9s} {:>9s} {:>9s}".format(
        "frm", "dx", "k", "dT", "sxx_exp", "sxx", "sd", "err/Y", "lat/Y", "ps_exp", "ps", "dps", "ps_kin",
        "KEint/E"))
    for r in rows:
        print("{:3s} {:4.1f} {:2d} {:4.0f} {:8.2f} {:8.2f} {:7.2f} {:+8.4f} {:9.2e} {:9.2e} {:9.2e} {:9.2e} {:9.2e} "
              "{:9.2e}".format(r["form"], r["dx_mm"], r["k"], r["dT_K"], r["sxx_expected_Pa"]/1e6,
                               r["sxx_mean_Pa"]/1e6, r["sxx_sd_Pa"]/1e6, r["sxx_err_over_Y"], r["lat_max_over_Y"],
                               r["ps_expected"], r["ps_mean"], r["ps_increment"], r["ps_kinematic"],
                               r["ke_internal_ratio"]))


INK, MUTED = "#1f2328", "#6e7781"
SHADES = {("rz", 2.2): "#c6dbef", ("rz", 1.1): "#2171b5", ("3d", 2.2): "#a1d99b", ("3d", 1.1): "#238b45"}
MARKS = {("rz", 2.2): "o", ("rz", 1.1): "s", ("3d", 2.2): "^", ("3d", 1.1): "D"}


def plots(rows, E, alpha):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    dT = np.linspace(0.0, 300.0, 301)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    axes[0].plot(dT, -np.minimum(E*alpha*dT, Y)/1e6, color=INK, lw=1.5, label="uniaxial: -min(E alpha dT, Y)")
    axes[1].plot(dT, np.maximum(alpha*dT - Y/E, 0.0), color=INK, lw=1.5, label="alpha dT - Y/E")
    for key in SHADES:
        rr = [r for r in rows if (r["form"], r["dx_mm"]) == key]
        if not rr:
            continue
        lab = "{} {} mm{}".format(key[0].upper(), key[1], " (collapsed)" if key in COLLAPSED else "")
        t = [r["dT_K"] for r in rr]
        axes[0].errorbar(t, [r["sxx_mean_Pa"]/1e6 for r in rr], yerr=[r["sxx_sd_Pa"]/1e6 for r in rr],
                         color=SHADES[key], marker=MARKS[key], ms=6, lw=1.5, capsize=2, label=lab)
        axes[1].plot(t, [r["ps_mean"] for r in rr], color=SHADES[key], marker=MARKS[key], ms=6, lw=1.5, label=lab)
    axes[0].set_xlabel("dT (K)")
    axes[0].set_ylabel("middle third: sigma_xx (MPa)")
    axes[1].set_xlabel("dT (K)")
    axes[1].set_ylabel("middle third: equivalent plastic strain")
    for ax in axes:
        ax.grid(color="#d0d7de", lw=0.5)
        ax.legend(fontsize=8, frameon=False)
    axes[1].set_ylim(-5e-4, 8e-3)
    fig.suptitle("Check 5: the heated bar held at both ends (Y = 100 MPa)", fontsize=11)
    fig.tight_layout()
    os.makedirs(os.path.join(REPO, OUT), exist_ok=True)
    fig.savefig(os.path.join(REPO, OUT, "bar.png"), dpi=130)
    plt.close(fig)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--only", choices=("fe", "run", "analyse"))
    p.add_argument("--max-load", type=float, default=2.0)
    p.add_argument("--force", action="store_true", help="start runs even above --max-load")
    a = p.parse_args(argv)
    os.chdir(REPO)
    if a.only in (None, "fe"):
        table = material.MaterialTable.load(os.path.join(REPO, PREPARED))
        E, nu, _ = constants(table)
        fe_reference(E, nu)
    if a.only in (None, "run"):
        build_bodies()
        for form, dx, procs in FORMS:
            print("{} {} mm: ".format(form, dx), end="", flush=True)
            print(launch(form, dx, procs, a.max_load, a.force), flush=True)
    if a.only in (None, "analyse"):
        analyse()
    return 0


if __name__ == "__main__":
    sys.exit(main())

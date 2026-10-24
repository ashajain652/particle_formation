"""Spheral-side driver of tests/test_spheral_frag_runner_packages.py (Spheral M2 plan, Task 4). Run through the launcher
from the repository:

    spheral_frag/container/spheral [-n N] tests/fixtures/spheral_frag/runner_packages_check.py MODE --form F ...

Every mode builds a body from an input .npz (x_ref, V, pid, dx; the lattice sphere of the test) with the thermal
material (PyThermalEOS + PyThermalStrength, law "elastic" unless --law), SolidSPH (S.SPH on the solid node list) and
CheapSynchronousRK2 as runner_material_check.py does, and writes a plain .npz from rank 0.

  temperature  --input IN --scale S --out OUT
               IN also holds T_table (n, n_frames) and times; the body starts at T_table[:, 0]; ImposedTemperature
               (monitor on) is advanced to every frame time in turn (the integrator's maxTime is the next frame
               time); after every step the gathered imposed temperature, u and the time.
  breathing    --input IN --zeta Z --steps N --eps E --out OUT
               the body at 300 K, s = 1, with the velocity field eps x (a breathing start), Damping(Z, omega_1); per
               step t, q = sum m v.x / sum m |x|^2 (the breathing mode's rate), the kinetic energy.
  loads        --input IN --tables uniform|fixture [--pressure P] [--fixture F --frames K0,K1] [--area A]
               [--exposure E] [--threshold G] [--decel] [--scale S] [--steps N] [--interval T] --out OUT
               SurfaceLoads (and Deceleration with --decel) on the body at 300 K; the uniform table is P at every
               inclination, no shear, base pressure P; the fixture tables are the frames' rows of
               tests/fixtures/spheral_frag/loads/loads_mvp.npz held over [0, T]. Writes the surface table, the
               sums, then N steps and the final momentum, velocities and the deceleration history.
  exposure     --input IN --out OUT
               surface detection and both exposure rules on IN's points (x, V, dx), without a Spheral body.
  cost         --body DIR --table TAB --fixture F --frames K0,K1 --steps N --warmup W [--max-load L] --out OUT
               a body stage's body at its frame-0 temperatures and the replay's scale for 10 crossings per 0.5 s;
               per configuration (bare; + ImposedTemperature; + Damping; + SurfaceLoads; + SurfaceLoads and
               Deceleration; all four; bare again) a fresh integrator, W warm-up and N timed steps (barrier to barrier), and
               the time spent inside each package's hooks (largest over the ranks), plus the packages'
               construction times."""
import argparse
import math
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "..")))
from spheral_frag import loads as core_loads, material  # noqa: E402
from spheral_frag.runner import nodes as rn  # noqa: E402
from spheral_frag.runner.eos import make_eos  # noqa: E402
from spheral_frag.runner.packages import damping as pd, get_vector, loads as pl, temperature as pt  # noqa: E402
from spheral_frag.runner.strength import make_strength  # noqa: E402

CFL = 0.25


def comm():
    from mpi4py import MPI
    return MPI.COMM_WORLD


def load(path):
    with np.load(path, allow_pickle=False) as z:
        return {k: z[k] for k in z.files}


def build(S, form, d, table, T0, scale, law):
    eos = make_eos(S, table, scale)
    strength = make_strength(S, eos, law)
    rho = scale*table.rho_free(T0)
    arrays = dict(x_ref=d["x_ref"], pid=d["pid"], rho=rho, m=rho*d["V"], u=table.enthalpy(T0))
    if "v" in d:
        arrays["v"] = d["v"]
    return rn.make_body(S, form, arrays, eos, strength, dx=float(d["dx"]))


def integrator_for(S, body, packages, t0=0.0, dt0=None):
    from SpheralController import SpheralController
    hydro = S.SPH(dataBase=body.db, W=body.kernel, cfl=CFL, compatibleEnergyEvolution=True,
                  densityUpdate=S.IntegrateDensity, XSPH=False, ASPH=False)
    integ = S.CheapSynchronousRK2Integrator(body.db)
    integ.appendPhysicsPackage(hydro)
    for p in packages:
        integ.appendPhysicsPackage(p)
    integ.dtMin, integ.dtMax, integ.verbose = 1e-15, 1.0, False
    ctl = SpheralController(integ, body.kernel, statsStep=10**9, printStep=10**9, redistributeStep=None,
                            restartStep=None, vizStep=None, vizTime=None, SPH=True, iterateInitialH=False)
    integ.currentTime = t0
    if dt0 is not None:
        integ.lastDt = dt0
    rn.prime_derivatives(body, integ)
    return hydro, integ, ctl


def step(integ, ctl, t_max=1e40):
    integ.step(t_max)
    ctl.totalSteps += 1
    ctl.doPeriodicWork()


def first_dt(table, scale, dx):
    return 0.1*CFL*dx/float(pd.longitudinal_speed(table, scale))


def momentum(body):
    m = np.array(body.nodes.mass().internalValues())
    v = get_vector(body.nodes.velocity(), body.dim)
    from spheral_frag.runner.packages import allsum
    return allsum((m[:, None]*v).sum(axis=0)), allsum(m.sum())


# ---------------------------------------------------------------------------------------------------------------

def temperature(S, a, table):
    d = load(a.input)
    Tt, times = d["T_table"], d["times"]
    body = build(S, a.form, d, table, Tt[:, 0], a.scale, a.law)
    hist = pt.TemperatureHistory(Tt, times, d["pid"])
    imp = pt.make_imposed_temperature(S, body, hist, table, float(times[0]), monitor=True)
    hydro, integ, ctl = integrator_for(S, body, [imp], float(times[0]), first_dt(table, a.scale, float(d["dx"])))
    rows = []

    def record():
        g = rn.gather_by_id(body.nodes, dict(pid=body.pid, T=imp.temperature, u=body.nodes.specificThermalEnergy()),
                            body.dim)
        if comm().rank == 0:
            rows.append((integ.currentTime, g["T"], g["u"], g["pid"]))

    record()
    for tf in times[1:]:
        while integ.currentTime < tf*(1.0 - 1e-15):
            step(integ, ctl, float(tf))
            record()
    n_eval, dmax, cmax = imp.monitor_summary()
    if comm().rank == 0:
        np.savez(a.out, t=np.array([r[0] for r in rows]), T=np.array([r[1] for r in rows]),
                 u=np.array([r[2] for r in rows]), pid=rows[0][3], times=times, T_table=Tt,
                 monitor=np.array([n_eval, dmax, cmax]),
                 monitor_rows=np.array(imp.monitor_rows), steps=ctl.totalSteps)
    return 0


def breathing(S, a, table):
    d = load(a.input)
    d["v"] = a.eps*d["x_ref"]
    n = len(d["pid"])
    body = build(S, a.form, d, table, np.full(n, 300.0), 1.0, a.law)
    D = pd.axial_extent(body)
    omega1 = pd.breathing_omega(table, 1.0, D)
    damp = pd.make_damping(S, body, a.zeta, omega1)
    hydro, integ, ctl = integrator_for(S, body, [damp], 0.0, first_dt(table, 1.0, float(d["dx"])))
    from spheral_frag.runner.packages import allsum
    rows = []

    def record():
        x = get_vector(body.nodes.positions(), body.dim)
        v = get_vector(body.nodes.velocity(), body.dim)
        m = np.array(body.nodes.mass().internalValues())
        s = allsum(np.array([np.sum(m*np.einsum("ij,ij->i", v, x)), np.sum(m*np.einsum("ij,ij->i", x, x)),
                             0.5*np.sum(m*np.einsum("ij,ij->i", v, v))]))
        rows.append((integ.currentTime, s[0]/s[1], s[2]))

    record()
    for _ in range(a.steps):
        step(integ, ctl)
        record()
    if comm().rank == 0:
        np.savez(a.out, rows=np.array(rows), omega1=omega1, D=D, zeta=a.zeta, eps=a.eps,
                 max_fraction=damp.max_fraction, c_l=float(pd.longitudinal_speed(table, 1.0)))
    return 0


def uniform_table(P):
    edges = core_loads.THETA_EDGES_DEG
    c = 0.5*(edges[1:] + edges[:-1])
    return core_loads.LoadTable(theta_deg=c, p=np.full(len(c), P), tau=np.zeros(len(c)), area=np.ones(len(c)),
                                theta_last_deg=float(c[-1]), p_stag=P/core_loads.BASE_FRACTION)


def fixture_tables(path, ks):
    z = load(path)
    rows = [int(np.flatnonzero(z["k"] == k)[0]) for k in ks]
    sub = {key: (z[key][rows] if key not in ("theta_deg", "edges_deg") else z[key]) for key in
           ("theta_deg", "edges_deg", "p", "tau", "area", "lee", "theta_last_deg", "p_stag")}
    return core_loads.unstack(sub), z["surface_area_m2"][rows], z["time_s"][rows]


def loads(S, a, table):
    d = load(a.input)
    n = len(d["pid"])
    body = build(S, a.form, d, table, np.full(n, 300.0), a.scale, a.law)
    if a.tables == "uniform":
        tables = [uniform_table(a.pressure)]
    else:
        tables, _, _ = fixture_tables(a.fixture, [int(k) for k in a.frames.split(",")])
    t0 = time.perf_counter()
    sl = pl.make_surface_loads(S, body, tables, (0.0, a.interval), a.area, threshold=a.threshold,
                               exposure=a.exposure)
    build_s = time.perf_counter() - t0
    packages = [sl]
    dec = None
    if a.decel:
        dec = pl.make_deceleration(S, sl, a.scale)
        packages.append(dec)
    st = sl.surface_table()
    # the packages' own balance at a time inside the interval: sum m a_loads against sum F, and with the frame's
    # deceleration sum m (a_loads - a_frame) (axial only in RZ: the rings' radial sums are not net forces)
    from spheral_frag.runner.packages import allsum
    tb = 0.3*a.interval
    ma = allsum((sl.local_mass[:, None]*sl.local_acceleration(tb)).sum(axis=0))
    F_t = sl.total_force(tb)
    a_frame = F_t/sl.mass
    bal_decel = allsum((sl.local_mass[:, None]*(sl.local_acceleration(tb) - a_frame)).sum(axis=0))
    hydro, integ, ctl = integrator_for(S, body, packages, 0.0, first_dt(table, a.scale, float(d["dx"])))
    p0, M = momentum(body)
    for _ in range(a.steps):
        step(integ, ctl, a.interval)
    p1, _ = momentum(body)
    g = rn.gather_by_id(body.nodes, dict(pid=body.pid, v=body.nodes.velocity(), x=body.nodes.positions(),
                                         m=body.nodes.mass()), body.dim)
    disp = sl.max_displacement()
    if comm().rank == 0:
        hist = dec.history() if dec is not None else np.zeros((0, 1 + 2*body.dim))
        np.savez(a.out, **{"st." + k: v for k, v in st.items()}, force_k=sl.force_k, mass=sl.mass,
                 share_factor=sl.share_factor, area_sum=sl.area_sum, n_surface=sl.n_surface,
                 n_exposed=sl.n_exposed, p0=p0, p1=p1, M=M, t=integ.currentTime, steps=ctl.totalSteps,
                 v=g["v"], x=g["x"], m=g["m"], pid=g["pid"], decel_history=hist, build_s=build_s,
                 max_displacement_dx=disp, processes=comm().size, sum_ma=ma, F_t=F_t, bal_decel=bal_decel)
    return 0


def exposure(S, a, table):
    """Surface detection and both exposure rules on IN's points (x, V, dx) alone, no Spheral body."""
    d = load(a.input)
    dx = float(d["dx"])
    x = d["x"]
    g, surf, normal = pl.detect_surface(a.form, x, d["V"], x, rn.NPERH*dx, a.threshold)
    out = dict(grad=g, surface=surf, normal=normal)
    for mode in pl.EXPOSURES:
        out["exposed." + mode] = pl.exposure(a.form, x, surf, dx, mode)
    np.savez(a.out, **out)
    return 0


def timed(cls, acc, key):
    """A subclass of a package class that adds the wall time of its hooks to acc[key]."""
    class Timed(cls):
        def evaluateDerivatives(self, *args):
            t0 = time.perf_counter()
            cls.evaluateDerivatives(self, *args)
            acc[key] += time.perf_counter() - t0

        def finalize(self, *args):
            t0 = time.perf_counter()
            if hasattr(cls, "finalize"):
                cls.finalize(self, *args)
            acc[key] += time.perf_counter() - t0
    return Timed


def cost(S, a, table):
    c = comm()
    load_avg = c.bcast(os.getloadavg() if c.rank == 0 else None, root=0)
    if load_avg[0] > a.max_load:
        if c.rank == 0:
            print(f"cost: load average {load_avg[0]:.2f} > {a.max_load}; not starting")
        return 2
    hist, b = pt.load_history(a.body)
    dx = float(b["dx_m"])
    T0 = np.asarray(hist.T[:, 0])
    form = a.form
    # the replay's scale for 10 crossings of the body's axial extent per 0.5 s (plan Task 5's formula)
    D = float(np.ptp(b["x_ref"][:, 0])) + dx
    c_l = float(pd.longitudinal_speed(table, 1.0))
    scale = (c_l*0.5/(10.0*D))**2
    d = dict(x_ref=b["x_ref"], pid=b["pid"], V=b["m"]/b["rho"], dx=np.array(dx))
    body = build(S, form, d, table, T0, scale, "flight")
    tables, areas, times = fixture_tables(a.fixture, [int(k) for k in a.frames.split(",")])
    omega1 = pd.breathing_omega(table, scale, pd.axial_extent(body))
    acc = {}
    construct = {}
    configs = [("bare", ()), ("temperature", ("T",)), ("damping", ("D",)), ("loads", ("L",)),
               ("loads+decel", ("L", "A")), ("all", ("T", "D", "L", "A")), ("bare_again", ())]
    out = {}
    for name, which in configs:
        for k in which:
            acc[k] = 0.0
        pk = []
        t0 = time.perf_counter()
        if "T" in which:
            pk.append(timed(pt.temperature_class(S), acc, "T")(body, hist, table, float(times[0])))
        if "D" in which:
            pk.append(timed(pd.damping_class(S), acc, "D")(body, pd.ZETA, omega1))
        sl = None
        if "L" in which:
            t1 = time.perf_counter()
            sl = timed(pl.loads_classes(S)[0], acc, "L")(body, tables, times, float(areas[0]))
            construct["loads"] = time.perf_counter() - t1
        if "A" in which:
            pk.append(timed(pl.loads_classes(S)[1], acc, "A")(sl, scale))
        if sl is not None:
            pk.insert(len(pk) - (1 if "A" in which else 0), sl)
        construct[name] = time.perf_counter() - t0
        hydro, integ, ctl = integrator_for(S, body, pk, float(times[0]), first_dt(table, scale, dx))
        for _ in range(a.warmup):
            step(integ, ctl, float(times[1]))
        for k in which:
            acc[k] = 0.0
        walls = []
        for _ in range(a.steps):
            c.Barrier()
            t0 = time.perf_counter()
            step(integ, ctl, float(times[1]))
            c.Barrier()
            walls.append(time.perf_counter() - t0)
        inside = {k: c.allreduce(acc[k], op=__import__("mpi4py").MPI.MAX)/a.steps for k in which}
        out[name] = (np.array(walls), inside)
        del integ, ctl, hydro
    n = c.allreduce(body.nodes.numInternalNodes)
    if c.rank == 0:
        save = dict(n=n, processes=c.size, load=np.array(load_avg), scale=scale, omega1=omega1)
        for name, (w, inside) in out.items():
            save[f"walls.{name}"] = w
            for k, v in inside.items():
                save[f"inside.{name}.{k}"] = v
            print(f"cost {form} {name:12s}: n={n} procs={c.size} median {np.median(w)*1e3:7.2f} ms  "
                  + "  ".join(f"{k} {v*1e3:.3f} ms" for k, v in inside.items()))
        for k, v in construct.items():
            save[f"construct.{k}"] = v
        print("construct", {k: round(v, 3) for k, v in construct.items()}, "load", load_avg)
        np.savez(a.out, **save)
    return 0


def main(argv):
    p = argparse.ArgumentParser()
    p.add_argument("mode", choices=("temperature", "breathing", "loads", "exposure", "cost"))
    p.add_argument("--form", choices=rn.FORMS, required=True)
    p.add_argument("--table", required=True)
    p.add_argument("--input")
    p.add_argument("--body")
    p.add_argument("--law", default="elastic")
    p.add_argument("--scale", type=float, default=1.0)
    p.add_argument("--zeta", type=float, default=pd.ZETA)
    p.add_argument("--eps", type=float, default=1.0)
    p.add_argument("--steps", type=int, default=50)
    p.add_argument("--warmup", type=int, default=5)
    p.add_argument("--tables", choices=("uniform", "fixture"), default="uniform")
    p.add_argument("--pressure", type=float, default=1.0e4)
    p.add_argument("--fixture", default="tests/fixtures/spheral_frag/loads/loads_mvp.npz")
    p.add_argument("--frames", default="182,183")
    p.add_argument("--area", type=float, default=4.0*math.pi*0.05**2)
    p.add_argument("--exposure", choices=pl.EXPOSURES, default="first_body")
    p.add_argument("--threshold", type=float, default=pl.THRESHOLD)
    p.add_argument("--decel", action="store_true")
    p.add_argument("--interval", type=float, default=0.5)
    p.add_argument("--max-load", type=float, default=2.0)
    p.add_argument("--out", required=True)
    a = p.parse_args(argv)
    S, _ = rn.spheral_module(a.form)
    table = material.MaterialTable.load(a.table)
    return {"temperature": temperature, "breathing": breathing, "loads": loads, "exposure": exposure,
            "cost": cost}[a.mode](S, a, table)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))

"""Cost of a SolidSPH step on a 100 mm AA7075 sphere, 3D or axisymmetric (Spheral M0 plan, Task 5).

The body is a sphere of radius --radius-mm (50) at rest on a lattice with a node at the centre, clipped to the
sphere: in 3D the full ball; in RZ (Spheral's axisymmetric form, x = z along the axis, y = r >= 0) the half-disk
of the meridional plane, with the radial rows at (j + 1/2) dx so that no particle sits on the axis. The node at the
centre puts the counts within 1.7 % of V/dx^3 (3D) and (pi R^2/2)/dx^2 (RZ) at every spacing the plan uses; a
cell-centred lattice is 5.4 % off at the 10 mm test spacing.

The material is for timing only. Asha chose (2026-10-05) to time SteinbergGuinanStrength, whose cost is what
production pays; the plan paired it with the Gruneisen equation of state, "both with the aluminium parameters of
Spheral's material library". Two things changed that (Asha, 2026-10-06):

  - Equation of state: Murnaghan with exponent 1, P = K (rho/rho0 - 1), K = E/(3(1 - 2 nu)) = 70.3 GPa for AA7075
    (E 71.7 GPa, nu 0.33), rho0 2,813 kg/m^3, atomic weight 26.98. Not Gruneisen: at 116c71f
    GruneisenEquationOfState::pressureAndDerivs returns dP/drho = rho0 C0^2 instead of C0^2 for rho <= rho0
    (introduced in c77981253, 2023-08-29; reaching the sound speed and bulk modulus through computeDPDrho since
    18a8420d1, 2023-09-26; still in develop on 2026-10-06). In SI its sound speed at rest and in tension is
    sqrt(rho0) = 53 times C0 and the damage model's longitudinal sound speed 44 times too high, so dt would be 44-53
    times too small. Murnaghan is checked correct (sound speed sqrt(K/rho0) on both sides of rho0, SI and CGS) and
    is the equation of state production uses on route 2 (spec section 8.2). The equation of state is a small part of
    a step's cost.
  - Steinberg-Guinan: Spheral's library has none for aluminium, so 6061-T6 from Steinberg, Cochran & Guinan (J. Appl. Phys. 51, 1498, 1980):
    G0 27.6 GPa, A = G'_p/G0 = 6.52 /Mbar, B = G'_T/G0 = 6.16e-4 /K, Y0 0.29 GPa, Ymax 0.68 GPa, beta 125, n 0.10,
    initial plastic strain 0; Yp 1e-3 as in every Spheral example. The cold-energy fit is zero and the melt-energy
    fit a constant rho0 cv (Tm0 - 300 K) with Tm0 1220 K and the Dulong-Petit cv of 26.98 g/mol: SCG's melt curve
    rises with compression, but a ninth-order polynomial costs the same whatever its coefficients. A, B, beta and n
    are non-zero, so Spheral takes the full pressure-, temperature- and hardening-dependent branch, not the
    constant-modulus shortcut.
  - Damage (--damage on): ProbabilisticDamageModel with the library's aluminium Weibull k and m (converted to SI by
    the factory), its documented defaults otherwise, and a fixed seed. In RZ the flaw count per particle follows
    Spheral's RZ node volume; that matters for fragmentation, not for cost.

summary.json says "material": "timing only".

Hydro: SolidSPH (the SPH factory on a solid node list), WendlandC4 kernel, nPerh 2.01, IntegrateDensity,
compatible energy, no XSPH, CheapSynchronousRK2, Courant 0.25, Spheral's own time-step vote. A uniform expansion
v = edot x (edot 10 /s, from the centre) is imposed at the start so that the strength (and the damage, when on) does
real work: the elastic stress grows from zero instead of staying zero.

Timing. Measured separately: imports (Spheral and its compiled modules), body generation and distribution over
the processes, set-up (physics packages, controller, the initial smoothing-length iteration), --warmup steps
(discarded, through the controller's advance), then each of --steps timed steps on its own, between MPI barriers,
as rank 0 sees it: the integrator step plus the controller's periodic-work check, i.e. the body of advance's loop.
No output, statistics, restarts or redistribution happen during the timed steps; the controller still runs Python
garbage collection every 100 steps, which production pays too. The numbers are native arm64: Spheral runs in the
arm64 image under Docker Desktop's Linux VM, with no Rosetta (the cluster will be x86-64).

The load average checked is the Linux VM's (/proc/loadavg in the container), read on rank 0 before Spheral is
imported; above --max-load (2.0) the run refuses to start (exit 2) unless --force. The host's load is invisible from
the container; a caller that measured it passes it with --host-load and it is recorded (run_benchmark.sh does).

Run through the launcher from the repository root, e.g.

    spheral_frag/container/spheral -n 18 spheral_frag/m0/bench.py --geometry 3d --dx-mm 2.2 --damage off

Writes OUT/<run name>/summary.json and steps.csv (per-step wall time and dt); the run name encodes the whole
configuration, and a run whose summary.json exists is skipped (exit 0). Exit 0 ok, 1 a Spheral failure, 2 bad
arguments or a machine that is not idle. --dt-votes prints each physics package's time-step vote and its reason
after the timed steps (how the Gruneisen sound speed was found).
"""
import argparse
import json
import math
import os
import platform
import resource
import sys
import time

SPHERAL_COMMIT = "116c71f"
PIN = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "container", "pin.json")

# Timing-only material, SI units (see the module docstring for the sources).
RHO0 = 2813.0                   # kg/m^3, AA7075
ETAMIN, ETAMAX = 0.2, 4.0
YOUNG, POISSON = 71.7e9, 0.33
MURNAGHAN = dict(n=1.0, K=YOUNG/(3.0*(1.0 - 2.0*POISSON)), atomicWeight=26.98)   # K = 70.3 GPa
SG = dict(G0=27.6e9, A=6.52e-11, B=6.16e-4, Y0=0.29e9, Ymax=0.68e9, Yp=1.0e-3, beta=125.0, gamma0=0.0, nhard=0.10)
TMELT0 = 1220.0                 # K, SCG 1980 for 6061-T6
CV_DULONG_PETIT = 3.0*8.314462618/26.98e-3   # J/(kg K), 924.5
EMELT_PER_VOLUME = RHO0*CV_DULONG_PETIT*(TMELT0 - 300.0)   # J/m^3, the constant melt-energy fit
EDOT = 10.0                     # 1/s, imposed uniform expansion
SEED = 48927592                 # ProbabilisticDamageModel's documented default seed, fixed
NPERH = 2.01
CFL = 0.25


def parse_args(argv):
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--geometry", choices=("3d", "rz"), required=True)
    p.add_argument("--dx-mm", type=float, required=True, help="lattice spacing in mm")
    p.add_argument("--damage", choices=("on", "off"), required=True)
    p.add_argument("--steps", type=int, default=100, help="timed steps (default 100)")
    p.add_argument("--warmup", type=int, default=10, help="discarded steps before timing (default 10)")
    p.add_argument("--out", default="spheral_output/m0/bench", help="parent of the run directory")
    p.add_argument("--radius-mm", type=float, default=50.0, help="sphere radius in mm (default 50: the 100 mm body)")
    p.add_argument("--max-load", type=float, default=2.0, help="refuse to start above this 1-min load average")
    p.add_argument("--force", action="store_true", help="start even above --max-load")
    p.add_argument("--host-load", default=None,
                   help="the host's load averages measured by the caller before launching (the container cannot see "
                        "them), recorded in summary.json, e.g. \"1.2 1.5 1.7\"")
    p.add_argument("--dt-votes", action="store_true",
                   help="after the timed steps, print every physics package's time-step vote and its reason")
    a = p.parse_args(argv)
    if a.dx_mm <= 0 or a.radius_mm <= 0 or a.steps < 1 or a.warmup < 0 or a.dx_mm > a.radius_mm:
        p.error("need 0 < dx <= radius, steps >= 1, warmup >= 0")
    return a


def mpi_size_from_env():
    """The process count before anything MPI is imported (mpirun sets it; a bare run is one process)."""
    for k in ("OMPI_COMM_WORLD_SIZE", "PMI_SIZE"):
        if k in os.environ:
            return int(os.environ[k])
    return 1


def run_name(a, nproc):
    return (f"bench_{a.geometry}_dx{a.dx_mm:.3f}mm_R{a.radius_mm:g}mm_n{nproc}_dmg{a.damage}"
            f"_steps{a.steps}_wu{a.warmup}_{SPHERAL_COMMIT}")


def stats(xs):
    xs = sorted(xs)
    n = len(xs)
    med = xs[n//2] if n % 2 else 0.5*(xs[n//2 - 1] + xs[n//2])
    return dict(mean=sum(xs)/n, median=med, min=xs[0], max=xs[-1])


def main(argv):
    try:
        a = parse_args(argv)
    except SystemExit as e:
        return 2 if e.code else 0

    nproc = mpi_size_from_env()
    rank = int(os.environ.get("OMPI_COMM_WORLD_RANK", os.environ.get("PMI_RANK", "0")))
    name = run_name(a, nproc)
    rundir = os.path.join(a.out, name)
    summary_path = os.path.join(rundir, "summary.json")
    if os.path.exists(summary_path):
        if rank == 0:
            print(f"bench: {name} exists, skipped")
        return 0

    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    load = os.getloadavg() if comm.rank == 0 else None
    load = comm.bcast(load, root=0)
    if load[0] > a.max_load and not a.force:
        if comm.rank == 0:
            print(f"bench: load average {load[0]:.2f} > {a.max_load} (the Linux VM's); not starting (--force overrides)")
        return 2

    # ---------------------------------------------------------------- imports
    t0 = time.perf_counter()
    import mpi
    if a.geometry == "3d":
        import Spheral3d as S
        from GenerateNodeDistribution3d import GenerateNodeDistribution3d
        from PeanoHilbertDistributeNodes import distributeNodes3d as distributeNodes
    else:
        import SpheralRZ as S
        from GenerateNodeDistribution2d import GenerateNodeDistribution2d, RZGenerator
        from PeanoHilbertDistributeNodes import distributeNodes2d as distributeNodes
    from SpheralController import SpheralController
    mpi.barrier()
    t_import = time.perf_counter() - t0

    R = 1e-3*a.radius_mm
    dx = 1e-3*a.dx_mm

    # ---------------------------------------------------------------- material
    t0 = time.perf_counter()
    units = S.MKS()
    eos = S.MurnaghanEquationOfState(RHO0, ETAMIN, ETAMAX, MURNAGHAN["n"], MURNAGHAN["K"], MURNAGHAN["atomicWeight"],
                                     units)
    zero = S.NinthOrderPolynomialFit(*([0.0]*10))
    melt = S.NinthOrderPolynomialFit(*([EMELT_PER_VOLUME] + [0.0]*9))
    strength = S.SteinbergGuinanStrength(eos, SG["G0"], SG["A"], SG["B"], SG["Y0"], SG["Ymax"], SG["Yp"],
                                         SG["beta"], SG["gamma0"], SG["nhard"], zero, melt)
    WT = S.TableKernel(S.WendlandC4Kernel(), 1000)
    big = 10.0*R
    box = (S.Vector(-big, -big, -big), S.Vector(big, big, big)) if a.geometry == "3d" else \
          (S.Vector(-big, -big), S.Vector(big, big))
    nodes = S.makeSolidNodeList("AA7075", eos, strength, nPerh=NPERH, kernelExtent=WT.kernelExtent,
                                hmin=1e-3*dx, hmax=100.0*dx, rhoMin=ETAMIN*RHO0, rhoMax=ETAMAX*RHO0,
                                xmin=box[0], xmax=box[1])
    t_material = time.perf_counter() - t0

    # ---------------------------------------------------------------- body
    t0 = time.perf_counter()
    k = int(math.floor(R/dx*(1.0 + 1e-12)))
    rclip = R*(1.0 + 1e-9)
    if a.geometry == "3d":
        n1 = 2*k + 1
        lo = -(k + 0.5)*dx
        gen = GenerateNodeDistribution3d(n1, n1, n1, RHO0, "lattice", xmin=(lo, lo, lo), xmax=(-lo, -lo, -lo),
                                         rmax=rclip, origin=(0.0, 0.0, 0.0), nNodePerh=NPERH, SPH=True)
        expected = (4.0/3.0)*math.pi*R**3/dx**3
    else:
        nz = 2*k + 1
        nr = int(math.ceil(R/dx - 0.5)) + 1
        zlo = -(k + 0.5)*dx
        gen = RZGenerator(GenerateNodeDistribution2d(nz, nr, RHO0, "lattice", xmin=(zlo, 0.0),
                                                     xmax=(-zlo, nr*dx), rmax=rclip, nNodePerh=NPERH, SPH=True))
        expected = 0.5*math.pi*R**2/dx**2
    distributeNodes((nodes, gen))
    nlocal = nodes.numInternalNodes
    ntotal = mpi.allreduce(nlocal, mpi.SUM)
    nmin_rank = mpi.allreduce(nlocal, mpi.MIN)
    nmax_rank = mpi.allreduce(nlocal, mpi.MAX)
    nodes.specificThermalEnergy(S.ScalarField("eps0", nodes, 0.0))
    pos = nodes.positions()
    vel = nodes.velocity()
    for i in range(nlocal):
        vel[i] = pos[i]*EDOT          # uniform expansion about the centre; in RZ (z, r) -> (edot z, edot r)
    mpi.barrier()
    t_body = time.perf_counter() - t0

    # ---------------------------------------------------------------- physics and controller
    t0 = time.perf_counter()
    db = S.DataBase()
    db.appendNodeList(nodes)
    hydro = S.SPH(dataBase=db, W=WT, cfl=CFL, compatibleEnergyEvolution=True, densityUpdate=S.IntegrateDensity,
                  XSPH=False, ASPH=False)
    integrator = S.CheapSynchronousRK2Integrator(db)
    integrator.appendPhysicsPackage(hydro)
    damage_model = None
    if a.damage == "on":
        damage_model = S.ProbabilisticDamageModel(materialName="aluminum", units=units, nodeList=nodes,
                                                  kernel=WT, seed=SEED)
        integrator.appendPhysicsPackage(damage_model)
    cs0 = math.sqrt((MURNAGHAN["K"] + (4.0/3.0)*SG["G0"])/RHO0)
    integrator.lastDt = 0.1*CFL*dx/cs0       # a small first guess; Spheral's vote takes over within a few steps
    integrator.dtMin = 1e-15
    integrator.dtMax = 1.0
    integrator.verbose = False
    controller = SpheralController(integrator, WT, statsStep=10**9, printStep=10**9, redistributeStep=None,
                                   restartStep=None, vizStep=None, vizTime=None, SPH=True)
    mpi.barrier()
    t_setup = time.perf_counter() - t0

    # ---------------------------------------------------------------- steps
    try:
        t0 = time.perf_counter()
        if a.warmup:
            controller.advance(1e40, a.warmup)
        mpi.barrier()
        t_warmup = time.perf_counter() - t0

        # The body of SpheralController.advance's loop, one step at a time. controller.step() is advance(1e40, 1),
        # which forces all the periodic work (a global conservation sum, a print, garbage collection) and prints a
        # timer report after every call; production's advance does that work only at its frequencies.
        walls, dts = [], []
        for _ in range(a.steps):
            mpi.barrier()
            ts = time.perf_counter()
            integrator.step(1e40)
            controller.totalSteps += 1
            controller.doPeriodicWork()
            mpi.barrier()
            walls.append(time.perf_counter() - ts)
            dts.append(integrator.lastDt)
    except Exception as e:  # noqa: BLE001 -- any Spheral failure is a model failure (exit 1)
        print(f"bench: Spheral failed on rank {mpi.rank}: {type(e).__name__}: {e}")
        return 1

    if a.dt_votes:
        st = S.State(db, integrator.physicsPackages())
        dv = S.StateDerivatives(db, integrator.physicsPackages())
        integrator.initializeDerivatives(controller.time(), integrator.lastDt, st, dv)
        integrator.evaluateDerivatives(controller.time(), integrator.lastDt, db, st, dv)
        for pkg in integrator.physicsPackages():
            v = pkg.dt(db, st, dv, controller.time())
            print("DTVOTE", type(pkg).__name__, v[0], v[1])
        print("DTVOTE integrator", integrator.lastDt, integrator.dtMin, integrator.dtMax, integrator.dtGrowth)
    rss_kb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss       # Linux: kilobytes
    rss_total = mpi.allreduce(rss_kb, mpi.SUM)*1024
    rss_max = mpi.allreduce(rss_kb, mpi.MAX)*1024
    damage_max = None
    if damage_model is not None:
        D = nodes.damage()
        dloc = max([D[i].eigenValues().maxElement() for i in range(nlocal)] or [0.0])
        damage_max = mpi.allreduce(dloc, mpi.MAX)

    if mpi.rank != 0:
        return 0

    w = stats(walls)
    pin = json.load(open(PIN))
    summary = dict(
        run=name, spheral_commit=SPHERAL_COMMIT, geometry=a.geometry,
        dx_mm=a.dx_mm, radius_mm=a.radius_mm, particles=ntotal,
        particles_expected=expected, particles_ratio=ntotal/expected,
        particles_per_process_min=nmin_rank, particles_per_process_max=nmax_rank,
        processes=mpi.procs, omp_num_threads=os.environ.get("OMP_NUM_THREADS"),
        damage=a.damage, steps=a.steps, warmup=a.warmup,
        import_s=t_import, material_s=t_material, body_generation_s=t_body, setup_s=t_setup,
        startup_s=t_import + t_material + t_setup, warmup_s=t_warmup,
        step_wall_s=w, step_wall_total_s=sum(walls),
        processor_s_per_particle_step=mpi.procs*w["median"]/ntotal,
        dt_s=stats(dts), time_reached_s=controller.time(),
        peak_rss_bytes_max_process=rss_max, peak_rss_bytes_total=rss_total,
        memory_bytes_per_particle=rss_total/ntotal,
        load_average_at_start=list(load), load_average_source="Linux VM (/proc/loadavg in the container), rank 0",
        host_load_average_at_start=[float(x) for x in a.host_load.split()] if a.host_load else None,
        forced=bool(a.force and load[0] > a.max_load),
        damage_max_eigenvalue_end=damage_max,
        machine=platform.machine(),
        platform_note="native arm64 Spheral in Docker Desktop's Linux VM; no Rosetta. The cluster will be x86-64.",
        material="timing only",
        material_detail=dict(
            eos="MurnaghanEquationOfState, exponent 1, K = E/(3(1 - 2 nu)) of AA7075 (not Gruneisen: its dP/drho "
                "for rho <= rho0 is rho0 times too large at 116c71f)",
            murnaghan=MURNAGHAN, rho0=RHO0, etamin=ETAMIN, etamax=ETAMAX,
            strength="SteinbergGuinanStrength, 6061-T6 Al from Steinberg, Cochran & Guinan 1980 (Spheral's "
                     "library has no Steinberg-Guinan set); cold fit 0, melt fit constant rho0 cv (Tm0 - 300 K)",
            steinberg_guinan=SG, Tmelt0_K=TMELT0, melt_energy_per_volume_J_m3=EMELT_PER_VOLUME,
            damage=("ProbabilisticDamageModel, library 'aluminum' kWeibull/mWeibull, documented defaults, seed "
                    f"{SEED}") if a.damage == "on" else "off",
            temperature="the strength model's temperature comes from Murnaghan's Dulong-Petit relation; it "
                        "shifts the shear modulus slightly and does not change the cost"),
        hydro=dict(scheme="SolidSPH (SPH factory)", kernel="WendlandC4 (TableKernel, 1000)", nPerh=NPERH,
                   densityUpdate="IntegrateDensity", compatibleEnergy=True, XSPH=False, ASPH=False,
                   integrator="CheapSynchronousRK2", cfl=CFL, imposed_expansion_rate_per_s=EDOT),
        pin=pin,
    )
    os.makedirs(rundir, exist_ok=True)
    with open(os.path.join(rundir, "steps.csv"), "w") as f:
        f.write("step,wall_s,dt_s\n")
        for i, (wt, d) in enumerate(zip(walls, dts)):
            f.write(f"{i + 1},{wt:.6e},{d:.6e}\n")
    with open(summary_path, "w") as f:
        json.dump(summary, f, indent=2)
        f.write("\n")
    print(f"bench: {name}: {ntotal} particles ({ntotal/expected - 1:+.2%} of the formula) on {mpi.procs} "
          f"processes, median step {w['median']:.4g} s, {summary['processor_s_per_particle_step']:.3g} "
          f"processor-s per particle-step, dt {summary['dt_s']['median']:.3g} s, "
          f"{summary['memory_bytes_per_particle']/1024:.1f} KiB per particle -> {summary_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))

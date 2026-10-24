"""Stress rotation under rigid-body rotation, checked numerically (Spheral M0 plan, Task 4).

A 2D square of AA7075-like solid carries a uniform deviatoric stress S0 (pure shear) and is spun about its centre
with v = omega z x (x - c) for a quarter turn. Under a correct Jaumann rate, dS/dt = W S - S W with W the spin, the
stress turns with the body: S(t) = R(theta) S0 R(theta)^T for the angle theta the body actually turned through, and
the von Mises stress is unchanged. This exercises SolidSPH's

    spinCorrection = (spin*Si - Si*spin).Symmetric()      (src/SPH/SolidSPH.cc, develop 116c71f, PR #520)

How rigidity is enforced. Nothing in a free body supplies the centripetal force of a rigid rotation, and the free
surfaces of a pre-stressed square are not traction-free, so a free run would fly apart radially and relax its edge
stress by elastic waves that cross the 20-particle square within the quarter turn: the stress would change for
reasons other than rotation. Instead a small Python physics package (RigidRotation) resets the velocity to the exact
rigid field v = omega z x (x - c), from the current positions, in postStateUpdate, which the integrator calls after
every stage's state update. The velocity field the stress update differentiates is then an exact spin (no
symmetric part, so no deformation rate enters S), and only the spin term moves S. The hydro's own accelerations are
computed and discarded. The positions are still integrated from v, so they drift outward by the integrator's
O((omega dt)^4) per step; that drift is reported as `radial_strain` (1.8e-6 at the default settings).

The time step is fixed at --courant-fraction (0.5) of Spheral's own first time-step vote (the sound-speed limit,
cs = sqrt(K/rho) = 5.0 km/s), and omega is chosen so that --steps of it make a quarter turn (omega dt = 0.021 rad at
75 steps). The comparison angle is the least-squares rotation of the particle positions about their centroid.

Measured 2026-10-05 at the defaults (20 x 20, 256 interior particles, 75 steps, SynchronousRK2): largest stress error
over the history 0.069 % of |S0|, von Mises drift 0.003 %; the error falls as 1/steps^2 (0.155 %, 0.039 %, 0.017 %
at 50, 100, 150 steps), so it is the integrator's phase error, not the spin term. CheapSynchronousRK2 gives 0.19 %
and 0.094 %: its first stage reuses the previous step's derivatives, which on the first step are the start-up ones,
evaluated before the spin was set, so its first step is forward Euler.

Why the check is made over the whole history, not only at the quarter turn. A deviatoric tensor in the plane
returns to itself under a half turn, so R(+pi/2) S0 R^T = R(-pi/2) S0 R^T: a wrong-sign spin term (the stress
turning by -theta while the body turns by +theta, a relative -2 theta) agrees with the right answer at exactly a
quarter turn. The two differ by up to 2 |S0| at an eighth of a turn. The stress is therefore compared with
R S0 R^T at every step, and the same comparison against the wrong-sign answer is reported to show the test's power.

Run through the launcher from the repository root:

    spheral_frag/container/spheral spheral_frag/m0/rigid_rotation.py [--nx 20] [--steps 75] [--shear-mpa 10]

Writes spheral_output/m0/rotation/<run name>/summary.json (and history.csv). Exit 0 when both checks pass within
--tol (default 1 %), 1 when a check fails, 2 on bad arguments.
"""
import argparse
import json
import math
import os
import sys
import time as walltime

import numpy as np

SPHERAL_COMMIT = "116c71f"

# AA7075-like solid, SI units.
RHO0 = 2813.0           # kg/m^3
YOUNG = 71.7e9          # Pa
POISSON = 0.33
K_BULK = YOUNG/(3.0*(1.0 - 2.0*POISSON))   # 70.3 GPa
G_SHEAR = YOUNG/(2.0*(1.0 + POISSON))      # 26.95 GPa
YIELD = 1.0e12          # Pa, far above any stress here: the material stays elastic


def parse_args(argv):
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--nx", type=int, default=20, help="particles per side (default 20)")
    p.add_argument("--side", type=float, default=0.1, help="side of the square in m (default 0.1)")
    p.add_argument("--steps", type=int, default=75, help="steps per quarter turn, 50-100 per the plan (default 75)")
    p.add_argument("--shear-mpa", type=float, default=10.0, help="initial pure shear S0_xy in MPa (default 10)")
    p.add_argument("--courant-fraction", type=float, default=0.5,
                   help="fixed step as a fraction of Spheral's own first time-step vote (default 0.5)")
    p.add_argument("--integrator", choices=("rk2", "cheaprk2"), default="rk2",
                   help="SynchronousRK2 (default) or CheapSynchronousRK2, whose first stage reuses the previous "
                        "step's derivatives; on the first step those are the start-up ones, evaluated before the "
                        "spin is set, so its first step is forward Euler")
    p.add_argument("--tol", type=float, default=0.01, help="relative tolerance of both checks (default 0.01)")
    p.add_argument("--outdir", default="spheral_output/m0/rotation", help="parent of the run directory")
    a = p.parse_args(argv)
    if a.nx < 8 or a.steps < 10 or a.side <= 0 or a.shear_mpa <= 0 or not 0 < a.courant_fraction <= 1 or a.tol <= 0:
        p.error("need nx >= 8, steps >= 10, side > 0, shear > 0, 0 < courant-fraction <= 1, tol > 0")
    return a


def rot(theta):
    c, s = math.cos(theta), math.sin(theta)
    return np.array([[c, -s], [s, c]])


def deviator3(s2):
    """2x2 in-plane deviatoric stress -> 3x3 with S_zz = -(S_xx + S_yy) (Spheral's 2D S is traceless in 3D)."""
    s = np.zeros((3, 3))
    s[:2, :2] = s2
    s[2, 2] = -(s2[0, 0] + s2[1, 1])
    return s


def von_mises(s2):
    s = deviator3(s2)
    return math.sqrt(1.5*np.sum(s*s))


def frob3(s2):
    return math.sqrt(np.sum(deviator3(s2)**2))


def main(argv):
    try:
        a = parse_args(argv)
    except SystemExit as e:
        return 2 if e.code else 0

    import mpi
    from SolidSpheral2d import (MKS, LinearPolynomialEquationOfState, ConstantStrength, makeSolidNodeList,
                                TableKernel, WendlandC4Kernel, DataBase, SPH, IntegrateDensity,
                                SynchronousRK2Integrator, CheapSynchronousRK2Integrator, Physics, HydroFieldNames,
                                ScalarField, Vector, SymTensor, State, StateDerivatives)
    from SpheralController import SpheralController
    from GenerateNodeDistribution2d import GenerateNodeDistribution2d
    from DistributeNodes import distributeNodes2d

    if mpi.procs != 1:
        print("rigid_rotation: run on one process (the checks gather particles locally)")
        return 2

    name = f"rotation_2d_n{a.nx}_L{a.side*1e3:.1f}mm_S{a.shear_mpa:g}MPa_steps{a.steps}_cf{a.courant_fraction:g}_{a.integrator}_{SPHERAL_COMMIT}"
    rundir = os.path.join(a.outdir, name)
    os.makedirs(rundir, exist_ok=True)
    t_start = walltime.time()

    units = MKS()
    eos = LinearPolynomialEquationOfState(RHO0, 0.5, 1.5,       # rho0, etamin, etamax
                                          0.0, K_BULK, 0.0, 0.0,  # P = a1 (rho/rho0 - 1): bulk modulus K
                                          0.0, 0.0, 0.0,          # no thermal (Gruneisen) terms
                                          27.0, units)
    strength = ConstantStrength(G_SHEAR, YIELD)
    WT = TableKernel(WendlandC4Kernel(), 1000)
    nPerh = 2.01
    dx = a.side/a.nx
    half = 0.5*a.side
    nodes = makeSolidNodeList("AA7075", eos, strength, nPerh=nPerh, kernelExtent=WT.kernelExtent,
                              hmin=1e-3*dx, hmax=100.0*dx, xmin=Vector(-10*a.side, -10*a.side),
                              xmax=Vector(10*a.side, 10*a.side))
    gen = GenerateNodeDistribution2d(a.nx, a.nx, RHO0, "lattice", xmin=(-half, -half), xmax=(half, half),
                                     nNodePerh=nPerh)
    distributeNodes2d((nodes, gen))
    n = nodes.numInternalNodes
    nodes.specificThermalEnergy(ScalarField("eps0", nodes, 0.0))

    centre = np.zeros(2)   # the lattice is symmetric about the origin
    x0 = np.array([[nodes.positions()[i].x, nodes.positions()[i].y] for i in range(n)])
    omega = 0.0            # set below once the time step is known
    tau = a.shear_mpa*1e6
    S0 = np.array([[0.0, tau], [tau, 0.0]])
    for i in range(n):
        nodes.deviatoricStress()[i] = SymTensor(S0[0, 0], S0[0, 1], S0[1, 0], S0[1, 1])

    # Interior: at least one kernel support (kernelExtent * h) inside every edge at the start.
    support = WT.kernelExtent*nPerh*dx
    interior = np.all(np.abs(x0) <= half - support + 1e-12*a.side, axis=1)
    if interior.sum() < 4:
        print(f"rigid_rotation: only {interior.sum()} interior particles at nx = {a.nx}; increase --nx")
        return 2

    class RigidRotation(Physics):
        """Holds v = omega z x (x - c) at every stage, from the current positions; votes no time step."""

        def __init__(self):
            Physics.__init__(self)
            self.omega = 0.0

        def evaluateDerivatives(self, t, dt, db, state, derivs):
            return

        def dt(self, db, state, derivs, t):
            return (1e100, "No vote")   # TimeStepType = pair<double, string>

        def registerState(self, db, state):
            return

        def registerDerivatives(self, db, derivs):
            return

        def label(self):
            return "RigidRotation"

        def impose(self, state):
            for pos, vel in zip(state.vectorFields(HydroFieldNames.position),
                                state.vectorFields(HydroFieldNames.velocity)):
                for i in range(pos.numInternalElements):
                    p = pos[i]
                    vel[i] = Vector(-self.omega*(p.y - centre[1]), self.omega*(p.x - centre[0]))

        def postStateUpdate(self, t, dt, db, state, derivs):
            self.impose(state)
            return True

    db = DataBase()
    db.appendNodeList(nodes)
    hydro = SPH(dataBase=db, W=WT, cfl=0.25, compatibleEnergyEvolution=True, XSPH=False,
                correctVelocityGradient=True, densityUpdate=IntegrateDensity)
    rigid = RigidRotation()
    integrator = (SynchronousRK2Integrator if a.integrator == "rk2" else CheapSynchronousRK2Integrator)(db)
    integrator.appendPhysicsPackage(hydro)
    integrator.appendPhysicsPackage(rigid)   # last, so its postStateUpdate has the final word on v
    integrator.verbose = False
    controller = SpheralController(integrator, WT, statsStep=10**9, printStep=10**9, restartStep=None,
                                   vizStep=None, SPH=True)

    # Spheral's own first time-step vote for the body at rest with S0 in it; the fixed step is a fraction of it.
    # The rigid rotation itself only adds the spin-rate vote, which is far weaker at omega dt ~ 0.02 rad.
    packages = integrator.physicsPackages()
    state = State(db, packages)
    derivs = StateDerivatives(db, packages)
    vote = hydro.dt(db, state, derivs, 0.0)
    dt_courant, dt_reason = vote[0], vote[1]
    dt = a.courant_fraction*dt_courant
    t_quarter = a.steps*dt
    omega = 0.5*math.pi/t_quarter
    rigid.omega = omega
    for i in range(n):
        nodes.velocity()[i] = Vector(-omega*(x0[i, 1] - centre[1]), omega*(x0[i, 0] - centre[0]))
    integrator.dtMin = dt
    integrator.dtMax = dt
    integrator.lastDt = dt
    integrator.dtGrowth = 1.0

    def snapshot():
        x = np.array([[nodes.positions()[i].x, nodes.positions()[i].y] for i in range(n)])
        S = np.array([[[nodes.deviatoricStress()[i].xx, nodes.deviatoricStress()[i].xy],
                       [nodes.deviatoricStress()[i].yx, nodes.deviatoricStress()[i].yy]] for i in range(n)])
        return x, S

    s0n = frob3(S0)
    vm0 = von_mises(S0)
    rows = []

    def measure(step):
        x, S = snapshot()
        xc0 = x0 - x0.mean(axis=0)
        xc = x - x.mean(axis=0)
        # Least-squares rotation angle of the whole body about its centroid, from the positions.
        theta = math.atan2(np.sum(xc0[:, 0]*xc[:, 1] - xc0[:, 1]*xc[:, 0]), np.sum(xc0*xc))
        r0 = np.linalg.norm(xc0[interior], axis=1)
        r = np.linalg.norm(xc[interior], axis=1)
        radial = float(np.max(np.abs(r/r0 - 1.0)))
        right = rot(theta) @ S0 @ rot(theta).T
        wrong = rot(-theta) @ S0 @ rot(-theta).T
        Si = S[interior]
        err = max(frob3(s - right) for s in Si)/s0n
        err_wrong = max(frob3(s - wrong) for s in Si)/s0n
        err_none = max(frob3(s - S0) for s in Si)/s0n
        vm = max(abs(von_mises(s)/vm0 - 1.0) for s in Si)
        rows.append(dict(step=step, time=controller.time(), theta=theta, theta_omega_t=omega*controller.time(),
                         err=err, err_wrong_sign=err_wrong, err_no_spin=err_none, von_mises_drift=vm,
                         radial_strain=radial, err_all_particles=max(frob3(s - right) for s in S)/s0n))

    measure(0)
    for k in range(1, a.steps + 1):
        controller.step()
        measure(k)

    with open(os.path.join(rundir, "history.csv"), "w") as f:
        keys = list(rows[0])
        f.write(",".join(keys) + "\n")
        for r in rows:
            f.write(",".join(f"{r[k]:.10g}" for k in keys) + "\n")

    last = rows[-1]
    worst = lambda k: max(r[k] for r in rows)
    summary = dict(
        run=name, spheral_commit=SPHERAL_COMMIT, hydro=f"SolidSPH (SPH factory), {type(integrator).__name__}, "
        "correctVelocityGradient, IntegrateDensity, WendlandC4, nPerh 2.01, XSPH off",
        rigid_rotation_enforcement="velocity reset to omega z x (x - c) from current positions in postStateUpdate "
        "after every integrator stage",
        material=dict(rho0=RHO0, E=YOUNG, nu=POISSON, K=K_BULK, G=G_SHEAR, yield_strength=YIELD,
                      eos="LinearPolynomial P = K (rho/rho0 - 1)", strength="ConstantStrength"),
        nx=a.nx, particles=n, interior_particles=int(interior.sum()), side_m=a.side, dx_m=dx,
        S0_xy_Pa=tau, S0_frobenius_Pa=s0n, von_mises0_Pa=vm0,
        dt_courant_vote_s=dt_courant, dt_courant_reason=dt_reason, dt_s=dt, omega_rad_s=omega,
        omega_dt_rad=omega*dt, steps=len(rows) - 1, final_time_s=last["time"],
        theta_final_rad=last["theta"], theta_final_deg=math.degrees(last["theta"]),
        theta_omega_t_final_rad=last["theta_omega_t"],
        tol=a.tol,
        max_err_over_history=worst("err"), err_final=last["err"],
        max_von_mises_drift_over_history=worst("von_mises_drift"), von_mises_drift_final=last["von_mises_drift"],
        max_radial_strain=worst("radial_strain"),
        max_err_all_particles_over_history=worst("err_all_particles"),
        wrong_sign_max_err_over_history=worst("err_wrong_sign"), wrong_sign_err_final=last["err_wrong_sign"],
        no_spin_max_err_over_history=worst("err_no_spin"),
        wall_time_s=walltime.time() - t_start,
    )
    summary["pass_rotation"] = summary["max_err_over_history"] <= a.tol
    summary["pass_von_mises"] = summary["max_von_mises_drift_over_history"] <= a.tol
    summary["wrong_sign_detected"] = summary["wrong_sign_max_err_over_history"] > a.tol
    summary["quarter_turn"] = abs(last["theta"] - 0.5*math.pi) < 0.02
    ok = summary["pass_rotation"] and summary["pass_von_mises"] and summary["quarter_turn"]
    summary["pass"] = ok
    with open(os.path.join(rundir, "summary.json"), "w") as f:
        json.dump(summary, f, indent=2)
        f.write("\n")
    print(json.dumps(summary, indent=2))
    print(f"rigid_rotation: {'PASS' if ok else 'FAIL'}  -> {os.path.join(rundir, 'summary.json')}")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))

"""Stage B: the box model of the plume, one stiff ODE system of extensive counts (molecules and atoms in the whole plume).

The plume is a horizontal line segment of length L (the trail) broadened into a Gaussian of mean-square radius r^2
(spec 4.2; plan decisions 1, 2, 12 and 19):

  * State. Total counts in the plume, not concentrations, so aluminium conservation is a plain sum (decision 1).
    Concentrations are counts / V with V = volume(r, L) = pi r^2 L / erf(L / 2r) (decision 2).
  * Radius. d(r^2)/dt = 4 D_total(r) with D for AlO in air (molecular plus eddy, transport.D_total); zero when
    Settings.freeze_radius. dV/dt = V dlnV_dr2(r, L) d(r^2)/dt.
  * Entrainment. Each ambient species gains c_ambient dV/dt; the aluminium-bearing species have no ambient source,
    so they are diluted only through V. Written for concentrations this is the spec's
    dc/dt = chemistry + (c_ambient - c) dlnV/dt.
  * Chemistry. Every box reaction runs on concentrations at the ambient temperature amb.T_K (held fixed for the run,
    so each rate constant is evaluated once, with a RateWarnings collector exposed as Box.warnings). Three-body steps
    use [M] = amb.n_total_m3 (decision 12): the plume's own total (the sum of its ambient-species concentrations)
    equals this to trace precision, since entrainment keeps every ambient species at its ambient value apart from
    the aluminium-driven consumption. Photolysis (R17, day_only) runs only when time_of_day == "day". The
    ambient species are dynamic: the aluminium reactions consume and produce them.
  * Particles. The layout carries the particle slices (bin numbers N, tallies by source molecule and route, O and H
    atoms per bin, smoke uptake by source, settled atoms by source and route); in this version (Task 7, gas phase
    only) their right-hand side is zero. Task 8 fills them.
  * Scaling (decision 19). The solver sees scaled variables: aluminium-bearing counts (and every particle slice)
    divided by the aluminium inventory (or by 1 when it is 0), each ambient count by c_ambient V0 (V0 when
    c_ambient is 0), r^2 by r0^2. BDF runs with rtol = Settings.rtol and atol = 1e-12 on these.

Two things are called "reactions": `Settings.reactions` is a tuple of reaction ids that restricts the network (None
keeps all of it), while the `reactions=` keyword of Box is a list of Reaction objects to draw from (default: the box
reactions of data/wake_column/reactions.csv).
"""
import math

import numpy as np
from scipy import sparse
from scipy.integrate import solve_ivp

from reentry_model.wake_column.particles import Bins, SmokeBackground, diameter, epstein_speed, volume_per_al
from reentry_model.wake_column.reactions import RateWarnings, load_reactions
from reentry_model.wake_column.species import AL_SPECIES, AMBIENT_SPECIES, AMU, load_species
from reentry_model.wake_column.transport import D_total, EddyProfile, dlnV_dr2, settling_rate, volume

GAS_SPECIES = AL_SPECIES + AMBIENT_SPECIES
SOURCES = AL_SPECIES                    # the molecule an aluminium atom arrived as
ROUTES = ("formation", "growth")        # how it entered a plume particle (smoke uptake has its own slice)
ATOL = 1e-12                            # absolute tolerance on the scaled variables (decision 19)
G0 = 9.80665                            # m/s^2 at the surface
R_EARTH_KM = 6371.0


class Layout:
    """Index slices of the state vector, in order: gas (the 7 aluminium-bearing species then the 9 ambient ones, in
    GAS_SPECIES order), r2 (one entry, m^2), N (K bins, particles), tally (K x 7 sources x 2 routes, Al atoms), O and H
    (K bins, atoms), smoke (7 sources, Al atoms taken up by background smoke), settled (7 sources x 2 routes, Al atoms
    that settled out)."""

    def __init__(self, K):
        self.K = int(K)
        self.gas_names = GAS_SPECIES
        self.n_al = len(AL_SPECIES)
        sizes = (("gas", len(GAS_SPECIES)), ("r2", 1), ("N", self.K), ("tally", self.K * len(SOURCES) * len(ROUTES)),
                 ("O", self.K), ("H", self.K), ("smoke", len(SOURCES)), ("settled", len(SOURCES) * len(ROUTES)))
        start = 0
        for name, n in sizes:
            setattr(self, name, slice(start, start + n))
            start += n
        self.size = start
        self.al = slice(self.gas.start, self.gas.start + self.n_al)             # aluminium-bearing gas species
        self.ambient = slice(self.gas.start + self.n_al, self.gas.stop)         # ambient gas species
        self._index = {name: self.gas.start + i for i, name in enumerate(GAS_SPECIES)}

    def index(self, name):
        """Position of gas species `name` in the state vector."""
        return self._index[name]

    def unpack(self, y):
        """Views of y by slice: gas (16,), r2 (float), N (K,), tally (K, 7, 2), O (K,), H (K,), smoke (7,),
        settled (7, 2)."""
        y = np.asarray(y)
        return {"gas": y[self.gas], "r2": float(y[self.r2][0]), "N": y[self.N],
                "tally": y[self.tally].reshape(self.K, len(SOURCES), len(ROUTES)), "O": y[self.O], "H": y[self.H],
                "smoke": y[self.smoke], "settled": y[self.settled].reshape(len(SOURCES), len(ROUTES))}


class Box:
    """The plume of one trail: initial state, right-hand side, Jacobian and solver (module docstring)."""

    def __init__(self, trail, hot, amb, time_of_day, settings, *, species=None, reactions=None, bins=None, smoke=None,
                 eddy=None, initial_bins=None):
        if time_of_day not in ("day", "night"):
            raise ValueError(f"time_of_day must be 'day' or 'night', got {time_of_day!r}")
        self.trail, self.hot, self.amb, self.time_of_day, self.settings = trail, hot, amb, time_of_day, settings
        self.species = species if species is not None else load_species()
        self.bins = bins if bins is not None else Bins(settings.bin_ratio, settings.d_max_nm)
        self.smoke = smoke if smoke is not None else SmokeBackground.load()
        self.eddy = eddy if eddy is not None else EddyProfile.load()
        self.layout = Layout(self.bins.K)
        self.initial_bins = dict(initial_bins or {})
        for k in self.initial_bins:
            if not 0 <= k < self.bins.K:
                raise ValueError(f"initial bin {k} outside 0..{self.bins.K - 1}")

        # geometry and spreading
        self.L = float(trail.L_m)
        self.r0 = float(hot.r0_m)
        if not (self.L > 0.0 and self.r0 > 0.0):
            raise ValueError(f"the box needs a positive trail length and initial radius (L = {self.L} m, "
                             f"r0 = {self.r0} m)")
        self.V0 = volume(self.r0, self.L)
        self.K_eddy = self.eddy.K(amb.z_km)
        self._alo = self.species["AlO"]

        # inventory: the vaporised aluminium plus any seeded particle atoms
        self.N_Al = float(trail.N_Al)
        self.inventory = self.N_Al + sum(N * self.bins.n_atoms[k] for k, N in self.initial_bins.items())

        # the network, with each rate constant evaluated once at the ambient temperature
        self._warn = RateWarnings()
        self.all_reactions = list(reactions) if reactions is not None else load_reactions()
        self.reactions = self._select(self.all_reactions)
        self._build_network()

        # scaling (decision 19)
        lay = self.layout
        self.c_amb = np.array([amb.n[s] for s in AMBIENT_SPECIES], dtype=float)
        s_al = max(self.inventory, 1.0)
        scale = np.full(lay.size, s_al)
        scale[lay.ambient] = np.where(self.c_amb > 0.0, self.c_amb, 1.0) * self.V0
        scale[lay.r2] = self.r0 ** 2
        self.scale = scale
        nb = lay.r2.stop                                    # the gas and r2 block is the first nb entries
        rows, cols = np.meshgrid(np.arange(nb), np.arange(nb), indexing="ij")
        self._jac_rows, self._jac_cols = rows.ravel(), cols.ravel()
        self._jac_scale = (scale[:nb, None] ** -1) * scale[None, :nb]

    # ---- set-up ---------------------------------------------------------------------------------------------------

    def _select(self, pool):
        """The active reactions: none without chemistry; Settings.reactions (rids) restricts; R17 by day only."""
        s = self.settings
        if not s.chemistry:
            return []
        if s.reactions is not None:
            known = {r.rid for r in pool}
            missing = [rid for rid in s.reactions if rid not in known]
            if missing:
                raise ValueError(f"Settings.reactions names unknown reaction ids {missing}")
            pool = [r for r in pool if r.rid in s.reactions]
        return [r for r in pool if not (r.day_only and self.time_of_day != "day")]

    def _build_network(self):
        """Index arrays for the reactions: two reactant slots per reaction (a missing one points at a constant 1
        appended to the gas vector), the effective rate constant k [M]^m, the volume exponent 1 - (number of
        reactant species) and the stoichiometric matrix."""
        ng = len(GAS_SPECIES)
        pos = {name: i for i, name in enumerate(GAS_SPECIES)}
        nr = len(self.reactions)
        self._ra = np.full(nr, ng, dtype=int)
        self._rb = np.full(nr, ng, dtype=int)
        self._keff = np.zeros(nr)
        self._vexp = np.zeros(nr)
        self._stoich = np.zeros((ng, nr))
        self.rate_constants = {}
        for j, rx in enumerate(self.reactions):
            k = rx.k(self.amb.T_K, self._warn)
            self.rate_constants[rx.rid] = k
            reac = [p for p in rx.reactants if p != "M"]
            n_m = len(rx.reactants) - len(reac)
            for p in reac + [p for p in rx.products if p != "M"]:
                if p not in pos:
                    raise ValueError(f"{rx.rid}: species {p!r} is not a box species")
            if not 1 <= len(reac) <= 2:
                raise ValueError(f"{rx.rid}: the box handles one or two reactant species, got {reac}")
            self._ra[j] = pos[reac[0]]
            if len(reac) == 2:
                self._rb[j] = pos[reac[1]]
            # [M] is the total air density amb.n_total_m3 (decision 12; equal to the plume's own sum of
            # ambient-species concentrations to trace precision, see the module docstring)
            self._keff[j] = k * self.amb.n_total_m3 ** n_m
            self._vexp[j] = 1.0 - len(reac)
            for p in reac:
                self._stoich[pos[p], j] -= 1.0
            for p in rx.products:
                if p != "M":
                    self._stoich[pos[p], j] += 1.0
            al_in = sum(self.species[p].n_Al for p in reac)
            al_out = sum(self.species[p].n_Al for p in rx.products if p != "M")
            if al_in != al_out:
                raise ValueError(f"{rx.rid} does not conserve aluminium ({al_in} in, {al_out} out)")

    # ---- public helpers -------------------------------------------------------------------------------------------

    @property
    def warnings(self):
        """Rate-clipping messages (one per reaction evaluated outside its validity range)."""
        return self._warn.messages

    def reaction(self, rid):
        """The Reaction object with id rid (from the pool, whether or not it is active in this run)."""
        for rx in self.all_reactions:
            if rx.rid == rid:
                return rx
        raise KeyError(f"no reaction {rid!r}")

    def g(self):
        """Gravitational acceleration (m/s^2) at the plume altitude."""
        return G0 * (R_EARTH_KM / (R_EARTH_KM + self.amb.z_km)) ** 2

    def rho_gas(self):
        """Ambient air mass density (kg/m^3) from the ambient species."""
        return sum(self.amb.n[s] * self.species[s].mass_amu for s in AMBIENT_SPECIES) * AMU

    def settling_rate_of_bin(self, k, r_m=None, h_per_al=0.0):
        """Loss rate (1/s) of bin k by settling out of a plume of radius r_m (default r0): 2 v_s / (pi r) with v_s
        the Epstein speed of a bin-k particle of hydrogen-to-aluminium ratio h_per_al (decision 10 sets its size and
        density: alumina Al2O3 mixed with Al(OH)3 in the share f = min(H/Al / 3, 1))."""
        r_m = self.r0 if r_m is None else float(r_m)
        n = float(self.bins.n_atoms[k])
        f = min(max(h_per_al / 3.0, 0.0), 1.0)
        m_al, m_o, m_h = (self.species[s].mass_amu for s in ("Al", "O", "H"))
        mass_per_al = ((1.0 - f) * (m_al + 1.5 * m_o) + f * (m_al + 3.0 * m_o + 3.0 * m_h)) * AMU
        rho_p = mass_per_al / volume_per_al(h_per_al)
        v_s = epstein_speed(0.5 * diameter(n, h_per_al), rho_p, self.rho_gas(), self.amb.T_K, self.g())
        return settling_rate(v_s, r_m)

    # ---- the ODE --------------------------------------------------------------------------------------------------

    def y0(self):
        """The scaled initial state: the hot-start fractions of the inventory as Al, AlO, OAlO; ambient counts
        c_ambient V0; r^2 = r0^2; seeded bins (initial_bins) tallied to (AlO, formation); everything else zero."""
        lay = self.layout
        x = np.zeros(lay.size)
        total = sum(self.hot.fractions.values())
        if self.N_Al > 0.0 and abs(total - 1.0) > 1e-9:
            raise ValueError(f"the hot-start fractions sum to {total!r}, not 1, so the box would not start from "
                             "the whole inventory")
        for name, frac in self.hot.fractions.items():
            if name not in AL_SPECIES:
                raise ValueError(f"hot-start species {name!r} is not an aluminium-bearing gas species")
            x[lay.index(name)] = frac * self.N_Al
        x[lay.ambient] = self.c_amb * self.V0
        x[lay.r2] = self.r0 ** 2
        tally = x[lay.tally].reshape(lay.K, len(SOURCES), len(ROUTES))
        for k, N in self.initial_bins.items():
            x[lay.N.start + k] = N
            tally[k, SOURCES.index("AlO"), ROUTES.index("formation")] = N * self.bins.n_atoms[k]
        x[lay.tally] = tally.ravel()
        return x / self.scale

    def _geometry(self, r2):
        """(V, dV/dt, d(r^2)/dt) at mean-square radius r2.

        r^2 only grows, so the true solution never has r2 < r0^2; a trial Newton iterate can (BDF taking a long step
        across the eddy switch has produced r2 < 0), and the floor at r0^2 keeps the right-hand side finite so the
        solver rejects the step and shortens it instead of failing on sqrt of a negative number."""
        r = math.sqrt(max(r2, self.r0 * self.r0))
        V = volume(r, self.L)
        if self.settings.freeze_radius:
            return V, 0.0, 0.0
        dr2 = 4.0 * D_total(r, self._alo, self.amb, self.K_eddy, self.settings)
        return V, V * dlnV_dr2(r, self.L) * dr2, dr2

    def _gas_rhs(self, gas, r2):
        """d(counts)/dt of the gas species and d(r^2)/dt, unscaled."""
        V, dV, dr2 = self._geometry(r2)
        ext = np.append(gas, 1.0)
        rates = self._keff * V ** self._vexp * ext[self._ra] * ext[self._rb]
        d = self._stoich @ rates
        d[len(AL_SPECIES):] += self.c_amb * dV
        return d, dr2

    def rhs(self, t, y):
        """Scaled right-hand side. The particle slices are zero in this version (Task 8 adds their terms)."""
        lay = self.layout
        x = y * self.scale
        out = np.zeros(lay.size)
        d, dr2 = self._gas_rhs(x[lay.gas], x[lay.r2][0])
        out[lay.gas] = d
        out[lay.r2] = dr2
        return out / self.scale

    def jac(self, t, y):
        """Scaled Jacobian, sparse: the analytic chemistry block of the gas, and the r^2 column (volume and spreading
        both depend on r^2) by a one-sided difference. Rows and columns of the particle slices are zero here."""
        lay = self.layout
        x = y * self.scale
        gas, r2 = x[lay.gas], x[lay.r2][0]
        ng = len(GAS_SPECIES)
        nb = lay.r2.stop
        V = volume(math.sqrt(max(r2, self.r0 * self.r0)), self.L)
        ext = np.append(gas, 1.0)
        coef = self._keff * V ** self._vexp
        drate = np.zeros((len(self.reactions), ng + 1))             # last column: the constant slot, dropped
        j = np.arange(len(self.reactions))
        np.add.at(drate, (j, self._ra), coef * ext[self._rb])
        np.add.at(drate, (j, self._rb), coef * ext[self._ra])
        J = np.zeros((nb, nb))
        J[:ng, :ng] = self._stoich @ drate[:, :ng]
        h = 1e-7 * r2
        d1, dr2_1 = self._gas_rhs(gas, r2 + h)
        d0, dr2_0 = self._gas_rhs(gas, r2)
        J[:ng, ng] = (d1 - d0) / h
        J[ng, ng] = (dr2_1 - dr2_0) / h
        J *= self._jac_scale
        return sparse.csc_matrix((J.ravel(), (self._jac_rows, self._jac_cols)), shape=(lay.size, lay.size))

    def time_grid(self):
        """Output times: t = 0, then log-spaced from t_first_s to t_end_s at points_per_decade."""
        s = self.settings
        if s.t_end_s <= s.t_first_s:
            return np.array([0.0, float(s.t_end_s)])
        n = max(int(round(math.log10(s.t_end_s / s.t_first_s) * s.points_per_decade)), 1) + 1
        grid = np.logspace(math.log10(s.t_first_s), math.log10(s.t_end_s), n)
        grid[-1] = s.t_end_s                    # exactly the end of the span
        return np.concatenate(([0.0], grid))

    def solve(self):
        """Integrate from 0 to t_end_s with BDF; returns a Solution holding the unscaled state at the output times."""
        grid = self.time_grid()
        res = solve_ivp(self.rhs, (0.0, grid[-1]), self.y0(), method="BDF", t_eval=grid, rtol=self.settings.rtol,
                        atol=ATOL, jac=self.jac)
        if not res.success:
            raise RuntimeError(f"box integration failed: {res.message}")
        return Solution(res.t, res.y.T * self.scale, self)


class Solution:
    """The box state at the output times: t (s), Y (unscaled state, one row per output time: counts, m^2 for r2,
    atoms), and the Box that produced it."""

    def __init__(self, t, Y, box):
        self.t, self.Y, self.box = np.asarray(t), np.asarray(Y), box

    def _u(self, i):
        return self.box.layout.unpack(self.Y[i])

    def r(self, i):
        """Plume radius (m) at output i."""
        return math.sqrt(self._u(i)["r2"])

    def V(self, i):
        """Plume volume (m^3) at output i."""
        return volume(self.r(i), self.box.L)

    def count(self, name, i):
        """Molecules of gas species `name` in the plume at output i."""
        return float(self.Y[i][self.box.layout.index(name)])

    def concentration(self, name, i):
        """Number density (m^-3) of gas species `name` at output i."""
        return self.count(name, i) / self.V(i)

    def al_budget(self, i):
        """Fractions of the aluminium inventory (vaporised plus any seeded particle atoms) in each aluminium-bearing
        gas species and in plume particles, background smoke and settled matter at output i (all zero when the
        inventory is zero)."""
        u = self._u(i)
        inv = self.box.inventory
        f = (1.0 / inv) if inv > 0.0 else 0.0
        out = {name: float(u["gas"][j]) * f for j, name in enumerate(AL_SPECIES)}
        out["particles"] = float(u["tally"].sum()) * f
        out["smoke"] = float(u["smoke"].sum()) * f
        out["settled"] = float(u["settled"].sum()) * f
        return out

    def bins(self, i):
        """The particle bins at output i: N (particles), n_atoms (Al atoms at each pivot), d_m (diameter, m), H_per_Al,
        O_per_Al, and tally[source][route] (Al atoms per bin)."""
        u = self._u(i)
        atoms = u["tally"].sum(axis=(1, 2))
        safe = np.where(atoms > 0.0, atoms, 1.0)
        h = np.where(atoms > 0.0, u["H"] / safe, 0.0)
        o = np.where(atoms > 0.0, u["O"] / safe, 0.0)
        n = self.box.bins.n_atoms
        return {"N": u["N"].copy(), "n_atoms": n.copy(),
                "d_m": np.array([diameter(n[k], h[k]) for k in range(len(n))]),
                "H_per_Al": h, "O_per_Al": o,
                "tally": {s: {r: u["tally"][:, a, b].copy() for b, r in enumerate(ROUTES)}
                          for a, s in enumerate(SOURCES)}}

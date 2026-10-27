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
  * Particles (spec 4.6; decisions 7-11). Sectional bins by aluminium atoms (particles.Bins, pivots n_k). Each bin
    carries its number N_k, its aluminium atoms tallied by source molecule and route (formation or growth), and its
    O and H atoms. Every rate runs on concentrations at the current V and is converted back to counts:
      - Formation: one labelled pathway per pair of aluminium-bearing species (all 28, particles.formation_pathways),
        event rate k_form c_A c_B V (1/2 for a self-pair), k_form = max(k_DD, k_HS) x spin x formation_eff[pid] x
        formation_scale (decision 7). Each event takes one A and one B and puts a 2-atom particle in bin 0, one Al
        tallied to (A, formation) and one to (B, formation), with the O and H the two molecules deliver (decision 9).
      - Growth: species s hits a bin-k particle at fm_kernel(r_s, m_s, r_k, m_k) c_s per particle, times
        sticking[s] x sticking_scale. The grown entity (n_k + 1 atoms) is placed by Bins.place; it carries the bin's
        tallies, O and H per particle plus one Al tallied to (s, growth) and the molecule's O and H; what each
        destination bin receives is in proportion to the aluminium atoms it receives (decision 8).
      - Oxygen uptake (decision 11): O, O2 and O3 collisions add one O atom each at the collision rate x
        oxidation_eff x g(O/Al), g = 1 below 1.485, 0 at 1.5 and above. The gas is not consumed (decision 9: the
        aluminium is a trace). Between 1.485 and 1.5, g is a cubic smoothstep, not the linear ramp decision 11
        names: the uptake is so fast (a 2-atom particle meets about 2e4 O2 per second at 76 km) that an
        under-oxidised bin sits pinned within about 1e-8 of 1.5, and the linear ramp's slope jumps there from
        -67 to 0. BDF keeps its Jacobian across steps, so a Jacobian taken on the ramp made Newton stall once
        oxygen-rich particles (from OAlO2 and Al(OH)2) pushed the bin past 1.5, and the stalled steps passed the
        root-mean-square convergence test of the 1000-component state: bin O/Al drifted to 2-3 at rtol 1e-6
        against 1.54 at rtol 1e-10. The smoothstep's slope is continuous and vanishes at both ends; with it the
        O/Al of every bin holding more than 1e-9 of the inventory agrees with the rtol 1e-10 solution to 1e-2
        at every absolute tolerance tried (76 km day and night, 60 km night, 90 km day), and to 2e-3 with the
        O tolerance below. The O content is a diagnostic (it changes no size, mass or rate other than its own),
        so it is integrated with its own absolute tolerance ATOL_O = 1e-18 on the scaled variable, which keeps
        its Newton updates visible in the convergence norm.
      - Background smoke: species s is lost at uptake_coefficient(s, z, T) x smoke_sticking per molecule (the smoke
        is a fixed background at the ambient Megner density, not a state); the atoms are tallied in smoke[s].
      - Coagulation: free-molecule kernel between every pair of bins (1/2 for a self-pair); the merged entity of
        n_i + n_j atoms is placed by Bins.place and carries both partners' tallies, O and H.
      - Settling: bin k is lost at 2 v_s / (pi r) with r the current plume radius and v_s the Epstein speed; its
        tallies move to settled[source, route] (settled atoms keep their route).
    Particle size and density follow decision 10 from each bin's H/Al (volume per Al atom mixed between alumina and
    Al(OH)3 by f = min(H/Al / 3, 1); mass per Al atom the matching mix of Al2O3 and Al(OH)3 formula masses, so the
    density is 3.0 g/cm3 for alumina and 2.42 g/cm3 for the hydroxide). The tracked O content is a composition
    diagnostic and does not change size or mass. The kernels therefore depend on composition and are evaluated
    (vectorised) at every right-hand-side call from per-run constants; split tables and every rate coefficient that
    does not depend on composition are computed once.
    Switches: Settings.particles=False removes every condensation term (formation, growth, oxygen uptake, smoke,
    coagulation, settling), leaving the gas phase alone; coagulation, settling and smoke each remove their own.
  * Scaling (decision 19). The solver sees scaled variables: aluminium-bearing counts, tallies, O, H, smoke and
    settled atoms divided by the aluminium inventory S (or by 1 when it is 0); the particle number of bin k by
    S / n_k, so a scaled N is the share of the inventory the bin holds and atol applies to its atoms; each ambient
    count by c_ambient V0 (V0 when c_ambient is 0); r^2 by r0^2. BDF runs with rtol = Settings.rtol and atol = 1e-12
    on these, except the O content of each bin, whose absolute tolerance is ATOL_O (oxygen uptake above).
  * Jacobian. Analytic for every term at frozen kernels (the kernels' weak dependence on composition through H/Al is
    left out, which only affects Newton convergence, not accuracy); the r^2 column is a one-sided difference of the
    whole right-hand side, since the volume, the spreading rate and the settling rate all depend on r^2.

Two things are called "reactions": `Settings.reactions` is a tuple of reaction ids that restricts the network (None
keeps all of it), while the `reactions=` keyword of Box is a list of Reaction objects to draw from (default: the box
reactions of data/wake_column/reactions.csv).
"""
import math

import numpy as np
from scipy import sparse
from scipy.integrate import solve_ivp

from reentry_model.wake_column.particles import (V_AL2O3, V_ALOH3, Bins, SmokeBackground, capture_rate,
                                                 diameter, epstein_speed, formation_pathways, spin_factor,
                                                 volume_per_al)
from reentry_model.wake_column.reactions import RateWarnings, load_reactions
from reentry_model.wake_column.species import AL_SPECIES, AMBIENT_SPECIES, AMU, ANGSTROM, load_species
from reentry_model.wake_column.transport import (AIR_MASS_AMU, K_B, D_total, EddyProfile, dlnV_dr2, settling_rate,
                                                 volume)

GAS_SPECIES = AL_SPECIES + AMBIENT_SPECIES
SOURCES = AL_SPECIES                    # the molecule an aluminium atom arrived as
ROUTES = ("formation", "growth")        # how it entered a plume particle (smoke uptake has its own slice)
PROVENANCE_ROUTES = ROUTES + ("smoke",)
OXIDANTS = ("O", "O2", "O3")            # decision 11
G_FULL, G_ZERO = 1.485, 1.5             # O/Al where oxygen uptake starts to slow and where it stops (decision 11)
OX_RAMP = "smoothstep"                  # shape of g between them (see oxidation_g and the module docstring)
NC = len(SOURCES) * len(ROUTES) + 2     # contents carried per bin besides N: 14 tallies, O, H
ATOL = 1e-12                            # absolute tolerance on the scaled variables (decision 19)
ATOL_O = 1e-18                          # absolute tolerance on the scaled O content of each bin (module docstring)
G0 = 9.80665                            # m/s^2 at the surface
R_EARTH_KM = 6371.0


def oxidation_g(x, ramp="smoothstep"):
    """(g, dg/dx) of the oxygen-uptake factor at O/Al ratio x (decision 11): 1 below G_FULL, 0 at G_ZERO and above.
    Between them `ramp` is "linear" (decision 11 as written) or "smoothstep" (3 s^2 - 2 s^3 in s = (G_ZERO - x) /
    (G_ZERO - G_FULL), whose slope is continuous and vanishes at both ends)."""
    s = np.clip((G_ZERO - np.asarray(x, dtype=float)) / (G_ZERO - G_FULL), 0.0, 1.0)
    inside = (s > 0.0) & (s < 1.0)
    if ramp == "linear":
        return s, np.where(inside, -1.0 / (G_ZERO - G_FULL), 0.0)
    if ramp == "smoothstep":
        return s * s * (3.0 - 2.0 * s), np.where(inside, -6.0 * s * (1.0 - s) / (G_ZERO - G_FULL), 0.0)
    raise ValueError(f"unknown oxidation ramp {ramp!r}")


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
        self._build_particles()

        # scaling (decision 19; the particle number of bin k by S / n_k, so atol bounds its atoms)
        lay = self.layout
        self.c_amb = np.array([amb.n[s] for s in AMBIENT_SPECIES], dtype=float)
        s_al = max(self.inventory, 1.0)
        scale = np.full(lay.size, s_al)
        scale[lay.ambient] = np.where(self.c_amb > 0.0, self.c_amb, 1.0) * self.V0
        scale[lay.r2] = self.r0 ** 2
        scale[lay.N] = s_al / self.bins.n_atoms
        self.scale = scale
        self._jac_scale = scale[None, :] / scale[:, None]

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

    def _placement(self, n_atoms):
        """[(bin, number fraction, share of the entity's atoms), ...] for an entity of n_atoms aluminium atoms, from
        Bins.place (so atoms above the last pivot are conserved); the shares sum to 1."""
        return [(d, F, F * self.bins.n_atoms[d] / n_atoms) for d, F in self.bins.place(n_atoms)]

    def _build_particles(self):
        """Every per-run constant of the particle terms (module docstring): formation and smoke uptake as extra
        columns of the network (rate k ext[a] ext[b] V^vexp, stoichiometry over the whole state), the content index
        map, the growth and coagulation placement tables, the molecular radii and masses for the kernels, and the
        settling constants."""
        s, lay, K = self.settings, self.layout, self.bins.K
        sp = self.species
        nsrc = len(SOURCES)
        self._part = bool(s.particles)
        self._coag = self._part and bool(s.coagulation)
        self._settle = self._part and bool(s.settling)
        self._smoke_on = self._part and bool(s.smoke)

        # state index of content c of bin k: tally (source a, route b) is c = 2a + b, then O (c = 14) and H (c = 15)
        idx = np.empty((K, NC), dtype=int)
        idx[:, :NC - 2] = lay.tally.start + np.arange(K * (NC - 2)).reshape(K, NC - 2)
        idx[:, NC - 2] = lay.O.start + np.arange(K)
        idx[:, NC - 1] = lay.H.start + np.arange(K)
        self._idxX = idx
        self._idxN = lay.N.start + np.arange(K)
        self._iform = ROUTES.index("formation")
        self._igrow = ROUTES.index("growth")

        # per-species efficiencies, checked against the names they refer to
        paths = formation_pathways()
        pids = {p.pid for p in paths}
        fe = dict(s.formation_eff)
        unknown = sorted(set(fe) - pids)
        if unknown:
            raise ValueError(f"Settings.formation_eff names unknown pathways {unknown}; known: {sorted(pids)}")
        st = dict(s.sticking)
        unknown = sorted(set(st) - set(AL_SPECIES))
        if unknown:
            raise ValueError(f"Settings.sticking names unknown species {unknown}; known: {list(AL_SPECIES)}")
        self.sticking = np.array([st.get(a, 1.0) * s.sticking_scale for a in AL_SPECIES])

        # what one molecule of each aluminium-bearing species delivers to a particle (decision 9), as contents
        dO = np.array([sp[a].deliver_O for a in AL_SPECIES], dtype=float)
        dH = np.array([sp[a].deliver_H for a in AL_SPECIES], dtype=float)
        self._delta = np.zeros((nsrc, NC))
        for a in range(nsrc):
            self._delta[a, 2 * a + self._igrow] = 1.0
        self._delta[:, NC - 2] = dO
        self._delta[:, NC - 1] = dH

        # formation and smoke uptake as extra network columns, with stoichiometry over the whole state
        T = self.amb.T_K
        ng = len(GAS_SPECIES)
        extra_ra, extra_rb, extra_k, extra_v, cols = [], [], [], [], []
        self.formation_rate_constants = {}
        for p in paths:
            k = (capture_rate(sp[p.a], sp[p.b], T) * spin_factor(sp[p.a], sp[p.b]) * fe.get(p.pid, 1.0)
                 * s.formation_scale)
            self.formation_rate_constants[p.pid] = k
            if not self._part:
                continue
            ia, ib = AL_SPECIES.index(p.a), AL_SPECIES.index(p.b)
            col = np.zeros(lay.size)
            col[lay.gas.start + ia] -= 1.0
            col[lay.gas.start + ib] -= 1.0
            content = np.zeros(NC)
            content[2 * ia + self._iform] += 1.0
            content[2 * ib + self._iform] += 1.0
            content[NC - 2] = dO[ia] + dO[ib]
            content[NC - 1] = dH[ia] + dH[ib]
            for d, F, W in self._placement(2.0):
                col[self._idxN[d]] += F
                col[idx[d]] += W * content
            extra_ra.append(ia); extra_rb.append(ib)
            extra_k.append(k * (0.5 if ia == ib else 1.0))        # self-pair: half the event rate (decision 7)
            extra_v.append(-1.0)                                   # k c_A c_B V = k N_A N_B / V
            cols.append(col)
        self.smoke_uptake = {}
        for a, name in enumerate(AL_SPECIES):
            u = self.smoke.uptake_coefficient(sp[name], self.amb.z_km, T) * s.smoke_sticking
            self.smoke_uptake[name] = u
            if not self._smoke_on:
                continue
            col = np.zeros(lay.size)
            col[lay.gas.start + a] = -1.0
            col[lay.smoke.start + a] = 1.0
            extra_ra.append(a); extra_rb.append(ng); extra_k.append(u); extra_v.append(0.0)
            cols.append(col)
        S = np.zeros((lay.size, len(self.reactions) + len(cols)))
        S[lay.gas, :len(self.reactions)] = self._stoich
        for j, col in enumerate(cols):
            S[:, len(self.reactions) + j] = col
        self._S = S
        self._ra = np.concatenate([self._ra, np.array(extra_ra, dtype=int)])
        self._rb = np.concatenate([self._rb, np.array(extra_rb, dtype=int)])
        self._keff = np.concatenate([self._keff, extra_k])
        self._vexp = np.concatenate([self._vexp, extra_v])

        # growth placement: the entity of n_k + 1 atoms from bin k (P_F number fractions, P_W atom shares)
        n = self.bins.n_atoms
        self._PF = np.zeros((K, K))
        self._PW = np.zeros((K, K))
        for k in range(K):
            for d, F, W in self._placement(n[k] + 1.0):
                self._PF[d, k] += F
                self._PW[d, k] += W
        # coagulation placement: entries (m, i, j, F, W) over ordered pairs, the entity of n_i + n_j atoms
        em, ei, ej, eF, eW = [], [], [], [], []
        for i in range(K):
            for j in range(K):
                for d, F, W in self._placement(n[i] + n[j]):
                    em.append(d); ei.append(i); ej.append(j); eF.append(F); eW.append(W)
        self._em, self._ei, self._ej = (np.array(v, dtype=int) for v in (em, ei, ej))
        self._eF, self._eW = np.array(eF), np.array(eW)
        self._emi = self._em * K + self._ei
        self._emj = self._em * K + self._ej

        # kernel constants: molecules (radius sigma/2, mass), per-Al particle masses (decision 10), sqrt(8 k T / pi)
        self._cT = math.sqrt(8.0 * K_B * T / math.pi)
        self._r_al = np.array([0.5 * sp[a].sigma_A * ANGSTROM for a in AL_SPECIES])
        self._m_al = np.array([sp[a].mass_amu * AMU for a in AL_SPECIES])
        self._ox_idx = np.array([lay.index(g) for g in OXIDANTS])
        self._r_ox = np.array([0.5 * sp[g].sigma_A * ANGSTROM for g in OXIDANTS])
        self._m_ox = np.array([sp[g].mass_amu * AMU for g in OXIDANTS])
        m_al, m_o, m_h = (sp[e].mass_amu * AMU for e in ("Al", "O", "H"))
        self._mpa_alumina = m_al + 1.5 * m_o                 # Al2O3 per Al atom
        self._mpa_hydroxide = m_al + 3.0 * m_o + 3.0 * m_h    # Al(OH)3 per Al atom
        self._settle_const = self.g() / (self.rho_gas() * math.sqrt(8.0 * K_B * T / (math.pi * AIR_MASS_AMU * AMU))
                                         * (1.0 + math.pi / 8.0))

    def _particle_props(self, f):
        """(radius m, mass kg, density kg/m3) of a particle at each pivot for hydroxide share f (decision 10)."""
        f = np.clip(f, 0.0, 1.0)
        vpa = (1.0 - f) * V_AL2O3 + f * V_ALOH3
        mpa = (1.0 - f) * self._mpa_alumina + f * self._mpa_hydroxide
        n = self.bins.n_atoms
        return (3.0 * n * vpa / (4.0 * math.pi)) ** (1.0 / 3.0), n * mpa, mpa / vpa

    def _kernel(self, r1, m1, r2, m2):
        """Free-molecule kernel pi (r1 + r2)^2 sqrt(8 k T / (pi mu)) (m3/s), broadcasting r1, m1 against r2, m2."""
        return math.pi * (r1 + r2) ** 2 * self._cT * np.sqrt(1.0 / m1 + 1.0 / m2)

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
        c_ambient V0; r^2 = r0^2; seeded bins (initial_bins) tallied to (AlO, formation) with one O per Al atom and
        no H (what AlO delivers, decision 9); everything else zero."""
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
            x[lay.O.start + k] = N * self.bins.n_atoms[k]          # AlO delivers one O and no H (decision 9)
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

    def _network_rates(self, gas, V):
        """Event rates (1/s, whole plume) of the network columns: reactions, formation pathways, smoke uptake."""
        ext = np.append(gas, 1.0)
        return self._keff * V ** self._vexp * ext[self._ra] * ext[self._rb]

    def _pstate(self, x, V, r):
        """The particle quantities the right-hand side and the Jacobian share, at unscaled state x, volume V and
        plume radius r: bin numbers and contents, composition, per-particle growth rates by species, oxygen
        uptake, and (when on) the coagulation kernel and the settling rates."""
        lay, n = self.layout, self.bins.n_atoms
        gas = x[lay.gas]
        N, X = x[self._idxN], x[self._idxX]
        al = N * n
        pos = al > 0.0
        f = np.zeros_like(al)
        np.divide(X[:, NC - 1], 3.0 * al, out=f, where=pos)
        rp, mp, rho_p = self._particle_props(f)
        p = {"N": N, "X": X, "gas": gas, "V": V}
        # growth: B[s, k] = kernel x sticking / V, so B[s, k] N_s is the rate at which one bin-k particle gains atoms
        B = self._kernel(self._r_al[:, None], self._m_al[:, None], rp[None, :], mp[None, :]) * self.sticking[:, None]
        p["B"] = B / V
        p["rho_ks"] = p["B"] * gas[:len(AL_SPECIES), None]
        p["rho"] = p["rho_ks"].sum(axis=0)
        # oxygen uptake (decision 11): Z per particle at g = 1, g(O/Al) and its slope
        kap = self._kernel(self._r_ox[:, None], self._m_ox[:, None], rp[None, :], mp[None, :])
        p["kap"] = kap * (self.settings.oxidation_eff / V)
        p["Z"] = gas[self._ox_idx] @ p["kap"]
        xo = np.full_like(al, G_ZERO)
        np.divide(X[:, NC - 2], al, out=xo, where=pos)
        p["xo"] = xo
        g, gd = oxidation_g(xo, OX_RAMP)
        p["gval"] = np.where(pos, g, 0.0)
        p["gder"] = np.where(pos, gd, 0.0)
        if self._coag:
            Kc = self._kernel(rp[:, None], mp[:, None], rp[None, :], mp[None, :]) / V
            kij = Kc[self._ei, self._ej]
            K = self.bins.K
            p["Kc"], p["kij"] = Kc, kij
            p["KN"] = Kc @ N
            p["M"] = np.bincount(self._emi, weights=self._eW * kij * N[self._ej], minlength=K * K).reshape(K, K)
            p["MF"] = np.bincount(self._emi, weights=self._eF * kij * N[self._ej], minlength=K * K).reshape(K, K)
        if self._settle:
            p["lam"] = 2.0 * rho_p * rp * self._settle_const / (math.pi * r)     # 2 v_s / (pi r), Epstein v_s
        return p

    def _rhs_x(self, x):
        """Unscaled right-hand side d(state)/dt at unscaled state x."""
        lay = self.layout
        gas, r2 = x[lay.gas], x[lay.r2][0]
        V, dV, dr2 = self._geometry(r2)
        out = self._S @ self._network_rates(gas, V)
        out[lay.ambient] += self.c_amb * dV
        out[lay.r2] = dr2
        if not self._part:
            return out
        p = self._pstate(x, V, math.sqrt(max(r2, self.r0 * self.r0)))
        N, X, rho = p["N"], p["X"], p["rho"]
        # growth: entities of n_k + 1 atoms carry the bin's contents per particle plus what the molecule delivers
        fN = rho * N
        dN = self._PF @ fN - fN
        D = p["rho_ks"].T @ self._delta
        dX = self._PW @ (rho[:, None] * X + N[:, None] * D) - rho[:, None] * X
        out[lay.al] -= p["rho_ks"] @ N
        # oxygen uptake
        dX[:, NC - 2] += N * p["Z"] * p["gval"]
        # coagulation: loss_i = content_i sum_j K_ij N_j / V; the merged contents are placed by atom share
        if self._coag:
            dN += 0.5 * (p["MF"] @ N) - N * p["KN"]
            dX += p["M"] @ X - X * p["KN"][:, None]
        # settling: the bin's tallies move to settled[source, route]
        if self._settle:
            lam = p["lam"]
            dN -= lam * N
            dX -= lam[:, None] * X
            out[lay.settled] += lam @ X[:, :NC - 2]
        out[self._idxN] += dN
        out[self._idxX] += dX
        return out

    def rhs(self, t, y):
        """Scaled right-hand side."""
        return self._rhs_x(y * self.scale) / self.scale

    def jac(self, t, y):
        """Scaled Jacobian (sparse): analytic for every term at frozen kernels, and the r^2 column by a one-sided
        difference of the whole right-hand side (module docstring)."""
        lay, K = self.layout, self.bins.K
        x = y * self.scale
        gas, r2 = x[lay.gas], x[lay.r2][0]
        ng = len(GAS_SPECIES)
        r = math.sqrt(max(r2, self.r0 * self.r0))
        V = volume(r, self.L)
        J = np.zeros((lay.size, lay.size))
        # network columns (reactions, formation, smoke): rate = coef ext[a] ext[b]
        ext = np.append(gas, 1.0)
        coef = self._keff * V ** self._vexp
        nr = len(self._keff)
        drate = np.zeros((nr, ng + 1))                              # last column: the constant slot, dropped
        j = np.arange(nr)
        np.add.at(drate, (j, self._ra), coef * ext[self._rb])
        np.add.at(drate, (j, self._rb), coef * ext[self._ra])
        J[:, :ng] = self._S @ drate[:, :ng]
        if self._part:
            p = self._pstate(x, V, r)
            N, X, rho, B = p["N"], p["X"], p["rho"], p["B"]
            iN, iX = self._idxN, self._idxX
            eye = np.eye(K)
            nal = len(AL_SPECIES)
            gal = lay.gas.start + np.arange(nal)
            # growth
            JNN = (self._PF - eye) * rho[None, :]
            JXX = (self._PW - eye) * rho[None, :]                   # the same for every content
            D = p["rho_ks"].T @ self._delta
            JXN = [self._PW * D[:, c][None, :] for c in range(NC)]
            J[np.ix_(iN, gal)] += (self._PF - eye) @ (B * N[None, :]).T
            T1 = np.einsum("kj,sj,jc->ksc", self._PW - eye, B, X)
            T2 = (self._PW @ (B * N[None, :]).T)[:, :, None] * self._delta[None, :, :]
            TG = T1 + T2
            for c in range(NC):
                J[np.ix_(iX[:, c], gal)] += TG[:, :, c]
            J[gal, gal] -= B @ N
            J[np.ix_(gal, iN)] -= B * gas[:nal, None]
            # oxygen uptake
            iO = iX[:, NC - 2]
            Z, gval, gder, xo = p["Z"], p["gval"], p["gder"], p["xo"]
            J[iO, iO] += Z * gder / self.bins.n_atoms
            J[iO, iN] += Z * (gval - xo * gder)
            J[np.ix_(iO, self._ox_idx)] += (N * gval)[:, None] * p["kap"].T
            # coagulation
            if self._coag:
                Kc, KN = p["Kc"], p["KN"]
                JNN = JNN + p["MF"] - np.diag(KN) - N[:, None] * Kc
                JXX = JXX + p["M"] - np.diag(KN)
                H = sparse.csr_matrix((self._eW * p["kij"], (self._emj, self._ei)), shape=(K * K, K))
                T = (H @ X).reshape(K, K, NC)                       # T[m, j, c]: d gain_X[m, c] / d N_j
                JXN = [JXN[c] + T[:, :, c] - X[:, c][:, None] * Kc for c in range(NC)]
            # settling
            if self._settle:
                lam = p["lam"]
                JNN = JNN - np.diag(lam)
                JXX = JXX - np.diag(lam)
                for c in range(NC - 2):
                    J[lay.settled.start + c, iX[:, c]] += lam
            J[np.ix_(iN, iN)] += JNN
            for c in range(NC):
                J[np.ix_(iX[:, c], iX[:, c])] += JXX
                J[np.ix_(iX[:, c], iN)] += JXN[c]
        # the r^2 column: volume, spreading and settling all depend on it
        ir2 = lay.r2.start
        h = 1e-7 * max(r2, self.r0 * self.r0)
        xp = x.copy()
        xp[ir2] += h
        J[:, ir2] = (self._rhs_x(xp) - self._rhs_x(x)) / h
        J *= self._jac_scale
        return sparse.csc_matrix(J)

    def time_grid(self):
        """Output times: t = 0, then log-spaced from t_first_s to t_end_s at points_per_decade."""
        s = self.settings
        if s.t_end_s <= s.t_first_s:
            return np.array([0.0, float(s.t_end_s)])
        n = max(int(round(math.log10(s.t_end_s / s.t_first_s) * s.points_per_decade)), 1) + 1
        grid = np.logspace(math.log10(s.t_first_s), math.log10(s.t_end_s), n)
        grid[-1] = s.t_end_s                    # exactly the end of the span
        return np.concatenate(([0.0], grid))

    def atol(self):
        """Absolute tolerances on the scaled state: ATOL everywhere, ATOL_O on the O content of each bin."""
        a = np.full(self.layout.size, ATOL)
        a[self._idxX[:, NC - 2]] = ATOL_O
        return a

    def solve(self):
        """Integrate from 0 to t_end_s with BDF; returns a Solution holding the unscaled state at the output times."""
        grid = self.time_grid()
        res = solve_ivp(self.rhs, (0.0, grid[-1]), self.y0(), method="BDF", t_eval=grid, rtol=self.settings.rtol,
                        atol=self.atol(), jac=self.jac)
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

    def provenance(self, i):
        """{(source, route): fraction of the inventory} at output i for every source molecule and every route
        ("formation", "growth", "smoke"; all 21 keys always present). Settled atoms count under the route by which
        they entered their particle, so the values sum to particles + smoke + settled of al_budget."""
        u = self._u(i)
        inv = self.box.inventory
        f = (1.0 / inv) if inv > 0.0 else 0.0
        in_plume = u["tally"].sum(axis=0) + u["settled"]            # (source, route)
        out = {}
        for a, s in enumerate(SOURCES):
            for b, route in enumerate(ROUTES):
                out[(s, route)] = float(in_plume[a, b]) * f
            out[(s, "smoke")] = float(u["smoke"][a]) * f
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

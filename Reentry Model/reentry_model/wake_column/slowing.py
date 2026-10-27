"""Stage A of the wake column: the slowing-down Monte Carlo of the hot aluminium vapour (spec section 4.4,
decision 15).

Aluminium atoms leave the droplet with its speed relative to the still air (1.5-7.6 km/s) plus an isotropic thermal
velocity at the droplet temperature, so their first collisions are hyperthermal and can open threshold channels that
are closed in ambient air (AlO + O2 -> OAlO + O above about 109 kJ/mol). Test particles are followed collision by
collision until they thermalise; their species and perpendicular distance from the droplet track at that moment give
the cylinder's starting composition and radius r0.

  * Free flight uses the null-collision method with the majorant rate sum_j n_j pi sigma_j^2 (|v| + 4 vbar_j) over
    the partners N2 (the remainder of the air), O2, O and H2O. A candidate partner is drawn by its share of the
    majorant, its velocity from the Maxwellian at the ambient temperature (or zero, `partners_at_rest`), and the
    collision is accepted with probability n_j pi sigma_j^2 |v - u| / w_j, where w_j = n_j pi sigma_j^2
    (|v| + 4 vbar_j) is that partner's term of the majorant; the product of the two draws gives partner j its true
    collision rate n_j pi sigma_j^2 |v - u|. (Dividing by the whole majorant instead would weight each partner by
    its share twice and undercount the minor partners O2, O and H2O.)
  * Collisions are hard-sphere and isotropic in the centre-of-mass frame. Before scattering, every reaction open to
    the tracked species and the partner reacts with the line-of-centers probability
    hot_reaction_scale * P0 (1 - E0/E) for E > E0 (E the centre-of-mass collision energy), the sum renormalised to 1
    when it exceeds 1. A product keeps the reactant's velocity and takes its own mass; the partner-side products
    are not tracked, and a reacting collision does not also scatter.
  * A particle is thermal once 1/2 m |v|^2 (in the air's frame) < thermal_cutoff * 3/2 k_B T_ambient; its species,
    sqrt(y^2 + z^2) (x is along the track), collision count and elapsed time are recorded then. Particles emitted
    below the cutoff stop at the droplet with no collision.

E0 and P0 are calibrated to the thermal rate expressions of the network: the line-of-centers thermal rate is exactly
P0 pi sigma12^2 vbar(T) exp(-E0 / k_B T), so E0 = k_B Ea and P0 = k(T_ref) exp(Ea/T_ref) / (pi sigma12^2 vbar(T_ref))
at the top of each reaction's validity range (decision 14), capped at 1. The threshold form (collisional
dissociation AlO + M) takes E0 and P0 from its row.
"""
from dataclasses import dataclass

import numpy as np
from scipy.constants import N_A, k as K_B

from reentry_model.wake_column.species import AMU, ANGSTROM

PARTNERS = ("N2", "O2", "O", "H2O")          # N2 stands for the remainder of the air (spec section 4.4)
TRACKED = ("Al", "AlO", "OAlO")              # the species a test particle can be
AL_IONISATION_J = 578e3 / N_A                # first ionisation energy of Al (5.986 eV, NIST ASD; spec 4.4), per atom
ALO_BOND_DEFAULT_J = 498e3 / N_A             # AlO bond energy (Gobrecht et al. 2022), used when no AlO + M row
MAX_STEPS = 200000                           # free-flight steps before the run is declared stuck


@dataclass(frozen=True)
class HotChannel:
    """One reaction of the hot stage seen from the tracked aluminium-bearing particle: `reactant` + `partner`
    ("M" = any partner) -> `products`; threshold E0_J and probability P0 above it, per collision (energies per
    molecule, J)."""
    rid: str
    reactant: str
    partner: str
    products: tuple
    E0_J: float
    P0: float

    def product(self, species):
        """The aluminium-bearing product, which the test particle becomes."""
        return next(p for p in self.products if p != "M" and species[p].n_Al > 0)


@dataclass(frozen=True)
class HotStartResult:
    fractions: dict                 # {"Al", "AlO", "OAlO"} -> share of the inventory thermalised as that species
    r0_m: float                     # weighted r.m.s. perpendicular distance of the thermalisation points, m
    mean_collisions: float          # real collisions per test particle before it thermalised
    mean_stop_time_s: float         # time from emission to thermalisation, s
    share_above_ionisation: float   # share of Al-atom collisions above the Al ionisation energy
    share_above_AlO_bond: float     # share of AlO collisions above the AlO bond energy
    per_sample: tuple               # one dict per emission sample with the same quantities, the sample's inputs and
                                    # the real-collision counts per partner (collisions_by_partner)


def _mu(ma_amu, mb_amu):
    return ma_amu * mb_amu / (ma_amu + mb_amu) * AMU


def _sigma12(a, b, species):
    return 0.5 * (species[a].sigma_A + species[b].sigma_A) * ANGSTROM


def _vbar(T, mu_kg):
    return np.sqrt(8.0 * K_B * T / (np.pi * mu_kg))


def calibrate(reactions, species):
    """The hot-stage channels of `reactions`: line-of-centers E0 and P0 for each bimolecular row, from the row's
    threshold for form "threshold". Rows with three or more reactants (three-body) and photolysis rows are skipped."""
    out = []
    for r in reactions:
        if r.order >= 3 or r.form == "photolysis":
            continue
        al = [s for s in r.reactants if s != "M" and species[s].n_Al > 0]
        if len(al) != 1:
            raise ValueError(f"{r.rid}: a hot-stage reaction needs exactly one aluminium-bearing reactant, "
                             f"got {r.reactants}")
        reactant = al[0]
        rest = list(r.reactants)
        rest.remove(reactant)
        partner = rest[0]
        if r.form == "threshold":
            E0 = r.E0_kJmol * 1e3 / N_A
            P0 = min(float(r.P0), 1.0)
        else:
            if partner == "M":
                raise ValueError(f"{r.rid}: a thermal rate with partner M cannot be calibrated per collision")
            Ea = r.Ea_K if np.isfinite(r.Ea_K) else 0.0
            T_ref = r.T_max
            mu = _mu(species[reactant].mass_amu, species[partner].mass_amu)
            sig = _sigma12(reactant, partner, species)
            E0 = Ea * K_B
            P0 = min(r.k(T_ref) * np.exp(Ea / T_ref) / (np.pi * sig ** 2 * _vbar(T_ref, mu)), 1.0)
        out.append(HotChannel(r.rid, reactant, partner, tuple(r.products), float(E0), float(P0)))
    return out


def loc_rate(ch, T, species, partner=None):
    """Analytic line-of-centers thermal rate P0 pi sigma12^2 vbar(T) exp(-E0 / k_B T) in m^3/s (`partner` names
    the collider when the channel's partner is M)."""
    b = partner or ch.partner
    mu = _mu(species[ch.reactant].mass_amu, species[b].mass_amu)
    return float(ch.P0 * np.pi * _sigma12(ch.reactant, b, species) ** 2 * _vbar(T, mu) * np.exp(-ch.E0_J / (K_B * T)))


def sample_rate(ch, T, species, n, rng, partner=None):
    """The same thermal average estimated from n relative velocities g drawn from the Maxwellian of the pair's
    reduced mass: mean of sigma_R(E) g with sigma_R = P0 pi sigma12^2 (1 - E0/E) for E = mu g^2 / 2 > E0. m^3/s."""
    b = partner or ch.partner
    mu = _mu(species[ch.reactant].mass_amu, species[b].mass_amu)
    g = np.linalg.norm(rng.normal(0.0, np.sqrt(K_B * T / mu), size=(n, 3)), axis=1)
    E = 0.5 * mu * g * g
    sig_R = ch.P0 * np.pi * _sigma12(ch.reactant, b, species) ** 2 * _p_shape(E, ch.E0_J)
    return float(np.mean(sig_R * g))


def _p_shape(E, E0):
    """(1 - E0/E) above the threshold and 0 below it; 1 everywhere for a barrierless channel (E0 = 0)."""
    if E0 <= 0.0:
        return np.ones_like(E)
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(E > E0, 1.0 - E0 / np.where(E > 0, E, 1.0), 0.0)


def _unit_vectors(rng, n):
    c = rng.uniform(-1.0, 1.0, n)
    phi = rng.uniform(0.0, 2.0 * np.pi, n)
    s = np.sqrt(1.0 - c * c)
    return np.column_stack((s * np.cos(phi), s * np.sin(phi), c))


def _partner_densities(a):
    """Densities [m^-3] of N2 (the remainder of the air), O2, O and H2O at ambient state a."""
    o2, o, h2o = (max(float(a.n[s]), 0.0) for s in ("O2", "O", "H2O"))
    return np.array([max(float(a.n_total_m3) - o2 - o - h2o, 0.0), o2, o, h2o])


def slow_down(samples, column, settings, species, channels, partners_at_rest=False, reactions_on=True):
    """Run Stage A over the emission samples and return the thermalised composition and r0 (HotStartResult)."""
    samples = tuple(samples)
    ns = len(samples)
    if ns == 0:
        raise ValueError("Stage A needs at least one emission sample (no aluminium vaporised?)")
    if settings.n_test_particles < ns:
        raise ValueError(f"n_test_particles = {settings.n_test_particles} is fewer than the {ns} emission samples")
    w = np.array([s.weight for s in samples], dtype=float)
    if not (np.all(w >= 0) and w.sum() > 0):
        raise ValueError(f"emission-sample weights must be non-negative with a positive sum, got {w}")

    rng = np.random.default_rng(settings.seed)
    nsp = len(TRACKED)
    m_sp = np.array([species[s].mass_amu for s in TRACKED]) * AMU
    m_pt = np.array([species[p].mass_amu for p in PARTNERS]) * AMU
    # pi sigma12^2 [m^2] for tracked species x partner
    xs = np.array([[np.pi * _sigma12(s, p, species) ** 2 for p in PARTNERS] for s in TRACKED])
    mu_sp_pt = m_sp[:, None] * m_pt[None, :] / (m_sp[:, None] + m_pt[None, :])

    # particles: equal split over the samples (the first ones take the remainder)
    counts = np.full(ns, settings.n_test_particles // ns)
    counts[: settings.n_test_particles % ns] += 1
    sid = np.repeat(np.arange(ns), counts)
    N = sid.size
    ambients = [column.at(s.h_km) for s in samples]
    dens = np.array([_partner_densities(a) for a in ambients])          # (ns, 4)
    for i, a in enumerate(ambients):
        if not dens[i].sum() > 0:
            raise ValueError(f"no collision partners at {samples[i].h_km} km: the air density is zero")
    T_amb = np.array([a.T_K for a in ambients])
    vbar_pt = np.sqrt(8.0 * K_B * T_amb[:, None] / (np.pi * m_pt[None, :]))   # (ns, 4)
    s_pt = np.sqrt(K_B * T_amb[:, None] / m_pt[None, :])                        # Maxwellian component spread

    # emission: isotropic thermal velocity, speed from the flux-weighted Maxwellian at T_d, plus V along +x
    V = np.array([0.0 if settings.hot_start == "droplet_temperature" else s.V_ms for s in samples])
    Td = np.array([s.T_d_K for s in samples])
    sp = np.zeros(N, dtype=int)                                         # every particle leaves as Al
    speed = np.sqrt(K_B * Td[sid] / m_sp[0]) * np.sqrt(2.0 * rng.gamma(2.0, 1.0, N))
    vel = speed[:, None] * _unit_vectors(rng, N)
    vel[:, 0] += V[sid]
    pos = np.zeros((N, 3))
    t = np.zeros(N)
    ncoll = np.zeros(N, dtype=np.int64)
    E_stop = settings.thermal_cutoff * 1.5 * K_B * T_amb                 # per sample

    # channels as arrays: reactant index, partner index (-1 = any), product index, E0, P0
    ch_list = [c for c in channels if c.reactant in TRACKED] if reactions_on else []
    ch_r = np.array([TRACKED.index(c.reactant) for c in ch_list], dtype=int)
    ch_p = np.array([-1 if c.partner == "M" else (PARTNERS.index(c.partner) if c.partner in PARTNERS else -2)
                     for c in ch_list], dtype=int)
    ch_prod = np.array([TRACKED.index(c.product(species)) for c in ch_list], dtype=int)
    scale = settings.hot_reaction_scale
    bond = next((c.E0_J for c in channels if c.reactant == "AlO" and c.partner == "M"), ALO_BOND_DEFAULT_J)

    # per-sample tallies
    n_al_coll = np.zeros(ns)
    n_al_ion = np.zeros(ns)
    n_alo_coll = np.zeros(ns)
    n_alo_bond = np.zeros(ns)
    loss_sum = np.zeros(ns)
    loss_n = np.zeros(ns)
    n_by_partner = np.zeros((ns, len(PARTNERS)))

    active = 0.5 * m_sp[sp] * np.einsum("ij,ij->i", vel, vel) >= E_stop[sid]
    idx = np.flatnonzero(active)
    steps = 0
    while idx.size:
        steps += 1
        if steps > MAX_STEPS:
            raise RuntimeError(f"Stage A did not thermalise {idx.size} test particles in {MAX_STEPS} steps")
        k = sid[idx]
        s_i = sp[idx]
        v = vel[idx]
        vmag = np.sqrt(np.einsum("ij,ij->i", v, v))
        nsig = dens[k] * xs[s_i]                                        # (n, 4) n_j pi sigma_j^2
        wj = nsig * (vmag[:, None] + vbar_pt[k])
        maj = wj.sum(axis=1)
        dt = -np.log1p(-rng.random(idx.size)) / maj                     # maj > 0: checked above, vbar > 0
        pos[idx] += v * dt[:, None]
        t[idx] += dt
        # candidate partner by its share of the majorant
        cum = np.cumsum(wj, axis=1)
        j = np.minimum((cum < (rng.random(idx.size) * maj)[:, None]).sum(axis=1), len(PARTNERS) - 1)
        if partners_at_rest:
            u = np.zeros_like(v)
        else:
            u = rng.normal(size=v.shape) * s_pt[k, j][:, None]
        g = v - u
        gmag = np.sqrt(np.einsum("ij,ij->i", g, g))
        # accept with the candidate's own majorant: n_j pi sigma_j^2 |v - u| / (n_j pi sigma_j^2 (|v| + 4 vbar_j)),
        # so that, with the candidate drawn by w_j / maj, partner j collides at the true rate n_j pi sigma_j^2 |g|
        rows = np.arange(idx.size)
        acc = rng.random(idx.size) * wj[rows, j] < nsig[rows, j] * gmag
        if not acc.any():
            continue
        a_idx, k, s_i, j, v, u, g, gmag = (x[acc] for x in (idx, k, s_i, j, v, u, g, gmag))
        ncoll[a_idx] += 1
        np.add.at(n_by_partner, (k, j), 1.0)
        mu = mu_sp_pt[s_i, j]
        E = 0.5 * mu * gmag * gmag
        is_al, is_alo = s_i == 0, s_i == 1
        np.add.at(n_al_coll, k[is_al], 1.0)
        np.add.at(n_al_ion, k[is_al & (E > AL_IONISATION_J)], 1.0)
        np.add.at(n_alo_coll, k[is_alo], 1.0)
        np.add.at(n_alo_bond, k[is_alo & (E > bond)], 1.0)

        # reactions, decided before scattering
        reacted = np.zeros(a_idx.size, dtype=bool)
        if len(ch_list):
            P = np.zeros((len(ch_list), a_idx.size))
            for c, ch in enumerate(ch_list):
                m = (s_i == ch_r[c]) & ((ch_p[c] == -1) | (j == ch_p[c]))
                if m.any():
                    P[c, m] = scale * ch.P0 * _p_shape(E[m], ch.E0_J)
            tot = P.sum(axis=0)
            over = tot > 1.0
            P[:, over] /= tot[over]
            cumP = np.cumsum(P, axis=0)
            r = rng.random(a_idx.size)
            hit = r < cumP[-1]
            choice = (cumP <= r[None, :]).sum(axis=0)
            reacted = hit
            if hit.any():
                sp[a_idx[hit]] = ch_prod[choice[hit]]
        # elastic hard-sphere scattering, isotropic in the centre-of-mass frame
        el = ~reacted
        if el.any():
            m1 = m_sp[s_i[el]]
            m2 = m_pt[j[el]]
            vcm = (m1[:, None] * v[el] + m2[:, None] * u[el]) / (m1 + m2)[:, None]
            v_new = vcm + (m2 / (m1 + m2) * gmag[el])[:, None] * _unit_vectors(rng, int(el.sum()))
            if partners_at_rest:
                is_n2 = j[el] == 0
                e0 = np.einsum("ij,ij->i", v[el], v[el])
                e1 = np.einsum("ij,ij->i", v_new, v_new)
                np.add.at(loss_sum, k[el][is_n2], (1.0 - e1 / e0)[is_n2])
                np.add.at(loss_n, k[el][is_n2], 1.0)
            vel[a_idx[el]] = v_new
        # thermalised?
        e_kin = 0.5 * m_sp[sp[a_idx]] * np.einsum("ij,ij->i", vel[a_idx], vel[a_idx])
        stopped = a_idx[e_kin < E_stop[sid[a_idx]]]
        if stopped.size:
            idx = np.setdiff1d(idx, stopped, assume_unique=True)

    # per-sample results and their weighted combination
    rperp2 = pos[:, 1] ** 2 + pos[:, 2] ** 2
    share = lambda a, b: float(a / b) if b > 0 else 0.0
    per = []
    F = np.zeros((ns, nsp))
    for i in range(ns):
        m = sid == i
        F[i] = np.bincount(sp[m], minlength=nsp) / m.sum()
        d = {"h_km": samples[i].h_km, "V_ms": float(V[i]), "T_d_K": samples[i].T_d_K, "T_amb_K": float(T_amb[i]),
             "weight": samples[i].weight, "n_particles": int(m.sum()),
             "fractions": {s: float(F[i, q]) for q, s in enumerate(TRACKED)},
             "r0_m": float(np.sqrt(rperp2[m].mean())), "mean_collisions": float(ncoll[m].mean()),
             "mean_stop_time_s": float(t[m].mean()),
             "share_above_ionisation": share(n_al_ion[i], n_al_coll[i]),
             "share_above_AlO_bond": share(n_alo_bond[i], n_alo_coll[i]),
             "collisions_by_partner": {p: int(n_by_partner[i, q]) for q, p in enumerate(PARTNERS)}}
        if partners_at_rest:
            d["mean_fractional_loss_vs_N2"] = share(loss_sum[i], loss_n[i]) if loss_n[i] else float("nan")
        per.append(d)

    wn = w / w.sum()
    fr = wn @ F
    fr = fr / fr.sum()
    pw = wn / counts                                                     # weight of one particle of each sample
    return HotStartResult(
        fractions={s: float(fr[q]) for q, s in enumerate(TRACKED)},
        r0_m=float(np.sqrt(sum(wn[i] * per[i]["r0_m"] ** 2 for i in range(ns)))),
        mean_collisions=float(wn @ [p["mean_collisions"] for p in per]),
        mean_stop_time_s=float(wn @ [p["mean_stop_time_s"] for p in per]),
        share_above_ionisation=share(pw @ n_al_ion, pw @ n_al_coll),
        share_above_AlO_bond=share(pw @ n_alo_bond, pw @ n_alo_coll),
        per_sample=tuple(per))

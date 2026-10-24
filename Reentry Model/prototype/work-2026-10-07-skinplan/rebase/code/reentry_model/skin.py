"""One-dimensional skins under the melting surface patches (Step 3 sub-plan 18; design:
docs/superpowers/specs/2026-10-07-melt-layer-skin-design.md, sections 5-7).

A skin belongs to one surface patch and has its area: a slab of `n_cells` equal cells along the patch's inward normal
that owns the outermost metal of that patch. Its state is each cell's metal mass `m` [kg] and metal enthalpy `H`
[J above material.T_REF], with the temperature `T` that follows from them. The liquid riding on the patch -- the film
and the deep account, which MeltingBody keeps as `m_f` and `m_d` -- shares the top cell's temperature: cell 0 and that
liquid are one node, of energy m_0 h(T_0) + L h_l(T_0). Every operation conserves mass and enthalpy exactly. Every
array is per skin, and all skins are solved together."""
import math
from dataclasses import dataclass

import numpy as np

from .thermal import SIGMA_SB

SKIN_PRESETS = {"production": 10.0e-6, "dev": 20.0e-6}  # m, the presets' cell sizes (spec section 11.1)
INTERFACE_TOLERANCE = 1.0e-3     # |mismatch| / (|E_top| + |E_base|) above which a macro step repeats (spec 6.2)
NEWTON_TOL = 1.0e-6              # K: a sub-step has converged when no cell moves by more than this
NEWTON_MAX = 40
LINE_SEARCH_MAX = 12           # halvings of a Newton step before it is taken anyway
DZ_FLOOR = 1.0e-12               # m: an emptied cell conducts as if this thin (it holds no heat)


@dataclass
class SkinSettings:
    thickness: float = 0.4e-3              # m, target thickness s
    cell: float = SKIN_PRESETS["production"]
    substep: float = 0.010                 # s, the skins' sub-step

    def __post_init__(self):
        if not (self.thickness > 0.0 and self.cell > 0.0 and self.substep > 0.0):
            raise ValueError("skin thickness, cell and substep must be positive")
        if 2.0 * self.cell > self.thickness * (1.0 + 1e-12):
            raise ValueError("a skin needs at least two cells: cell <= thickness / 2")

    @property
    def n_cells(self):
        return int(round(self.thickness / self.cell))


@dataclass
class SkinBook:
    """What a real pass did, per skin: energies [J] and masses [kg]."""
    E_top: np.ndarray        # convective heat taken in at the top
    E_rad: np.ndarray        # radiated from the top
    E_base: np.ndarray       # passed down into the 3D model through the base
    to_film: np.ndarray      # liquid handed to the film
    to_deep: np.ndarray      # liquid handed to the deep account
    drawn: np.ndarray        # metal drawn from the element below
    drawn_H: np.ndarray      # the enthalpy that metal brought
    L: np.ndarray            # the liquid on each top node after the pass


def thomas(sub, diag, sup, rhs):
    """Row-wise tridiagonal solve: sub and sup (n, k-1), diag and rhs (n, k). Returns x (n, k)."""
    n, k = diag.shape
    c = np.zeros((n, max(k - 1, 1)))
    d = np.zeros((n, k))
    if k > 1:
        c[:, 0] = sup[:, 0] / diag[:, 0]
    d[:, 0] = rhs[:, 0] / diag[:, 0]
    for i in range(1, k):
        denom = diag[:, i] - sub[:, i - 1] * c[:, i - 1]
        if i < k - 1:
            c[:, i] = sup[:, i] / denom
        d[:, i] = (rhs[:, i] - sub[:, i - 1] * d[:, i - 1]) / denom
    x = np.zeros((n, k))
    x[:, -1] = d[:, -1]
    for i in range(k - 2, -1, -1):
        x[:, i] = d[:, i] - c[:, i] * x[:, i + 1]
    return x


def cumulative(m, Q, x):
    """Per row: the integral from the top down to mass coordinate x (n, j) of a quantity Q (n, nc) carried uniformly
    within each cell of mass m (n, nc). Exact for piecewise-uniform cells; one sorted search for all rows."""
    n, nc = m.shape
    cm = np.concatenate([np.zeros((n, 1)), np.cumsum(m, axis=1)], axis=1)
    cQ = np.concatenate([np.zeros((n, 1)), np.cumsum(Q, axis=1)], axis=1)
    M = np.where(cm[:, -1] > 0.0, cm[:, -1], 1.0)[:, None]
    off = 2.0 * np.arange(n)[:, None]
    flat = (cm / M + off).ravel()
    i = np.searchsorted(flat, (x / M + off).ravel(), side="right").reshape(x.shape) - 1
    i = np.clip(i - (nc + 1) * np.arange(n)[:, None], 0, nc - 1)
    r = np.arange(n)[:, None]
    mi = m[r, i]
    frac = np.clip(np.divide(x - cm[r, i], mi, out=np.zeros_like(x, dtype=float), where=mi > 0.0), 0.0, 1.0)
    return cQ[r, i] + frac * Q[r, i]


class SkinField:
    """Every skin of the body, as arrays (module docstring)."""

    def __init__(self, material, settings=None):
        self.material = material
        self.settings = settings or SkinSettings()
        self.nc = self.settings.n_cells
        self.face_id = np.zeros(0, dtype=np.int64)
        self.area = np.zeros(0)
        self.m = np.zeros((0, self.nc))
        self.H = np.zeros((0, self.nc))
        self.T = np.zeros((0, self.nc))

    # -- books ------------------------------------------------------------------------------------------------------
    @property
    def n(self):
        return int(self.face_id.size)

    def mass(self):
        return float(self.m.sum())

    def energy(self):
        """Metal enthalpy held by the skins [J]. The liquid on their top nodes is MeltingBody.film_energy's."""
        return float(self.H.sum())

    def thickness(self):
        return self.m.sum(axis=1) / (self.material.rho * self.area) if self.n else np.zeros(0)

    def state(self):
        return tuple(a.copy() for a in (self.face_id, self.area, self.m, self.H, self.T))

    def set_state(self, state):
        self.face_id, self.area, self.m, self.H, self.T = (a.copy() for a in state)

    # -- membership -------------------------------------------------------------------------------------------------
    def add(self, face_id, area, T_profile, thickness):
        """Append skins with cell temperatures T_profile (k, nc) and thickness [m] (k,). Returns (mass, H) per new skin."""
        area = np.asarray(area, dtype=float)
        dz = np.asarray(thickness, dtype=float) / self.nc
        m = self.material.rho * (area * dz)[:, None] * np.ones((1, self.nc))
        T = np.asarray(T_profile, dtype=float).reshape(len(area), self.nc)
        H = m * self.material.enthalpy(T)
        self.face_id = np.concatenate([self.face_id, np.asarray(face_id, dtype=np.int64)])
        self.area = np.concatenate([self.area, area])
        self.m, self.H, self.T = np.vstack([self.m, m]), np.vstack([self.H, H]), np.vstack([self.T, T])
        return m.sum(axis=1), H.sum(axis=1)

    def keep(self, rows):
        rows = np.asarray(rows, dtype=np.int64)
        self.face_id, self.area = self.face_id[rows], self.area[rows]
        self.m, self.H, self.T = self.m[rows], self.H[rows], self.T[rows]

    def rows_of(self, face_ids):
        """The row of each face id's skin, -1 where it has none."""
        face_ids = np.asarray(face_ids, dtype=np.int64)
        if not self.n:
            return np.full(face_ids.shape, -1, dtype=np.int64)
        order = np.argsort(self.face_id)
        s = self.face_id[order]
        i = np.clip(np.searchsorted(s, face_ids), 0, len(s) - 1)
        return np.where(s[i] == face_ids, order[i], -1)

    # -- the top node -----------------------------------------------------------------------------------------------
    def top_energy(self, rows, L):
        """Energy of the top node of `rows` with liquid L [kg] on it [J]."""
        return self.H[rows, 0] + L * self.material.enthalpy_liquid(self.T[rows, 0])

    def equilibrate_top(self, rows, L, E_top):
        """Put cell 0 of `rows` and the liquid L [kg] on it at one temperature holding E_top [J] together. Returns T_0."""
        mat = self.material
        m0 = self.m[rows, 0]
        tot = m0 + L
        w = np.divide(L, tot, out=np.zeros_like(tot), where=tot > 0.0)
        h = np.divide(E_top, tot, out=np.zeros_like(tot), where=tot > 0.0)
        T0 = np.where(tot > 0.0, mat.temperature_from_enthalpy_mixed(h, w), self.T[rows, 0])
        self.T[rows, 0] = T0
        self.H[rows, 0] = m0 * mat.enthalpy(T0)
        return T0

    def book_top(self, rows, dE, L):
        """Add dE [J] to the top node of `rows` carrying liquid L [kg] and re-solve its temperature."""
        self.equilibrate_top(rows, L, self.top_energy(rows, L) + dE)

    def refreeze_top(self, rows, taken, L_after):
        """Liquid `taken` [kg] freezes onto cell 0 of `rows`; the node's energy is unchanged (spec 7.1 (6)). L_after is
        the liquid left on the node; the skin is remapped onto equal cells afterwards."""
        E = self.top_energy(rows, L_after + taken)
        self.m[rows, 0] += taken
        self.equilibrate_top(rows, L_after, E)
        self.remap(rows, L_after)

    # -- rezoning ---------------------------------------------------------------------------------------------------
    def _rezone(self, rows, start, D, hD, S=None, S_drawn=1.0):
        """Rows' columns below mass coordinate `start` [kg from the top], with D [kg] of metal added at the bottom at
        specific enthalpy hD [J/kg], redistributed onto nc equal cells: within an old cell the specific enthalpy is
        uniform, so the cumulative enthalpy is linear in the cumulative mass and the new cells' enthalpies are exact
        differences of it. S (an intensive per-cell field, optional) is remapped mass-weighted, the drawn metal
        carrying S_drawn. Returns (m, H, S) of the new cells."""
        m, H = self.m[rows], self.H[rows]
        nc = self.nc
        start, D, hD = (np.broadcast_to(np.asarray(v, dtype=float), (len(rows),)) for v in (start, D, hD))
        rest = m.sum(axis=1) - start
        M_new = rest + D
        y = M_new[:, None] * np.linspace(0.0, 1.0, nc + 1)[None, :]
        x = start[:, None] + np.minimum(y, rest[:, None])
        below = np.maximum(y - rest[:, None], 0.0)
        E_top = cumulative(m, H, start[:, None])
        cH = cumulative(m, H, x) - E_top + below * hD[:, None]
        cH[:, 0] = 0.0
        cH[:, -1] = H.sum(axis=1) - E_top[:, 0] + D * hD                      # exact totals
        m_new = np.repeat((M_new / nc)[:, None], nc, axis=1)
        S_new = None
        if S is not None:
            mS = m * S
            cS = cumulative(m, mS, x) - cumulative(m, mS, start[:, None]) + below * S_drawn
            S_new = np.divide(np.diff(cS, axis=1), m_new, out=np.ones_like(m_new) * S_drawn, where=m_new > 0.0)
        return m_new, np.diff(cH, axis=1), S_new

    def _set_cells(self, rows, m, H, L, E_liq_top):
        """Install cell masses and enthalpies on `rows`, the temperatures from them, the top node with its liquid."""
        mat = self.material
        self.m[rows], self.H[rows] = m, H
        h = np.divide(H, m, out=mat.enthalpy(self.T[rows]), where=m > 0.0)
        self.T[rows] = mat.temperature_from_enthalpy(h)
        self.equilibrate_top(rows, L, H[:, 0] + E_liq_top)

    def remap(self, rows, L):
        """Conservative remap of `rows` (any cell masses) onto nc equal cells; L [kg] is the liquid on their top nodes,
        whose energy stays with the top node."""
        rows = np.asarray(rows, dtype=np.int64)
        if not rows.size:
            return
        E_liq = L * self.material.enthalpy_liquid(self.T[rows, 0])
        m, H, _ = self._rezone(rows, 0.0, 0.0, 0.0)
        self._set_cells(rows, m, H, L, E_liq)

    def add_bottom(self, rows, mass, H, L):
        """Metal `mass` [kg] carrying `H` [J] joins the bottom of `rows` (distinct rows; an element's remainder at its
        death, spec section 9); L [kg] is the liquid on their top nodes."""
        rows = np.asarray(rows, dtype=np.int64)
        self.m[rows, -1] += mass
        self.H[rows, -1] += H
        self.remap(rows, L)

    # -- depths -----------------------------------------------------------------------------------------------------
    def depth_above(self, T_level):
        """Per skin: the depth [m] from the top to where the contiguous cells above T_level end, interpolated between
        cell centres, and whether every cell is above it (the layer then continues into the 3D model)."""
        if not self.n:
            return np.zeros(0), np.zeros(0, dtype=bool)
        k = np.cumprod(self.T > T_level, axis=1).sum(axis=1)
        dz = self.thickness() / self.nc
        r = np.arange(self.n)
        Ta = self.T[r, np.clip(k - 1, 0, self.nc - 1)]
        Tb = self.T[r, np.clip(k, 0, self.nc - 1)]
        frac = np.clip(np.divide(Ta - T_level, Ta - Tb, out=np.zeros_like(Ta), where=Ta > Tb), 0.0, 1.0)
        depth = np.where(k == 0, 0.0, np.where(k == self.nc, self.nc * dz, (k - 0.5 + frac) * dz))
        return depth, k == self.nc

    # -- hand-over --------------------------------------------------------------------------------------------------
    def hand_over(self, src, tgt_faces, weights, face_area):
        """Skins of vanished patches pass to the patches their deaths exposed, layer by layer (spec section 9).
        src (k,) rows; tgt_faces (k, j) face ids of the receiving patches (-1: none); weights (k, j), summing to 1 over
        the valid targets; face_area(face_ids) -> areas. Every source has equal cells, so cell-wise weighted sums do too.
        Returns (receiving rows, mass handed over). The caller removes `src` and any skin left without mass; the
        receiving skins' top nodes are not re-equilibrated with the liquid on them (the caller books that liquid's
        change of enthalpy, as the liquid's own hand-over does)."""
        src = np.asarray(src, dtype=np.int64)
        ok = (tgt_faces >= 0) & (weights > 0.0)
        if not src.size or not ok.any():
            return np.zeros(0, dtype=np.int64), 0.0
        faces = np.unique(tgt_faces[ok])
        new = faces[self.rows_of(faces) < 0]
        if new.size:
            k = new.size
            self.face_id = np.concatenate([self.face_id, new])
            self.area = np.concatenate([self.area, np.asarray(face_area(new), dtype=float)])
            self.m, self.H = np.vstack([self.m, np.zeros((k, self.nc))]), np.vstack([self.H, np.zeros((k, self.nc))])
            self.T = np.vstack([self.T, np.zeros((k, self.nc))])
        sj, tj = np.nonzero(ok)
        trows = self.rows_of(tgt_faces[sj, tj])
        w = weights[sj, tj][:, None]
        np.add.at(self.m, trows, w * self.m[src[sj]])
        np.add.at(self.H, trows, w * self.H[src[sj]])
        recv = np.unique(trows)
        full = recv[self.m[recv].sum(axis=1) > 0.0]
        if full.size:
            self.T[full] = self.material.temperature_from_enthalpy(self.H[full] / np.where(self.m[full] > 0.0, self.m[full], 1.0))
        return recv, float(self.m[src].sum())

    # -- the sub-step -----------------------------------------------------------------------------------------------
    def conduct(self, dt, q_top, T_amb, emissivity, T_b, L, T_surface=None, sens=None):
        """One implicit (backward-Euler) sub-step for every skin (spec section 5.4): the patch's convective flux q_top
        [W/m2] in at the top, radiation from the top cell, conduction through the cells (harmonic conductances), and the
        base held at T_b [K]. The liquid L [kg] on each top node shares cell 0's temperature and heat capacity. Newton
        on the enthalpy with the enthalpy-consistent update the 3D solver uses, damped if the change grows. T_surface
        (tests only) holds the top face at that temperature instead of the flux. `sens = (S, dEb)` carries dT/dT_b and
        the derivative of the base energy through the sub-step (the trial pass). Returns (E_top, E_rad, E_base, sens)
        [J per skin]; E_base is what the books leave, so they close exactly whatever the iteration error."""
        mat, n, nc = self.material, self.n, self.nc
        if n == 0:
            z = np.zeros(0)
            return z, z, z, sens
        A = self.area
        dz = np.maximum(self.m / (mat.rho * A[:, None]), DZ_FLOOR)
        Mt = self.m.copy()
        Mt[:, 0] += L
        w = np.zeros_like(Mt)
        w[:, 0] = np.divide(L, Mt[:, 0], out=np.zeros(n), where=Mt[:, 0] > 0.0)
        T_old = self.T.copy()
        E_old = Mt * mat.enthalpy_mixed(T_old, w)
        eps_s = np.asarray(emissivity, dtype=float) * SIGMA_SB
        q = np.asarray(q_top, dtype=float)
        T_b = np.asarray(T_b, dtype=float)

        def system(T):
            """Residual R [W] of every cell at T, the conductances, and the diagonal of the Newton matrix."""
            k = mat.k(T)
            G = A[:, None] / (0.5 * dz[:, :-1] / k[:, :-1] + 0.5 * dz[:, 1:] / k[:, 1:])
            Gb = A * k[:, -1] / (0.5 * dz[:, -1])
            flux = G * (T[:, :-1] - T[:, 1:])
            net = np.zeros((n, nc))
            diag = Mt * mat.cp_mixed(T, w) / dt
            if T_surface is None:
                net[:, 0] += q * A - eps_s * A * (T[:, 0] ** 4 - T_amb ** 4)
                diag[:, 0] += 4.0 * eps_s * A * T[:, 0] ** 3
            else:
                Gs = A * k[:, 0] / (0.5 * dz[:, 0])
                net[:, 0] += Gs * (T_surface - T[:, 0])
                diag[:, 0] += Gs
            net[:, :-1] -= flux
            net[:, 1:] += flux
            net[:, -1] -= Gb * (T[:, -1] - T_b)
            diag[:, :-1] += G
            diag[:, 1:] += G
            diag[:, -1] += Gb
            return (Mt * mat.enthalpy_mixed(T, w) - E_old) / dt - net, G, Gb, diag

        # Newton in the temperature with a backtracking line search per skin on the residual's norm. Not the 3D solver's
        # enthalpy-consistent update: in 10-micrometre cells conduction outweighs the heat capacity some five hundred
        # times over a 10 ms sub-step, so a cell must move with its neighbours, and inverting each cell's own h(T) moved
        # them apart by a tenth of the step wherever c_p,eff varies through the Scheil range -- the iteration stalled at
        # a 0.4 K change (measured while planning). Across a sharp latent ramp the line search does what that update did.
        T = T_old.copy()
        R, G, Gb, diag = system(T)
        for _ in range(NEWTON_MAX):
            dT = thomas(-G, diag, -G, -R)
            r0 = np.sqrt((R * R).sum(axis=1))
            lam = np.ones(n)
            for _ in range(LINE_SEARCH_MAX):
                T_try = T + lam[:, None] * dT
                R_try, G_try, Gb_try, diag_try = system(T_try)
                moved = np.abs(T_try - T).max(axis=1)
                worse = (np.sqrt((R_try * R_try).sum(axis=1)) > (1.0 - 1e-4 * lam) * r0) & (moved > NEWTON_TOL)
                if not worse.any():
                    break
                lam = np.where(worse, 0.5 * lam, lam)
            T, R, G, Gb, diag = T_try, R_try, G_try, Gb_try, diag_try
            if moved.max() <= NEWTON_TOL:
                break
        else:
            raise RuntimeError("a skin sub-step did not converge in {} iterations (last change {:.2e} K)".format(NEWTON_MAX, moved.max()))
        E_new = Mt * mat.enthalpy_mixed(T, w)
        if T_surface is None:
            E_top = q * A * dt
            E_rad = eps_s * A * (T[:, 0] ** 4 - T_amb ** 4) * dt
        else:
            E_top = A * mat.k(T[:, 0]) / (0.5 * dz[:, 0]) * (T_surface - T[:, 0]) * dt
            E_rad = np.zeros(n)
        E_base = E_top - E_rad - (E_new.sum(axis=1) - E_old.sum(axis=1))
        if sens is not None:
            S, dEb = sens
            rhs = Mt * mat.cp_mixed(T_old, w) / dt * S
            rhs[:, -1] += Gb
            S = thomas(-G, diag, -G, rhs)
            sens = (S, dEb + dt * Gb * (S[:, -1] - 1.0))
        self.T = T
        self.H = self.m * mat.enthalpy(T)
        return E_top, E_rad, E_base, sens

    # -- feed and draw ----------------------------------------------------------------------------------------------
    def feed_and_draw(self, L, room, avail, T_b, S=None):
        """The contiguous fully liquid cells at the top (T >= T_feed) leave the skin for the liquid on its patch -- the
        part within `room` [kg] (rho A delta_m, measured from the top, spec 5.4) to the film, the rest to the deep
        account -- and the skin draws metal from the element below, at most `avail` [kg], at its base temperature T_b,
        back to its target thickness. The remaining column and the drawn metal are rezoned onto equal cells, exactly in
        mass and enthalpy. S (the trial pass's dT/dT_b, optional) is rezoned with them, drawn metal carrying 1.
        Returns (to_film, to_deep, drawn, drawn_H, L_new, S_new) per skin."""
        mat, n = self.material, self.n
        if n == 0:
            z = np.zeros(0)
            return z, z, z, z, z, S
        contig = np.cumprod(self.T >= mat.T_feed, axis=1).astype(bool)
        R = (self.m * contig).sum(axis=1)
        E_rem = (self.H * contig).sum(axis=1)
        to_film = np.minimum(R, room)
        to_deep = R - to_film
        rest = self.m.sum(axis=1) - R
        D = np.minimum(np.clip(mat.rho * self.area * self.settings.thickness - rest, 0.0, None), np.maximum(avail, 0.0))
        hD = mat.enthalpy(np.asarray(T_b, dtype=float) * np.ones(n))
        E_liq = L * mat.enthalpy_liquid(self.T[:, 0])
        rows = np.arange(n)
        m_new, H_new, S_new = self._rezone(rows, R, D, hD, S, 1.0)
        L_new = L + R
        self._set_cells(rows, m_new, H_new, L_new, E_liq + E_rem)
        return to_film, to_deep, D, D * hD, L_new, S_new

    # -- the macro step ---------------------------------------------------------------------------------------------
    def advance(self, dt, q_top, T_amb, emissivity, T_b0, L, room, avail, T_b1=None):
        """The macro step's sub-steps (spec section 6). T_b1=None: the trial pass, on a copy, base held at T_b0, returning
        the base law (a, b) -- q_b(T) = a + b T [W/m2] of heat into the 3D model -- from the pass's base energy and its
        derivative with respect to the base temperature. Otherwise the real pass with the base held at T_b1, the 3D
        model's end-of-step facet temperature (what its backward-Euler step assumed for the boundary): feeds, draws and
        returns a SkinBook."""
        n = self.n
        N = max(1, int(math.ceil(dt / self.settings.substep - 1e-9)))
        h = dt / N
        trial = T_b1 is None
        saved = self.state() if trial else None
        T_b = np.asarray(T_b0 if trial else T_b1, dtype=float) * np.ones(n)
        L = np.array(L, dtype=float)
        avail = np.array(avail, dtype=float)
        book = SkinBook(*(np.zeros(n) for _ in range(7)), L)
        sens = (np.zeros((n, self.nc)), np.zeros(n)) if trial else None
        for _ in range(N):
            E_top, E_rad, E_base, sens = self.conduct(h, q_top, T_amb, emissivity, T_b, L, sens=sens)
            to_film, to_deep, drawn, drawn_H, L, S = self.feed_and_draw(L, room, avail, T_b, None if sens is None else sens[0])
            if sens is not None:
                sens = (S, sens[1])
            avail = avail - drawn
            book.E_top += E_top
            book.E_rad += E_rad
            book.E_base += E_base
            book.to_film += to_film
            book.to_deep += to_deep
            book.drawn += drawn
            book.drawn_H += drawn_H
        book.L = L
        if trial:
            self.set_state(saved)
            # the slope is clipped at zero: near the liquidus a warmer base can feed more and pass more heat down, and a
            # law rising with the base temperature would cost the 3D matrix its diagonal dominance; what the clip
            # misses is the interface mismatch's to book
            b = np.minimum(sens[1] / (self.area * dt), 0.0)
            return book.E_base / (self.area * dt) - b * T_b, b
        return book

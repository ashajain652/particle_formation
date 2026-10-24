"""Apply sub-plan 18 Task 7 (skins in MeltingBody) to the scratch copy's body.py."""
import sys

P = sys.argv[1]
s = open(P).read()


def sub(old, new, count=1):
    global s
    assert s.count(old) == count, (old[:80], s.count(old))
    s = s.replace(old, new)


# (a) imports, names, settings
sub("import numpy as np\n\n\n", "import numpy as np\n\nfrom .thermal import SIGMA_SB\n\n\n")
sub('REMOVAL_NAMES = ("girin", "instant")\n',
    'REMOVAL_NAMES = ("girin", "instant")\n'
    'SURFACE_MODEL_NAMES = ("elements", "skin")    # skin: a one-dimensional skin under every melting patch (sub-plan 18)\n')
sub("""                                      # spray (amendment of 2026-10-06)

    def __post_init__(self):
""", """                                      # spray (amendment of 2026-10-06)
    surface_model: str = "elements"   # skin: a one-dimensional skin of fine cells under every melting patch resolves the
                                      # melt layer and feeds the film continuously (sub-plan 18); elements: as before
    skin: object = None               # skin.SkinSettings for surface_model "skin"; None: its defaults

    def __post_init__(self):
        if self.surface_model not in SURFACE_MODEL_NAMES:
            raise ValueError("surface_model must be one of {}, got {!r}".format(SURFACE_MODEL_NAMES, self.surface_model))
""")

# (b) construction and the geometry refresh
sub("""        self._refresh_geometry()
        self.newtonian_drag_sphere = self.newtonian_drag()[0]      # the meshed sphere's own value: the normalisation
""", """        self._refresh_geometry()
        self.newtonian_drag_sphere = self.newtonian_drag()[0]      # the meshed sphere's own value: the normalisation
        self.skins = None                                          # sub-plan 18: the skins, with --surface-model skin
        if self.settings.surface_model == "skin" and self.settings.removal != "instant":
            from . import skin as skin_mod
            self.skins = skin_mod.SkinField(material, self.settings.skin or skin_mod.SkinSettings())
        self._draw_m, self._draw_H = np.zeros(mesh.n_elements), np.zeros(mesh.n_elements)   # metal the skins drew this step
        self.skins_created = self.skin_short_events = self.interface_repeats = 0
        self.skin_drawn_mass = self.skin_handover_mass = self.last_mismatch = self.skin_wall = 0.0
        self.last_skin_liquid = self.last_skin_nonrigid = None
        self._step_flow = None
        self._refresh_skin_index()
""")
sub("""        self.patch_of_face[self.surface.face_ids] = np.arange(self.surface.n_patches)
        self._fit_nose()
""", """        self.patch_of_face[self.surface.face_ids] = np.arange(self.surface.n_patches)
        self._fit_nose()
        self._refresh_skin_index()
""")

# (c) books: mass, energy, demise
sub("""    def mass(self, t):
        return float((self.phi * self.element_mass).sum() + self.m_f.sum() + self.m_d.sum())
""", """    def mass(self, t):
        m = float((self.phi * self.element_mass).sum() + self.m_f.sum() + self.m_d.sum())
        return m + self.skins.mass() if getattr(self, "skins", None) is not None else m
""")
sub("""    def energy(self):
        return self.solver.energy() + self.film_energy()
""", """    def energy(self):
        E = self.solver.energy() + self.film_energy()
        return E + self.skins.energy() if getattr(self, "skins", None) is not None else E
""")
sub("""        return self.consumed or float((self.phi * self.element_mass).sum()) < self.settings.demise_fraction * self.mass0
""", """        metal = float((self.phi * self.element_mass).sum()) + (self.skins.mass() if self.skins is not None else 0.0)
        return self.consumed or metal < self.settings.demise_fraction * self.mass0
""")

# (d) the step: dispatch, the element path unchanged, the skin path
sub("""    # -- the step --------------------------------------------------------------------------------------------------
    def advance(self, t, dt, loads, state=None):
""", """    # -- the step --------------------------------------------------------------------------------------------------
    def advance(self, t, dt, loads, state=None):
        if self.skins is None:
            return self._advance_elements(t, dt, loads, state)
        return self._advance_with_skins(t, dt, loads, state)

    def _advance_with_skins(self, t, dt, loads, state):
        \"\"\"The macro step with skins (sub-plan 18; spec sections 4 and 6): create skins; the trial pass gives each skin's
        base law; the 3D step takes that law on the skins' facets (no heating, no radiation there) with the deferred
        loads; the real pass, base held at the new facet temperature, feeds, draws and books; the mismatch goes to the 3D
        nodes, and a step whose mismatch exceeds the tolerance repeats once, re-linearised about that temperature. Then
        the rest of the melt step. Until the first skin exists this is the element model's step exactly.\"\"\"
        import dataclasses
        import time
        from . import skin as skin_mod
        from . import spray as spray_mod
        started = time.perf_counter()
        mat, liq, sk = self.material, self.liquid, self.skins
        self._create_skins(self.solver.temperature())
        if not sk.n:
            self._step_flow = None
            return self._advance_elements(t, dt, loads, state)
        on = self.skin_of_patch >= 0
        p = self.patch_of_skin
        flow = delta_m = None
        if state is not None:                     # evaluated before the passes: the feed's reach needs delta_m
            flow = self.flow.evaluate(state, self.theta, self.nose_radius(), liq.rho, self.surface_temperature())
            delta_m, _ = spray_mod.melt_layer(flow, liq)
        self._step_flow = (flow, delta_m)
        dm_p = delta_m[p] if delta_m is not None else np.full(sk.n, np.nan)
        room = np.where(np.isfinite(dm_p), mat.rho * sk.area * np.nan_to_num(dm_p), np.inf)
        avail = self._skin_available(p)
        L = (self.m_f + self.m_d)[p]
        q_top = loads.q_conv[p]
        eps = np.full(sk.n, self.emissivity)
        T_b0 = self.surface.facet_mean(self.solver.temperature())[p]
        q3 = np.where(on, 0.0, loads.q_conv)
        capacity = LOAD_DT_MAX * self.solver.nodal_capacity()
        applied = np.clip(self.pending_load, -capacity, capacity)
        solver_state, skin_state = self.solver.state(), sk.state()
        T_lin = T_b0
        for attempt in (0, 1):
            a, b = sk.advance(dt, q_top, self.T_ambient, eps, T_lin, L, room, avail)
            A_f, B_f = np.zeros(self.surface.n_patches), np.zeros(self.surface.n_patches)
            A_f[p], B_f[p] = a, b
            res = self.solver.step(dt, q3, self.T_ambient, nodal_load=applied / dt, interface=(on, A_f, B_f))
            T_b1 = self.surface.facet_mean(self.solver.temperature())[p]
            book = sk.advance(dt, q_top, self.T_ambient, eps, T_b0, L, room, avail, T_b1=T_b1)
            mismatch = book.E_base - (a + b * T_b1) * sk.area * dt
            scale = np.abs(book.E_top) + np.abs(book.E_base) + 1e-30
            if attempt == 0 and (np.abs(mismatch) > skin_mod.INTERFACE_TOLERANCE * scale).any():
                self.solver.set_state(solver_state)
                sk.set_state(skin_state)
                T_lin = T_b1
                self.interface_repeats += 1
                continue
            break
        self.pending_load = self.pending_load - applied
        full = np.zeros(self.surface.n_patches)
        full[p] = mismatch
        self._defer_to_facets(full)
        self.last_mismatch = float(np.abs(mismatch).sum())
        E_top, E_rad = float(book.E_top.sum()), float(book.E_rad.sum())
        self.integrated_heat += res.Q_conv * dt + E_top
        self.absorbed_heat += (res.Q_conv - res.Q_rad) * dt + E_top - E_rad
        self.radiated_heat += res.Q_rad * dt + E_rad
        self.iterations.append(res.iterations)
        self.last = dataclasses.replace(res, Q_conv=res.Q_conv + E_top / dt, Q_rad=res.Q_rad + E_rad / dt)
        self._applied_total = getattr(self, "_applied_total", 0.0) + res.Q_extra * dt
        self._dropped_total = getattr(self, "_dropped_total", 0.0) + res.Q_dropped * dt
        self.m_f[p] = self.m_f[p] + book.to_film
        self.m_d[p] = self.m_d[p] + book.to_deep
        e = self.surface.owner[p]
        np.add.at(self._draw_m, e, book.drawn)
        np.add.at(self._draw_H, e, book.drawn_H)
        self.skin_drawn_mass += float(book.drawn.sum())
        self.skin_short_events += int((sk.thickness() < 0.5 * sk.settings.thickness).sum())
        if book.to_film.sum() + book.to_deep.sum() > 0.0 and self.melt_onset is None:
            self.melt_onset = t
        self.melt_step(t, dt, state)
        self.skin_wall += time.perf_counter() - started

    def _advance_elements(self, t, dt, loads, state=None):
""")

# (e) the melt step
sub("""        flow = delta_m = None
        if state is not None and s.removal != "instant" and self.surface.n_patches and self.m_f.size:
            flow = self.flow.evaluate(state, self.theta, self.nose_radius(), liq.rho, self.surface_temperature())
""", """        flow = delta_m = None
        if self.skins is not None and self._step_flow is not None:      # with skins: the step's own evaluation
            flow, delta_m = self._step_flow
            self._step_flow = None
            if flow is not None:
                self.last_face_ids = self.surface.face_ids.copy()
        elif state is not None and s.removal != "instant" and self.surface.n_patches and self.m_f.size:
            flow = self.flow.evaluate(state, self.theta, self.nose_radius(), liq.rho, self.surface_temperature())
""")
sub("""        h_p = self.film_enthalpy(h_liq)                          # as film_energy() weighs the film: liquid, with L_f
""", """        h_p = self.film_enthalpy(h_liq)                          # as film_energy() weighs the film: liquid, with L_f
        if self.skins is not None and self._draw_m.any():        # the metal the skins drew leaves its elements
            self.phi = self.phi - self._draw_m / self.element_mass
            self._defer_to_elements(self._draw_m * h_e - self._draw_H)     # the element keeps what the drawn metal left
            self._draw_m[:], self._draw_H[:] = 0.0, 0.0
""")
sub("""        if reach is not None:
            f = f * reach
""", """        if reach is not None:
            f = f * reach
        if self.skins is not None and self.skins.n:              # a patch's skin feeds its film: its element does not
            skin_owner = np.zeros(self.mesh.n_elements, dtype=bool)
            skin_owner[self.surface.owner[self.skin_of_patch >= 0]] = True
            f = np.where(skin_owner, 0.0, f)
""")
sub("""        solid = (1.0 - mat.feed_fraction(self.film_temperature())) * (self.m_f + self.m_d) if s.removal != "instant" else np.zeros(0)
""", """        solid = (1.0 - mat.feed_fraction(self.film_temperature())) * (self.m_f + self.m_d) if s.removal != "instant" else np.zeros(0)
        skin_solid = None
        if self.skins is not None and self.skins.n and solid.size:      # on a skin it freezes into the skin's top cell
            on = self.skin_of_patch >= 0
            skin_solid, solid = np.where(on, solid, 0.0), np.where(on, 0.0, solid)
""")
sub("""            if s.deep_runoff and flow is not None:
                deep_liquid = self._deep_runoff(dt, flow, delta_m, reach, T, h_node, h_e, h_p)
            released = self._film_and_spray(t, dt, state, h_p, flow, delta_m)
        frozen = 0.0 if s.removal == "instant" else self._freeze_back(want, solid, h_e, h_p)
""", """            if s.deep_runoff and flow is not None:
                if self.skins is not None:                       # a skin's top moves when energy is booked on it
                    h_p = self.film_enthalpy(h_liq)
                deep_liquid = self._deep_runoff(dt, flow, delta_m, reach, T, h_node, h_e, h_p)
            if self.skins is not None:
                h_p = self.film_enthalpy(h_liq)
            released = self._film_and_spray(t, dt, state, h_p, flow, delta_m)
        frozen = 0.0 if s.removal == "instant" else self._freeze_back(want, solid, h_e, h_p)
        if skin_solid is not None:
            frozen += self._freeze_into_skins(skin_solid)
""")
sub("""        full = T[tets].mean(axis=1) >= mat.T_feed if s.molten_cascade and s.removal != "instant" else None
""", """        full = T[tets].mean(axis=1) >= mat.T_feed if s.molten_cascade and s.removal != "instant" and self.skins is None else None
""")
sub("""            else:
                extra, extra_h = np.zeros(self.mesh.n_elements), np.zeros(self.mesh.n_elements)
                extra[dead], extra_h[dead] = rest, rest * h_e[dead]
                delta, carried = self._add_to_film(extra, extra_h)
                self._defer_to_patches(carried - delta * h_p)
            self.phi[dead] = 0.0
            before = float(((self.m_f + self.m_d) * h_p).sum())
""", """            else:
                # with skins, an element under a skin passes its remainder to that skin's base (spec section 9)
                went = self._death_to_skins(dead, rest, h_e) if self.skins is not None and self.skins.n else np.zeros(dead.size, dtype=bool)
                extra, extra_h = np.zeros(self.mesh.n_elements), np.zeros(self.mesh.n_elements)
                extra[dead[~went]], extra_h[dead[~went]] = rest[~went], rest[~went] * h_e[dead[~went]]
                delta, carried = self._add_to_film(extra, extra_h)
                self._defer_to_patches(carried - delta * h_p)
            self.phi[dead] = 0.0
            if self.skins is not None:
                h_p = self.film_enthalpy(h_liq)
            before = float(((self.m_f + self.m_d) * h_p).sum())
""")

# (f) the film and spray: resolved depths, and the top's enthalpy after the runoff's booking
sub("""        nonrigid = self.nonrigid_depth() if s.rigid_substrate else None
        under = molten if nonrigid is None else nonrigid
""", """        nonrigid = self.nonrigid_depth() if s.rigid_substrate else None
        if self.skins is not None and self.skins.n:              # the skins resolve both depths (spec section 5.5)
            molten = self._compose_depth(mat.T_feed, molten)
            if nonrigid is not None:
                nonrigid = self._compose_depth(mat.T_rigid, nonrigid)
            self.last_skin_liquid, self.last_skin_nonrigid = molten, nonrigid
        under = molten if nonrigid is None else nonrigid
""")
sub("""            d = self.m_f - before
            self._spread_to_patches(-float((d * h_p).sum()), np.maximum(d, 0.0))
""", """            d = self.m_f - before
            self._spread_to_patches(-float((d * h_p).sum()), np.maximum(d, 0.0))
            if self.skins is not None and self.skins.n:          # the booking moved skin tops: the droplets leave at theirs
                h_p = self.film_enthalpy()
""")

# (g) deaths
sub("""            self.removed_enthalpy += float(((self.m_f + self.m_d) * self.film_enthalpy()).sum())
            self.m_f, self.m_d = np.zeros(0), np.zeros(0)
            return
""", """            self.removed_enthalpy += float(((self.m_f + self.m_d) * self.film_enthalpy()).sum())
            if self.skins is not None and self.skins.n:          # ... and so does the skins' metal
                self.removed_mass += self.skins.mass()
                self.removed_enthalpy += self.skins.energy()
                self.skins.keep(np.zeros(0, dtype=np.int64))
            self.m_f, self.m_d = np.zeros(0), np.zeros(0)
            self.skin_of_patch, self.patch_of_skin = np.zeros(0, dtype=np.int64), np.zeros(0, dtype=np.int64)
            return
""")
sub("""        self.m_d = self._hand_over(old_surface, old_m_d, keep, kept)   # the deep liquid beneath it goes the same way
""", """        self.m_d = self._hand_over(old_surface, old_m_d, keep, kept)   # the deep liquid beneath it goes the same way
        if self.skins is not None and self.skins.n:
            self._hand_over_skins(old_surface)
""")

# (h) readers and bookings that see the skin's top
sub("""    def film_temperature(self):
""", """    def surface_temperature(self):
        \"\"\"Wall temperature per patch: the skin's top cell where a patch has a skin (sub-plan 18), else the facet's.\"\"\"
        Tf = self.surface.facet_mean(self.solver.temperature())
        if getattr(self, "skins", None) is not None and self.skins.n:
            on = self.skin_of_patch >= 0
            Tf[on] = self.skins.T[self.skin_of_patch[on], 0]
        return Tf

    def radiated_power(self):
        P = self.solver.radiated_power(self.T_ambient)
        if getattr(self, "skins", None) is not None and self.skins.n:   # skins radiate from their tops, at patch area
            T0 = self.skins.T[:, 0]
            P += float((self.emissivity * SIGMA_SB * (T0 ** 4 - self.T_ambient ** 4) * self.skins.area).sum())
        return P

    def film_temperature(self):
""")
sub("""        if h_node is None:
            h_node = self.material.enthalpy_liquid(self.solver.temperature())
        return h_node[self.surface.faces].mean(axis=1) if self.m_f.size else np.zeros(0)
""", """        if h_node is None:
            h_node = self.material.enthalpy_liquid(self.solver.temperature())
        if not self.m_f.size:
            return np.zeros(0)
        h = h_node[self.surface.faces].mean(axis=1)
        if getattr(self, "skins", None) is not None and self.skins.n:   # liquid on a skin shares its top cell's temperature
            on = self.skin_of_patch >= 0
            h[on] = self.material.enthalpy_liquid(self.skins.T[self.skin_of_patch[on], 0])
        return h
""")
sub("""    def _defer_to_patches(self, energy):
        \"\"\"Book a per-patch energy [J] on the patches' nodes, to be applied over the next step.\"\"\"
        e = np.asarray(energy, dtype=float)
        if e.any():
            np.add.at(self.pending_load, self.surface.faces.ravel(), np.repeat(e / 3.0, 3))
""", """    def _defer_to_patches(self, energy):
        \"\"\"Book a per-patch energy [J]: on the nodes of a patch without a skin, to be applied over the next step; on a
        skin's top node at once (sub-plan 18), shared with the liquid that node carries.\"\"\"
        e = np.asarray(energy, dtype=float)
        if not e.any():
            return
        if getattr(self, "skins", None) is not None and self.skins.n:
            on = (self.skin_of_patch >= 0) & (e != 0.0)
            if on.any():
                p = np.flatnonzero(on)
                self.skins.book_top(self.skin_of_patch[p], e[p], (self.m_f + self.m_d)[p])
                e = np.where(self.skin_of_patch >= 0, 0.0, e)
        if e.any():
            np.add.at(self.pending_load, self.surface.faces.ravel(), np.repeat(e / 3.0, 3))
""")
sub("""        if not self.m_f.size:
            return np.zeros(self.mesh.n_nodes)
        return np.bincount(self.surface.faces.ravel(), np.repeat((self.m_f + self.m_d) / 3.0, 3), self.mesh.n_nodes)
""", """        if not self.m_f.size:
            return np.zeros(self.mesh.n_nodes)
        lm = self.m_f + self.m_d
        if getattr(self, "skins", None) is not None:              # liquid on a skin rides the skin's top, not the nodes
            lm = np.where(self.skin_of_patch >= 0, 0.0, lm)
        return np.bincount(self.surface.faces.ravel(), np.repeat(lm / 3.0, 3), self.mesh.n_nodes)
""")

# (i) the skins' own methods, after the liquid's hand-over
sub("""    # -- reporting -----------------------------------------------------------------------------------------------
    def on_current_surface(self, values, fill):
""", '''    # -- skins (sub-plan 18) -----------------------------------------------------------------------------------------
    def _refresh_skin_index(self):
        """skin_of_patch: the skin row under each patch of the current surface (-1: none); patch_of_skin the reverse
        (-1 for a skin whose patch has vanished, until the death loop hands it over)."""
        self.skin_of_patch = np.full(self.surface.n_patches, -1, dtype=np.int64)
        sk = getattr(self, "skins", None)
        if sk is None or not sk.n:
            self.patch_of_skin = np.zeros(0, dtype=np.int64)
            return
        self.patch_of_skin = self.patch_of_face[sk.face_id]
        ok = self.patch_of_skin >= 0
        self.skin_of_patch[self.patch_of_skin[ok]] = np.flatnonzero(ok)

    def _skin_available(self, patches):
        """Metal [kg] a skin on each of `patches` may take from its element: its share, by area among the element's
        patches, of what the element holds above PHI_DEATH."""
        e = self.surface.owner[patches]
        share = self.surface.areas[patches] / self.owner_area[e]
        return np.maximum(self.phi[e] - PHI_DEATH, 0.0) * self.element_mass[e] * share

    def _create_skins(self, T):
        """Patches whose surface has reached the solidus, or that hold liquid, get a skin (spec 5.2): a slab of the
        target thickness -- less where the element cannot give that above PHI_DEATH -- taken from the patch's element,
        its cells filled from the P1 field along the inward normal. The element's mass leaves at its own mean nodal
        enthalpy and the difference to what the cells hold is booked on its nodes; liquid already on the patch moves
        onto the skin's top with its energy. Returns how many skins were made."""
        sk, s, mat = self.skins, self.surface, self.material
        if sk is None or not s.n_patches:
            return 0
        Tf = s.facet_mean(T)
        idx = np.flatnonzero((self.skin_of_patch < 0) & ((Tf >= mat.T_solidus) | (self.m_f + self.m_d > 0.0)))
        if not idx.size:
            return 0
        thick = np.minimum(sk.settings.thickness, self._skin_available(idx) / (mat.rho * s.areas[idx]))
        idx, thick = idx[thick > 0.0], thick[thick > 0.0]
        if not idx.size:
            return 0
        e = s.owner[idx]
        z = (np.arange(sk.nc) + 0.5)[None, :] * (thick / sk.nc)[:, None]
        x = s.centroids[idx][:, None, :] - s.normals[idx][:, None, :] * z[:, :, None]
        P = self.mesh.points[self.mesh.tets[e]]
        inv = np.linalg.inv(np.stack([P[:, 1] - P[:, 0], P[:, 2] - P[:, 0], P[:, 3] - P[:, 0]], axis=2))
        lam = np.einsum("kij,kcj->kci", inv, x - P[:, None, 0])
        lam = np.concatenate([1.0 - lam.sum(axis=2, keepdims=True), lam], axis=2)
        lam = np.clip(lam, 0.0, None)
        lam = lam / lam.sum(axis=2, keepdims=True)            # a point past the element takes its nearest face's value
        Tp = np.einsum("kci,ki->kc", lam, T[self.mesh.tets[e]])
        L = (self.m_f + self.m_d)[idx]
        h_nodes = self.film_enthalpy()[idx]                   # the liquid's enthalpy on the 3D nodes, before the skins
        mass, H = sk.add(s.face_ids[idx], s.areas[idx], Tp, thick)
        h_e = mat.enthalpy(T)[self.mesh.tets[e]].mean(axis=1)
        np.subtract.at(self.phi, e, mass / self.element_mass[e])
        self._defer_to_elements(np.bincount(e, mass * h_e - H, self.mesh.n_elements))
        self._refresh_skin_index()
        if (L > 0.0).any():
            rows = self.skin_of_patch[idx]
            sk.book_top(rows, L * (h_nodes - mat.enthalpy_liquid(sk.T[rows, 0])), L)
        self.solver.set_fractions(self.phi)
        self.solver.set_film_mass(self._film_nodal())
        self.skins_created += int(idx.size)
        return int(idx.size)

    def _compose_depth(self, level, depth3d):
        """Per patch: the depth of the layer above `level` read from the skin where there is one, continued into the 3D
        model's own depth (measured from the skin's base) when the whole skin is above it (spec 5.5)."""
        out = np.array(depth3d, dtype=float, copy=True)
        if self.skins is None or not self.skins.n:
            return out
        d, through = self.skins.depth_above(level)
        on = self.skin_of_patch >= 0
        r = self.skin_of_patch[on]
        out[on] = d[r] + np.where(through[r], out[on], 0.0)
        return out

    def _defer_to_facets(self, energy):
        """Book a per-patch energy [J] on the 3D nodes of the patches' facets, skins or not (the interface mismatch)."""
        e = np.asarray(energy, dtype=float)
        if e.any():
            np.add.at(self.pending_load, self.surface.faces.ravel(), np.repeat(e / 3.0, 3))

    def _freeze_into_skins(self, solid):
        """Liquid on a skin patch whose top has fallen below the feed ramp refreezes into the skin's top cell, the node's
        energy unchanged (spec 7.1 (6)); the deep liquid, against the solid, goes first. Returns the mass [kg]."""
        p = np.flatnonzero(solid > 0.0)
        if not p.size:
            return 0.0
        from_deep = np.minimum(self.m_d[p], solid[p])
        from_film = np.minimum(self.m_f[p], solid[p] - from_deep)
        taken = from_deep + from_film
        self.m_d[p] = self.m_d[p] - from_deep
        self.m_f[p] = self.m_f[p] - from_film
        self.skins.refreeze_top(self.skin_of_patch[p], taken, (self.m_f + self.m_d)[p])
        return float(taken.sum())

    def _death_to_skins(self, dead, rest, h_e):
        """Dying elements under skins pass their remainder `rest` [kg] into those skins' bottoms, shared by area among
        the element's skin patches, at the element's own mean nodal enthalpy (spec 9). Returns which of `dead` did."""
        on = self.skin_of_patch >= 0
        p = np.flatnonzero(on & np.isin(self.surface.owner, dead))
        went = np.isin(dead, self.surface.owner[p])
        if p.size:
            e = self.surface.owner[p]
            skin_area = np.bincount(e, self.surface.areas[p], self.mesh.n_elements)
            full = np.zeros(self.mesh.n_elements)
            full[dead] = rest
            mass = full[e] * self.surface.areas[p] / skin_area[e]
            self.skins.add_bottom(self.skin_of_patch[p], mass, mass * h_e[e], (self.m_f + self.m_d)[p])
        return went

    def _hand_over_targets(self, old_surface, idx):
        """For vanished patches `idx` of `old_surface`: (targets (k, 4), patch indices of the current surface or -1;
        weights (k, 4)) -- the faces their own element exposed, weighted by area projected on the vanished patch's
        normal (raw area where none faces that way), or the nearest patch where the element exposed nothing. The rule
        of `_hand_over`, kept apart so that the liquid's hand-over keeps its exact arithmetic."""
        targets = self.patch_of_face[self.mesh.element_faces(old_surface.owner[idx])]
        ok = targets >= 0
        safe = np.where(ok, targets, 0)
        cosine = np.clip(np.einsum("nij,nj->ni", self.surface.normals[safe], old_surface.normals[idx]), 0.0, None)
        area = np.where(ok, self.surface.areas[safe] * cosine, 0.0)
        total = area.sum(axis=1)
        flat = ok & (total[:, None] <= 0.0)
        if flat.any():
            area = np.where(flat, self.surface.areas[safe], area)
            total = area.sum(axis=1)
        has = total > 0.0
        weights = np.where(has[:, None], area / np.where(has, total, 1.0)[:, None], 0.0)
        targets = np.where(ok & has[:, None], targets, -1)
        if (~has).any():
            orphan = np.flatnonzero(~has)
            targets[orphan, 0] = self._patch_tree.query(old_surface.centroids[idx[orphan]])[1]
            weights[orphan, 0] = 1.0
        return targets, weights

    def _hand_over_skins(self, old_surface):
        """Skins whose patch vanished pass to the faces the deaths exposed, layer by layer (spec 9); empty skins go."""
        sk = self.skins
        lost = np.flatnonzero(self.patch_of_face[sk.face_id] < 0)
        if lost.size:
            old_patch = np.full(len(self.patch_of_face), -1, dtype=np.int64)
            old_patch[old_surface.face_ids] = np.arange(old_surface.n_patches)
            targets, weights = self._hand_over_targets(old_surface, old_patch[sk.face_id[lost]])
            tgt_faces = np.where(targets >= 0, self.surface.face_ids[np.maximum(targets, 0)], -1)
            _, handed = sk.hand_over(lost, tgt_faces, weights, lambda f: self.surface.areas[self.patch_of_face[f]])
            self.skin_handover_mass += handed
        keep = np.ones(sk.n, dtype=bool)
        keep[lost] = False
        keep &= sk.m.sum(axis=1) > 0.0
        sk.keep(np.flatnonzero(keep))
        self._refresh_skin_index()

    # -- reporting -----------------------------------------------------------------------------------------------
    def on_current_surface(self, values, fill):
''')

open(P, "w").write(s)
print("body patched")

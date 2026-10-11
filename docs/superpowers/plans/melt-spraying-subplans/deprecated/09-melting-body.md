# Sub-plan: Task 9 — The melting body

> Extracted verbatim from `2026-09-20-melt-spraying.md` (current version, lines 4349–5703). Read `00-shared-context.md` first — this is the largest and highest-risk task in the plan and needs the most context, not the least. **Amended 2026-10-02: the deep runoff and the per-patch conjugate depth** (section below). **Amended 2026-10-03: the molten cascade** (the section after it). **Amended 2026-10-05 (runoff flux): what the body sees of the film's new flux** (the section after that). **Amended 2026-10-06: the thin branch needs a rigid substrate — the non-rigid depth and the regime layer** (the section after the runoff flux's). **Amended 2026-10-07: freeze-back never leaves a negative film** (the section after the rigid substrate's).


> **Amended 2026-09-27** by `docs/superpowers/specs/2026-09-27-surface-recession-remeshing-design.md`
> (measured facts 38–45 in `00-shared-context.md`). The block below takes precedence over the extracted body.

**Changes to Task 9.**

1. **`PHI_DEATH = 0.05` becomes `PHI_DEATH = 0.50`** (fact 44). A patch owner dies at half consumed, not at 95 %
   consumed. Two measured reasons: the frozen-geometry conduction error is lowest there (+1.97 K typical against
   +2.97 K at 0.05, +9.87 K against +14.83 K at peak flux, a shallow optimum at a half), and fact 42's quality
   envelope wants it if Design A is ever adopted. The `PHI_MIN = 1e-3` floor for interior elements is unchanged, as
   is fact 4's rule that the remainder joins the film and that deaths cascade within the step.
2. **Check the film can absorb it.** The film now receives half an element at a time rather than a twentieth.
   `film_blob_fraction` and `film_frozen_fraction` must stay within the ranges fact 25 records, and the stranded-film
   mass below the 0.47 % of initial mass measured there. This is an acceptance test, not an assumption.
3. **The shape consumers read the derived surface.** `nose_radius()`, `newtonian_drag()`, `reference_area()`,
   `drag_shape_factor()`, `equivalent_radius()`, `transverse_radius` and `mass_centre` take their geometry from
   `mesh.surface(derived=True)` and its `smoothed_normals()`, not from the raw element faces. `theta` and `t_hat`
   likewise. The **mass** account is untouched: it remains `Σ φ_e ρ V_e` over the **original** element volumes, and
   `φ_e` remains authoritative (spec §5, Design B).
4. **New attribute** `last_recession_residual` — the median absolute per-patch mismatch from Task 16's least-squares
   solve — and a `melt_stats()` entry for it, so how well the geometry tracks the mass is recorded per step rather
   than assumed. Fact 43: 0.03 % for a smooth recession field, 22 % for a patchy one, and which of those the real
   flight produces is unmeasured.
**Also carried here: measured fact 37 (2026-09-27), which the committed master plan is missing.** It was
regenerated into `prototype/plan3/plan.md` but never copied over `2026-09-20-melt-spraying.md`, so the extracted
body below predates it. `MeltingBody` gains the attribute `last_face_ids` and the method `on_current_surface`; add
both to the `Produces:` list, and note that `last_spray`/`last_flow` are now read through the latter:

```python
        self.last_face_ids = None                                          # in __init__, beside last_flow/last_spray
        self.last_face_ids = self.surface.face_ids.copy()                  # at step (iii), before the deaths of step (v)

    def on_current_surface(self, values, fill):
        """A per-patch array from the last flow and spray evaluation (`last_flow`, `last_spray`), carried onto the
        current surface. Element deaths later in the same melt step rebuild the surface, so the arrays no longer line
        up with it by index; patches are matched by face id instead, and a patch that was not on the surface at the
        evaluation (a face the deaths exposed) gets `fill`. Matching by index instead left every frame written after
        a step with a death -- every melting frame -- at the defaults: 0 of 104 frames of the 50 mm physics flight
        carried a wall pressure (measured 2026-09-27)."""
        n = self.surface.n_patches
        if values is None or self.last_face_ids is None or not self.last_face_ids.size:
            return np.full(n, fill, dtype=float)
        index = np.full(len(self.patch_of_face), -1, dtype=np.int64)
        index[self.last_face_ids] = np.arange(self.last_face_ids.size)
        j = index[self.surface.face_ids]
        return np.where(j >= 0, np.asarray(values, dtype=float)[np.maximum(j, 0)], fill)

```

Without it, every frame written after a step with a death — which is every melting step — carries the defaults:
**0 of 104** frames of the 50 mm physics flight had a wall pressure. Output only; the history columns move by at
most 1.0e-8 relative, against 1.6e-8 between two runs of one build.


5. **New `melt_stats()` entry** `conduction_dT_err_K` = max over live patch owners of `q L (1 − φ_e) / k`, the
   surface temperature error Design B accepts by conducting through material that is no longer there. Spec §10.3
   makes this the measured trigger for reconsidering Design A.
---

## Amendment of 2026-10-02 — the liquid below the conjugate depth runs off, and only the skin sprays

> Requested by Asha on 2026-10-02 for the large-fragment (Spheral) model that will read this model's frames. Measured
> in a throwaway copy of the prototype; the code below is the tested code, as a diff against the prototype's
> `body.py` (which is the extracted body below with the fact-37 code of the 2026-09-27 amendment applied). Facts
> 46–53 in `00-shared-context.md` carry the measurements. **Read facts 4, 5, 25, 26, 28 and 29 before reading the
> code**, and sub-plan 02's warning about the four enthalpy functions: the deep runoff is one more transfer of liquid
> between the elements, the patches and the outside world, and the two implementations rejected in facts 5 and 25 are
> as easy to reintroduce here as anywhere.

**The three-zone rule (decision of 2026-10-02).** Melted material belongs in one of three places. (1) The sprayable
skin: liquid above the 908 K liquidus within Girin's conjugate depth δ_m of the surface (215–303 µm on the 100 mm
flight, fact 28) stays in the film and is sprayed by Girin's mechanism. (2) The thin runoff layer: liquid below δ_m,
down to a film limit of about 2–3 mm, is not sprayed but runs off, driven by the pressure gradient along the surface
and the deceleration; it stays in this model. (3) Material thicker than the film limit goes to the large-fragment
model. Two changes to this task follow, both requested: **(A)** the conjugate depth of every patch is kept for the VTK
frames (`last_delta_m`, written by Task 10), and **(B)** the liquid of zone 2 runs off and is never sprayed.

**What the copy contained.** `prototype/proto3/` as it stood on 2026-10-02: Step 3 through the amendments of
2026-09-25, fact 37's `on_current_surface`, and sub-plan 01's dense band and derived-surface plumbing (`mesh.py`; the
melting runs mesh the 100 mm sphere with 177 363 tetrahedra and 18 830 surface patches). It does **not** contain sub-
plan 02's material amendments (every run here uses `AA7075_range`), this task's changes of 2026-09-27 (`PHI_DEATH`
0.50; the shape consumers reading the derived surface) or Tasks 16 and 17, which exist only as plans. The deep runoff
needs none of them. Two interactions are expected and are follow-ups: Task 17 must carry the new account across a
remesh (its amendment of this date), and Task 16's smoothed normals should remove most of the craters the deep liquid
now collects in (fact 47; Task 14's amendment asks for the re-measurement).

### The design: a deep account per patch, fed only by what flows

This is option (i) of the request — a separate, non-sprayable deep account per patch, moved by the pressure- and
deceleration-driven flux over the whole liquid depth — with one refinement: the account is fed only by what actually
flows out of the elements, not by the whole held inventory every step (the literal form, rejected below). One new
per-patch account, `m_d` [kg]: liquid below the conjugate depth that the runoff has moved onto a patch. One new
stage, `_deep_runoff`, between the feed (i) and the film step (ii) of `melt_step`:

1. **Where it acts.** Windward patches under Girin's closure (`np.isfinite(delta_m)`): only there does a conjugate
   depth exist, and so a layer below the shear's reach. Elsewhere — the Couette closure, the merged and free-molecular
   branches, the lee — nothing moves, and any deep liquid that arrives there (by runoff or at an element death) has
   reached the wall and becomes film at once.
2. **What it moves.** The deep liquid of a patch is the fully liquid inventory, min(f_feed φ_e, φ_e − PHI_MIN) ρ V_e,
   of the elements of its contiguous molten chain (the march of `molten_depth`, now `_molten_chain`, which returns the
   elements as well as the depth) that the feed gate held back because the shear does not reach them — shared among
   the patches whose chains pass through them in proportion to patch area — plus its own `m_d`.
3. **How it moves.** With h_D that liquid's thickness by mass and b the film's, by `film.deep_flux`
   G ((b + h_D)³ − b³)/(3 μ_l) (sub-plan 06's amendment; fact 29's pressure-driven part over the whole liquid depth),
   on the film's own linearly implicit upwind transport (`Runoff.transport`, fact 7). h_D is the mass actually there,
   not `molten_depth`'s geometric depth: that counts the owner element, whose liquid the feed has already moved into the
   film, and the whole volume of elements that are only partly liquid or partly consumed. On the 100 mm flight the
   geometric molten layer under the patches whose liquid is deeper than δ_m is a median 1.86 mm, while the liquid
   actually held below δ_m, spread over the same patches, is a few tenths of a millimetre (a median of about 0.3 mm,
   estimated with the mean windward facet area) — a factor of order a hundred in h³, and so in the flux.
4. **Who gives it up, who receives it.** A patch whose deep liquid falls gives it up in proportion from its chain's
   held elements and its `m_d`; an element is debited only by what actually flows out of it, never below `PHI_MIN`, so
   the melt pool keeps its place in the mesh and is still delivered to the film by the owner rule as the surface
   recedes to it. Whatever arrives goes to `m_d`. The arrivals are scaled to the departures (fact 48).
5. **Only the skin sprays.** `m_d` is never offered to the spray. It becomes film only from the top: after the
   transport the film is topped up from `m_d` to the conjugate depth (min(m_d, ρ_l A δ_m − m_f)), i.e. at most one
   skin's worth per macro step — the free surface receding into the deep liquid as the skin above it is stripped —
   and all at once where there is no conjugate depth.
6. **Everything else the film does, the deep account does with it.** It rides the patch's three nodes with the film
   (`_film_nodal`, so the solver carries its heat capacity and the mixed-enthalpy Newton update sees it), holds the
   same liquid enthalpy (`film_energy` sums both), freezes back first when the patch falls below the feed ramp (it
   lies against the solid), is handed to the faces a dying element exposes by the film's rule (`_hand_over`, factored
   out of `_kill`), leaves with the body if it is consumed, and counts in `mass()`, the layer the branch and
   Rayleigh–Taylor tests read, and `liquid_layer_depth()`.

### How every transfer is booked (facts 4, 5, 25 and 26)

- **Liquid leaving an element** is booked exactly as the feed books it. φ_e falls by the mass that left, which removes
  `out h_e` from the solver's energy (h_e the mean of the element's nodal enthalpies, `enthalpy` — the mixture
  function, because that is what the solver holds for the element), while the liquid takes away only `out h_out`, the
  enthalpy of the element's *molten part* (the f_feed-weighted mean of the nodal enthalpies). The difference
  `out (h_e − h_out)` stays with the element and is booked on its own nodes. Debiting only the destination — the first
  implementation fact 5 rejects — would leave the element holding melt it no longer has.
- **Liquid arriving on a patch** holds the patch's *liquid* enthalpy (`film_enthalpy` of `enthalpy_liquid`, never the
  mixture's — the second rejected implementation, fact 25), and the difference between what the departing liquid
  carried and what the arrivals now hold is booked once, at the arrivals, in proportion to the mass each received
  (`_spread_to_patches`) — never as two halves on two sets of nodes (fact 25's third rule).
- **Deep liquid becoming film, and freezing back.** The film and the deep account of one patch hold the same liquid
  enthalpy on the same nodes, so moving mass between them books nothing; the freeze-back draws the deep account first
  and books the transfer to the owner element exactly as before.
- **Nothing escapes the bounds of facts 4 and 26.** No element goes below `PHI_MIN`; every deferred energy goes
  through `pending_load` and is applied subject to `LOAD_DT_MAX`.
- The four enthalpy functions keep their roles (sub-plan 02's warning): `enthalpy`, the mixture, for what an element
  holds and gives up; `enthalpy_liquid` for both liquid accounts on a patch; and `enthalpy_mixed`, `cp_mixed` and
  `temperature_from_enthalpy_mixed` inside both solvers' Newton update, whose film weight now includes the deep account
  because it is liquid on those nodes. None of them changes.

### Alternatives considered, and why they were not taken

- **Option (i) in its literal form: feed the whole held inventory into the deep account every step.** Simplest to code
  (it is the feed into a second account), but it lifts the melt pool out of the mesh: the drained elements die as soon
  as the surface reaches them, so the surface recedes through the pool in one step, the pool is drawn above the new
  surface, and the owner elements' new melt — physically at the *bottom* of the liquid column — is fed to the skin
  *above* a non-sprayable pool.
- **(ii) Element-to-element transfer along the surface.** The receiving elements of a molten pool are full (φ = 1), so
  the inflow has nowhere to go in the mesh; the outflow opens partly empty interior elements, which `PHI_MIN` exists
  to prevent; and there is no element graph along the surface to transport on — the runoff graph is per patch.
- **(iii) Feed the liquid below δ_m into the film account with a non-sprayable flag.** The bookkeeping is the same as
  the chosen design's, but every consumer of the film — the spray branches, the film diagnostics, the blob and frozen
  fractions, the Rayleigh–Taylor criteria — would have to learn the flag. A separate account leaves them untouched.
- **(α) Deep liquid exists only within the in-situ pool, and surfaces into the film wherever it leaves it.** Considered
  and rejected because liquid reaching the edge of the pool would become film, and sprayable, however deep it lay —
  which is exactly what the rule forbids.
- **A two-directional partition by depth (the film also capped at δ_m, its excess becoming deep).** The literal reading
  of "only the skin within δ_m sprays", but it changes the feed and the film on every thick patch, including where no
  liquid lies below the conjugate depth, so it fails the requirement that those results stay unchanged.

Code (the tested change to `reentry_model/body.py`):

```diff
--- a/reentry_model/body.py
+++ b/reentry_model/body.py
@@ -156,6 +156,7 @@
     size_feedback: str = "current"    # current: Kn on the equivalent diameter of the remaining mass and the nose radius fitted
                                       # to the windward cap; initial: D0 and R0 throughout (SESAM's convention, for the devices)
     feed_depth: str = "conjugate"     # conjugate: only liquid the gas shear reaches leaves its element (fact 28); all: any depth
+    deep_runoff: bool = True          # the liquid below the conjugate depth runs off, never sprayed (amendment of 2026-10-02)
 
     def __post_init__(self):
         if self.removal not in REMOVAL_NAMES:
@@ -185,7 +186,11 @@
     (i) every element's liquid inventory f_feed(T_e) phi_e rho V_e becomes film on its patches (owners) or the
     nearest patch (interior elements keep PHI_MIN so no cavity can open), netted against (iv) below; what leaves is
     the element's molten part, at the enthalpy that part carries, so the melt front moves at the energy-limited rate
-    whatever the element size; (ii) surface flow, delta_m, lubrication, runoff transport; (iii) spraying and release;
+    whatever the element size; (ii-a) the deep runoff (`_deep_runoff`, amendment of 2026-10-02): the contiguous liquid
+    below the conjugate depth, which (i) holds in its elements, moves along the surface under the pressure gradient and
+    the deceleration into a per-patch deep account `m_d` that is never sprayed and becomes film only from the top, as
+    far as the film is thinner than the conjugate depth; (ii) surface flow, delta_m, lubrication, runoff transport;
+    (iii) spraying and release;
     (iv) re-solidification: the fraction 1 - f_feed(T_patch) of each patch's film returns to its owner element
     (raising phi_e, capped at 1) once the patch falls back through the ramp -- the mirror of the feed rule, netted
     with it so that no element both melts and freezes in one step; (v) patch owners at phi <= PHI_DEATH die: the
@@ -232,12 +237,16 @@
         self.phi = np.ones(mesh.n_elements)
         self.mass0 = float(self.element_mass.sum())
         self.m_f = np.zeros(self.surface.n_patches)
+        self.m_d = np.zeros(self.surface.n_patches)                         # kg, mobile liquid below the conjugate depth
         self.pending_load = np.zeros(mesh.n_nodes)                          # J, deferred melt energy for the next step
         self.sprayed_mass = self.runoff_mass = self.removed_mass = self.removed_enthalpy = self.frozen_mass = 0.0
+        self.deep_runoff_mass = 0.0                                         # kg the deep runoff took out of the elements
+        self.deep_surfaced_mass = 0.0                                       # kg of deep liquid that became film from the top
         self.n_released = 0.0
         self.source_rows = []
         self.hist_n, self.hist_m = np.zeros(spray_mod.N_BINS), np.zeros(spray_mod.N_BINS)
         self.last_flow = self.last_spray = None
+        self.last_delta_m = None                                           # the step's conjugate depth per patch [m]
         self.last_face_ids = None                                          # face ids of the surface they were evaluated on
         self.last_melt = {"n_dead": 0, "runoff_substeps": 0, "released_mass": 0.0, "n_released": 0.0, "feed_mass": 0.0}
         self.melt_onset = self.spray_onset = None
@@ -325,7 +334,7 @@
         return 2.0 * self.equivalent_radius()
 
     def mass(self, t):
-        return float((self.phi * self.element_mass).sum() + self.m_f.sum())
+        return float((self.phi * self.element_mass).sum() + self.m_f.sum() + self.m_d.sum())
 
     def reference_area(self):
         """Drag reference area: the hull's projected area, paired with the C_D of drag_shape_factor (the raw surface's
@@ -400,6 +409,8 @@
         if state is not None and s.removal != "instant" and self.surface.n_patches and self.m_f.size:
             flow = self.flow.evaluate(state, self.theta, self.nose_radius(), liq.rho, self.surface_temperature())
             delta_m, _ = spray_mod.melt_layer(flow, liq)
+            self.last_face_ids = self.surface.face_ids.copy()      # the surface delta_m belongs to (the deaths come later)
+        self.last_delta_m = delta_m                                # per patch, nan off Girin's closure; for the VTK frame
         h_node = mat.enthalpy(T)
         h_e = h_node[tets].mean(axis=1)                          # as solver.energy() weighs the solid
         h_liq = mat.enthalpy_liquid(T)
@@ -416,12 +427,14 @@
         # net effect, pinning the surface at T_feed and paying Newton iterations for it.
         fn = mat.feed_fraction(T)[tets]
         f = fn.mean(axis=1) * self.mesh.active
-        if s.feed_depth == "conjugate" and s.removal != "instant":
-            f = f * self._shear_reaches(delta_m)
+        reach = self._shear_reaches(delta_m) if s.feed_depth == "conjugate" and s.removal != "instant" else None
+        if reach is not None:
+            f = f * reach
         h_hot = (fn * h_node[tets]).mean(axis=1) * self.mesh.active     # what the molten part carries, per kg of element
         cap = np.where(self.owner_area > 0.0, self.phi, np.maximum(self.phi - PHI_MIN, 0.0))
         gross = np.minimum(f * self.phi, cap) * self.element_mass
-        solid = (1.0 - mat.feed_fraction(self.film_temperature())) * self.m_f if s.removal != "instant" else np.zeros(0)
+        # the liquid on each patch below the ramp re-solidifies: the film and the deep liquid beneath it alike
+        solid = (1.0 - mat.feed_fraction(self.film_temperature())) * (self.m_f + self.m_d) if s.removal != "instant" else np.zeros(0)
         want = np.bincount(self.surface.owner, solid, self.mesh.n_elements) if solid.size else np.zeros(self.mesh.n_elements)
         net = gross - want
         fed, want = np.maximum(net, 0.0), np.maximum(-net, 0.0)
@@ -436,6 +449,7 @@
         feed_mass = float(fed.sum())
         if feed_mass > 0.0 and self.melt_onset is None:
             self.melt_onset = t
+        deep_liquid = 0.0
         if s.removal == "instant":
             self.removed_mass += feed_mass
             self.removed_enthalpy += feed_mass * mat.h_liquid
@@ -445,6 +459,8 @@
             self._defer_to_elements(fed * h_e - fed_h)           # the element keeps only what the melt left behind
             delta, carried = self._add_to_film(fed, fed_h)
             self._defer_to_patches(carried - delta * h_p)        # ... and the melt arrives at its patch's temperature
+            if s.deep_runoff and flow is not None:
+                deep_liquid = self._deep_runoff(dt, flow, delta_m, reach, T, h_node, h_e, h_p)
             released = self._film_and_spray(t, dt, state, h_p, flow, delta_m)
         frozen = 0.0 if s.removal == "instant" else self._freeze_back(want, solid, h_e, h_p)
         # (v) death of consumed patch owners; a death exposes its neighbours, which die in turn if they are consumed
@@ -464,23 +480,26 @@
                 delta, carried = self._add_to_film(extra, extra_h)
                 self._defer_to_patches(carried - delta * h_p)
             self.phi[dead] = 0.0
-            before = float((self.m_f * h_p).sum())
-            self._kill(dead)                                     # the film of a dead patch moves to another patch...
+            before = float(((self.m_f + self.m_d) * h_p).sum())
+            self._kill(dead)                                     # the liquid of a dead patch moves to another patch...
             if self.m_f.size:                                    # ... which is at its own temperature
                 h_p = self.film_enthalpy(h_liq)
-                self._spread_to_patches(before - float((self.m_f * h_p).sum()), self.m_f)
+                self._spread_to_patches(before - float(((self.m_f + self.m_d) * h_p).sum()), self.m_f + self.m_d)
             n_dead += int(dead.size)
         self.solver.set_fractions(self.phi)
         self.solver.set_film_mass(self._film_nodal())
         self.frozen_mass += frozen
-        self.last_melt.update({"n_dead": n_dead, "released_mass": released, "feed_mass": feed_mass, "frozen_mass": frozen})
+        self.last_melt.update({"n_dead": n_dead, "released_mass": released, "feed_mass": feed_mass, "frozen_mass": frozen,
+                               "deep_liquid": deep_liquid})
 
     def _freeze_back(self, want, solid, h_e, h_p):
         """Give each element back the `want` kilograms of film that have fallen below the feed ramp, drawn from its
         own patches in proportion to what each wants to give (`solid`). An element can only recover film that is
         still there -- what ran off or sprayed away is gone -- and phi_e is capped at 1, because a patch may hold
         what its neighbours' runoff delivered and the mesh cannot grow a layer outside itself, so film with nowhere
-        to go simply stays film (`film_frozen_fraction` records how much)."""
+        to go simply stays film (`film_frozen_fraction` records how much). The deep liquid beneath the film (amendment
+        of 2026-10-02) lies against the solid, so it is drawn first; both hold the same liquid enthalpy on the same
+        nodes, so which of the two gives up the mass books nothing."""
         if not self.m_f.size or not want.any():
             return 0.0
         owner = self.surface.owner
@@ -488,12 +507,14 @@
         room = np.maximum(1.0 - self.phi, 0.0) * self.element_mass
         with np.errstate(divide="ignore", invalid="ignore"):
             share = np.where(have > 0.0, np.minimum(want, room) / np.where(have > 0.0, have, 1.0), 0.0)
-        taken = np.minimum(solid * share[owner], self.m_f)
+        taken = np.minimum(solid * share[owner], self.m_f + self.m_d)
         if not taken.any():
             return 0.0
         per_element = np.bincount(owner, taken, self.mesh.n_elements)
         self.phi = np.clip(self.phi + per_element / self.element_mass, 0.0, 1.0)
-        self.m_f = self.m_f - taken
+        from_deep = np.minimum(self.m_d, taken)
+        self.m_d = self.m_d - from_deep
+        self.m_f = self.m_f - (taken - from_deep)
         # the film arrives in the element it froze onto, which is colder than the surface it left
         self._defer_to_elements(np.bincount(owner, taken * h_p, self.mesh.n_elements) - per_element * h_e)
         return float(taken.sum())
@@ -522,13 +543,95 @@
 
     def _film_nodal(self):
         """The film's mass per node [kg]: each patch's film shared over its three nodes. The film's thermal state
-        lives on the nodes because that is where the solver's capacity and the solver's enthalpy live."""
+        lives on the nodes because that is where the solver's capacity and the solver's enthalpy live. The deep liquid
+        beneath the film rides the same nodes in the same way (amendment of 2026-10-02)."""
         if not self.m_f.size:
             return np.zeros(self.mesh.n_nodes)
-        return np.bincount(self.surface.faces.ravel(), np.repeat(self.m_f / 3.0, 3), self.mesh.n_nodes)
+        return np.bincount(self.surface.faces.ravel(), np.repeat((self.m_f + self.m_d) / 3.0, 3), self.mesh.n_nodes)
 
+    def _deep_runoff(self, dt, flow, delta_m, reach, T, h_node, h_e, h_p):
+        """(ii-a) The deep runoff (amendment of 2026-10-02): the contiguous liquid below the conjugate depth runs off
+        under the pressure gradient and the deceleration, and is never sprayed. Returns the deep liquid it saw [kg].
 
+        Plan fact 28(b) keeps melt in its element unless the gas shear can reach it, and fact 29 asked that the
+        pressure-gradient and deceleration-driven runoff integrate over the whole liquid depth rather than the film's.
+        Both hold here. The deep liquid of a patch is (a) the fully liquid inventory f_feed phi rho V of the elements of
+        its contiguous molten chain (`_molten_chain`, the march of `molten_depth`) that the shear does not reach -- the
+        ones the feed gate held back -- shared among the patches whose chains pass through them by patch area, plus
+        (b) the mobile deep liquid `m_d` the runoff has already delivered to it. With h_D that liquid's thickness and b
+        the film's, it moves by `film.deep_flux` G ((b + h_D)^3 - b^3) / (3 mu_l) on the same linearly implicit upwind
+        transport as the film, on windward patches under Girin's closure only: elsewhere there is no conjugate depth,
+        so there is no layer below the shear's reach. A patch that loses deep liquid gives it up in proportion from
+        its chain's elements and its `m_d`; what arrives goes to `m_d`. In situ the liquid stays where it is: an element
+        is debited only by what actually flows out of it, so the melt pool keeps its place in the mesh and is still
+        delivered to the film by the owner rule as the surface recedes to it.
+        Only the skin sprays. `m_d` is never offered to the spray: it becomes film only from the top, as far as the
+        film is thinner than the conjugate depth (the free surface has receded into it), and all at once where there
+        is no conjugate depth or on a leeward patch. Every transfer is booked by the rule of facts 5 and 25: liquid
+        leaves an element at the enthalpy its molten part carries, the element keeps the difference to its mean, the
+        liquid arrives at its destination patch's liquid enthalpy and the difference is released there; liquid moving
+        between `m_d` and the film of one patch moves at one temperature and books nothing. `deep_runoff_mass` counts
+        what left the elements, `deep_surfaced_mass` what became film from the top."""
+        liq, mat, mesh = self.liquid, self.material, self.mesh
+        n, areas = self.surface.n_patches, self.surface.areas
+        if not n:
+            return 0.0
+        girin = self.windward & np.isfinite(delta_m)          # a conjugate depth exists: liquid can lie below it
+        rows = e = share = None
+        D = self.m_d.copy()
+        if reach is not None:                                  # feed_depth "all" feeds every element: nothing is held back
+            _, chain = self._molten_chain()
+            rows, level = np.nonzero(chain >= 0)
+            e = chain[rows, level]
+            held = girin[rows] & (reach[e] <= 0.0)             # in a chain under a Girin patch, and the feed gate held it
+            rows, e = rows[held], e[held]
+            f = mat.feed_fraction(T)[mesh.tets].mean(axis=1) * mesh.active          # its fully liquid part (facts 5, 28)
+            liquid = np.minimum(f * self.phi, np.maximum(self.phi - PHI_MIN, 0.0)) * self.element_mass   # never below PHI_MIN
+            weight = np.bincount(e, areas[rows], mesh.n_elements)
+            share = liquid[e] * areas[rows] / np.where(weight[e] > 0.0, weight[e], 1.0)
+            D = D + np.bincount(rows, share, n)
+        seen = float(D.sum())
+        if self.runoff is not None and seen > 0.0:
+            b = self.m_f / (liq.rho * areas)
+            G = np.where(girin, flow.G, 0.0)
+            D_new, _, _ = self.runoff.transport(D, lambda hd: self._film_mod.deep_flux(G, b, b + hd, liq.mu),
+                                                    self.t_hat, liq.rho, areas, dt)
+            change = D_new - D
+            gain, loss = np.maximum(change, 0.0), np.maximum(-change, 0.0)
+            with np.errstate(divide="ignore", invalid="ignore"):
+                r = np.where(D > 0.0, np.minimum(loss / np.where(D > 0.0, D, 1.0), 1.0), 0.0)
+            from_m_d = r * self.m_d
+            carried, left = float((from_m_d * h_p).sum()), float(from_m_d.sum())
+            if e is not None and e.size:
+                out = np.bincount(e, r[rows] * share, mesh.n_elements)          # kg leaving each held element
+                fn = mat.feed_fraction(T)[mesh.tets]
+                fl = fn.mean(axis=1)
+                with np.errstate(divide="ignore", invalid="ignore"):            # per kg: what its molten part carries
+                    h_out = np.where(fl > 0.0, (fn * h_node[mesh.tets]).mean(axis=1) / np.where(fl > 0.0, fl, 1.0), h_e)
+                self.phi = self.phi - out / self.element_mass
+                self._defer_to_elements(out * h_e - out * h_out)               # the element keeps what the liquid left
+                carried += float((out * h_out).sum())
+                left += float(out.sum())
+                self.deep_runoff_mass += float(out.sum())                     # mobilised: what left the elements
+            # what arrives is exactly what left: the direct solve conserves only to its conditioning, 1e-12 to 2e-11 of
+            # the mass moved under these stiff coefficients (measured 2026-10-02), so the arrivals are scaled to the
+            # departures -- a correction of that size, which makes the books exact by construction
+            if gain.sum() > 0.0:
+                gain = gain * (left / float(gain.sum()))
+            self.m_d = self.m_d - from_m_d + gain
+            self._spread_to_patches(carried - float((gain * h_p).sum()), gain)  # ... which arrives at its patch's enthalpy
+        # the skin: deep liquid becomes film only from the top, as far as the film is thinner than the conjugate depth
+        with np.errstate(invalid="ignore"):
+            room = np.where(girin, np.maximum(liq.rho * areas * np.where(girin, delta_m, 0.0) - self.m_f, 0.0), np.inf)
+        up = np.minimum(self.m_d, room)
+        if up.any():
+            self.m_f = self.m_f + up
+            self.m_d = self.m_d - up
+            self.deep_surfaced_mass += float(up.sum())
+        return seen
 
+
+
     def _add_to_film(self, fed, carried):
         """Distribute each element's freed liquid to the film: owners to their own patches by area, interior elements
         to the NEAREST_PATCHES nearest patches by area. `carried` is the enthalpy that liquid takes with it [J per
@@ -567,10 +670,11 @@
         # tests below ask whether the gas shear reaches the bottom of the liquid, so they belong on this; the fluxes
         # and the release stay on the film, which is the mass that can actually move (decided 2026-09-22).
         molten = self.molten_depth()
+        deep = self.m_d / (liq.rho * areas)          # the deep liquid the runoff delivered: part of the layer, never film
         # (ii) lubrication and runoff
         n_sub, moved = 0, 0.0
         if self.runoff is not None:
-            q_of_b = lambda bb: self._film_mod.lubrication(flow.tau, flow.G, bb, delta_m, liq.mu, b_layer=bb + molten)[1]
+            q_of_b = lambda bb: self._film_mod.lubrication(flow.tau, flow.G, bb, delta_m, liq.mu, b_layer=bb + molten + deep)[1]
             before = self.m_f
             self.m_f, n_sub, moved = self.runoff.transport(self.m_f, q_of_b, self.t_hat, liq.rho, areas, dt)
             self.runoff_mass += moved                                              # mass that arrived on another patch
@@ -580,7 +684,7 @@
             d = self.m_f - before
             self._spread_to_patches(-float((d * h_p).sum()), np.maximum(d, 0.0))
         b = self.m_f / (liq.rho * areas)
-        layer = b + molten
+        layer = b + molten + deep
         v_s, q, _, thick = self._film_mod.lubrication(flow.tau, flow.G, b, delta_m, liq.mu, b_layer=layer)
         # the molten surface: the area over which a wave could form at all, either wetted by the film or molten in its
         # own right. Its contiguous extent is what every mode's wavelength is measured against (2026-09-24).
@@ -651,20 +755,26 @@
         return released
 
     def _kill(self, dead):
-        old_surface, old_m_f = self.surface, self.m_f
+        old_surface, old_m_f, old_m_d = self.surface, self.m_f, self.m_d
         gone, new = self.mesh.deactivate(dead)
         if self.mesh.n_active == 0:                              # nothing left: the run ends (demise) without a surface
-            self.consumed = True                                 # the film still on it leaves with the body
-            self.removed_mass += float(self.m_f.sum())
-            self.removed_enthalpy += float((self.m_f * self.film_enthalpy()).sum())
-            self.m_f = np.zeros(0)
+            self.consumed = True                                 # the liquid still on it leaves with the body
+            self.removed_mass += float(self.m_f.sum() + self.m_d.sum())
+            self.removed_enthalpy += float(((self.m_f + self.m_d) * self.film_enthalpy()).sum())
+            self.m_f, self.m_d = np.zeros(0), np.zeros(0)
             return
         self._refresh_geometry()
-        m_f = np.zeros(self.surface.n_patches)
         keep = self.patch_of_face[old_surface.face_ids]
         kept = keep >= 0
-        m_f[keep[kept]] += old_m_f[kept]
-        lost = ~kept & (old_m_f > 0.0)
+        self.m_f = self._hand_over(old_surface, old_m_f, keep, kept)
+        self.m_d = self._hand_over(old_surface, old_m_d, keep, kept)   # the deep liquid beneath it goes the same way
+
+    def _hand_over(self, old_surface, old, keep, kept):
+        """A per-patch liquid mass of the surface before the deaths, carried onto the current surface: kept patches keep
+        theirs; a vanished patch's goes to the faces its own element exposed (below)."""
+        m_f = np.zeros(self.surface.n_patches)
+        m_f[keep[kept]] += old[kept]
+        lost = ~kept & (old > 0.0)
         if lost.any():
             # The film on a vanished patch has not moved: the surface receded *into* it, so it now rests on the faces
             # of that same element which have just become boundary -- the face on its inward side, whose neighbour is
@@ -691,12 +801,12 @@
                 total = area.sum(axis=1)
             has = total > 0.0
             if has.any():
-                share = old_m_f[idx[has]][:, None] * area[has] / total[has][:, None]
+                share = old[idx[has]][:, None] * area[has] / total[has][:, None]
                 np.add.at(m_f, targets[has][ok[has]], share[ok[has]])
             if (~has).any():                                     # the element exposed nothing: its death opened a
                 orphan = idx[~has]                               # hole right through, so fall back to the nearest patch
-                np.add.at(m_f, self._patch_tree.query(old_surface.centroids[orphan])[1], old_m_f[orphan])
-        self.m_f = m_f
+                np.add.at(m_f, self._patch_tree.query(old_surface.centroids[orphan])[1], old[orphan])
+        return m_f
 
     # -- reporting -----------------------------------------------------------------------------------------------
     def on_current_surface(self, values, fill):
@@ -730,8 +840,9 @@
 
     def film_energy(self):
         """Enthalpy held by the film [J]. Identical to the nodal sum the solver carries (`set_film_mass` gives each
-        node sum_p m_p/3), so the film's sensible heat is inside 1^T M dT and the balance closes on it exactly."""
-        return float((self.m_f * self.film_enthalpy()).sum()) if self.m_f.size else 0.0
+        node sum_p m_p/3), so the film's sensible heat is inside 1^T M dT and the balance closes on it exactly. The
+        deep liquid beneath the film holds the same liquid enthalpy on the same nodes (amendment of 2026-10-02)."""
+        return float(((self.m_f + self.m_d) * self.film_enthalpy()).sum()) if self.m_f.size else 0.0
 
     def film_frozen_fraction(self):
         """Fraction of the film sitting on patches below the feed ramp: mass the enthalpy calls solid but the model
@@ -807,9 +918,16 @@
         This is the thickness that belongs in the thick/thin instability test against Girin's conjugate melt-layer
         depth delta_m: the film mass on a patch is the mobile inventory, not the depth of liquid beneath the wall.
         Measured on the 100 mm physics flight, the two differ by an order of magnitude and invert the branch choice."""
+        return self._molten_chain(Te, max_levels)[0]
+
+    def _molten_chain(self, Te=None, max_levels=8):
+        """`molten_depth`'s march, also returning the elements it passed through: (depth per patch [m], chain), the
+        chain an (n_patches, max_levels) array of element ids, level 0 the patch's owner and -1 past the chain's end.
+        The deep runoff needs the elements (amendment of 2026-10-02); the branch test needs only the depth."""
         mat, mesh = self.material, self.mesh
+        chain = np.full((self.surface.n_patches, max_levels), -1, dtype=np.int64)
         if not self.surface.n_patches or not mat.melts:
-            return np.zeros(self.surface.n_patches)
+            return np.zeros(self.surface.n_patches), chain
         if Te is None:
             Te = self.solver.temperature()[mesh.tets].mean(axis=1)
         molten = mesh.active & (Te >= mat.T_feed)
@@ -819,9 +937,10 @@
         alive = molten[cur]
         depth = np.zeros(self.surface.n_patches)
         rows = np.arange(len(cur))
-        for _ in range(max_levels):
+        for level in range(max_levels):
             if not alive.any():
                 break
+            chain[alive, level] = cur[alive]
             reach = np.einsum("ij,ij->i", centroid[cur] - self.surface.centroids, inward)      # projected depth so far
             depth = np.where(alive, np.maximum(depth, reach), depth)
             pair = mesh._face_elements[mesh.element_faces(cur)]                               # (n, 4, 2)
@@ -834,14 +953,16 @@
             alive = alive & good & molten[np.where(good, nxt, 0)]
             cur = np.where(good, nxt, cur)
         # the chain's last centroid sits half an element short of the far wall of that element: add that half
-        return np.maximum(depth, 0.0) * 1.5
+        return np.maximum(depth, 0.0) * 1.5, chain
 
     def liquid_layer_depth(self, Te=None):
-        """Depth of liquid at each patch [m]: its film plus the contiguous molten material beneath it. What the gas
-        shear sees, and so what decides whether the film is thick or thin against the conjugate melt-layer depth."""
+        """Depth of liquid at each patch [m]: its film plus the contiguous molten material beneath it, plus the deep
+        liquid the runoff has delivered there (amendment of 2026-10-02). What the gas shear sees, and so what decides
+        whether the film is thick or thin against the conjugate melt-layer depth."""
         if not self.m_f.size:
             return np.zeros(0)
-        return self.m_f / (self.liquid.rho * self.surface.areas) + self.molten_depth(Te)
+        a = self.liquid.rho * self.surface.areas
+        return self.m_f / a + self.molten_depth(Te) + self.m_d / a
 
     def film_thickness_max(self):
         """Thickest film [m] over the patches where a thickness means anything: at least a tenth of the median patch
@@ -858,6 +979,17 @@
         ok = (a >= 0.1 * np.median(a)) & (b <= np.sqrt(a))
         return float(b[ok].max()) if ok.any() else 0.0
 
+    def deep_blob_fraction(self):
+        """Fraction of the deep liquid (`m_d`) sitting deeper than its patch is wide (b > sqrt(A)), the deep account's
+        counterpart of `film_blob_fraction`. The deep runoff carries the liquid below the conjugate depth to where the
+        flow converges and nothing strips it there, so it piles up on the patches the patch graph cannot drain -- on
+        the 100 mm flight almost all of it (amendment of 2026-10-02). The mass is exact; the depth is not a film's."""
+        total = float(self.m_d.sum())
+        if total <= 0.0:
+            return 0.0
+        a = self.surface.areas
+        return float(self.m_d[self.m_d / (self.liquid.rho * a) > np.sqrt(a)].sum() / total)
+
     def film_blob_fraction(self):
         """Fraction of the film sitting deeper than its patch is wide (b > sqrt(A)) -- melt the lubrication film
         model cannot describe, because it no longer fits on the facet it is booked to. It is a geometry limit of the
@@ -933,6 +1065,9 @@
                 "molten_depth_mean_mm": lm.get("molten_depth_mean_mm", float("nan")),
                 "delta_m_mean_um": lm.get("delta_m_mean_um", float("nan")),
                 "thick_branch_fraction": lm.get("thick_branch_fraction", float("nan")),
+                "deep_liquid_kg": lm.get("deep_liquid", 0.0), "deep_mass_kg": float(self.m_d.sum()),
+                "deep_runoff_mass_kg": self.deep_runoff_mass, "deep_surfaced_mass_kg": self.deep_surfaced_mass,
+                "deep_blob_fraction": self.deep_blob_fraction(),
                 "removed_enthalpy_J": self.removed_enthalpy,
                 "film_thickness_max_mm": self.film_thickness_max() * 1e3, "film_thickness_mean_mm": self.film_thickness_mean() * 1e3,
                 "nose_radius_mm": self.nose_radius() * 1e3, "transverse_radius_mm": self.transverse_radius * 1e3,
```

Tests (the tested code, appended to `tests/test_reentry_model_melting.py`):

```diff
--- a/tests/test_reentry_model_melting.py
+++ b/tests/test_reentry_model_melting.py
@@ -275,3 +275,152 @@
     assert b.m_f[exposed] == pytest.approx(carried * w / w.sum())
     floor = exposed[int(np.argmax(cosine))]                         # the face most nearly parallel to the old patch
     assert b.m_f[floor] > 0.5 * carried and cosine.max() > 0.9      # takes the bulk of it
+
+
+# ---------------------------------------------------------------------------------------------------------------
+# Amendment of 2026-10-02: the liquid below the conjugate depth runs off, and only the skin sprays
+
+
+def girin_state(b, t=50.0):
+    """The 100 mm flight at 69.8 km (t = 50 s on the intact sphere's trajectory): the shock-layer branch, every windward
+    patch under Girin's closure with a conjugate depth of 0.17-0.30 mm, and G positive from the nose to about 80 deg
+    (measured 2026-10-02)."""
+    sim = simulator(b, t_max=90.0)
+    sim.advance(t)
+    return sim, sim.aero_state(sim.t, sim.y[:3], sim.y[3:])
+
+
+def molten_pool(b, T_pool=960.0, T_rest=850.0, depth=3e-3, theta_max=60.0):
+    """A molten pool: the nodes within `depth` of the surface and `theta_max` of the flight direction (+x) above the
+    liquidus, the rest of the body below it."""
+    p = b.mesh.points
+    r = np.linalg.norm(p, axis=1)
+    near = (r > 0.05 - depth) & (p[:, 0] > r * math.cos(math.radians(theta_max)))
+    b.solver.set_temperature(np.where(near, T_pool, T_rest))
+
+
+def flow_and_delta_m(b, a):
+    from reentry_model import spray
+    flow = b.flow.evaluate(a, b.theta, b.nose_radius(), b.liquid.rho, b.surface_temperature())
+    return flow, spray.melt_layer(flow, b.liquid)[0]
+
+
+def deep_stage_args(b, delta_m):
+    """What melt_step hands the deep stage: the shear-reach gate, the temperatures and the three enthalpies."""
+    T, mat = b.solver.temperature(), b.material
+    h_node = mat.enthalpy(T)
+    return b._shear_reaches(delta_m), T, h_node, h_node[b.mesh.tets].mean(axis=1), b.film_enthalpy(mat.enthalpy_liquid(T))
+
+
+def test_the_deep_runoff_moves_the_liquid_below_the_conjugate_depth_and_books_it_exactly(layered_mesh):
+    """Fact 28(b) holds melt in its element unless the gas shear reaches it; fact 29 asked that the pressure-gradient
+    and deceleration-driven flux integrate over the whole liquid depth. A pool 3 mm deep under the nose at 70 km: the
+    deep stage takes liquid only from the buried elements the feed gate held back (never from an owner, never from an
+    element the shear reaches, never below PHI_MIN), moves it outward as G pushes it there, and books mass and
+    enthalpy exactly -- what the elements lost is on the patches, and every enthalpy difference is in the deferred
+    loads. Only the skin is film afterwards: deep liquid became film only as far as the conjugate depth had room."""
+    from reentry_model import surface_flow as sf
+    b = melting_body(layered_mesh)
+    sim, a = girin_state(b)
+    molten_pool(b)
+    flow, delta_m = flow_and_delta_m(b, a)
+    assert flow.branch == sf.BRANCH_SHOCK_LAYER and np.isfinite(delta_m[b.windward]).all()
+    reach, T, h_node, h_e, h_p = deep_stage_args(b, delta_m)
+    phi0, mass0 = b.phi.copy(), b.mass(0.0)
+    E0 = b.energy() + b.pending_load.sum()
+    seen = b._deep_runoff(0.5, flow, delta_m, reach, T, h_node, h_e, h_p)
+    b.solver.set_fractions(b.phi)                                   # as melt_step does at its end
+    assert seen > 0.0 and b.deep_runoff_mass > 0.0                  # liquid lay below the conjugate depth, and it moved
+    assert b.deep_runoff_mass == pytest.approx(((phi0 - b.phi) * b.element_mass).sum(), rel=1e-12)   # what left the elements
+    assert b.deep_surfaced_mass == pytest.approx(b.m_f.sum(), rel=1e-12)                              # what became film
+    assert b.mass(0.0) == pytest.approx(mass0, rel=1e-12)
+    assert b.energy() + b.pending_load.sum() == pytest.approx(E0, rel=1e-12)
+    lost = phi0 - b.phi
+    owners = b.owner_area > 0.0
+    assert (lost >= 0.0).all() and lost.sum() > 0.0 and not lost[owners].any() and not lost[reach > 0.0].any()
+    assert (b.phi[lost > 0.0] >= body.PHI_MIN * (1.0 - 1e-12)).all()
+    left = lost * b.element_mass                                    # where it came from and where it went: outward
+    c = b.mesh.points[b.mesh.tets].mean(axis=1)
+    theta_from = np.arccos(c[:, 0] / np.linalg.norm(c, axis=1))
+    arrived = b.m_d + b.m_f
+    assert (arrived * b.theta).sum() / arrived.sum() > (left * theta_from).sum() / left.sum()
+    skin = b.liquid.rho * b.surface.areas * np.where(b.windward, delta_m, np.inf)
+    deep = b.m_d > 0.0
+    assert deep.any() and (b.m_f[deep] >= skin[deep] * (1.0 - 1e-12)).all() and (b.m_f <= skin * (1.0 + 1e-12)).all()
+    assert not b.m_d[~b.windward].any()                             # no conjugate depth on the lee: no deep liquid there
+
+
+def test_the_deep_liquid_is_never_sprayed_only_the_skin_is(layered_mesh):
+    """A millimetre of deep liquid on every windward patch and no film, over a body with nothing molten beneath its
+    surface: the deep stage tops the film up to the conjugate depth from the top of the deep liquid and leaves the
+    rest deep; spraying then strips only the film and leaves the deep liquid exactly as it was."""
+    b = melting_body(layered_mesh)
+    sim, a = girin_state(b)
+    b.solver.set_temperature(905.0)                                 # below T_feed: no element holds liquid
+    flow, delta_m = flow_and_delta_m(b, a)
+    liq = b.liquid
+    b.m_d = np.where(b.windward, 1e-3 * liq.rho * b.surface.areas, 0.0)
+    reach, T, h_node, h_e, h_p = deep_stage_args(b, delta_m)
+    total = b.m_d.sum()
+    seen = b._deep_runoff(0.5, flow, delta_m, reach, T, h_node, h_e, h_p)
+    assert seen == pytest.approx(total) and b.m_f.sum() + b.m_d.sum() == pytest.approx(total, rel=1e-12)
+    skin = liq.rho * b.surface.areas * np.where(b.windward, delta_m, np.inf)
+    deep = b.m_d > 0.0
+    assert deep.any() and b.m_f[deep] == pytest.approx(skin[deep], rel=1e-12)       # where any stays deep, the skin is full
+    film0, deep0 = b.m_f.copy(), b.m_d.copy()
+    assert b.deep_surfaced_mass == pytest.approx(film0.sum(), rel=1e-12) and b.deep_runoff_mass == 0.0   # none from elements
+    assert 0.0 <= b.deep_blob_fraction() <= 1.0
+    released = b._film_and_spray(sim.t, 0.5, a, h_p, flow, delta_m)
+    assert 0.0 < released <= film0.sum() * (1.0 + 1e-12)            # the skin sprays ...
+    assert np.array_equal(b.m_d, deep0)                             # ... and not a gram of the deep liquid
+    assert b.m_f.sum() + released == pytest.approx(film0.sum(), rel=1e-12)
+
+
+def test_nothing_changes_where_no_liquid_lies_below_the_conjugate_depth(layered_mesh):
+    """The deep runoff adds a transfer; it changes nothing else. At 70 km, with Girin's closure on every windward patch
+    and a surface just above the liquidus over a body that is not fully molten beneath its owner elements, there is
+    no liquid below the conjugate depth: a whole melt step is bit-identical with the deep runoff on and off. (Where the
+    closure is not Girin's there is no conjugate depth at all; the 50 mm flight, which never has it, is the
+    flight-level check, measured 2026-10-02.)"""
+    runs = []
+    for deep in (True, False):
+        b = melting_body(layered_mesh, deep_runoff=deep)
+        sim, a = girin_state(b)
+        boundary = np.zeros(b.mesh.n_nodes, dtype=bool)
+        boundary[np.unique(b.surface.faces)] = True
+        b.solver.set_temperature(np.where(boundary, 930.0, 880.0))  # owners molten, nothing beneath them
+        b.melt_step(sim.t, 0.5, a)
+        runs.append(b)
+    on, off = runs
+    assert on.last_melt["deep_liquid"] == 0.0 and on.deep_runoff_mass == 0.0 and on.sprayed_mass > 0.0
+    for name in ("phi", "m_f", "m_d", "pending_load"):
+        assert np.array_equal(getattr(on, name), getattr(off, name)), name
+    assert on.sprayed_mass == off.sprayed_mass and on.removed_enthalpy == off.removed_enthalpy
+    assert on.runoff_mass == off.runoff_mass and np.array_equal(on.last_delta_m, off.last_delta_m, equal_nan=True)
+
+
+def test_full_steps_with_the_deep_runoff_keep_the_books_exact(layered_mesh):
+    """Six coupled macro steps from 70 km with a molten pool under the nose and the physics loads: the deep runoff
+    acts (liquid below the conjugate depth is found and moved), elements die with deep liquid on their patches, and
+    the energy balance and the mass stay exact by the same measures as every other melting test."""
+    b = melting_body(layered_mesh)
+    sim, a = girin_state(b)
+    molten_pool(b)
+    b.energy0 = b.energy()
+    model = heating.PhysicsHeating()
+    seen = 0.0
+    for _ in range(6):
+        sim.advance(0.5)
+        a = sim.aero_state(sim.t, sim.y[:3], sim.y[3:])
+        loads = model.evaluate(a, b.theta, b.surface_temperature(), b.nose_radius(), T_mean=b.mean_temperature())
+        b.advance(sim.t, 0.5, loads, state=a)
+        seen += b.last_melt["deep_liquid"]
+        assert abs(b.energy_balance_residual()) < 1e-7
+        assert b.mass(0.0) + b.removed_mass == pytest.approx(b.mass0, rel=1e-12)
+        assert b.m_d.shape == (b.surface.n_patches,) and (b.m_d >= 0.0).all()
+        assert b.solver.film_mass.sum() == pytest.approx(b.m_f.sum() + b.m_d.sum())   # the solver carries both
+    assert seen > 0.0 and b.deep_runoff_mass > 0.0 and b.sprayed_mass > 0.0 and b.mesh.n_active < layered_mesh.n_elements
+    stats = b.melt_stats()
+    assert stats["deep_runoff_mass_kg"] == b.deep_runoff_mass and stats["deep_mass_kg"] == pytest.approx(b.m_d.sum())
+    assert stats["deep_surfaced_mass_kg"] == b.deep_surfaced_mass and 0.0 <= stats["deep_blob_fraction"] <= 1.0
+    assert stats["deep_liquid_kg"] == b.last_melt["deep_liquid"]
```

### Measured (2026-10-02)

- **Tests.** `tests/test_reentry_model_melting.py`: 14 passed, the ten above and the four new ones, each of
  which fails on the unamended body (it has no deep stage, no `deep_runoff` setting and no deep diagnostics to call).
  They pin that the deep stage takes liquid only
  from buried elements the feed gate held back, never from an owner or below `PHI_MIN`, books mass and enthalpy
  exactly (to 1e-12 of the body's mass and of its stored energy plus the deferred loads), moves the liquid the way G
  points and leaves film only up to the conjugate depth; that a supercritical skin over a millimetre of deep liquid
  sprays while the deep account is untouched to the bit; that a whole melt step is bit-identical with the deep runoff
  on and off when no liquid lies below the conjugate depth; and that six coupled steps over a molten pool keep the
  energy balance below 1e-7 and the mass to 1e-12. With the film, coupled and CLI tests of this amendment the whole
  unit tier gives 223 passed, 1 skipped, 2 failed and 3 errors — the prototype's own 216 passed plus the seven new
  tests, the failures and errors being the five known ones from the missing melting SESAM references (Task 11).
  `tests/test_reentry_model_fenicsx.py` passes in `fenicsx_env` (9 passed).
- **Unchanged where nothing lies below the conjugate depth** (fact 49): the 50 mm flight, which never has Girin's
  closure, is bit-identical with the deep runoff on, and the 100 mm flight is bit-identical with `--deep-runoff off`,
  both runs of each pair having numpy's random seed fixed (fact 52).
- **The 100 mm flight to 120 s** (fact 49): the deep runoff takes 0.285 kg out of the elements; sprayed mass rises
  7.1 % (1.018 to 1.091 kg) and the body at 120 s is 16 % lighter (0.454 to 0.381 kg), about fifty times the
  run-to-run spread; the median droplet radius does not move, the largest droplet grows 29 %, the front-surface
  Rayleigh–Taylor mode releases 37 % more, and the energy balance stays exact (−1.1e-10).
- **The whole 100 mm flight** (fact 49): at the default step 22.4 % instead of 30.8 % of the initial mass reaches the
  ground (0.330 against 0.453 kg), the deep liquid still present at 120 s surfacing and spraying by 229.5 s; every
  gram mobilised is accounted for at landing (285.1 g taken, 277.1 g became film, 8.0 g froze back).
- **Where the deep liquid goes** (fact 47): almost all of it into piles deeper than 2 mm and deeper than their facets
  are wide, in the craters at the eroding front — the patch graph's sinks.
- **Conservation** (fact 48): the stiff transport conserves only to 1e-14–2e-11 of what it moves per call, so the
  deep stage scales its arrivals to its departures.
- **Both backends** (fact 49): with the deep runoff active they agree to 1.3e-9 in mass and to 0.021 K, with both
  energy balances exact.
- **The macro step** (fact 50): the liquid held below the conjugate depth is a backlog of the melt step's order and
  vanishes as the step shrinks — a median 4.75 g, 0.23 g and 0.004 g at 0.5, 0.25 and 0.125 s — and the deep
  runoff's effect on the sprayed mass with it: +7.1 %, −0.3 % and −0.16 %, the last within the run-to-run spread.
  The same backlog sets the molten layer the branch test reads, so the unamended model's droplet population is not
  converged either (droplet count 1.72e7, 2.82e7 and 5.85e7). The default is the first open decision (fact 53).

### Not done here — for Asha to decide (fact 53)

Fact 53 sets them out with the options and a recommendation. First and most important, the time step: the liquid the
deep runoff moves is a backlog left by the order of the melt step, which vanishes as the step shrinks, and the same
backlog keeps the droplet population of the unamended model from converging (fact 50). The recommendation is to let the
surface recede through molten material by more than one element per step, as a Step 3 amendment of its own, and until
then to make `--deep-runoff off` the default; the amendment as built has it `on`. Then: what to do about the piles (keep
them as zone 3 for the large-fragment model to read, level them by surface tension, move liquid beyond the film limit to
a separate account, or re-measure after Task 16 first); how fast deep liquid should become skin; whether the
front-surface Rayleigh–Taylor mode may release it; whether anything should move under the Couette closure; whether the
shear passed down through the skin belongs in the deep flux; the film/deep split of the pressure-driven flux; the
laminar flux at Reynolds numbers of 10³–10⁴; and seeding numpy so that runs are reproducible.

## Amendment of 2026-10-03 — the molten cascade: the surface recedes through molten material within the step

> Asha's decision of 2026-10-02 on fact 53 (a): option (3), fix the cause of the molten backlog by letting the surface
> recede through more than one molten element per macro step. Measured in a throwaway copy of the prototype — the copy
> of the deep-runoff amendment above, verified before any change to be reproduced byte for byte by applying that
> amendment's nine diff blocks to a fresh `prototype/proto3/` — and the code below is the tested code, as a diff against
> that copy (the extracted body below with the fact-37 code of 2026-09-27 and the deep-runoff diff above applied).
> Facts 54–61 in `00-shared-context.md` carry the measurements. Read facts 4, 5, 25, 28 and 50 first. The copy keeps
> `PHI_DEATH = 0.05`: fact 44's 0.50 (this task's amendment of 2026-09-27) exists only as a plan, and the cascade reads
> the constant symbolically, so it works unchanged with either value (fact 61 (d) asks for the study to be repeated with
> 0.50).

**The cause, as fact 50 stated it.** The melt step runs after the conduction. Between two melt steps the heat flux
therefore melts the material beneath and beside the wall-owning element as well as the element itself; the feed hands
the owner's liquid to the film and the owner dies once it is consumed, but the element its death exposes was buried a
moment earlier, so the feed gate of fact 28(b) held its liquid in place, and it waited for the next step's feed. The
surface could recede through molten material by only one element per macro step, and whatever melted below the owner
waited as a backlog that shrank with the step.

**The design: one more rule inside the death loop of step (v).** After each pass of the existing death cascade
(fact 4), the elements that pass's deaths have just made wall owners are examined. Every one that is **fully molten**
and would survive the pass (φ > `PHI_DEATH`) is fed whole within the step (`_feed_exposed`): its remainder φρV leaves at
its mean nodal enthalpy h_e and arrives on the faces it now owns at their liquid enthalpy, the difference released
there; φ is set to zero. The next pass of the same loop kills it, `_kill` hands that liquid down to the faces it
exposes exactly as it hands down any dying patch's film, and the elements those deaths expose are examined in turn.
The loop ends when no death exposes a fully molten element. `MAX_CASCADE_PASSES = 32` caps the passes that feed in
one step; elements left by the cap wait for the next step's feed, as every exposed element did before, and
`cascade_capped_steps` counts the steps in which it bound. Two diagnostics are recorded: `cascade_passes` (the passes
that fed something this step) and `cascade_mass_kg` (the cumulative mass fed). `MeltSettings.molten_cascade` and
`--molten-cascade on|off` (default on) switch it; off reproduces the deep-runoff amendment bit for bit.

**What "fully molten" means, and why it is not "every node above the feed ramp".** The model already has a definition:
`molten_depth`, whose march sets the molten layer that the thick/thin branch test reads and the chain the deep runoff
drains, "stops at the first element that is not fully molten", meaning an element whose mean nodal temperature is below
T_feed (the top of the feed ramp, 910 K for `AA7075_range`). The cascade uses exactly that test, so it removes exactly
the layer fact 50 identified as the backlog, wherever a death exposes it. The first version of this amendment used the
stricter test that all four nodes be at or above T_feed, so that the owner feed of step (i) would take the whole
element by its own rule; measured on the 100 mm flight, **no element at the melt front ever passes that test**. Of the
1660–1920 elements in the molten chains below the wall-owning elements at 50.5, 60.5 and 70.5 s, not one has all four
nodes above 910 K: 43–45 % have one node below it and 50–51 % two, the coldest node sitting at a median of 897 K
(93 % liquid by the enthalpy, f_l = (897 − 750)/158), while the elements' mean temperatures run from 910 to 937 K
(median 918 K). Ninety-seven per cent of those cold nodes are shared with a wall-owning element, and about half of them
lie on the surface itself. The reason is the feed's own energy debit (fact 5): an element being fed loses the enthalpy
of its molten part and keeps the colder rest, which cools its nodes, so every node an element shares with an element
being fed sits on or just below the feed ramp. The strict cascade therefore fed 32 g in 120 s (in 129 steps, one or
two passes each) and left the backlog as it was — held liquid 4.41 g against 4.75 g, molten layer 1.21 mm against
1.26 mm, 16.7 million droplets against 17.2 million — so it was replaced by the model's own definition before any
other measurement was made. (Measured 2026-10-02 with the deep runoff off; with it on, 5.9 g fed and no change
beyond the run-to-run spread.)

**How it is booked (facts 4, 5 and 25).** The element leaves whole at its mean nodal enthalpy h_e, which is what the
solver holds for it, so nothing is booked on the element. Where every node is above T_feed this is the owner feed of
step (i) exactly — the molten part is then the whole element and its enthalpy is h_e. Where a node or two lies on or
just below the ramp it is the death rule of fact 4 applied to the whole element: the film holds the *liquid* enthalpy
(`enthalpy_liquid`, fact 25), so the residual latent heat of the cold nodes is paid by the faces the liquid lands on,
as a deferred load through `pending_load` and `LOAD_DT_MAX` (fact 26). Mass moves once, from φ to the film; energy is
booked once, as the difference, at the arrival. Neither of the two implementations rejected in facts 5 and 25 is
reintroduced: the element is debited (φ falls to zero) and the film is booked at the liquid enthalpy.

**What does not change.** The feed gate of fact 28(b) is respected by construction: the cascade feeds an element only
once a death has made it a wall owner, which the gate always admits, and it stops at the first exposed element that is
not fully molten, so molten material that cooler material separates from the wall is never fed. The deep runoff runs
before the spray on the liquid the gate holds; the cascade runs after it on what the step's deaths expose. The spray is
untouched: the cascade's liquid joins the film after the spray stage, like the death remainders of fact 4, and sprays
at the next step. Freeze-back runs before the death loop and is untouched. `removal = "instant"` feeds every element as
it melts, so nothing waits and the cascade is skipped there by construction. With no fully molten element exposed, a
melt step is bit-identical with the cascade on and off (unit test).

### Alternatives considered, and why they were not taken

- **Feed the newly exposed owners at the next step's start.** Before the next conduction step this is the same
  operation moved across the step boundary — the same temperature field and the same masses — so it only splits the
  death logic between two places and leaves the surface recorded at the end of each step (the history row, the frame,
  the drag area the trajectory flies with) one recession behind the mass account. Placed instead at the top of the next
  melt step, after the conduction and before the surface flow is evaluated, it would let that step's spray see the
  receded surface and its liquid at once; but then every owner the conduction left fully molten would also have to die
  before the flow is evaluated, moving those deaths from after the spray to before it on every melting step — a change
  wherever no molten element is ever exposed — or else the surface flow would have to be evaluated twice per step on
  two surfaces (it is the largest part of the melt step's cost, fact 23). Rejected.
- **An implicit recession**: re-solving the conduction on the receded surface within the step, or tracking the melt
  front inside the conduction solve. The latent heat of the exposed molten elements was already paid by this step's
  conduction, so moving their liquid needs no re-solve to keep mass and enthalpy exact; a re-solve per pass would cost a
  conduction step (0.6 s on the default mesh) per pass and change the operator splitting of every step, which is what a
  shorter macro step does anyway. Its one benefit — the rest of the step's heat reaching the receded surface directly —
  is the remaining splitting error, and the time-step study measures it. Rejected as not the smallest change.
- **Feed the whole molten chain under a dying owner at step (i), before the spray.** The spray would then see the
  liquid in the same step, but the chain's elements are still buried when they are fed, which is exactly what fact
  28(b)'s gate forbids, and it needs a prediction of a death that has not yet happened. Rejected.
- **Feed only elements whose every node is above T_feed** (the first version). Measured to do nothing on this flight,
  because the feed's own energy debit keeps every node an element shares with an element being fed on the ramp (the
  design section above). Replaced by the model's own definition of fully molten.
- **Feed partially molten exposed elements their molten fraction f φ only.** The element survives with (1 − f)φ, so the
  cascade stops at it anyway, and it changes results wherever an exposed element is molten in part — the case the
  amendment must leave alone. Rejected.
- **Kill the whole molten region below a dying owner in one pass.** The hand-over rule sends a vanished patch's liquid
  to the faces its own element exposed; when those faces vanish in the same batch they are not patches, and the liquid
  would fall back to the nearest patch — the concentrator fact 27 removed. Pass by pass, the rule stays exact. Rejected.

Code (the tested change to `reentry_model/body.py`):

```diff
--- a/reentry_model/body.py
+++ b/reentry_model/body.py
@@ -144,6 +144,8 @@
 PHI_DEATH = 0.05                 # a patch owner below this fraction dies (its remainder goes to the film): keeps the surface
                                  # nodes' thermal mass above 5 % of an element's, which the Newton iteration needs (measured
                                  # 2026-09-20: owners at 1e-3 made surface nodes swing by hundreds of K between iterates)
+MAX_CASCADE_PASSES = 32          # safety cap on the molten cascade's passes in one macro step; what it leaves waits for the
+                                 # next step, as every exposed element did before the cascade (amendment of 2026-10-03)
 REMOVAL_NAMES = ("girin", "instant")
 
 
@@ -157,6 +159,9 @@
                                       # to the windward cap; initial: D0 and R0 throughout (SESAM's convention, for the devices)
     feed_depth: str = "conjugate"     # conjugate: only liquid the gas shear reaches leaves its element (fact 28); all: any depth
     deep_runoff: bool = True          # the liquid below the conjugate depth runs off, never sprayed (amendment of 2026-10-02)
+    molten_cascade: bool = True       # a death that exposes a fully molten element (mean T >= T_feed) feeds it within the
+                                      # step, so the surface recedes through molten material by more than one element per
+                                      # macro step (amendment of 2026-10-03)
 
     def __post_init__(self):
         if self.removal not in REMOVAL_NAMES:
@@ -195,8 +200,12 @@
     (raising phi_e, capped at 1) once the patch falls back through the ramp -- the mirror of the feed rule, netted
     with it so that no element both melts and freezes in one step; (v) patch owners at phi <= PHI_DEATH die: the
     mesh's active set, the surface, the film (handed to the nearest surviving patch) and the solver's fractions are
-    refreshed. With `removal = "instant"` the liquid leaves the body at h_liquid instead, and the difference to the
-    element's own enthalpy is a nodal load on its nodes over the next step.
+    refreshed; with the molten cascade (`_feed_exposed`, amendment of 2026-10-03) an element a death exposes fully
+    molten (its mean temperature at or above T_feed, `molten_depth`'s test) is fed whole within the step, so that it
+    dies in turn and the surface recedes through molten material by as many elements as are molten rather than by one
+    per macro step. With `removal = "instant"`
+    the liquid leaves the body at h_liquid instead, and the difference to the element's own enthalpy is a nodal load on
+    its nodes over the next step.
 
     The film has the surface's own temperature, because it is thermally thin: q b / k_l = 0.22 K across a 10 um film
     at 2 MW/m2 (22 K even at 1 mm) and b^2/alpha = 0.3 ms against the 0.5 s macro step. Rather than give it an energy
@@ -242,6 +251,8 @@
         self.sprayed_mass = self.runoff_mass = self.removed_mass = self.removed_enthalpy = self.frozen_mass = 0.0
         self.deep_runoff_mass = 0.0                                         # kg the deep runoff took out of the elements
         self.deep_surfaced_mass = 0.0                                       # kg of deep liquid that became film from the top
+        self.cascade_mass = 0.0                                             # kg the molten cascade fed (2026-10-03)
+        self.cascade_capped_steps = 0                                       # steps whose cascade stopped at MAX_CASCADE_PASSES
         self.n_released = 0.0
         self.source_rows = []
         self.hist_n, self.hist_m = np.zeros(spray_mod.N_BINS), np.zeros(spray_mod.N_BINS)
@@ -463,8 +474,13 @@
                 deep_liquid = self._deep_runoff(dt, flow, delta_m, reach, T, h_node, h_e, h_p)
             released = self._film_and_spray(t, dt, state, h_p, flow, delta_m)
         frozen = 0.0 if s.removal == "instant" else self._freeze_back(want, solid, h_e, h_p)
-        # (v) death of consumed patch owners; a death exposes its neighbours, which die in turn if they are consumed
-        n_dead = 0
+        # (v) death of consumed patch owners; a death exposes its neighbours, which die in turn if they are consumed --
+        # and, with the molten cascade (amendment of 2026-10-03), are fed within the step if the death exposed them
+        # fully molten (`_feed_exposed`), so that they die in the next pass and expose the next
+        n_dead, passes, cascade, capped = 0, 0, 0.0, False
+        # fully molten as `molten_depth` has it -- the element's mean nodal temperature at or above T_feed -- so the cascade
+        # removes exactly the contiguous molten layer the branch test reads, wherever a death exposes it
+        full = T[tets].mean(axis=1) >= mat.T_feed if s.molten_cascade and s.removal != "instant" else None
         while not self.consumed:
             dead = np.flatnonzero((self.owner_area > 0.0) & self.mesh.active & (self.phi <= PHI_DEATH))
             if dead.size == 0:
@@ -481,16 +497,27 @@
                 self._defer_to_patches(carried - delta * h_p)
             self.phi[dead] = 0.0
             before = float(((self.m_f + self.m_d) * h_p).sum())
+            was_owner = self.owner_area > 0.0
             self._kill(dead)                                     # the liquid of a dead patch moves to another patch...
             if self.m_f.size:                                    # ... which is at its own temperature
                 h_p = self.film_enthalpy(h_liq)
                 self._spread_to_patches(before - float(((self.m_f + self.m_d) * h_p).sum()), self.m_f + self.m_d)
             n_dead += int(dead.size)
+            if full is not None and not self.consumed:          # the molten cascade: walls the deaths exposed fully molten
+                exposed = np.flatnonzero(full & (self.owner_area > 0.0) & ~was_owner & (self.phi > PHI_DEATH))
+                if exposed.size and passes < MAX_CASCADE_PASSES:
+                    cascade += self._feed_exposed(exposed, h_e, h_p)
+                    passes += 1
+                elif exposed.size:                               # the cap: these wait for the next step's feed
+                    capped = True
         self.solver.set_fractions(self.phi)
         self.solver.set_film_mass(self._film_nodal())
         self.frozen_mass += frozen
+        self.cascade_mass += cascade
+        self.cascade_capped_steps += int(capped)
         self.last_melt.update({"n_dead": n_dead, "released_mass": released, "feed_mass": feed_mass, "frozen_mass": frozen,
-                               "deep_liquid": deep_liquid})
+                               "deep_liquid": deep_liquid, "cascade_mass": cascade, "cascade_passes": passes,
+                               "cascade_capped": capped})
 
     def _freeze_back(self, want, solid, h_e, h_p):
         """Give each element back the `want` kilograms of film that have fallen below the feed ramp, drawn from its
@@ -519,6 +546,40 @@
         self._defer_to_elements(np.bincount(owner, taken * h_p, self.mesh.n_elements) - per_element * h_e)
         return float(taken.sum())
 
+    def _feed_exposed(self, exposed, h_e, h_p):
+        """(v-a) The molten cascade (amendment of 2026-10-03): feed the fully molten elements a death has just exposed,
+        within the step. Returns the mass fed [kg].
+
+        The melt step runs after the conduction, so between two melt steps the heat flux melts the material beneath
+        and beside the wall-owning element as well as the element itself. The feed hands the owner's liquid to the film
+        and the owner dies once consumed, but the element its death exposes was buried a moment earlier, so the feed
+        gate of fact 28(b) held its melt in place, and it used to wait for the next step's feed: the surface receded
+        through molten material by one element per macro step, and whatever melted below the owner waited as a backlog
+        that shrank with the step (plan fact 50). Exposed, the element is a wall owner -- in the sheared layer by
+        definition, which is what the gate admits -- so it is fed now.
+
+        "Fully molten" is `molten_depth`'s definition, the element's mean nodal temperature at or above T_feed, so the
+        cascade removes exactly the contiguous molten layer the branch test reads. It is not "every node above T_feed":
+        the feed's own energy debit keeps every element it is feeding on the feed ramp, so at the melt front no element
+        has all four nodes above it -- on the 100 mm flight not one of the 1700-1900 elements of the chains below the
+        wall-owning elements does, their coldest node sitting at a median 897 K, 93 % liquid by the enthalpy, and a
+        cascade on that test fed 32 g in 120 s (deep runoff off) and left the backlog as it was (measured 2026-10-02).
+        The element is fed whole: its remainder phi rho V leaves at its mean nodal enthalpy h_e and nothing is booked on
+        it, and it arrives on the faces it owns at their liquid enthalpy with the difference released there -- the owner
+        rule of step (i) where every node is above T_feed (the molten part is then all of it), and the death rule of fact
+        4 for the rest, whose residual latent heat the faces it lands on pay (facts 5 and 25). Left at phi = 0, it dies
+        in the next pass of the death loop and hands that liquid down to the faces it exposes, like any dying patch's
+        film; the cascade ends where the elements a death exposes are no longer fully molten -- the contiguity rule of
+        `molten_depth` -- so molten material that cooler material separates from the wall is never fed. Only elements
+        that would survive the pass are fed: one at or below PHI_DEATH dies anyway and leaves its remainder by the death
+        rule, as before."""
+        fed = np.zeros(self.mesh.n_elements)
+        fed[exposed] = self.phi[exposed] * self.element_mass[exposed]
+        self.phi[exposed] = 0.0
+        delta, carried = self._add_to_film(fed, fed * h_e)     # the whole element, at its mean nodal enthalpy
+        self._defer_to_patches(carried - delta * h_p)          # ... and the melt arrives at its patch's temperature
+        return float(fed.sum())
+
     def _defer_to_elements(self, energy, elements=None):
         """Book a per-element energy [J] on the elements' nodes, to be applied over the next step."""
         e = np.asarray(energy, dtype=float)
@@ -1068,6 +1129,7 @@
                 "deep_liquid_kg": lm.get("deep_liquid", 0.0), "deep_mass_kg": float(self.m_d.sum()),
                 "deep_runoff_mass_kg": self.deep_runoff_mass, "deep_surfaced_mass_kg": self.deep_surfaced_mass,
                 "deep_blob_fraction": self.deep_blob_fraction(),
+                "cascade_passes": float(lm.get("cascade_passes", 0)), "cascade_mass_kg": self.cascade_mass,
                 "removed_enthalpy_J": self.removed_enthalpy,
                 "film_thickness_max_mm": self.film_thickness_max() * 1e3, "film_thickness_mean_mm": self.film_thickness_mean() * 1e3,
                 "nose_radius_mm": self.nose_radius() * 1e3, "transverse_radius_mm": self.transverse_radius * 1e3,
```

**Tests** (the tested code below; each new test fails on the copy without the cascade, which has no `molten_cascade`
setting, no `_feed_exposed` and no cascade diagnostics):

1. *A molten column empties within one step*
   (`test_the_molten_cascade_empties_a_fully_molten_column_within_one_step`): a pool 4 mm deep under the nose at 960 K
   over a body at 850 K, one melt step. With the cascade the surface recedes through the whole pool in six passes —
   none of the 1 106 walls the step's deaths exposed is fully molten afterwards — every element the cascade fed was
   fully molten, and all of the 45.9 g it fed is in the film; 44 of the 7 842 fully molten elements remain, under
   owners the feed left partly consumed. Without it the 1 098 walls the first deaths exposed are all still fully molten
   and 6 764 fully molten elements remain. Mass to 1e-12 and enthalpy (stored plus queued) to
   1e-12.
2. *Nothing changes when the exposed elements are not fully molten*
   (`test_the_molten_cascade_changes_nothing_where_the_exposed_elements_are_not_fully_molten`): at 69.8 km under
   Girin's closure, with the surface nodes at 960 K over a body at 850 K and a millimetre of deep liquid on every
   windward patch, the wall-owning elements are fed and die, and every element their deaths expose has a mean
   temperature of at most 905 K; the deep runoff moves its liquid and the film sprays; the whole melt step is
   bit-identical with the cascade on and off (`phi`, `m_f`, `m_d`, `pending_load`, the active set, the sprayed, runoff
   and deep masses and the removed enthalpy).
3. *The gate* (`test_the_molten_cascade_feeds_only_the_wall_and_stops_at_material_that_is_not_fully_molten`): a fully
   molten skin 2 mm deep, an 880 K layer from 2 to 10 mm and a fully molten core below it, with the deep runoff on and
   off: every element the cascade fed was a wall owner at that moment and fully molten, the skin above the cool layer
   is gone, and every one of the core's elements keeps φ exactly.
4. *The deep liquid is still never sprayed* (`test_with_the_molten_cascade_the_deep_liquid_is_still_never_sprayed`):
   both mechanisms acting in one step, checked stage by stage — the spray leaves the deep account bit-identical, every
   death moves deep liquid without changing its total (to 1e-12), and each cascade feed lands in the film alone.
5. *Exact books over coupled steps* (`test_full_steps_with_the_molten_cascade_keep_the_books_exact`): six coupled steps
   from 69.8 km over a 6 mm pool with the physics loads, the cascade (at least two passes in a step), the deep runoff
   and the spray all acting: energy balance below 1e-7 of the absorbed heat and mass to 1e-12 at every step, the solver
   carrying the film and the deep liquid, the two history columns consistent with the body.
6. *Both backends* (`tests/test_reentry_model_fenicsx.py::test_the_molten_cascade_matches_the_skfem_backend`, in
   `fenicsx_env`): the same six steps on the coarse mesh, each backend seeded; identical passes per step (4, 3, 16, 4,
   5, 4) and active sets, mass to 1.3e-7, sprayed mass to 6.3e-7, cascade mass to 2.4e-7, the deep account to 3.7e-6,
   temperatures to 0.11 K, both balances exact (4.7e-10 and 2.4e-11); with the deep runoff off, 1e-10 and 4e-6 K.
7. *One existing test is pinned to the cascade off*: `test_film_temperature_freeze_back_and_the_netted_transfer` puts
   2.2 MW/m² on every facet — crater walls and lee included — of a body 50 K below its melting point. With the cascade
   the whole body melts within the 8 s of heating (consumed at step 16, the balance exact throughout), so there is no
   film left to freeze back, which is what that test is about. It now passes `molten_cascade=False`, with the reason in
   a comment; nothing else in it changes.

Tests (the tested code, appended to `tests/test_reentry_model_melting.py`, with the one-line change to the freeze-back
test; and the FEniCSx comparison appended to `tests/test_reentry_model_fenicsx.py`, run in `fenicsx_env`):

```diff
--- a/tests/test_reentry_model_melting.py
+++ b/tests/test_reentry_model_melting.py
@@ -210,7 +210,11 @@
     element once the surface falls through the feed ramp, and feed and freeze are netted, so neither runs while the
     other does. Mass and the coupled energy balance are exact throughout."""
     import types
-    b = melting_body(layered_mesh, name="AA7075")
+    # the molten cascade (amendment of 2026-10-03) stays off here: this device puts 2.2 MW/m2 on every facet, crater walls
+    # and lee included, of a body 50 K below its melting point, and once the surface may recede through molten material
+    # within the step the whole body melts in the 8 s of heating (measured: consumed at step 16, the books exact
+    # throughout), leaving no film to freeze back; the cascade's own books are tested with the amendment's tests below
+    b = melting_body(layered_mesh, name="AA7075", molten_cascade=False)
     b.solver.set_temperature(800.0)                                    # just below AA7075's melting point
     b.energy0 = b.energy()
     loads = lambda q: types.SimpleNamespace(q_conv=np.full(b.surface.n_patches, q))
@@ -424,3 +428,175 @@
     assert stats["deep_runoff_mass_kg"] == b.deep_runoff_mass and stats["deep_mass_kg"] == pytest.approx(b.m_d.sum())
     assert stats["deep_surfaced_mass_kg"] == b.deep_surfaced_mass and 0.0 <= stats["deep_blob_fraction"] <= 1.0
     assert stats["deep_liquid_kg"] == b.last_melt["deep_liquid"]
+
+
+# ---------------------------------------------------------------------------------------------------------------
+# Amendment of 2026-10-03: the molten cascade, the surface receding through molten material within the step
+
+
+def fully_molten(b):
+    """Per element: fully molten as `molten_depth` has it -- the mean nodal temperature at or above the top of the feed
+    ramp -- which is the test the molten cascade uses."""
+    return b.solver.temperature()[b.mesh.tets].mean(axis=1) >= b.material.T_feed
+
+
+def test_the_molten_cascade_empties_a_fully_molten_column_within_one_step(layered_mesh):
+    """Fact 50: the melt step runs after the conduction, so the material beneath the wall-owning element melts as well,
+    and the element the owner's death exposes was buried a moment earlier -- the feed gate held its melt in place --
+    so it waited, fully molten, for the next step's feed: the surface receded through molten material by one element
+    per macro step. The molten cascade feeds such an element within the step, so that it dies in turn and exposes the
+    next. A pool 4 mm deep under the nose, every node in it above the top of the feed ramp, over a body at 850 K: with
+    the cascade one melt step leaves no newly exposed wall-owning element fully molten -- the surface has receded
+    through the whole pool to material that is not -- every element the cascade fed was fully molten, and its whole
+    mass went to the film; without it the elements the first deaths exposed are still there, fully molten. Mass and
+    enthalpy books exact."""
+    out = {}
+    for cascade in (True, False):
+        b = melting_body(layered_mesh, molten_cascade=cascade)
+        molten_pool(b, depth=4e-3)
+        full, start = fully_molten(b), b.owner_area > 0.0
+        mass0, E0 = b.mass(0.0), b.energy() + b.pending_load.sum()
+        b.melt_step(0.0, 0.5, None)
+        assert b.mass(0.0) == pytest.approx(mass0, rel=1e-12)
+        assert b.energy() + b.pending_load.sum() == pytest.approx(E0, rel=1e-12)
+        out[cascade] = (b, full, start)
+    (on, full, start), (off, _, _) = out[True], out[False]
+    exposed = lambda b: (b.owner_area > 0.0) & b.mesh.active & ~start          # the wall this step's deaths exposed
+    assert (exposed(off) & full).sum() > 0                              # without: the exposed layer waits, fully molten
+    assert exposed(on).any() and not (exposed(on) & full).any()        # with: the surface receded through the pool
+    lm = on.last_melt
+    assert lm["cascade_passes"] >= 2 and not lm["cascade_capped"] and on.cascade_mass == lm["cascade_mass"] > 0.0
+    fed = off.phi - on.phi > 0.0                                        # what the cascade took beyond the step without it
+    assert fed.any() and full[fed].all() and not on.mesh.active[fed].any()     # fully molten elements only, now gone
+    assert on.cascade_mass == pytest.approx(((off.phi - on.phi) * on.element_mass)[fed].sum(), rel=1e-12)
+    assert on.m_f.sum() - off.m_f.sum() == pytest.approx(on.cascade_mass, rel=1e-9)   # ... all of it into the film
+    assert off.cascade_mass == 0.0 and off.last_melt["cascade_passes"] == 0 and not off.last_melt["cascade_capped"]
+    stats = on.melt_stats()
+    assert stats["cascade_passes"] == lm["cascade_passes"] and stats["cascade_mass_kg"] == on.cascade_mass
+
+
+def test_the_molten_cascade_changes_nothing_where_the_exposed_elements_are_not_fully_molten(layered_mesh):
+    """The cascade adds a feed and changes nothing else. At 69.8 km under Girin's closure, with the surface nodes at
+    960 K over a body at 850 K and a millimetre of deep liquid on every windward patch: the wall-owning elements (three
+    surface nodes of four) are molten, so the feed takes three quarters of what they hold and they die; every element
+    their deaths expose has at most two surface nodes, a mean temperature of at most 905 K, and so is not fully molten;
+    the deep runoff moves the deep liquid and the film sprays -- and the whole melt step is bit-identical with the
+    cascade on and off."""
+    runs = []
+    for cascade in (True, False):
+        b = melting_body(layered_mesh, molten_cascade=cascade)
+        sim, a = girin_state(b)
+        boundary = np.zeros(b.mesh.n_nodes, dtype=bool)
+        boundary[np.unique(b.surface.faces)] = True
+        b.solver.set_temperature(np.where(boundary, 960.0, 850.0))
+        owners = b.owner_area > 0.0
+        b.phi[owners] = 0.15                                             # the feed takes 3/4 of it: the owners die
+        b.solver.set_fractions(b.phi)
+        b.m_d = np.where(b.windward, 1e-3 * b.liquid.rho * b.surface.areas, 0.0)
+        b.melt_step(sim.t, 0.5, a)
+        runs.append(b)
+    on, off = runs
+    assert on.last_melt["n_dead"] > 0 and on.sprayed_mass > 0.0 and on.deep_surfaced_mass > 0.0
+    assert on.cascade_mass == 0.0 and on.last_melt["cascade_passes"] == 0 and not on.last_melt["cascade_capped"]
+    for name in ("phi", "m_f", "m_d", "pending_load"):
+        assert np.array_equal(getattr(on, name), getattr(off, name)), name
+    assert np.array_equal(on.mesh.active, off.mesh.active) and on.sprayed_mass == off.sprayed_mass
+    assert on.removed_enthalpy == off.removed_enthalpy and on.runoff_mass == off.runoff_mass
+    assert on.deep_runoff_mass == off.deep_runoff_mass and on.deep_surfaced_mass == off.deep_surfaced_mass
+
+
+def test_the_molten_cascade_feeds_only_the_wall_and_stops_at_material_that_is_not_fully_molten(layered_mesh, monkeypatch):
+    """The cascade respects the feed gate of fact 28(b) by construction: it feeds an element only once a death has made
+    it a wall owner -- in the sheared layer by definition, which is what the gate admits -- and it stops at the first
+    exposed element that is not fully molten, so molten material separated from the wall by material that is not is
+    never fed (the contiguity rule of `molten_depth`). A fully molten skin 2 mm deep, a layer at 880 K beneath it and a
+    fully molten core below that, at 69.8 km under Girin's closure, with the deep runoff on and off: one melt step feeds
+    the skin through the cascade -- every element it fed a wall owner at that moment and fully molten -- and leaves
+    every element of the buried core exactly as it was."""
+    for deep in (True, False):
+        b = melting_body(layered_mesh, deep_runoff=deep)
+        sim, a = girin_state(b)
+        r = np.linalg.norm(b.mesh.points, axis=1)
+        b.solver.set_temperature(np.where(r > 0.048, 960.0, np.where(r > 0.040, 880.0, 960.0)))
+        full = fully_molten(b)
+        core = (r[b.mesh.tets] < 0.040).all(axis=1)                     # every node in the buried molten core
+        calls, feed = [], b._feed_exposed
+
+        def spy(exposed, h_e, h_p, b=b, feed=feed, calls=calls):
+            calls.append(((b.owner_area[exposed] > 0.0).all(), full[exposed].all(), exposed.size))
+            return feed(exposed, h_e, h_p)
+
+        monkeypatch.setattr(b, "_feed_exposed", spy)
+        phi0 = b.phi.copy()
+        b.melt_step(sim.t, 0.5, a)
+        assert calls and all(owner and molten for owner, molten, _ in calls), deep
+        assert sum(n for _, _, n in calls) > 0 and b.cascade_mass > 0.0
+        assert core.sum() > 100 and full[core].all()
+        assert np.array_equal(b.phi[core], phi0[core]) and b.mesh.active[core].all(), deep   # the buried pool is untouched
+        assert not ((b.owner_area > 0.0) & b.mesh.active & full).any()  # and the skin above the 880 K layer is gone
+
+
+def test_with_the_molten_cascade_the_deep_liquid_is_still_never_sprayed(layered_mesh, monkeypatch):
+    """The cascade feeds the film, never the deep account, and it runs after the spray; so with both on, a step in which
+    the deep runoff moves the liquid below the conjugate depth and the cascade carries the surface through a molten
+    pool still sprays only the film. Checked stage by stage on a pool 4 mm deep at 69.8 km with a millimetre of deep
+    liquid already on every windward patch: the spray stage leaves the deep account bit-identical, every element death
+    moves deep liquid without changing its total, and each cascade feed lands in the film alone."""
+    b = melting_body(layered_mesh)
+    sim, a = girin_state(b)
+    molten_pool(b, depth=4e-3)
+    b.m_d = np.where(b.windward, 1e-3 * b.liquid.rho * b.surface.areas, 0.0)
+    seen = {"spray": [], "kill": [], "feed": []}
+    spray_stage, kill, feed = b._film_and_spray, b._kill, b._feed_exposed
+
+    def spy_spray(*args, **kw):
+        before = b.m_d.copy()
+        out = spray_stage(*args, **kw)
+        seen["spray"].append(np.array_equal(before, b.m_d))
+        return out
+
+    def spy_kill(dead):
+        before = float(b.m_d.sum())
+        kill(dead)
+        seen["kill"].append((before, float(b.m_d.sum())))
+
+    def spy_feed(exposed, h_e, h_p):
+        deep0, film0 = b.m_d.copy(), float(b.m_f.sum())
+        out = feed(exposed, h_e, h_p)
+        seen["feed"].append((np.array_equal(deep0, b.m_d), float(b.m_f.sum()) - film0, out))
+        return out
+
+    monkeypatch.setattr(b, "_film_and_spray", spy_spray)
+    monkeypatch.setattr(b, "_kill", spy_kill)
+    monkeypatch.setattr(b, "_feed_exposed", spy_feed)
+    b.melt_step(sim.t, 0.5, a)
+    assert seen["spray"] == [True] and b.sprayed_mass > 0.0             # the spray never touched the deep account
+    assert seen["kill"] and all(after == pytest.approx(before, rel=1e-12) for before, after in seen["kill"])
+    assert seen["feed"] and all(same and gain == pytest.approx(out, rel=1e-9) for same, gain, out in seen["feed"])
+    assert b.deep_runoff_mass > 0.0 and b.m_d.sum() > 0.0 and b.cascade_mass > 0.0
+
+
+def test_full_steps_with_the_molten_cascade_keep_the_books_exact(layered_mesh):
+    """Six coupled macro steps from 69.8 km with a pool 6 mm deep under the nose and the physics loads: the cascade acts
+    (deaths expose fully molten elements and the surface recedes through them within the step), the deep runoff and
+    the spray act with it, and the energy balance and the mass stay exact by the same measures as every other melting
+    test; the solver carries the film and the deep liquid, and the history's two new columns report the cascade."""
+    b = melting_body(layered_mesh)
+    sim, a = girin_state(b)
+    molten_pool(b, depth=6e-3)
+    b.energy0 = b.energy()
+    model = heating.PhysicsHeating()
+    passes = []
+    for _ in range(6):
+        sim.advance(0.5)
+        a = sim.aero_state(sim.t, sim.y[:3], sim.y[3:])
+        loads = model.evaluate(a, b.theta, b.surface_temperature(), b.nose_radius(), T_mean=b.mean_temperature())
+        b.advance(sim.t, 0.5, loads, state=a)
+        passes.append(b.last_melt["cascade_passes"])
+        assert abs(b.energy_balance_residual()) < 1e-7
+        assert b.mass(0.0) + b.removed_mass == pytest.approx(b.mass0, rel=1e-12)
+        assert b.solver.film_mass.sum() == pytest.approx(b.m_f.sum() + b.m_d.sum())
+    assert b.cascade_mass > 0.0 and max(passes) >= 2 and b.cascade_capped_steps == 0
+    assert b.sprayed_mass > 0.0 and b.deep_runoff_mass > 0.0
+    stats = b.melt_stats()
+    assert stats["cascade_mass_kg"] == b.cascade_mass and stats["cascade_passes"] == passes[-1]
```

```diff
--- a/tests/test_reentry_model_fenicsx.py
+++ b/tests/test_reentry_model_fenicsx.py
@@ -142,3 +142,44 @@
         s.set_temperature(450.0)
         assert s.facet_temperature() == pytest.approx(np.full(len(s.areas), 450.0)), name
         assert s.facet_temperature() == pytest.approx(s.facet_temperature(s.temperature())), name
+
+
+def test_the_molten_cascade_matches_the_skfem_backend(coarse_sphere_mesh):
+    """The molten cascade (amendment of 2026-10-03) in both backends: a pool 6 mm deep under the nose at 69.8 km and six
+    coupled steps with the physics loads -- deaths expose fully molten elements and the cascade feeds them within the
+    step, while the deep runoff and the spray act -- give the same cascade passes, the same active set, the same mass,
+    sprayed mass, cascade mass and deep account, and the same temperatures, with both energy balances exact. numpy's
+    generator is seeded for each backend because pyamg draws its starting vectors from it (plan fact 52). Measured on
+    2026-10-03: passes 4, 3, 16, 4, 5, 4 in both; mass 1.3e-7, sprayed 6.3e-7, cascade mass 2.4e-7, deep account 3.7e-6,
+    temperatures 0.11 K, balances 4.7e-10 and 2.4e-11 -- the stiff deep transport carries the backends' last-bit
+    differences further (plan fact 49); with the deep runoff off they agree to 1e-10 and 4e-6 K."""
+    pytest.importorskip("cantera")
+    from reentry_model import body, heating, mesh
+    from test_reentry_model_coupled import MASS_100MM, simulator
+    out = {}
+    for name in ("skfem", "fenicsx"):
+        np.random.seed(12345)
+        m = mesh.VolumeMesh(coarse_sphere_mesh.points, coarse_sphere_mesh.tets, dict(coarse_sphere_mesh.params))
+        b = body.MeltingBody(m, material.Material.from_drama_json("AA7075_range"), thermal.thermal_solver(name), MASS_100MM)
+        sim = simulator(b, t_max=90.0)
+        sim.advance(50.0)
+        p = b.mesh.points
+        r = np.linalg.norm(p, axis=1)
+        b.solver.set_temperature(np.where((r > 0.05 - 6e-3) & (p[:, 0] > 0.5 * r), 960.0, 850.0))
+        b.energy0 = b.energy()
+        model = heating.PhysicsHeating()
+        passes = []
+        for _ in range(6):
+            sim.advance(0.5)
+            a = sim.aero_state(sim.t, sim.y[:3], sim.y[3:])
+            b.advance(sim.t, 0.5, model.evaluate(a, b.theta, b.surface_temperature(), b.nose_radius(), T_mean=b.mean_temperature()),
+                      state=a)
+            passes.append(b.last_melt["cascade_passes"])
+        out[name] = (b, passes)
+    (s, ps), (f, pf) = out["skfem"], out["fenicsx"]
+    assert ps == pf and max(pf) >= 2 and f.cascade_mass > 0.0 and f.deep_runoff_mass > 0.0
+    assert np.array_equal(s.mesh.active, f.mesh.active)
+    assert f.cascade_mass == pytest.approx(s.cascade_mass, rel=1e-6) and f.mass(0.0) == pytest.approx(s.mass(0.0), rel=1e-6)
+    assert f.sprayed_mass == pytest.approx(s.sprayed_mass, rel=1e-5) and f.m_d.sum() == pytest.approx(s.m_d.sum(), rel=1e-4)
+    assert np.abs(f.solver.temperature() - s.solver.temperature()).max() < 0.5
+    assert abs(s.energy_balance_residual()) < 1e-8 and abs(f.energy_balance_residual()) < 1e-8
```

### Measured (2026-10-02/03)

- **Base verified.** A fresh copy of `prototype/proto3/` (unchanged since the deep-runoff amendment: all 114 files match
  that amendment's final manifest) with the nine diff blocks of 2026-10-02 applied is byte-identical to the tested copy
  of that amendment across the whole tree, and its unit tier gives that amendment's 223 passed, 1 skipped, 2 failed and
  3 errors.
- **Off is the old model, bit for bit** (fact 55): with `--molten-cascade off`, the 100 mm flight to 120 s (seeded)
  reproduces the deep-runoff amendment's run in all 241 rows of all 84 shared history columns, all 167 473 source rows
  and every result field but the run time.
- **Tests**: the seven above; the whole unit tier gives 229 passed, 1 skipped, 2 failed and 3 errors — the base's 223
  passed plus the six new tests (five melting tests and one bad-argument case), the failures and errors being the five
  known ones from the missing melting SESAM references (Task 11); `tests/test_reentry_model_fenicsx.py` 10 passed in
  `fenicsx_env`.
- **What the backlog is made of** (fact 54): the molten chains below the wall-owning elements run along the melt
  front, through elements molten by their mean temperature but never at all four nodes, because the feed's own energy
  debit keeps every node an element shares with an element being fed on the ramp. The cascade's first, strict form —
  every node above T_feed — fed 32 g in 120 s and changed nothing, so the model's own definition replaced it.
- **The default step** (fact 57): the cascade feeds 245 g over the flight (a median of 3 passes in 146 of 240 steps,
  never capped) and removes the backlog — liquid held below the conjugate depth a median 1.23 g instead of 4.75 g, an
  end-of-step backlog of 0.41 g — but the molten layer the branch test reads falls only from 1.26 to 0.97 mm, the
  thick-branch share from 12.4 % to 11.9 %, and the droplet population does not move: 18.1 against 17.2 million
  droplets (58.5 million at 0.125 s without the cascade), 180 µm median radius by number either way (76 µm at 0.125 s),
  thin-branch share 10.5 % against 9.4 % (38 %), Rayleigh–Taylor release 51.9 against 53.0 g (8.5 g), re-solidified
  6.2 g either way (48.4 g). Sprayed mass +3.6 %, mass at 120 s −8 %, balance exact (−6.6e-11).
- **The next cause** (fact 58): one molten wall-owning element counts as 0.67–1.16 mm of liquid in `molten_depth` on
  the production mesh — above the conjugate depth on every patch — so the branch test asks only whether the owner's mean
  temperature is above T_feed when the spray runs, which the length of the conduction step decides.
- **The deep runoff after the fix** (fact 59): 91 g mobilised instead of 285 g, a deep account peaking at 12 g instead
  of 117 g, and an effect on the sprayed mass of +0.5 % instead of +7.1 %.
- **The 50 mm flight** (fact 60): the cascade feeds 95.6 g; demise 8.0 s earlier (199.0 s, 70.7 km), sprayed mass
  +1.0 %, droplets −32 %, median radius by number +17 %, balance exact.
- **Both backends** (fact 56): identical passes and active sets, agreement to 1e-7 (deep runoff on) and 1e-10 (off).
- **Not measured** (fact 60): the 0.25 and 0.125 s flights with the cascade, the whole 100 mm flight;
  the machine ran on battery at 4–7 % and slept through most of the hours they need.

### Not done here — for Asha to decide (fact 61)

First, the branch test's liquid depth: the cascade removed the backlog, but the droplet population still depends on
the step because `molten_depth` counts whole elements read at the end of the conduction; the recommendation is a depth
by mass, as an amendment of its own with the time-step study repeated, and until then no droplet-population number from
the default step. Then: keep the cascade on and re-measure it at the smaller steps; keep the deep runoff on, its effect
now about half a per cent; measure again with fact 44's `PHI_DEATH = 0.50`; whether the cascade should run before the
deep stage; and seeding.

## Amendment of 2026-10-05 (runoff flux) — what the melting body sees of the film's new flux (no change to `body.py`)

> Part of the runoff-flux amendment: sub-plan 06's amendment of this date holds the change to `film.lubrication`, its
> design and the alternatives; facts 69–77 in `00-shared-context.md` carry the measurements. Measured in the same
> throwaway copy as sub-plan 06's amendment (the amendments of 2026-10-02, 2026-10-03 and 2026-10-05 (seeding) applied
> to `prototype/proto3/`, verified first to be reproduced byte for byte by their 20 diff blocks). The copy keeps
> `PHI_DEATH = 0.05`, as every amendment since 2026-10-02 has: fact 44's 0.50 exists only as a plan.

**No change to this task's code, and why none is needed.** The melting body calls `film.lubrication` twice in each melt
step, both times in `_film_and_spray`. The first call is inside the film's runoff transport (`q_of_b`, with the layer
`b_layer = bb + molten + deep` deciding the branch); that is the flux sub-plan 06 changes, and it is where the change
acts. The second call, after the transport, uses only the surface velocity V_s (handed to the spray, whose Weber number
reads it) and the branch flag (for `thick_branch_fraction`); neither changes. The deep stage `_deep_runoff` moves the
deep liquid with `film.deep_flux`, which is unchanged, and its emptying rate stays bounded as the deep liquid vanishes
(sub-plan 06's amendment shows why; fact 74 measures it in flight). One property of the existing code makes the
amendment of 2026-10-02's column identity hold wherever deep liquid sits: the deep stage tops the film up from the deep
account to the conjugate depth before the film moves, so a patch that carries deep liquid starts the film transport with
b ≥ δ_m, where the film's flux is exactly what it was. The change acts only on windward patches under Girin's closure
(the only place a conjugate depth exists) whose liquid layer is deeper than δ_m while their film is thinner, so a flight
that never has his closure is bit-identical with and without it — measured on the whole 50 mm flight — and the 100 mm
flight is bit-identical until its first step under his closure, at 49.5 s (fact 73).

**One melting test's tolerance moves, and why that is legitimate.**
`test_full_steps_with_the_molten_cascade_keep_the_books_exact` (this task's amendment of 2026-10-03) runs six coupled
steps over a 6 mm molten pool and holds the mass to 1e-12 of the body. With the new flux it failed by 1.6e-12. The drift
is the film transport's: fact 48 measured that its direct solve conserves the total only to its conditioning, and the
deep stage was then made exact by scaling its arrivals to its departures, while the film transport was left as it was
because correcting it would change every existing run in the last bits. This pool drives up to a quarter of a kilogram
of film per step through pairs of patches that drain into each other, with emptying rates that, multiplied by the
sub-step, reach about 10⁷ in both directions, and the cancellation in those pairs sets the drift. Measured over ten
states of numpy's generator (seeds 12345 and 1 to 9), the largest drift over the six steps, as a share of the body's
mass, is 6.8·10⁻¹³ to 2.8·10⁻¹² before the amendment (one state of ten already failing the 1e-12 check) and 9·10⁻¹⁴ to
4.4·10⁻¹² after it (six of ten failing), the largest single film transport changing the total by 4.0·10⁻¹² kg before and
6.4·10⁻¹² kg after it, 1.7·10⁻¹¹ and 2.6·10⁻¹¹ of the mass it moved. The check was therefore at the edge of the
transport's round-off before the amendment and passed in the seeded unit tier by chance; it now holds the mass to 1e-11
of the body, more than twice the largest drift measured. The same sweep on the deep-runoff test's 3 mm pool, which keeps
its 1e-12, drifts by at most 4.8·10⁻¹⁴ before and after. Correcting the film transport's books exactly is recorded as an
option for Asha (fact 77), not taken here, because it would move every run, the 50 mm flight included, in the last bits.

```diff
--- a/tests/test_reentry_model_melting.py
+++ b/tests/test_reentry_model_melting.py
@@ -580,7 +580,11 @@
     """Six coupled macro steps from 69.8 km with a pool 6 mm deep under the nose and the physics loads: the cascade acts
     (deaths expose fully molten elements and the surface recedes through them within the step), the deep runoff and
     the spray act with it, and the energy balance and the mass stay exact by the same measures as every other melting
-    test; the solver carries the film and the deep liquid, and the history's two new columns report the cascade."""
+    test; the solver carries the film and the deep liquid, and the history's two new columns report the cascade. The
+    mass is exact to the film transport's round-off, which its direct solve conserves only to its conditioning (plan
+    fact 48): this 6 mm pool moves up to a quarter of a kilogram of film per step through strongly coupled patch pairs,
+    and with the runoff flux of 2026-10-05 the six steps drift by up to 4.4e-12 of the body (measured over ten states
+    of numpy's generator; up to 2.8e-12 before that amendment), so the mass is held to 1e-11 here."""
     b = melting_body(layered_mesh)
     sim, a = girin_state(b)
     molten_pool(b, depth=6e-3)
@@ -594,7 +598,7 @@
         b.advance(sim.t, 0.5, loads, state=a)
         passes.append(b.last_melt["cascade_passes"])
         assert abs(b.energy_balance_residual()) < 1e-7
-        assert b.mass(0.0) + b.removed_mass == pytest.approx(b.mass0, rel=1e-12)
+        assert b.mass(0.0) + b.removed_mass == pytest.approx(b.mass0, rel=1e-11)
         assert b.solver.film_mass.sum() == pytest.approx(b.m_f.sum() + b.m_d.sum())
     assert b.cascade_mass > 0.0 and max(passes) >= 2 and b.cascade_capped_steps == 0
     assert b.sprayed_mass > 0.0 and b.deep_runoff_mass > 0.0
```

**The two backends, with the change active.** A new FEniCSx test puts the cascade test's 6 mm pool at 69.8 km through
ten steps of 0.05 s, where every film transport has thick patches carrying films thinner than δ_m, and asserts that the
change is exercised, that the largest emptying rate stays below 10⁹ per second (it reached 3.9·10²⁵ per second in this
setting before the change and is 2.7·10⁷ with it, set by the G b³ term of the deepest piles, not by a vanishing film)
and that the backends agree. The change also revises the reason facts 49 and 56 gave for the backends' looser agreement
when deep liquid is present: in the cascade test's own setting (six steps of 0.5 s) the passes are the same with and
without the change (4, 3, 16, 4, 5, 4), and the agreement improves from 1.3·10⁻⁷ to 7.6·10⁻¹⁰ in mass, 6.3·10⁻⁷ to
3.6·10⁻⁹ in sprayed mass, 2.4·10⁻⁷ to 1.0·10⁻¹⁰ in cascade mass, 3.7·10⁻⁶ to 6.9·10⁻¹¹ in the deep account and 0.11 K to
2.5·10⁻⁵ K in temperature. The deep transport is unchanged, so what carried the backends' last-bit differences further
was the film's runaway emptying rates on thick patches; the cascade test's docstring says so.

```diff
--- a/tests/test_reentry_model_fenicsx.py
+++ b/tests/test_reentry_model_fenicsx.py
@@ -152,7 +152,10 @@
     generator is seeded for each backend because pyamg draws its starting vectors from it (plan fact 52). Measured on
     2026-10-03: passes 4, 3, 16, 4, 5, 4 in both; mass 1.3e-7, sprayed 6.3e-7, cascade mass 2.4e-7, deep account 3.7e-6,
     temperatures 0.11 K, balances 4.7e-10 and 2.4e-11 -- the stiff deep transport carries the backends' last-bit
-    differences further (plan fact 49); with the deep runoff off they agree to 1e-10 and 4e-6 K."""
+    differences further (plan fact 49); with the deep runoff off they agree to 1e-10 and 4e-6 K. With the runoff flux
+    of 2026-10-05 the passes are the same and the backends agree to 7.6e-10 in mass, 3.6e-9 in sprayed mass, 1.0e-10 in
+    cascade mass, 6.9e-11 in the deep account and 2.5e-5 K: what carried the last-bit differences was the film's
+    runaway emptying rates on thick patches, not the deep transport (plan fact 72)."""
     pytest.importorskip("cantera")
     from reentry_model import body, heating, mesh
     from test_reentry_model_coupled import MASS_100MM, simulator
@@ -183,3 +186,60 @@
     assert f.sprayed_mass == pytest.approx(s.sprayed_mass, rel=1e-5) and f.m_d.sum() == pytest.approx(s.m_d.sum(), rel=1e-4)
     assert np.abs(f.solver.temperature() - s.solver.temperature()).max() < 0.5
     assert abs(s.energy_balance_residual()) < 1e-8 and abs(f.energy_balance_residual()) < 1e-8
+
+
+def test_the_thick_film_flux_matches_the_skfem_backend(coarse_sphere_mesh, monkeypatch):
+    """The runoff flux on a thick patch (amendment of 2026-10-05) in both backends: the 6 mm pool at 69.8 km of the
+    cascade test above, at ten steps of 0.05 s, where every film transport has thick patches carrying films thinner
+    than delta_m -- the case the amendment changes. Both backends give the same active set, mass, sprayed mass, runoff,
+    film, deep account and temperatures, with both energy balances exact, and the largest emptying rate of any edge
+    stays bounded: before the amendment it grew from sub-step to sub-step and reached 3.9e25 per second in this setting
+    (measured 2026-10-05); with it 2.7e7, set by the deepest piles' G b^3 / (3 mu) term, not by a vanishing film.
+    Measured: the same 7 983 active elements, mass 5.8e-10, sprayed mass 5.9e-9, runoff 7.5e-10, film 1.7e-8, deep
+    account 2.8e-10, temperatures 2.3e-4 K, balances 5.1e-9 and 1.1e-11."""
+    pytest.importorskip("cantera")
+    from reentry_model import body, film, heating, mesh
+    from test_reentry_model_coupled import MASS_100MM, simulator
+    seen = {"thin_on_thick": 0, "c_max": 0.0}
+    lub, coeff = film.lubrication, film.Runoff.edge_coefficients
+
+    def lubrication(tau, G, b, delta_m, mu_l, b_layer=None):
+        out = lub(tau, G, b, delta_m, mu_l, b_layer)
+        bb = np.asarray(b, dtype=float)
+        seen["thin_on_thick"] += int((out[3] & (bb > 0.0) & (bb < np.asarray(delta_m, dtype=float))).sum())
+        return out
+
+    def edge_coefficients(self, q, b, t_hat, areas):
+        c_ij, c_ji = coeff(self, q, b, t_hat, areas)
+        if len(c_ij):
+            seen["c_max"] = max(seen["c_max"], float(c_ij.max()), float(c_ji.max()))
+        return c_ij, c_ji
+
+    monkeypatch.setattr(film, "lubrication", lubrication)
+    monkeypatch.setattr(film.Runoff, "edge_coefficients", edge_coefficients)
+    out = {}
+    for name in ("skfem", "fenicsx"):
+        np.random.seed(12345)
+        m = mesh.VolumeMesh(coarse_sphere_mesh.points, coarse_sphere_mesh.tets, dict(coarse_sphere_mesh.params))
+        b = body.MeltingBody(m, material.Material.from_drama_json("AA7075_range"), thermal.thermal_solver(name), MASS_100MM)
+        sim = simulator(b, t_max=90.0)
+        sim.advance(50.0)
+        p = b.mesh.points
+        r = np.linalg.norm(p, axis=1)
+        b.solver.set_temperature(np.where((r > 0.05 - 6e-3) & (p[:, 0] > 0.5 * r), 960.0, 850.0))
+        b.energy0 = b.energy()
+        model = heating.PhysicsHeating()
+        for _ in range(10):
+            sim.advance(0.05)
+            a = sim.aero_state(sim.t, sim.y[:3], sim.y[3:])
+            b.advance(sim.t, 0.05, model.evaluate(a, b.theta, b.surface_temperature(), b.nose_radius(), T_mean=b.mean_temperature()),
+                      state=a)
+        out[name] = b
+    s, f = out["skfem"], out["fenicsx"]
+    assert seen["thin_on_thick"] > 0 and seen["c_max"] < 1e9 and f.runoff_mass > 0.0 and f.sprayed_mass > 0.0
+    assert np.array_equal(s.mesh.active, f.mesh.active)
+    assert f.mass(0.0) == pytest.approx(s.mass(0.0), rel=1e-6) and f.sprayed_mass == pytest.approx(s.sprayed_mass, rel=1e-5)
+    assert f.runoff_mass == pytest.approx(s.runoff_mass, rel=1e-5) and f.m_f.sum() == pytest.approx(s.m_f.sum(), rel=1e-4)
+    assert f.m_d.sum() == pytest.approx(s.m_d.sum(), rel=1e-4)
+    assert np.abs(f.solver.temperature() - s.solver.temperature()).max() < 0.5
+    assert abs(s.energy_balance_residual()) < 1e-8 and abs(f.energy_balance_residual()) < 1e-8
```

### Measured (2026-10-05)

- **Base verified.** A fresh copy of `prototype/proto3/` (unchanged since 2026-10-03: all 114 files match that
  amendment's final manifest) with the 20 diff blocks of 2026-10-02, 2026-10-03 and 2026-10-05 (seeding) applied is
  byte-identical to the tested copy of the seeding amendment across the whole tree.
- **Tests.** The whole unit tier gives 236 passed, 1 skipped, 2 failed and 3 errors — the base's 234 passed plus the two
  new film tests, the failures and errors being the five known ones from the missing melting SESAM references (Task 11);
  before the tolerance change above it gave 235 passed and one more failure, that test.
  `tests/test_reentry_model_fenicsx.py` in `fenicsx_env`: 11 passed (the ten of before and the new one).
- **The flights** (facts 73–75): the 50 mm flight bit-identical; the 100 mm flight at the default step bit-identical to
  49.5 s, its masses then moving by about the run-to-run scatter of fact 65, its droplet count and median radius by
  number within it, its front-surface Rayleigh–Taylor release (−16 %) outside it, and the film depth at release down by
  14 %; the switched-step flight at 0.0125 s, which crashed, now runs to 120 s with its emptying rates bounded, and so
  does 0.00625 s.

### Not done here — for Asha to decide (fact 77)

First the time step (fact 77 (a)): the switched-step series of fact 75 leaves the droplet population unconverged at
0.00625 s, each branch's droplets the same at every step and the branch split not, and the recommendation is the branch
test on a liquid depth by mass, with the series repeated. Then: whether to make the film transport's books exact as the
deep stage's are; the exit code a model failure is reported with; the re-solidification counter, which counts
freeze-and-re-melt cycles; and whether the step switch should become a model option.

## Amendment of 2026-10-06 — the thin branch needs a rigid substrate: the non-rigid depth and the regime layer

> Asha's request and decisions of 2026-10-06 (facts 78–87 in `00-shared-context.md` carry the measurements and the
> decisions in full). Measured in a throwaway copy of `prototype/proto3/` with the 24 diff blocks of the amendments of
> 2026-10-02, 2026-10-03 and both of 2026-10-05 applied in date order — every one exactly — and the copy verified before
> any change to reproduce fact 73's flights to the last digit. The code below is the tested code, as a diff against that
> copy. Sub-plan 07's amendment of this date holds the spray's half of the change, 06's what the film sees of it, 10's
> the history columns and the frame field, 13's the flag and 02's `Material.T_rigid`.

**Why.** The spray's regime test (fact 28(a), sub-plan 07) sends a patch to Girin's thin branch — his dominant ablation,
Girin & Kopyt's 1994 mode — when the liquid under it is shallower than the conjugate depth δ_m, on the premise that "the
rigid core still stabilises the disturbances". On an alloy with a 158 K melting range the material under the liquid is
mush, and mush more than half liquid is a slurry that flows (Li et al. 2014), not a rigid core: a film on it is deep melt
and belongs on the thick branch. Coherent mush, less than half liquid, carries load (Chen et al. 2016) and is a rigid
substrate. Fact 54 showed what lies under the film on the flights: the element directly beneath it always has a node in
the mushy range, at a median 897 K — 93 % liquid on the linear law.

**The change, in `body.py`.**

1. `MeltSettings.rigid_substrate: bool = True`, set by `--rigid-substrate on|off` (sub-plan 13).
2. `MeltingBody.nonrigid_depth(T=None, max_crossings=NONRIGID_MAX_CROSSINGS)` — the depth of non-rigid material under each
   patch (fact 79): a line from the patch centre along the inward normal through the elements it actually crosses, the
   P1 field's crossing of `Material.T_rigid` found exactly inside an element, stopping at the first rigid point, the
   active boundary or 64 elements, each element's stretch weighed by φ_e. `NONRIGID_MAX_CROSSINGS = 64`.
3. In `_film_and_spray`, with the switch on: the regime layer `b + nonrigid + deep` replaces `b + molten + deep` in both
   calls of `film.lubrication` (the runoff's flux and the spray's surface velocity and branch flag) and is passed to
   `spray.evaluate` as `regime_layer`; `on_slurry = nonrigid > 0` is passed as well, and the film mass on the windward
   patches it holds off Girin's closure is recorded. `layer = b + molten + deep` still goes to the spray as `b_layer`
   (the Rayleigh–Taylor criteria) and `molten` still defines the wave-fits region.
4. Four entries in `last_melt` and `melt_stats()`: `nonrigid_depth_mean_mm` (mean over the wet windward patches),
   `slurry_thick_fraction` (share of the wet windward patches under Girin's closure made thick only by the slurry),
   `rigid_thin_fraction` (share made thin by the non-rigid depth where the liquid layer said thick — added during the
   measurement, to report fact 58's count being removed) and `slurry_held_mass_kg` (the film held off the closure); all
   nan with the switch off. `last_nonrigid` keeps the step's depth for the surface frame (sub-plan 10).

With the switch off nothing of this runs and every flight is bit-identical to the copy before the change (fact 81).

**What does not change.** The feed, the feed gate of fact 28(b), the deep runoff (it still drains `_molten_chain`, the
fully molten chain), the molten cascade, freeze-back, deaths and every enthalpy booking: the change moves no mass and no
heat itself, it only decides which branch a patch takes. The energy balance stays exact on every flight (facts 82–84).

**Alternatives considered, and why they were not taken.**

- *Any slurry beneath the film triggers the thick branch.* Offered and declined (decision (2)): a slurry skin shallower
  than δ_m over coherent mush leaves rigid material inside the shear's reach, which is Girin's dominant-ablation case.
- *Judge the substrate element by element* (the owner's mean liquid fraction, or `molten_depth`'s march with a T_rigid
  test). Whole elements are 0.67–1.16 mm deep under the patches, all deeper than δ_m (fact 58), so every patch with a
  slurry owner would be thick however thin its slurry: the element-size artefact of fact 58 again. The line through the
  P1 field reads the depth inside the element.
- *Keep the liquid layer for `film.lubrication` and give only the spray Girin's surface velocity.* A patch would carry two
  surface velocities, the film flowing as over a rigid wall while spraying as deep melt; the regime layer is read by both
  instead, and on a patch made thick by slurry the film moves as the top b of the conjugate layer, as it has over fully
  molten elements since 2026-10-05 (sub-plan 06's amendment of this date).
- *Count the slurry in the Rayleigh–Taylor criteria.* They describe a liquid pool, the slurry is 700 to 4 000 times more
  viscous than the liquid (0.87–5.6 Pa s against 1.3 mPa s), and counting it would widen the applied front-surface mode
  over the nose; not asked for.
- *Keep the thin mode off Girin's closure.* Offered and declined (decision (3)). Facts 83 and 85 measure what the strict
  rule does; fact 87 (a) holds the follow-up.

**Two tests changed and one added.** `test_film_spraying_death_and_balances` starts the body at 880 K, which is 82 %
liquid on the linear law, so its held film piles up and the front-surface Rayleigh–Taylor mode sheds droplets above the
size histogram's 10 mm top edge (fact 85); its histogram check now compares the bins with the released droplets inside
their range and requires the rest to lie above it and within R/4. The 2026-10-05 two-backend test
(`test_the_thick_film_flux_matches_the_skfem_backend`) runs with the rule off: its interior at 850 K makes the whole body
slurry, and the backends then part at the spray floor (fact 86); it stays the record of the runoff-flux agreement, and
`test_the_rigid_substrate_matches_the_skfem_backend` checks the backends with the rule on, the interior at 820 K.

In `reentry_model/body.py`:

```diff
--- a/reentry_model/body.py
+++ b/reentry_model/body.py
@@ -146,6 +146,7 @@
                                  # 2026-09-20: owners at 1e-3 made surface nodes swing by hundreds of K between iterates)
 MAX_CASCADE_PASSES = 32          # safety cap on the molten cascade's passes in one macro step; what it leaves waits for the
                                  # next step, as every exposed element did before the cascade (amendment of 2026-10-03)
+NONRIGID_MAX_CROSSINGS = 64      # elements the non-rigid depth's line may cross before it stops (amendment of 2026-10-06)
 REMOVAL_NAMES = ("girin", "instant")
 
 
@@ -162,6 +163,9 @@
     molten_cascade: bool = True       # a death that exposes a fully molten element (mean T >= T_feed) feeds it within the
                                       # step, so the surface recedes through molten material by more than one element per
                                       # macro step (amendment of 2026-10-03)
+    rigid_substrate: bool = True      # the thin branch needs a rigid substrate: the regime test reads the film plus the slurry
+                                      # (above 50 % liquid) beneath it, and off Girin's closure a film on slurry does not
+                                      # spray (amendment of 2026-10-06)
 
     def __post_init__(self):
         if self.removal not in REMOVAL_NAMES:
@@ -258,6 +262,7 @@
         self.hist_n, self.hist_m = np.zeros(spray_mod.N_BINS), np.zeros(spray_mod.N_BINS)
         self.last_flow = self.last_spray = None
         self.last_delta_m = None                                           # the step's conjugate depth per patch [m]
+        self.last_nonrigid = None                                          # the step's non-rigid depth per patch [m]
         self.last_face_ids = None                                          # face ids of the surface they were evaluated on
         self.last_melt = {"n_dead": 0, "runoff_substeps": 0, "released_mass": 0.0, "n_released": 0.0, "feed_mass": 0.0}
         self.melt_onset = self.spray_onset = None
@@ -724,7 +729,7 @@
         liq, s, mat = self.liquid, self.settings, self.material
         from . import spray as spray_mod
         if state is None or self.m_f.sum() <= 0.0 or flow is None:
-            self.last_flow = self.last_spray = None
+            self.last_flow = self.last_spray = self.last_nonrigid = None
             return 0.0
         areas = self.surface.areas
         # the depth of liquid under each patch: its film plus the contiguous molten material beneath it. The branch
@@ -732,10 +737,16 @@
         # and the release stay on the film, which is the mass that can actually move (decided 2026-09-22).
         molten = self.molten_depth()
         deep = self.m_d / (liq.rho * areas)          # the deep liquid the runoff delivered: part of the layer, never film
+        # The thin branch needs a rigid substrate (amendment of 2026-10-06): the regime test -- lubrication's branch and
+        # the spray's -- reads what lies under the film down to rigid material, the slurry above 50 % liquid included,
+        # instead of the fully molten material alone. The Rayleigh-Taylor criteria and the molten region the wave-fits
+        # test measures against keep the liquid.
+        nonrigid = self.nonrigid_depth() if s.rigid_substrate else None
+        under = molten if nonrigid is None else nonrigid
         # (ii) lubrication and runoff
         n_sub, moved = 0, 0.0
         if self.runoff is not None:
-            q_of_b = lambda bb: self._film_mod.lubrication(flow.tau, flow.G, bb, delta_m, liq.mu, b_layer=bb + molten + deep)[1]
+            q_of_b = lambda bb: self._film_mod.lubrication(flow.tau, flow.G, bb, delta_m, liq.mu, b_layer=bb + under + deep)[1]
             before = self.m_f
             self.m_f, n_sub, moved = self.runoff.transport(self.m_f, q_of_b, self.t_hat, liq.rho, areas, dt)
             self.runoff_mass += moved                                              # mass that arrived on another patch
@@ -746,7 +757,8 @@
             self._spread_to_patches(-float((d * h_p).sum()), np.maximum(d, 0.0))
         b = self.m_f / (liq.rho * areas)
         layer = b + molten + deep
-        v_s, q, _, thick = self._film_mod.lubrication(flow.tau, flow.G, b, delta_m, liq.mu, b_layer=layer)
+        regime_layer = layer if nonrigid is None else b + nonrigid + deep
+        v_s, q, _, thick = self._film_mod.lubrication(flow.tau, flow.G, b, delta_m, liq.mu, b_layer=regime_layer)
         # the molten surface: the area over which a wave could form at all, either wetted by the film or molten in its
         # own right. Its contiguous extent is what every mode's wavelength is measured against (2026-09-24).
         wetted = (b >= spray_mod.B_MIN) | (molten > 0.0)
@@ -754,9 +766,13 @@
         # Girin & Kopyt's W sin(Theta): the deceleration normal to the film, which on this body is W cos(phi). It is the
         # whole deceleration at the stagnation point and vanishes at the equator, where W lies in the surface.
         w_n = flow.deceleration * np.maximum(np.cos(self.theta), 0.0)
-        # (iii) spraying
+        # (iii) spraying. Off Girin's closure a film whose base is slurry takes no shear mode; what it holds is recorded
+        on_slurry = None if nonrigid is None else nonrigid > 0.0
+        held = None if nonrigid is None else (self.windward & (b >= spray_mod.B_MIN) & ~np.isfinite(delta_m) & on_slurry)
+        held_mass = float("nan") if held is None else float(self.m_f[held].sum())
         res = self.spray.evaluate(flow, state, b, delta_m, v_s, self.windward, dt, areas, self.m_f,
-                                  self.transverse_radius, b_layer=layer, extent=extent, deceleration_n=w_n)
+                                  self.transverse_radius, b_layer=layer, extent=extent, deceleration_n=w_n,
+                                  regime_layer=None if nonrigid is None else regime_layer, on_slurry=on_slurry)
         released = float(res.dm.sum())
         if released > 0.0:
             if self.spray_onset is None:
@@ -772,8 +788,21 @@
             self.removed_mass += released
             self.n_released += float(res.dn.sum())
         self.m_f[self.m_f < 1e-30] = 0.0                       # no denormal films (they made 0/0 coefficients in the runoff)
-        self.last_flow, self.last_spray, self.last_b = flow, res, b
+        self.last_flow, self.last_spray, self.last_b, self.last_nonrigid = flow, res, b, nonrigid
         self.last_face_ids = self.surface.face_ids.copy()          # deaths later in this step rebuild the surface
+        # what the rigid-substrate test did this step, on the wet windward patches (thick_branch_fraction's set); the two
+        # fractions are over those under Girin's closure, where the depth decides: made thick only by the slurry, and made
+        # thin by the line's depth where the whole-element liquid layer said thick (fact 58's count, removed)
+        wet_w = self.windward & (self.m_f > 0.0)
+        girin_wet = wet_w & np.isfinite(delta_m)
+        if nonrigid is None or not girin_wet.any():
+            slurry_thick = rigid_thin = float("nan")
+        else:
+            with np.errstate(invalid="ignore"):
+                liquid_thick = layer > delta_m
+            slurry_thick = float((thick & ~liquid_thick)[girin_wet].mean())
+            rigid_thin = float((liquid_thick & ~thick)[girin_wet].mean())
+        nonrigid_mean = float(nonrigid[wet_w].mean() * 1e3) if nonrigid is not None and wet_w.any() else float("nan")
         # The Rayleigh-Taylor criterion is reported, never applied, so report it usefully: a body-level "any patch"
         # boolean says nothing about how much melt is involved or whether the unstable wave even fits on the nose.
         # `rt_mass_fraction` is the share of the film on unstable patches; `rt_wavelength_over_nose` is the shortest
@@ -812,7 +841,9 @@
             "molten_depth_mean_mm": float((self.m_f * molten).sum() / self.m_f.sum() * 1e3) if self.m_f.sum() > 0.0 else 0.0,
             "delta_m_mean_um": self._mean_conjugate_um(delta_m),
             "thick_branch_fraction": float(thick[self.windward & (self.m_f > 0.0)].mean())
-            if (self.windward & (self.m_f > 0.0)).any() else float("nan")})
+            if (self.windward & (self.m_f > 0.0)).any() else float("nan"),
+            "nonrigid_depth_mean_mm": nonrigid_mean, "slurry_thick_fraction": slurry_thick,
+            "rigid_thin_fraction": rigid_thin, "slurry_held_mass_kg": held_mass})
         return released
 
     def _kill(self, dead):
@@ -1024,6 +1055,57 @@
             return np.zeros(0)
         a = self.liquid.rho * self.surface.areas
         return self.m_f / a + self.molten_depth(Te) + self.m_d / a
+
+    def nonrigid_depth(self, T=None, max_crossings=NONRIGID_MAX_CROSSINGS):
+        """Depth of non-rigid material under each patch [m]: how far inward from the wall the material stays more than
+        half liquid (above `Material.T_rigid`), each element's stretch weighed by phi_e (amendment of 2026-10-06).
+
+        A line runs from the patch centre along the inward normal through the elements it actually crosses, leaving each
+        by the face its barycentric coordinates reach first. The P1 field is linear along the line inside an element, so
+        the point where it falls to T_rigid is found exactly instead of being counted in whole elements, and plan fact
+        58's element-size artefact does not arise. Only material continuous with the wall counts: the line stops at the
+        first point at or below T_rigid, at the active mesh's boundary, or after `max_crossings` elements. Each
+        element's stretch counts phi_e of its length, because what it has already fed to the film is in the film
+        account the layer adds on top (fact 58's double count). The thin branch is Girin's dominant ablation, where the
+        rigid core stabilises the disturbances, and slurry is no rigid core: this is the depth the regime test reads."""
+        mat, mesh, surf = self.material, self.mesh, self.surface
+        n = surf.n_patches
+        depth = np.zeros(n)
+        if not n or not mat.melts:
+            return depth
+        T = self.solver.temperature() if T is None else np.asarray(T, dtype=float)
+        T_r = mat.T_rigid
+        c, d = surf.centroids, -surf.normals
+        elem, s_in, live = surf.owner.copy(), np.zeros(n), np.arange(n)
+        for _ in range(max_crossings):
+            if live.size == 0:
+                break
+            e, s0 = elem[live], s_in[live]
+            x = mesh.points[mesh.tets[e]]                                                       # (k, 4, 3)
+            inv = np.linalg.inv(np.stack([x[:, 1] - x[:, 0], x[:, 2] - x[:, 0], x[:, 3] - x[:, 0]], axis=2))
+            grad = np.concatenate([-inv.sum(axis=1, keepdims=True), inv], axis=1)               # grad of lambda_0..3
+            lam = np.einsum("kij,kj->ki", inv, c[live] - x[:, 0])
+            lam = np.concatenate([1.0 - lam.sum(axis=1, keepdims=True), lam], axis=1)           # barycentrics of the centre
+            rate = np.einsum("kij,kj->ki", grad, d[live])                                       # d lambda / ds on the line
+            s_face = np.where(rate < 0.0, -lam / np.where(rate < 0.0, rate, -1.0), np.inf)       # where each face is reached
+            exit_face = s_face.argmin(axis=1)
+            rows = np.arange(len(e))
+            s1 = np.maximum(s_face[rows, exit_face], s0)
+            Tn = T[mesh.tets[e]]
+            T0, dT = (lam * Tn).sum(axis=1), (rate * Tn).sum(axis=1)                            # T on the line: T0 + s dT
+            Ta, Tb = T0 + s0 * dT, T0 + s1 * dT
+            hot = Ta > T_r
+            ends = hot & (Tb <= T_r)
+            with np.errstate(divide="ignore", invalid="ignore"):
+                s_end = np.where(ends, s0 + (Ta - T_r) / np.where(ends, Ta - Tb, 1.0) * (s1 - s0), s1)
+            depth[live] += np.where(hot, self.phi[e] * (s_end - s0), 0.0)
+            pair = mesh._face_elements[mesh.element_faces(e)[rows, exit_face]]                  # across the exit face
+            nxt = np.where(pair[:, 0] == e, pair[:, 1], pair[:, 0])
+            go = hot & ~ends & np.isfinite(s1) & (nxt >= 0)
+            go &= mesh.active[np.maximum(nxt, 0)]
+            elem[live[go]], s_in[live[go]] = nxt[go], s1[go]
+            live = live[go]
+        return depth
 
     def film_thickness_max(self):
         """Thickest film [m] over the patches where a thickness means anything: at least a tenth of the median patch
@@ -1126,6 +1208,10 @@
                 "molten_depth_mean_mm": lm.get("molten_depth_mean_mm", float("nan")),
                 "delta_m_mean_um": lm.get("delta_m_mean_um", float("nan")),
                 "thick_branch_fraction": lm.get("thick_branch_fraction", float("nan")),
+                "nonrigid_depth_mean_mm": lm.get("nonrigid_depth_mean_mm", float("nan")),
+                "slurry_thick_fraction": lm.get("slurry_thick_fraction", float("nan")),
+                "rigid_thin_fraction": lm.get("rigid_thin_fraction", float("nan")),
+                "slurry_held_mass_kg": lm.get("slurry_held_mass_kg", float("nan")),
                 "deep_liquid_kg": lm.get("deep_liquid", 0.0), "deep_mass_kg": float(self.m_d.sum()),
                 "deep_runoff_mass_kg": self.deep_runoff_mass, "deep_surfaced_mass_kg": self.deep_surfaced_mass,
                 "deep_blob_fraction": self.deep_blob_fraction(),
```

In `tests/test_reentry_model_melting.py` (the five new tests and the changed histogram check):

```diff
--- a/tests/test_reentry_model_melting.py
+++ b/tests/test_reentry_model_melting.py
@@ -95,7 +95,17 @@
     assert stats["film_T_max_K"] == pytest.approx(b.surface_temperature().max()) and stats["film_T_max_K"] > b.material.T_liquidus
     rows = np.array(b.source_rows)
     assert rows.shape[1] == 22 and np.all(rows[:, 14] > 0.0) and np.all(np.isfinite(rows[:, 12]))
-    assert b.hist_n.sum() == pytest.approx(b.n_released) and b.hist_m.sum() == pytest.approx(b.sprayed_mass)
+    # the histogram bins 1 um to 10 mm and drops the rest by design. With the rigid substrate on (amendment of 2026-10-06)
+    # this body, 880 K throughout, is slurry to the core: off Girin's closure its film is held from the shear modes,
+    # piles up, and the front-surface Rayleigh-Taylor mode sheds some of it above 10 mm (measured in two unseeded runs:
+    # 5 and 6 droplets, 5.2 and 6.1 % of the mass, 10.1-10.9 mm), so the bins hold exactly what was released inside
+    # their range and nothing is lost below it
+    from reentry_model import spray
+    r_rows, dn_rows, dm_rows = rows[:, 12], rows[:, 13], rows[:, 14]
+    inside = (r_rows >= spray.R_MIN) & (r_rows <= spray.R_MAX)
+    assert b.hist_n.sum() == pytest.approx(dn_rows[inside].sum()) and b.hist_m.sum() == pytest.approx(dm_rows[inside].sum())
+    assert dn_rows.sum() == pytest.approx(b.n_released) and dm_rows.sum() == pytest.approx(b.sprayed_mass)
+    assert (r_rows[~inside] > spray.R_MAX).all() and (r_rows <= 0.25 * 0.05).all()              # above the bins, below R/4
     assert not b.demised() and b.reference_area() < math.pi * 0.05 ** 2                            # the windward face has receded
 
 
@@ -604,3 +614,114 @@
     assert b.sprayed_mass > 0.0 and b.deep_runoff_mass > 0.0
     stats = b.melt_stats()
     assert stats["cascade_mass_kg"] == b.cascade_mass and stats["cascade_passes"] == passes[-1]
+
+
+# ---------------------------------------------------------------------------------------------------------------
+# Amendment of 2026-10-06: the thin branch needs a rigid substrate
+
+
+def linear_field(b, i, T_wall, gradient):
+    """Nodal temperatures falling linearly with depth below patch i's plane: T_wall on it, `gradient` [K/m] inward. P1
+    represents a linear field exactly, so the half-liquid point along the patch's inward normal is known to round-off."""
+    return T_wall + gradient * ((b.mesh.points - b.surface.centroids[i]) @ b.surface.normals[i])
+
+
+def test_the_nonrigid_depth_follows_the_temperature_field_to_the_half_liquid_point(layered_mesh):
+    """From the patch centre straight inward, through the elements the line actually crosses, to where the field falls
+    to T_rigid (829 K on AA7075_range): read inside the element, not counted in whole elements. Inside the second prism
+    layer and several elements down alike; a wall below T_rigid has no slurry under it at all."""
+    b = melting_body(layered_mesh)
+    i = b.i_stag
+    for depth in (0.4e-3, 3e-3):
+        b.solver.set_temperature(linear_field(b, i, 900.0, (900.0 - b.material.T_rigid) / depth))
+        assert b.nonrigid_depth()[i] == pytest.approx(depth, rel=1e-9)
+    b.solver.set_temperature(b.material.T_rigid - 1.0)                       # coherent mush or solid everywhere
+    assert not b.nonrigid_depth().any()
+
+
+def test_the_nonrigid_depth_weighs_each_element_by_what_is_left_of_it(layered_mesh):
+    """phi_e of an element's mass is still in it; the rest went to the film, which the layer counts already, so each
+    element's stretch of the line counts phi_e of its length (plan fact 58's double count, removed here)."""
+    b = melting_body(layered_mesh)
+    i = b.i_stag
+    b.solver.set_temperature(linear_field(b, i, 900.0, (900.0 - b.material.T_rigid) / 3e-3))
+    full = b.nonrigid_depth()[i]
+    b.phi[:] = 0.5
+    assert b.nonrigid_depth()[i] == pytest.approx(0.5 * full, rel=1e-12)
+
+
+def test_the_nonrigid_depth_stops_at_the_first_rigid_point(layered_mesh):
+    """Only slurry continuous with the wall counts: the line stops where the field first falls to T_rigid, whatever
+    lies deeper -- as `molten_depth` stops at the first element that is not fully molten."""
+    b = melting_body(layered_mesh)
+    i = b.i_stag
+    T = linear_field(b, i, 900.0, (900.0 - b.material.T_rigid) / 0.4e-3)
+    s = -(b.mesh.points - b.surface.centroids[i]) @ b.surface.normals[i]    # each node's depth below the patch's plane
+    b.solver.set_temperature(np.where(s > 1e-3, 900.0, T))                 # hot again below the two prism layers
+    assert b.nonrigid_depth()[i] == pytest.approx(0.4e-3, rel=1e-9)
+
+
+def slurry_under_the_surface(b, T_wall=905.0, depth=1e-3):
+    """The surface just below the feed ramp (nothing fully molten, so the liquid layer is the film alone) over material
+    that stays above T_rigid for `depth` below it: slurry, radially."""
+    d = 0.05 - np.linalg.norm(b.mesh.points, axis=1)
+    b.solver.set_temperature(T_wall - (T_wall - b.material.T_rigid) * d / depth)
+
+
+def test_a_film_on_slurry_deeper_than_the_conjugate_depth_is_thick_in_the_melt_step(layered_mesh):
+    """At 70 km, under Girin's closure on every windward patch (conjugate depth 0.17-0.30 mm), a 20 um film over a
+    millimetre of slurry: with the rigid substrate on, every wet windward patch is thick, and only because of the
+    slurry; with it off the liquid layer is the film alone and every such patch is thin."""
+    assert body.MeltSettings().rigid_substrate
+    seen = {}
+    for on in (True, False):
+        b = melting_body(layered_mesh, rigid_substrate=on)
+        sim, a = girin_state(b)
+        slurry_under_the_surface(b)
+        flow, delta_m = flow_and_delta_m(b, a)
+        assert np.isfinite(delta_m[b.windward]).all() and delta_m[b.windward].max() < 0.5e-3
+        b.m_f = np.where(b.windward, 2e-5 * b.liquid.rho * b.surface.areas, 0.0)
+        h_p = deep_stage_args(b, delta_m)[4]
+        b._film_and_spray(sim.t, 0.5, a, h_p, flow, delta_m)
+        seen[on] = b
+    on, off = seen[True], seen[False]
+    assert on.last_melt["thick_branch_fraction"] == 1.0 and off.last_melt["thick_branch_fraction"] == 0.0
+    assert on.last_melt["slurry_thick_fraction"] == 1.0 and on.last_melt["rigid_thin_fraction"] == 0.0
+    assert 0.9 < on.last_melt["nonrigid_depth_mean_mm"] < 1.1                 # the millimetre of slurry, read back
+    assert on.last_melt["slurry_held_mass_kg"] == 0.0                         # Girin's closure everywhere: nothing held
+    assert np.isnan(off.last_melt["slurry_thick_fraction"]) and np.isnan(off.last_melt["nonrigid_depth_mean_mm"])
+    from reentry_model import spray
+    wet = on.windward & (on.last_b >= spray.B_MIN)
+    assert not np.isin(on.last_spray.branch[wet], [spray.BRANCH_THIN, spray.BRANCH_RAREFIED]).any()
+    assert on.last_nonrigid is not None and on.last_nonrigid.shape == (on.surface.n_patches,) and off.last_nonrigid is None
+    for key in ("nonrigid_depth_mean_mm", "slurry_thick_fraction", "rigid_thin_fraction", "slurry_held_mass_kg"):
+        assert key in on.melt_stats()
+
+
+def test_off_girins_closure_a_film_on_slurry_does_not_spray_and_one_on_a_rigid_wall_does(layered_mesh):
+    """At 30 s of the 100 mm flight no patch has Girin's closure, so there is no conjugate depth. A film over slurry is
+    held from the shear modes, and the mass held is recorded; the same film over a rigid wall (surface below T_rigid)
+    takes the thin mode as before, and so does the film over slurry with the rigid substrate off."""
+    from reentry_model import spray
+    out = {}
+    for label, on, T_wall in (("slurry", True, 905.0), ("rigid", True, 820.0), ("off", False, 905.0)):
+        b = melting_body(layered_mesh, rigid_substrate=on)
+        sim, a = girin_state(b, t=30.0)
+        slurry_under_the_surface(b, T_wall=T_wall)
+        flow, delta_m = flow_and_delta_m(b, a)
+        assert not np.isfinite(delta_m).any()
+        b.m_f = np.where(b.windward, 2e-5 * b.liquid.rho * b.surface.areas, 0.0)
+        film = b.m_f.copy()
+        h_p = deep_stage_args(b, delta_m)[4]
+        b._film_and_spray(sim.t, 0.5, a, h_p, flow, delta_m)
+        out[label] = (b, film)
+    b, film = out["slurry"]
+    wet = b.windward & (b.last_b >= spray.B_MIN)
+    assert wet.any() and not np.isin(b.last_spray.branch[wet], [spray.BRANCH_THIN, spray.BRANCH_RAREFIED]).any()
+    assert b.last_melt["slurry_held_mass_kg"] == pytest.approx(b.m_f[wet].sum() + b.last_spray.dm[wet].sum(), rel=1e-12)
+    for label in ("rigid", "off"):
+        b, film = out[label]
+        wet = b.windward & (b.last_b >= spray.B_MIN)
+        assert np.isin(b.last_spray.branch[wet], [spray.BRANCH_THIN, spray.BRANCH_RAREFIED]).all()
+        assert b.sprayed_mass > 0.0
+    assert out["rigid"][0].last_melt["slurry_held_mass_kg"] == 0.0
```

In `tests/test_reentry_model_fenicsx.py` (the 2026-10-05 test with the rule off, and the new two-backend test):

```diff
--- a/tests/test_reentry_model_fenicsx.py
+++ b/tests/test_reentry_model_fenicsx.py
@@ -196,7 +196,12 @@
     stays bounded: before the amendment it grew from sub-step to sub-step and reached 3.9e25 per second in this setting
     (measured 2026-10-05); with it 2.7e7, set by the deepest piles' G b^3 / (3 mu) term, not by a vanishing film.
     Measured: the same 7 983 active elements, mass 5.8e-10, sprayed mass 5.9e-9, runoff 7.5e-10, film 1.7e-8, deep
-    account 2.8e-10, temperatures 2.3e-4 K, balances 5.1e-9 and 1.1e-11."""
+    account 2.8e-10, temperatures 2.3e-4 K, balances 5.1e-9 and 1.1e-11.
+    The rigid substrate (amendment of 2026-10-06) is off here: with the body's interior at 850 K, above AA7075_range's
+    829 K half-liquid point, the whole body is slurry, every wet windward patch takes the thick branch, and films of
+    about a micron that straddle the spray floor B_MIN are released whole by one backend and not by the other from the
+    sixth step on -- 8 to 29 patches, the runoff then differing by 1.6e-4 and the temperatures by 0.29 K (measured
+    2026-10-06). That is the model's threshold, not the backends'; the test after this one checks the two with it on."""
     pytest.importorskip("cantera")
     from reentry_model import body, film, heating, mesh
     from test_reentry_model_coupled import MASS_100MM, simulator
@@ -221,7 +226,8 @@
     for name in ("skfem", "fenicsx"):
         np.random.seed(12345)
         m = mesh.VolumeMesh(coarse_sphere_mesh.points, coarse_sphere_mesh.tets, dict(coarse_sphere_mesh.params))
-        b = body.MeltingBody(m, material.Material.from_drama_json("AA7075_range"), thermal.thermal_solver(name), MASS_100MM)
+        b = body.MeltingBody(m, material.Material.from_drama_json("AA7075_range"), thermal.thermal_solver(name), MASS_100MM,
+                             settings=body.MeltSettings(rigid_substrate=False))
         sim = simulator(b, t_max=90.0)
         sim.advance(50.0)
         p = b.mesh.points
@@ -243,3 +249,44 @@
     assert f.m_d.sum() == pytest.approx(s.m_d.sum(), rel=1e-4)
     assert np.abs(f.solver.temperature() - s.solver.temperature()).max() < 0.5
     assert abs(s.energy_balance_residual()) < 1e-8 and abs(f.energy_balance_residual()) < 1e-8
+
+
+def test_the_rigid_substrate_matches_the_skfem_backend(coarse_sphere_mesh):
+    """The rigid substrate (amendment of 2026-10-06) in both backends: the 6 mm pool at 69.8 km of the test above, with
+    the rest of the body at 820 K, below AA7075_range's 829 K half-liquid point, so slurry forms only where the pool
+    and the heating raise it -- ten steps of 0.05 s. The non-rigid depth reaches into the body, both backends make the
+    same branch decision on every patch at every step, and the active set, mass, sprayed mass, runoff, film and
+    temperatures agree with both energy balances exact. Measured: mass 6.7e-10, sprayed mass 1.0e-8, runoff 1.2e-8,
+    film 6.1e-8, temperatures 1.8e-4 K."""
+    pytest.importorskip("cantera")
+    from reentry_model import body, heating, mesh
+    from test_reentry_model_coupled import MASS_100MM, simulator
+    out = {}
+    for name in ("skfem", "fenicsx"):
+        np.random.seed(12345)
+        m = mesh.VolumeMesh(coarse_sphere_mesh.points, coarse_sphere_mesh.tets, dict(coarse_sphere_mesh.params))
+        b = body.MeltingBody(m, material.Material.from_drama_json("AA7075_range"), thermal.thermal_solver(name), MASS_100MM)
+        assert b.settings.rigid_substrate
+        sim = simulator(b, t_max=90.0)
+        sim.advance(50.0)
+        p = b.mesh.points
+        r = np.linalg.norm(p, axis=1)
+        b.solver.set_temperature(np.where((r > 0.05 - 6e-3) & (p[:, 0] > 0.5 * r), 960.0, 820.0))
+        b.energy0 = b.energy()
+        model = heating.PhysicsHeating()
+        branches = []
+        for _ in range(10):
+            sim.advance(0.05)
+            a = sim.aero_state(sim.t, sim.y[:3], sim.y[3:])
+            b.advance(sim.t, 0.05, model.evaluate(a, b.theta, b.surface_temperature(), b.nose_radius(), T_mean=b.mean_temperature()),
+                      state=a)
+            branches.append(None if b.last_spray is None else b.last_spray.branch.copy())
+        out[name] = (b, branches)
+    (s, bs), (f, bf) = out["skfem"], out["fenicsx"]
+    assert s.last_nonrigid is not None and s.last_nonrigid.max() > 6e-3 and f.sprayed_mass > 0.0
+    assert all((x is None and y is None) or np.array_equal(x, y) for x, y in zip(bs, bf))
+    assert np.array_equal(s.mesh.active, f.mesh.active)
+    assert f.mass(0.0) == pytest.approx(s.mass(0.0), rel=1e-6) and f.sprayed_mass == pytest.approx(s.sprayed_mass, rel=1e-5)
+    assert f.runoff_mass == pytest.approx(s.runoff_mass, rel=1e-5) and f.m_f.sum() == pytest.approx(s.m_f.sum(), rel=1e-4)
+    assert np.abs(f.solver.temperature() - s.solver.temperature()).max() < 0.5
+    assert abs(s.energy_balance_residual()) < 1e-8 and abs(f.energy_balance_residual()) < 1e-8
```

### Measured (2026-10-06)

- The non-rigid depth is exact to 1e-9 relative on a linear field, halves with φ and stops at the first rigid point
  (fact 79). Unit tier 251 passed, with only Task 11's five known reference failures; FEniCSx tier 12 passed (fact 80).
- With `--rigid-substrate off` both flights are bit-identical to the copy before the change (fact 81).
- 100 mm flight to 120 s, on against off: sprayed mass +1.2 %, the body at 120 s −3.0 %, droplets −32 % (17.89 to 12.10
  million, almost all of it the thin branch's), median radius by number +3.8 % and by mass +4.6 %; the thick branch's
  share of the sprayed mass 85.1 to 88.2 %, the thin branch's 10.4 to 2.5 %, the front-surface Rayleigh–Taylor mode's
  4.5 to 9.4 %; under Girin's closure the thick share of the wet windward patches rises from a median 27.4 to 62.8 %
  over a non-rigid depth of 8.1 mm (fact 82). The 50 mm flight, never under the closure: two-thirds of the sprayed
  mass leaves by the Rayleigh–Taylor mode and the mass median radius is 7.4 times larger (fact 83). Scheil's curve:
  the slurry band is 12.9 K instead of 79 K and the slurry-only thick share 6 % instead of 23 % (fact 84).
- The held film is shed by the front-surface Rayleigh–Taylor mode (fact 85); the two backends part at the spray floor
  when the whole body is slurry and agree to 1e-8 when it is not (fact 86). Energy balances exact throughout.

### Not done here — for Asha to decide (fact 87)

(a) Whether the front-surface Rayleigh–Taylor mode should also be held over slurry off the closure — decide this first;
(b) the material branch splits are quoted with; (c) the time-step series with the rule on; (d) sub-plans 14 and 15;
(e) the 64-element cap; (f) the φ weighting's assumption.

*Answered on 2026-10-07 (`00-shared-context.md` fact 88):* (a) it is not held; (b) `AA7075_scheil` is the melting
default; (c) the series repeated (facts 91–95); (d) sub-plans 14 and 15 amended.

## Amendment of 2026-10-07 — freeze-back never leaves a negative film

> Found by the time-step series of 2026-10-07 (facts 88–96 in `00-shared-context.md`). The code below is the tested
> code, as a diff against the copy of the rigid-substrate amendment of 2026-10-06 with sub-plan 13's Scheil default of
> this date applied.

**What failed.** The 0.00625 s run of the series stopped at its first fine step, at 49.50625 s, with the thermal
solver's input check "film mass must be one finite non-negative value per node" (exit code 2, fact 77 (c)). A probe that
checks the film and deep accounts after every stage of the melt step, re-running the same seeded flight, found the first
bad value right after `_freeze_back`: a film of −1.29e-25 kg on one windward patch, 51° from the stagnation point, whose
deep liquid freeze-back had just taken whole. Nothing was non-finite; the non-rigid depth, `film.lubrication`'s outputs
and the runoff transport were all clean up to that point.

**Why.** The deep-runoff amendment of 2026-10-02 caps what freeze-back takes from a patch at its film plus its deep
liquid, takes the deep liquid first and debits the film by the rest, computed as the capped amount less the deep part:
the film is left with m_f − ((m_f + m_d) − m_d). In floating point (m_f + m_d) − m_d is m_f rounded to the precision of
m_d, so when the deep liquid dwarfs the film, the film's part comes out one rounding unit above the film for about half
of all such pairs (measured on 200 000 random pairs: 99 860), and the film is left at about −1e-25 kg. The spray stage's
clean-up, which zeroes films below 1e-30 kg, runs before freeze-back, so the negative reaches the solver. It needs a
patch holding both film and deep liquid whose liquid all freezes back in one step. No earlier run stopped on it, and
none of the 2026-10-06 flights is changed by the fix: re-run on the fixed code, the 100 mm and 50 mm flights on the
linear range with the rule on and off and the 100 mm flight on Scheil's curve with the rule off are bit-identical to the
originals in every history column, result and source-table column (fact 89). The Scheil series met it at its first fine
step.

**The change.** Each account gives up at most what it holds, and the film's part is taken directly: `from_deep =
min(m_d, wanted)`, `from_film = min(m_f, wanted - from_deep)`, `taken = from_deep + from_film`. Each subtraction then
takes from a number no smaller than what it takes, which correctly rounded arithmetic never rounds below zero. Where the
cap does not bind, the amounts are those of before; where it binds they can differ in the last bit.

**Test.** `test_freeze_back_of_film_and_deep_liquid_leaves_neither_negative`: one patch with 2.34e-14 kg of film over
1.88e-9 kg of deep liquid — a pair the old arithmetic gets wrong, which the test asserts — all of it freezing back into
an owner element with room. Both accounts must be exactly zero afterwards and the owner must gain their mass. It failed
before the change with the film negative and passes after it. Unit tier 252 passed, with Task 11's five known reference
failures; FEniCSx tier 12 passed.

In `reentry_model/body.py`:

```diff
--- a/reentry_model/body.py
+++ b/reentry_model/body.py
@@ -539,14 +539,19 @@
         room = np.maximum(1.0 - self.phi, 0.0) * self.element_mass
         with np.errstate(divide="ignore", invalid="ignore"):
             share = np.where(have > 0.0, np.minimum(want, room) / np.where(have > 0.0, have, 1.0), 0.0)
-        taken = np.minimum(solid * share[owner], self.m_f + self.m_d)
+        # each account gives up at most what it holds, the film's part taken directly: as the capped sum less the deep
+        # part, (m_f + m_d) - m_d, it came out one rounding unit above m_f whenever m_d dwarfs m_f, leaving a negative
+        # film the solver rejects (amendment of 2026-10-07)
+        wanted = solid * share[owner]
+        from_deep = np.minimum(self.m_d, wanted)
+        from_film = np.minimum(self.m_f, wanted - from_deep)
+        taken = from_deep + from_film
         if not taken.any():
             return 0.0
         per_element = np.bincount(owner, taken, self.mesh.n_elements)
         self.phi = np.clip(self.phi + per_element / self.element_mass, 0.0, 1.0)
-        from_deep = np.minimum(self.m_d, taken)
         self.m_d = self.m_d - from_deep
-        self.m_f = self.m_f - (taken - from_deep)
+        self.m_f = self.m_f - from_film
         # the film arrives in the element it froze onto, which is colder than the surface it left
         self._defer_to_elements(np.bincount(owner, taken * h_p, self.mesh.n_elements) - per_element * h_e)
         return float(taken.sum())
```

In `tests/test_reentry_model_melting.py`:

```diff
--- a/tests/test_reentry_model_melting.py
+++ b/tests/test_reentry_model_melting.py
@@ -725,3 +725,33 @@
         assert np.isin(b.last_spray.branch[wet], [spray.BRANCH_THIN, spray.BRANCH_RAREFIED]).all()
         assert b.sprayed_mass > 0.0
     assert out["rigid"][0].last_melt["slurry_held_mass_kg"] == 0.0
+
+
+# ---------------------------------------------------------------------------------------------------------------
+# Amendment of 2026-10-07: freeze-back never leaves a negative film
+
+
+def test_freeze_back_of_film_and_deep_liquid_leaves_neither_negative(layered_mesh):
+    """Freeze-back draws a patch's deep liquid first and its film after, capped at what the two hold. Taking the cap as
+    their sum and the film's part as that sum less the deep liquid leaves the film one rounding unit below zero whenever
+    the deep liquid dwarfs the film -- for about half of all such pairs -- and the solver's film-mass check rejects the
+    negative: it stopped the 0.00625 s run of the 2026-10-07 time-step series at its first fine step, a film of -1.3e-25
+    kg on one patch. A patch whose film and deep liquid all freeze back must be left with exactly none of either, and its
+    owner element must gain exactly their mass."""
+    b = melting_body(layered_mesh)
+    b.solver.set_temperature(800.0)                                 # below the feed ramp: every patch's liquid freezes
+    i = b.i_stag
+    owner = b.surface.owner[i]
+    b.phi[owner] = 0.5                                              # room for it in the owner element
+    b.m_f[:], b.m_d[:] = 0.0, 0.0
+    b.m_f[i], b.m_d[i] = 2.3351899704880006e-14, 1.8789852661499484e-09   # (m_f + m_d) - m_d > m_f in floating point
+    assert b.m_f[i] - ((b.m_f[i] + b.m_d[i]) - b.m_d[i]) < 0.0      # the pair the old arithmetic got wrong
+    before = b.m_f[i] + b.m_d[i]
+    solid = (1.0 - b.material.feed_fraction(b.film_temperature())) * (b.m_f + b.m_d)
+    want = np.bincount(b.surface.owner, solid, b.mesh.n_elements)
+    _, _, _, h_e, h_p = deep_stage_args(b, np.full(b.surface.n_patches, np.nan))
+    frozen = b._freeze_back(want, solid, h_e, h_p)
+    assert (b.m_f >= 0.0).all() and (b.m_d >= 0.0).all()
+    assert b.m_f[i] == 0.0 and b.m_d[i] == 0.0
+    assert frozen == pytest.approx(before, rel=1e-15)
+    assert (b.phi[owner] - 0.5) * b.element_mass[owner] == pytest.approx(before, rel=1e-9)
```

---

**Depends on:** Tasks 1 through 7 — all of them. This is the integration point of the whole Step 3 plan.
**Produces, for later tasks:** the `MeltingBody` class: feed/freeze mass transfer between elements, the film's temperature and energy balance, spraying integration, the element-death cascade with hand-over of mass and heat to neighbors, nose-cap fitting, and deferred-load bookkeeping. Task 10's coupled loop is built directly on this class.
**Character:** physics and bookkeeping, the densest task in the plan — its code block runs to roughly 1,200 lines, larger than any other single task.
**Read before implementing:** the master plan's own text says, in its own words, "Read facts 25 and 26 before reading `melt_step`" — treat that as a direct instruction, not a suggestion. Facts 4, 5, 25, and 26 in the shared context describe two specific wrong implementations that were tried and rejected here during prototyping (debiting only the destination element; booking the film at mixture rather than liquid-only enthalpy), each of which silently reproduced a different bug. An agent working from this task's text alone, without that history, could plausibly re-derive and reintroduce either mistake while "simplifying" the bookkeeping.
**Recommendation carried over from the build-plan review:** do not split this task further across multiple sub-agents despite its size — the feed/freeze bookkeeping, film energy balance, and element-death hand-over are tightly coupled to each other, and a review checkpoint should happen after this task and before Task 10 begins.
**Refinement goal for the sub-agent:** turn the section below into a standalone implementation plan, and reproduce the "read facts 25/26 first" instruction and the two rejected-approach warnings directly in that plan's text so a future implementer doesn't have to go find them in the master document.

---

### Task 9: The melting body

**Files:**
- Modify (replace): `reentry_model/body.py`
- Modify: `reentry_model/trajectory.py` (three statements), `reentry_model/aero.py` (the Mach-1 step and the shape factor)
- Create: `reentry_model/data/atdb_disc.json` (extracted from DRAMA)
- Test: `tests/test_reentry_model_melting.py` (new); `tests/test_reentry_model_aero.py` (one assertion); `tests/test_reentry_model_fenicsx.py::test_melting_run_matches_the_skfem_backend` (from Task 3) now runs in `fenicsx_env`

**Interfaces:**
- Consumes: Tasks 1–7 (`mesh.deactivate/surface/active_nodes`, `Material.feed_fraction/liquid_fraction/enthalpy/h_liquid/liquid`, `thermal` `set_fractions/element_energies/step(nodal_load=)/pinned/temperature`, `surface_flow.SurfaceFlow`, `spray.SprayModel/melt_layer/source_rows/histogram/N_BINS`, `film.lubrication/Runoff`), `heating.HeatingResult`, `trajectory.AeroState`.
- Produces: `body.REMOVAL_NAMES = ("girin", "instant")`, `SIZE_FEEDBACK_NAMES = ("current", "initial")`, `NOSE_CAP_ANGLE = 30.0`, `NOSE_CAP_FACTOR = 1.67`, `CP_MAX_NEWTONIAN = 1.84`, `PHI_MIN = 1e-3`, `PHI_DEATH = 0.05`, `NEAREST_PATCHES = 4`, `LOAD_DT_MAX = 1000.0`, `FEED_DEPTH_NAMES = ("conjugate", "all")`; `fit_sphere(points) -> (centre, radius)`; `MeltSettings(removal, runoff, demise_fraction, particles, size_feedback="current", feed_depth="conjugate")`; `MeltingBody(mesh, material, solver, mass_kg, flow=None, spray_model=None, settings=None, T0, emissivity, T_ambient, v_hat)` with `.advance(t, dt, loads, state=None)`, `.melt_step(t, dt, state)`, `.mass(t)`, `.reference_area()`, `.reference_length()` (2 R_eq, or None with `initial`), `.nose_radius()` (the windward-cap fit, or R₀ with `initial`), `.newtonian_drag() -> (C_D, projected area)` (modified Newtonian over the windward convex hull), `.drag_shape_factor()` (that value over the meshed sphere's own, so exactly 1 while intact and 2.00 at the flat limit; 1 with `initial`), `.equivalent_radius()`, `.energy()` (FEM + film), `.film_temperature()` (the surface's own: the film is thermally thin), `.film_enthalpy(h_node=None)` (per patch, the mean of the *liquid* nodal enthalpies), `.film_energy()`, `.film_frozen_fraction()`, `.film_blob_fraction()` (the share of the film deeper than its patch is wide), `.molten_depth(Te=None, max_levels=8)` (the contiguous liquid depth inward from each patch, hopping element to element and stopping at the first one that is not molten), `.liquid_layer_depth()` (that depth plus the film, which is what the branch test and the Rayleigh-Taylor criterion see), `.region_extent(selected)` (each patch's own contiguous region as an equivalent diameter 2 sqrt(A/pi), which is the domain a wave has to fit inside -- the molten region, not the facet) and `.unstable_region_extent(unstable)` (its maximum, for the reported diagnostics; facts 30, 32 and 33), `.total_applied_load()`, `.total_dropped_load()` (the deferred-load accounts the energy balance closes against), `.mean_temperature()`, `.melt_front_depth()`, `.film_thickness_max/mean()`, `.demised()`, `.energy_balance_residual()`, `.melt_stats() -> dict` (the `coupled.MELT_COLUMNS` values plus `mass_kg`), attributes `phi, m_f, surface, theta, t_hat, windward, mass0, mass_centre, transverse_radius, fitted_nose_radius, cap_nose_radius, sprayed_mass, runoff_mass, removed_mass, removed_enthalpy, n_released, source_rows, hist_n, hist_m, melt_onset, spray_onset, consumed, last_flow, last_spray, last_b, last_melt, pending_load, liquid, flow, spray, runoff`; `Body.reference_area()`/`reference_length()`/`drag_shape_factor()` in the protocol (`ConstantBody`/`ThermalBody` return None; `ThermalBody.nose_radius()` returns the sphere's radius; `ThermalBody.advance` accepts `state=None`); `Simulator.aero_state` uses the body's reference area, Knudsen length and drag shape factor, and gives a consumed body no drag; `aero.drag_coefficient(kn, ma, tables, bridging, shape_factor=1.0)` scales the continuum entry only.

The film's temperature, its re-solidification and the way every transfer is booked are the subject of measured facts 25 and 26: the film has no energy equation of its own (its *mass* goes to the solver, which carries it on the boundary nodes with the liquid capacity), it holds `enthalpy_liquid` wherever it sits, the feed and the freeze are netted into one transfer per element, each transfer books only the enthalpy difference it carries and books it at the destination, and no deferred load may move a node more than `LOAD_DT_MAX` in one step. Read those two facts before reading `melt_step`.

- [ ] **Step 1: Write the failing tests**

Create `tests/test_reentry_model_melting.py`:

```python
"""body.MeltingBody: feed, instant removal in the lumped limit, film and spraying, element death with hand-over,
energy and mass balances, demise, and the shape feedback -- Knudsen length, nose-cap radius and drag (sections 8-11,
18)."""
import math
import os

import numpy as np
import pytest

from reentry_model import body, heating, material, mesh, thermal
from reentry_model.body import PHI_DEATH as PHI_DEATH_FOR_TEST
from test_reentry_model_coupled import MASS_100MM, simulator

pytest.importorskip("cantera")


@pytest.fixture(scope="module")
def layered_mesh(tmp_path_factory):
    return mesh.sphere_mesh(0.05, 4e-3, 20e-3, str(tmp_path_factory.mktemp("meshes")), layers=2)


def melting_body(the_mesh, name="AA7075_range", k_scale=1.0, **settings):
    mat = material.Material.from_drama_json(name)
    if k_scale != 1.0:
        mat.k_table = mat.k_table * k_scale
    m = mesh.VolumeMesh(the_mesh.points, the_mesh.tets, dict(the_mesh.params), the_mesh.element_layer.copy())
    return body.MeltingBody(m, mat, thermal.thermal_solver("skfem"), MASS_100MM, settings=body.MeltSettings(**settings))


def test_settings_and_construction(layered_mesh):
    with pytest.raises(ValueError):
        body.MeltSettings(removal="magic")
    with pytest.raises(ValueError):
        body.MeltSettings(demise_fraction=1.5)
    with pytest.raises(ValueError):
        body.MeltingBody(layered_mesh, material.Material.from_drama_json("AA7075_nomelt"), thermal.thermal_solver("skfem"), MASS_100MM)
    b = melting_body(layered_mesh)
    assert b.mass(0.0) == pytest.approx(b.mass0) and b.mass0 == pytest.approx(MASS_100MM, rel=3e-3)
    assert b.reference_area() == pytest.approx(math.pi * 0.05 ** 2, rel=2e-3) and b.equivalent_radius() == pytest.approx(0.05, rel=1e-3)
    assert b.m_f.shape == (b.surface.n_patches,) and b.m_f.sum() == 0.0 and not b.demised() and b.melt_front_depth() == 0.0
    stats = b.melt_stats()
    assert stats["mass_kg"] == b.mass(0.0) and stats["n_active_elements"] == layered_mesh.n_elements and stats["sprayed_mass_kg"] == 0.0


def test_instant_removal_in_the_lumped_limit_follows_q_over_l(layered_mesh):
    """k x 1e4, uniform 3e5 W/m2 on the whole surface, AA7075: once the body sits on the 850 K plateau the mass leaves
    at (Q_conv - Q_rad) / L_f with the geometry intact (SESAM's lumped law), the energy balance exact."""
    b = melting_body(layered_mesh, "AA7075", k_scale=1e4, removal="instant", runoff=False, size_feedback="initial")
    q = np.full(b.surface.n_patches, 3e5)
    loads = heating.HeatingResult(q, 3e5, 0.0, 0.0, 0.0)
    t, dt = 0.0, 0.5
    while b.melt_onset is None:
        t += dt
        b.advance(t, dt, loads)
    m1 = b.mass(0.0)
    for _ in range(20):
        t += dt
        b.advance(t, dt, loads)
    Q_net = b.last.Q_conv - b.last.Q_rad
    assert (m1 - b.mass(0.0)) / (20 * dt) == pytest.approx(Q_net / 4e5, rel=2e-2)
    assert b.removed_mass == pytest.approx(b.mass0 - b.mass(0.0)) and b.m_f.sum() == 0.0 and b.sprayed_mass == 0.0
    assert abs(b.energy_balance_residual()) < 1e-8 and b.mesh.n_active == layered_mesh.n_elements    # nothing dead yet: phi shrinks uniformly
    assert 0.5 < b.phi.min() < b.phi.max() < 1.0 and 848.0 < b.mean_temperature() < 852.0
    assert b.reference_area() == pytest.approx(math.pi * 0.05 ** 2, rel=2e-3)
    assert b.melt_onset is not None and b.spray_onset is None


def test_film_spraying_death_and_balances(layered_mesh):
    """Physics-mode loads at 71 km on the real body for 8 s: film forms, is stripped, the outer layer dies on the
    windward face with the film handed over, mass and energy balances hold, the source table grows."""
    b = melting_body(layered_mesh)
    sim = simulator(b, t_max=60.0)
    sim.advance(43.0)                                                                 # to 71 km on the constant-mass trajectory
    b.solver.set_temperature(880.0)                                                   # a hot body: the surface melts at once
    b.energy0 = b.energy()                                                            # the balance is reckoned from here
    model = heating.PhysicsHeating()
    for _ in range(16):
        sim.advance(0.5)
        a = sim.aero_state(sim.t, sim.y[:3], sim.y[3:])
        loads = model.evaluate(a, b.theta, b.surface_temperature(), 0.05, T_mean=b.mean_temperature())
        b.advance(sim.t, 0.5, loads, state=a)
    stats = b.melt_stats()
    assert b.sprayed_mass > 0.0 and b.n_released > 0.0 and len(b.source_rows) > 0 and stats["n_released"] == b.n_released
    assert b.removed_mass == pytest.approx(b.sprayed_mass) and b.mass(0.0) == pytest.approx(b.mass0 - b.sprayed_mass, rel=1e-6)
    assert b.mesh.n_active < layered_mesh.n_elements and b.mesh.element_layer[~b.mesh.active].max() <= 1     # only layer elements have died
    assert abs(b.energy_balance_residual()) < 1e-7
    assert b.surface.n_patches == b.m_f.size == len(b.solver.areas) and b.theta.size == b.surface.n_patches
    assert stats["theta_cr_deg"] < 45.0 and stats["spraying_area_m2"] > 0.0 and 20.0 < stats["r_median_um"] < 2000.0
    assert stats["closure_fraction_girin"] + stats["closure_fraction_couette_slip"] + stats["closure_fraction_couette_fm"] == pytest.approx(1.0)
    assert stats["kn_body"] > 0.0 and stats["kn_local_stag"] < stats["kn_body"] and stats["re_shock"] > 0.0      # the wall gas is compressed
    assert stats["flow_branch"] in (0.0, 1.0, 2.0) and stats["p_w_stag_Pa"] > stats["kn_body"] * 0.0
    assert stats["film_thickness_max_mm"] >= stats["film_thickness_mean_mm"] >= 0.0
    h_out = stats["removed_enthalpy_J"] / b.sprayed_mass                                           # the droplets leave at the film's
    assert b.material.h_liquid < h_out <= b.material.enthalpy(b.surface_temperature().max())       # own temperature, superheat and all
    assert stats["film_T_max_K"] == pytest.approx(b.surface_temperature().max()) and stats["film_T_max_K"] > b.material.T_liquidus
    rows = np.array(b.source_rows)
    assert rows.shape[1] == 22 and np.all(rows[:, 14] > 0.0) and np.all(np.isfinite(rows[:, 12]))
    assert b.hist_n.sum() == pytest.approx(b.n_released) and b.hist_m.sum() == pytest.approx(b.sprayed_mass)
    assert not b.demised() and b.reference_area() < math.pi * 0.05 ** 2                            # the windward face has receded


def test_a_dying_element_leaves_its_last_solid_on_the_face_it_exposes(layered_mesh):
    """An element is removed while it still holds up to PHI_DEATH of its material (typically far less). That remainder
    has to end up where the surface receded to -- on the faces the element itself has just exposed -- exactly like the
    film it was already carrying, and it does, because the remainder is credited to the element's own patches first and
    is then handed over with them (spec amendment 20). Held below the feed ramp so that nothing else melts this step."""
    b = melting_body(layered_mesh, name="AA7075")
    owner = int(b.surface.owner[0])
    b.m_f[:] = 0.0
    b.solver.set_temperature(b.material.T_feed - 30.0)             # below the feed ramp: no other element melts
    b.phi[owner] = 0.5 * PHI_DEATH_FOR_TEST                        # below the death threshold, with solid left
    b.solver.set_fractions(b.phi)
    rest = float(b.phi[owner] * b.element_mass[owner])
    faces = b.mesh.element_faces([owner])[0]
    mass_before = b.mass(0.0)
    b.melt_step(0.0, 0.5, None)
    exposed = np.array([q for q in b.patch_of_face[faces] if q >= 0])
    assert not b.mesh.active[owner] and exposed.size > 0            # it died and uncovered at least one face
    assert b.m_f.sum() == pytest.approx(rest) and b.m_f[exposed].sum() == pytest.approx(rest)
    assert b.mass(0.0) + b.removed_mass == pytest.approx(mass_before, rel=1e-12)


def test_demise_and_consumption(layered_mesh):
    b = melting_body(layered_mesh, "AA7075", k_scale=1e4, removal="instant", runoff=False, demise_fraction=0.5, size_feedback="initial")
    loads = heating.HeatingResult(np.full(b.surface.n_patches, 2e6), 2e6, 0.0, 0.0, 0.0)
    t = 0.0
    while not b.demised():
        t += 0.5
        b.advance(t, 0.5, loads)
        assert t < 200.0
    assert b.mass(0.0) < 0.5 * b.mass0 and not b.consumed and b.mesh.n_active > 0


def test_size_feedback_nose_fit_and_knudsen_length(layered_mesh):
    """Intact sphere: the windward-cap fit returns R0 and the Knudsen length is D0; a flattened front (the cap patches
    pushed onto the plane x = 0.8 R0, then 0.6 R0) fits a larger radius, capped at NOSE_CAP_FACTOR x the transverse
    radius; with size_feedback 'initial' both stay at their initial values. The trajectory's Kn follows the body's length."""
    b = melting_body(layered_mesh)
    assert b.nose_radius() == pytest.approx(0.05, rel=5e-3) and b.reference_length() == pytest.approx(0.1, rel=1e-3)
    assert b.transverse_radius == pytest.approx(0.05, rel=2e-3) and np.linalg.norm(b.mass_centre) < 1e-3
    centre, r = body.fit_sphere(b.surface.centroids)
    assert np.linalg.norm(centre) < 1e-4 and r == pytest.approx(0.05, rel=5e-3)
    original = b.surface.centroids.copy()
    flat = original.copy()
    flat[original[:, 0] > 0.8 * 0.05, 0] = 0.8 * 0.05
    b.surface.centroids = flat
    b._fit_nose()
    assert b.fitted_nose_radius > 0.06                                             # flat centre + curved shoulders: ~84 mm
    assert b.nose_radius() == pytest.approx(min(b.fitted_nose_radius, body.NOSE_CAP_FACTOR * b.transverse_radius))
    flat = original.copy()
    flat[original[:, 0] > 0.6 * 0.05, 0] = 0.6 * 0.05                               # a wider flat face: the fit exceeds the cap
    b.surface.centroids = flat
    b._fit_nose()
    assert b.fitted_nose_radius > body.NOSE_CAP_FACTOR * b.transverse_radius and b.nose_radius() == pytest.approx(body.NOSE_CAP_FACTOR * b.transverse_radius)
    fixed = melting_body(layered_mesh, size_feedback="initial")
    fixed.surface.centroids = flat
    fixed._fit_nose()
    assert fixed.nose_radius() == 0.05 and fixed.reference_length() is None
    with pytest.raises(ValueError):
        body.MeltSettings(size_feedback="magic")
    # the trajectory's Knudsen number uses the body's length: halve the body's mass and compare
    sim = simulator(b)
    a0 = sim.aero_state(0.0, sim.y[:3], sim.y[3:])
    b.phi *= 0.125                                                                  # equivalent diameter halves
    a1 = sim.aero_state(0.0, sim.y[:3], sim.y[3:])
    assert a1.kn == pytest.approx(2.0 * a0.kn, rel=1e-6) and a0.kn == pytest.approx(0.0298, rel=0.02)
    assert melt_stats_has_radii(b)


def melt_stats_has_radii(b):
    s = b.melt_stats()
    return s["nose_radius_mm"] > 0.0 and s["transverse_radius_mm"] > 0.0 and s["fitted_nose_radius_mm"] > 0.0


def test_drag_shape_factor_between_the_two_atdb_endpoints(layered_mesh):
    """The shape factor is exactly 1 for the intact sphere (so the ATDB sphere table is reproduced bit for bit) and
    reaches the ATDB disc entry at the flat-face limit: HTG's continuum database is modified Newtonian, its
    sphere/disc ratio being 0.49897 at every Mach number, and the same integral over a meshed flat plate gives 2.00."""
    import json
    from reentry_model import aero, mesh as mesh_mod
    from helpers import REPO_ROOT
    b = melting_body(layered_mesh)
    assert b.drag_shape_factor() == pytest.approx(1.0, abs=1e-12)
    cd, area = b.newtonian_drag()
    assert cd == pytest.approx(0.92, rel=0.02) and area == pytest.approx(math.pi * 0.05 ** 2, rel=0.02)
    sphere = json.load(open(os.path.join(REPO_ROOT, "reentry_model", "data", "atdb_sphere.json")))
    disc = json.load(open(os.path.join(REPO_ROOT, "reentry_model", "data", "atdb_disc.json")))
    ratios = np.array(sphere["cd_continuum"]) / np.array(disc["cd_continuum"])
    assert np.allclose(ratios, 0.49897, atol=1e-5) and disc["mach"] == sphere["mach"]
    assert np.allclose(np.array(disc["cd_free_molecular"]) / np.array(sphere["cd_free_molecular"]), 1.03, atol=0.04)
    # a flat plate normal to the flow: every windward element is square-on, so the integral is C_p,max itself
    plate = mesh_mod.box_mesh(0.004, 0.06, 0.06, 0.002)
    flat = body.MeltingBody(plate, material.Material.from_drama_json("AA7075_range"), thermal.thermal_solver("skfem"),
                            1.0, settings=body.MeltSettings())
    flat_cd = flat.newtonian_drag()[0]
    assert flat_cd == pytest.approx(body.CP_MAX_NEWTONIAN, rel=1e-6)
    assert flat_cd / b.newtonian_drag_sphere == pytest.approx(1.0 / 0.49897, rel=0.02)     # 2.00 vs ATDB's 2.0041
    fixed = melting_body(layered_mesh, size_feedback="initial")
    assert fixed.drag_shape_factor() == 1.0 and fixed.reference_area() == pytest.approx(math.pi * 0.05 ** 2, rel=2e-3)
    # the trajectory multiplies the continuum entry only
    tables, bridging = aero.SphereDragTables.from_json(), aero.SesamTable()
    assert aero.drag_coefficient(0.0, 10.0, tables, bridging, 2.0) == pytest.approx(2.0 * tables.cd_continuum(10.0))
    assert aero.drag_coefficient(1e6, 10.0, tables, bridging, 2.0) == pytest.approx(tables.cd_free_molecular(10.0))


def test_film_temperature_freeze_back_and_the_netted_transfer(layered_mesh):
    """2.2 MW/m2 for 6 s on a 100 mm body, then the heating off for 60 s, with no flow (so nothing runs off or
    sprays): the film takes the surface's temperature and its enthalpy is the mean of the nodes' -- not the enthalpy
    of the mean temperature, which differs by the whole latent heat on the ramp -- it freezes back onto its owner
    element once the surface falls through the feed ramp, and feed and freeze are netted, so neither runs while the
    other does. Mass and the coupled energy balance are exact throughout."""
    import types
    b = melting_body(layered_mesh, name="AA7075")
    b.solver.set_temperature(800.0)                                    # just below AA7075's melting point
    b.energy0 = b.energy()
    loads = lambda q: types.SimpleNamespace(q_conv=np.full(b.surface.n_patches, q))
    peak, tail = 0.0, 0.0
    for i in range(132):
        b.advance(i * 0.5, 0.5, loads(2.2e6 if i < 16 else 0.0))
        peak = max(peak, b.m_f.sum())
        if i >= 112:                                                   # the last 20 steps: nothing is being heated
            tail += b.last_melt["feed_mass"]
        assert abs(b.energy_balance_residual()) < 1e-6
        assert b.mass(0.0) + b.removed_mass == pytest.approx(b.mass0, rel=1e-12)
    assert tail / 20.0 < 1e-3 * peak   # netted: the feed decays away instead of cycling the film every step (run as
    #                                    two independent rates it sat at 3 % of the body per step with the film static)
    T, mat = b.solver.temperature(), b.material
    assert b.film_enthalpy() == pytest.approx(mat.enthalpy_liquid(T)[b.surface.faces].mean(axis=1))
    mix = mat.enthalpy(T)[b.surface.faces].mean(axis=1)                 # the film carries its latent heat wherever it
    assert (b.film_enthalpy() >= mix - 1e-6).all()                      # sits, so it is never below the mixture's
    assert (b.film_enthalpy() > mix + 0.2 * mat.latent_heat).any()      # and is well above it where the surface cooled
    assert b.film_energy() == pytest.approx((b.m_f * b.film_enthalpy()).sum())
    assert b.film_temperature() == pytest.approx(b.surface_temperature())
    assert peak > 0.0 and b.frozen_mass > 0.0 and b.phi.max() <= 1.0     # the film came back and phi never overflows
    assert b.solver.film_mass.sum() == pytest.approx(b.m_f.sum())        # the solver carries the same film
    assert b.melt_stats()["film_T_max_K"] == pytest.approx(b.surface_temperature().max())


def test_a_deferred_melt_load_never_moves_a_node_more_than_the_cap(layered_mesh):
    """A drained node has nothing to heat with the enthalpy its melt delivered: no node is moved more than
    LOAD_DT_MAX in one step, the remainder waits in `pending_load`, and the balance counts it either way."""
    import types
    b = melting_body(layered_mesh, name="AA7075")
    loads = lambda q: types.SimpleNamespace(q_conv=np.full(b.surface.n_patches, q))
    b.advance(0.0, 0.5, loads(0.0))
    b.pending_load[:] = 0.0
    b.pending_load[np.unique(b.surface.faces)] = 1.0e4                   # 10 kJ on every boundary node
    queued, T0 = b.pending_load.sum(), b.solver.temperature().copy()
    room = body.LOAD_DT_MAX * b.solver.nodal_capacity()
    b.advance(0.5, 0.5, loads(0.0))
    assert b.pending_load.sum() == pytest.approx(queued - np.minimum(1.0e4, room[np.unique(b.surface.faces)]).sum())
    assert (b.solver.temperature() - T0).max() < 1.05 * body.LOAD_DT_MAX
    assert b.pending_load.sum() + b.total_applied_load == pytest.approx(queued)     # nothing is lost on the way


def test_a_dying_element_hands_its_film_to_the_faces_it_exposes(layered_mesh):
    """When an element is removed the surface recedes into the film it carried, so that film belongs on the faces of
    that same element which have just become boundary -- its inward face and the walls of the pit it opens -- not on
    whichever surviving centroid happens to be nearest (which is a statement about the mesh's numbering, and which
    once put a third of a collapsing body's film onto one 0.6 mm2 facet)."""
    b = melting_body(layered_mesh, name="AA7075")
    owner = b.surface.owner[0]                                      # kill one patch owner, with film only on its patch
    b.m_f[:] = 0.0
    b.m_f[b.surface.owner == owner] = 1.0e-4
    carried, faces = float(b.m_f.sum()), b.mesh.element_faces([owner])[0]
    gone_normal = b.surface.normals[b.surface.owner == owner][0]
    b._kill(np.array([owner]))
    exposed = np.array([q for q in b.patch_of_face[faces] if q >= 0])
    assert exposed.size and b.m_f.sum() == pytest.approx(carried)    # nothing lost, and it all landed on that element's
    assert b.m_f[exposed].sum() == pytest.approx(carried)            # own newly exposed faces
    # shared by the area each exposed face presents to the vanished patch, i.e. its area projected on that patch's
    # outward normal: the floor the recession uncovered takes the film, the pit walls standing perpendicular take none
    cosine = np.clip(b.surface.normals[exposed] @ gone_normal, 0.0, None)
    w = b.surface.areas[exposed] * cosine
    assert b.m_f[exposed] == pytest.approx(carried * w / w.sum())
    floor = exposed[int(np.argmax(cosine))]                         # the face most nearly parallel to the old patch
    assert b.m_f[floor] > 0.5 * carried and cosine.max() > 0.9      # takes the bulk of it
```


- [ ] **Step 2: Run the tests to see them fail**

Run: `"$PY" -m pytest tests/test_reentry_model_melting.py -q`
Expected: AttributeError (`body.MeltSettings`).

- [ ] **Step 3: Replace `reentry_model/body.py`**

```python
"""The body the trajectory carries, as a stepper (spec section 4). ConstantBody keeps mass and one temperature
(Step 1 behaviour, `--thermal none`); ThermalBody wraps a mesh, a material and a conduction solver and advances the
temperature field with each macro step's per-patch convective loads, keeping the heat bookkeeping SESAM reports."""
import math
from dataclasses import dataclass
from typing import Protocol

import numpy as np


class Body(Protocol):
    def mass(self, t) -> float: ...
    def advance(self, t, dt, loads) -> None: ...          # loads: heating.HeatingResult applied over [t - dt, t]
    def reference_area(self): ...                        # m2 drag reference area, or None for the fixed pi D^2/4
    def reference_length(self): ...                      # m length of the body Knudsen number, or None for the fixed D
    def drag_shape_factor(self): ...                     # continuum C_D of the current shape / a sphere's (1 for a sphere)
    def surface_temperature(self) -> np.ndarray: ...     # K per patch
    def mean_temperature(self) -> float: ...             # K, energy-equivalent (spec 6.4)
    def energy(self) -> float: ...                       # J stored above the material's reference temperature
    def field(self): ...                                 # nodal temperatures, or None


def sphere_mass(diameter_m, density_kgm3):
    """Mass of a solid sphere [kg]."""
    return density_kgm3 * 4.0 / 3.0 * math.pi * (diameter_m / 2.0) ** 3


@dataclass
class ConstantBody:
    mass_kg: float
    temperature_K: float = 300.0

    def mass(self, t):
        return self.mass_kg

    def reference_area(self):
        return None

    def reference_length(self):
        return None

    def drag_shape_factor(self):
        return 1.0

    def advance(self, t, dt, loads):
        return None

    def surface_temperature(self):
        return np.array([self.temperature_K])

    def mean_temperature(self):
        return self.temperature_K

    def energy(self):
        return 0.0

    def field(self):
        return None


class ThermalBody:
    """Mesh + material + ThermalSolver. `v_hat` is the direction of motion in the body frame (fixed attitude: the
    stagnation patch is the one whose normal is along v_hat, spec 13.1)."""

    def __init__(self, mesh, material, solver, mass_kg, T0=300.0, emissivity=None, T_ambient=0.0, v_hat=(1.0, 0.0, 0.0)):
        self.mesh, self.material, self.solver, self.mass_kg = mesh, material, solver, mass_kg
        self.emissivity = material.emissivity if emissivity is None else float(emissivity)
        self.T_ambient = float(T_ambient)
        self.surface = mesh.surface()
        self.theta = self.surface.angles_to(v_hat)
        self.i_stag, self.i_back = self.surface.patch_toward(v_hat), self.surface.patch_toward(-np.asarray(v_hat, dtype=float))
        self.i_centre = mesh.centre_node()
        self.volume = mesh.volume()
        self.solver.setup(mesh, material, self.emissivity)
        self.solver.set_temperature(T0)
        self.energy0 = self.solver.energy()
        self.integrated_heat = 0.0       # J, integral of Q_conv dt (SESAM's integrated_heat_J)
        self.absorbed_heat = 0.0         # J, integral of (Q_conv - Q_rad) dt
        self.radiated_heat = 0.0         # J
        self.iterations = []             # Newton iterations per step
        self.last = None                 # thermal.StepResult of the last step

    def mass(self, t):
        return self.mass_kg

    def reference_area(self):
        return None

    def reference_length(self):
        """Length scale of the body Knudsen number, or None for the fixed initial diameter."""
        return None

    def drag_shape_factor(self):
        """Continuum drag of the body's shape relative to a sphere's: 1 for the sphere of Steps 1-2."""
        return 1.0

    def nose_radius(self):
        """Stagnation-point radius of curvature [m] for the heating and the surface flow: the sphere's."""
        return float(self.mesh.params.get("radius_m", 0.05))

    def advance(self, t, dt, loads, state=None):
        res = self.solver.step(dt, loads.q_conv, self.T_ambient)
        self.integrated_heat += res.Q_conv * dt
        self.absorbed_heat += (res.Q_conv - res.Q_rad) * dt
        self.radiated_heat += res.Q_rad * dt
        self.iterations.append(res.iterations)
        self.last = res

    def surface_temperature(self):
        return self.surface.facet_mean(self.solver.temperature())

    def mean_temperature(self):
        return float(self.material.temperature_from_enthalpy(self.energy() / (self.material.rho * self.volume)))

    def energy(self):
        return self.solver.energy()

    def field(self):
        return self.solver.temperature()

    def radiated_power(self):
        return self.solver.radiated_power(self.T_ambient)

    def energy_balance_residual(self):
        """(E - E0 - absorbed heat) / absorbed heat: zero for an exact discrete balance."""
        return (self.energy() - self.energy0 - self.absorbed_heat) / self.absorbed_heat if self.absorbed_heat else 0.0

    def surface_stats(self):
        T, Tf = self.solver.temperature(), self.surface_temperature()
        return {"surface_T_max_K": float(Tf.max()), "surface_T_min_K": float(Tf.min()),
                "surface_T_mean_K": float((Tf * self.surface.areas).sum() / self.surface.area),
                "T_stagnation_K": float(Tf[self.i_stag]), "T_back_K": float(Tf[self.i_back]), "T_centre_K": float(T[self.i_centre])}


SIZE_FEEDBACK_NAMES = ("current", "initial")
FEED_DEPTH_NAMES = ("conjugate", "all")     # how deep melt may leave an element for the film (Step 3, 2026-09-22)
CP_MAX_NEWTONIAN = 1.84          # only the ratio to the meshed sphere's own value is used, so this cancels
NOSE_CAP_ANGLE = 30.0            # deg: the windward cap fitted for the nose radius (depth (1 - cos 30 deg) R_t behind the front)
NOSE_CAP_FACTOR = 1.67           # a flat face of radius R_t heats like a sphere of 1.67 R_t (its stagnation velocity gradient is
                                 # ~0.6 x a sphere's of the same radius, Boison & Curtiss 1959): the cap on the fitted radius
PHI_MIN = 1.0e-3                 # element fraction kept by elements that do not own a patch (they cannot die: no cavities)
NEAREST_PATCHES = 4              # patches that receive an interior element's liquid (area-weighted)
LOAD_DT_MAX = 1000.0             # K, the most a deferred melt load may move one node in one macro step
PHI_DEATH = 0.05                 # a patch owner below this fraction dies (its remainder goes to the film): keeps the surface
                                 # nodes' thermal mass above 5 % of an element's, which the Newton iteration needs (measured
                                 # 2026-09-20: owners at 1e-3 made surface nodes swing by hundreds of K between iterates)
REMOVAL_NAMES = ("girin", "instant")


@dataclass
class MeltSettings:
    removal: str = "girin"           # girin: film + runoff + spraying; instant: liquid removed as it forms (verification device)
    runoff: bool = True
    demise_fraction: float = 0.01    # the run ends when the mass falls below this fraction of the initial mass
    particles: bool = True           # keep the source-table rows
    size_feedback: str = "current"    # current: Kn on the equivalent diameter of the remaining mass and the nose radius fitted
                                      # to the windward cap; initial: D0 and R0 throughout (SESAM's convention, for the devices)
    feed_depth: str = "conjugate"     # conjugate: only liquid the gas shear reaches leaves its element (fact 28); all: any depth

    def __post_init__(self):
        if self.removal not in REMOVAL_NAMES:
            raise ValueError("removal must be one of {}, got {!r}".format(REMOVAL_NAMES, self.removal))
        if self.size_feedback not in SIZE_FEEDBACK_NAMES:
            raise ValueError("size_feedback must be one of {}, got {!r}".format(SIZE_FEEDBACK_NAMES, self.size_feedback))
        if self.feed_depth not in FEED_DEPTH_NAMES:
            raise ValueError("feed_depth must be one of {}, got {!r}".format(FEED_DEPTH_NAMES, self.feed_depth))
        if not 0.0 < self.demise_fraction < 1.0:
            raise ValueError("demise_fraction must be within (0, 1)")


def fit_sphere(points):
    """Algebraic least-squares sphere through `points` (n, 3): (centre, radius)."""
    x = np.asarray(points, dtype=float)
    A = np.column_stack([2.0 * x, np.ones(len(x))])
    p = np.linalg.lstsq(A, (x * x).sum(axis=1), rcond=None)[0]
    c = p[:3]
    return c, float(math.sqrt(max(p[3] + c @ c, 0.0)))


class MeltingBody(ThermalBody):
    """ThermalBody with melting (spec Step 3 sections 6, 8-11): element fractions phi_e, the film per patch, the
    gas-side surface flow, spraying, element death and the mass/area/energy accounting.

    Per macro step (`advance`): the conduction step with the deferred melt loads of the previous step -> melt step:
    (i) every element's liquid inventory f_feed(T_e) phi_e rho V_e becomes film on its patches (owners) or the
    nearest patch (interior elements keep PHI_MIN so no cavity can open), netted against (iv) below; what leaves is
    the element's molten part, at the enthalpy that part carries, so the melt front moves at the energy-limited rate
    whatever the element size; (ii) surface flow, delta_m, lubrication, runoff transport; (iii) spraying and release;
    (iv) re-solidification: the fraction 1 - f_feed(T_patch) of each patch's film returns to its owner element
    (raising phi_e, capped at 1) once the patch falls back through the ramp -- the mirror of the feed rule, netted
    with it so that no element both melts and freezes in one step; (v) patch owners at phi <= PHI_DEATH die: the
    mesh's active set, the surface, the film (handed to the nearest surviving patch) and the solver's fractions are
    refreshed. With `removal = "instant"` the liquid leaves the body at h_liquid instead, and the difference to the
    element's own enthalpy is a nodal load on its nodes over the next step.

    The film has the surface's own temperature, because it is thermally thin: q b / k_l = 0.22 K across a 10 um film
    at 2 MW/m2 (22 K even at 1 mm) and b^2/alpha = 0.3 ms against the 0.5 s macro step. Rather than give it an energy
    equation, its *mass* is handed to the solver (`set_film_mass`), which weighs it with the liquid heat capacity, so
    the film rides the boundary nodes at the surface temperature by construction; it heats, cools and re-solidifies
    with the surface, and the droplets carry away whatever superheat the surface has. A film is liquid by
    construction, so it holds the *liquid* enthalpy h(T) + L_f (1 - f_l) wherever it sits (`film_enthalpy`): feeding
    it costs the latent heat its mass has not yet paid, which is what keeps a body that merely sits on the melting
    ramp from turning into film for free. `film_energy` sums that enthalpy over the patches -- identical to the nodal
    sum the solver's capacity matrix carries -- and every transfer -- feed, freeze-back, runoff between patches at
    different temperatures, the hand-over when a patch dies -- books the enthalpy difference it carries, at the
    patch it arrives on (`melt_step`), so the coupled balance stays exact (measured 1e-10). `removal = "instant"`
    removes the liquid inventory of every element as it forms, without film, runoff or spraying: the lumped-melting
    device that reproduces SESAM's Q/L_f law (SESAM hollows the sphere at fixed outer geometry, measured
    2026-09-20, so the geometry is kept until elements die).
    mass(t) = sum phi rho V + sum m_f; the drag reference area is the current surface's projection on the flight
    direction (pi R^2 while intact). Size feedback (decided 2026-09-21, replacing spec section 10's fixed R0): the body
    Knudsen number uses the equivalent diameter of the remaining mass, and the stagnation radius for the heating and
    the surface flow is fitted to the current windward cap -- a least-squares sphere through the patch centroids within
    (1 - cos 30 deg) R_t of the front-most point (R_t the transverse radius about the mass centre), bounded to
    [0.1 R_t, NOSE_CAP_FACTOR R_t] (the front erodes fastest and flattens, which lowers the stagnation heating as R^-1/2;
    a flat face heats like a sphere of 1.67 x its radius). The sphere drag tables are kept. `size_feedback = "initial"`
    keeps D0 and R0 (SESAM's convention, used by the verification devices)."""

    def __init__(self, mesh, material, solver, mass_kg, flow=None, spray_model=None, settings=None, T0=300.0,
                 emissivity=None, T_ambient=0.0, v_hat=(1.0, 0.0, 0.0)):
        if not material.melts or material.liquid is None:
            raise ValueError("MeltingBody needs a material with a latent heat and liquid properties (AA7075 or AA7075_range)")
        super().__init__(mesh, material, solver, mass_kg, T0, emissivity, T_ambient, v_hat)
        from . import film as film_mod, spray as spray_mod, surface_flow
        self.settings = settings or MeltSettings()
        self.flow = flow or surface_flow.SurfaceFlow()
        self.spray = spray_model or spray_mod.SprayModel(material.liquid)
        self._film_mod = film_mod
        self.v_hat = np.asarray(v_hat, dtype=float) / np.linalg.norm(v_hat)
        self.liquid = material.liquid
        self.element_mass = material.rho * mesh.element_volumes()          # kg at phi = 1
        self.phi = np.ones(mesh.n_elements)
        self.mass0 = float(self.element_mass.sum())
        self.m_f = np.zeros(self.surface.n_patches)
        self.pending_load = np.zeros(mesh.n_nodes)                          # J, deferred melt energy for the next step
        self.sprayed_mass = self.runoff_mass = self.removed_mass = self.removed_enthalpy = self.frozen_mass = 0.0
        self.n_released = 0.0
        self.source_rows = []
        self.hist_n, self.hist_m = np.zeros(spray_mod.N_BINS), np.zeros(spray_mod.N_BINS)
        self.last_flow = self.last_spray = None
        self.last_melt = {"n_dead": 0, "runoff_substeps": 0, "released_mass": 0.0, "n_released": 0.0, "feed_mass": 0.0}
        self.melt_onset = self.spray_onset = None
        self.consumed = False
        self._refresh_geometry()
        self.newtonian_drag_sphere = self.newtonian_drag()[0]      # the meshed sphere's own value: the normalisation

    # -- geometry -----------------------------------------------------------------------------------------------
    def _refresh_geometry(self):
        self.surface = self.mesh.surface()
        self.theta = self.surface.angles_to(self.v_hat)
        self.t_hat = self.surface.tangent_from(self.v_hat)
        self.windward = self.theta <= 0.5 * np.pi
        self.i_stag, self.i_back = self.surface.patch_toward(self.v_hat), self.surface.patch_toward(-self.v_hat)
        self.runoff = self._film_mod.Runoff(self.surface, self.mesh.points, self.windward) if self.settings.runoff else None
        owners = self.surface.owner
        self.owner_area = np.bincount(owners, self.surface.areas, self.mesh.n_elements)      # patch area per owner element
        from scipy.spatial import cKDTree
        self._patch_tree = cKDTree(self.surface.centroids)
        self.patch_of_face = np.full(len(self.mesh._face_nodes), -1, dtype=np.int64)     # face id -> patch index (-1: not a patch)
        self.patch_of_face[self.surface.face_ids] = np.arange(self.surface.n_patches)
        self._fit_nose()

    def _fit_nose(self):
        """Mass centre, transverse radius and the windward-cap radius of the current surface (module docstring)."""
        m, s, v = self.mesh, self.surface, self.v_hat
        w = self.phi * m.element_volumes()
        self.mass_centre = (w[:, None] * m.points[m.tets].mean(axis=1)).sum(axis=0) / max(w.sum(), 1e-300)
        X = m.points[np.unique(s.faces)] - self.mass_centre
        self.transverse_radius = float(np.sqrt(np.maximum(np.linalg.norm(X, axis=1) ** 2 - (X @ v) ** 2, 0.0)).max())
        x = (s.centroids - self.mass_centre) @ v
        band = x >= x.max() - (1.0 - math.cos(math.radians(NOSE_CAP_ANGLE))) * self.transverse_radius
        r0 = float(m.params.get("radius_m", 0.05))
        if band.sum() >= 4 and self.transverse_radius > 0.0:
            _, r = fit_sphere(s.centroids[band])
            self.fitted_nose_radius = r
            self.cap_nose_radius = float(min(max(r, 0.1 * self.transverse_radius), NOSE_CAP_FACTOR * self.transverse_radius))
        else:
            self.fitted_nose_radius = self.cap_nose_radius = r0

    def nose_radius(self):
        if self.settings.size_feedback == "initial":
            return float(self.mesh.params.get("radius_m", 0.05))
        return self.cap_nose_radius

    def newtonian_drag(self):
        """Modified-Newtonian drag coefficient of the current windward silhouette: sum C_p,max (n.v)^3 A / sum (n.v) A
        over the convex hull of the surface (0.92 for a sphere, 1.84 for a flat face, with C_p,max = 1.84).

        The hull, not the raw facets: the staircase left by element death scatters the normals and biases the integral
        low and non-monotonically (0.79-1.07 against the hull's smooth 0.92-1.73, measured 2026-09-21), and the flow
        sees the silhouette -- a shallow cavity recovers roughly the stagnation pressure at its mouth. It is therefore
        an upper bound where the face is cratered."""
        from scipy.spatial import ConvexHull
        points = self.mesh.points[np.unique(self.surface.faces)]
        try:
            hull = ConvexHull(points)
        except Exception:                                        # degenerate remnant: fall back to the sphere
            return CP_MAX_NEWTONIAN / 2.0, self.surface.projected_area(self.v_hat)
        tri = points[hull.simplices]
        n = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
        area = 0.5 * np.linalg.norm(n, axis=1)
        n = n / (2.0 * np.maximum(area, 1e-300))[:, None]
        inward = np.einsum("ij,ij->i", n, tri.mean(axis=1) - points.mean(axis=0)) < 0.0
        n[inward] *= -1.0
        c = n @ self.v_hat
        w = c > 0.0
        projected = float((area[w] * c[w]).sum())
        if projected <= 0.0:
            return CP_MAX_NEWTONIAN / 2.0, self.surface.projected_area(self.v_hat)
        return float((CP_MAX_NEWTONIAN * c[w] ** 3 * area[w]).sum() / projected), projected

    def drag_shape_factor(self):
        """The continuum drag of the current shape relative to a sphere's, for aero.drag_coefficient: exactly 1 while
        the body is intact (it is normalised by the meshed sphere's own Newtonian value, so the mesh discretisation
        error cancels and the ATDB sphere table is reproduced bit for bit) and 2.00 at the flat-face limit, where the
        ATDB disc entry sits. `size_feedback = "initial"` pins it to 1 (SESAM's convention)."""
        if self.settings.size_feedback == "initial":
            return 1.0
        return self.newtonian_drag()[0] / self.newtonian_drag_sphere

    def reference_length(self):
        if self.settings.size_feedback == "initial":
            return None
        return 2.0 * self.equivalent_radius()

    def mass(self, t):
        return float((self.phi * self.element_mass).sum() + self.m_f.sum())

    def reference_area(self):
        """Drag reference area: the hull's projected area, paired with the C_D of drag_shape_factor (the raw surface's
        projection counts forward-facing patches inside a crater and is 2-4 % larger)."""
        if self.settings.size_feedback == "initial":
            return self.surface.projected_area(self.v_hat)
        return self.newtonian_drag()[1]

    def equivalent_radius(self):
        return (3.0 * self.mass(0.0) / (4.0 * math.pi * self.material.rho)) ** (1.0 / 3.0)

    def energy(self):
        return self.solver.energy() + self.film_energy()

    def mean_temperature(self):
        m = self.mass(0.0)
        return float(self.material.temperature_from_enthalpy(self.energy() / m)) if m > 0.0 else 0.0

    def melt_front_depth(self):
        """Deepest molten point [m]: the distance from the original surface of the innermost element with f_l > 0."""
        f_l = self.material.liquid_fraction(self.solver.temperature())[self.mesh.tets].max(axis=1)
        molten = self.mesh.active & (f_l > 0.0)
        if not molten.any():
            return 0.0
        r = np.linalg.norm(self.mesh.points[self.mesh.tets[molten]].mean(axis=1), axis=1)
        return float(self.mesh.params.get("radius_m", r.max()) - r.min())

    # -- the step --------------------------------------------------------------------------------------------------
    def advance(self, t, dt, loads, state=None):
        # A deferred melt load is energy the transferred mass delivered to the nodes it arrived at; a node that has
        # since melted away has nothing to heat with it, so no node is asked to move more than LOAD_DT_MAX in one
        # step and the remainder waits (it is applied as soon as the node has the capacity, dropped into Q_dropped
        # if the node dies first, and counted either way -- the balance sees `pending_load` and `Q_dropped` alike).
        room = LOAD_DT_MAX * self.solver.nodal_capacity()
        applied = np.clip(self.pending_load, -room, room)
        self.pending_load = self.pending_load - applied
        res = self.solver.step(dt, loads.q_conv, self.T_ambient, nodal_load=applied / dt)
        self.integrated_heat += res.Q_conv * dt
        self.absorbed_heat += (res.Q_conv - res.Q_rad) * dt
        self.radiated_heat += res.Q_rad * dt
        self.iterations.append(res.iterations)
        self.last = res
        self._applied_total = getattr(self, "_applied_total", 0.0) + res.Q_extra * dt
        self._dropped_total = getattr(self, "_dropped_total", 0.0) + res.Q_dropped * dt
        self.melt_step(t, dt, state)

    def melt_step(self, t, dt, state):
        """Feed, film transport, spraying, re-solidification and element death (class docstring).

        Energy is moved between three places -- the finite-element solid, the film, and the outside world -- and every
        transfer is booked by the same rule: mass that moves at one temperature books nothing, and mass that arrives
        somewhere colder or hotter than it left books the difference it carries, as a deferred load on the nodes it
        arrives at. So a transfer of m kilograms from enthalpy h_src to h_dst changes the accounted energy by
        m (h_dst - h_src) and applies m (h_src - h_dst) to the destination: the balance closes exactly, and nothing
        larger than the difference ever touches the temperature field. Booking the two halves separately -- the
        solid's m h_src on its element's nodes, the film's m h_dst on its patch's nodes -- closes the balance just as
        exactly and wrecks the field: the same mass is spread by 1/4 over four nodes on one side and by 1/3 over
        three on the other, which left +-m h/12 on every surface node and swung the body to 1618 K and -1202 K in two
        macro steps (measured 2026-09-22). The enthalpies are the ones the two energy functionals actually use: the
        element's is the mean of h(T_i) over its four nodes, the film's the mean of the *liquid* h over its patch's
        three, and what leaves an element is its molten part, at the enthalpy that part carries."""
        mat, liq, s = self.material, self.liquid, self.settings
        from . import spray as spray_mod
        T = self.solver.temperature()
        tets = self.mesh.tets
        # The surface flow is evaluated once per step, here rather than inside the film step, because the feed needs
        # Girin's conjugate melt-layer depth: the gas shear penetrates the liquid only that far, so only liquid within
        # that depth of the wall can be carried away, and material deeper than it keeps its melt until the surface has
        # receded to it. Feeding from any depth let the interior drain through an intact skin -- by 82 s of the 100 mm
        # flight the core was being consumed faster than the outermost shell (measured 2026-09-22).
        flow = delta_m = None
        if state is not None and s.removal != "instant" and self.surface.n_patches and self.m_f.size:
            flow = self.flow.evaluate(state, self.theta, self.nose_radius(), liq.rho, self.surface_temperature())
            delta_m, _ = spray_mod.melt_layer(flow, liq)
        h_node = mat.enthalpy(T)
        h_e = h_node[tets].mean(axis=1)                          # as solver.energy() weighs the solid
        h_liq = mat.enthalpy_liquid(T)
        h_p = self.film_enthalpy(h_liq)                          # as film_energy() weighs the film: liquid, with L_f
        # (i) feed and (iv) re-solidification, netted. The feed hands the molten fraction of what each element still
        # holds to the film (f_feed of phi_e, so an element that is a quarter molten hands over a quarter of its
        # remainder); the mirror rule hands back the fraction of the film that has fallen below the ramp. What stops
        # the feed from eating a body that merely sits on the ramp is not the rule but the energy: film is liquid and
        # carries L_f wherever it sits (`Material.enthalpy_liquid`), so every kilogram fed debits the latent heat it
        # has not paid and the surface falls back to the solidus. Book the film at the mixture enthalpy instead and
        # melting is free on the ramp -- a 100 mm sphere turned entirely to film on a quarter of its latent heat
        # (measured 2026-09-22). The two directions are netted per element because they are one equilibrium seen
        # from opposite sides: run separately they cycled 3 % of the body's mass through the film every step with no
        # net effect, pinning the surface at T_feed and paying Newton iterations for it.
        fn = mat.feed_fraction(T)[tets]
        f = fn.mean(axis=1) * self.mesh.active
        if s.feed_depth == "conjugate" and s.removal != "instant":
            f = f * self._shear_reaches(delta_m)
        h_hot = (fn * h_node[tets]).mean(axis=1) * self.mesh.active     # what the molten part carries, per kg of element
        cap = np.where(self.owner_area > 0.0, self.phi, np.maximum(self.phi - PHI_MIN, 0.0))
        gross = np.minimum(f * self.phi, cap) * self.element_mass
        solid = (1.0 - mat.feed_fraction(self.film_temperature())) * self.m_f if s.removal != "instant" else np.zeros(0)
        want = np.bincount(self.surface.owner, solid, self.mesh.n_elements) if solid.size else np.zeros(self.mesh.n_elements)
        net = gross - want
        fed, want = np.maximum(net, 0.0), np.maximum(-net, 0.0)
        # the enthalpy the fed mass carries: the *molten* part of the element, not its mean. The mass that leaves sits
        # at the nodes that are above the ramp, so it carries h(T_i) there -- and what stays behind is the colder
        # rest, which is why the element must be debited the difference. Debit only the destination and the element
        # keeps a melt fraction it no longer has and feeds it again next step: the surface then melted three times
        # faster than the heat allowed and the body was gone in ten steps (measured 2026-09-22).
        with np.errstate(divide="ignore", invalid="ignore"):
            fed_h = np.where(f > 0.0, fed * h_hot / np.where(f > 0.0, f, 1.0), 0.0)
        self.phi = self.phi - fed / self.element_mass
        feed_mass = float(fed.sum())
        if feed_mass > 0.0 and self.melt_onset is None:
            self.melt_onset = t
        if s.removal == "instant":
            self.removed_mass += feed_mass
            self.removed_enthalpy += feed_mass * mat.h_liquid
            self._defer_to_elements(fed * (h_e - mat.h_liquid))  # it leaves the body at the liquidus, not at h_e
            released = feed_mass
        else:
            self._defer_to_elements(fed * h_e - fed_h)           # the element keeps only what the melt left behind
            delta, carried = self._add_to_film(fed, fed_h)
            self._defer_to_patches(carried - delta * h_p)        # ... and the melt arrives at its patch's temperature
            released = self._film_and_spray(t, dt, state, h_p, flow, delta_m)
        frozen = 0.0 if s.removal == "instant" else self._freeze_back(want, solid, h_e, h_p)
        # (v) death of consumed patch owners; a death exposes its neighbours, which die in turn if they are consumed
        n_dead = 0
        while not self.consumed:
            dead = np.flatnonzero((self.owner_area > 0.0) & self.mesh.active & (self.phi <= PHI_DEATH))
            if dead.size == 0:
                break
            rest = self.phi[dead] * self.element_mass[dead]
            if s.removal == "instant":
                self.removed_mass += float(rest.sum())
                self.removed_enthalpy += float(rest.sum()) * mat.h_liquid
                self._defer_to_elements(rest * (h_e[dead] - mat.h_liquid), dead)
            else:
                extra, extra_h = np.zeros(self.mesh.n_elements), np.zeros(self.mesh.n_elements)
                extra[dead], extra_h[dead] = rest, rest * h_e[dead]
                delta, carried = self._add_to_film(extra, extra_h)
                self._defer_to_patches(carried - delta * h_p)
            self.phi[dead] = 0.0
            before = float((self.m_f * h_p).sum())
            self._kill(dead)                                     # the film of a dead patch moves to another patch...
            if self.m_f.size:                                    # ... which is at its own temperature
                h_p = self.film_enthalpy(h_liq)
                self._spread_to_patches(before - float((self.m_f * h_p).sum()), self.m_f)
            n_dead += int(dead.size)
        self.solver.set_fractions(self.phi)
        self.solver.set_film_mass(self._film_nodal())
        self.frozen_mass += frozen
        self.last_melt.update({"n_dead": n_dead, "released_mass": released, "feed_mass": feed_mass, "frozen_mass": frozen})

    def _freeze_back(self, want, solid, h_e, h_p):
        """Give each element back the `want` kilograms of film that have fallen below the feed ramp, drawn from its
        own patches in proportion to what each wants to give (`solid`). An element can only recover film that is
        still there -- what ran off or sprayed away is gone -- and phi_e is capped at 1, because a patch may hold
        what its neighbours' runoff delivered and the mesh cannot grow a layer outside itself, so film with nowhere
        to go simply stays film (`film_frozen_fraction` records how much)."""
        if not self.m_f.size or not want.any():
            return 0.0
        owner = self.surface.owner
        have = np.bincount(owner, solid, self.mesh.n_elements)
        room = np.maximum(1.0 - self.phi, 0.0) * self.element_mass
        with np.errstate(divide="ignore", invalid="ignore"):
            share = np.where(have > 0.0, np.minimum(want, room) / np.where(have > 0.0, have, 1.0), 0.0)
        taken = np.minimum(solid * share[owner], self.m_f)
        if not taken.any():
            return 0.0
        per_element = np.bincount(owner, taken, self.mesh.n_elements)
        self.phi = np.clip(self.phi + per_element / self.element_mass, 0.0, 1.0)
        self.m_f = self.m_f - taken
        # the film arrives in the element it froze onto, which is colder than the surface it left
        self._defer_to_elements(np.bincount(owner, taken * h_p, self.mesh.n_elements) - per_element * h_e)
        return float(taken.sum())

    def _defer_to_elements(self, energy, elements=None):
        """Book a per-element energy [J] on the elements' nodes, to be applied over the next step."""
        e = np.asarray(energy, dtype=float)
        if not e.any():
            return
        nodes = self.mesh.tets if elements is None else self.mesh.tets[elements]
        np.add.at(self.pending_load, nodes.ravel(), np.repeat(e / 4.0, 4))

    def _spread_to_patches(self, energy, weights):
        """Book one energy [J] over the patches that received mass, in proportion to how much each received. Used
        where the transfer's per-edge detail is not carried back (runoff, and the hand-over when a patch dies): the
        total is exact and it lands on the arriving liquid, which is where the difference is released."""
        total = float(weights.sum())
        if total > 0.0 and energy != 0.0:
            self._defer_to_patches(energy * weights / total)

    def _defer_to_patches(self, energy):
        """Book a per-patch energy [J] on the patches' nodes, to be applied over the next step."""
        e = np.asarray(energy, dtype=float)
        if e.any():
            np.add.at(self.pending_load, self.surface.faces.ravel(), np.repeat(e / 3.0, 3))

    def _film_nodal(self):
        """The film's mass per node [kg]: each patch's film shared over its three nodes. The film's thermal state
        lives on the nodes because that is where the solver's capacity and the solver's enthalpy live."""
        if not self.m_f.size:
            return np.zeros(self.mesh.n_nodes)
        return np.bincount(self.surface.faces.ravel(), np.repeat(self.m_f / 3.0, 3), self.mesh.n_nodes)



    def _add_to_film(self, fed, carried):
        """Distribute each element's freed liquid to the film: owners to their own patches by area, interior elements
        to the NEAREST_PATCHES nearest patches by area. `carried` is the enthalpy that liquid takes with it [J per
        element]; it rides the same weights, so the caller knows what arrived on each patch. Returns the per-patch
        increment [kg] and the per-patch enthalpy carried in [J]."""
        before, got = self.m_f.copy(), np.zeros_like(self.m_f)
        k = np.flatnonzero(fed > 0.0)
        if k.size == 0:
            return np.zeros_like(self.m_f), got
        owner = self.owner_area[k] > 0.0
        # owners: shared by patch area; interior elements: the nearest patch
        if owner.any():
            ko = k[owner]
            share, share_h = np.zeros(self.mesh.n_elements), np.zeros(self.mesh.n_elements)
            share[ko], share_h[ko] = fed[ko] / self.owner_area[ko], carried[ko] / self.owner_area[ko]
            self.m_f += share[self.surface.owner] * self.surface.areas
            got += share_h[self.surface.owner] * self.surface.areas
        if (~owner).any():                                    # interior elements: the NEAREST_PATCHES nearest patches, by area
            ki = k[~owner]
            kk = min(NEAREST_PATCHES, self.surface.n_patches)
            _, near = self._patch_tree.query(self.mesh.points[self.mesh.tets[ki]].mean(axis=1), k=kk)
            near = near.reshape(len(ki), kk)
            w = self.surface.areas[near] / self.surface.areas[near].sum(axis=1, keepdims=True)
            np.add.at(self.m_f, near.ravel(), (fed[ki][:, None] * w).ravel())
            np.add.at(got, near.ravel(), (carried[ki][:, None] * w).ravel())
        return self.m_f - before, got

    def _film_and_spray(self, t, dt, state, h_p, flow=None, delta_m=None):
        liq, s, mat = self.liquid, self.settings, self.material
        from . import spray as spray_mod
        if state is None or self.m_f.sum() <= 0.0 or flow is None:
            self.last_flow = self.last_spray = None
            return 0.0
        areas = self.surface.areas
        # the depth of liquid under each patch: its film plus the contiguous molten material beneath it. The branch
        # tests below ask whether the gas shear reaches the bottom of the liquid, so they belong on this; the fluxes
        # and the release stay on the film, which is the mass that can actually move (decided 2026-09-22).
        molten = self.molten_depth()
        # (ii) lubrication and runoff
        n_sub, moved = 0, 0.0
        if self.runoff is not None:
            q_of_b = lambda bb: self._film_mod.lubrication(flow.tau, flow.G, bb, delta_m, liq.mu, b_layer=bb + molten)[1]
            before = self.m_f
            self.m_f, n_sub, moved = self.runoff.transport(self.m_f, q_of_b, self.t_hat, liq.rho, areas, dt)
            self.runoff_mass += moved                                              # mass that arrived on another patch
            # film that runs to a colder patch takes its enthalpy with it and arrives at that patch's temperature;
            # the difference is released where it lands (the per-edge detail is not carried back, so it is shared
            # over the patches that gained film, which is exact in total and second order in the attribution)
            d = self.m_f - before
            self._spread_to_patches(-float((d * h_p).sum()), np.maximum(d, 0.0))
        b = self.m_f / (liq.rho * areas)
        layer = b + molten
        v_s, q, _, thick = self._film_mod.lubrication(flow.tau, flow.G, b, delta_m, liq.mu, b_layer=layer)
        # the molten surface: the area over which a wave could form at all, either wetted by the film or molten in its
        # own right. Its contiguous extent is what every mode's wavelength is measured against (2026-09-24).
        wetted = (b >= spray_mod.B_MIN) | (molten > 0.0)
        extent = self.region_extent(wetted)
        # Girin & Kopyt's W sin(Theta): the deceleration normal to the film, which on this body is W cos(phi). It is the
        # whole deceleration at the stagnation point and vanishes at the equator, where W lies in the surface.
        w_n = flow.deceleration * np.maximum(np.cos(self.theta), 0.0)
        # (iii) spraying
        res = self.spray.evaluate(flow, state, b, delta_m, v_s, self.windward, dt, areas, self.m_f,
                                  self.transverse_radius, b_layer=layer, extent=extent, deceleration_n=w_n)
        released = float(res.dm.sum())
        if released > 0.0:
            if self.spray_onset is None:
                self.spray_onset = t
            if s.particles:
                self.source_rows += spray_mod.source_rows(t, state.h, state.V, self.theta, self.surface.centroids, self.t_hat, flow, res, b, liq, state)
            hn, hm = spray_mod.histogram(res.r, res.dn, res.dm)
            self.hist_n += hn
            self.hist_m += hm
            self.removed_enthalpy += float((res.dm * h_p).sum())   # the droplets keep the surface's superheat, and go
            self.m_f = self.m_f - res.dm
            self.sprayed_mass += released
            self.removed_mass += released
            self.n_released += float(res.dn.sum())
        self.m_f[self.m_f < 1e-30] = 0.0                       # no denormal films (they made 0/0 coefficients in the runoff)
        self.last_flow, self.last_spray, self.last_b = flow, res, b
        # The Rayleigh-Taylor criterion is reported, never applied, so report it usefully: a body-level "any patch"
        # boolean says nothing about how much melt is involved or whether the unstable wave even fits on the nose.
        # `rt_mass_fraction` is the share of the film on unstable patches; `rt_wavelength_over_nose` is the shortest
        # unstable wavelength against the nose diameter -- above 1 the mode cannot develop on the cap at all.
        rt_on, rt_lam, rt_tau = spray_mod.rayleigh_taylor(w_n, layer, liq)
        # the depth criterion above only says the pool is deep enough for the fastest mode to see it as deep
        # (W h^2 rho > 3 Sigma is h > lambda*/2pi); whether a wave fits across the pool is a second, lateral question,
        # answered against the molten region each patch belongs to rather than the body's largest one
        rt_bound, rt_bl, rt_bt = spray_mod.rayleigh_taylor_bounded(w_n, layer, extent, liq)
        rt_extent = float(extent[np.asarray(rt_on) & (self.m_f > 0.0)].max()) if (np.asarray(rt_on) & (self.m_f > 0.0)).any() else 0.0
        wet = self.m_f > 0.0
        rt_sel = np.asarray(rt_on) & wet
        rt_mass = float(self.m_f[rt_sel].sum() / self.m_f.sum()) if self.m_f.sum() > 0.0 else 0.0
        with np.errstate(invalid="ignore"):
            rt_fit = float(np.nanmin(rt_lam[rt_sel]) / (2.0 * self.nose_radius())) if rt_sel.any() else float("nan")
        r = res.r[res.dm > 0.0]
        self.last_melt.update({
            "runoff_substeps": n_sub, "n_released": float(res.dn.sum()), "regime_fractions": flow.closure_fractions(areas, self.windward),
            "kn_body": flow.kn_body, "re_shock": flow.re_shock, "flow_branch": float(flow.branch),
            "kn_local_stag": float(flow.kn_local[self.i_stag]), "p_w_stag": float(flow.p_w[self.i_stag]),
            "phi_sonic_deg": math.degrees(flow.phi_sonic) if np.isfinite(flow.phi_sonic) else float("nan"),
            "theta_cr_deg": float(np.degrees(self.theta[res.unstable].min())) if res.unstable.any() else float("nan"),
            "spraying_area_m2": float(areas[res.unstable].sum()), "rt_active": float(res.rt_active),
            "r_median_um": float(np.median(r) * 1e6) if r.size else float("nan"), "r_max_um": float(r.max() * 1e6) if r.size else float("nan"),
            "film_thickness_max_mm": float(b.max() * 1e3),
            "rt_mass_fraction": rt_mass, "rt_wavelength_over_nose": rt_fit,
            "rt_growth_ms": float(np.nanmedian(np.asarray(rt_tau)[rt_sel]) * 1e3) if rt_sel.any() else float("nan"),
            "rt_region_mm": rt_extent * 1e3,
            "rt_bounded_fraction": float(self.m_f[np.asarray(rt_bound) & (self.m_f > 0.0)].sum() / self.m_f.sum())
            if self.m_f.sum() > 0.0 else 0.0,
            "rt_bounded_growth_ms": float(np.nanmedian(np.asarray(rt_bt)[np.asarray(rt_bound) & (self.m_f > 0.0)]) * 1e3)
            if (np.asarray(rt_bound) & (self.m_f > 0.0)).any() else float("nan"),
            "spray_growth_ms": float(np.nanmedian(res.growth[res.dm > 0.0]) * 1e3)
            if res.growth is not None and (res.dm > 0.0).any() and np.isfinite(res.growth[res.dm > 0.0]).any() else float("nan"),
            "molten_depth_max_mm": float(molten.max() * 1e3),
            "molten_depth_mean_mm": float((self.m_f * molten).sum() / self.m_f.sum() * 1e3) if self.m_f.sum() > 0.0 else 0.0,
            "delta_m_mean_um": self._mean_conjugate_um(delta_m),
            "thick_branch_fraction": float(thick[self.windward & (self.m_f > 0.0)].mean())
            if (self.windward & (self.m_f > 0.0)).any() else float("nan")})
        return released

    def _kill(self, dead):
        old_surface, old_m_f = self.surface, self.m_f
        gone, new = self.mesh.deactivate(dead)
        if self.mesh.n_active == 0:                              # nothing left: the run ends (demise) without a surface
            self.consumed = True                                 # the film still on it leaves with the body
            self.removed_mass += float(self.m_f.sum())
            self.removed_enthalpy += float((self.m_f * self.film_enthalpy()).sum())
            self.m_f = np.zeros(0)
            return
        self._refresh_geometry()
        m_f = np.zeros(self.surface.n_patches)
        keep = self.patch_of_face[old_surface.face_ids]
        kept = keep >= 0
        m_f[keep[kept]] += old_m_f[kept]
        lost = ~kept & (old_m_f > 0.0)
        if lost.any():
            # The film on a vanished patch has not moved: the surface receded *into* it, so it now rests on the faces
            # of that same element which have just become boundary -- the face on its inward side, whose neighbour is
            # the next layer in, and the walls of the pit its lateral neighbours expose -- shared by their areas.
            # Handing it to the nearest surviving centroid instead sends it sideways to whatever triangle happens to
            # be closest, which put a third of a collapsing body's film onto one 0.6 mm2 face (measured 2026-09-22),
            # and is a rule about the mesh's numbering rather than about where the liquid is.
            idx = np.flatnonzero(lost)
            targets = self.patch_of_face[self.mesh.element_faces(old_surface.owner[idx])]   # (n, 4), -1: not a patch
            ok = targets >= 0
            safe = np.where(ok, targets, 0)
            # weight by each exposed face's area *projected on the vanished patch's own outward normal*, not by its raw
            # area: the face the recession uncovered underneath lies parallel to the patch that vanished (its new
            # outward normal points the same way), while the pit walls stand perpendicular to it. A prism element's
            # side faces can out-area its inward face -- 0.25 x 1.5 mm slivers against a 1 mm2 floor -- so raw areas
            # route most of the film sideways into walls no shear is pushing it towards, and concentrate it by the
            # element's aspect ratio (3 to 8 here). The projection is the footprint each face offers the liquid.
            cosine = np.clip(np.einsum("nij,nj->ni", self.surface.normals[safe], old_surface.normals[idx]), 0.0, None)
            area = np.where(ok, self.surface.areas[safe] * cosine, 0.0)
            total = area.sum(axis=1)
            flat = ok & (total[:, None] <= 0.0)                  # no face faces the right way: fall back to raw areas
            if flat.any():
                area = np.where(flat, self.surface.areas[safe], area)
                total = area.sum(axis=1)
            has = total > 0.0
            if has.any():
                share = old_m_f[idx[has]][:, None] * area[has] / total[has][:, None]
                np.add.at(m_f, targets[has][ok[has]], share[ok[has]])
            if (~has).any():                                     # the element exposed nothing: its death opened a
                orphan = idx[~has]                               # hole right through, so fall back to the nearest patch
                np.add.at(m_f, self._patch_tree.query(old_surface.centroids[orphan])[1], old_m_f[orphan])
        self.m_f = m_f

    # -- reporting -----------------------------------------------------------------------------------------------
    def film_temperature(self):
        """The film's temperature per patch: the surface's own (the film is thermally thin, class docstring)."""
        return self.surface_temperature()

    def film_enthalpy(self, h_node=None):
        """Specific enthalpy of each patch's film [J/kg]: the mean over the patch's three nodes of the *liquid*
        enthalpy h(T_i) + L_f (1 - f_l) -- the film is liquid by construction, so it carries its latent heat wherever
        it sits. Two things this must not be: the enthalpy of the mean temperature (the mean of the enthalpies and
        the enthalpy of the mean differ by the whole latent heat across the ramp, which left a -0.4 residual in the
        coupled balance), and the mixture enthalpy h(T) (which makes melting free on the ramp)."""
        if h_node is None:
            h_node = self.material.enthalpy_liquid(self.solver.temperature())
        return h_node[self.surface.faces].mean(axis=1) if self.m_f.size else np.zeros(0)

    def film_energy(self):
        """Enthalpy held by the film [J]. Identical to the nodal sum the solver carries (`set_film_mass` gives each
        node sum_p m_p/3), so the film's sensible heat is inside 1^T M dT and the balance closes on it exactly."""
        return float((self.m_f * self.film_enthalpy()).sum()) if self.m_f.size else 0.0

    def film_frozen_fraction(self):
        """Fraction of the film sitting on patches below the feed ramp: mass the enthalpy calls solid but the model
        still treats as liquid, because the mesh cannot grow a crust outside itself and its own element may be gone
        (section 9.5). It is the honest measure of that limitation, so it is recorded every step."""
        total = float(self.m_f.sum())
        if total <= 0.0:
            return 0.0
        return float(self.m_f[self.material.feed_fraction(self.film_temperature()) <= 0.0].sum() / total)

    @staticmethod
    def _mean_conjugate_um(delta_m):
        """Mean conjugate melt-layer depth [um] over the patches that have one; nan where the gate denies them all."""
        ok = np.isfinite(delta_m)
        return float(delta_m[ok].mean() * 1e6) if ok.any() else float("nan")

    def _shear_reaches(self, delta_m):
        """Per element, may its melt leave for the film this step? An element that owns part of the wall is in the
        sheared layer by definition and always may. A buried element may only if its centroid lies within Girin's
        conjugate melt-layer depth of the nearest patch, because that is how far the gas shear penetrates the liquid;
        deeper melt is stagnant and stays where it is until the surface reaches it. Where the wall Knudsen gate denies
        the Girin closure there is no conjugate depth, and then only wall-owning elements may feed, which is the
        conservative reading of the same statement."""
        owner = self.owner_area > 0.0
        if delta_m is None or not self.surface.n_patches:
            return owner.astype(float)
        dist, near = self._patch_tree.query(self.mesh.points[self.mesh.tets].mean(axis=1))
        d = delta_m[near]
        return (owner | (np.isfinite(d) & (dist <= d))).astype(float)

    def region_extent(self, selected):
        """Lateral extent [m] of *each* patch's own contiguous selected region, as that region's equivalent diameter
        2 sqrt(A/pi); zero where the patch is not selected.

        This is the domain a surface wave has to fit inside, and the region -- not the mesh facet -- is the physical
        boundary of it: melt is continuous across a facet edge, so a facet width is a bookkeeping artefact and using it
        would make the limit mesh-dependent (decided 2026-09-24). Contiguity is taken on the patch adjacency graph, so
        melt separated by solid does not add up to one wide pool. Each patch gets its own region rather than the
        body's largest, since a patch in an isolated puddle must not be credited with the width of the main pool."""
        n = self.surface.n_patches
        sel = np.asarray(selected, dtype=bool)
        out = np.zeros(n)
        if not n or not sel.any():
            return out
        from scipy.sparse import coo_matrix
        from scipy.sparse.csgraph import connected_components
        _, i, j = self.surface.edges()
        keep = sel[i] & sel[j]
        i, j = i[keep], j[keep]
        g = coo_matrix((np.ones(i.size * 2), (np.concatenate([i, j]), np.concatenate([j, i]))), shape=(n, n))
        count, label = connected_components(g, directed=False)
        area = np.bincount(label[sel], self.surface.areas[sel], count)
        out[sel] = 2.0 * np.sqrt(area[label[sel]] / np.pi)
        return out

    def unstable_region_extent(self, unstable):
        """Lateral extent [m] of the largest contiguous flagged region: the maximum of `region_extent`, kept for the
        reported diagnostics. A Rayleigh-Taylor wave has to fit inside the region where the melt is deep enough to
        support it, not on the nose as a whole: comparing its wavelength against the nose diameter asks whether the
        wave would fit on a surface most of which carries no pool (corrected 2026-09-23)."""
        e = self.region_extent(unstable)
        return float(e.max()) if e.size else 0.0

    def molten_depth(self, Te=None, max_levels=8):
        """Depth of *contiguous* fully molten material under each patch [m], measured inward from the patch.

        Marches inward through the mesh from the patch's owner element, hopping to the face-adjacent element furthest
        along the inward normal, and stops at the first element that is not fully molten (below `Material.T_feed`) or
        at the mesh boundary. Only contiguous liquid counts: molten material separated from the patch by solid is
        blocked from it, cannot be part of the layer the gas shear sees, and must not be credited to it -- accumulating
        every nearby molten element instead would count melt that cannot reach the wall (decided 2026-09-22).

        This is the thickness that belongs in the thick/thin instability test against Girin's conjugate melt-layer
        depth delta_m: the film mass on a patch is the mobile inventory, not the depth of liquid beneath the wall.
        Measured on the 100 mm physics flight, the two differ by an order of magnitude and invert the branch choice."""
        mat, mesh = self.material, self.mesh
        if not self.surface.n_patches or not mat.melts:
            return np.zeros(self.surface.n_patches)
        if Te is None:
            Te = self.solver.temperature()[mesh.tets].mean(axis=1)
        molten = mesh.active & (Te >= mat.T_feed)
        centroid = mesh.points[mesh.tets].mean(axis=1)
        inward = -self.surface.normals
        cur = self.surface.owner.copy()
        alive = molten[cur]
        depth = np.zeros(self.surface.n_patches)
        rows = np.arange(len(cur))
        for _ in range(max_levels):
            if not alive.any():
                break
            reach = np.einsum("ij,ij->i", centroid[cur] - self.surface.centroids, inward)      # projected depth so far
            depth = np.where(alive, np.maximum(depth, reach), depth)
            pair = mesh._face_elements[mesh.element_faces(cur)]                               # (n, 4, 2)
            nb = np.where(pair[:, :, 0] == cur[:, None], pair[:, :, 1], pair[:, :, 0])
            ok = nb >= 0
            step = np.einsum("nkj,nj->nk", centroid[np.where(ok, nb, 0)] - centroid[cur][:, None, :], inward)
            pick = np.where(ok, step, -np.inf).argmax(axis=1)
            nxt = nb[rows, pick]
            good = (nxt >= 0) & (step[rows, pick] > 0.0)
            alive = alive & good & molten[np.where(good, nxt, 0)]
            cur = np.where(good, nxt, cur)
        # the chain's last centroid sits half an element short of the far wall of that element: add that half
        return np.maximum(depth, 0.0) * 1.5

    def liquid_layer_depth(self, Te=None):
        """Depth of liquid at each patch [m]: its film plus the contiguous molten material beneath it. What the gas
        shear sees, and so what decides whether the film is thick or thin against the conjugate melt-layer depth."""
        if not self.m_f.size:
            return np.zeros(0)
        return self.m_f / (self.liquid.rho * self.surface.areas) + self.molten_depth(Te)

    def film_thickness_max(self):
        """Thickest film [m] over the patches where a thickness means anything: at least a tenth of the median patch
        area (slivers excluded) and no deeper than the patch is wide (b <= sqrt(A)). A film deeper than its patch is
        wide is not a film -- m_f/(rho_l A) is then a volume per area, not a depth, and the lubrication picture the
        number belongs to has already failed. Those patches are reported separately by `film_blob_fraction` rather
        than quietly averaged in: without the second test this diagnostic read 2431 mm on a 50 mm sphere in the last
        step of its flight, where thousands of dying patches handed their film to a few 0.6 mm2 survivors (measured
        2026-09-22; the mass was real, the depth was not)."""
        if not self.m_f.size:
            return 0.0
        a = self.surface.areas
        b = self.m_f / (self.liquid.rho * a)
        ok = (a >= 0.1 * np.median(a)) & (b <= np.sqrt(a))
        return float(b[ok].max()) if ok.any() else 0.0

    def film_blob_fraction(self):
        """Fraction of the film sitting deeper than its patch is wide (b > sqrt(A)) -- melt the lubrication film
        model cannot describe, because it no longer fits on the facet it is booked to. It is a geometry limit of the
        patch graph, not a mass error: the mass and the energy are still accounted exactly. Measured on the 50 mm
        physics flight it is zero until the last 5 % of the flight, where the body collapses from 3026 patches to 56."""
        total = float(self.m_f.sum())
        if total <= 0.0:
            return 0.0
        a = self.surface.areas
        return float(self.m_f[self.m_f / (self.liquid.rho * a) > np.sqrt(a)].sum() / total)

    def film_thickness_mean(self):
        """Film mass over the area it meaningfully wets [m]: patches at least B_MIN deep, the same floor below which
        spraying ignores a film. Counting every patch holding a nonzero number instead makes this diagnostic
        hypersensitive -- a patch with 1e-20 kg on it adds its whole area to the denominator, so a last-bit difference
        in the temperature field that nudges one element across the feed ramp moves the reported mean by 1-2 % while
        the film mass is bit-identical. That is how a reporting artefact was mistaken for a physics change
        (2026-09-23); the mass-based diagnostics (`film_blob_fraction`, `film_frozen_fraction`) never had the problem,
        because a patch holding nothing contributes nothing to a mass fraction."""
        from .spray import B_MIN
        if not self.m_f.size:                                    # no film account at all (e.g. --removal instant)
            return 0.0
        wet = self.m_f > self.liquid.rho * self.surface.areas * B_MIN
        return float(self.m_f[wet].sum() / (self.liquid.rho * self.surface.areas[wet].sum())) if wet.any() else 0.0

    def demised(self):
        """True once the body's material (film excluded) is below the demise fraction of the initial mass."""
        return self.consumed or float((self.phi * self.element_mass).sum()) < self.settings.demise_fraction * self.mass0

    def energy_balance_residual(self):
        """(E_body + E_film - E0 - absorbed heat + removed enthalpy + pending load + dropped load) / absorbed heat:
        zero for an exact discrete balance (the deferred loads cancel between the FEM and the removed material)."""
        if not self.absorbed_heat:
            return 0.0
        return (self.energy() - self.energy0 - self.absorbed_heat + self.removed_enthalpy + self.pending_load.sum()
                + self.total_dropped_load) / self.absorbed_heat

    @property
    def total_applied_load(self):
        return getattr(self, "_applied_total", 0.0)

    @property
    def total_dropped_load(self):
        return getattr(self, "_dropped_total", 0.0)

    def melt_stats(self):
        lm = self.last_melt
        fractions = lm.get("regime_fractions", [0.0, 0.0, 0.0])
        return {"mass_kg": self.mass(0.0), "film_mass_kg": float(self.m_f.sum()), "sprayed_mass_kg": self.sprayed_mass,
                "runoff_mass_kg": self.runoff_mass, "removed_mass_kg": self.removed_mass,
                "melt_front_depth_max_mm": self.melt_front_depth() * 1e3, "equivalent_radius_mm": self.equivalent_radius() * 1e3,
                "n_active_elements": float(self.mesh.n_active), "spraying_area_m2": lm.get("spraying_area_m2", 0.0),
                "theta_cr_deg": lm.get("theta_cr_deg", float("nan")), "n_released": self.n_released,
                "released_mass_kg": lm.get("released_mass", 0.0), "r_median_um": lm.get("r_median_um", float("nan")),
                "r_max_um": lm.get("r_max_um", float("nan")), "closure_fraction_girin": fractions[0],
                "closure_fraction_couette_slip": fractions[1], "closure_fraction_couette_fm": fractions[2], "rt_active": lm.get("rt_active", 0.0),
                "kn_body": lm.get("kn_body", float("nan")), "kn_local_stag": lm.get("kn_local_stag", float("nan")),
                "re_shock": lm.get("re_shock", float("nan")), "flow_branch": lm.get("flow_branch", float("nan")),
                "p_w_stag_Pa": lm.get("p_w_stag", float("nan")), "phi_sonic_deg": lm.get("phi_sonic_deg", float("nan")),
                "drag_shape_factor": self.drag_shape_factor(), "frozen_mass_kg": self.frozen_mass,
                "film_T_max_K": float(self.film_temperature().max()) if self.m_f.size else float("nan"),
                "film_T_mean_K": float((self.m_f * self.film_temperature()).sum() / self.m_f.sum()) if self.m_f.sum() > 0.0 else float("nan"),
                "film_frozen_fraction": self.film_frozen_fraction(), "unapplied_load_J": float(self.pending_load.sum()),
                "film_blob_fraction": self.film_blob_fraction(),
                "rt_mass_fraction": lm.get("rt_mass_fraction", float("nan")),
                "rt_growth_ms": lm.get("rt_growth_ms", float("nan")),
                "rt_region_mm": lm.get("rt_region_mm", float("nan")),
                "rt_bounded_fraction": lm.get("rt_bounded_fraction", float("nan")),
                "rt_bounded_growth_ms": lm.get("rt_bounded_growth_ms", float("nan")),
                "spray_growth_ms": lm.get("spray_growth_ms", float("nan")),
                "rt_wavelength_over_nose": lm.get("rt_wavelength_over_nose", float("nan")),
                "molten_depth_max_mm": lm.get("molten_depth_max_mm", float("nan")),
                "molten_depth_mean_mm": lm.get("molten_depth_mean_mm", float("nan")),
                "delta_m_mean_um": lm.get("delta_m_mean_um", float("nan")),
                "thick_branch_fraction": lm.get("thick_branch_fraction", float("nan")),
                "removed_enthalpy_J": self.removed_enthalpy,
                "film_thickness_max_mm": self.film_thickness_max() * 1e3, "film_thickness_mean_mm": self.film_thickness_mean() * 1e3,
                "nose_radius_mm": self.nose_radius() * 1e3, "transverse_radius_mm": self.transverse_radius * 1e3,
                "fitted_nose_radius_mm": self.fitted_nose_radius * 1e3,
                "n_dead_elements": float(self.mesh.n_elements - self.mesh.n_active)}
```


- [ ] **Step 4: Extract the flat-disc endpoint of the drag family**

```bash
"$PY" - <<'EOF'
import scipy.io as sio, numpy as np, json
c = sio.netcdf_file("/Applications/DRAMA-4.1.4/TOOLS/SARA/REENTRY/data/ATDB_CYLINDER.nc", "r", mmap=False)
ld = np.asarray(c.variables["log_LengthToDiameter"][:])
j = int(np.argmin(10.0 ** ld))                       # thinnest disc
ma = np.asarray(c.variables["MachNumber"][:])
get = lambda k: [round(float(v), 6) for v in np.asarray(c.variables[k][:])[0, :, j]]
doc = {
    "_provenance": ("Values of /Applications/DRAMA-4.1.4/TOOLS/SARA/REENTRY/data/ATDB_CYLINDER.nc (ESA DRAMA 4.1.4 aerothermal "
                    "database for the primitive CYLINDER, Hyperschall Technologie Goettingen GmbH) at angle of attack 0 (the flat "
                    "face normal to the flow) and the thinnest tabulated aspect ratio L/D = %.1e, read on 2026-09-22. This is the "
                    "flat-disc limit of the sphere-to-disc shape family: the face-on continuum C_D is independent of L/D over six "
                    "decades (1.824 at Ma 10 for every thickness), because hypersonic drag on a blunt body is pressure drag on the "
                    "frontal area while the side wall is parallel to the flow and the base sits in a near-vacuum wake. HTG's "
                    "continuum database is modified Newtonian: sphere/disc = 0.4990-0.4991 at every Mach number, i.e. exactly the "
                    "Newtonian ratio 1/2, so cd_continuum here IS C_p,max(Ma)." % (10.0 ** ld[j])),
    "mach": [float(v) for v in ma], "cd_free_molecular": get("FMF_CD"), "cd_continuum": get("CON_CD"),
    "heat_flux_factor_free_molecular": get("FMF_AvHeatFlux"), "heat_flux_factor_continuum": get("CON_AvHeatFlux"),
}
json.dump(doc, open("reentry_model/data/atdb_disc.json", "w"), indent=2)
print(doc["cd_continuum"])
EOF
```

Expected: `[1.801341, 1.824148, 1.828404, 1.829896, 1.830587, 1.830962]`, and the sphere table divided by it is 0.49897 at every Mach number.

In `reentry_model/aero.py` add, next to `DEFAULT_ATDB`:

```python
DISC_ATDB = os.path.join(DATA_DIR, "atdb_disc.json")     # the flat-disc limit (ATDB_CYLINDER at zero angle of attack)
```

and give `drag_coefficient` the shape factor:

```python
def drag_coefficient(kn, ma, tables, bridging, shape_factor=1.0):
    """SESAM's sphere C_D: continuum and free-molecular table values blended by f(Kn).

    `shape_factor` (Step 3) scales the continuum entry for a body that is no longer a sphere: it is the
    modified-Newtonian drag of the current windward shape divided by that of a sphere, so 1 for a sphere and 2.00 for
    a flat face (the ATDB disc entry is exactly 2 x the sphere entry at every Mach number -- HTG's continuum database
    is itself modified Newtonian). The free-molecular entry is left alone: face-on, a disc and a sphere differ by only
    3-4 % there (2.24 vs 2.15 at Ma 10), and a melting body is deep in the continuum by the time it flattens."""
    cd_c = tables.cd_continuum(ma) * shape_factor
    cd_fm = tables.cd_free_molecular(ma)
    f = bridging(kn) if kn > 0.0 else 0.0
    return cd_c + (cd_fm - cd_c) * f
```

- [ ] **Step 5: Edit `reentry_model/trajectory.py`**

In `Simulator.aero_state`, replace the three lines

```python
            kn = aero.knudsen(fs.rho, fs.m_bar, self.settings.diameter)
            cd = aero.drag_coefficient(kn, ma, self.tables, self.bridging)
            a_drag = -0.5 * fs.rho * V * v_rel * cd * self.area / self.body.mass(t)
```

with

```python
            kn = aero.knudsen(fs.rho, fs.m_bar, self.body.reference_length() or self.settings.diameter)   # a melting body's current size (Step 3)
            cd = aero.drag_coefficient(kn, ma, self.tables, self.bridging, self.body.drag_shape_factor())
            area = self.body.reference_area() or self.area          # a melting body's projected area (Step 3), else pi D^2/4
            m = self.body.mass(t)
            a_drag = -0.5 * fs.rho * V * v_rel * cd * area / m if m > 0.0 else np.zeros(3)      # a consumed body (Step 3) has no drag
```

- [ ] **Step 6: Smooth SESAM's Mach-1 drag step in `reentry_model/aero.py`**

A light melting remnant hovers at its terminal velocity near Ma 1, where the factor-2 step in `cd_continuum` stalls the adaptive integrator (1.3e5 RHS evaluations in one macro step, measured). After `DEFAULT_ATDB = os.path.join(DATA_DIR, "atdb_sphere.json")` add

```python
MACH_SWITCH_LO, MACH_SWITCH_HI = 0.98, 1.02      # SESAM halves C_D below Ma 1; the step is smoothed over this band (below)
```

replace `SphereDragTables.cd_continuum` with

```python
    def cd_continuum(self, ma):
        if ma < MACH_SWITCH_LO:
            return 0.5 * float(self.cd_c[0])
        if ma < MACH_SWITCH_HI:                                  # SESAM's factor-2 step at Ma 1, smoothed over +-2 % (module docstring)
            s = (ma - MACH_SWITCH_LO) / (MACH_SWITCH_HI - MACH_SWITCH_LO)
            return (0.5 + 0.5 * s * s * (3.0 - 2.0 * s)) * float(self.cd_c[0])
        if ma < self.mach[0]:
            return float(self.cd_c[0])
        return float(np.interp(ma, self.mach, self.cd_c))
```

and extend the class docstring's last sentence to: `simply clamped (Kn is negligible wherever Ma < 5). The factor-2 step is applied as a smooth (cubic) ramp over Ma 0.98-1.02: a discontinuous C_D stalls the adaptive integrator when a light body hovers at its transonic terminal velocity (a melting remnant, Step 3: 1.3e5 RHS evaluations in one macro step, measured 2026-09-21); the reference spheres cross Ma 1 in a fraction of a second, where the ramp changes nothing measurable (the drag test excludes |Ma - 1| <= 0.02 rows for SESAM's 3-decimal Mach column)."""`. In `tests/test_reentry_model_aero.py` replace the assertion `assert t.cd_continuum(3.0) == 0.898818 and t.cd_continuum(1.0) == 0.898818` with

```python
        assert t.cd_continuum(3.0) == 0.898818 and t.cd_continuum(1.02) == 0.898818
        assert t.cd_continuum(1.0) == pytest.approx(0.75 * 0.898818) and t.cd_continuum(0.98) == 0.5 * 0.898818   # the Ma-1 step smoothed over +-2 % (Step 3)
```

- [ ] **Step 7: Run the tests**

```bash
"$PY" -m pytest tests/test_reentry_model_melting.py tests/test_reentry_model_coupled.py tests/test_reentry_model_trajectory.py tests/test_reentry_model_aero.py -q
"$PY" -m pytest -m reference tests/test_reentry_model_reference.py -q                        # the Step 1 flights, ~2 min: unchanged by the ramp
FI_PROVIDER=tcp CC=/Users/ashajain/miniforge3/envs/fenicsx_env/bin/clang "$FX" -m pytest tests/test_reentry_model_fenicsx.py -q
```

Expected: 6 passed (melting), the Step 2 coupled/trajectory/aero tests unchanged, the eight Step 1 reference flights within their thresholds; in `fenicsx_env` 8 passed — the two backends give the same melting run (mass to 1e-6, temperatures to 0.5 K, the same dead elements).

- [ ] **Step 8: Commit**

```bash
git add reentry_model/body.py reentry_model/trajectory.py reentry_model/aero.py reentry_model/data/atdb_disc.json tests/test_reentry_model_melting.py tests/test_reentry_model_aero.py tests/test_reentry_model_data.py
git commit -m "Add the melting body: feed, film, spraying, element death, accounting and the projected area (Step 3 Task 9)

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---


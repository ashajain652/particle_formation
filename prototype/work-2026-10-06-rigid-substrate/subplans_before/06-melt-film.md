# Sub-plan: Task 6 — The melt film: lubrication and runoff

> Extracted verbatim from `2026-09-20-melt-spraying.md` (current version, lines 2778–3009). Read `00-shared-context.md` first. **Amended 2026-10-02: the deep runoff and the per-patch conjugate depth** (section below). **Amended 2026-10-05 (runoff flux): on a thick patch the film carries only its own part of the conjugate layer's flux** (the section after it).

## Amendment of 2026-10-02 — the flux of the liquid below the conjugate depth

**Why (decision of 2026-10-02).** Asha approved a three-zone rule for where melted material belongs, written for the
large-fragment (Spheral) model that will read this model's exported frames: (1) the sprayable skin, liquid above the
908 K liquidus within Girin's conjugate depth δ_m of the surface (215–303 µm on the 100 mm flight, fact 28), stays in
the film and is sprayed by Girin's mechanism; (2) the contiguous liquid below δ_m, down to a film limit of about
2–3 mm, is not sprayed but runs off, driven by the pressure gradient along the surface and the deceleration; (3)
material thicker than the film limit goes to the large-fragment model. Until now the middle zone never moved: fact
28(b)'s feed gate holds it in its elements, and fact 29 recorded its runoff flux as "one consequence not
implemented". This task supplies that flux; Task 9 moves the liquid with it (sub-plan 09's amendment of the same date
holds the design, the alternatives considered and the measurements; facts 46–53 in `00-shared-context.md`).

**The change: one new function, nothing else in this module.** `film.deep_flux(G, b, h, mu_l)` returns
G (h³ − b³)/(3 μ_l), the flux per unit width of the liquid lying between the film (thickness b, on top) and the full
contiguous liquid depth h, plus one paragraph of the module docstring. `lubrication` and `Runoff` are not touched, so
the film's own flux and its thin/thick branch test are exactly what they were on every patch.

**Why this form.** Fact 29 asked for "two depth limits in the same expression": the runoff flux of a liquid column is

    q = τ d²/(2 μ_l) + G h³/(3 μ_l),     d = min(δ_m, h),

where τ is the gas shear at the free surface [Pa], G the driving gradient along the surface (the gas-pressure gradient
plus the deceleration term −ρ_l a sin θ, [Pa/m]), μ_l the liquid viscosity, the first term the shear-driven part over
the conjugate depth only and the second the pressure- and deceleration-driven part over the whole liquid depth.
G h³/(3 μ_l) is the Poiseuille flux of a half-channel of depth h with no slip at the solid and no stress from this part
at the free surface — the gravity-driven film of classical lubrication theory with G in place of ρg.
`lubrication`'s thick branch already gives the film the first term and G b³/(3 μ_l) of the second, which is what fact
29 called under-integrated; `deep_flux` is the remainder, carried by the liquid beneath the film. It has no shear
term: the gas shear reaches only the conjugate depth, and that layer is the film's.

**One approximation, declared.** The split between the film and the liquid beneath it is not the exact split of the
half-channel profile. In that profile a film of thickness b riding on the deeper liquid moves faster than G b³/(3 μ_l)
gives it; the share of the column's pressure-driven flux that this simpler split moves to the liquid beneath instead
is about 1.5 b/h — measured by quadrature: 0.15 % for a 1 µm film on a 1 mm column, 1.5 % for 10 µm, 15 % for
100 µm and 35 % for a 250 µm film, i.e. a film as thick as the conjugate depth. Where the film is comparable to the
liquid beneath it the liquid beneath is therefore driven too hard: by a factor 2.8 when b equals the deep thickness.
The column's total is exact either way. The split was kept because it leaves `lubrication` untouched and so cannot
change any patch that has no liquid beneath its film, and because the deep transport is close to its steady state in
every macro step anyway (fact 47), so the factor moves how fast the deep liquid reaches where it collects, not where
that is. The exact split (the film gaining (G/μ_l) b s (2b + s)/2, with s the deep thickness) is recorded in fact 53
as an option.

Code (the tested change; apply to `reentry_model/film.py` as given in the body below, which is the prototype's file):

```diff
--- a/reentry_model/film.py
+++ b/reentry_model/film.py
@@ -10,6 +10,13 @@
 gradient for the instability is V_s / b or V_s / delta_m). Where delta_m is undefined (free-molecular patches) the
 film is thin-branch (Couette) by definition.
 
+The liquid below the film (amendment of 2026-10-02): where contiguous liquid lies beneath the film, the column's
+runoff flux needs two depth limits in one expression -- the shear-driven part over min(delta_m, h), the pressure-
+gradient and deceleration-driven part over the whole liquid depth h (plan fact 29):
+    q_column = tau d^2 / (2 mu_l) + G h^3 / (3 mu_l),    d = min(delta_m, h).
+The film keeps its own share above; `deep_flux` is the rest, G (h^3 - b^3) / (3 mu_l), carried by the liquid
+beneath it, which the shear does not reach.
+
 Runoff: for every edge shared by patches i and j, flux = q_donor * l_edge * (t_donor . n_edge)^+ with n_edge the
 in-plane edge normal pointing out of the donor; both directions are evaluated (a negative q reverses the flow toward
 the nose). The transport over a macro step is the linearly implicit upwind scheme (I + dt C) m_new = m_old with the
@@ -50,6 +57,19 @@
     return V, q, rate, thick
 
 
+def deep_flux(G, b, h, mu_l):
+    """Flux per unit width [m^2/s] of the liquid beneath the film: the pressure-gradient and deceleration-driven part
+    of a column of contiguous liquid of depth h whose top b is the film,
+        q_deep = G (h^3 - b^3) / (3 mu_l)            (zero where h <= b; signed like G, negative toward the nose).
+    G h^3 / (3 mu_l) is the Poiseuille flux of a half-channel of depth h -- no slip at the solid, no stress from this
+    part at the free surface -- and `lubrication` already gives the film G b^3 / (3 mu_l) of it, so with the film's
+    shear-driven share over the conjugate depth the column carries exactly plan fact 29's
+    tau d^2 / (2 mu_l) + G h^3 / (3 mu_l). There is no shear term here: the gas shear reaches only Girin's conjugate
+    depth, and that layer is the film's (amendment of 2026-10-02)."""
+    G, b, h = (np.asarray(x, dtype=float) for x in (G, b, h))
+    return G * (np.maximum(h, b) ** 3 - b ** 3) / (3.0 * mu_l)
+
+
 class Runoff:
     """Explicit upwind transport of the film mass on a SurfaceMesh (rebuilt whenever the surface changes)."""
 
```

Test (the tested code, appended to `tests/test_reentry_model_film.py`):

```diff
--- a/tests/test_reentry_model_film.py
+++ b/tests/test_reentry_model_film.py
@@ -73,3 +73,31 @@
     assert m.sum() == pytest.approx(m0.sum(), rel=1e-12) and (m >= 0.0).all() and moved > 0.0
     assert m[~windward].sum() == 0.0 and np.degrees(theta[m > 1e-9 * m.max()].max()) > 35.0        # spread outward (>= 4 patches), never leeward
     assert ro.transport(np.zeros(s.n_patches), q_of_b, t_hat, LIQ.rho, s.areas, 0.5)[1] == 0    # a dry surface costs nothing
+
+
+# ---------------------------------------------------------------------------------------------------------------
+# Amendment of 2026-10-02: the liquid below the conjugate depth runs off under the pressure gradient and deceleration
+
+
+def test_deep_flux_is_the_poiseuille_flux_and_completes_fact_29s_two_depth_limits():
+    """With no film on top, the liquid beneath it carries the body-force-driven lubrication flux of a layer of depth h,
+    q = G h^3 / (3 mu): the half-channel profile u(y) = (G / mu)(h y - y^2 / 2), no slip at the solid and no stress
+    from this part at the free surface, integrated over the depth. Under a film of thickness b it carries the rest of
+    that flux, so that the thick branch's film flux plus this is plan fact 29's column flux
+    tau d^2 / (2 mu) + G h^3 / (3 mu) with d = min(delta_m, h): the shear part over the conjugate depth only, the
+    pressure part over the whole liquid depth."""
+    mu, G, h = LIQ.mu, 4.0e4, 1.0e-3
+    y = np.linspace(0.0, h, 200001)
+    u = G / mu * (h * y - 0.5 * y * y)
+    assert film.deep_flux(G, 0.0, h, mu) == pytest.approx(np.trapezoid(u, y), rel=1e-9)
+    assert film.deep_flux(G, 0.0, h, mu) == pytest.approx(G * h ** 3 / (3.0 * mu), rel=1e-12)
+    tau, delta_m, b = 30.0, 2.5e-4, np.array([1e-5, 1e-4, 4e-4])       # films thinner and thicker than delta_m
+    V, q, _, thick = film.lubrication(np.full(3, tau), np.full(3, G), b, np.full(3, delta_m), mu, b_layer=np.full(3, h))
+    assert thick.all()                                                  # the column is deeper than delta_m
+    column = tau * delta_m ** 2 / (2.0 * mu) + G * h ** 3 / (3.0 * mu)
+    assert q + film.deep_flux(G, b, h, mu) == pytest.approx(np.full(3, column), rel=1e-12)
+    assert film.deep_flux(G, 2e-3, h, mu) == 0.0                        # no liquid beneath a film as deep as the column
+    assert film.deep_flux(-G, 0.0, h, mu) == pytest.approx(-G * h ** 3 / (3.0 * mu))   # the deceleration's pull: toward the nose
+    # a thin column (h <= delta_m) is all within the shear's reach: lubrication alone is the column, unchanged
+    _, q_thin, _, thick_thin = film.lubrication([tau], [G], [5e-5], [delta_m], mu, b_layer=[2e-4])
+    assert not thick_thin[0] and q_thin[0] == pytest.approx(tau * 5e-5 ** 2 / (2.0 * mu) + G * 5e-5 ** 3 / (3.0 * mu))
```

**Measured (2026-10-02, throwaway copy of the prototype).** `"$PY" -m pytest tests/test_reentry_model_film.py -q`:
4 passed (the three tests above and the new one). The new test fails on the unamended module (`AttributeError`, no
`deep_flux`) and passes with it. The quadrature of the half-channel velocity profile reproduces G h³/(3 μ_l) to
1e-9, and the thick-branch film flux plus `deep_flux` reproduces fact 29's column flux to 1e-12 for films of 10, 100
and 400 µm under a 1 mm column.

## Amendment of 2026-10-05 (runoff flux) — on a thick patch the film carries only its own part of the conjugate layer's flux

> Asha's decision of 2026-10-05, option (a): scale the film's flux on a thick patch by the film actually present. Facts
> 69–77 in `00-shared-context.md` carry the measurements: fact 69 the switched-step experiment that found the defect,
> fact 70 the crash and its cause, facts 71–77 the fix, its tests and what it does to the flights. Measured in a
> throwaway copy of `prototype/proto3/` with the amendments of 2026-10-02, 2026-10-03 and 2026-10-05 (seeding) applied —
> verified before any change to be reproduced byte for byte by their 20 diff blocks — and the code below is the tested
> code, as a diff against that copy. This sub-plan holds the change and its tests; sub-plan 09's amendment of this date
> holds what the melting body sees of it (no change to `body.py`), a melting test whose tolerance had to move and the
> two-backend test.

**Why.** The time-step study of 2026-10-05 (fact 69) ran the 100 mm flight with the macro step shortened once the body
Knudsen number falls below 0.01, the model's own continuum boundary. At a fine step of 0.0125 s the run failed at 61.6 s
of flight time with a NaN film mass on every patch (fact 70). The cause is in this module. `lubrication`'s thick branch
gave the film the flux of the whole conjugate layer, V_s δ_m / 2, whatever the film's own depth b. The branch is decided
on the liquid layer beneath the wall (fact 28: film, molten material and deep liquid together), not on the film, so a
patch can be thick while its film is a few microns, or a few times 10⁻⁴⁰ m, deep. `Runoff.edge_coefficients` divides the
flux by the film to get an emptying rate, c = q ℓ (t̂·n̂)⁺/(A b), which therefore grew like 1/b as the film vanished;
each linearly implicit sub-step then emptied such a film by a factor of about c Δt_s, which raised the next sub-step's
rate by the same factor, so the film's depth roughly squared from one sub-step to the next. Rates reached 4.4·10¹⁹⁷ per
second, and the direct solve reported "Matrix is exactly singular" and returned NaN, which spread to every patch. The
second new test below reproduces that failure, warning and NaN alike, in the way it most likely arose: where two such
patches drain into each other (their surface directions point at each other, as in the stagnation region and the craters
of the eroded front) the solve loses its pivots to cancellation — 1 + c Δt − c Δt is 0 in floating point once c Δt
exceeds about 10¹⁶. Which pair failed in the flight was not identified: the probe saved the state of the last non-finite
sub-step, by which time every patch was NaN. The same formula was also a physics error: it moved a film thinner than δ_m
on a thick patch as fast as a whole conjugate layer, overstating the runoff there.

**The change: one expression in `lubrication`.** The thick branch's shear-driven velocity is linear within Girin's
conjugate layer, V_s = τ δ_m / μ_l at the free surface and zero at depth δ_m, and the film is the top b of that layer.
The film's shear-driven flux is therefore the integral of that profile over its own depth,

    q_shear = ∫₀ᵇ V_s (1 − z/δ_m) dz = V_s b (1 − b/(2 δ_m))      for b < δ_m,
    q_shear = V_s δ_m / 2                                           for b ≥ δ_m (unchanged),

with z the depth below the free surface. It equals V_s δ_m / 2 at b = δ_m, so the flux is continuous; it tends to V_s b
as the film vanishes, so the emptying rate's shear part is at most V_s ℓ / A however thin the film, about 10⁴ per second
on a millimetre facet; and it is monotone in b. The pressure- and deceleration-driven part G b³/(3 μ_l) is added as
before. Nothing else in the module changes: the surface velocity V_s (which the spray's Weber number uses), the shear
rate and the branch flag are those of before, the thin branch is untouched, and the flux at b ≥ δ_m is computed by the
same expression as before, so it is bit-identical (unit test).

**What else was checked.** (1) *Every caller.* `lubrication` is called twice per melt step, both in
`MeltingBody._film_and_spray`: once inside the film's runoff transport, whose flux this changes, and once after it, of
which only V_s (for the spray) and the branch flag (for `thick_branch_fraction`) are used, and those do not change.
Nothing else in the package calls it (`spray.py` mentions it only in docstrings; `girin_case` does not use the film).
(2) *The deep runoff's flux.* `deep_flux` = G((b + h_D)³ − b³)/(3 μ_l) gives the deep liquid an emptying rate per unit
edge length of q_deep/(A h_D) = G (b² + b h_D + h_D²/3)/(μ_l A), at most |G| (b + h_D)²/(μ_l A), which stays bounded as
h_D vanishes because the flux vanishes with h_D. Measured: in the crashed run's probe at 0.0125 s, before the fix, the
deep transport's largest rate was 5.3·10⁵ per second while the film's reached 4.4·10¹⁹⁷; with the fix both stay bounded
(fact 74). It needs no change. (3) *The column flux of the amendment of 2026-10-02.* That amendment's test showed the
thick-branch film flux plus `deep_flux` equal to fact 29's column flux τ d²/(2 μ_l) + G h³/(3 μ_l), d = min(δ_m, h).
That identity now holds where the film fills the conjugate layer (b ≥ δ_m) and falls short by τ (δ_m − b)²/(2 μ_l) where
it does not. The shortfall is the shear-driven flux of the lower part of the conjugate layer, which in this model is
liquid held in the elements — most of it, by fact 58, the wall-owning element that `molten_depth` counts whole although
its liquid has already been fed into the film — and no account moves it. Where deep liquid lies under a patch the
identity holds by construction at the start of every film transport, because the deep stage tops the film up to δ_m
before the film moves (sub-plan 09's amendment of 2026-10-02). The test is changed to say so, and `deep_flux`'s
docstring with it.

**Alternatives considered, and why they were not taken.**

- *A numerical floor on the film depth in `edge_coefficients`* (b replaced by max(b, b_min)). It stops the overflow but
  treats the symptom: the rate is then V_s δ_m ℓ/(2 A b_min), which a small enough b_min leaves as large as needed to
  break the solve and a large one makes a tuning constant with no physical meaning; and a film of a micron would still
  move like a whole conjugate layer.
- *A cap on the emptying rate* (c ≤ c_max, or c Δt_s ≤ K). The same objection, and it would also slow every patch whose
  legitimate rate exceeds the cap, the deep piles of fact 47 among them.
- *Scaling the whole layer's flux linearly by the share of the layer present,* q = (V_s δ_m/2)(b/δ_m) = V_s b/2. Bounded
  and continuous, but it gives a film at the top of the layer the layer's mean velocity V_s/2, so a thin skin would move
  at half the surface speed that the spray's Weber number credits it with.
- *Treating the film as resting on a stationary base,* the thin branch's Couette flux τ b²/(2 μ_l). Also bounded and, as
  it happens, continuous at b = δ_m, but it contradicts the thick branch's premise that the conjugate layer beneath the
  skin is liquid and sheared, and with it the surface velocity V_s = τ δ_m/μ_l that the spray uses; at b = δ_m/10 it
  moves the film nineteen times more slowly than the profile does.
- *Moving the lower part of the conjugate layer as well,* so that the column identity of 2026-10-02 holds everywhere.
  That liquid is held in its elements by the feed gate of fact 28(b), and moving it is the decision fact 53 (f) left
  open ("the shear passed down from above"); it is not part of a fix.
- *A switch to restore the old flux* (the pattern of `--deep-runoff` and `--molten-cascade`). Not added: the old flux is
  a defect that crashes the transport and overstates the runoff, not a modelling alternative; every earlier run remains
  reproducible from the amendments that made it.

Code (the tested change to `reentry_model/film.py`):

```diff
--- a/reentry_model/film.py
+++ b/reentry_model/film.py
@@ -10,6 +10,16 @@
 gradient for the instability is V_s / b or V_s / delta_m). Where delta_m is undefined (free-molecular patches) the
 film is thin-branch (Couette) by definition.
 
+A thick patch whose film is thinner than delta_m (amendment of 2026-10-05): the branch is decided on the liquid layer
+beneath the wall (plan fact 28), so a patch is thick wherever film, molten material and deep liquid together are
+deeper than delta_m, however little of it is film. The shear-driven velocity is linear within the conjugate layer,
+V_s at the surface and zero at depth delta_m, and the film is the top b of that layer, so the film carries
+    q = V_s b (1 - b / (2 delta_m)) + G b^3 / (3 mu_l)     for b < delta_m,
+which is V_s delta_m / 2 at b = delta_m (continuous) and V_s b as b -> 0, so the runoff's emptying rate q / (A b) stays
+below V_s l / A. Handing the whole layer's V_s delta_m / 2 to a vanishing film made that rate grow like 1/b, emptied
+the film by a factor c dt in one implicit sub-step, squared its depth from sub-step to sub-step and broke the direct
+solve (NaN film on every patch at a 0.0125 s macro step); it also moved such a film as fast as a whole conjugate layer.
+
 The liquid below the film (amendment of 2026-10-02): where contiguous liquid lies beneath the film, the column's
 runoff flux needs two depth limits in one expression -- the shear-driven part over min(delta_m, h), the pressure-
 gradient and deceleration-driven part over the whole liquid depth h (plan fact 29):
@@ -39,7 +49,9 @@
     branch is decided on the layer, because whether the gas shear penetrates the whole liquid or only its top delta_m
     is a property of the liquid's depth, not of how much of it the melt bookkeeping has mobilised; the fluxes stay on
     the film. Deciding the branch on the film instead put 79-99.8 % of the patches on the thin branch where the layer
-    says 94-100 % are thick (measured on the 100 mm physics flight, 2026-09-22).
+    says 94-100 % are thick (measured on the 100 mm physics flight, 2026-09-22). On a thick patch the film is the top
+    b of the conjugate layer and carries only its own part of that layer's shear-driven flux, V_s b (1 - b/(2 delta_m))
+    while b < delta_m (module docstring, amendment of 2026-10-05); V_s, the shear rate and the branch are unchanged.
     """
     tau, G, b = (np.asarray(x, dtype=float) for x in (tau, G, b))
     delta_m = np.asarray(delta_m, dtype=float)
@@ -49,7 +61,13 @@
     q_thin = tau * b * b / (2.0 * mu_l) + G * b ** 3 / (3.0 * mu_l)
     d = np.where(thick, delta_m, 0.0)
     V_thick = tau * d / mu_l
-    q_thick = V_thick * d / 2.0 + G * b ** 3 / (3.0 * mu_l)
+    # the film is the top b of the conjugate layer, whose velocity falls linearly from V_thick at the surface to zero
+    # at depth delta_m, so it carries V b (1 - b / (2 delta_m)) while it is thinner than delta_m and the whole layer's
+    # V delta_m / 2 once it fills it (amendment of 2026-10-05: the whole layer's flux handed to a vanishing film made
+    # the runoff's emptying rate grow like 1/b and the transport break down)
+    with np.errstate(divide="ignore", invalid="ignore"):
+        q_shear = np.where(b < d, V_thick * b * (1.0 - b / (2.0 * np.where(thick, d, 1.0))), V_thick * d / 2.0)
+    q_thick = q_shear + G * b ** 3 / (3.0 * mu_l)
     V = np.where(thick, V_thick, V_thin)
     q = np.where(thick, q_thick, q_thin)
     with np.errstate(divide="ignore", invalid="ignore"):
@@ -65,7 +83,11 @@
     part at the free surface -- and `lubrication` already gives the film G b^3 / (3 mu_l) of it, so with the film's
     shear-driven share over the conjugate depth the column carries exactly plan fact 29's
     tau d^2 / (2 mu_l) + G h^3 / (3 mu_l). There is no shear term here: the gas shear reaches only Girin's conjugate
-    depth, and that layer is the film's (amendment of 2026-10-02)."""
+    depth, and that layer is the film's (amendment of 2026-10-02). Since the amendment of 2026-10-05 the film carries
+    only its own part of the conjugate layer's shear flux, so the column total is fact 29's where the film fills the
+    conjugate layer (b >= delta_m) and short by tau (delta_m - b)^2 / (2 mu_l) where it does not: that is the shear
+    flux of the part of the conjugate layer held in the elements, which no account moves. The emptying rate this flux
+    gives the deep liquid, q_deep / (A h_D) <= G (b + h_D)^2 / (mu_l A) per unit edge length, stays bounded as h_D -> 0."""
     G, b, h = (np.asarray(x, dtype=float) for x in (G, b, h))
     return G * (np.maximum(h, b) ** 3 - b ** 3) / (3.0 * mu_l)
 
```

Tests (the tested change to `tests/test_reentry_model_film.py`: the test of 2026-10-02 changed as described above, and
two new tests appended):

```diff
--- a/tests/test_reentry_model_film.py
+++ b/tests/test_reentry_model_film.py
@@ -85,7 +85,9 @@
     from this part at the free surface, integrated over the depth. Under a film of thickness b it carries the rest of
     that flux, so that the thick branch's film flux plus this is plan fact 29's column flux
     tau d^2 / (2 mu) + G h^3 / (3 mu) with d = min(delta_m, h): the shear part over the conjugate depth only, the
-    pressure part over the whole liquid depth."""
+    pressure part over the whole liquid depth. Since the amendment of 2026-10-05 that holds where the film fills the
+    conjugate layer; a thinner film carries only its own part of the layer's shear flux, so the column is short by the
+    shear flux of the rest of the layer, tau (delta_m - b)^2 / (2 mu), which lies in the elements and no account moves."""
     mu, G, h = LIQ.mu, 4.0e4, 1.0e-3
     y = np.linspace(0.0, h, 200001)
     u = G / mu * (h * y - 0.5 * y * y)
@@ -95,9 +97,101 @@
     V, q, _, thick = film.lubrication(np.full(3, tau), np.full(3, G), b, np.full(3, delta_m), mu, b_layer=np.full(3, h))
     assert thick.all()                                                  # the column is deeper than delta_m
     column = tau * delta_m ** 2 / (2.0 * mu) + G * h ** 3 / (3.0 * mu)
-    assert q + film.deep_flux(G, b, h, mu) == pytest.approx(np.full(3, column), rel=1e-12)
+    held = tau * (delta_m - np.minimum(b, delta_m)) ** 2 / (2.0 * mu)  # the shear flux of the layer below a thin film
+    assert q + film.deep_flux(G, b, h, mu) == pytest.approx(np.full(3, column) - held, rel=1e-12)
+    assert held[2] == 0.0 and (held[:2] > 0.0).all()                    # exact where the film fills the conjugate layer
     assert film.deep_flux(G, 2e-3, h, mu) == 0.0                        # no liquid beneath a film as deep as the column
     assert film.deep_flux(-G, 0.0, h, mu) == pytest.approx(-G * h ** 3 / (3.0 * mu))   # the deceleration's pull: toward the nose
     # a thin column (h <= delta_m) is all within the shear's reach: lubrication alone is the column, unchanged
     _, q_thin, _, thick_thin = film.lubrication([tau], [G], [5e-5], [delta_m], mu, b_layer=[2e-4])
     assert not thick_thin[0] and q_thin[0] == pytest.approx(tau * 5e-5 ** 2 / (2.0 * mu) + G * 5e-5 ** 3 / (3.0 * mu))
+
+
+# ---------------------------------------------------------------------------------------------------------------
+# Amendment of 2026-10-05: on a thick patch the film carries only its own part of the conjugate layer's flux
+
+
+def test_a_thick_patch_moves_only_the_film_present():
+    """The thick branch's shear-driven velocity is linear within Girin's conjugate layer -- V at the surface, zero at
+    depth delta_m -- and the film is the top b of that layer, so a film thinner than delta_m carries
+    V b (1 - b / (2 delta_m)): the integral of the profile over its own depth. That goes to zero with the film, keeps
+    q / b below V however thin the film, and equals the whole layer's V delta_m / 2 at b = delta_m. The surface
+    velocity, the shear rate and the branch are those of before, and so is everything at b >= delta_m and on thin
+    patches, bit for bit. Before the amendment a vanishing film on a thick patch was handed the whole layer's
+    V delta_m / 2, so q / b grew like 1/b (the crash of the 0.0125 s time-step study, plan facts 69-70)."""
+    mu, tau, G, dm = LIQ.mu, 30.0, 4.0e4, 2.5e-4
+    V = tau * dm / mu
+    b = np.array([5e-324, 1e-300, 1e-200, 1e-100, 1e-30, 1e-12, 1e-9, 1e-6, 1e-5, 1e-4, 2.4e-4])
+    n = b.size
+    with np.errstate(over="ignore"):        # the thin branch's V / b, evaluated everywhere and discarded on thick patches
+        V_s, q, rate, thick = film.lubrication(np.full(n, tau), np.zeros(n), b, np.full(n, dm), mu, b_layer=np.full(n, 2e-3))
+    assert thick.all() and np.array_equal(V_s, np.full(n, V)) and np.array_equal(rate, np.full(n, V / dm))
+    assert (q > 0.0).all() and (q <= V * b).all()                       # vanishes with the film; q / b <= V
+    assert q[1:] / b[1:] == pytest.approx(V * (1.0 - b[1:] / (2.0 * dm)), rel=1e-12)
+    for bb in (1e-5, 1e-4, 2.4e-4):                                     # the profile's integral over the top b
+        z = np.linspace(0.0, bb, 100001)
+        q1 = film.lubrication([tau], [0.0], [bb], [dm], mu, b_layer=[2e-3])[1][0]
+        assert q1 == pytest.approx(np.trapezoid(V * (1.0 - z / dm), z), rel=1e-9)
+    # continuous at b = delta_m, with the pressure and deceleration part G b^3 / (3 mu) on top as before
+    q_at = [film.lubrication([tau], [G], [x], [dm], mu, b_layer=[2e-3])[1][0] for x in (dm * (1 - 1e-9), dm, dm * (1 + 1e-9))]
+    assert q_at[0] == pytest.approx(q_at[1], rel=1e-8) and q_at[2] == pytest.approx(q_at[1], rel=1e-8)
+    # unchanged, bit for bit: a film at least delta_m deep on a thick patch, and every thin patch
+    bt = np.array([dm, 2.6e-4, 1e-3, 5e-3])
+    tt, gt, dt_ = np.full(4, tau), np.full(4, G), np.full(4, dm)
+    V_old = tt * dt_ / mu
+    assert np.array_equal(film.lubrication(tt, gt, bt, dt_, mu, b_layer=bt + 1e-3)[1], V_old * dt_ / 2.0 + gt * bt ** 3 / (3.0 * mu))
+    for bn, dm_thin, extra in ((np.array([1e-9, 5e-5, 2e-4, dm]), dt_, 0.0),               # the layer within delta_m
+                               (np.array([1e-9, 5e-5, 1e-3, 5e-3]), np.full(4, np.nan), 1e-3)):  # no conjugate depth
+        V_n, q_n, _, thick_n = film.lubrication(tt, gt, bn, dm_thin, mu, b_layer=bn + extra)
+        assert not thick_n.any()
+        assert np.array_equal(q_n, tt * bn * bn / (2.0 * mu) + gt * bn ** 3 / (3.0 * mu))
+        assert np.array_equal(V_n, tt * bn / mu + gt * bn * bn / (2.0 * mu))
+
+
+def test_runoff_stays_finite_on_nearly_dry_thick_patches():
+    """Two thick patches whose surface directions point at each other drain into each other -- as they do in the
+    stagnation region and in the craters of the eroded front. With a film of 1e-40 kg on each, the thick branch's flux
+    before the amendment gave them emptying rates of 1e37 per second, the linearly implicit system lost its pivots to
+    cancellation (1 + c dt - c dt = 0) and spsolve returned NaN for every patch: the failure of the 0.0125 s flight,
+    "Matrix is exactly singular" included. Now the rate is bounded by V l / A, and the transport stays finite,
+    non-negative and conservative for nearly dry films all the way down to the smallest positive double. The deep
+    liquid's flux needed no change: its rate q_deep / (A h_D) <= |G| (b + h_D)^2 / (mu A) stays bounded as h_D -> 0."""
+    import warnings
+    from scipy.sparse.linalg import MatrixRankWarning
+    pts, s = strip_surface()
+    ro = film.Runoff(s, pts)
+    n = s.n_patches
+    t_hat = np.tile([1.0, 0.0, 0.0], (n, 1))
+    t_hat[1::2] = [-1.0, 0.0, 0.0]                                      # the two triangles of each square face each other
+    tau, dm, G = 30.0, 2.5e-4, np.full(n, -4.0e4)
+    molten = np.full(n, 2e-3)                                           # 2 mm of molten material under every patch
+    q_of_b = lambda b: film.lubrication(np.full(n, tau), G, b, np.full(n, dm), LIQ.mu, b_layer=b + molten)[1]
+    V = tau * dm / LIQ.mu
+    bound = (V + 4.0e4 * (1e-5 * 1.001) ** 2 / (3.0 * LIQ.mu)) * ro.length.max() / s.areas.min()
+    for tiny in (1e-20, 1e-40, 1e-100, 1e-200, 1e-300, np.finfo(float).tiny, 5e-324):
+        m = np.zeros(n)
+        m[10:14] = 1e-5 * LIQ.rho * s.areas[10:14]                      # ten microns of film on two squares
+        m[6] = m[7] = tiny                                              # and a nearly dry pair [kg]
+        b = m / (LIQ.rho * s.areas)
+        with warnings.catch_warnings(), np.errstate(over="ignore"):    # (the thin branch's V / b on a 5e-324 film)
+            assert np.concatenate(ro.edge_coefficients(q_of_b(b), b, t_hat, s.areas)).max() <= bound
+            warnings.simplefilter("error", MatrixRankWarning)           # no singular matrix
+            out = m.copy()
+            for _ in range(8):                                          # eight macro steps of 0.0125 s, 32 sub-steps
+                out, k, moved = ro.transport(out, q_of_b, t_hat, LIQ.rho, s.areas, 0.0125)
+        assert np.isfinite(out).all() and (out >= 0.0).all() and k == 4
+        assert out.sum() == pytest.approx(m.sum(), rel=1e-12, abs=0.0)
+    # the deep liquid beneath a film, nearly gone: a bounded rate, a finite transport
+    b_film = np.full(n, 1e-4)
+    deep = lambda hd: film.deep_flux(G, b_film, b_film + hd, LIQ.mu)
+    for tiny in (1e-40, 1e-200, 5e-324):
+        D = np.zeros(n)
+        D[6] = D[7] = tiny
+        D[10:14] = 1e-4 * LIQ.rho * s.areas[10:14]
+        hd = D / (LIQ.rho * s.areas)
+        c = np.concatenate(ro.edge_coefficients(deep(hd), hd, t_hat, s.areas))
+        assert c.max() <= 4.0e4 * (2e-4 * 1.001) ** 2 / LIQ.mu * ro.length.max() / s.areas.min()
+        with warnings.catch_warnings():
+            warnings.simplefilter("error", MatrixRankWarning)
+            out, _, _ = ro.transport(D, deep, t_hat, LIQ.rho, s.areas, 0.0125)
+        assert np.isfinite(out).all() and (out >= 0.0).all() and out.sum() == pytest.approx(D.sum(), rel=1e-12, abs=0.0)
```

**Measured (2026-10-05, throwaway copy).** `"$PY" -m pytest tests/test_reentry_model_film.py -q`: 6 passed (the four
tests of before, one of them changed, and the two new ones). Against the unamended module all three fail: the changed
test because the column identity is exact there for every film, the first new one because the old flux exceeds V_s b for
every film thinner than δ_m, and the second because the old emptying rate of a pair of nearly dry patches is 1.7·10¹⁷
per second against a bound of 1.6·10⁴. Run without that bound, the second test's transport reproduces the flight's
failure on the unamended module: with films between 10⁻²⁵ and 10⁻³⁰⁰ kg on the mutually draining pair, eight macro steps
of 0.0125 s end in "Matrix is exactly singular" and a NaN mass on every patch for six of the seven film masses tried
(10⁻⁶⁰ kg survived by chance), while with the fix the transport stays finite and conserves the mass to 2.5·10⁻¹⁴ for
films down to 5·10⁻³²⁴ kg, the smallest positive double, its largest rate 1.15·10⁴ per second, which is V_s ℓ (t̂·n̂)/A
for that strip. The flights are in facts 73–75: the 50 mm flight, which never has Girin's closure, is bit-identical; the
100 mm flight at the default step is bit-identical until the first step under his closure (49.5 s) and then moves its
masses by about the run-to-run scatter; and the switched-step flights at 0.0125 s, which crashed, and at 0.00625 s now
run to 120 s.

---

**Depends on:** Task 1 (surface mesh edges, normals, centroids).
**Produces, for later tasks:** the melt-film lubrication and runoff solver that Task 9 calls every step.
**Character:** numerics — small and self-contained.
**Read before implementing:** Measured fact 7 in the shared context: an explicit CFL-limited time-stepping scheme needed ten thousand to a hundred thousand sub-steps per macro step during prototyping and was rejected as too slow. The interface must use the linearly-implicit, unconditionally-stable scheme described below — do not substitute a more "obvious" explicit scheme.
**Refinement goal for the sub-agent:** turn the section below into a standalone implementation plan.

---

### Task 6: The melt film — lubrication and runoff

**Files:**
- Create: `reentry_model/film.py`
- Test: `tests/test_reentry_model_film.py`

**Interfaces:**
- Consumes: `SurfaceMesh.edges/normals/centroids/areas` (Task 1).
- Produces: `lubrication(tau, G, b, delta_m, mu_l, b_layer=None) -> (V_s, q, shear_rate, thick)` -- `b` is the mobile film, `b_layer` the contiguous liquid depth the branch test uses when given (fact 28); `Runoff(surface, points, windward=None)` with `.i, .j, .length, .n_i, .n_j, .n_patches`, `.edge_coefficients(q, b, t_hat, areas) -> (c_ij, c_ji)`, `.transport(m_f, q_of_thickness, t_hat, rho_l, areas, dt, substeps=4) -> (m_f, n_solves, moved_kg)`.

- [ ] **Step 1: Write the failing tests**

Create `tests/test_reentry_model_film.py`:

```python
"""film.py: the lubrication branches, the implicit upwind runoff (conservation, positivity, the strip steady state)."""
import numpy as np
import pytest

from reentry_model import film, material, mesh

LIQ = material.LiquidProperties(2400.0, 1.3e-3, 0.86)


def test_lubrication_branches():
    tau, G, mu = 30.0, 5e4, LIQ.mu
    V, q, rate, thick = film.lubrication([tau, tau, tau], [G, G, G], [5e-5, 5e-4, 5e-4], [1e-4, 1e-4, np.nan], mu)
    b, d = 5e-5, 1e-4
    assert not thick[0] and V[0] == pytest.approx(tau * b / mu + G * b * b / (2 * mu)) and q[0] == pytest.approx(tau * b ** 2 / (2 * mu) + G * b ** 3 / (3 * mu))
    assert rate[0] == pytest.approx(V[0] / b)
    b = 5e-4
    assert thick[1] and V[1] == pytest.approx(tau * d / mu) and q[1] == pytest.approx(V[1] * d / 2 + G * b ** 3 / (3 * mu)) and rate[1] == pytest.approx(V[1] / d)
    assert not thick[2] and V[2] == pytest.approx(tau * b / mu + G * b * b / (2 * mu))               # no delta_m: thin branch
    assert film.lubrication([tau], [-1e6], [1e-3], [np.nan], mu)[1][0] < 0.0                         # a strong adverse gradient reverses q


def strip_surface(n=20, w=1e-3):
    """A planar strip of n squares (2 triangles each) along +x, width w, in the plane z = 0 with the normal +z."""
    x = np.arange(n + 1) * w
    pts = np.array([[xi, yi, 0.0] for xi in x for yi in (0.0, w)])
    faces = []
    for k in range(n):
        a, b, c, d = 2 * k, 2 * k + 1, 2 * k + 2, 2 * k + 3
        faces += [[a, c, b], [b, c, d]]
    faces = np.array(faces)
    p = pts[faces]
    normals = np.cross(p[:, 1] - p[:, 0], p[:, 2] - p[:, 0])
    areas = 0.5 * np.linalg.norm(normals, axis=1)
    normals = normals / (2 * areas)[:, None]
    return pts, mesh.SurfaceMesh(faces, p.mean(axis=1), normals, areas, np.arange(len(faces)), np.arange(len(faces)))


def test_runoff_conserves_mass_and_reaches_the_strip_steady_state():
    """Constant shear on a strip fed at its first patch: in the steady state the thin-film flux through every edge
    equals the feed, so b = sqrt(2 mu S / (rho tau w)) on the interior patches (closed form) within 0.1 %."""
    pts, s = strip_surface()
    ro = film.Runoff(s, pts)
    assert len(ro.i) == 2 * (s.n_patches // 2) - 1                                     # 20 diagonals + 19 shared verticals
    t_hat = np.tile([1.0, 0.0, 0.0], (s.n_patches, 1))
    tau, w = 30.0, 1e-3
    S = 2e-6                                                                            # kg/s fed into patch 0
    q_of_b = lambda b: film.lubrication(np.full(s.n_patches, tau), np.zeros(s.n_patches), b, np.full(s.n_patches, np.nan), LIQ.mu)[1]
    m = np.zeros(s.n_patches)
    dt = 0.05
    for k in range(400):
        m[0] += S * dt
        total = m.sum()
        m, n, moved = ro.transport(m, q_of_b, t_hat, LIQ.rho, s.areas, dt, substeps=1)   # one implicit step per feed pulse
        assert m.sum() == pytest.approx(total, rel=1e-12) and (m >= 0.0).all() and n == 1
        m[-2:] = 0.0                                                                    # the strip's end is stripped (a sink)
    b = m / (LIQ.rho * s.areas)
    b_exact = np.sqrt(2.0 * LIQ.mu * S / (LIQ.rho * tau * w))
    assert np.abs(b[4:-4] / b_exact - 1.0).max() < 1e-3 and moved > 0.0


def test_runoff_on_the_sphere_stops_at_the_equator(coarse_sphere_mesh):
    s = coarse_sphere_mesh.surface()
    v = np.array([1.0, 0.0, 0.0])
    theta, t_hat = s.angles_to(v), s.tangent_from(v)
    windward = theta <= np.pi / 2
    ro = film.Runoff(s, coarse_sphere_mesh.points, windward)
    assert windward[ro.i].all() and windward[ro.j].all()
    m0 = np.where(theta < 0.5, 1e-4 * s.areas * LIQ.rho, 0.0)
    tau = np.where(windward, 30.0 * np.sin(theta), 0.0)
    G = np.where(windward, 5e4 * np.sin(theta) * np.cos(theta), 0.0)
    q_of_b = lambda b: film.lubrication(tau, G, b, np.full(b.size, np.nan), LIQ.mu)[1]
    m, n, moved = ro.transport(m0, q_of_b, t_hat, LIQ.rho, s.areas, 0.5)
    assert m.sum() == pytest.approx(m0.sum(), rel=1e-12) and (m >= 0.0).all() and moved > 0.0
    assert m[~windward].sum() == 0.0 and np.degrees(theta[m > 1e-9 * m.max()].max()) > 35.0        # spread outward (>= 4 patches), never leeward
    assert ro.transport(np.zeros(s.n_patches), q_of_b, t_hat, LIQ.rho, s.areas, 0.5)[1] == 0    # a dry surface costs nothing
```


- [ ] **Step 2: Run the tests to see them fail**

Run: `"$PY" -m pytest tests/test_reentry_model_film.py -q`
Expected: ImportError (`reentry_model.film`).

- [ ] **Step 3: Create `reentry_model/film.py`**

```python
"""The melt film on the surface patches: lubrication velocity and flux, and the explicit upwind runoff transport on
the patch graph (spec Step 3 section 8).

State: m_f per patch [kg], thickness b = m_f / (rho_l A). Lubrication solution with the patch's shear tau, driving
gradient G, thickness b, melt boundary-layer thickness delta_m (from Girin's Eq. 2, spray.melt_layer) and the
liquid viscosity mu_l:
    thin (b <= delta_m):  V_s = tau b / mu_l + G b^2 / (2 mu_l),   q = tau b^2 / (2 mu_l) + G b^3 / (3 mu_l)
    thick (b > delta_m):  V_s = tau delta_m / mu_l,               q = V_s delta_m / 2 + G b^3 / (3 mu_l)
(q per unit width [m^2/s], signed along the surface direction t away from the stagnation point; the velocity
gradient for the instability is V_s / b or V_s / delta_m). Where delta_m is undefined (free-molecular patches) the
film is thin-branch (Couette) by definition.

Runoff: for every edge shared by patches i and j, flux = q_donor * l_edge * (t_donor . n_edge)^+ with n_edge the
in-plane edge normal pointing out of the donor; both directions are evaluated (a negative q reverses the flow toward
the nose). The transport over a macro step is the linearly implicit upwind scheme (I + dt C) m_new = m_old with the
edge coefficients c = q l (t . n)^+ / (A b) [1/s] taken at the start of the step: unconditionally stable, positive,
and mass-conserving to round-off (every edge flux leaves one patch and enters another); its steady state is the
exact nonlinear one. The spec's explicit sub-stepped
scheme was replaced on 2026-09-20: films driven by the free-molecular shear near the equator move at ~10 m/s and
cross the hemisphere many times per 0.5 s step, so an explicit CFL needed 1e4-1e5 sub-steps per macro step."""
from dataclasses import dataclass

import numpy as np


def lubrication(tau, G, b, delta_m, mu_l, b_layer=None):
    """(V_s, q, shear_rate, thick) per patch for the thin/thick lubrication branches (module docstring); `thick` is the
    branch mask, true where the sheared layer does not reach the bottom of the liquid.

    `b` is the mobile film -- the mass that is actually available to move -- while `b_layer`, if given, is the depth of
    liquid beneath the wall (film plus the contiguous molten material under it, `MeltingBody.liquid_layer_depth`). The
    branch is decided on the layer, because whether the gas shear penetrates the whole liquid or only its top delta_m
    is a property of the liquid's depth, not of how much of it the melt bookkeeping has mobilised; the fluxes stay on
    the film. Deciding the branch on the film instead put 79-99.8 % of the patches on the thin branch where the layer
    says 94-100 % are thick (measured on the 100 mm physics flight, 2026-09-22).
    """
    tau, G, b = (np.asarray(x, dtype=float) for x in (tau, G, b))
    delta_m = np.asarray(delta_m, dtype=float)
    layer = b if b_layer is None else np.maximum(np.asarray(b_layer, dtype=float), b)
    thick = np.isfinite(delta_m) & (layer > delta_m)
    V_thin = tau * b / mu_l + G * b * b / (2.0 * mu_l)
    q_thin = tau * b * b / (2.0 * mu_l) + G * b ** 3 / (3.0 * mu_l)
    d = np.where(thick, delta_m, 0.0)
    V_thick = tau * d / mu_l
    q_thick = V_thick * d / 2.0 + G * b ** 3 / (3.0 * mu_l)
    V = np.where(thick, V_thick, V_thin)
    q = np.where(thick, q_thick, q_thin)
    with np.errstate(divide="ignore", invalid="ignore"):
        rate = np.where(thick, V / np.where(thick, d, 1.0), np.where(b > 0.0, V / np.where(b > 0.0, b, 1.0), 0.0))
    return V, q, rate, thick


class Runoff:
    """Explicit upwind transport of the film mass on a SurfaceMesh (rebuilt whenever the surface changes)."""

    def __init__(self, surface, points, windward=None):
        e, i, j = surface.edges()
        if windward is not None:                               # the flow stops at the equator: no edge leads into a leeward patch
            keep = windward[i] & windward[j]
            e, i, j = e[keep], i[keep], j[keep]
        x1, x2 = points[e[:, 0]], points[e[:, 1]]
        self.length = np.linalg.norm(x2 - x1, axis=1)
        mid = 0.5 * (x1 + x2)
        self.i, self.j = i, j
        self.n_i = self._edge_normal(x2 - x1, surface.normals[i], mid - surface.centroids[i])
        self.n_j = self._edge_normal(x2 - x1, surface.normals[j], mid - surface.centroids[j])
        self.n_patches = surface.n_patches

    @staticmethod
    def _edge_normal(edge, normal, outward):
        n = np.cross(edge, normal)
        n = n / np.maximum(np.linalg.norm(n, axis=1), 1e-300)[:, None]
        sign = np.sign(np.einsum("ij,ij->i", n, outward))
        return n * np.where(sign == 0.0, 1.0, sign)[:, None]

    def edge_coefficients(self, q, b, t_hat, areas):
        """Emptying-rate coefficients c_ij, c_ji [1/s] of every edge: the flux i -> j is m_i c_ij with
        c_ij = q_i l (t_i . n_ij)^+ / (A_i b_i) (zero where the patch is dry); q, b per patch."""
        with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
            rate = np.where(b > 0.0, q / (areas * np.where(b > 0.0, b, 1.0)), 0.0)          # 1/m per unit edge length
        rate = np.where(np.isfinite(rate), rate, 0.0)                                       # a degenerate (sliver) patch moves nothing
        c_ij = np.maximum(0.0, rate[self.i] * self.length * np.einsum("ij,ij->i", t_hat[self.i], self.n_i))
        c_ji = np.maximum(0.0, rate[self.j] * self.length * np.einsum("ij,ij->i", t_hat[self.j], self.n_j))
        return c_ij, c_ji

    def transport(self, m_f, q_of_thickness, t_hat, rho_l, areas, dt, substeps=4):
        """Advance m_f over dt in `substeps` linearly implicit upwind steps (I + dt_s C) m_new = m_old, the edge
        coefficients C taken at the start of each sub-step; each M-matrix system is solved directly, so the scheme is
        unconditionally stable, positive and conservative to round-off (the columns of C sum to zero), and its
        steady state C(m) m = 0 is the exact nonlinear one. Dry patches have no coefficient, so a wetting front
        advances one patch per sub-step (a documented limit; the film that matters is stripped where it forms). A
        Picard iteration on the fully implicit form does not contract when dt x c >> 1, which is the case for micron
        films at 30 Pa on millimetre patches (measured 2026-09-20). `q_of_thickness(b)` returns the signed flux per
        patch for thickness b. Returns (m_f, n_solves, mass that arrived on another patch [kg])."""
        import scipy.sparse as sp
        import scipy.sparse.linalg as spla
        m = np.array(m_f, dtype=float)
        if not (m > 0.0).any() or len(self.i) == 0:
            return m, 0, 0.0
        moved, dt_s, n = 0.0, dt / substeps, 0
        for n in range(1, substeps + 1):
            b = m / (rho_l * areas)
            c_ij, c_ji = self.edge_coefficients(q_of_thickness(b), b, t_hat, areas)
            out = np.bincount(self.i, c_ij, self.n_patches) + np.bincount(self.j, c_ji, self.n_patches)
            rows = np.concatenate([np.arange(self.n_patches), self.j, self.i])
            cols = np.concatenate([np.arange(self.n_patches), self.i, self.j])
            vals = np.concatenate([1.0 + dt_s * out, -dt_s * c_ij, -dt_s * c_ji])
            A = sp.csc_matrix((vals, (rows, cols)), shape=(self.n_patches, self.n_patches))
            m_new = np.maximum(spla.spsolve(A, m), 0.0)
            moved += float(np.maximum(m_new - m, 0.0).sum())
            m = m_new
        return m, n, moved
```


- [ ] **Step 4: Run the tests**

Run: `"$PY" -m pytest tests/test_reentry_model_film.py -q`
Expected: 3 passed (the strip's interior reaches the closed-form thickness to 1e-3 exactly at the scheme's steady state).

- [ ] **Step 5: Commit**

```bash
git add reentry_model/film.py tests/test_reentry_model_film.py
git commit -m "Add the melt film: lubrication branches and the linearly implicit upwind runoff (Step 3 Task 6)

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---


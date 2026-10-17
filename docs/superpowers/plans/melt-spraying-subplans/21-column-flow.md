# Sub-plan 21 — flow in the columns (stage 3 of the layered-skin design) — scope

> Read `00-shared-context.md` first, then the spec. **This sub-plan is a scoped outline, not yet a code-level plan.**
> Its tasks are written out with tested code after sub-plan 20 is built and its run time measured (spec §11.2).

**Goal:** The flowing layer moves: every skin cell carries a velocity in its column's tangent plane, coupled through the
depth by the viscosity at each interface, pushed by the gas pressure gradient, the deceleration and the two thickness
pressures, dragged by the gas shear at the top and held at the rigid front; mass, momentum and heat cross between
neighbouring columns layer by layer.

**Spec:** §6.3 (heat: the sideways conduction made implicit, amended 2026-10-17), §6.4 (flow), §6.5 (sideways exchange
and the arrival rule), §11 (build stage 3; level 1's flow and Rayleigh–Taylor checks; level 2 with flow on), §12 (the
risks this stage carries: capillary stiffness at the grid scale, the flowing-layer coordinate, momentum between tangent
planes).

**Builds on:** sub-plan 20 — the columns, the packed cells (with their unused velocity), the surface gradient and
curvature operators, and the base operations.

## Amendment of 2026-10-17 — the sideways conduction becomes implicit (Asha's decision)

> Design of record: spec §6.3 as amended on 2026-10-17; facts 147 and 148. **Write and build this task first, before
> the flow**: it is heat only, so the plates verify it with flow off before anything else changes.

**Why.** Sub-plan 20 conducts heat between neighbouring columns explicitly (`skin.sideways_heat`), sub-cycled inside each
sub-step to stay below the explicit stability limit, and the smallest columns set that limit: on the 100 mm sphere's
2 mm surface the columns' dual areas range from 1.26 to 4.59 mm² (median 3.31), so a cycle is 3.8 ms and every 10 ms
sub-step takes three, where the median column alone would need one. On the level-3 flight of the 100 mm sphere (no
melting; 14.6 h, 72 s a macro step) the sideways conduction took three fifths of a step, against an eighth for the
conduction through the depth that the skins exist to resolve (fact 147). Option A (fact 148) removed the work the
sub-cycles repeated without changing a result; what remains is the sub-cycling itself, which an explicit step cannot do
without. Task 4's transport of mass, enthalpy and momentum, layer by layer on the same depth matching, is already
linearly implicit (spec §6.5); the conduction joins it.

**What.** Backward Euler for the sideways conduction, so a sub-step has no sideways stability limit and no sub-cycles.
The code-level plan chooses, with a measurement of each on the 100 mm sphere's first six steps, between:

- (a) one implicit system per sub-step: the columns' depth tridiagonals coupled by the edges' bin exchanges, solved
  inside the existing Newton iteration by a Krylov method preconditioned by the depth solve (Thomas on every column).
  The sub-step stays one backward-Euler step.
- (b) a split sub-step: the depth solve, then an implicit sideways solve on the depth grid's bins (a graph Laplacian over
  the columns in each bin, with the exchange's k w l weights), the bins' heat given back to the cells by the existing
  overlaps. Simpler, but it adds a splitting error of first order in the sub-step.

Either keeps the exchange antisymmetric edge by edge, so the books close to round-off, and exact for a field linear
along the surface; either then carries Task 4's implicit transport.

**Tests (level 1).**
- Sub-plan 20's four still hold:
  - the heat sums to zero between columns of any depth;
  - it is exact for a field linear along the surface;
  - the steady two-dimensional temperature between held boundary columns is reached;
  - the checkerboard over 0.2 s, about 30 times the old limit on 1 mm columns, stays within its initial bounds, now in
    one step and monotone.
- New: a warm spot spreading sideways in an insulated slab, against the two-dimensional heat equation's exact solution.
  The check runs at sub-steps of 1, 3 and 10 times the old explicit limit, and the error must fall in proportion to the
  sub-step. The books close to round-off.
- Sub-plan 20's pin of the explicit heat to the build's (`test_the_sub_cycles_share_their_geometry_and_give_the_heat_as_built_to_the_last_bit`)
  is retired with the explicit step.

**Verification.**
- Levels 2 and 2b: the cold and hot plates, flow off, at both macro steps, pass as before. The changes from sub-plan
  20's results are reported, because the plates' heating falls off along a half cosine and so carries sideways heat.
- Level 3's 100 mm flight under SESAM-equivalent heating is flown and timed again. That one sideways solve a sub-step
  brings a flight to roughly 5–6 hours is an estimate, not a measurement. The solve must cost less than the three
  explicit cycles it replaces, or it is not worth having.

## Amendment of 2026-10-12 — port first, move into `shared/` afterwards (Asha's decision; supersedes the paths of the 2026-10-11 amendment below)

> Design of record: `docs/superpowers/specs/levelset/2026-10-11-level-set-deforming-body-design.md` §4.2, §4.4 and §14 as revised on 2026-10-12. **When this sub-plan's code is ported from
> the prototype into the package, use the module paths the prototype uses (`reentry_model/<module>.py`), exactly as the
> body and the earlier amendments say, and ignore the `reentry_model/shared/...` paths of the 2026-10-11 amendment
> below.** There is no `reentry_model/shared/` package on `main`. After the skyscrapers port has landed on `main`, one
> combined refactor (the level-set design's stage 0 + S, with its own plan) moves the architecture-independent modules
> into `reentry_model/shared/`, keeps every old path working as an alias of the same module object, and checks that
> every flight's history is unchanged. The 2026-10-11 amendment stays as the record of **which** of this sub-plan's
> modules will be shared; only **when** and **where they are created** is superseded. Helpful, not required: where it
> costs nothing, keep the shared physics free of skyscrapers-only imports, so that the later move stays mechanical.

**Shared later from this sub-plan:** the slurry viscosity law it reads; the column flow is skyscrapers' own.

## Amendment of 2026-10-11 — the shared physics package (Asha's instruction)

> Design of record: `docs/superpowers/specs/levelset/2026-10-11-level-set-deforming-body-design.md` §4 and §14 (stages 0 and S). The repository now has two model architectures,
> **skyscrapers** (this plan's) and **level set**, chosen by `reentry_model run --architecture skyscrapers|levelset`
> (default `skyscrapers`, whose runs, run names and results are unchanged). Code that both architectures use lives in
> `reentry_model/shared/`; each moved module's old path stays importable as a one-line re-export, so nothing that
> imports the old path breaks. **This amendment changes where files live, not what they do.** It takes precedence over
> the file paths in the body below and in earlier amendments. Test files keep their names and places and import the
> module from `reentry_model.shared`. The skins are never shared between the architectures. Read the matching amendment
> at the end of `00-shared-context.md` first.

**What changes in this sub-plan.**

- The viscosity in the cells comes from **`reentry_model/shared/rheology.py`** (S1); the column flow is skyscrapers' own.

## Tasks (to be written in full)

0. **What sub-plan 20 left to this stage** (its rulings 7 and 9, 2026-10-12): the cells' velocity, stored and carried
   through `Skins`' split, merge, take and add operations with exact momentum; the surface gradient and curvature
   operators on the corner graph (`columns.Columns`), tested against a sphere's exact curvature and a linear field's
   gradient; and the interface rule Asha chooses from sub-plan 20's Task 10 measurements.
1. **Viscosity per interface**: liquid 1.3 mPa·s; Li et al. (2014) slurry law with the bridge near the liquidus;
   fixed-point iteration on the shear rate (at most five passes, held at 367 1/s); rigid cells anchored. Tests against
   the law's tabulated values.
2. **Momentum through the depth**: implicit tridiagonal viscous coupling per column, no slip at the base, the gas shear
   on the top cell. Tests: shear alone gives a linear profile; the push alone a half-parabola; liquid over slurry the
   joined profiles with the kink at the viscosity jump; Stokes' first problem.
3. **The push**: −∇_s p_w − ρ_l a_t + ρ_l a_n ∇_s H − σ ∇_s κ on every layer, with the two thickness terms implicit
   (linearised in the column heights). Tests: a bump on a flat layer under deceleration grows at the thin-film
   Rayleigh–Taylor linear rate; bumps below the surface-tension cutoff decay.
4. **Sideways exchange** (in the solve of the implicit sideways conduction, amendment of 2026-10-17): implicit upwind
   transport of mass, enthalpy and momentum by fraction of the flowing depth; the
   arrival rule (same fraction where the receiver flows from its top down, on top of a rigid top otherwise); velocities
   projected between tangent planes. Tests: a translating pulse keeps its heat and arrives where it should; books to
   round-off; runoff landing on a solid skin freezes from below and the 3D surface moves out by the deposit's volume.
5. **The droplet momentum sink** placeholder at the top cells for stage 4, tested with a prescribed release.
6. **Verification**: level 2's plate (spec §11.1, amended 2026-10-12: the box's 90-degree edges fold the surface; sub-plan
   20's `analysis/skin_plate.py`) with flow on and a push along the face (mass, momentum and energy exact; sensible
   behaviour); the momentum-balance residual and the roughness history columns; the cost of the sideways solve.

## Inputs to settle before the code-level plan is written

- **The skins' interface** (sub-plan 20's Results; facts 118–125). Amended after the build (Asha, 2026-10-12): damped
  repetition, the profile valuation and the parts make the coupling converge and level 2 pass at 0.5 s, but a melting
  step on the 100 mm sphere costs up to 1 041 s (fact 125), and buried superheated cells remain at 6 MW/m² (fact 124).
  The onset excess at 0.1–0.25 s steps was cured by the melting part cap (fact 132); levels 2 and 2b now verify with the
  skins at the solid's density (`swell=False`), and compared like that the plates melt 2–6 % late on 2 mm elements, the
  3D model's resolution in depth (facts 133–137), now resolved by the fine band, 0.7 mm elements to 4 mm under the
  2 mm surface, which makes a melting step on the 100 mm sphere 45 % dearer (facts 138–140). To settle: whether flow's
  level-2 run uses the same setting, since with it the flowing layer's thickness is that of solid metal, not liquid.
- **The skins' 0.8 mm minimum** (option B′, Asha, 2026-10-16; facts 141–146): where the non-rigid layer is thinner than
  0.8 mm a skin's bottom is rigid metal, so the column flow must hold rigid cells still and count only the cells above
  T_rigid in the flowing depth; the finest cells stay in the top 0.4 mm (`SkinSettings.fine_zone`), and a column at its
  minimum holds 59 cells rather than 40, which the flow's per-cell work will feel. Flow makes each pass dearer,
  so the pass count is the first thing to settle: long parts during melting were rejected and needed many passes
  (fact 126), and parts of at most 0.0625 s while melting now cut a melting step to 122–226 s (fact 129). Absorption
  is limited only by the move allowance and the hot plate (level 2b) checks it; coarse surface meshes are not supported
  under strong melting, so flow's verification should use 2 mm surfaces (facts 130–131).
- Sub-plan 20's measured run time (fact 120 under the one-repeat rule; fact 125 under the amended default: 60–93 s a
  macro step solid, 123–1 041 s melting on the 100 mm sphere, almost all in the skins), which was to decide whether the
  sideways solve can be implicit over every cell or needs sub-cycling. Decided (Asha, 2026-10-17): implicit, the
  amendment above, from the level-3 flight's profile (facts 147 and 148).
- Whether the implicit thickness pressures are enough to keep grid-scale growth out (level 1's growth test), given
  that the capillary length (about 1.8 mm) is close to the column spacing.

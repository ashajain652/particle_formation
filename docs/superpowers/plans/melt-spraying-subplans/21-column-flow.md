# Sub-plan 21 — flow in the columns (stage 3 of the layered-skin design) — scope

> Read `00-shared-context.md` first, then the spec. **This sub-plan is a scoped outline, not yet a code-level plan.**
> Its tasks are written out with tested code after sub-plan 20 is built and its run time measured (spec §11.2).

**Goal:** The flowing layer moves: every skin cell carries a velocity in its column's tangent plane, coupled through the
depth by the viscosity at each interface, pushed by the gas pressure gradient, the deceleration and the two thickness
pressures, dragged by the gas shear at the top and held at the rigid front; mass, momentum and heat cross between
neighbouring columns layer by layer.

**Spec:** §6.4 (flow), §6.5 (sideways exchange and the arrival rule), §11 (build stage 3; level 1's flow and
Rayleigh–Taylor checks; level 2 with flow on), §12 (the risks this stage carries: capillary stiffness at the grid scale,
the flowing-layer coordinate, momentum between tangent planes).

**Builds on:** sub-plan 20 — the columns, the packed cells (with their unused velocity), the surface gradient and
curvature operators, and the base operations.

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
4. **Sideways exchange**: implicit upwind transport of mass, enthalpy and momentum by fraction of the flowing depth; the
   arrival rule (same fraction where the receiver flows from its top down, on top of a rigid top otherwise); velocities
   projected between tangent planes. Tests: a translating pulse keeps its heat and arrives where it should; books to
   round-off; runoff landing on a solid skin freezes from below and the 3D surface moves out by the deposit's volume.
5. **The droplet momentum sink** placeholder at the top cells for stage 4, tested with a prescribed release.
6. **Verification**: level 2's plate (spec §11.1, amended 2026-10-12: the box's 90-degree edges fold the surface; sub-plan
   20's `analysis/skin_plate.py`) with flow on and a push along the face (mass, momentum and energy exact; sensible
   behaviour); the momentum-balance residual and the roughness history columns; the cost of the sideways solve.

## Inputs to settle before the code-level plan is written

- **The skins' interface at 0.5 s macro steps** (sub-plan 20's Results; facts 118–119, built 2026-10-12). At the
  flight's step the deferred mismatch holds half the absorbed heat and the plate never melts, and the base
  reconciliation can bury superheated cells. Stage 2 is verified at 0.05 s only. Asha chooses the cure first; flow
  would build on the same interface.
- Sub-plan 20's measured run time (fact 120: 62–70 s a macro step solid, 82–133 s melting on the 100 mm sphere, 98–99 %
  in the skins), which decides whether the sideways solve can be implicit over every cell or needs sub-cycling.
- Whether the implicit thickness pressures are enough to keep grid-scale growth out (level 1's growth test), given
  that the capillary length (about 1.8 mm) is close to the column spacing.

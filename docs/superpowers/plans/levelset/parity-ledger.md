# Level set: reconciliation ledger with the skyscrapers

**Purpose.** This ledger records every change to the skyscrapers model (the layered-skin surface model) that the level-set
architecture must answer, from the parity target written into the design on 2026-11-10 (§20.10 of
`docs/superpowers/specs/levelset/2026-10-11-level-set-deforming-body-design.md`) onward. When the stage 3 plan is
written, this ledger becomes the pinned parity target, so nothing has to be rebuilt from the skyscrapers' history.

**How to use it.** Add a row whenever a skyscrapers change lands (a commit in the build, an approved design, or a decision
by Asha). Each row names the change, where it is recorded, the route by which it reaches the level set, the stage that
needs it, and its status. The routes follow the design's isolation rule (§4.1):

- **shared**: architecture-independent physics. It moves into `reentry_model/shared/` in stage 0 + S, and from then on
  both architectures import it unchanged.
- **skin amendment**: skin physics. By Asha's rule of 2026-10-11 it is never shared, so the level set gains it only
  through a written amendment to the design.
- **replaced**: the level set answers the same need with its own mechanism. The difference is deliberate and is
  reported in the parity comparisons, not removed.
- **not applicable**: the change fixes something that exists only in the column geometry of the skyscrapers.
- **decision needed**: the route is not obvious and is Asha's to choose.

Statuses: *open* (route not yet settled), *decided* (route settled, not built), *built* (in the level-set code), *pinned*
(part of the parity target fixed for stage 3).

Seeded 2026-11-11 from §20.10 and from the skyscrapers' build history and designs up to commit ef5e826 (2026-11-10). On
2026-11-11 the skyscrapers' planning copy (`prototype/work-2026-10-28-skin-spraying-plan/`, in the `subplan20` worktree)
also held uncommitted edits to `cli.py`, `skinbody.py`, `spray.py`, three test files and `part_22.md`. They were not
read, and must be checked and given rows when they are committed.

## Related work already under way

The shared-package refactor (stage 0 + S) has started on branch `levelset-stage0` (head da4c549, 2026-11-11). It has
already moved the body seam, the Step 2 thermal protocol, the run command's flight arguments and `heating.py` into
`reentry_model/shared/`. A row routed **shared** is closed once its module is in that package on the branch that stage 3
builds on.

## The ledger

### 1. The windward loads from a fitted sphere (sub-plan 22 Task 5g)

- **Skyscrapers change.** Until 2026-11-10 the windward loads used the sphere's original radius whatever its shape. These
  loads are the stagnation-point heating (Fay–Riddell, Sutton–Graves or DKR, each scaling as one over the square root of
  the radius) and the surface flow's pressure gradient, shear, shock Reynolds number and wall Knudsen number. Now a
  sphere is fitted every macro step to the tops of the skin columns whose normals lie within 30 degrees of the direction
  of motion. The flat-face bound (1.67 times the windward transverse radius, Boison and Curtiss 1959) applies on the flat
  side, by Asha's amendment of 2026-11-10 (option B). This is on by default (`--size-feedback current`).
- **Recorded in.** `docs/superpowers/specs/2026-11-10-fitted-nose-radius-design.md` and the plan
  `docs/superpowers/plans/2026-11-10-fitted-nose-radius-build.md`, both in the `subplan20` worktree.
- **Route: shared plus a level-set fit (decided by Asha, 2026-11-11).** The heating and surface-flow formulas that take
  the radius are shared. The fit itself belongs to the level set. Asha decided that the level set fits a sphere to the
  windward side of its gas–liquid level set φ (the outer surface, on top of any melt). The windward side is the part of
  φ's zero contour that lies within 90 degrees above and below the body's mass centroid.
  - *Interpretation, to be confirmed when it is planned:* a point of the surface counts as windward when the angle
    between the flight direction (pointing upstream) and the line from the mass centroid to the point is at most 90
    degrees. That is the half of the surface facing the flow, cut by the plane through the mass centroid normal to the
    flight direction.
  - *Two differences from the skyscrapers, both deliberate:* the skyscrapers select by the surface normal within 30
    degrees of the flow; the level set selects by position relative to the mass centroid within 90 degrees. Whether the
    flat-face bound of 1.67 times the transverse radius also applies is not yet decided for the level set.
- **Stage.** 3, the first stage that applies gas loads.
- **Status.** Decided (route and selection rule). Open: whether the flat-face bound applies, and whether the fit is
  weighted by patch area.

### 2. New defaults for melting skins (sub-plan 22 Task 5h)

- **Skyscrapers change.**
  - Re-frozen metal is no longer handed down into the core (`--hand-down off`).
  - The number of parts a macro step is split into is set by melting alone (`--part-count-from melting`).
  - The momentum double count with the hand-down off was fixed (ae91993).
- **Recorded in.** `docs/superpowers/specs/2026-11-08-hand-down-and-part-count-design.md`; build commits 968ccdd and
  ae91993.
- **Route: not applicable.** The hand-down and the part count belong to the columns. The level set evolves the skin base
  β and the molten–rigid surface ψ instead, and sets its own sub-steps (§10 of the design).
- **Stage.** None.
- **Status.** Decided.

### 3. The Rayleigh–Taylor droplet without the body-radius cap (sub-plan 22 Task 5h)

- **Skyscrapers change.** The Rayleigh–Taylor droplet size is no longer capped at a quarter of the body's radius.
- **Recorded in.** Build commit 968ccdd.
- **Route: shared**, as part of the spray model.
- **Stage.** 3.
- **Status.** Open until the spray model is in the shared package.

### 4. Knudsen-number gates and spraying in the slip band (sub-plan 22 Task 5i)

- **Skyscrapers change.**
  - Girin and Kopyt's thin branch sprays only in the continuum regime (`--spray-gates 2026-11-10`).
  - With `--slip-band on`, both of Girin's branches may also spray in the slip band. That is where the body's Knudsen
    number (the gas's mean free path over the body's size) lies between 0.01 and 0.1 and the patch's own wall Knudsen
    number is below 0.1.
  - In the slip band, the melt layer's depth comes from a momentum balance. This follows the plan of 2026-11-06, by Asha's
    amendment of 2026-11-10.
- **Recorded in.** Build commit b91d69f; `docs/superpowers/plans/2026-11-06-slip-band-spraying.md`; the fitted-radius
  design's amendment.
- **Route.** The gates are **shared** with the spray model. The momentum-balance depth reads a skin's layer, so it is a
  **decision needed**: the level set's skins can supply the same depth through an amendment, or the depth can come from
  the level set's own resolved layer.
- **Stage.** 3.
- **Status.** Open.

### 5. The leeward test by the mass centroid (sub-plan 22 Task 5j)

- **Skyscrapers change.** A patch is leeward when it lies behind the body's mass centroid (`--leeward-test centroid`, the
  default).
- **Recorded in.** Build commit ef5e826.
- **Route: replaced.** The level set uses its ray-cast shadow test (§20.8). The centroid test is kept, through the shared
  package, as the one compatibility switch `--leeward centroid` for the parity runs (§20.14 item 3, decided by Asha
  2026-11-10).
- **Stage.** 3.
- **Status.** Decided.

### 6. Deposit spreading

- **Skyscrapers change.** Re-frozen deposits spread to neighbouring cells to keep the surface's angles reasonable.
- **Recorded in.** `docs/superpowers/specs/2026-11-06-deposit-spreading-design.md` (built, not yet flown).
- **Route: decision needed.** This is skin physics, so it is a skin amendment under the rule. A level-set surface carries
  deposits as metal under φ, though, and may not need a spreading rule at all.
- **Stage.** 3.
- **Status.** Open.

### 7. Re-frozen metal kept in the skins, and sideways heat between cells

- **Skyscrapers change.** Proposed only: re-frozen metal stays in the skins, and heat flows sideways between cells that
  touch.
- **Recorded in.** `docs/superpowers/specs/2026-11-08-refrozen-metal-in-skins-design.md` (draft, nothing built).
- **Route.** Not applicable as drafted. The level set already keeps re-frozen metal as metal under φ, and its 3D thermal
  solve (stage 2) carries sideways heat. Revisit if the draft is built.
- **Stage.** None.
- **Status.** Watching.

### 8. The viscous-film gate (branch `we-critical-308`)

- **Skyscrapers change.** The thin-film spray gate gains a viscous film model, with We_cr = 3.08 (fa02a4c, a54520c).
- **Recorded in.** Branch `we-critical-308`, not merged.
- **Route: shared**, with the spray model, if it is merged before the target is pinned (§20.10).
- **Stage.** 3.
- **Status.** Watching (unmerged).

### 9. The rest of the 2026-11-10 parity target

The items below were already listed in §20.10 and have not changed since: layered skins with a velocity per cell, the
10 µm fine zone and graded cells; runoff over the whole surface with the arrival rule and deposits freezing on cold
metal; Girin's thick and thin branches with the conjugate depth; the slurry spray switch.

- **Route.** Skin rules by amendment. Girin's spray model is shared.
- **Stage.** 3.
- **Status.** Pinned when stage 3's plan is written.

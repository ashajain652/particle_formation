# Handoff: the level-set architecture with two tracked interfaces (2026-11-10)

Written by the skyscrapers session (sub-plan 22) for a new session that Asha is starting on the level-set architecture.
Read this first. Paths are relative to a checkout of this repository; the files marked *subplan20 only* exist so far only
in the worktree `.claude/worktrees/subplan20` (uncommitted there), so read them by absolute path:
`/Users/ashajain/Documents/University Documents /MIT Graduate Work/Research/Space Sustainability/Particle Wake Evolution/.claude/worktrees/subplan20/<path>`.

## 1. What Asha wants

Asha, 2026-11-10, in her words: "I think I want to use a level-set method to track the surface layer with greater
refinement which the skyscrapers currently approximate. … I think the level set will enable an arbitrary shape to be
tracked for both the liquid-gas interface and the molten-rigid interface. Then we can discretize the 'molten' space in
between these two regions with our denser skyscraper-like cells." And: "I want the shape of the melt to be driven by
physics, not our gridding or column assumptions."

**Your first job is the design, not code**: amend the existing level-set design for this two-interface idea with Asha
(brainstorming), and get her approval before any plan or code. Asha's reporting rules are in her global `CLAUDE.md`
(full sentences, every number explained, decisions with options and a recommendation). Use they/them for Asha.

## 2. What already exists

- **The design**: `docs/superpowers/specs/levelset/2026-10-11-level-set-deforming-body-design.md` (draft for Asha's
  review). One level set on a fixed tetrahedral box mesh; its zero contour is the bulk metal's surface and the **base** of
  thin skins carried on surface markers (a 0.4 mm metal zone plus a liquid zone, §9); a sharp cut-cell thermal solve (§7);
  an incompressible Navier-Stokes flow solve in the metal where metal more than half liquid forms a region at least three
  cells thick, with the half-liquid surface (895.1 K) as a no-slip wall (§8; diffuse interface first). So in that design
  the molten-rigid boundary is an isotherm of the thermal field, not a tracked surface, and deep molten metal goes into
  the volume flow solve rather than into columns. Decisions of record with dates in §19, open items in §17.
- **Item 8 of §17** (*subplan20 only*, added 2026-11-10, uncommitted): the capillary curvature is to come from the level
  set, kappa = div(grad phi / |grad phi|), not the skyscrapers' thin-film rule; copy it into your branch.
- **Asha's isolation rule (2026-10-11)**: the existing approach is called **skyscrapers**; level-set code goes in
  `reentry_model/levelset/`, code shared by both in `reentry_model/shared/` (old paths as re-exports), a run flag
  `--architecture skyscrapers|levelset`, and the level set's docs, plans, prototype and tests in `levelset` folders. The
  skins are **not** shared between the architectures (Asha, 2026-10-11).
- **Status**: stage 0 (the shared-package refactor and the flag) was built on branch `levelset-stage0` (worktree
  `.claude/worktrees/levelset-stage0`, head da4c549) and never merged. On 2026-10-12 Asha changed the order: the
  skyscrapers' Step 3 code is ported first, then one combined refactor (stage 0 plus the shared-physics sub-plans S1-S6,
  index `docs/superpowers/plans/shared/README.md`) whose plan is not yet written; da4c549's tested pieces may be reused,
  the branch itself is not to be merged. Stage 0's plan: `docs/superpowers/plans/levelset/00-stage0-isolation-and-shared-package.md`.
- **The October assessment** (why a level set, and what it costs): Asha's memory note `level-set-assessment.md` in
  `/Users/ashajain/.claude/projects/-Users-ashajain-Documents-University-Documents--MIT-Graduate-Work-Research-Space-Sustainability-Particle-Wake-Evolution/memory/`
  (Chen et al. 2024; resolving the windward film is unaffordable, so the film stays sub-grid; cost tiers; the background
  mesh about 3.6 million tetrahedra in a 150 mm box at 2 mm, at most about 0.55 million holding metal, about 25 times more
  at 1 mm, estimated).

## 3. Why it matters now: what the skyscraper flights showed (measured 2026-11-08 to 2026-11-10)

- **An overhang stopped the 100 mm sphere.** The flight with the part count set by melting stopped at 39.69 s when a
  surface corner's move "degrades the mesh even after a rebuild". Re-frozen metal handed down into the core built a ridge
  that travelled back over the shoulder (0.5 mm high at 20 degrees from the nose at 25 s, about 6.8 mm high at about 75
  degrees by 39.5 s); its back face became vertical by 37 s and leaned over its own foot at 39.0 s. A column is a single
  height along its base normal and cannot represent that. Figure and script (*subplan20 only*):
  `prototype/work-2026-10-28-skin-spraying/measure/task5_compare/d100__parts-melting_fold_cross_section.png` and
  `plot_fold_cross_section.py` beside it.
- **The 50 mm sphere melted almost through.** By 187.0 s the core held 10.7-14.6 % of its volume at 171 s, its nose had
  receded 4.5-5.4 mm past the original centre (mesh radius 24.2 mm), and the melt sat in a rear pool: about half of the
  2 367 columns, all 120-178 degrees from the nose, were more than 20 mm deep, the deepest 51.5-66.5 mm, fully liquid.
  Fixed-area columns put that pool's surface about twice as far out as the same volume would sit as a spherical shell
  (about 22-24 mm deep). Spray came from the ragged, misplaced column tops, and the skin-core interface check
  (`jump_faults`) rose from zero to 190 000-334 000 (cause still under investigation).
- **Three column assumptions** set the geometry rather than physics: a fixed cross-section (the base patch's area), a
  single-valued height along the base normal, and a thin-film capillary curvature (base curvature minus the Laplacian of
  the thicknesses; Asha's option A of 2026-10-21, chosen over the top surface's mesh curvature, which drew metal onto
  pile-ups at corners touching obtuse triangles).

## 4. Open questions to settle with Asha before designing further

1. **The cells between the two surfaces.** Straight columns along one surface's normal cross each other or miss the other
   surface where the molten region is thick compared with its curvature (the 50 mm rear pool was as deep as the core was
   wide) or where a surface folds; columns following the distance field from the rigid surface widen correctly but
   collide where columns from opposite sides meet. The skyscrapers session's view: the existing design's split (columns
   where the molten layer is thin, a volume solve on the background mesh where it is thick) is probably still needed.
2. **What the molten-rigid interface is**: a second tracked level set, or the half-liquid isotherm of the thermal field as
   the existing design has it, and how the two must stay consistent.
3. **Leeward mass loss.** Nothing removes metal from the rear: the ring-and-shedding plan
   `docs/superpowers/plans/2026-09-27-deformation-ring-shedding.md` (a rear cap or ring; capillary pinch-off, tearing of
   semi-solid metal, forcing by the wake) is not built. Without it a rear pool grows without limit in either architecture.
4. **Loads on a deformed shape.** The gas loads come from each surface patch's own angle to the flight direction, with no
   shadowing; sub-plan 22 Task 5j is adding a leeward test by the body's mass centroid to the skyscrapers. The design's
   §8.4 mentions a shelter test; §17 item 3 keeps a CFD hook for later.

## 5. Shared physics changing this week in the skyscrapers (sub-plan 22 Tasks 5g-5j)

Plan `docs/superpowers/plans/2026-11-10-fitted-nose-radius-build.md`, spec
`docs/superpowers/specs/2026-11-10-fitted-nose-radius-design.md` (the later sections *subplan20 only* until committed):
the windward loads from a sphere fitted to the gas-molten interface (5g, built, 248c617); new defaults (hand-down off,
part count from melting) and the Rayleigh-Taylor droplet without the quarter-body-radius cap (5h); the thin branch gated
to the continuum and Girin's spraying in the slip band (body and wall Knudsen numbers below 0.1, conjugate depth from a
momentum balance), following `docs/superpowers/plans/2026-11-06-slip-band-spraying.md` (5i); the leeward test by the mass
centroid (5j). These are spraying and loads physics the level-set architecture will want through the shared package.

## 6. Working rules

- Work in your own worktree on a new branch from `step3-rigid-substrate-scheil` at or after 248c617. Never edit the
  `subplan20` worktree, the planning copy `prototype/work-2026-10-28-skin-spraying-plan/`, the build repository
  `prototype/work-2026-10-28-skin-spraying/`, or anything under `code/reentry_model_output/`: the skyscrapers build and a
  flight run from there.
- Never touch the DRAMA databases, the wrapper scripts `sphere_reentry.py` and `sphere_sweep.py`, the package at the
  repository root, `prototype/proto3/`, `plan3`, or `.claude/worktrees/we-critical-308`. Other sessions' files are read
  only: `docs/superpowers/plans/melt-spraying-subplans/23-notes.md`, `docs/superpowers/plans/2026-11-06-slip-band-spraying.md`,
  `docs/superpowers/specs/2026-11-08-refrozen-metal-in-skins-design.md`.
- The laptop is shared with a running flight, a long diagnostic run and this week's builds (each unit tier takes about
  37 minutes): run tests only on mains power (`pmset -g batt`), keep them focused, and plan long runs for the Mac Studio.
- Commit only when Asha asks; never push without asking; stage named paths only, never `git add -A`; never bare
  `git stash`. Paths contain a space: quote them. Interpreter `/Users/ashajain/miniforge3/envs/drama_env/bin/python`.

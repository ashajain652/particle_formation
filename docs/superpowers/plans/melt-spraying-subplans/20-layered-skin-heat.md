# Sub-plan 20 — skins on every corner, heat only (stage 2 of the layered-skin design) — scope

> Read `00-shared-context.md` first, then the spec. **This sub-plan is a scoped outline, not yet a code-level plan.**
> Its tasks are written out with tested code, as sub-plan 19's are, after sub-plan 19 is built and measured: the spec
> gates each build stage on the one before (§11.2), and three of this stage's inputs come from stage 1's measurements
> (listed at the end). That is the same order sub-plan 18 followed for its later stages.

**Goal:** Every corner of the 3D mesh's surface carries a layered skin from time zero that holds the outermost metal —
its minimum mass, and everything more than half liquid — in graded cells, coupled to the moving 3D mesh at its base
node, with liquid removed as it forms in place of flow and spraying.

**Spec:** `docs/superpowers/specs/2026-10-11-layered-skin-moving-mesh-design.md` §5.1–5.2, §6.1–6.3, §6.6–6.7, §7, §9
(the rows this stage retires), §10, §11 (build stage 2; level 1's heat checks; level 2 with flow off; level 3's SESAM
verification).

**Builds on:** sub-plan 19 — `MovingMesh.move_with_policy` (each column's draw, absorption and hand-down become the corner
displacement along its normal; `MoveResult.swept_volume` and `swept_enthalpy` are what each column takes from or gives to
the 3D model), `MovingMesh.frame()` (column directions and dual areas), `MovingMesh.smoothing()`, and both backends'
`update_geometry`. The stage-1 skin's conduction solver (`skin.py` in `prototype/work-2026-10-08-skin/code/`, Newton in
temperature with a line search) and its conservative cell operations are reused.

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

- The slurry viscosity law comes from **`reentry_model/shared/rheology.py`** (stage S's **S1**), `material` from
  `reentry_model.shared`. `columns.py` and the rewritten `skin.py` are skyscrapers' own and are **never shared** with
  the level set.

## Tasks (to be written in full)

0. **Volume-neutral smoothing** (`reentry_model/movemesh.py`, decided 2026-10-12): `MovingMesh.smoothing` returns the
   tangential step plus the normal correction that cancels each corner's swept volume, found with `swept_volumes` (one
   correction step, iterated if the remainder is not negligible). Tests: a smoothing-only move on the 100 mm sphere
   books no corner more than a small fraction of its earlier sweep (fact 105: up to about 1 µm of equivalent depth),
   the books stay exact, and sub-plan 19's recession measurement is repeated to show the fold limit and the rebuild
   record are not made worse.
1. **Columns** (`reentry_model/columns.py`): dual cells, areas, face lengths, neighbours, smoothed normals, the surface
   gradient and curvature operators on the corner graph. Tests against a sphere's exact areas, normals and curvature.
2. **Cells with variable counts** (`reentry_model/skin.py`, rewritten): packed per-column arrays of mass, enthalpy and
   velocity (velocity unused until stage 3); target sizes (10 µm fine zone over the minimum, 10 % growth below); split
   above 1.5 × target, merge below 0.5 ×, exact in mass, enthalpy and momentum. Tests: conservation to round-off; the top
   cells untouched by a split or merge deeper down.
3. **The minimum mass** (§6.1): the larger of 400 µm of solid and δ_m plus 50 µm as liquid, per column. Tests: a
   column at its minimum, above it, and with δ_m deeper than 469 µm.
4. **Heat through the depth** (§6.3): the stage-1 implicit solve on variable cell counts, convective heating and
   radiation at the top cell, the base law at the bottom. Tests: steady ablation with liquid removed as it forms
   (recession speed q/(ρ (L + c_p ΔT)) and the exponential profile).
5. **Heat sideways** (§6.3): explicit conduction between neighbouring columns' cells at the same depth below the top,
   stable below Δx²/(4α) (about 25 ms at 2 mm). Tests: a known two-dimensional steady temperature; energy closed.
6. **The base** (§6.6): draw below the minimum, absorb metal turning non-rigid at the base, hand rigid bottom cells
   beyond the minimum to the 3D model — each converted to the corner's displacement for `move_with_policy`; the
   temperature-jump fault (5 K, provisional). Tests: each operation's books with the moving mesh; a drained column
   drawing back to its minimum.
7. **The interface** (§7.2): the trial and real passes at the base nodes, the affine law applied as a nodal load over
   the dual area, the mismatch booked, one repeat above 10⁻³ of the throughput; 3D radiation off on the base surface.
   Tests: both backends agree with the interface active.
8. **The body and the CLI**: `--surface-model skin` with skins on every corner from t = 0 (no creation events),
   `--skin-min-thickness`, `--skin-cell`, `--skin-growth`, `--skin-substep`, `--skin-preset`, `--deform-depth`; run
   names; history columns of §10.2 that exist at this stage. `--surface-model elements` stays bit-identical.
9. **Verification and measurement**: level 2's box with flow off (recession within 5 % at 1.0, 1.5 and 2.0 s of the
   stage-1 fine box and of a new one-dimensional 1 µm reference; mass 1e-12, energy 1e-8); level 3's Step 2 verification
   against SESAM with skins on every corner; **the skins' run time** on the box and a short 100 mm flight — the main
   unknown of the design, measured here before the flow is built on it.

## Inputs to settle before the code-level plan is written

- **The fold limit** found in sub-plan 19 (the fixed 2 mm surface triangulation fails once neighbouring triangles fold
  by more than about 75–95 degrees). Decided 2026-10-11 (Asha, option 1): this stage measures the fold under
  heating-driven recession — the box and the short 100 mm flight of Task 9 report the largest fold every step — and
  local surface refinement is planned only if real flights approach the limit.
- **The tangential smoothing is volume-neutral per corner** (decided 2026-10-12, Asha; spec §5.1): each corner is moved
  along its normal by minus the volume its smoothing swept over its dual area, so a column is never credited with
  volume the smoothing alone swept. This is Task 0 below.
- **The enthalpy of metal coming in on an outward move** (hand-down, freezing): `MovingMesh.move` books it at the
  corner's enthalpy; this stage must pass the skin's bottom-cell enthalpy and put any difference on the corner node.
- The rebuild cost and residual measured in sub-plan 19 (1.6–2.1 s; 0.07–0.27 %, one-signed, offset kept to the zone;
  fact 102): record the residual against each step's heat input in the history, and set the threshold for the
  intersection-based transfer in those terms.
- The mesh work per macro step measured in sub-plan 19 (0.26–0.43 s against 0.2 s of conduction; fact 106), which decides
  whether the zone motion needs optimising before skins are added.
- Integration points the moving mesh leaves: `ThermalBody.i_centre` and the centre temperature after a rebuild (no
  centre node is embedded; the nearest node moved 4.4 mm from the origin in one check), `params["path"]` still naming
  the original mesh file, and the assumption that every element is active.

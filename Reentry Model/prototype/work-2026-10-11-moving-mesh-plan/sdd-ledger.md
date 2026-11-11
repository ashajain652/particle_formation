# SDD ledger — plan: docs/superpowers/plans/melt-spraying-subplans/19-moving-mesh.md

Spec: docs/superpowers/specs/2026-10-11-layered-skin-moving-mesh-design.md (committed cc15b3f, updated 16ba4ab).
Work copy: prototype/work-2026-10-12-moving-mesh/code (nested git repo, branch moving-mesh from 49d88e0; ignored by the
outer repo via .gitignore line added in Task 0 -- that outer edit is uncommitted until the end).
Design ledger (other plan, read-only here): .superpowers/sdd/18-melt-layer-skin/progress.md.

## Pre-flight scan (2026-10-12)

| Pair / task | Produces vs consumes | Found |
|---|---|---|
| T1 -> T3 | corner_frame, orient_outward (T3 tests) | consistent signatures |
| T1 -> T4 | signed_volumes, radius_ratio, nodal_gradient(…, elements), lumped_masses, corner_frame, orient_outward, tangential_smoothing(…, edges), surface_edges, triangle_quality, max_dihedral | consistent; T1 test helper `inward` calls MovingMesh.frame (T4) but no T1 test calls it |
| T2 -> T4 | solver.update_geometry(elements) | consistent (T4 move calls update_geometry(touched)) |
| T2 -> T5 | skfem setup resets AMG; fenicsx setup on new mesh | consistent; T5 rebuild calls solver.setup |
| T3 -> T4 | ZoneMap(points, corners, depth).displacement | consistent |
| T4 -> T5 | MovingMesh.corners/faces_corner/_set_surface/compression/steps_since_rebuild | consistent; T5A is a method appended inside the class, T5B module-level after it |
| T4/T5 -> T6 | move, rebuild, compression, min_gap | consistent; T6 methods inserted after rebuild, before the rebuilding section |
| T1-T7 tests file | appended by section markers "(Task N)" | generator verified reassembly byte-exact |
| fenicsx tests T2/T5/T7 | three functions appended in order | generator verified |
| T7 -> fenicsx T7 | run_moving_front imported from the movemesh test module | consistent |
| T1 | tests vs code | 5 tests, all against T1 functions |
| T2 | tests vs diffs | diffs apply (patch-checked); AMG test has no assert (passes by not raising) |
| T3 | tests vs code | consistent |
| T4 | tests vs code | consistent; tolerances from measured values |
| T5 | tests vs code | validity-gate test parametrised x3 -> +6 tests |
| T6 | tests vs code | consistent |
| T7 | tests only | consistent |
| T8 | harness + docs edits outside the work copy | documents only, allowed by the plan |

Ruling: keep `test_a_solver_set_up_on_a_new_mesh_drops_its_old_amg_hierarchy` without an explicit assert — it is a
regression test whose failure mode is the ValueError the old code raised on the second step; completing the step is
the assertion — costs a weak test if wrong (an explicit `assert result.iterations > 0` could be added later).
Ruling: Task 0 is done by the controller (workspace setup: clone, ignore line, baseline), not by a subagent — it is the
worktree step of this skill, not implementation — costs nothing if wrong.

## Progress
Task 0: complete (work copy at 49d88e0 on branch moving-mesh; baseline 303 passed, 1 skipped, 2 failed + 3 errors known; FEniCSx 14 passed)
Ruling: files end with a single newline instead of the blocks' trailing blank lines; each later append must leave PEP 8's two blank lines between top-level definitions, and the controller checks the final files against prototype/work-2026-10-11-moving-mesh-plan/code byte for byte after Task 7 — costs a whitespace fix if wrong.
Task 1: minor (deferred): orient_outward's flip branch, nodal_gradient(elements=), tangential_smoothing(edges=), lumped_masses not directly tested in Task 1 (plan-mandated; later tasks exercise some).
Task 1: minor (deferred): max_dihedral divides by zero on a zero-area triangle (plan-mandated).
Task 1: minor (deferred): imports used only by later tasks (time, dataclass, mesh_mod; erfc and helpers in tests).
Task 1: complete (commits 49d88e0..bf88a6d, review clean)
Task 2: Ruling: the skfem geometry test recomputes all elements, so the subset path of update_geometry is untested (plan-mandated, Important) — strengthen it: move only the outer nodes, pass just the elements they touch, compare with a fresh solver — costs a slightly longer test if wrong; the work copy then deviates from the plan copy in this test, and the plan copy and sub-plan 19 are synced to the built code at the end.
Task 2: Ruling: FenicsxThermalSolver.cell_volumes goes stale after update_geometry (plan-mandated, Important) — refresh it in update_geometry and assert it in the FEniCSx moved-geometry test — costs two lines if wrong; same sync note.
Task 2: minor (deferred): skfem AMG hierarchy kept across update_geometry until amg_rebuild_every (iteration count only); geometry.x writable relies on dolfinx 0.11; input_global_indices map serial only.
Task 2: fix round 1/5 (2 addressed, 0 open — subset geometry test; FEniCSx cell volumes refreshed via _cell_volumes helper; commits 4a9f78f..c682560)
Task 2: minor (deferred): fenicsx cell_nodes/element_of_cell computed at setup only (fine: connectivity unchanged by moves).
Task 2: complete (commits bf88a6d..c682560, review clean after 1 fix round)
Task 3: minor (deferred): redundant setdiff1d in zone test; compression test checks the median only; ZoneMap.idx/K/w undocumented; empty-interior and fewer-than-12-corners paths untested (all plan-mandated).
Task 3: complete (commits c682560..2985482, review clean)
Task 4: Ruling: plan defect — MovingMesh.__init__ calls default_sizes, which the plan defines only in Task 5, so Task 4 alone raised NameError; the implementer copied default_sizes into the Task 4 section. Keep it for now; Task 5 removes that copy and adds the rebuild section verbatim (where default_sizes belongs), so the module ends identical to the tested copy; the plan text is corrected when the plan copy and sub-plan 19 are synced at the end — costs a misplaced function if wrong.
Task 4: Ruling: reviewer's Important finding (plan-mandated) is real — scaling the corners' first-order swept volumes to the exact total blows them up where they cancel (tangential move: scale factor ~ -7e14), and stage 2 credits each column with its share. Replace with exact per-triangle prism volumes, a third to each corner, plus a bounded remainder spread in proportion to |dV| (tested in the scratch copy: remainder 5e-7 of the total; inward per-corner within 1.0 % of dual area x displacement; tangential per-corner <= 0.6 % of A|u|; moving-front error unchanged 0.0010; a dual-area-times-|u| spread was tried first and broke the moving front — it puts volume on sliding side faces). New tests: per-corner checks in the books and tangential tests, deep nodes unchanged, and an exact-prism box test. Diff: .superpowers/sdd/19-moving-mesh/task-4-fix1.diff — costs one function and three tests if wrong; plan copy synced at the end.
Task 4: fix round 1/5 (2 addressed, 0 open — exact prism swept volumes per corner; per-corner tests; commits 37496d5..164724e)
Task 4: minor (deferred): tangential test bound strict with |tangent|; box test relies on exact zeros; remainder dropped if all dV are zero (cannot occur for a real move); refusal test does not say which gate fired; compression/steps_since_rebuild untested; state half-changed if update_geometry raises.
Task 4: complete (commits 2985482..164724e, review clean after 1 fix round)
Task 5: minor (deferred): rebuild residual test one-sided (abs missing); rebuild tests do not check the preceding move was not refused; no check that the new boundary is exactly the corner set; locate's batched solve raises LinAlgError on a singular candidate and extrapolates unbounded; gmsh errors surface as gmsh exceptions, not RemeshError; default_sizes untested on a real sphere; FEniCSx rebuild test's step asserts nothing; check_surface ignores orientation consistency (all plan-mandated).
Task 5: minor (deferred): implementer commits in the nested work repo carry a "Claude Haiku 5.5" Co-Authored-By trailer (nested repo is outer-ignored; harmless).
Task 5: complete (commits 164724e..31001e9, review clean)
Task 6: Ruling: reviewer's Important finding (plan-mandated) is real — the minimum-gap test is satisfied without the gap clause (threshold and gap met on the same step). Replace it: depth 4 mm, 0.3 mm moves, min_gap 9, ratio_min 1e-6; rebuild expected at call index 9 (index 7 without the gap, both measured in the scratch copy) — costs one test if wrong.
Task 6: fix round 1/5 (1 addressed, 0 open — minimum-gap test now discriminates; commits b8358ca..30463a9)
Task 6: minor (deferred): a compression rebuild lags one call behind the threshold; one call can rebuild twice (compression then quality) and last_rebuild keeps only the second reason; compression accumulates inward motion only; min_gap counts non-zero moves, not macro steps; no test of compression-then-refused path.
Task 6: complete (commits 31001e9..30463a9, review clean after 1 fix round)
Task 7: minor (deferred): front test repeats constants run_moving_front holds; FEniCSx checked only against skfem; sys.path insert inside the FEniCSx test; side-plane property not asserted; single-resolution tolerance.
Task 7: complete (commits 30463a9..82a8867, review clean; measured front error 0.001035, backend gap 4.3e-11 K)
Ruling: Task 8's document steps (facts in 00-shared-context, the results section of sub-plan 19) are done by the controller after the implementer's measurements, because sub-plan 19 is generated from part_19.md and must not be hand-edited; the implementer does the whitespace fix, the harness, the bit-identity flight and the measurement runs, and writes the numbers to its report — costs a slower record if wrong.
Task 8: minor (deferred): harness failure handling covers RuntimeError only and writes JSON at the end; rebuild counter includes the failed move's rebuild; no reproducibility metadata (gmsh threads, commit, power) in the JSON; docstring cosmetics; three blank lines in two places and a trailing blank line (matching the tested copy); report wording on the fold maximum.
Task 8: complete (commits 82a8867..3b3034f, review clean; documents to be written by the controller)
Final review (opus, 49d88e0..3b3034f): Ready with fixes; no Critical.
Ruling: I1 (swept slab booked at the surface enthalpy, a first-order one-signed bias, 0.42 % of swept) — fix now: inward-moving corners book the mean of the old and the remapped pre-offset enthalpy — costs a small change to the books if wrong.
Ruling: I2 (tangential smoothing sweeps second-order volume of either sign per corner; stage 2 cannot credit negative volume to an empty skin) — a spec-level choice (§5.1/5.2 invariants), left for Asha before sub-plan 20 with the reviewer's recommendation (b): make smoothing volume-neutral per corner — costs nothing in stage 1 (no columns yet) if deferred.
Ruling: I3 (rebuild residual one-signed, 99-186 J per rebuild, offset taken from the cold core) — fix now: apply the rebuild offset over the new deforming zone only, and keep a running total `rebuild_residual_total` — costs a different spread of a 0.1-0.3 % residual if wrong; supermesh threshold left to sub-plan 20.
Ruling: I4 (incoming metal on an outward move booked at the corner enthalpy) — deferred to sub-plan 20's interface design, recorded there — costs nothing in stage 1.
Ruling: fix now from the deferred minors: the one-sided rebuild-residual assertion (abs) and asserting the preceding moves were not refused; a runtime check that the rebuilt boundary is exactly the corner set (RemeshError); gmsh generation errors re-raised as RemeshError; the stale module docstring (swept-volume method; one-thread measurements). Leave: locate misses (0.17 K), the wasted second rebuild, performance, Taubin wording (docs), stage-2 integration hazards (to sub-plan 20).
Final fix wave: commit 42da6c3 (F1-F6; movemesh 32, FEniCSx 17, unit tier 335).
Ruling: revert F1 — booking inward slabs at the mean of old and remapped enthalpy shrank the books' correction (-16.2 J -> +0.29 J per 0.3 mm move) but made the exact-solution moving-front test nine times worse (0.0010 -> 0.0093 of the rise): the surface value is consistent with the lumped quadrature, and the correction is not an accuracy measure; the spec also books metal crossing the rigid front at the base node's temperature (§6.6). Keep the surface value, restore the front tolerance 0.005, drop F1's test; done by the same fixer within this wave — costs the reviewer's 0.4 % one-signed correction if wrong (recorded for sub-plan 20).
Final fix wave: complete (commits 3b3034f..d0421e2 plus a18b531 relocation; re-review: all findings addressed; F1 reverted by ruling).
Documents (controller): plan copy synced byte-identical to the build; part_19.md updated (counts, post-review changes, Results section) and sub-plan 19 regenerated (generator checks pass); facts 102-107 and the index line in 00-shared-context; sub-plan 20 inputs; spec §5.2, §5.3, §12 aligned with the build. Outer repo changes uncommitted.
Asha 2026-10-12: keep the branch moving-mesh; volume-neutral smoothing per corner (spec §5.1, sub-plan 20 Task 0); commit everything; bring the other sub-plans up to date or deprecate them.

# Sub-plan 21 — flow in the columns (stage 3 of the layered-skin design)

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or
> superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.
> Read `00-shared-context.md` first, then the spec. **Written at code level: Tasks 1–5, 6a, 6b, 6c, 6d, 7 and 8
> (2026-10-18/24).** Task 9, the base where runoff lands (Asha, 2026-10-25), is outlined at the end and not yet
> written. The outline points to the tasks, and the inputs it listed are settled, each where it says.
> Task 8's coupled pushed plate ran to 1.0 s, not 2.0 s (its ruling 8: about 6.5 hours alone to 2.0 s).

**Goal:** The flowing layer moves: every skin cell carries a velocity in its column's tangent plane, coupled through the
depth by the viscosity at each interface, pushed by the gas pressure gradient, the deceleration and the two thickness
pressures, dragged by the gas shear at the top and held at the rigid front; mass, momentum and heat cross between
neighbouring columns layer by layer.

**Architecture:** Everything happens inside the skin body's existing sub-steps; the 3D model, the moving mesh and the
coupled run are unchanged. Task 1 first makes the sideways conduction implicit: one coupled solve per Newton iteration
of each sub-step, the columns' tridiagonal depth solves as its preconditioner. The flow tasks then give the cells a
velocity, couple it through the depth by an implicit tridiagonal viscous solve per column, and carry mass, enthalpy and
momentum between neighbouring columns in an implicit upwind transport of their own, after the momentum step (Task 6a;
ruling, 2026-10-21).

**Tech Stack:** Python 3.12 (`drama_env`), numpy, scipy (sparse matrices; the conjugate gradients are written out),
scikit-fem for the 3D model, gmsh for its rebuilds.

**Spec:** `docs/superpowers/specs/2026-10-11-layered-skin-moving-mesh-design.md` — §6.3 (heat: the sideways conduction
made implicit, amended 2026-10-17), §6.4 (flow), §6.5 (sideways exchange and the arrival rule), §11 (build stage 3;
level 1's flow and Rayleigh–Taylor checks; level 2 with flow on), §12 (the risks this stage carries: capillary stiffness
at the grid scale, the flowing-layer coordinate, momentum between tangent planes).

**Builds on:** sub-plan 20 — the columns, the packed cells (with their unused velocity), the surface gradient and
curvature operators, and the base operations. Build stage 2 is complete: levels 1, 2 and 3 pass (fact 149).

**How this plan's code is checked.** Every code and diff block below was written and tested in a planning copy,
`prototype/work-2026-10-18-column-flow-plan/code/`, made from sub-plan 20's tested copy,
`prototype/work-2026-10-12-skin-heat-plan/code/`, as committed at 1d01be5 (identical to the skins build's f11e326,
option A included). This document is generated from it by `prototype/work-2026-10-18-column-flow-plan/make_plan21.py`,
which replays the tasks written so far on a copy of the base — every diff through `patch`, each from the state the
tasks before it left (kept, for files a later task changes again, in `states/task<N>/` as the build repository
committed them) — and checks that the result is the tested copy byte for byte. Never edit this document by hand: change `part_21.md` or the tested copy, and
regenerate.

## Global Constraints

- Work in a fresh work copy, `prototype/work-2026-10-18-column-flow/code/`: a git repository of its own, ignored by the
  outer one, cloned from the skins build (`prototype/work-2026-10-12-skin-heat/code/`, branch `skin-heat` at f11e326).
  Commit there after every task. Never edit `prototype/proto3/`, the package `reentry_model/` at the repo root, the
  specs or the generated plans. Step 3 code stays in the prototype until the skyscrapers architecture is verified
  (Asha, 2026-10-12).
- Interpreters: `PY=/Users/ashajain/miniforge3/envs/drama_env/bin/python`;
  `FX=/Users/ashajain/miniforge3/envs/fenicsx_env/bin/python`, run with
  `FI_PROVIDER=tcp CC=/Users/ashajain/miniforge3/envs/fenicsx_env/bin/clang`.
- Unit tier, run in the work copy: `"$PY" -m pytest -m "not drama and not reference" -q --ignore=tests/test_integration_drama.py --ignore=tests/test_sphere_reentry.py --ignore=tests/test_sphere_sweep.py`.
  Baseline at f11e326: 405 passed, 1 skipped. Record the counts after every task; no failure is acceptable.
- FEniCSx tier: `FI_PROVIDER=tcp CC=... "$FX" -m pytest tests/test_reentry_model_fenicsx.py -q`. Baseline 19 passed.
- `--surface-model elements` (the default) is unchanged by this sub-plan.
- Books: in every skin test and run, the mass is closed to 1e-12 and the energy balance residual is below 1e-8.
- Exit codes: a Newton or coupled solve that does not converge raises `RuntimeError` (exit code 1 at the command line);
  bad settings exit 2.
- Long runs under `caffeinate -i`, on mains power, the power source logged at start and end (`pmset -g batt`).
- Quote every path: the repository's path contains a space (`University Documents /MIT`).

## Review Focus

1. **A coupled solve that cannot converge** — it stops the run (exit code 1), never a sub-step taken on an unconverged
   Newton step. Pinned by `test_a_coupled_solve_that_cannot_converge_raises` (Task 1).
2. **Cells that change between sub-steps** (draws, removal, a regrid) — a sub-step's sideways map is rebuilt, never the
   last one's geometry reused. Pinned by
   `test_a_sideways_map_reuses_the_last_geometry_only_while_the_cells_keep_their_thicknesses` (Task 1).
3. **Loosely converged coupled solves** — the inexact Newton steps reach the answer tight ones give. Pinned by
   `test_the_inexact_coupled_solves_reach_the_tight_solve_s_answer` (Task 1).
4. **Neighbouring columns of very different depths** — a melting column beside one at its minimum: the exchange still
   sums to zero and stays symmetric. Pinned by `test_the_sideways_map_sums_to_zero_between_columns_of_any_depth` and
   `test_the_sideways_map_is_symmetric_and_negative_semi_definite` (Task 1).
5. **A sub-step far beyond the old explicit limit** (a finer surface, smaller columns) — bounded, the heat kept. Pinned
   by `test_a_sub_step_far_beyond_the_explicit_limit_stays_bounded_and_keeps_its_heat` (Task 1).
6. **Metal crossing many columns in one sub-step** (Courant numbers near 100 in flight) — positive, never overshooting,
   the books exact. Pinned by `test_courant_numbers_far_above_one_stay_positive_bounded_and_conservative` (Task 6a).
7. **Metal arriving at rigid metal** — it lands on top, never inside, and rigid metal joined to the substrate never
   moves. Pinned by `test_metal_arriving_on_a_rigid_top_lands_on_it_and_the_receiver_s_own_cells_are_untouched`,
   `test_a_floating_crust_carried_across_lands_on_top_and_the_slurry_beneath_it_joins_by_fraction` and
   `test_rigid_metal_joined_to_the_substrate_never_leaves_whatever_its_velocity` (Task 6a).
8. **The transport's upwind order** comes from how scipy numbers strong components, checked every call — were it ever
   wrong, the same answer, more slowly. Pinned by `test_the_upwind_order_changes_the_cost_never_the_answer` (Task 6a).
9. **Surface tension at obtuse corners** — the top surface's own curvature ran backwards there; the thin film's must
   rise with each corner's own height and its product with the faces' exchange have no growing mode. Pinned by
   `test_a_pile_up_lowers_its_own_curvature_at_an_obtuse_corner_of_the_top_but_raises_the_film_s` and
   `test_the_faces_exchange_times_the_film_operator_has_no_growing_mode` (Task 6b).
10. **The predicted heights against the transport's** — the thickness terms are implicit in the heights the thickness
    step predicts, which the transport makes to 0.62 % of the change on the test sphere but only to about 4 % at a
    melting state's rim (6.5 µm of 0.15–0.17 mm, at a Courant number of 1.8; not pinned by a test). Measured in 6c on
    the sphere's draining front: 11–16 % of the largest change at 1 ms, at a few columns of the pile at the shoulder
    (90–91°, 1 mm deep: 3.7·10⁻⁵ m predicted, 1–6·10⁻⁶ m made). Found in 6d: not the mixing's density (the pile is all
    liquid) but where the transport puts arriving metal, at the same fraction of the receiver's flowing depth it left,
    so the donor's fast top layers arrive in the pile's fast top layers and pass on within the sub-step, while 6b's
    column model spreads a column's gain over its whole moving thickness; the miss falls only once the top cells'
    Courant number is well below one (0.16, 0.22 and 0.08 of the largest change at 1, 0.25 and 0.05 ms). Since 6d the
    column model also does not know that landed metal flows on: at the front's new tops it misses by up to 0.59 of the
    change in a sub-step at 1 ms (0.27 at 0.25 ms), where the tops are rigid and only the deceleration's thickness term,
    a few pascals across a face, depends on it. At 10 ms, with the linearisations at their cap, 28 times the change, and
    the correction then adds kinetic energy (Tasks 6c, 6d). Where it matters, a column model weighting each column's
    implicit departures by where its arrivals land, or a second linearisation at the transport's own heights, would
    close it.
11. **Rayleigh–Taylor growth under-resolved in time** — at dt·s near or above one, backward Euler damps a growing mode
    spuriously instead of following it. `ThicknessBook.growth` reports it every sub-step; no test can pin what the
    flow's sub-step must avoid (Task 6b). At 6c's 1 ms it is at most 0.0067 (3.4 mm of liquid at 48 m/s²); since 6d the
    adaptive count keeps it near FLOW_GROWTH, 0.02, where backward Euler's growth rate is within about 1 %, and the
    history's `flow_growth_max` reports each macro step's.
12. **Quiescent layers ringing under surface tension** — two-way flows through the depth join most moving cells into one
    cycle of the transport (Task 6b). Measured in 6c on the 100 mm sphere: 10 315 of 141 000 moving cells in one cycle
    within 12 sub-steps of 1 ms, and with every column melted flows meeting at the rear stagnation point make cycles of
    18 000–57 000 cells, 35–651 s a transport. The iteration from the upwind solve is built in 6d (`_solve`, every cycle
    over 500 cells iterated, the rest factored): 0.10 and 0.28 s against 6c's LU's 0.20 and 0.76 s on the quiescent
    melt, 0.16 and 0.25 s against 4.65 and 25.3 s with every column melted; a cost still, not an error (a quiescent flow
    sub-step 1.8 s against 2.9 s).
13. **A runoff front crossing rigid metal** — 6a's arrival rule held metal landing on a rigid top until the next flow
    sub-step, a front advancing at most one column a flow sub-step (2 m/s on 2 mm columns at 1 ms), which set 6c's ten
    flow sub-steps. Lifted in 6d: the landed metal flows on within the same solve, coasting at the velocity it arrived
    with until the next momentum step drags it on the rigid top beneath. A front of thin deposits therefore runs ahead
    by the drag it would have felt within one flow sub-step: on a strip at 5 ms flow sub-steps the implicit tail reached
    the strip's end and the deposit's centroid was 2 mm ahead of the 0.25 ms run's. At the Courant numbers the adaptive
    count keeps, about 6, the sphere's runoff over its rigid rear is within 6.2 % of the 0.25 ms run's in every band of
    15 degrees. Pinned by
    `test_metal_landing_on_rigid_tops_flows_on_within_the_sub_step_at_the_velocity_it_arrived_with` and
    `test_a_runoff_front_crosses_rigid_metal_faster_than_a_column_a_flow_sub_step_and_converges_to_a_fine_one`.
14. **The flow's sub-step is chosen** — since 6d each heat sub-step from the last one's largest Courant number and
    growth, against 6 and 0.02, between 1 and 40 flow sub-steps (10 ms to 0.25 ms), on one state's study at 93 s (Task
    6d). Watch `flow_cycles_max` at 40 (the Courant number then above its target: `flow_courant_max`) and
    `flow_passes_max` at 8: the thickness step's linearisations at their cap marked the unstable 10 ms runs and reached
    it at 1.5 and 2 ms. A flow that speeds up within a heat sub-step takes one heat sub-step beyond the target before
    the count follows; the first from rest takes 6c's ten. On Task 8's pushed plate (45 um of liquid at 3-4 m/s running
    onto cold metal) the chosen count, at Courant 5.9, leaves where the metal went in 50 ms 25 % (L1) off the 0.25 ms
    run's and its kinetic energy 9 % high, against 2.5 % and 5 % on 6d's sphere: two thirds of it the deposit's coasting
    front (Review Focus 13), 13.7 of 338 mg landing 10-35 mm further; 11 % at 0.5 ms flow sub-steps.
15. **The flow's cost** — 6c: a heat sub-step on the heated melting state cost 48 s (54–59 s re-measured beside 6d under
    load), ten flow sub-steps each about as much as its conduction, the transport three fifths. 6d: 43 s at ten flow
    sub-steps and 30–33 s at the nine the adaptive count chose, the transport under half; one flow sub-step 3.7 s
    against 7.2 s, but on a driven state with small cycles no cheaper (2.67 against 2.63 s), the landed metal's estimate
    and the iteration of a 1 700-cell cycle cancelling the momentum step's saving. A melting part's trial and real
    passes took 124 and 115 s against 6c's 1 113 and 924 s on the heated state, 14 s without the flow: the flow, not the
    trial pass, is the cost, about an hour a melting macro step. The momentum step still runs Task 4's five shear-rate
    passes in a heated state (Review Focus 18). Task 8, on the level-2 plate 170 mm across (16 500 columns) alone: a
    10 ms skin sub-step of the pushed skins 6.7 s with the flow (five flow sub-steps) against 3.5 s without, the
    transport a third of the flow; the trial pass the same. A flow sub-step costs about 0.25 s there however few cells
    flow. The regression's melting macro steps, the liquid released as it forms, cost a quarter to two thirds more than
    the flow-off plate's (under load); the pushed plate's first melting macro step took 5 229 s alone, for Review Focus
    17 and 24.
16. **The transport's iterative solve** — stops at a normwise backward error of 1e-13 or after 40 corrections, then
    takes the LU; every cell is installed from the solve's fluxes, so the books do not depend on the tolerance and the
    masses stay positive at any (the iterates rise to the masses from below). Pinned by
    `test_the_iterative_transport_solve_gives_the_direct_solve_s_answer_where_flows_meet`,
    `test_the_books_close_exactly_and_the_masses_stay_positive_however_loosely_the_iteration_stops` and
    `test_an_iteration_that_cannot_finish_falls_back_to_the_direct_solve`. The 500-cell threshold suits the quiescent
    melt's many mid-sized cycles (0.88 s a transport against 1.45 s at 2 000 cells) and costs a driven state's lone
    1 700-cell cycle 0.1 s (Task 6d).
17. **The real pass's flow counts** — taken from its trial pass, so that the base law is made by the pass it predicts;
    the skins' memory of the flow is part of their state (`Skins.state()` has seven arrays since 6d; code that unpacks
    six must take the first six). Pinned by
    `test_the_real_pass_takes_the_trial_pass_s_counts_and_the_memory_is_part_of_the_skins_state` and
    `test_a_skin_body_s_flow_takes_the_counts_each_trial_pass_chose_in_its_real_pass` (Task 6d). Found in Task 8: the
    trial pass holds its base at the start-of-step temperature, so at a melting onset it can melt less and flow slower
    than its real pass, which then runs at the trial's counts above the Courant target -- one flow sub-step a skin
    sub-step at Courant 18.8 on the pushed plate's third melting part, where 12.6 put the runoff 16 % off on 6d's
    sphere; and a real pass whose cells flow only because its base was held at the 3D model's first answer replays the
    ten its trial pass chose with no memory, 380 flow sub-steps in which nothing moved, 96 s of a 124 s pass. Counts
    taken from the real pass's own memory where it exceeds the target would end both; Task 9 (Asha, 2026-10-25).
18. **Task 4's shear-rate fixed point** — five passes wherever bridge or slurry cells carry a stress on the falling
    branch of Li et al.'s flow curve (between its maximum at about 228 1/s and the 367 1/s at which the law is held):
    there the fixed point has an unstable middle root and drifts away from it, 1.1–1.8 times a pass, to the held rate or
    the lower branch; 60 passes still leave a few cells moving. The cap costs 5.3·10⁻⁴ of the fastest speed at worst and
    1.4·10⁻⁶ on average on the heated state; a liquid layer ends in 1–4 passes. Stopping on the velocities rather than
    the viscosities (Task 4's design), or holding Li's rate at the curve's maximum (spec 6.4's rheology), would end it
    sooner, each pass about 0.14 s of a heated flow sub-step: for Asha to decide (Task 6d).
19. **What a release may take** — liquid first, the exposed slurry only with the flag, never rigid metal, at most what
    the cells hold; a column it empties leaves the flow sub-steps after it and the draw working. Pinned by
    `test_a_release_takes_liquid_first_from_the_top_cells_slurry_only_with_the_flag_and_never_rigid_metal` and
    `test_a_release_that_empties_columns_leaves_the_flow_and_the_draw_working_and_the_books_closed` (Task 7).
20. **The droplet-momentum ratio on unsheared tops** — fact 29's Pi is against the gas shear, which the surface flow
    makes zero on the leeward face, where 6d's runoff runs over the rigid rear at metres a second: Pi and its flag are
    taken over the sheared tops, so a release there (Girin's modes release only windward) is in neither, and the
    history's Pi measures a developing layer below its steady speed early in melting. Pinned by
    `test_the_droplet_momentum_ratio_is_fact_29_s_v_melt_delta_over_nu_on_a_sheared_layer` (Task 7).
21. **Runoff over cold metal freezes only at the skin sub-step's conduction** — the flow sub-steps conduct no heat, so
    liquid landing on cold rigid metal stays liquid and flows on until the next conduction. On Task 8's pushed plate a
    film of 10-20 um moves at 0.6-1.2 m/s on a rigid top, 6-12 mm in a 10 ms skin sub-step, where it would really freeze
    within some 0.03 ms: with 10 ms skin sub-steps 11 % of the deposit lands 5-10 mm further than with 5 or 2.5 ms, which
    agree. Not pinned by a test; on the sphere's rigid rear the deposits would spread that much further than physics does.
    Task 9 (Asha, 2026-10-25).
22. **The release's threshold against stage 2's** — Task 7's release takes liquid at the liquidus, 908 K, stage 2's
    removal at T_feed, 910 K, which keeps a liquid layer about 2 K k_l / q thick on the front: the plates' exact recession
    moves by 1.5-6.1 % (cold) and 5.4-17.6 % (hot), the hot plate's non-rigid depth by 9 %. Level 2's regression with
    the flow on is therefore judged against references with the liquid leaving at the liquidus (`skin_plate.references`,
    `ablation_1d.py --remove-at liquidus`), and meets them as the flow-off plate meets stage 2's (Task 8, item 3).
    Pinned by the references' values in `analysis/skin_plate.py`; the regression itself is an analysis run.
23. **The depth-uniform pressure at short wavelengths** (spec 14 assumption 1) — measured in Task 8
    (`analysis/flow_wavelengths.py`): the share of the flowing layers' volume flux at kh > 0.56 is about 1 % or less on
    the element model's 1.2 mm layers at 105 s and of order a quarter (5-36 % over the first 0.1 s, 78 % at 0.4 s of a
    flux six times smaller) on its 3.4 mm slurry layers at 120 s, where the columns' rates are on average 5-8 % too fast
    at 0.1 s by the spec's inertial factor and 94-201 % by the viscous factor the slurry obeys; on states built from the
    element model's mean depths, moved by flow alone. Whether that is significant, and so whether a further task replaces
    the depth-uniform pressure by a non-hydrostatic one on the cells (Asha, 2026-10-22), is decided once sub-plan 22 has
    measured the share again on the flow-on flight's own states (Asha, 2026-10-25, option a). The measure's definition is pinned by
    `test_a_single_sinusoid_on_a_flat_layer_has_the_share_and_excess_of_its_own_wavenumber` and
    `test_the_share_of_two_waves_is_the_short_one_s_power_over_both_and_a_vector_field_counts_every_component`.
24. **Runoff deposits and the part cap** — deposits make skins thicker than their minimum, their rigid bottoms are
    handed down and the bases move, so a melting part of 0.0625 s can be rejected for its moves after its damped passes
    (fact 126's cost, which the part cap removed): on the coupled pushed plate every part from the fifth was rejected
    after 5-7 passes, run again as three, each of those rejected and run again as two of 0.0104 s, 880-950 s of processor
    time for each 0.0625 s of melting against about 90 s, and its first melting macro step took 5 229 s alone (Task 8,
    item 6). Not pinned; the history's `step_parts` and the interface repeats show it. Task 9 (Asha, 2026-10-25).
25. **Base faults where runoff deposits** — on the coupled pushed plate 1 761 base faults in its first melting macro
    step and 90 columns around the heat's downstream rim and the deposit with a bottom cell 5-400 K above its base node,
    up to 1 253 K, none on the flow-off or released plates and none on the pushed skins alone: the base's reconciliation
    (spec 6.6) putting the move's sweep less the skins' exchange into small bottom cells of columns handing deposits down
    beside columns drawing -- spec 12's buried superheated cells back (assumed from the books, not traced). The energy
    books still close (6.5e-11). Not pinned; for the base's own design to answer before a flow-on flight: Task 9 (Asha,
    2026-10-25, option a, over leaving it to sub-plan 22).

## Amendment of 2026-10-17 — the sideways conduction becomes implicit (Asha's decision)

> Design of record: spec §6.3 as amended on 2026-10-17; facts 147 and 148. **Write and build this task first, before
> the flow**: it is heat only, so the plates verify it with flow off before anything else changes. **Decided
> (Asha, 2026-10-18): design (a)**, measured against (b) on the 100 mm sphere (fact 150); written out as Task 1 below,
> where the steady two-dimensional test is replaced by the map's symmetry and definiteness and the warm spot's
> convergence, which test the same exchange more directly.

**Why.** Sub-plan 20 conducts heat between neighbouring columns explicitly (`skin.sideways_heat`), sub-cycled inside each
sub-step to stay below the explicit stability limit, and the smallest columns set that limit: on the 100 mm sphere's
2 mm surface the columns' dual areas range from 1.26 to 4.59 mm² (median 3.31), so a cycle is 3.8 ms and every 10 ms
sub-step takes three, where the median column alone would need one. On the level-3 flight of the 100 mm sphere (no
melting; 14.6 h, 72 s a macro step) the sideways conduction took three fifths of a step, against an eighth for the
conduction through the depth that the skins exist to resolve (fact 147). Option A (fact 148) removed the work the
sub-cycles repeated without changing a result; what remains is the sub-cycling itself, which an explicit step cannot do
without. Task 6's transport of mass, enthalpy and momentum, layer by layer on the same depth matching, is already
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
along the surface; either then carries Task 6's implicit transport.

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

---

### Task 1: Implicit sideways conduction (design (a))

**Files:**
- Modify: `reentry_model/columns.py` (`Columns.incidence`), `reentry_model/skin.py` (the coupled sub-step and the
  `Sideways` map, which replaces the explicit exchange), `reentry_model/skinbody.py` (the hook that makes the map)
- Test: `tests/test_reentry_model_skin.py` (the sideways tests, Task 5 of sub-plan 20, rewritten)

**Interfaces:**
- Consumes: sub-plan 20's `skin._sideways_geometry`, `skin._interp_weights`, `skin._interp_apply`, `skin.depth_grid`,
  `skin.thomas`, `Columns.edges` and `Columns.conduct_weight`.
- Produces: `skin.Sideways(sk, columns, previous=None)` with `.heat(x) -> (n, K)` [W per cell], linear in the cell
  temperatures x, and `.bins(x) -> (n, bins)`; `skin._conduct(sk, dt, q_top, A_top, T_amb, eps, T_b, sideways=None)`;
  `skin.advance(..., side=...)` with `side(sk, previous) -> Sideways`; `Columns.incidence` (scipy CSR, columns by
  edges); `skin.KRYLOV_TOL`, `skin.KRYLOV_MAX` and the diagnostic counter `skin.KRYLOV`.
- Retires: `skin.sideways_heat`, `skin._sideways_once`, `skin.SIDEWAYS_CFL`, and sub-plan 20's pin of the explicit heat
  to its build (`test_the_sub_cycles_share_their_geometry_and_give_the_heat_as_built_to_the_last_bit`).

The sideways exchange keeps sub-plan 20's construction: each column's cells averaged onto the common depth grid, every
surface edge exchanging k w l (T_j − T_i) in every bin both columns reach, each column's bin heat given back to its
cells by their overlap with the bin. Written as a map on the cells, S = Pᵀ L P (P the averaging onto the bins, L the
edges' conductances, Pᵀ the return to the cells), it is symmetric and negative semi-definite, so it can join the
implicit sub-step: every Newton step solves J dT = −R with J the columns' tridiagonal depth Jacobians less S, by
conjugate gradients preconditioned by the columns' own tridiagonal solves (`_coupled_solve`). Three choices keep it
cheap, each measured on the 100 mm sphere's first six steps (fact 150):

- **Inexact Newton steps.** Each coupled solve stops at `KRYLOV_TOL` = 10⁻³ of its right-hand side and the Newton
  iteration corrects the rest; started from the columns' own solve, it takes about one conjugate-gradient iteration
  (1.13 measured). Solving each step to 10⁻⁶ instead moves no column of the history by more than 6 parts in a billion
  and costs 61 % more.
- **The map is linear**, so the line search's first try takes the sideways heat at T + dT from the solve's own
  products instead of evaluating it again; only a column whose step is cut back costs an evaluation.
- **The geometry is reused** from the last sub-step while the cells keep their thicknesses — in a body that does not
  melt, every sub-step; k is taken afresh at each sub-step's start.

The sensitivity the trial pass tracks for the base law (dT/dT_b) stays through the depth only: the law each column gives
the 3D model is its own, and the neighbours' share is the interface mismatch's to book. Against (b) — the same
exchange split from the depth solve and solved on the bins alone — (a) has no splitting error and costs more: 260 s for
the six steps against 200 s for (b) and 254 s for sub-plan 20's explicit exchange after option A, its processor time
(224 s) a little below the explicit exchange's (235 s).

- [ ] **Step 1: Replace the sideways tests** in `tests/test_reentry_model_skin.py` (everything after `grid_skins`, the
  helper the tests share):

```diff
--- a/tests/test_reentry_model_skin.py
+++ b/tests/test_reentry_model_skin.py
@@ -351,117 +351,141 @@
     return sk, col, pts, edge
 
 
-def test_sideways_heat_sums_to_zero_between_columns_of_any_depth(coarse_sphere_mesh):
+SPOT_SPACE_ERROR = 1.5            # K: the warm spot's error on 2 mm columns from the grid alone (1.29 K measured, 2026-10-18)
+
+
+def random_sphere_skins(coarse_sphere_mesh, seed=4):
+    """Skins on the 4 mm test sphere's corners, each 1-3 times its minimum deep, every cell at its own temperature."""
     from test_reentry_model_columns import sphere_columns
     mat = scheil()
     col, _ = sphere_columns(coarse_sphere_mesh)
     sk = skin.Skins(mat, col.area, skin.SkinSettings())
-    rng = np.random.default_rng(4)
+    rng = np.random.default_rng(seed)
     sk.fill(skin.minimum_mass(mat, sk.area, sk.settings) * rng.uniform(1.0, 3.0, col.n), rng.uniform(500.0, 900.0, col.n))
     sk.regrid()
     sk.T = np.where(sk.valid, sk.T + rng.normal(0.0, 20.0, sk.T.shape), sk.T)          # every cell different
     sk.H = sk.m * mat.enthalpy(sk.T)
-    Q = skin.sideways_heat(sk, col, 0.01)
+    return sk, col, rng
+
+
+def test_the_sideways_map_sums_to_zero_between_columns_of_any_depth(coarse_sphere_mesh):
+    sk, col, _ = random_sphere_skins(coarse_sphere_mesh)
+    Q = skin.Sideways(sk, col).heat(sk.T)
     assert abs(Q.sum()) < 1e-12 * np.abs(Q).sum()
     assert np.all(Q[~sk.valid] == 0.0)
 
 
-def test_sideways_heat_is_exact_for_a_field_linear_along_the_surface():
+def test_the_sideways_map_is_exact_for_a_field_linear_along_the_surface():
     mat = sharp()                                                    # constant k: a linear field carries no net flux
     sk, col, pts, edge = grid_skins(mat)
     T = 600.0 + 2e4 * pts[:, 0] - 1e4 * pts[:, 1]                       # 20 K/mm and -10 K/mm, uniform in depth
     sk.T = np.where(sk.valid, T[:, None], sk.T)
     sk.H = sk.m * mat.enthalpy(sk.T)
-    Q = skin.sideways_heat(sk, col, 1e-4)                            # one cycle: the edge columns do not move first
+    Q = skin.Sideways(sk, col).heat(sk.T)
     assert np.abs(Q[~edge]).max() < 1e-9 * np.abs(Q[edge]).max()
 
 
-def test_conduction_between_columns_reaches_the_known_steady_temperature():
-    """Columns on a flat plate, insulated top and bottom, the boundary columns held at a linear field: the steady state
-    is that linear field everywhere (Laplace's equation, with a constant conductivity)."""
-    mat = sharp()
-    sk, col, pts, edge = grid_skins(mat, n=6)
-    T_edge = 600.0 + 3e4 * pts[:, 0] + 1e4 * pts[:, 1]
-    # 4 ms keeps one sideways cycle per step even at the grid's corner columns (a sixth of the interior's area): with
-    # two, the held columns would drift between cycles and the interior settle 1.8 K off (found while planning)
-    for _ in range(200):                                             # 0.8 s: the slowest mode decays e^-22 over it
-        sk.T[edge] = np.where(sk.valid[edge], T_edge[edge][:, None], 300.0)
-        sk.H = sk.m * mat.enthalpy(sk.T)
-        skin._conduct(sk, 0.004, np.zeros(sk.n), sk.area, 0.0, 0.0, None, skin.sideways_heat(sk, col, 0.004))
-    err = np.abs(sk.T[~edge][sk.valid[~edge]] - np.repeat(T_edge[~edge], sk.count[~edge]))
-    assert err.max() < 1e-9 * (T_edge.max() - T_edge.min())           # 3e-13 K measured: exact
-
-
-def test_a_sub_step_beyond_the_explicit_limit_is_sub_cycled_and_stays_bounded():
+def test_the_sideways_map_is_symmetric_and_negative_semi_definite(coarse_sphere_mesh):
+    """What the conjugate gradients of the implicit sub-step rely on: y . S x = x . S y, x . S x <= 0, and no heat at a
+    temperature uniform over the body."""
+    sk, col, rng = random_sphere_skins(coarse_sphere_mesh)
+    S = skin.Sideways(sk, col)
+    x, y = (np.where(sk.valid, rng.normal(size=sk.T.shape), 0.0) for _ in range(2))
+    yx, xy = float((y * S.heat(x)).sum()), float((x * S.heat(y)).sum())
+    assert abs(yx - xy) < 1e-10 * abs(yx)
+    assert float((x * S.heat(x)).sum()) < 0.0
+    uniform = np.where(sk.valid, 700.0, 0.0)
+    assert np.abs(S.heat(uniform)).max() < 1e-9 * np.abs(S.heat(x)).max() * 700.0
+
+
+def test_a_coupled_sub_step_without_sideways_differences_is_the_uncoupled_one():
+    """Identical columns heated alike carry nothing sideways: the coupled sub-step gives the uncoupled one's
+    temperatures, and its sideways books are zero."""
+    mat = scheil()
+    a, col, _, _ = grid_skins(mat, n=4)
+    b, _, _, _ = grid_skins(mat, n=4)
+    q = np.full(a.n, 2e6)
+    ra = skin._conduct(a, 0.01, q, a.area, 0.0, 0.15, np.full(a.n, 700.0), skin.Sideways(a, col))
+    rb = skin._conduct(b, 0.01, q, b.area, 0.0, 0.15, np.full(b.n, 700.0))
+    assert np.abs(a.T - b.T)[a.valid].max() < 1e-9
+    assert np.abs(ra[3]).max() < 1e-9 * np.abs(ra[0]).max()
+
+
+def test_a_sub_step_far_beyond_the_explicit_limit_stays_bounded_and_keeps_its_heat():
+    """A checkerboard 400 K / 800 K on 1 mm columns over one 0.2 s sub-step, about 30 times the explicit exchange's
+    stability limit there, insulated: implicit, it neither overshoots nor loses heat, and it smooths the board (to
+    622-629 K, measured 2026-10-18). The heat is kept to the Newton iteration's tolerance (3.5e-13 measured)."""
     mat = scheil()
     sk, col, pts, _ = grid_skins(mat, n=6)
-    hot = (np.rint(pts[:, 0] / 1e-3) + np.rint(pts[:, 1] / 1e-3)) % 2 == 0          # a checkerboard
+    hot = (np.rint(pts[:, 0] / 1e-3) + np.rint(pts[:, 1] / 1e-3)) % 2 == 0
     sk.T = np.where(sk.valid, np.where(hot, 800.0, 400.0)[:, None], sk.T)
     sk.H = sk.m * mat.enthalpy(sk.T)
     E0 = sk.energy()
-    Q = skin.sideways_heat(sk, col, 0.2)                             # 0.2 s: about 30 times the limit on 1 mm columns
-    sk.H = sk.H + Q * 0.2
-    T = mat.temperature_from_enthalpy(sk.H[sk.valid] / sk.m[sk.valid])
-    assert sk.energy() == pytest.approx(E0, rel=1e-13)
+    skin._conduct(sk, 0.2, np.zeros(sk.n), sk.area, 0.0, 0.0, None, skin.Sideways(sk, col))
+    T = sk.T[sk.valid]
+    assert sk.energy() == pytest.approx(E0, rel=1e-11)
     assert T.min() >= 400.0 - 1e-6 and T.max() <= 800.0 + 1e-6
-
-
-def _sideways_as_built(sk, columns, dt):
-    """The sideways heat as sub-plan 20 built it (before 2026-10-17): the depth grid, the interpolation weights and the
-    temperatures made again in every sub-cycle, the last included."""
-    def interp(x, xp, fp):
-        n, q = xp.shape
-        off = (np.arange(n) * skin._SPAN)[:, None]
-        i = np.searchsorted((xp + off).ravel(), (x + off).ravel(), side="right").reshape(x.shape) - 1
-        i = np.clip(i - (np.arange(n) * q)[:, None], 0, q - 2)
-        r = np.arange(n)[:, None]
-        x0, x1, f0, f1 = xp[r, i], xp[r, i + 1], fp[r, i], fp[r, i + 1]
-        t = np.clip(np.divide(x - x0, x1 - x0, out=np.zeros_like(x, dtype=float), where=x1 > x0), 0.0, 1.0)
-        return f0 + t * (f1 - f0)
-
-    def once(T, dz):
-        n = T.shape[0]
-        zb = np.concatenate([np.zeros((n, 1)), np.cumsum(dz, axis=1)], axis=1)
-        H = zb[:, -1]
-        zc = np.minimum(skin.depth_grid(sk.settings, float(H.max()))[None, :], H[:, None])
-        Cg = interp(zc, zb, np.concatenate([np.zeros((n, 1)), np.cumsum(T * dz, axis=1)], axis=1))
-        length = np.diff(zc, axis=1)
-        Tbar = np.divide(np.diff(Cg, axis=1), length, out=np.zeros_like(length), where=length > 0.0)
-        i, j = columns.edges[:, 0], columns.edges[:, 1]
-        k = sk.props.k(Tbar)
-        flux = (0.5 * (k[i] + k[j]) * columns.conduct_weight[:, None] * np.minimum(length[i], length[j])
-                * (Tbar[j] - Tbar[i]))
-        Qg = np.zeros_like(length)
-        np.add.at(Qg, i, flux)
-        np.add.at(Qg, j, -flux)
-        return np.diff(interp(zb, zc, np.concatenate([np.zeros((n, 1)), np.cumsum(Qg, axis=1)], axis=1)), axis=1)
-
-    mat, valid, dz = sk.material, sk.valid, sk.dz()
-    wsum = np.bincount(columns.edges.ravel(), weights=np.repeat(columns.conduct_weight, 2), minlength=sk.n)
-    rate = float((mat.k_table.max() * wsum / (mat.rho * mat.cp_table.min() * sk.area)).max())
-    cycles = max(1, int(np.ceil(dt * rate / skin.SIDEWAYS_CFL)))
-    T, H, total = sk.T.copy(), sk.H.copy(), np.zeros_like(sk.T)
-    for _ in range(cycles):
-        Q = np.where(valid, once(T, dz), 0.0)
-        total += Q / cycles
-        if cycles > 1:
-            H = H + Q * dt / cycles
-            T = np.where(valid, mat.temperature_from_enthalpy(np.divide(H, sk.m, out=np.zeros_like(H), where=sk.m > 0.0)), T)
-    return total, cycles
-
-
-def test_the_sub_cycles_share_their_geometry_and_give_the_heat_as_built_to_the_last_bit():
-    """Asha, 2026-10-17: the sub-cycles of a sideways step share its depth grid and interpolation weights, and the
-    temperatures after the last are not computed -- a third of the sideways heat's time on the 100 mm sphere. Columns of
-    different depths, every cell at its own temperature, over several sub-cycles: the heat is the build's exactly."""
-    mat = scheil()
-    sk, col, pts, _ = grid_skins(mat, n=6)
-    rng = np.random.default_rng(7)
-    sk.fill(skin.minimum_mass(mat, sk.area, sk.settings) * rng.uniform(1.0, 2.5, col.n), rng.uniform(450.0, 850.0, col.n))
-    sk.regrid()
-    sk.T = np.where(sk.valid, sk.T + rng.normal(0.0, 15.0, sk.T.shape), sk.T)
-    sk.H = sk.m * mat.enthalpy(sk.T)
-    for dt in (1e-4, 0.05):
-        reference, cycles = _sideways_as_built(sk, col, dt)
-        np.testing.assert_array_equal(skin.sideways_heat(sk, col, dt), reference)
-    assert cycles >= 3                                               # the second dt is sub-cycled
+    assert T.max() - T.min() < 20.0
+
+
+def test_a_warm_spot_spreads_as_the_heat_equation_says_and_converges_with_the_sub_step():
+    """A Gaussian warm spot (100 K, 4 mm) in a flat array of 2 mm columns, uniform in depth and insulated, spreads over
+    0.2 s as the two-dimensional heat equation's exact solution, to the grid's own spatial error; against the same
+    grid at a 2.5 ms sub-step the time error falls in proportion to the sub-step (backward Euler, first order: the
+    errors at 40, 20 and 10 ms stand as 37.5 : 17.5 : 7.5, their distances from the reference's 2.5 ms; measured 2.79,
+    1.32 and 0.57 K, 2026-10-18)."""
+    mat = sharp()
+    alpha = 100.0 / (2800.0 * 1000.0)
+    t0, T_end = (4e-3) ** 2 / (4.0 * alpha), 0.2
+
+    def run(h):
+        sk, col, pts, _ = grid_skins(mat, n=21, h=2e-3, T=500.0)
+        c = pts[:, :2] - pts[:, :2].mean(axis=0)
+        r2 = (c * c).sum(axis=1)
+        sk.T = np.where(sk.valid, (500.0 + 100.0 * np.exp(-r2 / (4.0 * alpha * t0)))[:, None], sk.T)
+        sk.H = sk.m * mat.enthalpy(sk.T)
+        for _ in range(int(round(T_end / h))):
+            skin._conduct(sk, h, np.zeros(sk.n), sk.area, 0.0, 0.0, None, skin.Sideways(sk, col))
+        exact = 500.0 + 100.0 * t0 / (t0 + T_end) * np.exp(-r2 / (4.0 * alpha * (t0 + T_end)))
+        return sk.T[:, 0], exact, r2 < (12e-3) ** 2
+
+    ref, exact, inner = run(0.0025)
+    err = {h: np.abs(run(h)[0] - ref)[inner].max() for h in (0.04, 0.02, 0.01)}
+    assert err[0.04] / err[0.02] == pytest.approx(37.5 / 17.5, rel=0.1)
+    assert err[0.02] / err[0.01] == pytest.approx(17.5 / 7.5, rel=0.1)
+    assert np.abs(ref - exact)[inner].max() < SPOT_SPACE_ERROR
+
+
+def test_a_sideways_map_reuses_the_last_geometry_only_while_the_cells_keep_their_thicknesses(coarse_sphere_mesh):
+    sk, col, _ = random_sphere_skins(coarse_sphere_mesh)
+    first = skin.Sideways(sk, col)
+    assert skin.Sideways(sk, col, first).geometry is first.geometry
+    sk.add_bottom(np.full(sk.n, 1e-7), np.full(sk.n, float(sk.material.enthalpy(700.0))), np.zeros(sk.n))
+    sk.regrid()
+    again, fresh = skin.Sideways(sk, col, first), skin.Sideways(sk, col)
+    assert again.geometry is not first.geometry
+    np.testing.assert_array_equal(again.heat(sk.T), fresh.heat(sk.T))
+
+
+def test_the_inexact_coupled_solves_reach_the_tight_solve_s_answer(coarse_sphere_mesh, monkeypatch):
+    """Columns of different depths melting under 3 MW/m2 on the front, every cell at its own temperature: each Newton
+    step solved to KRYLOV_TOL gives the temperatures that steps solved to 1e-10 give, to the Newton tolerance."""
+    out = []
+    for tol in (skin.KRYLOV_TOL, 1e-10):
+        monkeypatch.setattr(skin, "KRYLOV_TOL", tol)
+        sk, col, _ = random_sphere_skins(coarse_sphere_mesh)
+        q = 3e6 * np.clip(col.normal[:, 0], 0.0, None)
+        books = skin._conduct(sk, 0.01, q, sk.area, 0.0, 0.15, np.full(sk.n, 850.0), skin.Sideways(sk, col))
+        assert abs(books[3].sum()) < 1e-12 * np.abs(books[3]).sum()            # the sideways books sum to zero
+        out.append(sk.T[sk.valid].copy())
+    assert np.abs(out[0] - out[1]).max() < 1e-5
+
+
+def test_a_coupled_solve_that_cannot_converge_raises(coarse_sphere_mesh, monkeypatch):
+    """A Newton step whose conjugate gradients run out of iterations stops the run (exit code 1 from the command
+    line), never a sub-step taken on an unconverged step."""
+    monkeypatch.setattr(skin, "KRYLOV_MAX", 0)
+    monkeypatch.setattr(skin, "KRYLOV_TOL", 1e-14)
+    sk, col, _ = random_sphere_skins(coarse_sphere_mesh)
+    with pytest.raises(RuntimeError, match="coupled sideways solve did not converge"):
+        skin._conduct(sk, 0.01, np.zeros(sk.n), sk.area, 0.0, 0.0, np.full(sk.n, 700.0), skin.Sideways(sk, col))
```

- [ ] **Step 2: Run them to see them fail.** `"$PY" -m pytest tests/test_reentry_model_skin.py -q -k "sideways or
  coupled or warm_spot or explicit_limit"`. Expected: `AttributeError: module 'reentry_model.skin' has no attribute
  'Sideways'`.

- [ ] **Step 3: The edges' incidence matrix** in `reentry_model/columns.py`:

```diff
--- a/reentry_model/columns.py
+++ b/reentry_model/columns.py
@@ -10,12 +10,14 @@
 with a and b the angles opposite the edge in its two triangles: the flux of a piecewise-linear field through the two
 dual segments across the edge is exactly w (T_j - T_i) per unit conductivity and depth, so the discrete operator is
 exact for a field linear along the surface. A weight is clipped at zero where the opposite angles are obtuse enough
-to make it negative: a negative conductance would break the explicit step's positivity, and it is rare (3 of 6 756
-edges on the 4 mm test sphere, 4 of 28 245 on the default 2 mm sphere, measured 2026-10-12).
+to make it negative: a negative conductance would break the explicit step's positivity (and, since sub-plan 21 makes the
+sideways conduction implicit, the definiteness its conjugate gradients need), and it is rare (3 of 6 756 edges on the
+4 mm test sphere, 4 of 28 245 on the default 2 mm sphere, measured 2026-10-12).
 The surface gradient and curvature operators the flow needs are added in sub-plan 21, the first stage that uses them."""
 from dataclasses import dataclass
 
 import numpy as np
+import scipy.sparse as sp
 
 from .movemesh import corner_frame
 
@@ -82,6 +84,17 @@
         np.add.at(out, j, -d)
         return out
 
+    @property
+    def incidence(self):
+        """B (n columns x E edges), +1 at each edge's first column and -1 at its second, so B @ f gives every column the
+        sum over its edges of f, with the first column's sign. Built once: the connectivity is fixed."""
+        if getattr(self, "_incidence", None) is None:
+            i, j = self.edges[:, 0], self.edges[:, 1]
+            e = np.arange(len(i))
+            self._incidence = sp.csr_matrix((np.r_[np.ones(len(i)), -np.ones(len(i))], (np.r_[i, j], np.r_[e, e])),
+                                            shape=(self.n, len(i)))
+        return self._incidence
+
 
 @dataclass
 class CornerSurface:
```

- [ ] **Step 4: The coupled sub-step and the sideways map** in `reentry_model/skin.py`:

```diff
--- a/reentry_model/skin.py
+++ b/reentry_model/skin.py
@@ -47,7 +47,10 @@
 NEWTON_MAX = 40
 LINE_SEARCH_MAX = 12              # halvings of a Newton step before it is taken anyway
 DZ_FLOOR = 1.0e-12                # m: the conduction length of a vanishing cell
-SIDEWAYS_CFL = 0.9                # explicit sideways conduction: sub-cycles keep dt k sum(w) / (rho c_p A) below this
+KRYLOV_TOL = 1.0e-3               # each Newton step's coupled solve (sideways conduction in the implicit sub-step) stops
+                                  # at this residual relative to its right-hand side: the Newton iteration corrects the rest
+KRYLOV_MAX = 200                  # conjugate-gradient iterations of one coupled solve at most
+KRYLOV = {"solves": 0, "iterations": 0}   # the coupled solves made and their iterations, for measuring the cost
 TABLE_DT = 0.01                   # K, spacing of the skins' property tables (PropertyTable)
 TABLE_RANGE = (150.0, 3000.0)     # K, their extent; beyond it the end intervals are extended linearly
 PROFILE_POINTS = 8                # Gauss-Legendre points of an enthalpy averaged over a depth of the 3D profile
@@ -459,14 +462,53 @@
     drawn_H: np.ndarray      # the enthalpy it brought
 
 
-def _conduct(sk, dt, q_top, A_top, T_amb, eps, T_b, Q_side=None):
+def _coupled_solve(G, diag, rhs, valid, sideways):
+    """One Newton step of the coupled sub-step: J x = rhs with J the columns' tridiagonal Jacobians less the sideways
+    map (`Sideways.heat`, linear), symmetric positive definite, by conjugate gradients preconditioned by the columns' own
+    tridiagonal solves (`thomas`) and started from them, to KRYLOV_TOL of the right-hand side. Returns x and the sideways
+    heat at x, built from the products the iteration made (the map is linear), which the line search reuses."""
+    def depth(x):
+        y = diag * x
+        y[:, :-1] -= G * x[:, 1:]
+        y[:, 1:] -= G * x[:, :-1]
+        return y
+
+    def solve_columns(r):
+        return np.where(valid, thomas(-G, diag, -G, r), 0.0)
+
+    x = solve_columns(rhs)
+    Hx = sideways.heat(x)
+    r = np.where(valid, rhs - depth(x) + Hx, 0.0)
+    target = KRYLOV_TOL * math.sqrt(float((rhs * rhs).sum()))
+    z = solve_columns(r)
+    p = z.copy()
+    rz = float((r * z).sum())
+    KRYLOV["solves"] += 1
+    for it in range(KRYLOV_MAX + 1):
+        if math.sqrt(float((r * r).sum())) <= target:
+            KRYLOV["iterations"] += it
+            return x, Hx
+        Hp = sideways.heat(p)
+        Jp = np.where(valid, depth(p) - Hp, 0.0)
+        a = rz / float((p * Jp).sum())
+        x, Hx, r = x + a * p, Hx + a * Hp, r - a * Jp
+        z = solve_columns(r)
+        rz, rz_old = float((r * z).sum()), rz
+        p = z + (rz / rz_old) * p
+    raise RuntimeError("a skin sub-step's coupled sideways solve did not converge in {} iterations".format(KRYLOV_MAX))
+
+
+def _conduct(sk, dt, q_top, A_top, T_amb, eps, T_b, sideways=None):
     """One implicit (backward-Euler) sub-step for every column: the convective flux q_top [W/m2] over the top area
     A_top in at the top cell, radiation from it over the same area, conduction through the cells (harmonic
-    conductances), the explicit sideways heat Q_side [W per cell], and the base held at T_b [K] through the bottom
-    half-cell (T_b None: an adiabatic base, for the tests). Newton in the temperature with a backtracking line search
-    per column (stage 1's solver, sub-plan 18). Returns (E_top, E_rad, E_base, E_side) [J per column] and, when the
-    skins track S, updates it to dT/dT_b and returns the base energy's derivative with respect to T_b as a fifth item.
-    E_base is what the books leave, so they close exactly whatever the iteration error."""
+    conductances), the sideways conduction `sideways` (a `Sideways` map made at the sub-step's start, implicit with the
+    rest: sub-plan 21, design (a)), and the base held at T_b [K] through the bottom half-cell (T_b None: an adiabatic
+    base, for the tests). Newton in the temperature with a backtracking line search per column (stage 1's solver,
+    sub-plan 18); with the sideways map the columns are coupled and each Newton step is solved by `_coupled_solve`.
+    Returns (E_top, E_rad, E_base, E_side) [J per column] and, when the skins track S, updates it to dT/dT_b and returns
+    the base energy's derivative with respect to T_b as a fifth item -- through the depth only: the base law each column
+    gives the 3D model is its own, and what the neighbours' share misses is the interface mismatch's to book. E_base is
+    what the books leave, so they close exactly whatever the iteration error."""
     mat = sk.props
     n, K = sk.m.shape
     valid = sk.valid
@@ -482,14 +524,17 @@
     q = np.asarray(q_top, dtype=float) * np.asarray(A_top, dtype=float)
     base = T_b is not None
     T_b = np.asarray(T_b if base else 0.0, dtype=float) * np.ones(n)
-    Qs = np.zeros((n, K)) if Q_side is None else np.where(valid, Q_side, 0.0)
-
-    def system(T):
+    side = (lambda T: np.zeros((n, K))) if sideways is None else sideways.heat
+
+    def system(T, Q=None):
+        """The residual, the conductances and the Jacobian's diagonal at T, and the sideways heat Q at T (given when the
+        caller has it already)."""
+        Q = side(T) if Q is None else Q
         h, cp, k = mat.conduction(T)
         G = np.where(link, A[:, None] / (0.5 * dz[:, :-1] / k[:, :-1] + 0.5 * dz[:, 1:] / k[:, 1:]), 0.0)
         Gb = A * k[rows, last] / (0.5 * dz[rows, last]) if base else np.zeros(n)
         flux = G * (T[:, :-1] - T[:, 1:])
-        net = Qs.copy()
+        net = Q.copy()
         net[:, :-1] -= flux
         net[:, 1:] += flux
         net[:, 0] += q - eps_s * (T[:, 0] ** 4 - T_amb ** 4)
@@ -500,27 +545,31 @@
         diag[:, 1:] += G
         diag[rows, last] += Gb
         R = np.where(valid, (m * h - E_old) / dt - net, 0.0)
-        return R, G, Gb, np.where(valid, diag, 1.0)
+        return R, G, Gb, np.where(valid, diag, 1.0), Q
 
     # Newton in the temperature with a backtracking line search per column on the residual's norm: in 10-micrometre
     # cells conduction outweighs the heat capacity some five hundred times over a 10 ms sub-step, so a cell must move
     # with its neighbours, and the 3D solver's per-node enthalpy-consistent update stalled here (sub-plan 18)
     T = T_old.copy()
-    R, G, Gb, diag = system(T)
+    R, G, Gb, diag, Q = system(T)
     moved = np.zeros(n)
     for _ in range(NEWTON_MAX):
-        dT = thomas(-G, diag, -G, -R)
+        if sideways is None:
+            dT, Q_dT = thomas(-G, diag, -G, -R), None
+        else:
+            dT, Q_dT = _coupled_solve(G, diag, -R, valid, sideways)
         r0 = np.sqrt((R * R).sum(axis=1))
         lam = np.ones(n)
         for _ in range(LINE_SEARCH_MAX):
             T_try = T + lam[:, None] * dT
-            R_try, G_try, Gb_try, diag_try = system(T_try)
+            full = Q_dT is not None and bool(np.all(lam == 1.0))         # the map is linear: Q(T + dT) = Q(T) + Q(dT)
+            R_try, G_try, Gb_try, diag_try, Q_try = system(T_try, Q + Q_dT if full else None)
             moved = np.abs(T_try - T).max(axis=1)
             worse = (np.sqrt((R_try * R_try).sum(axis=1)) > (1.0 - 1e-4 * lam) * r0) & (moved > NEWTON_TOL)
             if not worse.any():
                 break
             lam = np.where(worse, 0.5 * lam, lam)
-        T, R, G, Gb, diag = T_try, R_try, G_try, Gb_try, diag_try
+        T, R, G, Gb, diag, Q = T_try, R_try, G_try, Gb_try, diag_try, Q_try
         if moved.max() <= NEWTON_TOL:
             break
     else:
@@ -530,7 +579,7 @@
     E_new = m * mat.enthalpy(T)
     E_top = q * dt
     E_rad = eps_s * (T[:, 0] ** 4 - T_amb ** 4) * dt
-    E_side = Qs.sum(axis=1) * dt
+    E_side = Q.sum(axis=1) * dt
     E_base = E_top - E_rad + E_side - (E_new.sum(axis=1) - E_old.sum(axis=1))
     out = [E_top, E_rad, E_base, E_side]
     if sk.S is not None:
@@ -544,8 +593,10 @@
 
 
 def advance(sk, dt, q_top, A_top, T_amb, eps, T_b0, m_min, T_b1=None, side=None, remove=True, g_base=None):
-    """The macro step's sub-steps for every column (spec 7.2). Each sub-step: the sideways heat `side(sk, h)` [W per
-    cell] if given (Task 5), the implicit conduction, then -- with `remove` and a melting material -- the contiguous top
+    """The macro step's sub-steps for every column (spec 7.2). Each sub-step: the implicit conduction, through the depth
+    and, with `side`, sideways (`side(sk, previous)` makes the sub-step's `Sideways` map from the last one, Task 5),
+    then -- with `remove` and a
+    melting material -- the contiguous top
     cells at or above T_feed leave, a column below its minimum mass m_min [kg] draws back to it, and the cells are split
     and merged. The draw takes rigid metal (spec 6.6): at the base temperature, but never above T_rigid -- a 3D base
     node can overshoot it within a macro step (on 4 mm elements under 6 MW/m2 it reached 930 K in one 0.5 s step, found
@@ -573,9 +624,10 @@
     s_b = (T_b < mat.T_rigid).astype(float)           # drawn metal follows the base only where the draw is not capped
     drawn_depth = np.zeros(sk.n)                      # m, how far down the 3D profile this pass has drawn
     melting = remove and mat.melts
+    sideways = None
     for _ in range(N):
-        Q_side = side(sk, h) if side is not None else None
-        res = _conduct(sk, h, q_top, A_top, T_amb, eps, T_b, Q_side)
+        sideways = side(sk, sideways) if side is not None else None
+        res = _conduct(sk, h, q_top, A_top, T_amb, eps, T_b, sideways)
         book.E_top += res[0]
         book.E_rad += res[1]
         book.E_base += res[2]
@@ -613,8 +665,8 @@
 
 def _interp_weights(x, xp):
     """Row-wise linear interpolation at x (n, p) of a function given at non-decreasing xp (n, q), held at the ends, as
-    (rows, index, fraction) for `_interp_apply`: they depend on the abscissae only, so the sub-cycles of a sideways step,
-    whose cells keep their thicknesses, share them."""
+    (rows, index, fraction) for `_interp_apply`: they depend on the abscissae only, so a sub-step's sideways map makes
+    them once and applies them at every iteration of its solve."""
     n, q = xp.shape
     off = (np.arange(n) * _SPAN)[:, None]
     i = np.searchsorted((xp + off).ravel(), (x + off).ravel(), side="right").reshape(x.shape) - 1
@@ -642,9 +694,9 @@
 
 
 def _sideways_geometry(sk, columns, dz):
-    """What every sub-cycle of a sideways step shares, from the cells' thicknesses dz alone: the cells' boundaries from
-    the top, each column's share of each bin of the depth grid, the depth each edge's two columns share in each bin,
-    and the interpolation weights from cells to bins and back."""
+    """What a sub-step's sideways map is built on, from the cells' thicknesses dz alone: the cells' boundaries from the
+    top, each column's share of each bin of the depth grid, the depth each edge's two columns share in each bin, and the
+    interpolation weights from cells to bins and back."""
     n = dz.shape[0]
     zb = np.concatenate([np.zeros((n, 1)), np.cumsum(dz, axis=1)], axis=1)          # cell boundaries, from the top
     H = zb[:, -1]
@@ -656,54 +708,43 @@
             "to_bins": _interp_weights(zc, zb), "to_cells": _interp_weights(zb, zc)}
 
 
-def _sideways_once(sk, columns, T, dz, geometry=None):
-    """Sideways heat [W per cell] at temperatures T: each column's cells averaged onto the depth grid, every surface
-    edge exchanging k w l (T_j - T_i) in every bin both columns reach (k the mean of the two bins' conductivities, w the
-    edge's cotangent weight, l the depth the two columns share in the bin), and each column's bin heat given back to its
-    cells in proportion to their overlap with the bin. Antisymmetric edge by edge, so the total is zero to round-off;
-    exact for a field linear along the surface where neighbouring columns have the same cells. `geometry`: the step's
-    `_sideways_geometry`, made here if not given."""
-    n = T.shape[0]
-    g = _sideways_geometry(sk, columns, dz) if geometry is None else geometry
-    C = np.concatenate([np.zeros((n, 1)), np.cumsum(T * dz, axis=1)], axis=1)      # integral of T from the top
-    Cg = _interp_apply(g["to_bins"], C)
-    length = g["length"]
-    Tbar = np.divide(np.diff(Cg, axis=1), length, out=np.zeros_like(length), where=length > 0.0)
-    i, j = columns.edges[:, 0], columns.edges[:, 1]
-    k = sk.props.k(Tbar)
-    flux = 0.5 * (k[i] + k[j]) * columns.conduct_weight[:, None] * g["shared"] * (Tbar[j] - Tbar[i])
-    Qg = np.zeros_like(length)
-    np.add.at(Qg, i, flux)
-    np.add.at(Qg, j, -flux)
-    Qcum = np.concatenate([np.zeros((n, 1)), np.cumsum(Qg, axis=1)], axis=1)
-    return np.diff(_interp_apply(g["to_cells"], Qcum), axis=1)
-
-
-def sideways_heat(sk, columns, dt):
-    """The sideways conduction of one sub-step dt [s] as an explicit heat source [W per cell], the time average over
-    as many sub-cycles as keep dt k_max sum(w) / (rho c_p,min A) below SIDEWAYS_CFL (spec 6.3); the cells' temperatures
-    move between sub-cycles, the source being what `_conduct` then adds. The smallest columns set the count: on the
-    100 mm sphere's 2 mm surface, dual areas of 1.26-4.59 mm2 (median 3.31) give a cycle of 3.8 ms, three to a 10 ms
-    sub-step, where the median column alone would need one (measured 2026-10-17). The sub-cycles share the step's
-    geometry, and the temperatures after the last are not computed (Asha, 2026-10-17: remaking the one in every
-    sub-cycle and computing the other took a third of the sideways heat's time; the result is the same to the last
-    bit)."""
-    mat = sk.material
-    valid = sk.valid
-    dz = sk.dz()
-    wsum = np.bincount(columns.edges.ravel(), weights=np.repeat(columns.conduct_weight, 2), minlength=sk.n)
-    rate = float((mat.k_table.max() * wsum / (mat.rho * mat.cp_table.min() * sk.area)).max())
-    cycles = max(1, int(math.ceil(dt * rate / SIDEWAYS_CFL)))
-    geometry = _sideways_geometry(sk, columns, dz)
-    T, H = sk.T.copy(), sk.H.copy()
-    total = np.zeros_like(T)
-    for c in range(cycles):
-        Q = np.where(valid, _sideways_once(sk, columns, T, dz, geometry), 0.0)
-        total += Q / cycles
-        if c < cycles - 1:
-            H = H + Q * dt / cycles
-            T = np.where(valid, mat.temperature_from_enthalpy(np.divide(H, sk.m, out=np.zeros_like(H), where=sk.m > 0.0)), T)
-    return total
+class Sideways:
+    """The sideways conduction of one sub-step as a linear map on the cells' temperatures (spec 6.3; implicit since
+    sub-plan 21, design (a), Asha 2026-10-18). Each column's cells are averaged onto the common depth grid, every surface
+    edge exchanges c (x_j - x_i) in every bin both columns reach, with c = k w l (k the mean of the two bins'
+    conductivities at the sub-step's start, w the edge's cotangent weight, l the depth the two columns share in the bin),
+    and each column's bin heat goes back to its cells in proportion to their overlap with the bin. `heat(x)` [W per cell]
+    is antisymmetric edge by edge, so it sums to zero to round-off, and exact for a field linear along the surface where
+    neighbouring columns have the same cells. As a map it is symmetric and negative semi-definite -- the bins' averaging
+    and the heat's return to the cells are each other's transposes -- so `_conduct` takes it into its implicit solve
+    (`_coupled_solve`), and the sub-step has no sideways stability limit: the explicit exchange it replaces needed three
+    sub-cycles a 10 ms sub-step on the 100 mm sphere's smallest columns (facts 147-148)."""
+
+    def __init__(self, sk, columns, previous=None):
+        """`previous`: the last sub-step's map, whose geometry is reused while the cells keep their thicknesses (in a
+        body that does not melt, all of its sub-steps)."""
+        self.valid, self.dz = sk.valid, sk.dz()
+        same = previous is not None and previous.dz.shape == self.dz.shape and np.array_equal(previous.dz, self.dz)
+        self.geometry = previous.geometry if same else _sideways_geometry(sk, columns, self.dz)
+        self.i, self.j = columns.edges[:, 0], columns.edges[:, 1]
+        k = sk.props.k(self.bins(sk.T))
+        self.c = 0.5 * (k[self.i] + k[self.j]) * columns.conduct_weight[:, None] * self.geometry["shared"]
+        self.B = columns.incidence
+
+    def bins(self, x):
+        """The mean of x (per cell) over each bin of the depth grid, (n, bins); zero in bins beyond a column."""
+        n = x.shape[0]
+        g = self.geometry
+        C = np.concatenate([np.zeros((n, 1)), np.cumsum(x * self.dz, axis=1)], axis=1)      # integral from the top
+        return np.divide(np.diff(_interp_apply(g["to_bins"], C), axis=1), g["length"], out=np.zeros_like(g["length"]),
+                         where=g["length"] > 0.0)
+
+    def heat(self, x):
+        """The heat [W per cell] the neighbouring columns conduct into each cell at cell temperatures x (linear in x)."""
+        xb = self.bins(x)
+        Qg = self.B @ (self.c * (xb[self.j] - xb[self.i]))
+        Qcum = np.concatenate([np.zeros((x.shape[0], 1)), np.cumsum(Qg, axis=1)], axis=1)
+        return np.where(self.valid, np.diff(_interp_apply(self.geometry["to_cells"], Qcum), axis=1), 0.0)
 
 
 # -- the base (Task 7) -------------------------------------------------------------------------------------------------
```

- [ ] **Step 5: The skin body's hook** in `reentry_model/skinbody.py`:

```diff
--- a/reentry_model/skinbody.py
+++ b/reentry_model/skinbody.py
@@ -208,7 +208,7 @@
         A_top = self.surface.areas
         eps = np.full(sk.n, self.emissivity)
         m_min = skin_mod.minimum_mass(mat, sk.area, self.settings)
-        side = lambda s, h: skin_mod.sideways_heat(s, self.columns, h)
+        side = lambda s, previous: skin_mod.Sideways(s, self.columns, previous)   # implicit (sub-plan 21)
         nodes = mm.corners
         T_b0 = self.solver.temperature()[nodes]
         g0 = self._normal_gradient()                         # the profile's slope below every base, for this part
```

- [ ] **Step 6: Run them to see them pass.** Expected: `tests/test_reentry_model_skin.py` and
  `tests/test_reentry_model_columns.py` 31 passed (the warm spot's grid error 1.29 K; its time errors at 40, 20 and
  10 ms 2.79, 1.32 and 0.57 K against a 2.5 ms run, as first order predicts; the checkerboard smoothed to 622–629 K in
  one 0.2 s sub-step with its heat kept to 3.5·10⁻¹³). Then the unit tier: 409 passed, 1 skipped; and the FEniCSx tier:
  19 passed.

- [ ] **Step 7: Level 2 with the implicit exchange.** Run the four plates, flow off, side by side:
  `"$PY" analysis/skin_plate.py --case cold|hot --dt 0.05|0.5 --out <json>`. Each must pass as it did with the
  explicit exchange (spec §11.1: recession within 5 % of the references at 1.5 and 2.0 s, and at 1.0 s for the cold
  plate; for the hot plate the flowing depth within 10 % and the base within 1 K of T_rigid; books to 10⁻¹² and 10⁻⁸).

  Measured (2026-10-18, the four run side by side, each once): all four pass, and the implicit exchange moves no
  judged number by more than 0.2 points from sub-plan 20's explicit exchange (B′ with the look-ahead, facts 145-146).
  The table gives the centre's recession against the one-dimensional reference at 1.0, 1.5 and 2.0 s (negative is
  late), with the explicit exchange's in brackets:

| Plate, macro step | 1.0 s | 1.5 s | 2.0 s | Other criteria |
|---|---|---|---|---|
| Cold, 0.05 s | +3.7 % [+3.6] | +0.2 % [0.0] | −0.9 % [−0.9] | fine box +1.7, −0.5, −0.6 % [+1.6, −0.6, −0.6] |
| Cold, 0.5 s | +3.2 % [+3.2] | +0.5 % [+0.4] | −0.9 % [−0.9] | fine box +1.2, −0.2, −0.7 % [+1.2, −0.2, −0.6] |
| Hot, 0.05 s | −3.1 % (not judged) [−3.1] | −1.5 % [−1.5] | −0.8 % [−0.8] | flowing depth within 0.9 %, base at most 0.08 K above T_rigid |
| Hot, 0.5 s | −9.0 % (not judged) [−9.0] | −3.0 % [−2.8] | −1.1 % [−1.1] | flowing depth within 1.1 %, base at most 0.09 K above T_rigid |

  Mass closed to 4·10⁻¹⁶ and energy to 2.2·10⁻¹⁰ (the cold plate at 0.05 s) or better; no base faults; 84, 76, 12 and
  12 interface repeats. The plates' heated region tapers along a half cosine from 32 to 45 mm, so the implicit exchange
  carries heat at its rim, but the judged centre is one-dimensional: these runs show the change leaves level 2 where
  it was, not that the exchange is more accurate. Their wall times (598–1 546 s) were taken four at a time and are not
  comparable with sub-plan 20's, timed alone.

- [ ] **Step 8: Commit.**

```bash
git add reentry_model/columns.py reentry_model/skin.py reentry_model/skinbody.py tests/test_reentry_model_skin.py
git commit -m "Skins: the sideways conduction is implicit, solved with each sub-step's depth conduction by conjugate gradients preconditioned by the columns' own solves (sub-plan 21 Task 1, design (a))"
```

---

### Task 2: The cells' velocity and the surface operators

**Files:**
- Modify: `reentry_model/skin.py` (`Skins.U`, `Skins.momentum_out`, `Skins.column_momentum`, `Skins.momentum`; every
  cell operation carries the velocity), `reentry_model/columns.py` (`triangle_geometry`, `Columns.gradient`,
  `Columns.curvature`)
- Test: `tests/test_reentry_model_skin.py`, `tests/test_reentry_model_columns.py`

**Interfaces:**
- Consumes: Task 1's `skin.py` and `columns.py`.
- Produces: `Skins.U` (n, K, 3) [m/s]; `Skins.momentum_out` (n, 3) [kg m/s], what the last `take_top` or
  `take_bottom` took away; `Skins.column_momentum() -> (n, 3)`, `Skins.momentum() -> (3,)`;
  `Skins.add_bottom(dm, h, s=1.0, u=None)` (u the arriving metal's velocity, at rest by default);
  `Skins.state()` now six arrays, the velocity last; `columns.triangle_geometry(points, faces) -> (area, unit_normal,
  cot)` with cot (nt, 3) the cotangent at each corner; `Columns.gradient(f, points=None) -> (n, 3)`;
  `Columns.curvature(points=None) -> (n,)` [1/m].

What sub-plan 20 left to this stage (its ruling 7): the velocity the flow will give each cell, and the two surface
operators the push reads. (Its ruling 9's interface rule was settled after that build: damped repetition, Asha,
2026-10-12.) Nothing here moves a cell yet — every velocity stays zero until Task 4 — so this task changes no result: the
100 mm sphere's first six steps give Task 1's history in every column (processor time 221 s against 224 s).

**The velocity** is a vector per cell, stored in three components so that Task 6 can project it between neighbouring
columns' tangent planes. It is part of the cells' state, saved and restored with them, and every cell operation carries
it so that each column's momentum is exact: a split's halves keep their cell's velocity, a merged cell takes the
mass-weighted mean, a partly taken cell keeps its own, and metal drawn from the 3D model arrives at rest (it is rigid)
unless a velocity is given. What leaves at the top or the bottom takes its momentum with it, recorded in
`momentum_out` for the momentum books of Task 8.

**The gradient** of a per-corner field averages the gradients of its linear interpolation on the corner's triangles, by
area, and projects the result onto the corner's tangent plane: exact for a field linear on a flat surface, within
0.09 % of a linear field's tangential part on average on the 4 mm test sphere (0.8 % at its one irregular corner) and
0.02 % on the 2 mm sphere. **The curvature** is the mean curvature, the sum of the principal curvatures, from the
cotangent Laplacian of the position. Sub-plan 20 measured it over the barycentric dual area, right on average but
scattered 7.6 % from corner to corner on the test sphere, which the push's −σ ∇_s κ would turn into flow; over the
mixed Voronoi area (Meyer et al. 2003) the scatter falls to 0.45 % on the 2 mm sphere, and what is left is at the
corners touching an obtuse triangle, where that area is only approximate (33 of the test sphere's 2 254 corners, off by
up to 176 %; 48 of the 2 mm sphere's 9 417, up to 23 %). Every other corner is within 1.1·10⁻⁴ of 2/R, so a corner
touching an obtuse triangle takes the mean of its neighbours that touch none, or of those filled in the pass before:
every corner of both spheres is then within 8.4·10⁻⁵ (test sphere) and 1.3·10⁻⁵ (2 mm) of 2/R (measured 2026-10-18,
`prototype/work-2026-10-18-column-flow/measure/curvature.py`).

- [ ] **Step 1: Append the velocity tests** to `tests/test_reentry_model_skin.py`:

```diff
--- a/tests/test_reentry_model_skin.py
+++ b/tests/test_reentry_model_skin.py
@@ -489,3 +489,58 @@
     sk, col, _ = random_sphere_skins(coarse_sphere_mesh)
     with pytest.raises(RuntimeError, match="coupled sideways solve did not converge"):
         skin._conduct(sk, 0.01, np.zeros(sk.n), sk.area, 0.0, 0.0, np.full(sk.n, 700.0), skin.Sideways(sk, col))
+
+
+# -- the cells' velocity (sub-plan 21 Task 2) ---------------------------------------------------------------------------
+def moving_skins(seed=5):
+    """Five columns filled to two to four times their minimum in fine cells (so the first regrid merges the deep ones),
+    every cell at its own velocity in metres a second."""
+    mat = scheil()
+    s = skin.SkinSettings()
+    sk = skin.Skins(mat, np.full(5, 4e-6), s)
+    rng = np.random.default_rng(seed)
+    sk.fill(skin.minimum_mass(mat, sk.area, s) * rng.uniform(2.0, 4.0, 5), rng.uniform(600.0, 850.0, 5))
+    sk.U = np.where(sk.valid[..., None], rng.normal(0.0, 0.3, sk.m.shape + (3,)), 0.0)
+    return sk, rng
+
+
+def test_a_filled_skin_is_at_rest_and_its_velocity_is_part_of_its_state():
+    mat = scheil()
+    sk = filled(mat)
+    assert np.all(sk.U == 0.0) and sk.U.shape == sk.m.shape + (3,)
+    moving, _ = moving_skins()
+    saved = moving.state()
+    P = moving.column_momentum()
+    moving.U = np.zeros_like(moving.U)
+    moving.set_state(saved)
+    np.testing.assert_array_equal(moving.column_momentum(), P)
+
+
+def test_split_and_merge_keep_every_column_s_momentum():
+    sk, _ = moving_skins()
+    P0, count0 = sk.column_momentum(), sk.count.copy()
+    sk.regrid()                                                       # merges the deep fine cells
+    assert np.all(sk.count < count0)
+    np.testing.assert_allclose(sk.column_momentum(), P0, rtol=1e-12, atol=1e-18)
+    sk.add_bottom(np.full(sk.n, 0.5 * sk.column_mass().min()), np.full(sk.n, float(sk.material.enthalpy(700.0))))
+    P1, count1 = sk.column_momentum(), sk.count.copy()
+    sk.regrid()                                                       # splits the bottom cell the draw thickened
+    assert np.all(sk.count > count1)
+    np.testing.assert_allclose(sk.column_momentum(), P1, rtol=1e-12, atol=1e-18)
+    np.testing.assert_allclose(P1, P0, rtol=1e-12, atol=1e-18)        # the drawn metal arrived at rest
+
+
+def test_what_leaves_at_the_top_or_the_bottom_takes_its_momentum_and_what_arrives_brings_its_own():
+    sk, rng = moving_skins()
+    sk.regrid()
+    P = sk.column_momentum()
+    mass, _ = sk.take_top(np.array([0, 1, 2, 3, 1]))
+    np.testing.assert_allclose(sk.column_momentum() + sk.momentum_out, P, rtol=1e-12, atol=1e-18)
+    assert np.all((np.linalg.norm(sk.momentum_out, axis=1) > 0.0) == (mass > 0.0))
+    P = sk.column_momentum()
+    sk.take_bottom(0.3 * sk.column_mass())                            # whole cells and a part of one
+    np.testing.assert_allclose(sk.column_momentum() + sk.momentum_out, P, rtol=1e-12, atol=1e-18)
+    P = sk.column_momentum()
+    dm, u = np.full(sk.n, 1e-7), rng.normal(0.0, 0.2, (sk.n, 3))
+    sk.add_bottom(dm, np.full(sk.n, float(sk.material.enthalpy(700.0))), u=u)
+    np.testing.assert_allclose(sk.column_momentum(), P + dm[:, None] * u, rtol=1e-12, atol=1e-18)
```

- [ ] **Step 2: Append the surface operators' tests** to `tests/test_reentry_model_columns.py`:

```diff
--- a/tests/test_reentry_model_columns.py
+++ b/tests/test_reentry_model_columns.py
@@ -87,3 +87,46 @@
     v = np.array([-1.0, 0.0, 0.0])
     assert top.angles_to(v)[top.patch_toward(v)] < np.radians(5.0)
     assert top.projected_area(v) == pytest.approx(np.pi * (0.045 + 1e-3) ** 2, rel=0.01)
+
+
+# -- the surface operators (sub-plan 21 Task 2) -------------------------------------------------------------------------
+def obtuse_corners(col, points=None):
+    p = col.points if points is None else points
+    _, _, cot = columns.triangle_geometry(p, col.faces)
+    out = np.zeros(col.n, dtype=bool)
+    for v in range(3):
+        np.logical_or.at(out, col.faces[:, v], (cot < 0.0).any(axis=1))
+    return out
+
+
+def test_the_gradient_is_exact_for_a_field_linear_on_a_flat_surface():
+    pts, faces, inner = flat_grid()
+    col = columns.Columns(pts, faces)
+    g = col.gradient(2.0 + 300.0 * pts[:, 0] - 700.0 * pts[:, 1])
+    assert np.abs(g - np.array([300.0, -700.0, 0.0])).max() < 1e-9          # every corner, the edges' included
+
+
+def test_on_a_sphere_the_gradient_of_a_linear_field_is_its_tangential_part(coarse_sphere_mesh):
+    """On the 4 mm test sphere within 0.09 % of a linear field's tangential part on average and 0.8 % at worst, at the
+    irregular corner gmsh left (0.04 % at worst on the 2 mm sphere; measured 2026-10-18)."""
+    col, _ = sphere_columns(coarse_sphere_mesh)
+    a = np.array([3.0, -1.0, 2.0])
+    exact = a - (col.normal @ a)[:, None] * col.normal
+    err = np.linalg.norm(col.gradient(col.points @ a) - exact, axis=1) / np.linalg.norm(a)
+    assert err.mean() < 0.002 and err.max() < 0.01
+
+
+def test_the_curvature_of_a_sphere_is_two_over_its_radius_at_every_corner(coarse_sphere_mesh):
+    """Every corner within 2e-4 of 2/R, those touching an obtuse triangle included (33 of the 2 254, which the mixed
+    Voronoi area alone left off by up to 176 %; 8.4e-5 measured at worst, 2026-10-18), and the same on a sphere 2 %
+    larger whose corners are given."""
+    col, _ = sphere_columns(coarse_sphere_mesh)
+    assert obtuse_corners(col).sum() > 0                                     # the fallback is exercised
+    R = np.linalg.norm(col.points, axis=1)
+    assert np.abs(col.curvature() * R / 2.0 - 1.0).max() < 2e-4
+    assert np.abs(col.curvature(1.02 * col.points) * 1.02 * R / 2.0 - 1.0).max() < 2e-4
+
+
+def test_the_curvature_of_a_plane_is_zero():
+    pts, faces, inner = flat_grid()
+    assert np.abs(columns.Columns(pts, faces).curvature()[inner]).max() < 1e-9
```

- [ ] **Step 3: Run them to see them fail.** `"$PY" -m pytest tests/test_reentry_model_skin.py
  tests/test_reentry_model_columns.py -q`. Expected: `AttributeError: 'Skins' object has no attribute 'U'` and
  `AttributeError: module 'reentry_model.columns' has no attribute 'triangle_geometry'`.

- [ ] **Step 4: The velocity in the cells** in `reentry_model/skin.py`:

```diff
--- a/reentry_model/skin.py
+++ b/reentry_model/skin.py
@@ -143,7 +143,10 @@
 class Skins:
     """The skins of every column, as padded arrays (module docstring). `area` [m2] is each column's area, the dual
     area of its base corner; the body updates it when the corners move (the cells keep their mass). `S`, when not None,
-    is an intensive per-cell field every cell operation carries (the trial pass's dT/dT_b)."""
+    is an intensive per-cell field every cell operation carries (the trial pass's dT/dT_b). `U` (n, K, 3) is every
+    cell's velocity [m/s] (sub-plan 21 Task 2): part of the cells' state, carried by every operation with its momentum
+    exact -- split halves keep it, merged cells take the mass-weighted mean, metal drawn from the 3D model arrives at
+    rest -- and `momentum_out` (n, 3) [kg m/s] is what the last `take_top` or `take_bottom` took away with the metal."""
 
     def __init__(self, material, area, settings=None):
         self.material = material
@@ -153,6 +156,8 @@
         self.m, self.H, self.T = np.zeros((n, 1)), np.zeros((n, 1)), np.full((n, 1), 300.0)
         self.count = np.zeros(n, dtype=np.int64)
         self.S = None
+        self.U = np.zeros((n, 1, 3))
+        self.momentum_out = np.zeros((n, 3))
         self.rho_l = material.liquid.rho if material.liquid is not None and self.settings.swell else material.rho
         self.props = PropertyTable(material)
 
@@ -174,11 +179,18 @@
     def column_mass(self):
         return self.m.sum(axis=1)
 
+    def column_momentum(self):
+        """Each column's momentum [kg m/s], (n, 3)."""
+        return (self.m[..., None] * self.U).sum(axis=1)
+
+    def momentum(self):
+        return self.column_momentum().sum(axis=0)
+
     def state(self):
-        return tuple(a.copy() for a in (self.area, self.m, self.H, self.T, self.count))
+        return tuple(a.copy() for a in (self.area, self.m, self.H, self.T, self.count, self.U))
 
     def set_state(self, state):
-        self.area, self.m, self.H, self.T, self.count = (a.copy() for a in state)
+        self.area, self.m, self.H, self.T, self.count, self.U = (a.copy() for a in state)
         self.S = None
 
     # -- geometry of the cells --------------------------------------------------------------------------------------
@@ -202,10 +214,10 @@
         return s.cell + (s.growth - 1.0) * np.maximum(z - s.fine_zone, 0.0)
 
     # -- installing cells -------------------------------------------------------------------------------------------
-    def _install(self, m, H, count, S=None, T=None):
+    def _install(self, m, H, count, S=None, T=None, U=None):
         """Install padded cell arrays. Temperatures follow from the enthalpies, except where `T` is given for cells
         whose specific enthalpy the operation did not change (a shifted, split or partly taken cell); nan in T marks
-        the cells to recompute."""
+        the cells to recompute. `U` is the cells' velocity (at rest if not given)."""
         valid = np.arange(m.shape[1])[None, :] < count[:, None]
         self.m, self.H, self.count = np.where(valid, m, 0.0), np.where(valid, H, 0.0), np.asarray(count, dtype=np.int64)
         redo = valid if T is None else valid & np.isnan(T)
@@ -214,6 +226,7 @@
             self.T[redo] = self.material.temperature_from_enthalpy(self.H[redo] / self.m[redo])
         if self.S is not None:
             self.S = np.where(valid, S, 0.0)
+        self.U = np.zeros(self.m.shape + (3,)) if U is None else np.where(valid[..., None], U, 0.0)
 
     def fill(self, mass, T):
         """Every column filled with `mass` [kg] (per column) of metal at temperature T [K] (per column), in equal cells of
@@ -235,6 +248,7 @@
             self.m, self.H, self.T = pad(self.m, 0.0), pad(self.H, 0.0), pad(self.T, 300.0)
             if self.S is not None:
                 self.S = pad(self.S, 0.0)
+            self.U = np.concatenate([self.U, np.zeros((self.n, extra, 3))], axis=1)
 
     # -- split and merge --------------------------------------------------------------------------------------------
     def regrid(self):
@@ -254,7 +268,7 @@
                 self._split(np.flatnonzero(big.any(axis=1)), big)
         raise RuntimeError("the skin cells did not settle in {} split-and-merge passes".format(REGRID_PASSES))
 
-    def _put_rows(self, rows, m, H, T, S, count):
+    def _put_rows(self, rows, m, H, T, S, count, U):
         """Write rebuilt cell arrays for `rows` back (temperatures of nan are recomputed from the enthalpies), padding
         or trimming the arrays to the longest column."""
         self._width(m.shape[1])
@@ -270,11 +284,14 @@
         self.count[rows] = count
         if self.S is not None:
             self.S[rows] = np.where(valid, pad(S), 0.0)
+        U = np.concatenate([U, np.zeros((len(rows), K - U.shape[1], 3))], axis=1) if U.shape[1] < K else U
+        self.U[rows] = np.where(valid[..., None], U, 0.0)
         width = int(self.count.max()) if self.n else 1
         if width < K:
             self.m, self.H, self.T = self.m[:, :width], self.H[:, :width], self.T[:, :width]
             if self.S is not None:
                 self.S = self.S[:, :width]
+            self.U = self.U[:, :width]
 
     def _merge(self, rows, small):
         """Each run of small cells joins the cell beneath that ends it; a run at the bottom with no cell to end it joins
@@ -303,10 +320,13 @@
             mS = np.zeros(n * K)
             np.add.at(mS, flat, (m0 * self.S[rows])[valid])
             S = np.divide(mS, m, out=np.zeros_like(m), where=m > 0.0).reshape(n, K)
-        self._put_rows(rows, m.reshape(n, K), H.reshape(n, K), T.reshape(n, K), S, groups)
+        mU = np.zeros((n * K, 3))
+        np.add.at(mU, flat, (m0[..., None] * self.U[rows])[valid])
+        U = np.divide(mU, m[:, None], out=np.zeros_like(mU), where=m[:, None] > 0.0).reshape(n, K, 3)
+        self._put_rows(rows, m.reshape(n, K), H.reshape(n, K), T.reshape(n, K), S, groups, U)
 
     def _split(self, rows, big):
-        """Every big cell in `rows` becomes two equal halves of its mass and enthalpy (and its T and S)."""
+        """Every big cell in `rows` becomes two equal halves of its mass and enthalpy (and its T, S and velocity)."""
         valid = self.valid[rows]
         reps = np.where(valid, 1 + big[rows], 0)
         n, K = reps.shape
@@ -316,6 +336,7 @@
         r_all = np.repeat(np.arange(n)[:, None], K, axis=1)
         m, H, T = np.zeros((n, K2)), np.zeros((n, K2)), np.zeros((n, K2))
         S = np.zeros((n, K2)) if self.S is not None else None
+        U = np.zeros((n, K2, 3))
         share = 1.0 / reps.clip(1)
         for part in (0, 1):
             sel = valid & (part < reps)
@@ -325,27 +346,32 @@
             T[r, c] = self.T[rows][sel]                                    # halves keep the specific enthalpy
             if S is not None:
                 S[r, c] = self.S[rows][sel]
-        self._put_rows(rows, m, H, T, S, count)
+            U[r, c] = self.U[rows][sel]
+        self._put_rows(rows, m, H, T, S, count, U)
 
     # -- the top and the bottom -------------------------------------------------------------------------------------
     def take_top(self, k):
-        """Remove the top k[i] cells of every column; returns (mass, enthalpy) removed per column."""
+        """Remove the top k[i] cells of every column; returns (mass, enthalpy) removed per column, and leaves the momentum
+        that went with them in `momentum_out`."""
         k = np.asarray(k, dtype=np.int64)
         n, K = self.m.shape
         top = np.arange(K)[None, :] < k[:, None]
         mass, H = (self.m * top).sum(axis=1), (self.H * top).sum(axis=1)
+        self.momentum_out = ((self.m * top)[..., None] * self.U).sum(axis=1)
         if k.any():
             src = np.minimum(np.arange(K)[None, :] + k[:, None], K - 1)
             shift = lambda a: np.take_along_axis(a, src, axis=1)
             S = shift(self.S) if self.S is not None else None
-            self._install(shift(self.m), shift(self.H), self.count - k, S, shift(self.T))
+            U = np.take_along_axis(self.U, src[..., None], axis=1)
+            self._install(shift(self.m), shift(self.H), self.count - k, S, shift(self.T), U)
         return mass, H
 
-    def add_bottom(self, dm, h, s=1.0):
+    def add_bottom(self, dm, h, s=1.0, u=None):
         """Metal dm >= 0 [kg] at specific enthalpy h [J/kg] joins the bottom cell of every column (an emptied column gets
         a new cell); returns the enthalpy added. With S tracked, the metal carries the sensitivity s per column
         (default 1: it is at the base temperature and follows it; 0 for metal whose enthalpy does not depend on the
-        base temperature, such as rigid metal drawn at T_rigid's enthalpy from a base above T_rigid)."""
+        base temperature, such as rigid metal drawn at T_rigid's enthalpy from a base above T_rigid). It arrives with
+        the velocity u (n, 3) [m/s], at rest by default: metal from the 3D model is rigid."""
         dm, h = np.broadcast_to(np.asarray(dm, dtype=float), (self.n,)), np.broadcast_to(np.asarray(h, dtype=float), (self.n,))
         s = np.broadcast_to(np.asarray(s, dtype=float), (self.n,))
         some = dm > 0.0
@@ -359,16 +385,20 @@
         self.H[rows, b] += dm[rows] * h[rows]
         if self.S is not None:
             self.S[rows, b] = (m_old * self.S[rows, b] + dm[rows] * s[rows]) / self.m[rows, b]
+        u = np.zeros((self.n, 3)) if u is None else np.broadcast_to(np.asarray(u, dtype=float), (self.n, 3))
+        self.U[rows, b] = (m_old[:, None] * self.U[rows, b] + dm[rows, None] * u[rows]) / self.m[rows, b][:, None]
         self.count[rows] = np.maximum(self.count[rows], 1)
         self.T[rows, b] = self.material.temperature_from_enthalpy(self.H[rows, b] / self.m[rows, b])
         return np.where(some, dm * h, 0.0)
 
     def take_bottom(self, dm):
-        """Remove dm >= 0 [kg] per column from the bottom up, each cell's metal at its own specific enthalpy (a partly
-        taken cell keeps its specific enthalpy); returns the enthalpy removed per column."""
+        """Remove dm >= 0 [kg] per column from the bottom up, each cell's metal at its own specific enthalpy and velocity
+        (a partly taken cell keeps both); returns the enthalpy removed per column, and leaves the momentum that went in
+        `momentum_out`."""
         dm = np.broadcast_to(np.asarray(dm, dtype=float), (self.n,))
         if np.any(dm > self.column_mass() * (1.0 + 1e-12) + 1e-300):
             raise ValueError("cannot take more metal from a skin's bottom than the skin holds")
+        self.momentum_out = np.zeros((self.n, 3))
         if not np.any(dm > 0.0):
             return np.zeros(self.n)
         n, K = self.m.shape
@@ -377,10 +407,11 @@
         taken = np.clip(dm[:, None] - below, 0.0, self.m)                 # how much of each cell goes
         frac = np.divide(taken, self.m, out=np.zeros_like(taken), where=self.m > 0.0)
         H_taken = (self.H * frac).sum(axis=1)
+        self.momentum_out = (taken[..., None] * self.U).sum(axis=1)
         m, H = self.m - taken, self.H * (1.0 - frac)
         keep = valid & (m > 0.0)
         count = keep.sum(axis=1)                                          # the kept cells are the top ones
-        self._install(m, H, count, self.S, self.T)                        # a partly taken cell keeps its T
+        self._install(m, H, count, self.S, self.T, self.U)                # a partly taken cell keeps its T and velocity
         return H_taken
 
     def exchange_bottom(self, dm, dH):
```

- [ ] **Step 5: The surface operators** in `reentry_model/columns.py`:

```diff
--- a/reentry_model/columns.py
+++ b/reentry_model/columns.py
@@ -13,7 +13,17 @@
 to make it negative: a negative conductance would break the explicit step's positivity (and, since sub-plan 21 makes the
 sideways conduction implicit, the definiteness its conjugate gradients need), and it is rare (3 of 6 756 edges on the
 4 mm test sphere, 4 of 28 245 on the default 2 mm sphere, measured 2026-10-12).
-The surface gradient and curvature operators the flow needs are added in sub-plan 21, the first stage that uses them."""
+
+The flow (sub-plan 21) reads two surface operators on the same graph. `gradient` averages the gradients of a field's
+linear interpolation over each corner's triangles, by area, and projects the result onto the corner's tangent plane.
+`curvature` is the mean curvature (the sum of the principal curvatures) from the cotangent Laplacian of the position
+over each corner's mixed Voronoi area (Meyer et al. 2003). That area is exact only where no touching triangle is
+obtuse, so a corner touching one takes the mean of its neighbours that touch none (33 of the 4 mm test sphere's 2 254
+corners and 48 of the 2 mm sphere's 9 417, which the mixed area alone left off by 15 and 3 % on average and by up to 176
+and 23 %). On these 50 mm-radius spheres the curvature is then within 8.4·10^-5 of 2/R at every corner (1.3·10^-5 on
+the 2 mm sphere), where the barycentric area had scattered 7.6 and 3.9 %, and the gradient of a linear field within
+0.09 % of its tangential part on average (0.02 % on the 2 mm sphere), 0.8 % at the test sphere's one irregular corner
+(measured 2026-10-18)."""
 from dataclasses import dataclass
 
 import numpy as np
@@ -69,6 +79,62 @@
         self.conduct_weight = np.maximum(self.cot_weight, 0.0)
         self.face_length = face
 
+    def gradient(self, f, points=None):
+        """The tangential gradient at every corner of a per-corner field f (n,) [f per m] on the surface through
+        `points` (default: the columns' own corners): the gradients of f's linear interpolation on the touching
+        triangles, averaged by area, projected onto the corner's tangent plane. Exact for a field linear on a flat
+        surface."""
+        p = self.points if points is None else np.asarray(points, dtype=float)
+        normal = self.normal if points is None else corner_frame(p, self.faces, np.arange(self.n))[0]
+        A, unit, _ = triangle_geometry(p, self.faces)
+        a, b, c = (p[self.faces[:, v]] for v in range(3))
+        fa, fb, fc = (np.asarray(f, dtype=float)[self.faces[:, v]] for v in range(3))
+        g = (fa[:, None] * np.cross(unit, c - b) + fb[:, None] * np.cross(unit, a - c)
+             + fc[:, None] * np.cross(unit, b - a)) / (2.0 * A[:, None])
+        G, W = np.zeros((self.n, 3)), np.zeros(self.n)
+        for v in range(3):
+            np.add.at(G, self.faces[:, v], g * A[:, None])
+            np.add.at(W, self.faces[:, v], A)
+        G /= W[:, None]
+        return G - np.einsum("ij,ij->i", G, normal)[:, None] * normal
+
+    def curvature(self, points=None):
+        """The mean curvature, the sum of the principal curvatures [1/m], at every corner of the surface through
+        `points` (default: the columns' own corners), positive where the surface bends away from its outward normal (2/R
+        on a sphere): the cotangent Laplacian of the position along the corner's normal over its mixed Voronoi area. A
+        corner touching an obtuse triangle takes the mean of its neighbours that touch none, or, where every neighbour
+        touches one, of those that took a value in the pass before (module docstring)."""
+        p = self.points if points is None else np.asarray(points, dtype=float)
+        normal = self.normal if points is None else corner_frame(p, self.faces, np.arange(self.n))[0]
+        A, _, cot = triangle_geometry(p, self.faces)
+        lap, area = np.zeros((self.n, 3)), np.zeros(self.n)
+        obtuse = np.zeros(self.n, dtype=bool)
+        blunt = (cot < 0.0).any(axis=1)
+        for v in range(3):
+            i, j, k = self.faces[:, v], self.faces[:, (v + 1) % 3], self.faces[:, (v + 2) % 3]
+            cj, ck = cot[:, (v + 1) % 3], cot[:, (v + 2) % 3]
+            dj, dk = p[j] - p[i], p[k] - p[i]
+            np.add.at(lap, i, 0.5 * (ck[:, None] * dj + cj[:, None] * dk))
+            voronoi = 0.125 * (np.einsum("ij,ij->i", dj, dj) * ck + np.einsum("ij,ij->i", dk, dk) * cj)
+            np.add.at(area, i, np.where(cot[:, v] < 0.0, 0.5 * A, np.where(blunt, 0.25 * A, voronoi)))
+            np.logical_or.at(obtuse, i, blunt)
+        kappa = -np.einsum("ij,ij->i", lap, normal) / area
+        i, j = self.edges[:, 0], self.edges[:, 1]
+        known = ~obtuse
+        while not known.all():                     # outward from the corners whose own value is exact, a ring a pass
+            good = known.astype(float)
+            total, count = np.zeros(self.n), np.zeros(self.n)
+            np.add.at(total, i, kappa[j] * good[j])
+            np.add.at(total, j, kappa[i] * good[i])
+            np.add.at(count, i, good[j])
+            np.add.at(count, j, good[i])
+            fill = ~known & (count > 0.0)
+            if not fill.any():
+                break
+            kappa[fill] = total[fill] / count[fill]
+            known = known | fill
+        return kappa
+
     def top_points(self, thickness):
         """The top surface: each corner moved out along its normal by its column's thickness [m]."""
         return self.points + np.asarray(thickness, dtype=float)[:, None] * self.normal
@@ -96,6 +162,17 @@
         return self._incidence
 
 
+def triangle_geometry(points, faces):
+    """Per triangle: its area [m2], its unit normal (outward for outward-oriented faces) and the cotangent of the angle
+    at each of its three corners, (nt, 3)."""
+    a, b, c = (points[faces[:, v]] for v in range(3))
+    n2 = np.cross(b - a, c - a)
+    twice = np.linalg.norm(n2, axis=1)
+    legs = ((b - a, c - a), (c - b, a - b), (a - c, b - c))
+    cot = np.stack([np.einsum("ij,ij->i", u, w) / np.linalg.norm(np.cross(u, w), axis=1) for u, w in legs], axis=1)
+    return 0.5 * twice, n2 / twice[:, None], cot
+
+
 @dataclass
 class CornerSurface:
     """The gas-facing surface of a corner-column body as the coupled run reads it: one 'patch' per corner (its top dual
```

- [ ] **Step 6: Run them to see them pass.** Expected: `tests/test_reentry_model_skin.py` 29 passed and
  `tests/test_reentry_model_columns.py` 9 passed; the unit tier 416 passed, 1 skipped. The first six steps of the
  100 mm level-3 flight (`--heating sesam --surface-model skin --t-max 3`) give Task 1's history in every column.

- [ ] **Step 7: Commit.**

```bash
git add reentry_model/skin.py reentry_model/columns.py tests/test_reentry_model_skin.py tests/test_reentry_model_columns.py
git commit -m "Skins: every cell carries a velocity with its momentum exact through split, merge, take and add, and the corner graph gains the surface gradient and the mean curvature over mixed Voronoi areas (sub-plan 21 Task 2)"
```

---

### Task 3: Viscosity per interface

**Files:**
- Create: `reentry_model/rheology.py`, `tests/test_reentry_model_rheology.py`
- Modify: `reentry_model/skin.py` (`cell_viscosity`, `interface_drag`), `tests/test_reentry_model_skin.py`

**Interfaces:**
- Consumes: Task 2's `skin.py`; `material.Material.melts`, `.liquid.mu` and `.liquid_fraction`; `Skins.props`.
- Produces: `rheology.li_viscosity(f_s, gamma)` and `rheology.viscosity(T, f_s, gamma, mu_liquid, f_s_bridge)` [Pa s],
  with `rheology.GAMMA_MAX` (367 1/s), `BRIDGE_LOW` and `BRIDGE_HIGH` (906.4 and 908 K) and `RIGID_SOLID_FRACTION`
  (0.5); `skin.cell_viscosity(sk, gamma) -> (n, K)` [Pa s], infinite in rigid and missing cells;
  `skin.interface_drag(sk, gamma) -> (n, K-1)` [Pa s/m], gamma the interfaces' shear rates [1/s].

The viscosity law is the Step 4 plan's (§4), as spec §2 records it: liquid aluminium's 1.3 mPa·s at and above the
liquidus; a log-linear bridge over 906.4–908 K from Li et al.'s value where Scheil's curve stands at 10 % solid (9.75 %
measured); Li Yageng et al. (2014), η = [0.871 − 0.00849 γ̇^0.74924] exp(3.7311 f_s) Pa·s, from there down to half
solid, the shear rate held at no more than 367 1/s (its bracket would reach zero at 483 1/s); and infinite at and above
half solid, where a cell is rigid and anchored. At rest that is 0.87 Pa·s at f_s = 0 and 5.6 Pa·s at f_s = 0.5, 670 to
4 330 times the liquid. It lives in `rheology.py`, which imports neither architecture's modules, so that the later move
into `reentry_model/shared/` stays mechanical (the 2026-10-12 amendment); the Step 4 plan's Chen branch above half
solid is not needed where rigid cells are anchored, and is left to that plan.

**Ruling (2026-10-19):** the drag between neighbouring cells puts their two half-cells in series,
D = 1 / (dz_i / 2η_i + dz_j / 2η_j) per unit area, each cell's viscosity taken at the interface's shear rate — the
spec's "viscosity at their interface" (§6.4) in the form the stress's continuity across the interface requires, the way
`_conduct` combines conductivities. Across a jump from liquid to slurry the liquid's half dominates and the profile
Task 4 builds on it has the exact kink; a rigid cell resists without limit and adds nothing to the sum, so a flowing
cell drags on the rigid cell beneath over its own half-thickness alone — no slip at the rigid front. The fixed-point
iteration on the shear rate (at most five passes) is Task 4's: the rate comes from the velocities it solves for.
(Corrected in Task 4, 2026-10-20: each cell's viscosity is taken at its own shear rate, not the interface's, which
would give a slurry half-cell under liquid the liquid's rate.)

- [ ] **Step 1: Write the rheology tests**, `tests/test_reentry_model_rheology.py`:

```python
"""The viscosity of liquid and slurry AA7075 (sub-plan 21 Task 3; the Step 4 plan's section 4): Li et al.'s slurry law,
its shear rate held at 367 1/s, the bridge to the liquid at the liquidus, and rigid metal at half solid."""
import numpy as np
import pytest

from reentry_model import material, rheology

MU_L = 1.3e-3


def scheil():
    return material.Material.from_drama_json(material.MATERIAL_NAMES["AA7075_scheil"])


def test_at_rest_li_et_al_runs_from_0_87_to_5_6_pa_s_670_to_4300_times_the_liquid():
    assert rheology.li_viscosity(0.0, 0.0) == pytest.approx(0.871)
    assert rheology.li_viscosity(0.5, 0.0) == pytest.approx(0.871 * np.exp(0.5 * 3.7311))         # 5.63 Pa s
    assert rheology.li_viscosity(0.0, 0.0) / MU_L == pytest.approx(670.0)
    assert rheology.li_viscosity(0.5, 0.0) / MU_L == pytest.approx(4328.0, rel=1e-3)
    # the published form at a measured point: f_s = 0.1, 61 1/s
    assert rheology.li_viscosity(0.1, 61.0) == pytest.approx((0.871 - 0.00849 * 61.0 ** 0.74924) * np.exp(0.37311))


def test_the_shear_rate_is_held_at_367_per_second_where_the_fit_still_holds():
    f_s = np.linspace(0.1, 0.5, 5)
    np.testing.assert_array_equal(rheology.li_viscosity(f_s, 450.0), rheology.li_viscosity(f_s, 367.0))
    assert np.all(rheology.li_viscosity(f_s, 367.0) > 0.0)                    # the bracket's zero is at 483 1/s
    rates = np.array([0.0, 61.0, 150.0, 367.0])
    assert np.all(np.diff(rheology.li_viscosity(0.3, rates)) < 0.0)           # shear-thinning


def test_the_bridge_runs_log_linearly_from_the_slurry_to_the_liquid_at_the_liquidus():
    f_b = 0.0975
    top = rheology.li_viscosity(f_b, 100.0)
    eta = lambda T, f_s: float(rheology.viscosity(T, f_s, 100.0, MU_L, f_b))
    assert eta(906.4, f_b) == pytest.approx(top)
    assert eta(907.2, 0.05) == pytest.approx(np.sqrt(top * MU_L))             # half way: the geometric mean
    assert eta(908.0, 0.0) == MU_L and eta(950.0, 0.0) == MU_L


def test_metal_half_solid_or_more_is_rigid():
    assert np.isinf(rheology.viscosity(895.0, 0.5, 50.0, MU_L, 0.0975))
    assert np.isinf(rheology.viscosity(850.0, 0.9, 50.0, MU_L, 0.0975))
    assert rheology.viscosity(895.2, 0.4999, 50.0, MU_L, 0.0975) == pytest.approx(rheology.li_viscosity(0.4999, 50.0))


def test_on_aa7075_s_scheil_curve_the_viscosity_falls_without_a_break_from_the_rigid_front_to_the_liquid():
    mat = scheil()
    f_b = 1.0 - float(mat.liquid_fraction(np.array([rheology.BRIDGE_LOW]))[0])
    assert f_b == pytest.approx(0.0975, abs=1e-4)                              # Scheil's curve at 906.4 K
    T = np.linspace(mat.T_rigid + 1e-3, 910.0, 4001)
    eta = rheology.viscosity(T, 1.0 - mat.liquid_fraction(T), 100.0, mat.liquid.mu, f_b)
    assert np.all(np.isfinite(eta)) and np.all(np.diff(eta) <= 1e-12 * eta[:-1])
    for T_edge in (rheology.BRIDGE_LOW, rheology.BRIDGE_HIGH):                 # no jump where the bridge begins or ends
        Ts = np.array([T_edge - 1e-9, T_edge + 1e-9])
        lo, hi = rheology.viscosity(Ts, 1.0 - mat.liquid_fraction(Ts), 100.0, mat.liquid.mu, f_b)
        assert hi == pytest.approx(lo, rel=1e-6)
```

- [ ] **Step 2: Append the drag tests** to `tests/test_reentry_model_skin.py`:

```diff
--- a/tests/test_reentry_model_skin.py
+++ b/tests/test_reentry_model_skin.py
@@ -544,3 +544,46 @@
     dm, u = np.full(sk.n, 1e-7), rng.normal(0.0, 0.2, (sk.n, 3))
     sk.add_bottom(dm, np.full(sk.n, float(sk.material.enthalpy(700.0))), u=u)
     np.testing.assert_allclose(sk.column_momentum(), P + dm[:, None] * u, rtol=1e-12, atol=1e-18)
+
+
+# -- viscosity (sub-plan 21 Task 3) -------------------------------------------------------------------------------------
+def column_at(T_cells, mat=None):
+    """One column of 10 um cells at the given temperatures, top first."""
+    mat = mat or scheil()
+    T_cells = np.asarray(T_cells, dtype=float)
+    sk = skin.Skins(mat, np.array([4e-6]), skin.SkinSettings(**CELLS_04))
+    m = np.full((1, len(T_cells)), mat.rho * 4e-6 * 1e-5)
+    sk._install(m, m * mat.enthalpy(T_cells)[None, :], np.array([len(T_cells)]))
+    return sk
+
+
+def test_two_liquid_cells_drag_on_each_other_with_the_liquid_s_viscosity_over_their_spacing():
+    sk = column_at([950.0, 940.0, 930.0])
+    D = skin.interface_drag(sk, 200.0)
+    dz = sk.dz()[0]
+    np.testing.assert_allclose(D[0], 1.3e-3 / (0.5 * (dz[:-1] + dz[1:])), rtol=1e-12)
+
+
+def test_across_a_viscosity_jump_the_two_half_cells_act_in_series():
+    """Liquid over slurry: the drag is 1 / (dz_l / 2 mu_l + dz_s / 2 eta_s), which the liquid's half dominates -- the
+    stress is continuous through the jump, and the profile Task 4 builds on it has the exact kink."""
+    sk = column_at([920.0, 900.0])
+    eta = skin.cell_viscosity(sk, 150.0)[0]
+    dz = sk.dz()[0]
+    assert eta[0] == 1.3e-3 and eta[1] > 500 * eta[0]
+    assert skin.interface_drag(sk, 150.0)[0, 0] == pytest.approx(1.0 / (0.5 * dz[0] / eta[0] + 0.5 * dz[1] / eta[1]))
+
+
+def test_a_flowing_cell_drags_on_the_rigid_cell_beneath_over_its_own_half_thickness_alone():
+    sk = column_at([905.0, 890.0, 880.0])                         # slurry, then two rigid cells (T_rigid 895.1 K)
+    eta = skin.cell_viscosity(sk, 50.0)[0]
+    assert np.isfinite(eta[0]) and np.all(np.isinf(eta[1:]))
+    D = skin.interface_drag(sk, 50.0)[0]
+    assert D[0] == pytest.approx(eta[0] / (0.5 * sk.dz()[0, 0]))  # no slip at the rigid front
+    assert D[1] == 0.0                                            # rigid on rigid: nothing moves
+
+
+def test_metal_that_does_not_melt_never_flows():
+    mat = material_mod.Material.from_drama_json("AA7075_nomelt")
+    sk = column_at([1500.0, 1400.0], mat)
+    assert np.all(np.isinf(skin.cell_viscosity(sk, 10.0))) and np.all(skin.interface_drag(sk, 10.0) == 0.0)
```

- [ ] **Step 3: Run them to see them fail.** `"$PY" -m pytest tests/test_reentry_model_rheology.py
  tests/test_reentry_model_skin.py -q`. Expected: `ModuleNotFoundError: No module named 'reentry_model.rheology'`.

- [ ] **Step 4: Create `reentry_model/rheology.py`:**

```python
"""The viscosity of liquid and slurry AA7075 (the Step 4 plan's section 4, docs/superpowers/plans/2026-09-27-deformation-
ring-shedding.md; the layered-skin design section 6.4). Shared physics: it imports neither architecture's modules, so
that the later move into reentry_model/shared/ stays mechanical (sub-plan 21's amendment of 2026-10-12).

- Liquid, at and above the liquidus: the material's own viscosity (1.3 mPa s for AA7075).
- Bridge, BRIDGE_LOW to BRIDGE_HIGH (906.4-908 K): log-linear in the temperature from Li et al.'s value at BRIDGE_LOW,
  where Scheil's curve stands at 10 % solid, to the liquid at the liquidus (labelled: no data).
- Slurry, from the bridge down to half solid: Li Yageng et al. (2014, China Foundry 11, 79-84),
  eta = [0.871 - 0.00849 gamma^0.74924] exp(3.7311 f_s) Pa s, measured at f_s = 0.1-0.5 and shear rates gamma of
  61-490 1/s on continuously sheared slurry; gamma is held at no more than GAMMA_MAX = 367 1/s, the highest rate at
  which the fit still matches the data (its bracket would reach zero at 483 1/s). Between zero and 61 1/s it is
  extrapolated, and there it is the runny side of the data (1.4-1.9 times low at 61 1/s; the Step 4 plan's fact 10).
  From 0.87 Pa s (f_s = 0) to 5.6 Pa s (f_s = 0.5) at rest, 670 to 4 300 times the liquid.
- Rigid, at or above half solid: infinite. A rigid cell does not flow; it is anchored (spec 6.4)."""
import numpy as np

LI_A = 0.871                       # Pa s
LI_B = 0.00849                     # Pa s^(1 + LI_N)
LI_N = 0.74924
LI_C = 3.7311
GAMMA_MAX = 367.0                  # 1/s: Li et al.'s law is held at this shear rate
BRIDGE_LOW, BRIDGE_HIGH = 906.4, 908.0     # K: the log-linear bridge from the slurry to the liquid
RIGID_SOLID_FRACTION = 0.5         # at or above it a cell is rigid (material.RIGID_LIQUID_FRACTION's complement)


def li_viscosity(f_s, gamma):
    """Li et al.'s slurry viscosity [Pa s] at solid fraction f_s and shear rate gamma [1/s], the rate held within
    [0, GAMMA_MAX]."""
    g = np.clip(np.asarray(gamma, dtype=float), 0.0, GAMMA_MAX)
    return (LI_A - LI_B * g ** LI_N) * np.exp(LI_C * np.asarray(f_s, dtype=float))


def viscosity(T, f_s, gamma, mu_liquid, f_s_bridge):
    """The viscosity [Pa s] of metal at temperature T [K] and solid fraction f_s, sheared at gamma [1/s] (module
    docstring); `mu_liquid` [Pa s] is the liquid's and `f_s_bridge` the material's solid fraction at BRIDGE_LOW (0.098 on
    AA7075's Scheil curve). Infinite where the metal is rigid."""
    T, f_s = np.asarray(T, dtype=float), np.asarray(f_s, dtype=float)
    top = li_viscosity(f_s_bridge, gamma)
    w = np.clip((T - BRIDGE_LOW) / (BRIDGE_HIGH - BRIDGE_LOW), 0.0, 1.0)
    bridge = np.exp((1.0 - w) * np.log(top) + w * np.log(mu_liquid))
    eta = np.where(T >= BRIDGE_HIGH, mu_liquid, np.where(T > BRIDGE_LOW, bridge, li_viscosity(f_s, gamma)))
    return np.where(f_s >= RIGID_SOLID_FRACTION, np.inf, eta)
```

- [ ] **Step 5: The cells' viscosity and the drag between them** in `reentry_model/skin.py`:

```diff
--- a/reentry_model/skin.py
+++ b/reentry_model/skin.py
@@ -32,6 +32,7 @@
 
 import numpy as np
 
+from . import rheology
 from .thermal import SIGMA_SB
 
 SKIN_PRESETS = {"production": 10.0e-6, "dev": 20.0e-6}   # m, the fine cell of each preset (spec section 6.2)
@@ -778,6 +779,39 @@
         return np.where(self.valid, np.diff(_interp_apply(self.geometry["to_cells"], Qcum), axis=1), 0.0)
 
 
+# -- viscosity (sub-plan 21 Task 3) ------------------------------------------------------------------------------------
+def cell_viscosity(sk, gamma):
+    """Every cell's viscosity [Pa s], (n, K), sheared at gamma [1/s] (per cell, or one for all): `rheology.viscosity` at
+    the cell's temperature and Scheil solid fraction; infinite in a rigid cell (half solid or more) and beyond each
+    column's cells, and everywhere for a material that does not melt."""
+    mat = sk.material
+    if not mat.melts or mat.liquid is None:
+        return np.full(sk.m.shape, np.inf)
+    f_s = 1.0 - sk.props.liquid_fraction(sk.T)
+    f_s_bridge = 1.0 - float(mat.liquid_fraction(np.array([rheology.BRIDGE_LOW]))[0])
+    eta = rheology.viscosity(sk.T, f_s, gamma, mat.liquid.mu, f_s_bridge)
+    return np.where(sk.valid, eta, np.inf)
+
+
+def interface_drag(sk, gamma):
+    """Per column, the viscous drag per unit area between neighbouring cells' centres [Pa s/m], (n, K-1):
+    1 / (dz_i / (2 eta_i) + dz_j / (2 eta_j)), both cells' viscosities taken at the interface's shear rate gamma [1/s]
+    (n, K-1) -- the two half-cells in series, as the shear stress's continuity across the interface requires, the way
+    `_conduct` combines conductivities. A rigid cell resists without limit and adds nothing to the sum, so a flowing cell
+    drags on the rigid cell beneath it over its own half-thickness alone: no slip at the rigid front. Zero between two
+    rigid cells and across a missing one."""
+    if sk.m.shape[1] < 2:
+        return np.zeros((sk.n, 0))
+    g = np.asarray(gamma, dtype=float) * np.ones((sk.n, sk.m.shape[1] - 1))
+    eta_up = cell_viscosity(sk, np.pad(g, ((0, 0), (0, 1)), mode="edge"))[:, :-1]
+    eta_dn = cell_viscosity(sk, np.pad(g, ((0, 0), (1, 0)), mode="edge"))[:, 1:]
+    dz = sk.dz()
+    resist = 0.5 * dz[:, :-1] / eta_up + 0.5 * dz[:, 1:] / eta_dn
+    link = sk.valid[:, 1:]
+    flows = np.isfinite(eta_up) | np.isfinite(eta_dn)
+    return np.where(link & flows, 1.0 / np.where(link & flows, resist, 1.0), 0.0)
+
+
 # -- the base (Task 7) -------------------------------------------------------------------------------------------------
 def hand_down(sk, m_min, T_rigid):
     """Rigid cells at the bottom of each column (at or below T_rigid, contiguous from the bottom) pass to the 3D model,
```

- [ ] **Step 6: Run them to see them pass.** Expected: `tests/test_reentry_model_rheology.py` 5 passed and
  `tests/test_reentry_model_skin.py` 33 passed; the unit tier 425 passed, 1 skipped. No result changes: nothing reads
  the viscosity until Task 4.

- [ ] **Step 7: Commit.**

```bash
git add reentry_model/rheology.py tests/test_reentry_model_rheology.py reentry_model/skin.py tests/test_reentry_model_skin.py
git commit -m "Skins: the viscosity of every cell from Li et al.'s slurry law with its bridge to the liquid, and the drag between neighbouring cells with their half-cells in series (sub-plan 21 Task 3)"
```

---

### Task 4: Momentum through the depth

**Files:**
- Modify: `reentry_model/skin.py` (`series_drag`, and `interface_drag` at each cell's own shear rate; `SHEAR_PASSES`,
  `SHEAR_TOL`, `CRUST_VISCOSITY`, `MomentumBook`, `anchored_cells`, `momentum_step`), `tests/test_reentry_model_skin.py`

**Interfaces:**
- Consumes: Task 2's `Skins.U` and `Skins.column_momentum`; Task 3's `cell_viscosity`; `thomas`.
- Produces: `skin.series_drag(sk, eta) -> (n, K-1)` [Pa s/m] from the cells' viscosities; `skin.interface_drag(sk,
  gamma)` now with gamma each cell's own shear rate; `skin.momentum_step(sk, dt, push, shear, A_top) -> MomentumBook`
  with push (n, 3) [N/m3] the same on every layer of a column, shear (n, 3) [Pa] on the top cell and A_top (n,) [m2];
  `MomentumBook(push, shear, wall, passes)`, the impulses [N s] (n, 3); `skin.anchored_cells(sk, flowing) -> (n, K)`
  bool, the rigid cells joined to the substrate; `skin.SHEAR_PASSES` (5), `skin.SHEAR_TOL` (10⁻³) and
  `skin.CRUST_VISCOSITY` (10⁹ Pa·s).

Every column's momentum through its depth in one implicit sub-step (spec §6.4), each velocity component apart: the
cells' mass over the sub-step, the drag between neighbours from Task 3, the push on every flowing cell's volume and the
gas shear on the top cell's top area — a tridiagonal system per column and component, solved by `thomas` as the heat
is. The rigid metal joined to the substrate — each column's run of rigid cells up from its bottom cell, which rests on
the rigid 3D model — is held at rest (its momentum, if it had any, is booked to the wall); a moving cell drags on it,
and a flowing bottom cell on the rigid 3D metal beneath the base, with no slip over its own half-thickness. **A crust
floats** (Asha, 2026-10-20, option B, revising the spec's assumption that a crust is anchored): rigid cells with flowing
metal between them and the substrate move as one piece, given a viscosity of 10⁹ Pa·s, dragged by the gas shear on
their top, pushed like the rest and dragged by the slurry beneath. The columns carry no sideways stress in a solid, so
a crust over many columns moves as column-sized plates — as a warm crust, barely solid near half liquid, would — and a
crust held at its rim is not held; anchoring a crust only where it touches rigid metal at the surface (option C) is
kept as a refinement should Task 8 find crusts over slurry common.
The shear rates come from a fixed point started at the sub-step's velocities: each pass solves the system with the
cells' viscosities, then takes every flowing cell's shear rate as the mean over its two faces of the stress there over
its viscosity, and stops after five passes or when no viscosity moves by more than 10⁻³. **Each cell at its own rate**
(found while writing this task): Task 3 took both half-cells at an interface at that interface's rate — the velocity
difference over the distance between the centres — which at a jump from liquid to slurry mixes the liquid's
770 1/s with the slurry's 0.3 1/s and gives the slurry's half-cell five times too low a viscosity; at its own rate the
two-layer profile is exact. The impulses (push, shear, the wall's drag and held momentum) add up to the momentum
change to round-off, for Task 8's momentum books. Nothing calls `momentum_step` yet: moving the metal is Task 6's, so
no result changes.

Measured against exact solutions (2026-10-19/20): shear alone, Couette's line to 1.1·10⁻¹⁰ in one pass; the push
alone, the half-parabola to 2.8·10⁻⁴ of its top speed at 30 cells, 6.9·10⁻⁵ at 60 and 1.7·10⁻⁵ at 120 (second
order); 1 Pa on liquid over slurry at 900 K (37 % solid), the slurry at 3.46 Pa·s shearing at 0.29 1/s, both
profiles to 1.1·10⁻⁵ in two passes; Stokes' first problem (1 Pa suddenly on 0.9 mm of liquid, 18.5 ms, a diffusion
length of 100 µm), the top cell to 1.3·10⁻³, 6.6·10⁻⁴ and 3.3·10⁻⁴ at 100, 200 and 400 steps (first order, backward
Euler), the momentum the shear's impulse exactly; a crust of three rigid cells on 15 of slurry at 900 K under 1 Pa of
gas shear rides the top of the slurry's Couette line as one piece, to 10⁻⁴.

The layer's response develops over time, not at once (Asha's question, 2026-10-20): the push acts on every layer from
the moment it is applied — the pressure crosses a millimetre of metal in about 10⁻⁷ s and the deceleration acts on all
of it — but the velocity spreads from the anchored front, and the gas shear's from the top, over the viscous time
h²/ν: 0.018 s for 0.1 mm of liquid, 0.46 s for 0.5 mm, 1.85 s for 1 mm, and 0.2 ms for 0.5 mm of slurry at 3.5 Pa·s,
whose high viscosity carries momentum fastest. The momentum sub-step keeps the cells' inertia, so it follows that
development. Started from rest and compared over 1 s with 0.1 ms sub-steps on the same cells, the 10 ms sub-steps were
off by 15 % of the speed at 1 s in the first sub-step for 0.1 mm of liquid under 10 kPa/m, which settles within about
two sub-steps, and by 0.04 % after 0.1 s; by at most 1.0 % for 0.5 mm of liquid under the same push and 1.9 % under
10 Pa of gas shear; and 1.1 % in the first sub-step for 0.5 mm of slurry at 900 K. In flight the forcing changes over
seconds, so only the thinnest liquid layers' first sub-step after a sudden onset lags.

- [ ] **Step 1: Append the momentum tests** to `tests/test_reentry_model_skin.py`:

```diff
--- a/tests/test_reentry_model_skin.py
+++ b/tests/test_reentry_model_skin.py
@@ -4,7 +4,7 @@
 import pytest
 
 from reentry_model import material as material_mod
-from reentry_model import skin
+from reentry_model import rheology, skin
 
 
 def scheil():
@@ -587,3 +587,145 @@
     mat = material_mod.Material.from_drama_json("AA7075_nomelt")
     sk = column_at([1500.0, 1400.0], mat)
     assert np.all(np.isinf(skin.cell_viscosity(sk, 10.0))) and np.all(skin.interface_drag(sk, 10.0) == 0.0)
+
+
+# -- momentum through the depth (sub-plan 21 Task 4) --------------------------------------------------------------------
+MU_L = 1.3e-3
+STEADY = 1e9                       # s: a backward-Euler step this long gives the steady profile
+
+
+def heights(sk, k):
+    """The first k cells' centres' height [m] above the rigid front beneath them."""
+    dz = sk.dz()[0][:k]
+    return dz.sum() - (np.cumsum(dz) - 0.5 * dz)
+
+
+def books_close(sk, P0, book):
+    """The impulses add up to the momentum change, to round-off of the largest of them (in a steady step they are
+    large and nearly cancel)."""
+    parts = (book.push, book.shear, book.wall, P0, sk.column_momentum())
+    scale = max(np.abs(x).max() for x in parts) or 1e-30
+    return np.abs(sk.column_momentum() - P0 - book.push - book.shear - book.wall).max() <= 1e-12 * scale
+
+
+def test_shear_alone_gives_a_linear_profile():
+    """Liquid sheared at its top by 1 Pa over rigid metal: u = tau s / mu at every cell centre, s the height above the
+    rigid front (Couette), exactly; the drag to the rigid front over the bottom liquid cell's half-thickness is no slip."""
+    sk = column_at([950.0] * 30 + [880.0] * 10)
+    P0 = sk.column_momentum()
+    book = skin.momentum_step(sk, STEADY, np.zeros(3), np.array([1.0, 0.0, 0.0]), sk.area)
+    np.testing.assert_allclose(sk.U[0, :30, 0], heights(sk, 30) / MU_L, rtol=1e-9)
+    assert np.all(sk.U[0, :, 1:] == 0.0) and np.all(sk.U[0, 30:] == 0.0) and book.passes == 1
+    assert books_close(sk, P0, book)
+
+
+def test_a_push_alone_gives_a_half_parabola_converging_at_second_order():
+    """A push of 10 kPa/m on a free liquid layer over rigid metal: u = (G / mu)(h s - s^2 / 2). Measured 2026-10-19:
+    2.8e-4 of the top's speed off at 30 cells, 6.9e-5 at 60."""
+    def error(N):
+        sk = column_at([950.0] * N + [880.0] * 10)
+        skin.momentum_step(sk, STEADY, np.array([1e4, 0.0, 0.0]), np.zeros(3), sk.area)
+        s, h = heights(sk, N), sk.dz()[0][:N].sum()
+        exact = 1e4 / MU_L * (h * s - 0.5 * s * s)
+        return np.abs(sk.U[0, :N, 0] - exact).max() / exact.max()
+    e30, e60 = error(30), error(60)
+    assert e30 < 5e-4 and e30 / e60 == pytest.approx(4.0, rel=0.1)
+
+
+def test_liquid_over_slurry_joins_two_linear_profiles_with_the_kink_at_the_viscosity_jump():
+    """1 Pa of gas shear on 15 liquid cells over 15 of slurry at 900 K (37 % solid) over rigid metal: the stress is 1 Pa
+    all through, so the slurry shears at 1 Pa over its own viscosity at that rate (3.46 Pa s, 0.29 1/s) and the liquid at
+    1 Pa / 1.3 mPa s, the profiles joining at the jump. Each cell at its own shear rate makes it exact; the fixed point
+    takes two passes."""
+    sk = column_at([950.0] * 15 + [900.0] * 15 + [880.0] * 10)
+    mat = sk.material
+    f_s = 1.0 - float(mat.liquid_fraction(np.array([900.0]))[0])
+    f_b = 1.0 - float(mat.liquid_fraction(np.array([rheology.BRIDGE_LOW]))[0])
+    gamma = 0.0
+    for _ in range(50):
+        gamma = 1.0 / float(rheology.viscosity(900.0, f_s, gamma, MU_L, f_b))
+    book = skin.momentum_step(sk, STEADY, np.zeros(3), np.array([1.0, 0.0, 0.0]), sk.area)
+    s = heights(sk, 30)
+    h_s = sk.dz()[0][15:30].sum()
+    np.testing.assert_allclose(sk.U[0, 15:30, 0], gamma * s[15:], rtol=1e-4)
+    np.testing.assert_allclose(sk.U[0, :15, 0], gamma * h_s + (s[:15] - h_s) / MU_L, rtol=1e-6)
+    assert book.passes <= skin.SHEAR_PASSES
+
+
+def test_a_layer_spun_up_from_rest_by_a_sudden_shear_follows_stokes_first_problem():
+    """1 Pa suddenly applied at the top of 0.9 mm of liquid at rest, for 18.5 ms (a diffusion length sqrt(nu t) of
+    100 um): u(y, t) = (2 tau / mu) sqrt(nu t) ierfc(y / 2 sqrt(nu t)), y the depth. The top cell within 1e-3 at 200
+    steps, the error halving with the step (backward Euler; 1.3e-3, 6.6e-4, 3.3e-4 at 100, 200, 400 steps, measured
+    2026-10-19); the momentum is the shear's impulse exactly."""
+    from scipy.special import erfc
+    nu, T_end = MU_L / 2400.0, 0.0185
+    L = np.sqrt(nu * T_end)
+
+    def run(steps):
+        sk = column_at([950.0] * 80 + [880.0] * 5)
+        for _ in range(steps):
+            skin.momentum_step(sk, T_end / steps, np.zeros(3), np.array([1.0, 0.0, 0.0]), sk.area)
+        dz = sk.dz()[0][:80]
+        y = np.cumsum(dz) - 0.5 * dz
+        exact = 2.0 / MU_L * L * (np.exp(-(y / (2 * L)) ** 2) / np.sqrt(np.pi) - y / (2 * L) * erfc(y / (2 * L)))
+        return sk, abs(sk.U[0, 0, 0] / exact[0] - 1.0), np.abs(sk.U[0, :80, 0] - exact).max() / exact[0]
+    sk, top200, profile = run(200)
+    assert top200 < 1e-3 and profile < 1e-3
+    assert run(100)[1] / top200 == pytest.approx(2.0, rel=0.1)
+    assert sk.column_momentum()[0, 0] == pytest.approx(1.0 * sk.area[0] * T_end, rel=1e-12)
+
+
+def slurry_rate(sk, T, tau):
+    """The steady shear rate [1/s] of slurry at temperature T under the stress tau [Pa]: tau over its own viscosity at
+    that rate, by fixed point."""
+    mat = sk.material
+    f_s = 1.0 - float(mat.liquid_fraction(np.array([T]))[0])
+    f_b = 1.0 - float(mat.liquid_fraction(np.array([rheology.BRIDGE_LOW]))[0])
+    gamma = 0.0
+    for _ in range(50):
+        gamma = tau / float(rheology.viscosity(T, f_s, gamma, MU_L, f_b))
+    return gamma
+
+
+def test_a_floating_crust_slides_on_the_slurry_beneath_it_as_one_piece():
+    """Asha, 2026-10-20, option B: a crust (three rigid cells at 880 K) on 15 cells of slurry at 900 K over the rigid
+    metal joined to the substrate, 1 Pa of gas shear on its top: the stress is 1 Pa all through, the slurry shears in
+    Couette's line from the anchored front, and the crust rides its top as one piece, at the slurry's speed there."""
+    sk = column_at([880.0] * 3 + [900.0] * 15 + [880.0] * 5)
+    gamma = slurry_rate(sk, 900.0, 1.0)
+    book = skin.momentum_step(sk, STEADY, np.zeros(3), np.array([1.0, 0.0, 0.0]), sk.area)
+    s = heights(sk, 18)
+    h_s = sk.dz()[0][3:18].sum()
+    np.testing.assert_allclose(sk.U[0, 3:18, 0], gamma * s[3:], rtol=1e-4)
+    np.testing.assert_allclose(sk.U[0, :3, 0], gamma * h_s, rtol=1e-4)
+    assert np.ptp(sk.U[0, :3, 0]) < 1e-7 * sk.U[0, 0, 0]                    # one piece
+    assert np.all(sk.U[0, 18:] == 0.0) and np.all(book.shear[:, 0] > 0.0)
+
+
+def test_the_rigid_metal_joined_to_the_substrate_is_held_and_the_momentum_books_close():
+    """A crust over slurry and liquid over the anchored rigid metal, every cell moving at the start and pushed along the
+    surface: the anchored cells stop -- their momentum booked to the wall -- while the crust moves and takes the gas
+    shear, and the impulses add up to the momentum change."""
+    sk = column_at([880.0] + [920.0] * 10 + [902.0] * 10 + [880.0] * 5)
+    rng = np.random.default_rng(3)
+    sk.U = rng.normal(0.0, 0.05, sk.U.shape)
+    P0 = sk.column_momentum()
+    flowing = sk.valid & np.isfinite(skin.cell_viscosity(sk, 0.0))
+    assert np.array_equal(skin.anchored_cells(sk, flowing)[0], np.arange(26) >= 21)
+    book = skin.momentum_step(sk, 1e-3, np.array([2e4, -1e4, 0.0]), np.array([3.0, 0.0, 0.0]), sk.area)
+    assert np.all(sk.U[0, 21:] == 0.0) and np.all(sk.U[0, :21] != 0.0)
+    assert book.shear[0, 0] == pytest.approx(3.0 * sk.area[0] * 1e-3)
+    assert books_close(sk, P0, book)
+
+
+def test_a_column_whose_bottom_cell_flows_anchors_nothing():
+    """Slurry down to the base under a crust: no rigid cell is joined to the substrate, so the crust floats and the
+    slurry's bottom cell meets the rigid 3D metal beneath the base with no slip."""
+    sk = column_at([880.0] * 2 + [900.0] * 10)
+    flowing = sk.valid & np.isfinite(skin.cell_viscosity(sk, 0.0))
+    assert not skin.anchored_cells(sk, flowing).any()
+    gamma = slurry_rate(sk, 900.0, 1.0)
+    skin.momentum_step(sk, STEADY, np.zeros(3), np.array([1.0, 0.0, 0.0]), sk.area)
+    dz = sk.dz()[0]
+    s = dz[:12].sum() - (np.cumsum(dz[:12]) - 0.5 * dz[:12])                # heights above the base
+    np.testing.assert_allclose(sk.U[0, 2:12, 0], gamma * s[2:], rtol=1e-4)
```

- [ ] **Step 2: Run them to see them fail.** `"$PY" -m pytest tests/test_reentry_model_skin.py -q`. Expected:
  `AttributeError: module 'reentry_model.skin' has no attribute 'momentum_step'`.

- [ ] **Step 3: The momentum sub-step, and the drag at each cell's own shear rate,** in `reentry_model/skin.py`:

```diff
--- a/reentry_model/skin.py
+++ b/reentry_model/skin.py
@@ -793,23 +793,141 @@
     return np.where(sk.valid, eta, np.inf)
 
 
-def interface_drag(sk, gamma):
-    """Per column, the viscous drag per unit area between neighbouring cells' centres [Pa s/m], (n, K-1):
-    1 / (dz_i / (2 eta_i) + dz_j / (2 eta_j)), both cells' viscosities taken at the interface's shear rate gamma [1/s]
-    (n, K-1) -- the two half-cells in series, as the shear stress's continuity across the interface requires, the way
-    `_conduct` combines conductivities. A rigid cell resists without limit and adds nothing to the sum, so a flowing cell
-    drags on the rigid cell beneath it over its own half-thickness alone: no slip at the rigid front. Zero between two
-    rigid cells and across a missing one."""
+def series_drag(sk, eta):
+    """Per column, the viscous drag per unit area between neighbouring cells' centres [Pa s/m], (n, K-1), from the
+    cells' viscosities eta [Pa s] (n, K): 1 / (dz_i / (2 eta_i) + dz_j / (2 eta_j)), the two half-cells in series, as
+    the shear stress's continuity across the interface requires, the way `_conduct` combines conductivities. A rigid
+    cell resists without limit and adds nothing to the sum, so a flowing cell drags on the rigid cell beneath it over
+    its own half-thickness alone: no slip at the rigid front. Zero between two rigid cells and across a missing one."""
     if sk.m.shape[1] < 2:
         return np.zeros((sk.n, 0))
-    g = np.asarray(gamma, dtype=float) * np.ones((sk.n, sk.m.shape[1] - 1))
-    eta_up = cell_viscosity(sk, np.pad(g, ((0, 0), (0, 1)), mode="edge"))[:, :-1]
-    eta_dn = cell_viscosity(sk, np.pad(g, ((0, 0), (1, 0)), mode="edge"))[:, 1:]
     dz = sk.dz()
-    resist = 0.5 * dz[:, :-1] / eta_up + 0.5 * dz[:, 1:] / eta_dn
-    link = sk.valid[:, 1:]
-    flows = np.isfinite(eta_up) | np.isfinite(eta_dn)
-    return np.where(link & flows, 1.0 / np.where(link & flows, resist, 1.0), 0.0)
+    resist = 0.5 * dz[:, :-1] / eta[:, :-1] + 0.5 * dz[:, 1:] / eta[:, 1:]
+    link = sk.valid[:, 1:] & (np.isfinite(eta[:, :-1]) | np.isfinite(eta[:, 1:]))
+    return np.where(link, 1.0 / np.where(link, resist, 1.0), 0.0)
+
+
+def interface_drag(sk, gamma):
+    """`series_drag` with every cell's viscosity at its own shear rate gamma [1/s] (per cell, or one for all). Each
+    cell, not each interface: across a jump from liquid to slurry the liquid shears some thousand times faster than the
+    slurry, and a rate taken from the velocity difference between the two centres would give the slurry's half-cell the
+    liquid's rate (found while writing Task 4)."""
+    return series_drag(sk, cell_viscosity(sk, gamma))
+
+
+# -- momentum through the depth (sub-plan 21 Task 4) --------------------------------------------------------------------
+SHEAR_PASSES = 5                  # passes of the fixed point on the cells' shear rates at most (spec 6.4)
+SHEAR_TOL = 1.0e-3                # a pass that moves no flowing cell's viscosity by more than this, relatively, ends it
+CRUST_VISCOSITY = 1.0e9           # Pa s: a floating crust's, so that it moves as one piece (some 1e8 times the slurry's)
+
+
+@dataclass
+class MomentumBook:
+    """What a momentum sub-step did, per column: impulses [N s], (n, 3), which add up to the momentum change."""
+    push: np.ndarray         # of the push on the moving cells
+    shear: np.ndarray        # of the gas shear on the top cell, if it moves
+    wall: np.ndarray         # of the anchored rigid metal -- the skin's rigid cells joined to the substrate and the 3D
+                             # model beneath the base -- on the moving cells, and the momentum a cell loses on joining it
+    passes: int              # the fixed point's passes
+
+
+def anchored_cells(sk, flowing):
+    """The rigid cells joined to the substrate, (n, K): each column's run of rigid cells from its bottom cell up, which
+    rests on the rigid 3D model. A column whose bottom cell flows has none."""
+    n, K = sk.m.shape
+    j = np.arange(K)[None, :]
+    from_bottom = np.clip(sk.count[:, None] - 1 - j, 0, K - 1)                 # the j-th cell from the bottom
+    rigid_up = np.take_along_axis(sk.valid & ~flowing, from_bottom, axis=1) & (j < sk.count[:, None])
+    held = np.cumprod(rigid_up, axis=1).sum(axis=1)
+    return sk.valid & (j >= (sk.count - held)[:, None])
+
+
+def momentum_step(sk, dt, push, shear, A_top):
+    """One implicit sub-step dt [s] of every column's momentum through its depth (spec 6.4), each velocity component
+    apart:
+
+        m_i (u_i' - u_i) / dt = A [D_(i-1) (u'_(i-1) - u'_i) + D_i (u'_(i+1) - u'_i)] + G V_i + tau A_top [i the top],
+
+    A the column's area, D the drag between neighbouring cells (`series_drag`), G the push per unit volume [N/m3,
+    (n, 3)] -- every layer of a column feels the same (spec 6.4) -- V_i the cell's volume and tau the gas shear [Pa,
+    (n, 3)] on the top cell over its top area A_top [m2, (n,)]. The rigid metal joined to the substrate (`anchored_cells`)
+    is held at rest; a moving cell drags on it, and a flowing bottom cell on the rigid 3D metal beneath the base, with no
+    slip over its own half-thickness. A crust -- rigid cells with flowing metal between them and the substrate -- floats:
+    it moves as one piece (CRUST_VISCOSITY), dragged by the gas shear on its top, pushed like the rest and dragged by the
+    slurry beneath (Asha, 2026-10-20, option B; the column model has no sideways stress in a solid, so a crust over many
+    columns moves as column-sized plates, as a warm, barely solid crust would). Each flowing cell's viscosity is at its
+    own shear rate, the mean over its two faces of the stress there over its viscosity, found by fixed-point iteration
+    from the start-of-sub-step velocities: at most SHEAR_PASSES passes, ending when no flowing cell's viscosity moves by
+    more than SHEAR_TOL. Updates sk.U; returns a MomentumBook whose impulses add up to the columns' momentum change to
+    round-off."""
+    n, K = sk.m.shape
+    valid = sk.valid
+    rows = np.arange(n)
+    last = np.maximum(sk.count - 1, 0)
+    dz, A, m = sk.dz(), sk.area, sk.m
+    push = np.broadcast_to(np.asarray(push, dtype=float), (n, 3))
+    shear = np.broadcast_to(np.asarray(shear, dtype=float), (n, 3))
+    A_top = np.broadcast_to(np.asarray(A_top, dtype=float), (n,))
+    flowing = valid & np.isfinite(cell_viscosity(sk, 0.0))
+    anchored = anchored_cells(sk, flowing)
+    moving = valid & ~anchored                                            # the flowing cells and the floating crust
+    crust = moving & ~flowing
+    P0 = sk.column_momentum()
+    held = ((m * anchored)[..., None] * sk.U).sum(axis=1)                 # what anchored cells carried, now at rest
+    U_old = np.where(moving[..., None], sk.U, 0.0)
+    top_moves = moving[:, 0]
+    force = np.where(moving[..., None], push[:, None, :] * (dz * A[:, None])[..., None], 0.0)
+    force[:, 0] += np.where(top_moves[:, None], shear * A_top[:, None], 0.0)
+    bottom_flows = flowing[rows, last] & (sk.count > 0)
+
+    def viscosities(gamma):
+        return np.where(crust, CRUST_VISCOSITY, cell_viscosity(sk, gamma))
+
+    def drags(eta):
+        Db = np.where(bottom_flows, eta[rows, last] / (0.5 * np.maximum(dz[rows, last], DZ_FLOOR)), 0.0)
+        return series_drag(sk, eta), np.where(np.isfinite(Db), Db, 0.0)
+
+    def solve(eta):
+        D, Db = drags(eta)
+        diag = np.where(valid, m / dt, 1.0)
+        diag[:, :-1] += A[:, None] * D
+        diag[:, 1:] += A[:, None] * D
+        diag[rows, last] += A * Db
+        off = -A[:, None] * D
+        sub, sup = np.where(moving[:, 1:], off, 0.0), np.where(moving[:, :-1], off, 0.0)
+        diag = np.where(moving, diag, 1.0)
+        rhs = np.where(moving[..., None], (m / dt)[..., None] * U_old + force, 0.0)
+        U = np.stack([thomas(sub, diag, sup, rhs[..., c]) for c in range(3)], axis=-1)
+        return np.where(moving[..., None], U, 0.0), D, Db
+
+    def rates(U, D, Db, eta):
+        """Every flowing cell's shear rate: the mean over its faces of the stress there over its viscosity."""
+        stress = np.zeros((n, K + 1))                                     # at every face, the top face first
+        stress[:, 1:K] = D * np.linalg.norm(U[:, :-1] - U[:, 1:], axis=2)
+        stress[:, 0] = np.where(top_moves, np.linalg.norm(shear, axis=1), 0.0)
+        stress[rows, last + 1] = Db * np.linalg.norm(U[rows, last], axis=1)
+        face = 0.5 * (stress[:, :-1] + stress[:, 1:])
+        return np.where(flowing, np.divide(face, eta, out=np.zeros_like(face), where=np.isfinite(eta)), 0.0)
+
+    eta = viscosities(0.0)
+    eta = viscosities(rates(U_old, *drags(eta), eta))                       # at the start-of-sub-step rates
+    for passes in range(1, SHEAR_PASSES + 1):
+        U, D, Db = solve(eta)
+        eta_new = viscosities(rates(U, D, Db, eta))
+        change = np.abs(eta_new[flowing] - eta[flowing]) / eta[flowing]
+        if change.size == 0 or change.max() <= SHEAR_TOL:
+            break
+        eta = eta_new
+    wall = np.zeros((n, 3))
+    if K > 1:
+        down = moving[:, :-1] & anchored[:, 1:]                           # a moving cell on the anchored metal beneath
+        wall -= ((A[:, None] * D * down)[..., None] * U[:, :-1]).sum(axis=1)
+    wall -= (A * Db)[:, None] * U[rows, last]
+    book = MomentumBook(push=dt * np.where(moving[..., None], push[:, None, :] * (dz * A[:, None])[..., None], 0.0).sum(axis=1),
+                        shear=dt * np.where(top_moves[:, None], shear * A_top[:, None], 0.0),
+                        wall=dt * wall - held, passes=passes)
+    sk.U = U
+    return book
 
 
 # -- the base (Task 7) -------------------------------------------------------------------------------------------------
```

- [ ] **Step 4: Run them to see them pass.** Expected: `tests/test_reentry_model_skin.py` 40 passed (with Task 3's
  drag tests, which hold at one rate for all cells); the unit tier 432 passed, 1 skipped.

- [ ] **Step 5: Commit.**

```bash
git add reentry_model/skin.py tests/test_reentry_model_skin.py
git commit -m "Skins: every column's momentum through its depth in one implicit sub-step, each cell's viscosity at its own shear rate by fixed point, the rigid metal joined to the substrate held with no slip and a crust floating as one piece, the books exact (sub-plan 21 Task 4)"
```

---

### Task 5: The push

**Files:**
- Modify: `reentry_model/surface_flow.py` (`SurfaceFlowResult.dp_ds`), `reentry_model/skin.py` (`along_meridian`,
  `push`)
- Test: `tests/test_reentry_model_surface_flow.py`, `tests/test_reentry_model_skin.py`

**Interfaces:**
- Consumes: Task 2's `Columns.gradient`, `Columns.curvature` and `Columns.top_points`; the surface flow's
  `SurfaceFlowResult`.
- Produces: `SurfaceFlowResult.dp_ds` (n,) [Pa/m], the wall pressure's gradient along the meridian on the windward
  face, zero leeward; `skin.along_meridian(normals, v_hat, magnitude) -> (n, 3)`, vectors in each corner's tangent plane
  away from the stagnation point (for the pressure term here and the gas shear, `SurfaceFlowResult.tau`, in Task 6);
  `skin.push(columns, thickness, dp_ds, g, v_hat, rho_l, sigma) -> (G (n, 3) [N/m3], terms)`, the terms a dict of
  `pressure`, `deceleration`, `thickness` and `capillary`.

The push on every layer of a column (spec §6.4), G = −∇_s p_w + ρ_l g_t + ρ_l (g·n) ∇_s H − σ ∇_s κ, each term in the
column's own tangent plane, where Task 4 solves for the velocities. g is the acceleration the metal feels in the body's
frame — minus the body's own; for a body slowed by drag, the deceleration along the flight direction — so it pulls the
metal toward the nose along the surface and out of the front face: the spec's −ρ_l a_t and ρ_l a_n (it writes the
first with the body's acceleration and the second with the deceleration's pull), and the element model's
−ρ_l W sin θ on the windward face. Unlike the element model it acts on the leeward face too, where it pulls the metal
toward the equator. H is each column's thickness and κ the mean curvature of the top surface the thicknesses make.
(Task 6b replaces this κ by its thin-film form, on the top surface's cotangents, and applies surface tension only where a
column's top flows: Task 2's obtuse-corner rule made this term run backwards at those corners, Asha, 2026-10-21.)

**The pressure term takes the surface flow's gradient, not the corners'.** Differentiating each corner's own p_w would
carry the scatter of the corners' normals — 0.06° off radial on average on the 2 mm sphere, 0.87° at worst — amplified
some 25 times, the radius over the spacing: a Newtonian pressure's push came out 1.1 % off on average and 13 % at worst
that way (1.5 % and 10 % on the 4 mm test sphere), against 0.1 % and 1.2 % with each corner's exact angle; normals
taken from the curvature vector, or averaged over the neighbours, did no better (measured 2026-10-20,
`prototype/work-2026-10-18-column-flow/measure/`). The surface flow already differentiates p_w(θ) on its own fine bins
along the meridian for the element model's G; it now also returns that gradient, `dp_ds`, on its own, with every
field it returned before unchanged to the last bit at 20, 80 and 200 s of the 100 mm flight, and the push takes
−dp_ds along the meridian away from the stagnation point.

**Sizes at peak heating** on the 100 mm sphere (estimated, not run: 93 s, 60 km, a dynamic pressure of 6.1 kPa from the
level-3 physics flight, a stagnation pressure of 1.84 times it and a drag coefficient of 0.92): the pressure term peaks
at p_stag / R, about 2.2·10⁵ Pa/m at 45° from the nose; the deceleration, 30 m/s², gives at most ρ_l g = 7.2·10⁴ Pa/m
at 90°, pulling the other way. Their ratio, 1.84 × 4ρ_s / (3 ρ_l C_D) ≈ 3.1, is the same for a solid sphere of any size
at any point of its flight, since both scale with the dynamic pressure. The thickness and capillary terms matter at
short wavelengths: a 0.1 mm bump 4 mm long gives about 3·10⁵ Pa/m of capillary push, which is why Task 6 takes both
implicitly; the sphere's own mesh, with the curvature exact to 10⁻⁴, gives at most 11 Pa/m of spurious capillary push
(0.3 on average).

The two thickness terms are evaluated at the thicknesses given; their implicit treatment, linearised in the column
heights, belongs to Task 6's solve, and so do the Rayleigh–Taylor growth and the surface-tension cutoff's decay, which
need the thickness to evolve (Asha, 2026-10-20). Nothing calls `push` yet: the skin body evaluates the surface flow on
its corner surface and drives `momentum_step` with the push and the gas shear when the flow joins its sub-steps
(Task 6), so no result changes.

- [ ] **Step 1: The surface flow's test** — append to `tests/test_reentry_model_surface_flow.py`:

```diff
--- a/tests/test_reentry_model_surface_flow.py
+++ b/tests/test_reentry_model_surface_flow.py
@@ -126,3 +126,15 @@
     assert r.branch == sf.BRANCH_FREE_MOLECULAR and np.all(np.isnan(r.delta_a)) and np.all(r.u_e == 0.0)
     assert np.isnan(r.p_stag) and np.all(r.p_w == r.p_inf) and np.all(r.closure == sf.CLOSURE_COUETTE)
     assert r.tau[1] == pytest.approx(a.freestream.rho * a.V ** 2 * math.sin(theta[1]) * math.cos(theta[1])) and r.tau[3] == 0.0
+
+
+def test_g_is_the_pressure_gradient_along_the_meridian_and_the_deceleration():
+    """Sub-plan 21 Task 5: the surface flow hands the skins the pressure part of G on its own, dp_ds; on the windward
+    face G = -dp_ds - rho_l W sin(theta), with p_w falling from the stagnation point, and leeward both are zero."""
+    state = aero_state_at(80.0)
+    theta = np.radians(np.linspace(0.0, 150.0, 61))
+    r = sf.SurfaceFlow().evaluate(state, theta, 0.05, 2400.0)
+    windward = theta <= 0.5 * math.pi
+    np.testing.assert_allclose(r.G[windward], -r.dp_ds[windward] - 2400.0 * r.deceleration * np.sin(theta[windward]),
+                               rtol=1e-12, atol=1e-9)
+    assert np.all(r.dp_ds[windward][1:-1] < 0.0) and np.all(r.dp_ds[~windward] == 0.0) and np.all(r.G[~windward] == 0.0)
```

- [ ] **Step 2: The push's tests** — append to `tests/test_reentry_model_skin.py`:

```diff
--- a/tests/test_reentry_model_skin.py
+++ b/tests/test_reentry_model_skin.py
@@ -729,3 +729,89 @@
     dz = sk.dz()[0]
     s = dz[:12].sum() - (np.cumsum(dz[:12]) - 0.5 * dz[:12])                # heights above the base
     np.testing.assert_allclose(sk.U[0, 2:12, 0], gamma * s[2:], rtol=1e-4)
+
+
+# -- the push (sub-plan 21 Task 5) --------------------------------------------------------------------------------------
+V_HAT = np.array([1.0, 0.0, 0.0])
+RHO_L, SIGMA = 2400.0, 0.8
+
+
+def regular_grid(n_side=81, h=1e-3):
+    """Columns on a regular flat square grid of n_side corners a side, h apart, in z = 0, normals +z."""
+    from reentry_model import columns
+    x, y = np.meshgrid(np.arange(n_side) * h, np.arange(n_side) * h, indexing="ij")
+    pts = np.column_stack([x.ravel(), y.ravel(), np.zeros(x.size)])
+    idx = np.arange(n_side * n_side).reshape(n_side, n_side)
+    a, b, c, d = idx[:-1, :-1].ravel(), idx[1:, :-1].ravel(), idx[1:, 1:].ravel(), idx[:-1, 1:].ravel()
+    return columns.Columns(pts, np.vstack([np.column_stack([a, b, c]), np.column_stack([a, c, d])])), pts
+
+
+def test_along_the_meridian_is_in_the_tangent_plane_away_from_the_stagnation_point(coarse_sphere_mesh):
+    from test_reentry_model_columns import sphere_columns
+    col, _ = sphere_columns(coarse_sphere_mesh)
+    s = skin.along_meridian(col.normal, V_HAT, np.ones(col.n))
+    side = np.linalg.norm(V_HAT - (col.normal @ V_HAT)[:, None] * col.normal, axis=1) > 1e-6
+    assert np.abs(np.einsum("ij,ij->i", s, col.normal)).max() < 1e-12
+    np.testing.assert_allclose(np.linalg.norm(s[side], axis=1), 1.0, rtol=1e-12)
+    assert np.all(s[side] @ V_HAT < 0.0)                                       # away from the stagnation point
+    poles = skin.along_meridian(np.array([V_HAT, -V_HAT]), V_HAT, np.ones(2))
+    assert np.all(poles == 0.0)
+
+
+def test_the_deceleration_pulls_the_metal_toward_the_nose_everywhere(coarse_sphere_mesh):
+    """30 m/s2 along the flight direction on the 4 mm sphere: rho_l g sin(theta) toward the nose -- the element model's
+    -rho_l W sin(theta) on the windward face, and on the leeward face too, where it pulls the metal toward the equator."""
+    from test_reentry_model_columns import sphere_columns
+    col, _ = sphere_columns(coarse_sphere_mesh)
+    G, terms = skin.push(col, np.full(col.n, 8e-4), np.zeros(col.n), 30.0 * V_HAT, V_HAT, RHO_L, SIGMA)
+    sin = np.linalg.norm(V_HAT - (col.normal @ V_HAT)[:, None] * col.normal, axis=1)
+    away = skin.along_meridian(col.normal, V_HAT, np.ones(col.n))
+    np.testing.assert_allclose(np.einsum("ij,ij->i", terms["deceleration"], away), -RHO_L * 30.0 * sin, rtol=1e-9, atol=1e-9)
+    assert np.all(terms["pressure"] == 0.0) and np.abs(terms["thickness"]).max() < 1e-9          # uniform: round-off
+
+
+def test_the_pressure_term_is_the_surface_flow_s_gradient_away_from_the_stagnation_point(coarse_sphere_mesh):
+    """A Newtonian wall pressure falling from the stagnation point, dp/ds = -2 p0 sin cos / R: its push points away
+    from the stagnation point at the surface flow's magnitude -- the input's own, not re-differentiated."""
+    from test_reentry_model_columns import sphere_columns
+    col, _ = sphere_columns(coarse_sphere_mesh)
+    cos = np.clip(col.normal @ V_HAT, 0.0, 1.0)
+    dp_ds = -2.0 * 11e3 * np.sqrt(1.0 - cos ** 2) * cos / 0.05
+    _, terms = skin.push(col, np.full(col.n, 8e-4), dp_ds, np.zeros(3), V_HAT, RHO_L, SIGMA)
+    away = skin.along_meridian(col.normal, V_HAT, np.ones(col.n))
+    np.testing.assert_allclose(np.einsum("ij,ij->i", terms["pressure"], away), -dp_ds, rtol=1e-12, atol=1e-9)
+
+
+def test_a_thicker_column_under_a_pull_out_of_the_surface_draws_the_metal_toward_it():
+    """A flat layer thickening by 1 mm a metre under a pull of 100 m/s2 out of its surface: rho_l (g . n) grad H."""
+    col, pts = regular_grid(21)
+    H = 5e-4 + 1e-3 * pts[:, 0]
+    _, terms = skin.push(col, H, np.zeros(col.n), np.array([0.0, 0.0, 100.0]), V_HAT, RHO_L, SIGMA)
+    np.testing.assert_allclose(terms["thickness"], np.tile([RHO_L * 100.0 * 1e-3, 0.0, 0.0], (col.n, 1)), atol=1e-9)
+    assert np.abs(terms["deceleration"]).max() == 0.0                          # the pull lies along the normal
+
+
+def test_a_bump_s_surface_tension_pushes_the_metal_from_its_crest_converging_at_second_order():
+    """A cosine bump of 20 um on a 0.5 mm layer: kappa = -H'' for small slopes, so -sigma grad kappa = sigma H''' =
+    sigma eps k^3 sin(kx). Measured 2026-10-20: 2.45 % of its peak off at 20 spacings a wavelength, 0.62 % at 40."""
+    def error(lam):
+        col, pts = regular_grid()
+        k, eps = 2.0 * np.pi / lam, 2e-5
+        _, terms = skin.push(col, 5e-4 + eps * np.cos(k * pts[:, 0]), np.zeros(col.n), np.zeros(3), V_HAT, RHO_L, SIGMA)
+        exact = SIGMA * eps * k ** 3 * np.sin(k * pts[:, 0])
+        inner = np.all((pts[:, :2] > 16e-3) & (pts[:, :2] < 64e-3), axis=1)
+        assert np.abs(terms["capillary"][inner, 1:]).max() < 1e-9 * np.abs(exact).max()
+        return np.abs(terms["capillary"][inner, 0] - exact[inner]).max() / np.abs(exact).max()
+    e20, e40 = error(20e-3), error(40e-3)
+    assert e40 < 0.01 and e20 / e40 == pytest.approx(4.0, rel=0.1)
+
+
+def test_a_sphere_s_own_mesh_drives_no_capillary_flow_to_speak_of(coarse_sphere_mesh):
+    """Uniform 0.8 mm skins on the 4 mm sphere: the curvature is uniform, so the capillary push is the mesh's error
+    alone -- 11 Pa/m at worst, 0.28 on average (measured 2026-10-20), against 2e5 Pa/m from the pressure at peak
+    heating (fact 153)."""
+    from test_reentry_model_columns import sphere_columns
+    col, _ = sphere_columns(coarse_sphere_mesh)
+    _, terms = skin.push(col, np.full(col.n, 8e-4), np.zeros(col.n), np.zeros(3), V_HAT, RHO_L, SIGMA)
+    size = np.linalg.norm(terms["capillary"], axis=1)
+    assert size.max() < 20.0 and size.mean() < 0.5
```

- [ ] **Step 3: Run them to see them fail.** `"$PY" -m pytest tests/test_reentry_model_surface_flow.py
  tests/test_reentry_model_skin.py -q`. Expected: `AttributeError: 'SurfaceFlowResult' object has no attribute 'dp_ds'`
  and `AttributeError: module 'reentry_model.skin' has no attribute 'along_meridian'`.

- [ ] **Step 4: The pressure gradient on its own** in `reentry_model/surface_flow.py`:

```diff
--- a/reentry_model/surface_flow.py
+++ b/reentry_model/surface_flow.py
@@ -213,6 +213,8 @@
     tau_continuum: np.ndarray
     tau_fm: np.ndarray
     G: np.ndarray
+    dp_ds: np.ndarray = None      # the wall pressure's gradient along the meridian [Pa/m], from the fine bins (the
+                                  # pressure part of G; the skins' push takes it on its own, sub-plan 21 Task 5)
 
     @property
     def branch_name(self):
@@ -283,7 +285,7 @@
             return SurfaceFlowResult(BRANCH_FREE_MOLECULAR, state.kn, float("nan"), state.ma, float("nan"), p_inf,
                                      float("nan"), decel, self.gamma_pm, np.full(n, p_inf), zero, zero, zero, zero, zero,
                                      zero, nan, nan, np.full(n, np.inf), np.full(n, CLOSURE_COUETTE), np.ones(n, dtype=bool),
-                                     zero, zero, zero, zero)
+                                     zero, zero, zero, zero, zero)
         G_inertial = -rho_liquid * decel * sin
         if branch == BRANCH_FREE_MOLECULAR:
             # Schaaf-Chambre from the freestream directly: no compression model, no edge state, no boundary layer
@@ -292,7 +294,7 @@
                                      self.gamma_pm, np.full(n, p_inf), zero, zero, np.full(n, fs.rho), np.full(n, fs.T),
                                      zero, zero, nan, np.full(n, lam_inf), np.full(n, lam_inf / radius),
                                      np.full(n, CLOSURE_COUETTE), np.ones(n, dtype=bool), tau_fm, zero, tau_fm,
-                                     np.where(windward, G_inertial, 0.0))
+                                     np.where(windward, G_inertial, 0.0), zero)
         # shock-layer construction: the reservoir (also the bridge's continuum endpoint in the merged branch)
         stag = self.air.stagnation(fs.rho, fs.T, V)
         re_shock = float(fs.rho * V * radius / stag.mu)
@@ -328,9 +330,10 @@
         closure = np.where(windward, closure, CLOSURE_COUETTE)
         dp_ds = np.gradient(p_b, s_b)                              # from the actual p_w(theta), not a Newtonian formula
         G = np.where(windward, -at(dp_ds) + G_inertial, 0.0)
+        dp_ds_w = np.where(windward, at(dp_ds), 0.0)
         mach_e = np.where(a_e > 0.0, u_e / np.maximum(a_e, 1e-300), 0.0)
         u_e, rho_e, mu_e = (np.where(np.isfinite(x), x, 0.0) for x in (u_e, rho_e, mu_e))
         return SurfaceFlowResult(branch, float(state.kn), re_shock, float(state.ma), float(stag.p), p_inf, float(phi_star),
                                  decel, self.gamma_pm, p_w, u_e, u_e * slip, rho_e, T_e, mu_e, mach_e,
                                  np.where(windward, delta_a, np.nan), lam_w, np.where(windward, kn_local, np.inf),
-                                 closure, flagged | ~windward, tau, np.where(windward, tau_c, 0.0), tau_fm, G)
+                                 closure, flagged | ~windward, tau, np.where(windward, tau_c, 0.0), tau_fm, G, dp_ds_w)
```

- [ ] **Step 5: The push** in `reentry_model/skin.py`:

```diff
--- a/reentry_model/skin.py
+++ b/reentry_model/skin.py
@@ -930,6 +930,50 @@
     return book
 
 
+# -- the push (sub-plan 21 Task 5) -------------------------------------------------------------------------------------
+def along_meridian(normals, v_hat, magnitude):
+    """Vectors of the given magnitudes (n,) in each corner's tangent plane, pointing along the meridian away from the
+    stagnation point -- the direction of minus the flight direction's tangential part, the way the gas flows and the
+    surface flow's s runs; zero where that part vanishes (the stagnation and rear points)."""
+    v = np.asarray(v_hat, dtype=float) / np.linalg.norm(v_hat)
+    away = -(v[None, :] - (normals @ v)[:, None] * normals)
+    size = np.linalg.norm(away, axis=1)
+    unit = np.divide(away, size[:, None], out=np.zeros_like(away), where=size[:, None] > 1e-12)
+    return np.asarray(magnitude, dtype=float)[:, None] * unit
+
+
+def push(columns, thickness, dp_ds, g, v_hat, rho_l, sigma):
+    """The push along the surface on every layer of every column [N/m3], (n, 3), and its four terms (spec 6.4), each in
+    the column's own tangent plane -- the base's, where `momentum_step` solves for the velocities:
+
+        pressure      -dp_w/ds along the meridian away from the stagnation point, dp_ds [Pa/m, (n,)] the surface flow's
+                      (`SurfaceFlowResult.dp_ds`), differentiated on its own fine bins of the angle. The gradient of the
+                      corners' own p_w would carry the scatter of their normals' angles, amplified some 25 times (the
+                      radius over the spacing) by the differentiation: on the 2 mm sphere a Newtonian pressure's push
+                      came out 1.1 % off on average and 13 % at worst that way, from normals 0.06 deg off radial on
+                      average (measured 2026-10-20).
+        deceleration  rho_l g_t, g [m/s2, (3,) or (n, 3)] the acceleration the metal feels in the body's frame, minus the
+                      body's own -- for a body slowed by drag, the deceleration along the flight direction v_hat -- so
+                      that it pulls the metal forward, toward the nose along the surface and out of the front face;
+                      the spec's -rho_l a_t, and the element model's -rho_l W sin(theta) on the windward face.
+        thickness     rho_l (g . n) grad_s H, H [m, (n,)] the columns' thicknesses: a thicker column under a pull out
+                      of the surface has the lower pressure beneath it and draws metal in (the front face's
+                      Rayleigh-Taylor mechanism).
+        capillary     -sigma grad_s kappa, kappa the mean curvature (`Columns.curvature`) of the top surface the
+                      thicknesses make (`Columns.top_points`).
+
+    The two thickness terms are evaluated at the given H; Task 6 takes them implicitly in the column heights."""
+    H = np.asarray(thickness, dtype=float)
+    normals = columns.normal
+    g = np.broadcast_to(np.asarray(g, dtype=float), (columns.n, 3))
+    g_n = np.einsum("ij,ij->i", g, normals)
+    terms = {"pressure": -along_meridian(normals, v_hat, dp_ds),
+             "deceleration": rho_l * (g - g_n[:, None] * normals),
+             "thickness": rho_l * g_n[:, None] * columns.gradient(H),
+             "capillary": -sigma * columns.gradient(columns.curvature(columns.top_points(H)))}
+    return sum(terms.values()), terms
+
+
 # -- the base (Task 7) -------------------------------------------------------------------------------------------------
 def hand_down(sk, m_min, T_rigid):
     """Rigid cells at the bottom of each column (at or below T_rigid, contiguous from the bottom) pass to the 3D model,
```

- [ ] **Step 6: Run them to see them pass.** Expected: `tests/test_reentry_model_surface_flow.py` 8 passed and
  `tests/test_reentry_model_skin.py` 46 passed; the unit tier 439 passed, 1 skipped.

- [ ] **Step 7: Commit.**

```bash
git add reentry_model/surface_flow.py tests/test_reentry_model_surface_flow.py reentry_model/skin.py tests/test_reentry_model_skin.py
git commit -m "Skins: the push along the surface on every layer -- the surface flow's pressure gradient along the meridian, the deceleration, the thickness and the capillary terms -- in each column's tangent plane (sub-plan 21 Task 5)"
```

---

### Task 6a: Sideways transport at given velocities

**Files:**
- Modify: `reentry_model/columns.py` (`Columns.face_vector`), `reentry_model/skin.py` (`SLIVER`, `TransportBook`,
  `tangent_frames`, `_fractions`, `_pieces`, `_upwind_order`, `_factor`, `transport`)
- Test: `tests/test_reentry_model_columns.py`, `tests/test_reentry_model_skin.py`

**Interfaces:**
- Consumes: Task 2's `Skins.U`, Task 3's `cell_viscosity`, Task 4's `anchored_cells`; `Columns.edges`, `.normal`,
  `.area`; `Skins.dz`, `Skins._install` and `Skins.regrid`.
- Produces: `Columns.face_vector` (E, 3) [m], each dual face's normal times its length, from the edge's first corner
  to its second; `skin.transport(sk, columns, dt, regrid=True) -> TransportBook`, which moves the cells at their
  velocities `sk.U` and, with `regrid`, splits and merges them; `TransportBook(mass, enthalpy, momentum, wall, courant)`,
  per column what it gained from its neighbours less what it gave them — mass (n,) [kg], enthalpy (n,) [J], momentum
  (n, 3) [kg m/s] as it left the donors — the wall's impulse (n, 3) [N s], the normal parts removed, and the largest
  cell's outflow rate times dt; `skin.tangent_frames(normals) -> (e1, e2)`, each (n, 3); `skin.SLIVER` (10⁻¹²); and,
  private, `_fractions(depth) -> (n, K+1)`, `_pieces(sig, tau) -> (row, cell of sig, cell of tau, length)`,
  `_upwind_order(src, dst, n) -> (order, upwind)` and `_factor(A, upwind)`.

The metal moves sideways at the velocities the cells already have (spec §6.5); in the tests they are prescribed, and
the skin body drives them with Task 4's `momentum_step` in 6c. The transport is its own step in each sub-step, after
the momentum step (ruling, 2026-10-21); the conduction stays in Task 1's solve. Across each face, cell k of the donor
column a sends c m_k′ toward its neighbour b, with c = (u_k·f)⁺ / A_a [1/s] — the velocity's component across the face
times the face's length, per unit of the column's area — and m_k′ the cell's mass at the sub-step's end; the metal
carries the cell's specific enthalpy and its velocity. Every cell then satisfies m_k′ (1 + dt Σc) = m_k + dt (what
arrives): (I + dt C) m′ = m, the film's runoff scheme (`film.Runoff`) on the cells. The matrix is an M-matrix whose
columns sum to one, so the masses stay positive, the books close to round-off and nothing limits dt; the enthalpy and
the trial pass's sensitivity S solve the same system, and so does the momentum, with every arrival projected onto its
receiver's plane. Only the cells that send or receive are solved for: with every velocity zero, or only rigid metal
joined to the substrate given one, nothing changes at all.

**Where the metal goes** (spec §6.5, §14 assumptions 3 and 4). A flowing donor cell's span of fractions of the donor's
flowing depth joins the receiver's cells between the same fractions of the receiver's flowing layer — the run of
flowing cells from its top — each in proportion to the overlap. Where the receiver's top is rigid (anchored metal, or
a crust over flowing metal) everything arriving lands on top of its top cell, as a new cell; so does a floating
crust's metal, wherever it goes. Rigid metal joined to the substrate neither leaves nor receives, whatever velocity it
holds. After the exchange the cells are split and merged as everywhere else (spec §6.2): a deposit thinner than half
its target merges into the cell beneath it, as any thin cell does.

**The faces.** A velocity crosses each dual face at the face's vector, its two segments' normals times their lengths
(`Columns.face_vector`), not the edge's direction times the face's length: a column's faces close exactly on a plane,
so a uniform flow carries as much into every column as out of it, where the edge's direction gathered or spread a
column's metal at up to 539 1/s under 1 m/s on the columns' tests' jittered plane (2.9·10⁻¹³ with the faces' vectors).
On the 4 mm test sphere the faces close in the tangent plane to 3.7 % of a face length at worst and 0.03 % on average;
the rest is along the normal, the curvature, which a tangential velocity does not see.

**The momentum is projected at every arrival.** Each cell's momentum is held in two components in its own plane
(`tangent_frames`), and every arrival enters its receiver through the 2×2 projection between the two planes, so a
second system of twice the size is solved with the same structure; any normal part a cell held is booked to the wall
too. Projecting once at the sub-step's end instead — the Cartesian components carried by the scalar system, each cell
then projected — costs nothing extra but carries metal that crosses several columns within a sub-step straight on, off
the curved surface: on the 100 mm sphere a uniform meridian flow kept 99.7 % of its kinetic energy in one sub-step at a
Courant number of 1 against 99.8 % projected at arrival, 90 % against 96 % at 9.6, and 24 % against 49 % at 96, where
the two differed by up to 1.8 times the largest cell momentum. Courant numbers near 100 are what the flight will see:
0.5 mm of liquid under the 2.2·10⁵ Pa/m of Task 5's pressure push moves at about 20 m/s at its top (G h²/2μ), 100 columns
of 2 mm in a 10 ms sub-step.

**The solve, measured** (2026-10-21, `prototype/work-2026-10-18-column-flow/measure/task6a/`, on the 100 mm sphere's
skins, 9 417 columns of 10 µm cells, flowing at 0.1 m/s along the meridian, every cell liquid and moving: 584 000
unknowns, 2.3 million arrivals). A sparse direct solve in the cells' own order took 19.6 s to factor; in SuperLU's
fill-reducing order (COLAMD) 0.66 s; minimum-degree 54 s. Taken in upwind order, every donor before its receivers, the
system is triangular but for its cycles where flows meet (here at most 7 cells, near the rear stagnation point), and
SuperLU with no pivoting factors it in 0.06 s with no fill at all — and without pivoting an M-matrix's elimination keeps
every sign, so the masses come out non-negative to the last bit. The order is scipy's strong components, 0.06 s, which
its algorithm (Pearce's) numbers in reverse topological order; that is how scipy numbers them, not what it promises,
so `_upwind_order` checks it every call and `_factor` falls back to COLAMD with pivoting, as exact and slower, should
it fail. Not built (estimated): an upwind sweep written in numpy would need the levels of that order, a hundred or more
relaxations over the 2.3 million arrivals; Gauss–Seidel or a Krylov iteration gains about ν/(1 + ν) a sweep, so far
above a Courant number of one it needs many, and its books close only to its tolerance. The pieces where a donor's and
a receiver's partitions cut each other come from a rank merge (`_pieces`, 0.13 s against 0.22 s for a row-wise sort).

Whole transport sub-steps (without the regrid, which the sub-step does anyway): 0.95 s with all 584 000 cells moving,
of which the momentum's system 0.30 s and the pieces 0.13 s; 0.27 s in a melting state — the front hemisphere's top
0.3 mm liquid over rigid metal, 124 000 of 553 000 cells moving along a half-parabola in depth — the same at Courant
numbers of 1, 9.6 and 96. On the same machine a coupled conduction sub-step on the sphere's skins at 700 K takes
0.6–0.8 s and a momentum sub-step on the melting state 0.54 s; the regrid after the transport 0.09–0.25 s, more at high
Courant numbers, where more cells change size.

**Rulings** (2026-10-21, within spec §6.5 and the rulings of 2026-10-20/21):

- *The face's vector, not the edge's direction times the face's length* — the face is two segments in two triangles,
  and only their vectors close around a column — if wrong, nothing: on symmetric triangles the two are the same.
- *The donor's flowing depth is all its flowing cells, top down; a moving rigid cell is no part of it and lands on top*
  — the spec's "flowing layer"; where flowing metal lies under a floating rigid layer (melt beneath a refrozen deposit),
  it is matched by its share of the flowing metal and the rigid layer lands on the receiver's top — if wrong, such a
  layer arrives above metal that was above it; rare, and transient.
- *A donor cell's pieces are normalised to sum to one, and pieces under 10⁻¹² of the cell are dropped* — the pieces'
  lengths are differences of cumulative fractions, and their rounding left a cell's outflow distributed 5.5·10⁻¹⁴ short
  or over, which grew to 7.6·10⁻¹³ of the mass over ten sub-steps; the slivers are the rounding of a boundary the two
  partitions share — if wrong, a piece of 10⁻¹² of a cell is moved to its neighbour.
- *The diagonal is the sum of the arrivals made* — every cell sends exactly what its arrivals take, so the books close
  whatever the pieces — if wrong, nothing: it equals dt Σc to round-off.
- *The momentum is projected at every arrival* (above) — the spec's text, and, at high Courant numbers, nearer the
  parallel transport it stands for — if wrong, 0.30 s of the worst case's 0.95 s is spent for nothing.
- *The regrid is optional* (`regrid=True` by default) — 6c calls the transport before the sub-step's own regrid and
  passes False, so the cells are split and merged once — if wrong, a sub-step regrids twice.

What the tests measure (2026-10-21): a pulse half again as thick as the liquid around it and 50 K warmer, 3 mm wide, on
1 mm columns at 0.1 m/s for 0.1 s, keeps its mass and heat, its centroids move the 10 mm to 6·10⁻⁸, 10⁻⁵ and
8.3·10⁻⁴ mm at Courant numbers 0.5, 1 and 2 (the scheme's long upwind tail reaching the grid's ends), and its variance
grows by the one-dimensional implicit upwind scheme's h² ν (1 + ν) a sub-step, less 4.2, 3.2 and 2.4 % — the numerical
diffusion u h (1 + ν)/2; depth fractions and both arrival rules exactly; on the 4 mm sphere's mixed skins every book to
4·10⁻¹⁵ of its scale, the wall taking 2.5 % of the largest column's momentum; a 10 s sub-step (Courant numbers to
875) positive, within the old enthalpies and speeds, conservative and regridded; on the sphere a velocity carried into
three neighbours 4.5–4.9° away loses 2.3–7.7 % of its momentum to the wall.

- [ ] **Step 1: The faces' test** — append to `tests/test_reentry_model_columns.py`:

```diff
--- a/tests/test_reentry_model_columns.py
+++ b/tests/test_reentry_model_columns.py
@@ -130,3 +130,31 @@
 def test_the_curvature_of_a_plane_is_zero():
     pts, faces, inner = flat_grid()
     assert np.abs(columns.Columns(pts, faces).curvature()[inner]).max() < 1e-9
+
+
+# -- the dual faces as vectors (sub-plan 21 Task 6a) --------------------------------------------------------------------
+def test_the_dual_faces_close_every_column_so_a_uniform_flow_neither_gathers_nor_spreads_metal(coarse_sphere_mesh):
+    """`face_vector`: each dual face's two segments, each its normal in its triangle times its length, summed and
+    pointing from the edge's first corner to its second -- what a velocity crossing the face is dotted with. Around
+    every interior corner of a jittered plane the faces close to round-off, so a uniform flow carries as much into each
+    column as out of it; the edge's direction times the face's length does not close where the triangles are not
+    symmetric about their edges, and on this plane a uniform 1 m/s would gather or spread a column's metal at up to
+    539 1/s (measured 2026-10-21; 2.9e-13 with the faces' vectors). On the test sphere the faces close in each column's
+    tangent plane to 3.7 % of the mean face length at worst and 0.03 % on average; what is left is along the normal,
+    the surface's curvature (24 % of a face length on average)."""
+    pts, faces, inner = flat_grid()
+    col = columns.Columns(pts, faces)
+    i, j = col.edges[:, 0], col.edges[:, 1]
+    S = col.face_vector
+    assert np.all(np.einsum("ij,ij->i", S, pts[j] - pts[i]) > 0.0)                    # from the first corner out
+    assert np.all(np.linalg.norm(S, axis=1) <= col.face_length * (1.0 + 1e-12))
+    u = np.array([1.0, 0.3, 0.0])
+    assert np.abs((col.incidence @ (S @ u))[inner] / col.area[inner]).max() < 1e-9
+    e = pts[j] - pts[i]
+    naive = e / np.linalg.norm(e, axis=1)[:, None] * col.face_length[:, None]
+    assert np.abs((col.incidence @ (naive @ u))[inner] / col.area[inner]).max() > 100.0
+    col, _ = sphere_columns(coarse_sphere_mesh)
+    closure = col.incidence @ col.face_vector
+    tangential = closure - np.einsum("ij,ij->i", closure, col.normal)[:, None] * col.normal
+    assert np.linalg.norm(tangential, axis=1).max() < 0.075 * col.face_length.mean()
+    assert np.linalg.norm(tangential, axis=1).mean() < 5e-4 * col.face_length.mean()
```

- [ ] **Step 2: The transport's tests** — append to `tests/test_reentry_model_skin.py`:

```diff
--- a/tests/test_reentry_model_skin.py
+++ b/tests/test_reentry_model_skin.py
@@ -815,3 +815,296 @@
     _, terms = skin.push(col, np.full(col.n, 8e-4), np.zeros(col.n), np.zeros(3), V_HAT, RHO_L, SIGMA)
     size = np.linalg.norm(terms["capillary"], axis=1)
     assert size.max() < 20.0 and size.mean() < 0.5
+
+
+# -- sideways transport (sub-plan 21 Task 6a) ---------------------------------------------------------------------------
+LIQUID, SLURRY, RIGID = 950.0, 902.0, 880.0       # K: all liquid; 30 % solid, flowing; rigid (T_rigid 895.1 K)
+
+
+def strip_grid(nx, ny, h=1e-3):
+    """Columns on a flat regular grid of nx by ny corners, h apart (regular_grid's triangles), normals +z."""
+    from reentry_model import columns
+    x, y = np.meshgrid(np.arange(nx) * h, np.arange(ny) * h, indexing="ij")
+    pts = np.column_stack([x.ravel(), y.ravel(), np.zeros(x.size)])
+    idx = np.arange(nx * ny).reshape(nx, ny)
+    a, b, c, d = idx[:-1, :-1].ravel(), idx[1:, :-1].ravel(), idx[1:, 1:].ravel(), idx[:-1, 1:].ravel()
+    return columns.Columns(pts, np.vstack([np.column_stack([a, b, c]), np.column_stack([a, c, d])])), pts
+
+
+def stacked(col, temps):
+    """Skins on the columns `col`, column i holding a cell of 10 um of solid metal's mass at each temperature of
+    temps[i], top first."""
+    mat = scheil()
+    sk = skin.Skins(mat, col.area, skin.SkinSettings(**CELLS_04))
+    count = np.array([len(t) for t in temps])
+    K = int(count.max())
+    T = np.array([list(t) + [300.0] * (K - len(t)) for t in temps])
+    m = np.where(np.arange(K)[None, :] < count[:, None], mat.rho * col.area[:, None] * 1e-5, 0.0)
+    sk._install(m, m * mat.enthalpy(T), count)
+    return sk
+
+
+def moving_from(sk, columns, donor, u, cells=None):
+    """The velocity u (3,) in the donor column's cells (all, or those listed), every other cell at rest; returns the
+    donor's outflow rate across each edge per unit mass [1/s] (E,) toward each of its neighbours, and the neighbours."""
+    sk.U = np.zeros_like(sk.U)
+    cells = np.arange(sk.count[donor]) if cells is None else np.asarray(cells)
+    sk.U[donor, cells] = u
+    i, j = columns.edges[:, 0], columns.edges[:, 1]
+    flux = np.where(i == donor, columns.face_vector @ u, np.where(j == donor, -(columns.face_vector @ u), 0.0))
+    other = np.where(i == donor, j, i)
+    out = flux > 0.0
+    return np.maximum(flux[out], 0.0) / sk.area[donor], other[out]
+
+
+def test_with_every_velocity_zero_the_transport_changes_nothing(coarse_sphere_mesh):
+    sk, col, _ = random_sphere_skins(coarse_sphere_mesh)
+    before = sk.state()
+    book = skin.transport(sk, col, 0.01)
+    assert all(np.array_equal(x, y) for x, y in zip(before, sk.state()))
+    assert not (book.mass.any() or book.enthalpy.any() or book.momentum.any() or book.wall.any()) and book.courant == 0.0
+
+
+def test_rigid_metal_joined_to_the_substrate_never_leaves_whatever_its_velocity():
+    """Anchored cells hold zero velocity after every momentum step (Task 4); given one anyway they do not move."""
+    col, _ = strip_grid(3, 3)
+    sk = stacked(col, [[RIGID] * 6] * col.n)
+    sk.U = np.where(sk.valid[..., None], np.array([0.3, -0.2, 0.0]), 0.0)
+    before = sk.state()
+    skin.transport(sk, col, 0.01)
+    assert all(np.array_equal(x, y) for x, y in zip(before, sk.state()))
+
+
+def test_a_warm_thick_pulse_carried_by_a_uniform_flow_arrives_where_the_flow_takes_it():
+    """A pulse half again as thick as the layer around it and 50 K warmer (a Gaussian 3 mm wide) in liquid on 1 mm
+    columns, every cell moving at 0.1 m/s along x for 0.1 s, against the same grid without the pulse (whose open ends
+    drain and pile up alike): its mass and heat are kept and their centroids move by u t, 10 mm (to 6e-8, 1e-5 and
+    8.3e-4 mm at Courant numbers nu = u dt / h of 0.5, 1 and 2, measured 2026-10-21, as the scheme's long upwind tail
+    reaches the grid's edges), while the implicit upwind scheme spreads them as the one-dimensional scheme does, by
+    h^2 nu (1 + nu) of variance a sub-step -- the numerical diffusion u h (1 + nu) / 2 -- less 4.2, 3.2 and 2.4 %. The
+    books close to round-off and every velocity stays u."""
+    h, u, x0, w = 1e-3, 0.1, 20e-3, 3e-3
+
+    def run(dt, pulse):
+        col, pts = strip_grid(61, 9, h)
+        g = np.exp(-0.5 * ((pts[:, 0] - x0) / w) ** 2) if pulse else np.zeros(col.n)
+        s = skin.SkinSettings(cell=20e-6)
+        sk = skin.Skins(scheil(), col.area, s)
+        sk.fill(skin.minimum_mass(sk.material, sk.area, s) * (1.0 + 0.5 * g), LIQUID + 50.0 * g)
+        sk.regrid()
+        sk.U = np.where(sk.valid[..., None], np.array([u, 0.0, 0.0]), 0.0)
+        M0, E0 = sk.mass(), sk.energy()
+        start = (sk.column_mass(), sk.H.sum(axis=1))
+        for _ in range(int(round(0.1 / dt))):
+            skin.transport(sk, col, dt)
+        assert sk.mass() == pytest.approx(M0, rel=1e-13) and sk.energy() == pytest.approx(E0, rel=1e-13)
+        assert np.abs(sk.U[sk.valid] - np.array([u, 0.0, 0.0])).max() < 1e-12
+        return start, (sk.column_mass(), sk.H.sum(axis=1)), pts
+
+    def moments(f, x):
+        mean = (x * f).sum() / f.sum()
+        return mean, ((x - mean) ** 2 * f).sum() / f.sum(), f.sum()
+
+    for dt in (0.005, 0.01, 0.02):
+        (m0, H0), (m1, H1), pts = run(dt, True)
+        (r0, R0), (r1, R1), _ = run(dt, False)
+        inner = (pts[:, 1] > 0.0) & (pts[:, 1] < pts[:, 1].max())
+        nu = u * dt / h
+        for before, after in (((m0 - r0), (m1 - r1)), ((H0 - R0), (H1 - R1))):
+            a, b = moments(before[inner], pts[inner, 0]), moments(after[inner], pts[inner, 0])
+            assert b[2] == pytest.approx(a[2], rel=1.5e-5)                          # 6.2e-6 at nu = 2
+            assert b[0] - a[0] == pytest.approx(u * 0.1, abs=2e-6)
+            assert b[1] - a[1] == pytest.approx(h * h * nu * (1.0 + nu) * 0.1 / dt, rel=0.09)
+
+
+def test_depth_is_matched_by_fraction_of_the_flowing_layer():
+    """Spec 6.5: metal leaving the donor between fractions s1 and s2 of its flowing depth joins each receiver between
+    the same fractions of the receiver's, which flows from its top down. The donor's ten liquid cells over rigid metal
+    move only in their first and third (fractions 0-0.1 and 0.2-0.3); a receiver of 20 liquid cells takes them in its
+    cells 0-1 and 4-5, a receiver of 4 in its cells 0 and 1 (0.2-0.3 straddles its 0.25), each in proportion to the
+    overlap; nothing reaches rigid metal."""
+    col, _ = strip_grid(3, 3)
+    temps = [[LIQUID] * 20 + [RIGID] * 5 if k % 2 else [LIQUID] * 4 + [RIGID] * 10 for k in range(col.n)]
+    temps[4] = [LIQUID] * 10 + [RIGID] * 5
+    sk = stacked(col, temps)
+    dt, u = 0.01, np.array([0.08, 0.05, 0.0])
+    rate, receivers = moving_from(sk, col, 4, u, cells=[0, 2])
+    before = sk.m.copy()
+    book = skin.transport(sk, col, dt, regrid=False)
+    m_left = before[4, [0, 2]] / (1.0 + dt * rate.sum())               # the donor's cells, implicitly
+    for r, c in zip(receivers, rate):
+        gain = sk.m[r] - before[r]
+        sent = dt * c * m_left                                        # its first and third cells' share
+        if sk.count[r] == 25:
+            np.testing.assert_allclose(gain[[0, 1, 4, 5]], [sent[0] / 2, sent[0] / 2, sent[1] / 2, sent[1] / 2],
+                                       rtol=1e-12)
+            assert np.all(gain[2:4] == 0.0) and np.all(gain[6:] == 0.0)
+        else:
+            np.testing.assert_allclose(gain[:2], [sent[0] + sent[1] / 2, sent[1] / 2], rtol=1e-12)
+            assert np.all(gain[2:] == 0.0)
+        assert book.mass[r] == pytest.approx(sent.sum(), rel=1e-12)
+    assert len(receivers) >= 2 and {25, 14} <= set(sk.count[receivers])
+
+
+def test_metal_arriving_on_a_rigid_top_lands_on_it_and_the_receiver_s_own_cells_are_untouched():
+    """Spec 6.5's arrival rule: where the receiver's top is rigid -- anchored solid, or a crust over flowing metal --
+    the arriving metal lands on top of its top cell, as a new cell holding exactly what arrived; nothing passes
+    through rigid metal."""
+    col, _ = strip_grid(3, 3)
+    temps = [[RIGID] * 12 if k % 2 else [RIGID] * 2 + [LIQUID] * 6 + [RIGID] * 4 for k in range(col.n)]
+    temps[4] = [LIQUID] * 8 + [RIGID] * 4
+    sk = stacked(col, temps)
+    rate, receivers = moving_from(sk, col, 4, np.array([0.1, 0.04, 0.0]), cells=np.arange(8))
+    before = sk.state()
+    book = skin.transport(sk, col, 0.01, regrid=False)
+    h_donor = sk.material.enthalpy(LIQUID)
+    for r in range(col.n):
+        if r == 4:
+            continue
+        if r in receivers:
+            assert sk.count[r] == before[4][r] + 1
+            assert sk.m[r, 0] == pytest.approx(book.mass[r], rel=1e-12)
+            assert sk.H[r, 0] == pytest.approx(book.mass[r] * h_donor, rel=1e-12)
+            np.testing.assert_allclose(sk.U[r, 0], [0.1, 0.04, 0.0], rtol=1e-12)
+            shifted = slice(1, sk.count[r])
+        else:
+            assert book.mass[r] == 0.0
+            shifted = slice(0, sk.count[r])
+        for now, then in ((sk.m, before[1]), (sk.H, before[2]), (sk.T, before[3]), (sk.U, before[5])):
+            assert np.array_equal(now[r, shifted], then[r, :before[4][r]])
+    assert any(r % 2 for r in receivers) and any(r % 2 == 0 for r in receivers)          # solid tops and crusts
+
+
+def test_a_floating_crust_carried_across_lands_on_top_and_the_slurry_beneath_it_joins_by_fraction():
+    """Asha, 2026-10-20 (Task 4): a crust over flowing metal floats and moves as one piece. Carried into a neighbour it
+    lands on top of the neighbour's top; the slurry beneath it, ten cells, joins each receiver's 20 liquid cells two by
+    two; the donor's anchored cells, given the same velocity, stay where they are."""
+    col, _ = strip_grid(3, 3)
+    temps = [[LIQUID] * 20 + [RIGID] * 5] * col.n
+    temps[4] = [RIGID] * 2 + [SLURRY] * 10 + [RIGID] * 5
+    sk = stacked(col, temps)
+    dt = 0.01
+    rate, receivers = moving_from(sk, col, 4, np.array([0.1, 0.04, 0.0]))       # every cell, the anchored too
+    before = sk.state()
+    book = skin.transport(sk, col, dt, regrid=False)
+    left = before[1][4, :12] / (1.0 + dt * rate.sum())
+    np.testing.assert_allclose(sk.m[4, :12], left, rtol=1e-12)
+    for now, then in ((sk.m, before[1]), (sk.H, before[2]), (sk.T, before[3]), (sk.U, before[5])):
+        assert np.array_equal(now[4, 12:17], then[4, 12:17])                    # the anchored cells
+    h_crust, h_slurry = sk.material.enthalpy(RIGID), sk.material.enthalpy(SLURRY)
+    for r, c in zip(receivers, rate):
+        assert sk.count[r] == 26
+        assert sk.m[r, 0] == pytest.approx(dt * c * left[:2].sum(), rel=1e-12)
+        assert sk.H[r, 0] == pytest.approx(sk.m[r, 0] * h_crust, rel=1e-12)
+        gain = sk.m[r, 1:21] - before[1][r, :20]
+        np.testing.assert_allclose(gain, np.repeat(dt * c * left[2:12] / 2.0, 2), rtol=1e-12)
+        H_gain = sk.H[r, 1:21] - before[2][r, :20]
+        np.testing.assert_allclose(H_gain, gain * h_slurry, rtol=1e-9)
+        assert np.array_equal(sk.m[r, 21:], before[1][r, 20:25])
+        assert book.mass[r] == pytest.approx(dt * c * left.sum(), rel=1e-12)
+
+
+def mixed_sphere_skins(coarse_sphere_mesh, seed=6, speed=0.05):
+    """Skins on the 4 mm test sphere, 1-3 times their minimum deep, every cell at its own temperature between 860 and
+    960 K -- liquid, slurry, crusts and anchored rigid metal in every arrangement -- and moving at its own velocity in
+    its column's tangent plane, about `speed` m/s."""
+    sk, col, rng = random_sphere_skins(coarse_sphere_mesh, seed)
+    sk.T = np.where(sk.valid, rng.uniform(860.0, 960.0, sk.T.shape), sk.T)
+    sk.H = sk.m * sk.material.enthalpy(sk.T)
+    U = rng.normal(0.0, speed, sk.U.shape)
+    U -= np.einsum("ikx,ix->ik", U, col.normal)[..., None] * col.normal[:, None, :]
+    sk.U = np.where(sk.valid[..., None], U, 0.0)
+    return sk, col, rng
+
+
+def test_the_books_close_on_skins_of_every_kind_on_the_sphere(coarse_sphere_mesh):
+    """Mass, enthalpy and the cells' carried sensitivity are conserved to round-off, each column's change is its book;
+    the momentum of each column changes by what was carried in less what was carried out plus the wall's impulse, the
+    normal parts the projection removed, which lie along the column's normal; every velocity ends in its column's
+    plane, and the anchored cells are untouched. Measured 2026-10-21 at Courant numbers up to 0.87: every book to
+    4e-15 of its scale, the wall's impulse 2.5 % of the largest column's momentum."""
+    sk, col, rng = mixed_sphere_skins(coarse_sphere_mesh)
+    sk.S = np.where(sk.valid, rng.uniform(0.0, 1.0, sk.m.shape), 0.0)
+    flowing = sk.valid & np.isfinite(skin.cell_viscosity(sk, 0.0))
+    anchored = skin.anchored_cells(sk, flowing)
+    assert flowing.sum() > 1000 and anchored.sum() > 1000 and (sk.valid & ~flowing & ~anchored).sum() > 100
+    before = sk.state()
+    mS0, P0, M0, E0 = (sk.m * sk.S).sum(), sk.column_momentum(), sk.column_mass(), sk.H.sum(axis=1)
+    book = skin.transport(sk, col, 0.01, regrid=False)
+    assert book.courant > 0.1
+    np.testing.assert_allclose(sk.column_mass() - M0, book.mass, rtol=0, atol=1e-13 * M0.max())
+    np.testing.assert_allclose(sk.H.sum(axis=1) - E0, book.enthalpy, rtol=0, atol=1e-13 * np.abs(E0).max())
+    assert abs(book.mass.sum()) < 1e-13 * M0.sum() and abs(book.enthalpy.sum()) < 1e-13 * np.abs(E0).sum()
+    assert (sk.m * sk.S).sum() == pytest.approx(mS0, rel=1e-13)
+    scale = np.abs(P0).max()
+    np.testing.assert_allclose(sk.column_momentum() - P0, book.momentum + book.wall, rtol=0, atol=1e-13 * scale)
+    assert np.abs(book.momentum.sum(axis=0)).max() < 1e-13 * scale
+    assert np.linalg.norm(np.cross(book.wall, col.normal), axis=1).max() < 1e-13 * scale
+    assert np.abs(book.wall).max() > 1e-6 * scale                                  # the planes differ: it is not idle
+    assert np.abs(np.einsum("ikx,ix->ik", sk.U, col.normal)).max() < 1e-12 * np.abs(sk.U).max()
+    top = (sk.count - before[4])[:, None]                                         # a new top cell shifts the rest
+    k = np.arange(anchored.shape[1])[None, :]
+    rows, cells = np.nonzero(anchored)
+    for now, then in ((sk.m, before[1]), (sk.H, before[2]), (sk.T, before[3])):
+        assert np.array_equal(now[rows, cells + top[rows, 0]], then[rows, cells])
+
+
+def test_courant_numbers_far_above_one_stay_positive_bounded_and_conservative(coarse_sphere_mesh):
+    """A 10 s sub-step on the mixed sphere skins (Courant numbers up to 875): every cell keeps a positive mass, every
+    specific enthalpy stays within the old ones' range and every speed below the old fastest (the scheme mixes, it
+    never overshoots), the books close and the cells end within their split and merge bounds."""
+    sk, col, _ = mixed_sphere_skins(coarse_sphere_mesh)
+    v = sk.valid
+    h0, M0, E0 = sk.H[v] / sk.m[v], sk.mass(), sk.energy()
+    fastest = np.linalg.norm(sk.U, axis=2).max()
+    book = skin.transport(sk, col, 10.0)
+    v = sk.valid
+    assert book.courant > 100.0
+    assert np.all(sk.m[v] > 0.0)
+    h = sk.H[v] / sk.m[v]
+    assert h.min() >= h0.min() - 1e-9 * abs(h0.min()) and h.max() <= h0.max() + 1e-9 * abs(h0.max())
+    assert np.linalg.norm(sk.U, axis=2).max() <= fastest * (1.0 + 1e-12)
+    assert sk.mass() == pytest.approx(M0, rel=1e-13) and sk.energy() == pytest.approx(E0, rel=1e-13)
+    dz, t = sk.dz(), sk.targets()
+    assert np.all(dz[v] >= skin.MERGE_BELOW * t[v]) and np.all(dz[v] <= skin.SPLIT_ABOVE * t[v])
+
+
+def test_a_velocity_carried_into_a_column_is_projected_onto_its_plane(coarse_sphere_mesh):
+    """One column on the 4 mm test sphere, liquid throughout, moving at 1 m/s in its plane into its neighbours, which are
+    liquid at rest: each neighbour gains the projection onto its own plane of the momentum that left toward it, and
+    the wall takes the normal part -- the planes of neighbouring 4 mm columns on a 50 mm sphere are 4.6 degrees apart
+    on average, the three receivers' 4.5-4.9 degrees from the donor's, and the normal parts 2.3-7.7 % of the momentum
+    carried (measured 2026-10-21)."""
+    from test_reentry_model_columns import sphere_columns
+    col, _ = sphere_columns(coarse_sphere_mesh)
+    sk = skin.Skins(scheil(), col.area, skin.SkinSettings())
+    sk.fill(skin.minimum_mass(sk.material, sk.area, sk.settings), LIQUID)
+    sk.regrid()
+    donor = int(np.argmin(np.abs(col.normal[:, 2])))                              # on the equator of z
+    u = np.cross(col.normal[donor], [0.0, 0.0, 1.0])
+    rate, receivers = moving_from(sk, col, donor, u / np.linalg.norm(u))
+    P0 = sk.column_momentum()
+    book = skin.transport(sk, col, 0.01, regrid=False)
+    n = col.normal[receivers]
+    carried = book.momentum[receivers]
+    inplane = carried - np.einsum("ij,ij->i", carried, n)[:, None] * n
+    np.testing.assert_allclose(sk.column_momentum()[receivers] - P0[receivers], inplane, rtol=0,
+                               atol=1e-12 * np.abs(carried).max())
+    np.testing.assert_allclose(book.wall[receivers], inplane - carried, rtol=0, atol=1e-12 * np.abs(carried).max())
+    normal_part = np.linalg.norm(carried - inplane, axis=1) / np.linalg.norm(carried, axis=1)
+    assert normal_part.max() > 0.02 and normal_part.max() < 0.16
+
+
+def test_the_upwind_order_changes_the_cost_never_the_answer(coarse_sphere_mesh, monkeypatch):
+    """The solve takes the cells in upwind order -- scipy's strong components, which come out in reverse topological
+    order -- so that its sparse LU has almost no fill; were that order ever not upwind, the LU falls back to its own
+    fill-reducing order and gives the same answer, more slowly."""
+    first, _, _ = mixed_sphere_skins(coarse_sphere_mesh)
+    second, col, _ = mixed_sphere_skins(coarse_sphere_mesh)
+    a = skin.transport(first, col, 0.05, regrid=False)
+    rng = np.random.default_rng(0)
+    monkeypatch.setattr(skin, "_upwind_order", lambda src, dst, n: (rng.permutation(n), False))
+    b = skin.transport(second, col, 0.05, regrid=False)
+    for x, y in zip(first.state(), second.state()):
+        np.testing.assert_allclose(x, y, rtol=1e-11, atol=1e-11 * np.abs(y).max())
+    np.testing.assert_allclose(a.wall, b.wall, rtol=0, atol=1e-11 * np.abs(a.momentum).max())
```

- [ ] **Step 3: Run them to see them fail.** `"$PY" -m pytest tests/test_reentry_model_columns.py
  tests/test_reentry_model_skin.py -q`. Expected: 11 failed — `AttributeError: 'Columns' object has no attribute
  'face_vector'` (5) and `AttributeError: module 'reentry_model.skin' has no attribute 'transport'` (6).

- [ ] **Step 4: The faces' vectors** in `reentry_model/columns.py`:

```diff
--- a/reentry_model/columns.py
+++ b/reentry_model/columns.py
@@ -23,7 +23,14 @@
 and 23 %). On these 50 mm-radius spheres the curvature is then within 8.4·10^-5 of 2/R at every corner (1.3·10^-5 on
 the 2 mm sphere), where the barycentric area had scattered 7.6 and 3.9 %, and the gradient of a linear field within
 0.09 % of its tangential part on average (0.02 % on the 2 mm sphere), 0.8 % at the test sphere's one irregular corner
-(measured 2026-10-18)."""
+(measured 2026-10-18).
+
+The flow's sideways transport (sub-plan 21 Task 6a) crosses each dual face at its `face_vector`: the face's two
+segments, each its normal in its own triangle times its length, summed and pointing from the edge's first corner to
+its second, so that a velocity u carries u . face_vector per unit depth across the face. A column's faces close
+exactly on a plane, so a uniform flow carries as much into every column as out of it; the edge's direction times the
+face's length closes only where the triangles are symmetric about their edges, and on the tests' jittered plane of 1 mm
+squares a uniform 1 m/s would gather or spread a column's metal at up to 539 1/s that way (measured 2026-10-21)."""
 from dataclasses import dataclass
 
 import numpy as np
@@ -69,15 +76,21 @@
         mid = 0.5 * (p[i] + p[j])
         cot = np.zeros(len(self.edges))
         face = np.zeros(len(self.edges))
+        vector = np.zeros((len(self.edges), 3))
+        _, unit, _ = triangle_geometry(p, self.faces)
         for s in range(2):
             ok = self._tri[:, s] >= 0
             k = self._opp[ok, s]
             u, v = p[i[ok]] - p[k], p[j[ok]] - p[k]
             cot[ok] += np.einsum("ij,ij->i", u, v) / np.linalg.norm(np.cross(u, v), axis=1)
-            face[ok] += np.linalg.norm(mid[ok] - self.centroid[self._tri[ok, s]], axis=1)
+            segment = self.centroid[self._tri[ok, s]] - mid[ok]
+            face[ok] += np.linalg.norm(segment, axis=1)
+            across = np.cross(segment, unit[self._tri[ok, s]])             # the segment's normal times its length
+            vector[ok] += across * np.sign(np.einsum("ij,ij->i", across, p[j[ok]] - p[i[ok]]))[:, None]
         self.cot_weight = 0.5 * cot
         self.conduct_weight = np.maximum(self.cot_weight, 0.0)
         self.face_length = face
+        self.face_vector = vector
 
     def gradient(self, f, points=None):
         """The tangential gradient at every corner of a per-corner field f (n,) [f per m] on the surface through
```

- [ ] **Step 5: The transport** in `reentry_model/skin.py`:

```diff
--- a/reentry_model/skin.py
+++ b/reentry_model/skin.py
@@ -31,6 +31,9 @@
 from dataclasses import dataclass
 
 import numpy as np
+import scipy.sparse as sp
+import scipy.sparse.linalg as spla
+from scipy.sparse.csgraph import connected_components
 
 from . import rheology
 from .thermal import SIGMA_SB
@@ -974,6 +977,227 @@
     return sum(terms.values()), terms
 
 
+# -- sideways transport (sub-plan 21 Task 6a) ---------------------------------------------------------------------------
+SLIVER = 1.0e-12                  # a piece of a donor cell shorter than this share of it is the round-off of a boundary
+                                  # the donor's and the receiver's partitions share, and is left out
+
+
+@dataclass
+class TransportBook:
+    """What a transport sub-step carried between the columns, per column: each column's gain from its neighbours less
+    what it gave them. The mass, the enthalpy and the momentum as it left the donors each sum to zero over the body to
+    round-off; every arrival's velocity is projected onto its receiver's plane, and the normal parts removed are the
+    wall's impulse, so that each column's momentum changes by momentum + wall."""
+    mass: np.ndarray         # kg (n,)
+    enthalpy: np.ndarray     # J (n,)
+    momentum: np.ndarray     # kg m/s (n, 3), as it left the donors, before the projection
+    wall: np.ndarray         # N s (n, 3), the normal parts removed, each along its own column's normal
+    courant: float           # the largest outflow rate of any cell [1/s] times the sub-step
+
+
+def tangent_frames(normals):
+    """Two orthonormal vectors (e1, e2), each (n, 3), spanning every column's tangent plane, e1 x e2 along the normal:
+    e1 the global x axis's part in the plane, or y's where the normal lies within 26 degrees of x."""
+    ref = np.where((np.abs(normals[:, 0]) < 0.9)[:, None], np.array([1.0, 0.0, 0.0]), np.array([0.0, 1.0, 0.0]))
+    e1 = ref - np.einsum("ij,ij->i", ref, normals)[:, None] * normals
+    e1 /= np.linalg.norm(e1, axis=1)[:, None]
+    return e1, np.cross(normals, e1)
+
+
+def _fractions(depth):
+    """Per column, its cells' boundaries (n, K+1) as fractions of the depth they hold (`depth` (n, K) each cell's share,
+    zero for a cell that does not count), from 0 at the top to 1 at the last that counts and held there; 1 throughout
+    where none counts."""
+    z = np.concatenate([np.zeros((depth.shape[0], 1)), np.cumsum(depth, axis=1)], axis=1)
+    total = z[:, -1:]
+    return np.where(total > 0.0, z / np.where(total > 0.0, total, 1.0), 1.0)
+
+
+def _pieces(sig, tau):
+    """The pieces two partitions of [0, 1] cut each other into, row by row: sig and tau (D, K+1) non-decreasing
+    boundaries. Returns, for every piece of positive length, its row, its cell of sig, its cell of tau and its length.
+    A rank merge: each boundary's place in the merged row is its own index plus the number of the other partition's
+    boundaries before it (tau's first where two are equal)."""
+    D, B = sig.shape
+    rows = np.arange(D)[:, None]
+    off = (np.arange(D) * 4.0)[:, None]                      # boundaries lie in [0, 1]: rows 4 apart never mix
+    before = lambda xp, x, side: (np.searchsorted((xp + off).ravel(), (x + off).ravel(), side=side).reshape(x.shape)
+                                  - rows * B)
+    at_s = np.arange(B)[None, :] + before(tau, sig, "right")
+    at_t = np.arange(B)[None, :] + before(sig, tau, "left")
+    merged = np.empty((D, 2 * B))
+    merged[rows, at_s], merged[rows, at_t] = sig, tau
+    from_sig = np.zeros((D, 2 * B), dtype=np.int64)
+    from_sig[rows, at_s] = 1
+    seen = np.cumsum(from_sig, axis=1)[:, :-1]                # sig's boundaries at or before each piece's start
+    length = np.diff(merged, axis=1)
+    r, p = np.nonzero(length > 0.0)
+    return r, seen[r, p] - 1, p - seen[r, p], length[r, p]
+
+
+def _upwind_order(src, dst, n):
+    """An order of n unknowns in which every donor (src) comes before its receivers (dst), except within a cycle of
+    transfers, which stays together; and whether it holds. The order is scipy's strong components, which its (Pearce's)
+    algorithm numbers in reverse topological order -- how scipy numbers them, not what it promises, so it is checked."""
+    graph = sp.csr_matrix((np.ones(len(src)), (src, dst)), shape=(n, n))
+    _, label = connected_components(graph, directed=True, connection="strong")
+    return np.argsort(-label, kind="stable"), bool(np.all(label[src] >= label[dst]))
+
+
+def _factor(A, upwind):
+    """A sparse LU of A, an M-matrix or one by blocks, in the given order with no pivoting where the order is upwind (the
+    matrix is then triangular but for its cycles: no fill, every sign kept), else in SuperLU's own fill-reducing order."""
+    if upwind:
+        return spla.splu(A, permc_spec="NATURAL", diag_pivot_thresh=0.0, options=dict(SymmetricMode=True))
+    return spla.splu(A)
+
+
+def transport(sk, columns, dt, regrid=True):
+    """One linearly implicit upwind sub-step dt [s] of the metal moving sideways between neighbouring columns at the
+    cells' velocities sk.U (spec 6.5; Asha, 2026-10-21: its own step, after the momentum step). Updates the cells and
+    returns a TransportBook.
+
+    Across each face, cell k of the donor column a sends c m_k' toward the receiver b, c = (u_k . f)+ / A_a [1/s], f the
+    face's vector (`Columns.face_vector`) and m_k' the cell's mass at the sub-step's end; it carries the cell's specific
+    enthalpy and velocity. Every cell is then m_k' (1 + dt sum c) = m_k + dt (what arrives): (I + dt C) m' = m, the
+    film's runoff scheme (`film.Runoff`) on the cells -- an M-matrix whose columns sum to one, so the masses stay
+    positive, the books close to round-off and nothing limits dt. The enthalpy, the carried sensitivity S and the
+    momentum solve the same system, the momentum with every arrival projected onto its receiver's tangent plane, in two
+    components per cell (`tangent_frames`), the normal part booked to the wall (spec 6.4, section 14 assumption 2).
+
+    Where metal joins the receiver (spec 6.5, section 14 assumption 3): a flowing donor cell's span of fractions of the
+    donor's flowing depth (its flowing cells, top down) joins the receiver's cells between the same fractions of its
+    flowing layer, the run of flowing cells from its top, each in proportion to the overlap. Where the receiver's top is
+    rigid -- anchored metal, or a crust -- everything lands on top of it, as a new top cell; so does a floating crust's
+    metal (Task 4, Asha 2026-10-20) wherever it goes, being no part of any flowing depth. Rigid metal joined to the
+    substrate (`anchored_cells`) neither leaves nor receives, whatever velocity it holds: nothing passes through rigid
+    metal. Only the cells that send or receive are solved for; with every velocity zero nothing changes.
+
+    The cells are taken in upwind order (`_upwind_order`), in which the system is triangular but for its cycles where
+    flows meet, so its sparse LU has almost no fill. With `regrid` the cells are then split and merged (spec 6.2); a
+    sub-step that regrids anyway passes False.
+
+    Measured 2026-10-21 on the 100 mm sphere's skins (9 417 columns, 10 um cells), flowing along the meridian away from
+    the stagnation point: 0.27 s a sub-step with the front hemisphere's top 0.3 mm flowing (124 000 of 553 000 cells
+    moving), the same at Courant numbers of 1 and 96, and 0.95 s with all 584 000 cells moving, against 0.6-0.8 s for a
+    coupled conduction sub-step on the sphere's skins at 700 K. With every cell moving, the mass's LU in the cells' own
+    order took 19.6 s and in SuperLU's fill-reducing order 0.66 s, against 0.06 s upwind. Projected once at the
+    sub-step's end instead of at every arrival, a
+    uniform meridian flow kept 90 % of its kinetic energy in one sub-step at Courant number 9.6 against 96 %, and 24 %
+    against 49 % at 96: metal crossing many columns in a sub-step would carry its velocity straight on, off the curved
+    surface."""
+    n, K = sk.m.shape
+    valid = sk.valid
+    flowing = valid & np.isfinite(cell_viscosity(sk, 0.0))
+    moving = valid & ~anchored_cells(sk, flowing) & (sk.U != 0.0).any(axis=2)
+    book = TransportBook(np.zeros(n), np.zeros(n), np.zeros((n, 3)), np.zeros((n, 3)), 0.0)
+    if not moving.any():
+        return book
+    # every face a moving column has, each way, and every cell's outflow rate across it [1/s]
+    i, j = columns.edges[:, 0], columns.edges[:, 1]
+    a, b, f = np.r_[i, j], np.r_[j, i], np.r_[columns.face_vector, -columns.face_vector]
+    keep = np.flatnonzero(moving.any(axis=1)[a])
+    a, b, f = a[keep], b[keep], f[keep]
+    c = np.where(moving[a], np.maximum(np.einsum("dkx,dx->dk", sk.U[a], f), 0.0), 0.0) / sk.area[a][:, None]
+    keep = np.flatnonzero((c > 0.0).any(axis=1))
+    if not len(keep):
+        return book
+    a, b, c = a[keep], b[keep], c[keep]
+    D = len(a)
+    # where it joins: by fraction of the flowing depth, or on top
+    dz = sk.dz()
+    run = np.cumprod(flowing, axis=1).sum(axis=1)             # the receiver's flowing layer: its top run of flowing cells
+    sig = _fractions(np.where(flowing, dz, 0.0))
+    tau = _fractions(np.where(np.arange(K)[None, :] < run[:, None], dz, 0.0))
+    joins = np.where(flowing[a] & (run[b] > 0)[:, None], c, 0.0)
+    rows = np.flatnonzero((joins > 0.0).any(axis=1))
+    r, kd, kr, length = _pieces(sig[a[rows]], tau[b[rows]])
+    some = joins[rows[r], kd] > 0.0
+    edge, kd, kr, length = rows[r[some]], kd[some], kr[some], length[some]
+    group = edge * K + kd
+    span = np.bincount(group, length, D * K)[group]
+    some = length > SLIVER * span                             # not the round-off of a boundary the two share
+    edge, kd, kr, length, group = edge[some], kd[some], kr[some], length[some], group[some]
+    share = length / np.bincount(group, length, D * K)[group]  # each donor cell's pieces sum to one
+    te, tk = np.nonzero((c > 0.0) & ~(flowing[a] & (run[b] > 0)[:, None]))
+    # every arrival: its face, the receiving cell (n K + b: a new top cell of column b), the sending cell, dt c share
+    e_edge = np.r_[edge, te]
+    e_row = np.r_[b[edge] * K + kr, n * K + b[te]]
+    e_col = np.r_[a[edge] * K + kd, a[te] * K + tk]
+    e_val = dt * np.r_[c[edge, kd] * share, c[te, tk]]
+    # the cells that send or receive, in upwind order; each sends what its arrivals take, so the books close exactly
+    active = np.zeros(n * K + n, dtype=bool)
+    active[e_row], active[e_col] = True, True
+    act = np.flatnonzero(active)
+    N = len(act)
+    local = np.full(n * K + n, -1, dtype=np.int64)
+    local[act] = np.arange(N)
+    order, upwind = _upwind_order(local[e_col], local[e_row], N)
+    at = np.empty(N, dtype=np.int64)
+    at[order] = np.arange(N)
+    to, frm = at[local[e_row]], at[local[e_col]]
+    diag = 1.0 + np.bincount(local[e_col], e_val, N)
+    # mass, enthalpy and the carried sensitivity
+    lu = _factor(sp.csc_matrix((np.r_[diag, -e_val], (np.r_[at, to], np.r_[at, frm])), shape=(N, N)), upwind)
+    cell = lambda v: np.r_[v.reshape((n * K,) + v.shape[2:]), np.zeros((n,) + v.shape[2:])][act]
+    rhs = [cell(sk.m), cell(sk.H)] + ([cell(sk.m * sk.S)] if sk.S is not None else [])
+    x = np.empty((N, len(rhs)))
+    x[order] = lu.solve(np.column_stack(rhs)[order])
+    # the momentum, each arrival projected onto its receiver's plane: two components per cell in its own plane's frame
+    owner = np.where(act < n * K, act // K, act - n * K)
+    e1, e2 = tangent_frames(columns.normal)
+    P = cell(sk.m[..., None] * sk.U)
+    p = np.column_stack([np.einsum("ij,ij->i", P, e1[owner]), np.einsum("ij,ij->i", P, e2[owner])])
+    held = P - p[:, :1] * e1[owner] - p[:, 1:] * e2[owner]   # any normal part a cell held: the wall's
+    R = np.stack([np.stack([np.einsum("ij,ij->i", e[b], g[a]) for g in (e1, e2)], axis=1) for e in (e1, e2)], axis=1)
+    Re = R[e_edge]                                            # [arrival, receiver's component, donor's component]
+    rows2 = [2 * at, 2 * at + 1] + [2 * to + q for q in (0, 0, 1, 1)]
+    cols2 = [2 * at, 2 * at + 1] + [2 * frm + q for q in (0, 1, 0, 1)]
+    vals2 = [diag, diag] + [-e_val * Re[:, q, s] for q, s in ((0, 0), (0, 1), (1, 0), (1, 1))]
+    lu2 = _factor(sp.csc_matrix((np.concatenate(vals2), (np.concatenate(rows2), np.concatenate(cols2))),
+                                shape=(2 * N, 2 * N)), upwind)
+    rhs2 = np.empty(2 * N)
+    rhs2[2 * at], rhs2[2 * at + 1] = p[:, 0], p[:, 1]
+    y = lu2.solve(rhs2)
+    p = np.column_stack([y[2 * at], y[2 * at + 1]])
+    P = p[:, :1] * e1[owner] + p[:, 1:] * e2[owner]
+    # the books, face by face: what every arrival took from its sending cell at the sub-step's end
+    moved = e_val[:, None] * np.column_stack([x[:, :2], p])[local[e_col]]
+    sent = np.stack([np.bincount(e_edge, moved[:, v], D) for v in range(4)], axis=1)
+    carried = sent[:, 2:3] * e1[a] + sent[:, 3:4] * e2[a]
+    lost = np.einsum("ij,ij->i", carried, columns.normal[b])[:, None] * columns.normal[b]
+    net = lambda v: np.bincount(b, v, n) - np.bincount(a, v, n)
+    book = TransportBook(net(sent[:, 0]), net(sent[:, 1]), np.stack([net(carried[:, v]) for v in range(3)], axis=1),
+                         -np.stack([np.bincount(b, lost[:, v], n) + np.bincount(owner, held[:, v], n)
+                                    for v in range(3)], axis=1), float((diag - 1.0).max()))
+    # install: new values for the cells solved for, a new top cell where metal landed on top
+    m, H, T, U = sk.m.copy(), sk.H.copy(), sk.T.copy(), sk.U.copy()
+    S = sk.S.copy() if sk.S is not None else None
+    old = act < n * K
+    m.ravel()[act[old]], H.ravel()[act[old]], T.ravel()[act[old]] = x[old, 0], x[old, 1], np.nan
+    U.reshape(-1, 3)[act[old]] = P[old] / x[old, :1]
+    if S is not None:
+        S.ravel()[act[old]] = x[old, 2] / x[old, 0]
+    count = sk.count.copy()
+    tops = act[~old] - n * K
+    if len(tops):
+        new = np.zeros(n, dtype=bool)
+        new[tops] = True
+        src = np.clip(np.arange(K + 1)[None, :] - new[:, None], 0, K)
+        down = lambda v: np.take_along_axis(np.concatenate([v, np.zeros((n, 1) + v.shape[2:])], axis=1),
+                                            src.reshape(src.shape + (1,) * (v.ndim - 2)), axis=1)
+        m, H, T, U = down(m), down(H), down(T), down(U)
+        m[tops, 0], H[tops, 0], T[tops, 0], U[tops, 0] = x[~old, 0], x[~old, 1], np.nan, P[~old] / x[~old, :1]
+        if S is not None:
+            S = down(S)
+            S[tops, 0] = x[~old, 2] / x[~old, 0]
+        count = count + new
+    sk._install(m, H, count, S, T, U)
+    if regrid:
+        sk.regrid()
+    return book
+
+
 # -- the base (Task 7) -------------------------------------------------------------------------------------------------
 def hand_down(sk, m_min, T_rigid):
     """Rigid cells at the bottom of each column (at or below T_rigid, contiguous from the bottom) pass to the 3D model,
```

- [ ] **Step 6: Run them to see them pass.** Expected: `tests/test_reentry_model_columns.py` 10 passed and
  `tests/test_reentry_model_skin.py` 56 passed; the unit tier 450 passed, 1 skipped. Nothing calls `transport` yet: the flow
  joins the skin body's sub-steps in 6c, so no result changes.

- [ ] **Step 7: Commit.**

```bash
git add reentry_model/columns.py tests/test_reentry_model_columns.py reentry_model/skin.py tests/test_reentry_model_skin.py
git commit -m "Skins: the metal moves sideways between neighbouring columns at the cells' velocities, linearly implicit and upwind, by fraction of the flowing depth and on top of a rigid top, each arrival's velocity projected onto its receiver's plane (sub-plan 21 Task 6a)"
```

---

### Task 6b: The thickness terms implicit in the column heights

**Files:**
- Modify: `reentry_model/columns.py` (`Columns.film_operator`, `Columns.film_curvature`), `reentry_model/skin.py`
  (`MomentumBook.response` and `.wall_response`; `momentum_step(..., response=False)`, each floating crust one unknown,
  the unused `P0` removed; `push(..., wet=None)` with the capillary term in its thin-film form; `transport(...,
  faces=None)`; `THICKNESS_PASSES`, `THICKNESS_TOL`, `flowing_tops`, `_walled`, `_capillary`, `FacePush`,
  `ThicknessBook`, `thickness_step`)
- Test: `tests/test_reentry_model_columns.py`, `tests/test_reentry_model_skin.py` (and Task 5's two capillary tests,
  re-measured)

**Interfaces:**
- Consumes: Task 2's `Columns.gradient`, `.curvature` (on the base alone now), `.top_points`, `.incidence`,
  `.conduct_weight`; Task 4's `momentum_step`, `anchored_cells`, `cell_viscosity`; Task 5's `push`; 6a's `transport` and
  `Columns.face_vector`.
- Produces: `Columns.film_operator(thickness) -> (n, n)` sparse [1/m2], minus the cotangent Laplacian of the top surface
  over each column's area; `Columns.film_curvature(thickness) -> (n,)` [1/m], the base's curvature plus that operator
  times the thicknesses; `skin.flowing_tops(sk) -> (n,)` bool, where a column's top cell flows; `push(columns,
  thickness, dp_ds, g, v_hat, rho_l, sigma, wet=None)`, its capillary term −σ ∇_s κ with κ the film's curvature, zero
  where `wet` is False; `momentum_step(sk, dt, push, shear, A_top, response=False)`, with `response` its book carrying
  `response` (n, K) [m/s per N/m3], each cell's velocity per unit push, and `wall_response` (n,) [m3], the anchored
  metal's drag on it; `transport(sk, columns, dt, regrid=True, faces=None)`, `faces` a `FacePush`;
  `skin.thickness_step(sk, columns, dt, book, applied, g, rho_l, sigma) -> ThicknessBook(push, wall, heights, faces,
  growth, passes)`, push and wall (n, 3) [N s], heights (n,) [m] the thicknesses predicted at the end of the transport,
  faces a `FacePush(response (n, K), across (2E,) [Pa])` for the transport, growth the largest dt·s of the
  Rayleigh–Taylor growth over the columns, passes the linearisations made; `skin.THICKNESS_PASSES` (8),
  `skin.THICKNESS_TOL` (10⁻³).
- **What 6c calls, in this order, every flow sub-step:**
  1. `wet = flowing_tops(sk)`; `G, terms = push(columns, sk.thickness(), dp_ds, g, v_hat, rho_l, sigma, wet=wet)`;
  2. `mbook = momentum_step(sk, dt, G, shear, A_top, response=True)`;
  3. `tbook = thickness_step(sk, columns, dt, mbook, terms['thickness'] + terms['capillary'], g, rho_l, sigma)`;
  4. `trbook = transport(sk, columns, dt, regrid=False, faces=tbook.faces)`;
  5. the sub-step's own regrid.
  The column's momentum changes by `mbook.push + mbook.shear + mbook.wall + tbook.push + tbook.wall + trbook.momentum +
  trbook.wall`. `tbook.growth` must stay well under one: backward Euler follows a Rayleigh–Taylor growth only where the
  sub-step resolves it.

**Decided before this task was written (Asha, 2026-10-21; spec §6.4, §14 assumption 1, §16).** A first attempt stopped
on a finding: Task 5's capillary term made surface tension run backwards. It took κ as Task 2's curvature of the top
surface, and Task 2 gives a corner touching an obtuse triangle its neighbours' mean curvature. Such a corner cannot see
its own height: raising it makes its neighbours concave and lowers its own curvature, so a pile-up lowers its pressure
and draws metal in. Raising each of the test sphere's 33 such corners by 1 µm lowered its own curvature at 32 (by up to
0.083 1/m) and left it unchanged at the other, where every other corner's rose by at least 0.21 1/m.

The faces' exchange times that term's linearisation has eigenvalues down to −1.48·10¹⁰ 1/m⁴ on the test sphere and
−2.40·10¹¹ 1/m⁴ on the 2 mm sphere (48 such corners). That is a growth of about 80 and 1 300 1/s for 0.3 mm of liquid,
unstable at any sub-step and in any time scheme. A uniform 0.3 mm liquid layer with 1 µm of noise on the test sphere,
its velocities reset after each 10 ms sub-step, roughened from 0.48 to 4.0 µm in 30 sub-steps.

Asha's decisions:
- **(option A)** Surface tension's pressure takes the thin-film form on the top surface's geometry. The pressure is
  (1/A_i) dE/dH_i, where E is the thin film's surface energy on the top surface and A_i is the column's area, the one
  that turns the column's volume into its height in the skins.
- Surface tension acts only where a column's top cell flows.
- The depth-uniform pressure stays, stated as a limitation.

Option B, the top's full curvature over each column's own area with no obtuse rule, was also stable but left a uniform
layer rippled by the mesh: 7.7 µm on the test sphere. Keeping Task 5's term, explicitly or with the energy-consistent
increment, stayed unstable.

**The thin film's curvature.** κ_i = κ_base,i + (J H)_i, where:
- κ_base is Task 2's curvature of the base, its obtuse-corner rule acting on the fixed base alone;
- J = −L_top / A is minus the cotangent Laplacian of the top surface over each column's area (`Columns.film_operator`).

A·J is symmetric and non-negative, and so is the faces' exchange B diag(w) Bᵀ. Their product therefore has no growing
mode: on the test sphere with a rough top (0–3.4 mm) every eigenvalue is real and non-negative to round-off. A cosine's
curvature on a plane is the compact second difference of the heights (to 4·10⁻⁷, from its own top's cotangents). A
uniform layer's on the sphere is the base's exactly, so a uniform 0.8 mm skin on the test sphere now gets at most
0.32 Pa/m of capillary push (0.031 on average) from the mesh, where the top's own curvature gave 11 Pa/m (0.28).

A cotangent depends on angles alone, so the top's cotangents differ from the base's only where the thickness varies:
the two operators on a smooth bump at 3.4 mm differ by 0.35 % on average on the test sphere. With the base's area, which
the energy argument requires, a thick layer's bump curvature stays about 2H/R above the top surface's own (14 % at
3.4 mm on the 50 mm sphere): the columns hold their volume over the base's area. Asha kept the prism columns and
stated the error (spec §14 assumption 10, 2026-10-22), over wedge columns whose area grows with height.

**Where surface tension acts.** Surface tension acts only where a column's top cell flows (`flowing_tops`). An anchored
rigid top and a floating crust feel no capillary push, and no capillary push crosses a face with a non-flowing top on
either side. A rigid neighbour's height enters a flowing column's curvature through the top's Laplacian, as a wall's
would. In the column push's gradient, a wall takes its flowing neighbours' mean curvature (`_walled`), so it pushes
nothing of its own. The energy argument holds with the wall's heights as data: the capillary flux runs only between
flowing tops.

**The scheme** (spec §6.4's last paragraph: "inside the transport solve", read since 2026-10-21 as within the sub-step's
transport step, so the thickness terms see the heights at the end of the transport). The momentum step is linear in
the push for fixed viscosities. With `response` it also solves its last pass's system for a unit push alone: one more
tridiagonal solve, giving each cell's response r and the anchored metal's drag on it.

The two terms are the gradient of P = σκ − ρ_l (g·n) H. Across each dual face P pushes by its two-point difference,
φ = w (P_a − P_b), where:
- w is the face's clipped cotangent weight (Task 1's);
- the capillary part acts only between flowing tops;
- the deceleration's part uses the mean of the two columns' g·n.

φ is linear in the heights through J. Each sub-step then does three things.

1. **It predicts the heights at the end of the transport, H + dH**, from the transport linearised on the columns. A
   column sends across a face what its cells carry, at their velocity without the thickness push plus their response
   times the face's push. This is implicit in its moving metal, as the transport is, and linearised in the face's push
   through the mobility Σ r dz of the cells that cross. Newton's method on that piecewise linear flux re-solves at the
   new face pushes until the flux the last solve missed would move no column by more than THICKNESS_TOL of the largest
   change. It took 1–3 solves on the test sphere and 3 on the melting state.
2. **It corrects every moving cell's velocity** by its response times the change of `push`'s two terms at the predicted
   heights. The cells' momentum follows the push of Task 5's form, and the impulse and the wall's share are booked.
3. **It hands `transport` each face's own push** (`FacePush`). A cell's velocity across a face becomes its velocity
   without the thickness push plus its response times that face's push.

Positivity is 6a's: every rate stays clipped at zero, so the masses stay positive at any velocity.

**A floating crust is one unknown.** In `momentum_step` each run of crust cells is now one unknown of the tridiagonal
system, a rigid body: its mass is the cells', and its drags are those at its top and bottom faces.
- Task 4 coupled the crust's cells by 10⁹ Pa s. Those rows, some 10¹² times the cells' inertia, amplified the solver's
  round-off: on mixed sphere skins the books closed to 1.2·10⁻⁵ of the largest impulse. They now close to
  6.6·10⁻¹⁴ (1 ms) and 6.2·10⁻¹³ (10 ms).
- The velocities moved by up to 1.1·10⁻³ of the largest at 10 ms. This is the coupled crust's own round-off: it moved
  by 1.1·10⁻⁴, 1.1·10⁻⁵ and 9.6·10⁻⁷ at 10⁸, 10⁷ and 10⁶ Pa s, and by 4.2·10⁻⁷ and 4.2·10⁻⁶ at 10⁵ and 10⁴, where the
  crust's flexibility shows.
- On Task 4's own crust column the change is 9.8·10⁻⁸ of the top speed, and its books close to 8.9·10⁻¹⁵ instead of
  9.9·10⁻⁸.
- Columns without a crust are unchanged to the last bit.

**Measured** (2026-10-21, `prototype/work-2026-10-18-column-flow/measure/task6b/`; the tests' numbers).

*Rayleigh–Taylor growth* (spec §11.1 level 1). The layer is 1 mm of slurry at 902 K (30 % solid, 2.68 Pa s at rest)
on a strip of 1 mm columns two corners wide, pulled out of its surface at 100 m/s²: a capillary cutoff of 11.5 mm. This
is the lubrication regime the formula assumes, the viscous time ρh²/μ being 0.9 ms against a growth time of 0.63 s.
- A cosine 24 mm long grows at the thin-film rate h³k²(ρg − σk²)/(3μ) = 1.578 1/s to within −0.41, −0.22 and +0.42 % at
  sub-steps of 1, 2.5 and 10 ms. The grid's compact difference accounts for −0.40 %, and backward Euler adds +dt·s/2.
- A cosine 8 mm long, below the cutoff, decays at −19.44 1/s to within 1.2 % at 1 ms. At 10 ms it decays 8.7 % slower,
  backward Euler's own damping, whose rate it meets within 0.10 % (0.26 % at 1 ms).

*The two-column mode.* 0.3 mm of liquid on 2 mm columns, alternating by 1 µm, falls to the following fraction of itself
in one sub-step, and never above 0.009 of itself after that:

| Sub-step | Under 30 m/s² | Without the pull |
|---|---|---|
| 10 ms | 0.112 | 0.117 |
| 100 ms | 5.5·10⁻⁴ | −6.3·10⁻⁵ |
| 1 s | −0.051 | −0.056 |

The table shows the grid-scale mode damped within one sub-step at every sub-step from 10 ms to 1 s, with or without the
pull.

*The faces' own push against the columns' alone.* A bump 1 µm high on one column of 1.1 mm of slurry relaxes to
12.5 % in one 10 ms sub-step and 3.0 % in five. With the columns' pushes alone it stays at 99.8 %: they are centred
differences, zero on the bump's column by symmetry, and the upwind transport takes metal out of a column only at its
own velocity.

*A uniform layer.* 0.3 mm of liquid on the test sphere stays uniform to 26 nm, with speeds of at most 2.2 µm/s, over 20
sub-steps of 10 ms.

*Prediction and Newton passes.* With 1 µm of roughness on 0.3 mm of liquid on the test sphere:
- the transport makes the predicted heights to 0.62 % of the largest change;
- the Newton passes were 1, 3, 3 and 3;
- the roughness falls to 34 % in four sub-steps.

*Books.* On mixed sphere skins, every momentum book closes to 7.8·10⁻¹⁴ (1 ms) and 8.8·10⁻¹³ (10 ms) of the largest
impulse, and mass and enthalpy are exact. With σ = 0 and g = 0 the whole step equals the momentum step followed by the
transport, to the last bit.

*Cost* on the 100 mm sphere's melting state. The melted face is the one where g·n > 0 under the flight direction given
to `push`: its top 0.3 mm is liquid, with 141 000 cells moving at 0.1 m/s along the meridian, under 30 m/s² and
σ = 0.8 N/m. Each figure was measured three times.
- The thickness step takes 0.29–0.30 s at 3 passes.
- The response adds 0.03–0.05 s to the momentum step's 0.24–0.25 s.
- The faces' push adds nothing measurable to the transport's 0.31–0.37 s.
- Together 6b adds about 0.33 s a sub-step.
- At the melt's rim, the transport made the predicted heights to 6.5 µm of a change of 0.15–0.17 mm.

**Rulings** (2026-10-21; each is what, then why, then the cost if wrong):
- *Correction form.* The momentum step is given `push`'s whole G, and the thickness step adds the change. Why: the
  slurry's shear-rate fixed point then sees the whole push, and the curvature's own nonlinearity enters at the start of
  the sub-step. Cost if wrong: none for fixed viscosities, as both forms give the same velocities.
- *The faces' own push in the transport, the stored momentum from the column push* (the coordinator's ruling, the
  column form of Rhie–Chow interpolation). Why: the bump above. Cost if wrong: the mass and the momentum are carried at
  slightly different velocities at the grid scale.
- *The heights predicted from the transport linearised on the columns, by Newton's method on the crossing cells, stopped
  on the flux missed* (the coordinator's ruling). Why: the transport is implicit in the cells' masses, which a column
  model of it can follow but not reproduce. Cost if wrong: where the cells of a column carry very different shares, the
  prediction differs from the transport by a few per cent of the change, and the thickness terms are that much
  explicit.
- *The face's push by Task 1's clipped cotangent weights, the deceleration's term with the two columns' mean g·n.* Why:
  the energy argument needs a non-negative exchange, and the spec's form is ρ_l (g·n) ∇H. Cost if wrong: 4 of the
  28 245 faces on the 2 mm sphere are clipped.
- *A wall takes its flowing neighbours' mean curvature in the column push's gradient.* Why: the face carries no
  capillary push into the wall, and the column's mean of its face velocities is then about half the push from its
  liquid side. Cost if wrong: through the cells' inertia, a velocity toward a wall carries metal onto its top a sub-step
  later.
- *A crust as one unknown* (the coordinator's ruling). Why: section above. Cost if wrong: none, as it is the rigid limit
  of Task 4's coupling.
- *The tests on strips two corners wide* (the coordinator's ruling). Why: with the film's curvature an open edge is a
  mirror for surface tension. Cost if wrong: the strips' ends disturb the cosines, so the tests measure the middle wave.

- [ ] **Step 1: The film curvature's tests** — append to `tests/test_reentry_model_columns.py`:

```diff
--- a/tests/test_reentry_model_columns.py
+++ b/tests/test_reentry_model_columns.py
@@ -158,3 +158,65 @@
     tangential = closure - np.einsum("ij,ij->i", closure, col.normal)[:, None] * col.normal
     assert np.linalg.norm(tangential, axis=1).max() < 0.075 * col.face_length.mean()
     assert np.linalg.norm(tangential, axis=1).mean() < 5e-4 * col.face_length.mean()
+
+
+# -- the thin film's curvature (sub-plan 21 Task 6b) --------------------------------------------------------------------
+def test_the_film_curvature_is_the_base_s_plus_the_top_s_cotangent_laplacian_of_the_heights(coarse_sphere_mesh):
+    """On a plane a cosine's curvature is the compact second difference of the heights (to 4e-7, the cotangents of its
+    own tilted top), and on the test sphere a uniform layer's is the base's, Task 2's, exactly. The cotangents are the
+    top surface's, but a cotangent depends on angles alone: on a sphere a uniform layer's top has the base's angles but
+    for the corners' normals not being quite radial (the two operators on a smooth bump differ by 0.35 % on average at
+    3.4 mm, measured 2026-10-21), and the area is the column's, the base's, which the energy argument requires."""
+    pts, faces, inner = flat_grid(n=24, jitter=0.0)
+    col = columns.Columns(pts, faces)
+    H = 1e-6 * np.cos(2.0 * np.pi * pts[:, 0] / 8e-3)
+    k2 = (2.0 * np.sin(np.pi * 1e-3 / 8e-3) / 1e-3) ** 2
+    np.testing.assert_allclose(col.film_curvature(H)[inner], k2 * H[inner], rtol=0, atol=1e-6 * k2 * 1e-6)
+    col, _ = sphere_columns(coarse_sphere_mesh)
+    uniform = np.full(col.n, 3.4e-3)
+    np.testing.assert_allclose(col.film_curvature(uniform), col.curvature(), rtol=0, atol=1e-9)
+    bump = 1e-6 * np.cos(col.points[:, 0] / 0.02)
+    top, base = col.film_operator(uniform) @ bump, col.film_operator(np.zeros(col.n)) @ bump
+    big = np.abs(base) > 0.5 * np.abs(base).max()
+    assert np.mean(top[big] / base[big]) == pytest.approx(1.0, abs=0.01)
+
+
+def test_a_pile_up_lowers_its_own_curvature_at_an_obtuse_corner_of_the_top_but_raises_the_film_s(coarse_sphere_mesh):
+    """Why surface tension takes the thin film's curvature (Task 6b; Asha, 2026-10-21): a corner touching an obtuse
+    triangle takes its neighbours' mean curvature (Task 2), so raising it makes its neighbours concave and lowers its
+    own -- a pile-up there would draw metal in, surface tension running backwards: raising each of the test sphere's 33
+    such corners by 1 um lowers its own curvature at 32 and leaves it unchanged at the other, by up to 0.083 1/m, where
+    every other corner's rises by at least 0.21 1/m. The film's curvature rises at every one of them, by 0.19-0.75 1/m
+    (measured 2026-10-21)."""
+    col, _ = sphere_columns(coarse_sphere_mesh)
+    H = np.full(col.n, 3e-4)
+    _, _, cot = columns.triangle_geometry(col.top_points(H), col.faces)
+    obtuse = np.zeros(col.n, dtype=bool)
+    for v in range(3):
+        np.logical_or.at(obtuse, col.faces[:, v], (cot < 0.0).any(axis=1))
+    corners = np.flatnonzero(obtuse)
+    assert len(corners) > 20
+    own_top, own_film = [], []
+    for c in corners:
+        bumped = H.copy()
+        bumped[c] += 1e-6
+        own_top.append(col.curvature(col.top_points(bumped))[c] - col.curvature(col.top_points(H))[c])
+        own_film.append(col.film_curvature(bumped)[c] - col.film_curvature(H)[c])
+    own_top, own_film = np.array(own_top), np.array(own_film)
+    assert np.all(own_top <= 0.0) and (own_top < 0.0).sum() >= 30 and np.all(own_film > 0.0)
+
+
+def test_the_faces_exchange_times_the_film_operator_has_no_growing_mode(coarse_sphere_mesh):
+    """The film's curvature is (1 / A_i) dE / dH_i for the thin film's surface energy E, A_i the column's area that turns
+    its volume into its height (Asha, 2026-10-21): the faces' exchange B diag(w) B^T / A times the film operator is a
+    product of two symmetric non-negative matrices once scaled by the areas, so no pattern of heights grows by itself --
+    on the test sphere with a rough top (0-3.4 mm) every eigenvalue is real and at least -1e-14 of the largest (measured
+    3.3e-5 against 7.7e11 for a uniform top, 2026-10-21; Task 2's curvature of the top gave -1.5e10, a growth of about
+    80 1/s for 0.3 mm of liquid)."""
+    import scipy.linalg as la
+    col, _ = sphere_columns(coarse_sphere_mesh)
+    H = np.random.default_rng(3).uniform(0.0, 3.4e-3, col.n)
+    B = col.incidence.toarray()
+    K = (B * col.conduct_weight[None, :]) @ B.T / col.area[:, None]
+    ev = la.eigvals(K @ col.film_operator(H).toarray())
+    assert np.abs(ev.imag).max() <= 1e-9 * ev.real.max() and ev.real.min() > -1e-14 * ev.real.max()
```

- [ ] **Step 2: The thickness terms' tests**, and Task 5's two capillary tests re-measured, in
  `tests/test_reentry_model_skin.py`:

```diff
--- a/tests/test_reentry_model_skin.py
+++ b/tests/test_reentry_model_skin.py
@@ -793,7 +793,8 @@
 
 def test_a_bump_s_surface_tension_pushes_the_metal_from_its_crest_converging_at_second_order():
     """A cosine bump of 20 um on a 0.5 mm layer: kappa = -H'' for small slopes, so -sigma grad kappa = sigma H''' =
-    sigma eps k^3 sin(kx). Measured 2026-10-20: 2.45 % of its peak off at 20 spacings a wavelength, 0.62 % at 40."""
+    sigma eps k^3 sin(kx). Measured 2026-10-20: 2.45 % of its peak off at 20 spacings a wavelength, 0.62 % at 40; the
+    same to three digits with the thin film's curvature (Task 6b, 2026-10-21)."""
     def error(lam):
         col, pts = regular_grid()
         k, eps = 2.0 * np.pi / lam, 2e-5
@@ -808,13 +809,14 @@
 
 def test_a_sphere_s_own_mesh_drives_no_capillary_flow_to_speak_of(coarse_sphere_mesh):
     """Uniform 0.8 mm skins on the 4 mm sphere: the curvature is uniform, so the capillary push is the mesh's error
-    alone -- 11 Pa/m at worst, 0.28 on average (measured 2026-10-20), against 2e5 Pa/m from the pressure at peak
-    heating (fact 153)."""
+    alone -- with the thin film's curvature (Task 6b) the base's, Task 2's operator: 0.32 Pa/m at worst, 0.031 on
+    average (measured 2026-10-21; 11 and 0.28 from the top surface's own curvature, 2026-10-20), against 2e5 Pa/m from
+    the pressure at peak heating (fact 153)."""
     from test_reentry_model_columns import sphere_columns
     col, _ = sphere_columns(coarse_sphere_mesh)
     _, terms = skin.push(col, np.full(col.n, 8e-4), np.zeros(col.n), np.zeros(3), V_HAT, RHO_L, SIGMA)
     size = np.linalg.norm(terms["capillary"], axis=1)
-    assert size.max() < 20.0 and size.mean() < 0.5
+    assert size.max() < 0.7 and size.mean() < 0.07
 
 
 # -- sideways transport (sub-plan 21 Task 6a) ---------------------------------------------------------------------------
@@ -1108,3 +1110,281 @@
     for x, y in zip(first.state(), second.state()):
         np.testing.assert_allclose(x, y, rtol=1e-11, atol=1e-11 * np.abs(y).max())
     np.testing.assert_allclose(a.wall, b.wall, rtol=0, atol=1e-11 * np.abs(a.momentum).max())
+
+
+# -- the thickness terms, implicit (sub-plan 21 Task 6b) ----------------------------------------------------------------
+def film_strip(nx, h, h0, T, ncells, eps=0.0, lam=None):
+    """Columns on a strip two corners wide (`strip_grid`), each column ncells cells at T holding h0 of thickness --
+    perturbed by eps cos(2 pi x / lam) -- over five 10 um cells of rigid metal, which the cells' cosine leaves flat."""
+    col, pts = strip_grid(nx, 2, h)
+    mat = scheil()
+    sk = skin.Skins(mat, col.area, skin.SkinSettings(**CELLS_04))
+    top = h0 + (eps * np.cos(2.0 * np.pi * pts[:, 0] / lam) if eps else 0.0)
+    m = np.zeros((col.n, ncells + 5))
+    m[:, :ncells] = (top * float(sk.density(np.array([T]))[0]) * col.area / ncells)[:, None]
+    m[:, ncells:] = (mat.rho * col.area * 1e-5)[:, None]
+    sk._install(m, m * mat.enthalpy(np.r_[np.full(ncells, T), np.full(5, RIGID)])[None, :], np.full(col.n, ncells + 5))
+    return sk, col, pts
+
+
+def flow_step(sk, col, dt, g, sigma=SIGMA, dp_ds=None):
+    """One flow sub-step in the order the skin body is to call it (6c): the push, surface tension where a column's top
+    flows; the momentum step with its response; the thickness step; the transport with the faces' own push."""
+    G, terms = skin.push(col, sk.thickness(), np.zeros(col.n) if dp_ds is None else dp_ds, g, V_HAT, RHO_L, sigma,
+                         wet=skin.flowing_tops(sk))
+    mb = skin.momentum_step(sk, dt, G, np.zeros(3), sk.area, response=True)
+    tb = skin.thickness_step(sk, col, dt, mb, terms["thickness"] + terms["capillary"], g, RHO_L, sigma)
+    return mb, tb, skin.transport(sk, col, dt, regrid=False, faces=tb.faces)
+
+
+def thin_film_rate(T, h0, k, g_out):
+    """The thin-film linear rate s(k) = h0^3 k^2 (rho g - sigma k^2) / (3 mu) [1/s], mu the slurry's at rest."""
+    mat = scheil()
+    mu = float(rheology.li_viscosity(1.0 - float(mat.liquid_fraction(np.array([T]))[0]), 0.0))
+    return h0 ** 3 * k ** 2 * (RHO_L * g_out - SIGMA * k ** 2) / (3.0 * mu)
+
+
+def cosine_rate(nx, h, lam, g_out, dt, span, waves):
+    """The rate [1/s] at which a cosine 1 um high on 1 mm of slurry at 902 K grows or decays on a strip, measured over
+    the second half of `span` seconds on the middle wave."""
+    sk, col, pts = film_strip(nx, h, 1e-3, SLURRY, 25, eps=1e-6, lam=lam)
+    x = pts[:, 0]
+    mid = (x > lam - 1e-9) & (x < (waves - 1) * lam + 1e-9)
+    c = np.cos(2.0 * np.pi * x / lam)[mid]
+    amplitude = lambda: float(((sk.thickness()[mid] - sk.thickness()[mid].mean()) * c).sum() / (c * c).sum())
+    N = int(round(span / dt))
+    a = []
+    for k in range(N):
+        flow_step(sk, col, dt, np.array([0.0, 0.0, g_out]))
+        if k + 1 in (N // 2, N):
+            a.append(amplitude())
+    return np.log(a[1] / a[0]) / ((N - N // 2) * dt)
+
+
+def test_a_cosine_on_a_slurry_layer_pulled_out_of_its_surface_grows_at_the_thin_film_rayleigh_taylor_rate():
+    """Spec 11.1 level 1. 1 mm of slurry at 902 K (30 % solid, 2.68 Pa s at rest) on a strip of 1 mm columns, pulled out
+    of its surface at 100 m/s2 (a capillary cutoff 2 pi sqrt(sigma / rho g) of 11.5 mm), a cosine 24 mm long and 1 um
+    high: the lubrication regime the formula assumes, the layer's viscous time rho h0^2 / mu 0.9 ms against a growth time
+    of 0.63 s. s(k) = h0^3 k^2 (rho g - sigma k^2) / (3 mu) = 1.578 1/s, met within -0.41, -0.22 and +0.42 % at sub-steps
+    of 1, 2.5 and 10 ms (measured 2026-10-21): the grid's own -0.40 % (its compact second difference) and backward
+    Euler's growth, +dt s / 2."""
+    s = thin_film_rate(SLURRY, 1e-3, 2.0 * np.pi / 24e-3, 100.0)
+    for dt in (1e-3, 2.5e-3, 1e-2):
+        assert cosine_rate(73, 1e-3, 24e-3, 100.0, dt, 0.63, 3) == pytest.approx(s, rel=0.01)
+
+
+def test_a_cosine_shorter_than_the_cutoff_decays_at_the_thin_film_rate():
+    """The same layer and pull, a cosine 8 mm long, below the 11.5 mm cutoff, on 0.25 mm columns: s = -19.44 1/s, met
+    within 1.2 % at 1 ms sub-steps; at 10 ms the rate is 8.7 % slower, backward Euler's own damping, -ln(1 + dt |s|) / dt,
+    which it meets within 0.10 % (0.26 % at 1 ms; measured 2026-10-21)."""
+    s = thin_film_rate(SLURRY, 1e-3, 2.0 * np.pi / 8e-3, 100.0)
+    rates = {dt: cosine_rate(97, 0.25e-3, 8e-3, 100.0, dt, 0.06, 3) for dt in (1e-3, 1e-2)}
+    assert rates[1e-3] == pytest.approx(s, rel=0.025)
+    for dt, rate in rates.items():
+        assert rate == pytest.approx(-np.log(1.0 - dt * s) / dt, rel=0.006)
+
+
+def test_the_two_column_mode_of_liquid_on_2_mm_columns_decays_at_any_sub_step():
+    """0.3 mm of liquid aluminium (viscous time 0.17 s: inertial at every sub-step here) on 2 mm columns, pulled out of
+    its surface at 30 m/s2 (the 100 mm sphere's peak deceleration) or not at all, its thickness alternating by 1 um from
+    column to column: surface tension relaxes it within each sub-step of 10 ms to 1 s -- to 0.112, 5.5e-4 and -0.051 of
+    itself in one sub-step of 10 ms, 100 ms and 1 s under the pull (0.117, -6.3e-5 and -0.056 without), never above
+    0.009 of itself after that (measured 2026-10-21)."""
+    for g_out in (30.0, 0.0):
+        for dt in (1e-2, 1e-1, 1.0):
+            sk, col, pts = film_strip(33, 2e-3, 0.3e-3, LIQUID, 15)
+            i = np.rint(pts[:, 0] / 2e-3).astype(int)
+            sign = np.where(i % 2 == 0, 1.0, -1.0)
+            m, H = sk.m.copy(), sk.H.copy()
+            m[:, :15] *= (1.0 + 1e-6 / 0.3e-3 * sign)[:, None]
+            H[:, :15] *= (1.0 + 1e-6 / 0.3e-3 * sign)[:, None]
+            sk._install(m, H, sk.count, T=sk.T)
+            inner = (i > 3) & (i < 29)
+            mode = lambda: float(((sk.thickness() - sk.thickness()[inner].mean()) * sign)[inner].mean())
+            a = [mode()]
+            for _ in range(3 if g_out else 6):
+                flow_step(sk, col, dt, np.array([0.0, 0.0, g_out]))
+                a.append(mode())
+            assert abs(a[1]) < 0.25 * a[0] and max(abs(x) for x in a[2:]) < 0.02 * a[0]
+
+
+def test_with_neither_surface_tension_nor_a_pull_into_the_surface_the_step_is_the_momentum_step_and_the_transport(
+        coarse_sphere_mesh):
+    """Spec 6.4's thickness terms vanish with sigma = 0 and g = 0: the momentum step with its response, the thickness
+    step and the transport with the faces' push leave exactly what the momentum step and the transport leave, to the last
+    bit, on mixed sphere skins pushed by a Newtonian pressure gradient."""
+    from test_reentry_model_columns import sphere_columns
+    col, _ = sphere_columns(coarse_sphere_mesh)
+    cos = np.clip(col.normal @ V_HAT, 0.0, 1.0)
+    dp_ds = -2.0 * 11e3 * np.sqrt(1.0 - cos ** 2) * cos / 0.05
+    first, _, _ = mixed_sphere_skins(coarse_sphere_mesh)
+    second, _, _ = mixed_sphere_skins(coarse_sphere_mesh)
+    G, _ = skin.push(col, first.thickness(), dp_ds, np.zeros(3), V_HAT, RHO_L, 0.0)
+    a = skin.momentum_step(first, 0.01, G, np.zeros(3), first.area)
+    ta = skin.transport(first, col, 0.01, regrid=False)
+    b, tb, tr = flow_step(second, col, 0.01, np.zeros(3), sigma=0.0, dp_ds=dp_ds)
+    assert all(np.array_equal(x, y) for x, y in zip(first.state(), second.state()))
+    assert all(np.array_equal(getattr(a, k), getattr(b, k)) for k in ("push", "shear", "wall"))
+    assert all(np.array_equal(getattr(ta, k), getattr(tr, k)) for k in ("mass", "enthalpy", "momentum", "wall"))
+    assert not (tb.push.any() or tb.wall.any())
+
+
+def test_the_books_close_with_the_thickness_impulse_and_the_masses_stay_positive(coarse_sphere_mesh):
+    """Mixed sphere skins -- liquid, slurry, crusts and anchored metal, each column 1-3 times its minimum deep -- under a
+    Newtonian pressure gradient, 30 m/s2 of deceleration and surface tension: mass and enthalpy are kept, and the
+    momentum changes by the push's, the shear's and the wall's impulses of the momentum step, the thickness step's push
+    and wall, and the transport's carried momentum and wall, to round-off (7.8e-14 and 8.8e-13 of the largest impulse at
+    1 and 10 ms, measured 2026-10-21); every cell keeps a positive mass."""
+    from test_reentry_model_columns import sphere_columns
+    col, _ = sphere_columns(coarse_sphere_mesh)
+    cos = np.clip(col.normal @ V_HAT, 0.0, 1.0)
+    dp_ds = -2.0 * 11e3 * np.sqrt(1.0 - cos ** 2) * cos / 0.05
+    for dt in (1e-3, 1e-2):
+        sk, _, _ = mixed_sphere_skins(coarse_sphere_mesh)
+        P0, M0, E0 = sk.column_momentum(), sk.mass(), sk.energy()
+        mb, tb, tr = flow_step(sk, col, dt, 30.0 * V_HAT, dp_ds=dp_ds)
+        parts = (mb.push, mb.shear, mb.wall, tb.push, tb.wall, tr.momentum, tr.wall)
+        scale = max(np.abs(x).max() for x in parts + (P0,))
+        miss = np.abs(sk.column_momentum() - P0 - sum(parts)).max() / scale
+        assert miss < 2e-12
+        assert sk.mass() == pytest.approx(M0, rel=1e-13) and sk.energy() == pytest.approx(E0, rel=1e-13)
+        assert np.all(sk.m[sk.valid] > 0.0) and np.abs(tb.push).max() > 1e-3 * np.abs(mb.push).max()
+
+
+def test_surface_tension_acts_only_where_a_column_s_top_flows():
+    """Asha, 2026-10-21: liquid 0.29 mm deep on a strip beside a wall of rigid metal 0.4 mm high, and a crust over
+    slurry in the middle. The wall's height enters the curvature of the liquid beside it, as a wall's would; neither the
+    wall nor the crust feels a capillary push, and surface tension moves no metal into either (what lands there is the
+    round-off of a velocity the faces' push cancels, 1.5e-15 of what moves); the wall's cells stay as they were and
+    the crust and the slurry beneath it, pushed by nothing, stay at rest."""
+    col, pts = strip_grid(21, 2)
+    x = np.rint(pts[:, 0] / 1e-3).astype(int)
+    temps = [[RIGID] * 40 if k >= 15 else [RIGID] * 2 + [SLURRY] * 20 + [RIGID] * 5 if k == 7
+             else [LIQUID] * 25 + [RIGID] * 5 for k in x]
+    sk = stacked(col, temps)
+    wet = skin.flowing_tops(sk)
+    assert np.array_equal(wet, (x < 15) & (x != 7))
+    H = sk.thickness()
+    kappa = col.film_curvature(H)
+    beside = np.flatnonzero(x == 14)
+    np.testing.assert_allclose(kappa[beside], (col.film_operator(H) @ H)[beside], rtol=1e-12)
+    assert np.all(kappa[beside] < -50.0)                                    # a wall higher than the liquid: concave
+    G, terms = skin.push(col, H, np.zeros(col.n), np.zeros(3), V_HAT, RHO_L, SIGMA, wet=wet)
+    assert np.all(terms["capillary"][~wet] == 0.0) and np.abs(terms["capillary"][wet]).max() > 1e3
+    before = sk.state()
+    mb, tb, tr = flow_step(sk, col, 0.01, np.zeros(3))
+    assert np.abs(tr.mass[~wet]).max() < 1e-14 * np.abs(tr.mass).max()
+    top = sk.count - before[4]                                              # round-off slivers landed on rigid tops
+    for c in np.flatnonzero(~wet):
+        for now, then in ((sk.m, before[1]), (sk.H, before[2])):
+            assert np.array_equal(now[c, top[c]:sk.count[c]], then[c, :before[4][c]])
+    for c in np.flatnonzero(x == 7):
+        assert np.all(sk.U[c, top[c]:sk.count[c]] == 0.0)
+
+
+def test_a_uniform_layer_on_the_test_sphere_stays_uniform(coarse_sphere_mesh):
+    """0.3 mm of liquid over rigid metal on every column of the 4 mm test sphere, surface tension only: a sphere's
+    curvature is the same everywhere, so the layer stays uniform but for the base curvature's error (Task 2's, within
+    8.4e-5 of 2 / R): 26 nm and 2.2 um/s at most after 20 sub-steps of 10 ms (measured 2026-10-21)."""
+    from test_reentry_model_columns import sphere_columns
+    col, _ = sphere_columns(coarse_sphere_mesh)
+    sk = stacked(col, [[LIQUID] * 30 + [RIGID] * 5] * col.n)
+    H0 = sk.thickness()
+    for _ in range(20):
+        flow_step(sk, col, 0.01, np.zeros(3))
+    assert np.abs(sk.thickness() - H0).max() < 6e-8 and np.abs(sk.U).max() < 5e-6
+
+
+def test_a_one_column_bump_relaxes_through_the_faces_own_push():
+    """A bump 1 um high on one column of 1.1 mm of slurry on 1 mm columns, surface tension only: the faces' own push --
+    a two-point difference across each face, which sees a disturbance one column wide -- relaxes it, to 12.5 % of itself
+    in one 10 ms sub-step and 3.0 % in five (measured 2026-10-21). The columns' pushes are centred differences, zero on
+    the bump's own column by symmetry, and with them alone (no faces' push) the transport, which takes metal out of a
+    column only at its own velocity, leaves the bump at 99.8 % of its height (the column form of Rhie-Chow
+    interpolation)."""
+    out = {}
+    for faces in (True, False):
+        col, pts = strip_grid(21, 21)
+        sk = stacked(col, [[SLURRY] * 100 + [RIGID] * 5] * col.n)
+        c = int(np.argmin(np.linalg.norm(pts - pts.mean(axis=0), axis=1)))
+        base = sk.thickness()[c]
+        m, H = sk.m.copy(), sk.H.copy()
+        m[c, :100] *= 1.0 + 1e-6 / (base - 5e-5)
+        H[c, :100] *= 1.0 + 1e-6 / (base - 5e-5)
+        sk._install(m, H, sk.count, T=sk.T)
+        heights = []
+        for _ in range(5):
+            G, terms = skin.push(col, sk.thickness(), np.zeros(col.n), np.zeros(3), V_HAT, RHO_L, SIGMA,
+                                 wet=skin.flowing_tops(sk))
+            mb = skin.momentum_step(sk, 0.01, G, np.zeros(3), sk.area, response=True)
+            applied = terms["thickness"] + terms["capillary"]
+            tb = skin.thickness_step(sk, col, 0.01, mb, applied, np.zeros(3), RHO_L, SIGMA)
+            skin.transport(sk, col, 0.01, regrid=False, faces=tb.faces if faces else None)
+            heights.append(sk.thickness()[c] - base)
+        out[faces] = np.array(heights) / 1e-6
+    assert out[True][0] < 0.25 and out[True][-1] < 0.06 and np.all(out[False] > 0.99)
+
+
+def test_a_floating_crust_is_one_unknown_and_the_momentum_books_close_on_mixed_skins(coarse_sphere_mesh):
+    """Each floating crust is one unknown in the momentum step's system, a rigid body, where 6a's coupled its cells by
+    a viscosity of 1e9 Pa s, whose rows some 1e12 times the cells' inertia amplified the solver's round-off: on mixed
+    sphere skins the books closed to 1.2e-5 of the largest impulse and now to round-off; the velocities moved by 1.1e-3
+    of the largest at 10 ms, as the coupled crust's moved by 1.1e-5 at 1e7 Pa s and 4.2e-7 at 1e5 (measured
+    2026-10-21), so the change is that round-off, and columns without a crust are unchanged to the last bit. The books
+    close to 6.6e-14 and 6.2e-13 of the largest impulse at 1 and 10 ms."""
+    from test_reentry_model_columns import sphere_columns
+    col, _ = sphere_columns(coarse_sphere_mesh)
+    for dt in (1e-3, 1e-2):
+        sk, _, _ = mixed_sphere_skins(coarse_sphere_mesh)
+        flowing = sk.valid & np.isfinite(skin.cell_viscosity(sk, 0.0))
+        crust = sk.valid & ~flowing & ~skin.anchored_cells(sk, flowing)
+        assert crust.sum() > 1000
+        P0 = sk.column_momentum()
+        book = skin.momentum_step(sk, dt, np.array([3e4, -1e4, 2e3]), np.zeros(3), sk.area)
+        parts = (book.push, book.shear, book.wall, P0, sk.column_momentum())
+        miss = (np.abs(sk.column_momentum() - P0 - book.push - book.shear - book.wall).max()
+                / max(np.abs(x).max() for x in parts))
+        assert miss < 1.5e-12
+        same = crust[:, 1:] & crust[:, :-1]
+        assert np.array_equal(sk.U[:, 1:][same], sk.U[:, :-1][same])
+
+
+def test_the_momentum_response_is_what_a_change_of_the_push_gives():
+    """With the viscosities held -- liquid aluminium's here -- the momentum step is linear in the push: a change dG of
+    the push changes every cell's velocity by its response times dG, and the wall's impulse by -dt wall_response dG."""
+    sk = column_at([950.0] * 30 + [880.0] * 10)
+    U0 = np.random.default_rng(2).normal(0.0, 0.01, sk.U.shape) * sk.valid[..., None]
+    G, dG = np.array([2e4, -1e4, 0.0]), np.array([300.0, 500.0, 0.0])
+    sk.U = U0.copy()
+    a = skin.momentum_step(sk, 0.01, G, np.zeros(3), sk.area, response=True)
+    Ua = sk.U.copy()
+    sk.U = U0.copy()
+    b = skin.momentum_step(sk, 0.01, G + dG, np.zeros(3), sk.area)
+    np.testing.assert_allclose(sk.U - Ua, a.response[..., None] * dG, rtol=0, atol=1e-12 * np.abs(Ua).max())
+    np.testing.assert_allclose(b.wall - a.wall, -0.01 * a.wall_response[:, None] * dG, rtol=1e-9, atol=1e-20)
+    assert np.all(a.response[0, :30] > 0.0) and np.all(a.response[0, 30:] == 0.0)
+
+
+def test_the_predicted_heights_are_the_transport_s_and_the_linearisations_stop_on_the_flux(coarse_sphere_mesh):
+    """Liquid aluminium 0.29 mm deep with 1 um of random roughness on the test sphere, surface tension only: the
+    roughness decays, the heights the thickness step predicts are those the transport makes, and the linearisations of
+    the crossing cells stop within a few passes, on the flux they miss: over four sub-steps of 10 ms the roughness falls
+    to 34 % of itself, the transport makes the predicted heights to 0.62 % of the largest change, and the passes are
+    1, 3, 3 and 3 (measured 2026-10-21)."""
+    from test_reentry_model_columns import sphere_columns
+    col, _ = sphere_columns(coarse_sphere_mesh)
+    sk = stacked(col, [[LIQUID] * 25 + [RIGID] * 5] * col.n)
+    rng = np.random.default_rng(1)
+    f = (1.0 + 1e-6 / 2.9e-4 * rng.normal(size=col.n))[:, None]
+    m, H = sk.m.copy(), sk.H.copy()
+    m[:, :25] *= f
+    H[:, :25] *= f
+    sk._install(m, H, sk.count, T=sk.T)
+    spread0 = np.std(sk.thickness())
+    worst, passes = 0.0, []
+    for _ in range(4):
+        H0 = sk.thickness()
+        mb, tb, tr = flow_step(sk, col, 0.01, np.zeros(3))
+        worst = max(worst, np.abs(tb.heights - sk.thickness()).max() / np.abs(sk.thickness() - H0).max())
+        passes.append(tb.passes)
+    assert worst < 0.0125 and max(passes) <= 4 and np.std(sk.thickness()) < 0.5 * spread0
```

- [ ] **Step 3: Run them to see them fail.** `"$PY" -m pytest tests/test_reentry_model_columns.py
  tests/test_reentry_model_skin.py -q`. Expected: 15 failed:
  - `AttributeError: 'Columns' object has no attribute 'film_curvature'` (2);
  - `AttributeError: 'Columns' object has no attribute 'film_operator'` (1);
  - `AttributeError: module 'reentry_model.skin' has no attribute 'flowing_tops'` (9);
  - `TypeError: momentum_step() got an unexpected keyword argument 'response'` (1);
  - the crust's books, `assert ... 1.32e-05 < 1.5e-12` (1);
  - Task 5's sphere test, `assert (... 10.99 < 0.7)` (1).

- [ ] **Step 4: The thin film's curvature** in `reentry_model/columns.py`:

```diff
--- a/reentry_model/columns.py
+++ b/reentry_model/columns.py
@@ -30,7 +30,21 @@
 its second, so that a velocity u carries u . face_vector per unit depth across the face. A column's faces close
 exactly on a plane, so a uniform flow carries as much into every column as out of it; the edge's direction times the
 face's length closes only where the triangles are symmetric about their edges, and on the tests' jittered plane of 1 mm
-squares a uniform 1 m/s would gather or spread a column's metal at up to 539 1/s that way (measured 2026-10-21)."""
+squares a uniform 1 m/s would gather or spread a column's metal at up to 539 1/s that way (measured 2026-10-21).
+
+Surface tension's pressure (sub-plan 21 Task 6b; Asha, 2026-10-21, option A) takes the curvature in its thin-film form,
+`film_curvature`: the base's own curvature (`curvature`, its obtuse-corner rule acting on the fixed base alone) plus
+`film_operator` times the thicknesses, minus the cotangent Laplacian of the top surface over each column's area. That
+is (1/A_i) dE/dH_i for the thin film's surface energy E measured on the top surface, A_i the column's area, which turns
+its volume into its height in the skins: the operator times the areas is symmetric and non-negative, as the faces'
+exchange is, so their product has no growing mode (on the test sphere with a rough top every eigenvalue is real and
+non-negative to round-off). The top's own curvature has one: a corner touching an obtuse triangle takes its neighbours'
+mean, cannot see its own height, and a pile-up on it lowers its pressure and draws metal in -- the faces' exchange times
+its linearisation has eigenvalues down to -1.5e10 1/m4 on the test sphere and -2.4e11 on the 2 mm sphere, a growth of
+about 1 300 1/s for 0.3 mm of liquid on the latter (measured 2026-10-21). A cotangent depends on angles alone, so the
+top's differ from the base's only where the thickness varies; with the base's area, a thick layer's bump curvature stays
+about 2 H / R above the top surface's own (14 % at 3.4 mm on the 50 mm sphere), as the columns hold their volume over the
+base's area."""
 from dataclasses import dataclass
 
 import numpy as np
@@ -152,6 +166,27 @@
         """The top surface: each corner moved out along its normal by its column's thickness [m]."""
         return self.points + np.asarray(thickness, dtype=float)[:, None] * self.normal
 
+    def film_operator(self, thickness):
+        """Minus the cotangent Laplacian of the top surface the thicknesses [m] make, over each column's area: a sparse
+        (n, n) matrix [1/m2], symmetric and non-negative once multiplied by the areas, whose product with the
+        thicknesses is the thin film's curvature less the base's (module docstring). The cotangents are the top's."""
+        _, _, cot = triangle_geometry(self.top_points(thickness), self.faces)
+        rows, cols, vals = [], [], []
+        for v in range(3):
+            i, j, k = self.faces[:, v], self.faces[:, (v + 1) % 3], self.faces[:, (v + 2) % 3]
+            cj, ck = cot[:, (v + 1) % 3], cot[:, (v + 2) % 3]
+            rows += [i, i, i]
+            cols += [j, k, i]
+            vals += [-0.5 * ck, -0.5 * cj, 0.5 * (ck + cj)]
+        L = sp.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(self.n, self.n))
+        return sp.diags(1.0 / self.area) @ L
+
+    def film_curvature(self, thickness):
+        """The thin film's mean curvature [1/m] at every column: the base's (`curvature`, whose obtuse-corner rule acts
+        on the base alone) plus `film_operator` times the thicknesses (module docstring)."""
+        H = np.asarray(thickness, dtype=float)
+        return self.curvature() + self.film_operator(H) @ H
+
     def laplacian(self, f, weights=None):
         """sum over neighbours j of w_ij (f_j - f_i) per corner (f per corner, any trailing shape)."""
         w = self.cot_weight if weights is None else weights
```

- [ ] **Step 5: The momentum's response, each crust one unknown, the push's thin-film capillary term, the faces' push in
  the transport and the thickness step** in `reentry_model/skin.py`:

```diff
--- a/reentry_model/skin.py
+++ b/reentry_model/skin.py
@@ -821,7 +821,8 @@
 # -- momentum through the depth (sub-plan 21 Task 4) --------------------------------------------------------------------
 SHEAR_PASSES = 5                  # passes of the fixed point on the cells' shear rates at most (spec 6.4)
 SHEAR_TOL = 1.0e-3                # a pass that moves no flowing cell's viscosity by more than this, relatively, ends it
-CRUST_VISCOSITY = 1.0e9           # Pa s: a floating crust's, so that it moves as one piece (some 1e8 times the slurry's)
+CRUST_VISCOSITY = 1.0e9           # Pa s: a floating crust's half-cell at its interface with flowing metal (some 1e8 times
+                                  # the slurry's, adding nothing to the drag); the crust itself is one unknown (Task 6b)
 
 
 @dataclass
@@ -832,6 +833,10 @@
     wall: np.ndarray         # of the anchored rigid metal -- the skin's rigid cells joined to the substrate and the 3D
                              # model beneath the base -- on the moving cells, and the momentum a cell loses on joining it
     passes: int              # the fixed point's passes
+    response: np.ndarray = None       # with response=True, (n, K) [m/s per N/m3]: every cell's velocity per unit push,
+                                      # the last pass's system solved for the push alone (zero where a cell is held)
+    wall_response: np.ndarray = None  # with response=True, (n,) [m3]: the anchored metal's drag on that response, so
+                                      # that a change dG of the push over the sub-step gives the wall -dt wall_response dG
 
 
 def anchored_cells(sk, flowing):
@@ -845,7 +850,7 @@
     return sk.valid & (j >= (sk.count - held)[:, None])
 
 
-def momentum_step(sk, dt, push, shear, A_top):
+def momentum_step(sk, dt, push, shear, A_top, response=False):
     """One implicit sub-step dt [s] of every column's momentum through its depth (spec 6.4), each velocity component
     apart:
 
@@ -856,13 +861,21 @@
     (n, 3)] on the top cell over its top area A_top [m2, (n,)]. The rigid metal joined to the substrate (`anchored_cells`)
     is held at rest; a moving cell drags on it, and a flowing bottom cell on the rigid 3D metal beneath the base, with no
     slip over its own half-thickness. A crust -- rigid cells with flowing metal between them and the substrate -- floats:
-    it moves as one piece (CRUST_VISCOSITY), dragged by the gas shear on its top, pushed like the rest and dragged by the
-    slurry beneath (Asha, 2026-10-20, option B; the column model has no sideways stress in a solid, so a crust over many
-    columns moves as column-sized plates, as a warm, barely solid crust would). Each flowing cell's viscosity is at its
-    own shear rate, the mean over its two faces of the stress there over its viscosity, found by fixed-point iteration
-    from the start-of-sub-step velocities: at most SHEAR_PASSES passes, ending when no flowing cell's viscosity moves by
-    more than SHEAR_TOL. Updates sk.U; returns a MomentumBook whose impulses add up to the columns' momentum change to
-    round-off."""
+    it moves as one piece, one unknown of the system (its mass the cells', its drags those at its top and bottom faces),
+    dragged by the gas shear on its top, pushed like the rest and dragged by the slurry beneath (Asha, 2026-10-20,
+    option B; the column model has no sideways stress in a solid, so a crust over many columns moves as column-sized
+    plates, as a warm, barely solid crust would). Task 4 coupled the crust's cells by a viscosity of 1e9 Pa s, whose rows,
+    some 1e12 times the cells' inertia, amplified the solver's round-off: on mixed sphere skins the books closed to
+    1.2e-5 of the largest impulse, now to 6e-14 (1 ms) and 6e-13 (10 ms), and the velocities moved by up to 1.1e-3 of the
+    largest -- the coupled crust's own round-off, which fell as 1.1e-4, 1.1e-5 and 9.6e-7 at 1e8, 1e7 and 1e6 Pa s;
+    columns without a crust are unchanged to the last bit (sub-plan 21 Task 6b, measured 2026-10-21). Each flowing
+    cell's viscosity is at its own shear rate, the mean over its two faces of the stress there over its viscosity, found
+    by fixed-point iteration from the start-of-sub-step velocities: at most SHEAR_PASSES passes, ending when no flowing
+    cell's viscosity moves by more than SHEAR_TOL. Updates sk.U; returns a MomentumBook whose impulses add up to the
+    columns' momentum change to round-off. With `response`, the book also carries the velocities a unit push alone
+    gives through the last pass's system and the wall's drag on them: with the viscosities held, the sub-step is linear
+    in the push, so a change of the push changes every velocity by its response times the change (sub-plan 21 Task
+    6b)."""
     n, K = sk.m.shape
     valid = sk.valid
     rows = np.arange(n)
@@ -875,13 +888,20 @@
     anchored = anchored_cells(sk, flowing)
     moving = valid & ~anchored                                            # the flowing cells and the floating crust
     crust = moving & ~flowing
-    P0 = sk.column_momentum()
     held = ((m * anchored)[..., None] * sk.U).sum(axis=1)                 # what anchored cells carried, now at rest
     U_old = np.where(moving[..., None], sk.U, 0.0)
     top_moves = moving[:, 0]
     force = np.where(moving[..., None], push[:, None, :] * (dz * A[:, None])[..., None], 0.0)
     force[:, 0] += np.where(top_moves[:, None], shear * A_top[:, None], 0.0)
     bottom_flows = flowing[rows, last] & (sk.count > 0)
+    # the unknowns: every cell, but a crust's cells one unknown together (sub-plan 21 Task 6b)
+    same = np.zeros((n, K), dtype=bool)
+    same[:, 1:] = crust[:, 1:] & crust[:, :-1]
+    unit = np.cumsum(~same, axis=1) - 1
+    Kr = int(unit.max()) + 1
+    lump = lambda x: np.bincount((rows[:, None] * Kr + unit).ravel(), x.ravel(), n * Kr).reshape(n, Kr)
+    moving_r, valid_r, last_r = lump(moving.astype(float)) > 0.0, lump(valid.astype(float)) > 0.0, unit[rows, last]
+    cross_r, cross_k = np.nonzero(~same[:, 1:])                          # the interfaces between unknowns
 
     def viscosities(gamma):
         return np.where(crust, CRUST_VISCOSITY, cell_viscosity(sk, gamma))
@@ -890,18 +910,26 @@
         Db = np.where(bottom_flows, eta[rows, last] / (0.5 * np.maximum(dz[rows, last], DZ_FLOOR)), 0.0)
         return series_drag(sk, eta), np.where(np.isfinite(Db), Db, 0.0)
 
+    def system(D, Db):
+        Dr = np.zeros((n, Kr - 1))
+        Dr[cross_r, unit[cross_r, cross_k]] = D[cross_r, cross_k]
+        diag = np.where(valid_r, lump(m) / dt, 1.0)
+        diag[:, :-1] += A[:, None] * Dr
+        diag[:, 1:] += A[:, None] * Dr
+        diag[rows, last_r] += A * Db
+        off = -A[:, None] * Dr
+        sub, sup = np.where(moving_r[:, 1:], off, 0.0), np.where(moving_r[:, :-1], off, 0.0)
+        return sub, np.where(moving_r, diag, 1.0), sup
+
+    def unknowns(sub, diag, sup, rhs):
+        """Each cell's value from the system on the unknowns, its right-hand side summed over each unknown's cells."""
+        return np.where(moving, thomas(sub, diag, sup, lump(np.where(moving, rhs, 0.0)))[rows[:, None], unit], 0.0)
+
     def solve(eta):
         D, Db = drags(eta)
-        diag = np.where(valid, m / dt, 1.0)
-        diag[:, :-1] += A[:, None] * D
-        diag[:, 1:] += A[:, None] * D
-        diag[rows, last] += A * Db
-        off = -A[:, None] * D
-        sub, sup = np.where(moving[:, 1:], off, 0.0), np.where(moving[:, :-1], off, 0.0)
-        diag = np.where(moving, diag, 1.0)
-        rhs = np.where(moving[..., None], (m / dt)[..., None] * U_old + force, 0.0)
-        U = np.stack([thomas(sub, diag, sup, rhs[..., c]) for c in range(3)], axis=-1)
-        return np.where(moving[..., None], U, 0.0), D, Db
+        sub, diag, sup = system(D, Db)
+        rhs = (m / dt)[..., None] * U_old + force
+        return np.stack([unknowns(sub, diag, sup, rhs[..., c]) for c in range(3)], axis=-1), D, Db
 
     def rates(U, D, Db, eta):
         """Every flowing cell's shear rate: the mean over its faces of the stress there over its viscosity."""
@@ -929,6 +957,12 @@
     book = MomentumBook(push=dt * np.where(moving[..., None], push[:, None, :] * (dz * A[:, None])[..., None], 0.0).sum(axis=1),
                         shear=dt * np.where(top_moves[:, None], shear * A_top[:, None], 0.0),
                         wall=dt * wall - held, passes=passes)
+    if response:                                                          # the last pass's system, for a unit push alone
+        r = unknowns(*system(D, Db), dz * A[:, None])
+        book.response = r
+        book.wall_response = A * Db * r[rows, last]
+        if K > 1:
+            book.wall_response += (A[:, None] * D * (moving[:, :-1] & anchored[:, 1:]) * r[:, :-1]).sum(axis=1)
     sk.U = U
     return book
 
@@ -945,7 +979,7 @@
     return np.asarray(magnitude, dtype=float)[:, None] * unit
 
 
-def push(columns, thickness, dp_ds, g, v_hat, rho_l, sigma):
+def push(columns, thickness, dp_ds, g, v_hat, rho_l, sigma, wet=None):
     """The push along the surface on every layer of every column [N/m3], (n, 3), and its four terms (spec 6.4), each in
     the column's own tangent plane -- the base's, where `momentum_step` solves for the velocities:
 
@@ -962,10 +996,16 @@
         thickness     rho_l (g . n) grad_s H, H [m, (n,)] the columns' thicknesses: a thicker column under a pull out
                       of the surface has the lower pressure beneath it and draws metal in (the front face's
                       Rayleigh-Taylor mechanism).
-        capillary     -sigma grad_s kappa, kappa the mean curvature (`Columns.curvature`) of the top surface the
-                      thicknesses make (`Columns.top_points`).
-
-    The two thickness terms are evaluated at the given H; Task 6 takes them implicitly in the column heights."""
+        capillary     -sigma grad_s kappa, kappa the thin film's mean curvature (`Columns.film_curvature`): the base's
+                      plus minus the cotangent Laplacian of the thicknesses on the top surface over each column's area
+                      (sub-plan 21 Task 6b; Asha, 2026-10-21, option A) -- not the top surface's own curvature, whose
+                      corners touching an obtuse triangle take their neighbours' mean, cannot see their own height, and
+                      make surface tension draw metal onto a pile-up there. It acts only where a column's top flows
+                      (`wet`, (n,) bool, `flowing_tops`; None for everywhere): not on an anchored rigid top nor on a
+                      floating crust; a rigid neighbour's height enters a flowing column's curvature as a wall's would,
+                      and in the gradient a wall takes its flowing neighbours' mean curvature (`_walled`).
+
+    The two thickness terms are evaluated at the given H; `thickness_step` takes them implicitly in the column heights."""
     H = np.asarray(thickness, dtype=float)
     normals = columns.normal
     g = np.broadcast_to(np.asarray(g, dtype=float), (columns.n, 3))
@@ -973,7 +1013,7 @@
     terms = {"pressure": -along_meridian(normals, v_hat, dp_ds),
              "deceleration": rho_l * (g - g_n[:, None] * normals),
              "thickness": rho_l * g_n[:, None] * columns.gradient(H),
-             "capillary": -sigma * columns.gradient(columns.curvature(columns.top_points(H)))}
+             "capillary": _capillary(columns, columns.film_curvature(H), wet, sigma)}
     return sum(terms.values()), terms
 
 
@@ -1052,10 +1092,12 @@
     return spla.splu(A)
 
 
-def transport(sk, columns, dt, regrid=True):
+def transport(sk, columns, dt, regrid=True, faces=None):
     """One linearly implicit upwind sub-step dt [s] of the metal moving sideways between neighbouring columns at the
     cells' velocities sk.U (spec 6.5; Asha, 2026-10-21: its own step, after the momentum step). Updates the cells and
-    returns a TransportBook.
+    returns a TransportBook. `faces`, a FacePush from `thickness_step` (sub-plan 21 Task 6b), adds to every moving cell's
+    velocity across each face its response times the face's own thickness push, less its column's; None leaves the
+    velocities as they are.
 
     Across each face, cell k of the donor column a sends c m_k' toward the receiver b, c = (u_k . f)+ / A_a [1/s], f the
     face's vector (`Columns.face_vector`) and m_k' the cell's mass at the sub-step's end; it carries the cell's specific
@@ -1090,6 +1132,8 @@
     valid = sk.valid
     flowing = valid & np.isfinite(cell_viscosity(sk, 0.0))
     moving = valid & ~anchored_cells(sk, flowing) & (sk.U != 0.0).any(axis=2)
+    if faces is not None:
+        moving |= valid & ~anchored_cells(sk, flowing) & (faces.response != 0.0)
     book = TransportBook(np.zeros(n), np.zeros(n), np.zeros((n, 3)), np.zeros((n, 3)), 0.0)
     if not moving.any():
         return book
@@ -1098,7 +1142,10 @@
     a, b, f = np.r_[i, j], np.r_[j, i], np.r_[columns.face_vector, -columns.face_vector]
     keep = np.flatnonzero(moving.any(axis=1)[a])
     a, b, f = a[keep], b[keep], f[keep]
-    c = np.where(moving[a], np.maximum(np.einsum("dkx,dx->dk", sk.U[a], f), 0.0), 0.0) / sk.area[a][:, None]
+    across = np.einsum("dkx,dx->dk", sk.U[a], f)
+    if faces is not None:
+        across = across + faces.response[a] * faces.across[keep][:, None]
+    c = np.where(moving[a], np.maximum(across, 0.0), 0.0) / sk.area[a][:, None]
     keep = np.flatnonzero((c > 0.0).any(axis=1))
     if not len(keep):
         return book
@@ -1198,6 +1245,155 @@
     return book
 
 
+# -- the thickness terms, implicit (sub-plan 21 Task 6b) ----------------------------------------------------------------
+THICKNESS_PASSES = 8              # linearisations of the transport's fluxes at most
+THICKNESS_TOL = 1.0e-3            # the linearisations stop once the flux the last one missed would move no column by
+                                  # more than this share of the largest height change it predicts
+
+
+def flowing_tops(sk):
+    """Whether each column's top cell flows, (n,) bool: where surface tension acts (Asha, 2026-10-21), not on an
+    anchored rigid top nor on a floating crust."""
+    return (sk.count > 0) & np.isfinite(cell_viscosity(sk, 0.0)[:, 0])
+
+
+def _walled(columns, kappa, wet):
+    """The curvatures `kappa` (n,) where a column's top flows; at a column whose top does not -- a wall -- the mean of
+    its flowing neighbours', so that the gradient beside a wall carries no capillary push of the wall's own."""
+    if wet is None or wet.all():
+        return kappa
+    i, j = columns.edges[:, 0], columns.edges[:, 1]
+    total = np.bincount(i, kappa[j] * wet[j], columns.n) + np.bincount(j, kappa[i] * wet[i], columns.n)
+    count = np.bincount(i, wet[j], columns.n) + np.bincount(j, wet[i], columns.n)
+    return np.where(wet, kappa, np.where(count > 0.0, total / np.maximum(count, 1.0), kappa))
+
+
+def _capillary(columns, kappa, wet, sigma):
+    """-sigma grad_s kappa where a column's top flows, zero elsewhere (`push`'s capillary term)."""
+    G = -sigma * columns.gradient(_walled(columns, kappa, wet))
+    return G if wet is None else np.where(wet[:, None], G, 0.0)
+
+
+@dataclass
+class FacePush:
+    """What `transport` adds to the cells' velocities across the faces."""
+    response: np.ndarray     # (n, K) [m/s per N/m3]: every cell's velocity per unit push (MomentumBook.response)
+    across: np.ndarray       # (2E,) [Pa]: for each face, from its edge's first column to its second and then the other
+                             # way (the order `transport` builds), the face's own thickness push times the face, less
+                             # the sending column's push dotted with the face's vector
+
+
+@dataclass
+class ThicknessBook:
+    """What the implicit thickness terms did in a sub-step, per column."""
+    push: np.ndarray         # N s (n, 3): the impulse of their change on the moving cells
+    wall: np.ndarray         # N s (n, 3): the anchored metal's drag on the velocities that change gave
+    heights: np.ndarray      # m (n,): the thicknesses predicted at the end of the transport
+    faces: FacePush          # for `transport`
+    growth: float            # the largest Rayleigh-Taylor growth in a sub-step, dt s, over the columns
+    passes: int              # the linearisations made
+
+
+def thickness_step(sk, columns, dt, book, applied, g, rho_l, sigma):
+    """The push's thickness terms, rho_l (g . n) grad_s H and -sigma grad_s kappa, implicit in the column heights at
+    the end of the sub-step's transport (spec 6.4's last paragraph; Asha, 2026-10-21: the transport is its own step
+    after the momentum step). Called between `momentum_step`, given `push`'s G with `wet=flowing_tops(sk)` and
+    response=True (its book is `book`), and `transport`, given the returned book's `faces`; `applied` (n, 3) [N/m3] is the
+    part of that push the two terms made, terms['thickness'] + terms['capillary']. Corrects sk.U; returns a
+    ThicknessBook.
+
+    The two terms are the gradient of P = sigma kappa - rho_l (g . n) H, kappa the thin film's curvature
+    (`Columns.film_curvature`), surface tension only where a column's top flows (`flowing_tops`). Across each dual face
+    P pushes by its two-point difference, phi = w (P_a - P_b) [Pa] with w the face's clipped cotangent weight
+    (`Columns.conduct_weight`; the capillary part only between two flowing tops, the deceleration's with the mean of the
+    two columns' g . n), linear in the heights through `Columns.film_operator`. Each sub-step:
+
+      1. predicts the heights at the end of the transport, H + dH, from the transport linearised on the columns: a column
+         sends across a face what its cells carry at their velocity without the thickness push plus their response
+         times the face's push -- implicitly in its moving metal, as the transport does -- linearised in the face's push
+         through the mobility of the cells that cross, sum r dz. Newton's method on that piecewise linear flux: the
+         system is solved again at the new face pushes until the flux the last solve missed would move no column by more
+         than THICKNESS_TOL of the largest change (1-3 solves on the test sphere and 3 on the 100 mm sphere's melting
+         state, measured 2026-10-21);
+      2. corrects every moving cell's velocity by its response times the change of the column's push -- `push`'s terms at
+         the predicted heights, less `applied` -- so that the cells' momentum follows the push of Task 5's form, and books
+         the impulse of that change and the anchored metal's drag on it (`MomentumBook.wall_response`);
+      3. hands `transport` each face's own push (ruling, 2026-10-21; the column form of Rhie-Chow interpolation): a
+         cell's velocity across a face becomes its velocity without the thickness push plus its response times that
+         face's push. The two-point difference sees a disturbance one column wide; the columns' pushes are centred
+         differences, blind to it, and with them alone a bump 1 um high on one column of 1.1 mm of slurry stayed at
+         99.8 % of its height while the faces' push took it to 12.5 % in one 10 ms sub-step.
+
+    Measured 2026-10-21. A cosine 24 mm long on 1 mm of slurry at 902 K pulled out of its surface at 100 m/s2 grows at
+    the thin-film rate h^3 k^2 (rho g - sigma k^2) / 3 mu within 0.42 % at sub-steps of 1-10 ms; one of 8 mm, below the
+    11.5 mm cutoff, decays at it within 1.2 % at 1 ms and within 0.26 % of backward Euler's own rate at 1-10 ms; a mode
+    alternating from column to column on 0.3 mm of liquid on 2 mm columns falls to at most 0.12 of itself in one sub-step
+    of 10 ms to 1 s. On the 100 mm sphere's melting state (the face where g . n > 0 melted 0.3 mm deep, 141 000 cells
+    moving) the step takes 0.30 s and the response 0.03-0.05 s more in the momentum step, against 0.31-0.37 s for the
+    transport. The growth reported is the largest dt M k^2 (rho_l g . n - sigma k^2) over the columns, M = sum r dz their
+    mobility and k at its fastest up to the columns' pi / sqrt(A): backward Euler follows the Rayleigh-Taylor growth only
+    where it is well under one."""
+    n, K = sk.m.shape
+    r = book.response
+    dz, A = sk.dz(), sk.area
+    H = dz.sum(axis=1)
+    flowing = sk.valid & np.isfinite(cell_viscosity(sk, 0.0))
+    moving = sk.valid & ~anchored_cells(sk, flowing)
+    wet = flowing_tops(sk)
+    normals = columns.normal
+    g = np.broadcast_to(np.asarray(g, dtype=float), (n, 3))
+    g_n = np.einsum("ij,ij->i", g, normals)
+    i, j = columns.edges[:, 0], columns.edges[:, 1]
+    w, g_e, both = columns.conduct_weight, 0.5 * (g_n[i] + g_n[j]), (wet[i] & wet[j]).astype(float)
+    kappa, J = columns.film_curvature(H), columns.film_operator(H)
+    phi0 = w * (sigma * both * (kappa[i] - kappa[j]) - rho_l * g_e * (H[i] - H[j]))      # each face's, first to second
+    Bt = columns.incidence.T.tocsr()
+    dphi = sp.diags(w) @ (sp.diags(sigma * both) @ (Bt @ J) - sp.diags(rho_l * g_e) @ Bt)  # its change per change of H
+    # every face each way whose sender moves; its cells' velocity across it without the thickness push [m2/s]
+    a, b, f = np.r_[i, j], np.r_[j, i], np.r_[columns.face_vector, -columns.face_vector]
+    keep = np.flatnonzero(moving.any(axis=1)[a])
+    a, b, f, phi_d = a[keep], b[keep], f[keep], np.r_[phi0, -phi0][keep]
+    dphi_d = sp.vstack([dphi, -dphi]).tocsr()[keep]
+    un = np.where(moving[a], np.einsum("dkx,dx->dk", (sk.U - r[..., None] * applied[:, None, :])[a], f), 0.0)
+    ra, da = r[a], np.where(moving[a], dz[a], 0.0)
+    flux = lambda phi: (np.maximum(un + ra * phi[:, None], 0.0) * da).sum(axis=1)       # m3/s
+    thick = np.where(moving, dz, 0.0).sum(axis=1)
+    D = len(a)
+    C = sp.csr_matrix((np.r_[np.ones(D), -np.ones(D)], (np.r_[a, b], np.r_[np.arange(D), np.arange(D)])), shape=(n, D))
+    active = np.zeros(n, dtype=bool)
+    active[a], active[b] = True, True
+    act = np.flatnonzero(active)
+    dH, phi, passes = np.zeros(n), phi_d, 0
+    F = flux(phi)
+    for passes in range(1, (THICKNESS_PASSES if D else 0) + 1):
+        M = (ra * da * (un + ra * phi[:, None] > 0.0)).sum(axis=1)                       # m5/(N s): the crossing cells'
+        Pa = sp.csr_matrix((F / np.where(thick[a] > 0.0, thick[a], 1.0), (np.arange(D), a)), shape=(D, n))
+        N = (sp.diags(A) + dt * C @ (Pa + sp.diags(M) @ dphi_d)).tocsr()[act][:, act]
+        dH = np.zeros(n)
+        dH[act] = spla.splu(N.tocsc()).solve(-dt * (C @ (F + M * (phi_d - phi)))[act])
+        new = phi_d + dphi_d @ dH
+        F_new = flux(new)
+        missed = dt * np.abs(C @ (F_new - F - M * (new - phi))) / A
+        phi, F = new, F_new
+        if missed.max() <= THICKNESS_TOL * np.abs(dH).max():
+            break
+    # the cells' velocities by the change of the columns' push; the faces' own push for the transport
+    G = rho_l * g_n[:, None] * columns.gradient(H + dH) + _capillary(columns, kappa + J @ dH, wet, sigma)
+    dG = G - applied
+    sk.U = np.where(moving[..., None], sk.U + r[..., None] * dG[:, None, :], sk.U)
+    across = np.zeros(2 * len(i))
+    across[keep] = phi - np.einsum("dx,dx->d", G[a], f)
+    # the books, and the fastest Rayleigh-Taylor growth in a sub-step
+    V = np.where(moving, dz, 0.0).sum(axis=1) * A
+    s_tension = sigma * wet
+    k2 = np.pi ** 2 / A
+    fastest = np.maximum(rho_l * g_n, 0.0) / (2.0 * np.maximum(s_tension, 1e-300))         # k^2 of the fastest growth
+    k2 = np.where(s_tension > 0.0, np.minimum(k2, fastest), k2)
+    s = (r * np.where(moving, dz, 0.0)).sum(axis=1) * k2 * (rho_l * g_n - s_tension * k2)
+    return ThicknessBook(push=dt * V[:, None] * dG, wall=-dt * book.wall_response[:, None] * dG, heights=H + dH,
+                         faces=FacePush(r, across), growth=float(max(dt * s.max(), 0.0)) if n else 0.0, passes=passes)
+
+
 # -- the base (Task 7) -------------------------------------------------------------------------------------------------
 def hand_down(sk, m_min, T_rigid):
     """Rigid cells at the bottom of each column (at or below T_rigid, contiguous from the bottom) pass to the 3D model,
```

- [ ] **Step 6: Run them to see them pass.** Expected: `tests/test_reentry_model_columns.py` 13 passed and
  `tests/test_reentry_model_skin.py` 67 passed; the unit tier 464 passed, 1 skipped. Nothing calls `thickness_step`
  yet: the flow joins the skin body's sub-steps in 6c, so no flight's result changes.

- [ ] **Step 7: Commit.**

```bash
git add reentry_model/columns.py tests/test_reentry_model_columns.py reentry_model/skin.py tests/test_reentry_model_skin.py
git commit -m "Skins: the push's thickness terms implicit in the column heights at the end of the transport -- surface tension in its thin-film form, only where a column's top flows, each face's own push in the transport and the cells' velocities corrected by their response -- and each floating crust one unknown of the momentum step (sub-plan 21 Task 6b)"
```

---

### Task 6c: The flow in the skin body

**Files:**
- Modify: `reentry_model/skin.py` (`FLOW_CYCLES`, `FLOW_SECONDS`, `Flow`, `FlowBook`, `flow_substeps`;
  `SkinBook.momentum_out` and `.flow`; `advance(..., flow=None)`), `reentry_model/skinbody.py` (`FLOW_DEFAULT`,
  `FlowForcing`; `SkinBody(..., flow, flow_cycles, surface_flow)`, `flow_forcing`, `advance(..., forcing=None)`, the
  momentum books and `momentum_residual`, the flow's history columns), `reentry_model/coupled.py` (`SKINBODY_COLUMNS`,
  `skin_results`, the frames' `top_speed`), `reentry_model/cli.py` (`--skin-flow`, `--skin-flow-cycles`,
  `skin_name_suffix(..., flow, flow_cycles)`, their checks, `build_thermal`)
- Test: `tests/test_reentry_model_skin.py`, `tests/test_reentry_model_skinbody.py`, `tests/test_reentry_model_cli.py`
  (the skin run's settings gain the flow's two keys)

**Interfaces:**
- Consumes: Task 5's `push` and `along_meridian`, `SurfaceFlowResult.dp_ds` and `.tau`; Task 4's `momentum_step`; 6a's
  `transport`; 6b's `flowing_tops`, `momentum_step(..., response=True)`, `thickness_step` and `FacePush`;
  `Skins.momentum_out`; `SurfaceFlow.evaluate`.
- Produces: `skin.Flow(columns, dp_ds, shear, g, v_hat, rho_l, sigma, A_top, cycles=FLOW_CYCLES)`, what drives the flow
  over a part (dp_ds (n,) [Pa/m], shear (n, 3) [Pa], g (3,) [m/s2], A_top (n,) [m2]); `skin.flow_substeps(sk, flow,
  dt) -> FlowBook`; `FlowBook(push, shear, wall, carried, courant, growth, passes, sub_steps)`, the impulses (n, 3)
  [N s] and the carried momentum [kg m/s], with `.zero(n)` and `.add(other)`; `skin.FLOW_CYCLES` (10) and the
  diagnostic `skin.FLOW_SECONDS`; `skin.advance(..., flow=None)`, its book gaining `momentum_out` (n, 3) [kg m/s], what
  the removed liquid took, and `flow`, the sub-steps' FlowBook; `skinbody.FlowForcing(dp_ds, tau, g)` with
  `.toward(other, frac)`; `SkinBody(..., flow=FLOW_DEFAULT, flow_cycles=None, surface_flow=None)` with `.flow_on`,
  `.flow_cycles`, `.flow_forcing(state)`, `.advance(t, dt, loads, state=None, forcing=None)`, `.momentum_residual()`,
  `.impulse`, `.momentum_lost`, `.impulse_size`, `.last_flow`; `skinbody.FLOW_DEFAULT` (False); history columns
  `momentum_residual`, `speed_max_ms`, `roughness_max`, `flow_courant_max`, `flow_growth_max`, `flow_passes_max`,
  `wall_flow_s`; results `momentum_residual`, `flow_courant_max`; the frames' `top_speed`; `run --skin-flow on|off`
  (run-name suffix `_skinflow-on`) and `--skin-flow-cycles <n>` (`_flowcycles<n>` when not FLOW_CYCLES).

The flow joins the skin body. Every macro step, the surface flow is evaluated on the corners where the heating is —
from the step's aero state, the top surface's angles to the flight direction, the nose radius and the tops'
temperatures at the step's start — for its pressure gradient along the meridian (Task 5's `dp_ds`) and its gas shear
`tau`; with the deceleration, which the metal feels as g, the drag's deceleration along the flight direction (minus the
body's own acceleration, Task 5), that is the step's `FlowForcing`. Each part takes it interpolated linearly in time
between the last macro step's and this one's at the part's middle, as it takes the heating (spec 7.1), the shear laid
along the meridian in each column's plane. Every sub-step of both passes then runs, after its conduction, `cycles` flow
sub-steps of its length over `cycles`, each in 6b's order — `push` with surface tension where a top flows,
`momentum_step` with its response, `thickness_step`, `transport` with the faces' push and no regrid — the cells split
and merged between them; the removal (with `remove`), the draw and the sub-step's own regrid follow the last (spec 7.1:
heat, flow, sideways exchange and spray). The trial pass runs the flow as the real pass does (spec 7.2's sub-steps are
spec 7.1's), the sensitivity carried by the transport. Where no cell flows nothing is computed: a velocity left in rigid
cells is booked to the wall and cleared, as the momentum step would.

**The momentum books.** The flow's impulses — the momentum and thickness steps' push, the gas shear, the walls (the
anchored metal's drag, the momentum cells held on joining it, the transport's projections) and what the transport
carried, which sums to zero over the body — add up in each part's book; what leaves with removed liquid, handed-down
cells or a reconciliation at the bottom is booked as lost; metal drawn arrives at rest. `momentum_residual` is the
skins' momentum less its value at time zero, the impulses and the losses, over the impulses' total size: zero for exact
books.

**The flow's sub-step** (option c, Asha, 2026-10-21, chosen here from measurements;
`prototype/work-2026-10-18-column-flow/measure/task6c/`, 2026-10-22, the laptop under a load average of 5 to 140 from
other sessions, so the times are processor times and indicative). The 100 mm sphere's 9 417 columns of 10 µm cells, its
front hemisphere's top 0.3 mm liquid at 950 K over rigid metal at 850 K (6b's destabilising face), at rest at t = 0 and
driven by the surface flow of the 100 mm flight at 93 s — a stagnation pressure of 11.9 kPa, -dp/ds up to
3.4·10⁵ Pa/m at 45°, gas shear up to 87 Pa at 49°, 29.8 m/s² of deceleration — with surface tension, for 50 ms of
flow alone (no conduction), at flow sub-steps of 10 ms to 0.25 ms:

| Flow sub-step | Courant, largest (median) | KE the transport destroys / KE made | KE at 50 ms / 0.25 ms run | Thickness change, L1 error / 0.25 ms run | Volume past 105° / 0.25 ms run | Linearisations, most |
|---|---|---|---|---|---|---|
| 10 ms | 1 266 (244) | 0.59 | 4.81 | 95 % | 0 | 8 (the cap) |
| 5 ms | 26 (23) | 0.62 | 0.31 | 85 % | 0.11 | 7 |
| 2.5 ms | 12.6 (11) | 0.52 | 0.36 | 59 % | 0.89 | 5 |
| 2 ms | 9.9 (8.4) | 0.43 | 0.48 | 44 % | 0.96 | 3 |
| 1 ms | 4.8 (4.1) | 0.27 | 1.04 | 3.9 % | 1.007 | 4 |
| 0.5 ms | 2.4 (2.0) | 0.28 | 1.02 | 1.4 % | 1.005 | 2 |
| 0.25 ms | 1.2 (1.0) | 0.27 | 1 | — | 1 | 2 |

The table shows the result converging at 1 ms and not above it. In 50 ms the liquid drains off the front at a few metres
a second and crosses the shoulder onto the rigid rear: the volume that crossed 85° is within 1.1 % of the 0.25 ms run's
(4 589 mm³) at 5 ms and below (3 % short at 10 ms), but the runoff over the rear is not. At 1, 0.5 and 0.25 ms it
reaches the rear stagnation point, every band of 15° from 90° to 180° within 7 % (1 ms) and 4 % (0.5 ms) of the finest
run's; at 2, 2.5, 5 and 10 ms it stops at 150°, 135°, 120° and 105°. The cause is 6a's arrival rule, not the scheme's
smearing: metal landing on a rigid top becomes a new top cell after the solve, so it flows on only in the next flow
sub-step, and a front crossing rigid metal advances at most one column a flow sub-step — 0.2 m/s on 2 mm columns at 10
ms, 2 m/s at 1 ms. That predicts the fronts seen: 25, 20, 10 and 5 columns of 2 mm in 50 ms are 57°, 46°, 23° and 11° of
the 50 mm radius beyond the shoulder at 2, 2.5, 5 and 10 ms. Inside the flowing layer the scheme's own smearing at 1 ms
is already at its floor: the transport destroys 27 % of the kinetic energy the push and the shear make, the same at 0.5
and 0.25 ms, so that is the 2 mm grid's upwind diffusion and the projections between planes, not the time step's. At 10
ms the thickness step's linearisations stop at their cap of 8 with the predicted heights 13 times the change off, and
its correction adds kinetic energy: the fastest cell is at 22 m/s at 30 ms against 4.9 at 1 ms. With every column's top
melted (no front; 9 417 columns, 282 000 cells moving), 10 ms flow sub-steps blow up — the fastest cell 72 m/s at the
third sub-step, 2 563 m/s at the fourth, the predicted heights 16 and 91 times the change off — and 5 ms keeps the
linearisations at their cap with the predictions about one change off; at 1 ms it is stable (Courant numbers to 4.7, two
linearisations, the predictions within 0.24 of the change over 10 ms), but its transport took 35-114 s a sub-step, the
flows meeting at the rear stagnation point joining 18 000-23 000 cells into one cycle (below).

*Rayleigh-Taylor growth* (`ThicknessBook.growth`, dt s, one flow sub-step from rest on the sphere's front melted 0.3,
1.0 and 3.4 mm deep, the flowing depths of the element model's flight at about 85, 105 and 120 s, under the forcing at
93 s, 29.8 m/s², and 115 s, 48.4 m/s²):

| Liquid depth, deceleration | 10 ms | 5 ms | 2 ms | 1 ms |
|---|---|---|---|---|
| 0.3 mm, 29.8 m/s² | 0.018 | 0.0048 | 8.3·10⁻⁴ | 2.1·10⁻⁴ |
| 0.3 mm, 48.4 m/s² | 0.047 | 0.013 | 0.0022 | 5.6·10⁻⁴ |
| 1.0 mm, 29.8 m/s² | 0.072 | 0.018 | 0.0030 | 7.6·10⁻⁴ |
| 1.0 mm, 48.4 m/s² | 0.19 | 0.049 | 0.0079 | 0.0020 |
| 3.4 mm, 29.8 m/s² | 0.25 | 0.063 | 0.010 | 0.0025 |
| 3.4 mm, 48.4 m/s² | 0.67 | 0.17 | 0.027 | 0.0067 |

The table shows the growth falling as the square of the sub-step — the layer's response in a sub-step is limited by
its inertia — and well under one at 5 ms and below in every case; at 10 ms the deep layers late in flight reach 0.67,
where backward Euler's growth in a sub-step, 1/(1 - dt s) = 3.0, is half again the exact e^(dt s) = 1.95.

**Ruling: ten flow sub-steps a heat sub-step (`FLOW_CYCLES` = 10, 1 ms under the 10 ms heat sub-step)** — at 1 ms the
measured state is within 3.9 % of the 0.25 ms run in where the metal went and 4 % in its kinetic energy, the runoff
reaches where the finest run's does, the growth is at most 0.0067 even for 3.4 mm of liquid at 48 m/s², the Courant
numbers are at most 4.8 and the kinetic energy the transport destroys is the grid's own 27 %; 2 ms holds the runoff back
by up to 30° of the rear and loses half the kinetic energy by 50 ms, and 10 ms is unstable in a fully melted layer.
Accepted: Courant numbers to about 5, where the implicit upwind scheme's time smearing no longer shows against the
grid's, and the grid's 27 % (the 2 mm surface's; only finer columns reduce it) — if wrong, the flow costs ten times what
one sub-step would: a 10 ms heat sub-step of the heated melting state takes 48 s with ten against 8 s with one and 4.5 s
for its conduction alone (*Cost*, below). The sub-step is fixed, not adapted: the Courant numbers grow with the layer's
speed, which grows with its depth and the push late in flight, and a front faster than 2 m/s on 2 mm columns is held
back at 1 ms too (Review Focus 13 and 14).

**Ruling: the flow is off by default (`FLOW_DEFAULT`, `--skin-flow off`), and with it on the liquid stays** — stage 2's
skins remove the liquid as it forms "instead of flow and spray" (spec 11.2); with the flow on and no spraying until
stage 4, the liquid stays and flows (`remove` False), and a flight would only pile it up — if wrong, a stage-3 flight
must ask for the flow. Flow off is the heat-only skins: the skin body's tests pass unchanged, and the 3 MW/m² removal
run of `test_heating_with_removal_keeps_mass_and_energy_exact` gives Task 6b's code's every array to the last bit
(2026-10-22, numpy's generator seeded as a run seeds it; unseeded, pyamg's starting vectors make any two runs differ in
round-off).

**Ruling: on the command line the flow takes Girin's removal** — `--skin-flow on` needs `--surface-model skin` and
`--melt on`, and `--removal girin` (the default): the liquid stays, and sub-plan 22's spraying will remove it on the
same command line; `--removal instant` is refused with the flow on, being the heat-only skins' stand-in for flow and
spray (spec 11.2), as `--removal girin` is refused with it off. `--skin-flow-cycles` needs the flow on and is an integer
of at least 1; the run name gains `_skinflow-on`, and `_flowcycles<n>` only for other than ten (spec 10.1's convention)
— if wrong, a stage-3 run with the flow on takes the girin removal's continuum step (`--dt-continuum`, 0.0125 s once
the body's Knudsen number is below 0.01) and its name says `_melt-girin` although nothing sprays until sub-plan 22;
`--dt-continuum off` keeps the 0.5 s step.

**Ruling: the forcing interpolated in time as the heating is, the shear in each column's plane** — spec 7.1's
treatment of the heating, at no cost; the push's density is the skins' liquid density (`Skins.rho_l`, the solid's
with `--skin-swell off`, as the cells' masses have it) — if wrong, the push differs by the forcing's change over half a
macro step.

**Ruling: the trial pass runs the flow** (spec 7.1 and 7.2: a sub-step is heat, flow, sideways exchange and spray) — it
changes the base law: on the 100 mm sphere's heated melting state (`heated_state`: the front's skins at the conduction
profile of 3 MW/m², 0.30 mm of liquid over 0.55 mm of slurry, the base held at the profile's 874 K, the flow developed
over 20 ms under the 93 s forcing), one part of 0.0625 s with the flow and without it gives base laws that differ, at
the base temperature, by 6.0 % of the energy through the columns' tops and bases in total; 3 428 of the 9 417 columns
differ by more than the interface tolerance of 10⁻³ of their own (the median of the flowing columns 7·10⁻⁴), and the
slope by up to 46 %. A trial pass without the flow would leave the real pass's base heat that far from its law: repeated
by the damped repetition, and what it leaves deferred to the nodes. The real pass's base heat meets the law made with
the flow to 9·10⁻¹⁵. It costs what the flow costs: the trial pass took 1 216 s of processor time with ten flow sub-steps
a heat sub-step against 17 s without the flow. If wrong, half the flow's cost is spent on the law.

**Ruling: the momentum residual over the impulses' total size**, the push's, the shear's and the walls' summed in size
over every column and part since time zero, as the energy residual is over the heat absorbed — if wrong, only the
scale of a number that should be round-off.

*Cost* (processor time, on `heated_state` with 310 000 cells flowing, the flow developed over 20 ms; each measured twice,
under other sessions' load, so indicative). A 10 ms heat sub-step:

| Flow sub-steps | Whole | Conduction | Push | Momentum | Thickness | Transport | Regrid | Rest (draw, regrid) |
|---|---|---|---|---|---|---|---|---|
| 1 (10 ms) | 7.9-8.5 s | 4.1-4.6 | 0.05 | 1.0 | 1.0 | 1.4 | — | 0.4 |
| 10 (1 ms) | 47.8-48.3 s | 4.5 | 0.5 | 9.4 | 3.2-3.5 | 28-29 | 1.4-1.5 | 0.3 |

The table shows the transport taking three fifths of a heat sub-step at ten flow sub-steps, 2.9 s each, against 1.4 s
for the one 10 ms sub-step, whose Courant number of 490 and eight linearisations (the cap) mark it unresolved; the
momentum step's shear-rate fixed point runs its five passes every sub-step, the slurry where liquid meets rigid metal. A
macro step while melting is eight parts of seven sub-steps, each run by the trial and the real pass, 1.8 times on
average (fact 126): about 200 heat sub-steps, which on this state is about 2.7 h with ten flow sub-steps, 27 minutes
with one and 16 minutes with the flow off — the flow, not the heat, is now the skins' cost. (Fact 126 measured a melting
skin pass at 0.65 s a 10 ms sub-step; these skins are deeper, 1.6 mm, and the laptop was loaded, so the ratios are the
finding, not the seconds.) The remedy for the transport's large cycles below would take its share down; the momentum
step's solve over every column, flowing or not, is the next.

*A quiescent melted region* (Review Focus 12). The sphere's melting state at rest under surface tension alone, no gas
forcing and no deceleration: the melt's rim at the shoulder rounds at 5-17 mm/s (Courant numbers 0.01-0.2), and the
transport's cycles of meeting flows grow from 7 cells to 10 315 of the 141 000 moving within 12 flow sub-steps of 1 ms
(to 7 022 within 6 of 10 ms), its time from 1.1 to 4.0-4.8 s a sub-step. Profiled, 2.3 s of a 3.35 s transport is
SuperLU's factorisation of its two systems, whose lower triangles fill 5.3 and 5.8 times. Under the 93 s forcing the
largest cycle was 1 700 cells after 20 sub-steps of 1 ms and the factorisations took 0.06 and 0.15 s, filling 1.1 times;
with every column melted, flows meeting at the rear stagnation point joined 57 000 cells into one cycle within 9
sub-steps of 5 ms, and a transport took 651 s. *A remedy measured, not built:* solve each system by the iteration x ← x
+ L⁻¹ (b − A x) started from the upwind solve, L its lower triangle in upwind order — the whole system but the cycles'
back-couplings, which an upwind order makes triangular with no fill. On the quiescent state it took 2-3 corrections to
reach 10⁻¹⁰-10⁻¹³ of the right-hand side, 0.05 s against the LU's 0.30 s for the mass and 0.16 s against 1.23 s for the
momentum, the masses positive and the books missing 2·10⁻¹⁷ of the total at 10⁻¹³; under the forcing (small cycles) 4-5
corrections, 0.10 s against the LU's 0.06 s; with every column melted (a cycle of 18 379 cells) 7-9 corrections, 0.19 s
against 7.4 s and 0.41 s against 41.8 s. Installing each cell's mass from the iterate's fluxes, m + dt (arrivals −
departures), would close the books exactly by construction. It changes 6a's solver, so it is reported, not built (Review
Focus 12).

**Measured against the tests** (2026-10-22): one flow sub-step is the push, the momentum step, the thickness step and the
transport called by hand, to the last bit, at one cycle and three, the books to 2.1·10⁻¹² and 9.4·10⁻¹³ of the largest
impulse (6b's on these skins); four flow sub-steps in each 10 ms heat sub-step give sixteen heat sub-steps of 2.5 ms
with one each, to the last bit, with the conduction off; a pass with the conduction on, removal and the flow closes its
momentum to 4.9·10⁻¹⁶; on the 4 mm test sphere under 3 MW/m² from 880 K for two 0.25 s steps with the flow on the mass
is kept to 1.1·10⁻¹⁶, the energy to 9.8·10⁻¹³ and the momentum to 2.8·10⁻¹⁷, the fastest cell 1.6 m/s, Courant numbers
to 6.9.

- [ ] **Step 1: The flow sub-steps' tests** — append to `tests/test_reentry_model_skin.py`:

```diff
--- a/tests/test_reentry_model_skin.py
+++ b/tests/test_reentry_model_skin.py
@@ -1388,3 +1388,171 @@
         worst = max(worst, np.abs(tb.heights - sk.thickness()).max() / np.abs(sk.thickness() - H0).max())
         passes.append(tb.passes)
     assert worst < 0.0125 and max(passes) <= 4 and np.std(sk.thickness()) < 0.5 * spread0
+
+
+# -- the flow in the sub-steps (sub-plan 21 Task 6c) --------------------------------------------------------------------
+def sphere_flow(col, cycles=1, tau=50.0, A_top=None):
+    """A Flow on the test sphere's columns: a Newtonian pressure falling away from the stagnation point (11 kPa of
+    stagnation pressure), `tau` Pa of gas shear along the meridian, 30 m/s2 of deceleration along V_HAT, surface
+    tension."""
+    cos = np.clip(col.normal @ V_HAT, 0.0, 1.0)
+    dp_ds = -2.0 * 11e3 * np.sqrt(1.0 - cos ** 2) * cos / 0.05
+    shear = skin.along_meridian(col.normal, V_HAT, np.where(cos > 0.0, tau, 0.0))
+    return skin.Flow(col, dp_ds, shear, 30.0 * V_HAT, V_HAT, RHO_L, SIGMA, col.area if A_top is None else A_top, cycles)
+
+
+def by_hand(sk, flow, dt):
+    """One flow sub-step called piece by piece, in 6b's order."""
+    G, terms = skin.push(flow.columns, sk.thickness(), flow.dp_ds, flow.g, flow.v_hat, flow.rho_l, flow.sigma,
+                         wet=skin.flowing_tops(sk))
+    mb = skin.momentum_step(sk, dt, G, flow.shear, flow.A_top, response=True)
+    tb = skin.thickness_step(sk, flow.columns, dt, mb, terms["thickness"] + terms["capillary"], flow.g, flow.rho_l,
+                             flow.sigma)
+    return mb, tb, skin.transport(sk, flow.columns, dt, regrid=False, faces=tb.faces)
+
+
+def test_a_flow_sub_step_is_the_push_the_momentum_the_thickness_step_and_the_transport_in_that_order(coarse_sphere_mesh):
+    """6b's ruling: every flow sub-step calls `push` with surface tension where a column's top flows, `momentum_step`
+    with its response, `thickness_step` and `transport` with the faces' push, in that order; `flow_substeps` does exactly
+    that, to the last bit, and books every impulse: the columns' momentum changes by the push's (the thickness step's
+    included), the gas shear's, the walls' (the anchored metal's, the thickness step's and the transport's projections)
+    and what the transport carried, to round-off (2.1e-12 of the largest impulse at one cycle and 9.4e-13 at three,
+    measured 2026-10-22, as 6b's books on these skins). With `cycles` the sub-step runs as that many flow sub-steps,
+    the cells split and merged between them; the heat sub-step's own regrid follows the last."""
+    from test_reentry_model_columns import sphere_columns
+    col, _ = sphere_columns(coarse_sphere_mesh)
+    for cycles in (1, 3):
+        first, _, _ = mixed_sphere_skins(coarse_sphere_mesh)
+        second, _, _ = mixed_sphere_skins(coarse_sphere_mesh)
+        flow = sphere_flow(col, cycles, A_top=col.area * 1.01)
+        P0 = first.column_momentum()
+        book = skin.flow_substeps(first, flow, 0.01)
+        parts = []
+        for k in range(cycles):
+            if k:
+                second.regrid()
+            parts.append(by_hand(second, flow, 0.01 / cycles))
+        assert all(np.array_equal(x, y) for x, y in zip(first.state(), second.state()))
+        assert book.sub_steps == cycles and book.courant == max(tr.courant for _, _, tr in parts) > 0.0
+        assert book.growth == max(tb.growth for _, tb, _ in parts) and book.passes == max(tb.passes for _, tb, _ in parts)
+        total = lambda f: sum(f(*p) for p in parts)
+        np.testing.assert_allclose(book.push, total(lambda mb, tb, tr: mb.push + tb.push), rtol=1e-14, atol=0.0)
+        np.testing.assert_allclose(book.shear, total(lambda mb, tb, tr: mb.shear), rtol=1e-14, atol=0.0)
+        np.testing.assert_allclose(book.wall, total(lambda mb, tb, tr: mb.wall + tb.wall + tr.wall), rtol=1e-14, atol=1e-30)
+        np.testing.assert_allclose(book.carried, total(lambda mb, tb, tr: tr.momentum), rtol=1e-14, atol=1e-30)
+        impulse = book.push + book.shear + book.wall + book.carried
+        scale = max(np.abs(x).max() for x in (book.push, book.shear, book.wall, P0))
+        assert np.abs(first.column_momentum() - P0 - impulse).max() < 5e-12 * scale
+
+
+def test_where_nothing_flows_the_flow_changes_no_cell_and_books_what_the_cells_held(coarse_sphere_mesh):
+    """Skins all below T_rigid: nothing flows, so the flow sub-step changes no mass, enthalpy or temperature, to the
+    last bit; any velocity left in the rigid cells is booked to the wall, as `momentum_step` books an anchored cell's."""
+    from test_reentry_model_columns import sphere_columns
+    col, _ = sphere_columns(coarse_sphere_mesh)
+    sk, _, rng = random_sphere_skins(coarse_sphere_mesh)
+    sk.T = np.minimum(sk.T, 850.0)
+    sk.H = sk.m * sk.material.enthalpy(sk.T)
+    sk.U = np.where(sk.valid[..., None], rng.normal(0.0, 0.01, sk.U.shape), 0.0)
+    before, P0 = sk.state(), sk.column_momentum()
+    book = skin.flow_substeps(sk, sphere_flow(col, 2), 0.01)
+    assert all(np.array_equal(x, y) for x, y in zip(before[:5], sk.state()[:5])) and not sk.U.any()
+    np.testing.assert_array_equal(book.wall, -P0)
+    assert not (book.push.any() or book.shear.any() or book.carried.any()) and book.sub_steps == 0
+
+
+def strip_layer(nx=24, h=1e-3, liquid=20, bump=0.0, T=LIQUID):
+    """Flat columns (`strip_grid`, two corners wide), each `liquid` cells of liquid aluminium at T over five of rigid
+    metal, the liquid's mass raised by `bump` (a fraction) on the columns 4-7 mm along x."""
+    col, pts = strip_grid(nx, 2, h)
+    sk = stacked(col, [[T] * liquid + [RIGID] * 5] * col.n)
+    if bump:
+        f = 1.0 + bump * ((pts[:, 0] > 3.5e-3) & (pts[:, 0] < 7.5e-3))
+        m, H = sk.m.copy(), sk.H.copy()
+        m[:, :liquid] *= f[:, None]
+        H[:, :liquid] *= f[:, None]
+        sk._install(m, H, sk.count, T=sk.T)
+    return sk, col, pts
+
+
+def strip_flow(col, dp_ds=-2e4, g=np.zeros(3), sigma=SIGMA, cycles=1, tau=0.0):
+    """A Flow on flat columns (normals +z), the flight direction tilted from -z toward -x so that the meridian runs
+    along +x: the pressure falling at `dp_ds` along it pushes the metal toward +x, the gas shear tau too."""
+    v_hat = np.array([-0.6, 0.0, -0.8])
+    shear = skin.along_meridian(col.normal, v_hat, np.full(col.n, tau))
+    return skin.Flow(col, np.full(col.n, dp_ds), shear, g, v_hat, RHO_L, sigma, col.area, cycles)
+
+
+def no_conduction(sk, dt, q_top, A_top, T_amb, eps, T_b, sideways=None):
+    z = np.zeros(sk.n)
+    return (z, z, z, z) + ((z,) if sk.S is not None else ())
+
+
+def test_n_flow_sub_steps_in_a_heat_sub_step_are_n_heat_sub_steps_of_their_length_with_the_conduction_off(monkeypatch):
+    """The flow sub-cycles within a heat sub-step (Asha, 2026-10-21, option c): with the conduction off, four flow
+    sub-steps in each 10 ms heat sub-step give exactly what sixteen heat sub-steps of 2.5 ms with one each give, to the
+    last bit -- a bump on liquid aluminium on a strip of 1 mm columns, pushed along it, pulled out of its surface and
+    smoothed by surface tension; the books agree to round-off. With the conduction on the two differ by the splitting
+    between heat and flow, of first order in the heat sub-step."""
+    monkeypatch.setattr(skin, "_conduct", no_conduction)
+    out = []
+    for substep, cycles in ((0.01, 4), (0.0025, 1)):
+        sk, col, _ = strip_layer(bump=0.2)
+        sk.settings = skin.SkinSettings(**CELLS_04, substep=substep)
+        T_b = np.full(sk.n, RIGID)
+        book = skin.advance(sk, 0.04, np.zeros(sk.n), sk.area, 0.0, 0.0, T_b, np.zeros(sk.n), T_b1=T_b, remove=False,
+                            flow=strip_flow(col, g=np.array([0.0, 0.0, 30.0]), cycles=cycles))
+        out.append((sk, book))
+    (a, A), (b, B) = out
+    assert all(np.array_equal(x, y) for x, y in zip(a.state(), b.state()))
+    assert A.flow.sub_steps == B.flow.sub_steps == 16 and A.flow.courant == B.flow.courant > 0.0
+    for key in ("push", "shear", "wall", "carried"):
+        x, y = getattr(A.flow, key), getattr(B.flow, key)
+        np.testing.assert_allclose(x, y, rtol=0.0, atol=1e-13 * np.abs(A.flow.push).max())
+
+
+def test_a_uniform_push_moves_a_flat_layer_as_the_momentum_step_and_the_transport_alone(monkeypatch):
+    """With neither surface tension nor a pull out of the surface, a heat sub-step's flow is the momentum step and the
+    transport alone (6b), and so is `advance`'s with the conduction off: a uniform push of 20 kPa/m along a strip of
+    liquid aluminium 0.23 mm deep on 1 mm columns moves the metal along it as those two steps and the regrid do, to the
+    last bit; the mass is kept and the first columns drain into the last."""
+    monkeypatch.setattr(skin, "_conduct", no_conduction)
+    first, col, pts = strip_layer()
+    second, _, _ = strip_layer()
+    first.settings = second.settings = skin.SkinSettings(**CELLS_04, substep=0.01)
+    T_b, M0, H0 = np.full(first.n, RIGID), first.mass(), first.thickness()
+    flow = strip_flow(col, sigma=0.0)
+    skin.advance(first, 0.04, np.zeros(first.n), first.area, 0.0, 0.0, T_b, np.zeros(first.n), T_b1=T_b, remove=False,
+                 flow=flow)
+    G, _ = skin.push(col, second.thickness(), flow.dp_ds, flow.g, flow.v_hat, RHO_L, 0.0)
+    np.testing.assert_allclose(G, np.tile([2e4, 0.0, 0.0], (col.n, 1)), rtol=1e-12, atol=1e-9)
+    for _ in range(4):
+        skin.momentum_step(second, 0.01, G, flow.shear, col.area)
+        skin.transport(second, col, 0.01, regrid=False)
+        second.regrid()
+    assert all(np.array_equal(x, y) for x, y in zip(first.state(), second.state()))
+    x = pts[:, 0]
+    assert first.mass() == pytest.approx(M0, rel=1e-13)
+    assert first.thickness()[x < 1e-9].max() < H0.max() and first.thickness()[x > x.max() - 1e-9].min() > H0.max()
+
+
+def test_the_books_of_a_pass_with_flow_close_exactly():
+    """A pass of 0.1 s with the conduction on, liquid removed as it forms and the flow on: liquid aluminium at 909 K, just
+    below T_feed, on a strip heated at 30 MW/m2, pushed along it and sheared: the mass and the energy close as without
+    flow, and each column's momentum changes by the flow's impulses less what the removed liquid took with it, to
+    round-off (4.9e-16 of the largest impulse, measured 2026-10-22)."""
+    sk, col, _ = strip_layer(nx=12, T=909.0)
+    sk.settings = skin.SkinSettings(**CELLS_04, substep=0.01)
+    T_b = np.full(sk.n, RIGID)
+    m_min = skin.minimum_mass(sk.material, sk.area, sk.settings)
+    E0, M0, P0 = sk.energy(), sk.mass(), sk.column_momentum()
+    q = np.full(sk.n, 3e7)
+    book = skin.advance(sk, 0.1, q, sk.area, 300.0, 0.3, T_b, m_min, T_b1=T_b, flow=strip_flow(col, tau=50.0))
+    assert book.removed.sum() > 0.0 and np.abs(book.momentum_out).max() > 0.0 and book.flow.courant > 0.0
+    assert sk.mass() == pytest.approx(M0 - book.removed.sum() + book.drawn.sum(), rel=1e-13)
+    gained = book.E_top - book.E_rad - book.E_base + book.E_side - book.removed_H + book.drawn_H
+    assert sk.energy() == pytest.approx(E0 + gained.sum(), rel=1e-12)
+    f = book.flow
+    impulse = f.push + f.shear + f.wall + f.carried - book.momentum_out
+    scale = max(np.abs(x).max() for x in (f.push, f.shear, f.wall, book.momentum_out))
+    assert np.abs(sk.column_momentum() - P0 - impulse).max() < 1e-13 * scale
```

- [ ] **Step 2: The skin body's tests** — append to `tests/test_reentry_model_skinbody.py`:

```diff
--- a/tests/test_reentry_model_skinbody.py
+++ b/tests/test_reentry_model_skinbody.py
@@ -361,3 +361,107 @@
         b.advance(0.5 * (k + 1), 0.5, front_loads(b, 3e6))
         assert b.mass() == fresh()
     assert b.mesh.volume() != V0                                     # the corners moved
+
+
+# -- the flow in the skin body (sub-plan 21 Task 6c) -------------------------------------------------------------------
+def newtonian(b, p_stag=11e3, tau=50.0, decel=30.0):
+    """A prescribed gas-side forcing on the body's corners: a Newtonian pressure falling away from the stagnation point,
+    `tau` Pa of gas shear on the windward face and the deceleration along the flight direction."""
+    cos = np.clip(b.columns.normal @ b.v_hat, 0.0, 1.0)
+    return skinbody.FlowForcing(-2.0 * p_stag * np.sqrt(1.0 - cos ** 2) * cos / 0.05, np.where(cos > 0.0, tau, 0.0),
+                                decel * b.v_hat)
+
+
+def same_body(a, b):
+    """Whether two skin bodies hold the same state, to the last bit."""
+    return (np.array_equal(a.solver.temperature(), b.solver.temperature()) and np.array_equal(a.mesh.points, b.mesh.points)
+            and all(np.array_equal(x, y) for x, y in zip(a.skins.state(), b.skins.state()))
+            and np.array_equal(a.pending, b.pending))
+
+
+def test_the_flow_is_off_by_default_and_turned_on_where_nothing_flows_it_changes_nothing(coarse_sphere_mesh):
+    """The flow is off by default (the ruling of 2026-10-22: until stage 4's spraying gives the liquid a sink), and the
+    heat-only skins are then unchanged: their tests above pass as they did, and the 3 MW/m2 removal run of
+    `test_heating_with_removal_keeps_mass_and_energy_exact` reproduces Task 6b's code to the last bit, every array of
+    the 3D model and the skins (measured 2026-10-22, numpy's generator seeded as a run seeds it: pyamg draws its
+    starting vectors from it, so two unseeded runs differ in round-off). Turned on in a body where nothing flows -- 1
+    MW/m2 on a body at 600 K for two 0.5 s steps -- the flow is driven and its books kept but no cell moves, and the
+    body is the flow-off body to the last bit."""
+    bodies = []
+    for flow in (False, True):
+        np.random.seed(12345)
+        b = make(coarse_sphere_mesh, T0=600.0, flow=flow, remove=not flow)
+        for k in range(2):
+            b.advance(0.5 * (k + 1), 0.5, front_loads(b, 1e6), forcing=newtonian(b))
+        bodies.append(b)
+    off, on = bodies
+    assert off.flow_on is False and skinbody.FLOW_DEFAULT is False and on.flow_on is True
+    assert make(coarse_sphere_mesh, T0=600.0).flow_on is False
+    assert same_body(off, on) and np.all(on.skins.T[on.skins.valid] < on.material.T_rigid)
+    a, c = off.skin_stats(), on.skin_stats()
+    assert all(a[k] == c[k] for k in a if not k.startswith("wall_"))
+    assert c["speed_max_ms"] == 0.0 and c["momentum_residual"] == 0.0 and c["flow_courant_max"] == 0.0
+
+
+def test_a_skin_body_with_flow_keeps_its_mass_and_energy_closes_its_momentum_and_reports_the_flow(coarse_sphere_mesh):
+    """3 MW/m2 on the front of the 4 mm test sphere from 880 K for two 0.25 s steps, flow on with a Newtonian pressure
+    gradient, 50 Pa of gas shear and 30 m/s2 of deceleration, the liquid kept, one flow sub-step a heat sub-step (the
+    books, not the flow's accuracy, are tested; the default's ten cost ten times the flow): the mass is the starting
+    mass to 1e-13, the energy books close to 1e-8 of the heat absorbed, and the momentum balance -- the skins' momentum
+    against the impulses of the push, the shear and the walls, and what left with metal handed down -- closes to
+    round-off. The history's flow columns are the state's: the largest speed in any cell, the corner-to-corner roughness
+    of the top, the flow sub-steps' largest Courant number, Rayleigh-Taylor growth and linearisations. Measured
+    2026-10-22: mass to 1.1e-16, energy 9.8e-13, momentum 2.8e-17; the fastest cell 1.6 m/s, the flowing layer up to
+    1.38 mm deep, Courant numbers to 6.9, growth 2.1e-3, roughness 0.14."""
+    b = make(coarse_sphere_mesh, flow=True, remove=False, flow_cycles=1)
+    M0 = b.mass()
+    for k in range(2):
+        b.advance(0.25 * (k + 1), 0.25, front_loads(b, 3e6), forcing=newtonian(b))
+    stats = b.skin_stats()
+    assert b.removed_mass == 0.0 and b.mass() == pytest.approx(M0, rel=1e-13)
+    assert abs(b.energy_balance_residual()) < 1e-8
+    assert stats["momentum_residual"] < 1e-12 and b.impulse_size > 0.0
+    sk = b.skins
+    speed = np.linalg.norm(sk.U, axis=2)
+    assert stats["speed_max_ms"] == float(speed.max()) > 0.0
+    i, j = b.columns.edges[:, 0], b.columns.edges[:, 1]
+    H = sk.thickness()
+    assert stats["roughness_max"] == float((np.abs(H[i] - H[j]) / b.columns.edge_length).max()) > 0.0
+    assert stats["flow_courant_max"] > 0.0 and 0.0 <= stats["flow_growth_max"] < 1.0 and stats["flow_passes_max"] >= 1
+    assert stats["wall_flow_s"] > 0.0 and stats["flowing_depth_max_mm"] > 0.0
+
+
+def test_the_flow_is_driven_by_the_surface_flow_on_the_corners_taken_at_each_part_s_middle(coarse_sphere_mesh, monkeypatch):
+    """Every macro step the skin body evaluates the surface flow on its corners, as the heating is, from the step's aero
+    state and the tops' temperatures: the pressure gradient along the meridian, the gas shear along it on the top cells,
+    and the deceleration, which the metal feels along the flight direction. Each part takes them interpolated linearly
+    between the last macro step's and this one's at its middle, as it takes the heating."""
+    pytest.importorskip("cantera")
+    from test_reentry_model_surface_flow import aero_state_at
+    from reentry_model import surface_flow
+    b = make(coarse_sphere_mesh, T0=600.0, flow=True, remove=False)
+    seen = []
+    advance = skin.advance
+
+    def recorded(sk, dt, *a, **k):
+        seen.append((b.columns.normal.copy(), k["flow"]))              # the part's normals: the corners move after it
+        return advance(sk, dt, *a, **k)
+    monkeypatch.setattr(skinbody.skin_mod, "advance", recorded)
+    evaluate = lambda a: surface_flow.SurfaceFlow().evaluate(a, b.theta, b.nose_radius(), b.material.liquid.rho,
+                                                             b.surface_temperature())
+    states = (aero_state_at(90.0), aero_state_at(92.0))
+    first = evaluate(states[0])
+    b.advance(0.5, 0.5, front_loads(b, 1e6), state=states[0])
+    assert all(np.array_equal(f.dp_ds, first.dp_ds) for _, f in seen) and len(seen) >= 2
+    normal, f = seen[0]
+    np.testing.assert_allclose(f.shear, skin.along_meridian(normal, b.v_hat, first.tau), rtol=1e-12, atol=1e-12)
+    np.testing.assert_allclose(f.g, first.deceleration * b.v_hat, rtol=1e-12)
+    assert f.rho_l == b.skins.rho_l and f.sigma == b.material.liquid.sigma and f.cycles == skin.FLOW_CYCLES
+    second = evaluate(states[1])
+    seen.clear()
+    b._hint = 2                                                          # two parts, to see the interpolation
+    b.advance(1.0, 0.5, front_loads(b, 1e6), state=states[1])
+    assert b.step_parts == 2
+    for frac in (0.25, 0.75):                                           # the two parts' middles
+        dp = first.dp_ds + (second.dp_ds - first.dp_ds) * frac
+        assert any(np.allclose(f.dp_ds, dp, rtol=1e-12, atol=1e-9) for _, f in seen)
```

- [ ] **Step 3: The command line's tests** — in `tests/test_reentry_model_cli.py` (the skin run's settings gain the
  flow's two keys):

```diff
--- a/tests/test_reentry_model_cli.py
+++ b/tests/test_reentry_model_cli.py
@@ -379,7 +379,7 @@
     assert s["surface_model"] == "skin" and s["prism_layers"] == 0
     assert s["skin"] == {"min_thickness_mm": 0.8, "fine_zone_mm": 0.4, "cell_um": 20.0, "growth": 1.1, "substep_ms": 10.0, "preset": "dev",
                          "deform_depth_mm": 15.0, "n_columns": r["skin_columns"], "fine_band_mm": 4.0, "fine_size_mm": 2.0,
-                         "swell": True}
+                         "swell": True, "flow": "off", "flow_cycles": cli.skin.FLOW_CYCLES}
     assert r["skin_columns"] == s["n_patches"] and abs(r["energy_balance_residual"]) < 1e-8
     rows = list(csv.DictReader(open(path[:-5] + ".csv")))
     assert all(k in rows[0] for k in coupled.SKINBODY_COLUMNS) and float(rows[-1]["skin_cells"]) > 0.0
@@ -399,7 +399,7 @@
     MELT + ["--removal", "instant", "--surface-model", "skin", "--deform-depth", "0"],
     MELT + ["--removal", "instant", "--surface-model", "skin", "--fine-band", "-1"],
     MELT + ["--removal", "instant", "--surface-model", "skin", "--fine-size", "0"],
-    MELT + ["--surface-model", "skin"],                                     # Girin's removal needs flow and spraying
+    MELT + ["--surface-model", "skin"],                                     # Girin's removal needs the flow (and spraying)
     MELT + ["--removal", "instant", "--surface-model", "skin", "--animate"],
     MELT + ["--removal", "instant", "--surface-model", "skin", "--prism-layers", "4"],     # the rebuild cannot make prisms
     ["run", "--surface-model", "skin"],                                     # no thermal model
@@ -422,3 +422,55 @@
     assert "--prism-layers" in err and "rebuild cannot reproduce prism layers" in err
     given = cli.build_parser().parse_args(FEM + ["--melt", "on", "--removal", "instant", "--surface-model", "skin"])
     assert given.prism_layers is None                                           # not given: nothing to refuse
+
+
+# ---------------------------------------------------------------------------------------------------------------
+# Sub-plan 21 Task 6c: the flow in the skins, its setting and its sub-cycles
+
+def test_the_flow_settings_reach_the_run_name_only_when_not_default():
+    from reentry_model import skin
+    suffix = lambda **k: cli.skin_name_suffix(0.8, None, 1.10, 10.0, "production", 15.0, **k)
+    assert suffix() == suffix(flow=False) == ""
+    assert suffix(flow=True) == suffix(flow=True, flow_cycles=skin.FLOW_CYCLES) == "_skinflow-on"
+    assert suffix(flow=True, flow_cycles=skin.FLOW_CYCLES + 1) == "_skinflow-on_flowcycles{}".format(skin.FLOW_CYCLES + 1)
+    assert cli.skin_name_suffix(0.8, None, 1.10, 10.0, "dev", 15.0, swell=False, flow=True) == "_skin-dev_noswell_skinflow-on"
+    args = cli.build_parser().parse_args(["run", "--diameter", "100", "--velocity", "7.5", "--altitude", "77.5"])
+    assert args.skin_flow == "off" and args.skin_flow_cycles is None
+
+
+def test_a_skin_run_with_flow_keeps_its_liquid_and_writes_the_flow_s_columns_and_frames(tmp_path):
+    """0.5 s from 71 km with a warm body, the layered skins (dev preset) and the flow on at two flow sub-steps a skin
+    sub-step: the liquid stays and flows -- no removal, the spraying that will take it is stage 4's -- the run name ends
+    in _skinflow-on_flowcycles2 after the removal's _melt-girin, the settings record the flow and its sub-cycles, and the
+    history's flow columns and the frames' top speed are written."""
+    import glob
+    import pyvista as pv
+    argv = MELT + ["--altitude", "71", "--velocity", "7.24", "--temperature", "880", "--t-max", "0.5",
+                   "--surface-model", "skin", "--skin-preset", "dev", "--frames-every", "1", "--fine-size", "2",
+                   "--skin-flow", "on", "--skin-flow-cycles", "2", "--outdir", str(tmp_path)]
+    assert cli.main(argv) == 0
+    (path,) = glob.glob(str(tmp_path / "*.json"))
+    doc = json.load(open(path))
+    assert doc["run_name"].endswith("_melt-girin_surface-skin_skin-dev_fine4x2mm_skinflow-on_flowcycles2")
+    s, r = doc["settings"], doc["results"]
+    assert s["skin"]["flow"] == "on" and s["skin"]["flow_cycles"] == 2
+    assert r["removed_mass_kg"] == 0.0 and abs(r["energy_balance_residual"]) < 1e-8 and r["momentum_residual"] < 1e-12
+    rows = list(csv.DictReader(open(path[:-5] + ".csv")))
+    assert all(k in rows[0] for k in coupled.SKINBODY_COLUMNS)
+    surf = sorted(glob.glob(str(tmp_path / "**" / "surface_*.vtp"), recursive=True))
+    assert len(surf) == 2 and "top_speed" in pv.read(surf[-1]).point_data
+
+
+@pytest.mark.parametrize("argv", [
+    MELT + ["--surface-model", "skin", "--skin-flow", "maybe"],
+    MELT + ["--surface-model", "skin", "--skin-flow", "on", "--skin-flow-cycles", "0"],
+    MELT + ["--surface-model", "skin", "--skin-flow", "on", "--skin-flow-cycles", "2.5"],
+    MELT + ["--removal", "instant", "--surface-model", "skin", "--skin-flow-cycles", "4"],     # sub-cycles of no flow
+    MELT + ["--removal", "instant", "--surface-model", "skin", "--skin-flow", "on"],          # removal is instead of flow
+    MELT + ["--skin-flow", "on"],                                                             # the element model's runoff
+    FEM + ["--surface-model", "skin", "--skin-flow", "on"],                                   # nothing melts
+])
+def test_bad_flow_arguments_exit_2(argv, tmp_path):
+    with pytest.raises(SystemExit) as exc:
+        cli.main(argv + ["--outdir", str(tmp_path)])
+    assert exc.value.code == 2
```

- [ ] **Step 4: Run them to see them fail.** `"$PY" -m pytest tests/test_reentry_model_skin.py
  tests/test_reentry_model_skinbody.py tests/test_reentry_model_cli.py -q`. Expected: 11 failed, 143 passed —
  - `AttributeError: module 'reentry_model.skin' has no attribute 'Flow'` (4) and `... 'flow_substeps'` (1);
  - `TypeError: SkinBody.__init__() got an unexpected keyword argument 'flow'` (3);
  - `TypeError: skin_name_suffix() got an unexpected keyword argument 'flow'` (1);
  - `SystemExit: 2`, `unrecognized arguments: --skin-flow on --skin-flow-cycles 2` (1);
  - `AttributeError: module 'reentry_model.skin' has no attribute 'FLOW_CYCLES'`, the skin run's settings (1).

  The seven exit-code cases pass already, their flags being unknown to the parser.

- [ ] **Step 5: The flow sub-steps and the pass's books** in `reentry_model/skin.py`:

```diff
--- a/reentry_model/skin.py
+++ b/reentry_model/skin.py
@@ -28,6 +28,7 @@
 (`skinbody.SkinBody`) turns the drawn, absorbed and handed-down metal into the corners' displacements and reconciles
 each column with what the moving mesh actually swept."""
 import math
+import time
 from dataclasses import dataclass
 
 import numpy as np
@@ -495,6 +496,8 @@
     removed_H: np.ndarray    # the enthalpy it took
     drawn: np.ndarray        # rigid metal drawn from the 3D model in the sub-steps
     drawn_H: np.ndarray      # the enthalpy it brought
+    momentum_out: np.ndarray = None   # kg m/s (n, 3): what the removed liquid took with it (sub-plan 21 Task 6c)
+    flow: "FlowBook" = None           # the flow sub-steps' books, with the flow on (sub-plan 21 Task 6c)
 
 
 def _coupled_solve(G, diag, rhs, valid, sideways):
@@ -627,16 +630,17 @@
     return tuple(out)
 
 
-def advance(sk, dt, q_top, A_top, T_amb, eps, T_b0, m_min, T_b1=None, side=None, remove=True, g_base=None):
+def advance(sk, dt, q_top, A_top, T_amb, eps, T_b0, m_min, T_b1=None, side=None, remove=True, g_base=None, flow=None):
     """The macro step's sub-steps for every column (spec 7.2). Each sub-step: the implicit conduction, through the depth
     and, with `side`, sideways (`side(sk, previous)` makes the sub-step's `Sideways` map from the last one, Task 5),
-    then -- with `remove` and a
-    melting material -- the contiguous top
-    cells at or above T_feed leave, a column below its minimum mass m_min [kg] draws back to it, and the cells are split
-    and merged. The draw takes rigid metal (spec 6.6): at the base temperature, but never above T_rigid -- a 3D base
-    node can overshoot it within a macro step (on 4 mm elements under 6 MW/m2 it reached 930 K in one 0.5 s step, found
-    while planning), and drawing at its temperature would draw liquid that the next sub-step removes, without limit;
-    what the 3D model actually loses is reconciled by the body after the base moves. With `g_base` [K/m per column], the
+    then, with `flow` (a `Flow`), the flow's sub-steps (`flow_substeps`: heat, then flow, sideways exchange and -- in
+    stage 4 -- spray, spec 7.1; sub-plan 21 Task 6c), then -- with `remove` and a melting material -- the contiguous top
+    cells at or above T_feed leave, taking their momentum with them, a column below its minimum mass m_min [kg] draws
+    back to it, and the cells are split and merged. The draw takes rigid metal (spec 6.6): at the base temperature, but
+    never above T_rigid -- a 3D base node can overshoot it within a macro step (on 4 mm elements under 6 MW/m2 it
+    reached 930 K in one 0.5 s step, found while planning), and drawing at its temperature would draw liquid that the
+    next sub-step removes, without limit; what the 3D model actually loses is reconciled by the body after the base
+    moves. With `g_base` [K/m per column], the
     3D model's temperature gradient along the outward normal below the base, each draw takes the next layer down: its
     enthalpy is h(min(T, T_rigid)) averaged along T(s) = T_b - g_base s over the depth it takes (`profile_enthalpy`).
 
@@ -644,14 +648,17 @@
     (a, b) per column, E_base / dt = a + b T_b [W into the 3D model], from the pass's base energy and its derivative with
     respect to the base temperature, the slope clipped at zero (stage 1's choice: a law rising with T_b would cost the 3D
     matrix its diagonal dominance; what the clip misses is the interface mismatch's to book). Otherwise the real pass,
-    base held at T_b1, the 3D model's end-of-step base temperature; returns a SkinBook."""
+    base held at T_b1, the 3D model's end-of-step base temperature; returns a SkinBook. The trial pass runs the flow too
+    (spec 7.2's sub-steps are spec 7.1's: heat, flow, sideways exchange and spray), the sensitivity S carried with the
+    metal the transport moves."""
     mat = sk.material
     N = max(1, int(math.ceil(dt / sk.settings.substep - 1e-9)))
     h = dt / N
     trial = T_b1 is None
     saved = sk.state() if trial else None
     T_b = np.asarray(T_b0 if trial else T_b1, dtype=float) * np.ones(sk.n)
-    book = SkinBook(*(np.zeros(sk.n) for _ in range(8)))
+    book = SkinBook(*(np.zeros(sk.n) for _ in range(8)), momentum_out=np.zeros((sk.n, 3)),
+                    flow=FlowBook.zero(sk.n) if flow is not None else None)
     dEb = np.zeros(sk.n)
     if trial:
         sk.S = np.zeros_like(sk.m)
@@ -669,11 +676,14 @@
         book.E_side += res[3]
         if trial:
             dEb += res[4]
+        if flow is not None:
+            book.flow.add(flow_substeps(sk, flow, h))
         if melting:
             k = np.cumprod(sk.valid & (sk.T >= mat.T_feed), axis=1).sum(axis=1)
             gone, gone_H = sk.take_top(k)
             book.removed += gone
             book.removed_H += gone_H
+            book.momentum_out += sk.momentum_out
         short = np.maximum(m_min - sk.column_mass(), 0.0)
         if g_base is not None and np.any(short > 0.0):
             rows = np.flatnonzero(short > 0.0)
@@ -1394,6 +1404,107 @@
                          faces=FacePush(r, across), growth=float(max(dt * s.max(), 0.0)) if n else 0.0, passes=passes)
 
 
+# -- the flow in the sub-steps (sub-plan 21 Task 6c) --------------------------------------------------------------------
+FLOW_CYCLES = 10                  # flow sub-steps per heat sub-step: 1 ms under the 10 ms heat sub-step (Task 6c's
+                                  # ruling, 2026-10-22). On the 100 mm sphere's front melted 0.3 mm deep under the 93 s
+                                  # forcing, 50 ms of flow at 1 ms match 0.25 ms to 3.9 % in where the metal went and
+                                  # 4 % in kinetic energy, the Courant numbers at most 4.8 and the transport destroying
+                                  # 27 % of the kinetic energy made, as at 0.25 ms (the grid's own upwind smearing); at
+                                  # 2-10 ms 44-95 % and 0.3-4.8 times, because a runoff front crossing rigid metal
+                                  # advances at most one column a flow sub-step (metal landing on a rigid top joins the
+                                  # flow at the next), 2 m/s on 2 mm columns at 1 ms. The Rayleigh-Taylor growth dt s
+                                  # stays under 0.007 at 1 ms on 3.4 mm of liquid at 48 m/s2 (0.67 at 10 ms).
+FLOW_SECONDS = {"push": 0.0, "momentum": 0.0, "thickness": 0.0, "transport": 0.0, "regrid": 0.0}   # wall time, for
+                                  # measuring the cost (as KRYLOV counts the coupled solves)
+
+
+@dataclass
+class Flow:
+    """What drives the skins' flow over a part of a macro step, held over it (spec 6.4): `dp_ds` [Pa/m, (n,)] the
+    surface flow's wall-pressure gradient along the meridian, `shear` [Pa, (n, 3)] the gas shear on the top cells in
+    each column's plane (`along_meridian` of the surface flow's tau), `g` [m/s2, (3,)] the acceleration the metal feels
+    in the body's frame, minus the body's own, `v_hat` the flight direction, `rho_l` [kg/m3] and `sigma` [N/m] the
+    liquid's density and surface tension, `A_top` [m2, (n,)] the top areas the shear acts over, and `cycles`, the flow
+    sub-steps in each heat sub-step."""
+    columns: object
+    dp_ds: np.ndarray
+    shear: np.ndarray
+    g: np.ndarray
+    v_hat: np.ndarray
+    rho_l: float
+    sigma: float
+    A_top: np.ndarray
+    cycles: int = FLOW_CYCLES
+
+
+@dataclass
+class FlowBook:
+    """What the flow sub-steps did, per column: impulses [N s] and momentum [kg m/s], (n, 3), which add up to the
+    columns' momentum change, and the sub-steps' worst numbers."""
+    push: np.ndarray         # the push's impulse on the moving cells, its thickness terms' implicit change included
+    shear: np.ndarray        # the gas shear's on the top cells
+    wall: np.ndarray         # the anchored metal's drag (the momentum and thickness steps'), the momentum cells held on
+                             # joining it, and the normal parts the transport's projections removed
+    carried: np.ndarray      # what the transport carried in less what it sent, as it left the donors: zero over the body
+    courant: float = 0.0     # the transport's largest outflow rate of any cell times its sub-step
+    growth: float = 0.0      # the thickness step's largest Rayleigh-Taylor growth dt s
+    passes: int = 0          # the thickness step's most linearisations
+    sub_steps: int = 0       # the flow sub-steps taken (none where nothing flows)
+
+    @classmethod
+    def zero(cls, n):
+        return cls(*(np.zeros((n, 3)) for _ in range(4)))
+
+    def add(self, other):
+        """Add another book's impulses and keep the worst numbers of both."""
+        for key in ("push", "shear", "wall", "carried"):
+            setattr(self, key, getattr(self, key) + getattr(other, key))
+        self.courant, self.growth = max(self.courant, other.courant), max(self.growth, other.growth)
+        self.passes, self.sub_steps = max(self.passes, other.passes), self.sub_steps + other.sub_steps
+
+
+def flow_substeps(sk, flow, dt):
+    """`flow.cycles` flow sub-steps of dt / cycles [s] each, every one in the order of 6b's ruling: the push, surface
+    tension only where a column's top flows (`push`, `flowing_tops`); the momentum through the depth with its response
+    (`momentum_step`); the thickness terms implicit in the heights at the end of the transport (`thickness_step`); the
+    transport with the faces' own push (`transport`, regrid=False). The cells are split and merged between the flow
+    sub-steps, and the heat sub-step's own regrid follows the last (spec 6.2: once after each sub-step's exchanges).
+    Updates the cells; returns a FlowBook. Where nothing flows nothing moves and no sub-step is taken: any velocity left
+    in the rigid cells is booked to the wall and cleared, as `momentum_step` books an anchored cell's.
+
+    The flow's sub-step is chosen from measurements on the 100 mm sphere's melting state (sub-plan 21 Task 6c, Asha
+    2026-10-21, option c; FLOW_CYCLES)."""
+    n = sk.n
+    book = FlowBook.zero(n)
+    if not (sk.valid & np.isfinite(cell_viscosity(sk, 0.0))).any():
+        book.wall -= sk.column_momentum()
+        sk.U = np.zeros_like(sk.U)
+        return book
+    h = dt / flow.cycles
+    col = flow.columns
+    clock = lambda key, started: FLOW_SECONDS.__setitem__(key, FLOW_SECONDS[key] + time.perf_counter() - started)
+    for k in range(flow.cycles):
+        if k:
+            started = time.perf_counter()
+            sk.regrid()
+            clock("regrid", started)
+        started = time.perf_counter()
+        G, terms = push(col, sk.thickness(), flow.dp_ds, flow.g, flow.v_hat, flow.rho_l, flow.sigma, wet=flowing_tops(sk))
+        clock("push", started)
+        started = time.perf_counter()
+        mb = momentum_step(sk, h, G, flow.shear, flow.A_top, response=True)
+        clock("momentum", started)
+        started = time.perf_counter()
+        tb = thickness_step(sk, col, h, mb, terms["thickness"] + terms["capillary"], flow.g, flow.rho_l, flow.sigma)
+        clock("thickness", started)
+        started = time.perf_counter()
+        tr = transport(sk, col, h, regrid=False, faces=tb.faces)
+        clock("transport", started)
+        book.add(FlowBook(mb.push + tb.push, mb.shear, mb.wall + tb.wall + tr.wall, tr.momentum, tr.courant, tb.growth,
+                          tb.passes, 1))
+    return book
+
+
 # -- the base (Task 7) -------------------------------------------------------------------------------------------------
 def hand_down(sk, m_min, T_rigid):
     """Rigid cells at the bottom of each column (at or below T_rigid, contiguous from the bottom) pass to the 3D model,
```

- [ ] **Step 6: The flow in the skin body** in `reentry_model/skinbody.py`:

```diff
--- a/reentry_model/skinbody.py
+++ b/reentry_model/skinbody.py
@@ -38,7 +38,15 @@
 step receives the trapezoid of the two (Asha, 2026-10-12).
 
 The body's surface for the coupled run is the top surface, one 'patch' per corner (`columns.CornerSurface`), so the
-heating, the history and the trajectory's loop are unchanged."""
+heating, the history and the trajectory's loop are unchanged.
+
+The flow (sub-plan 21 Task 6c; `flow`, off by default): every macro step the surface flow is evaluated on the corners
+as the heating is, from the step's aero state and the tops' temperatures -- its wall-pressure gradient along the
+meridian and its gas shear -- with the deceleration, which the metal feels along the flight direction; each part takes
+them interpolated in time at its middle, as it takes the heating, and every sub-step of both passes runs the flow after
+its conduction (`skin.flow_substeps`), `flow_cycles` flow sub-steps to each. With the flow on the liquid stays and flows
+(`remove` False in a run; spraying is stage 4's). The momentum books: the skins' momentum changes by the flow's impulses
+less what leaves with removed or handed-down metal, and `momentum_residual` reports what they miss."""
 import dataclasses
 import math
 import time
@@ -50,6 +58,9 @@
 from . import skin as skin_mod
 from .body import ThermalBody
 from .thermal import SIGMA_SB
+
+FLOW_DEFAULT = False              # the flow is off until stage 4's spraying gives the liquid a sink: with it on and nothing
+                                  # sprayed, a flight would only pile the liquid up (Task 6c's ruling, 2026-10-22)
 
 GRADIENT_FLOOR = 1.0              # K/m: the smallest gradient along the normal that is divided by. A corner above T_rigid
                                   # whose gradient is below it (a flat or reversed field) is given this one, so its depth
@@ -88,6 +99,22 @@
                                   # (measured 2026-10-14). The element model's runs keep their solver's own tolerance.
 
 
+@dataclasses.dataclass
+class FlowForcing:
+    """The gas side's push on the skins' metal at a macro step, per corner: `dp_ds` [Pa/m, (n,)] the wall pressure's
+    gradient along the meridian and `tau` [Pa, (n,)] the gas shear (the surface flow's), and `g` [m/s2, (3,)] the
+    acceleration the metal feels in the body's frame -- minus the body's own: the deceleration along the flight
+    direction."""
+    dp_ds: np.ndarray
+    tau: np.ndarray
+    g: np.ndarray
+
+    def toward(self, other, frac):
+        """This forcing moved the fraction `frac` of the way to `other`, linearly."""
+        return FlowForcing(*(a + (b - a) * frac for a, b in ((self.dp_ds, other.dp_ds), (self.tau, other.tau),
+                                                             (self.g, other.g))))
+
+
 def jump_faults(T_bottom, T_node, T_deep, crossed, tol=None):
     """The columns whose skin's bottom cell lies more than `tol` (JUMP_TOLERANCE) outside the temperatures of the metal
     that crossed their base (spec 6.6): between the base node's and the 3D profile's at the depth swept, `T_deep` (equal
@@ -101,13 +128,20 @@
 
 class SkinBody(ThermalBody):
     """A ThermalBody whose surface is resolved by corner columns with layered skins on a moving mesh (module docstring).
-    `remove`: liquid leaves as it forms (stage 2); `demise_fraction`: the run ends below this share of the initial mass."""
+    `remove`: liquid leaves as it forms (stage 2); `demise_fraction`: the run ends below this share of the initial mass;
+    `flow`: the metal flows (sub-plan 21), `flow_cycles` flow sub-steps in each heat sub-step (None: skin.FLOW_CYCLES),
+    driven by `surface_flow` (a `surface_flow.SurfaceFlow`, its defaults if None) from each macro step's aero state."""
 
     def __init__(self, mesh, material, solver, mass_kg, settings=None, T0=300.0, emissivity=None, T_ambient=0.0,
                  v_hat=(1.0, 0.0, 0.0), deform_depth=movemesh.DEFORM_DEPTH, remove=True, demise_fraction=0.01,
-                 repeat_rule=DEFAULT_REPEAT_RULE, part_max=PART_MAX):
+                 repeat_rule=DEFAULT_REPEAT_RULE, part_max=PART_MAX, flow=FLOW_DEFAULT, flow_cycles=None,
+                 surface_flow=None):
         if repeat_rule not in REPEAT_RULES:
             raise ValueError("repeat_rule must be one of {}, got {!r}".format(REPEAT_RULES, repeat_rule))
+        self.flow_cycles = skin_mod.FLOW_CYCLES if flow_cycles is None else flow_cycles
+        if not (isinstance(self.flow_cycles, (int, np.integer)) and self.flow_cycles >= 1):
+            raise ValueError("flow_cycles must be an integer of at least 1, got {!r}".format(flow_cycles))
+        self.flow_on, self.surface_flow = bool(flow), surface_flow
         self.repeat_rule = repeat_rule
         self.part_max = part_max                             # s, or None: no cap on the parts while melting
         self.material, self.solver, self.mass_kg = material, solver, mass_kg
@@ -138,11 +172,18 @@
         self._dt_last = None                                 # s, and its length: the look-ahead of the part cap
         self.capped_ahead = 0                                # macro steps capped because a top was about to be non-rigid
         self._heating = None                                 # this macro step's (t0, dt, q at its start, q at its end)
-        self.wall = {"skins": 0.0, "3d": 0.0, "move": 0.0, "rebuild": 0.0}
+        self.wall = {"skins": 0.0, "3d": 0.0, "move": 0.0, "rebuild": 0.0, "flow": 0.0}
+        self.last_flow = None                                # the surface flow's last evaluation (SurfaceFlowResult)
+        self._forcing = self._forcing_prev = None            # this macro step's (start, end) FlowForcing; the last end
+        self.impulse, self.momentum_lost = np.zeros(3), np.zeros(3)   # the momentum books since time zero [N s, kg m/s]
+        self.impulse_size = 0.0                              # N s: the push's, shear's and walls' impulses, in size
+        self.flow_courant = self.flow_growth = 0.0           # this macro step's worst flow sub-step
+        self.flow_passes = 0
         self._volume = None                                  # m3, the 3D model's volume since the last move (mass)
         self.iterations, self.last, self.last_move = [], None, None
         self._fill(T0)
         self.energy0 = self.energy()
+        self.momentum0 = self.skins.momentum()
         self._top0, self._normal0 = self.surface.points.copy(), self.surface.normals.copy()
 
     # -- set-up -----------------------------------------------------------------------------------------------------
@@ -175,10 +216,18 @@
         return self.moving.mesh
 
     # -- the macro step ---------------------------------------------------------------------------------------------
-    def advance(self, t, dt, loads, state=None):
+    def advance(self, t, dt, loads, state=None, forcing=None):
         """The macro step from t - dt to t (module docstring) as equal parts, as many as the last macro step's largest
-        move called for."""
+        move called for. With the flow on, the gas side's push is the surface flow's on the corners from the aero state
+        `state` (`flow_forcing`), or `forcing`, a FlowForcing given instead (a prescribed push, in the tests)."""
         n, self._ratio, self.step_parts = self._hint, 0.0, 0
+        if self.flow_on:
+            end = forcing if forcing is not None else self.flow_forcing(state)
+            last = self._forcing_prev
+            start = last if last is not None and last.dp_ds.shape == end.dp_ds.shape else end
+            self._forcing, self._forcing_prev = (start, end), end
+            self.flow_courant = self.flow_growth = 0.0
+            self.flow_passes = 0
         sk, mat = self.skins, self.material
         top = np.where(sk.valid[:, 0], sk.T[:, 0], np.nan)
         if self.part_max is not None and mat.melts:
@@ -199,13 +248,36 @@
         self._hint = max(1, int(math.ceil(MOVE_HINT_MARGIN * self._ratio)))
         self._q_prev = q_end.copy()
 
+    def flow_forcing(self, state):
+        """The gas side's push on the metal at this macro step (FlowForcing): the surface flow evaluated on the corners
+        from the aero state, as the heating is -- the top surface's angles to the flight direction, the nose radius and
+        the tops' temperatures at the step's start -- and the deceleration along the flight direction."""
+        if state is None:
+            raise ValueError("a skin body with the flow on needs the macro step's aero state, or a forcing given")
+        if self.surface_flow is None:
+            from . import surface_flow
+            self.surface_flow = surface_flow.SurfaceFlow()
+        mat = self.material
+        rho_l = mat.liquid.rho if mat.liquid is not None else mat.rho
+        res = self.surface_flow.evaluate(state, self.theta, self.nose_radius(), rho_l, self.surface_temperature())
+        self.last_flow = res
+        return FlowForcing(np.asarray(res.dp_ds, dtype=float), np.asarray(res.tau, dtype=float),
+                           res.deceleration * self.v_hat)
+
     def _advance_part(self, t, dt, dt_macro, loads, depth):
         """One part of a macro step, ending at t: the passes, then -- if its draws and hand-downs would move a base more
         than MOVE_MAX of the mean surface edge -- the part again as shorter parts, else the books and the move."""
         mat, sk, mm = self.material, self.skins, self.moving
         t0, dt_m, q_start, q_end = self._heating
-        q = q_start + (q_end - q_start) * (((t - 0.5 * dt) - t0) / dt_m)     # the heating at the part's middle
+        frac = ((t - 0.5 * dt) - t0) / dt_m
+        q = q_start + (q_end - q_start) * frac                              # the heating at the part's middle
         A_top = self.surface.areas
+        flow = None
+        if self.flow_on:                                                    # the gas side's push at the part's middle
+            f = self._forcing[0].toward(self._forcing[1], frac)
+            sigma = mat.liquid.sigma if mat.liquid is not None else 0.0
+            flow = skin_mod.Flow(self.columns, f.dp_ds, skin_mod.along_meridian(self.columns.normal, self.v_hat, f.tau),
+                                 f.g, self.v_hat, sk.rho_l, sigma, A_top, self.flow_cycles)
         eps = np.full(sk.n, self.emissivity)
         m_min = skin_mod.minimum_mass(mat, sk.area, self.settings)
         side = lambda s, previous: skin_mod.Sideways(s, self.columns, previous)   # implicit (sub-plan 21)
@@ -218,10 +290,11 @@
         T_lin = T_b0
         damped = self.repeat_rule == "damped"
         omega, r_prev = DAMPED_OMEGA0, None
+        flow_seconds = sum(skin_mod.FLOW_SECONDS.values())
         for attempt in range(DAMPED_MAX_PASSES if damped else 2):
             started = time.perf_counter()
             a, b = skin_mod.advance(sk, dt, q, A_top, self.T_ambient, eps, T_lin, m_min, side=side, remove=self.remove,
-                                    g_base=g0)
+                                    g_base=g0, flow=flow)
             self.wall["skins"] += time.perf_counter() - started
             started = time.perf_counter()
             res = self.solver.step(dt, np.zeros(len(self.solver.areas)), self.T_ambient, nodal_load=load,
@@ -230,7 +303,7 @@
             T_b1 = self.solver.temperature()[nodes]
             started = time.perf_counter()
             book = skin_mod.advance(sk, dt, q, A_top, self.T_ambient, eps, T_b0, m_min, T_b1=T_b1, side=side,
-                                    remove=self.remove, g_base=g0)
+                                    remove=self.remove, g_base=g0, flow=flow)
             self.wall["skins"] += time.perf_counter() - started
             mismatch = book.E_base - (a + b * T_b1) * dt
             scale = np.abs(book.E_top) + np.abs(book.E_base) + 1e-30
@@ -264,6 +337,7 @@
                 self.interface_repeats += 1
                 continue
             break
+        self.wall["flow"] += sum(skin_mod.FLOW_SECONDS.values()) - flow_seconds
         ratio = self._move_ratio(book, m_min)
         if ratio > 1.0 and depth < SPLIT_DEPTH_MAX:          # too far for one part: run it again as shorter parts
             self.solver.set_state(solver_state)
@@ -284,6 +358,13 @@
         self.column_removed += book.removed
         self.removed_enthalpy += float(book.removed_H.sum())
         self.drawn_mass += float(book.drawn.sum())
+        self.momentum_lost += book.momentum_out.sum(axis=0)
+        if book.flow is not None:                            # the flow's impulses (the transport's carried sum to zero)
+            fb = book.flow
+            self.impulse += (fb.push + fb.shear + fb.wall + fb.carried).sum(axis=0)
+            self.impulse_size += float(sum(np.linalg.norm(x, axis=1).sum() for x in (fb.push, fb.shear, fb.wall)))
+            self.flow_courant, self.flow_growth = max(self.flow_courant, fb.courant), max(self.flow_growth, fb.growth)
+            self.flow_passes = max(self.flow_passes, fb.passes)
         self.iterations.append(res.iterations)
         self.last = dataclasses.replace(res, Q_conv=E_top / dt, Q_rad=E_rad / dt)    # the gas-side powers, for the history
         self._move_base(book, m_min, g0)
@@ -343,6 +424,8 @@
         started = time.perf_counter()
         T_node = self.solver.temperature()[mm.corners]
         handed, handed_H = skin_mod.hand_down(sk, m_min, mat.T_rigid) if mat.melts else (np.zeros(sk.n), np.zeros(sk.n))
+        if mat.melts:
+            self.momentum_lost += sk.momentum_out.sum(axis=0)                     # what the handed-down cells held
         moved = (book.drawn - handed) / (mat.rho * col.area)                      # m inward by drawing and handing down
         depth = self._absorb_depth(np.maximum(MOVE_MAX * col.edge_length.mean() - moved, 0.0))
         inward = book.drawn / mat.rho + depth * col.area - handed / mat.rho          # m3 per corner
@@ -357,6 +440,7 @@
         T_deep = T_node - g0 * (book.drawn / (mat.rho * col.area) + depth)   # the profile at the depth swept
         self.jump_faults += int(jump_faults(T_bottom, T_node, T_deep, crossed).sum())
         sk.exchange_bottom(mat.rho * r.swept_volume - (book.drawn - handed), r.swept_enthalpy - (book.drawn_H - handed_H))
+        self.momentum_lost += sk.momentum_out.sum(axis=0)                         # what a reconciliation took away
         sk.regrid()
         self.absorbed_mass += float((mat.rho * depth * col.area).sum())
         self.handed_mass += float(handed.sum())
@@ -384,6 +468,13 @@
             return 0.0
         return (self.energy() - self.energy0 - self.absorbed_heat + self.removed_enthalpy) / self.absorbed_heat
 
+    def momentum_residual(self):
+        """|P - P0 - impulses + momentum lost| over the impulses' total size: zero for exact books (spec 10.2)."""
+        if self.impulse_size <= 0.0:
+            return 0.0
+        miss = self.skins.momentum() - self.momentum0 - self.impulse + self.momentum_lost
+        return float(np.linalg.norm(miss) / self.impulse_size)
+
     def mean_temperature(self):
         return float(self.material.temperature_from_enthalpy(self.energy() / self.mass()))
 
@@ -420,6 +511,7 @@
         flowing = sk.depth_above(mat.T_rigid) if mat.melts else np.zeros(sk.n)
         liquid = sk.depth_above(mat.T_liquidus) if mat.melts else np.zeros(sk.n)
         r = self.last_move
+        i, j = self.columns.edges[:, 0], self.columns.edges[:, 1]
         return {"skin_cells": int(sk.count.sum()), "skin_thickness_max_mm": float(th.max() * 1e3),
                 "skin_above_min": int((sk.column_mass() > m_min * (1.0 + 1e-9)).sum()),
                 "flowing_depth_mean_mm": float(flowing.mean() * 1e3), "flowing_depth_max_mm": float(flowing.max() * 1e3),
@@ -437,7 +529,15 @@
                 "removed_mass_kg": self.removed_mass, "removed_enthalpy_J": self.removed_enthalpy,
                 "base_recession_max_mm": float(self._recession().max() * 1e3),
                 "wall_skins_s": self.wall["skins"], "wall_3d_s": self.wall["3d"], "wall_move_s": self.wall["move"],
-                "wall_rebuild_s": self.wall["rebuild"], "step_parts": self.step_parts}
+                "wall_rebuild_s": self.wall["rebuild"], "step_parts": self.step_parts,
+                # the flow (sub-plan 21 Task 6c): its momentum books, the fastest cell, the top's roughness -- the largest
+                # height difference between neighbours over their spacing -- and the flow sub-steps' worst Courant
+                # number, Rayleigh-Taylor growth and linearisations in this macro step
+                "momentum_residual": self.momentum_residual(),
+                "speed_max_ms": float(np.linalg.norm(sk.U, axis=2).max()) if sk.U.size else 0.0,
+                "roughness_max": float((np.abs(th[i] - th[j]) / self.columns.edge_length).max()),
+                "flow_courant_max": self.flow_courant, "flow_growth_max": self.flow_growth,
+                "flow_passes_max": self.flow_passes, "wall_flow_s": self.wall["flow"]}
 
     def _recession(self):
         """How far each corner's top has moved inward along its starting direction since time zero [m]."""
```

- [ ] **Step 7: The history's flow columns, the results and the frames' top speed** in `reentry_model/coupled.py`:

```diff
--- a/reentry_model/coupled.py
+++ b/reentry_model/coupled.py
@@ -44,7 +44,10 @@
                     "handed_mass_kg", "shortfall_kg", "jump_faults", "absorb_capped", "interface_mismatch_J",
                     "interface_repeats", "repeat_columns_over", "mismatch_worst_ratio", "mismatch_total_ratio",
                     "removed_mass_kg", "removed_enthalpy_J", "base_recession_max_mm",
-                    "wall_skins_s", "wall_3d_s", "wall_move_s", "wall_rebuild_s", "step_parts"]
+                    "wall_skins_s", "wall_3d_s", "wall_move_s", "wall_rebuild_s", "step_parts",
+                    # the flow (sub-plan 21 Task 6c, spec 10.2)
+                    "momentum_residual", "speed_max_ms", "roughness_max", "flow_courant_max", "flow_growth_max",
+                    "flow_passes_max", "wall_flow_s"]
 PVD_TEMPLATE = '<?xml version="1.0"?>\n<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">\n<Collection>\n{}</Collection>\n</VTKFile>\n'
 
 
@@ -179,7 +182,9 @@
                "demise_time_s": float(last["time_s"][-1]) if history.end_reason == "demise" else None,
                "surface_fold_max_deg": float(np.max(last["surface_fold_deg"])),
                "surface_quality_min": float(np.min(last["surface_quality_min"])),
-               "zone_ratio_min": float(np.min(last["zone_ratio_min"]))}
+               "zone_ratio_min": float(np.min(last["zone_ratio_min"])),
+               "momentum_residual": float(last["momentum_residual"][-1]),
+               "flow_courant_max": float(np.max(last["flow_courant_max"]))}
         out.update({k: float(last[k][-1]) for k in ("removed_mass_kg", "removed_enthalpy_J", "drawn_mass_kg", "absorbed_mass_kg",
                                                   "handed_mass_kg", "rebuild_residual_J", "wall_skins_s", "wall_3d_s",
                                                   "wall_move_s", "wall_rebuild_s")})
@@ -296,7 +301,8 @@
 def write_skin_frame(output_dir, k, body, loads):
     """A frame of the layered-skin body (sub-plan 20, spec 10.2): field_<k>.vtu, the 3D model's nodal T on its current
     (moving) mesh; surface_<k>.vtp, the top surface with per-corner point data -- top temperature, convective flux,
-    liquid and flowing depths, skin thickness and cells, and the base's position (the corner on the 3D mesh)."""
+    liquid and flowing depths, skin thickness and cells, the base's position (the corner on the 3D mesh) and the top
+    cell's speed (sub-plan 21 Task 6c)."""
     import pyvista as pv
     mesh, sk, mat = body.mesh, body.skins, body.material
     cells = np.hstack([np.full((mesh.n_elements, 1), 4), mesh.tets]).ravel()
@@ -312,6 +318,7 @@
     poly.point_data["liquid_depth"] = sk.depth_above(mat.T_liquidus) if mat.melts else np.zeros(sk.n)
     poly.point_data["flowing_depth"] = sk.depth_above(mat.T_rigid) if mat.melts else np.zeros(sk.n)
     poly.point_data["base_position"] = body.columns.points
+    poly.point_data["top_speed"] = np.where(sk.count > 0, np.linalg.norm(sk.U[:, 0], axis=1), 0.0)
     poly.save(os.path.join(output_dir, "surface_{}.vtp".format(k)))
 
 
```

- [ ] **Step 8: The settings** in `reentry_model/cli.py`:

```diff
--- a/reentry_model/cli.py
+++ b/reentry_model/cli.py
@@ -123,13 +123,14 @@
 
 
 def skin_name_suffix(min_thickness_mm, cell_um, growth, substep_ms, preset, deform_mm, swell=True,
-                     fine_band_mm=FINE_BAND_MM, fine_size_mm=FINE_SIZE_MM, fine_zone_mm=0.4):
+                     fine_band_mm=FINE_BAND_MM, fine_size_mm=FINE_SIZE_MM, fine_zone_mm=0.4, flow=False, flow_cycles=None):
     """The run-name suffix of a skin run's non-default settings (layered-skin design section 10.1): `_skin-dev` for the
     dev preset, `_skinmin<mm>mm` for another minimum thickness, `_skinfine<mm>mm` for another fine zone (both changed by
     option B', 2026-10-16: 0.8 and 0.4 mm), `_cell<um>um` for an explicit cell other than the
     preset's, `_growth<ratio>`, `_sub<ms>ms` and `_deform<mm>mm` for another growth, sub-step or deforming depth,
-    `_noswell` for skins held at the solid's density, and `_fine<mm>x<mm>mm` for another fine band (`_fine0` for none);
-    empty for the defaults."""
+    `_noswell` for skins held at the solid's density, `_fine<mm>x<mm>mm` for another fine band (`_fine0` for none), and,
+    with the flow on (sub-plan 21 Task 6c), `_skinflow-on` followed by `_flowcycles<n>` for other than skin.FLOW_CYCLES
+    flow sub-steps in each heat sub-step; empty for the defaults."""
     out = "_skin-dev" if preset == "dev" else ""
     if min_thickness_mm != 0.8:
         out += "_skinmin{:g}mm".format(min_thickness_mm)
@@ -149,6 +150,10 @@
         out += "_fine0"
     elif (fine_band_mm, fine_size_mm) != (FINE_BAND_MM, FINE_SIZE_MM):
         out += "_fine{:g}x{:g}mm".format(fine_band_mm, fine_size_mm)
+    if flow:
+        out += "_skinflow-on"
+        if flow_cycles is not None and flow_cycles != skin.FLOW_CYCLES:
+            out += "_flowcycles{}".format(flow_cycles)
     return out
 
 
@@ -302,6 +307,13 @@
                     help="on: each skin cell at its own phase's density, so melting swells the skin (default). off: every "
                          "cell at the solid's density, the verification setting of levels 2 and 2b, whose references hold "
                          "one density")
+    me.add_argument("--skin-flow", choices=("on", "off"), default="off",
+                    help="with --surface-model skin and --melt on: on, the skins' liquid and slurry flow over the surface "
+                         "(sub-plan 21) and stay -- the spraying that will take them is sub-plan 22's, so give --removal "
+                         "girin, whose continuum step (--dt-continuum) follows; off: the heat-only skins, liquid removed "
+                         "as it forms with --removal instant (default)")
+    me.add_argument("--skin-flow-cycles", type=int, default=None,
+                    help="with --skin-flow on: flow sub-steps in each skin sub-step (default {})".format(skin.FLOW_CYCLES))
     me.add_argument("--skin-preset", choices=tuple(skin.SKIN_PRESETS), default="production",
                     help="skin cell-size preset: production 10 um, dev 20 um (default %(default)s)")
     me.add_argument("--deform-depth", type=float, default=15.0,
@@ -366,9 +378,13 @@
         mat.k_table = mat.k_table * args.k_scale                  # verification device (near-isothermal body), not physical
     solver = thermal.thermal_solver(args.thermal_solver, linear_solver=args.linear_solver, lumped_mass=not args.consistent_mass)
     if skins:
+        flow = args.skin_flow == "on"                       # the liquid stays and flows; instant removal is instead of it
+        flow_model = surface_flow.SurfaceFlow(rarefied_shear=args.rarefied_shear, gamma_pm=args.gamma_pm,
+                                              kn_body_shock=args.kn_body_shock) if flow else None
         the_body = skinbody.SkinBody(the_mesh, mat, solver, mass, skin_settings(args), T0=args.temperature,
                                      emissivity=args.emissivity, T_ambient=args.t_ambient, deform_depth=args.deform_depth * 1e-3,
-                                     remove=melting, demise_fraction=args.demise_fraction)
+                                     remove=melting and not flow, demise_fraction=args.demise_fraction, flow=flow,
+                                     flow_cycles=args.skin_flow_cycles, surface_flow=flow_model)
     elif melting:
         flow = surface_flow.SurfaceFlow(rarefied_shear=args.rarefied_shear, gamma_pm=args.gamma_pm, kn_body_shock=args.kn_body_shock)
         spray_model = spray.SprayModel(mat.liquid, k_r=args.kr, k_t=args.kt, we_critical=args.we_critical,
@@ -399,7 +415,7 @@
                         "growth": args.skin_growth, "substep_ms": args.skin_substep, "preset": args.skin_preset,
                         "deform_depth_mm": args.deform_depth, "n_columns": the_body.columns.n,
                         "fine_band_mm": args.fine_band, "fine_size_mm": min(args.fine_size, h_surface * 1e3),
-                        "swell": args.skin_swell == "on"}
+                        "swell": args.skin_swell == "on", "flow": args.skin_flow, "flow_cycles": the_body.flow_cycles}
     if args.heating == "physics":
         info.update({"stagnation": args.stagnation, "bridging_heat": args.bridging_heat, "matting_n": args.matting_n,
                      "accommodation": args.accommodation, "catalycity": args.catalycity})
@@ -477,8 +493,18 @@
         parser.error("--skin-cell must be at most half of --skin-min-thickness and of --skin-fine-zone")
     if args.surface_model == "skin" and args.thermal != "fem":
         parser.error("--surface-model skin needs --thermal fem")
-    if args.surface_model == "skin" and args.melt == "on" and args.removal != "instant":
-        parser.error("--surface-model skin with --melt on needs --removal instant (sub-plan 20: no flow or spraying yet)")
+    if args.skin_flow == "on" and not (args.surface_model == "skin" and args.melt == "on"):
+        parser.error("--skin-flow on needs --surface-model skin and --melt on: it moves the skins' molten metal")
+    if args.skin_flow_cycles is not None and args.skin_flow != "on":
+        parser.error("--skin-flow-cycles needs --skin-flow on")
+    if args.skin_flow_cycles is not None and args.skin_flow_cycles < 1:
+        parser.error("--skin-flow-cycles must be at least 1")
+    if args.surface_model == "skin" and args.melt == "on" and args.skin_flow == "on" and args.removal != "girin":
+        parser.error("--skin-flow on keeps the liquid, which only spraying removes: it needs --removal girin, not instant "
+                     "removal, which stands in for flow and spraying in the heat-only skins (spec 11.2)")
+    if args.surface_model == "skin" and args.melt == "on" and args.skin_flow == "off" and args.removal != "instant":
+        parser.error("--surface-model skin with --melt on needs --removal instant, or --skin-flow on with Girin's removal "
+                     "(sub-plan 21: the skins' spraying is sub-plan 22's)")
     if args.surface_model == "skin" and (args.animate or args.stills):
         parser.error("--animate/--stills read the element model's frames; with --surface-model skin use --frames-every")
     if args.surface_model == "skin" and args.prism_layers is not None and args.prism_layers > 0:
@@ -501,7 +527,8 @@
                                        skin_suffix=skin_name_suffix(args.skin_min_thickness, args.skin_cell, args.skin_growth,
                                                                     args.skin_substep, args.skin_preset, args.deform_depth,
                                                                     args.skin_swell == "on", args.fine_band, args.fine_size,
-                                                                    args.skin_fine_zone))
+                                                                    args.skin_fine_zone, args.skin_flow == "on",
+                                                                    args.skin_flow_cycles))
     run_dir = os.path.join(args.outdir, name)
     os.makedirs(args.outdir, exist_ok=True)
     thermal_info, the_body = {}, body.ConstantBody(mass, args.temperature)
```

- [ ] **Step 9: Run them to see them pass.** Expected: `tests/test_reentry_model_skin.py` 72 passed,
  `tests/test_reentry_model_skinbody.py` 16 passed (2 deselected, the reference tier's), `tests/test_reentry_model_cli.py`
  66 passed, `tests/test_reentry_model_coupled.py` 9 and `tests/test_reentry_model_columns.py` 13 as before (176 in
  860 s); the unit tier 481 passed, 1 skipped (1 015 s). The flow is off by default, so no flight's result changes:
  the heat-only skins reproduce Task 6b's to the last bit.

- [ ] **Step 10: Commit.**

```bash
git add reentry_model/skin.py reentry_model/skinbody.py reentry_model/coupled.py reentry_model/cli.py tests/test_reentry_model_skin.py tests/test_reentry_model_skinbody.py tests/test_reentry_model_cli.py
git commit -m "Skins: the flow joins the skin body -- the surface flow's pressure gradient and shear on the corners and the deceleration drive ten flow sub-steps of push, momentum, thickness and transport in every 1 ms of each heat sub-step, in the trial and real passes, the momentum books in the history; off by default until spraying, the liquid kept with it on (sub-plan 21 Task 6c)"
```

---

### Task 6d: Making the flow affordable

**Files:**
- Modify: `reentry_model/skin.py` (`TRANSPORT_TOL`, `TRANSPORT_CORRECTIONS`, `TRANSPORT_CYCLE`, `TRANSPORT`,
  `LANDED_PASSES`, `LANDED_TOL`; `_upwind_order` returning the cycles' sizes, `_solve`, `_landed_velocity`, `transport`
  installing the cells from its fluxes and letting landed metal flow on; `Skins.flow_memory`, `.state()`,
  `.set_state()`, `.rows()`; `momentum_step` on the moving columns, `_momentum_columns`; `FLOW_CYCLES_MIN`,
  `FLOW_CYCLES_MAX`, `FLOW_COURANT`, `FLOW_GROWTH`, `FLOW_RELAX`, `flow_cycles`; `Flow.cycles` None by default and
  `Flow.plan`; `FlowBook.cycles`; `flow_substeps(..., cycles=None)`; `advance` giving the real pass its trial pass's
  counts), `reentry_model/skinbody.py` (`flow_cycles` None by default; the history's `flow_cycles_max`),
  `reentry_model/coupled.py` (`SKINBODY_COLUMNS`), `reentry_model/cli.py` (`--skin-flow-cycles` a fixed count, named
  `_flowcycles<n>` whenever given; the settings' `flow_cycles`)
- Test: `tests/test_reentry_model_skin.py` (eight tests; 6a's arrival test, whose metal now flows on, and its
  upwind-order test's stand-in order, which now returns the cycles' sizes), `tests/test_reentry_model_skinbody.py` (one
  test; two 6c tests: the fixed count's history column and the chosen count's default),
  `tests/test_reentry_model_cli.py` (the run name and the settings of a chosen count)

**Interfaces:**
- Consumes: 6a's `transport`, `_upwind_order`, `_factor`, `_fractions`; 6b's `thickness_step`; 6c's `Flow`, `FlowBook`,
  `flow_substeps`, `advance(..., flow)`, `SkinBody(..., flow, flow_cycles)`; Task 4's `momentum_step`.
- Produces: `skin.transport(sk, columns, dt, regrid=True, faces=None)` as before, each cell installed from the solve's
  fluxes and flowing metal landing on a rigid top flowing on in the same solve; `skin._solve(A, rhs, upwind, big=None)`;
  `skin._upwind_order(src, dst, n) -> (order, upwind, size)`, size (n,) each unknown's cycle;
  `skin._landed_velocity(columns, area, run, dt, base_m, base_p, stay) -> (n, 3)` [m/s]; `skin.TRANSPORT_TOL` (1e-13),
  `TRANSPORT_CORRECTIONS` (40), `TRANSPORT_CYCLE` (500 cells), `TRANSPORT` (the solves, corrections and fallbacks),
  `LANDED_PASSES` (200), `LANDED_TOL` (1e-6); `Skins.flow_memory` (3,) [Courant, growth, s], part of `state()` (seven
  arrays; `set_state` takes six as before, the memory then empty); `Skins.rows(idx)`; `skin._momentum_columns(sk, dt,
  push, shear, A_top, response)`, `momentum_step` solving only the columns where something moves;
  `skin.flow_cycles(memory, dt) -> int`, `FLOW_CYCLES` (10, now the count with no memory), `FLOW_CYCLES_MIN` (1),
  `FLOW_CYCLES_MAX` (40), `FLOW_COURANT` (6), `FLOW_GROWTH` (0.02), `FLOW_RELAX` (2); `Flow(..., cycles=None,
  plan=None)`, None choosing the count each heat sub-step; `FlowBook.cycles`, the most flow sub-steps a heat sub-step
  took; `flow_substeps(sk, flow, dt, cycles=None)`; `SkinBody(..., flow_cycles=None)`, None choosing; the history column
  `flow_cycles_max`; `run --skin-flow-cycles <n>` fixing the count (`_flowcycles<n>`), the settings' `flow_cycles`
  "adaptive" or n.

6c measured the flow at about 48 s of processor time a melting 10 ms heat sub-step against 4.5 s for the heat alone, the
transport three fifths of it, and the trial pass with the flow at 1 216 s against 17 s without it for one 0.0625 s part:
hours a melting macro step (Asha, 2026-10-22, option a: make it affordable before Tasks 7 and 8). The flow's physics and
exactness stay: the books close to round-off, the masses stay positive, and results change only where an accuracy
improvement below says so. Measured 2026-10-22 in `prototype/work-2026-10-18-column-flow/measure/task6d/`, on 6c's
states of the 100 mm sphere (9 417 columns of 10 um cells) made once with 6c's code (`states.py`): `heated` (the front's
skins at the conduction profile of 3 MW/m2, the flow developed over 20 ms), `driven` (the front melted 0.3 mm deep after
20 flow sub-steps of 1 ms under the 93 s forcing), `quiet` (the same after 10 under surface tension alone) and `all`
(every column melted, 5 driven); processor times, on the laptop under other sessions' load (load averages of 4 to 80),
so the ratios are the finding.

**1. The iterative transport solve** (6c's remedy, built). In upwind order the transport's systems are triangular but
for the back-couplings of their cycles, where flows meet, and their LU fills only within those: nothing for a cycle of a
few cells, five or six times over for one of ten thousand. `_solve` factors the matrix without the back-couplings among
the unknowns of cycles larger than TRANSPORT_CYCLE (500 cells) — the lower triangle and the small cycles, with almost no
fill — and iterates x_(k+1) = M^-1 (b - U x_k) from M^-1 b, U those back-couplings, whose residual U (x_k - x_(k+1))
costs no product with the whole matrix, until no row's residual exceeds TRANSPORT_TOL (1e-13) of the largest |b| + |A|
|x| in its column (a normwise backward error), at most TRANSPORT_CORRECTIONS (40) corrections; then the LU of the whole
system, as 6a's. **The books are exact by construction**: every cell is installed from the solve's fluxes, what it had
less what it sent plus what arrived, each arrival dt c times its donor's solved value, so each column's change is its
book whatever the solve's error. For the masses M is an M-matrix and -U >= 0, a regular splitting: the iterates rise
monotonically from below, every residual is non-negative, and the installed masses, x + (b - A x), are positive.

| State (largest cycle) | 6c's LU: mass, momentum | 6d: mass, momentum | Corrections | Off the LU |
|---|---|---|---|---|
| every column melted (18 379 cells) | 4.65 s, 25.3 s | 0.16 s, 0.25 s | 9, 9 | 1.5e-13 |
| quiescent melt (10 190) | 0.20 s, 0.76 s | 0.10 s, 0.28 s | 3, 3 | 9e-15 |
| heated melting state (11 109) | 0.43 s, 2.04 s | 0.14 s, 0.21 s | 5, 5 | 9e-15 |
| driven at 93 s (1 700) | 0.03 s, 0.07 s | 0.07 s, 0.10 s | 4, 4 | 1.4e-13 |

The table shows the large cycles' cost gone — 2 to 100 times faster where 6c's LU filled — at the price of a few
hundredths of a second where the largest cycle is small enough for the LU (the driven state's, iterated at 500 cells,
filled 1.1 times); each answer within 1.5e-13 of the LU's, the masses non-negative, the total mass changed by 2e-16.
Whole flow sub-steps of 1 ms, 6c's code and 6d's one after the other from each state (`solve.py`, under loads of 23-42,
so the splits, wall time as 6c's `FLOW_SECONDS` are, run above the processor totals): 2.91 against 1.84 s on the
quiescent melt (the transport 2.10 against 1.19 s), 47.3 against 4.4 s with every column melted, 7.2 against 3.7 s on
the heated state (the transport 5.1 against 2.0 s), and 2.63 against 2.67 s driven at 93 s, where the cycles are small
and the landed metal's estimate (0.07 s), the iterated 1 700-cell cycle (about 0.1 s more than its LU) and the
installation from the fluxes cancel the momentum step's saving (1.04 against 0.78 s).

**2. The front limit lifted.** 6a's rule made metal landing on a rigid top a new top cell that moved only in the next
flow sub-step, so a front crossing rigid metal advanced at most one column a flow sub-step — 2 m/s on 2 mm columns at 1
ms, the reason for 6c's ten flow sub-steps. The landed metal still lands on top (spec 6.5), but its new top cell is now
a flowing donor of the same solve: it sends across its column's faces at c = (u . f)+ / A, by fraction into a flowing
receiver's whole flowing layer and on top of a rigid one, whose new cell flows on in turn. Its velocity u must be known
before the solve for the step to stay linear: **it coasts at the velocity it arrived with** — the mass-weighted mean of
what lands there, each arrival's velocity projected onto the column's plane — and the drag of the rigid top beneath it
reaches it at the next momentum step, as every other cell moves at the velocity its momentum step gave it.
`_landed_velocity` estimates the weights before the solve (each donor cell's share dt c / (1 + dt sum c), and a landed
column sending on its estimate likewise), pass by pass a column further, until no velocity moves by more than 1e-6 of
the fastest, at most 200 passes; landed metal estimated below SLIVER (1e-12) of the largest landing stays where it lands
(the far end of the implicit scheme's geometric tail; in the one-cell variant below its products underflowed to a cell
of no mass). The weights average the velocities only: the solve carries the masses, the enthalpy and the momentum, the
arrivals' projected, so the system stays an M-matrix whose columns sum to one, positive and conservative, and the
momentum is booked as before. A new top cell that crust metal reached stays (rigid metal on a rigid top is anchored).

Measured: liquid moving at 0.5 m/s onto rigid metal on a strip of 1 mm columns crosses all 15 rigid columns in one 10 ms
sub-step (a Courant number of 5), each holding 0.831-0.836 of the last — the one-dimensional implicit upwind scheme's nu
/ (1 + nu) = 5/6 — at its arrival velocity to 2e-16. A sheared layer running onto rigid metal (20 Pa, 20 ms, the fastest
cell 1.4 m/s; `strip_runoff`): where the metal went beyond the edge, against the 0.25 ms run, at flow sub-steps of 5, 2,
1 and 0.5 ms:

| Flow sub-step | 6c (6a's rule) | Coasting (built) | One cell on its rigid top |
|---|---|---|---|
| 5 ms | 1.25 | 0.435 | 0.36 |
| 2 ms | 0.28 | 0.209 | 0.25 |
| 1 ms | 0.195 | 0.135 | 0.174 |
| 0.5 ms | 0.10 | 0.068 | 0.104 |

The table shows coasting the closest at the sub-steps the adaptive count chooses, and first order in the sub-step. The
alternative measured, the velocity a momentum step would give the landed metal as one cell on its rigid top (its wall
drag over its own half-thickness, the push and the shear), holds a front back almost as 6a's rule did, because the metal
landing within one sub-step is a thin precursor of the deposits the fine run builds up over many (about 145 um here,
coasting for tens of milliseconds); on the sphere it matched coasting within 0.012 at 1-2.5 ms (0.025, 0.081 and 0.171
against 0.022, 0.081 and 0.159) and stayed stable at 10 ms flow sub-steps, where coasting does not, but the adaptive
count does not choose those for a fast flow.

**3. The flow sub-step again** (6c's study, `substep.py` and `analyse_substep.py`: the 100 mm sphere's front melted 0.3
mm deep over rigid metal at 850 K, from rest, the forcing of the flight at 93 s, surface tension, 50 ms of flow alone,
against the 0.25 ms run):

| Flow sub-step | Courant, largest (median) | Where the metal went, L1 (6c) | KE at 50 ms / 0.25 ms run (6c) | Growth dt s | Linearisations, most |
|---|---|---|---|---|---|
| 10 ms | 4 143 (297) | 0.665 (0.953) | 2.98 (4.81) | 0.018 | 8 (the cap) |
| 5 ms | 26.4 (24.9) | 0.432 (0.854) | 0.90 (0.31) | 0.0048 | 8 |
| 2.5 ms | 12.6 (11.0) | 0.159 (0.586) | 0.95 (0.36) | 0.0013 | 3 |
| 2 ms | 9.9 (8.4) | 0.081 (0.439) | 1.06 (0.48) | 8.3e-4 | 8 |
| 1.5 ms | 7.3 (6.2) | 0.036 | 1.12 | 4.7e-4 | 8 |
| 1.25 ms | 6.1 (5.1) | 0.025 | 1.05 | 3.3e-4 | 4 |
| 1 ms | 4.8 (4.1) | 0.022 (0.039) | 1.03 (1.04) | 2.1e-4 | 3 |
| 0.5 ms | 2.4 (2.0) | 0.014 (0.014) | 1.01 (1.02) | 5.4e-5 | 2 |

The table shows the runoff converging in the flow sub-step from about 2 ms where 6c's converged only at 1 ms: at 2 ms
the metal is 8.1 % off the finest run's places against 44 %, and over the rigid rear the bands of 15 degrees from 90 to
150 degrees are within 7 % of the finest run's, 150-165 degrees 13 % high and 165-180 degrees 32 % short (6c's stopped
at 150 degrees); at 1.25 ms every band is within 6.2 %. The transport destroys 0.28-0.30 of the kinetic energy made at
1-2 ms, the grid's own 0.27 (0.43 at 2 ms in 6c). The growth is the same as 6c's: it does not depend on the front. At 10
ms the thickness step's linearisations stop at their cap, its predictions 28 times the change off, and the fastest cell
reaches 69 m/s at 40 ms: 6b's thickness step is unstable there, with the front limit or without it. At 1.5 and 2 ms its
linearisations reached their cap in a sub-step too, the kinetic energy 6-12 % high.

**4. Adaptive flow sub-steps.** The skins remember the last heat sub-step's flow (`Skins.flow_memory`: its largest
Courant number C and growth g over its flow sub-steps, and their length h_0, part of the skins' state so that a trial
pass, a repeated part and a shortened part all start from it). `flow_cycles` takes the next flow sub-step as the longest
of FLOW_RELAX h_0, h_0 FLOW_COURANT / C and h_0 (FLOW_GROWTH / g)^(1/2) — the Courant number grows in proportion to the
flow sub-step and the growth as its square (6c's table) — and the count as the heat sub-step over it, rounded up,
between FLOW_CYCLES_MIN (1: 10 ms) and FLOW_CYCLES_MAX (40: 0.25 ms); FLOW_CYCLES (10, 6c's count) where the skins
remember no flow, at the start or after a heat sub-step in which nothing flowed. The trial pass leaves the counts it
chose in the flow's plan and the real pass after it takes the same, so that the base law is made by the pass it
predicts; a real pass with no trial pass before it chooses its own. `--skin-flow-cycles <n>` still fixes the count, and
since a fixed count is no longer the default its run name always says so (`_flowcycles<n>`, spec 10.1's convention).
Measured on a bump of liquid aluminium on 1 mm columns: a still layer under surface tension relaxes from 10 flow
sub-steps to the floor (10, 5, 3, 2, 1, 1, ...), one sheared at 50 Pa and pushed at 20 kPa/m settles at 6 (10, 7, 9, 8,
7, 6, ...); on the heated state the trial pass of one part chose 8, 7, 5, 3, 2, 2 and 2 as its flow slowed.

**5. The momentum step on the moving columns.** 6c's step solved every column, flowing or not. It now solves the columns
with a moving cell (`Skins.rows`, a copy of their cells) and gives the others rest, the momentum their anchored cells
held booked to the wall, as the full solve did: every velocity, impulse and response the same to the last bit (each
column's system is its own). On the heated state 6 450 of 9 417 columns move: 1.34-1.39 s a momentum step against
1.90-1.96 s; on the driven state 7 108: 1.13-1.15 s against 1.44-1.46 s.

**6. What this task found in Tasks 4 and 6b.**

*Task 4's shear-rate fixed point, five passes in every melting state 6c measured* (`shear.py`). It is neither an
unreachable tolerance nor an oscillation. Li et al.'s fit makes the slurry's flow curve, stress against shear rate, γ
η(γ), rise to a maximum at about 228 1/s and fall to the 367 1/s at which the law is held, beyond which the stress rises
again with the held viscosity: a stress between the two has three rates, and the middle one, on the falling branch, is
unstable under the fixed point (its slope -γ η'/η exceeds one). Cells near it drift away by 1.1-1.8 times a pass, to the
held rate or the lower branch. On the heated state the cells still moving by more than 1e-3 after five passes are 1 571
of the 35 419 in the bridge (906.4-908 K, Li's law at 10 % solid in the bridge's log-linear blend) and 19 of the 146 155
in the slurry, none of the 128 341 liquid ones; the worst one's shear rate walks 247, 260, 274, 293, 320 and 369 1/s
(beyond the 367 1/s hold) over the first six passes. 60 passes still leave 6 cells moving. What the cap costs: after
five passes the velocities are 5.3e-4 of the fastest from the 60-pass result at worst and 1.4e-6 on average over the
moving cells (1.5e-3 and 8.8e-6 after one). On the driven state, liquid but for 8 slurry cells, the fixed point ends at
4 passes, and all liquid at 1. Not changed (Task 4's design): stopping on the velocities instead of the viscosities, or
holding Li's rate at the curve's maximum (a rheology change, spec 6.4), would end it sooner; each pass is about 0.14 s
of a heated flow sub-step's 1.0-1.4 s momentum step (Review Focus 18).

*6b's predicted heights at the shoulder's pile, 11-16 % of the change off at 1 ms* (`diag_pred.py`). Not the mixing's
density: the pile is liquid at 950 K throughout, its density change 1e-19 m. It is where the transport puts arriving
metal: at the same fraction of the receiver's flowing depth it left the donor's, so the donor's fast top layers arrive
in the pile column's fast top layers and pass on through it within the sub-step (cells 0-5 of column 9249 gained, 6-7
lost), while 6b's column model spreads a column's gain over its whole moving thickness at its mean rate, F / thick, and
so predicts the gain staying: 3.9e-5 m predicted against 1.8e-6 m made. The miss shrinks only once the top cells'
Courant number is well below one, 0.16, 0.22 and 0.08 of the largest change at 1, 0.25 and 0.05 ms: 6b's ruling's "where
the cells of a column carry very different shares", located. Not changed: a column model that weighted each column's
implicit departures by where its arrivals land would follow it, but the transport's fraction matching is what a column
model cannot reproduce (6b's ruling). Since 6d the column model also does not know that landed metal flows on, so at the
front's new tops it misses by up to 0.59 of the largest change in a sub-step (0.27 at 0.25 ms); there the tops are
rigid, no surface tension acts, and only the deceleration's thickness term, a few pascals across a face, depends on it
(Review Focus 10).

**7. The cost** (`cost.py` and `trial.py`, processor time, each case twice).

One 10 ms heat sub-step on the heated state, 6c's code and 6d's in turn from the same cells (the chosen count's after
one heat sub-step had set the skins' memory), under a load of 5-15:

| Code, flow sub-steps | Whole | Conduction | Push | Momentum | Thickness | Transport | Regrid |
|---|---|---|---|---|---|---|---|
| 6c, ten | 54-59 s | 5.0-5.9 | 0.6-0.7 | 11.5-12.5 | 3.6-3.8 | 31-33 | 1.9-2.0 |
| 6d, ten (fixed) | 43 s | 5.5-6.4 | 0.9 | 10.1-10.8 | 4.5-4.7 | 18.7-18.8 | 2.3-2.4 |
| 6d, chosen: nine | 30-33 s | 1.6 | 0.7 | 7.4-8.3 | 3.8-4.3 | 14.9-16.5 | 1.6-1.8 |

The table shows a heat sub-step 1.3 times cheaper at 6c's count and 1.6-2.0 times at the chosen one, the transport down
from three fifths of it to under half; the flow is still most of a heat sub-step (the chosen count's is the next 10 ms,
whose conduction took 1.6 s against 5-6 s). The quiescent and all-melted states freeze within a heat sub-step without
heating (their Courant numbers fell to 0.01), so their flowing cost is the flow sub-step's above.

*A melting part with its trial pass* (`trial.py`: one 0.0625 s part of seven heat sub-steps on the heated state, its
memory set as above): 6c's code (from the warm state) took 1 113 s for the trial pass with the flow and 924 s for the
real pass, 17 s for the trial pass without it; 6d's (10 ms on) 124 s and 115 s, choosing 8, 7, 5, 3, 2, 2 and 2 flow
sub-steps as the flow slowed, the real pass taking the same and meeting the law to 2.1e-15, and 14 s without the flow,
whose base law is 1.05 % of the throughput off the flow's (3 178 columns beyond the interface tolerance; 6.0 % from 6c's
state). A part is about 8.5 times cheaper, and the trial pass no longer dominates: it costs about what the real pass
does, and the flow is most of both. A melting macro step, eight such parts each run by both passes 1.8 times on average
(fact 126), is then about 8 x 1.8 x 239 s, an hour, against about 8 hours by 6c's code on this state (6c's own estimate,
2.7 h, took one heat sub-step for all, and its transport grows dearer as the flows meet). Options for a cheaper trial
pass, measured, not built (the flow stays in it, spec 16, 2026-10-22): at a fixed two flow sub-steps a heat sub-step it
took 64 s, its base law 0.97 % of the throughput from the chosen count's (2 597 columns beyond 1e-3); at one, 50 s and
2.0 % (worse than no flow, one 10 ms flow sub-step being unresolved); at four, 122 s and 0.52 % (both under a load of
80). The chosen count is near two here already, so the trial pass's own count saves little on this state.

**Rulings** (2026-10-22, within spec 6.4-6.5, 7.1-7.2 and the rulings of 6a-6c; each is what, then why, then the cost
if wrong):
- *The large cycles iterated, the small factored, at 500 cells* — the LU fills only within cycles: one of 1 700 cells
  filled 1.1 times (0.03 s), ten thousand 5-6 times (0.2-25 s) — if wrong, a few hundredths of a second a flow sub-step
  where a cycle of 500-2 000 cells would have factored faster.
- *The tolerance a normwise backward error of 1e-13, at most 40 corrections, then the LU* — what 6c's remedy reached in
  2-9 corrections, above the round-off floor of the residual at high Courant numbers; the cells are installed from the
  fluxes, so the books do not depend on it — if wrong, a solve stopped loosely moves metal by its error (1.9e-4 of the
  largest mass at a tolerance of 1e-3), the books still exact.
- *The landed metal coasts at the velocity it arrived with* — the transport moves every cell at the velocity its
  momentum step left it, and the landed metal's momentum is its arrivals'; measured closer to the fine runs than the
  one-cell estimate at the sub-steps the adaptive count chooses — if wrong, a front of thin deposits over rigid metal
  runs ahead by the drag it would have felt in one flow sub-step (the strip's front at 5 ms 18 columns past the fine
  run's tail of 1e-3 of its deposit, its centroid 2 mm ahead).
- *The adaptive count from the last heat sub-step's largest Courant number and growth, the targets 6 and 0.02* — at
  Courant 6 the runoff is within 2.5 % of the finest run and its kinetic energy 5 %, the linearisations at 4, where 7.3
  reached their cap; dt s = 0.02 keeps backward Euler's growth rate within about 1 % — if wrong, a part whose flow
  speeds up within a heat sub-step takes one flow sub-step beyond the target before the next count follows (the first
  from rest is 6c's ten).
- *At most twice the last flow sub-step, from 10 ms to 0.25 ms* — a flow sub-step from rest or a slow flow grows over a
  few heat sub-steps (10 to 1 in four), and 10 ms flow sub-steps are only taken where the flow is slow (the thickness
  step is unstable at 10 ms in a fast one); 0.25 ms was 6c's finest — if wrong, a flow needing more than 40 flow
  sub-steps runs above the Courant target (the history's `flow_courant_max` and `flow_cycles_max` say so).
- *The real pass takes its trial pass's counts; the memory is part of the skins' state* — so that the base law is made
  by the pass it predicts (required in writing this task); the memory rolls back with the cells on a trial pass, a
  damped repeat and a split part — if wrong, nothing: a real pass's counts could otherwise differ only where a Courant
  number crosses a rounding of the count.
- *The momentum step on the columns with a moving cell* — the others' cells are anchored and end at rest, as the full
  solve gave, to the last bit — if wrong, nothing.

**Measured against the tests** (2026-10-22): the iterative solve on the mixed sphere skins' cycles (up to 35 cells,
every cycle iterated) gives the direct solve's state to 4.5e-14 of each array's largest value in 44 corrections, the
column books to 4e-16 and the total to 2e-18 of the mass; stopped at a tolerance of 1e-3 its masses are 1.9e-4 off with
the books as exact and every mass positive; with no corrections allowed both systems fall back to the LU, to the last
bit; the landed metal and the strip's runoff as above; `flow_cycles` against its targets, floor and cap; the trial and
real passes' counts on a strip and in the skin body; the momentum step on the moving columns to the last bit.

- [ ] **Step 1: The tests** — in `tests/test_reentry_model_skin.py` (6a's arrival test now lets the landed metal flow
  on; its upwind-order test's stand-in order returns the cycles' sizes):

```diff
--- a/tests/test_reentry_model_skin.py
+++ b/tests/test_reentry_model_skin.py
@@ -950,8 +950,10 @@
 
 def test_metal_arriving_on_a_rigid_top_lands_on_it_and_the_receiver_s_own_cells_are_untouched():
     """Spec 6.5's arrival rule: where the receiver's top is rigid -- anchored solid, or a crust over flowing metal --
-    the arriving metal lands on top of its top cell, as a new cell holding exactly what arrived; nothing passes
-    through rigid metal."""
+    the arriving metal lands on top of its top cell, as a new cell holding exactly what stayed there; nothing passes
+    through rigid metal. Since Task 6d the landed metal flows on within the same sub-step, at the velocity it arrived
+    with, onto the rigid tops beyond: columns the donor does not touch gain a new top cell too, and every new top cell
+    holds the donor's metal at its specific enthalpy and velocity."""
     col, _ = strip_grid(3, 3)
     temps = [[RIGID] * 12 if k % 2 else [RIGID] * 2 + [LIQUID] * 6 + [RIGID] * 4 for k in range(col.n)]
     temps[4] = [LIQUID] * 8 + [RIGID] * 4
@@ -963,18 +965,21 @@
     for r in range(col.n):
         if r == 4:
             continue
-        if r in receivers:
+        if book.mass[r] > 0.0:
             assert sk.count[r] == before[4][r] + 1
             assert sk.m[r, 0] == pytest.approx(book.mass[r], rel=1e-12)
             assert sk.H[r, 0] == pytest.approx(book.mass[r] * h_donor, rel=1e-12)
             np.testing.assert_allclose(sk.U[r, 0], [0.1, 0.04, 0.0], rtol=1e-12)
             shifted = slice(1, sk.count[r])
         else:
-            assert book.mass[r] == 0.0
+            assert book.mass[r] == 0.0 and sk.count[r] == before[4][r]
             shifted = slice(0, sk.count[r])
         for now, then in ((sk.m, before[1]), (sk.H, before[2]), (sk.T, before[3]), (sk.U, before[5])):
             assert np.array_equal(now[r, shifted], then[r, :before[4][r]])
+    assert all(book.mass[r] > 0.0 for r in receivers)
     assert any(r % 2 for r in receivers) and any(r % 2 == 0 for r in receivers)          # solid tops and crusts
+    beyond = [r for r in range(col.n) if r != 4 and r not in receivers and book.mass[r] > 0.0]
+    assert beyond                                                                     # it flowed on (Task 6d)
 
 
 def test_a_floating_crust_carried_across_lands_on_top_and_the_slurry_beneath_it_joins_by_fraction():
@@ -1105,7 +1110,7 @@
     second, col, _ = mixed_sphere_skins(coarse_sphere_mesh)
     a = skin.transport(first, col, 0.05, regrid=False)
     rng = np.random.default_rng(0)
-    monkeypatch.setattr(skin, "_upwind_order", lambda src, dst, n: (rng.permutation(n), False))
+    monkeypatch.setattr(skin, "_upwind_order", lambda src, dst, n: (rng.permutation(n), False, np.ones(n)))
     b = skin.transport(second, col, 0.05, regrid=False)
     for x, y in zip(first.state(), second.state()):
         np.testing.assert_allclose(x, y, rtol=1e-11, atol=1e-11 * np.abs(y).max())
@@ -1255,8 +1260,8 @@
     """Asha, 2026-10-21: liquid 0.29 mm deep on a strip beside a wall of rigid metal 0.4 mm high, and a crust over
     slurry in the middle. The wall's height enters the curvature of the liquid beside it, as a wall's would; neither the
     wall nor the crust feels a capillary push, and surface tension moves no metal into either (what lands there is the
-    round-off of a velocity the faces' push cancels, 1.5e-15 of what moves); the wall's cells stay as they were and
-    the crust and the slurry beneath it, pushed by nothing, stay at rest."""
+    round-off of a velocity the faces' push cancels, 1.5e-15 of what moves; 1.4e-15 with Task 6d's transport); the
+    wall's cells stay as they were and the crust and the slurry beneath it, pushed by nothing, stay at rest."""
     col, pts = strip_grid(21, 2)
     x = np.rint(pts[:, 0] / 1e-3).astype(int)
     temps = [[RIGID] * 40 if k >= 15 else [RIGID] * 2 + [SLURRY] * 20 + [RIGID] * 5 if k == 7
@@ -1417,7 +1422,7 @@
     that, to the last bit, and books every impulse: the columns' momentum changes by the push's (the thickness step's
     included), the gas shear's, the walls' (the anchored metal's, the thickness step's and the transport's projections)
     and what the transport carried, to round-off (2.1e-12 of the largest impulse at one cycle and 9.4e-13 at three,
-    measured 2026-10-22, as 6b's books on these skins). With `cycles` the sub-step runs as that many flow sub-steps,
+    measured 2026-10-22, as 6b's books on these skins; 9.6e-13 at three with Task 6d's transport). With `cycles` the sub-step runs as that many flow sub-steps,
     the cells split and merged between them; the heat sub-step's own regrid follows the last."""
     from test_reentry_model_columns import sphere_columns
     col, _ = sphere_columns(coarse_sphere_mesh)
@@ -1540,7 +1545,7 @@
     """A pass of 0.1 s with the conduction on, liquid removed as it forms and the flow on: liquid aluminium at 909 K, just
     below T_feed, on a strip heated at 30 MW/m2, pushed along it and sheared: the mass and the energy close as without
     flow, and each column's momentum changes by the flow's impulses less what the removed liquid took with it, to
-    round-off (4.9e-16 of the largest impulse, measured 2026-10-22)."""
+    round-off (4.9e-16 of the largest impulse, measured 2026-10-22; 1.9e-16 with Task 6d's transport)."""
     sk, col, _ = strip_layer(nx=12, T=909.0)
     sk.settings = skin.SkinSettings(**CELLS_04, substep=0.01)
     T_b = np.full(sk.n, RIGID)
@@ -1556,3 +1561,180 @@
     impulse = f.push + f.shear + f.wall + f.carried - book.momentum_out
     scale = max(np.abs(x).max() for x in (f.push, f.shear, f.wall, book.momentum_out))
     assert np.abs(sk.column_momentum() - P0 - impulse).max() < 1e-13 * scale
+
+
+# -- making the flow affordable (sub-plan 21 Task 6d) ------------------------------------------------------------------
+def transported(coarse_sphere_mesh, monkeypatch, dt=0.05, **settings):
+    """The mixed sphere skins after one transport sub-step of dt at the given module settings, with the column masses
+    before it, the book and the solves' counts."""
+    for key, value in settings.items():
+        monkeypatch.setattr(skin, key, value)
+    sk, col, _ = mixed_sphere_skins(coarse_sphere_mesh)
+    M0, counts = sk.column_mass(), dict(skin.TRANSPORT)
+    book = skin.transport(sk, col, dt, regrid=False)
+    return sk, M0, book, {k: skin.TRANSPORT[k] - counts[k] for k in counts}
+
+
+def test_the_iterative_transport_solve_gives_the_direct_solve_s_answer_where_flows_meet(coarse_sphere_mesh, monkeypatch):
+    """The mixed sphere skins' random velocities make flows meet in cycles (of up to 35 cells here) at a Courant number of
+    4.4. Iterating every cycle (TRANSPORT_CYCLE 1) gives what factoring them all gives to 4.5e-14 of each array's largest
+    value, in 44 corrections over the two systems (measured 2026-10-22); the books close exactly by construction -- each
+    column's mass changes by its book to 4e-16 of the largest column, the books sum to 2e-18 of the mass -- and every
+    mass is positive."""
+    direct, M0, _, d = transported(coarse_sphere_mesh, monkeypatch, TRANSPORT_CYCLE=10 ** 9)
+    monkeypatch.undo()
+    sk, M0, book, c = transported(coarse_sphere_mesh, monkeypatch, TRANSPORT_CYCLE=1)
+    assert d["corrections"] == 0 and c["corrections"] > 2 and c["fallbacks"] == 0 and book.courant > 4.0
+    for x, y in zip(sk.state(), direct.state()):
+        np.testing.assert_allclose(x, y, rtol=0, atol=1e-13 * np.abs(y).max())
+    np.testing.assert_allclose(sk.column_mass() - M0, book.mass, rtol=0, atol=1e-15 * M0.max())
+    assert abs(book.mass.sum()) < 1e-16 * M0.sum() and np.all(sk.m[sk.valid] > 0.0)
+
+
+def test_the_books_close_exactly_and_the_masses_stay_positive_however_loosely_the_iteration_stops(coarse_sphere_mesh,
+                                                                                                monkeypatch):
+    """Each cell is installed from the solve's fluxes -- what it had, less what it sent, plus what arrived -- so the
+    books close whatever the iteration's error, and the iteration rises to the masses from below (a regular splitting of
+    an M-matrix), so they stay positive: stopped at 1e-3 of the right-hand side, the masses are 1.9e-4 of the largest
+    off the direct solve's, while the books close to 4e-16 and every mass is positive (measured 2026-10-22)."""
+    direct, _, _, _ = transported(coarse_sphere_mesh, monkeypatch, TRANSPORT_CYCLE=10 ** 9)
+    monkeypatch.undo()
+    sk, M0, book, c = transported(coarse_sphere_mesh, monkeypatch, TRANSPORT_CYCLE=1, TRANSPORT_TOL=1e-3)
+    off = np.abs(sk.m - direct.m).max() / direct.m.max()
+    assert 1e-6 < off < 2e-3 and c["fallbacks"] == 0
+    np.testing.assert_allclose(sk.column_mass() - M0, book.mass, rtol=0, atol=1e-15 * M0.max())
+    assert abs(book.mass.sum()) < 1e-16 * M0.sum() and np.all(sk.m[sk.valid] > 0.0)
+
+
+def test_an_iteration_that_cannot_finish_falls_back_to_the_direct_solve(coarse_sphere_mesh, monkeypatch):
+    """With no corrections allowed the iteration cannot reach its tolerance where cycles are iterated, and both systems
+    fall back to the LU, which gives what factoring every cycle gives, to the last bit."""
+    direct, _, _, _ = transported(coarse_sphere_mesh, monkeypatch, TRANSPORT_CYCLE=10 ** 9)
+    monkeypatch.undo()
+    sk, _, _, c = transported(coarse_sphere_mesh, monkeypatch, TRANSPORT_CYCLE=1, TRANSPORT_CORRECTIONS=0)
+    assert c["fallbacks"] == 2 and all(np.array_equal(x, y) for x, y in zip(sk.state(), direct.state()))
+
+
+def test_metal_landing_on_rigid_tops_flows_on_within_the_sub_step_at_the_velocity_it_arrived_with():
+    """Liquid on the first five columns of a strip of 1 mm columns moving at 0.5 m/s toward rigid metal, one 10 ms
+    sub-step (a Courant number of 5 along the strip): the metal landing on the first rigid top flows on in the same
+    solve (Task 6d; 6a's rule held it there until the next sub-step), over every rigid column to the strip's end, each
+    holding 0.831-0.836 of the last -- the one-dimensional implicit upwind scheme's nu / (1 + nu) = 5/6 for a pulse
+    entering columns at rest -- at the velocity it arrived with, to 2e-16; the books close and no rigid cell changes
+    (measured 2026-10-22)."""
+    col, pts = strip_grid(20, 2)
+    x = np.rint(pts[:, 0] / 1e-3).astype(int)
+    sk = stacked(col, [[LIQUID] * 10 + [RIGID] * 5 if xi < 5 else [RIGID] * 15 for xi in x])
+    u = np.array([0.5, 0.0, 0.0])
+    sk.U = np.where((sk.valid & (sk.T > 900.0))[..., None], u, 0.0)
+    before, M0 = sk.state(), sk.column_mass()
+    book = skin.transport(sk, col, 0.01, regrid=False)
+    per = np.array([book.mass[x == k].sum() for k in range(20)])
+    assert np.all(per[5:] > 0.0)                                         # 15 columns in one sub-step
+    np.testing.assert_allclose(per[6:13] / per[5:12], 5.0 / 6.0, rtol=0.01)
+    new = sk.count > before[4]
+    assert np.array_equal(new, x >= 5) and np.abs(sk.U[new, 0] - u).max() < 1e-12
+    assert abs(book.mass.sum()) < 1e-16 * M0.sum()
+    for c in np.flatnonzero(new):
+        for now, then in ((sk.m, before[1]), (sk.H, before[2]), (sk.T, before[3])):
+            assert np.array_equal(now[c, 1:sk.count[c]], then[c, :before[4][c]])
+
+
+def strip_runoff(cycles, span=0.02, tau=20.0, edge=10):
+    """Liquid aluminium 0.2 mm deep (20 cells at 950 K over five rigid ones) on the first `edge` columns of a strip of
+    1 mm columns two corners wide, rigid metal (25 cells) beyond, sheared along it by tau Pa for `span` s, the
+    conduction off, `cycles` flow sub-steps in each 10 ms heat sub-step: the volume [m3] each column beyond the edge
+    gained."""
+    col, pts = strip_grid(40, 2)
+    x = np.rint(pts[:, 0] / 1e-3).astype(int)
+    sk = stacked(col, [[LIQUID] * 20 + [RIGID] * 5 if xi < edge else [RIGID] * 25 for xi in x])
+    sk.settings = skin.SkinSettings(**CELLS_04, substep=0.01)
+    T_b, H0 = np.full(sk.n, RIGID), sk.thickness()
+    skin.advance(sk, span, np.zeros(sk.n), sk.area, 0.0, 0.0, T_b, np.zeros(sk.n), T_b1=T_b, remove=False,
+                 flow=strip_flow(col, dp_ds=0.0, tau=tau, cycles=cycles))
+    dV = (sk.thickness() - H0) * col.area
+    return np.array([dV[x == k].sum() for k in range(edge, 40)])
+
+
+def test_a_runoff_front_crosses_rigid_metal_faster_than_a_column_a_flow_sub_step_and_converges_to_a_fine_one(monkeypatch):
+    """The front of a sheared liquid layer running onto rigid metal (`strip_runoff`, 20 ms, the fastest cell about
+    1.4 m/s): with 6a's arrival rule the front advanced at most one column a flow sub-step, and with flow sub-steps of
+    5 ms reached 3 of the 10 columns the 0.25 ms run's covers; now metal flows on over the rigid tops within the
+    sub-step. Against the 0.25 ms run, where the metal went beyond the edge (L1) is 0.435, 0.135 and 0.068 off at 5, 1
+    and 0.5 ms (1.25, 0.195 and 0.10 with 6a's rule), first order in the sub-step, and the share beyond four columns
+    from the edge -- out of reach in four flow sub-steps of 5 ms under 6a's rule -- is there at 5 ms (measured
+    2026-10-22)."""
+    monkeypatch.setattr(skin, "_conduct", no_conduction)
+    ref = strip_runoff(40)
+    err = {c: np.abs(strip_runoff(c) - ref).sum() / np.abs(ref).sum() for c in (2, 10, 20)}
+    coarse = strip_runoff(2)
+    assert coarse[5:].sum() > 0.5 * ref[5:].sum() > 0.0
+    assert err[2] < 0.6 and err[10] < 0.2 and err[20] < 0.6 * err[10]
+
+
+def test_the_adaptive_count_meets_its_targets_within_its_floor_and_cap():
+    """`flow_cycles` (Task 6d): from no memory, FLOW_CYCLES; the next flow sub-step meets the Courant target in
+    proportion to the sub-step and the growth target in proportion to its square, grows by at most FLOW_RELAX, and the
+    count stays within its floor and cap; the same memory gives the same count."""
+    dt, C, g = 0.01, skin.FLOW_COURANT, skin.FLOW_GROWTH
+    assert skin.flow_cycles(np.zeros(3), dt) == skin.FLOW_CYCLES
+    assert skin.flow_cycles([2.0 * C, 0.0, 1e-3], dt) == 20                  # twice the target: half the sub-step
+    assert skin.flow_cycles([0.5 * C, 0.0, 1e-3], dt) == 5                   # half: twice it, the most it may grow
+    assert skin.flow_cycles([0.1 * C, 0.0, 1e-3], dt) == 5
+    assert skin.flow_cycles([0.0, 4.0 * g, 1e-3], dt) == 20                  # four times the growth: half the sub-step
+    assert skin.flow_cycles([0.0, 0.0, 1e-2], dt) == skin.FLOW_CYCLES_MIN == 1
+    assert skin.flow_cycles([1e4 * C, 0.0, 1e-3], dt) == skin.FLOW_CYCLES_MAX
+    assert skin.flow_cycles([1.3 * C, 0.3 * g, 2e-3], dt) == skin.flow_cycles([1.3 * C, 0.3 * g, 2e-3], dt) == 7
+
+
+def test_the_real_pass_takes_the_trial_pass_s_counts_and_the_memory_is_part_of_the_skins_state(monkeypatch):
+    """A bump on liquid aluminium on flat columns, its flow's count chosen each heat sub-step: a still layer under
+    surface tension relaxes from FLOW_CYCLES to the floor (10, 5, 3, 2, 1, ...), one sheared at 50 Pa and pushed at
+    20 kPa/m settles at 6 (10, 7, 9, 8, 7, 6, ..., measured 2026-10-22). The trial pass leaves the skins and their memory
+    as they were and its counts in the flow's plan; the real pass after it, at a base 5 K warmer, takes the same and
+    then forgets the plan, and the skins remember its last heat sub-step."""
+    monkeypatch.setattr(skin, "_conduct", no_conduction)
+    seen = []
+    substeps = skin.flow_substeps
+
+    def recorded(sk, flow, dt, cycles=None):
+        book = substeps(sk, flow, dt, cycles)
+        seen.append(book.cycles)
+        return book
+    monkeypatch.setattr(skin, "flow_substeps", recorded)
+    for tau, dp_ds, settled in ((0.0, 0.0, 1), (50.0, -2e4, 6)):
+        sk, col, _ = strip_layer(bump=0.2)
+        sk.settings = skin.SkinSettings(**CELLS_04, substep=0.01)
+        T_b = np.full(sk.n, RIGID)
+        flow = strip_flow(col, dp_ds=dp_ds, tau=tau, cycles=None)
+        before = sk.state()
+        seen.clear()
+        skin.advance(sk, 0.1, np.zeros(sk.n), sk.area, 0.0, 0.0, T_b, np.zeros(sk.n), remove=False, flow=flow)
+        trial = list(seen)
+        assert all(np.array_equal(x, y) for x, y in zip(sk.state(), before)) and flow.plan == trial
+        assert trial[0] == skin.FLOW_CYCLES and trial[-1] == settled and len(trial) == 10
+        seen.clear()
+        skin.advance(sk, 0.1, np.zeros(sk.n), sk.area, 0.0, 0.0, T_b, np.zeros(sk.n), T_b1=T_b + 5.0, remove=False,
+                     flow=flow)
+        assert seen == trial and flow.plan is None
+        assert sk.flow_memory[2] == pytest.approx(0.01 / settled) and sk.flow_memory[0] > 0.0
+
+
+def test_the_momentum_step_on_the_moving_columns_alone_is_the_solve_over_every_column_to_the_last_bit(
+        coarse_sphere_mesh):
+    """Task 6d: the momentum step solves only the columns where something moves (6c's solved every column). On mixed
+    sphere skins with two fifths of the columns cold, 1 354 of 2 254 columns solved, every velocity, impulse and
+    response is the solve over every column's, to the last bit, in the same five passes."""
+    sk, col, rng = mixed_sphere_skins(coarse_sphere_mesh)
+    cold = rng.uniform(size=col.n) < 0.4
+    sk.T = np.where(cold[:, None] & sk.valid, 870.0, sk.T)
+    sk.H = sk.m * sk.material.enthalpy(sk.T)
+    other = skin.Skins(sk.material, sk.area, sk.settings)
+    other.set_state(sk.state())
+    G = np.array([3e4, -1e4, 2e3])
+    a = skin.momentum_step(sk, 0.01, G, np.zeros(3), sk.area, response=True)
+    b = skin._momentum_columns(other, 0.01, np.broadcast_to(G, (col.n, 3)), np.zeros((col.n, 3)), other.area, True)
+    flowing = sk.valid & np.isfinite(skin.cell_viscosity(sk, 0.0))
+    solved = (sk.valid & ~skin.anchored_cells(sk, flowing)).any(axis=1)
+    assert 0 < solved.sum() < col.n and np.array_equal(sk.U, other.U) and a.passes == b.passes
+    assert all(np.array_equal(getattr(a, k), getattr(b, k)) for k in ("push", "shear", "wall", "response", "wall_response"))
```

- [ ] **Step 2: The skin body's tests** — in `tests/test_reentry_model_skinbody.py`:

```diff
--- a/tests/test_reentry_model_skinbody.py
+++ b/tests/test_reentry_model_skinbody.py
@@ -428,6 +428,7 @@
     H = sk.thickness()
     assert stats["roughness_max"] == float((np.abs(H[i] - H[j]) / b.columns.edge_length).max()) > 0.0
     assert stats["flow_courant_max"] > 0.0 and 0.0 <= stats["flow_growth_max"] < 1.0 and stats["flow_passes_max"] >= 1
+    assert stats["flow_cycles_max"] == 1
     assert stats["wall_flow_s"] > 0.0 and stats["flowing_depth_max_mm"] > 0.0
 
 
@@ -456,7 +457,7 @@
     normal, f = seen[0]
     np.testing.assert_allclose(f.shear, skin.along_meridian(normal, b.v_hat, first.tau), rtol=1e-12, atol=1e-12)
     np.testing.assert_allclose(f.g, first.deceleration * b.v_hat, rtol=1e-12)
-    assert f.rho_l == b.skins.rho_l and f.sigma == b.material.liquid.sigma and f.cycles == skin.FLOW_CYCLES
+    assert f.rho_l == b.skins.rho_l and f.sigma == b.material.liquid.sigma and f.cycles is None    # chosen (Task 6d)
     second = evaluate(states[1])
     seen.clear()
     b._hint = 2                                                          # two parts, to see the interpolation
@@ -465,3 +466,37 @@
     for frac in (0.25, 0.75):                                           # the two parts' middles
         dp = first.dp_ds + (second.dp_ds - first.dp_ds) * frac
         assert any(np.allclose(f.dp_ds, dp, rtol=1e-12, atol=1e-9) for _, f in seen)
+
+
+def test_a_skin_body_s_flow_takes_the_counts_each_trial_pass_chose_in_its_real_pass(coarse_sphere_mesh, monkeypatch):
+    """Sub-plan 21 Task 6d: with the flow on, the flow sub-steps in each heat sub-step are chosen from the last one's
+    Courant number and growth (the default), and the real pass of every part takes the counts its trial pass chose, so
+    that the base law is made by the pass it predicts; the history reports the most a heat sub-step took. 3 MW/m2 on the
+    front of the 4 mm test sphere from 880 K for a 0.25 s step, as above: the books close as with a fixed count."""
+    calls = []
+    advance, substeps = skin.advance, skin.flow_substeps
+
+    def recorded_advance(sk, dt, *a, **k):
+        calls.append(("trial" if k.get("T_b1") is None else "real", []))
+        return advance(sk, dt, *a, **k)
+
+    def recorded_substeps(sk, flow, dt, cycles=None):
+        book = substeps(sk, flow, dt, cycles)
+        calls[-1][1].append((cycles, book.cycles))                    # the count handed down, and the one used
+        return book
+    monkeypatch.setattr(skin, "advance", recorded_advance)
+    monkeypatch.setattr(skin, "flow_substeps", recorded_substeps)
+    b = make(coarse_sphere_mesh, flow=True, remove=False)
+    assert b.flow_cycles is None
+    M0 = b.mass()
+    b.advance(0.25, 0.25, front_loads(b, 3e6), forcing=newtonian(b))
+    kinds = [k for k, _ in calls]
+    assert kinds == ["trial", "real"] * (len(calls) // 2) and len(calls) >= 2
+    chosen = [[c for c, _ in cs] for _, cs in calls]
+    assert all(chosen[i] == chosen[i + 1] for i in range(0, len(calls), 2))
+    used = [u for _, cs in calls for _, u in cs if u]                 # where something flowed
+    assert used and max(used) <= skin.FLOW_CYCLES_MAX and min(used) >= skin.FLOW_CYCLES_MIN
+    stats = b.skin_stats()
+    assert max(u for _, u in calls[-1][1]) <= stats["flow_cycles_max"] <= max(used) and stats["flow_courant_max"] > 0.0
+    assert b.mass() == pytest.approx(M0, rel=1e-13) and abs(b.energy_balance_residual()) < 1e-8
+    assert stats["momentum_residual"] < 1e-12
```

- [ ] **Step 3: The command line's tests** — in `tests/test_reentry_model_cli.py`:

```diff
--- a/tests/test_reentry_model_cli.py
+++ b/tests/test_reentry_model_cli.py
@@ -379,7 +379,7 @@
     assert s["surface_model"] == "skin" and s["prism_layers"] == 0
     assert s["skin"] == {"min_thickness_mm": 0.8, "fine_zone_mm": 0.4, "cell_um": 20.0, "growth": 1.1, "substep_ms": 10.0, "preset": "dev",
                          "deform_depth_mm": 15.0, "n_columns": r["skin_columns"], "fine_band_mm": 4.0, "fine_size_mm": 2.0,
-                         "swell": True, "flow": "off", "flow_cycles": cli.skin.FLOW_CYCLES}
+                         "swell": True, "flow": "off", "flow_cycles": "adaptive"}
     assert r["skin_columns"] == s["n_patches"] and abs(r["energy_balance_residual"]) < 1e-8
     rows = list(csv.DictReader(open(path[:-5] + ".csv")))
     assert all(k in rows[0] for k in coupled.SKINBODY_COLUMNS) and float(rows[-1]["skin_cells"]) > 0.0
@@ -428,11 +428,13 @@
 # Sub-plan 21 Task 6c: the flow in the skins, its setting and its sub-cycles
 
 def test_the_flow_settings_reach_the_run_name_only_when_not_default():
+    """Task 6d: the flow sub-steps are chosen each skin sub-step by default, so any fixed count is named."""
     from reentry_model import skin
     suffix = lambda **k: cli.skin_name_suffix(0.8, None, 1.10, 10.0, "production", 15.0, **k)
     assert suffix() == suffix(flow=False) == ""
-    assert suffix(flow=True) == suffix(flow=True, flow_cycles=skin.FLOW_CYCLES) == "_skinflow-on"
-    assert suffix(flow=True, flow_cycles=skin.FLOW_CYCLES + 1) == "_skinflow-on_flowcycles{}".format(skin.FLOW_CYCLES + 1)
+    assert suffix(flow=True) == suffix(flow=True, flow_cycles=None) == "_skinflow-on"
+    assert suffix(flow=True, flow_cycles=skin.FLOW_CYCLES) == "_skinflow-on_flowcycles{}".format(skin.FLOW_CYCLES)
+    assert suffix(flow=True, flow_cycles=2) == "_skinflow-on_flowcycles2"
     assert cli.skin_name_suffix(0.8, None, 1.10, 10.0, "dev", 15.0, swell=False, flow=True) == "_skin-dev_noswell_skinflow-on"
     args = cli.build_parser().parse_args(["run", "--diameter", "100", "--velocity", "7.5", "--altitude", "77.5"])
     assert args.skin_flow == "off" and args.skin_flow_cycles is None
```

- [ ] **Step 4: Run them to see them fail.** `"$PY" -m pytest tests/test_reentry_model_skin.py
  tests/test_reentry_model_skinbody.py tests/test_reentry_model_cli.py -q`. Expected: 15 failed, 148 passed (in
  1 135 s):
  - `assert []`, 6a's arrival test: no column the donor does not touch gained metal (it does not flow on yet);
  - `ValueError: too many values to unpack (expected 2)`, 6a's upwind-order test: the stand-in order returns the
    cycles' sizes, which 6c's transport does not take;
  - `AttributeError: <module 'reentry_model.skin' ...> has no attribute 'TRANSPORT_CYCLE'` (3), the iterative solve's
    tests;
  - `assert np.False_`, the landed metal reaching all 15 rigid columns, and `assert 1.25e-24 > (0.5 * 5.30e-10)`, the
    strip's runoff beyond four columns at 5 ms;
  - `AttributeError: module 'reentry_model.skin' has no attribute 'FLOW_COURANT'`, `TypeError: flow_substeps() takes 3
    positional arguments but 4 were given` and `AttributeError: module 'reentry_model.skin' has no attribute
    '_momentum_columns'`;
  - `KeyError: 'flow_cycles_max'` and `assert 10 is None` (2), the skin body's;
  - the skin run's settings (`'flow_cycles': 10` against `'adaptive'`) and `assert '_skinflow-on' ==
    '_skinflow-on_flowcycles10'`, the command line's.

- [ ] **Step 5: The iterative solve, the landed metal flowing on, the momentum step on the moving columns and the
  adaptive count** in `reentry_model/skin.py`:

```diff
--- a/reentry_model/skin.py
+++ b/reentry_model/skin.py
@@ -151,7 +151,10 @@
     is an intensive per-cell field every cell operation carries (the trial pass's dT/dT_b). `U` (n, K, 3) is every
     cell's velocity [m/s] (sub-plan 21 Task 2): part of the cells' state, carried by every operation with its momentum
     exact -- split halves keep it, merged cells take the mass-weighted mean, metal drawn from the 3D model arrives at
-    rest -- and `momentum_out` (n, 3) [kg m/s] is what the last `take_top` or `take_bottom` took away with the metal."""
+    rest -- and `momentum_out` (n, 3) [kg m/s] is what the last `take_top` or `take_bottom` took away with the metal.
+    `flow_memory` (3,) is what the adaptive count of flow sub-steps remembers of the last heat sub-step's flow (sub-plan
+    21 Task 6d, `flow_cycles`): its largest Courant number and Rayleigh-Taylor growth and its flow sub-step's length [s],
+    zero before any flow; part of the state, so that a trial pass, a repeated or a shortened part starts from it."""
 
     def __init__(self, material, area, settings=None):
         self.material = material
@@ -163,6 +166,7 @@
         self.S = None
         self.U = np.zeros((n, 1, 3))
         self.momentum_out = np.zeros((n, 3))
+        self.flow_memory = np.zeros(3)
         self.rho_l = material.liquid.rho if material.liquid is not None and self.settings.swell else material.rho
         self.props = PropertyTable(material)
 
@@ -192,10 +196,23 @@
         return self.column_momentum().sum(axis=0)
 
     def state(self):
-        return tuple(a.copy() for a in (self.area, self.m, self.H, self.T, self.count, self.U))
+        return tuple(a.copy() for a in (self.area, self.m, self.H, self.T, self.count, self.U, self.flow_memory))
+
+    def rows(self, idx):
+        """Skins holding copies of the columns idx alone, with this one's material, settings and property tables: for
+        a solve that needs only some columns (sub-plan 21 Task 6d)."""
+        sub = object.__new__(Skins)
+        sub.material, sub.settings, sub.rho_l, sub.props = self.material, self.settings, self.rho_l, self.props
+        sub.area, sub.m, sub.H, sub.T, sub.count, sub.U = (a[idx] for a in (self.area, self.m, self.H, self.T, self.count,
+                                                                           self.U))
+        sub.S = None if self.S is None else self.S[idx]
+        sub.momentum_out, sub.flow_memory = np.zeros((len(sub.area), 3)), self.flow_memory.copy()
+        return sub
 
     def set_state(self, state):
-        self.area, self.m, self.H, self.T, self.count, self.U = (a.copy() for a in state)
+        """Restore a `state()`; one of the six arrays before Task 6d's memory starts the flow with none."""
+        self.area, self.m, self.H, self.T, self.count, self.U = (a.copy() for a in state[:6])
+        self.flow_memory = np.array(state[6], dtype=float) if len(state) > 6 else np.zeros(3)
         self.S = None
 
     # -- geometry of the cells --------------------------------------------------------------------------------------
@@ -650,7 +667,9 @@
     matrix its diagonal dominance; what the clip misses is the interface mismatch's to book). Otherwise the real pass,
     base held at T_b1, the 3D model's end-of-step base temperature; returns a SkinBook. The trial pass runs the flow too
     (spec 7.2's sub-steps are spec 7.1's: heat, flow, sideways exchange and spray), the sensitivity S carried with the
-    metal the transport moves."""
+    metal the transport moves. Where the flow's count is chosen (`flow.cycles` None, Task 6d), the trial pass leaves the
+    counts it chose in `flow.plan`, and the real pass after it takes the same, so that the base law is made by the pass
+    it predicts; a real pass with no trial before it chooses its own."""
     mat = sk.material
     N = max(1, int(math.ceil(dt / sk.settings.substep - 1e-9)))
     h = dt / N
@@ -667,7 +686,11 @@
     drawn_depth = np.zeros(sk.n)                      # m, how far down the 3D profile this pass has drawn
     melting = remove and mat.melts
     sideways = None
-    for _ in range(N):
+    plan = replay = None
+    if flow is not None and flow.cycles is None:                  # the trial pass's counts, for its real pass
+        replay = not trial and flow.plan is not None and len(flow.plan) == N
+        plan = list(flow.plan) if replay else []
+    for i in range(N):
         sideways = side(sk, sideways) if side is not None else None
         res = _conduct(sk, h, q_top, A_top, T_amb, eps, T_b, sideways)
         book.E_top += res[0]
@@ -677,7 +700,9 @@
         if trial:
             dEb += res[4]
         if flow is not None:
-            book.flow.add(flow_substeps(sk, flow, h))
+            if plan is not None and not replay:
+                plan.append(flow_cycles(sk.flow_memory, h))
+            book.flow.add(flow_substeps(sk, flow, h, cycles=plan[i] if plan is not None else None))
         if melting:
             k = np.cumprod(sk.valid & (sk.T >= mat.T_feed), axis=1).sum(axis=1)
             gone, gone_H = sk.take_top(k)
@@ -696,6 +721,8 @@
         book.drawn += short
         book.drawn_H += sk.add_bottom(short, h_b, s_b)
         sk.regrid()
+    if plan is not None:
+        flow.plan = plan if trial else None
     if trial:
         sk.set_state(saved)
         b = np.minimum(dEb / dt, 0.0)
@@ -885,15 +912,40 @@
     columns' momentum change to round-off. With `response`, the book also carries the velocities a unit push alone
     gives through the last pass's system and the wall's drag on them: with the viscosities held, the sub-step is linear
     in the push, so a change of the push changes every velocity by its response times the change (sub-plan 21 Task
-    6b)."""
+    6b).
+
+    Only the columns where something moves are solved (sub-plan 21 Task 6d): in the others every cell is anchored, so it
+    ends at rest and the momentum it held is the wall's, as the solve over every column gave -- to the last bit, each
+    column's system being its own (6c's step solved all 9 417 columns of the 100 mm sphere where half of them flowed)."""
+    n = sk.n
+    push = np.broadcast_to(np.asarray(push, dtype=float), (n, 3))
+    shear = np.broadcast_to(np.asarray(shear, dtype=float), (n, 3))
+    A_top = np.broadcast_to(np.asarray(A_top, dtype=float), (n,))
+    flowing = sk.valid & np.isfinite(cell_viscosity(sk, 0.0))
+    live = np.flatnonzero((sk.valid & ~anchored_cells(sk, flowing)).any(axis=1))
+    if len(live) == n:
+        return _momentum_columns(sk, dt, push, shear, A_top, response)
+    book = MomentumBook(np.zeros((n, 3)), np.zeros((n, 3)), -sk.column_momentum(), 0,
+                        np.zeros(sk.m.shape) if response else None, np.zeros(n) if response else None)
+    U = np.zeros_like(sk.U)
+    if len(live):
+        sub = sk.rows(live)
+        part = _momentum_columns(sub, dt, push[live], shear[live], A_top[live], response)
+        U[live] = sub.U
+        book.push[live], book.shear[live], book.wall[live], book.passes = part.push, part.shear, part.wall, part.passes
+        if response:
+            book.response[live], book.wall_response[live] = part.response, part.wall_response
+    sk.U = U
+    return book
+
+
+def _momentum_columns(sk, dt, push, shear, A_top, response):
+    """`momentum_step` on every column of sk, push (n, 3), shear (n, 3) and A_top (n,) given whole."""
     n, K = sk.m.shape
     valid = sk.valid
     rows = np.arange(n)
     last = np.maximum(sk.count - 1, 0)
     dz, A, m = sk.dz(), sk.area, sk.m
-    push = np.broadcast_to(np.asarray(push, dtype=float), (n, 3))
-    shear = np.broadcast_to(np.asarray(shear, dtype=float), (n, 3))
-    A_top = np.broadcast_to(np.asarray(A_top, dtype=float), (n,))
     flowing = valid & np.isfinite(cell_viscosity(sk, 0.0))
     anchored = anchored_cells(sk, flowing)
     moving = valid & ~anchored                                            # the flowing cells and the floating crust
@@ -1030,6 +1082,15 @@
 # -- sideways transport (sub-plan 21 Task 6a) ---------------------------------------------------------------------------
 SLIVER = 1.0e-12                  # a piece of a donor cell shorter than this share of it is the round-off of a boundary
                                   # the donor's and the receiver's partitions share, and is left out
+TRANSPORT_TOL = 1.0e-13           # the iterative solve stops once no row's residual exceeds this share of the largest
+                                  # |b| + |A| |x| of its right-hand side's column (a normwise backward error; Task 6d)
+TRANSPORT_CORRECTIONS = 40        # corrections of the iterative solve at most, before the direct LU takes over
+TRANSPORT_CYCLE = 500             # cells: a cycle of meeting flows larger than this is iterated, not factored (Task 6d;
+                                  # one of 1 700 cells filled its triangles 1.1 times, one of 10 000 5-6 times)
+TRANSPORT = {"solves": 0, "corrections": 0, "fallbacks": 0}   # the transport's systems solved, the iteration's
+                                  # corrections and the solves the direct LU finished, for measuring the cost
+LANDED_PASSES = 200               # passes at most of the estimate of the landed metal's velocities over the columns
+LANDED_TOL = 1.0e-6               # it stops once no estimated velocity moves by more than this share of the fastest
 
 
 @dataclass
@@ -1087,11 +1148,12 @@
 
 def _upwind_order(src, dst, n):
     """An order of n unknowns in which every donor (src) comes before its receivers (dst), except within a cycle of
-    transfers, which stays together; and whether it holds. The order is scipy's strong components, which its (Pearce's)
-    algorithm numbers in reverse topological order -- how scipy numbers them, not what it promises, so it is checked."""
+    transfers, which stays together; whether it holds; and the size of every unknown's cycle (its strong component; 1
+    outside any cycle). The order is scipy's strong components, which its (Pearce's) algorithm numbers in reverse
+    topological order -- how scipy numbers them, not what it promises, so it is checked."""
     graph = sp.csr_matrix((np.ones(len(src)), (src, dst)), shape=(n, n))
     _, label = connected_components(graph, directed=True, connection="strong")
-    return np.argsort(-label, kind="stable"), bool(np.all(label[src] >= label[dst]))
+    return np.argsort(-label, kind="stable"), bool(np.all(label[src] >= label[dst])), np.bincount(label)[label]
 
 
 def _factor(A, upwind):
@@ -1102,6 +1164,91 @@
     return spla.splu(A)
 
 
+def _solve(A, rhs, upwind, big=None):
+    """x with A x = rhs, rhs (N,) or (N, p), for the transport's systems (sub-plan 21 Task 6d). In upwind order A is
+    lower triangular but for the back-couplings of its cycles, where flows meet, and its LU fills only within them:
+    a cycle of a few cells costs nothing, but one where meeting flows join thousands of cells filled the triangles 5-6
+    times, and the LU then took most of a flow sub-step (6c). `big` (N,) marks the unknowns in such cycles, of more than
+    TRANSPORT_CYCLE cells: the back-couplings among them, U, are left out of the factored matrix M = A - U, which fills
+    only within the small cycles, and x_(k+1) = M^-1 (b - U x_k) from x_0 = M^-1 b, whose residual b - A x_(k+1) is
+    U (x_k - x_(k+1)), until no row's residual exceeds TRANSPORT_TOL of the largest |b| + |A| |x| in its column, at most
+    TRANSPORT_CORRECTIONS corrections; then the LU of A itself. Where no cycle is large, or `big` is None, M is A and x
+    its solve, as 6a's; where the order is not upwind, the LU in SuperLU's own order (`_factor`).
+
+    For the masses M is an M-matrix and N = -U >= 0 (a regular splitting), so M^-1 >= 0 and from x_0 >= 0 the iterates
+    rise monotonically to the solution, every residual b - A x being non-negative: the masses installed from an
+    iterate's fluxes, x + (b - A x), are positive (`transport`).
+
+    Measured 2026-10-22 on 6c's states of the 100 mm sphere, one flow sub-step of 1 ms, the mass's system and the
+    momentum's: with every column melted (one cycle of 18 379 cells) 6c's LU took 4.65 and 25.3 s, filling its
+    triangles 6.0 and 6.3 times, and this solve 0.16 and 0.25 s in 9 corrections each; on a quiescent melt (10 190
+    cells) 0.20 and 0.76 s against 0.10 and 0.28 s (3); on the heated melting state (11 109) 0.43 and 2.04 s against
+    0.14 and 0.21 s (5); driven at 93 s (1 700 cells, filling 1.1 times) 0.03 and 0.07 s against 0.07 and 0.10 s (4).
+    Each within 2e-13 of the LU's answer, the masses non-negative."""
+    TRANSPORT["solves"] += 1
+    if not upwind:
+        return _factor(A, upwind).solve(rhs)
+    U = sp.triu(A, k=1, format="coo")
+    keep = np.zeros(U.nnz, dtype=bool) if big is None else big[U.row] & big[U.col]
+    U = sp.csr_matrix((U.data[keep], (U.row[keep], U.col[keep])), shape=A.shape)
+    if U.nnz == 0:
+        return _factor(A, upwind).solve(rhs)
+    lu = _factor((A - U).tocsc(), upwind)
+    x = lu.solve(rhs)
+    Ux = U @ x
+    limit = TRANSPORT_TOL * (np.abs(rhs) + abs(A) @ np.abs(x)).max(axis=0)
+    r = -Ux
+    for k in range(TRANSPORT_CORRECTIONS + 1):
+        if np.all(np.abs(r).max(axis=0) <= limit):
+            TRANSPORT["corrections"] += k
+            return x
+        if k == TRANSPORT_CORRECTIONS:
+            break
+        x = lu.solve(rhs - Ux)
+        Ux, r = U @ x, Ux
+        r = r - Ux
+    TRANSPORT["fallbacks"] += 1
+    return _factor(A, upwind).solve(rhs)
+
+
+def _landed_velocity(columns, area, run, dt, base_m, base_p, stay):
+    """The velocity [m/s, (n, 3)] at which the metal landing on each column's rigid top flows on within the transport's
+    sub-step (sub-plan 21 Task 6d): the mass-weighted mean of what lands there, each arrival's velocity projected onto
+    the column's plane -- the landed metal coasts at the velocity it arrived with until the next momentum step drags it
+    on the rigid top -- zero where nothing lands, where the top flows (`run` > 0: the metal joins it by fraction) or
+    where a crust landed too (`stay`: rigid metal on a rigid top is anchored). What lands from the cells, base_m [kg,
+    (n,)] and base_p [kg m/s, (n, 3)] with its velocities projected, is known before the solve; what the landed metal
+    carries on onto further rigid tops is estimated column by column -- a column sends dt c m / (1 + dt sum c) of what
+    landed, c = (u . f)+ / A at its estimated velocity u -- pass by pass, a pass reaching one column further, until no
+    velocity moves by more than LANDED_TOL of the fastest, at most LANDED_PASSES passes; metal estimated at less than
+    SLIVER of the largest landing stays where it lands (the far end of the implicit scheme's geometric tail). These are
+    estimates for the rates alone: the solve carries the masses, the enthalpy and the momentum exactly."""
+    n = len(area)
+    i, j = columns.edges[:, 0], columns.edges[:, 1]
+    src, dst, f = np.r_[i, j], np.r_[j, i], np.r_[columns.face_vector, -columns.face_vector]
+    normal = columns.normal
+    can = (run == 0) & ~stay
+    onto = run[dst] == 0
+    u = np.zeros((n, 3))
+    m_est, p_est = base_m, base_p
+    for _ in range(LANDED_PASSES):
+        on = can & (m_est > SLIVER * m_est.max())
+        u_new = np.where(on[:, None], p_est / np.where(on, m_est, 1.0)[:, None], 0.0)
+        change = np.abs(u_new - u).max()
+        u = u_new
+        if change <= LANDED_TOL * np.abs(u).max():
+            break
+        d = np.flatnonzero(on[src])
+        rate = np.maximum(np.einsum("dx,dx->d", u[src[d]], f[d]), 0.0) / area[src[d]]
+        out = np.bincount(src[d], rate, n)
+        sent = dt * rate * m_est[src[d]] / (1.0 + dt * out[src[d]])
+        d, sent = d[onto[d]], sent[onto[d]]
+        v = u[src[d]] - np.einsum("dx,dx->d", u[src[d]], normal[dst[d]])[:, None] * normal[dst[d]]
+        m_est = base_m + np.bincount(dst[d], sent, n)
+        p_est = base_p + np.stack([np.bincount(dst[d], sent * v[:, x], n) for x in range(3)], axis=1)
+    return u
+
+
 def transport(sk, columns, dt, regrid=True, faces=None):
     """One linearly implicit upwind sub-step dt [s] of the metal moving sideways between neighbouring columns at the
     cells' velocities sk.U (spec 6.5; Asha, 2026-10-21: its own step, after the momentum step). Updates the cells and
@@ -1113,9 +1260,11 @@
     face's vector (`Columns.face_vector`) and m_k' the cell's mass at the sub-step's end; it carries the cell's specific
     enthalpy and velocity. Every cell is then m_k' (1 + dt sum c) = m_k + dt (what arrives): (I + dt C) m' = m, the
     film's runoff scheme (`film.Runoff`) on the cells -- an M-matrix whose columns sum to one, so the masses stay
-    positive, the books close to round-off and nothing limits dt. The enthalpy, the carried sensitivity S and the
-    momentum solve the same system, the momentum with every arrival projected onto its receiver's tangent plane, in two
-    components per cell (`tangent_frames`), the normal part booked to the wall (spec 6.4, section 14 assumption 2).
+    positive and nothing limits dt. The enthalpy, the carried sensitivity S and the momentum solve the same system, the
+    momentum with every arrival projected onto its receiver's tangent plane, in two components per cell
+    (`tangent_frames`), the normal part booked to the wall (spec 6.4, section 14 assumption 2). Every cell is installed
+    from the solve's fluxes -- what it had, less what it sent, plus what arrived, each arrival taking dt c times its
+    donor's solved value -- so the books close exactly whatever the solve's own error (Task 6d).
 
     Where metal joins the receiver (spec 6.5, section 14 assumption 3): a flowing donor cell's span of fractions of the
     donor's flowing depth (its flowing cells, top down) joins the receiver's cells between the same fractions of its
@@ -1125,9 +1274,18 @@
     substrate (`anchored_cells`) neither leaves nor receives, whatever velocity it holds: nothing passes through rigid
     metal. Only the cells that send or receive are solved for; with every velocity zero nothing changes.
 
+    Flowing metal that lands on a rigid top flows on in the same solve (Task 6d): its new top cell is a flowing donor
+    of the same system, sending across its column's faces at the velocity it arrived with (`_landed_velocity`: the
+    arrivals' mean, projected onto the column's plane, fixed before the solve so that the step stays linear), by
+    fraction into a flowing receiver's whole flowing layer and on top of a rigid one, whose new cell flows on in turn.
+    6a's rule held it there until the next sub-step, so a front crossing rigid metal advanced at most one column a
+    sub-step: 2 m/s on 2 mm columns at 1 ms, which set 6c's flow sub-step. Its momentum is the solve's, the arrivals'
+    projected; the drag of the rigid top beneath reaches it at the next momentum step. A new top cell that crust metal
+    reached stays where it lands (rigid metal on a rigid top is anchored).
+
     The cells are taken in upwind order (`_upwind_order`), in which the system is triangular but for its cycles where
-    flows meet, so its sparse LU has almost no fill. With `regrid` the cells are then split and merged (spec 6.2); a
-    sub-step that regrids anyway passes False.
+    flows meet, and solved by `_solve`: the small cycles factored, those of more than TRANSPORT_CYCLE cells iterated.
+    With `regrid` the cells are then split and merged (spec 6.2); a sub-step that regrids anyway passes False.
 
     Measured 2026-10-21 on the 100 mm sphere's skins (9 417 columns, 10 um cells), flowing along the meridian away from
     the stagnation point: 0.27 s a sub-step with the front hemisphere's top 0.3 mm flowing (124 000 of 553 000 cells
@@ -1182,6 +1340,42 @@
     e_row = np.r_[b[edge] * K + kr, n * K + b[te]]
     e_col = np.r_[a[edge] * K + kd, a[te] * K + tk]
     e_val = dt * np.r_[c[edge, kd] * share, c[te, tk]]
+    # metal landing on a rigid top flows on in the same solve (Task 6d): the new top cell sends across its column's
+    # faces at the velocity of what lands there, by fraction into a flowing receiver and on top of a rigid one
+    lands = flowing[a[te], tk] & (run[b[te]] == 0)
+    stay = np.zeros(n, dtype=bool)
+    stay[b[te[~flowing[a[te], tk]]]] = True                  # a crust landed: rigid metal on a rigid top stays
+    out = np.bincount(e_col, e_val, n * K)[a[te] * K + tk]
+    w = (e_val[len(edge):] * sk.m[a[te], tk] / (1.0 + out))[lands]   # what lands from each cell, estimated
+    rec = b[te[lands]]
+    v = sk.U[a[te], tk][lands]
+    v = v - np.einsum("dx,dx->d", v, columns.normal[rec])[:, None] * columns.normal[rec]
+    u_on = _landed_velocity(columns, sk.area, run, dt, np.bincount(rec, w, n),
+                            np.stack([np.bincount(rec, w * v[:, x], n) for x in range(3)], axis=1), stay)
+    s_all, r_all, f_all = np.r_[i, j], np.r_[j, i], np.r_[columns.face_vector, -columns.face_vector]
+    on = np.flatnonzero((u_on != 0.0).any(axis=1)[s_all])
+    c_on = np.maximum(np.einsum("dx,dx->d", u_on[s_all[on]], f_all[on]), 0.0) / sk.area[s_all[on]]
+    on, c_on = on[c_on > 0.0], c_on[c_on > 0.0]
+    reach = np.zeros(n, dtype=bool)                           # the new top cells metal reaches: only those send
+    reach[rec] = True
+    while len(on):
+        more = reach[s_all[on]] & (run[r_all[on]] == 0) & ~reach[r_all[on]]
+        if not more.any():
+            break
+        reach[r_all[on][more]] = True
+    on, c_on = on[reach[s_all[on]]], c_on[reach[s_all[on]]]
+    if len(on):
+        s, r = s_all[on], r_all[on]
+        into = np.flatnonzero(run[r] > 0)                     # by fraction of the receiver's flowing layer, whole
+        piece = np.diff(tau[r[into]], axis=1) * (np.arange(K)[None, :] < run[r[into]][:, None])
+        q, kq = np.nonzero(piece > 0.0)
+        onto = np.flatnonzero(run[r] == 0)                    # on top of a rigid receiver, its new top cell
+        e_edge = np.r_[e_edge, D + into[q], D + onto]
+        e_row = np.r_[e_row, r[into[q]] * K + kq, n * K + r[onto]]
+        e_col = np.r_[e_col, n * K + s[into[q]], n * K + s[onto]]
+        e_val = np.r_[e_val, dt * c_on[into[q]] * piece[q, kq], dt * c_on[onto]]
+        a, b = np.r_[a, s], np.r_[b, r]
+        D = len(a)
     # the cells that send or receive, in upwind order; each sends what its arrivals take, so the books close exactly
     active = np.zeros(n * K + n, dtype=bool)
     active[e_row], active[e_col] = True, True
@@ -1189,38 +1383,44 @@
     N = len(act)
     local = np.full(n * K + n, -1, dtype=np.int64)
     local[act] = np.arange(N)
-    order, upwind = _upwind_order(local[e_col], local[e_row], N)
+    order, upwind, size = _upwind_order(local[e_col], local[e_row], N)
+    big = size[order] > TRANSPORT_CYCLE                       # the unknowns in large cycles, in the solve's order
     at = np.empty(N, dtype=np.int64)
     at[order] = np.arange(N)
     to, frm = at[local[e_row]], at[local[e_col]]
-    diag = 1.0 + np.bincount(local[e_col], e_val, N)
-    # mass, enthalpy and the carried sensitivity
-    lu = _factor(sp.csc_matrix((np.r_[diag, -e_val], (np.r_[at, to], np.r_[at, frm])), shape=(N, N)), upwind)
+    lr, lc = local[e_row], local[e_col]
+    diag = 1.0 + np.bincount(lc, e_val, N)
+    # mass, enthalpy and the carried sensitivity, installed from the solve's fluxes: each cell keeps what it had, less
+    # what it sent, plus what arrived, so the books close exactly whatever the solve's own error (Task 6d)
     cell = lambda v: np.r_[v.reshape((n * K,) + v.shape[2:]), np.zeros((n,) + v.shape[2:])][act]
-    rhs = [cell(sk.m), cell(sk.H)] + ([cell(sk.m * sk.S)] if sk.S is not None else [])
-    x = np.empty((N, len(rhs)))
-    x[order] = lu.solve(np.column_stack(rhs)[order])
+    b0 = np.column_stack([cell(sk.m), cell(sk.H)] + ([cell(sk.m * sk.S)] if sk.S is not None else []))
+    x = np.empty_like(b0)
+    x[order] = _solve(sp.csc_matrix((np.r_[diag, -e_val], (np.r_[at, to], np.r_[at, frm])), shape=(N, N)), b0[order],
+                      upwind, big)
+    moved = e_val[:, None] * x[lc]
+    x = b0 + np.stack([np.bincount(lr, moved[:, v], N) - np.bincount(lc, moved[:, v], N) for v in range(x.shape[1])],
+                      axis=1)
     # the momentum, each arrival projected onto its receiver's plane: two components per cell in its own plane's frame
     owner = np.where(act < n * K, act // K, act - n * K)
     e1, e2 = tangent_frames(columns.normal)
     P = cell(sk.m[..., None] * sk.U)
-    p = np.column_stack([np.einsum("ij,ij->i", P, e1[owner]), np.einsum("ij,ij->i", P, e2[owner])])
-    held = P - p[:, :1] * e1[owner] - p[:, 1:] * e2[owner]   # any normal part a cell held: the wall's
+    p0 = np.column_stack([np.einsum("ij,ij->i", P, e1[owner]), np.einsum("ij,ij->i", P, e2[owner])])
+    held = P - p0[:, :1] * e1[owner] - p0[:, 1:] * e2[owner]  # any normal part a cell held: the wall's
     R = np.stack([np.stack([np.einsum("ij,ij->i", e[b], g[a]) for g in (e1, e2)], axis=1) for e in (e1, e2)], axis=1)
     Re = R[e_edge]                                            # [arrival, receiver's component, donor's component]
     rows2 = [2 * at, 2 * at + 1] + [2 * to + q for q in (0, 0, 1, 1)]
     cols2 = [2 * at, 2 * at + 1] + [2 * frm + q for q in (0, 1, 0, 1)]
     vals2 = [diag, diag] + [-e_val * Re[:, q, s] for q, s in ((0, 0), (0, 1), (1, 0), (1, 1))]
-    lu2 = _factor(sp.csc_matrix((np.concatenate(vals2), (np.concatenate(rows2), np.concatenate(cols2))),
-                                shape=(2 * N, 2 * N)), upwind)
     rhs2 = np.empty(2 * N)
-    rhs2[2 * at], rhs2[2 * at + 1] = p[:, 0], p[:, 1]
-    y = lu2.solve(rhs2)
-    p = np.column_stack([y[2 * at], y[2 * at + 1]])
+    rhs2[2 * at], rhs2[2 * at + 1] = p0[:, 0], p0[:, 1]
+    y = _solve(sp.csc_matrix((np.concatenate(vals2), (np.concatenate(rows2), np.concatenate(cols2))),
+                             shape=(2 * N, 2 * N)), rhs2, upwind, np.repeat(big, 2))
+    sent_p = e_val[:, None] * np.column_stack([y[2 * at], y[2 * at + 1]])[lc]
+    got = np.einsum("eqs,es->eq", Re, sent_p)
+    p = p0 + np.stack([np.bincount(lr, got[:, v], N) - np.bincount(lc, sent_p[:, v], N) for v in range(2)], axis=1)
     P = p[:, :1] * e1[owner] + p[:, 1:] * e2[owner]
     # the books, face by face: what every arrival took from its sending cell at the sub-step's end
-    moved = e_val[:, None] * np.column_stack([x[:, :2], p])[local[e_col]]
-    sent = np.stack([np.bincount(e_edge, moved[:, v], D) for v in range(4)], axis=1)
+    sent = np.stack([np.bincount(e_edge, v, D) for v in (moved[:, 0], moved[:, 1], sent_p[:, 0], sent_p[:, 1])], axis=1)
     carried = sent[:, 2:3] * e1[a] + sent[:, 3:4] * e2[a]
     lost = np.einsum("ij,ij->i", carried, columns.normal[b])[:, None] * columns.normal[b]
     net = lambda v: np.bincount(b, v, n) - np.bincount(a, v, n)
@@ -1405,17 +1605,41 @@
 
 
 # -- the flow in the sub-steps (sub-plan 21 Task 6c) --------------------------------------------------------------------
-FLOW_CYCLES = 10                  # flow sub-steps per heat sub-step: 1 ms under the 10 ms heat sub-step (Task 6c's
-                                  # ruling, 2026-10-22). On the 100 mm sphere's front melted 0.3 mm deep under the 93 s
-                                  # forcing, 50 ms of flow at 1 ms match 0.25 ms to 3.9 % in where the metal went and
-                                  # 4 % in kinetic energy, the Courant numbers at most 4.8 and the transport destroying
-                                  # 27 % of the kinetic energy made, as at 0.25 ms (the grid's own upwind smearing); at
-                                  # 2-10 ms 44-95 % and 0.3-4.8 times, because a runoff front crossing rigid metal
-                                  # advances at most one column a flow sub-step (metal landing on a rigid top joins the
-                                  # flow at the next), 2 m/s on 2 mm columns at 1 ms. The Rayleigh-Taylor growth dt s
-                                  # stays under 0.007 at 1 ms on 3.4 mm of liquid at 48 m/s2 (0.67 at 10 ms).
+FLOW_CYCLES = 10                  # flow sub-steps in a heat sub-step where the skins remember no flow (6c's fixed count,
+                                  # 1 ms under the 10 ms heat sub-step, chosen on the 100 mm sphere's runoff when a front
+                                  # crossing rigid metal advanced one column a flow sub-step; sub-plan 21 Task 6c)
+FLOW_CYCLES_MIN, FLOW_CYCLES_MAX = 1, 40   # the adaptive count's floor and cap (Task 6d): 10 ms and 0.25 ms flow
+                                  # sub-steps under the 10 ms heat sub-step
+FLOW_COURANT = 6.0                # the adaptive count's target: the largest Courant number of a flow sub-step (Task 6d).
+                                  # On the 100 mm sphere's runoff at 93 s, with the front limit lifted, 50 ms of flow at
+                                  # 1.25 ms (Courant 6.1) put the metal within 2.5 % of the 0.25 ms run's places and
+                                  # its kinetic energy within 5 %, the thickness step's linearisations at 4 at most; at
+                                  # 1.5 ms (7.3) 3.6 % and 12 %, the linearisations at their cap of 8; at 1 ms 2.2 %
+FLOW_GROWTH = 0.02                # and the largest Rayleigh-Taylor growth dt s of one: backward Euler's growth rate is
+                                  # then within about dt s / 2 = 1 % (6b's cosine met it within 0.42 % at 0.016)
+FLOW_RELAX = 2.0                  # a flow sub-step at most this many times the last heat sub-step's
 FLOW_SECONDS = {"push": 0.0, "momentum": 0.0, "thickness": 0.0, "transport": 0.0, "regrid": 0.0}   # wall time, for
                                   # measuring the cost (as KRYLOV counts the coupled solves)
+
+
+def flow_cycles(memory, dt):
+    """The number of flow sub-steps in a heat sub-step dt [s] (sub-plan 21 Task 6d), from the skins' `flow_memory` --
+    the last heat sub-step's largest Courant number C and Rayleigh-Taylor growth g over its flow sub-steps, and their
+    length h_0 -- by the scalings measured on the 100 mm sphere: the Courant number in proportion to the flow sub-step,
+    the growth to its square. The flow sub-step is the longest of FLOW_RELAX h_0, h_0 FLOW_COURANT / C and
+    h_0 (FLOW_GROWTH / g)^(1/2), so that the next one meets both targets and grows by at most FLOW_RELAX a heat
+    sub-step, and the count dt over it, rounded up, within FLOW_CYCLES_MIN and FLOW_CYCLES_MAX; FLOW_CYCLES where the
+    skins remember no flow (at the start, or when nothing flowed in the last heat sub-step). Deterministic: the same
+    memory and dt give the same count."""
+    courant, growth, length = (float(x) for x in memory)
+    if length <= 0.0:
+        return FLOW_CYCLES
+    h = FLOW_RELAX * length
+    if courant > 0.0:
+        h = min(h, length * FLOW_COURANT / courant)
+    if growth > 0.0:
+        h = min(h, length * math.sqrt(FLOW_GROWTH / growth))
+    return int(min(max(math.ceil(dt / h - 1e-9), FLOW_CYCLES_MIN), FLOW_CYCLES_MAX))
 
 
 @dataclass
@@ -1425,7 +1649,8 @@
     each column's plane (`along_meridian` of the surface flow's tau), `g` [m/s2, (3,)] the acceleration the metal feels
     in the body's frame, minus the body's own, `v_hat` the flight direction, `rho_l` [kg/m3] and `sigma` [N/m] the
     liquid's density and surface tension, `A_top` [m2, (n,)] the top areas the shear acts over, and `cycles`, the flow
-    sub-steps in each heat sub-step."""
+    sub-steps in each heat sub-step: None (the default since Task 6d) chooses them each heat sub-step (`flow_cycles`),
+    a number fixes them. `plan` is what a trial pass chose, for its real pass (`advance`)."""
     columns: object
     dp_ds: np.ndarray
     shear: np.ndarray
@@ -1434,7 +1659,8 @@
     rho_l: float
     sigma: float
     A_top: np.ndarray
-    cycles: int = FLOW_CYCLES
+    cycles: int = None
+    plan: list = None
 
 
 @dataclass
@@ -1450,6 +1676,7 @@
     growth: float = 0.0      # the thickness step's largest Rayleigh-Taylor growth dt s
     passes: int = 0          # the thickness step's most linearisations
     sub_steps: int = 0       # the flow sub-steps taken (none where nothing flows)
+    cycles: int = 0          # the most flow sub-steps in one heat sub-step (Task 6d)
 
     @classmethod
     def zero(cls, n):
@@ -1461,29 +1688,38 @@
             setattr(self, key, getattr(self, key) + getattr(other, key))
         self.courant, self.growth = max(self.courant, other.courant), max(self.growth, other.growth)
         self.passes, self.sub_steps = max(self.passes, other.passes), self.sub_steps + other.sub_steps
-
-
-def flow_substeps(sk, flow, dt):
-    """`flow.cycles` flow sub-steps of dt / cycles [s] each, every one in the order of 6b's ruling: the push, surface
-    tension only where a column's top flows (`push`, `flowing_tops`); the momentum through the depth with its response
-    (`momentum_step`); the thickness terms implicit in the heights at the end of the transport (`thickness_step`); the
-    transport with the faces' own push (`transport`, regrid=False). The cells are split and merged between the flow
-    sub-steps, and the heat sub-step's own regrid follows the last (spec 6.2: once after each sub-step's exchanges).
-    Updates the cells; returns a FlowBook. Where nothing flows nothing moves and no sub-step is taken: any velocity left
-    in the rigid cells is booked to the wall and cleared, as `momentum_step` books an anchored cell's.
-
-    The flow's sub-step is chosen from measurements on the 100 mm sphere's melting state (sub-plan 21 Task 6c, Asha
-    2026-10-21, option c; FLOW_CYCLES)."""
+        self.cycles = max(self.cycles, other.cycles)
+
+
+def flow_substeps(sk, flow, dt, cycles=None):
+    """`cycles` flow sub-steps of dt / cycles [s] each -- by default `flow.cycles`, or, where that is None, as many as
+    `flow_cycles` chooses from the skins' `flow_memory`, which they then update (sub-plan 21 Task 6d) -- every one in the
+    order of 6b's ruling: the push, surface tension only where a column's top flows (`push`, `flowing_tops`); the
+    momentum through the depth with its response (`momentum_step`); the thickness terms implicit in the heights at the
+    end of the transport (`thickness_step`); the transport with the faces' own push (`transport`, regrid=False). The
+    cells are split and merged between the flow sub-steps, and the heat sub-step's own regrid follows the last (spec
+    6.2: once after each sub-step's exchanges). Updates the cells; returns a FlowBook. Where nothing flows nothing moves
+    and no sub-step is taken: any velocity left in the rigid cells is booked to the wall and cleared, as `momentum_step`
+    books an anchored cell's, and the skins forget the flow.
+
+    The flow's sub-step was fixed in 6c at 1 ms (FLOW_CYCLES) and is chosen since 6d from the last heat sub-step's
+    Courant number and growth, against targets set from the sub-step study on the 100 mm sphere (FLOW_COURANT,
+    FLOW_GROWTH)."""
     n = sk.n
     book = FlowBook.zero(n)
+    adaptive = flow.cycles is None
     if not (sk.valid & np.isfinite(cell_viscosity(sk, 0.0))).any():
         book.wall -= sk.column_momentum()
         sk.U = np.zeros_like(sk.U)
+        if adaptive:
+            sk.flow_memory = np.zeros(3)
         return book
-    h = dt / flow.cycles
+    if cycles is None:
+        cycles = flow_cycles(sk.flow_memory, dt) if adaptive else flow.cycles
+    h = dt / cycles
     col = flow.columns
     clock = lambda key, started: FLOW_SECONDS.__setitem__(key, FLOW_SECONDS[key] + time.perf_counter() - started)
-    for k in range(flow.cycles):
+    for k in range(cycles):
         if k:
             started = time.perf_counter()
             sk.regrid()
@@ -1502,6 +1738,9 @@
         clock("transport", started)
         book.add(FlowBook(mb.push + tb.push, mb.shear, mb.wall + tb.wall + tr.wall, tr.momentum, tr.courant, tb.growth,
                           tb.passes, 1))
+    book.cycles = cycles
+    if adaptive:
+        sk.flow_memory = np.array([book.courant, book.growth, h])
     return book
 
 
```

- [ ] **Step 6: The chosen count in the skin body** in `reentry_model/skinbody.py`:

```diff
--- a/reentry_model/skinbody.py
+++ b/reentry_model/skinbody.py
@@ -44,7 +44,9 @@
 as the heating is, from the step's aero state and the tops' temperatures -- its wall-pressure gradient along the
 meridian and its gas shear -- with the deceleration, which the metal feels along the flight direction; each part takes
 them interpolated in time at its middle, as it takes the heating, and every sub-step of both passes runs the flow after
-its conduction (`skin.flow_substeps`), `flow_cycles` flow sub-steps to each. With the flow on the liquid stays and flows
+its conduction (`skin.flow_substeps`), as many flow sub-steps to each as `skin.flow_cycles` chooses from the last
+heat sub-step's Courant number and growth (sub-plan 21 Task 6d; `flow_cycles` fixes them). With the flow on the liquid
+stays and flows
 (`remove` False in a run; spraying is stage 4's). The momentum books: the skins' momentum changes by the flow's impulses
 less what leaves with removed or handed-down metal, and `momentum_residual` reports what they miss."""
 import dataclasses
@@ -129,8 +131,9 @@
 class SkinBody(ThermalBody):
     """A ThermalBody whose surface is resolved by corner columns with layered skins on a moving mesh (module docstring).
     `remove`: liquid leaves as it forms (stage 2); `demise_fraction`: the run ends below this share of the initial mass;
-    `flow`: the metal flows (sub-plan 21), `flow_cycles` flow sub-steps in each heat sub-step (None: skin.FLOW_CYCLES),
-    driven by `surface_flow` (a `surface_flow.SurfaceFlow`, its defaults if None) from each macro step's aero state."""
+    `flow`: the metal flows (sub-plan 21), `flow_cycles` flow sub-steps in each heat sub-step (None: as many as
+    `skin.flow_cycles` chooses each heat sub-step, Task 6d), driven by `surface_flow` (a `surface_flow.SurfaceFlow`, its
+    defaults if None) from each macro step's aero state."""
 
     def __init__(self, mesh, material, solver, mass_kg, settings=None, T0=300.0, emissivity=None, T_ambient=0.0,
                  v_hat=(1.0, 0.0, 0.0), deform_depth=movemesh.DEFORM_DEPTH, remove=True, demise_fraction=0.01,
@@ -138,9 +141,9 @@
                  surface_flow=None):
         if repeat_rule not in REPEAT_RULES:
             raise ValueError("repeat_rule must be one of {}, got {!r}".format(REPEAT_RULES, repeat_rule))
-        self.flow_cycles = skin_mod.FLOW_CYCLES if flow_cycles is None else flow_cycles
-        if not (isinstance(self.flow_cycles, (int, np.integer)) and self.flow_cycles >= 1):
-            raise ValueError("flow_cycles must be an integer of at least 1, got {!r}".format(flow_cycles))
+        self.flow_cycles = flow_cycles
+        if not (flow_cycles is None or isinstance(flow_cycles, (int, np.integer)) and flow_cycles >= 1):
+            raise ValueError("flow_cycles must be None or an integer of at least 1, got {!r}".format(flow_cycles))
         self.flow_on, self.surface_flow = bool(flow), surface_flow
         self.repeat_rule = repeat_rule
         self.part_max = part_max                             # s, or None: no cap on the parts while melting
@@ -178,7 +181,7 @@
         self.impulse, self.momentum_lost = np.zeros(3), np.zeros(3)   # the momentum books since time zero [N s, kg m/s]
         self.impulse_size = 0.0                              # N s: the push's, shear's and walls' impulses, in size
         self.flow_courant = self.flow_growth = 0.0           # this macro step's worst flow sub-step
-        self.flow_passes = 0
+        self.flow_passes = self.flow_cycles_max = 0          # and the most flow sub-steps in a heat sub-step
         self._volume = None                                  # m3, the 3D model's volume since the last move (mass)
         self.iterations, self.last, self.last_move = [], None, None
         self._fill(T0)
@@ -227,7 +230,7 @@
             start = last if last is not None and last.dp_ds.shape == end.dp_ds.shape else end
             self._forcing, self._forcing_prev = (start, end), end
             self.flow_courant = self.flow_growth = 0.0
-            self.flow_passes = 0
+            self.flow_passes = self.flow_cycles_max = 0
         sk, mat = self.skins, self.material
         top = np.where(sk.valid[:, 0], sk.T[:, 0], np.nan)
         if self.part_max is not None and mat.melts:
@@ -364,7 +367,7 @@
             self.impulse += (fb.push + fb.shear + fb.wall + fb.carried).sum(axis=0)
             self.impulse_size += float(sum(np.linalg.norm(x, axis=1).sum() for x in (fb.push, fb.shear, fb.wall)))
             self.flow_courant, self.flow_growth = max(self.flow_courant, fb.courant), max(self.flow_growth, fb.growth)
-            self.flow_passes = max(self.flow_passes, fb.passes)
+            self.flow_passes, self.flow_cycles_max = max(self.flow_passes, fb.passes), max(self.flow_cycles_max, fb.cycles)
         self.iterations.append(res.iterations)
         self.last = dataclasses.replace(res, Q_conv=E_top / dt, Q_rad=E_rad / dt)    # the gas-side powers, for the history
         self._move_base(book, m_min, g0)
@@ -532,12 +535,14 @@
                 "wall_rebuild_s": self.wall["rebuild"], "step_parts": self.step_parts,
                 # the flow (sub-plan 21 Task 6c): its momentum books, the fastest cell, the top's roughness -- the largest
                 # height difference between neighbours over their spacing -- and the flow sub-steps' worst Courant
-                # number, Rayleigh-Taylor growth and linearisations in this macro step
+                # number, Rayleigh-Taylor growth and linearisations in this macro step, and the most of them a heat
+                # sub-step took (Task 6d)
                 "momentum_residual": self.momentum_residual(),
                 "speed_max_ms": float(np.linalg.norm(sk.U, axis=2).max()) if sk.U.size else 0.0,
                 "roughness_max": float((np.abs(th[i] - th[j]) / self.columns.edge_length).max()),
                 "flow_courant_max": self.flow_courant, "flow_growth_max": self.flow_growth,
-                "flow_passes_max": self.flow_passes, "wall_flow_s": self.wall["flow"]}
+                "flow_passes_max": self.flow_passes, "flow_cycles_max": self.flow_cycles_max,
+                "wall_flow_s": self.wall["flow"]}
 
     def _recession(self):
         """How far each corner's top has moved inward along its starting direction since time zero [m]."""
```

- [ ] **Step 7: The history's count** in `reentry_model/coupled.py`:

```diff
--- a/reentry_model/coupled.py
+++ b/reentry_model/coupled.py
@@ -47,7 +47,7 @@
                     "wall_skins_s", "wall_3d_s", "wall_move_s", "wall_rebuild_s", "step_parts",
                     # the flow (sub-plan 21 Task 6c, spec 10.2)
                     "momentum_residual", "speed_max_ms", "roughness_max", "flow_courant_max", "flow_growth_max",
-                    "flow_passes_max", "wall_flow_s"]
+                    "flow_passes_max", "flow_cycles_max", "wall_flow_s"]
 PVD_TEMPLATE = '<?xml version="1.0"?>\n<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">\n<Collection>\n{}</Collection>\n</VTKFile>\n'
 
 
```

- [ ] **Step 8: The count's setting** in `reentry_model/cli.py`:

```diff
--- a/reentry_model/cli.py
+++ b/reentry_model/cli.py
@@ -129,8 +129,8 @@
     option B', 2026-10-16: 0.8 and 0.4 mm), `_cell<um>um` for an explicit cell other than the
     preset's, `_growth<ratio>`, `_sub<ms>ms` and `_deform<mm>mm` for another growth, sub-step or deforming depth,
     `_noswell` for skins held at the solid's density, `_fine<mm>x<mm>mm` for another fine band (`_fine0` for none), and,
-    with the flow on (sub-plan 21 Task 6c), `_skinflow-on` followed by `_flowcycles<n>` for other than skin.FLOW_CYCLES
-    flow sub-steps in each heat sub-step; empty for the defaults."""
+    with the flow on (sub-plan 21 Task 6c), `_skinflow-on` followed by `_flowcycles<n>` where the flow sub-steps in each
+    heat sub-step are fixed at n rather than chosen each heat sub-step (Task 6d); empty for the defaults."""
     out = "_skin-dev" if preset == "dev" else ""
     if min_thickness_mm != 0.8:
         out += "_skinmin{:g}mm".format(min_thickness_mm)
@@ -152,7 +152,7 @@
         out += "_fine{:g}x{:g}mm".format(fine_band_mm, fine_size_mm)
     if flow:
         out += "_skinflow-on"
-        if flow_cycles is not None and flow_cycles != skin.FLOW_CYCLES:
+        if flow_cycles is not None:
             out += "_flowcycles{}".format(flow_cycles)
     return out
 
@@ -313,7 +313,9 @@
                          "girin, whose continuum step (--dt-continuum) follows; off: the heat-only skins, liquid removed "
                          "as it forms with --removal instant (default)")
     me.add_argument("--skin-flow-cycles", type=int, default=None,
-                    help="with --skin-flow on: flow sub-steps in each skin sub-step (default {})".format(skin.FLOW_CYCLES))
+                    help="with --skin-flow on: flow sub-steps in each skin sub-step, fixed (default: chosen each skin "
+                         "sub-step from the flow's Courant number and growth, {} to {})".format(skin.FLOW_CYCLES_MIN,
+                                                                                            skin.FLOW_CYCLES_MAX))
     me.add_argument("--skin-preset", choices=tuple(skin.SKIN_PRESETS), default="production",
                     help="skin cell-size preset: production 10 um, dev 20 um (default %(default)s)")
     me.add_argument("--deform-depth", type=float, default=15.0,
@@ -415,7 +417,8 @@
                         "growth": args.skin_growth, "substep_ms": args.skin_substep, "preset": args.skin_preset,
                         "deform_depth_mm": args.deform_depth, "n_columns": the_body.columns.n,
                         "fine_band_mm": args.fine_band, "fine_size_mm": min(args.fine_size, h_surface * 1e3),
-                        "swell": args.skin_swell == "on", "flow": args.skin_flow, "flow_cycles": the_body.flow_cycles}
+                        "swell": args.skin_swell == "on", "flow": args.skin_flow,
+                        "flow_cycles": "adaptive" if the_body.flow_cycles is None else the_body.flow_cycles}
     if args.heating == "physics":
         info.update({"stagnation": args.stagnation, "bridging_heat": args.bridging_heat, "matting_n": args.matting_n,
                      "accommodation": args.accommodation, "catalycity": args.catalycity})
```

- [ ] **Step 9: Run them to see them pass.** Expected: `tests/test_reentry_model_skin.py` 80 passed,
  `tests/test_reentry_model_skinbody.py` 17 passed (2 deselected, the reference tier's),
  `tests/test_reentry_model_cli.py` 66 passed, `tests/test_reentry_model_coupled.py` 9 and
  `tests/test_reentry_model_columns.py` 13 as before (185 in 1 269 s); the unit tier 490 passed, 1 skipped (1 411 s).
  The flow is off by default, so no flight's result changes; with it on, a front over rigid metal now flows on within
  the flow sub-step and the count is chosen.

- [ ] **Step 10: Commit.**

```bash
git add reentry_model/skin.py reentry_model/skinbody.py reentry_model/coupled.py reentry_model/cli.py tests/test_reentry_model_skin.py tests/test_reentry_model_skinbody.py tests/test_reentry_model_cli.py
git commit -m "Skins: the flow made affordable -- the transport's large cycles iterated from the upwind solve with the cells installed from its fluxes, so the books are exact by construction; metal landing on a rigid top flows on in the same solve at the velocity it arrived with, lifting the one-column-a-flow-sub-step front limit; the flow sub-steps chosen each heat sub-step from the last one's Courant number and growth, the real pass taking its trial pass's counts; the momentum step on the moving columns only (sub-plan 21 Task 6d)"
```

---

### Task 7: The droplet momentum sink at the top cells (a placeholder for stage 4's spraying)

**Files:**
- Modify: `reentry_model/skin.py` (`Skins.take_top_mass`; `PI_FLAG`, `ReleaseBook`, `release_top`,
  `prescribed_release`, `droplet_momentum_ratio`; `Flow.release` and `Flow.slurry`; `FlowBook.released`,
  `.released_H`, `.released_p` and `.released_solid`, with `.zero` and `.add`; `flow_substeps` releasing in every flow
  sub-step; `FLOW_SECONDS["release"]`), `reentry_model/skinbody.py` (`SkinBody(..., release=None,
  slurry_spray=False)`, `.released_mass`, `.released_enthalpy`, `.droplet_momentum`, `.droplet_shear`, the books, the
  history's three columns), `reentry_model/coupled.py` (`SKINBODY_COLUMNS`, `skin_results`)
- Test: `tests/test_reentry_model_skin.py` (seven tests), `tests/test_reentry_model_skinbody.py` (one test)

**Interfaces:**
- Consumes: Task 2's `Skins.U`, `Skins.take_top` and `momentum_out`; Task 3's `cell_viscosity`; 6c's `Flow`,
  `FlowBook`, `flow_substeps`, `advance(..., flow)` and the skin body's momentum books (`momentum_lost`,
  `momentum_residual`); 6d's chosen counts, which the real pass takes from its trial pass.
- Produces: `Skins.take_top_mass(dm, admit) -> (mass, enthalpy, solid)` [kg, J, kg per column], `momentum_out` set;
  `skin.release_top(sk, dm, slurry=False) -> ReleaseBook(mass, enthalpy, momentum, solid)` ((n,) [kg], (n,) [J],
  (n, 3) [kg m/s], (n,) [kg]); `skin.prescribed_release(rate)`, a release at a given rate [kg/(m2 s)];
  `skin.droplet_momentum_ratio(released, shear) -> (Pi, each (n,), flagged)` and `skin.PI_FLAG` (0.2);
  `Flow(..., release=None, slurry=False)`; `FlowBook.released` (n,) [kg], `.released_H` (n,) [J], `.released_p`
  (n, 3) [kg m/s] and `.released_solid` (n,) [kg]; `SkinBody(..., release=None, slurry_spray=False)` with
  `.released_mass` [kg] and `.released_enthalpy` [J] since time zero and `.droplet_momentum` (n, 3) [kg m/s] and
  `.droplet_shear` (n, 3) [N s] over the last macro step; history columns `released_mass_kg`,
  `droplet_momentum_ratio` and `droplet_momentum_flagged`; results `released_mass_kg`, `droplet_momentum_ratio_max`
  and `droplet_momentum_flagged_steps`.
- **What sub-plan 22 supplies and calls:**
  1. *A release*, the one thing it adds: a callable `release(sk) -> (n,)` [kg/(m2 s)], Girin's stripping rate per unit
     of each column's top area (`SprayResult.mdot`, per unit wetted area), from the skins as they are when it is asked
     -- b the flowing depth `sk.depth_above(T_rigid)` and the liquid depth `sk.depth_above(T_liquidus)` (spec 6.7) --
     and the part's gas state, which it holds; given as `SkinBody(..., flow=True, release=..., slurry_spray=...)`, which
     hands it to every part's `Flow`. It is asked in every flow sub-step of the trial pass, of every damped repeat and
     of parts run again, so it keeps no records of its own.
  2. *Its cap and its cascade*: `release_top` takes rate A_top h from the top cells and never more than they may give
     (fact 29's inventory cap). The cascade rule's own rates, cell by cell (the 2026-10-07 design 8.1), are sub-plan
     22's: `release_top(sk, dm, slurry)` takes whatever dm that rule computes.
  3. *Its records*: the pass's FlowBook carries what left each column, summed over the flow sub-steps -- spec 8's
     `h_kg` is `released_H / released` and `solid_fraction` `released_solid / released` -- and the skin body adds the
     accepted passes' only. Records per mode need the mode, which the book does not carry yet: sub-plan 22 extends it.

Stage 4's spraying (spec 8; sub-plan 22) takes Girin's release from each column's top cells. This task builds what
the release does to the cells and the flow, so that sub-plan 22 only has to supply the rate: the release from the top
cells, the momentum it carries away (fact 29's sink), and the droplet-momentum ratio Pi that spec 10.2 puts in the
history.

**The release.** `release_top` takes the mass asked from each column's top down through the cells it may take: the run
of fully liquid cells from the top, at or above the liquidus, and, with `slurry` (spec 8's `--slurry-spray on`), every
flowing cell down to the first rigid one. At most what they hold leaves. A slurry top without the flag, a floating crust
and anchored metal release nothing: nothing passes through rigid metal (spec 6.5). Each cell's metal leaves at its own
specific enthalpy and velocity; a partly taken cell keeps both, its temperature and its sensitivity, and the cells
beneath are untouched, so the thinned top cell merges at the next regrid, as spec 6.2 has spray do. A column whose
cells are all taken is emptied: the flow sub-steps after it land the runoff on it as on any rigid top, and the heat
sub-step's draw restores its minimum.

**The release speed** (spec 6.4: "The droplets carry momentum away from the top cells at their release speed (fact
29's sink)"). Fact 29's sink is mdot_strip V_s, V_s the sheared layer's surface speed: its Pi = v_melt delta_m / nu is
mdot V_s / tau_gas with V_s = tau delta_m / mu, the top of a Couette layer delta_m deep. The columns resolve that
profile, so V_s is the top cells' own velocity, and the droplets, made of the top cells' metal, leave with it. Removing
the cells with their momentum is then the sink, with nothing to add to the books. A droplet speed other than the cells'
-- Girin's closure V_s, say -- would leave the difference times the released mass to be booked on the metal left
behind or supplied by the gas, which neither the spec nor fact 29 asks for (ruling 1).

**Where it runs: in every flow sub-step, after the transport.** Spec 7.1's sub-step is heat, flow, sideways exchange
and spray, and since 6c the last three run in every flow sub-step: the rate is asked of the skins as the transport left
them, and rate A_top h is released. Not once a heat sub-step: the top cells cross a column for each unit of their
Courant number, which 6d's count keeps near 6 a flow sub-step, so they cross some 25 columns of 2 mm in a 10 ms heat
sub-step at the 5 m/s of the 100 mm sphere's driven state, and up to 240 at the count's cap of 40. A release made once
would spray metal from where the runoff had taken it rather than where its rate was found, and fact 29's operator split
between the runoff and the stripping would be first order in the heat sub-step rather than in the flow's. The
trial pass runs the real pass's flow sub-steps (6c) at the real pass's counts (6d), so both release at the same
moments, and the base law includes the release: on a heated strip it moves the law by 0.20 % of the throughput. With the
flow off nothing asks for a release: the heat-only skins' instant removal stays as it is, and a skin body refuses a
release without the flow.

**Pi** (spec 10.2; fact 29). For each column over a macro step, Pi_i = |p_i| / |J_i|: p_i the momentum the release
took from its top cells (`FlowBook.released_p`) and J_i the gas shear's impulse on its top (`FlowBook.shear`), both
integrated over the top and the step, so Pi_i is mdot V_s / tau_gas. The history's Pi is the sum of |p_i| over the sum
of |J_i| over the sheared columns: the share of the gas shear's momentum that left with droplets, each column's ratio
weighted by its shear. `droplet_momentum_flagged` counts the sheared columns whose own Pi_i exceeds fact 29's 0.2. A top
the gas does not shear has no Pi: on the 100 mm sphere's driven state 2 403 liquid tops run over the rigid rear at
2.5-4.7 m/s (10th to 90th percentile) under no gas shear, and a release prescribed there made Pi_i infinite on all 4 712
unsheared columns; Girin's modes release only on the windward face (`SprayModel.evaluate`'s `windward`). By hand, on
liquid aluminium 0.234 mm deep at its steady Couette profile under 20 Pa, released at the top: Pi_i =
mdot (H - dz_0 / 2) / mu, fact 29's v_melt delta_m / nu, 0.42 at 1 mm/s and 0.042 at 0.1 mm/s.

**Measured on the 100 mm sphere** (2026-10-24, `prototype/work-2026-10-18-column-flow/measure/task7/sphere_pi.py`;
processor times under other sessions' load). 6d's driven state -- the front's top 0.3 mm liquid, driven from rest for
20 ms by the forcing of the flight at 93 s; 160 000 cells flowing, 7 108 liquid tops -- released at rho_l v_melt over the
windward face for ten more flow sub-steps of 1 ms:

| v_melt | Released | Pi | Flagged, of 4 705 sheared columns | Pi_i median, 90th percentile, largest |
|---|---|---|---|---|
| 1 mm/s | 0.376 g | 0.145 | 1 146 | 0.136, 0.235, 0.263 |
| 0.1 mm/s | 0.038 g | 0.015 | 0 | 0.014, 0.024, 0.027 |

The table shows Pi at a third of fact 29's estimates, 0.48 and 0.048, which took the film at its steady Couette speed
tau delta_m / mu. The 0.3 mm layer, its viscous time h^2 / nu 0.17 s, is far from steady after 30 ms -- its top at a
median 4.3 m/s under 64 Pa, against 15 m/s for the steady profile -- so these numbers measure the developing flow, and a
flight's steadier layers will carry more. At 1 mm/s a quarter of the sheared columns already pass fact 29's flag. The
release costs 0.054 s of processor time a call on the 9 417 columns, against about 2 s for a 1 ms flow sub-step of this
state (20.0 s for ten without the release, 21.0 s with it).

**No release, no change.** With no release (the default) every flow sub-step is 6d's, and so is every result:
against 6d's code, the flow-off removal run of `test_heating_with_removal_keeps_mass_and_energy_exact` (numpy's
generator seeded), a skin body with the flow on over two 0.25 s steps under 3 MW/m2 and two heat sub-steps' flow on the
mixed sphere skins give the same 109 arrays -- the 3D field, the mesh, the skins' state, every history column but the
wall times and the three new ones, and the flow's books -- to the last bit (`measure/task7/bitid.py`). A release of
nothing, a zero rate or a rate only where a top is not liquid, does the same.

**Rulings** (2026-10-24, within spec 6.4, 6.5, 7.1, 8 and 10.2 and fact 29; each is what, then why, then the cost if
wrong):
1. *The droplets leave at their cells' own velocity.* Why: fact 29's V_s is the sheared layer's surface speed, which
   the columns resolve, and the spec has the droplets' momentum leave the top cells (6.4, 8); removing the cells with
   their momentum is then the sink, and the books close with nothing new. Cost if wrong (the droplets leaving at
   Girin's closure speed instead): the difference times the released mass is an impulse on the metal left behind or
   the gas's, and Pi is off by the ratio of the two speeds; sub-plan 22 computes both.
2. *The release in every flow sub-step, after the transport.* Why: above. Cost if wrong: the rate is asked up to 40
   times a heat sub-step instead of once, 0.05 s each against a 2 s flow sub-step.
3. *Liquid is a cell at or above the liquidus.* Why: it is fully liquid by the material's own law, the history's liquid
   depth (spec 6.7) and the viscosity's liquid branch; not stage 2's T_feed, 2 K higher, the element model's film ramp.
   Cost if wrong: cells within those 2 K wait for the slurry flag.
4. *The run from the top, cell by cell as the cells lie: liquid without the flag, every flowing cell with it, stopping
   at the first rigid cell.* Why: spec 8's "liquid first" and "exposed slurry", the 2026-10-07 design's cascade, which
   works down to the first cell at or below T_rigid, and spec 6.5's rigid metal. Cost if wrong: with the flag, a slurry
   cell over liquid goes before the liquid it covers, which sub-plan 22's cascade may order otherwise through its dm.
5. *The rate per unit of the top area, capped by what the cells may give.* Why: Girin's mdot is per unit wetted area
   (`SprayModel.evaluate` releases min(mdot A dt, m_f)), the heating's q is per top area, and fact 29's checklist caps
   each removal by the mass present. Cost if wrong: the release differs by the top area's excess over the base area,
   about 2H/R (spec 14 assumption 10).
6. *Pi over the sheared tops, weighted by their shear; the flag per sheared column.* Why: above. Cost if wrong: a
   column's large Pi_i is hidden in the mean, but counted by the flag, and a release from an unsheared top is in
   neither.
7. *The release a callable of the skins that keeps no records; a release needs the flow.* Why: the trial pass, the
   damped repeats and parts run again all ask for it, and with the flow off the heat-only skins' removal is stage 2's
   stand-in for flow and spray (spec 11.2). Cost if wrong: nothing for this task.
8. *The droplets' momentum is booked as lost with the metal, as the removed liquid's is.* Why: it leaves with the
   metal; `momentum_residual`'s scale stays the impulses' size. Cost if wrong: only that scale.

**Measured against the tests** (2026-10-24): what the release takes from six columns of every arrangement, to 1.9e-16
(the books) and 2.2e-16 (the cells left); on the mixed sphere skins each column short of exactly its book, to 1.6e-16,
2.3e-16 and 2.0e-16 of the largest column's mass, enthalpy and momentum, the rigid metal untouched to 2.2e-16, and 6.4 %
of the released mass solid with the flag; the flow sub-steps with a release equal to the four steps and `release_top`
by hand to the last bit, the momentum books to 2.1e-12 and 1.1e-12 of the largest impulse at one and three cycles;
columns emptied in the first of four flow sub-steps, the pass's books to 4.1e-16; no release, or nothing to release, to
the last bit; Pi by hand to 2.2e-15; the trial and real passes releasing the same to the last bit, the law met to
4.0e-17 of the throughput; the skin body over two 0.25 s steps releasing 0.53 g, its mass to the last bit, its energy
to 9.8e-13 and its momentum to 2.3e-17, Pi 0.0048.

- [ ] **Step 1: The release's tests** — append to `tests/test_reentry_model_skin.py`:

```diff
--- a/tests/test_reentry_model_skin.py
+++ b/tests/test_reentry_model_skin.py
@@ -1738,3 +1738,253 @@
     solved = (sk.valid & ~skin.anchored_cells(sk, flowing)).any(axis=1)
     assert 0 < solved.sum() < col.n and np.array_equal(sk.U, other.U) and a.passes == b.passes
     assert all(np.array_equal(getattr(a, k), getattr(b, k)) for k in ("push", "shear", "wall", "response", "wall_response"))
+
+
+# -- the droplet momentum sink at the top cells (sub-plan 21 Task 7) ----------------------------------------------------
+def test_a_release_takes_liquid_first_from_the_top_cells_slurry_only_with_the_flag_and_never_rigid_metal():
+    """Spec 8: the release leaves from the top cells, liquid first; exposed slurry only under --slurry-spray on (here
+    `slurry`); nothing passes through rigid metal (spec 6.5). Six columns of 10 um cells, every cell moving: liquid over
+    slurry over rigid metal; slurry exposed over liquid; a floating crust over liquid; rigid metal alone; two liquid cells
+    over rigid metal; and liquid to the bottom. Asked for two and a half cells each, the release takes them from the top
+    of the liquid, at most the liquid there is, and with the flag the exposed slurry and the liquid it covers; asked for
+    everything, it takes the run of liquid cells from the top, or with the flag every flowing cell down to the first
+    rigid one, and a column liquid to its bottom is emptied. The metal leaves at its own cells' specific enthalpy and
+    velocity, with its solid by their liquid fractions; a partly taken cell keeps its temperature, specific enthalpy and
+    velocity, and every cell beneath is untouched (the books to 1.9e-16 and the cells to 2.2e-16, measured
+    2026-10-24)."""
+    col, _ = strip_grid(3, 2)
+    temps = [[LIQUID] * 3 + [SLURRY] * 3 + [RIGID] * 4, [SLURRY] * 2 + [LIQUID] * 3 + [RIGID] * 4,
+             [RIGID] * 2 + [LIQUID] * 4 + [RIGID] * 4, [RIGID] * 8, [LIQUID] * 2 + [RIGID] * 6, [LIQUID] * 7]
+    taken = {(False, 2.5): [2.5, 0, 0, 0, 2, 2.5], (True, 2.5): [2.5, 2.5, 0, 0, 2, 2.5],
+             (False, 1e3): [3, 0, 0, 0, 2, 7], (True, 1e3): [6, 5, 0, 0, 2, 7]}       # cells taken from each column's top
+    rng = np.random.default_rng(7)
+    for (slurry, ask), cells in taken.items():
+        sk = stacked(col, temps)
+        sk.U = np.where(sk.valid[..., None], rng.normal(0.0, 0.2, sk.U.shape), 0.0)
+        area, m, H, T, count, U, _ = sk.state()
+        P0, f_l = sk.column_momentum(), sk.props.liquid_fraction(T)
+        book = skin.release_top(sk, ask * m[:, 0], slurry=slurry)
+        for c, e in enumerate(cells):
+            whole, part = int(e), e - int(e)
+            t = np.zeros(count[c])
+            t[:whole] = m[c, :whole]
+            t[whole:whole + 1] = part * m[c, whole:whole + 1]
+            assert book.mass[c] == pytest.approx(t.sum(), rel=4e-16, abs=0.0)
+            assert book.enthalpy[c] == pytest.approx((t * H[c, :count[c]] / m[c, :count[c]]).sum(), rel=4e-16, abs=0.0)
+            assert book.solid[c] == pytest.approx((t * (1.0 - f_l[c, :count[c]])).sum(), rel=4e-16, abs=0.0)
+            tU = t[:, None] * U[c, :count[c]]
+            np.testing.assert_allclose(book.momentum[c], tU.sum(axis=0), rtol=0.0, atol=4e-16 * np.abs(tU).sum())
+            assert sk.count[c] == count[c] - whole
+            left = slice(whole, count[c])
+            assert np.array_equal(sk.T[c, :sk.count[c]], T[c, left]) and np.array_equal(sk.U[c, :sk.count[c]], U[c, left])
+            np.testing.assert_allclose(sk.m[c, :sk.count[c]], m[c, left] - t[left], rtol=4e-16, atol=0.0)
+            np.testing.assert_allclose(sk.H[c, :sk.count[c]], H[c, left] * (1.0 - t[left] / m[c, left]), rtol=5e-16,
+                                       atol=0.0)
+        assert np.all(book.solid[[0, 4, 5]] == 0.0) if not slurry else book.solid[1] > 0.0
+        np.testing.assert_array_equal(book.momentum, sk.momentum_out)
+        np.testing.assert_allclose(sk.column_momentum() + book.momentum, P0, rtol=1e-13, atol=1e-30)
+    assert sk.count[5] == 0                                                         # liquid to its bottom: emptied
+
+
+def test_what_a_release_takes_leaves_each_column_short_of_exactly_its_book(coarse_sphere_mesh):
+    """On the test sphere's mixed skins -- liquid, slurry, crusts and anchored rigid metal in every arrangement, every
+    cell moving and carrying a sensitivity -- a release of up to four top cells' mass per column, without and with the
+    slurry: each column's mass, enthalpy and momentum fall by exactly what the book says left (to 1.6e-16, 2.3e-16 and
+    2.0e-16 of the largest column's, measured 2026-10-24), the cells left keep their sensitivity, and no rigid cell is
+    touched -- the rigid metal in every column is as it was, to 2.2e-16 of the largest; liquid alone has no solid, and
+    with the slurry 6.4 % of the released mass is solid."""
+    for slurry in (False, True):
+        sk, col, rng = mixed_sphere_skins(coarse_sphere_mesh)
+        sk.S = np.where(sk.valid, rng.uniform(0.0, 1.0, sk.m.shape), 0.0)
+        M0, H0, P0 = sk.column_mass(), sk.H.sum(axis=1), sk.column_momentum()
+        rigid = lambda s: np.where(s.valid & ~np.isfinite(skin.cell_viscosity(s, 0.0)), s.m, 0.0).sum(axis=1)
+        R0, S0 = rigid(sk), [sk.S[c, :sk.count[c]].copy() for c in range(col.n)]
+        count0 = sk.count.copy()
+        dm = rng.uniform(0.0, 4.0, col.n) * sk.m[:, 0]
+        book = skin.release_top(sk, dm, slurry=slurry)
+        assert book.mass.sum() > 0.0 and np.all(book.mass <= dm * (1.0 + 1e-15))
+        np.testing.assert_allclose(sk.column_mass(), M0 - book.mass, rtol=0.0, atol=4e-16 * M0.max())
+        np.testing.assert_allclose(sk.H.sum(axis=1), H0 - book.enthalpy, rtol=0.0, atol=5e-16 * np.abs(H0).max())
+        np.testing.assert_allclose(sk.column_momentum(), P0 - book.momentum, rtol=0.0, atol=4e-16 * np.abs(P0).max())
+        np.testing.assert_allclose(rigid(sk), R0, rtol=0.0, atol=5e-16 * R0.max())
+        for c in range(col.n):
+            assert np.array_equal(sk.S[c, :sk.count[c]], S0[c][count0[c] - sk.count[c]:])
+        assert (book.solid > 0.0).any() == slurry and np.all(book.solid <= book.mass)
+
+
+def released(skins, rate, A_top, dt, slurry=False):
+    """A release of `rate` [kg/(m2 s)] from the top cells over dt, by hand."""
+    return skin.release_top(skins, rate * A_top * dt, slurry=slurry)
+
+
+def test_a_flow_sub_step_releases_after_the_transport_and_its_books_close_with_the_droplets_momentum(coarse_sphere_mesh):
+    """Spec 7.1's order within each flow sub-step -- the flow, the sideways exchange, the spray: `flow_substeps` with a
+    release calls `push`, `momentum_step`, `thickness_step`, `transport` and then `release_top` with the rate the
+    release gives at that moment times the top area and the flow sub-step, to the last bit, at one cycle and three, and
+    books what left: the columns' momentum changes by the push's, the shear's, the walls' and what the transport carried,
+    less the momentum the droplets took at their cells' own velocities (fact 29's sink), to 2.1e-12 and 1.1e-12 of the
+    largest impulse; the mass and the enthalpy fall by what was released, to 1.2e-16 and 4.1e-16 (measured 2026-10-24,
+    the books of 6c's test on these skins). The slurry is released too: the transport mixes these skins' random cells,
+    at a Courant number of 4.4, into slurry at their tops (3 of the 1 177 liquid tops are left after one flow sub-step),
+    so that a release of liquid alone would take almost nothing here."""
+    from test_reentry_model_columns import sphere_columns
+    col, _ = sphere_columns(coarse_sphere_mesh)
+    rate = np.random.default_rng(9).uniform(0.0, 8.0, col.n)
+    for cycles, slurry in ((1, True), (3, True)):
+        first, _, _ = mixed_sphere_skins(coarse_sphere_mesh)
+        second, _, _ = mixed_sphere_skins(coarse_sphere_mesh)
+        flow = sphere_flow(col, cycles, A_top=col.area * 1.01)
+        flow.release, flow.slurry = skin.prescribed_release(rate), slurry
+        P0, M0, E0 = first.column_momentum(), first.mass(), first.energy()
+        book = skin.flow_substeps(first, flow, 0.01)
+        parts = []
+        for k in range(cycles):
+            if k:
+                second.regrid()
+            parts.append(by_hand(second, flow, 0.01 / cycles)
+                         + (released(second, rate, flow.A_top, 0.01 / cycles, slurry),))
+        assert all(np.array_equal(x, y) for x, y in zip(first.state(), second.state()))
+        for key, part in (("released", "mass"), ("released_H", "enthalpy"), ("released_p", "momentum"),
+                          ("released_solid", "solid")):
+            np.testing.assert_allclose(getattr(book, key), sum(getattr(p[3], part) for p in parts), rtol=1e-14, atol=0.0)
+        assert book.released.sum() > 0.0 and book.released_solid.sum() > 0.0
+        assert first.mass() == pytest.approx(M0 - book.released.sum(), rel=3e-16)
+        assert first.energy() == pytest.approx(E0 - book.released_H.sum(), rel=8e-16)
+        impulse = book.push + book.shear + book.wall + book.carried - book.released_p
+        scale = max(np.abs(x).max() for x in (book.push, book.shear, book.wall, book.released_p, P0))
+        assert np.abs(first.column_momentum() - P0 - impulse).max() < 5e-12 * scale
+
+
+def test_a_release_that_empties_columns_leaves_the_flow_and_the_draw_working_and_the_books_closed(monkeypatch):
+    """Columns liquid to their bottom, asked for far more than they hold, are emptied in the first of four flow
+    sub-steps; the flow sub-steps after it land the runoff on them as on any rigid top, the heat sub-step's draw restores
+    their minimum, and the pass's mass, energy and momentum close with the release, to 4.1e-16, 2.2e-16 and 4.1e-16
+    (measured 2026-10-24; the conduction off, so that the layer stays liquid). A negative rate is refused."""
+    monkeypatch.setattr(skin, "_conduct", no_conduction)
+    col, pts = strip_grid(12, 2)
+    x = np.rint(pts[:, 0] / 1e-3).astype(int)
+    deep = (x == 3) | (x == 4)
+    sk = stacked(col, [[LIQUID] * 6 if d else [LIQUID] * 10 + [RIGID] * 5 for d in deep])
+    sk.settings = skin.SkinSettings(**CELLS_04, substep=0.01)
+    T_b, m_min = np.full(sk.n, RIGID), skin.minimum_mass(sk.material, sk.area, sk.settings)
+    flow = strip_flow(col, tau=50.0, cycles=4)
+    emptied = []
+    take = skin.release_top
+
+    def recorded(skins, dm, slurry=False):
+        book = take(skins, dm, slurry)
+        emptied.append(skins.count[deep].copy())
+        return book
+    monkeypatch.setattr(skin, "release_top", recorded)
+    flow.release = skin.prescribed_release(np.where(deep, 1e3, 0.5))
+    E0, M0, P0 = sk.energy(), sk.mass(), sk.column_momentum()
+    book = skin.advance(sk, 0.03, np.full(sk.n, 3e6), sk.area, 300.0, 0.3, T_b, m_min, T_b1=T_b, remove=False, flow=flow)
+    f = book.flow
+    assert np.all(emptied[0] == 0) and np.all(book.drawn[deep] > 0.0) and np.all(sk.column_mass() >= m_min * (1 - 1e-12))
+    assert sk.mass() == pytest.approx(M0 - f.released.sum() + book.drawn.sum(), rel=1e-15)
+    gained = book.E_top - book.E_rad - book.E_base + book.E_side - f.released_H + book.drawn_H
+    assert sk.energy() == pytest.approx(E0 + gained.sum(), rel=5e-16)
+    scale = max(np.abs(v).max() for v in (f.push, f.shear, f.wall, f.released_p))
+    assert np.abs(sk.column_momentum() - P0 - (f.push + f.shear + f.wall + f.carried - f.released_p)).max() < 1e-15 * scale
+    flow.release = lambda skins: np.full(skins.n, -1.0)
+    with pytest.raises(ValueError):
+        skin.flow_substeps(sk, flow, 0.01)
+
+
+def test_with_no_release_or_nothing_it_may_take_the_flow_is_unchanged_to_the_last_bit(coarse_sphere_mesh):
+    """With no release (the default) the flow sub-steps are Task 6d's; a release of nothing -- a zero rate, or a rate
+    only where a column's top cell is not liquid at the moment it is asked, without the slurry -- gives the same cells
+    and the same books to the last bit, the released metal zero. (Asked once from the starting tops instead, the second
+    rate releases metal: the transport lands liquid on some of those tops within the flow sub-step.)"""
+    from test_reentry_model_columns import sphere_columns
+    col, _ = sphere_columns(coarse_sphere_mesh)
+    dry = lambda s: np.where(s.T[:, 0] < s.material.T_liquidus, 5.0, 0.0)              # no liquid at the top
+    probe, _, _ = mixed_sphere_skins(coarse_sphere_mesh)
+    assert 0.0 < (dry(probe) > 0.0).mean() < 1.0
+    out = []
+    for release in (None, skin.prescribed_release(0.0), dry):
+        sk, _, _ = mixed_sphere_skins(coarse_sphere_mesh)
+        flow = sphere_flow(col, 2)
+        flow.release = release
+        out.append((sk, skin.flow_substeps(sk, flow, 0.01)))
+    (a, A) = out[0]
+    for b, B in out[1:]:
+        assert all(np.array_equal(x, y) for x, y in zip(a.state(), b.state()))
+        for key in ("push", "shear", "wall", "carried", "courant", "growth", "passes", "sub_steps", "cycles"):
+            assert np.array_equal(getattr(A, key), getattr(B, key))
+        assert not (B.released.any() or B.released_H.any() or B.released_p.any() or B.released_solid.any())
+    assert not A.released.any() and A.released_p.shape == (col.n, 3)
+
+
+def test_the_droplet_momentum_ratio_is_fact_29_s_v_melt_delta_over_nu_on_a_sheared_layer():
+    """Fact 29: Pi = mdot V_s / tau_gas, the droplets' momentum flux against the gas shear driving the film. By hand:
+    liquid aluminium 0.234 mm deep on a strip of 1 mm columns, sheared at 20 Pa along it at its steady Couette profile
+    (u = tau s / mu at every cell's centre, s its height above the rigid front: the momentum step keeps it and the
+    transport carries it along the strip unchanged), released at mdot = rho_l v_melt from the top in one 1 ms flow
+    sub-step. The top cell leaves at u_0 = tau (H - dz_0 / 2) / mu, so every column's Pi_i = mdot u_0 / tau =
+    mdot (H - dz_0 / 2) / mu, fact 29's v_melt delta_m / nu with delta_m the top cell's height: 0.42 at v_melt = 1 mm/s,
+    flagged above 0.2, and 0.042 at 0.1 mm/s, to 2.2e-15 (measured 2026-10-24); the body's Pi, the shear-weighted mean,
+    is the same; without a release Pi is zero; a top the gas does not shear has no Pi of its own (infinite) and is
+    neither in the body's Pi nor flagged."""
+    for v_melt, flagged in ((1e-3, True), (1e-4, False)):
+        sk, col, _ = strip_layer(nx=12)
+        flow = strip_flow(col, dp_ds=0.0, sigma=0.0, tau=20.0)
+        dz = sk.dz()[:, :20]
+        H, dz0 = dz[0].sum(), dz[0, 0]
+        sk.U[:, :20, 0] = 20.0 * (H - (np.cumsum(dz, axis=1) - 0.5 * dz)) / MU_L       # the steady Couette profile
+        mdot = RHO_L * v_melt
+        flow.release = skin.prescribed_release(mdot)
+        book = skin.flow_substeps(sk, flow, 1e-3, cycles=1)
+        Pi, each, count = skin.droplet_momentum_ratio(book.released_p, book.shear)
+        expected = mdot * (H - 0.5 * dz0) / MU_L
+        np.testing.assert_allclose(each, expected, rtol=5e-15)
+        assert Pi == pytest.approx(expected, rel=5e-15)
+        assert expected == pytest.approx(v_melt * (H - 0.5 * dz0) / (MU_L / RHO_L))              # v_melt delta_m / nu
+        assert count == (col.n if flagged else 0) and (expected > skin.PI_FLAG) == flagged
+    Pi, each, count = skin.droplet_momentum_ratio(np.zeros((3, 3)), np.ones((3, 3)))
+    assert Pi == 0.0 and not each.any() and count == 0                                 # nothing released
+    Pi, each, count = skin.droplet_momentum_ratio(np.array([[1e-9, 0, 0], [0, 0, 0], [1e-9, 0, 0]]),
+                                                  np.array([[0, 0, 0], [0, 0, 0], [0, 1e-8, 0]]))
+    assert Pi == pytest.approx(0.1) and each[0] == np.inf and each[1] == 0.0 and count == 0   # an unsheared top: no Pi
+
+
+def test_the_trial_pass_releases_what_its_real_pass_releases_and_its_base_law_includes_the_release(monkeypatch):
+    """The trial pass runs the flow's sub-steps as the real pass does (6c), the release in them: liquid aluminium on a
+    strip of 1 mm columns heated at 30 MW/m2 (its top stays liquid over the base at 880 K), pushed and sheared, released
+    at 3 kg/(m2 s), the conduction on. A real pass at the trial pass's base temperature releases the same metal in
+    every flow sub-step, to the last bit, and its base heat meets the trial pass's law to 4.0e-17 of the throughput;
+    the trial pass leaves the skins as they were; and the law with the release differs from the law without by 0.20 %
+    of the throughput (measured 2026-10-24)."""
+    seen = []
+    take = skin.release_top
+
+    def recorded(sk, dm, slurry=False):
+        book = take(sk, dm, slurry)
+        seen.append((np.copy(dm), book))
+        return book
+    monkeypatch.setattr(skin, "release_top", recorded)
+    laws = {}
+    for rate in (None, 3.0):
+        sk, col, _ = strip_layer(nx=12)
+        sk.settings = skin.SkinSettings(**CELLS_04, substep=0.01)
+        T_b, m_min, q = np.full(sk.n, RIGID), np.zeros(sk.n), np.full(sk.n, 3e7)
+        flow = strip_flow(col, tau=50.0, cycles=2)
+        flow.release = None if rate is None else skin.prescribed_release(rate)
+        before = sk.state()
+        seen.clear()
+        laws[rate] = skin.advance(sk, 0.05, q, sk.area, 300.0, 0.3, T_b, m_min, remove=False, flow=flow)
+        trial = list(seen)
+        assert all(np.array_equal(x, y) for x, y in zip(sk.state(), before))
+        seen.clear()
+        book = skin.advance(sk, 0.05, q, sk.area, 300.0, 0.3, T_b, m_min, T_b1=T_b, remove=False, flow=flow)
+        assert len(seen) == len(trial) == (10 if rate else 0)
+        for (dm, b), (dm_t, b_t) in zip(seen, trial):
+            assert np.array_equal(dm, dm_t) and all(np.array_equal(getattr(b, k), getattr(b_t, k))
+                                                    for k in ("mass", "enthalpy", "momentum", "solid"))
+        a, s = laws[rate]
+        throughput = np.abs(book.E_top).sum() + np.abs(book.E_base).sum()
+        assert np.abs(book.E_base / 0.05 - (a + s * T_b)).sum() * 0.05 < 1e-16 * throughput
+        assert rate is None or book.flow.released.sum() > 0.0
+    (a0, s0), (a1, s1) = laws[None], laws[3.0]
+    change = np.abs((a1 + s1 * T_b) - (a0 + s0 * T_b)).sum() * 0.05 / throughput
+    assert 1e-3 < change < 4e-3
```

- [ ] **Step 2: The skin body's test** — append to `tests/test_reentry_model_skinbody.py`:

```diff
--- a/tests/test_reentry_model_skinbody.py
+++ b/tests/test_reentry_model_skinbody.py
@@ -500,3 +500,31 @@
     assert max(u for _, u in calls[-1][1]) <= stats["flow_cycles_max"] <= max(used) and stats["flow_courant_max"] > 0.0
     assert b.mass() == pytest.approx(M0, rel=1e-13) and abs(b.energy_balance_residual()) < 1e-8
     assert stats["momentum_residual"] < 1e-12
+
+
+# -- the droplet momentum sink (sub-plan 21 Task 7) ---------------------------------------------------------------------
+def test_a_skin_body_with_a_prescribed_release_keeps_its_books_and_reports_the_droplet_momentum_ratio(coarse_sphere_mesh):
+    """Sub-plan 21 Task 7, the placeholder for stage 4's spraying: a release of 2 kg/(m2 s) from every column's top, liquid
+    first, as sub-plan 22's Girin rate will be given, on the 4 mm test sphere heated at 3 MW/m2 from 880 K with the flow on
+    (as above: two 0.25 s steps, one flow sub-step a heat sub-step). The tops melt in the second step and the release
+    takes the liquid as it forms, 0.53 g. It leaves the body: the mass left and the mass released are the starting mass
+    to the last bit, the energy books close with the released enthalpy to 9.8e-13 of the heat absorbed, and the
+    momentum balance, the droplets' momentum lost with the metal, to 2.3e-17. The history reports the released mass,
+    fact 29's droplet-momentum ratio Pi over the macro step, which is its definition on the step's books (0.0048 here,
+    with the tops barely liquid and slow), and the columns flagged above 0.2, none (measured 2026-10-24). A release
+    needs the flow: with it off the heat-only skins remove the liquid as it forms."""
+    with pytest.raises(ValueError):
+        make(coarse_sphere_mesh, release=skin.prescribed_release(2.0))
+    b = make(coarse_sphere_mesh, flow=True, remove=False, flow_cycles=1, release=skin.prescribed_release(2.0))
+    M0 = b.mass()
+    for k in range(2):
+        b.advance(0.25 * (k + 1), 0.25, front_loads(b, 3e6), forcing=newtonian(b))
+    stats = b.skin_stats()
+    assert b.removed_mass == 0.0 and b.released_mass > 0.0 and stats["released_mass_kg"] == b.released_mass
+    assert b.mass() + b.released_mass == pytest.approx(M0, rel=1e-13)
+    assert abs(b.energy_balance_residual()) < 1e-8 and b.released_enthalpy > 0.0
+    assert stats["momentum_residual"] < 1e-12 and np.abs(b.droplet_momentum).max() > 0.0
+    Pi, _, flagged = skin.droplet_momentum_ratio(b.droplet_momentum, b.droplet_shear)
+    assert stats["droplet_momentum_ratio"] == Pi > 0.0 and stats["droplet_momentum_flagged"] == flagged
+    off = make(coarse_sphere_mesh, flow=True, remove=False, flow_cycles=1)
+    assert off.skin_stats()["droplet_momentum_ratio"] == 0.0 and off.skin_stats()["released_mass_kg"] == 0.0
```

- [ ] **Step 3: Run them to see them fail.** `"$PY" -m pytest tests/test_reentry_model_skin.py
  tests/test_reentry_model_skinbody.py -q -k "release or droplet or nothing_it_may_take or empties"`. Expected: 8
  failed, 99 deselected --
  - `AttributeError: module 'reentry_model.skin' has no attribute 'release_top'` (4: the take, the books, the emptied
    columns and the trial pass);
  - `AttributeError: module 'reentry_model.skin' has no attribute 'prescribed_release'` (4: the flow sub-steps, no
    release, Pi and the skin body).

- [ ] **Step 4: The release from the top cells, its books in the flow and Pi** in `reentry_model/skin.py`:

```diff
--- a/reentry_model/skin.py
+++ b/reentry_model/skin.py
@@ -151,7 +151,8 @@
     is an intensive per-cell field every cell operation carries (the trial pass's dT/dT_b). `U` (n, K, 3) is every
     cell's velocity [m/s] (sub-plan 21 Task 2): part of the cells' state, carried by every operation with its momentum
     exact -- split halves keep it, merged cells take the mass-weighted mean, metal drawn from the 3D model arrives at
-    rest -- and `momentum_out` (n, 3) [kg m/s] is what the last `take_top` or `take_bottom` took away with the metal.
+    rest -- and `momentum_out` (n, 3) [kg m/s] is what the last `take_top`, `take_top_mass` or `take_bottom` took away
+    with the metal.
     `flow_memory` (3,) is what the adaptive count of flow sub-steps remembers of the last heat sub-step's flow (sub-plan
     21 Task 6d, `flow_cycles`): its largest Courant number and Rayleigh-Taylor growth and its flow sub-step's length [s],
     zero before any flow; part of the state, so that a trial pass, a repeated or a shortened part starts from it."""
@@ -387,6 +388,30 @@
             U = np.take_along_axis(self.U, src[..., None], axis=1)
             self._install(shift(self.m), shift(self.H), self.count - k, S, shift(self.T), U)
         return mass, H
+
+    def take_top_mass(self, dm, admit):
+        """Remove dm >= 0 [kg] per column from its top down through the run of admissible cells from its top (`admit`,
+        (n, K) bool), at most what that run holds, each cell's metal at its own specific enthalpy and velocity: a partly
+        taken cell keeps both, its temperature and its sensitivity, and the cells beneath are untouched (sub-plan 21 Task
+        7: the release of spec 8, which thins the top cell until it merges). Returns (mass, enthalpy, solid) taken per
+        column, solid the mass of solid in it by the cells' liquid fractions, and leaves the momentum that went in
+        `momentum_out`."""
+        dm = np.broadcast_to(np.asarray(dm, dtype=float), (self.n,))
+        if not np.all(np.isfinite(dm) & (dm >= 0.0)):
+            raise ValueError("a release from a skin's top must be a finite, non-negative mass")
+        run = np.cumprod(self.valid & admit, axis=1).astype(bool)
+        above = np.concatenate([np.zeros((self.n, 1)), np.cumsum(self.m, axis=1)[:, :-1]], axis=1)   # mass above each
+        taken = np.where(run, np.clip(dm[:, None] - above, 0.0, self.m), 0.0)
+        frac = np.divide(taken, self.m, out=np.zeros_like(taken), where=self.m > 0.0)
+        mass, H = taken.sum(axis=1), (self.H * frac).sum(axis=1)
+        solid = (taken * (1.0 - self.props.liquid_fraction(self.T))).sum(axis=1)
+        momentum = (taken[..., None] * self.U).sum(axis=1)
+        whole = np.cumprod(run & (taken >= self.m), axis=1).sum(axis=1)        # the cells wholly taken, from the top
+        part = (taken > 0.0) & (np.arange(self.m.shape[1])[None, :] >= whole[:, None])
+        self.m, self.H = np.where(part, self.m - taken, self.m), np.where(part, self.H * (1.0 - frac), self.H)
+        self.take_top(whole)
+        self.momentum_out = momentum
+        return mass, H, solid
 
     def add_bottom(self, dm, h, s=1.0, u=None):
         """Metal dm >= 0 [kg] at specific enthalpy h [J/kg] joins the bottom cell of every column (an emptied column gets
@@ -650,9 +675,10 @@
 def advance(sk, dt, q_top, A_top, T_amb, eps, T_b0, m_min, T_b1=None, side=None, remove=True, g_base=None, flow=None):
     """The macro step's sub-steps for every column (spec 7.2). Each sub-step: the implicit conduction, through the depth
     and, with `side`, sideways (`side(sk, previous)` makes the sub-step's `Sideways` map from the last one, Task 5),
-    then, with `flow` (a `Flow`), the flow's sub-steps (`flow_substeps`: heat, then flow, sideways exchange and -- in
-    stage 4 -- spray, spec 7.1; sub-plan 21 Task 6c), then -- with `remove` and a melting material -- the contiguous top
-    cells at or above T_feed leave, taking their momentum with them, a column below its minimum mass m_min [kg] draws
+    then, with `flow` (a `Flow`), the flow's sub-steps (`flow_substeps`: heat, then flow, sideways exchange and spray,
+    spec 7.1; sub-plan 21 Task 6c, the spray a release from the top cells with `flow.release` since Task 7), then --
+    with `remove` and a melting material -- the contiguous top cells at or above T_feed leave, taking their momentum
+    with them, a column below its minimum mass m_min [kg] draws
     back to it, and the cells are split and merged. The draw takes rigid metal (spec 6.6): at the base temperature, but
     never above T_rigid -- a 3D base node can overshoot it within a macro step (on 4 mm elements under 6 MW/m2 it
     reached 930 K in one 0.5 s step, found while planning), and drawing at its temperature would draw liquid that the
@@ -1618,8 +1644,8 @@
 FLOW_GROWTH = 0.02                # and the largest Rayleigh-Taylor growth dt s of one: backward Euler's growth rate is
                                   # then within about dt s / 2 = 1 % (6b's cosine met it within 0.42 % at 0.016)
 FLOW_RELAX = 2.0                  # a flow sub-step at most this many times the last heat sub-step's
-FLOW_SECONDS = {"push": 0.0, "momentum": 0.0, "thickness": 0.0, "transport": 0.0, "regrid": 0.0}   # wall time, for
-                                  # measuring the cost (as KRYLOV counts the coupled solves)
+FLOW_SECONDS = {"push": 0.0, "momentum": 0.0, "thickness": 0.0, "transport": 0.0, "regrid": 0.0, "release": 0.0}
+                                  # wall time, for measuring the cost (as KRYLOV counts the coupled solves)
 
 
 def flow_cycles(memory, dt):
@@ -1650,7 +1676,15 @@
     in the body's frame, minus the body's own, `v_hat` the flight direction, `rho_l` [kg/m3] and `sigma` [N/m] the
     liquid's density and surface tension, `A_top` [m2, (n,)] the top areas the shear acts over, and `cycles`, the flow
     sub-steps in each heat sub-step: None (the default since Task 6d) chooses them each heat sub-step (`flow_cycles`),
-    a number fixes them. `plan` is what a trial pass chose, for its real pass (`advance`)."""
+    a number fixes them. `plan` is what a trial pass chose, for its real pass (`advance`).
+
+    `release` is the spray's placeholder (sub-plan 21 Task 7), what sub-plan 22's spraying supplies: None releases
+    nothing; otherwise a callable `release(sk) -> (n,)` [kg/(m2 s)], the release rate per unit of each column's top area
+    A_top at the skins' present state, evaluated in every flow sub-step after the transport and taken from the top cells
+    by `release_top` -- liquid first, the exposed slurry too where `slurry` (spec 8's `--slurry-spray on`, a setting of
+    sub-plan 22's). It is called in the trial pass, in every repeat of a part and in parts run again, so it keeps no
+    records of its own: what was released is in the pass's FlowBook (`released`, `released_H`, `released_p`,
+    `released_solid`)."""
     columns: object
     dp_ds: np.ndarray
     shear: np.ndarray
@@ -1661,12 +1695,15 @@
     A_top: np.ndarray
     cycles: int = None
     plan: list = None
+    release: object = None
+    slurry: bool = False
 
 
 @dataclass
 class FlowBook:
     """What the flow sub-steps did, per column: impulses [N s] and momentum [kg m/s], (n, 3), which add up to the
-    columns' momentum change, and the sub-steps' worst numbers."""
+    columns' momentum change less what the release took (`released_p`), the sub-steps' worst numbers, and the metal
+    released from the top cells (Task 7; None in a single sub-step's book that released nothing)."""
     push: np.ndarray         # the push's impulse on the moving cells, its thickness terms' implicit change included
     shear: np.ndarray        # the gas shear's on the top cells
     wall: np.ndarray         # the anchored metal's drag (the momentum and thickness steps'), the momentum cells held on
@@ -1677,15 +1714,24 @@
     passes: int = 0          # the thickness step's most linearisations
     sub_steps: int = 0       # the flow sub-steps taken (none where nothing flows)
     cycles: int = 0          # the most flow sub-steps in one heat sub-step (Task 6d)
+    released: np.ndarray = None        # kg (n,): the metal released from the top cells (Task 7)
+    released_H: np.ndarray = None      # J (n,): its enthalpy, at its cells' own specific enthalpies
+    released_p: np.ndarray = None      # kg m/s (n, 3): its momentum, at its cells' own velocities -- fact 29's sink
+    released_solid: np.ndarray = None  # kg (n,): the solid in it, by its cells' liquid fractions
 
     @classmethod
     def zero(cls, n):
-        return cls(*(np.zeros((n, 3)) for _ in range(4)))
+        return cls(*(np.zeros((n, 3)) for _ in range(4)), released=np.zeros(n), released_H=np.zeros(n),
+                   released_p=np.zeros((n, 3)), released_solid=np.zeros(n))
 
     def add(self, other):
-        """Add another book's impulses and keep the worst numbers of both."""
+        """Add another book's impulses and releases and keep the worst numbers of both."""
         for key in ("push", "shear", "wall", "carried"):
             setattr(self, key, getattr(self, key) + getattr(other, key))
+        for key in ("released", "released_H", "released_p", "released_solid"):
+            mine, theirs = getattr(self, key), getattr(other, key)
+            if theirs is not None:
+                setattr(self, key, theirs.copy() if mine is None else mine + theirs)
         self.courant, self.growth = max(self.courant, other.courant), max(self.growth, other.growth)
         self.passes, self.sub_steps = max(self.passes, other.passes), self.sub_steps + other.sub_steps
         self.cycles = max(self.cycles, other.cycles)
@@ -1696,15 +1742,25 @@
     `flow_cycles` chooses from the skins' `flow_memory`, which they then update (sub-plan 21 Task 6d) -- every one in the
     order of 6b's ruling: the push, surface tension only where a column's top flows (`push`, `flowing_tops`); the
     momentum through the depth with its response (`momentum_step`); the thickness terms implicit in the heights at the
-    end of the transport (`thickness_step`); the transport with the faces' own push (`transport`, regrid=False). The
-    cells are split and merged between the flow sub-steps, and the heat sub-step's own regrid follows the last (spec
-    6.2: once after each sub-step's exchanges). Updates the cells; returns a FlowBook. Where nothing flows nothing moves
-    and no sub-step is taken: any velocity left in the rigid cells is booked to the wall and cleared, as `momentum_step`
-    books an anchored cell's, and the skins forget the flow.
+    end of the transport (`thickness_step`); the transport with the faces' own push (`transport`, regrid=False); then,
+    with `flow.release`, the spray's placeholder (sub-plan 21 Task 7): the release the rate gives at that moment over
+    the flow sub-step, rate A_top h, from the top cells (`release_top`) -- spec 7.1's heat, flow, sideways exchange and
+    spray, the last three in every flow sub-step. The cells are split and merged between the flow sub-steps, and the
+    heat sub-step's own regrid follows the last (spec 6.2: once after each sub-step's exchanges). Updates the cells;
+    returns a FlowBook. Where nothing flows nothing moves, nothing can be released and no sub-step is taken: any
+    velocity left in the rigid cells is booked to the wall and cleared, as `momentum_step` books an anchored cell's, and
+    the skins forget the flow.
 
     The flow's sub-step was fixed in 6c at 1 ms (FLOW_CYCLES) and is chosen since 6d from the last heat sub-step's
     Courant number and growth, against targets set from the sub-step study on the 100 mm sphere (FLOW_COURANT,
-    FLOW_GROWTH)."""
+    FLOW_GROWTH). The release runs in every flow sub-step, not once a heat sub-step (Task 7's ruling): the top cells
+    cross a column for each unit of their Courant number, which the count keeps near 6 a flow sub-step, so they cross
+    some 25 columns of 2 mm in a heat sub-step at the 100 mm sphere's 5 m/s and up to 240 at the count's cap; a release
+    once a heat sub-step would spray metal from where the runoff has taken it rather than from where its rate was
+    found, and fact 29's split between the runoff and the stripping would be first order in the heat sub-step rather
+    than in the flow's. Its cost: 0.05 s of processor time a release on the 100 mm
+    sphere's 9 417 columns, against about 2 s for a 1 ms flow sub-step of its driven melting state (measured
+    2026-10-24)."""
     n = sk.n
     book = FlowBook.zero(n)
     adaptive = flow.cycles is None
@@ -1736,12 +1792,90 @@
         started = time.perf_counter()
         tr = transport(sk, col, h, regrid=False, faces=tb.faces)
         clock("transport", started)
-        book.add(FlowBook(mb.push + tb.push, mb.shear, mb.wall + tb.wall + tr.wall, tr.momentum, tr.courant, tb.growth,
-                          tb.passes, 1))
+        step = FlowBook(mb.push + tb.push, mb.shear, mb.wall + tb.wall + tr.wall, tr.momentum, tr.courant, tb.growth,
+                        tb.passes, 1)
+        if flow.release is not None:
+            started = time.perf_counter()
+            rate = np.broadcast_to(np.asarray(flow.release(sk), dtype=float), (n,))
+            if not np.all(np.isfinite(rate) & (rate >= 0.0)):
+                raise ValueError("a release rate must be finite and non-negative")
+            rb = release_top(sk, rate * flow.A_top * h, flow.slurry)
+            step.released, step.released_H = rb.mass, rb.enthalpy
+            step.released_p, step.released_solid = rb.momentum, rb.solid
+            clock("release", started)
+        book.add(step)
     book.cycles = cycles
     if adaptive:
         sk.flow_memory = np.array([book.courant, book.growth, h])
     return book
+
+
+# -- the droplet momentum sink at the top cells (sub-plan 21 Task 7) ----------------------------------------------------
+PI_FLAG = 0.2                     # fact 29: above this the droplets' momentum is no longer small against the gas shear
+                                  # that drives the film, and the sink must be solved with the profile (as here it is)
+
+
+@dataclass
+class ReleaseBook:
+    """What a release took from the top cells, per column."""
+    mass: np.ndarray         # kg (n,)
+    enthalpy: np.ndarray     # J (n,), at the cells' own specific enthalpies
+    momentum: np.ndarray     # kg m/s (n, 3), at the cells' own velocities: the droplets' momentum, fact 29's sink
+    solid: np.ndarray        # kg (n,): the solid in it, by the cells' liquid fractions (the droplets' solid fraction)
+
+
+def release_top(sk, dm, slurry=False):
+    """Release dm [kg, (n,)] per column from its top cells (spec 8; sub-plan 21 Task 7, the placeholder for stage 4's
+    spraying, whose rate sub-plan 22 supplies): liquid first -- the run of fully liquid cells from the top, at or above
+    the liquidus -- and, with `slurry` (spec 8's `--slurry-spray on`), the exposed slurry too, every flowing cell down
+    to the first rigid one; never rigid metal (spec 6.5: nothing passes through it), so a top that is slurry without
+    the flag, a floating crust or anchored metal releases nothing. At most what those cells hold leaves; the order
+    within them, and each one's rate, is the cascade rule's (sub-plan 22). The metal leaves at its cells' own specific
+    enthalpy and velocity (`Skins.take_top_mass`): the droplets carry momentum away from the top cells at their release
+    speed (spec 6.4), the speed of the metal they are made of, the top of the resolved profile -- fact 29's sink
+    mdot V_s with V_s the sheared layer's surface speed -- so each column's mass, enthalpy and momentum fall by exactly
+    what the book says left. Updates the cells; returns a ReleaseBook."""
+    flowing = sk.valid & np.isfinite(cell_viscosity(sk, 0.0))
+    admit = flowing if slurry else flowing & (sk.T >= sk.material.T_liquidus)
+    mass, H, solid = sk.take_top_mass(dm, admit)
+    return ReleaseBook(mass, H, sk.momentum_out.copy(), solid)
+
+
+def prescribed_release(rate):
+    """A release (`Flow.release`) at the rate given [kg/(m2 s), per column or one for all] whatever the skins' state:
+    the tests' device, and any prescribed spraying's, which sub-plan 22 replaces with Girin's rate per column."""
+    rate = np.array(rate, dtype=float)
+    return lambda sk: np.broadcast_to(rate, (sk.n,))
+
+
+def droplet_momentum_ratio(released, shear):
+    """Fact 29's droplet-momentum ratio Pi = mdot V_s / tau_gas: the momentum the droplets take from the film against the
+    gas shear that drives it, from the books of a time (a macro step in the history): released (n, 3) [kg m/s] the
+    momentum the release took from each column's top cells (`FlowBook.released_p`), shear (n, 3) [N s] the gas shear's
+    impulse on its top over the same time (`FlowBook.shear`), both integrated over the column's top and the time.
+    Returns (Pi, each, flagged): each column's Pi_i = |released_i| / |shear_i|, zero where nothing was released and
+    infinite where metal left a top the gas did not shear; Pi = sum |released_i| / sum |shear_i| over the sheared
+    columns, the share of the gas shear's momentum that left with droplets -- each column's ratio weighted by its shear,
+    so that columns barely sheared do not dominate it -- zero without a release; and flagged, the number of sheared
+    columns above PI_FLAG, where fact 29 calls the stripping and the runoff no longer separable. Pi is defined against
+    the shear, so a top the gas does not shear has none and is neither in Pi nor flagged: Girin's modes release only on
+    the windward face, where the gas shears, but a prescribed release on the leeward runoff (moving at metres a second
+    on the 100 mm sphere's rigid rear, under no shear) would otherwise make Pi infinite. On a layer at its steady
+    Couette profile released at the top, Pi_i = mdot (H - dz_0 / 2) / mu: fact 29's v_melt delta_m / nu with delta_m the
+    top cell's height above the rigid front.
+
+    Measured 2026-10-24 on the 100 mm sphere's driven melting state (the front's top 0.3 mm liquid, driven from rest by
+    the flight's forcing at 93 s for 20 ms), released at rho_l v_melt over the windward face for 10 ms more: Pi = 0.145
+    at v_melt = 1 mm/s, 1 146 of the 4 705 sheared columns above 0.2 (the largest 0.26), and 0.015 at 0.1 mm/s, none.
+    Fact 29 estimated 0.48 and 0.048, taking the film at its steady Couette speed tau delta_m / mu; the 0.3 mm layer,
+    whose viscous time h^2 / nu is 0.17 s, has not reached it in 30 ms (its top at about 4.3 m/s under 64 Pa, against
+    15 m/s)."""
+    p, J = np.linalg.norm(released, axis=1), np.linalg.norm(shear, axis=1)
+    sheared = J > 0.0
+    with np.errstate(divide="ignore", invalid="ignore"):
+        each = np.where(p > 0.0, p / J, 0.0)
+    Pi = float(p[sheared].sum() / J.sum()) if sheared.any() else 0.0
+    return Pi, each, int((each[sheared] > PI_FLAG).sum())
 
 
 # -- the base (Task 7) -------------------------------------------------------------------------------------------------
```

- [ ] **Step 5: The release in the skin body** in `reentry_model/skinbody.py`:

```diff
--- a/reentry_model/skinbody.py
+++ b/reentry_model/skinbody.py
@@ -48,7 +48,13 @@
 heat sub-step's Courant number and growth (sub-plan 21 Task 6d; `flow_cycles` fixes them). With the flow on the liquid
 stays and flows
 (`remove` False in a run; spraying is stage 4's). The momentum books: the skins' momentum changes by the flow's impulses
-less what leaves with removed or handed-down metal, and `momentum_residual` reports what they miss."""
+less what leaves with removed or handed-down metal, and `momentum_residual` reports what they miss.
+
+The droplet momentum sink (sub-plan 21 Task 7; `release`, None by default): the placeholder for stage 4's spraying, a
+release rate per column (`skin.Flow.release`, what sub-plan 22 supplies from Girin's model) taken in every flow
+sub-step from the top cells, liquid first, the exposed slurry too with `slurry_spray`. The released metal leaves the
+body with its cells' enthalpy and momentum, booked as the removed liquid is, and the history reports the released mass
+and fact 29's droplet-momentum ratio Pi over each macro step (`skin.droplet_momentum_ratio`)."""
 import dataclasses
 import math
 import time
@@ -133,18 +139,24 @@
     `remove`: liquid leaves as it forms (stage 2); `demise_fraction`: the run ends below this share of the initial mass;
     `flow`: the metal flows (sub-plan 21), `flow_cycles` flow sub-steps in each heat sub-step (None: as many as
     `skin.flow_cycles` chooses each heat sub-step, Task 6d), driven by `surface_flow` (a `surface_flow.SurfaceFlow`, its
-    defaults if None) from each macro step's aero state."""
+    defaults if None) from each macro step's aero state; `release`, with the flow on, a release from the top cells in
+    every flow sub-step (`skin.Flow.release`: a callable of the skins giving a rate [kg/(m2 s)] per unit top area per
+    column; sub-plan 21 Task 7), the exposed slurry released too where `slurry_spray`."""
 
     def __init__(self, mesh, material, solver, mass_kg, settings=None, T0=300.0, emissivity=None, T_ambient=0.0,
                  v_hat=(1.0, 0.0, 0.0), deform_depth=movemesh.DEFORM_DEPTH, remove=True, demise_fraction=0.01,
                  repeat_rule=DEFAULT_REPEAT_RULE, part_max=PART_MAX, flow=FLOW_DEFAULT, flow_cycles=None,
-                 surface_flow=None):
+                 surface_flow=None, release=None, slurry_spray=False):
         if repeat_rule not in REPEAT_RULES:
             raise ValueError("repeat_rule must be one of {}, got {!r}".format(REPEAT_RULES, repeat_rule))
         self.flow_cycles = flow_cycles
         if not (flow_cycles is None or isinstance(flow_cycles, (int, np.integer)) and flow_cycles >= 1):
             raise ValueError("flow_cycles must be None or an integer of at least 1, got {!r}".format(flow_cycles))
         self.flow_on, self.surface_flow = bool(flow), surface_flow
+        if release is not None and not self.flow_on:
+            raise ValueError("a release from the top cells needs the flow on: without it the heat-only skins remove the "
+                             "liquid as it forms (stage 2's stand-in for flow and spray)")
+        self.release, self.slurry_spray = release, bool(slurry_spray)
         self.repeat_rule = repeat_rule
         self.part_max = part_max                             # s, or None: no cap on the parts while melting
         self.material, self.solver, self.mass_kg = material, solver, mass_kg
@@ -180,6 +192,9 @@
         self._forcing = self._forcing_prev = None            # this macro step's (start, end) FlowForcing; the last end
         self.impulse, self.momentum_lost = np.zeros(3), np.zeros(3)   # the momentum books since time zero [N s, kg m/s]
         self.impulse_size = 0.0                              # N s: the push's, shear's and walls' impulses, in size
+        self.released_mass = self.released_enthalpy = 0.0    # kg, J: released from the top cells since time zero (Task 7)
+        self.droplet_momentum = np.zeros((self.columns.n, 3))   # kg m/s: what the droplets took from each column, and
+        self.droplet_shear = np.zeros((self.columns.n, 3))      # N s: the gas shear's impulse on its top, this macro step
         self.flow_courant = self.flow_growth = 0.0           # this macro step's worst flow sub-step
         self.flow_passes = self.flow_cycles_max = 0          # and the most flow sub-steps in a heat sub-step
         self._volume = None                                  # m3, the 3D model's volume since the last move (mass)
@@ -231,6 +246,7 @@
             self._forcing, self._forcing_prev = (start, end), end
             self.flow_courant = self.flow_growth = 0.0
             self.flow_passes = self.flow_cycles_max = 0
+            self.droplet_momentum, self.droplet_shear = np.zeros((self.columns.n, 3)), np.zeros((self.columns.n, 3))
         sk, mat = self.skins, self.material
         top = np.where(sk.valid[:, 0], sk.T[:, 0], np.nan)
         if self.part_max is not None and mat.melts:
@@ -280,7 +296,8 @@
             f = self._forcing[0].toward(self._forcing[1], frac)
             sigma = mat.liquid.sigma if mat.liquid is not None else 0.0
             flow = skin_mod.Flow(self.columns, f.dp_ds, skin_mod.along_meridian(self.columns.normal, self.v_hat, f.tau),
-                                 f.g, self.v_hat, sk.rho_l, sigma, A_top, self.flow_cycles)
+                                 f.g, self.v_hat, sk.rho_l, sigma, A_top, self.flow_cycles, release=self.release,
+                                 slurry=self.slurry_spray)
         eps = np.full(sk.n, self.emissivity)
         m_min = skin_mod.minimum_mass(mat, sk.area, self.settings)
         side = lambda s, previous: skin_mod.Sideways(s, self.columns, previous)   # implicit (sub-plan 21)
@@ -368,6 +385,11 @@
             self.impulse_size += float(sum(np.linalg.norm(x, axis=1).sum() for x in (fb.push, fb.shear, fb.wall)))
             self.flow_courant, self.flow_growth = max(self.flow_courant, fb.courant), max(self.flow_growth, fb.growth)
             self.flow_passes, self.flow_cycles_max = max(self.flow_passes, fb.passes), max(self.flow_cycles_max, fb.cycles)
+            self.released_mass += float(fb.released.sum())   # the release from the top cells (Task 7): it leaves the body
+            self.released_enthalpy += float(fb.released_H.sum())
+            self.momentum_lost += fb.released_p.sum(axis=0)  # the droplets' momentum, at their cells' own velocities
+            self.droplet_momentum += fb.released_p
+            self.droplet_shear += fb.shear
         self.iterations.append(res.iterations)
         self.last = dataclasses.replace(res, Q_conv=E_top / dt, Q_rad=E_rad / dt)    # the gas-side powers, for the history
         self._move_base(book, m_min, g0)
@@ -466,10 +488,11 @@
         return self.solver.energy() + self.skins.energy() + float(self.pending.sum())
 
     def energy_balance_residual(self):
-        """(E - E0 - absorbed heat + removed enthalpy) / absorbed heat: zero for exact books."""
+        """(E - E0 - absorbed heat + removed and released enthalpy) / absorbed heat: zero for exact books."""
         if not self.absorbed_heat:
             return 0.0
-        return (self.energy() - self.energy0 - self.absorbed_heat + self.removed_enthalpy) / self.absorbed_heat
+        return ((self.energy() - self.energy0 - self.absorbed_heat + self.removed_enthalpy + self.released_enthalpy)
+                / self.absorbed_heat)
 
     def momentum_residual(self):
         """|P - P0 - impulses + momentum lost| over the impulses' total size: zero for exact books (spec 10.2)."""
@@ -515,6 +538,7 @@
         liquid = sk.depth_above(mat.T_liquidus) if mat.melts else np.zeros(sk.n)
         r = self.last_move
         i, j = self.columns.edges[:, 0], self.columns.edges[:, 1]
+        Pi, _, flagged = skin_mod.droplet_momentum_ratio(self.droplet_momentum, self.droplet_shear)
         return {"skin_cells": int(sk.count.sum()), "skin_thickness_max_mm": float(th.max() * 1e3),
                 "skin_above_min": int((sk.column_mass() > m_min * (1.0 + 1e-9)).sum()),
                 "flowing_depth_mean_mm": float(flowing.mean() * 1e3), "flowing_depth_max_mm": float(flowing.max() * 1e3),
@@ -542,7 +566,11 @@
                 "roughness_max": float((np.abs(th[i] - th[j]) / self.columns.edge_length).max()),
                 "flow_courant_max": self.flow_courant, "flow_growth_max": self.flow_growth,
                 "flow_passes_max": self.flow_passes, "flow_cycles_max": self.flow_cycles_max,
-                "wall_flow_s": self.wall["flow"]}
+                "wall_flow_s": self.wall["flow"],
+                # the droplet momentum sink (Task 7): the mass released from the top cells since time zero, fact 29's
+                # droplet-momentum ratio Pi over this macro step and the columns whose own Pi exceeds skin.PI_FLAG
+                "released_mass_kg": self.released_mass, "droplet_momentum_ratio": Pi,
+                "droplet_momentum_flagged": flagged}
 
     def _recession(self):
         """How far each corner's top has moved inward along its starting direction since time zero [m]."""
```

- [ ] **Step 6: The history's columns and the results** in `reentry_model/coupled.py`:

```diff
--- a/reentry_model/coupled.py
+++ b/reentry_model/coupled.py
@@ -47,7 +47,9 @@
                     "wall_skins_s", "wall_3d_s", "wall_move_s", "wall_rebuild_s", "step_parts",
                     # the flow (sub-plan 21 Task 6c, spec 10.2)
                     "momentum_residual", "speed_max_ms", "roughness_max", "flow_courant_max", "flow_growth_max",
-                    "flow_passes_max", "flow_cycles_max", "wall_flow_s"]
+                    "flow_passes_max", "flow_cycles_max", "wall_flow_s",
+                    # the droplet momentum sink (sub-plan 21 Task 7, spec 10.2): fact 29's Pi and its flag
+                    "released_mass_kg", "droplet_momentum_ratio", "droplet_momentum_flagged"]
 PVD_TEMPLATE = '<?xml version="1.0"?>\n<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">\n<Collection>\n{}</Collection>\n</VTKFile>\n'
 
 
@@ -184,7 +186,10 @@
                "surface_quality_min": float(np.min(last["surface_quality_min"])),
                "zone_ratio_min": float(np.min(last["zone_ratio_min"])),
                "momentum_residual": float(last["momentum_residual"][-1]),
-               "flow_courant_max": float(np.max(last["flow_courant_max"]))}
+               "flow_courant_max": float(np.max(last["flow_courant_max"])),
+               "released_mass_kg": float(last["released_mass_kg"][-1]),
+               "droplet_momentum_ratio_max": float(np.max(last["droplet_momentum_ratio"])),
+               "droplet_momentum_flagged_steps": int((last["droplet_momentum_flagged"] > 0).sum())}
         out.update({k: float(last[k][-1]) for k in ("removed_mass_kg", "removed_enthalpy_J", "drawn_mass_kg", "absorbed_mass_kg",
                                                   "handed_mass_kg", "rebuild_residual_J", "wall_skins_s", "wall_3d_s",
                                                   "wall_move_s", "wall_rebuild_s")})
```

- [ ] **Step 7: Run them to see them pass.** Expected: `tests/test_reentry_model_skin.py` 87 passed,
  `tests/test_reentry_model_skinbody.py` 18 passed (2 deselected, the reference tier's), `tests/test_reentry_model_cli.py`
  66 passed, `tests/test_reentry_model_coupled.py` 9 and `tests/test_reentry_model_columns.py` 13 as before (193, 803 s
  within the unit tier, the skin body's 595 s of it); the unit tier 498 passed, 1 skipped (1 040 s). Without a release
  nothing changes: the flow's results are 6d's to the last bit.

- [ ] **Step 8: Commit.**

```bash
git add reentry_model/skin.py reentry_model/skinbody.py reentry_model/coupled.py tests/test_reentry_model_skin.py tests/test_reentry_model_skinbody.py
git commit -m "Skins: the droplet momentum sink at the top cells, the placeholder for stage 4's spraying -- a release rate per column taken in every flow sub-step after the transport from the top cells, liquid first and the exposed slurry with the flag, never rigid metal, the metal leaving at its cells' own enthalpy and velocity so that every book closes; fact 29's droplet-momentum ratio and its flag in the history (sub-plan 21 Task 7)"
```

---

### Task 8: Verification — level 2's plate with the flow on, the flow's cost, and the share of the flow at short wavelengths

**Files:**
- Modify: `analysis/skin_plate.py` (`PUSH_ALPHA`, `PUSHES`, `RELEASES`, `FACE_TILT`, `REFERENCE_1D_LIQUIDUS`,
  `HOT_LIQUIDUS`; `references`, `inclined`, `face_columns`, `plate_forcing`, `column_record`, `save_state`;
  `run(..., flow, push, release, flow_cycles, record, save_at, save_path)`; `--flow`, `--push`, `--release`,
  `--flow-cycles`, `--columns`, `--save-state`; the history's flow columns), `analysis/ablation_1d.py` (`REMOVE_AT`,
  `run(..., remove_at)`, `--remove-at`)
- Create: `analysis/flow_wavelengths.py`
- Test: `tests/test_reentry_model_skinbody.py` (three tests), `tests/test_flow_wavelengths.py` (new, five tests)

**Interfaces:**
- Consumes: 6c's `SkinBody(..., flow, flow_cycles)`, `FlowForcing`, `advance(..., forcing)`, `momentum_residual` and the
  history's flow columns; 6d's chosen counts; Task 7's `prescribed_release` and `SkinBody(..., release)`; Task 5's
  `along_meridian`; sub-plan 20's plate (`run`, `projected_flux`, `solid_top`, its references); `Columns.edges`,
  `.conduct_weight`, `.area`.
- Produces: `skin_plate.run(..., flow=False, push=None, release=None, flow_cycles=None, record=False, save_at=None,
  save_path=None)` -> (body, rows, centre columns, failure[, records]), refusing a push, a release or a count with the
  flow off; `skin_plate.PUSHES` ({"none": None, "flight": (45.0, -3.424e5, 80.29, 29.77)}: deg, Pa/m, Pa, m/s2),
  `RELEASES` ({"none": None, "all": 1e9} kg/(m2 s)), `PUSH_ALPHA` (45), `FACE_TILT` (5 deg); `inclined(alpha_deg) ->
  (3,)`, `face_columns(body) -> (n,) bool`, `plate_forcing(body, face, push) -> FlowForcing`, `column_record(body) ->
  dict` (mass, thickness, top, T_top, flowing, liquid, top_speed, flux), `save_state(path, body, q)`, `references(case,
  at_liquidus=False)`, `REFERENCE_1D_LIQUIDUS`, `HOT_LIQUIDUS`; `ablation_1d.run(..., remove_at="feed")`, `REMOVE_AT`;
  `flow_wavelengths.wavelength_share(columns, field, depth, region=None, kh_short=KH_SHORT, basis=None) -> dict` (share,
  inertial, viscous, and per mode k, h, kh, power), `modes(columns, region=None) -> (k, phi, idx)`,
  `inertial_factor(kh)`, `viscous_factor(kh)`, `column_flux(sk)`, `kh_at(factor, target)`, `KH_SHORT` (0.56), `THICK`
  (1 mm).

Spec 11.1's level 2 ends: "With flow on and a push along the face: mass, momentum and energy close exactly and the
result is physically sensible (no exact reference exists)." This task runs sub-plan 20's plate with the flow on
(`analysis/skin_plate.py --flow on`), measures the flow's cost there, and measures the share of the flight's flow at
wavelengths where the depth-uniform pressure is not exact (spec 14 assumption 1; Asha, 2026-10-21, option a). Nothing in
`reentry_model/` changes. Measured 2026-10-24 in `prototype/work-2026-10-18-column-flow/measure/task8/`, each run under
`caffeinate -i` with its power source logged. Until 18:29 the laptop was on mains power and the plates ran three to five
side by side with the short-wave segments (load averages 5-9 on its eight cores: their wall times are under load); from
18:29 to about 20:10 it was on battery power, when the first coupled run of the pushed plate was stopped and the pushed
plate was measured on its skins alone, each run under 30 minutes and alone (item 2; the cost split, item 6); back on mains
power the coupled pushed plate ran alone to 1.0 s (22:11-23:48, detached from the session after a session's exit had
killed an earlier attempt 1 h 33 min in).

**1. The push along the face.** On the plate the flight direction is the face's normal, (-1, 0, 0), so
`along_meridian` gives no direction on the face: neither the wall pressure's gradient nor the gas shear can push along
it, and the deceleration pulls straight out of it. What the flight direction feeds in the skin body was checked: the
surface flow's angles (`theta`, read only by `flow_forcing`, which a given forcing bypasses), the stagnation and rear
patches of the surface statistics (not read by the plate), the shear's and the pressure term's directions in every
column's plane and the deceleration (`FlowForcing.g`); the plate's heating is its own projected flux. So the push tilts
the flight direction: the face becomes a patch of the 100 mm sphere's shoulder `PUSH_ALPHA` = 45 degrees from the
stagnation point (`inclined`: v_hat = (-cos 45, sin 45, 0), the stagnation point off the plate toward +y, the face's
meridian running to -y), and every column of the face (`face_columns`: corner at x = 0, normal within 5 degrees of the
face's; the edges' columns, leaning 45 degrees, are no part of a shoulder) takes the surface flow of that point at 93 s
of the 100 mm flight (`measure/task8/forcing45.py`): dp/ds = -3.42e5 Pa/m and 80.3 Pa of gas shear down the meridian,
and the deceleration of 29.8 m/s2 along the tilted flight direction -- 21.1 m/s2 out of the face (the front face's
Rayleigh-Taylor pull, the thickness term) and 21.1 m/s2 along it toward the nose, which with the liquid's 2 400 kg/m3
takes 5.1e4 Pa/m off the gas's push: a net 2.9e5 Pa/m down the face. Flight-like by construction: 45 degrees is where
6c's forcing has its largest pressure gradient (3.4e5 Pa/m between 44 and 46.5 degrees; the shear is largest, 87 Pa, at
49 degrees); at 115 s the same point has 5.5e5 Pa/m and 95 Pa. A forcing built directly -- g laid along the face, say --
would have exercised none of the meridian's geometry the flight uses.

**2. The plate pushed along its face.** The pushed level-2 plate (cold, 6 MW/m2, the skins' own density law) runs on a
plate 170 mm across instead of 110 (`--width-mm 170`; ruling): runoff leaving the heat freezes only where the next
conduction catches it, and the 110 mm plate has 10 mm of cold face beyond the heat -- on a reduced plate with 2 mm the
runoff reached the edge in the first step and ran over it (the unit test's plate). The coupled run, every pass logged
(`measure/task8/run_logged.py`; `runs/push_w170_to1`, 15 812 columns), ran to 1.0 s alone in 5 747 s (ruling 8: to 2.0 s
it would take about 6.5 hours, item 6), so it covers the first 0.3-0.4 s of melting. Its books close: mass to 4.4e-16,
energy to 6.5e-11 of the heat absorbed (two rebuilds), momentum to 1.2e-17 of the impulses. At 1.0 s, against the flow-off
plate of the same width at 1.0 s (`runs/off_w170_to1`), the face's columns downstream of the heat gained 0.22 g more and
those upstream 0.02 g (12 times less); the face's flux, 0.16 m3/s, points 2.2 degrees from -y; the liquid on 437 face
columns (33-83 um deep, median 54 um, viscous time 5.5 ms) has its top cell at 0.75 of the Couette-Poiseuille speed
above its base cell (10th-90th percentile 0.65-0.83); a deposit 60, 54 and 27 um thick on average 0-5, 5-10 and 10-15 mm
past the heat (within 20 mm of the centre line), none 15-35 mm past it, 136 of its 550 tops still above T_rigid; the
fastest cell 4.8 m/s, Courant numbers to 9.2 (the real passes at their trial passes' counts, item 6), growth at most
3.2e-4; the roughness 0.061, at the deposit; the grid scale's share of the face's thickness variance 7.4e-4; off the face
nothing moves and nothing is above T_rigid -- the columns there differ from the flow-off plate's only by their base
exchange (0.58 g, mostly at the box's far corners, 24 um on the face's own edge row 40 mm past the heat with no deposit
between). The centre recedes 0.795 mm (0.857 mm with the flow off; the push carries the liquid off). *One check fails*:
the history counts 1 761 base faults in the macro step from 0.5 to 1.0 s, none before and none on the flow-off and
regression plates: at 1.0 s 90 columns of the face around the heat's downstream rim and the deposit (39-61 mm from the
centre, y = -21 to -61 mm) have a bottom cell 5 to 400 K above its base node -- up to 1 253 K, which no heat reaches --
their masses from 6 % below to 5 % above their minimum. The
skins alone (below) have no cell above 910 K, so it is the coupling: the reconciliation at the base, which puts what the
move swept less what the skins exchanged into the bottom cell (spec 6.6), here into the small bottom cells of columns
whose deposits are being handed down while their neighbours draw -- the buried superheated cells of spec 12, which
absorption keeping up had removed, back where runoff deposits (assumed from the books, not traced cell by cell; Review
Focus 25). The pushed plate's physics over a longer push were measured on the plate's own skins
(`measure/task8/plate_skins.py`): the 170 mm plate with the flow off to 1.0 s (19 min alone; its centre
0.857 mm receded, 0.39 mm of slurry under the heat), its skins, columns, heating, base temperatures and the 3D gradient
kept (`--save-state 1.0`), then pushed for 0.2 s by `skin.advance`'s real pass in 10 ms parts with the base held at its
temperatures, the conduction on, the flow sub-steps chosen. Over those 0.2 s the books of the skins alone close to
round-off: mass to 3.5e-16 of the skins' (4.47 g drawn from below as the melted metal ran off), energy to 1.4e-15 of
the heat taken in, momentum to 4.4e-18 of the impulses' size. On the coupled body the books with the flow on are the
regression's (item 3: mass to 2.2e-16, energy 5.2e-13, momentum 1.3e-16 over 2 s), the reduced pushed plate's (the
unit test: 2.2e-16, 1.3e-12, 3.1e-17) and the pushed plate's above.

*Physically sensible* (`measure/task8/analyse_states.py`, `grid_roughness.py`; thresholds a ruling). After 0.2 s of
push:
- *the runoff goes the way the push points*: the face's columns downstream of the heat (y < 0) gained 4.37 g, those
  upstream 0.04 g (108 times less), and the face's volume flux, 9.3e-2 m3/s, points 0.8 degrees from -y (thresholds:
  ten times, 10 degrees);
- *the liquid's speed is the Couette-Poiseuille profile's*: on the 115 face columns holding at least three liquid
  cells (liquid 32-60 um deep, median 45 um, its viscous time h^2 / nu 3.7 ms, so near steady within a 10 ms heat
  sub-step), the top cell's speed over the liquid's base cell's, along the push, is 0.78 of tau h / mu + G h^2 / (2 mu)
  between their centres (10th-90th percentile 0.70-0.84; G = 2.93e5 Pa/m) -- the top at a median 3.4 m/s against 3.8.
  The 22 % short of steady is the transport's own loss of the kinetic energy it carries, 27-30 % at Courant numbers near
  6 (6c, 6d) (threshold: within 30 %);
- *metal leaving the heat freezes where it lands on cold metal*: the liquid reaches 2 mm past the heat's end (y = -47
  mm), and the deposit beyond it is a ridge 1.5 mm thick on average 0-5 mm past the heat, 0.7 mm at 5-10 mm, 0.14 mm at
  10-15 mm and 20 um at 15-20 mm, ending 18 mm past the heat (gain above 1 % of its largest), 108 of its 2 560 tops still
  above T_rigid (the bands within 20 mm of the centre line). Where the speed and the freezing time say: the deposit's
  thin leading edge, some 10-20 um of liquid on a rigid top, moves at about tau h / mu, 0.6-1.2 m/s, until the next
  conduction freezes it -- 6-12 mm in a 10 ms heat sub-step -- while such a film on 700 K metal would freeze in some
  0.03 ms (estimated, d^2 / alpha (L / c dT)^2). The heat sub-step sets it (item 4):
  halving it to 5 ms moves 11 % of the deposit 5-10 mm back toward the heat, and 2.5 ms changes nothing more;
- *no growth of grid-scale roughness*: the roughness history column rises from 0.018 to 0.63 at y = -49 mm -- the
  ridge's flank, a thickness of 2.5 mm over a few columns -- while the share of the face's thickness variance at the
  grid's shortest wavelengths (k above 0.8 of the largest, under about 5 mm) stays 0.14-0.16 % between 0.2 and 0.25 s
  as the whole field's rms grows from 0.39 to 0.50 mm: the grid scale grows with the ridge, not by itself (threshold:
  that share below 1 % and not growing against the rest);
- *nothing at the plate's edges*: the columns off the face gained 0.06 g in all, the largest 0.4 mg at a far corner of
  the box, 119 mm from the face's centre (a column there drawing toward its minimum), their fastest top 3.6e-7 m/s, none above
  T_rigid (threshold: no metal arriving, nothing moving).

**3. The regression against the exact references.** With no push and a release taking the liquid as it forms
(`--release all`: Task 7's `prescribed_release` at 1e9 kg/(m2 s), some 1e5 times what a column holds in a 1 ms flow
sub-step), the flow-on plate is the flow-off plate's physics with the flow running: surface tension still acts where a
top flows, and the release takes the liquid. It cannot take it as instant removal does: the release takes a cell at or
above the liquidus, 908 K (Task 7's ruling 3), instant removal at T_feed, 910 K. Stage 2's T_feed keeps a liquid layer
on the front across which the last 2 K fall, about 2 K k_l / q thick, and the melting starts that much later: the
one-dimensional reference with its front cells removed at the liquidus (`ablation_1d.py --remove-at liquidus`, 1 um,
0.5 ms) recedes further by 0.048, 0.060 and 0.079 mm at 1.0, 1.5 and 2.0 s on the cold plate (6.1, 2.2 and 1.5 %) and by
0.139, 0.149 and 0.158 mm on the hot plate (17.6, 8.2 and 5.4 %), its non-rigid depth 0.13-0.14 mm (9 %) thinner --
beyond level 2's 5 % on the hot plate at 1.5 and 2.0 s. So the regression is judged against the references with the
liquid leaving at the liquidus (`references(case, at_liquidus=True)`): the one-dimensional reference so computed, and
the fine box moved by the one-dimensional difference (it has no liquidus run of its own). Both plates at the flight's
0.5 s macro step (where both pass with the flow off), the skins at the solid's density, development cells, the fine
band, to 2.0 s:

| Plate | 1.0 s (hot: onset, not judged) | 1.5 s | 2.0 s | Hot: flowing depth, base above T_rigid |
|---|---|---|---|---|
| cold, flow off, against T_feed's: 1D; fine box | +3.2 %; +1.2 % | +0.5 %; -0.2 % | -0.9 %; -0.6 % | |
| cold, flow on and released, against the liquidus's: 1D; fine box | +2.4 %; +0.6 % | +0.3 %; -0.3 % | -0.8 %; -0.6 % | |
| hot, flow off, against T_feed's | -9.0 % | -3.0 % | -1.1 % | within 1.1 %; at most 0.09 K |
| hot, flow on and released, against the liquidus's | -6.4 % | -1.6 % | -1.1 % | within 1.2 %; at most 0.10 K |

The table shows the flow-on plates meeting level 2 against the liquidus's references as closely as the flow-off
plates meet it against stage 2's: every judged recession within 1.6 % (5 % allowed), the hot plate's flowing depth
within 1.2 % of the exact non-rigid depth (10 %) and its base within 0.10 K of T_rigid (1 K); against stage 2's
references the flow-on plates are 6.5 % ahead at 1.5 s (hot) and 8.7 % at 1.0 s (cold), the threshold's effect, not the
flow's. The flow-off plates reproduce Task 1's level-2 numbers (fact 150) to the printed digit. Books with the flow on:
mass to 2.2e-16, energy to 5.2e-13 (cold) and 2.5e-13 (hot), momentum to 1.4e-17 and 1.3e-16; no base fault; 49.5 g
(cold) and 32.7 g (hot) released; the fastest cell 9 mm/s (cold) and 4.8 cm/s (hot, the slurry at the melt's rim,
where the crater's curved top gives surface tension some 5 000 Pa/m, estimated). The momentum residual and the
roughness are the history's: the roughness of the hot plate's rim, where the 1.5-2.0 mm skins end over the heat's
taper, is 0.25, 0.38 and 0.45 at 1.0, 1.5 and 2.0 s with the flow on against 0.35, 0.47 and 0.60 with it off (measured
from the column records over the grid's sides and diagonals, `measure/task8/roughness.py`): surface tension smooths the
rim; on the cold plate both are the box corners' 0.007-0.013.

**4. The flow sub-step and the heat sub-step** (6c's and 6d's study on the plate, `plate_skins.py substep`: 50 ms more
from the pushed skins at 0.2 s, the conduction on, against the 0.25 ms run):

| Flow sub-step | Largest Courant | Where the metal went, L1 / the finest run's change | Of it beyond the heat | KE at 50 ms / finest | Processor time |
|---|---|---|---|---|---|
| 0.25 ms (40) | 0.71 | -- | -- | 1 | 134 s |
| 0.5 ms (20) | 1.43 | 0.112 | 0.087 | 0.997 | 77 s |
| 1 ms (10) | 2.90 | 0.196 | 0.153 | 1.028 | 51 s |
| chosen: 2 ms (5) | 5.91 | 0.254 | 0.171 | 1.090 | 35 s |

The table shows the plate's runoff at the count the skins choose (Courant 5.9, at the target of 6; the growth at most
3.7e-4 and one linearisation, so the Courant number sets the count) 25 % off the finest run's places in L1 and its
kinetic energy 9 % high, against 2.5 % and 5 % on 6d's sphere: two thirds of the difference
is in the deposit beyond the heat, where coasting landed metal runs ahead by the drag it would have felt in one flow
sub-step (6d's ruling; Review Focus 13): at the chosen count 13.7 mg of the 338 mg deposited in 50 ms within 20 mm of
the centre line lies 10-35 mm past the heat, against 0.5 mg in the finest run. The metal drawn from below as the melt
ran off agrees within 0.6 % (1.054 against 1.060 g). *The heat sub-step*, at 0.5 ms flow sub-steps: of about 345 mg
deposited within 20 mm of the centre line in 50 ms, 39 mg lies 5-10 mm past the heat with 10 ms heat sub-steps, 2.1 mg
with 5 ms and 1.8 mg with 2.5 ms, and 81, 55 and 52 deposit tops are still above T_rigid at the end: the flow sub-steps conduct
no heat, so liquid landing on cold metal stays liquid until the next conduction; at the 10 ms heat sub-step 11 % of
the deposit lands 5-10 mm further than it does once the heat sub-step is 5 ms or less (Review Focus 21).

**5. The swelling setting** ("Inputs to settle"). Ruling: the regression keeps the skins at the solid's density
(`--swell off`) -- it is judged against references that hold one density (spec 11.1, Asha, 2026-10-14) -- and the pushed
plate runs with the skins' own law (`--swell on`), the flight's. Measured on the pushed skins over 50 ms at the chosen
count, the same state read at both densities: at the solid's density the liquid is 45 um deep where it is 58 um
swollen, the face's volume flux 0.63 of the swollen layer's (Couette's h^2 and Poiseuille's h^3 at 0.85 of the depth
give 0.72 and 0.61), its kinetic energy 0.56, its fastest top 3.0 m/s against 3.9, and the metal deposited beyond the
heat in 50 ms 80 mg against 338 mg -- if wrong, a push case at one density moves a third less metal.

**6. The cost of the flow on the plate.** *Per skin sub-step, alone* (`plate_skins.py cost`, the pushed skins at 0.2 s,
16 500 columns, two repetitions each, processor time):

| One 10 ms skin sub-step | Whole | Conduction (with the implicit sideways solve) | Push | Momentum | Thickness | Transport | Release | Regrid | Rest |
|---|---|---|---|---|---|---|---|---|---|
| flow off (removal) | 3.48-3.55 s | 3.27-3.35 | -- | -- | -- | -- | -- | -- | 0.20 |
| flow on, real pass, 5 flow sub-steps | 6.66-6.77 s | 3.30-3.33 | 0.35-0.37 | 0.52 | 0.77 | 1.09-1.11 | 0 | 0.40-0.43 | 0.24-0.25 |
| flow on, trial pass | 6.46-6.55 s | 3.17-3.20 | 0.32-0.33 | 0.51 | 0.75-0.76 | 1.10-1.13 | 0 | 0.35-0.39 | 0.23-0.27 |

The table shows the flow doubling a skin sub-step on the pushed plate (1.9 times), the transport the largest part of
it, a third; the trial pass costs what the real pass does (6d). The push case releases nothing; in the cold
regression's last macro step the release took 30 s of the flow's 277 s (Task 7: 0.054 s a call on the sphere).

*Per macro step*, from the runs (wall time under load, so indicative; the flow's split from `skin.FLOW_SECONDS`):

| 0.5 s macro step | cold, flow off | cold, flow on and released | hot, flow off | hot, flow on and released |
|---|---|---|---|---|
| 0.5-1.0 s | 547 s | 741 s (flow 177) | 216 s | 292 s (flow 120) |
| 1.0-1.5 s | 625 s | 1 008 s (flow 262) | 204 s | 262 s (flow 97) |
| 1.5-2.0 s | 796 s | 1 019 s (flow 277) | 187 s | 349 s (flow 135) |

The table shows the regression's flow -- slurry at millimetres to centimetres a second, one flow sub-step a skin
sub-step after the first -- adding a quarter to two thirds to a melting macro step; of the flow's 277 s in the cold
plate's last step, transport 106, thickness 64, momentum 52, release 30, push 26. *The pushed plate* (alone, every pass
logged, `measure/task8/push_cost.py`): its first macro step 518 s, its second, where melting starts, 5 229 s -- skins
4 172 s (of it flow 1 888: transport 706, thickness 495, momentum 287, push 231, regrid 170), 3D solves 785 s, moves 261 s
-- in 312 passes, against 547 s for the flow-off 110 mm plate's second step under load. Three reasons: (a) the parts are
rejected for their moves -- the runoff's deposits are handed down to the 3D model and move bases, so the 0.0625 s parts
the melting cap makes are rejected after their damped passes as long parts were before the cap (fact 126): the first
three parts took 4 passes each, the fourth 9 and was run again as two; from the fifth every part was rejected after 5-7
passes, run again as three parts of 0.0208 s, each of those rejected after 3 passes and run again as two of 0.0104 s --
880-950 s of processor time for each 0.0625 s of melting against about 90 s for a part's trial and real pass that is not
rejected; (b) a real pass whose trial pass, its base held at the start-of-step temperature, melted less runs
at the trial's counts: one flow sub-step a skin sub-step at a Courant number of 18.8 on the third part, three times the
target (6d's ruling; Review Focus 17); (c) in the first macro step, before any metal flowed in an accepted pass, the
first real pass held a few cells above T_rigid at the 3D model's first answer for the base and ran 380 flow sub-steps
in which nothing moved (Courant 0), 96 s of its 124 s: a flow sub-step costs about 0.25 s on this plate however few
cells flow, and the real pass replayed the ten its trial pass chose with no memory. 13 of the run's 156 real passes ran
above the Courant target, the largest at 18.8. Estimated alone from these parts: a melting macro step of the pushed plate
about 8 x 1 000 s, 2.2-2.5 hours, and the plate to 2.0 s about 6.5 hours (ruling 8).

**7. The share of the flight's flow at short wavelengths** (spec 14 assumption 1; Asha, 2026-10-21, option a). The
columns push every layer of a column alike, so a disturbance whose wavelength is not long against the flowing depth,
kh not small, responds too fast: by sqrt(kh / tanh kh) in an inertial layer (spec 14.1), and in a viscous one by the
lubrication rate h^3 k^2 (rho g + sigma k^2) / (3 mu) -- what columns with a resolved viscous profile and a depth-uniform
pressure give -- over the Stokes rate of a viscous layer on a rigid base, (rho g + sigma k^2)(sinh 2kh - 2kh) / (2 mu k
(cosh 2kh + 2 (kh)^2 + 1)): 1.04 and 1.45 at kh = 0.5, 1.15 and 2.77 at kh = 1. Which applies is set by the layer's
viscous time h^2 / nu: 2.7 s for 1.2 mm of liquid (inertial against any disturbance the flight drives), 0.01-0.04 s
for 1-3.4 mm of slurry at 1-5 Pa s (viscous). The element model's thick layers late in flight are slurry under a thin
film (its mean film 0.05 mm at 105 and 120 s, `slurry_thick_fraction` 0.41 and 0.54), so the viscous factor is the one
that applies there, and it is not the spec's.

*The measure* (`analysis/flow_wavelengths.py`; ruling): the share of the flowing layers' volume flux q -- each column's
sum of its cells' thicknesses times their velocities [m2/s] -- carried at kh above KH_SHORT = 0.56 (where the inertial
factor reaches 1.05; the viscous one does at 0.17), over the columns whose flowing depth is at least 1 mm, and the
flux-weighted mean excess of either factor, sum P_n (F(k_n h_n) - 1) / sum P_n. Wavelengths are the modes of the
columns' cotangent Laplacian over that region (free at its border, the columns' areas the mass), each at its own mean
depth h_n = sum_j A_j phi_nj^2 h_j; a field's power in mode n is (phi_n^T A f)^2, summed over its Cartesian components.
On a flat grid of squares a sampled cosine is one mode exactly, of wavenumber (2 / dx) sin(kappa dx / 2) (the test's
cylinder). Why the flux: the columns' error is in how fast they move metal in response to disturbances along the
surface, and the flux is what moves it -- where the runoff takes it and how deep the flowing layer is that Girin's
spraying reads (spec 6.7, 8); a share s means that at most s of where the metal goes is subject to the overstated
rates, and the mean excess says by how much. The top's height is measured the same way (its variance), for the layer's
shape.

*The states* (`measure/task8/short_waves.py`). A full flow-on flight is out of reach (item 6), so the 100 mm sphere's
9 417 columns of production cells are given the flowing depth the element model's 0.0125 s flight had on its wet
windward patches -- a mean of 1.225 mm at 105 s and 3.384 mm at 120 s -- in two shapes over the angle from the
stagnation point (uniform tapering as cos^2 to nothing between 70 and 90 degrees, `flat`, or h0 cos theta, `cos`),
scaled to that mean; as `slurry`, the element model's mean film (0.056 and 0.045 mm) liquid at T_feed over slurry
falling linearly from the liquidus to T_rigid, or as `liquid`, all of it at 950 K (the inertial case); at rest. A segment
of 0.1 s of flow alone under the flight's forcing at that time (the surface flow's dp/ds and shear, the deceleration of
39.3 and 53.2 m/s2, surface tension; the flow sub-steps chosen) moves it, and the measure is read after 10, 50 and
100 ms (processor time 216-372 s a segment, 749 s for the liquid):

| State | Columns over 1 mm at 100 ms (mean depth) | Flux at kh > 0.56 at 10, 50, 100 ms | Mean excess at 100 ms, inertial; viscous | Flux at kh > 0.17 at 100 ms |
|---|---|---|---|---|
| 105 s, flat, slurry | 2 300 (1.44 mm) | 0.2 %, 0.1 %, 0.1 % | 0.1 %; 1.4 % | 5.6 % |
| 105 s, cos, slurry | 3 252 (1.36 mm) | 1.0 %, 1.2 %, 0.4 % | 0.2 %; 1.7 % | 4.4 % |
| 120 s, flat, slurry | 4 108 (3.73 mm) | 5.3 %, 26 %, 26 % | 7.8 %; 201 % | 75 % |
| 120 s, cos, slurry | 4 424 (3.37 mm) | 6.6 %, 8.9 %, 36 % | 5.0 %; 94 % | 71 % |
| 120 s, flat, liquid | 6 649 (2.48 mm) | 4.3 %, 2.9 %, 5.1 % | 0.9 %; 15 % | 14 % |

The table shows the share negligible at the element model's depth at 105 s and not at 120 s. At 1.2 mm at most 1.2 %
of the flux is at kh > 0.56 within 0.1 s, and the columns' rates are on average 0.2 % too fast if the layer is inertial
and 1.7 % if viscous. At 3.4 mm a quarter to a third of the slurry's flux is at kh > 0.56 after 50-100 ms, and on
average the columns' rates are 5-8 % too fast if it were inertial and 94-201 % if viscous -- the slurry is. Where it is
(`measure/task8/sw_where.py`, the 100 ms states): the slurry drains off the stagnation region at up to 5 m/s in the
first 10 ms (its liquid film first) and piles at the shoulder, 8.6 mm deep on average between 75 and 90 degrees (`flat`)
and 7.7 mm (`cos`), the stagnation region left at 1.2-1.5 mm; the short-wave flux lies on that pile's fronts, between
45 and 75 degrees and past 90, at wavelengths of 18-36 mm (kh 0.6-1.5), smooth from column to column (the short part's
correlation between neighbours 0.85-0.90: a front, not a column-scale ripple). Where the 3.4 mm layer is liquid it runs
off faster and spreads (6 649 columns over 1 mm at 100 ms) and 5 % of its flux is short. Run on to 0.4 s (the `cos`
slurry state, back on mains power), the slurry settles into the pile, 9.6 mm deep on average between 75 and 90 degrees,
its fastest cell slowing from 0.76 to 0.16 m/s and the flux's rms falling six times (from 4.5e-4 to 7e-5 m2/s): the share
at kh > 0.56 rises to 57, 70 and 78 % at 0.2, 0.3 and 0.4 s, the inertial excess to 51-75 % and the viscous one to factors
of 45-87, carried by the deep pile's creep at kh of 3-5, now less smooth from column to column (correlation 0.23). A
rerun of the `flat` slurry segment gave every number to the last bit.

*Uncertainty and what the states leave out.* Each state was measured once; the two shapes of the same mean depth give
26 and 36 % at 100 ms and 26 and 9 % at 50 ms, and the share grows with time without settling: at 3.4 mm it is "of order
a quarter, between 5 and 36 % over the first 0.1 s", rising to three quarters by 0.4 s of a flux six times smaller as
the slurry comes to rest in a pile; at 1.2 mm "about 1 % or less". The states leave out the melting
that feeds the layer and the heating (flow alone, no conduction: the slurry crossing the shoulder onto the rear does not
freeze, so the pile builds where in flight it would partly freeze); spraying (Girin's release would take the film and,
with the slurry flag, slurry); the flow's history before 105 and 120 s (the states start at rest from the element
model's mean depth; its distribution over the surface is not in its outputs and is assumed here); and the receding
body's shape. Above all the depths are the element model's, whose runoff is not the columns': the columns drained the
3.4 mm slurry at metres a second within tens of milliseconds, so a flow-on flight may not hold layers that deep where the
push is strong, and may hold deeper piles where it is weak. *What it would change*: a non-hydrostatic pressure (spec
14.1's remedy; Asha, 2026-10-22) would slow exactly this short-wave part -- the advance of steep fronts and piles, and
Rayleigh-Taylor growth at kh near 1 -- by the factors above at each wavelength, and leave the long-wave runoff, three
quarters of the flux at 120 s and 99 % at 105 s, unchanged; on the 3.4 mm slurry states the flux-weighted overstatement
is a factor of about 2-3, on the 1.2 mm ones about 1.02. **Decided (Asha, 2026-10-25, option a, over building the
non-hydrostatic pressure now):** the share is measured again on the flow-on flight's own states in sub-plan 22 before
deciding whether it is significant (Review Focus 23).

**Rulings** (2026-10-24, within spec 6.4-6.7, 7.1, 11.1 and 14 and the rulings of Tasks 1-7; each is what, then why,
then the cost if wrong):
1. *The push: the face as a patch of the 100 mm sphere's shoulder at 45 degrees at 93 s, the flight direction tilted,
   the forcing on the face's columns only.* Why: item 1; it drives every term of the push through the code the flight
   uses, at the flight's largest pressure gradient of that moment; the edges' columns lean 45 degrees and belong to no
   shoulder. Cost if wrong: a push not met in flight is not exercised -- the late flight's are 1.6 times larger -- and
   runoff reaching an edge piles where the gas stops pushing.
2. *The pushed plate 170 mm across.* Why: item 2, so that the runoff's deposit is not cut off by the edge. Cost if wrong:
   2.2 times the columns of the flow-off part of a skin sub-step (the conduction's 3.3 s of the 6.7 s).
3. *The regression judged against references with the liquid leaving at the liquidus.* Why: item 3; Task 7's release
   takes liquid at the liquidus, and the references' T_feed moves the hot plate's exact recession by 5-18 %. Cost if
   wrong: none for the flow; the comparison then measures the threshold.
4. *The release that takes the liquid as it forms is 1e9 kg/(m2 s).* Why: some 1e5 times what a column holds in a 1 ms
   flow sub-step, so `release_top` takes every liquid cell. Cost if wrong: nothing, the take being capped.
5. *The swelling: the regression at the solid's density, the pushed plate at the skins' own law.* Why: item 5. Cost if
   wrong: a third of the runoff.
6. *The sensibility checks and thresholds* as listed in item 2. Why: each is what an exact solution would fix if one
   existed -- the direction, the steady profile, the freezing, a mode the grid cannot hold, the boundary. Cost if wrong:
   a check passed loosely; every measured value is in item 2.
7. *The short-wave measure: the flux's share at kh > 0.56 over flowing depths of at least 1 mm, the surface Laplacian's
   modes, both factors' mean excess.* Why: item 7. Cost if wrong: the height's share, the other candidate, is 0.20-0.65
   at 120 s and 0.05-0.38 at 105 s, but much of it is the states' own shape (0.42 before any flow in the 120 s `flat`
   state); the kinetic energy's share weighs the fast thin film, which drained within 10 ms.
8. *The coupled pushed plate to 1.0 s, its longer push on its skins alone.* Why: from the passes of the run a
   session's exit killed (its first seven melting parts, 5 030 s of wall time alone; 880-950 s of processor time for each
   0.0625 s once deposits split the parts), a melting macro step costs 2.2-2.5 hours alone and the plate to 2.0 s about
   6.5 hours, over the six the controller set; to 1.0 s it took 1 h 36 min. The skins carry the flow's physics, the
   coupled run its books and the coupling's response. Cost if wrong: the coupled plate's last second -- the deposit
   growing, the base faults' course -- is not measured; the skins alone cover 0.2 s of push after 1.0 s.

**What this task found in earlier tasks.**
- *Task 7's threshold* (ruling 3, "cells within those 2 K wait for the slurry flag"): releasing at the liquidus instead
  of T_feed takes away the liquid layer stage 2's removal keeps, about 2 K k_l / q thick, and the plates' exact recession
  moves by 1.5-6.1 % (cold) and 5.4-17.6 % (hot) at 1.0-2.0 s, the non-rigid depth by 9 % (item 3). Not changed: the
  liquidus is the material's own (Task 7's reason); the references with it are now in `skin_plate.py`.
- *6d's counts and the coupling*: a real pass runs at the counts its trial pass chose, which with its base held at the
  start-of-step temperature can melt less and flow slower, so at the onset the real pass ran at a Courant number of 18.8,
  three times the target; and a real pass whose cells flow only because of its base runs the trial's no-memory ten
  flow sub-steps at 0.25 s each on the 170 mm plate however little flows (item 6).
- *The parts and the deposits*: runoff landing beyond the heat makes skins thicker than their minimum, the bottoms are
  handed down to the 3D model and bases move, so melting parts of 0.0625 s are rejected for their moves after their
  damped passes and run again as parts of 0.0104 s (item 6): the cost the part cap removed (fact 126) returns where runoff
  deposits, ten times a part's cost.
- *The base's reconciliation and the deposits*: in the same columns 1 761 base faults and bottom cells up to 1 253 K
  in one macro step (item 2), the buried superheated cells of spec 12 back; the books stay exact.
- *6d's coasting*: on the plate two thirds of the chosen count's difference from the finest run is the deposit's front
  running ahead (item 4), as 6d's ruling said it would.

**Measured against the tests** (2026-10-24): the flow-off plate through Task 8's `run` is the body built and stepped
by hand as sub-plan 20's script did, to the last bit, and sub-plan 20's script on the planning copy's package gives
Task 8's every array and history column to the last bit over three steps (`measure/task8/bitid_flowoff.py`); the
reduced pushed plate's books to 2.2e-16 (mass), 1.3e-12 (energy) and 3.1-4.7e-17 (momentum), 464 times more metal
gained downstream than upstream, the face's flux 5.0 degrees from -y; with nothing flowing the flow-on plate is the
flow-off plate to the last bit, and with slurry tops and no push it differs by 6.8e-7 of the largest thickness and
6.3e-4 K after 0.2 s, at 1.6e-5 m/s (`measure/task8/measure_tests.py`); the measure on the cylinder's sinusoids to
1e-15 and 3e-15 of its definition.

- [ ] **Step 1: The plate's tests** -- append to `tests/test_reentry_model_skinbody.py`:

```diff
--- a/tests/test_reentry_model_skinbody.py
+++ b/tests/test_reentry_model_skinbody.py
@@ -528,3 +528,95 @@
     assert stats["droplet_momentum_ratio"] == Pi > 0.0 and stats["droplet_momentum_flagged"] == flagged
     off = make(coarse_sphere_mesh, flow=True, remove=False, flow_cycles=1)
     assert off.skin_stats()["droplet_momentum_ratio"] == 0.0 and off.skin_stats()["released_mass_kg"] == 0.0
+
+
+# -- level 2 with the flow on: the plate (sub-plan 21 Task 8) -----------------------------------------------------------
+def skin_plate_module():
+    import os
+    import sys
+    here = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "analysis")
+    if here not in sys.path:
+        sys.path.insert(0, here)
+    import skin_plate
+    return skin_plate
+
+
+# a reduced plate for the unit tier: 24 mm across, heated within 7 mm of its centre, the flux gone at 10 mm; 2 mm
+# elements and no fine band (box_mesh is structured: no meshing), the skins' development cells
+SMALL_PLATE = dict(verbose=False, width=24e-3, r_flat=7e-3, r_zero=10e-3, fine=None)
+
+
+def test_the_plate_with_the_flow_off_is_sub_plan_20_s_plate():
+    """`analysis/skin_plate.py` with the flow off builds and steps the plate as sub-plan 20's script did, so its level-2
+    results stand: a reduced plate from 880 K under 6 MW/m2 for two 0.1 s steps, which melts, removes liquid and hands
+    metal down, is the body built and advanced by hand to the last bit (numpy's generator seeded: pyamg draws from it).
+    Sub-plan 20's script itself, run on the planning copy's package, gives Task 8's every array and history column to the
+    last bit on such a plate over three steps (measured 2026-10-24, measure/task8/bitid_flowoff.py). A push, a release or
+    a flow sub-step count with the flow off is refused."""
+    sp_ = skin_plate_module()
+    np.random.seed(12345)
+    body, rows, _, failure = sp_.run(dt=0.1, t_end=0.2, T_start=880.0, q0=6e6, **SMALL_PLATE)
+    np.random.seed(12345)
+    mat = material.Material.from_drama_json(material.MATERIAL_NAMES["AA7075_scheil"])
+    m = mesh.box_mesh(sp_.LENGTH, 24e-3, 24e-3, 2e-3)
+    m.points[:, 1:] -= 12e-3
+    hand = skinbody.SkinBody(m, mat, thermal_solver("skfem"), mat.rho * m.volume(),
+                             skin.SkinSettings(cell=skin.SKIN_PRESETS["dev"], swell=False), T0=880.0, emissivity=0.0,
+                             v_hat=(-1.0, 0.0, 0.0), deform_depth=4e-3, repeat_rule=skinbody.DEFAULT_REPEAT_RULE)
+    p = hand._top0
+    radius = np.where(np.abs(p[:, 0]) < 0.05e-3, np.hypot(p[:, 1], p[:, 2]), np.inf)
+    for k in range(2):
+        q = sp_.projected_flux(hand, radius, 7e-3, 10e-3, 6e6)
+        hand.advance(0.1 * (k + 1), 0.1, heating.HeatingResult(q, 6e6, 0.0, 0.0, 0.0))
+    assert failure is None and body.removed_mass > 0.0 and same_body(body, hand)
+    for kw in ({"push": sp_.PUSHES["flight"]}, {"release": 1.0}, {"flow_cycles": 2}):
+        with pytest.raises(ValueError):
+            sp_.run(t_end=0.0, **SMALL_PLATE, **kw)
+
+
+def test_a_plate_pushed_along_its_face_closes_its_books_and_its_metal_runs_the_way_the_push_points():
+    """Spec 11.1's level 2 with the flow on: the reduced plate from 880 K under 6 MW/m2 for two 0.1 s steps, its face a
+    patch of the 100 mm sphere's shoulder 45 degrees from the stagnation point at 93 s (`skin_plate.PUSHES['flight']`:
+    dp/ds -3.42e5 Pa/m and 80 Pa of shear down the face's meridian, -y, and 29.8 m/s2 of deceleration along the tilted
+    flight direction), the flow sub-steps chosen, the liquid kept. The books close to round-off -- the mass to 1e-13, the
+    energy to 1e-8 of the heat absorbed, the momentum to 1e-12 of the impulses -- and the metal runs down the face: on the
+    face the columns downstream of the heat (y < 0) gain more than a hundred times what those upstream do, and the face's
+    volume flux points within 10 degrees of -y. Measured 2026-10-24 (measure/task8/measure_tests.py): mass 2.2e-16,
+    energy 1.3e-12, momentum 3.1e-17 (4.7e-17 in a second run: pyamg's starting vectors are not seeded here); gained
+    downstream 3.3e-5 kg, upstream 7.2e-8 kg (464 times); the face's flux 5.0 degrees from -y (on the plate's crater
+    wall the face's columns tilt), the fastest cell 2.4 m/s, 10 and 4 flow sub-steps a heat sub-step at most. On this
+    plate the runoff reaches the edge 2 mm beyond the heat within the first step and runs over it: with the side's
+    columns the whole flux is 19 degrees from -y."""
+    sp_ = skin_plate_module()
+    body, rows, _, failure, records = sp_.run(dt=0.1, t_end=0.2, T_start=880.0, q0=6e6, flow=True,
+                                              push=sp_.PUSHES["flight"], record=True, **SMALL_PLATE)
+    assert failure is None and body.flow_on and body.removed_mass == 0.0 and body.released_mass == 0.0
+    assert max(abs(x["mass_error"]) for x in rows) < 1e-13
+    assert max(abs(x["energy_residual"]) for x in rows) < 1e-8
+    assert max(x["momentum_residual"] for x in rows) < 1e-12 and body.impulse_size > 0.0
+    face = sp_.face_columns(body)
+    y = body._top0[:, 1]
+    gained = np.maximum(records[-1]["mass"] - records[0]["mass"], 0.0)
+    assert gained[face & (y < 0.0)].sum() > 100.0 * gained[face & (y > 0.0)].sum() > 0.0
+    flux = records[-1]["flux"][face].sum(axis=0)
+    assert -flux[1] / np.linalg.norm(flux) > np.cos(np.radians(10.0))
+    assert rows[-1]["speed_max_ms"] > 0.0 and rows[-1]["flow_cycles_max"] >= 1
+
+
+def test_a_plate_with_the_flow_on_and_nothing_flowing_is_the_flow_off_plate():
+    """With no push and no metal above T_rigid -- the reduced plate from 700 K under 1 MW/m2 for two 0.1 s steps -- the
+    flow-on plate is the flow-off plate to the last bit (numpy's generator seeded): no cell flows, so no flow sub-step
+    is taken, nothing moves and nothing is released. Once the top turns to slurry the two part, even with no push and
+    nothing liquid: surface tension acts where a top flows (6b's ruling) and moves the slurry wherever the tops' heights
+    differ -- on the plate where the draws and hand-downs have left neighbouring columns a little apart -- so that
+    starting at 880 K under 2 MW/m2, 45 tops of slurry and none liquid, the flow-on plate's skins differ from the
+    flow-off plate's by 6.8e-7 of their largest thickness after 0.2 s and its 3D field by 6.3e-4 K, at speeds of at most
+    1.6e-5 m/s (measured 2026-10-24, measure/task8/measure_tests.py)."""
+    sp_ = skin_plate_module()
+    bodies = []
+    for flow in (False, True):
+        np.random.seed(12345)
+        bodies.append(sp_.run(dt=0.1, t_end=0.2, T_start=700.0, q0=1e6, flow=flow, **SMALL_PLATE)[0])
+    off, on = bodies
+    assert on.flow_on and same_body(off, on) and np.all(on.skins.T[on.skins.valid] < on.material.T_rigid)
+    assert on.skin_stats()["speed_max_ms"] == 0.0 and on.momentum_residual() == 0.0
```

- [ ] **Step 2: The short-wave measure's tests** -- create `tests/test_flow_wavelengths.py`:

```python
"""The share of the columns' flow at short wavelengths (sub-plan 21 Task 8, `analysis/flow_wavelengths.py`): the
measure of the depth-uniform pressure's limitation (layered-skin design section 14, assumption 1)."""
import os
import sys

import numpy as np
import pytest

from reentry_model import columns

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "analysis"))
import flow_wavelengths as fw  # noqa: E402


def cylinder(N=40, M=4, dx=1e-3):
    """An open cylinder of N x M square facets of side dx, each split in two: flat to itself and periodic around, so
    the cotangent Laplacian is the five-point one (weights 1 on the squares' sides, 0 on their diagonals, half on the
    open ends' edges, whose columns have half the area) and a sampled cosine around it is one of its modes exactly."""
    R = dx / (2.0 * np.sin(np.pi / N))
    a = 2.0 * np.pi * np.arange(N) / N
    points = np.array([(R * np.cos(t), R * np.sin(t), r * dx) for r in range(M + 1) for t in a])
    at = lambda r, j: r * N + j % N
    faces = [f for r in range(M) for j in range(N)
             for f in ((at(r, j), at(r, j + 1), at(r + 1, j + 1)), (at(r, j), at(r + 1, j + 1), at(r + 1, j)))]
    return columns.Columns(points, np.array(faces))


def test_the_overstatement_factors_are_the_inertial_and_viscous_layers_rate_ratios():
    """sqrt(kh / tanh kh) (spec 14.1) and the lubrication rate over the Stokes rate of a viscous layer on a rigid base,
    by hand at kh = 1: sqrt(1 / tanh 1) = 1.14588 and (2/3)(cosh 2 + 3) / (sinh 2 - 2) = 2.77112; both 1 at kh = 0, the
    viscous one's series 1 + 9/5 (kh)^2 meeting the formula at kh = 0.05 to 3e-7; 1.05 at kh = 0.560 (KH_SHORT) and
    0.167."""
    assert fw.inertial_factor(1.0) == pytest.approx(np.sqrt(1.0 / np.tanh(1.0)), rel=1e-15)
    assert fw.viscous_factor(1.0) == pytest.approx((2.0 / 3.0) * (np.cosh(2.0) + 3.0) / (np.sinh(2.0) - 2.0), rel=1e-14)
    assert fw.inertial_factor(0.0) == 1.0 and fw.viscous_factor(0.0) == 1.0
    assert fw.viscous_factor(0.05) == pytest.approx(fw.viscous_factor(0.0500001), abs=1e-6)
    assert fw.inertial_factor(fw.KH_SHORT) == pytest.approx(1.05, abs=1e-3)
    assert fw.kh_at(fw.viscous_factor, 1.05) == pytest.approx(0.167, abs=1e-3)


@pytest.mark.parametrize("m", [1, 4, 8])
def test_a_single_sinusoid_on_a_flat_layer_has_the_share_and_excess_of_its_own_wavenumber(m):
    """cos(2 pi m j / N) around the cylinder of 40 facets of 1 mm on a layer 1.2 mm deep: all its power is in its own
    mode, of the five-point Laplacian's wavenumber k = (2 / dx) sin(pi m / N) -- 156.9, 618.0 and 1 175.6 1/m, kh 0.19,
    0.74 and 1.41 -- so its short-wave share is 0 below KH_SHORT and 1 above, and its mean excesses are the factors'
    at that kh less one, to 1e-15 and 3e-15 (measured 2026-10-24); its power is its area-weighted square, half the
    cylinder's 160 mm2."""
    c = cylinder()
    j = np.arange(c.n) % 40
    res = fw.wavelength_share(c, np.cos(2.0 * np.pi * m * j / 40), np.full(c.n, 1.2e-3))
    kh = 2.0 / 1e-3 * np.sin(np.pi * m / 40) * 1.2e-3
    assert res["share"] == pytest.approx(1.0 if kh > fw.KH_SHORT else 0.0, abs=1e-15)    # 2.5e-30 leaks at m = 1
    assert res["inertial"] == pytest.approx(float(fw.inertial_factor(kh)) - 1.0, rel=1e-12, abs=1e-15)
    assert res["viscous"] == pytest.approx(float(fw.viscous_factor(kh)) - 1.0, rel=1e-12, abs=1e-15)
    assert res["power"].sum() == pytest.approx(80e-6, rel=1e-12)
    assert res["k"][np.argmax(res["power"])] == pytest.approx(kh / 1.2e-3, rel=1e-12)


def test_the_share_of_two_waves_is_the_short_one_s_power_over_both_and_a_vector_field_counts_every_component():
    """a cos of m = 1 (kh 0.19) plus b cos of m = 8 (kh 1.41): the share is b^2 / (a^2 + b^2) and the inertial excess
    the factors' power-weighted mean, by hand; a field of three Cartesian components sums their powers, so the same
    waves along x and z give the same share; and over a region the modes are the region's own, its columns' depths
    weighted by each mode's square."""
    c = cylinder()
    j = np.arange(c.n) % 40
    lo, hi = np.cos(2.0 * np.pi * j / 40), np.cos(2.0 * np.pi * 8 * j / 40)
    a, b = 2.0, 0.5
    res = fw.wavelength_share(c, a * lo + b * hi, np.full(c.n, 1.2e-3))
    kl, kh_ = (2e3 * np.sin(np.pi * m / 40) * 1.2e-3 for m in (1, 8))
    assert res["share"] == pytest.approx(b ** 2 / (a ** 2 + b ** 2), rel=1e-12)
    mean = (a ** 2 * (fw.inertial_factor(kl) - 1.0) + b ** 2 * (fw.inertial_factor(kh_) - 1.0)) / (a ** 2 + b ** 2)
    assert res["inertial"] == pytest.approx(float(mean), rel=1e-12)
    vec = np.stack([a * lo, np.zeros(c.n), b * hi], axis=1)
    assert fw.wavelength_share(c, vec, np.full(c.n, 1.2e-3))["share"] == pytest.approx(res["share"], rel=1e-12)
    ring = np.arange(c.n) < 120                                       # the first three rings: a cylinder 2 mm long
    depth = np.where(np.arange(c.n) < 80, 1e-3, 2e-3)
    sub = fw.wavelength_share(c, hi, depth, region=ring)
    k, phi, idx = fw.modes(c, ring)
    assert len(idx) == 120 and sub["h"] == pytest.approx((c.area[idx][:, None] * phi ** 2 * depth[idx][:, None]).sum(0))
```

- [ ] **Step 3: Run them to see them fail.** `"$PY" -m pytest tests/test_reentry_model_skinbody.py
  tests/test_flow_wavelengths.py -q -k "plate and not hot_plate or wavelength or overstatement or sinusoid or two_waves"`.
  Expected: `ModuleNotFoundError: No module named 'flow_wavelengths'` collecting `tests/test_flow_wavelengths.py`; and,
  `tests/test_reentry_model_skinbody.py` alone, 3 failed, 20 deselected --
  - `AttributeError: module 'skin_plate' has no attribute 'PUSHES'` (2: the flow-off plate's refusals, after its
    comparison with the hand-built body has passed on sub-plan 20's script, and the pushed plate);
  - `TypeError: run() got an unexpected keyword argument 'flow'` (the plate with nothing flowing).

- [ ] **Step 4: The flow on the plate** in `analysis/skin_plate.py`:

```diff
--- a/analysis/skin_plate.py
+++ b/analysis/skin_plate.py
@@ -38,8 +38,36 @@
 numerical delay of the same size. `--swell on` runs the skins' own law; the swelling itself is verified at level 1
 (`tests/test_reentry_model_skin.py`, against the one-dimensional reference with the skins' density law).
 
+With the flow on (`--flow on`, sub-plan 21 Task 8; spec 11.1's last sentence: "with flow on and a push along the
+face: mass, momentum and energy close exactly and the result is physically sensible") the liquid stays in the skins
+and flows; nothing is removed as it forms. Two settings drive it:
+
+- `--push flight`: the face as a patch of the 100 mm sphere's shoulder PUSH_ALPHA (45 degrees) from the stagnation
+  point at 93 s of its flight, where the surface flow's pressure gradient is largest. The flight direction is tilted by
+  that angle about the z axis (`inclined`), so that along the face the meridian runs to -y, and every column of the
+  face (`face_columns`: its corners at x = 0 whose normal is the face's, the edges' columns not) takes that point's
+  forcing (`plate_forcing`): dp/ds -3.42e5 Pa/m and 80.3 Pa of gas shear down the meridian, and the deceleration of
+  29.8 m/s2 along the flight direction -- 21.1 m/s2 out of the face, the front face's Rayleigh-Taylor pull, and 21.1
+  m/s2 along it toward the nose, against the gas. On the plate's own flight direction, normal to the face, neither the
+  pressure gradient nor the shear has a direction along it and the deceleration pulls straight out. The flight
+  direction feeds nothing else here: the heating is the plate's own projected flux, and the stagnation and rear
+  patches it picks are not read. `--push none` (the default) drives nothing but surface tension.
+- `--release all`: Task 7's release from the top cells at RELEASE_ALL, far more than any column holds in the shortest
+  flow sub-step, so that the liquid leaves as it forms -- at the liquidus (Task 7's ruling 3) rather than at T_feed,
+  2 K higher, where instant removal takes it. With no push this is the flow-off plate's physics with the flow on, the
+  regression against the references.
+
+The flow sub-steps are chosen by the skins (Task 6d) or fixed (`--flow-cycles`). With the flow on the history gains the
+skin body's flow columns (the momentum residual, the fastest cell, the roughness, the flow sub-steps' worst Courant
+number, growth and linearisations, the released mass), and `--columns PATH.npz` keeps every column's mass, thickness,
+top, top temperature, depths, top speed and volume flux at every macro step, for the sensibility checks;
+`--save-state T` keeps the skins, the columns, the heating and the base at the first macro step ending at or after T
+(`<out>.state.npz`), for the flow sub-step study and the cost split on the skins alone.
+
     python analysis/skin_plate.py [--case cold|hot] [--dt 0.05] [--preset dev|production] [--h-mm 2] [--t-end 2]
                                   [--swell off|on] [--fine-mm 4] [--fine-size-mm 0.7] [--out PATH.json]
+                                  [--flow off|on] [--push none|flight] [--release none|all] [--flow-cycles N]
+                                  [--columns PATH.npz] [--save-state T]
 """
 import argparse
 import json
@@ -59,9 +87,25 @@
 # R_FLAT - R_READ is half again the diffusion length 2 sqrt(alpha t), about 19 mm at 2 s: with full flux only to 12 mm
 # the centre lost heat sideways and fell 9 % behind the exact solution by 2 s (measured while planning)
 REFERENCE_1D = {1.0: 0.7890, 1.5: 2.7120, 2.0: 5.2190}       # mm, analysis/ablation_1d.py at 1 um and 0.5 ms (2026-10-12)
+# The references with the liquid leaving at the liquidus, as the release from the top cells takes it (sub-plan 21 Task 7's
+# ruling 3), rather than at T_feed (stage 2's removal): analysis/ablation_1d.py --remove-at liquidus at 1 um and 0.5 ms
+# (2026-10-24). The cold plate's recession [mm], and the hot plate's recession and non-rigid depth [mm]; the fine box
+# is moved by the one-dimensional difference (`references`)
+REFERENCE_1D_LIQUIDUS = {1.0: 0.8370, 1.5: 2.7720, 2.0: 5.2980}
+HOT_LIQUIDUS = {1.0: (0.9280, 1.3530), 1.5: (1.9660, 1.4330), 2.0: (3.0750, 1.5370)}
 CASES = {"cold": (T0, Q0), "hot": (880.0, 2.0e6)}            # K, W/m2: the level-2 plate and the hot plate (2026-10-14)
 HOT_FLOWING_TOL, HOT_BASE_TOL_K = 0.10, 1.0                    # the hot plate's flowing depth (relative) and base excess
 FINE = (4.0e-3, 0.7e-3)                                         # m: the fine band's depth and size, as a skin run's
+# The flow (sub-plan 21 Task 8): the push along the face, the face as a patch of the 100 mm sphere's shoulder at
+# PUSH_ALPHA from the stagnation point at 93 s of the flight -- the surface flow there (dp/ds [Pa/m], gas shear [Pa]) and
+# the deceleration [m/s2], evaluated 2026-10-24 (measure/task8/forcing45.py): the largest pressure gradient of 6c's
+# forcing, 3.4e5 Pa/m between 44 and 46.5 degrees, the shear 87 Pa at its largest at 49 degrees
+PUSH_ALPHA = 45.0                                               # deg
+PUSHES = {"none": None, "flight": (PUSH_ALPHA, -3.424e5, 80.29, 29.77)}
+RELEASES = {"none": None, "all": 1.0e9}  # kg/(m2 s): a top area of 4 mm2 then releases 4 kg in a 1 ms flow sub-step,
+                                         # some 1e5 times what a column holds: the liquid leaves as it forms
+FACE_TILT = 5.0                          # deg: a column of the face has its corner at x = 0 and its normal within this
+                                         # of the face's; the edges' columns, tilted 45 degrees, take no gas forcing
 
 
 def hot_reference():
@@ -76,6 +120,19 @@
     here = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
     doc = json.load(open(os.path.join(here, "data", "reference_values", "box_q6p_fine.json")))
     return {t: float(np.interp(t, doc["t_s"], doc["r_total_mm"])) for t in (1.0, 1.5, 2.0)}
+
+
+def references(case, at_liquidus=False):
+    """The references a plate is judged against at 1.0, 1.5 and 2.0 s: (one-dimensional, fine box) recession [mm] for the
+    cold plate, {t: (recession, non-rigid depth)} [mm] and None for the hot one; `at_liquidus`: with the liquid leaving
+    at the liquidus (REFERENCE_1D_LIQUIDUS, HOT_LIQUIDUS), the fine box moved by the one-dimensional difference -- the
+    liquid layer that stage 2's T_feed keeps on the front."""
+    if case == "hot":
+        return (HOT_LIQUIDUS if at_liquidus else hot_reference()), None
+    fine = fine_box()
+    if not at_liquidus:
+        return REFERENCE_1D, fine
+    return REFERENCE_1D_LIQUIDUS, {t: fine[t] + REFERENCE_1D_LIQUIDUS[t] - REFERENCE_1D[t] for t in fine}
 
 
 def taper(r, r_flat=R_FLAT, r_zero=R_ZERO):
@@ -105,31 +162,100 @@
     return c.points + (body.skins.column_mass() / (body.material.rho * c.area))[:, None] * c.normal
 
 
+def inclined(alpha_deg):
+    """The flight direction tilted by alpha [deg] from the plate's own (-1, 0, 0) about the z axis: the heated face then
+    lies alpha from the stagnation point, which is off the plate toward +y, and its meridian runs to -y."""
+    a = np.radians(alpha_deg)
+    return (-np.cos(a), np.sin(a), 0.0)
+
+
+def face_columns(body, tilt=FACE_TILT):
+    """The heated face's columns, (n,) bool: corners at x = 0 at time zero whose normal is within `tilt` [deg] of the
+    face's (-1, 0, 0); the columns on the plate's edges, whose normals lean 45 degrees, are not the face's."""
+    return (np.abs(body._top0[:, 0]) < 0.05e-3) & (-body._normal0[:, 0] >= np.cos(np.radians(tilt)))
+
+
+def plate_forcing(body, face, push):
+    """The gas side's push on the plate (skinbody.FlowForcing): with `push` (alpha [deg], dp/ds [Pa/m], gas shear [Pa],
+    deceleration [m/s2]; PUSHES) the surface flow of that point of the shoulder on every column of the face and the
+    deceleration along the body's flight direction, which `inclined(alpha)` gave it; None: nothing (surface tension
+    alone drives the flow)."""
+    n = body.columns.n
+    if push is None:
+        return skinbody.FlowForcing(np.zeros(n), np.zeros(n), np.zeros(3))
+    _, dp_ds, tau, decel = push
+    return skinbody.FlowForcing(np.where(face, dp_ds, 0.0), np.where(face, tau, 0.0), decel * body.v_hat)
+
+
+def column_record(body):
+    """Every column's state for the sensibility checks: its mass [kg], thickness [m], top [m, (n, 3)], top temperature
+    [K], flowing and liquid depths [m], top speed [m/s] and volume flux [m2/s, (n, 3)] -- the sum of its cells'
+    thicknesses times their velocities, what crosses a unit width of the surface."""
+    sk, mat = body.skins, body.material
+    dz = sk.dz()
+    return {"mass": sk.column_mass(), "thickness": dz.sum(axis=1), "top": body.surface.points.copy(),
+            "T_top": sk.T[:, 0].copy(), "flowing": sk.depth_above(mat.T_rigid), "liquid": sk.depth_above(mat.T_liquidus),
+            "top_speed": np.linalg.norm(sk.U[:, 0], axis=1), "flux": (dz[..., None] * sk.U).sum(axis=1)}
+
+
+def save_state(path, body, q):
+    """The skins, their columns (base corners, triangles, top areas), the heating q [W/m2], the base temperatures and the
+    3D model's gradient below them: what `skin.advance` needs to run a part on the skins alone."""
+    np.savez(path, **dict(zip(("area", "m", "H", "T", "count", "U", "flow_memory"), body.skins.state())),
+             points=body.columns.points, faces=body.columns.faces, A_top=body.surface.areas, q=q,
+             T_b=body.solver.temperature()[body.moving.corners], g_base=body._normal_gradient(), v_hat=body.v_hat,
+             rho_l=body.skins.rho_l, face=face_columns(body))
+
+
 def run(dt=0.05, preset="dev", deform=4e-3, t_end=2.0, solver="skfem", verbose=True, h=2e-3, width=WIDTH, r_flat=R_FLAT,
-        r_zero=R_ZERO, repeat_rule=skinbody.DEFAULT_REPEAT_RULE, T_start=T0, q0=Q0, swell=False, fine=FINE):
+        r_zero=R_ZERO, repeat_rule=skinbody.DEFAULT_REPEAT_RULE, T_start=T0, q0=Q0, swell=False, fine=FINE, flow=False,
+        push=None, release=None, flow_cycles=None, record=False, save_at=None, save_path=None):
+    """The plate (module docstring) to t_end at macro steps dt; returns (body, rows, centre columns, failure) and, with
+    `record`, every macro step's `column_record` as a fifth item. `flow` turns the flow on (the liquid kept), driven by
+    `push` (PUSHES) with the release `release` [kg/(m2 s), RELEASES] and `flow_cycles` (None: chosen); `save_at` [s]
+    keeps the skins' state at the first macro step ending at or after it in `save_path`."""
+    if not flow and (push is not None or release is not None or flow_cycles is not None):
+        raise ValueError("a push, a release and the flow sub-steps need the flow on")
     mat = material.Material.from_drama_json(material.MATERIAL_NAMES["AA7075_scheil"])
     m = mesh.box_mesh(LENGTH, width, width, h)
     m.points[:, 1:] -= 0.5 * width                                 # the face's centre on the x axis
     if fine is not None and fine[0] > 0.0:                         # the fine band (depth, size) m, as a skin run's
         m = movemesh.refined_mesh(m, fine[0], min(fine[1], h))
-    body = skinbody.SkinBody(m, mat, thermal_solver(solver), mat.rho * m.volume(),
-                             skin.SkinSettings(cell=skin.SKIN_PRESETS[preset], swell=swell), T0=T_start, emissivity=0.0,
-                             v_hat=(-1.0, 0.0, 0.0), deform_depth=deform, repeat_rule=repeat_rule)
+    settings = skin.SkinSettings(cell=skin.SKIN_PRESETS[preset], swell=swell)
+    if flow:                                                       # the liquid stays and flows (sub-plan 21 Task 8)
+        body = skinbody.SkinBody(m, mat, thermal_solver(solver), mat.rho * m.volume(), settings, T0=T_start,
+                                 emissivity=0.0, v_hat=(-1.0, 0.0, 0.0) if push is None else inclined(push[0]),
+                                 deform_depth=deform, repeat_rule=repeat_rule, remove=False, flow=True,
+                                 flow_cycles=flow_cycles, release=None if release is None else skin.prescribed_release(release))
+    else:
+        body = skinbody.SkinBody(m, mat, thermal_solver(solver), mat.rho * m.volume(), settings, T0=T_start,
+                                 emissivity=0.0, v_hat=(-1.0, 0.0, 0.0), deform_depth=deform, repeat_rule=repeat_rule)
     p = body._top0                                                 # the top surface at time zero: the plate's own
     front = np.abs(p[:, 0]) < 0.05e-3
     radius = np.where(front, np.hypot(p[:, 1], p[:, 2]), np.inf)
     centre = radius <= R_READ
+    face = face_columns(body)
     solid0, area0 = solid_top(body), body.surface.areas.copy()
     M0, rows, t, started, failure = body.mass(), [], 0.0, time.perf_counter(), None
+    records = [column_record(body)] if record else None
     for k in range(int(round(t_end / dt))):
         try:
             q = projected_flux(body, radius, r_flat, r_zero, q0)
-            body.advance(t + dt, dt, heating.HeatingResult(q, q0, 0.0, 0.0, 0.0))
+            if flow:
+                body.advance(t + dt, dt, heating.HeatingResult(q, q0, 0.0, 0.0, 0.0),
+                             forcing=plate_forcing(body, face, push))
+            else:
+                body.advance(t + dt, dt, heating.HeatingResult(q, q0, 0.0, 0.0, 0.0))
         except RuntimeError as exc:                        # a move no rebuild can absorb: record where, and stop
             failure = {"t": t + dt, "error": str(exc), "fold_deg": body.moving.surface_fold()}
             print(json.dumps(failure), flush=True)
             break
         t += dt
+        if save_at is not None and t >= save_at - 1e-9 and save_path is not None:
+            save_state(save_path, body, q)
+            save_at = None
+        if record:
+            records.append(column_record(body))
         st = body.skin_stats()
         rows.append({"t": t, "recession_mm": float(((solid_top(body) - solid0) @ np.array([1.0, 0.0, 0.0]))[centre].mean() * 1e3),
                      "top_recession_mm": float(body._recession()[centre].mean() * 1e3),
@@ -138,7 +264,7 @@
                      "T_top_centre": float(body.surface_temperature()[centre].mean()),
                      "flowing_mm": float(body.skins.depth_above(mat.T_rigid)[centre].mean() * 1e3),
                      "base_excess_K": float((body.solver.temperature()[body.moving.corners][centre] - mat.T_rigid).max()),
-                     "mass_error": (body.mass() + body.removed_mass) / M0 - 1.0,
+                     "mass_error": (body.mass() + body.removed_mass + body.released_mass) / M0 - 1.0,
                      "energy_residual": body.energy_balance_residual(), "fold_deg": st["surface_fold_deg"],
                      "quality": st["surface_quality_min"], "zone_ratio": st["zone_ratio_min"], "rebuilds": st["rebuilds"],
                      "faults": st["jump_faults"], "repeats": st["interface_repeats"], "cells": st["skin_cells"],
@@ -146,8 +272,20 @@
                      "total_ratio": st["mismatch_total_ratio"],
                      "absorbed_kg": st["absorbed_mass_kg"], "handed_kg": st["handed_mass_kg"],
                      "wall_s": time.perf_counter() - started})
+        if flow:                                           # the skin body's flow columns (sub-plan 21 Tasks 6c-7)
+            sk = body.skins
+            rows[-1].update({key: st[key] for key in ("momentum_residual", "speed_max_ms", "roughness_max",
+                                                      "flow_courant_max", "flow_growth_max", "flow_passes_max",
+                                                      "flow_cycles_max", "released_mass_kg", "wall_flow_s",
+                                                      "liquid_depth_max_mm", "flowing_depth_max_mm")})
+            rows[-1].update({"liquid_mm": float(sk.depth_above(mat.T_liquidus)[centre].mean() * 1e3),
+                             "kinetic_J": 0.5 * float((sk.m * (sk.U ** 2).sum(axis=2)).sum()),
+                             "wall_3d_s": st["wall_3d_s"], "wall_skins_s": st["wall_skins_s"],
+                             "wall_move_s": st["wall_move_s"], "flow_seconds": dict(skin.FLOW_SECONDS)})
         if verbose and (k + 1) % max(1, int(round(0.25 / dt))) == 0:
             print(json.dumps({key: (round(v, 6) if isinstance(v, float) else v) for key, v in rows[-1].items()}), flush=True)
+    if record:
+        return body, rows, int(centre.sum()), failure, records
     return body, rows, int(centre.sum()), failure
 
 
@@ -170,25 +308,43 @@
                    help="off: the skins at the solid's density, as the references hold their metal (default); on: the "
                         "skins' own law, each cell at its phase's density")
     a.add_argument("--out", default="skin_plate.json")
+    a.add_argument("--flow", choices=("off", "on"), default="off",
+                   help="on: the skins' liquid stays and flows (sub-plan 21 Task 8); off: removed as it forms (default)")
+    a.add_argument("--push", choices=tuple(PUSHES), default="none",
+                   help="with --flow on, flight: the face as the 100 mm sphere's shoulder at 45 degrees at 93 s")
+    a.add_argument("--release", choices=tuple(RELEASES), default="none",
+                   help="with --flow on, all: the liquid released from the top cells as it forms")
+    a.add_argument("--flow-cycles", type=int, default=None, help="with --flow on: flow sub-steps fixed (default chosen)")
+    a.add_argument("--columns", default=None, help="keep every column's state at every macro step in this npz")
+    a.add_argument("--save-state", type=float, default=None, help="keep the skins' state at this time [s]")
     args = a.parse_args()
+    if args.flow == "off" and (args.push != "none" or args.release != "none" or args.flow_cycles is not None):
+        a.error("--push, --release and --flow-cycles need --flow on")
     T_start, q0 = CASES[args.case]
-    body, rows, n_centre, failure = run(args.dt, args.preset, args.deform_mm * 1e-3, args.t_end, args.solver,
-                                        h=args.h_mm * 1e-3, width=args.width_mm * 1e-3, r_flat=args.r_flat_mm * 1e-3,
-                                        r_zero=args.r_zero_mm * 1e-3, repeat_rule=args.repeat_rule, T_start=T_start, q0=q0,
-                                        swell=args.swell == "on", fine=(args.fine_mm * 1e-3, args.fine_size_mm * 1e-3))
+    out = run(args.dt, args.preset, args.deform_mm * 1e-3, args.t_end, args.solver,
+              h=args.h_mm * 1e-3, width=args.width_mm * 1e-3, r_flat=args.r_flat_mm * 1e-3,
+              r_zero=args.r_zero_mm * 1e-3, repeat_rule=args.repeat_rule, T_start=T_start, q0=q0,
+              swell=args.swell == "on", fine=(args.fine_mm * 1e-3, args.fine_size_mm * 1e-3), flow=args.flow == "on",
+              push=PUSHES[args.push], release=RELEASES[args.release], flow_cycles=args.flow_cycles,
+              record=args.columns is not None, save_at=args.save_state,
+              save_path=os.path.splitext(args.out)[0] + ".state.npz")
+    body, rows, n_centre, failure = out[:4]
+    if args.columns is not None:
+        np.savez(args.columns, t=np.r_[0.0, [x["t"] for x in rows]], face=face_columns(body), top0=body._top0,
+                 area=body.columns.area, **{key: np.array([rec[key] for rec in out[4]]) for key in out[4][0]})
     at = lambda t_, key: float(np.interp(t_, [x["t"] for x in rows], [x[key] for x in rows]))
     compare, verdict = {}, []
-    fine = fine_box() if args.case == "cold" else None
-    hot = hot_reference() if args.case == "hot" else None
+    ref, fine = references(args.case, at_liquidus=args.release == "all")     # the release takes liquid at the liquidus
+    hot = ref if args.case == "hot" else None
     for t_ in (1.0, 1.5, 2.0):
         if not rows or t_ > rows[-1]["t"] + 1e-9:
             continue
         r, g = at(t_, "recession_mm"), at(t_, "removed_recession_mm")
         if fine is not None:
-            compare[str(t_)] = {"centre_mm": r, "removed_mm": g, "fine_mm": fine[t_], "ref1d_mm": REFERENCE_1D[t_],
-                                "rel_fine": r / fine[t_] - 1.0, "rel_1d": r / REFERENCE_1D[t_] - 1.0,
-                                "removed_rel_1d": g / REFERENCE_1D[t_] - 1.0}
-            verdict.append(abs(r / fine[t_] - 1.0) <= 0.05 and abs(r / REFERENCE_1D[t_] - 1.0) <= 0.05)
+            compare[str(t_)] = {"centre_mm": r, "removed_mm": g, "fine_mm": fine[t_], "ref1d_mm": ref[t_],
+                                "rel_fine": r / fine[t_] - 1.0, "rel_1d": r / ref[t_] - 1.0,
+                                "removed_rel_1d": g / ref[t_] - 1.0}
+            verdict.append(abs(r / fine[t_] - 1.0) <= 0.05 and abs(r / ref[t_] - 1.0) <= 0.05)
         else:
             ref_r, ref_nr = hot[t_]
             flow = at(t_, "flowing_mm")
@@ -197,7 +353,9 @@
                                 "base_excess_K": at(t_, "base_excess_K")}
             verdict.append((t_ < 1.25 or abs(r / ref_r - 1.0) <= 0.05) and abs(flow / ref_nr - 1.0) <= HOT_FLOWING_TOL
                            and at(t_, "base_excess_K") <= HOT_BASE_TOL_K)
-    summary = {"case": args.case, "T0_K": T_start, "q_Wm2": q0, "passes": bool(verdict) and all(verdict) and failure is None,
+    judged = args.push == "none"                     # a push along the face has no exact reference (spec 11.1)
+    summary = {"case": args.case, "T0_K": T_start, "q_Wm2": q0,
+               "passes": judged and bool(verdict) and all(verdict) and failure is None,
                "failure": failure, "completed_s": rows[-1]["t"] if rows else 0.0, "h_mm": args.h_mm, "dt": args.dt,
                "width_mm": args.width_mm, "r_flat_mm": args.r_flat_mm, "r_zero_mm": args.r_zero_mm,
                "repeat_rule": args.repeat_rule, "swell": args.swell, "fine_mm": args.fine_mm,
@@ -209,6 +367,13 @@
                "rebuilds": rows[-1]["rebuilds"] if rows else 0, "faults": rows[-1]["faults"] if rows else 0,
                "repeats": rows[-1]["repeats"] if rows else 0, "wall_s": rows[-1]["wall_s"] if rows else 0.0,
                "wall_split": body.wall}
+    if args.flow == "on":
+        summary.update({"flow": "on", "push": args.push, "release": args.release, "flow_cycles": args.flow_cycles,
+                        "judged": judged, "released_kg": body.released_mass,
+                        "worst_momentum_residual": max(abs(x["momentum_residual"]) for x in rows) if rows else 0.0,
+                        "roughness_max": max(x["roughness_max"] for x in rows) if rows else 0.0,
+                        "speed_max_ms": max(x["speed_max_ms"] for x in rows) if rows else 0.0,
+                        "flow_seconds": dict(skin.FLOW_SECONDS)})
     print(json.dumps(summary, indent=1))
     json.dump({"summary": summary, "rows": rows}, open(args.out, "w"), indent=1)
 
```

- [ ] **Step 5: The one-dimensional reference with the liquid leaving at the liquidus** in `analysis/ablation_1d.py`:

```diff
--- a/analysis/ablation_1d.py
+++ b/analysis/ablation_1d.py
@@ -13,8 +13,15 @@
 swelling; `nonrigid` swells only the cells at or above T_rigid, the metal the skins hold above a 3D model of one
 density. Lengths are reported as solid metal (cells times `cell`).
 
+`--remove-at liquidus` (sub-plan 21 Task 8, 2026-10-24): the front cells leave at or above the liquidus (908 K) instead
+of T_feed (910 K), as sub-plan 21 Task 7's release from the top cells takes the liquid. Stage 2's T_feed keeps a liquid
+layer on the front across which the last 2 K fall, about 2 K k_l / q thick (0.03 mm at 6 MW/m2, 0.09 mm at 2 MW/m2),
+and the slab's melting starts that much later: at the liquidus the cold plate's slab recedes 0.048, 0.060 and 0.079 mm
+further by 1.0, 1.5 and 2.0 s (6.1, 2.2 and 1.5 %) and the hot plate's 0.139, 0.149 and 0.158 mm (17.6, 8.2 and 5.4 %),
+its non-rigid depth 0.13-0.14 mm (9 %) thinner (1 um, 0.5 ms).
+
     python analysis/ablation_1d.py [--cell-um 1] [--dt-ms 0.5] [--q 6e6] [--t0 700] [--length-mm 8] [--t-end 2]
-                                   [--swell off|all|nonrigid] [--out PATH]
+                                   [--swell off|all|nonrigid] [--remove-at feed|liquidus] [--out PATH]
 """
 import argparse
 import json
@@ -29,6 +36,7 @@
 
 
 SWELL_MODES = ("off", "all", "nonrigid")
+REMOVE_AT = ("feed", "liquidus")      # the front cells leave at or above T_feed (stage 2's) or the liquidus (Task 7's)
 
 
 def thickness(mat, T, cell, swell):
@@ -40,11 +48,16 @@
     return dz if swell == "all" else np.where(T >= mat.T_rigid, dz, cell)
 
 
-def run(mat, q=6.0e6, T0=700.0, length=8.0e-3, cell=1.0e-6, dt=5.0e-4, t_end=2.0, every=0.05, swell="off"):
+def run(mat, q=6.0e6, T0=700.0, length=8.0e-3, cell=1.0e-6, dt=5.0e-4, t_end=2.0, every=0.05, swell="off",
+        remove_at="feed"):
     """Recession [m] of the slab's front and its non-rigid depth [m] (front to T_rigid), both as solid metal, at
-    multiples of `every` [s], as (times, recession, non-rigid depth) arrays."""
+    multiples of `every` [s], as (times, recession, non-rigid depth) arrays; the front cells leave at or above T_feed or,
+    with `remove_at` "liquidus", the liquidus."""
     if swell not in SWELL_MODES:
         raise ValueError("swell must be one of {}, got {!r}".format(SWELL_MODES, swell))
+    if remove_at not in REMOVE_AT:
+        raise ValueError("remove_at must be one of {}, got {!r}".format(REMOVE_AT, remove_at))
+    T_remove = mat.T_feed if remove_at == "feed" else mat.T_liquidus
     n = int(round(length / cell))
     T = np.full(n, T0)
     front = 0                                            # index of the first cell still in the slab
@@ -80,7 +93,7 @@
         else:
             raise RuntimeError("the 1D reference did not converge at step {}".format(step))
         T[front:] = Tk
-        liquid = np.cumprod(Tk >= mat.T_feed)
+        liquid = np.cumprod(Tk >= T_remove)
         front += int(liquid.sum())
         if step % record == 0:
             times.append(step * dt)
@@ -99,12 +112,14 @@
     p.add_argument("--length-mm", type=float, default=8.0)
     p.add_argument("--t-end", type=float, default=2.0)
     p.add_argument("--swell", choices=SWELL_MODES, default="off")
+    p.add_argument("--remove-at", choices=REMOVE_AT, default="feed")
     p.add_argument("--out", default="ablation_1d.json")
     a = p.parse_args()
     mat = material.Material.from_drama_json(material.MATERIAL_NAMES["AA7075_scheil"])
-    t, r, nr = run(mat, a.q, a.t0, a.length_mm * 1e-3, a.cell_um * 1e-6, a.dt_ms * 1e-3, a.t_end, swell=a.swell)
+    t, r, nr = run(mat, a.q, a.t0, a.length_mm * 1e-3, a.cell_um * 1e-6, a.dt_ms * 1e-3, a.t_end, swell=a.swell,
+                   remove_at=a.remove_at)
     doc = {"q_Wm2": a.q, "T0_K": a.t0, "length_mm": a.length_mm, "cell_um": a.cell_um, "dt_ms": a.dt_ms, "swell": a.swell,
-           "material": mat.name, "t_s": [round(x, 6) for x in t.tolist()],
+           "remove_at": a.remove_at, "material": mat.name, "t_s": [round(x, 6) for x in t.tolist()],
            "recession_mm": [round(x, 6) for x in (r * 1e3).tolist()], "nonrigid_mm": [round(x, 6) for x in (nr * 1e3).tolist()]}
     json.dump(doc, open(a.out, "w"), indent=0)
     for t_ in (1.0, 1.5, 2.0):
```

- [ ] **Step 6: The short-wave measure** -- create `analysis/flow_wavelengths.py`:

```python
"""Sub-plan 21 Task 8: the share of the columns' flow at short wavelengths, for the depth-uniform pressure's limitation
(layered-skin design section 14, assumption 1; Asha, 2026-10-21, option a).

The columns give every layer of a column the same push along the surface (spec 6.4): the pressure varies with depth only
through the deceleration's component into the surface. That is exact for a disturbance whose wavelength is long against
the flowing depth h, kh small (k = 2 pi / lambda). Where kh is not small a disturbance really stirs only the top 1/k or
so of the layer, and the columns' response to it is too fast: its rate -- the growth of a Rayleigh-Taylor mode, the
decay of a capillary ripple, the speed of a wave -- is overstated by

    inertial layer (liquid; its viscous time h^2 / nu is 2.7 s at 1.2 mm):   sqrt(kh / tanh kh)       (spec 14.1)
    viscous layer (slurry; 0.01-0.04 s at 1-3.4 mm and 1-5 Pa s):            (2/3) (kh)^3 (cosh 2kh + 2 (kh)^2 + 1)
                                                                              / (sinh 2kh - 2kh)

the second the ratio of the lubrication rate h^3 k^2 (rho g + sigma k^2) / (3 mu), which is what columns with a viscous
profile and a depth-uniform pressure give, to the Stokes-flow rate of a viscous layer on a rigid base,
(rho g + sigma k^2) (sinh 2kh - 2kh) / (2 mu k (cosh 2kh + 2 (kh)^2 + 1)). At kh = 0.5 they are 1.04 and 1.45, at kh = 1
1.15 and 2.77. Both are 1 for a long wave: the runoff along the meridian, which varies over tens of millimetres, is
where they are a few per cent at most.

The measure (ruling, 2026-10-24): the share of the flowing layers' volume flux q -- each column's sum of its cells'
thicknesses times their velocities [m2/s], what crosses a unit width of the surface -- carried at wavelengths where kh
exceeds KH_SHORT, over the columns whose flowing depth is at least THICK (1 mm). The flux is the measure because it is
what the columns compute wrongly where kh is not small, and what moves the metal: where the runoff takes it and how deep
the flowing layer is that Girin's spraying reads (spec 6.7, 8). A short-wave share of the flux s means that at most the
fraction s of where the metal goes is subject to the overstated rates; the flux-weighted mean excess of the factors
above over the modes, sum P_n (F(k_n h_n) - 1) / sum P_n, says by how much on average. The top surface's height is
measured the same way (its variance at short wavelengths: the shape of the layer, which the rates move). KH_SHORT is
where the inertial factor reaches 1.05 (kh = 0.56); the viscous factor does at kh = 0.17.

Wavelengths on the columns are the modes of their surface Laplacian over the region (`modes`): L phi = k^2 M phi, L the
cotangent Laplacian on the region's edges (each edge's clipped weight, `Columns.conduct_weight`, with free ends at the
region's border), M the columns' areas, so that each mode is a wave of wavenumber k on the surface -- on a flat grid of
squares a sampled cosine of wavenumber kappa is one exactly, with k = (2 / dx) sin(kappa dx / 2). A field f (per column,
or per column and Cartesian component) has the power P_n = sum over components of (phi_n^T M f)^2 in mode n, and each
mode its own depth, h_n = sum_j A_j phi_nj^2 h_j (the modes are M-orthonormal)."""
import numpy as np
import scipy.linalg as sla
import scipy.sparse as sp

KH_SHORT = 0.56            # kh above which the columns overstate an inertial layer's rates by more than 5 %
THICK = 1.0e-3             # m: the flowing layers measured (spec 14.1: "in layers over about 1 mm")


def inertial_factor(kh):
    """sqrt(kh / tanh kh): how much the columns overstate an inertial layer's rates at kh (1 at kh = 0)."""
    kh = np.asarray(kh, dtype=float)
    safe = np.where(kh > 1e-8, kh, 1.0)
    return np.where(kh > 1e-8, np.sqrt(safe / np.tanh(safe)), 1.0)


def viscous_factor(kh):
    """The lubrication rate over the Stokes rate of a viscous layer on a rigid base at kh (module docstring; 1 at kh = 0):
    how much the columns overstate a viscous layer's rates. Its series 1 + (9/5) (kh)^2 below kh = 0.05, where the
    difference sinh 2kh - 2kh loses its digits."""
    kh = np.asarray(kh, dtype=float)
    x = np.where(kh > 0.05, kh, 1.0)
    exact = (2.0 / 3.0) * x ** 3 * (np.cosh(2 * x) + 2 * x ** 2 + 1.0) / (np.sinh(2 * x) - 2 * x)
    return np.where(kh > 0.05, exact, 1.0 + 1.8 * kh ** 2)


def modes(columns, region=None):
    """The surface Laplacian's modes over the columns `region` ((n,) bool, all if None): (k [1/m, ascending], phi
    (m, m) M-orthonormal by columns, idx) for the m columns idx of the region, from L phi = k^2 M phi -- L the cotangent
    Laplacian on the edges with both ends in the region, its clipped weights, M the columns' areas (module docstring).
    Dense: about 30 s for 4 000 columns."""
    region = np.ones(columns.n, dtype=bool) if region is None else np.asarray(region, dtype=bool)
    idx = np.flatnonzero(region)
    local = np.full(columns.n, -1)
    local[idx] = np.arange(len(idx))
    i, j = columns.edges[:, 0], columns.edges[:, 1]
    keep = region[i] & region[j]
    a, b, w = local[i[keep]], local[j[keep]], columns.conduct_weight[keep]
    m = len(idx)
    L = sp.csr_matrix((np.r_[w, w, -w, -w], (np.r_[a, b, a, b], np.r_[a, b, b, a])), shape=(m, m)).toarray()
    lam, phi = sla.eigh(L, np.diag(columns.area[idx]))
    return np.sqrt(np.clip(lam, 0.0, None)), phi, idx


def wavelength_share(columns, field, depth, region=None, kh_short=KH_SHORT, basis=None):
    """The share of `field` ((n,) or (n, 3), per column) carried at kh above `kh_short` over the columns `region` (all if
    None), and the field-weighted mean excess of the columns' rates over the exact ones (module docstring), each mode at
    its own depth from `depth` [m, (n,)]. `basis`: `modes(columns, region)` made already. Returns a dict: share,
    inertial and viscous (the mean excesses, F - 1), and per mode k [1/m], h [m], kh and power."""
    k, phi, idx = modes(columns, region) if basis is None else basis
    f = np.asarray(field, dtype=float)[idx]
    f = f[:, None] if f.ndim == 1 else f
    A = columns.area[idx]
    power = ((phi.T @ (A[:, None] * f)) ** 2).sum(axis=1)
    h = (A[:, None] * phi ** 2 * np.asarray(depth, dtype=float)[idx][:, None]).sum(axis=0)
    kh = k * h
    total = power.sum()
    if total <= 0.0:
        return {"share": 0.0, "inertial": 0.0, "viscous": 0.0, "k": k, "h": h, "kh": kh, "power": power}
    return {"share": float(power[kh > kh_short].sum() / total),
            "inertial": float((power * (inertial_factor(kh) - 1.0)).sum() / total),
            "viscous": float((power * (viscous_factor(kh) - 1.0)).sum() / total),
            "k": k, "h": h, "kh": kh, "power": power}


def column_flux(sk):
    """Each column's volume flux [m2/s, (n, 3)]: its cells' thicknesses times their velocities, summed."""
    return (sk.dz()[..., None] * sk.U).sum(axis=1)


def kh_at(factor, target, lo=1e-6, hi=10.0):
    """The kh at which `factor` (inertial_factor or viscous_factor, increasing) reaches `target`, by bisection."""
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if float(factor(mid)) < target else (lo, mid)
    return 0.5 * (lo + hi)


if __name__ == "__main__":
    print("inertial factor 1.05 at kh = %.3f; viscous 1.05 at kh = %.3f" % (kh_at(inertial_factor, 1.05),
                                                                         kh_at(viscous_factor, 1.05)))
    for kh in (0.1, 0.25, 0.5, 1.0, 2.0):
        print("kh %.2f: inertial %.4f, viscous %.4f" % (kh, inertial_factor(kh), viscous_factor(kh)))
```

- [ ] **Step 7: Run them to see them pass.** Expected: `tests/test_reentry_model_skinbody.py` 21 passed (2 deselected,
  the reference tier's; the three new plates 69 s), `tests/test_flow_wavelengths.py` 5 passed (0.4 s),
  `tests/test_reentry_model_skin.py` 87 as before; the unit tier 506 passed, 1 skipped (2 205 s under the plates' load;
  498 before). Nothing in `reentry_model/` changes, so every flight is Task 7's.

- [ ] **Step 8: The analysis runs** (`PLATE="analysis/skin_plate.py"`, each on mains power under `caffeinate -i`):
  - the regression: `"$PY" "$PLATE" --case cold --dt 0.5 --flow on --release all --swell off --out reg_cold.json`
    and `--case hot` likewise (item 3; 50 and 25 minutes under load), against `--case cold --dt 0.5` and `--case hot
    --dt 0.5` with the flow off;
  - the pushed plate: `"$PY" "$PLATE" --case cold --dt 0.5 --flow on --push flight --swell on --width-mm 170 --t-end 1.0
    --save-state 1.0 --columns push.npz --out push.json` (item 2: 1 h 36 min alone; its coupled books, sensibility and
    per-macro-step cost), the flow-off plate of that width to 1.0 s with its state kept (`--swell on --width-mm 170
    --t-end 1.0 --save-state 1.0`, 19 min), and the skins-alone studies from that state as in
    `measure/task8/plate_skins.py` (items 2, 4, 5, 6);
  - the short-wave segments: `measure/task8/short_waves.py 105|120 flat|cos slurry|liquid 0.1` (item 7).

- [ ] **Step 9: Commit.**

```bash
git add analysis/skin_plate.py analysis/ablation_1d.py analysis/flow_wavelengths.py tests/test_reentry_model_skinbody.py tests/test_flow_wavelengths.py
git commit -m "Skins: level 2 with the flow on -- the plate pushed along its face as a patch of the 100 mm sphere's shoulder at 45 degrees at 93 s, its books exact and its runoff, Couette-Poiseuille speed, freezing deposit, roughness and edges sensible; the release regression against references with the liquid leaving at the liquidus; the flow's cost on the plate; the share of the flow at short wavelengths, analysis/flow_wavelengths.py, about 1 % at 1.2 mm and a quarter at 3.4 mm of slurry (sub-plan 21 Task 8)"
```

---

## The outline (every item written out above)

6. **Sideways exchange**, split in three (Asha, 2026-10-21) and a fourth (2026-10-22), each with its own tests and
   commit; 6a, the transport at given velocities, 6b, the thickness terms implicit, 6c, the flow in the skin body, and
   6d, making it affordable, are written out above as Tasks 6a, 6b, 6c and 6d.
   - **6d. Making the flow affordable** (Asha, 2026-10-22, option a, before Tasks 7 and 8) is written out above as Task
     6d.
7. **The droplet momentum sink** placeholder at the top cells for stage 4 is written out above as Task 7.
8. **Verification** — level 2's plate with the flow on and a push along the face, the flow's cost, and the share of the
   flight's flow at short wavelengths — is written out above as Task 8. Whether that share is significant, and so
   whether a further task replaces the depth-uniform pressure by a non-hydrostatic pressure solved on the cells (Asha,
   2026-10-22), is decided once sub-plan 22 has measured it again on the flow-on flight's own states (Asha, 2026-10-25,
   option a; Review Focus 23).
9. **The base where runoff lands** (Asha, 2026-10-25, option a, over leaving it to sub-plan 22): on Task 8's coupled
   pushed plate, where runoff deposits on cold metal and is handed down, the base's reconciliation leaves faulted bottom
   cells up to 1 253 K and the part cap rejects parts down to 0.0104 s (Review Focus 24 and 25). Task 9 finds the cause
   and fixes it before any flow-on flight, and with it the heat the runoff loses to cold metal within the flow sub-steps
   (Review Focus 21) and the real pass's flow counts above the Courant target at a melting onset (Review Focus 17). Not
   yet written.

## Inputs that were to be settled (all settled)

- **The skins' interface** (sub-plan 20's Results; facts 118–125). Amended after the build (Asha, 2026-10-12): damped
  repetition, the profile valuation and the parts make the coupling converge and level 2 pass at 0.5 s, but a melting
  step on the 100 mm sphere costs up to 1 041 s (fact 125), and buried superheated cells remain at 6 MW/m² (fact 124).
  The onset excess at 0.1–0.25 s steps was cured by the melting part cap (fact 132); levels 2 and 2b now verify with the
  skins at the solid's density (`swell=False`), and compared like that the plates melt 2–6 % late on 2 mm elements, the
  3D model's resolution in depth (facts 133–137), now resolved by the fine band, 0.7 mm elements to 4 mm under the
  2 mm surface, which makes a melting step on the 100 mm sphere 45 % dearer (facts 138–140). To settle: whether flow's
  level-2 run uses the same setting, since with it the flowing layer's thickness is that of solid metal, not liquid.
  **Settled in Task 8 (item 5, ruling 5):** the regression against the references keeps the solid's density, the pushed
  plate runs with the skins' own law; at the solid's density the pushed plate's runoff is a third smaller.
- **The skins' 0.8 mm minimum** (option B′, Asha, 2026-10-16; facts 141–146): where the non-rigid layer is thinner than
  0.8 mm a skin's bottom is rigid metal, so the column flow must hold rigid cells still and count only the cells above
  T_rigid in the flowing depth; the finest cells stay in the top 0.4 mm (`SkinSettings.fine_zone`), and a column at its
  minimum holds 59 cells rather than 40, which the flow's per-cell work will feel. Flow makes each pass dearer,
  so the pass count is the first thing to settle: long parts during melting were rejected and needed many passes
  (fact 126), and parts of at most 0.0625 s while melting now cut a melting step to 122–226 s (fact 129). Absorption
  is limited only by the move allowance and the hot plate (level 2b) checks it; coarse surface meshes are not supported
  under strong melting, so flow's verification should use 2 mm surfaces (facts 130–131). **Settled:** rigid cells are
  held and only cells above T_rigid flow (Task 4's anchored cells, Task 6a's flowing depth); the pass count stays the
  part cap's with the flow in the trial pass (Task 6c's ruling; the real pass takes its trial pass's counts, 6d); the
  per-cell work was made affordable in 6d; Task 8's plates use 2 mm surfaces, and its pushed plate shows the parts
  rejected for their moves again where runoff deposits are handed down (Task 8, item 6).
- Sub-plan 20's measured run time (fact 120 under the one-repeat rule; fact 125 under the amended default: 60–93 s a
  macro step solid, 123–1 041 s melting on the 100 mm sphere, almost all in the skins), which was to decide whether the
  sideways solve can be implicit over every cell or needs sub-cycling. Decided (Asha, 2026-10-17): implicit, the
  amendment above, from the level-3 flight's profile (facts 147 and 148); design (a) (Asha, 2026-10-18; fact 150).
- Whether the implicit thickness pressures are enough to keep grid-scale growth out (level 1's growth test), given
  that the capillary length (about 1.8 mm) is close to the column spacing. **Settled:** in Task 6b, with surface
  tension in its thin-film form (Asha, 2026-10-21), the faces' exchange times the film operator has no growing mode and
  a mode alternating from column to column falls to at most 0.12 of itself in a sub-step of 10 ms to 1 s; on Task 8's
  pushed plate the grid scale's share of the face's thickness variance stays at 0.14-0.16 % while a deposit ridge grows.

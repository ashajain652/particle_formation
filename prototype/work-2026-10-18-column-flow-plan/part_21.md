# Sub-plan 21 — flow in the columns (stage 3 of the layered-skin design)

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or
> superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.
> Read `00-shared-context.md` first, then the spec. **Written at code level: Tasks 1–5, 6a, 6b, 6c, 6d, 7, 8 and 9
> (2026-10-18/27).** The outline points to the tasks, and the inputs it listed are settled, each where it says.
> Task 8's coupled pushed plate ran to 1.0 s, not 2.0 s (its ruling 8: about 6.5 hours alone to 2.0 s); Task 9 ran it
> again with its fixes, closing its books to 1.5 s and failing at 1.60 s where its deposits folded the base surface,
> and the pushed plate is judged at 1.5 s, the deposit's folding a watched risk for sub-plan 22 (Asha, 2026-10-27;
> Task 9, item 9).

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
    front (Review Focus 13), 13.7 of 338 mg landing 10-35 mm further; 11 % at 0.5 ms flow sub-steps. Since Task 9 the real pass's count
    follows its own flow where that is faster than its trial pass's (Review Focus 17): on the pushed plate's first
    melting macro step 9 real passes ran above a Courant number of 6.5, the largest at 9.2 (13 and 18.8 before), the
    heat sub-steps above it the first of a pass or one after the flow sped up (13 of 523).
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
    17 and 24. Task 9: that macro step took 3 878 s (one
    run beside it): no part rejected, but each 0.0625 s part still five damped passes of 45-77 s, the real passes now
    at the Courant target with three to five flow sub-steps a heat sub-step where Task 8's early parts ran one; the
    landing conduction 51 s of its 1 650 s of flow (3 %); the real pass whose cells flowed only because of its base
    took 54 flow sub-steps (24 s, 0.44 s each) where Task 8's replayed 380 (126 s).
16. **The transport's iterative solve** — stops at a normwise backward error of 1e-13 or after 40 corrections, then
    takes the LU; every cell is installed from the solve's fluxes, so the books do not depend on the tolerance and the
    masses stay positive at any (the iterates rise to the masses from below). Pinned by
    `test_the_iterative_transport_solve_gives_the_direct_solve_s_answer_where_flows_meet`,
    `test_the_books_close_exactly_and_the_masses_stay_positive_however_loosely_the_iteration_stops` and
    `test_an_iteration_that_cannot_finish_falls_back_to_the_direct_solve`. The 500-cell threshold suits the quiescent
    melt's many mid-sized cycles (0.88 s a transport against 1.45 s at 2 000 cells) and costs a driven state's lone
    1 700-cell cycle 0.1 s (Task 6d).
17. **The real pass's flow counts** — 6d took them from its trial pass, so that the base law is made by the pass it
    predicts; the skins' memory of the flow is part of their state (`Skins.state()` has seven arrays since 6d; code
    that unpacks six must take the first six). Task 8 found the trial pass, its base held at the start-of-step
    temperature, melting less and flowing slower at a melting onset, so that the real pass ran at the trial's counts
    above the Courant target (one flow sub-step a heat sub-step at 18.8 on the pushed plate's third melting part, where
    12.6 put the runoff 16 % off on 6d's sphere), and a real pass whose cells flowed only because of its base replayed
    the ten flow sub-steps its trial pass chose with no memory (380 in which nothing moved, 96 s of a 124 s pass). Since
    Task 9 (spec 16, 2026-10-26) the real pass takes in each heat sub-step the larger of its trial pass's count and the
    count its own memory chooses, and its own where the trial took none (0 in the plan); where the two passes flow
    alike nothing changes. On the pushed plate's first melting macro step 9 real passes ran above a Courant number of
    6.5, the largest at 9.2 (13 and 18.8 before), and the first macro step's real pass took 54 flow sub-steps (24 s)
    where it replayed 380 (126 s). Pinned by
    `test_a_real_pass_takes_the_larger_of_its_trial_pass_s_count_and_its_own_and_its_own_where_the_trial_flowed_nothing`
    and `test_a_skin_body_s_real_pass_takes_at_least_its_trial_pass_s_counts_and_its_own_where_the_trial_flowed_nothing`;
    6d's `test_the_real_pass_takes_the_trial_pass_s_counts_and_the_memory_is_part_of_the_skins_state` holds where the
    passes flow alike.
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
21. **Runoff over cold metal freezes within the flow sub-steps** (Task 9) — the flow sub-steps conducted no heat, so
    landed metal stayed liquid and coasted on until the next heat sub-step's conduction: on Task 8's pushed plate a film
    of 10-20 um at 0.6-1.2 m/s ran 6-12 mm in a 10 ms heat sub-step, where it freezes within 0.005-0.03 ms. Now the
    columns metal landed on top of have their top 0.4 mm conducted adiabatically in every flow sub-step from the landing
    on (`landing_conduct`): at 10 ms heat sub-steps the share of the pushed plate's deposit landing 5-10 mm past the heat
    falls from 9.9 % to 3.0 % (0.9 and 0.7 % at 5 and 2.5 ms), its centroid from 1.84 to 1.05 mm past the heat's end (0.68
    and 0.65 mm), at 2-3 % of the flow's cost. What remains at 10 ms is the flowing layer's own cooling at the heat's rim,
    which only the heat sub-step's conduction does (not measured further); the block's own cells, conducted with their
    column already, are conducted again within the block (ruling 5 of Task 9). Pinned by
    `test_metal_landing_on_cold_rigid_metal_loses_its_heat_to_it_within_the_flow_sub_steps` and
    `test_runoff_onto_cold_metal_freezes_where_it_lands_and_the_books_close`.
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
24. **Runoff deposits and the part cap** — Task 8: deposits handed down as they grew moved bases by 1.1-3.2 move
    allowances a 0.0625 s part (a rigid deposit grows as conduction into cold metal allows, 0.6-1.75 mm in 0.0625 s on
    850-735 K metal), so every melting part was rejected for its moves after its damped passes, down to 0.0104 s (fact
    126's cost, which the part cap removed): the allowance doing its job. Since Task 9 (spec 6.6, option C) metal more
    than 5 K warmer than the base node waits in the skin until it cools, so a deposit moves its base only as fast as its
    heat leaves it: on Task 8's pushed plate to 1.5 s no part was rejected (the largest move 0.8 allowances), and the
    first melting macro step took 3 878 s against 5 229 s. What it costs: the skins hold the warm deposit, 12.0 g in
    632 columns at 1.5 s, up to 4.05 mm, the plate's cells 2.5 % more (Task 9, item 6). But the limit releases a held
    cell whole once it is near its node's temperature, and the deep cells of a thick skin are as thick as the
    allowance: run on to 2.0 s, the plate split its parts in vain from 1.56 s (the same 1.026 allowances at 62.5 to 7.8
    ms), kept 8 of its last 15 over the allowance and handed the 12.4 g held down in 0.039 s before its mesh failed at
    1.60 s. Since 2026-10-27 (spec 6.6, 7.1) a part hands down at most what moves its base by the allowance, a cell
    going in part, so a hand-down never splits a part (Task 9, item 9). Run on with the cap from 1.5 s the plate split no part (parts of 0.05 s) and still failed at 1.80 s, on the deposit's
    folding (item 9). Pinned by
    `test_hand_down_passes_only_bottom_metal_within_the_fault_tolerance_of_the_base_node`,
    `test_hand_down_takes_at_most_its_cap_a_cell_going_in_part` and
    `test_a_deposit_thicker_than_the_move_allowance_comes_down_over_several_parts_and_splits_none`; the parts and the held
    mass are the logged runs' (`measure/task9/run_logged9.py`, `item4/run_fold.py`).
25. **Base faults where runoff deposits** — found in Task 9: sub-plan 20's moving mesh booked an outward sweep at the
    base node's enthalpy while the skins handed down warmer metal, and the reconciliation put the difference into the
    skin's new bottom cell, often a sliver of a cell (3 709 K before a merge, bottom cells kept up to 1 111 K on a reduced
    pushed plate, 1 253 K on Task 8's); Task 8's reading, the move's sweep less the skins' exchange beside columns
    drawing, had the place right and the quantity wrong. Since Task 9 (spec 6.6, options A and C) handed metal is booked
    at its own enthalpy and is within 5 K of its node: no bottom cell is heated by the reconciliation, and on Task 8's
    pushed plate to 1.5 s the hottest cell is 911.8 K (1 253 K before).
    The counter now judges the bottom cell after the reconciliation against the metal that crossed, for a hand-down
    the range from the handed cells' coldest and warmest to the cell it uncovered, as the hand-down left it (K2, in
    the wording Asha adopted on 2026-10-26): the skin's own warmer interior that a hand-down uncovers is no fault, a
    reconciliation that heats it is. Judged by the handed cells alone it had counted the uncovered cell wherever the
    skin's own gradient made that cell more than 5 K warmer than the warmest handed one, 1 900 of the 1 942 faults on
    the pushed plate to 1.5 s (Task 9, item 7). Pinned by
    `test_a_warm_deposit_goes_down_only_near_the_node_s_temperature_at_its_own_enthalpy_and_heats_no_cell`,
    `test_a_hand_down_s_bottom_cell_is_judged_against_the_handed_cells_and_the_cell_it_uncovered`,
    `test_the_cell_a_hand_down_uncovers_is_no_fault_unless_the_reconciliation_heats_it` and
    `test_an_outward_sweep_is_booked_at_the_enthalpy_the_caller_gives_for_it`.
26. **The metal the fault tolerance holds in the skins** (Task 9) — option C keeps rigid metal more than 5 K warmer
    than its base node in the skin until it has cooled, with no bound (Asha, 2026-10-26), and since 2026-10-27 a part
    hands down at most its move allowance, the rest waiting too: on Task 8's pushed plate it held 12.0 g in 632 columns
    at 1.5 s, up to 4.05 mm deep, rising at 11 g/s then and 7 g/s over the next part, to 12.4 g (Task 9, items 6 and 9) --
    over the 30 s or so of melting the 100 mm sphere sees, about 0.3 kg of its 1.47 kg if its runoff deposited as the
    plate's does (an extrapolation). The metal is still the body's and the tops are where they would be; what grows is
    the skins' depth and cells (the plate's cells 2.5 % more by 1.5 s). The skin body records it at every macro step
    (`held_mass_kg`, `held_depth_max_mm`; a run's results `held_mass_max_kg`, `held_depth_max_mm`), and **sub-plan
    22's flow-on flight reports it**, with the fold of the base surface the deposit makes (spec 12, 2026-10-27); a
    bound is set only if the skins grow too thick there. With the cap it came down at about 26 g/s, from 12.0 g at 1.5 s to 10.6 g at 1.75 s, where without it the whole
    12.4 g came down in 0.039 s. Pinned by
    `test_a_warm_deposit_goes_down_only_near_the_node_s_temperature_at_its_own_enthalpy_and_heats_no_cell` (the
    steep skins' whole excess held) and `test_hand_down_passes_only_bottom_metal_within_the_fault_tolerance_of_the_base_node`
    (`handable`, what a hand-down would take).

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

{{DIFF_TEST_SKIN_T1}}

- [ ] **Step 2: Run them to see them fail.** `"$PY" -m pytest tests/test_reentry_model_skin.py -q -k "sideways or
  coupled or warm_spot or explicit_limit"`. Expected: `AttributeError: module 'reentry_model.skin' has no attribute
  'Sideways'`.

- [ ] **Step 3: The edges' incidence matrix** in `reentry_model/columns.py`:

{{DIFF_COLUMNS_T1}}

- [ ] **Step 4: The coupled sub-step and the sideways map** in `reentry_model/skin.py`:

{{DIFF_SKIN_T1}}

- [ ] **Step 5: The skin body's hook** in `reentry_model/skinbody.py`:

{{DIFF_SKINBODY_T1}}

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

{{DIFF_TEST_SKIN_T2}}

- [ ] **Step 2: Append the surface operators' tests** to `tests/test_reentry_model_columns.py`:

{{DIFF_TEST_COLUMNS_T2}}

- [ ] **Step 3: Run them to see them fail.** `"$PY" -m pytest tests/test_reentry_model_skin.py
  tests/test_reentry_model_columns.py -q`. Expected: `AttributeError: 'Skins' object has no attribute 'U'` and
  `AttributeError: module 'reentry_model.columns' has no attribute 'triangle_geometry'`.

- [ ] **Step 4: The velocity in the cells** in `reentry_model/skin.py`:

{{DIFF_SKIN_T2}}

- [ ] **Step 5: The surface operators** in `reentry_model/columns.py`:

{{DIFF_COLUMNS_T2}}

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

{{TEST_RHEOLOGY_T3}}

- [ ] **Step 2: Append the drag tests** to `tests/test_reentry_model_skin.py`:

{{DIFF_TEST_SKIN_T3}}

- [ ] **Step 3: Run them to see them fail.** `"$PY" -m pytest tests/test_reentry_model_rheology.py
  tests/test_reentry_model_skin.py -q`. Expected: `ModuleNotFoundError: No module named 'reentry_model.rheology'`.

- [ ] **Step 4: Create `reentry_model/rheology.py`:**

{{RHEOLOGY_T3}}

- [ ] **Step 5: The cells' viscosity and the drag between them** in `reentry_model/skin.py`:

{{DIFF_SKIN_T3}}

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

{{DIFF_TEST_SKIN_T4}}

- [ ] **Step 2: Run them to see them fail.** `"$PY" -m pytest tests/test_reentry_model_skin.py -q`. Expected:
  `AttributeError: module 'reentry_model.skin' has no attribute 'momentum_step'`.

- [ ] **Step 3: The momentum sub-step, and the drag at each cell's own shear rate,** in `reentry_model/skin.py`:

{{DIFF_SKIN_T4}}

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

{{DIFF_TEST_SURFACE_FLOW_T5}}

- [ ] **Step 2: The push's tests** — append to `tests/test_reentry_model_skin.py`:

{{DIFF_TEST_SKIN_T5}}

- [ ] **Step 3: Run them to see them fail.** `"$PY" -m pytest tests/test_reentry_model_surface_flow.py
  tests/test_reentry_model_skin.py -q`. Expected: `AttributeError: 'SurfaceFlowResult' object has no attribute 'dp_ds'`
  and `AttributeError: module 'reentry_model.skin' has no attribute 'along_meridian'`.

- [ ] **Step 4: The pressure gradient on its own** in `reentry_model/surface_flow.py`:

{{DIFF_SURFACE_FLOW_T5}}

- [ ] **Step 5: The push** in `reentry_model/skin.py`:

{{DIFF_SKIN_T5}}

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

{{DIFF_TEST_COLUMNS_T6A}}

- [ ] **Step 2: The transport's tests** — append to `tests/test_reentry_model_skin.py`:

{{DIFF_TEST_SKIN_T6A}}

- [ ] **Step 3: Run them to see them fail.** `"$PY" -m pytest tests/test_reentry_model_columns.py
  tests/test_reentry_model_skin.py -q`. Expected: 11 failed — `AttributeError: 'Columns' object has no attribute
  'face_vector'` (5) and `AttributeError: module 'reentry_model.skin' has no attribute 'transport'` (6).

- [ ] **Step 4: The faces' vectors** in `reentry_model/columns.py`:

{{DIFF_COLUMNS_T6A}}

- [ ] **Step 5: The transport** in `reentry_model/skin.py`:

{{DIFF_SKIN_T6A}}

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

{{DIFF_TEST_COLUMNS_T6B}}

- [ ] **Step 2: The thickness terms' tests**, and Task 5's two capillary tests re-measured, in
  `tests/test_reentry_model_skin.py`:

{{DIFF_TEST_SKIN_T6B}}

- [ ] **Step 3: Run them to see them fail.** `"$PY" -m pytest tests/test_reentry_model_columns.py
  tests/test_reentry_model_skin.py -q`. Expected: 15 failed:
  - `AttributeError: 'Columns' object has no attribute 'film_curvature'` (2);
  - `AttributeError: 'Columns' object has no attribute 'film_operator'` (1);
  - `AttributeError: module 'reentry_model.skin' has no attribute 'flowing_tops'` (9);
  - `TypeError: momentum_step() got an unexpected keyword argument 'response'` (1);
  - the crust's books, `assert ... 1.32e-05 < 1.5e-12` (1);
  - Task 5's sphere test, `assert (... 10.99 < 0.7)` (1).

- [ ] **Step 4: The thin film's curvature** in `reentry_model/columns.py`:

{{DIFF_COLUMNS_T6B}}

- [ ] **Step 5: The momentum's response, each crust one unknown, the push's thin-film capillary term, the faces' push in
  the transport and the thickness step** in `reentry_model/skin.py`:

{{DIFF_SKIN_T6B}}

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

{{DIFF_TEST_SKIN_T6C}}

- [ ] **Step 2: The skin body's tests** — append to `tests/test_reentry_model_skinbody.py`:

{{DIFF_TEST_SKINBODY_T6C}}

- [ ] **Step 3: The command line's tests** — in `tests/test_reentry_model_cli.py` (the skin run's settings gain the
  flow's two keys):

{{DIFF_TEST_CLI_T6C}}

- [ ] **Step 4: Run them to see them fail.** `"$PY" -m pytest tests/test_reentry_model_skin.py
  tests/test_reentry_model_skinbody.py tests/test_reentry_model_cli.py -q`. Expected: 11 failed, 143 passed —
  - `AttributeError: module 'reentry_model.skin' has no attribute 'Flow'` (4) and `... 'flow_substeps'` (1);
  - `TypeError: SkinBody.__init__() got an unexpected keyword argument 'flow'` (3);
  - `TypeError: skin_name_suffix() got an unexpected keyword argument 'flow'` (1);
  - `SystemExit: 2`, `unrecognized arguments: --skin-flow on --skin-flow-cycles 2` (1);
  - `AttributeError: module 'reentry_model.skin' has no attribute 'FLOW_CYCLES'`, the skin run's settings (1).

  The seven exit-code cases pass already, their flags being unknown to the parser.

- [ ] **Step 5: The flow sub-steps and the pass's books** in `reentry_model/skin.py`:

{{DIFF_SKIN_T6C}}

- [ ] **Step 6: The flow in the skin body** in `reentry_model/skinbody.py`:

{{DIFF_SKINBODY_T6C}}

- [ ] **Step 7: The history's flow columns, the results and the frames' top speed** in `reentry_model/coupled.py`:

{{DIFF_COUPLED_T6C}}

- [ ] **Step 8: The settings** in `reentry_model/cli.py`:

{{DIFF_CLI_T6C}}

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

{{DIFF_TEST_SKIN_T6D}}

- [ ] **Step 2: The skin body's tests** — in `tests/test_reentry_model_skinbody.py`:

{{DIFF_TEST_SKINBODY_T6D}}

- [ ] **Step 3: The command line's tests** — in `tests/test_reentry_model_cli.py`:

{{DIFF_TEST_CLI_T6D}}

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

{{DIFF_SKIN_T6D}}

- [ ] **Step 6: The chosen count in the skin body** in `reentry_model/skinbody.py`:

{{DIFF_SKINBODY_T6D}}

- [ ] **Step 7: The history's count** in `reentry_model/coupled.py`:

{{DIFF_COUPLED_T6D}}

- [ ] **Step 8: The count's setting** in `reentry_model/cli.py`:

{{DIFF_CLI_T6D}}

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

{{DIFF_TEST_SKIN_T7}}

- [ ] **Step 2: The skin body's test** — append to `tests/test_reentry_model_skinbody.py`:

{{DIFF_TEST_SKINBODY_T7}}

- [ ] **Step 3: Run them to see them fail.** `"$PY" -m pytest tests/test_reentry_model_skin.py
  tests/test_reentry_model_skinbody.py -q -k "release or droplet or nothing_it_may_take or empties"`. Expected: 8
  failed, 99 deselected --
  - `AttributeError: module 'reentry_model.skin' has no attribute 'release_top'` (4: the take, the books, the emptied
    columns and the trial pass);
  - `AttributeError: module 'reentry_model.skin' has no attribute 'prescribed_release'` (4: the flow sub-steps, no
    release, Pi and the skin body).

- [ ] **Step 4: The release from the top cells, its books in the flow and Pi** in `reentry_model/skin.py`:

{{DIFF_SKIN_T7}}

- [ ] **Step 5: The release in the skin body** in `reentry_model/skinbody.py`:

{{DIFF_SKINBODY_T7}}

- [ ] **Step 6: The history's columns and the results** in `reentry_model/coupled.py`:

{{DIFF_COUPLED_T7}}

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

{{DIFF_TEST_SKINBODY_T8}}

- [ ] **Step 2: The short-wave measure's tests** -- create `tests/test_flow_wavelengths.py`:

{{TEST_FLOW_WAVELENGTHS_T8}}

- [ ] **Step 3: Run them to see them fail.** `"$PY" -m pytest tests/test_reentry_model_skinbody.py
  tests/test_flow_wavelengths.py -q -k "plate and not hot_plate or wavelength or overstatement or sinusoid or two_waves"`.
  Expected: `ModuleNotFoundError: No module named 'flow_wavelengths'` collecting `tests/test_flow_wavelengths.py`; and,
  `tests/test_reentry_model_skinbody.py` alone, 3 failed, 20 deselected --
  - `AttributeError: module 'skin_plate' has no attribute 'PUSHES'` (2: the flow-off plate's refusals, after its
    comparison with the hand-built body has passed on sub-plan 20's script, and the pushed plate);
  - `TypeError: run() got an unexpected keyword argument 'flow'` (the plate with nothing flowing).

- [ ] **Step 4: The flow on the plate** in `analysis/skin_plate.py`:

{{DIFF_SKIN_PLATE_T8}}

- [ ] **Step 5: The one-dimensional reference with the liquid leaving at the liquidus** in `analysis/ablation_1d.py`:

{{DIFF_ABLATION_1D_T8}}

- [ ] **Step 6: The short-wave measure** -- create `analysis/flow_wavelengths.py`:

{{FLOW_WAVELENGTHS_T8}}

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

### Task 9: The base where runoff lands — handed-down metal at its own enthalpy and only near the node's temperature, the runoff's heat lost within the flow sub-steps, and the real pass's flow counts at a melting onset

**Files:**
- Modify: `reentry_model/movemesh.py` (`MovingMesh.move(..., h_outward=None)` and `.move_with_policy(...,
  h_outward=None)`, an outward sweep booked at the enthalpy given; the module docstring), `reentry_model/skin.py`
  (`hand_down(..., T_limit=None, span=False, cap=None)` and `handable`, which share `_bottom_run`; `TransportBook.landed`;
  `LANDING_DEPTH`, `landing_conduct`; `flow_substeps` conducting what landed after each transport;
  `FLOW_SECONDS["landing"]`; `advance`'s real pass taking the larger of its trial pass's count and its own),
  `reentry_model/skinbody.py` (`jump_faults(..., handed=None)`; `HAND_DOWN_SHARE` and `SkinBody._hand_cap`, a part's
  hand-down capped at its move allowance (item 9); `SkinBody._move_ratio` and `._move_base`, which records the metal the
  limit and the cap hold, `held_mass` and `held_depth_max`, in `skin_stats`; the module docstring),
  `reentry_model/coupled.py` (two history columns and two results for the held metal)
- Test: `tests/test_reentry_model_movemesh.py` (one test), `tests/test_reentry_model_skinbody.py` (seven tests, two of
  them the cap's; 6d's counts test rewritten for the new rule; the parts test's recorded move, Task 8's pushed plate's
  threshold and the flow-off plate's recorded hand-down changed), `tests/test_reentry_model_skin.py` (three tests; the two flow sub-step tests built by hand take the landing
  conduction), `tests/test_reentry_model_cli.py` (a skin run's results carry the held metal)

**Interfaces:**
- Consumes: sub-plan 20's `MovingMesh.move_with_policy(u, h_inward)`, `skin.hand_down`, `SkinBody._move_base`,
  `_move_ratio`, `jump_faults`, `JUMP_TOLERANCE`; Task 1's `_conduct`; 6a's and 6d's `transport` (its new top cells where
  metal lands on a rigid top); 6d's `flow_cycles`, `Flow.plan`, `advance`'s trial and real passes; Task 7's release
  after the transport.
- Produces: `MovingMesh.move(u, h_inward=None, h_outward=None)` and `.move_with_policy(...)`, `h_outward` (n_corners,)
  [J/kg], nan where none is given; `skin.hand_down(sk, m_min, T_rigid, T_limit=None, span=False, cap=None) -> (mass,
  enthalpy)` or, with `span`, `(mass, enthalpy, coldest, warmest)` [kg, J, K, K per column; nan where nothing went],
  `T_limit` [K per column], `cap` [kg per column], at most that much taken, a cell going in part;
  `skinbody.HAND_DOWN_SHARE` (1 - 1e-6 of the move allowance) and `SkinBody._hand_cap() -> mass` [kg per column]; `skin.handable(sk, m_min, T_rigid, T_limit=None) -> mass` [kg per column], what `hand_down` would take,
  taking nothing; `skinbody.jump_faults(T_bottom, T_node, T_deep, crossed, tol=None, handed=None)`, `handed` = (mask,
  coldest, warmest[, uncovered]), the last the bottom cell's temperature as the hand-down left it, `crossed` now the
  columns that drew or absorbed; `SkinBody.held_mass` [kg] and `.held_depth_max` [m], what the limit and the cap hold
  at the last base move, as the history's `held_mass_kg` and `held_depth_max_mm` and the results' `held_mass_max_kg` and
  `held_depth_max_mm`; `TransportBook.landed` (n,) bool;
  `skin.landing_conduct(sk, landed, dt, depth=None)`, `skin.LANDING_DEPTH` (0.4 mm); `FLOW_SECONDS["landing"]`;
  `Flow.plan` holding 0 where the trial pass took no flow sub-step. Nothing else changes in an interface: the command
  line is as it was, and the history gains the two columns of the held metal, every other keeping its name
  (`jump_faults` now counts by spec 6.6 as amended).

Task 8's coupled pushed plate closed its books, but where runoff left the heat and landed on cold rigid metal its first
melting macro step counted 1 761 base faults, left bottom cells up to 1 253 K -- spec 12's buried superheated cells,
which sub-plan 20 had removed -- and had every part from the fifth rejected for its moves, down to 0.0104 s (5 229 s for
the macro step); the flow sub-steps conducted no heat, so landed metal ran on 6-12 mm before freezing; and a real pass
took its trial pass's flow counts even where it flowed faster or where its trial pass had flowed nothing (Review Focus
17, 21, 24, 25). Asha, 2026-10-25 (option a): find the cause before any flow-on flight and fix it, with the other two.
Investigated and measured 2026-10-25/26 in `prototype/work-2026-10-18-column-flow/measure/task9/`; the cause was
reported first, and Asha decided the fix on 2026-10-26 (spec 6.6 and 16 as amended that day): options A and C and
the counter K2 below, and the counts rule of item 4. The pushed plate's run to 2.0 s then failed at 1.60 s; its cause,
reported in a third hand-back (item 9), led Asha on 2026-10-27 to cap a part's hand-down at its move allowance and to
judge the pushed plate at 1.5 s, the deposit's folding a watched risk for sub-plan 22's flight (option 1, then option
2; spec 6.6, 7.1, 12 and 16 as amended that day).

**1. The cause, traced** (`measure/task9/trace_move.py`, `analyse_trace.py`; ruling 1). *The reproduction:* Task 8's
pushed plate reduced to 60 mm across (2 282 columns against 15 812), heated within 15 mm and to nothing at 21 mm, from
700 K under 6 MW/m2, the flight push, the skins' own density law, 0.5 s macro steps, no fine band; every kept part's
base operations recorded column by column (262 s of processor time to 1.0 s, alone). It gives faults of the same kind
and size: every fault a column handing metal down on the face just downstream of the heat (y = -12 to -22 mm), 199 by
1.0 s (Task 8: 1 761 on seven times the columns); bottom cells after a move up to 266 K above their node, 1 111 K
(Task 8: 5-400 K, 1 253 K); and parts rejected for hand-downs at the heat's downstream rim, split down to 0.0156 s
(Task 8: 0.0104 s). Task 8's state at 1.0 s fits: its 90 hot bottom cells are all on the face, the cell above each a
median 10 K over its node, the bottom cells 0.5-1.4 times their target, as a hand-down and a regrid leave them.

*One faulted column through a part* (column 201, the part ending at 0.8125 s): through the real pass's conduction, flow
sub-steps, arrivals, release, draws and regrids its bottom cell stays within 1-5 K of the node (811.8 K over 810.8 K);
the hand-down then takes 0.316 mm of its bottom, whose mean is 818.0 K -- the skin above its base is 35-67 K/mm warmer
than the base node where a deposit's heat is conducting down, against the 3D model's 13-28 K/mm below it -- and leaves
of the cell it took in part a sliver of 8.5e-9 kg (1.4 % of a cell) at 824.3 K; the moving mesh books the outward
sweep at the node's enthalpy, h(810.8 K), extrapolates its own gradient over the new metal and closes its books by its
uniform offset; the reconciliation (`exchange_bottom`) gives the skin the sweep less what it exchanged -- in mass
2.0e-9 kg, but in enthalpy 0.038 J beyond the node's value, the handed slab's excess over the node -- and the sliver
reaches 3 709 K; the regrid merges it, 859.6 K, 49 K above the node. Where the cell the hand-down left is too big to
merge (60 % of its target) the heat stays: 1 111 K. Over each part the enthalpy the reconciliation returned beyond the
node's equals the handed metal's excess over its node to 0.1 % (0.16-1.9 J a part, against 0.0005-0.005 J before any
deposit formed); the faulted columns drew nothing, and the sweep's mass matched the hand-down to within 2 % of a cell.
*The cause:* sub-plan 20's moving mesh booked an outward sweep at the corner's own enthalpy while the skins hand down
their bottom cells at theirs, which spec 6.6 says the 3D model gains; harmless while hand-downs were thin and at the
node's temperature, but runoff thickens the skins by 0.2-1 mm a part, so the slab handed spans their warmer interior,
and the difference landed in one bottom cell, often a sliver -- spec 12's mechanism of "metal crossing the base at
different temperatures", returned through hand-downs. Task 8's reading (Review Focus 25, "the move's sweep less the
skins' exchange into small bottom cells of columns handing deposits down beside columns drawing") had the place right
and the quantity wrong: it is the valuation, not the mass, and the neighbours' draws play no part. The fault counter
fired one step earlier for another reason: it judged a hand-down against the node alone, so it counted the warmer cell
a hand-down uncovers (6-30 K above the node) rather than the superheating, which it never saw (199 faults, against 270
bottom cells the reconciliation heated by more than 5 K). Confirmed by
`test_a_warm_deposit_goes_down_only_near_the_node_s_temperature_at_its_own_enthalpy_and_heats_no_cell`, which fails
before the fix with a cell at 1 614 K in a column whose hottest metal was 890 K.

*One rejected part* (0.9375-1.0 s, 0.0625 s): 12 bases beyond the allowance (0.538 mm), the worst column 141, 1.2 mm
past the heat's end, by its hand-down of 1.055 mm (1.96 allowances); the largest draw anywhere 0.19 mm. Every rejection
on the reduced plate was a hand-down at the heat's downstream rim (y = -20 to -22 mm).

*Are the rejections the allowance doing its job?* Yes. A rigid deposit freezing onto cold metal under a sustained
supply grows as conduction into the metal allows: 2 e (T_rigid - T_s) sqrt(t / pi) / (rho dh), e = sqrt(k rho c) the
metal's effusivity and dh = h(908 K) - h(T_rigid) = 2.1e5 J/kg the heat that makes liquid rigid -- 1.75, 1.16 and 0.61
mm in 0.0625 s on metal at 735, 800 and 850 K (a mean of 28, 18.5 and 9.7 mm/s). The runoff supplies about 5-9e-5
m2/s per unit width (Task 8's 4.37 g in 0.2 s over the downstream rim's 90-140 mm), 28-43 mm/s into one 2 mm column,
so the first cold columns grow at the conduction limit -- the reduced plate's front column 1.06 mm in 0.0625 s on
810-860 K metal -- and the rest spreads over a few columns. Handed down as it grows, a deposit moves its base out by
1.1-3.2 quarter-edges a 0.0625 s part, which needs parts of 0.02-0.06 s; it was the allowance doing its job, at ten
times a part's cost (fact 126's waste, each part converging 5-7 damped passes before being rejected). Freezing within
the flow sub-steps (item 3) pulls back the deposit's thin leading edge, not its front column, already at the
conduction limit, so it changes the part length little (estimated). Held in the skin until it cools (option C), the
deposit no longer moves the base faster than its heat leaves it.

**2. The fix** (Asha, 2026-10-26; spec 6.6 as amended). *Option A:* the moving mesh books an outward sweep at the
handed metal's own enthalpy, per corner the mean of what was handed there (`MovingMesh.move(..., h_outward)`, nan where
none was, those corners' sweeps at their own enthalpy as before); the 3D model gains exactly what the skins handed, and
the reconciliation returns only the sweep's mass mismatch, at that enthalpy. *Option C:* `hand_down(..., T_limit)` passes
down only rigid bottom cells no warmer than the base node plus the fault tolerance (`T_node + JUMP_TOLERANCE`, the node
after the part's 3D solve; ruling 2); warmer metal waits in the skin until it has cooled. Where the limit ends the take
at a cell's top, the round-off that the take leaves of that cell (a remainder below SLIVER of it) goes too (ruling 3):
left, a remainder of 1e-22 kg took the reconciliation's round-off to thousands of kelvin and the counter read it.
*Counter K2:* `jump_faults(..., handed=(mask, coldest, warmest, uncovered))` judges the bottom cell, now read after
the reconciliation (before the regrid), against the metal that crossed -- for a draw or an absorption the node's and
the profile's temperatures at the depth swept, as before, for a hand-down the range from the handed cells' coldest and
warmest (`hand_down(..., span=True)`) to the cell it uncovers, as the hand-down left it (Asha, 2026-10-26, adopting
item 7's wording after the build); a column that did both is a fault only outside both. On the reduced plate to 1.0 s
(prototypes in `measure/task9/proto_A`, `proto_C`, `proto_AC`; the final code in `count_faults.py`):

| | before | A alone | C alone | A and C, as built (with items 3 and 4) |
|---|---|---|---|---|
| bottom cells the reconciliation heated by more than 5 K | 270 | 0 | 141 | 0 |
| largest excess of a bottom cell over its node after a move | 266 K | 32 K | 44 K | 11 K |
| parts attempted, rejected for their moves | 17, 4 (to 0.0156 s) | 18, 4 | 9, 0 | 9, 0 |
| processor time to 1.0 s (alone) | 262 s | 239 s | 173 s | 169 s |

The table shows A removing the superheating and C the rejections, and only both together removing both; the counter
as amended then counts 17 faults where it counted 199 (item 7). Mass, energy and momentum close to 9e-16, 5e-11 and
2e-17 in every variant.

**3. The runoff's heat lost within the flow sub-steps** (`skin.landing_conduct`; ruling 5). The flow sub-steps conducted
no heat, so metal landing on a rigid top stayed liquid and coasted on (6d's rule) until the next heat sub-step's
conduction caught it: on Task 8's pushed plate a film of 10-20 um at 0.6-1.2 m/s ran 6-12 mm over metal at about 700 K
in a 10 ms heat sub-step, where it freezes in 0.005-0.03 ms (a contact estimate: the heat that makes 10-20 um of liquid
at 908 K rigid, rho_l d dh with dh = 2.1e5 J/kg, against 2 e (T_c - T_s) sqrt(t / pi) into metal at T_s = 735-850 K
from the contact temperature T_c, their mean), 3-36 um of travel. Now the transport marks the columns metal landed on
top of (`TransportBook.landed`, its new top cells), and from the flow sub-step it landed in to the heat sub-step's end
each such column's top block -- the cells within LANDING_DEPTH, 0.4 mm, of its top, the landed metal and the cold metal
beneath -- is conducted for the flow sub-step with no heat across the block's top or bottom (`_conduct` on the blocks
alone), after the transport and before the release. The block is adiabatic, so the books stay exact (the Newton
iteration's error, about 1e-10 J a cell, kept by its bottom cell) and the base law and its sensitivity are untouched;
the film freezes where it lands, coasting at most a flow sub-step, 0.3-0.6 mm at 0.5 ms and 1.2-2.4 mm at the 2 ms the
pushed plate's count chooses. Measured: a 20 um film at 950 K on 0.8 mm of metal at 700 K, one 0.5 ms flow sub-step,
ends at 763-765 K, the block's energy to 1.1e-16; on a strip, liquid 0.2 mm deep sheared onto metal at 700 K for one
10 ms heat sub-step, the heat sub-step's conduction off, is rigid beyond the edge's first column, the deposit's centroid
0.30 columns past the edge against 1.63 without (the unit tests, 2026-10-26).

*Where the pushed plate's deposit lands against the heat sub-step* (Task 8's study, item 4: 50 ms from its pushed skins
at 0.2 s, the conduction on, the base held; `measure/task9/plate_skins9.py`, `compare_substeps9.py`): of the metal
deposited beyond the heat within 20 mm of the centre line, the share landing 5-10 mm past the heat's end, the deposit's
centroid past it, the tops beyond the heat still above T_rigid at the end, and the processor time (under a load of two
or three runs):

| Heat sub-step, flow sub-steps | 5-10 mm past, without; with | Centroid past, without; with | Tops above T_rigid | Processor time (landing) |
|---|---|---|---|---|
| 10 ms, 0.5 ms (20) | 39.0 mg, 9.9 %; 11.4 mg, 3.0 % | 1.84; 1.05 mm | 71; 45 | 71; 86 s (1.1 s) |
| 5 ms, 0.5 ms (10) | 2.1 mg, 0.6 %; 3.1 mg, 0.9 % | 0.69; 0.68 mm | 51; 40 | 89; 118 s (1.1 s) |
| 2.5 ms, 0.5 ms (5) | 1.8 mg, 0.5 %; 2.4 mg, 0.7 % | 0.67; 0.65 mm | 48; 37 | 108; 151 s (1.0 s) |
| 10 ms, chosen (5, 2 ms) | 44.0 mg, 11.2 %; 21.1 mg, 5.8 % | 1.99; 1.34 mm | 74; 34 | 45; 50 s (0.7 s) |

The table shows the landing conduction pulling the deposit back toward where the shorter heat sub-steps put it: at 10
ms the share landing 5-10 mm past the heat falls from 9.9 % to 3.0 % (0.9 and 0.7 % at 5 and 2.5 ms) and its centroid
from 1.84 to 1.05 mm past the heat's end (0.68 and 0.65 mm), half the tops beyond the heat that were still above
T_rigid are rigid, and at 5 and 2.5 ms, which agreed before, it changes the deposit by under 0.1 mm. What remains at 10
ms is not landed metal: the liquid still flowing over rigid metal at the heat's rim cools only in each heat sub-step's
conduction (not measured further). With the conduction off, Task 9's code gives Task 8's 10 ms run to the last bit.
At the depths of ruling 5's alternatives, 0.2 and 0.8 mm, 3.3 and 2.7 % land 5-10 mm past, the centroid at 1.11 and
0.97 mm. The landing conduction costs 1.0-1.1 s of a flow's 50-77 s (2 %) at 0.5 ms flow sub-steps and 0.7 s of 23 s
(3 %) at the chosen count.

**4. The real pass's flow counts** (spec 16, 2026-10-26, amending 6d's "the real pass taking its trial pass's counts";
ruling 6). The trial pass holds its base at the start-of-step temperature, so at a melting onset it can melt less and
flow slower than its real pass, which then ran at the trial's counts: on Task 8's pushed plate one flow sub-step a heat
sub-step at a Courant number of 18.8 (target 6), and a real pass whose cells flowed only because its base was held at
the 3D model's first answer replayed the ten flow sub-steps its trial pass chose with no memory, 380 in which nothing
moved, 96 s of a 124 s pass. Now the trial pass leaves in its plan the count it took in each heat sub-step, 0 where
nothing flowed (`flow_substeps` took no sub-step), and the real pass takes the larger of that and the count its own
memory chooses, and its own where the plan holds 0: where the two passes flow alike the counts are the same and the base
law is made by the pass it predicts (6d's reason), and a real pass flowing faster follows its own flow from its next
heat sub-step, as any pass's count follows a flow that speeds up. On the pushed plate's first melting macro step (item
5) 9 real passes ran above a Courant number of 6.5, the largest at 9.2 (Task 8: 13, and 18.8); of their 523 heat
sub-steps 13 were above 6.5, 5 of them the first of a pass. In its first macro step the real pass whose cells flowed
only because of its base took 54 flow sub-steps (24 s of processor time) where Task 8's replayed 380 (126 s): ten in
the heat sub-step its cells began to flow, as a flow from rest takes, then 5, 3, 2 and one a heat sub-step, at Courant
numbers of 0.001-0.004. On the reduced pushed plate to 1.0 s, 10 of 71 real passes ran above 6.5 before, the largest
at 24.0, and 6 of 38 after, the largest at 9.2 (`measure/task9/runs/red60_*.passes.jsonl`).

**5. Task 8's coupled pushed plate again** (its settings: cold, 6 MW/m2, the flight push, the skins' own density law,
170 mm across, 0.5 s macro steps; `measure/task9/run_logged9.py`, every pass and part logged; `analyse_push9.py`). It
ran detached on mains power from 14:15, alone but for the unit tests of this task for most of its first melting macro
step (one process beside it, so its times are under a light load), then beside the other runs of this task. Its first
melting macro step (0.5-1.0 s), against Task 8's run of it:

| 0.5-1.0 s | Task 8 | Task 9 |
|---|---|---|
| parts of the macro step; rejected for their moves | every part from the fifth rejected, run again down to 0.0104 s | 8 of 0.0625 s; none (the largest move 0.53 allowances) |
| skin passes (real); interface repeats | 303 (152); 107 | 75 (38); 29 |
| processor time in the skins' passes; the flow's; the landing conduction's | 4 209 s; 1 888 s; -- | 3 215 s; 1 650 s; 51 s |
| wall time of the macro step | 5 229 s (alone) | 3 878 s (one run beside it) |
| real passes above a Courant number of 6.5; the largest | 13; 18.8 | 9; 9.2 |
| base faults counted | 1 761 (the old counter) | 456 (K2, every one the cell a hand-down uncovered) |
| bottom cells' largest excess over their node after a move; the hottest cell | 400 K; 1 253 K | 20.2 K; 910.8 K |
| metal handed down; held in the skins by the limit at 1.0 s | 10.5 g; -- | 5.4 g; 4.7 g in 381 columns, 2.0 mm at most |
| cells at 1.0 s | 523 355 | 527 340 |

The table shows the fixes doing what they were for: no part rejected, no bottom cell heated beyond its neighbours (the
hottest cell 910.8 K, liquid at the top, against 1 253 K at a bottom), the real passes at their target. The macro step
is not the ten times cheaper that removing the rejections promised: each 0.0625 s part still takes four or five damped
passes (4, 4, 4 and 9 in Task 8's first parts), and its passes cost 19-77 s each against 7.5-40 s, because the real
passes now run the flow at the Courant target -- three to five flow sub-steps a heat sub-step where Task 8's early parts
ran one at a Courant number of 18.8 -- and the skins hold the deposit (4 000 cells more). In the first macro step the
real pass whose cells flowed only because of its base took 54 flow sub-steps (24 s) where Task 8's replayed 380 (126 s).
*From 1.0 to 1.5 s* (beside the other runs of this task; Task 8 did not run it): 7 118 s of wall time, eight parts of
0.0625 s, none rejected (the largest move 0.8 allowances), four damped passes each, 64 passes, 6 339 s of processor time
in them (the flow 3 797 s, the landing conduction 77 s), the real passes at six or seven flow sub-steps a heat sub-step
(three to five in the macro step before) and a Courant number of at most 6.1, 1 486 faults counted (1 444 the uncovered
cell), 10.0 g handed down; the centre receded to 2.756 mm. In the last part (1.4375-1.5 s) the bottom cells' largest
excess over their node rose to 68.6 K, the counter counted 33 faults that were not uncovered cells, the sharpest fold of
the base surface went from 85.7 degrees (1.0 s) to 95.3 (98.0 at 0.5 s, the plate's edges; spec 12's fold limit is 75-95
degrees), the third rebuild had taken place and a column held 6.1 mm of metal above T_rigid: the runoff had spilled
over the plate's square edge, examined with the run's continuation to 2.0 s in item 9. Books: the
history prints the mass error, the energy residual and the momentum residual as 0.0 to its six decimals at 1.0 and 1.5 s
(the run's summary, with them in full, is written only at its end); the reduced plate's are 6.7e-16, 3.4e-11 and
1.5e-17. At 1.0 s the centre has receded 0.852 mm (Task 8 0.795 mm; the flow-off plate of that width 0.857 mm): the
limit's doing (on the reduced pushed plate the centre recedes 0.704 mm before this task, 0.729 mm after it, 0.700 mm
with the limit undone and 0.729 mm with any other change undone), though at 1.0 s the limit holds nothing in the fully
heated region's skins (1.15 g in the taper's 153 columns, 3.25 g in 208 beyond the heat): it acts through the heat's rim
and the flow there, not traced further.

**6. The metal option C holds in the skins** (`run_logged9.py`: after every kept part, what a hand-down with no limit
would take less what went; ruling 7). The hold is unbounded (Asha, 2026-10-26): the skin body records it at every base
move (`held_mass`, `held_depth_max`; the history's `held_mass_kg` and `held_depth_max_mm`, a run's results' largest),
by `handable` with no limit less what the hand-down took, and sub-plan 22's flow-on flight reports it; a bound is set
only if the skins grow too thick there. On the reduced pushed plate to 1.0 s the body's record equals the copy's
measure to the last bit in every part (0.31 g in 52 columns at 1.0 s, the deepest 1.37 mm; `runs/red60_k2`). On the
pushed plate it grows as long as runoff keeps arriving, and had not levelled off by 1.5 s:

| End of part | 0.75 s | 0.875 s | 1.0 s | 1.125 s | 1.25 s | 1.375 s | 1.5 s | 1.5625 s | 1.6016 s |
|---|---|---|---|---|---|---|---|---|---|
| held [g] | 0.48 | 2.34 | 4.68 | 6.95 | 9.13 | 10.83 | 11.96 | 12.38 | 2.14 |
| columns holding; the deepest held [mm] | 177; 0.67 | 286; 1.56 | 381; 2.01 | 462; 2.36 | 542; 2.70 | 582; 2.84 | 632; 4.05 | 651; 4.04 | 316; 2.35 |
| cells (523 777 before any melting) | 523 382 | 527 019 | 527 340 | 530 182 | 531 814 | 534 653 | 536 754 | 537 569 | 530 119 |

The table shows the held metal rising at about 17 g/s from 0.75 to 1.25 s and at 11 g/s over the next quarter-second,
while 10 g was handed down from 1.0 to 1.5 s (15.4 g in all by 1.5 s): it slows as the oldest deposits cool to their
node and go down. The last two columns are the run to 2.0 s (item 9; its earlier columns agree with these to a few
hundredths of a gram): 7 g/s over the next part, levelling off, and then the limit released it in whole cells, 10.2 g
in 0.039 s, until the mesh failed; with the cap (item 9) it fell from 12.0 g at 1.5 s to 10.6 g at 1.75 s, about 26 g/s, the deepest held 3.6-3.9 mm. Where it sits at 1.0 s: 1.15 g in 153 columns of the heat's taper, 3.25 g in
208 beyond the heat downstream, none in the fully heated region or upstream. A holding column's skin carries 49 cells
on average at 1.0 s (at most 62) against 33 at its minimum (development cells; graded cells grow 10 % a cell, so a
skin's cells grow with the logarithm of its thickness: with production cells a 4 mm skin holds about 76 against 59),
and the plate's cells grew 2.5 % by 1.5 s, which a skin sub-step's conduction, a third to a half of its cost, follows in
proportion -- about 1 % of a skin sub-step. *Over a flight* (estimated, not measured): extrapolated at the last
quarter-second's 11 g/s, about 0.3 kg would be held in the skins after the 30 s or so of melting the 100 mm sphere sees
-- a fifth of its 1.47 kg -- if its runoff deposited as this plate's does; how much reaches the cold rear depends on the
spraying (sub-plan 22). The metal is still the body's, in the skins rather than the 3D model, and the skins' tops are
where they would be; what grows is their depth and their cells, a few millimetres and tens of cells a deposit column
after a second.

**7. The counter as amended, and the cell a hand-down uncovers.** Read after the reconciliation and judged against the
handed cells for a hand-down, the counter counts what the reconciliation did wrong where the old counter counted the
skin's own interior: on the reduced pushed plate to 1.0 s before the fixes it would have counted the 270 bottom cells
the reconciliation heated by more than 5 K, where the old counter counted 199 uncovered cells. After the fixes, judged
by the handed cells alone, it counted 17 there (none to 0.75 s, then 2, 1, 4 and 10 a part), and on Task 8's pushed
plate 456 to 1.0 s and 1 942 to 1.5 s. None of the reduced plate's was a heated cell, nor were 1 900 of the pushed
plate's (the other 42 came in its last four parts): every one was the cell the hand-down uncovered, within 1 K of its
temperature before the reconciliation (the column `faults_uncovered` of `run_logged9.py`). The limit ends a hand-down
at the first cell more than 5 K warmer than the node, and where the skin's own gradient makes that cell, now the bottom
one, more than 5 K warmer than the warmest cell handed, judging by the handed cells alone counts it. Asha adopted on
2026-10-26 the wording that counts only the coupling's faults -- "... for a hand-down the handed cells' and the cell it
uncovers', as the hand-down left it" -- built as ruling 4 says, the uncovered cell read straight after the hand-down
and the bottom cell still after the reconciliation. On the reduced pushed plate to 1.0 s it counts none of the 17
(`runs/red60_k2`); `test_the_cell_a_hand_down_uncovers_is_no_fault_unless_the_reconciliation_heats_it` pins that the
cell a hand-down uncovers is no fault and that the same cell heated by a reconciliation returning 0.05 J more is. On
the pushed plate: item 5.

**8. Flow off and nothing landing.** With the flow off only option A moves a result, by booking hand-downs at their own
enthalpy; option C holds nothing back where the skins' bottoms sit at their nodes' temperature, which
`test_on_the_flow_off_plate_the_fault_tolerance_holds_no_hand_down_back` pins on the reduced flow-off plate, and the
level-2 plates with C alone are the same to the last bit, every array and history column, to 1.0 s (cold and hot,
seeded; `measure/task9/flowoff_state.py`, `TASK9_NO_A=1`). With A, at 1.0 s: the centre's recession moves by -1.1e-8
mm (cold) and +1.3e-8 mm (hot); 38 and 37 of the 6 932 columns change their mean specific enthalpy by more than 100
J/kg (at most 7 279 and 1 090 J/kg), their top cells by more than 0.1 K in 9 and 1 columns (at most 1.1 and 0.12 K);
on the cold plate 6 columns hold one cell more or fewer (a split or merge on the other side of its threshold), and the
one rebuild at 0.5-1.0 s meshes the volume differently (88 898 against 88 849 nodes; 86 686 against 86 664 hot), the
mesher answering round-off with a different mesh (fact 124). On the reduced flow-off plate of the unit tests (880 K,
6 MW/m2, 0.2 s), cells move by at most 0.048 K, the 3D field by 0.0042 K, the masses and the recession not at all. The
flow-off plates over 2.0 s still pass levels 2 and 2b, against Task 8's runs of them (unseeded both):

| Flow off | 1.0 s | 1.5 s | 2.0 s |
|---|---|---|---|
| cold, against the 1D reference; the fine box: Task 8 | +3.2 %; +1.2 % | +0.5 %; -0.2 % | -0.9 %; -0.6 % |
| cold: Task 9 | +3.2 %; +1.2 % | +0.4 %; -0.2 % | -0.8 %; -0.6 % |
| hot, against its 1D reference: Task 8 (flowing depth against the exact; base over T_rigid) | -9.0 % (onset) | -3.0 % (+1.1 %; 0.012 K) | -1.1 % (-0.1 %; 0.069 K) |
| hot: Task 9 | -9.0 % (onset) | -3.0 % (+1.1 %; 0.038 K) | -1.1 % (-0.4 %; 0.070 K) |

The table shows the flow-off plates meeting levels 2 and 2b as before: the cold plate's centre moves by -0.0008 mm at
1.5 s and +0.0038 mm at 2.0 s (of 2.72 and 5.17 mm), the hot plate's by -0.00002 and -0.0004 mm (of 1.76 and 2.88 mm),
no base fault, mass and energy to 2.2e-16 and 5.2e-13 (cold) and 4.5e-13 (hot). No test asserted a flow-off result of an
earlier build to the last bit; each such claim was a measurement in a docstring (Task 6c's
`test_the_flow_is_off_by_default...`, Task 8's `test_the_plate_with_the_flow_off_is_sub_plan_20_s_plate`), whose
assertions compare two bodies of the same code and stand. *With the flow on and nothing landing on cold metal*: the
landing conduction acts only where metal landed on a rigid top, so it changes nothing else; the real pass's counts
change only where its own memory asks for more flow sub-steps than its trial pass took or where its trial pass flowed
nothing; A changes every result in which metal is handed down (by the same order as with the flow off where the handed
metal is near its node's temperature, which C now ensures); C changes the results where a hand-down would have taken
metal more than 5 K warmer than its node -- deposits, and skins thickened by arriving melt over a steep profile; and the
counter counts differently wherever metal is handed down.

*The regression with the flow on* (Task 8's item 3: no push, the liquid released as it forms, the skins at the solid's
density, against the references whose liquid leaves at the liquidus):

| Flow on and released | 1.0 s | 1.5 s | 2.0 s | Hot: flowing depth, base over T_rigid |
|---|---|---|---|---|
| cold, 1D; fine box: Task 8 | +2.4 %; +0.6 % | +0.3 %; -0.3 % | -0.8 %; -0.6 % | |
| cold: Task 9 | +2.4 %; +0.6 % | +0.5 %; -0.1 % | -1.2 %; -1.0 % | |
| hot: Task 8 | -6.4 % (onset) | -1.6 % | -1.1 % | within 1.2 %; at most 0.10 K |
| hot: Task 9 | -6.7 % (onset) | -3.4 % | -1.5 % | within 0.9 %; at most 0.09 K |

The table shows both plates still meeting level 2 against the liquidus's references (5 % allowed), the hot plate 1.8
points further behind at 1.5 s. No single change of this task accounts for that: undoing each in turn on the hot plate
to 1.5 s, its centre recedes 1.935 mm in Task 8, 1.899 mm with all of them, and 1.922, 1.916 and 1.911 mm with the
booking (A), the limit (C) and the landing conduction undone, 1.898 mm with the counts rule undone
(`measure/task9/plate_toggle.py`): metal is handed down there (2.1-2.7 g), some of it warmer than its node, and the
rim's slurry lands on rigid tops, and each moves the coupled run by a fraction of a per cent. Books: mass to 2.2e-16
and 4.4e-16, energy to 5.2e-13 and 4.3e-13, momentum to 7.4e-18 and 4.9e-17 (cold, hot); no base fault; released by
2.0 s 49.47 g (cold, as Task 8's) and 32.72 g (hot; 32.70 g).

**9. The pushed plate to 2.0 s, and the cap on a part's hand-down** (2026-10-27; `measure/task9/item4/`, the third and
fourth hand-backs in `measure/task9/report.md`). Run again to 2.0 s on mains power (`run_logged9.py`, as item 5 but with
the end's examination logged), the plate closed its books to 1.5 s and failed inside the 1.5-2.0 s macro step: "a corner
move of up to 0.000399 m degrades the mesh even after a rebuild", the base surface's sharpest fold 158.7 degrees. *What
folds* (`run_fold.py`, `folds_state.py`): the deposit. Runoff freezing on cold metal is handed down, so the base grows out
under it, unevenly at the scale of one triangle -- a one-column ridge along the flow in the heating's taper (at a radius
of about 36-40 mm, corners 1-1.5 mm proud of their neighbours 2 mm away, the triangles either side tilted about 35
degrees), the heated face's own sharpest fold rising from 24 degrees at 1.0 s to 69 at 1.5 s, about 90 degrees a
second -- and, from about 1.3 s, the runoff spilling over the plate's square y = -85 mm edge, across the 8 mm side face
(skins up to 7.9 mm, 6.1 mm of it above T_rigid, at 1.5 s) to the back edge, where the surface starts at 90 degrees.
Spec 12's limit for neighbouring triangles is about 75-95 degrees. *What made it fast*: option C's limit releases a held
cell whole once it is within 5 K of its node, and the deep cells of a thick skin are large (a cell's target is a tenth
of its depth below the fine zone): at 1.5 s the lowest held cells were 0.18 mm of solid metal at the median and up to
0.79 mm (20 of 613 holding columns over the 0.556 mm allowance, 163 over half of it), a median 1.8 K above the limit. A
release does not shrink with the part, and each outward move raises the base node by the 3D gradient times the move,
which admits the next cell at the next part (offline, with no time passing, 0.94, 0.55 and 0.36 g in the first three
rounds; `ratchet.py`, `cells_held.py`). From 1.5625 s the part asked for 1.026 allowances at 62.5, 31, 16 and 7.8 ms;
8 of the last 15 parts were kept at the split limit (1.95 ms) over the allowance (1.015-1.19), each handing down
0.65-1.38 g whatever its length, against 0.8-1.7 g per 62.5 ms part before: the 12.38 g held came down to 2.14 g in
0.039 s, about 260 g/s, and the folds at the side edge rose from 117 to 159 degrees. *Not Task 9 alone* (the reduced
pushed plate to 2.0 s with each rule for the hand-down, every move, fold and refused move logged): Task 8's own code
folded the face's ridge past 100 degrees by 1.28 s, with moves refused and rescued by rebuilds from 1.13 s; Task 9
without the limit failed at 1.92 s on the deposit's rim beside the spilled edge; Task 9 as built folded to 131 degrees
by 1.81 s; with the release capped at the allowance (a diagnostic) it reached 2.0 s splitting no part after 1.3 s, its
face's fold past 90 degrees at about 1.85 s; on a 90 mm plate whose runoff did not reach the edge before about 1.7 s
the face's ridge alone passed 90 degrees at about 1.70 s. *Decided* (Asha, 2026-10-27, option 1 then option 2, over
option 2 alone, a hand-down limited by its neighbours' depth or local refinement (option 3, the fallback) and a test
plate without square edges (option 4); spec 6.6, 7.1, 12, 16): a part hands down at most what moves its base by the move
allowance (`hand_down(..., cap)`, `SkinBody._hand_cap`, HAND_DOWN_SHARE of the allowance; ruling 8), a cell going in
part, the rest waiting, so a hand-down never splits a part; and the pushed plate is judged at 1.5 s, the deposit's
folding a watched risk for sub-plan 22's flight. *The capped full plate*
(`item4/run_fold.py --resume`, from the planning copy's run again to 1.5 s, pickled there; `runs/cap9_w170`). The cap
would not have bitten before 1.5 s: the largest hand-down in that run's 17 parts was 0.532 mm, 0.957 of the cap, so the
capped code takes every hand-down to the last bit as it did (`item4/cap_check.py`), and the resumed body's first part
is the uninterrupted run's to the last digit logged. From 1.5 s the macro step ran in parts of 0.05 s, none split (no
move ratio above the cap's 0.999999); the cap held metal back from 1.60 s, in 1, 10, 11, 17 and 28 columns a part --
at the spilled edge, its side and the face's deposit, which would have handed down up to 0.64-1.14 mm -- the parts
handed down 1.6-2.1 g each, and the held metal fell from 12.0 g at 1.5 s to 10.6 g at 1.75 s, about 26 g/s against 260
g/s without the cap. The surface still folded: its sharpest fold 120 degrees at 1.55 s, 131 at 1.65 s on the spilled
front edge and 127.5 at 1.75 s on the side's back edge (triangles tilted 73 and 159 degrees from the face's direction);
the face's own 74, 78, 90, 98 and 99 degrees at a radius of about 40 mm, 36 then 110 face corners over 45 degrees; the
surface quality from 0.73 to 0.58. A move refused at the face's centre at 1.70 s was rescued by a rebuild; in the part
to 1.80 s the move was refused at the face's deposit 50 mm out and, after the rebuild, at corners of the spilled edge,
its side and its back edge: "a corner move of up to 0.000603 m degrades the mesh even after a rebuild" -- 0.2 s later
than without the cap, on the deposit's folding, as option 2 expects. Books at every kept part to 1.75 s (the replay run again with
them logged, `runs/cap9b_w170`, the same to the last digit logged): mass to 2.2e-16, energy to 7.1e-11, momentum to
4.9e-18; the counter counted 7 faults in the five parts kept (30 to 1.5 s; not examined -- the reduced plate's were the
side face's draws). The five parts cost 2 738 s of processor time, about 550 s a 0.05 s part (on battery power, alone).

*The regression* (`item4/cap_state.py`: the reduced plate to 2.0 s, seeded, with the cap and with it lifted, compared bit
for bit by `compare_states.py`). With the flow off (700 K, 6 MW/m2, the liquid removed as it forms) the cap never holds
metal back and the two runs are the same to the last bit, every array and history column. With the flow on and the
flight push they are the same in every history column at 0.5 and 1.0 s; the cap first holds metal back at 1.4375 s, in
one column at the spilled edge by 0.014 mm, then in each of the next ten parts (from 1.5 s parts of 0.05 s), in at most
8 columns and by at most 0.24 mm, mostly at the edge (38 of 62 holdings, the face 12, the side 12). At 2.0 s, with the
cap against without: the sharpest fold 103 against 132 degrees, the surface quality 0.57 against 0.48, 2 rebuilds
against 8, 77 interface repeats against 105, 6.16 g handed down against 7.92 g, 0.28 g held against 0.01 g, the
centre receded 4.600 against 4.587 mm; the counter counts 6 faults against none, all on the spilled side face at 1.4375-1.5 s, in columns that
drew while their slurry drained (one also handing metal down), their bottom cells 9-17 K above their nodes -- the kind the full plate counted at
1.5 s -- and none in a column whose hand-down the cap held back (`item4/fault_trace.py`). Mass, energy and momentum
close to 9e-16, 2e-11 and 2e-17 in both.

**Rulings** (2026-10-25/27, within spec 6.5, 6.6 and 7.1-7.2 as amended 2026-10-26/27 and the rulings of Tasks 1-8; each is
what, then why, then the cost if wrong):
1. *The reproduction: Task 8's pushed plate reduced to 60 mm, heated within 15 mm and to nothing at 21 mm, no fine band,
   to 1.0 s.* Why: item 1; it gives faults, bottom cells and rejections of the same kind and of comparable size at a
   twentieth of the cost, and its base operations can be recorded column by column. Cost if wrong: a mechanism of the
   full plate's alone missed -- the full plate re-run (item 5) shows the same fixes act there; its runoff reaches the
   reduced plate's edge after 1.3 s, so it is no measure of the plate's later steps.
2. *The limit on hand-downs is the base node's temperature after the part's 3D solve plus JUMP_TOLERANCE (5 K), the
   fault tolerance spec 6.6 names, in the part's move ratio and its move alike.* Why: that node is the surface the
   handed metal joins, and the move ratio must see the hand-down the move will make. Cost if wrong: hand-downs judged
   against a temperature the 3D model no longer has -- held a part longer or passed a part sooner.
3. *Where the limit ends a hand-down at a cell's top, the take's round-off remainder of that cell (below SLIVER of it)
   goes too.* Why: `take_bottom` leaves 1e-22 kg of the last cell when the take is the sum of whole cells, and the
   reconciliation's round-off then raised such a remainder to thousands of kelvin, which the counter, now reading
   after the reconciliation, counted. Applied only where the limit binds, so no hand-down it does not bind changes.
   Cost if wrong: none measurable (a round-off of mass).
4. *The counter reads the bottom cell after the reconciliation and before the regrid, and a hand-down's range runs
   from the coldest and warmest of the cells it took from, whole or in part, to the cell it uncovered, read straight
   after the hand-down.* Why: spec 6.6 as amended (K2, with the uncovered cell: Asha, 2026-10-26); before the regrid,
   so that a sliver the reconciliation heated is seen before a merge hides it; the uncovered cell as the hand-down left
   it, so that the skin's own interior is no fault and what the reconciliation then does to it is judged. Cost if
   wrong: a reconciliation that heats an uncovered cell by less than the gap between it and the handed cells goes
   uncounted -- with A that gap is the skin's own gradient, and the reconciliation returns only the sweep's mass
   mismatch at the handed metal's enthalpy.
5. *The landing conduction: in every flow sub-step after the transport, the columns that metal landed on top of in this
   heat sub-step, from the flow sub-step it landed in, have the cells within LANDING_DEPTH = 0.4 mm of their top (at least
   two) conducted through the depth for the flow sub-step, adiabatically -- no heat across the block's top, whose
   heating and radiation are the heat sub-step's, nor across its bottom -- the Newton iteration's error kept by the
   block's bottom cell; before the release.* Why: item 3; the block holds the landed metal and as much of the cold metal
   beneath as heat reaches in a flow sub-step (sqrt(alpha h), 0.33 mm at 2 ms), so the film freezes where it lands, the
   base law is untouched and the books stay exact; adiabatic, its diffusion counted twice only within the block (the
   heat sub-step's conduction already took its cells). Cost if wrong: the measured sensitivity (item 3): 0.2 and 0.8 mm
   put 3.3 and 2.7 % of the deposit 5-10 mm past the heat against 3.0 % (its centroid 1.11 and 0.97 mm past against
   1.05 mm); a block too shallow bottles the film's heat and lets it flow on, a deep one doubles the
   diffusion of more of the column.
6. *A trial pass's flow sub-step count is 0 in its plan where its heat sub-step took none (nothing flowed), and the real
   pass then takes its own memory's count; elsewhere the larger of the trial's and its own.* Why: spec 16, 2026-10-26
   (item 4). Cost if wrong: none for the base law where the two passes flow alike (the counts are equal there); where
   they do not, the base law is made by a pass that flowed differently, which the damped repetition corrects.
7. *C's held metal measured as what a hand-down with no limit would take less what went, after every kept part; in
   the skin body as `handable` with no limit, just before the hand-down, less what the hand-down took.* Why: it is the
   metal C keeps in the skins and nothing else: a skin is above its minimum also by its flowing layer, by the sweep's
   mismatch and by its area's change, none of them C's; `handable` computes the take a hand-down with no limit would
   make from the same cells, so it is that measure exactly, without the copy of the skins or the take that
   `run_logged9.py` made. Cost if wrong: none for the model; a reading.
8. *A part's hand-down is capped at HAND_DOWN_SHARE = 1 - 1e-6 of the move allowance at every column, whatever the
   column drew, and the held metal is what a hand-down with neither the limit nor the cap would take less what went.*
   Why: item 9 (Asha, 2026-10-27); a cap at the whole allowance split parts by round-off (a move ratio of 1.000 plus a
   few 1e-16 in the diagnostic run, `item4/run_variant.py`), since the ratio divides the capped mass by the area again;
   a draw and a hand-down in one part move a base opposite ways, so each within the allowance keeps the net move within
   it; and metal waiting for the cap is held in the skins as surely as the limit's. Cost if wrong: a column that drew a
   whole allowance could hand down only one more in the part, a part's delay.

**What this task found in earlier tasks.**
- *Sub-plan 20's moving mesh* booked an outward sweep at the corner's own enthalpy, where spec 6.6 said the 3D model
  gains the handed cells' (item 1); harmless while hand-downs were thin and at the node's temperature, it made the
  superheated bottom cells where runoff deposits. Changed (A).
- *Sub-plan 20's fault counter* judged a hand-down against the node alone, before the reconciliation, so it counted the
  warmer cell a hand-down uncovers and never saw what the reconciliation did (item 1). Changed (K2, the uncovered cell
  in a hand-down's range; item 7).
- *`Skins.take_bottom`* leaves a round-off remainder of about 1e-22 kg of the last cell when the take is the sum of
  whole cells; the regrid merged it before, the counter now reads it (ruling 3). Folded where the limit ends a
  hand-down; left elsewhere, where the regrid still merges it.
- *Task 8's reading of the faults* (Review Focus 25): the reconciliation, yes, but its enthalpy, not the sweep's mass,
  and no neighbour's draw (item 1).
- *6d's counts*: the real pass at its trial pass's counts above the Courant target at an onset, and replaying counts
  chosen where nothing flowed (item 4). Changed (spec 16, 2026-10-26).
- *Task 8's heat sub-step study*: the 11 % of the deposit landing 5-10 mm further at 10 ms is 9.9 % of what lands beyond
  the heat within 20 mm of the centre line (Task 8 counted it against 345 mg deposited, this task against the 395 mg
  beyond the heat) (item 3).
- *The damped repetition where deposits form*: each 0.0625 s part of the pushed plate takes four or five passes,
  against fact 126's 1.8 on the sphere, which is now most of a melting macro step's cost (item 5; not investigated).
- *The fold and the warm bottom cells at 1.5 s*: in the pushed plate's last part the base surface's sharpest fold rose
  to 95.3 degrees and the bottom cells' largest excess over their node to 68.6 K, with 33 faults that were not uncovered
  cells (item 5): the runoff spilled over the plate's square edge, the warm bottom cells and the faults are the
  slurry's on its side face, and the fold is the deposit's (item 9).
- *Option C's limit* releases a held cell whole once it is near its node's temperature, a release that does not
  shrink with the part, so the part cap's splitting failed for it and the held metal came down at once (item 9).
  Changed (the cap, spec 6.6 and 7.1, 2026-10-27).

- [ ] **Step 1: The moving mesh's test** -- in `tests/test_reentry_model_movemesh.py`:

{{DIFF_TEST_MOVEMESH_T9}}

- [ ] **Step 2: The base's tests and the skin body's counts** -- in `tests/test_reentry_model_skinbody.py` (the hand-down
  limited by the node's temperature and what it would take, and capped, a cell going in part; the warm deposit that was
  the cause and the metal held; a deposit thicker than the move allowance coming down over several parts, splitting
  none; the counter judging a hand-down by the handed cells and the cell it uncovers, in its arrays and on a skin body;
  the real pass's counts in the skin body, which replaces 6d's test of the trial's counts; and the flow-off plate's
  hand-downs, which neither the limit nor the cap holds back; the parts test's recorded move passes the new valuation
  on, and the pushed plate's threshold is Task 8's ruling 6's ten, its docstring saying why):

{{DIFF_TEST_SKINBODY_T9}}

- [ ] **Step 3: The landing conduction's and the counts' tests** -- in `tests/test_reentry_model_skin.py` (the two flow
  sub-step tests built by hand take the landing conduction after the transport):

{{DIFF_TEST_SKIN_T9}}

- [ ] **Step 4: The skin run's results** -- in `tests/test_reentry_model_cli.py`:

{{DIFF_TEST_CLI_T9}}

- [ ] **Step 5: Run them to see them fail.** `"$PY" -m pytest tests/test_reentry_model_movemesh.py
  tests/test_reentry_model_skinbody.py tests/test_reentry_model_skin.py tests/test_reentry_model_cli.py -q -m "not
  reference" -k "outward_sweep_is_booked or only_bottom_metal_within or warm_deposit_goes_down or judged_against_the_handed
  or uncovers_is_no_fault or real_pass_takes_at_least or holds_no_hand_down_back or landing_on_cold_rigid or
  runoff_onto_cold_metal or larger_of_its_trial or in_that_order or releases_after_the_transport or too_far or
  writes_its_columns_frames or cap_a_cell or thicker_than_the_move_allowance"`. Expected: 16 failed --
  - `TypeError: MovingMesh.move() got an unexpected keyword argument 'h_outward'`;
  - `AttributeError: module 'reentry_model.skin' has no attribute 'handable'`, `TypeError: hand_down() got an unexpected
    keyword argument 'cap'` (the flow-off plate's) and `TypeError: hand_down() takes 3 positional arguments but 4 were
    given` (the cap's);
  - `TypeError: MovingMesh.move_with_policy() takes from 2 to 3 positional arguments but 4 were given` (the deposit
    thicker than the allowance; with Task 9 built but the cap lifted, `assert (2.513518154313309 < 1.0)`, the first
    part asking for 2.51 allowances at every length from 0.5 s to 0.69 ms, measured 2026-10-27);
  - `AssertionError: 1614.0 K in a column whose hottest metal was 890.0 K` (the cause);
  - `TypeError: jump_faults() got an unexpected keyword argument 'handed'`;
  - `AssertionError: assert (np.float64(99.8087120132434) < 0.5)` (the uncovered cell on a skin body: with no limit the
    warm metal goes down too, and the bottom cell the reconciliation leaves is 97-100 K above the 760 K expected);
  - `KeyError: 'held_mass_max_kg'` (the skin run's results);
  - `assert [10, 10, 10, 10, 10, 5, ...] == []` (the skin body's real pass chose no count of its own);
  - `AttributeError: module 'reentry_model.skin' has no attribute 'landing_conduct'`;
  - `assert [np.int64(10), np.int64(11), np.int64(12), np.int64(13), np.int64(14)] in ([], [10])` (landed metal still
    above T_rigid five columns past the edge);
  - `assert [1, 1, 1, 1, 1, 1, ...] == [10, 7, 9, 8, 7, 6, ...]` (a real pass at its trial pass's one flow sub-step);
  - `AttributeError: 'TransportBook' object has no attribute 'landed'` (2: the flow sub-steps built by hand);
  - `TypeError: MovingMesh.move_with_policy() takes from 2 to 3 positional arguments but 4 were given` (the parts test's
    recorded move).

- [ ] **Step 6: Handed-down metal at its own enthalpy** in `reentry_model/movemesh.py`:

{{DIFF_MOVEMESH_T9}}

- [ ] **Step 7: The limited hand-down, the landing conduction and the real pass's counts** in `reentry_model/skin.py`:

{{DIFF_SKIN_T9}}

- [ ] **Step 8: The base's operations, the held metal and the counter** in `reentry_model/skinbody.py`:

{{DIFF_SKINBODY_T9}}

- [ ] **Step 9: The held metal's history columns and results** in `reentry_model/coupled.py`:

{{DIFF_COUPLED_T9}}

- [ ] **Step 10: Run them to see them pass.** Expected: `tests/test_reentry_model_movemesh.py` 39 passed,
  `tests/test_reentry_model_skinbody.py` 28 passed (2 deselected, the reference tier's; 21 before),
  `tests/test_reentry_model_skin.py` 90 passed (87 before), `tests/test_reentry_model_cli.py` 66 passed (one
  assertion added); the unit tier 517 passed, 1 skipped (1 276 s alone before the cap's two tests, 2 025 s under a load
  of about ten processes after; 506 before). Every skin flight with the flow off moves only by
  the booking of handed-down metal (item 8); with the flow on, by all four changes.

- [ ] **Step 11: The analysis runs** (each on mains power under `caffeinate -i`, the scripts in `measure/task9/`, the
  long ones detached by `launch9.sh`, which logs the power source at their start and end):
  - the cause, on the planning copy before this task and on the prototypes (`TASK9_CODE=.../proto_A`, `proto_C`,
    `proto_AC`): `"$PY" trace_move.py 60 1.0` (262 s alone) and `"$PY" analyse_trace.py runs/trace_w60`; the failing
    test first written as `test_task9_cause.py`; the fixed code on the same plate, `"$PY" count_faults.py 60 1.0
    runs/fix_w60.json` (169 s alone; item 2);
  - the landing conduction: `run_substeps9.sh` (`plate_skins9.py substep`, 50 ms from Task 8's pushed skins at 0.2 s,
    heat sub-steps of 10, 5 and 2.5 ms at 0.5 ms flow sub-steps, depths of 0.2 and 0.8 mm, the chosen count, with the
    landing conduction and without) and `"$PY" compare_substeps9.py` (item 3; 45-151 s a run under load);
  - the counts and the parts on the reduced plate before and after: `run_logged9.py --case cold --dt 0.5 --t-end 1.0
    --flow on --push flight --swell on --width-mm 60 --r-flat-mm 15 --r-zero-mm 21 --fine-mm 0` (item 4), and again
    under the adopted fault wording, `--out red60_k2.json` (item 7; 202 s alone);
  - the pushed plate: `launch9.sh push9_w170 run_logged9.py --case cold --dt 0.5 --flow on --push flight --swell on
    --width-mm 170 --t-end 2.0 --save-state 1.0 --columns push9_w170.npz` and `"$PY" analyse_push9.py` (items 5 and 6;
    to 1.5 s in 3 h 12 min of the run's own time, stopped there when the laptop went onto battery power; its outputs
    kept as `runs/push9_w170_to1.5.*`, its logging script as `run_logged9_to1.5.py`); then, with the logging extended
    for the examination of its end (the faults' split, the warmest bottom cell over its node, the fold's place, the
    body's held metal), the same to 2.0 s with `--save-state 1.5`, and `"$PY" analyse_rerun9.py runs/push9_w170`;
  - the failure at 1.60 s and the cap (item 9): `launch9.sh fold9_w170 item4/run_fold.py --pickle-at 1.5 ...` (the
    pushed plate again, every part tried, kept and refused logged, the whole body pickled at 1.5 s) and the same on Task
    8's code (`TASK9_CODE=.../item4/task8tree/...`, `fold8_w170`); the reduced plate to 2.0 s with each rule for the
    hand-down (`run_fold.py`, and `item4/run_variant.py noC|capped`; `red8_*`, `red9_*`, `rednoC_nof`, `redcap_nof`,
    `red9w90_nof`), read with `"$PY" item4/analyse_fold.py runs/NAME`; the held cells and the per-move release on the
    1.5 s state (`item4/cells_held.py`, `item4/ratchet.py`); the cap checked against the run as built to 1.5 s
    (`"$PY" item4/cap_check.py runs/fold9_w170 1.5`) and the replay from the pickle with the cap, `launch9.sh cap9_w170
    item4/run_fold.py --resume fold9_w170.body.pkl --case cold --dt 0.5 --push flight --width-mm 170 --t-end 2.0`; the
    regression, the reduced plate to 2.0 s with the cap and with it lifted (`CAP_OFF=1`), flow on and off, compared bit
    for bit (`item4/cap_state.py`, `compare_states.py`);
  - flow off: `run_flowoff9.sh` (the level-2 plates to 1.0 s, seeded, with the fix and with the limit alone,
    `TASK9_NO_A=1`, against the same before this task) and `"$PY" compare_states.py`; `launch9.sh off_cold9
    run_logged9.py --case cold --dt 0.5` and `off_hot9 ... --case hot` to 2.0 s (item 8; 39 and 13 minutes under load);
  - the regression: `launch9.sh reg_cold9 run_logged9.py --case cold --dt 0.5 --flow on --release all --swell off` and
    `reg_hot9 ... --case hot ...` (51 and 17 minutes under load), and on the hot plate each change undone in turn,
    `"$PY" plate_toggle.py noA|noC|noland|oldcounts --case hot --dt 0.5 --t-end 1.5 --flow on --release all --swell
    off`, and on the reduced pushed plate (`run_toggles60.sh`) (item 8).

- [ ] **Step 12: Commit.**

```bash
git add reentry_model/movemesh.py reentry_model/skin.py reentry_model/skinbody.py reentry_model/coupled.py tests/test_reentry_model_movemesh.py tests/test_reentry_model_skin.py tests/test_reentry_model_skinbody.py tests/test_reentry_model_cli.py
git commit -m "Skins: the base where runoff lands -- handed-down metal booked at its own enthalpy and only metal within the fault tolerance of the base node handed down, so that no bottom cell is superheated and no part is rejected for a deposit's moves (spec 6.6 as amended 2026-10-26, options A and C); a fault judged against the metal that crossed, for a hand-down the handed cells and the cell it uncovers; the metal the fault tolerance holds recorded in the history; metal landing on a rigid top losing its heat to the metal beneath within the flow sub-steps; and the real pass taking the larger of its trial pass's flow count and its own; and a part's hand-down capped at its move allowance, a cell going in part, so that a hand-down never splits a part (spec 6.6 and 7.1 as amended 2026-10-27) (sub-plan 21 Task 9)"
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
9. **The base where runoff lands** (Asha, 2026-10-25, option a, over leaving it to sub-plan 22; its fix decided
   2026-10-26, spec 6.6 and 16) -- the cause of the base faults and rejected parts where runoff lands on cold metal,
   handed-down metal booked at its own enthalpy and only metal near the node's temperature handed down, the fault judged
   against the metal that crossed, the runoff's heat lost within the flow sub-steps and the real pass's flow counts --
   is written out above as Task 9.

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

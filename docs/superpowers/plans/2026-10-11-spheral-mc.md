# Spheral MC: Coupling the re-entry model and Spheral at the surface — Implementation Plan

> **For agentic workers:** execute task by task with one sub-agent per task; the controller reviews each task's
> diff and numbers before the next starts. Steps use checkbox (`- [ ]`) syntax for tracking.

> **Rescoped at review (Asha, 2026-10-11; "Decisions at review" below).** This version of Spheral captures large
> fragmentation and necking; it does not lose mass through liquid spraying. The continuous spray sink and everything
> built only for it — release tables by inclination, the sink package, the physical-mass field, continuum-surface-force
> shares, slabs and the synthetic release, check 15, the sink on the flight — moved to **"Deferred: sub-grid Girin
> spraying (later milestone)"** at the end of this plan, with their facts and decisions kept. MC now has four tasks.

**Goal:** the permanent particle IDs survive deletion, restarts, the switch to a window and MPI redistribution (spec
§6.2, §14.1, the ID part of check 14), so that a fragment's particles can be traced to where they started in the sphere
and followed from one check to the next; the window driver's core, which M3's separation windows and M4's slurry and
detachment windows extend, runs from a replay checkpoint through M2's switch with its accounts closed; and the release
the frames carry, together with the release on faces a step's deaths removed, closes against the finite-element
history's sprayed mass (the re-entry model's film and mass account; M1 decision 16, Step 3 fact 107 (a)). Mass leaves
Spheral in this version only by the replay's whole-particle centre-crossing removal (M2, D11).

**Done when** (spec §14, amended below): the ID part of check 14 (spec §12) passes with SolidSPH in RZ at 2.2 and
1.1 mm and in 3D at 2.2 mm, through centre-crossing removals, restarts, the switch and a window, on 1, 4 and 18
processes (3D on 4 and 18); Spheral's accounts close to round-off through `analyse` in the replay and in a window; the
frames' release plus the removed faces' release closes against the history's sprayed mass within the residual Task 1
measures, with the M1 release thresholds on the closed sum; every acceptance threshold is measured and recorded in the
README, `docs/spheral_frag_assumptions.md` and the tests; MC's measured facts and the rescoping are written into the
spec (Task 4).

**Answers:** whether the permanent IDs survive deletion, restarts, the switch and redistribution on any process count;
whether the frames' release closes the history's sprayed mass once the removed faces' release is added, and what the
remaining positive excess is (fact 11); what a window's restart costs, which sets the window's sub-interval. The spec's
original MC answers — whether mass can leave Spheral continuously, the deletion floor, whether the smoothing length
shrinks with the mass — go with the deferred sink.

**If it fails:** IDs lost or duplicated → a defect in the runner, fixed in Task 3 (fact 9 shows Spheral itself carries
registered fields through deletion and redistribution); the flight rerun with export addition 07 not byte-identical →
stop and report, D7's fallback (`--removed-release history`) closes the total and the frames stay as they are; a body
out of balance after the switch → M2 Task 11's finding, reported before Task 2 builds on it. Whole-particle deletion
(spec §19) is no longer MC's fallback: it belongs to the deferred sink.

**Spec:** `docs/superpowers/specs/2026-10-02-spheral-large-fragments-design.md` — §1 (purposes and routes), §2 (facts:
particle identity and mass, traps), §6.2 (permanent IDs), §9.1 (the replay; centre crossing), §9.3 (windows, the
switch), §10 (the three zones), §11.3 (accounts), §12 check 14, **§14.1 (milestone MC)**, §15 (risk "Particle identity
and mass"), §17 (the release gap); for the deferred work also §8.1–8.2, §12 check 15 and §19. Builds on the M2 plan
(`2026-10-09-spheral-m2.md`: its facts, D1–D6 and "Decisions at review" 1–9, cited as *M2 decision n*), M1 decision 16
(the release on removed faces deferred to MC) and Step 3 shared context fact 107 (a)
(`docs/superpowers/plans/melt-spraying-subplans/00-shared-context.md`). Ordinary task-by-task plan with code shapes
and tests, not generated from a prototype (Asha, 2026-10-01).

**Amendment to spec §14.1 and §12 (M2 decision 8, Asha 2026-10-11).** The scheme is SolidSPH: SolidCRKSPHRZ fails the
patch test at `116c71f` and SolidCRKSPH in 3D is no better overall (check 7). §14.1's condition 2 ("mass-over-density
volumes under CRKSPH") therefore does not apply, and condition 1 (integrated density) already holds in the replay
(`run.make_hydro`: `densityUpdate=S.IntegrateDensity` for both schemes). Where spec §12 and §14 say "for each scheme
still in contention", MC runs **SolidSPH** only. `--scheme crksph` stays a bracket for a later Spheral version, and MC
does not test it.

**Not in MC:** the continuous spray sink and everything built only for it (deferred, end of this plan); the windows
proper — triggers, end conditions (clear, settled), separation, slurry and detachment windows, chaining (M3, M4); the
softened configuration, including the Murnaghan equation of state (M4; deferred task S4 Step 3 records its shape);
the heating package and free energy in windows (M4); the fragment finder in a run (M3); the hand-back (M6); the
cluster. MC builds the window driver's core (Task 2) on M2 Task 11's switch, which M3 and M4 extend.

## Decisions at review

Asha answered the plan's seven questions on 2026-10-11 (the questions as asked are kept at the end of this plan):

1. **The release on removed faces (Q1, D7): option (4) accepted.** Export addition 07 writes a side file
   (`release_removed.npz`) from a byte-identical rerun of the reconstruction's flight; option (2), from the history,
   is the fallback and a bracket. It stays in MC (Task 1) because it closes the re-entry model's own film and mass
   account — the frames' release against the history's sprayed mass, the 12 % gap of spec §17 — whether or not
   Spheral ever sprays; the per-face detail is what the deferred release tables (S1) will read.
2. **Which particles sink (Q2, D5) — answered with the scope.** Asha's words: *"this early version of spheral is not
   being used to lose mass through liquid spraying. We want to capture large fragmentation and necking processes. We
   will later implement Girin spraying as a sub-grid parameterization."* Consequences:
   - **No continuous spray sink in this version of Spheral.** Spec §14.1's decision of 2026-10-02 ("spray mass leaves
     continuously") is deferred with it, not reversed: when spraying is taken up, D5's question (the frame's release
     or the particle's temperature decides which particles sink) is open again, neither accepted nor rejected.
   - **Spheral's purpose in this version is large fragmentation and necking** (spec §1's routes 1–3: tearing, loss of
     strength, detachment of pieces that melting has nearly cut off).
   - **Recession of the body stays M2's mechanism:** in the replay a particle is removed whole when the finite-element
     surface passes its reference position's centre (M2 decision 2, D11, Q5), booked as `centre_crossing` with its
     body mass and held temperature's enthalpy. That is how necks thin in Spheral. Windows remove no mass in MC (new
     question 1).
   - **Girin spraying comes later as a sub-grid parameterization.** The sink package, the physical-mass field, the
     csf shares, the release tables by inclination, slabs and the synthetic release, check 15 and the sink on the
     flight (old Tasks 2, 3, 4, 6, 7, 9) and decisions D2, D3, D4, D5, D6, D8, D9, D10, D12 move to "Deferred: sub-grid
     Girin spraying", content and probe facts kept. Where that work sits is new question 2.
   - **Kept in MC:** Task 1 (above), the window driver's core without sink, floor, physical-mass field, free energy and
     softened equation of state (Task 2, was Task 5), the ID part of check 14 with centre crossings and restarts in
     place of sink steps and floor deletions (Task 3, was Task 8), documentation (Task 4, was Task 10).
3. **The mapping by inclination, normalised per step (Q3, D6): accepted.** Applies when the deferred work is taken up.
4. **The area shares under the sink (Q4, D4): accepted** — build the `csf` shares and choose by check 6 and check 15.
   Deferred with the sink (S3); the replay and MC's windows keep the patch shares (M2 decision 9). S3 records why the
   csf rule's surface re-detection is also a candidate for M3's and M4's moving surfaces.
5. **The replay keeps centre-crossing removal (Q5, D11): accepted.** Spec §14.1's open item "The replay" is closed.
   With the sink deferred, centre crossing is the only way mass leaves Spheral in this version, apart from the
   confirmed fragments M3 removes.
6. **Check 15's extensions (Q6): accepted.** Apply to check 15 when it is run (S5).
7. **Floor removals at sub-boundaries (Q7, D8): accepted.** Apply when the sink is built (S2); Task 2's sub-boundary
   hook is where they will act.

**Spec amendments these imply** (written by Task 4, not now):

- **§3 (scope):** the first version (M0–M4, MC) no longer includes "the spray mass sink"; permanent particle IDs stay.
  The sink joins the "built only when triggered" list with the sub-grid Girin spraying (new question 2).
- **§4.3 (package layout):** `packages/`'s "spray sink" marked as the later milestone's.
- **§10:** "In windows, through a mass sink …" replaced by this version's rule (windows remove no sprayed mass; per
  new question 1, the omitted recession recorded per window); the sink paragraph moves to the later milestone's
  section with D5's open question; the M1 note's last sentence (the sink's `release_rate` inherits the 12 % deficit)
  replaced by Task 1's closure.
- **§12:** check 14's ID part as MC runs it (centre crossings, restarts, the switch, a window, redistribution; no sink
  or floor); check 15 marked for the later milestone, with Q6's extensions.
- **§14, the milestone table:** MC's row — Builds: the permanent IDs verified through the replay, the switch and a
  window; the window driver's core; the release on removed faces (export addition 07) and the film account's closure.
  Answers: as this plan's "Answers". Done when: as this plan's "Done when". M4's row loses "MC's spray sink in the
  windows". A row or note for the later sub-grid spraying milestone (new question 2). "If an early milestone fails":
  MC's clause ("the continuous sink misbehaves — whole-particle deletion") moves to the later milestone. The
  flowchart's MC node ("mass leaves continuously") becomes "IDs, window core, release closure".
- **§14.1:** the decision "spray mass leaves continuously" (2026-10-02) recorded as deferred (Asha, 2026-10-11, her
  words); "The sink", the four conditions and the smoothing-length open item move to the later milestone's section
  with MC's facts 1–8 and 14–16 (condition 1 holds as built, condition 2 void by M2 decision 8); "Open within MC: the
  replay" closed by Q5; "Checks" becomes check 14's ID part; "If MC fails" rewritten as this plan's.
- **§2, fact "Particle identity and mass":** amended with facts 1–10 (measured, no longer "untested").
- **§15, risk "Particle identity and mass":** the mass part moves to the later milestone; the identity part is
  answered by Task 3.
- **§17, the release gap:** closed by Task 1 (or bounded by its fallback).
- **§18:** "the spray mass sink from the frame's rate … the sink's deletion floor" moves to the later milestone;
  centre-crossing removal stays.
- **§19:** the whole-particle-deletion row becomes the deferred sink's fallback, not MC's; the row "Droplet spraying
  as a sub-grid model in Spheral now — deferred to M9" gains Asha's 2026-10-11 statement.
- **§20:** a 2026-10-11 entry — the rescoping in Asha's words, Q1 and Q3–Q7 accepted, SolidSPH (M2 decision 8).

## Facts established while writing this plan (2026-10-11)

Facts 9–13, 17 and 18 bear on MC as rescoped. Facts 1–8 and 14–16 were measured for the sink; they are kept for the
deferred work (S1–S6), and fact 16 also for M3's moving-surface loads (S3's note).

Measured in the pinned image (`particle-formation/spheral:116c71f-arm64`, through `spheral_frag/container/spheral`)
with throwaway probes in `spheral_output/mc_probe/` (git-ignored: `sink_probe.py`, `release_gap.py`,
`skin_temperature.py`, `csf_area.py`), read from Spheral's source in the image (`/home/spheral/workspace/src`) or from
the code on `main`. The Spheral probes use a lattice slab at dx = 2 mm (3D 8 × 10 × 10 = 800 particles; RZ 8 axial ×
12 radial rows = 96 rings), AA7075 at 300 K through the runner's own `PyThermalEOS`, `PyThermalStrength` (`elastic`)
and `nodes.make_body`, SolidSPH with compatible energy as `run.make_hydro`, CheapSynchronousRK2; the top layer
(largest x) loses about half its mass over 200 steps.

1. **Spheral's mass has no update policy.** `SPHBase::registerState` enrols the mass field bare
   (`state.enroll(mass)`, `SPH/SPHBase.cc` line 193); SolidSPH adds none. A package that changes the mass field
   between steps is therefore not overwritten.
2. **Lowering the mass in a package's `finalize` works and closes to round-off.** A Python `Physics` package that
   subtracts q_i dt from the top layer's masses in `finalize` after every step (as `ImposedTemperature` sets u):
   3D, mass lost 6.27541903623708e-2 of the total against q t = 6.27541903623705e-2, closure (remaining + booked −
   initial) −1.9e-15 of the total; RZ 1.7e-16. With integrated density the sinking particles' density is unchanged
   bitwise (3D) or to 1.8e-15 (RZ), the layer below unchanged, and the stress-free slab stays at rest (|v| 0.0 and
   |S| 0.0 in 3D; in RZ 3.9e-12 m/s and 1e-4 Pa, against 2.6e-12 m/s and 1.8e-5 Pa without the sink).
3. **Spheral's increment policy cannot carry the mass.** Enrolling the mass with `ScalarIncrementState` and filling a
   derivative field "delta mass" fails on the first step: "IncrementState ERROR: unable to find unique match for
   derivative field key delta mass" (`DataBase/IncrementStateInline.hh` line 79), because the policy matches every
   derivative key that *starts with* its prefix + field name, and "delta mass density" also does. A wildcard policy
   would add the density's rate to the mass. The sink therefore changes the mass in `finalize` (D2).
4. **In RZ the density is derived from the ring's two masses.** `AxisymmetricMassDensityPolicy` sets ρ = m / (2π r A)
   with A = massRZ / ρ_RZ (`SPH/SPHRZUtilities.hh`, `computeMassDensityFromArealMassDensity`). Scaling m and massRZ by
   the same factor keeps ρ (fact 2); lowering m alone does not: in the same 200 steps ρ fell by 8–27 % at the top,
   the slab reached 452 m/s and |S| 8.6 GPa.
5. **The source of the RZ masses' last-bit moves** (`runner/nodes.py`: "the source was not located") is
   `SolidSPHRZ::applyGhostBoundaries` and `enforceBoundaries`, which divide every mass by r before the boundaries and
   multiply it back after (`SPH/SolidSPHRZ.cc`, read from the code). Accounts that must close exactly therefore keep
   their own mass, not Spheral's: the replay already books body.npz's masses (fact 17), and so does MC's window
   (Task 2); the deferred sink would need D12.
6. **Summed density fails at the free surface with or without a sink.** With `RigorousSumDensity` the same slab,
   *without* any sink, reached ρ −1.3 to +24 % of ρ_free at the top, 54 m/s and |S| 10 GPa in 200 steps (with the sink
   309 m/s and 10.7 GPa). Integrated density (spec §14.1 condition 1) is required anyway; M2 already uses it.
7. **The smoothing length does not respond to the mass.** The SPH factory's default `HUpdate = IdealH` enrols H with a
   `ReplaceBoundedState` policy, the ideal H computed from the neighbours' positions (`SmoothingScale/
   SmoothingScaleBase.cc` lines 69–77). With and without the sink the final H fields are identical bitwise (max
   relative difference 0.0). On the slab H itself moved from 2.01 dx to 2.42–3.11 dx near the free surfaces in both
   runs (IdealH answering the missing neighbours, not the mass).
8. **H can be shrunk with the mass from Python, but not naively.** The ideal H is a derivative field "new H" that a
   package appended after the hydro may rewrite in `evaluateDerivatives`; only the rewritten particles' H changes
   (all others bitwise equal). Multiplying it by (m₀/m)^(1/dim) at every step compounds through IdealH's own
   iteration: at m/m₀ = 0.545 the top layer's h ended at 0.56–0.66 of the unshrunk run's, against the intended
   0.817, and the shorter h shortened Spheral's step (200 steps covered 1.64e-5 s against 1.81e-5 s). A
   non-compounding rule is needed (D3).
9. **MPI redistribution carries the IDs.** On 4 processes, after deleting half of rank 0's particles,
   `PeanoHilbertOrderRedistributeNodes(kernelExtent).redistributeNodes(db, [])` moved particles between ranks (3D
   counts 100/200/200/200 → 175 each; RZ 12/24/24/24 → 21 each), and the registered `pid` IntField, the `x_ref`
   VectorField, the mass and a user ScalarField followed every particle bitwise. The replay itself redistributes at
   every frame boundary through `distributeNodes` (M2 Task 1, 4 → 1 processes bitwise); in-run redistribution
   (`redistributeStep`) is off.
10. **`deleteNodes(vector_of_ULL)` mid-run** with SolidSPH alone: stepping continues on the remaining particles (795 of
    800, no motion). The runner's own packages cache their local rows and refuse a changed count
    (`packages.LocalCount`), so deletion and redistribution happen only at restarts unless the packages are rebuilt.
11. **The release gap on the frames Spheral uses** (the reconstruction's 100 mm Scheil flight, M2 fact 8;
    `data/spheral/m1_frames.csv`): Σ `release_rate` · A over the flight is 1.025 kg against the history's sprayed
    1.169 kg (−12.3 %); per spraying frame (207, k 49–426) −29.5 % to +3.8 %. The positive excess is small and rare:
    13 frames, at most +3.8 % (k 396) and +1.4 % (k 426), the other eleven below 1.5e-3. Shifting the frames by one
    step against the history makes it far worse (−75 % to +360 %), so on these frames the excess is not a one-step
    offset (fact 107 (a)'s guess); fact 107's −40 % to +13 % were on the newer frame-export flight.
12. **The removed faces' release cannot be recovered from the frames.** Estimating it as the previous frame's
    `release_rate` · A on the faces that disappear between frames k − 1 and k gives 0.224 kg over the flight against
    a gap of 0.144 kg (correlation 0.90 per step; after adding it, −11 % to +36 % per step and +6.9 % in total). A
    step's dying faces do not spray at their previous rate; only the export knows what they released (D7).
13. **Spray is the recession on this flight.** From the history, cumulative sprayed mass over cumulative solid-mass
    loss (mass − film − deep) is 0.965 at 50 s and 0.996 at 100 s, and 0.96–1.03 over the ten steps before 50 and
    100 s. In the replay, centre crossing therefore removes material at the spray rate within a few percent (D11);
    in a window that removes nothing, the omitted recession is the frame's release over ρ_l (new question 1).
14. **Spray rates per area** (every third frame, k 49–426): the area-mean over a frame's spraying patches has median
    3.3 and maximum 5.7 kg/m²/s (k 196: 2.4 mm/s of liquid at 2,400 kg/m³); a frame's largest single patch rate has
    median 57 and maximum 355 kg/m²/s (k 151). Per patch the release is spiky; by inclination it is smooth. At the
    area-mean maximum a 50 ms stiff window omits 0.12 mm of recession (0.05 of a 2.2 mm spacing), a 1 s soft window
    2.4 mm (about one spacing).
15. **Most of the spray comes from patches whose material at a particle's centre is below the liquidus.** At frames
    60, 100, 150 and 182, the share of the frame's released mass whose patch has T ≥ 908 K at a depth d along
    `n_derived`: d = 0.25 mm 1.3/80/83/94 %; 0.55 mm (a 1.1 mm particle's centre) 1.1/47/71/81 %; 1.1 mm (a 2.2 mm
    particle's centre) 0.3/11/41/32 %; 2.2 mm 0/0.8/23/8 %. Median T at 1.1 mm: 890–907 K (mush). The skin that sprays
    is thinner than a particle (spec §10), so §10's "exposed surface particles above the liquidus" would drop most of
    the spray (D5).
16. **A continuum-surface-force share follows a shrinking surface.** On lattice spheres of radius 50 mm, Σᵢ |∇c|ᵢ Vᵢ
    (|∇c| from `packages/loads.colour_gradient` at h = 2.01 dx) is 0.459 (5 mm) and 0.485 (2.2 mm) of 4πR² — the free
    surface's half-kernel, a constant to calibrate — and its share in 30° inclination bands lies within 0.465–0.509 of
    the exact bands at 2.2 mm (−4.0 to +5.0 % about the mean; M2 Task 4's volume shares misplace 10–13 %, M2
    decision 9's patch shares are exact). When the outer layer's volumes drop to a quarter, the total moves by −3.3 %
    and the outer layer's share falls from 0.378 to 0.096: the load passes to the layer below by itself (D4).
17. **What exists.** `record.py` already has `REASON_SINK` and `REASON_FLOOR` and splits the film account by reason;
    `analyse` closes the accounts to 1e-12 m₀. The replay's record checks book body.npz's masses, not Spheral's
    (`Replay.write_check`), and its removals book body.npz's mass and the held temperature's enthalpy. The permanent
    IDs (`pid`, `x_ref`) are built (M2 Task 1, M2 decision 6) and carried through every checkpoint. M2 Tasks 1–9 are
    built (decision 8 SolidSPH confirmed, decision 9 patch shares the default); Task 10 (check 13) is being run by
    another session on 2026-10-11; Tasks 11 (the switch) and 12 (documentation) are not yet. MC leaves
    `REASON_SINK` and `REASON_FLOOR` in place, unused, for the deferred work.
18. **The finite-element run is seeded** (seed 12345, spec §2 "Frames as built"); Step 3 fact 105 found coupled runs
    with and without frames identical in every history column once seeded, so a rerun that only adds a write-only
    output can be checked byte for byte (Task 1).

## Global constraints

- **Interpreters.** `PY=/Users/ashajain/miniforge3/envs/drama_env/bin/python` runs the core, the finite-element
  prototype, the body stage and the analysis. Everything that imports Spheral runs in the container through
  `spheral_frag/container/spheral [-n N]`. Spheral's Python never imports `reentry_model`.
- **Three stages, files only** (spec §4.2). The prepared frames and `prepare.json` are not rewritten; Task 1 adds a
  side file beside the finite-element run, so every M2 body and run keeps its provenance.
- **The frames** are the reconstruction's MVP frames (M2 fact 8).
- **Package layout** as M2: `spheral_frag/runner/run.py`'s `window` mode; tests
  `tests/test_spheral_frag_runner_{window,check14}.py` (marker `spheral`). No `sink.py` or `release.py` in MC.
- **Units** SI; release rates kg/m²/s at the particles, kg/m² per macro step in the frames (M1's contract).
- **Outputs** git-ignored: runs `naming.run_dir(spheral_output/runs, name)` (the nested layout of M2 Task 7b), check
  studies `spheral_output/mc/<check>/`, probes `spheral_output/mc_probe/`. Committed: code, tests, fixtures under
  `tests/fixtures/spheral_frag/`, measured tables `data/spheral/mc_*.csv`, documentation.
- **Run names encode the whole configuration** (CLAUDE.md): MC's new brackets (D13) are suffixes only when they
  differ from their defaults, so every existing replay name stays unchanged; a run whose `summary.json` exists is
  skipped.
- **Exit codes** as in `reentry_model`: 0 ok; 1 a model failure (Spheral error, accounts not closing, a slowness or
  monitor failure); 2 bad arguments, a missing input, or Spheral not available.
- **Thresholds.** Each from an analytic expectation stated in its task or from a measured run rounded up to one
  significant figure, written in the test with a comment naming its source, changed only together with the README and
  `docs/spheral_frag_assumptions.md`. No silent tolerance changes.
- **Machine.** Runs on more than 4 processes check the Linux VM's load as `bench.py` does and refuse above 2.0 unless
  `--force`; long runs (Task 3's 3D segment) run one at a time; another session may be using the machine.
- **Nothing for the deferred sink is built in MC**, not even as a disabled option: no `--sink` flag, no `m_phys`, no
  sink brackets. Task 2 leaves a named, empty sub-boundary hook list and nothing else.
- **Commits on `main`,** one per task, scope prefix and a sentence describing the behaviour, ending with the
  attribution line; `git add` only the task's files (other sessions commit concurrently). Push only when Asha says so.

## Review focus

1. **Every gram is booked once.** In the replay, centre crossings go to the film account exactly once with body.npz's
   mass and the held temperature's enthalpy (M2); a window removes nothing, its mass constant to round-off, and its
   record checks book the body's masses per ID, not Spheral's last-bit-moved masses (fact 5). → Tasks 2, 3.
2. **The frames' provenance is untouched.** Task 1's rerun is byte-identical in every vtu, vtp and the history before
   the side file is copied; otherwise D7's fallback and a report. → Task 1.
3. **The closure is honest.** The release thresholds move to the closed sum only with the measured residual, and the
   positive excess's cause is found or stated as not found with its numbers. → Task 1.
4. **Every ID is accounted for.** The set of IDs at every check is the same on every process count, each removed ID
   appears in exactly one removal row, and checkpoints read across process counts bitwise. → Task 3.
5. **The window core does not pre-empt M3 or M4.** It holds the frame's temperatures and loads, has no triggers or end
   conditions beyond its length, and no softened or free-energy configuration; its sub-boundary hooks are empty. → Task 2.
6. **Nothing of the deferred sink leaks into MC** (Global constraints), and every replay name is unchanged. → Task 2.

## Plan decisions

D1, D7, D11 and D13 apply to MC. D2–D6, D8–D10 and D12 belong to the deferred sink and are kept, with their text, in
"Deferred: sub-grid Girin spraying"; their numbers are kept so that the facts and the questions still point at them.

- **D1, the scheme** (M2 decision 8, an amendment above): SolidSPH only. Check 14's ID part in RZ at 2.2 and 1.1 mm,
  3D at 2.2 mm.
- **D7, the release on removed faces: fact 107 (a)'s option (4), written beside the frames** (Q1, accepted). A
  write-only export addition 07 in the reconstruction's prototype records, per macro step, every face that step's
  deaths removed with the release it had before removal (centroid, `n_derived`, area, θ, mass, film temperature), and
  the spray step's own Σ dm, in one side file per run (`release_removed.npz`), leaving every `field_<k>.vtu`,
  `surface_<k>.vtp` and the history byte-identical. The flight is rerun, compared byte for byte with the existing one,
  and the side file copied in (Task 1). Why (4): option (1) puts the mass on faces of other areas and inclinations,
  chained through the cascade; (2) gives the total but no inclination; the frames alone cannot reconstruct it (fact
  12); (4) is a superset of (2) at the same cost (one rerun), and the deferred release tables (S1) need its
  inclinations. Fallback if the rerun is not byte-identical, or for flights without the side file: option (2) from
  the history itself, the per-step difference between its sprayed increment and Σ `release_rate` · A
  (`--removed-release history`, labelled), which needs no rerun.
- **D11, the replay keeps centre-crossing removal** (spec §14.1, "Open within MC"; Q5, accepted). Fact 13: on this
  flight the spray rate and the recession agree within a few percent, so centre crossing in the replay removes the
  sprayed material at the right rate on average. With the sink deferred it is the only recession in this version.
- **D13, run-name brackets of the window mode** (defaults in brackets, suffixes only when different): a window's name
  is its source's name plus `__window-k<K>-<length>ms`, the source being the replay run it switches from
  (`--from-replay`) or the body it starts stress-free from (`--frame` alone, marked `-free`), so windows from
  different replays never collide; then `sub_interval` (Task 2's measured default) and M2's `loads`, `pressure` and
  `area_shares` brackets as the replay emits them. No replay name changes.

---

### Task 1: The release on removed faces — export addition 07 and the history's closure

**Files:**
- Modify: `prototype/work-2026-10-08-spheral-mvp/code/reentry_model/coupled.py` and `body.py` (addition 07,
  write-only), its tests; `prototype/work-2026-10-08-spheral-mvp/README.md` (addition 07, the new package SHA-256);
  `prototype/work-2026-10-08-spheral-mvp/rebuild/post/` (07 as a diff, as 01–06)
- Modify: `spheral_frag/contract.py` (an optional run-level item `release_removed`), `analysis/spheral_m1_flight.py`
  (columns `release_removed_kg`, `release_spray_step_kg`), `data/spheral/m1_frames.csv`,
  `tests/test_spheral_frag_flight.py` (the release thresholds)

Unchanged by the rescoping (Q1): the closure belongs to the re-entry model's film and mass account. If Task 3 or
M2 Task 10 is committing in `analysis/`, `data/` or `tests/` at the time, the controller sequences Step 4's commit
after theirs.

- [ ] **Step 1: Addition 07** (D7). Where a step's deaths remove surface faces (`MeltingBody._hand_over` and the
  cascade's passes), record for each removed face the release it carried in that step before removal: centroid,
  `n_derived` (as 04 writes it), area, θ against `v_hat_body`, mass released (kg), film temperature; and the spray
  step's own Σ dm over all faces. Write them per run as `release_removed.npz` (arrays concatenated over steps, with a
  frame index per row, plus the per-step scalar), nothing in the vtu, vtp or history. Test in the prototype: a hot
  6 s flight with deaths every step, run with and without 07, gives byte-identical vtu, vtp and history; Σ over a
  step of (frame faces' release · A + removed faces' mass) equals the step's Σ dm to round-off.
- [ ] **Step 2: Rerun** the reconstruction's 100 mm Scheil flight with 07 (its recorded flags, seed 12345, frames every
  step) into a scratch output directory; compare the history CSV and every `field_<k>.vtu` and `surface_<k>.vtp`
  byte for byte with `reentry_model_output/<the M2 flight>/`. Identical: copy `release_removed.npz` into the existing
  run directory and record the new package SHA-256 beside the old one (the frames' provenance is unchanged). Not
  identical: stop and report (D7's fallback then applies and the frames stay as they are).
- [ ] **Step 3: The positive excess** (fact 11: 13 frames, at most +3.8 %). Compare per step the frame's Σ
  `release_rate` · A, the spray step's own Σ dm and the history's sprayed increment; test the candidate that the
  frame's patch areas (at write time, after the step's recession moved nodes) differ from the areas the spray step
  divided by. Record the cause found, or that it is not found, with the numbers.
- [ ] **Step 4: Closure.** Per spraying frame, (Σ `release_rate` · A + Σ removed) / sprayed increment − 1, and over the
  flight; the M1 table gains the two columns; `release_step_low`, `release_step_high` and `release_total` become
  thresholds on the closed sum (the measured worst case rounded up), the open-sum ones kept as diagnostics. Commit
  with the README's and the assumptions' release entries.

### Task 2: The window driver's core (was Task 5)

**Files:**
- Modify: `spheral_frag/runner/run.py` (`window` mode), `spheral_frag/naming.py` (mode `window`, D13's brackets),
  `spheral_frag/analyse.py` (window runs: their record checks and accounts)
- Create: `tests/test_spheral_frag_runner_window.py` (marker `spheral`)

Depends on M2 Task 11 (the `switch` conversion and the first real-time run); if Task 11 is not committed when this
task starts, the controller asks before building the conversion here. What it keeps from the old Task 5 is what M3's
separation windows and M4's slurry and detachment windows need from every window: a start from a replay checkpoint or
a stress-free body, real time, the frame's loads and temperatures held, restarts at sub-boundaries with checkpoints
and resume, the record and the accounts. What it drops went with the sink (`--sink`, `--floor`, `--sink-h`,
`--release-*`, `m_phys`, floor deletions), with M4 (`--eos murnaghan`, `--energy free`, heating), or with M3 (the
fragment finder, the clear and settled tests).

- [ ] **Step 1: The loop.** `run.py window --body DIR (--from-replay RUN --frame K | --frame K) --length S
  [--sub-interval S] [--loads on|off] [--pressure P] [--area-shares patch|volume] [--out DIR] [--resume]`. Start: a
  replay checkpoint converted by M2 Task 11's switch (masses, densities and free densities divided by s, velocities
  zero, u unchanged), or the body's particles stress-free at rest. Real time from the frame's time, s = 1,
  temperatures held at frame k (`ImposedTemperature`). Each sub-interval: a fresh body from the last checkpoint (the
  restart, M2 Task 1's proven route; the restart is also where particles are redistributed, fact 9), packages
  (imposed temperature, damping off, loads with frame k's table held or a prescribed uniform pressure, the
  deceleration, Hold for held particles as the replay), advance, then at the sub-boundary: the sub-boundary hooks (a
  named list, empty in MC: M3's fragment finder and end tests and, later, the deferred floor deletions act there), a
  record check booking the body's masses per ID (as `Replay.write_check`), a row in `subintervals.csv` (KE, max |v|,
  max surface stress, the stepping cost, the restart's own cost, and the recession the window omits: frame k's
  Σ `release_rate` · A × elapsed time / (ρ_l × the frame's surface area), in spacings — new question 1), the
  checkpoint. `summary.json` last.
- [ ] **Step 2: Accounts.** `analyse` reads window runs: the record's per-ID masses, the switch's conversion
  recorded in `meta.json`; it closes the accounts to 1e-12 m₀ (the replay's centre crossings before the switch carried
  in from the source run, none added by the window) and reports the film account by reason.
- [ ] **Step 3: Tests:** a window on the RZ 2.2 mm body of frame 0 (stress-free, frame 182's loads held, as check 6)
  and one switched from a replay checkpoint, stopped after 3 of 6 sub-intervals and resumed, equals the uninterrupted
  run bitwise; the body's mass per ID is constant through the window and the accounts close to 1e-14 m₀; the run name
  parses back; a replay's name is unchanged by the new brackets. Measure a restart's cost at RZ 2.2 and 1.1 mm and 3D
  2.2 mm (18 processes); the default sub-interval is the shortest at which a restart costs at most 5 % of the
  sub-interval's stepping, rounded up to one significant figure, recorded with the numbers. Commit.

### Task 3: The ID part of check 14 (was Task 8)

**Files:**
- Create: `tests/test_spheral_frag_runner_check14.py` (marker `spheral`), `analysis/spheral_mc_check14.py`,
  `data/spheral/mc_check14.csv`
- Modify: `spheral_frag/runner/run.py` only if a defect is found

The old Task 8 passed the IDs through sink steps and floor deletions on a slab; with the sink deferred, the deletions
are the replay's centre crossings and the window's restarts, on the flight's own bodies.

- [ ] **Step 1: Redistribution** (fact 9) as a test: on 4 and 18 processes, after unbalancing deletions,
  Spheral's Peano–Hilbert redistribution carries `pid`, `x_ref`, the mass and every per-particle field the
  checkpoints hold, bitwise.
- [ ] **Step 2: Through the replay, the switch and a window.** RZ 2.2 and 1.1 mm: the MVP flight replayed from k = 0
  to 100 (centre crossings begin once the film appears at k ≈ 48; RZ 2.2 mm about 11 min at M2's 6.6 s per
  interval) on 1, 4 and 18 processes; 3D 2.2 mm: from M2 Task 10's k = 47 checkpoint to k = 100 on 4 and 18 processes
  (about 2.6 h each at 178 s per interval, one at a time; 1 process would take about 18 times longer). Each then
  switched at k = 100 into a 5 ms window of 5 sub-intervals (Task 2). At every check: the set of IDs is the same for
  every process count; each removed ID appears in exactly one removal row, reason `centre_crossing`, at the same
  frame on every process count; the masses per ID agree to round-off across process counts (not bitwise: the
  reductions differ); a checkpoint written on N processes is read on M (all pairs of 1, 4, 18; 3D 4 and 18) bitwise.
  The number of crossings per segment is recorded in `data/spheral/mc_check14.csv`.
- [ ] **Step 3:** In-run redistribution (`redistributeStep`) stays off: the packages refuse a changed local count
  (fact 10), and every sub-boundary and frame boundary redistributes through the restart. Recorded with the numbers.
  Commit.

### Task 4: Documentation and the spec amendment (was Task 10)

**Files:**
- Modify: `README.md` (the large-fragment section: the window command, MC's verification table — check 14's ID
  part, the release closure, the restart cost), `docs/spheral_frag_assumptions.md` (centre crossing as the only
  recession in this version; the window's held temperatures and loads and the omitted recession; the release
  closure; each marked verified, analytic, measured or assumed), the spec (every item of "Spec amendments these
  imply" above, MC as built and measured, D1's amendment, D7's decision), this plan's "Decisions at review" (the
  answers to the new questions).
- [ ] One commit, `docs: Spheral MC as built and measured: ...`.

## Order and parallelism

1. **Task 1** (the prototype and drama_env) can start at once; it needs no Spheral. Its Step 4 touches `analysis/`,
   `data/` and `tests/`, so its commit is sequenced after any M2 commit in those directories.
2. **Task 2** after M2 Task 11.
3. **Task 3's Step 1** at once; its Step 2 after Task 2 and after M2 Task 10 has kept its 3D k = 47 checkpoint.
4. **Task 4** last.

Wall-time estimate: Task 1's flight rerun about an hour; Task 2's tests and restart measurements under an hour; Task
3's RZ replays tens of minutes per process count, its 3D segments about 5 h in all.

## Questions for Asha

### New questions (2026-10-11, after the rescoping)

1. **Recession in windows.** With the sink deferred, a window removes no mass, while the finite element keeps
   receding. In a stiff window (≤ 50 ms) the omitted recession is at most 0.12 mm (fact 14's largest area-mean rate,
   0.05 of a 2.2 mm spacing); in a soft slurry window up to 1 s it reaches about 2.4 mm, about one spacing, and in a
   detachment window it is the neck's own thinning. Options: (a) remove nothing in windows, record the omitted
   recession per sub-interval (Task 2) and flag a window where it exceeds half a spacing; (b) centre crossing in
   windows against the frames' surfaces, interpolated in time — sound only while particles stay near their
   reference positions, which a slurry or separating window does not; (c) keep long windows short by chaining them
   from the replay. **Recommended: (a) in MC**, and the choice for long windows made in M4's plan with the measured
   omission in hand; necking itself is followed in the replay, where centre crossing thins the neck at the right rate
   (fact 13).
2. **Where the later sub-grid Girin spraying lives, and which rate drives it.** Two readings of "later implement Girin
   spraying as a sub-grid parameterization": (a) the deferred sink driven by the finite element's own Girin release
   (the frames' release tables, S1–S6 as written); (b) Spheral's own Girin closure on a film tracked below each
   surface particle (spec §13.6, M9; §19's row "droplet spraying as a sub-grid model in Spheral" already defers that
   to M9). Both need the same mechanics — mass leaving particles continuously, integrated density, both RZ masses,
   the floor, check 15 — and differ only in the rate's source. **Recommended:** a lettered milestone of its own ("MS,
   sub-grid spraying"), planned when the scoping verdict (M4) shows that spray-driven mass loss changes the
   fragments, starting from the deferred section with reading (a), the finite element's verified Girin release; (b)
   stays M9's and reuses MS's sink. Spec §14 gains the row (Task 4).

### Questions as asked on 2026-10-11 (answered; see "Decisions at review")

1. **The release on removed faces** (D7, fact 107 (a), M1 decision 16). Recommended: option (4), a write-only side
   file (`release_removed.npz`) from export addition 07, with the reconstruction's flight rerun and compared byte for
   byte, so the frames and every M2 run keep their provenance; the history-derived option (2) as the fallback and a
   bracket. The frames alone cannot reconstruct the removed release (fact 12). Or: option (2) only, no rerun.
   *Answered: yes (option 4).*
2. **Which particles sink** (D5). Recommended: every exposed surface particle in an inclination bin with release,
   dropping §10's and §14.1's "above the liquidus" — at a 2.2 mm particle's centre only 0.3–41 % of the sprayed mass
   comes from material above the liquidus (fact 15). Amends §10 and §14.1. *Answered with the scope: no spray sink in
   this version (Decisions at review 2); D5 open when the deferred work is taken up.*
3. **The mapping and its normalisation** (D6). Recommended: by inclination as the spec says, normalised per step to the
   frame's total release (a labelled device, factor recorded; S1 measures how far it is from 1); per-patch ownership
   as a bracket. *Answered: yes; deferred.*
4. **The area shares under the sink** (D4). Recommended: build the `csf` shares and choose between them and the patch
   shares in windows by check 6 (S3) and check 15 (S5); the replay keeps the patch shares (M2 decision 9).
   *Answered: yes; deferred.*
5. **The replay keeps centre-crossing removal** (D11, spec §14.1's open item). Recommended; the spray and the recession
   agree within a few percent on this flight (fact 13). *Answered: yes.*
6. **Check 15's extensions** (S5). Recommended: the top pressure on a held cup, the shares compared, the flight's
   own rates, and a mushy slab with free and imposed energy, beside the spec's core matrix (32 runs); and the softened
   slab without viscous rows (D10), since those are M4's. *Answered: yes; deferred.*
7. **Floor removals at sub-boundaries** (D8): the sink stops at the floor and the particle is deleted at the next
   restart, at most one sub-interval late. Recommended. *Answered: yes; deferred.*

---

## Deferred: sub-grid Girin spraying (later milestone)

> **Not built in this version of Spheral** (Asha, 2026-10-11: *"this early version of spheral is not being used to lose
> mass through liquid spraying. We want to capture large fragmentation and necking processes. We will later implement
> Girin spraying as a sub-grid parameterization."*). Kept as measured and designed on 2026-10-11 so that the later
> milestone (new question 2) starts from it. Q3, Q4, Q6 and Q7 were accepted and apply when it is taken up; D5 (Q2)
> is open again then. Facts 1–8 and 14–16 are its measured basis. Task references below use the deferred numbering
> S1–S6; "Task 1" and "Task 2" are MC's (the release side file and the window driver, which the sink extends with its
> brackets, its removal reasons and its floor deletions at the sub-boundary hook). Before it is executed it is
> re-planned against the code as it then stands (the window driver, M3's and M4's packages).

**Its goal, as written for MC:** mass that the finite-element model sprays leaves Spheral's particles continuously, step
by step, through a spray-sink package, without false pressure, motion or interpolation error; every sink step, floor
removal and centre-crossing removal is booked with its mass and enthalpy into the film account. **Done when:** check 15
(spec §12, with Q6's extensions) passes with SolidSPH in RZ and 3D at 2.2 and 1.1 mm; the film account's sink share
closes against the finite element's sprayed mass within Task 1's residual. **Answers:** whether mass can leave Spheral
continuously without false pressure, motion or interpolation error; the deletion floor (a quarter or a half of the
starting mass); whether the smoothing length shrinks with the mass. **If it fails:** false tension or motion that
neither the floor, the smoothing-length rule nor the area-share rule removes → whole-particle deletion (spec §19), with
fragment masses reported net of each particle's owed mass; S5 decides this.

**Its review focus:** the sink's decrements and the floor's remainders booked once, with D9's enthalpy, from the sink's
own mass field (D12); the release tables closing against the history per step (Task 1's residual) and nothing sinking
in the replay (D11), so spray is never counted twice; RZ moving m and massRZ by one factor (fact 4); no
mass-proportional quantity left behind (loads, deceleration, damping, shares see the current mass); nothing tuned to
pass check 15; every device a setting with a recorded default and bracket.

### Its decisions

- **D2, the mass changes in the sink package's `finalize`** (facts 1–3): after every step each sinking particle's
  mass drops by its rate times the step, in RZ m and massRZ by the same factor (fact 4); integrated density keeps
  its density, so it shrinks in volume. Not an update policy: Spheral's increment policy cannot address the mass
  (fact 3), and a Python update policy is untried and costs per call (spec §2). The mass is constant within a step
  (first order in dt, at most a few 1e-5 of a particle's mass per step at the rates of fact 14).
- **D3, the smoothing length: held by default, shrunk as the bracket.** Held: Spheral's IdealH as it is (fact 7).
  Shrunk: the sink rewrites the sinking particles' ideal H so that h = h_ideal,unshrunk (m/m₀)^(1/dim), where the
  unshrunk ideal is recovered each step by dividing out the factor applied the step before (fact 8's compounding);
  checked on the slab to 1 % of the intended ratio. Spec §14.1 keeps whichever gives the smaller surface error in
  check 15; held is the default until then because shrinking shortens the step and leaves the particle fewer
  neighbours.
- **D4, the area shares under the sink** (fact 16; Q4 accepted). A sinking particle keeps internal forces
  proportional to its mass (SPH pair forces carry m_i m_j), so a fixed patch share (M2 decision 9) overloads it as it
  empties: its pressure share p A / m grows while what it can pass on shrinks. A third share rule, `csf`: aᵢ = |∇c|ᵢ
  Vᵢ with the current V = m/ρ, scaled by one factor so that the shares sum to the frame's surface area, recomputed
  every `--surface-every` steps within a sub-interval (default 50; S3 measures the cost). Check 6 is rerun with it
  (S3) and check 15 compares `patch` and `csf` under the sink (S5). Recommended: `csf` in windows with the sink if
  check 6's errors with it stay within 1.5 times the patch shares'; the replay keeps `patch`.
- **D5, which particles sink: the frame's release decides, not the particle's temperature** (fact 15; would amend §10
  and §14.1's "above the liquidus"). Every exposed surface particle (`first_body` exposure, M2 Task 4) whose
  inclination falls in a bin with release sinks at that bin's rate. Spray comes from a skin thinner than a particle,
  so the particle's own temperature says nothing about it; Spheral never decides what is sprayed (§10). *Open: Q2 was
  answered with the scope, not on this point.*
- **D6, the release by inclination, normalised per step** (Q3 accepted). Prepare-side release tables per frame k: for
  each 1° bin of θ the released mass per step (the frame's patches plus the removed faces of D7) over the frame's patch
  area in the bin, as kg/m²/s, linear in time between frames k and k + 1 as the loads. The particle rate is
  qᵢ = r(θᵢ) aᵢ. Because the particles' shares by bin need not equal the frame's areas by bin, Σ qᵢ differs from the
  frame's total; with `--release-normalise on` (default, a labelled device) the rates are scaled by one factor per
  step so that Σ qᵢ equals the frame's total release rate, the factor recorded. Per-patch ownership (each patch's
  release to the particle that owns its area) is a bracket, `--release-map patch`; it is exact per step but spiky
  (fact 14: a 2.2 mm particle under the median per-frame maximum, 57 kg/m²/s, empties in about 0.1 s).
- **D8, floor removals at sub-boundaries** (Q7 accepted). Below the floor (default ¼ of the particle's starting mass,
  bracket ½) the sink stops lowering a particle and books what it still owes; at the next sub-boundary (a full
  restart, Task 2's hook) the particle is deleted and its remaining mass goes to the film account with its own
  enthalpy (reason floor), the owed part with the sprayed enthalpy (reason sink). The delay is at most one
  sub-interval (Task 2 measures a restart's cost). Live deletion with rebuilt packages (fact 10) is the alternative if
  restarts prove too expensive.
- **D9, the energy of what leaves.** Removed mass carries h(T_film) of its θ bin (the frame's mass-weighted film
  temperature of the releasing patches and removed faces, through §8's table: spec §14.1 condition 4). With imposed
  temperatures (the replay's `ImposedTemperature`, the windows' default) the particle's u stays imposed and only the
  accounts see h_spray; with `--energy free` the particle keeps the rest of its energy, u' = (m u − dm h_spray) /
  (m − dm), and the slab's energy closes (check 15).
- **D10, check 15's configurations.** Stiff: the thermal EOS at real stiffness (K held at its solidus value, 70.3 GPa,
  M1 decision 2). Softened: Spheral's `MurnaghanEquationOfState` with exponent 1, K = 30 MPa and ρ₀ = ρ_free(T) of the
  slab (spec §8.2; its sound speed checked correct in M0), the thermal strength (no shear above 50 % liquid), no
  viscous rows (M4's). Labelled: the softened slab is inviscid but for Spheral's artificial viscosity.
- **D12, a physical-mass field.** The sink keeps a registered ScalarField `m_phys` per particle (kg, the ring's in RZ),
  lowered by exactly the booked decrements; Spheral's mass is set from it (times the scale factor, 1 in windows). Record
  checks and accounts read `m_phys`, so they close exactly despite fact 5.
- **Its run-name brackets** (D13's part for the sink; defaults in brackets): `sink` (off; on in windows that spray),
  `floor` (0.25), `sink_h` (held | shrink), `release_map` (theta | patch), `release_normalise` (on),
  `removed_release` (side | history | none), `energy` (imposed | free), `eos` (thermal | murnaghan:<Pa>), and the
  slab's geometry name.

### S1: Release tables by inclination (was Task 2; core, drama_env)

**Files:**
- Create: `spheral_frag/release.py` (numpy only at module level), the subcommand `spheral_frag release`,
  `tests/test_spheral_frag_release.py`, `data/spheral/mc_release.csv`
- Modify: `spheral_frag/__main__.py`, `spheral_frag/naming.py` (the release brackets)

- [ ] **Step 1: The tables** (D6, D7):

  ```python
  def release_table(frame, removed, bins_deg=1.0, dt_step=0.5) -> dict   # per bin: rate_kgm2s, area_m2, mass_kg,
                                                                          # removed_kg, T_film_K, h_spray_Jkg
  def write_tables(prepared_dir, fe_run=None, mode="side") -> str        # <prepared>/release.npz + release.json
  def evaluate(tables, k, k_next, w, theta_deg) -> (rate, h_spray)       # linear in time, as loads.evaluate
  ```

  θ of the frame's patches as the load tables bin them (`patch_theta`); removed faces by their exported θ; a bin
  with removed mass but no frame area passes its mass to the nearest bin with area (counted). `mode="history"` (D7's
  fallback): the per-step deficit against the history spread in proportion to the bins. `release.json` records the
  inputs' SHA-256 (prepare.json, the side file), the mode and the per-frame closure. The prepared frames and
  `prepare.json` are not touched.
- [ ] **Step 2: Closure on the flight.** Σ over bins of mass equals the frame's Σ `release_rate` · A + Σ removed to
  1e-12 (exact by construction) and the history's sprayed increment to Task 1's residual; per frame into
  `data/spheral/mc_release.csv` with the area-mean and largest rates (fact 14's numbers in both modes).
- [ ] **Step 3: The mapping mismatch** (D6's normalisation). On bodies of the flight at k = 100, 150 and 182 (the M2
  body stage, RZ 2.2 and 1.1 mm, 3D 2.2 mm), with the surface, normals and `patch` and `csf` shares computed as the
  loads package does (numpy, `packages/loads.py` functions): Σᵢ r(θᵢ) aᵢ against the frame's total. Record the factor
  the normalisation would apply; it sets whether `--release-normalise on` stays the default.
- [ ] **Step 4: Tests** (synthetic frames, `tests/spheral_frag_synthetic.py`): a uniform release on a sphere gives a
  flat table; removed faces land in their bins; the closure is exact; the history mode reproduces the history's total
  per step; `evaluate` is bitwise at frame times. Commit.

### S2: The sink package and the physical-mass field (was Task 3)

**Files:**
- Create: `spheral_frag/runner/packages/sink.py`, `tests/test_spheral_frag_runner_sink.py` (marker `spheral`)
- Modify: `spheral_frag/runner/packages/__init__.py` (docstring), `spheral_frag/runner/nodes.py` (`m_phys` and
  `m_start` in the state and the checkpoints, schema version up; old checkpoints read with m_phys = m / s)

- [ ] **Step 1: `SpraySink`** (D2, D5, D6, D8, D9, D12). Built per sub-interval like the other packages, from the
  loads package's surface, exposure, inclinations and shares and the release tables of frames k and k + 1 (or a
  synthetic uniform table for the slab, S4). In `finalize` after every step: qᵢ from the table at t (normalised
  per D6), dmᵢ = min(qᵢ dt, m_physᵢ − floor · m_startᵢ); m_phys −= dm; Spheral's m = s m_phys, in RZ massRZ scaled by
  the same ratio; u per D9; per-particle accumulators of booked mass, enthalpy (dm h_spray) and owed mass (the part
  the floor held back). `dt()` votes no limit (asserted: qᵢ dt ≪ mᵢ). Exposes `booked()` for the driver.
- [ ] **Step 2: The smoothing-length bracket** (D3). `sink_h="shrink"`: in `evaluateDerivatives` (appended after
  the hydro) the sinking particles' "new H" is multiplied by (m_startᵢ / m_physᵢ)^(1/dim) divided by the factor
  applied at the previous step, so the factor does not compound (fact 8).
- [ ] **Step 3: Tests** (`spheral` marker; the probe slab of the facts, 1 and 4 processes, RZ and 3D): Σ m_phys +
  booked = initial to 1e-15 of the total; Spheral's density of sinking particles unchanged to 1e-14; a stress-free
  slab stays at rest; RZ with m only (a test-only switch) reproduces fact 4's failure, so the test catches the
  omission; `held` leaves H bitwise equal to a run without the sink; `shrink` gives h / h_held within 1 % of
  (m/m_start)^(1/dim) after 200 steps; the floor stops the mass at floor · m_start and books the owed rest; `free`
  energy closes Σ m u + Σ booked enthalpy to round-off on a slab at rest. Measure the package's cost per step (3D
  2.2 mm, 18 processes; RZ 1.1 mm). Then extend Task 2's window driver: `--sink`, `--floor`, `--sink-h`,
  `--release-map`, `--release-normalise`, `--energy` and their brackets; floor deletions at the sub-boundary hook;
  the sink's bookings as `record.Removal` rows (reasons sink and floor; one per particle per sub-interval); record
  checks and `analyse` from `m_phys`; `subintervals.csv` gains sink and floor mass and enthalpy and the normalisation
  factor; Task 2's resume test repeated with a synthetic release. Commit.

### S3: CSF area shares, re-detected within a sub-interval; check 6 with them (was Task 4)

**Files:**
- Modify: `spheral_frag/runner/packages/loads.py` (`area_shares="csf"`, `surface_every`), `spheral_frag/naming.py`
  (the `area_shares` bracket gains `csf`), `tests/test_spheral_frag_runner_packages.py`
- Modify: `analysis/spheral_m2_check6.py` (a csf row set), `data/spheral/m2_check6.csv`

**A candidate for M3 and M4, independent of spraying (reasoning, 2026-10-11).** M2 decision 9's patch shares give each
frame-k patch's area to the surface particle nearest it *by reference position*. That holds while the particle surface
is the frame's: the replay, and MC's windows, which hold frame k's loads on frame k's shape for a short time and never
change topology. In M3's separation and M4's detachment windows a piece opens a new face, moves and turns; in a slurry
window the surface flows; necking narrows a section the frames may describe only coarsely. The reference positions
then no longer say which particle faces the flow, so loads in those windows need the surface, normals, exposure and θ
re-detected from current positions (this task's `surface_every`) and shares that follow the current surface — which
the csf rule does (fact 16: within −4 to +5 % by band, the load passing inward as a layer thins). Two parts of this
task would change for that use: the shares should be calibrated by the measured constant (fact 16's 0.459–0.485 of
the true area at each spacing), not normalised to the frame's area, because a splitting body's area grows by its new
faces; and whether gas loads act on a fresh crack face at all is a modelling choice M3 must make. Check 6 on an
undeformed sphere cannot choose a rule for moving surfaces, so nothing of this task is built in MC; M3's plan should
consider the re-detection and the csf rule with a necking or splitting test case to choose by.

- [ ] **Step 1: The rule** (D4, fact 16). aᵢ = |∇c|ᵢ Vᵢ with Vᵢ = m_physᵢ / ρᵢ, one factor so that Σ aᵢ over the
  surface particles equals the frame's area; with `surface_every` N > 0 the surface, normals, exposure, θ and shares
  are recomputed from the current positions and volumes every N steps (gathered as at construction), the forces
  blended in time as now. N = 0 keeps M2's once-per-interval behaviour (the replay's default).
- [ ] **Step 2: Tests:** on the 2.2 mm lattice sphere the csf shares sum to 4πR² exactly after calibration and their
  30° bands lie within fact 16's −4 to +5 %; halving the outer layer's volumes moves share inward and keeps the total
  within 5 %; the recomputation leaves an unchanged body's forces bitwise equal; cost per recomputation at 3D 2.2 mm
  (18 processes) and RZ 1.1 mm.
- [ ] **Step 3: Check 6 with csf shares** (RZ 2.2/1.1 mm, 3D 3.0/2.2 mm, as M2 Task 8): drag against the table's
  integral and the interior and core errors against the finite-element reference, beside the patch and volume rows.
  Commit with the rows; D4's recommendation is confirmed or withdrawn on these numbers.

### S4: Slabs, the synthetic release and the softened configuration (was Task 6)

**Files:**
- Modify: `spheral_frag/body.py` (`--geometry slab`), `spheral_frag/runner/eos.py` (`make_eos` for
  `murnaghan:<Pa>`), `spheral_frag/release.py` (a synthetic uniform table), their tests

Step 3 (the softened equation of state) is M4's in any case (its slurry configuration); if M4 is planned first it
builds it there and S4 reuses it. The slab of Step 1 may also serve M4's check 9 (a sheared or draining layer).

- [ ] **Step 1: The slab** (check 15). A layer of thickness H on a held cup: 3D a W × W × H block (W = 20 mm,
  H = 10 mm) whose bottom two layers and side walls two layers thick are held (`hold`, as check 5's bar); RZ a disk
  of radius W/2 on the axis with the same base and rim held. Cell-centred lattice, no frame clipping, uniform
  temperature (`--synthetic uniform`, the slab's T a body option); named `_geometry-slab-<W>x<H>mm`.
- [ ] **Step 2: The synthetic release:** a uniform rate R (kg/m²/s) and film temperature on every exposed particle
  facing +x (θ < 45°), as a one-frame table, for the slab only.
- [ ] **Step 3: The softened EOS** (D10): Spheral's Murnaghan, exponent 1, K from the option, ρ₀ = ρ_free(T_slab);
  the thermal strength unchanged; its sound speed checked against √(K/ρ₀) both sides of ρ₀ (M0's check). Tests;
  commit.

### S5: Check 15 — the spray sink on a slab (was Task 7)

**Files:**
- Create: `analysis/spheral_mc_check15.py`, `data/spheral/mc_check15.csv`, `tests/test_spheral_frag_runner_check15.py`
  (marker `spheral`)

- [ ] **Step 1: The cases** (Q6 accepted). The core matrix of spec §12 check 15: forms RZ and 3D; dx 2.2 and 1.1 mm;
  stiff (liquid at 920 K, real stiffness, length 2 ms) and softened (D10, length 0.2 s); floors ¼ and ½; `sink_h`
  held and shrink — 32 runs, each with the top pressure on (16.6 kPa, the flight's peak dynamic pressure, spec §2) and
  the shares `csf`. The rate R recedes the surface by 2 spacings of the coarser dx over the run: R = 2 ρ_l · 2.2 mm /
  length (5,280 kg/m²/s stiff, 52.8 soft; a 2.2 mm particle then empties in about 2,400 (stiff) and 5,000 (soft) of
  its own sound-crossing times, so the loss stays quasi-static). Then at dx 2.2 mm, RZ and 3D, held h and the default
  floor: the shares `patch` against `csf`; the pressure off; the flight's own rates (fact 14: 5.7 kg/m²/s soft over
  1 s, 355 stiff); a mushy slab (895 K, half liquid, spraying at a 920 K film: D9's enthalpy difference) with
  `--energy free` and `imposed`.
- [ ] **Step 2: Measure** per run (spec §12 check 15): mass closure (remaining + booked − initial) and energy closure
  (free energy: Σ m (u + v²/2) + booked enthalpy − initial − the load's work); the stress in the surface layer and the
  layer below against the applied pressure (false tension or compression, in units of p and of the elastic scale);
  velocities against the noise of the same run without the sink; the surface position where the kernel-estimated
  density (Σ m_j W, diagnostic only) falls to half its interior value, against R t / ρ_l; the neighbour count and h of
  sinking particles; the floor removals' timing; cost per step.
- [ ] **Step 3: Decide and record** (spec §14.1): the floor (¼ or ½), the smoothing-length rule (the smaller surface
  error), the share rule (D4), each with the alternatives' numbers; thresholds from the runs. **If the sink fails** —
  false tension or motion that neither floor, h rule nor share rule removes — stop and report to Asha (spec §14.1, the
  fallback of §19). Commit with the table and plots in `spheral_output/mc/check15/`. Then repeat Task 3's Step 2
  window with the sink and floor on, so that the IDs are also checked through sink steps and floor deletions.

### S6: The sink on the flight's frames (was Task 9)

**Files:**
- Create: `analysis/spheral_mc_flight.py`, `data/spheral/mc_flight.csv`

- [ ] **Step 1: Replay** the MVP flight's RZ 2.2 mm body from k = 0 to 182 at the defaults (about 20 min at M2's 6.6 s
  per interval), and its 3D 2.2 mm body to k = 100 (about 5 h at 178 s per interval), keeping their checkpoints
  (Task 3's replays to k = 100 can be continued).
- [ ] **Step 2: Windows** at k = 100, 150 and 182 (RZ 2.2 and 1.1 mm via a 1.1 mm replay to 182, 3D 2.2 mm at 100):
  stiff, 5 ms, the frame's loads held, deceleration on, the sink from S1's tables with S5's chosen floor, h rule
  and shares; and one softened window at k = 150 (RZ 2.2 mm, 0.2 s).
- [ ] **Step 3: Measure:** the sink's total against the frame's release rate times the window length (normalised: equal
  by construction, the factor recorded; unnormalised: the mismatch of S1 Step 3), the film account by reason,
  `analyse`'s closure, the stress change at the surface against the window without the sink, and the switch's
  balance (M2 Task 11's measure) with the sink on. Commit.

### Its order and cost (as estimated for MC)

S1 alongside S2 and S3 (S1's Step 3 needs S3's csf shares as numpy functions); S2 after Task 2; S4 after S2 (its Step
3 at once); S5 after S3 and S4, deciding whether the sink works; S6 after S1 and S5. The stiff slab runs 2 ms at
0.05–0.1 µs per step (20,000–40,000 steps; 3D 1.1 mm about 3,000 particles, minutes to half an hour each), the softened
ones 0.2 s at about 2.5–5 µs (40,000–80,000 steps); check 15's 32 core runs and about 12 variations take a day or two;
S6's replays about 6 h, its windows (stiff 5 ms: 50,000–100,000 steps) hours.

# Spheral MC: Coupling the re-entry model and Spheral at the surface — Implementation Plan

> **For agentic workers:** execute task by task with one sub-agent per task; the controller reviews each task's
> diff and numbers before the next starts. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** mass that the finite-element model sprays leaves Spheral's particles continuously, step by step, through a
spray-sink package, without false pressure, motion or interpolation error; every sink step, floor removal and
centre-crossing removal is booked with its mass and enthalpy into the film account; the release the frames carry
closes against the finite-element history's sprayed mass; and the permanent particle IDs survive deletion, restarts,
the switch to a window and MPI redistribution (spec §14.1).

**Done when** (spec §14, amended below): the ID part of check 14 and check 15 (spec §12) pass with SolidSPH in RZ and
3D at two spacings each (2.2 and 1.1 mm); Spheral's accounts close to round-off through `analyse`; the film account's
sink share closes against the finite element's sprayed mass within the residual Task 1 measures; every acceptance
threshold is measured and recorded in the README, `docs/spheral_frag_assumptions.md` and the tests; MC's measured
facts are written into the spec.

**Answers** (spec §14): whether mass can leave Spheral continuously without false pressure, motion or interpolation
error (check 15); the deletion floor (a quarter or a half of the starting mass); whether the smoothing length shrinks
with the mass (check 15, both ways).

**If it fails** (spec §14, §14.1): false tension or motion that neither the floor, the smoothing-length rule nor the
area-share rule removes → the windows fall back to whole-particle deletion (spec §19), with fragment masses reported
net of each particle's owed mass. Task 7 decides this; the controller stops and reports to Asha before Task 8 if it
fails, and the fallback is planned then, not here.

**Spec:** `docs/superpowers/specs/2026-10-02-spheral-large-fragments-design.md` — §2 (facts: particle identity and
mass, traps), §6.2 (permanent IDs), §8.1–8.2 (energy table, pressure law; the softened configuration), §9.1 (the
replay; centre crossing), §9.3 (windows, the switch, softening rules), §10 (the three zones, the sink as specified),
§11.3 (accounts), §12 checks 14 and 15, **§14.1 (milestone MC)**, §15 (risk "Particle identity and mass"), §17 (the
release gap), §19 (whole-particle deletion as the fallback). Builds on the M2 plan (`2026-10-09-spheral-m2.md`: its
facts, D1–D6 and "Decisions at review" 1–9, cited as *M2 decision n*), M1 decision 16 (the release on removed faces
deferred to MC) and Step 3 shared context fact 107 (a)
(`docs/superpowers/plans/melt-spraying-subplans/00-shared-context.md`). Ordinary task-by-task plan with code shapes
and tests, not generated from a prototype (Asha, 2026-10-01).

**Amendment to spec §14.1 and §12 check 15 (M2 decision 8, Asha 2026-10-11).** The scheme is SolidSPH: SolidCRKSPHRZ
fails the patch test at `116c71f` and SolidCRKSPH in 3D is no better overall (check 7). §14.1's condition 2
("mass-over-density volumes under CRKSPH") therefore does not apply, and condition 1 (integrated density) already holds
in the replay (`run.make_hydro`: `densityUpdate=S.IntegrateDensity` for both schemes). Where spec §12 says "for each
scheme still in contention" and §14 "at two spacings for each scheme still in contention", MC runs **SolidSPH in RZ
and 3D at two spacings each** (2.2 and 1.1 mm). `--scheme crksph` stays a bracket for a later Spheral version, and MC
does not test it.

**Not in MC:** the windows proper — triggers, separation, slurry and detachment windows, chaining (M3, M4); the
viscous slurry rows and the rest of the softened configuration (M4: check 15's softened slab here is inviscid but for
Spheral's artificial viscosity, D10); the heating package (M4); the fragment finder in a run (M3); the hand-back
(M6); the cluster. MC builds the window driver's core (Task 5) on M2 Task 11's switch, which M3 and M4 extend.

## Facts established while writing this plan (2026-10-11)

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
   their own mass, not Spheral's (D12).
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
    100 s: a sink at the spray rate in a window stands in for the finite element's recession within a few percent.
14. **Spray rates per area** (every third frame, k 49–426): the area-mean over a frame's spraying patches has median
    3.3 and maximum 5.7 kg/m²/s (k 196: 2.4 mm/s of liquid at 2,400 kg/m³); a frame's largest single patch rate has
    median 57 and maximum 355 kg/m²/s (k 151). Per patch the release is spiky; by inclination it is smooth.
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
    decision 9's patch shares are exact). When the outer layer's volumes drop to a quarter, the total moves by −3.3 % and the outer layer's share
    falls from 0.378 to 0.096: the load passes to the layer below by itself (D4).
17. **What exists.** `record.py` already has `REASON_SINK` and `REASON_FLOOR` and splits the film account by reason;
    `analyse` closes the accounts to 1e-12 m₀. The replay's record checks book body.npz's masses, not Spheral's
    (`Replay.write_check`), and its removals book body.npz's mass and the held temperature's enthalpy. M2 Tasks 1–9
    are built (decision 8 SolidSPH confirmed, decision 9 patch shares the default, being committed by another
    session); Tasks 10 (check 13), 11 (the switch) and 12 (documentation) are not yet.
18. **The finite-element run is seeded** (seed 12345, spec §2 "Frames as built"); Step 3 fact 105 found coupled runs
    with and without frames identical in every history column once seeded, so a rerun that only adds a write-only
    output can be checked byte for byte (Task 1).

## Global constraints

- **Interpreters.** `PY=/Users/ashajain/miniforge3/envs/drama_env/bin/python` runs the core, the prepare-side release
  tables, the finite-element prototype, the body stage and the analysis. Everything that imports Spheral runs in the
  container through `spheral_frag/container/spheral [-n N]`. Spheral's Python never imports `reentry_model`.
- **Three stages, files only** (spec §4.2). The release tables are a prepare-side product written beside the prepared
  frames, never computed in the runner.
- **The frames** are the reconstruction's MVP frames (M2 fact 8). The prepared frames and `prepare.json` are not
  rewritten; MC adds files beside them (Task 2), so every M2 body and run keeps its provenance.
- **Package layout** as M2: `spheral_frag/runner/packages/sink.py`, `runner/run.py`'s `window` mode, core additions
  `spheral_frag/release.py`; tests `tests/test_spheral_frag_release.py` (drama_env) and
  `tests/test_spheral_frag_runner_{sink,window,check14,check15}.py` (marker `spheral`).
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
  `--force`; long studies (Tasks 7, 9) run one at a time; another session may be using the machine.
- **Commits on `main`,** one per task, scope prefix and a sentence describing the behaviour, ending with the
  attribution line; `git add` only the task's files (other sessions commit concurrently). Push only when Asha says so.

## Review focus

1. **Every gram is booked once.** The sink's decrements, the floor's remainders and the centre crossings each go to
   the film account exactly once, with the enthalpy D9 assigns; Spheral's accounts close to round-off from the
   sink's own mass field (D12), not from Spheral's last-bit-moved masses (fact 5). → Tasks 3, 5, 7, 9.
2. **The sink removes what the finite element sprayed, no more.** The release tables close against the history's
   sprayed mass per step (Task 1's residual); in the replay nothing sinks (centre crossing removes, D11), so spray is
   never counted twice. → Tasks 1, 2, 9.
3. **RZ moves both masses.** m and massRZ change by one factor (fact 4); a test fails otherwise. → Task 3.
4. **No mass-proportional quantity is left behind.** Anything that divides by a mass (loads, deceleration, damping,
   the CSF shares) sees the sinking particle's current mass; the area share of a sinking particle shrinks with it or
   passes to the layer below (D4). → Tasks 3, 4, 7.
5. **Nothing is tuned to pass check 15.** The floor, the smoothing-length rule and the share rule are chosen by the
   measured surface error, with the alternatives' numbers recorded. → Task 7.
6. **Every device is a setting with a recorded default and bracket**: the floor, the smoothing-length rule, the share
   rule, the normalisation of the release table, the sub-interval. → Tasks 2–7.

## Plan decisions (Asha to confirm; see Questions)

- **D1, the scheme** (M2 decision 8, an amendment above): SolidSPH only, RZ and 3D at 2.2 and 1.1 mm.
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
- **D4, the area shares under the sink** (fact 16). A sinking particle keeps internal forces proportional to its
  mass (SPH pair forces carry m_i m_j), so a fixed patch share (M2 decision 9) overloads it as it empties: its
  pressure share p A / m grows while what it can pass on shrinks. MC adds a third share rule, `csf`: aᵢ = |∇c|ᵢ Vᵢ
  with the current V = m/ρ, scaled by one factor so that the shares sum to the frame's surface area, recomputed every
  `--surface-every` steps within a sub-interval (default 50; Task 4 measures the cost). Check 6 is rerun with it
  (Task 4) and check 15 compares `patch` and `csf` under the sink (Task 7). Recommended: `csf` in windows with the
  sink if check 6's errors with it stay within 1.5 times the patch shares'; the replay keeps `patch`.
- **D5, which particles sink: the frame's release decides, not the particle's temperature** (fact 15; amends §10 and
  §14.1's "above the liquidus"). Every exposed surface particle (`first_body` exposure, M2 Task 4) whose inclination
  falls in a bin with release sinks at that bin's rate. Spray comes from a skin thinner than a particle, so the
  particle's own temperature says nothing about it; Spheral never decides what is sprayed (§10).
- **D6, the release by inclination, normalised per step.** Prepare-side release tables per frame k: for each 1° bin of
  θ the released mass per step (the frame's patches plus the removed faces of D7) over the frame's patch area in the
  bin, as kg/m²/s, linear in time between frames k and k + 1 as the loads. The particle rate is qᵢ = r(θᵢ) aᵢ. Because
  the particles' shares by bin need not equal the frame's areas by bin, Σ qᵢ differs from the frame's total; with
  `--release-normalise on` (default, a labelled device) the rates are scaled by one factor per step so that Σ qᵢ
  equals the frame's total release rate, the factor recorded. Per-patch ownership (each patch's release to the
  particle that owns its area) is a bracket, `--release-map patch`; it is exact per step but spiky (fact 14: a
  2.2 mm particle under the median per-frame maximum, 57 kg/m²/s, empties in about 0.1 s).
- **D7, the release on removed faces: fact 107 (a)'s option (4), written beside the frames.** A write-only export
  addition 07 in the reconstruction's prototype records, per macro step, every face that step's deaths removed with
  the release it had before removal (centroid, `n_derived`, area, θ, mass, film temperature), and the spray step's
  own Σ dm, in one side file per run (`release_removed.npz`), leaving every `field_<k>.vtu`, `surface_<k>.vtp` and the
  history byte-identical. The flight is rerun, compared byte for byte with the existing one, and the side file
  copied in (Task 1). Why (4): option (1) puts the mass on faces of other areas and inclinations, chained through the
  cascade; (2) gives the total but no inclination; the frames alone cannot reconstruct it (fact 12); (4) is a
  superset of (2) at the same cost (one rerun). Fallback if the rerun is not byte-identical, or for flights without
  the side file: option (2) from the history itself, the per-step difference between its sprayed increment and Σ
  `release_rate` · A spread over the table in proportion (`--removed-release history`, labelled), which needs no
  rerun.
- **D8, floor removals at sub-boundaries.** Below the floor (default ¼ of the particle's starting mass, bracket ½)
  the sink stops lowering a particle and books what it still owes; at the next sub-boundary (a full restart, M2
  Task 1's proven route) the particle is deleted and its remaining mass goes to the film account with its own
  enthalpy (reason floor), the owed part with the sprayed enthalpy (reason sink). The delay is at most one
  sub-interval (`--sub-interval`, default 1 % of the window length; Task 5 measures a restart's cost). Live deletion
  with rebuilt packages (fact 10) is the alternative if restarts prove too expensive.
- **D9, the energy of what leaves.** Removed mass carries h(T_film) of its θ bin (the frame's mass-weighted film
  temperature of the releasing patches and removed faces, through §8's table: spec §14.1 condition 4). With imposed
  temperatures (the replay's `ImposedTemperature`, the windows' default) the particle's u stays imposed and only the
  accounts see h_spray; with `--energy free` the particle keeps the rest of its energy, u' = (m u − dm h_spray) /
  (m − dm), and the slab's energy closes (check 15).
- **D10, check 15's configurations.** Stiff: the thermal EOS at real stiffness (K held at its solidus value, 70.3 GPa,
  M1 decision 2). Softened: Spheral's `MurnaghanEquationOfState` with exponent 1, K = 30 MPa and ρ₀ = ρ_free(T) of the
  slab (spec §8.2; its sound speed checked correct in M0), the thermal strength (no shear above 50 % liquid), no
  viscous rows (M4's). Labelled: MC's softened slab is inviscid but for Spheral's artificial viscosity.
- **D11, the replay keeps centre-crossing removal** (spec §14.1, "Open within MC", recommended there): the sink acts
  in windows only. Fact 13 shows that the spray rate and the recession agree within a few percent on this flight, so
  windows and replay remove the same material at the same rate.
- **D12, a physical-mass field.** The sink keeps a registered ScalarField `m_phys` per particle (kg, the ring's in RZ),
  lowered by exactly the booked decrements; Spheral's mass is set from it (times the scale factor, 1 in windows). Record
  checks and accounts read `m_phys`, so they close exactly despite fact 5.
- **D13, run-name brackets** (defaults in brackets, suffixes only when different): `sink` (off; on in windows that
  spray), `floor` (0.25), `sink_h` (held | shrink), `release_map` (theta | patch), `release_normalise` (on),
  `removed_release` (side | history | none), `energy` (imposed | free), `sub_interval` (window length / 100),
  `eos` (thermal | murnaghan:<Pa>), and the slab's geometry name; the window mode's name carries its frame and length.

---

### Task 1: The release on removed faces — export addition 07 and the history's closure

**Files:**
- Modify: `prototype/work-2026-10-08-spheral-mvp/code/reentry_model/coupled.py` and `body.py` (addition 07,
  write-only), its tests; `prototype/work-2026-10-08-spheral-mvp/README.md` (addition 07, the new package SHA-256);
  `prototype/work-2026-10-08-spheral-mvp/rebuild/post/` (07 as a diff, as 01–06)
- Modify: `spheral_frag/contract.py` (an optional run-level item `release_removed`), `analysis/spheral_m1_flight.py`
  (columns `release_removed_kg`, `release_spray_step_kg`), `data/spheral/m1_frames.csv`,
  `tests/test_spheral_frag_flight.py` (the release thresholds)

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

### Task 2: Release tables by inclination (core, drama_env)

**Files:**
- Create: `spheral_frag/release.py` (numpy only at module level), the subcommand `spheral_frag release`,
  `tests/test_spheral_frag_release.py`, `data/spheral/mc_release.csv`
- Modify: `spheral_frag/__main__.py`, `spheral_frag/naming.py` (D13's release brackets)

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

### Task 3: The sink package and the physical-mass field

**Files:**
- Create: `spheral_frag/runner/packages/sink.py`, `tests/test_spheral_frag_runner_sink.py` (marker `spheral`)
- Modify: `spheral_frag/runner/packages/__init__.py` (docstring), `spheral_frag/runner/nodes.py` (`m_phys` and
  `m_start` in the state and the checkpoints, schema version up; old checkpoints read with m_phys = m / s)

- [ ] **Step 1: `SpraySink`** (D2, D5, D6, D8, D9, D12). Built per sub-interval like the other packages, from the
  loads package's surface, exposure, inclinations and shares and the release tables of frames k and k + 1 (or a
  synthetic uniform table for the slab, Task 6). In `finalize` after every step: qᵢ from the table at t (normalised
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
  2.2 mm, 18 processes; RZ 1.1 mm). Commit.

### Task 4: CSF area shares, re-detected within a sub-interval; check 6 with them

**Files:**
- Modify: `spheral_frag/runner/packages/loads.py` (`area_shares="csf"`, `surface_every`), `spheral_frag/naming.py`
  (the `area_shares` bracket gains `csf`), `tests/test_spheral_frag_runner_packages.py`
- Modify: `analysis/spheral_m2_check6.py` (a csf row set), `data/spheral/m2_check6.csv`

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

### Task 5: The window driver's core

**Files:**
- Modify: `spheral_frag/runner/run.py` (`window` mode), `spheral_frag/naming.py` (mode `window`, D13's brackets),
  `spheral_frag/analyse.py` (sink and floor removals, `m_phys` in the checks)
- Create: `tests/test_spheral_frag_runner_window.py` (marker `spheral`)

Depends on M2 Task 11 (the `switch` conversion and the first real-time run); if Task 11 is not committed when this
task starts, the controller asks before building the conversion here.

- [ ] **Step 1: The loop.** `run.py window --body DIR (--from-replay RUN --frame K | --frame K) --length S
  [--sub-interval S] [--sink on|off] [--floor 0.25] [--sink-h held|shrink] [--release-map theta|patch]
  [--release-normalise on|off] [--energy imposed|free] [--area-shares patch|csf] [--surface-every N] [--eos
  thermal|murnaghan:<Pa>] [--loads on|off] [--pressure P] [--out DIR] [--resume]`. Start: a replay checkpoint
  converted by M2 Task 11's switch (masses, densities and free densities divided by s, velocities zero, u unchanged),
  or the body's particles stress-free at rest. Real time from the frame's time, s = 1, temperatures held at frame k
  (or free energy). Each sub-interval: a fresh body from the last checkpoint (the restart), packages (imposed
  temperature unless free, damping off, loads with frame k's table held or a prescribed uniform pressure, the
  deceleration, Hold for held particles, the sink), advance, then at the sub-boundary: floor deletions (D8), the
  sink's bookings as `record.Removal` rows (one per particle per sub-interval, reasons sink and floor), a record
  check from `m_phys`, a row in `subintervals.csv` (sink and floor mass and enthalpy, the normalisation factor,
  KE, max |v|, max surface stress, cost, the restart's own cost), the checkpoint. `summary.json` last.
- [ ] **Step 2: Accounts.** `analyse` reads `m_phys` and the sink and floor rows; it closes the accounts to 1e-12 m₀
  and reports the film account by reason (centre crossing, sink, floor) and its enthalpy.
- [ ] **Step 3: Tests:** a window on the probe slab with a synthetic release stopped after 3 of 6 sub-intervals and
  resumed equals the uninterrupted run bitwise; accounts close to 1e-14 m₀; the run name parses back; a replay's
  name is unchanged by the new brackets. Measure a restart's cost at the slab sizes of Task 7 and at 3D 2.2 mm on
  the sphere (18 processes): this sets D8's default sub-interval. Commit.

### Task 6: Slabs, the synthetic release and the softened configuration

**Files:**
- Modify: `spheral_frag/body.py` (`--geometry slab`), `spheral_frag/runner/eos.py` (`make_eos` for
  `murnaghan:<Pa>`), `spheral_frag/release.py` (a synthetic uniform table), their tests

- [ ] **Step 1: The slab** (check 15). A layer of thickness H on a held cup: 3D a W × W × H block (W = 20 mm,
  H = 10 mm) whose bottom two layers and side walls two layers thick are held (`hold`, as check 5's bar); RZ a disk
  of radius W/2 on the axis with the same base and rim held. Cell-centred lattice, no frame clipping, uniform
  temperature (`--synthetic uniform`, the slab's T a body option); named `_geometry-slab-<W>x<H>mm`.
- [ ] **Step 2: The synthetic release:** a uniform rate R (kg/m²/s) and film temperature on every exposed particle
  facing +x (θ < 45°), as a one-frame table, for the slab only.
- [ ] **Step 3: The softened EOS** (D10): Spheral's Murnaghan, exponent 1, K from the option, ρ₀ = ρ_free(T_slab);
  the thermal strength unchanged; its sound speed checked against √(K/ρ₀) both sides of ρ₀ (M0's check). Tests;
  commit.

### Task 7: Check 15 — the spray sink on a slab

**Files:**
- Create: `analysis/spheral_mc_check15.py`, `data/spheral/mc_check15.csv`, `tests/test_spheral_frag_runner_check15.py`
  (marker `spheral`)

- [ ] **Step 1: The cases.** The core matrix of spec §12 check 15: forms RZ and 3D; dx 2.2 and 1.1 mm; stiff
  (liquid at 920 K, real stiffness, length 2 ms) and softened (D10, length 0.2 s); floors ¼ and ½; `sink_h` held and
  shrink — 32 runs, each with the top pressure on (16.6 kPa, the flight's peak dynamic pressure, spec §2) and the
  shares `csf`. The rate R recedes the surface by 2 spacings of the coarser dx over the run: R = 2 ρ_l · 2.2 mm /
  length (5,280 kg/m²/s stiff, 52.8 soft; a 2.2 mm particle then empties in about 2,400 (stiff) and 5,000 (soft) of its
  own sound-crossing times, so the loss stays quasi-static). Then at dx 2.2 mm, RZ and 3D, held h and the default floor: the shares `patch` against `csf`; the pressure
  off; the flight's own rates (fact 14: 5.7 kg/m²/s soft over 1 s, 355 stiff); a mushy slab (895 K, half liquid,
  spraying at a 920 K film: D9's enthalpy difference) with `--energy free` and `imposed`.
- [ ] **Step 2: Measure** per run (spec §12 check 15): mass closure (remaining + booked − initial) and energy closure
  (free energy: Σ m (u + v²/2) + booked enthalpy − initial − the load's work); the stress in the surface layer and the
  layer below against the applied pressure (false tension or compression, in units of p and of the elastic scale);
  velocities against the noise of the same run without the sink; the surface position where the kernel-estimated
  density (Σ m_j W, diagnostic only) falls to half its interior value, against R t / ρ_l; the neighbour count and h of
  sinking particles; the floor removals' timing; cost per step.
- [ ] **Step 3: Decide and record** (spec §14.1): the floor (¼ or ½), the smoothing-length rule (the smaller surface
  error), the share rule (D4), each with the alternatives' numbers; thresholds from the runs. **If the sink fails** —
  false tension or motion that neither floor, h rule nor share rule removes — stop and report to Asha (spec §14.1, the
  fallback of §19). Commit with the table and plots in `spheral_output/mc/check15/`.

### Task 8: The ID part of check 14

**Files:**
- Create: `tests/test_spheral_frag_runner_check14.py` (marker `spheral`), `data/spheral/mc_check14.csv`
- Modify: `spheral_frag/runner/run.py` only if a defect is found

- [ ] **Step 1: Redistribution** (fact 9) as a test: on 4 and 18 processes, after unbalancing deletions, Spheral's
  Peano–Hilbert redistribution carries `pid`, `x_ref`, `m_phys`, `m_start` and the sink's accumulators bitwise.
- [ ] **Step 2: Through a window:** the slab and a flight window (Task 9's RZ 2.2 mm case) on 1, 4 and 18 processes
  through a switch, sink steps, floor deletions and restarts: the set of IDs at every check is the same for every
  process count, each removed ID appears in exactly one removal row of the right reason, the masses per ID agree to
  round-off across process counts (not bitwise: the reductions differ), and a checkpoint written on N processes is
  read on M (all pairs of 1, 4, 18) bitwise.
- [ ] **Step 3:** In-run redistribution (`redistributeStep`) stays off: the packages refuse a changed local count, and
  every sub-boundary redistributes through the restart. Recorded with the numbers. Commit.

### Task 9: The sink on the flight's frames

**Files:**
- Create: `analysis/spheral_mc_flight.py`, `data/spheral/mc_flight.csv`

- [ ] **Step 1: Replay** the MVP flight's RZ 2.2 mm body from k = 0 to 182 at the defaults (about 20 min at M2's 6.6 s
  per interval), and its 3D 2.2 mm body to k = 100 (about 5 h at 178 s per interval), keeping their checkpoints.
- [ ] **Step 2: Windows** at k = 100, 150 and 182 (RZ 2.2 and 1.1 mm via a 1.1 mm replay to 182, 3D 2.2 mm at 100):
  stiff, 5 ms, the frame's loads held, deceleration on, the sink from Task 2's tables with Task 7's chosen floor, h rule
  and shares; and one softened window at k = 150 (RZ 2.2 mm, 0.2 s).
- [ ] **Step 3: Measure:** the sink's total against the frame's release rate times the window length (normalised: equal
  by construction, the factor recorded; unnormalised: the mismatch of Task 2 Step 3), the film account by reason,
  `analyse`'s closure, the stress change at the surface against the window without the sink, and the switch's
  balance (M2 Task 11's measure) with the sink on. Commit.

### Task 10: Documentation and the spec amendment

**Files:**
- Modify: `README.md` (the large-fragment section: the window command, the sink, the release tables, MC's
  verification table), `docs/spheral_frag_assumptions.md` (the sink, its devices D2–D9 each marked verified,
  analytic, measured or assumed; the release closure), the spec (§§2, 10, 11.3, 12 checks 14–15, 14, 14.1, 15, 17, 18,
  19, 20: MC as built and measured, D1's amendment, D5's and D7's decisions), this plan's "Decisions at review".
- [ ] One commit, `docs: Spheral MC as built and measured: ...`.

## Order and parallelism

1. **Tasks 1 and 2** (drama_env and the prototype) run in parallel with **Tasks 3 and 4** (the runner); Task 2's
   Step 3 needs Task 4's csf shares as numpy functions, its other steps do not.
2. **Task 5** after Task 3 and M2 Task 11. **Task 6** after Task 5 (its Step 3 can start at once).
3. **Task 7** after Tasks 4 and 6. It decides whether the sink works: stop and report if it fails.
4. **Task 8** after Task 5 (its Step 1 at once). **Task 9** after Tasks 2 and 7.
5. **Task 10** last.

Wall-time estimate at the measured step (spec §2; the slab's cost per particle as the sphere's): the stiff slab runs
2 ms at 0.05–0.1 µs per step (20,000–40,000 steps; 3D 1.1 mm about 3,000 particles, minutes to half an hour each),
the softened ones 0.2 s at about 2.5–5 µs (40,000–80,000 steps); check 15's 32 core runs and about 12 variations
take a day or two; Task 9's replays about 6 h, its windows (stiff 5 ms: 50,000–100,000 steps) hours; Task 1's flight
rerun about an hour and its re-prepare-free release tables minutes.

## Questions for Asha

1. **The release on removed faces** (D7, fact 107 (a), M1 decision 16). Recommended: option (4), a write-only side
   file (`release_removed.npz`) from export addition 07, with the reconstruction's flight rerun and compared byte for
   byte, so the frames and every M2 run keep their provenance; the history-derived option (2) as the fallback and a
   bracket. The frames alone cannot reconstruct the removed release (fact 12). Or: option (2) only, no rerun.
2. **Which particles sink** (D5). Recommended: every exposed surface particle in an inclination bin with release,
   dropping §10's and §14.1's "above the liquidus" — at a 2.2 mm particle's centre only 0.3–41 % of the sprayed mass
   comes from material above the liquidus (fact 15). Amends §10 and §14.1.
3. **The mapping and its normalisation** (D6). Recommended: by inclination as the spec says, normalised per step to the
   frame's total release (a labelled device, factor recorded; Task 2 measures how far it is from 1); per-patch
   ownership as a bracket.
4. **The area shares under the sink** (D4). Recommended: build the `csf` shares and choose between them and the patch
   shares in windows by check 6 (Task 4) and check 15 (Task 7); the replay keeps the patch shares (M2 decision 9).
5. **The replay keeps centre-crossing removal** (D11, spec §14.1's open item). Recommended; the spray and the recession
   agree within a few percent on this flight (fact 13).
6. **Check 15's extensions** (Task 7). Recommended: the top pressure on a held cup, the shares compared, the flight's
   own rates, and a mushy slab with free and imposed energy, beside the spec's core matrix (32 runs); and the softened
   slab without viscous rows (D10), since those are M4's.
7. **Floor removals at sub-boundaries** (D8): the sink stops at the floor and the particle is deleted at the next
   restart, at most one sub-interval late. Recommended.

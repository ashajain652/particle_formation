# Level-set architecture: to-do list

Started 2026-11-14 at Asha's request. This file collects the work still to do on the level-set architecture, in the
order it is expected to be done. It points to the design
(`docs/superpowers/specs/levelset/2026-10-11-level-set-deforming-body-design.md`, called "the design" below) and to the
stage plans in this folder, and it does not repeat them. The design governs where they disagree. Mark an item done here
with its date and the commit that finished it.

## 1. Finish stage 3 as amended by design section 21

Build the one-skin-depth geometry with plan `03b-stage3-one-skin-depth.md`. In it, every skin is 1.2 mm deep. The gas
surface φ is the evolved level set, and the skin base β is derived 1.2 mm inside it. The columns come from one division
of the layer, and the plan ends with the merge of main, the composite rheology and the kinked-profile spraying.

Status: started 2026-11-14 on branch `levelset-stage3`, run subagent-driven.

## 2. Settle the calls in design section 21 that are marked for Asha to check

Section 21 marks four of Claude's calls for Asha to check:
- the default skin depth of 1.2 mm (§21.3);
- superseding her hand-down rule of 2026-11-10 (§21.3);
- the volume keeping the solid's density, so that a molten cell returned to it shrinks by about 15 % (§21.6);
- the amended spraying line, where a skin whose flowing layer reaches its base uses the base as a floor (§21.7).

None of them blocks the build. The plates of plan 03b's Task 7 test the first call: if the melting front crossing the
skin base makes them fail, a larger depth is Asha's decision then (§21.12).

## 3. The volume melt flow (design §21.9, with §8 and §20.7)

Solve the flow of all molten metal inside β with incompressible Navier–Stokes, as §8 sets out. The rigid wall is ψ's
zero contour, where the metal is half liquid. The free surface is β, where the flow couples to the skins. It gets its
own design after item 1 is built. Asha brought it forward from stage 4.

Cost estimate: of order 30 to 170 hours a flight on the laptop while thick melt persists (§20.13). This is *estimated*,
very roughly, from about 10⁴ flow steps of 10 ms on about a million cells, and it is to be measured. The parallel plan
for it is §20.15's third layer: PETSc's solvers with several MPI ranks on the Mac Studio.

Item 5 below asks whether this solve should reach further than the volume.

## 4. Finer markers near edges (design §21.10), and thin metal (design §21.11)

- **Finer markers near edges.** Space the markers by φ's curvature, down to about 0.4 mm near edges and 2 mm elsewhere.
  Two prerequisites come first:
  - runoff's cost must grow no faster than the number of markers (1293 s of a 1528 s step at 0.4 mm, *measured*
    once, 2026-11-13);
  - the part limit must be untied from the marker spacing.
- **Thin metal.** Design the treatment of metal thinner than twice the skin depth. Until then the run stops there
  (§21.11).

## 5. Assess, and if justified add, Navier–Stokes solves of the melt flow (Asha, 2026-11-14)

**The request.** Asha asked whether the model solves the Navier–Stokes equations for the melt flow. She wants to know
whether the level-set architecture should add real Navier–Stokes solves of the melt. Two forms are to be weighed:
- full Navier–Stokes;
- a thin-film form that keeps the terms the skyscrapers' skins drop.

The item must record what each would cost and how it relates to the future coupled-CFD gas loads.

### What the skins solve today, in both architectures

Checked in the code on 2026-11-14. Each skin, a column of cells through the flowing layer, solves a one-dimensional,
implicit, time-dependent momentum equation through its depth: `momentum_step` in `levelset/skin.py` (line 980) and in
`skyscrapers/skin.py`. In that equation:
- each cell's inertia balances the viscous drag from the cells above and below it;
- a push per unit volume drives the flow. The push has four terms:
  - the gas pressure gradient;
  - the deceleration;
  - the melt's own hydrostatic pressure, ρ_l (g · n) ∇_s H, where ρ_l is the liquid's density, g · n the
    deceleration's component out of the surface, and ∇_s H the slope of the layer's thickness along the surface;
  - surface tension's capillary pressure, −σ ∇_s κ, where σ is the surface tension and κ the surface's curvature;
- the gas shear acts on the top cell.

The two thickness terms are solved implicitly in the column heights:
- the skyscrapers': `push` and `thickness_step` in `skyscrapers/skin.py`, lines 1182 and 1647;
- the level set's: `_push` and `thickness_step` in `levelset/runoff.py`, lines 239 and 282.

Runoff then moves metal between neighbouring columns, upwind and implicit. Each cell's momentum travels with its metal
(`levelset/runoff.py` lines 525 and 637–663, and `skyscrapers/skin.py`'s `transport`). That is a first-order form of
the along-surface inertia.

This corrects the answer Asha was given on 2026-11-14, which said that surface tension and the melt's own hydrostatic
pressure are absent and that along-surface inertia is dropped. Relative to full Navier–Stokes, the skins really drop
three things:
- the viscous stress between neighbouring columns (sideways shear);
- any pressure set by the melt's own motion through incompressibility (the pressure is taken as hydrostatic plus
  capillary across the layer);
- solving velocity and mass transport together (they are split into a momentum step, a thickness step and a transport
  step).

The element model is a steady closed-form lubrication solution (`shared/lubrication.py`). Girin's spraying criterion
alone comes from linearised Navier–Stokes.

### Where the level set already plans Navier–Stokes

The design already puts full incompressible Navier–Stokes in the metal volume (§8.1, §20.7, §21.9; item 3 above). So
the open question is not whether the level set gets a Navier–Stokes solve. It is whether that solve, or a fuller
thin-film form, should also govern the 1.2 mm layer that the skins hold.

The thin-film form holds while two conditions are met:
- the layer is much thinner than the radius of curvature of the surface beneath it;
- the layer's top slopes gently (§17 item 8).

On the skyscrapers' 100 mm flight to 30 s, every flowing layer was 0.84 to 0.93 mm deep (*measured* once,
2026-11-12, §21.3), well inside those conditions. The 50 mm comparison flights broke them: by 187.0 s the rear pool
stood 51.5 to 66.5 mm deep on a base about 20 mm from the centre (§17 item 8). In the level set, that pool belongs to
the volume flow (item 3), not to the skins.

### The two forms

**(a) Full Navier–Stokes in the layer as well.** The skins would be replaced by the volume flow everywhere.
- **The obstacle is resolution.** Girin's spraying reads the flowing layer's velocity profile at 10 µm cells.
- **Isotropic cells.** Isotropic 10 µm cells through a 1.2 mm layer on the 100 mm sphere would number about 4 × 10¹⁰:
  3.1 × 10⁻² m² of surface times 1.2 mm of depth, divided by (10 µm)³. That is out of reach (*estimated*).
- **Anisotropic cells.** Thin prism layers in the volume mesh were rejected in §20.6. The skyscrapers' prism layers
  folded at 90-degree edges, and the band would be remeshed every step.

So full Navier–Stokes in the layer is not a practical option at the resolution spraying needs. Full Navier–Stokes
stays where the design puts it: in the volume, at the volume's 0.7 mm fine size, coupled to the skins at β.

**(b) A fuller thin-film form in the skins.** Keep the skins and add what they drop.
- **Sideways viscous stress** between columns. Plan 03b's division gives the shared boundary between every pair of
  neighbouring columns at every depth. Plan 03b already uses those faces for sideways heat conduction, so viscous stress
  can use them too. Each pair of cells at the same depth exchanges momentum at a rate set by the viscosity, the shared
  face's area and the distance between the columns. This has the same form and the same implicit solve as the
  sideways conduction of plan 03b's Task 4.
- **A non-hydrostatic pressure correction.** This would come from a depth-integrated or second-order thin-film model
  (for example of the integral boundary-layer kind), where the layer's thickness changes over lengths not much longer
  than its depth.
- **Coupled transport.** Velocity and mass transport would be coupled in one sub-step instead of split. This is the
  least valuable of the three while the sub-steps are short.
- **Cost** (*estimated*): the sideways stress makes each momentum sub-step one coupled solve over all skins, as
  sideways heat makes each heat sub-step (§21.12). That cost will be known once plan 03b's Task 4 measures the sideways
  heat solve. It should be far cheaper than item 3's volume flow.

### How to assess it

Before building either form, measure how large the dropped terms are where they could matter.

1. **A diagnostic.** Record a cheap diagnostic on existing and new flights, computed from the solved velocities without
   feeding back into the flow. For each column, it gives the sideways viscous stress and the non-hydrostatic pressure
   gradient, each as a share of the push and the drag that are kept. Report it where spraying happens (the slip band,
   the shoulder, rims and the ring).
2. **A benchmark.** Run a case where thin-film theory is known to fail, such as a layer flowing over a sharp edge or a
   step, or a pool spreading, through the skins. Compare it against the volume Navier–Stokes solve once item 3 exists.
3. **The decision.** Add form (b) if the dropped terms change the push or the drag by more than a tolerance to be
   agreed with Asha (a few per cent) where spraying happens. If they do not, record the diagnostic as the evidence that
   the thin-film form suffices.

### Relation to the future coupled-CFD gas loads (Asha, 2026-11-12)

- **The same seam.** Whatever form the melt solve takes, it receives the gas pressure, shear and heat flux on its
  surface as loads. Today these come from the analytical laws in `shared/surface_flow.py` and `shared/heating.py`,
  through the level set's `loads.py`. A coupled CFD solution would replace those loads at that seam. So the two pieces
  of work meet only at the gas-loads interface. That interface should stay one replaceable interface in
  `reentry_model/shared/`, as the note of 2026-11-12 asks, so either can be built first.
- **Separate solvers.** The gas and the melt cannot share one Navier–Stokes solve. The design's melt solver is
  incompressible, and the hypersonic gas is outside its reach (§8.1). Chen et al. (2024) solve gas and metal together
  only because their gas is slow. A coupled model would therefore alternate:
  - the CFD solves the gas around the current shape;
  - the melt and thermal solve advance under those loads;
  - the new shape goes back to the CFD.

  A hypersonic CFD solve per shape is likely to take hours, so the exchange would happen every several macro steps,
  or when the shape has changed by more than a set amount, not every step (*estimated*). The level set's φ is a
  watertight surface, which helps: a CFD mesh can be generated from it directly.
- **The CFD would make a fuller melt solve worth more.** The analytical loads are smooth fits, and they are least
  certain behind the shoulder and at edges. A CFD solution would give sharp local pressure and shear changes at edges
  and in separated flow, which is where the thin-film assumptions are weakest. While the loads remain analytical, the
  accuracy of a fuller melt solve there is capped by the loads themselves. So item 5's assessment should be repeated
  once CFD loads exist.

Status: not started. It waits for item 1, since its diagnostic needs the division's shared faces. The benchmark waits
for item 3.

## 6. Parity flights and stage 4's measurements (design §20.11 items 4 and 5)

- **The flow solve's fine-size convergence:** a pool at fine sizes of 0.7, 0.5 and 0.35 mm.
- **Parity flights:** the 50 mm and 100 mm flights against the skyscrapers', two seeds each, judged within the scatter
  between seeds up to the time the skyscrapers' column assumptions first break.
- **Measured cost:** the flow solve's real cost, which replaces item 3's estimate.

## 7. The run driver

`reentry_model run --architecture levelset` still exits with code 2, because the level set has no run driver
(`Reentry Model/CLAUDE.md`). Build it once item 1 lands, so that level-set flights run from the same command as the
skyscrapers'.

## 8. Step 2: the leeward mass-loss mechanisms (design §20.9) and the tunnel configuration (§11)

- **The leeward mechanisms:**
  - collection into a ring or rear cap;
  - resolved pinch-off into fragment records;
  - sub-grid rim shedding;
  - the mean wake load on protrusions.
- **Deferred:** tearing, which needs a viscoplastic strain solve and a critical strain, and the wake's fluctuating load.
- **The tunnel configuration** waits for the tunnel record (§17 item 2).

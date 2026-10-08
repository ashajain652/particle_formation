# `spheral_frag` core (M1) — assumptions and methods

State of the Spheral-free core on `main` as of 2026-10-08: spec
`superpowers/specs/2026-10-02-spheral-large-fragments-design.md`, plan `superpowers/plans/2026-10-07-spheral-m1.md`
(its "Decisions at review" 1–14 are cited as *decision n*). Referenced from `model_assumptions.md` (decision 8). Each
item is marked **(verified)** where a test compares it with the finite-element model's own numbers or an exact answer,
**(analytic)** where the reference is a closed form, **(measured)** where it was measured on the real flight without a
reference to agree with, and **(assumed)** where it is a choice or a placeholder. The numbers are in the README's
verification table (section "Large fragments — `spheral_frag` core (M1)") and in `tests/test_spheral_frag_*.py`.

**The input it was measured on.** One real flight: the 100 mm `AA7075_scheil` physics flight of 2026-10-08
(`reentry_model_output/model_d100.00mm_v07.50000kms_h077.500km_us76_sesam-table_none_fem-physics_melt-girin_material-AA7075_scheil_deeprunoff-off/`),
1,255 frames at the 0.5 s step, seed 12345, molten cascade on, deep runoff off. It was written by the reconstructed
Step 3 prototype `prototype/work-2026-10-08-spheral-mvp/code/` (package SHA-256 9629f39c…), not by `reentry_model`,
which does not melt yet: an MVP input to build the Spheral pipeline on (decision 13), whose physics is not the Step 3
line of record (no runoff-flux fix, no rigid-substrate rule) and whose 0.5 s step leaves the molten backlog of spec §2
in the frames.

## 1. The frame contract and the import

- **One table names every finite-element field** (`spheral_frag/contract.py`, `FE_FIELDS`, 47 items by node,
  tetrahedron, surface patch, frame, history column and run key); no other module spells a finite-element name. Every
  item is present under its name on the flight and marked `confirmed` (decision 14). **(verified)**
- **Units read from the export, not from names.** `release_rate` is kg/m² *per macro step*; `deep_thickness` is a mass
  per area written as a thickness, m_d / (ρ_l A); `p_w` = 0 means "not evaluated", which on the flight is frame 0 alone
  (`flow_eval` 0). The faces a step's element deaths exposed carry the record's flow evaluation (`flow_eval` 2, up to
  1,517 faces). **(verified)**
- **The export's winding is not trusted.** Up to 1,662 faces of a frame are wound inward in `surface_<k>.vtp`; every one
  has exactly one active owner tetrahedron, and `prepare` rewinds each face away from its owner. After that the surface
  is closed as an oriented surface on every frame (0 open directed edges) and its divergence volume equals the
  tetrahedra's to 2.2e-16. **(verified)** It is not manifold: up to 288 non-manifold edges (on 1,204 frames) and 15
  vertices (on 131), where the staircase surface pinches after element deaths. What Spheral's polyhedron builder needs
  is M2's question. **(measured)**
- **Mass** is Σ φρV over the active tetrahedra plus the film and deep accounts, and matches the history's `mass_kg` to
  4.2e-9 relative, within the CSV's 9 significant digits (threshold 5e-9). The frame's shape does not shrink with φ, so
  φ must be read for the mass to be right. **(verified)**
- **The flight direction** is the run JSON's `settings.v_hat_body`, (1, 0, 0). **(verified)**
- **The lee base pressure's reference** is the history's `p_w_stag_step_Pa`, every step's own stagnation wall pressure;
  `p_w_stag_Pa` is NaN until the first film (0–23.5 s) while those frames already carry loads (decision 14). The
  largest `p_w` on a frame equals it to a median 6e-10, except on the steps where the stagnation patch died (worst
  13 %, k 64). **(verified)**
- **The release on removed faces is lost.** Σ release_rate · A over a frame's surface is 12 % below the history's
  sprayed mass over the flight (1.025 of 1.169 kg; −30 % to +4 % per step): what was released from the faces that
  step's deaths removed is not carried onto the frame. A sink that reads `release_rate` (spec §10, milestone MC)
  inherits that deficit until the export carries it. **(measured)**

## 2. The prepared frame

- **Stored reduced on mesh versions** (decision 12): one mesh file holds the finite-element nodes and tetrahedra; each
  frame stores a bit mask of its active tetrahedra, φ where it is not 1, the nodes whose position moved, its faces and
  fields, deflated. Nodes, tetrahedra and patch geometry are rebuilt on loading, bitwise against the finite-element
  files. The flight is 791 MB in total, with a median frame of 0.51 MB and one mesh version; preparing it took
  1,192 s. **(verified; measured)**
- **Not copied:** the informational fields `q_rad`, `T_patch`, `r_droplet`, `we_s`, `kn_local`, `flow_eval`, which the
  finite-element run keeps. **(assumed)**

## 3. The material table

- **The finite element's own arithmetic, not a fit.** The table stores the finite-element material's enthalpy nodes
  (piecewise quadratic h(T) with its quadratic-formula inverse) and liquid-fraction nodes, and repeats its operations
  one for one with numpy alone. At 100,543 temperatures (uniform 250–1,500 K plus every node ± 1e-9 K), h, f_l and c_p
  are bitwise and the inverse is within 1.6e-12 K. On every frame of the flight, nodal `liquid_fraction` equals the
  table's f_l(T) bitwise. **(verified)**
- **Zero of energy at 300 K** (spec §8.1), against the finite element's 293 K: h = h_FE(T) − h_FE(300 K), with the
  offset 6,152.3 J/kg subtracted after the finite element's own value, so differences are exact. **(verified)**
- **np.interp's last operation** is fused (FMA) by drama_env's numpy and not by Spheral's numpy 1.26: 6 f_l and 38 c_p
  values of 10,000 differ by one ulp. The table records which numpy wrote it and reproduces that rounding, so Spheral's
  Python matches drama_env bitwise. **(verified)**
- **Viscosity:**
  - Slurry follows Li et al. (2014), η = [0.871 − 0.00849 γ̇^0.74924] exp(3.7311 f_s) Pa·s for 0.1 ≤ f_s ≤ 0.5, with
    the shear rate held at no more than 367 s⁻¹.
  - A log-linear bridge in T runs from f_l = 0.9 (906.37 K) to the liquid's 1.3 mPa·s at 908 K.
  - Below 50 % liquid the viscosity is infinite.

  The solid and mush rows and the tearing law are M2–M3's. **(assumed; values from the Step 4 plan §4)**
- **Mechanical columns are labelled placeholders** (decision 2), listed as `provisional` in every table:
  - E is 71.7 GPa up to the 748 K solidus, then falls linearly to zero at 50 % liquid (895.10 K).
  - ν = 0.33.
  - K is held at its solidus value, 70.29 GPa, through mush, slurry and liquid.
  - α = 23.4 × 10⁻⁶ K⁻¹.
  - ρ_free mixes the expanded solid with the finite element's 2,400 kg/m³ liquid (decision 3). That implies 12.3 % on
    melting; the spec's 6.5 % is recorded but unused.

  **(assumed)**

## 4. Thickness and layer depths

- **Thickness** is the distance along the inward normal from each patch centroid to the first other surface
  triangle. The ray test is Möller–Trumbore without back-face culling, so it does not depend on the winding of the
  triangle it hits. On frame 0 every patch is within 5.8e-4 of 2R (faceting). No ray escapes on any frame. The
  thinnest patch over the flight is 0.21 mm (k 166). **(verified on a sphere, a dumbbell neck and a slab; measured on
  the flight)**
- **Depths are marched along the derived surface's normals** `n_derived` (decision 9; Taubin-smoothed, unit, and within
  90° of the facet normal on every patch).
  - The marching samples the P1 f_l every 0.05 mm and places the crossing by linear interpolation. Slurry is
    contiguous f_l > 0.5; liquid is f_l ≥ 1 − 1e-9.
  - On the slab the depths are exact to 5e-18 m. On a coarse sphere they are within one element of the continuum.
    **(analytic)**
  - The facet-normal depths stay a diagnostic: on frame 100, 43 % of the rays along facet normals left the staircase
    surface before f_l fell to 0.5.
- **On the flight**, slurry reaches the body's far side on 250 frames (25.5–150 s; deepest 81 mm), and liquid is at
  most 2.25 mm deep. A layer reaching the far side is read as the whole remaining thickness. **(measured)**

## 5. The three zones

- **The rule** (spec §10), per patch, with layer = film + contiguous slurry depth:
  - zone 1 (sprayable skin): layer ≤ δ_m;
  - zone 2 (thin runoff): δ_m < layer ≤ the film limit, or 0 < layer ≤ the film limit where δ_m is NaN;
  - zone 3 (bulk, Spheral's): layer > the film limit.

  The film limit is 2 mm, with a 3 mm bracket. **(analytic on the slab)**
- **Both readings of a layer thicker than the film limit** are reported (decision 6): all of it is Spheral's, or only the
  part below the limit; the choice is made before M4. `deep_thickness` stays out of the layer; `include_deep` adds it
  and changes nothing else. **(assumed)**
- **On the flight**, zone 3 exists on 336 frames (24.5–192 s), at most 152 cm² at the 2 mm limit and 144 cm² at 3 mm.
  Deep runoff is off, so `deep_thickness` changes nothing. These are 0.5 s-step frames, so part of that layer is the
  molten backlog (spec §2). **(measured)**

## 6. Loads

- **The windward table** bins the frame's own `p_w` and `tau` area-weighted into 1° bins of θ, the angle between the
  outward facet normal and v_hat. A patch is valid where p_w > 0. Binning changes the drag by at most 7.9e-4 on the
  flight. **(verified)**
- **The lee extension** (decision 5; labelled "Step 4 §5's lee model, to be replaced by its module"):
  - pressure: beyond the last bin with loads, 3 % of `p_w_stag_step_Pa`, held to 180°;
  - shear: half the last bin's value to the separation angle (180°), zero beyond.

  The brackets are 1 and 5 %, 150°, and 0 or 1 times the last shear. On the flight the frames carry loads to 179.5°,
  so the extension never acts. **(assumed)**
- **The drag compared with the history** reads the table at each facet's derived-surface inclination θ(n_derived)
  (decision 11). Pressure acts on A cos θ; shear acts along the derived tangent on A (n_facet · n_derived). Summing the
  staircase facets' own loads instead loses 26–37 %.
- **The history's drag is not met within spec §7.3's 5 %** (Review focus 6: reported, not fitted). D_smooth / D_hist − 1
  is:

  | Phase | Times | D_smooth / D_hist − 1 |
  |---|---|---|
  | Before the first film (no flow branch recorded) | 0.5–23.5 s | +10.3 to +10.6 % |
  | Merged branch | 24–49 s | +4.7 to +10.4 % |
  | Shock layer | 49.5–192 s | −35.9 to +9.1 %, worst just before 192 s |
  | Free-molecular branch, which the finite element also uses at Mach ≤ 1 | 192.5–626.5 s | constant −11.0 % |

  The history's drag comes from the SESAM-table drag coefficient, and the frame's loads come from Step 3's surface-flow
  model; the +10 % already present on the intact sphere, before anything has melted, points to that difference rather than
  to the core. **(measured)**

## 7. The fragment record and the debris log (definitions chosen where the spec leaves them open)

- **Dust** is a particle with largest principal damage ≥ 0.99. It attaches to the group of its nearest intact particle
  within 1.5 mean smoothing lengths; otherwise it goes to the debris log.
- **Groups:** the most massive group is the main body. A group of fewer than 30 intact particles is debris.
- **Clear** (spec §9.3) means farther than 3 mean smoothing lengths and moving away.
- **Release** is the first check at which the group was separate, decided once it is clear. Groups are followed across
  checks by permanent particle IDs. **(assumed; tested)**
- **3D shape:**
  - Principal lengths are 2 √(5 λᵢ) of the mass-weighted covariance (exact for a uniform ellipsoid; 0.09 mm on a
    10/6/4 mm lattice ellipsoid at 0.5 mm spacing).
  - The equivalent diameter is (6V/π)^⅓, with V = Σ m/ρ.
  - The area is the ellipsoid's, by Thomsen's approximation (at most 1.061 %).

  **(analytic)**
- **Axisymmetric shape:**
  - A group is a ring when every member's r exceeds its smoothing length; otherwise it is a cap.
  - A ring's lengths are its cross-section's (2 √(4 λ), exact for a uniform ellipse) plus 2πR. Its area is by Pappus
    with Ramanujan's perimeter.
  - A cap's lengths come from its moments as a body of revolution.

  **(analytic)**
- **Phase state:** fluid if every particle is more than half liquid; mixed if some are; solid or mush otherwise.
- **Weber, Ohnesorge and the breakup flag** apply to fluid fragments only:
  - σ = 0.80 N/m;
  - the viscosity is Li et al.'s at the mass-weighted temperature in the zero-shear limit (decision 7);
  - breakup when We > 12 (1 + 1.077 Oh^1.6) (Pilch and Erdman 1987).

  **(assumed)**
- **Origin** (route and mechanism) comes from the run's mode until M3's runner records each group's own. **(assumed)**
- **The debris log's size** is an upper bound: the group's largest extent plus one spacing. **(assumed)**

## 8. Mass accounts

- The starting mass equals main body + resolved fragments (with attached dust) + unresolved debris + film account, at
  every check, summed exactly (`math.fsum`). Each removal books its mass and its own enthalpy.
- `analyse` exits 1 beyond 1e-12 m₀. Measured closures:
  - synthetic history: 2e-18 m₀;
  - synthetic runs on the real flight's frames 100–110 and 800–810 (3D, dx 2.2 mm): 7e-18 and 4e-17 m₀.

  **(verified)**
- Those runs are `fake_run`'s stand-in for the runner (M3): a lattice with centre-crossing removals and a carved ring.
  They test the bookkeeping, not fragmentation.

## 9. What M1 does not settle

- The drag gap of §6, and the lost release of §1, are properties of the input frames.
- The non-manifold surface is M2's problem.
- The molten backlog is in the 0.5 s frames; M3 and M4 need frames free of it (spec §13.1).
- The mechanical columns are placeholders.
- The frames are to be rewritten once Step 3 adopts the export additions 02, 03, 04 and 06 (decision 13);
  `prepare`'s contract checks then confirm the new frames.

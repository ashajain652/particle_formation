# Viscous Rayleigh–Taylor instability on liquid over slurry (parked 2026-10-29)

This note records a goal Asha set on 2026-10-29: to treat the Rayleigh–Taylor (RT) instability of the melt layer with
its viscosity, rather than with the inviscid rules the model uses now. Nothing here is implemented. The literature
below is recalled from memory, not read, and must be checked against the sources before it is cited. The one number
worked out here is a back-of-envelope estimate, made once.

## How this came up

1. **The model's RT rules are inviscid.** `reentry_model/spray.py` (`rayleigh_taylor`, `rayleigh_taylor_bounded`)
   decides whether the mode grows from the deceleration, the density, surface tension and the depth of the pool. The
   depth enters only through the rule that a wave fits if its wavelength is shorter than about 2π times that depth.
   Viscosity does not appear.

2. **On the layered skins the depth is ambiguous where liquid lies on slurry.** Slurry is the partly solid metal
   between the liquidus and the rigid point. Sub-plan 22 Task 1 (2026-10-28) found that at the nose, under the ablating
   flight's peak deceleration of 103 m/s² (10.5 g), with 0.3 mm of liquid over 3.1 mm of slurry, the choice decides the
   outcome. Read over the whole flowing depth (liquid plus slurry), the RT mode takes 79 columns 1–10 degrees from the
   stagnation point, at a wavelength of 19.7 mm, with droplets of 3.05 mm. Read over the liquid alone, it takes none.
   Details are in `prototype/work-2026-10-28-skin-spraying/measure/task1/report.md` (section 3.4, decision 1) in the
   `subplan20` worktree.

3. **Neither reading includes viscosity.** Counting slurry as depth treats it as if it flowed as freely as the liquid.
   Excluding it assumes it takes no part at all. The slurry's viscosity is 0.87–5.6 Pa·s (Li et al.'s law), about a
   thousand times the liquid's (about 0.001 Pa·s). Over the inviscid mode's growth time of about 7 ms, viscous effects
   reach 2.1–5.4 mm into it, which is the whole layer. So the slurry layer is viscous through its depth, and the
   inviscid rules do not describe it.

## A first estimate: the slurry layer may still be unstable

It is tempting to conclude that viscosity simply stops the mode. A rough estimate says otherwise, so the question is
open. For a thin viscous layer hanging under a wall, the lubrication (thin-film) limit gives a growth rate

s(k) = (h³ / 3μ) · (ρ g k² − γ k⁴),

where:
- h is the layer's depth;
- μ is its dynamic viscosity;
- ρ is its density;
- g is the deceleration pressing the layer against the body;
- γ is the surface tension;
- k is the wavenumber.

The fastest wave has k² = ρg / (2γ), and its growth rate is s_max = h³ (ρg)² / (12 μ γ).

The estimate takes h = 3.1 mm, ρ = 2 500 kg/m³, g = 103 m/s² and γ = 0.9 N/m. The fastest wavelength is then about
17 mm, close to the 19.7 mm the inviscid rule picked. Its growth time, the inverse of s_max, is:
- about 5 ms at the low end of the slurry viscosity (0.87 Pa·s);
- about 31 ms at the high end (5.6 Pa·s).

So a viscous layer this deep may grow nearly as fast as the inviscid mode, or several times slower. Either way it is
not obviously stable.

**The estimate's caveats:**
- The lubrication limit wants the wavelength much longer than the depth; here kh is about 1.2, so it is marginal.
- It ignores inertia, the gas shear and the flow carrying the layer away from the nose.
- It uses the slurry's viscosity at rest. The model's slurry (`rheology.py`, Li et al. 2014) is a viscous fluid with
  no strength, but its viscosity falls as it is sheared, which a growing wave would do. Like the model, it has no yield
  stress: the model's slurry flows at any stress up to half solid and is rigid from there. Real semi-solid aluminium
  near that point may hold a yield stress. That claim is recalled, not checked, and it is not part of the model.
- It treats the liquid and slurry as one layer, not two layers of very different viscosity.

## What to work out

1. **Whether it matters.** From sub-plan 22's flights, record the mass of liquid and slurry on columns where the two
   depth readings disagree, the time those columns spend under strong deceleration, and how fast the flow carries
   their metal away. Compare that residence time with the viscous growth time. If the metal leaves long before a wave
   could grow, record it as a stated limitation and stop.
2. **The sources (to find and check).**
   - Chandrasekhar's *Hydrodynamic and Hydromagnetic Stability* (1961), the chapter on RT instability with viscosity:
     the general viscous dispersion relation.
   - Bellman and Pennington (1954, Q. Appl. Math.): viscosity and surface tension together.
   - Fermigier, Limat, Wesfreid, Boudinet and Quilliet (1992, J. Fluid Mech. 236): RT instability of a thin viscous
     layer under a plate, the lubrication limit used above, and the pattern it forms.
   - Work on two-layer RT, with a less viscous layer over a more viscous one: our case of liquid over slurry.
   - Only if the semi-solid literature shows that slurry near half solid holds a yield stress, which the model does
     not include: work on RT in a yield-stress fluid.
3. **The model.** Pick the simplest rule that brackets the physics:
   - (a) a viscous growth rate for the RT branch from the full dispersion relation, on the layer's own viscosity
     profile;
   - (b) a damping gate like the thin branch's viscous film (`spray.viscous_film`, Γτ* ≥ 1 blocks), applied to RT;
   - (c) a two-layer relation for liquid over slurry.

   Each needs a droplet-release rule after the instability grows. That rule is not obvious for a viscous layer, which
   may form long drips rather than droplets.
4. **Validation.** Check the chosen rule against a published experiment or exact solution before using it in a flight.
5. **Asha decides,** with the mass at stake in hand. Record the decision under the Step 3 conventions: a dated
   amendment and facts in `00-shared-context.md`, and the code in a work copy with its tests.

## Related

- `2026-10-22-viscous-film-instabilities.md` in this folder: the same question for the thin shear branch.
- Sub-plan 22 (`docs/superpowers/plans/melt-spraying-subplans/22-skin-spraying-flight.md`, branch
  `step3-rigid-substrate-scheil`): Task 1c's decision 1 on which depth the RT rules read is the interim choice this
  milestone would replace.

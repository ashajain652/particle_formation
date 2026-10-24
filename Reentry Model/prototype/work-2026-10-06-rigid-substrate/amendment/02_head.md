## Amendment of 2026-10-06 — `Material.T_rigid`, the temperature at which the material is half liquid

> Part of the rigid-substrate amendment (sub-plan 09's amendment of this date holds the design; facts 78–87 in
> `00-shared-context.md` the measurements). The code below is the tested code, as a diff against this sub-plan's
> `material.py` of 2026-09-28 (the copy the measurements used carried it for the Scheil runs).

`RIGID_LIQUID_FRACTION = 0.5` and the property `Material.T_rigid`: the temperature at which the material's own liquid
fraction is one half, found by bisection (the liquid fraction is monotonic), infinite for a material that does not melt.
Below it mush is coherent and carries load, a rigid substrate for the melt film; above it mush is a slurry that flows —
Chen et al. 2016's semi-solid law ends at 50 % liquid and Li et al. 2014's slurry data begin there, the boundary Step 4
and the large-fragment design use. Values (fact 78): 829.0 K for `AA7075_range`, 895.10 K for `AA7075_scheil`, 850.0 K
for `AA7075`. Nothing else in the material changes; the existing tests are unchanged and the new one is appended after
the empirical-data test.


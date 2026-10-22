# Sphere drag across rarefaction and compressibility: sources

Collected 2026-10-22 while assessing a drag law for molten AA7075 droplets in the transitional regime (Knudsen
number 0.1-20 with speed ratio above 0.05), the corner the droplet drag law in
`prototype/work-2026-10-14-particle-settling/particle_trajectory.py` flags as uncalibrated. Bibliographic details
were confirmed against Crossref (api.crossref.org) unless marked otherwise.

## Primary sources

1. **Loth, E. (2008).** Compressibility and Rarefaction Effects on Drag of a Spherical Particle. *AIAA Journal*
   46(9), 2219-2228. DOI [10.2514/1.28943](https://doi.org/10.2514/1.28943).
   The correlation itself: rarefaction-dominated branch for particle Reynolds number below 45 (Cunningham-corrected
   Schiller-Naumann blended with a Reynolds-bridged Schaaf & Chambre free-molecule sphere on Ma^4/(1+Ma^4)), and a
   compression-dominated branch above 45. Abstract: [ADS](https://ui.adsabs.harvard.edu/abs/2008AIAAJ..46.2219L/abstract),
   [Semantic Scholar](https://www.semanticscholar.org/paper/Compressibility-and-Rarefaction-Effects-on-Drag-of-Loth/9f900a8d473f5880f1cf4e11a2a22d73665a2215),
   [ResearchGate](https://www.researchgate.net/publication/233845201_Compressibility_and_Rarefaction_Effects_on_Drag_of_a_Spherical_Particle).

2. **Harrison, A. K. (2021).** Comment on "Compressibility and Rarefaction Effects on Drag of a Spherical Particle".
   *AIAA Journal* 59(8), 3288-3289. DOI [10.2514/1.J060681](https://doi.org/10.2514/1.J060681).
   Reports two errors in the 2008 paper's rarefaction-regime correlation. Accepted manuscript (Los Alamos,
   LA-UR-21-21429): [OSTI](https://www.osti.gov/servlets/purl/1812671). Not yet read: the errors are not
   identified in our notes; read this before implementing Loth (2008).

3. **Loth, E., Daspit, J. T., Jeong, M., Nagata, T., Nonomura, T. (2021).** Supersonic and Hypersonic Drag
   Coefficients for a Sphere. *AIAA Journal* 59(8), 3261-3274. DOI
   [10.2514/1.J060153](https://doi.org/10.2514/1.J060153) (publisher page:
   [arc.aiaa.org](https://arc.aiaa.org/doi/10.2514/1.J060153)).
   The updated model; corrects the 2008 correlation's behaviour near Reynolds number 45 (as reported by source 5).

4. **Singh, N., et al. General Drag Coefficient for Flow over Spherical Particles.** *AIAA Journal*. DOI
   [10.2514/1.J060648](https://doi.org/10.2514/1.J060648)
   (publisher page: [arc.aiaa.org](https://arc.aiaa.org/doi/10.2514/1.J060648); author PDF:
   [Stanford](https://web.stanford.edu/group/ihmegroup/cgi-bin/MatthiasIhme/wp-content/papercite-data/pdf/singh2020drag.pdf);
   [ResearchGate](https://www.researchgate.net/publication/354973102_General_Drag_Coefficient_for_Flow_over_Spherical_Particles)).
   A newer general correlation reported to do better than Loth's at all Mach numbers. Author list, volume and
   pages not yet confirmed (Crossref lookup was rate-limited).

## Secondary source used

5. **arXiv:2510.00479**, appendix B.3, "On the joint estimation of flow fields and particle properties from
   Lagrangian data": [arxiv.org/html/2510.00479](https://arxiv.org/html/2510.00479).
   Restates Loth (2008) in full (equations B.7-B.10). Its free-molecule term prints erf(-s^2) where the physics
   and Schaaf & Chambre require exp(-s^2), and its C_M expressions are typeset ambiguously; check against
   sources 1-3 before use. Also states that Loth's 2008 model "was not empirically supported near Re_p = 45" and
   that Loth et al. (2021) corrects this.

## BibTeX

```bibtex
@article{loth2008compressibility,
  author  = {Loth, E.},
  title   = {Compressibility and Rarefaction Effects on Drag of a Spherical Particle},
  journal = {AIAA Journal},
  volume  = {46},
  number  = {9},
  pages   = {2219--2228},
  year    = {2008},
  doi     = {10.2514/1.28943}
}

@article{harrison2021comment,
  author  = {Harrison, Alan K.},
  title   = {Comment on ``Compressibility and Rarefaction Effects on Drag of a Spherical Particle''},
  journal = {AIAA Journal},
  volume  = {59},
  number  = {8},
  pages   = {3288--3289},
  year    = {2021},
  doi     = {10.2514/1.J060681}
}

@article{loth2021supersonic,
  author  = {Loth, Eric and Daspit, John Tyler and Jeong, Michael and Nagata, Takayuki and Nonomura, Taku},
  title   = {Supersonic and Hypersonic Drag Coefficients for a Sphere},
  journal = {AIAA Journal},
  volume  = {59},
  number  = {8},
  pages   = {3261--3274},
  year    = {2021},
  doi     = {10.2514/1.J060153}
}
```

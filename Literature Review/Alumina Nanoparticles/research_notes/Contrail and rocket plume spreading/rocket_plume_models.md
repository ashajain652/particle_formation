# Rocket plume dispersion modelling and measurements (stratosphere, mesosphere, lower thermosphere): what is comparable to a 50-90 km aluminium vapour trail

Coverage note: about 14 tool calls. Most primary papers were paywalled, so many findings come from abstracts, catalogue records and the local PDFs. Anything not seen in a source is listed under Gaps, not asserted.

## 1. Plume models for stratospheric rocket exhaust: dispersion equations and diffusivities

### Takeaway
The only quantified stratospheric plume-growth data found are Beiting's (Aerospace Corp.) imaging of Shuttle and Titan IV plumes at 18, 24 and 30 km, where plume diameter grows linearly with time (4-9 m/s), not as the square root of time that a constant-diffusivity cylinder gives. I could not retrieve the equations or K values used by Ross, Danilin, Jackman, Lohn, Popp, Newman, Dallas, Larson, Ryan, Maloney or Brown; the sources reached only describe them at abstract level.

### Cited Findings
- Beiting (Aerospace Corp. technical report TR-99(1306)-1, April 1999): plume expansion measured for nine Shuttle and Titan IV launches at 18, 24 and 30 km, from electronic images of polarized near-infrared sunlight scattered by exhaust particles; diameters "increase linearly with time", tracked up to 50 min after the vehicle reached altitude (measured). — [OSTI 20005912](https://www.osti.gov/biblio/20005912)
- Same report: mean diameter expansion rates 4.3 +/- 1.0 m/s at 18 km, 6.8 +/- 1.9 m/s at 24 km, 8.7 +/- 2.5 m/s at 30 km; launch-to-launch scatter up to a factor of 1.6 (18 km), 2.2 (24 km), 2.7 (30 km); expansion rate rose with altitude at Cape Canaveral (six cases) but fell between 24 and 30 km for the one Vandenberg case; no correlation with wind speed or wind shear. — [OSTI 20005912](https://www.osti.gov/biblio/20005912)
- Same report: data were compared with "several models for diffusivity" and used to update a particle model of solid-motor exhaust; the models and comparison results are not given on the abstract page. — [OSTI 20005912](https://www.osti.gov/biblio/20005912)
- Dao et al. 1996 (SPIE 2833), scanning lidar on Titan IV K-21 (6 Nov 1995), STS-76 (22 Mar 1996) and Titan IV K-16 (24 Apr 1996): thin particle layers at 18-43 km persisting up to 3 h after launch, with vertical thickness under 250 m in over 700 profiles (measured). The plume spreads into a horizontally extended, very thin sheet, i.e. strongly anisotropic. — [OSTI 492298](https://www.osti.gov/biblio/492298); [SPIE record](https://biomedicaloptics.spiedigitallibrary.org/conference-proceedings-of-spie/2833/0000/Lidar-measurement-of-TITAN-IV-rocket-exhaust-plumes-at-stratospheric/10.1117/12.258156.full) (search-result summary of the 700-profile and 250 m figures)
- Ross et al. 1997 (Nature 390, 62-64): after two 1996 Titan IV launches, ozone fell to near zero across regions 4-8 km wide within 30 min of launch and recovered to ambient within about 30 min more (measured, aircraft in situ). — [Springer/Nature](https://link.springer.com/article/10.1038/36318), figures via search summary
- A chemical-kinetics plus dispersion model of a Titan IV plume predicted an ozone hole about 10 km in diameter at 20 km altitude and about 28 km at 30 km altitude, lasting hours to a day, using a worst-case and a "best guess" dispersion rate (modelled; the paper is the OSTI 7070382/45889 family, not individually verified). — [search summary, OSTI 492298 neighbours](https://www.osti.gov/biblio/7070382)
- Pellett and Staton (NASA Langley) used a Gaussian multilayer diffusion model for nine solid-motor exhaust cloud cases, but in the troposphere near Cape Canaveral (not stratosphere). — [NTRS 19820006799](https://ntrs.nasa.gov/api/citations/19820006799/downloads/19820006799.pdf), via search summary
- A TRW study (J. Spacecraft and Rockets, 1994) found plume dispersion rate has a very significant effect on local ozone loss. — search summary of [OSTI 45889](https://www.osti.gov/biblio/45889)
- Danilin et al. 2003 used an AER plume model for HNO3 and H2O uptake on alumina in the Athena II plume (WB-57F at about 18 km); Popp et al. 2002 gave an NOy emission index of 2.7 +/- 0.6 g NO2 per kg fuel; Schmid et al. 2003 measured alumina 0.004-1.2 um in the Athena II wake. — [Danilin](https://earthref.org/ERR/71262), [Popp](https://earthref.org/ERR/80492), [Schmid](https://earthref.org/ERR/72383/)
- Brown et al. 2024 (Earth and Space Science): 2019 stratospheric launch emissions about 0.91 Gg alumina and 0.50 Gg reactive chlorine; this is an emissions inventory. — [DOAJ](https://doaj.org/article/0ab52facb92f488da067ab8c10224ce9)

### Inferences
- Linear-in-time diameter growth (4-9 m/s) means a diameter of order 10 km after about 20-40 min. A constant-K Gaussian cylinder gives sigma proportional to sqrt(t), so the model's growth law differs in form; the measurements imply scale-dependent (or shear-driven, along-plume) growth in the stratosphere. This is the closest analogue to the model's K_R(r) law, but it is at 18-30 km where K and N differ from 50-90 km.
- The Titan IV lidar thickness (under 250 m) against widths of kilometres shows horizontal and vertical growth are very unequal in the stratosphere, so an isotropic K is questionable there. Whether this holds at 50-90 km is unknown.
- Converting Beiting's 4-9 m/s into an effective K requires an assumed law (for example K about 0.5 * d(sigma^2)/dt); I did not do this because I could not see how the report defines "diameter".

### Gaps
- Equations and K values (horizontal vs vertical) used by Ross, Danilin, Jackman, Lohn, Newman, Dallas, Larson 2017, Ryan 2022 and Maloney: not retrieved (paywalled or not indexed by search).
- Beiting's comparison with diffusivity models; full text unavailable. Larson 2017 and Ryan 2022 were not found by name in any search.
- Popp 2002 and Danilin 2003 plume widths against age: not in any snippet seen.
- Plume sizes for Athena II at 18 km; WB-57F plume age windows are noted in a search summary only (encounters 5-125 min after launch at 17-19.5 km, Titan IV and Shuttle).

## 2. Measured plume growth from in-situ or imaging

### Takeaway
Measured stratospheric widths: diameter 4-9 m/s growth (Beiting), ozone holes 4-8 km wide at under 30 min (Ross), layers under 250 m thick at 18-43 km (Dao). No measured widths for Falcon 9 mesospheric plumes were found.

### Cited Findings
- See Section 1 for Beiting, Dao and Ross numbers (all measured at 17-43 km).
- WB-57F flew through the Athena II plume six times at about 18 km. — [NOAA CSL](https://csl.noaa.gov/groups/csl2/instruments/palms/imagemaps/accent19990924.html)
- Falcon 9 upper stage, uncontrolled reentry 19 Feb 2025, about 03:42 UTC at about 100 km off western Ireland (52 N, 12 W). Lithium lidar at Kuehlungsborn, Germany saw the plume about 20 h later after advection of about 1600 km. The plume arrived about 00:21 UTC and lasted to the end of sounding at 00:48 UTC, i.e. about 27-40 min of passage (the paper gives both figures in different places). Peak density 31 +/- 8.3 atoms per cm^3 at 96.1 km versus 3.0 +/- 1.1 atoms per cm^3 background; vertical extent 94.5-96.8 km (about 2.3 km thick after 20 h; lidar resolution 200 m with 620 m smoothing) (measured). — local PDF `Wing_2026_Lithium_Plume_Falcon9_Reentry.pdf` ([doi 10.1038/s43247-025-03154-8](https://doi.org/10.1038/s43247-025-03154-8))
- Same paper: about 30 kg lithium in the Al-Li 2198 tank walls, against natural meteoric input of about 80 g per day (estimated). The paper does not give a horizontal width or a diffusivity; it tests only trajectory consistency using UA-ICON winds. — same PDF

### Inferences
- If the Falcon 9 lithium layer was about 2.3 km thick after 20 h at 96 km, this is the nearest measurement to the model's 12-hour end state (radius about 1 km, i.e. diameter about 2 km at 12 h). The thickness is consistent with the model order of magnitude, but lidar smoothing (0.62 km) and trajectory uncertainty mean the match is weak. This is my comparison, not the authors'. Note the lithium source was a reentry at about 100 km, not 50-90 km.
- Horizontal width is only constrained by the passage time (27-40 min at advection speed of about 1600 km / 20 h = about 22 m/s gives roughly 35-50 km along-wind), but wind at the plume level at observation time was not extracted here, so I treat this as unverified arithmetic.

### Gaps
- Space Shuttle SRM and Titan IV widths versus age beyond Beiting's summary; Falcon 9 launch plumes in the mesosphere; rocket-induced ozone holes in the mesosphere: no sources found.

## 3. Mesospheric plumes (shuttle main-engine water, noctilucent clouds, shear spreading)

### Takeaway
Shuttle water plumes released at 100-115 km spread to 3000-4000 km in 18 h, far faster than molecular or modest eddy diffusion would give, and the eddy diffusivity needed to carry water down to 90 km for ice formation was 400-1000 m^2/s, 40-500 times the 2-10 m^2/s of Megner et al. These are above 90 km and so not directly the 50-90 km regime.

### Cited Findings
- STS-135 (8 Jul 2011): three main engines released about 350 t of water vapour at 100-115 km; plume spread horizontally in all directions over 3000-4000 km in 18 h; reached northern Europe after 21 h and produced polar mesospheric clouds brighter than over 99% of those seen in that region (measured by satellite and ground). A 1-D cloud model initialised by a plume diffusion model matched brightness with an eddy diffusion coefficient of 400-1000 m^2/s, which carries water down to about 90 km (fitted to observations). The authors state the horizontal spreading is not well understood. — [Chalmers record for Stevens et al., JGR 117, D19206, 2012](https://research.chalmers.se/en/publication/166331)
- Stevens et al. 2003 (GRL), STS-85 (Aug 1997): bright OH seen by MAHRSI 1-2 days after launch at the Arctic above 85 km; water reached the Arctic in a little over a day, faster than models predicted. — [ScienceDaily](https://www.sciencedaily.com/releases/2003/06/030604084634.htm); [NRL summary](https://webarchive.library.unt.edu/web/20130315182456/http://www.nrl.navy.mil/content_images/07Atmospheric_Stevens.pdf)
- Roughly half of the shuttle main-tank water is injected at 103-114 km (as stated in press coverage). — [ScienceDaily](https://www.sciencedaily.com/releases/2003/06/030604084634.htm)
- Wing 2026 lists the Falcon 9 trajectory-consistent transport as "advection" only; no diffusivity is derived. Concurrent wind shear at the lidar site was described as relatively weak (qualitative). — local PDF `Wing_2026_Lithium_Plume_Falcon9_Reentry.pdf`

### Inferences
- The rapid dispersal of shuttle plumes (many hundreds of km per day) is dominated by large-scale horizontal winds and tidal or gravity-wave shear, not local diffusion. A 12-hour 1 km cylinder with no shear omits the process that made the dominant observed spreading at 100-115 km. Whether shear matters as strongly at 50-90 km is not shown by these sources.

### Gaps
- Kelley, Niciejewski and Meier papers on plume spreading; Niciejewski 2011 (STS-107 plume to Antarctica); measured shear values from the plumes: not retrieved. Meier Lyman-alpha spreading rates: not found.

## 4. Chemical tracer trails (TMA, sodium, lithium) and shear

### Takeaway
Trail expansion between about 85 and 105 km gives effective diffusivities larger than typical eddy values (Bishop et al. 2004, TOMEX), and chemical-release winds show strong shears near and below the 100-110 km wind maximum; no number was retrieved.

### Cited Findings
- Bishop, Larsen, Hecht, Liu, Gardner 2004, JGR 109, D02S03 (TOMEX, White Sands, 26 Oct 2000): TMA trails give horizontal winds, turbulence and diffusivity at about 85-140 km; diffusivities from trail expansion rates were larger than typically assumed eddy coefficients and consistent with molecular diffusion above 110 km (measured; abstract-level, from search result). — [core.ac.uk PDF (not accessible to me)](https://core.ac.uk/download/217151685.pdf)
- Larsen 2002 (JGR, doi 10.1029/2001JA000218): over 400 chemical tracer wind profiles since 1958; wind maximum at 100-110 km at low and mid latitudes; winds above 100 m/s in more than 60% of observations; shears strongest on the bottomside of the maximum and often unstable or near instability. — [EarthRef record / search summary](https://earthref.org/ERR/72620)
- A sodium-lidar study at 80-105 km reported shears above 50 m/s per km (different dataset, not Larsen). — via search summary of [CSU lidar result](https://core.ac.uk/download/51491302.pdf)

### Gaps
- Numerical effective diffusivity values by altitude, shear magnitudes at 50-90 km, and how trails are sheared into ribbons: the full papers could not be read. Chemical release trails below 90 km at all are scarce in what I found.

## 5. Gas-particle chemistry in plumes and dilution

### Takeaway
Plume chemistry studies treat dilution as an input rate (a dispersion parameter), and the one condensation-in-spreading-plume case found is HNO3/H2O uptake on alumina (Danilin 2003). I found no model of alumina or metal vapour condensation in a spreading plume that is comparable to the model.

### Cited Findings
- Danilin et al. 2003: AER plume model used for HNO3 and H2O uptake on Athena II alumina, with surface residence times for HNO3 of 5-52 min and uncertainty of at least a factor of 2 (search summary). — [EarthRef](https://earthref.org/ERR/71262)
- Ferreira et al. 2024 (GRL, local PDF): says earlier work found turbulent diffusion efficient for clusters at 90-95 km (Hunten et al. 1980) but negligible below 80 km, and gives molecular-dynamics alumina cluster growth; it does not model a spreading plume. — local PDF `Ferreira_2024_Satellite_Demise_Ozone_Depletion_MD.pdf`
- Ross et al. 1997: the limited extent and duration of observed ozone loss (4-8 km, under 1 h) argues solid-motor chemistry is not globally significant. — [Springer](https://link.springer.com/article/10.1038/36318)

### Gaps
- Danilin's dilution rule; any alumina nucleation-in-plume model. Not found.

## 6. Reentry plumes: dilution assumptions

### Takeaway
None of the reentry-metal papers checked model plume dilution; they use global grids or observe a plume.

### Cited Findings
- Maloney et al. 2025 (JGR Atmospheres): WACCM6 with CARMA, 10 Gg/yr Al2O3 injected uniformly between 60 and 70 km; accumulation in the polar regions on timescales of months to 2 years; about 1.5 K warming; no interactive ozone chemistry. No plume dilution treatment mentioned in sources seen (abstract-level). — [EGU25-3866](https://meetingorganizer.copernicus.org/EGU25/EGU25-3866.html)
- Murphy 2023 (stratospheric aerosol measurement paper, local PDF), Barker 2024 and Schulz 2021 (local PDFs): searches of the full text for dilution, plume, eddy diffusion and shear terms with numbers returned only instrument-dilution and afterburning (Barker) context; they do not assume a plume dilution. Barker 2024 uses a global 3-D emissions inventory. — local PDFs `Murphy_2023_...`, `Barker_2024_...`, `Schulz_2021_...`
- Ferreira 2024: see Section 5, no dilution. Wing 2026: see Section 2, advection only.

### Inferences
- The model's explicit near-field cylinder with K_eddy 2-10 m^2/s appears to be unmatched in the reentry literature found; the closest published choices are grid-box instant mixing (Maloney, Barker) and tracer advection (Wing).

### Gaps
- Wing 2026 plume dilution of lithium (30 kg injected against measured layer) was not computed. Murphy, Schulz, Barker details were checked only by keyword search of numbered sentences, so a section I missed is possible.

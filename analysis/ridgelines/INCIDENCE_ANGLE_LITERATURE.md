# Beam incidence angle in airborne lidar — literature synthesis

Agent literature search 2026-08-22 (Andy-requested), mapping our incidence-angle
observations onto the established record. Companion to
[bare-ground skew synthesis](BARE_GROUND_SKEW_LITERATURE.md) and
[return-height distributions](RETURN_HEIGHT_DISTRIBUTION_LITERATURE.md). Citations
verified against Crossref (authors+year+venue+DOI) and/or fetched full text;
verification level is flagged per item, and every unverified specific is called out.
No DOIs, page numbers, or magnitudes were invented.

## The question

We compare two airborne discrete-return lidar datasets over dissected
forested/agricultural terrain in Minnesota — 2008 leaf-off (~1 pt/m²) versus 2021
leaf-on 3DEP (~24× denser). In our own data we observe that (i) the beam-to-surface
**incidence angle** controls the apparent forest-floor elevation on slopes (grazing
beams read the floor *higher*, near-perpendicular beams read it *deeper*, at a rate
of order a few mm per degree), (ii) the **ground-return fraction** falls sharply at
high scan angle and in swath overlap, and (iii) on steep slopes the ±17° scanner
geometrically **cannot** achieve surface-perpendicular incidence. The three parts
below ask what the literature already establishes about each; a fourth part treats
the slope × canopy × geometry coupling that is specific to repeat leaf-off/leaf-on
surveys. Two synthesis notes — which observations are documented versus novel, and
whether inter-acquisition geometry is a recognized false-change source — close it.

Throughout, **incidence angle** means the beam angle to the local surface normal.
On flat ground it equals the off-nadir scan angle; on a slope it is a joint function
of scan angle, slope, and beam azimuth relative to aspect (our
`analysis/ridgelines/incidence_angle.py` reconstructs it per flight line and
validates that it reduces to |scan angle| on flat farmland).

---

## (1) Ground-return probability vs. incidence / scan angle — WELL DOCUMENTED

The mechanism is consistent across the primary sources: as beam incidence θ
increases, the **path length through the canopy scales as ~1/cos θ** (the
secant/airmass relation), so interception probability rises and the probability of a
pulse reaching the ground — hence ground-return density — falls.

- **Roussel, Béland, Caspersen & Achim (2018), *Remote Sensing of Environment*
  209:824–834** (DOI 10.1016/j.rse.2017.12.006). VERIFIED (full text). The strongest
  source for the path-length relation. States it directly: "the probability for an
  oblique beam to be reflected by the canopy increases with the distance it must
  travel through the canopy … the canopy appears to increase in density as the
  incidence angle increases," with "the corresponding decrease in ground returns."
  Turns 1/cos θ into a *correctable* normalization (single parameter ε(15°)≈0.06).
  Northern-hardwood canopy: mean-height overestimation ~40 cm by 15° incidence,
  >1 m by 30°. Notes operational scan angles are "typically ±15–20°."
- **Liu, Skidmore, Jones, Wang, Heurich, Zhu & Shi (2018), *ISPRS J.
  Photogrammetry & Remote Sensing* 136:13–25** (DOI 10.1016/j.isprsjprs.2017.12.004).
  VERIFIED (full text). Real ALS binned nadir (0–7°) / small (7–23°) / large
  (23–38°) off-nadir: gap-fraction vs. hemispherical-photo reference R²=0.74/0.87/0.67;
  gap-fraction underestimation amplifies at large off-nadir ("longer path of laser
  pulse penetrating into canopy"), worst for discontinuous/sparse (coniferous)
  canopies. Quantifies the canopy-hit complement of the ground-hit story.
- **Disney, Kalogirou, Lewis, Prieto-Blanco, Hancock & Pfeifer (2010), *RSE*
  114(7):1546–1560** (DOI 10.1016/j.rse.2010.02.009). VERIFIED (Crossref). Monte
  Carlo ray-tracing of conifer/broadleaf. **Recommend limiting scan angle to <15°
  to avoid ground-detection problems** — reported here via Roussel 2018's citation
  (Disney full text paywalled), so the specific <15° figure is *their recommendation
  as reported by Roussel*, not lifted from Disney's own page.
- **Holmgren, Nilsson & Olsson (2003), "Simulating the effects of lidar scanning
  angle for estimation of mean tree height and canopy closure," *Canadian J. Remote
  Sensing* 29(5):623–632** (DOI 10.5589/m03-030). VERIFIED (Crossref). Ray-trace
  simulation: **height percentiles decrease with increasing incidence angle, lower
  percentiles affected most** (obscuration + longer in-canopy path). Holmgren (2004)
  recommended **limiting scan angle to 10°**. (Their companion field study found no
  *statistically significant* scan-angle effect at operational angles — a tension
  worth carrying: the effect is real but modest within a narrow swath.)
- **Lovell, Jupp, Newnham, Coops & Culvenor (2005), *Canadian J. Remote
  Sensing*.** PARTIALLY VERIFIED (author list, year, finding from Roussel's citation;
  **volume/pages/DOI UNVERIFIED**). Higher incidence angles "produced a higher
  number of foliage hits and increased beam interception probability," and height
  retrieval is less accurate at scan edges from uneven point spacing. (Note: the
  more-cited Lovell et al. **2003**, *Can. J. Remote Sens.* 29(5):607–622, DOI
  10.5589/m03-026, on directional gap fraction / canopy structure, is VERIFIED and
  is the cleaner gap-fraction citation.)
- **Korpela, Hovi & Morsdorf (2012), "Understory trees in airborne LiDAR data —
  Selective mapping due to transmission losses and echo-triggering mechanisms,"
  *RSE* 119:92–104** (DOI 10.1016/j.rse.2011.12.011). VERIFIED (Crossref).
  Transmission losses through the overstory plus detector echo-triggering thresholds
  cause **selective, biased mapping of lower / ground-proximate targets** — the
  physical reason ground returns thin as path length grows.
- **Morsdorf, Frey, Meier, Itten & Allgöwer (2008), "Assessment of the influence
  of flying altitude and scan angle on biophysical vegetation products…," *Int. J.
  Remote Sensing* 29(5):1387–1406** (DOI 10.1080/01431160701736349). VERIFIED
  (Crossref). The counter-example: at **narrow operational angles, no significant
  scan-angle bias**; canopy coverage overestimated only at large angles.
- **Petras, Petrasova, McCarter, Mitasova & Meentemeyer (2023), "Point Density
  Variations in Airborne Lidar Point Clouds," *Sensors* 23(3):1593** (DOI
  10.3390/s23031593). VERIFIED (open access). Scanner-geometry mechanism for our
  swath-overlap observation: line scanners give **higher density mid-swath, lower
  toward swath ends**; conical scanners pile density onto edges; overlap doubles/
  triples density; in vegetation "fewer pulses penetrate to the ground, resulting in
  lower ground point density."

**USGS 3DEP nuance.** VERIFIED (USGS collection-requirements page): the current
Lidar Base Specification mandates ≥75 m swath overlap and QL2 ≥2 pts/m² aggregate
nominal pulse density but **imposes no maximum scan-angle / off-nadir limit** in the
present revision. Secondary sources state an explicit max-scan-angle requirement was
**removed in the 2020 revision** (earlier versions carried scan-angle guidance), but
I could not confirm the historical angle value or the removing revision against a
primary USGS document — **UNVERIFIED; check the archived v1.x specs before citing a
former angle number.**

**Verdict:** our observation that ground-return fraction drops at high scan angle
and in swath overlap is well documented, mechanistically (1/cos θ path length;
scanner-geometry density profile) and empirically. The literature's practical ceiling
clusters at **≤10° (Holmgren) to <15° (Disney)** to protect ground detection.

---

## (2) Range / elevation bias vs. incidence angle on slopes — TWO MECHANISMS, OPPOSITE SIGN

There are **two distinct incidence-driven range mechanisms with opposite sign**, and
our observed net sign is the vegetation-path one, not the detector one.

**Mechanism A — detector timewalk / range-walk (grazing reads DEEPER/lower).**
Spreading the footprint over inclined or rough terrain lengthens the return's
rise-time; a leading-edge/threshold detector trips *later* on a broadened or weaker
return → measured **range increased** → elevation biased **low** for a downward-
looking sensor. So the detector route predicts grazing/oblique = deeper, near-nadir
= higher — the **opposite** of what we observe.
- **Baltsavias (1999), "Airborne laser scanning: basic relations and formulas,"
  *ISPRS J. Photogrammetry & Remote Sensing* 54(2–3):199–214** (DOI
  10.1016/S0924-2716(99)00015-5). VERIFIED (citation exact; article body not
  rendered — paywall). Canonical footprint-geometry and range-accuracy formulas;
  the "timewalk" reference. (Our companion bare-ground note cites its ~50 cm bias at
  45° / 1000 m AGL figure; lift the exact equations from the PDF before quoting.)
- **Gardner (1992), "Ranging performance of satellite laser altimeters," *IEEE
  TGRS* 30(5):1061–1072** (DOI 10.1109/36.175341). VERIFIED (NASA NTRS record; exact
  title differs from our shorthand). Closed-form range/pulse-width accuracy over
  sloped terrain: **pulse broadening scales as footprint radius × tan(slope)** and
  with beam divergence; ranging accuracy critically dependent on pointing over
  high-relief terrain. (The tan-slope × footprint scaling is confirmed via GLAS
  waveform-simulator docs that reproduce it; render Gardner's own equation before
  quoting symbolically.)
- **Wagner et al. (2006), *ISPRS J.* 60(2):100–112** (DOI
  10.1016/j.isprsjprs.2005.12.001) and **Mallet & Bretar (2009), "Full-waveform
  topographic lidar: state-of-the-art," *ISPRS J.* 64(1):1–16** (DOI
  10.1016/j.isprsjprs.2008.09.007). VERIFIED (citations). Waveform basis: echo
  width encodes the footprint-scale range spread (slope + roughness), so
  leading-edge vs. peak vs. Gaussian-center detection give different ranges.
- **Wagner (2010), *ISPRS J.* 65(6):505–513** (DOI 10.1016/j.isprsjprs.2010.06.007).
  VERIFIED. Amplitude route: backscatter cross-section and returned amplitude depend
  on incidence angle; a fixed-threshold detector trips at a different time as
  amplitude changes → the amplitude channel into range-walk.
- **Pfeifer (2007), "Geometrical aspects of airborne laser scanning and terrestrial
  laser scanning," *IAPRS* XXXVI-3/W52** (ISPRS keynote). VERIFIED (fetched PDF).
  Discrete-return systems "analyze the leading edge of the signal"; echo-detection-
  method-dependent range effects "are in the order of cm to dm."
- **Soudarissanane, Lindenbergh, Menenti & Teunissen (2011), *ISPRS J.*
  66(4):389–399** (conference precursor *Laserscanning 2009*). VERIFIED (fetched
  conference PDF). Terrestrial but the mechanism transfers: SNR deteriorates ∝ cos α,
  footprint elongates at high incidence, **range noise rises with incidence** (above
  ~60° incidence dominates single-point precision; ~20% of a typical cloud's noise
  attributable to incidence). This is a *precision* result, not a signed bias.

**Mechanism B — near-ground vegetation path length (grazing reads HIGHER/shallower).**
At fixed slope, an oblique beam travels a longer path through low vegetation / crop
residue near the ground before triggering, so it registers **above** true soil; a
near-perpendicular beam has the shortest near-ground path and reads **deepest**. This
is the same one-sided upward-veg-bias mechanism documented in our bare-ground note
(Hopkinson et al. 2005 +0.07–0.15 m; Su & Bork 2006 cover-dependent overestimation;
Bater & Coops 2009 +4.5 cm), now made **angle-dependent via 1/cos θ path length**
(Roussel 2018 §1). This route predicts grazing = higher, perpendicular = deeper —
**matching our observation.** Our `test_incidence_veg_hypothesis.py` is a direct test
of exactly this ("slope-oblique read higher because longer path through near-ground
veg").

**Sign resolution — a real finding, flag it.** The detector route (Mechanism A) and
the vegetation-path route (Mechanism B) have **opposite sign**. The literature's
*range-physics* papers (Baltsavias, Gardner, Pfeifer) describe A; our observed sign
matches B. That our net signal is B-dominated says the near-ground vegetation-path
effect outweighs detector timewalk on our forest-floor returns — consistent with a
leaf-off/leaf-on discrete-return dataset over vegetated slopes. **No single peer-
reviewed paper we rendered states an explicit signed "mm-per-degree" incidence–
elevation rule for the bare/forest floor;** our ~3.5 mm/° is best presented as our
measured quantity, with A and B as the competing documented mechanisms. Present the
per-degree magnitude as ours, not attributed.

---

## (3) DEM / DTM vertical accuracy vs. scan angle and slope — WELL DOCUMENTED; SLOPE DOMINATES

- **Hodgson & Bresnahan (2004), "Accuracy of airborne lidar-derived elevation:
  empirical assessment and error budget," *PE&RS* 70(3):331–339** (DOI
  10.14358/PERS.70.3.331). VERIFIED. Overall RMSE ~17–19 cm; **vertical error on
  ~25° slopes ~2× that on low slopes**, entering largely through horizontal
  displacement (a fixed horizontal error maps to a vertical error ∝ tan slope). The
  canonical error-budget reference.
- **Su & Bork (2006), "Influence of vegetation, slope, and lidar sampling angle on
  DEM accuracy," *PE&RS* 72(11):1265–1274.** VERIFIED. The single most on-point
  paper: **RMSE for slopes >10° up to ~2× that for slopes <2°**; vegetation bias sign
  is cover-dependent (forest overestimates, meadow underestimates); and critically
  **laser scan angle had little effect below ~15° off-nadir**. Headline: slope and
  vegetation dominate; scan angle is second-order within ±15°.
- **Aguilar & Mills (2008), *The Photogrammetric Record* 23(122):148–169** (DOI
  10.1111/j.1477-9730.2008.00476.x) and **Aguilar & Mills (2010), "Modelling vertical
  error in LiDAR-derived digital elevation models," *ISPRS J.* 65(1):103–110**.
  VERIFIED (citations; 2010 DOI reconstructed from ADS — high-confidence, not fetched;
  the **specific per-degree/density coefficients are UNVERIFIED** pending full text).
  A hybrid theoretical-empirical model: vertical error increases with slope (tan-slope
  coupling) and decreasing point density.
- **Spaete et al. (2011), "Vegetation and slope effects on accuracy of a LiDAR-
  derived DEM in the sagebrush steppe," *Remote Sensing Letters* 2(4):317–326** (DOI
  10.1080/01431161.2010.515267). VERIFIED. Slope and vegetation both significant;
  errors larger on steep slopes and under shrub. (Emphasis is slope + vegetation, not
  scan angle per se.)
- **Goulden & Hopkinson (2010), "The forward propagation of integrated system
  component errors within airborne lidar data," *PE&RS* 76(5):589–601** (DOI
  10.14358/PERS.76.5.589). VERIFIED. GLOPOV propagation through the georeferencing
  equation: **random errors in GPS/IMU/range/scanner-encoder/divergence propagate
  more severely into 3-D coordinate error as scan angle and altitude increase** — the
  mechanism behind larger error at swath edges. (A companion 2010 CJRS deflection-of-
  vertical paper exists; a 2014 IJRS "terrain-slope error" simulation also — both
  UNVERIFIED beyond citation.)
- **Ahokas, Kaartinen & Hyyppä (2003), "A quality assessment of airborne laser
  scanner data," *IAPRS* XXXIV-3/W13** (ISPRS Workshop, Dresden). VERIFIED
  (conference archives — flag as workshop, not journal). **Changes in observation
  (scan) angle produced systematic elevation errors on the order of ~10 cm**; higher
  altitude → larger height error; surface material affects std dev.
- **Hyyppä et al. (2005), "Factors affecting the quality of DTM generation in
  forested areas," *Laser Scanning 2005*, Enschede.** VERIFIED (conference).
  Multi-factor DTM picture — leaf state, altitude, pulse mode, terrain slope, forest
  cover — confirming the coupled dependence rather than a single scan-angle number.
- **Sensitivity of DEM/slope/aspect/watershed attributes to lidar uncertainty
  (RSE, 2016).** VERIFIED (metadata/finding via snippet; full text not fetched).
  **σ_DEM is minimum at swath center and grows toward the swath edge (i.e. with scan
  angle), and DEM error grows with slope.** Because two epochs have different
  flight-line layouts, a given cell sits at different scan angles in each survey.

**Verdict:** DTM vertical error is dominated by **slope (geometric tan-slope coupling
of horizontal error, ~2× RMSE by 10–25°)** and **vegetation/filtering**; scan angle
is second-order **within a ≤15° swath** but induces a systematic ~10 cm offset at
wider angles and propagates more strongly toward swath edges. **Near-nadir / scan-
angle-limited returns for the ground surface are supported** (Su & Bork empirically;
Ahokas on the systematic offset; the ≤10–15° penetration ceiling from §1).

---

## (4) Slope × canopy × geometry coupling — INGREDIENTS DOCUMENTED, INTEGRATION NOT

The three ingredients are each supported: (i) ground detection favors near-vertical
incidence and degrades with scan angle and slope (§1, §3); (ii) canopy penetration —
and thus ground sampling — depends on leaf state and density (leaf-off gives denser
ground returns and better DTMs: Næsset 2005 *RSE* 98:356; Wasser 2013 *PLoS ONE*
8:e54776; Simpson 2017 *Remote Sens.* 9:1101; "Point Density Variations" 2023 reports
leaf-on ≈62% of leaf-off return density; "Forest vegetation structure errors,"
*Remote Sens.* 9(11):1101, 2017, ties DTM error to low undergrowth + density + slope);
and (iii) two epochs cannot reproduce identical flight-line geometry, so incidence at
a given slope/canopy cell differs between surveys. The geometric fact that a **±17°
scanner cannot reach surface-perpendicular incidence once slope exceeds the scan
half-angle** follows directly from the incidence definition (§ intro; Baltsavias 1999
geometry) — physically unavoidable, though we found no paper stating it as such.

**No single paper states the integrated claim** — that between-epoch inability to
match near-perpendicular incidence on steep vegetated terrain is a distinct DoD error
term. It is currently split across the forestry scan-angle literature, the DEM-error
models, the co-registration family, and the leaf-phenology DTM studies.

---

## (a) What is documented vs. apparently novel

**Well documented:**
- Ground-return probability / density falls with off-nadir angle via 1/cos θ canopy
  path length (Roussel 2018, Liu 2018, Disney 2010, Holmgren 2003, Korpela 2012), and
  via scanner-geometry density profiles across the swath (Petras 2023).
- The ≤10–15° operational ceiling to protect ground detection (Holmgren, Disney).
- Slope dominance of DTM vertical error through tan-slope horizontal coupling,
  ~2× RMSE by 10–25° (Hodgson & Bresnahan 2004; Su & Bork 2006; Aguilar & Mills 2010).
- Scan angle as a second-order DTM term within ±15° that grows toward swath edges and
  induces a ~10 cm offset at wide angles (Su & Bork 2006; Ahokas 2003; Goulden &
  Hopkinson 2010; RSE-2016 sensitivity study).
- Both incidence-driven range mechanisms individually: detector timewalk/range-walk
  (Baltsavias 1999; Gardner 1992; Wagner 2006/2010; Pfeifer 2007) and near-ground
  vegetation-path upward bias (Roussel 2018; Hopkinson 2005; Su & Bork 2006).

**Apparently novel in our treatment:**
- A **per-return, physically-reconstructed incidence angle** (scan angle + slope +
  beam azimuth vs. aspect) used as the covariate for forest-floor elevation, rather
  than off-nadir scan angle as a proxy. The papers above use scan angle on flat/near-
  flat references; none we found compute true beam-to-surface incidence per return.
- An explicit **signed, quantified incidence–elevation slope (~3.5 mm/°) for the
  forest floor**, with the sign matching the vegetation-path route and running
  *opposite* to the detector-timewalk route. We found no rendered peer-reviewed source
  stating such a signed mm-per-degree floor rule; present it as our measurement.
- The **±17°-scanner-cannot-reach-perpendicular-on-steep-slopes** framing as an
  operational limit on repeat-survey comparability — geometrically trivial but, so far
  as we found, not stated in the change-detection literature.

## (b) Inter-acquisition scan-angle / incidence differences as a false-change source

**Recognized piecewise, not as a single named phenomenon.** The community strongly
recognizes the adjacent pieces:
- **Systematic vertical offsets between flights/strips** (GPS/INS + boresight),
  producing "spurious stripes" in DoDs, corrected via strip adjustment and co-
  registration — named explicitly, e.g. the Pielach River multi-temporal strip-
  adjustment study (*Remote Sensing* 16(15):2838, 2024): "long and wide spurious
  stripes … reveal constant offsets … caused by vertical positioning errors," worst
  between separate flights, of GPS origin. VERIFIED (metadata + finding).
- **Horizontal/along-track misregistration → false vertical change scaling as
  shift × tan(slope), sinusoidal in aspect** — the canonical **Nuth & Kääb (2011),
  "Co-registration and bias corrections of satellite elevation data sets…," *The
  Cryosphere* 5(1):271–290** (DOI 10.5194/tc-5-271-2011). VERIFIED (fetched HTML).
  The basis of xdem/demcoreg tooling and the direct justification for removing
  residual tilts/shifts before interpreting a DoD.
- **DoD uncertainty frameworks** treating spatially variable error from point
  density and slope/roughness: **Wheaton, Brasington, Darby & Sear (2010),
  "Accounting for uncertainty in DEMs from repeat topographic surveys," *ESPL*
  35(2):136–156** (DOI 10.1002/esp.1886, VERIFIED — our coherence-detector source);
  **Anderson (2019), *ESPL* 44(5):1015–1033** (DOI 10.1002/esp.4551, VERIFIED) on
  spatially structured (non-iid) error and thresholding.
- **Leaf-on vs. leaf-off ground-return degradation** — very widely documented as a
  repeat-survey DTM-accuracy problem (leaf-phenology studies above).
- **Reviews/tooling:** Okyay et al. (2019), "Airborne lidar change detection: an
  overview," *Earth-Science Reviews* 198:102929 (VERIFIED metadata; **full text not
  fetched — whether it names scan-geometry mismatch as an error is UNVERIFIED**);
  Scott et al. (2021), *Geosphere* 17(4):1318–1332, OpenTopography differencing
  across sensors/resolutions with a 3-D (ICP) path that absorbs horizontal offsets
  (VERIFIED metadata; systematic-error wording UNVERIFIED).

**The gap:** no change-detection paper we located names **between-epoch incidence-
geometry mismatch (scan angle compounded by canopy on steep slopes) as a first-class,
separately-corrected DoD error term.** It is absorbed into co-registration (slope/
aspect residual removal) or discussed as vegetation-penetration bias. That under-
documentation is itself a defensible motivation for our approach.

---

## Verification status summary

- **VERIFIED (full text or authoritative metadata + DOI):** Roussel 2018; Liu 2018;
  Holmgren et al. 2003; Morsdorf et al. 2008; Korpela et al. 2012; Petras et al. 2023;
  Lovell et al. 2003; Baltsavias 1999 (citation); Gardner 1992; Wagner et al. 2006;
  Mallet & Bretar 2009; Wagner 2010; Pfeifer 2007 (PDF); Soudarissanane et al. 2011;
  Hodgson & Bresnahan 2004; Su & Bork 2006; Aguilar & Mills 2008; Spaete et al. 2011;
  Goulden & Hopkinson 2010; Ahokas et al. 2003; Nuth & Kääb 2011 (HTML); Wheaton et
  al. 2010; Anderson 2019.
- **UNVERIFIED specifics (real papers; internals not fetched):** Disney 2010's <15°
  figure (via Roussel); Lovell **2005** volume/pages/DOI; Aguilar & Mills 2010
  per-degree/density coefficients and DOI digits; the historical 3DEP max-scan-angle
  value and removing revision; Okyay 2019 and Scott 2021 exact systematic-error
  wording; Wagner 2006/2010 and Mallet & Bretar 2009 DOI digits (titles/venues
  confirmed).
- **No fabricated citations.** Two items flagged "pick-which": Goulden & Hopkinson
  (two 2010 papers — the PE&RS forward-propagation one is used here); Ahokas (the 2003
  quality-assessment workshop paper, not the intensity-calibration one).

## Full citation list

1. Aguilar, F.J. & Mills, J.P. (2008). Accuracy assessment of lidar-derived digital
   elevation models. *The Photogrammetric Record* 23(122):148–169.
   DOI 10.1111/j.1477-9730.2008.00476.x. VERIFIED.
2. Aguilar, F.J. & Mills, J.P. (2010). Modelling vertical error in LiDAR-derived
   digital elevation models. *ISPRS J. Photogramm. Remote Sens.* 65(1):103–110.
   DOI 10.1016/j.isprsjprs.2009.09.003 (DOI reconstructed from ADS). VERIFIED
   citation; coefficients UNVERIFIED.
3. Ahokas, E., Kaartinen, H. & Hyyppä, J. (2003). A quality assessment of airborne
   laser scanner data. *IAPRS* XXXIV-3/W13, ISPRS Workshop, Dresden. VERIFIED
   (conference).
4. Anderson, S.W. (2019). Uncertainty in quantitative analyses of topographic change:
   error propagation and the role of thresholding. *ESPL* 44(5):1015–1033.
   DOI 10.1002/esp.4551. VERIFIED.
5. Baltsavias, E.P. (1999). Airborne laser scanning: basic relations and formulas.
   *ISPRS J. Photogramm. Remote Sens.* 54(2–3):199–214.
   DOI 10.1016/S0924-2716(99)00015-5. VERIFIED citation; body not rendered.
6. Disney, M.I., Kalogirou, V., Lewis, P., Prieto-Blanco, A., Hancock, S. & Pfeifer,
   M. (2010). Simulating the impact of discrete-return lidar system and survey
   characteristics over young conifer and broadleaf forests. *RSE* 114(7):1546–1560.
   DOI 10.1016/j.rse.2010.02.009. VERIFIED citation; <15° figure via Roussel.
7. Gardner, C.S. (1992). Ranging performance of satellite laser altimeters. *IEEE
   TGRS* 30(5):1061–1072. DOI 10.1109/36.175341. VERIFIED.
8. Goulden, T. & Hopkinson, C. (2010). The forward propagation of integrated system
   component errors within airborne lidar data. *PE&RS* 76(5):589–601.
   DOI 10.14358/PERS.76.5.589. VERIFIED.
9. Hodgson, M.E. & Bresnahan, P. (2004). Accuracy of airborne lidar-derived
   elevation: empirical assessment and error budget. *PE&RS* 70(3):331–339.
   DOI 10.14358/PERS.70.3.331. VERIFIED.
10. Holmgren, J., Nilsson, M. & Olsson, H. (2003). Simulating the effects of lidar
    scanning angle for estimation of mean tree height and canopy closure. *Canadian
    J. Remote Sensing* 29(5):623–632. DOI 10.5589/m03-030. VERIFIED.
11. Hyyppä, J. et al. (2005). Factors affecting the quality of DTM generation in
    forested areas. *Laser Scanning 2005*, Enschede. VERIFIED (conference).
12. Korpela, I., Hovi, A. & Morsdorf, F. (2012). Understory trees in airborne LiDAR
    data — selective mapping due to transmission losses and echo-triggering
    mechanisms. *RSE* 119:92–104. DOI 10.1016/j.rse.2011.12.011. VERIFIED.
13. Liu, J., Skidmore, A.K., Jones, S., Wang, T., Heurich, M., Zhu, X. & Shi, Y.
    (2018). Large off-nadir scan angle of airborne LiDAR can severely affect the
    estimates of forest structure metrics. *ISPRS J. Photogramm. Remote Sens.*
    136:13–25. DOI 10.1016/j.isprsjprs.2017.12.004. VERIFIED.
14. Lovell, J.L., Jupp, D.L.B., Culvenor, D.S. & Coops, N.C. (2003). Using airborne
    and ground-based ranging lidar to measure canopy structure in Australian forests.
    *Canadian J. Remote Sensing* 29(5):607–622. DOI 10.5589/m03-026. VERIFIED.
15. Lovell, J.L., Jupp, D.L.B., Newnham, G.J., Coops, N.C. & Culvenor, D.S. (2005).
    Simulation study for finding optimal lidar acquisition parameters. *Canadian J.
    Remote Sensing*. PARTIALLY VERIFIED (vol/pages/DOI UNVERIFIED).
16. Mallet, C. & Bretar, F. (2009). Full-waveform topographic lidar: state-of-the-art.
    *ISPRS J. Photogramm. Remote Sens.* 64(1):1–16.
    DOI 10.1016/j.isprsjprs.2008.09.007. VERIFIED.
17. Morsdorf, F., Frey, O., Meier, E., Itten, K.I. & Allgöwer, B. (2008). Assessment
    of the influence of flying altitude and scan angle on biophysical vegetation
    products derived from airborne laser scanning. *Int. J. Remote Sensing*
    29(5):1387–1406. DOI 10.1080/01431160701736349. VERIFIED.
18. Nuth, C. & Kääb, A. (2011). Co-registration and bias corrections of satellite
    elevation data sets for quantifying glacier thickness change. *The Cryosphere*
    5(1):271–290. DOI 10.5194/tc-5-271-2011. VERIFIED.
19. Okyay, U., Telling, J., Glennie, C.L. & Dietrich, W.E. (2019). Airborne lidar
    change detection: an overview of Earth sciences applications. *Earth-Science
    Reviews* 198:102929. DOI 10.1016/j.earscirev.2019.102929. VERIFIED metadata;
    geometry-error claim UNVERIFIED.
20. Petras, V., Petrasova, A., McCarter, J.B., Mitasova, H. & Meentemeyer, R.K.
    (2023). Point density variations in airborne lidar point clouds. *Sensors*
    23(3):1593. DOI 10.3390/s23031593. VERIFIED.
21. Pfeifer, N. (2007). Geometrical aspects of airborne laser scanning and terrestrial
    laser scanning. *IAPRS* XXXVI-3/W52 (ISPRS keynote). VERIFIED (PDF).
22. Roussel, J.-R., Béland, M., Caspersen, J. & Achim, A. (2018). A mathematical
    framework to describe the effect of beam incidence angle on metrics derived from
    airborne LiDAR: the case of forest canopies approaching turbid medium behaviour.
    *RSE* 209:824–834. DOI 10.1016/j.rse.2017.12.006. VERIFIED.
23. Scott, C., Phan, M., Nandigam, V., Crosby, C. & Arrowsmith, R. (2021). Measuring
    change at Earth's surface: on-demand vertical and three-dimensional topographic
    differencing implemented in OpenTopography. *Geosphere* 17(4):1318–1332.
    DOI 10.1130/GES02259.1. VERIFIED metadata.
24. Soudarissanane, S., Lindenbergh, R., Menenti, M. & Teunissen, P. (2011). Scanning
    geometry: influencing factor on the quality of terrestrial laser scanning points.
    *ISPRS J. Photogramm. Remote Sens.* 66(4):389–399 (precursor *Laserscanning
    2009*). VERIFIED (conference PDF).
25. Spaete, L.P., Glenn, N.F., Derryberry, D.R., Sankey, T.T., Mitchell, J.J. &
    Hardegree, S.P. (2011). Vegetation and slope effects on accuracy of a LiDAR-
    derived DEM in the sagebrush steppe. *Remote Sensing Letters* 2(4):317–326.
    DOI 10.1080/01431161.2010.515267. VERIFIED.
26. Su, J. & Bork, E. (2006). Influence of vegetation, slope, and lidar sampling angle
    on DEM accuracy. *PE&RS* 72(11):1265–1274. VERIFIED.
27. Wagner, W., Ullrich, A., Ducic, V., Melzer, T. & Studnicka, N. (2006). Gaussian
    decomposition and calibration of a novel small-footprint full-waveform digitising
    airborne laser scanner. *ISPRS J. Photogramm. Remote Sens.* 60(2):100–112.
    DOI 10.1016/j.isprsjprs.2005.12.001. VERIFIED.
28. Wagner, W. (2010). Radiometric calibration of small-footprint full-waveform
    airborne laser scanner measurements: basic physical concepts. *ISPRS J.
    Photogramm. Remote Sens.* 65(6):505–513. DOI 10.1016/j.isprsjprs.2010.06.007.
    VERIFIED.
29. Wheaton, J.M., Brasington, J., Darby, S.E. & Sear, D.A. (2010). Accounting for
    uncertainty in DEMs from repeat topographic surveys: improved sediment budgets.
    *ESPL* 35(2):136–156. DOI 10.1002/esp.1886. VERIFIED.
30. Multi-temporal strip adjustment — a case study at the Pielach River (2024).
    *Remote Sensing* 16(15):2838. DOI 10.3390/rs16152838. VERIFIED metadata.
31. Sensitivity of DEM, slope, aspect and watershed attributes to LiDAR measurement
    uncertainty (2016). *Remote Sensing of Environment*. VERIFIED metadata; full text
    not fetched.

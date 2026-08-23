# Ground-elevation bias vs. beam-to-slope incidence — vegetation-path mechanism

Agent literature search 2026-08-15 (Andy-requested). Tightly scoped companion to
[INCIDENCE_ANGLE_LITERATURE.md](INCIDENCE_ANGLE_LITERATURE.md) (broad scan-angle
review), [BARE_GROUND_SKEW_LITERATURE.md](BARE_GROUND_SKEW_LITERATURE.md), and
[RETURN_HEIGHT_DISTRIBUTION_LITERATURE.md](RETURN_HEIGHT_DISTRIBUTION_LITERATURE.md).
This file does **one** thing: test the specific, *signed* mechanism below against the
record, without re-treading the general scan-angle material. Citations verified per
item against fetched full text or authoritative metadata; every unverified specific is
flagged. No DOIs, magnitudes, or attributions were invented. Two attribution errors
carried in the companion files are corrected here (see the boxed notes).

## The mechanism under test

In airborne lidar over vegetated terrain, a pulse striking the surface more nearly
**perpendicular to the local slope** travels a *shorter* path through the near-ground
vegetation layer (understory, low branches, litter, crop residue), so it is more
likely to penetrate to true bare ground and records a **lower (deeper), truer** ground
elevation. A pulse at *oblique* incidence to the surface travels a *longer* near-ground
path (path length ~ layer_thickness / cos(incidence)), is more likely to stop **on** low
vegetation, and records a **higher (shallower)** apparent ground elevation. Net: as
beam-to-**local-surface** incidence rises, measured ground elevation is biased **upward**;
near-slope-perpendicular beams give the truest (lowest) ground. Because a bounded
scanner (ours ±17°) cannot achieve near-perpendicular incidence once slope exceeds the
scan half-angle, the upward bias grows with slope. Our per-return gen1 (2008
discrete-return, dissected forested MN) measurement: near-slope-perpendicular ground
returns read the floor consistently deeper than the full mixed-angle set, the gap
widening ~6 mm (0–8° slope) to ~35 mm (16–24°).

Throughout, **incidence** = beam angle to the *local surface normal* (scan angle, slope,
and beam azimuth-vs-aspect jointly), not off-nadir scan angle. This distinction is the
crux: the published physics is almost entirely framed on scan angle against a flat
reference.

---

## (1) Path-length-through-vegetation control on reaching ground (∝ 1/cos θ) — WELL DOCUMENTED (for probability/density; on scan angle)

The 1/cos θ near-ground path-length control on whether a pulse reaches the ground is
established, but as a control on ground-return **probability/density**, framed on the
**off-nadir/beam angle**, not on the local slope.

- **Roussel, Béland, Caspersen & Achim (2018), *Remote Sensing of Environment*
  209:824–834** (DOI 10.1016/j.rse.2017.12.006). VERIFIED — full author post-print read.
  The strongest path-length source. Verbatim: *"When a beam arrives at an angle of θ
  degrees, its travel distance through each layer is 1/cos(θ) times longer than that of
  a vertical beam. Thus, the probability of interacting with canopy elements increases
  with the incidence angle …"*, and increasing incidence *"[shifts] the expected height
  distribution upwards (increasing the number of canopy while decreasing the number of
  ground points)."* Beer–Lambert / turbid-medium interception; the model normalizes the
  point vertical distribution *"as if all data had been sampled at-nadir."* Northern-
  hardwood canopy: mean-height overestimation ~40 cm by 15°, >1 m by 30°.
- **MacArthur & Horn (1969), *Ecology* 50:802–804** and **Ni-Meister, Jupp & Dubayah
  (2001), *IEEE TGRS* 39(9):1943–1958** (GORT). VERIFIED (companion RETURN_HEIGHT file).
  The Poisson gap-probability basis: P_gap(z) = exp(−G·LAI(z)/cos θ) — the cos θ in the
  denominator *is* the slant-path lengthening. Ni-Meister makes the exponential
  extinction-with-depth explicit. These give the mechanism its analytic form.
- **Liu et al. (2018), *ISPRS J. P&RS* 136:13–25** (DOI 10.1016/j.isprsjprs.2017.12.004);
  **Korpela, Hovi & Morsdorf (2012), *RSE* 119:92–104** (DOI 10.1016/j.rse.2011.12.011).
  VERIFIED. Empirical/analytic support that the canopy-hit complement grows and ground/
  low-target sampling thins as the in-canopy path lengthens (Liu is CANOPY gap-fraction/
  LAI despite a DTM-sounding title — do not cite it as a ground result).

> **CORRECTION — carried error in the companion file.** The parameter **ε(15°) ≈ 0.06**
> in Roussel (2018) is the canopy **absorption** parameter (fraction of energy absorbed,
> the model's single empirical constant), **not** a ground-return-probability
> normalization. The companion INCIDENCE_ANGLE file attaches ε(15°)≈0.06 to the "1/cos θ
> ground-return decrease as correctable"; that is a conflation. The 1/cos θ term drives
> the *canopy interaction* rise; the points-per-pulse decrease with angle Roussel reports
> as a **separate, mechanistically uncertain, minor** add-on (*"To our knowledge, this
> phenomenon has not yet been reported"*). Do not cite ε as a ground normalization.

**Verdict (1):** the path-length (1/cos θ) control on reaching ground is well documented
as a **probability/density** effect on **scan/beam angle**. Its transposition to
**local-slope** incidence is geometrically identical but not itself published.

---

## (2) The signed ELEVATION bias (oblique reads HIGHER) — DOCUMENTED on FLAT ground; sign OPPOSITE to detector timewalk

The *elevation* consequence — pulses stopping on low vegetation register ground **high** —
is well documented with magnitudes, and its sign is the opposite of detector range-walk.

- **Hopkinson, Chasmer, Zsigovics, Creed, Sitar, Treitz & Maher (2005), *IAPRS*
  XXXVI-8/W2:108–113.** VERIFIED — primary PDF read. Ground-classified returns lay
  **above** GPS in vegetation: aquatic **+0.15 m** (σ 0.22), all-vegetated mean
  **+0.07 m** (±0.16). Verbatim mechanism: *"In the absence of strong laser backscatter
  from the true ground surface, laser returns will therefore tend to be biased upwards
  into the more highly reflective foliage."* They **sampled flat areas**, cleanly
  isolating penetration from the geometric tan-slope term — but consequently do **not**
  condition it on slope.
- **Ewald (2013), M.S. thesis, Oregon State University** ("Where's the Ground Surface?
  Elevation Bias in LIDAR-derived DEMs Due to Dense Vegetation in Oregon Tidal Marshes").
  VERIFIED — primary PDF read. *Thesis, not peer-reviewed — flag.* DEM ground
  **positively biased** under dense vegetation: mean **+4.5 cm** (1.4 m cell), 10–30 cm
  typical, up to **+36.6 cm** (*Carex obnupta*) and **+48.8 cm** (*Carex lyngbyei*);
  open-terrain RMSE 4.5 cm with no consistent bias. Reviews the *Spartina* tidal-marsh
  literature (Schmid 2011; Chassereau 2011; Hladik & Alber 2012; Wang 2009) reporting
  **+10 to +45 cm** lidar-above-GPS from non-penetration (verified as summarized by
  Ewald, not at source).
- **Töyrä, Pietroniro, Hopkinson & Kalbfleisch (2003), *Can. J. Remote Sensing*
  29(6):679–690.** Ground overestimation **+0.07 m** (graminoid) to **+0.15 m** (willow
  scrub). VERIFIED only *as quoted in* Hopkinson (2005); original not opened.
- **Su & Bork (2006), *PE&RS* 72(11):1265–1274.** VERIFIED. Cover-dependent **sign**:
  *"elevations were over-estimated in forest but under-estimated within meadow habitats"*
  — the upward bias is specifically the tall/dense-cover (foliage-hit) case.

**Detector timewalk / range-walk — OPPOSITE SIGN, correctly excluded.**
- **Baltsavias (1999), *ISPRS J. P&RS* 54(2–3):199–214** (DOI 10.1016/S0924-2716(99)00015-5);
  **Gardner (1992), *IEEE TGRS* 30(5):1061–1072** (DOI 10.1109/36.175341). VERIFIED
  (citations). Spreading the footprint over inclined/rough terrain lengthens rise-time;
  a leading-edge detector trips *later* → range *increased* → elevation biased **low**
  (grazing reads *deeper*). Pulse broadening ∝ footprint × tan(slope) (Gardner). This is
  the **opposite** sign to our observation, and confirms the two mechanisms are distinct.

**Sign resolution:** our observed sign (oblique/high-slope-incidence reads *higher*)
matches the vegetation-path route, not the detector route — so on our forest-floor
returns the near-ground-vegetation effect dominates detector timewalk.

> **CORRECTION — carried error in the companion/skew files.** The **+4.5 cm** figure is
> **Ewald (2013)**, *not* Bater & Coops (2009); and the **+0.07/+0.15 m graminoid/scrub**
> pair is **Töyrä et al. (2003)** (quoted in Hopkinson), not Hopkinson's own class result
> (his own is aquatic +0.15 m, all-veg mean +0.07 m). Re-attribute accordingly.

**Verdict (2):** the signed upward ground bias from stopping on low vegetation is
**documented, with magnitudes**, but on **flat / near-flat** terrain, and framed by
cover class — **not** as a function of beam-to-slope **incidence angle**.

---

## (3) Explicit relation to SLOPE — TWO mechanisms; the geometric one is documented, the penetration one is entangled

Two slope→ground-error mechanisms exist and must be kept apart.

- **Geometric tan-slope (documented, clean, quantitative).** A fixed horizontal error
  maps to vertical error ∝ tan(slope). **Hodgson & Bresnahan (2004), *PE&RS*
  70(3):331–339** (DOI 10.14358/PERS.70.3.331): elevation error on ~25° slopes ~**2×**
  that on ~1.5° slopes, entering through horizontal displacement (verified via
  concordant secondaries; a citing paper restates *"100 cm horizontal error on a 10°
  slope → up to 18 cm elevation error"*, i.e. 100·tan10° = 17.6 cm). **Su & Bork (2006)**:
  RMSE for slopes >10° roughly **double** that for <2°. This is **not** the vegetation
  route — it is geometry, and it is unsigned/scatter-like.
- **Vegetation-penetration-on-slope (documented as physics, NOT cleanly isolated).**
  **Mohd Salleh, Ismail & Abdul Rahman (2015), *ISPRS Annals* II-2/W2:183–189**. VERIFIED
  — primary PDF read. The **single closest paper**, naming *both* mechanisms side by
  side: the tan-slope geometric term (citing Hodgson & Bresnahan), then *"high density
  of vegetation also reduce the number of LiDAR ground points due to less signal can
  penetrate the canopy and reach to the ground along with the slope"* — and, citing
  **Lewis & Hancock (2007)**, that on a slope a ground return *"is at a higher altitude,"*
  complicating filtering. Slope correlates with DTM RMSE r≈0.87–0.99, comparable to
  canopy cover. But slope and canopy are **entangled** in their RMSE — no isolated,
  penetration-only, *signed* slope-bias magnitude, and RMSE is not a signed bias.
- **Goodwin et al. (2007)** (cited by Roussel): larger incidence *"produced a higher
  number of foliage hits and increased beam interception probability"* — the foliage-hit
  half of the mechanism. VERIFIED as cited, original not opened.

**Verdict (3):** slope is explicitly tied to ground error, but the documented, clean,
quantitative slope mechanism is the **geometric tan-slope** one. The
slope→higher-incidence→longer-near-ground-path→**signed upward** penetration bias is
stated qualitatively (Mohd Salleh 2015; Lewis & Hancock 2007) but **never isolated with
vegetation held fixed**, and never as a per-return incidence quantity.

---

## (4) Recommendations to select near-nadir / normalize by path length — GROUND recs are thin and partly NEGATIVE; a 1/cos ground normalization is a GAP

> **CORRECTION — carried error in the companion file.** **Disney et al. (2010)** and
> **Holmgren (2003/2004)** are **CANOPY-metric** papers; their "<15°" / "≤10°" scan-angle
> limits protect **tree-height / canopy** estimation (crown shadowing), *not* ground/DTM
> detection. The "<15° to avoid ground-detection problems" framing traces to an
> untraceable paraphrase and should **not** be cited for a ground argument. The Holmgren
> "≤10°" verbatim could not be confirmed from the Holmgren papers themselves (paywalled);
> it is safer attributed to the Nordic acquisition-guideline literature. Re-cite the
> ground/DTM recommendations to the papers below.

- **Ahokas, Yu, Oksanen, Hyyppä, Kaartinen & Hyyppä (2005), "Optimization of the
  Scanning Angle for Countrywide Laser Scanning," *ISPRS Laser Scanning 2005*, Enschede.**
  VERIFIED — full text read. The correct **ground/DTM** near-nadir citation: *"at
  scanning angles of more than 10 degrees off-nadir, the amount of shadowed area
  increases substantially, i.e., the number of measured ground hits decreases and gaps
  in the DTM occur more frequently"*; scan angles "up to 15 degrees seems to be usable."
- **Estornell, Ruiz, Velázquez-Martí & Hermosilla (2011), *Int. J. Digital Earth*
  4(6):521–538.** PARTIALLY VERIFIED. Steep shrub: *"the likelihood of obtaining laser
  pulse returns from the ground can increase with narrow scan angles."*
- **Roussel et al. (2018)** proposes *correcting* (normalize-to-nadir) rather than
  discarding oblique data — but for **canopy** metrics, not ground/DTM.
- **Honest counter-evidence (near-nadir hypothesis is not universally supported for
  elevation).** **Su & Bork (2006)**: sampling angle had *little* effect on DEM error
  below ~15°; slope and vegetation dominated. **Ahokas et al. (2011), *Remote Sensing*
  3(7):1365–1379** (DOI 10.3390/rs3071365): over 0–15° in boreal forest, scan angle had
  **no significant** effect on ground hits/transmittance. Our per-return signal runs
  partly against these pooled scan-angle nulls — which is expected, since they use
  scan angle on flat-ish references, not per-return incidence to the local slope.
- **Per-return local-incidence geometry exists — but built for intensity, not elevation.**
  **Tan, Cheng & Zhao (2021), *Remote Sensing* 13(3):511** (DOI 10.3390/rs13030511)
  computes the beam-to-surface-normal incidence from POS + point coordinates (slope +
  aspect + scan geometry, exactly our construction) and warns it can *over*correct on
  steep slopes — but the target is **intensity** for land-cover classification.
  **Höfle & Pfeifer (2007), *ISPRS J. P&RS* 62(6):415–433** and **Soudarissanane et al.
  (2011), *ISPRS J. P&RS* 66(4):389–399** use incidence (cos α) for intensity radiometry
  and TLS range precision respectively — geometry identical, target not ground DTM.

**Verdict (4):** ground/DTM near-nadir recommendations exist (Ahokas 2005; Estornell
2011) but are qualitative and partly contradicted at pooled scan angles (Su & Bork;
Ahokas 2011). A **1/cos-path-length normalization of ground returns for DTM bias**, and
**per-return local-incidence filtering of ground returns**, are **genuine gaps**: Roussel
did 1/cos for canopy, Tan built the slope+aspect+scan incidence for intensity, but no one
has combined them for bare earth.

---

## VERDICT — is the signed mechanism documented, or is our measurement a contribution?

**The signed "perpendicular-reaches-ground → lower; oblique-stops-on-veg → higher"
elevation-bias mechanism is documented only in *disassembled* form. Its parts are each
well supported; the assembled, slope-resolved, per-return-incidence claim is not stated
in the peer-reviewed record, making our signed measurement a genuine contribution.**

Supported parts:
- 1/cos θ near-ground path length controls reaching ground (Roussel 2018; Ni-Meister
  2001; MacArthur–Horn 1969) — as **probability/density**, on **scan angle**.
- Stopping on low vegetation biases ground **upward**, with magnitudes +4.5 to +48.8 cm
  (Hopkinson 2005; Ewald 2013; Töyrä 2003; *Spartina* literature) — on **flat** ground,
  by **cover class**.
- The sign is opposite to detector timewalk (Baltsavias 1999; Gardner 1992), so the two
  are cleanly distinguished.
- Slope degrades ground DTMs and thins ground returns; both mechanisms named together
  once (Mohd Salleh 2015; Lewis & Hancock 2007 for the on-slope "higher altitude"
  ground return).

Not documented (our contribution):
1. The bias tied to beam-to-**local-slope** incidence (not off-nadir scan angle), via a
   **per-return, physically reconstructed** incidence (scan + slope + aspect).
2. A **signed, slope-resolved magnitude** for the *bare/forest floor* (~6 mm at 0–8° to
   ~35 mm at 16–24°) — no rendered source gives a signed per-slope floor bias of this
   kind; the closest (Hopkinson, Ewald) are flat and cover-classed.
3. Selecting **near-slope-perpendicular** returns (or 1/cos-normalizing) to *deepen and
   de-bias* the ground surface — no ground/DTM paper does this; the machinery exists only
   for intensity (Tan 2021) or canopy (Roussel 2018).
4. The **±17°-scanner-cannot-reach-perpendicular-on-steep-slopes** ceiling as an
   operational limit on repeat-survey ground comparability — geometrically trivial, not
   stated in the change-detection literature.

Present items 1–4 as our findings, not attributed. The pooled scan-angle nulls (Su &
Bork; Ahokas 2011) are not a refutation — they test scan angle on near-flat references,
whereas our covariate is per-return incidence to the local slope; but they should be
cited as the tension our per-return, slope-resolved treatment resolves.

---

## Verification status summary

- **VERIFIED (full text read):** Roussel et al. 2018; Hopkinson et al. 2005; Ewald 2013
  (thesis); Mohd Salleh et al. 2015; Ahokas et al. 2005 (workshop); Tan et al. 2021.
- **VERIFIED (authoritative metadata/DOI; body via concordant secondaries):** Su & Bork
  2006; Hodgson & Bresnahan 2004; Ni-Meister et al. 2001; MacArthur & Horn 1969;
  Baltsavias 1999; Gardner 1992; Liu et al. 2018; Korpela et al. 2012; Ahokas et al.
  2011; Estornell et al. 2011; Höfle & Pfeifer 2007; Soudarissanane et al. 2011.
- **VERIFIED only as quoted in another source (original not opened):** Töyrä et al. 2003
  (in Hopkinson 2005); Goodwin et al. 2007, Lewis & Hancock 2007 (in Mohd Salleh 2015 /
  Roussel 2018); the *Spartina*-marsh +10–45 cm figures (in Ewald 2013).
- **CORRECTED attributions (do not repeat the companion-file errors):** +4.5 cm = Ewald
  2013, not Bater & Coops 2009; +0.07/+0.15 m graminoid/scrub = Töyrä 2003, not
  Hopkinson's own; ε(15°)≈0.06 = Roussel's absorption parameter, not a ground
  normalization; Disney 2010 / Holmgren 2003–2004 are canopy, not ground/DTM.
- **No fabricated citations or magnitudes.** Unverified specifics are flagged inline.

## Full citation list

1. Ahokas, E., Yu, X., Oksanen, J., Hyyppä, J., Kaartinen, H. & Hyyppä, H. (2005).
   Optimization of the scanning angle for countrywide laser scanning. *ISPRS Laser
   Scanning 2005*, Enschede. VERIFIED (workshop, full text).
2. Ahokas, E., Hyyppä, J., Yu, X. & Holopainen, M. (2011). Transmittance of airborne
   laser scanning pulses for boreal forest elevation modeling. *Remote Sensing*
   3(7):1365–1379. DOI 10.3390/rs3071365. VERIFIED (partial-negative result).
3. Baltsavias, E.P. (1999). Airborne laser scanning: basic relations and formulas.
   *ISPRS J. Photogramm. Remote Sens.* 54(2–3):199–214.
   DOI 10.1016/S0924-2716(99)00015-5. VERIFIED (citation; timewalk, opposite sign).
4. Estornell, J., Ruiz, L.A., Velázquez-Martí, B. & Hermosilla, T. (2011). Analysis of
   the factors affecting LiDAR DTM accuracy in a steep shrub area. *Int. J. Digital
   Earth* 4(6):521–538. PARTIALLY VERIFIED.
5. Ewald, M.J. (2013). Where's the ground surface? Elevation bias in LIDAR-derived
   DEMs due to dense vegetation in Oregon tidal marshes. M.S. thesis, Oregon State
   University. VERIFIED (full text). *Thesis — not peer-reviewed.*
6. Gardner, C.S. (1992). Ranging performance of satellite laser altimeters. *IEEE TGRS*
   30(5):1061–1072. DOI 10.1109/36.175341. VERIFIED (opposite-sign timewalk).
7. Goodwin, N.R., Coops, N.C. & Culvenor, D.S. (2007). Development of a simulation model
   to predict LiDAR interception in forested environments. *RSE* 111(4):481–492 (venue
   from Roussel's citation — CONFIRM DOI before quoting). VERIFIED as cited only.
8. Hodgson, M.E. & Bresnahan, P. (2004). Accuracy of airborne lidar-derived elevation:
   empirical assessment and error budget. *PE&RS* 70(3):331–339.
   DOI 10.14358/PERS.70.3.331. VERIFIED (geometric tan-slope mechanism).
9. Höfle, B. & Pfeifer, N. (2007). Correction of laser scanning intensity data: data and
   model-driven approaches. *ISPRS J. Photogramm. Remote Sens.* 62(6):415–433. VERIFIED
   (intensity, not elevation).
10. Hopkinson, C., Chasmer, L.E., Zsigovics, G., Creed, I.F., Sitar, M., Treitz, P. &
    Maher, R.V. (2005). Errors in LiDAR ground elevation and wetland vegetation height
    estimates. *IAPRS* XXXVI-8/W2:108–113. VERIFIED (full text; upward bias +0.07/+0.15 m,
    flat terrain). *Companion journal article: Hopkinson et al. 2005, Can. J. Remote
    Sens. 31(2):191–206.*
11. Korpela, I., Hovi, A. & Morsdorf, F. (2012). Understory trees in airborne LiDAR data
    — selective mapping due to transmission losses and echo-triggering mechanisms. *RSE*
    119:92–104. DOI 10.1016/j.rse.2011.12.011. VERIFIED.
12. Lewis, P. & Hancock, S. (2007). LiDAR for vegetation applications. UCL (cited in Mohd
    Salleh 2015 for on-slope ground-return-at-higher-altitude geometry). VERIFIED as
    cited only; original not opened.
13. Liu, J., Skidmore, A.K., Jones, S., Wang, T., Heurich, M., Zhu, X. & Shi, Y. (2018).
    Large off-nadir scan angle of airborne LiDAR can severely affect the estimates of
    forest structure metrics. *ISPRS J. Photogramm. Remote Sens.* 136:13–25.
    DOI 10.1016/j.isprsjprs.2017.12.004. VERIFIED (CANOPY, despite DTM-sounding title).
14. MacArthur, R.H. & Horn, H.S. (1969). Foliage profile by vertical measurements.
    *Ecology* 50(5):802–804. VERIFIED (gap-probability basis).
15. Mohd Salleh, M.R., Ismail, Z. & Abdul Rahman, M.Z. (2015). Accuracy assessment of
    LiDAR-derived DTM with different slope and canopy cover in tropical forest region.
    *ISPRS Annals* II-2/W2:183–189. VERIFIED (full text; names both mechanisms).
16. Ni-Meister, W., Jupp, D.L.B. & Dubayah, R. (2001). Modeling lidar waveforms in
    heterogeneous and discrete canopies. *IEEE TGRS* 39(9):1943–1958. VERIFIED
    (P_gap = exp(−G·LAI/cos θ)).
17. Roussel, J.-R., Béland, M., Caspersen, J. & Achim, A. (2018). A mathematical
    framework to describe the effect of beam incidence angle on metrics derived from
    airborne LiDAR: the case of forest canopies approaching turbid medium behaviour.
    *RSE* 209:824–834. DOI 10.1016/j.rse.2017.12.006. VERIFIED (full text). *ε(15°)≈0.06
    is the absorption parameter, NOT a ground normalization.*
18. Soudarissanane, S., Lindenbergh, R., Menenti, M. & Teunissen, P. (2011). Scanning
    geometry: influencing factor on the quality of terrestrial laser scanning points.
    *ISPRS J. Photogramm. Remote Sens.* 66(4):389–399. VERIFIED (TLS range precision).
19. Su, J. & Bork, E. (2006). Influence of vegetation, slope, and lidar sampling angle on
    DEM accuracy. *PE&RS* 72(11):1265–1274. VERIFIED (cover-dependent sign; scan-angle
    null below ~15°).
20. Tan, K., Cheng, X. & Zhao, C. (2021). Airborne LiDAR intensity correction based on a
    new method for incidence angle correction for improving land-cover classification.
    *Remote Sensing* 13(3):511. DOI 10.3390/rs13030511. VERIFIED (full local-incidence
    geometry — but for INTENSITY, not elevation).
21. Töyrä, J., Pietroniro, A., Hopkinson, C. & Kalbfleisch, W. (2003). Assessment of
    airborne scanning laser altimetry in a deltaic wetland environment. *Can. J. Remote
    Sensing* 29(6):679–690. Ground overestimate +0.07 m (graminoid) to +0.15 m (scrub).
    VERIFIED only as quoted in Hopkinson 2005; original not opened.

### Not to be cited for this argument (canopy, not ground/DTM)
- Disney et al. (2010), *RSE* 114(7):1546–1560 — canopy-height scan-angle simulation;
  the "<15° to protect ground detection" framing is a paraphrase, not the paper.
- Holmgren, Nilsson & Olsson (2003), *Can. J. Remote Sens.* 29(5):623–632; Holmgren
  (2004), *Scand. J. For. Res.* 19(6):543–553 — canopy tree-height/closure; the "≤10°"
  is Nordic acquisition guideline, unverified verbatim from these papers.

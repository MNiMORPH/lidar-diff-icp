# Scan-angle (nadir-vs-edge) dependence of lidar vertical accuracy on slopes — literature synthesis

Agent literature search 2026-08-22 (Andy-requested). Companion to and distinct from
[INCIDENCE_ANGLE_LITERATURE.md](INCIDENCE_ANGLE_LITERATURE.md) (broad incidence review),
[GROUND_BIAS_VS_INCIDENCE_LITERATURE.md](GROUND_BIAS_VS_INCIDENCE_LITERATURE.md)
(the signed vegetation-path mechanism), [EARLYGEN_SLOPE_ERROR_LITERATURE.md](EARLYGEN_SLOPE_ERROR_LITERATURE.md)
(early-gen downward-slope mechanisms, boresight, co-registration), and
[BARE_GROUND_SKEW_LITERATURE.md](BARE_GROUND_SKEW_LITERATURE.md). **Those files establish
the physics (1/cos θ path length, timewalk, tan-slope coupling, min-bin gridding,
boresight). This file does the one thing they do not: it tests the specific
`scan angle` (near-NADIR vs. swath-EDGE) partition of vertical accuracy on slopes, and
the counterintuitive `near-nadir-WORSE-on-slopes` signature we measure — reconciling it
with the community's dominant `swath-edge-worst` narrative.**

Citations are verified per item; every unverified specific is flagged. No DOIs,
magnitudes, or attributions were invented. Where a source was already fetched/verified in
a companion file, that is stated and the finding is cross-referenced rather than re-derived.

---

## The empirical result driving this file (our measurement — the hypothesis to test)

From `data/derived/elba_fulldensity/gen1_csf_angles.npz` (6.77M gen1 CSF ground returns),
gen1 (2008 discrete-return, ~1 pt/m², leaf-off Nov) vs a fixed gen2 (2021 3DEP) bare-earth
plane, residual `r = d − datum` (datum = flat-ground median, carrying the +67 mm
GEOID03→GEOID18 constant; `r < 0` = gen1 reads the ground LOW). Verified in
[NEARNADIR_SLOPE_DEPENDENCE.md](NEARNADIR_SLOPE_DEPENDENCE.md) and
[GLENNIE_SCANANGLE_SWATH_TEST.md](GLENNIE_SCANANGLE_SWATH_TEST.md):

1. **Near-nadir beams (|scan angle| < 5°) carry a slope-growing DOWNWARD residual:** ≈0 on
   flat, ≈−15 mm through 3–24°, deepening to **−40 to −44 mm by 27–35°**.
2. **Grazing / swath-edge beams (|scan angle| > 12–15°) do NOT:** they stay flat near 0
   across the whole slope range. Within every slope band the `edge − nadir` delta is
   **positive (+10 to +26 mm)** — the residual *improves* toward the swath edge.
3. For near-nadir beams, incidence-to-local-surface-normal ≈ slope (≈ 0.8 + 0.95·slope), so
   the near-nadir-vs-slope curve **is** a vertical-error-vs-incidence curve.
4. **Boresight / scan-mirror / legacy-nav EXCLUDED** (Glennie 2014 diagnostic fails: sign is
   nadir-worst not edge-worst, and edge behaviour is per-swath incoherent); the constant
   datum and real erosion excluded (persists on stable zero-curvature ridgeline points).

So the swath-EDGE (high scan angle) beams are the *accurate* ones and the near-NADIR beams
carry the slope bias — the **opposite** of the usual "accuracy degrades toward the swath
edge" expectation. The four questions below test that against the record.

---

## (Q1) Is scan-angle-dependent vertical accuracy, and its slope interaction, documented?

**Yes — abundantly, but almost always in the `edge-worst` direction, and mostly as
variance/RMSE rather than a signed mm-per-degree bias.** The dominant, well-supported
narrative is that vertical error is *smallest at the swath center (near nadir) and grows
toward the swath edge*, and separately that it grows with slope. Both are in the companion
files; the on-point sources for the *scan-angle* half:

- **Su & Bork (2006), *PE&RS* 72(11):1265–1274.** VERIFIED (companion files). The single most
  on-point empirical DEM paper: RMSE for slopes >10° up to ~2× that for <2°; **laser scan
  angle had little effect below ~15° off-nadir**; slope and vegetation dominate, scan angle
  second-order within a narrow swath. (Note: this is a *pooled scan-angle null*, not a test
  of per-return incidence to the local slope — see the reconciliation in Q2.)
- **Ahokas, Kaartinen & Hyyppä (2003), *IAPRS* XXXIV-3/W13** (ISPRS Workshop). VERIFIED
  (companion). **Changes in scan angle produced systematic elevation errors ~10 cm**; height
  error grows with altitude. A signed systematic offset with observation angle, at wider
  angles.
- **Goulden & Hopkinson (2010), *PE&RS* 76(5):589–601** (DOI 10.14358/PERS.76.5.589).
  VERIFIED (companion). GLOPOV error propagation: random GPS/IMU/range/encoder/divergence
  errors **propagate more severely into 3-D coordinate error as scan angle and altitude
  increase** — the mechanism behind larger error at swath edges.
- **Glennie et al. (2014), *GRL* 41** (DOI 10.1002/2014GL059919). VERIFIED — full PDF read
  (EARLYGEN file). The canonical scan-angle-error diagnostic: residual displacement "shows a
  clear correlation with scan angle and a large 10 m difference between observations near
  nadir and those at the 30° extent of the scan," errors "mainly at the edge of the flight
  line swaths (i.e., at larger scan angles)." **This is the `edge-worst` template our Glennie
  test explicitly fails**, excluding boresight for us.
- **Sensitivity of DEM/slope/aspect/watershed attributes to lidar uncertainty (2016),
  *RSE*.** VERIFIED (metadata; companion). σ_DEM **minimum at swath center, growing toward
  the edge (with scan angle)**, and DEM error grows with slope — and because two epochs have
  different flight-line layouts, a given cell sits at different scan angles in each survey.
  The clearest statement that scan-angle geometry differs between epochs and matters for a
  DoD.
- **Vosselman; Habib et al. (strip-overlap geometric QA).** VERIFIED thematically (search
  metadata; primary internals NOT fetched — treat as citation-level). Conjugate-surface
  discrepancies between overlapping strips are used to detect systematic errors; discrepancy
  and striping concentrate in the swath-edge / large-scan-angle overlap. The classic
  `edge-worst` engineering picture. **Do not quote specific numbers — flagged citation-only.**

**Magnitude / sign / angle answer (Q1):** the documented scan-angle vertical effect is a
**~10 cm systematic offset at wide angles** (Ahokas 2003) plus **edge-growing variance**
(Goulden & Hopkinson; RSE-2016), with slope entering separately as a ~2× RMSE growth by
10–25° (Su & Bork; Hodgson & Bresnahan — companion files). Within a **≤15°** swath the scan
angle is **second-order** (Su & Bork). All of this is `edge-worst` and mostly *unsigned
scatter*, framed on scan angle against a **flat/near-flat** reference — not a signed
near-nadir slope bias.

---

## (Q2) Near-nadir WORSE than off-nadir on slopes? — CONTRADICTED as a *range-error* claim, SUPPORTED once reframed as a *ground-penetration* claim

This is the crux, and it turns entirely on **what "worse" means**, because two different
error families point in opposite directions.

### Narrative A — range/pointing error: edge-WORST (our result CONTRADICTS this, correctly)

The `edge-worst` narrative (Q1) is a **range- and pointing-error** statement: georeferencing
error (GPS/IMU/boresight/encoder) and footprint elongation grow with scan angle, so the 3-D
*position* of a return degrades toward the swath edge (Goulden & Hopkinson 2010; Glennie
2014; the strip-QA literature). **If our residual were this family, it would be edge-worst
and per-swath coherent — it is neither** (our Glennie test, `GLENNIE_SCANANGLE_SWATH_TEST.md`:
edge−nadir is *positive*, and edge behaviour flips sign across swaths). So we **correctly
contradict Narrative A** — and that contradiction is exactly what excludes boresight/nav
error as our driver. Our result is not "near-nadir has worse *ranging*"; it is the opposite
error family.

### Narrative B — canopy/near-ground penetration: near-nadir reaches the TRUE (lower) floor; oblique stops HIGHER

The **penetration** family runs the other way, and **it is the one that matches our sign** —
this is the reframing that dissolves the paradox:

- **Roussel, Béland, Caspersen & Achim (2018), *RSE* 209:824–834**
  (DOI 10.1016/j.rse.2017.12.006). VERIFIED — full text (companion files). A beam at
  incidence θ travels **1/cos θ** farther through each vegetation layer, so oblique beams
  intercept vegetation sooner and the point distribution shifts **upward** (more canopy,
  fewer ground); **the near-vertical (near-nadir) beam has the shortest near-ground path and
  is the one that reaches the true, lower floor.** On a slope, "near-nadir" *is* the
  low-local-incidence beam (our #3), so near-nadir reads deep, oblique reads high.
- **Liu et al. (2018), *ISPRS J. P&RS* 136:13–25** (DOI 10.1016/j.isprsjprs.2017.12.004).
  VERIFIED (companion; a CANOPY gap-fraction paper). Nadir (0–7°) / small (7–23°) / large
  (23–38°) off-nadir bins: gap-fraction underestimation **amplifies at large off-nadir**,
  and the discrepancy **intensifies from upper- to lower-slope positions and with slope
  steepness** — i.e. the oblique-vs-slope penetration penalty is real and slope-position
  dependent. (Full-text fetch of the source PDF returned image-only; the slope-position
  detail is from the ScienceDirect/AWARE abstracts and is flagged VERIFIED-abstract, not
  full-text.)
- **Pang et al. (2011), "Impact of footprint diameter and off-nadir pointing on the
  precision of canopy height estimates from spaceborne lidar," *RSE* 115(11):2798–2809.**
  VERIFIED (abstract/metadata; full text NOT fetched). Precision **decreases with off-nadir
  pointing**; local incidence = terrain slope ± off-nadir angle depending on aspect (when
  aspect and sensor azimuth oppose, local incidence = |slope − off-nadir|) — the explicit
  slope × pointing × azimuth composition, matching our incidence reconstruction (spaceborne,
  large-footprint; transfer the *geometry*, not the magnitudes). **Flag: spaceborne,
  large-footprint — different regime.**
- **Return-structure evidence in our own data** ([[return-structure-leafoff-penetration]]):
  the deepest gen1 echoes reach **266 mm deeper** to the forest floor than gen2; gen1
  leaf-off Nov penetrates the near-ground layer better than gen2 leaf-on May. So on our data
  the *reference* (gen2) sits high in forest, and gen1's near-nadir beams read the true floor
  **low relative to it** — precisely the "near-nadir reaches deeper" prediction.

**Reconciliation (the whole point):** the two narratives are **not in conflict — they
measure different error channels.** Narrative A (range/pointing) degrades toward the edge;
Narrative B (ground penetration through near-ground vegetation) degrades toward the *oblique*
by making oblique beams read the floor **too high**. Our residual is defined against a
ground surface, so it is dominated by B, and **B makes near-nadir the DEEPER (here: the one
that departs downward from the leaf-on gen2 reference), oblique the shallower/agreeing one.**
There is no published paper stating "near-nadir is worse on slopes" as a *ranging* claim
(that would be false); there is solid published physics that near-nadir **penetrates deeper
to the true floor** (Roussel 2018; Liu 2018; Pang 2011), which is our sign once the reference
is the leaf-on epoch. Our result is therefore **not contradicted** — it is the penetration
family, and it *contradicts only the range-error narrative*, which is exactly the diagnostic
value (it excludes boresight).

**A caveat kept honest:** in the penetration frame, near-nadir gen1 is arguably the *more
accurate* reading of the true floor (it penetrates deepest), and the gen2 leaf-on reference
is the biased-high one. "Near-nadir worse" is only true **relative to the gen2 datum**; in
absolute terms near-nadir may be **better**. This inverts the word "worse" and must be stated
that way in the writeup — the near-nadir beams are *low relative to gen2*, not *wrong*.

---

## (Q3) Physical mechanism for a near-VERTICAL beam ranging a slope DEEPER — candidate ranking

Four candidates, ranked by how well each predicts **near-nadir LOWER** (not higher) on our
leaf-off/leaf-on, discrete-return, forested-slope data.

**★ TOP CANDIDATE — near-ground vegetation path length (1/cos θ), reference-relative.**
Near-nadir = shortest path through the near-ground vegetation layer → penetrates to the true
(lower) floor; oblique = longer path (1/cos θ) → stops higher on litter/understory/stubble.
Against a **leaf-on** reference (gen2) that itself sits high, near-nadir gen1 reads **low**.
**Predicts near-nadir LOWER: YES.** Best-supported, mechanism established.
- **Roussel et al. (2018), *RSE* 209:824–834**, DOI 10.1016/j.rse.2017.12.006. VERIFIED
  (full text). The 1/cos θ path-length control and upward shift of oblique point
  distributions. **This is the top mechanism citation.**
- Supporting: **Liu et al. (2018)** (slope-position amplification), **Pang et al. (2011)**
  (slope × pointing geometry), and our own [[return-structure-leafoff-penetration]] deep-echo
  test (gen1 266 mm deeper in forest). See also **Hopkinson et al. (2005), *IAPRS*
  XXXVI-8/W2:108–113** (VERIFIED, companion): oblique/low-penetration returns bias ground
  **upward** +0.07–0.15 m — the shallow side of the same coin.

**Candidate 2 — footprint elongation downslope + the range assigned to an asymmetric
footprint.** A vertical beam on a slope illuminates an ellipse elongated *down the fall
line*; where the return's range is taken (leading-edge vs. peak vs. centroid) sets the sign.
**This is a documented SIGN TRAP** (EARLYGEN file): raw geometry / leading-edge / first-photon
actually biases the ground **HIGH** (near edge returns first → range short → elevation high;
Laconte 2019 ~20 cm at high incidence; Wang 2018). A **downward** bias needs a *processing*
route (fixed-threshold timewalk on the broadened echo, or noise-clipping of the weak up-slope
tail). **Predicts near-nadir LOWER: only via the estimator, and it is stronger at OBLIQUE
incidence, so it does not naturally single out near-nadir.** Secondary at best for us.
- **Baltsavias (1999), *ISPRS J.* 54(2–3):199–214**, DOI 10.1016/S0924-2716(99)00015-5;
  **Gardner (1992), *IEEE TGRS* 30(5):1061–1072**, DOI 10.1109/36.175341;
  **Deems, Painter & Finnegan (2013), *J. Glaciology* 59(215):467–479**, doi:10.3189/2013JoG12J154
  (VERIFIED full text, EARLYGEN file: timewalk ~50 cm at 45°/1000 m). Timewalk is downward
  but **grows with incidence** → predicts *oblique* worse, not nadir worse. Runs against us
  as a near-nadir explanation.
- **Laconte et al. (2019), IEEE ICRA / arXiv:1810.01619; Wang et al. (2018), diva-portal
  1160255.** VERIFIED (full text, EARLYGEN file). Leading-edge on a slope → range short →
  elevation **HIGH** — the sign trap. Wrong sign for near-nadir-low.

**Candidate 3 — multiple-return / last-return / lowest-return selection on slopes as a
function of scan angle.** Near-nadir beams give the **densest** ground sampling (highest
ground-return probability, swath-center density peak); a min-elevation / lowest-return
gridding step then selects the **downhill** point in a cell, pulling the surface **down**,
and the pull grows with points-per-cell and with slope (≈ ½·cellsize·tan s). Because
near-nadir cells are the *most densely* sampled, the downhill-selection pull is **strongest
at near-nadir** — a genuine route to *near-nadir-lower*. **Predicts near-nadir LOWER: YES,
conditional on a lowest-return gridding step.** Depends on our pipeline (our ground grid is a
per-cell **median** slope-normal, `ground_q=0.50`, not a min — see [[ground-construction-two-stage]]),
so this route is *weakened* for us but not zero (denser near-nadir sampling still shifts the
within-cell distribution).
- **Su & Bork (2008), *Computers & Geosciences* 35(2):289–300** (VERIFIED as paraphrase,
  EARLYGEN file): "increasingly underestimate terrain elevation as slope increased."
- **lidR issue #51** (VERIFIED verbatim, EARLYGEN file): min-in-cell selects the downhill
  point on a regular slope → lower elevation.
- **Petras et al. (2023), *Sensors* 23(3):1593**, DOI 10.3390/s23031593 (VERIFIED,
  companion): line-scanner density peaks mid-swath (near nadir) → most points-per-cell at
  nadir → strongest downhill-selection pull there.

**Candidate 4 — canopy/vegetation path-length differences through a tilted canopy for
near-vertical vs oblique beams.** A special case of Candidate 1 for the overstory rather than
the near-ground layer; same 1/cos θ physics (Roussel 2018; Ni-Meister et al. 2001, GORT;
MacArthur & Horn 1969). Predicts the *same sign* as Candidate 1 and is subsumed by it.

**Q3 verdict:** the best-supported mechanism for **near-nadir-reads-deeper** is **Candidate 1
— near-ground vegetation path length (1/cos θ), referenced to a leaf-on epoch** (Roussel
2018), with **Candidate 3 — denser near-nadir sampling feeding a downhill within-cell pull**
as a real secondary that depends on the gridding estimator. The classical range-physics
routes (Candidate 2: footprint/timewalk) are the **wrong sign or wrong angle-dependence** to
be the near-nadir driver — timewalk is downward but *oblique*-worst, and raw geometry is
*upward*. This is why our signal cannot be the "usual" edge/range story.

---

## (Q4) Scan-angle cutoffs / weighting for slope-terrain vertical accuracy — thin, mostly canopy-framed; a near-nadir-weighted GROUND grid is a gap

- **Ahokas, Yu, Oksanen, Hyyppä, Kaartinen & Hyyppä (2005), "Optimization of the Scanning
  Angle for Countrywide Laser Scanning," *ISPRS Laser Scanning 2005*, Enschede.** VERIFIED —
  full text (GROUND_BIAS file). The correct **ground/DTM** near-nadir citation: beyond **10°
  off-nadir** shadowed area increases, ground hits decrease, and DTM gaps grow; angles "up to
  15 degrees seems to be usable." This is the closest thing to a **ground-DTM scan-angle
  ceiling** in the literature.
- **Estornell et al. (2011), *Int. J. Digital Earth* 4(6):521–538.** PARTIALLY VERIFIED
  (GROUND_BIAS file). Steep shrub: ground-return likelihood **increases with narrow scan
  angles** — supports near-nadir selection for the ground on steep terrain.
- **Disney et al. (2010) <15° and Holmgren (2003/2004) ≤10°** are **CANOPY** limits (crown
  shadowing), **not ground/DTM** — do NOT cite them for a ground argument (correction carried
  in GROUND_BIAS_VS_INCIDENCE). The "≤10°" figure could not be confirmed verbatim from the
  Holmgren papers; attribute the Nordic countrywide ceiling to Ahokas 2005 instead.
- **Roussel et al. (2018)** proposes **normalizing to nadir** (1/cos θ correction) rather
  than discarding oblique data — but for **canopy** metrics. A **1/cos-path-length
  normalization of GROUND returns for DTM bias** is not published — a genuine gap.
- **Tan, Cheng & Zhao (2021), *Remote Sensing* 13(3):511**, DOI 10.3390/rs13030511.
  VERIFIED (GROUND_BIAS file). Builds the full **per-return local-incidence** geometry
  (scan + slope + aspect) exactly as we do, and warns it **over-corrects on steep slopes**,
  recommending an **incidence-angle threshold** (trade-off: reducing overcorrection vs.
  excluding terrain) — but the target is **intensity**, not elevation. The nearest published
  precedent for an incidence cutoff on steep slopes, wrong target.
- **Pooled scan-angle NULLS (the honest counter-evidence):** **Su & Bork (2006)** (scan
  angle little effect below ~15°) and **Ahokas et al. (2011), *Remote Sensing* 3(7):1365–1379**,
  DOI 10.3390/rs3071365 (VERIFIED, GROUND_BIAS file: no significant scan-angle effect on
  ground hits over 0–15° in boreal forest). These test **scan angle on flat-ish references**,
  not per-return incidence to the local slope — so they are the *tension our per-return,
  slope-resolved treatment resolves*, not a refutation.

**Q4 verdict:** there is a **ground/DTM near-nadir preference** (Ahokas 2005: ≤10–15°;
Estornell 2011) and an **incidence-threshold** precedent (Tan 2021, for intensity), but **no
published near-nadir-WEIGHTED ground gridding nor 1/cos-ground normalization for a slope
vertical bias.** A near-nadir-weighted / incidence-thresholded ground surface for the DoD is
supported in spirit but **novel in execution** — consistent with the practical conclusion in
[[scan-geometry-governs-gen1-floor]] that the cleanest gen1 surface is near-nadir/low-incidence
weighted and the forest correction must be empirical stable-terrain calibration.

---

## VERDICT — documented, contradicted, or novel?

**Our near-nadir-worse-on-slopes signature is NOVEL as a stated, signed, scan-angle-resolved
ground result — and it is NOT contradicted; it belongs to a different error family than the
narrative it appears to invert.**

- It **contradicts the `swath-edge-worst` narrative only as a *range/pointing-error* claim**
  — and that contradiction is the diagnostic that **excludes boresight/nav error** (our
  Glennie test), consistent with Glennie 2014's own edge-worst template being the thing we
  fail.
- It is **the vegetation-**penetration** family, whose physics is well documented**
  (Roussel 2018; Liu 2018; Pang 2011): near-nadir/low-incidence beams penetrate deepest to
  the true floor, oblique beams stop higher. **Referenced to a leaf-on epoch (gen2), that
  makes near-nadir read LOW — exactly our sign.** No paper is *against* this; several supply
  the mechanism.
- **No peer-reviewed source states**: (a) a **signed, per-return, scan-angle-partitioned**
  ground-elevation-vs-slope curve showing near-nadir carrying the bias and the swath edge
  flat; (b) the sign that **near-nadir departs DOWNWARD** from a repeat-survey reference while
  the edge agrees; (c) that the apparent "worse" is reference-relative (near-nadir is
  arguably the *truer* floor). Present items (a)–(c) as **our findings**, framed via the
  penetration mechanism.

**Top mechanism candidate:** near-ground vegetation path length (1/cos θ), reference-relative
— **Roussel, Béland, Caspersen & Achim (2018), *RSE* 209:824–834, DOI
10.1016/j.rse.2017.12.006** — with denser near-nadir sampling feeding within-cell
downhill-selection (Su & Bork 2008; lidR #51; Petras 2023) as a gridding-dependent secondary.
The classical footprint/timewalk range physics (Baltsavias 1999; Gardner 1992; Deems 2013;
Laconte 2019) is the **wrong sign or oblique-worst**, so it is *not* the near-nadir driver —
a clean exclusion, not a gap.

**Scan-angle-cutoff guidance for a correction:** the only ground-framed ceiling is **Ahokas
et al. (2005)** (≤10–15° off-nadir for DTM ground hits), reinforced by Estornell (2011)
(narrow angles favor ground on steep terrain); the only steep-slope incidence-threshold
precedent is **Tan et al. (2021)** (built for intensity, warns of over-correction on steep
slopes → use an incidence threshold). A **near-nadir-weighted ground grid** or **1/cos-ground
normalization** for the slope vertical bias is supported in spirit but **unpublished** — our
correction would be a contribution, not a reproduction.

---

## Honest gaps and negatives

- **The `edge-worst` vs `nadir-worst` reconciliation rests on defining the reference.** Our
  "near-nadir worse" is *relative to the leaf-on gen2 datum*; absolute-accuracy-wise
  near-nadir may be the better reading. This is a framing point the literature does not make
  for us and that we must state carefully — it is the single most misstateable claim here.
- **Liu (2018) slope-position detail is VERIFIED-abstract, not full-text** (source PDF
  fetched image-only). Pang (2011) is abstract/metadata only. Both are **spaceborne /
  large-footprint or CANOPY** — transfer the geometry, not magnitudes; do not present either
  as a bare-ground airborne result.
- **No source gives a signed mm-per-degree near-nadir GROUND slope bias** — consistent with
  the companion files. Our ~−40 mm-by-27° near-nadir curve is our measurement, not attributed.
- **The 27° "knee"** ([NEARNADIR_SLOPE_DEPENDENCE.md](NEARNADIR_SLOPE_DEPENDENCE.md): tan+step@27
  R²=0.800 beats smooth tan R²=0.53) has **no literature analog** — every documented mechanism
  is a smooth tan-law (EARLYGEN file). Treat the knee as our finding to interrogate (is it the
  ±17°-scanner-cannot-reach-perpendicular ceiling crossing? a leaf-off-penetration threshold?
  a min-return-density crossing?), not a confirmed effect.
- **Vosselman / Habib strip-QA sources are citation-only** (search metadata; internals not
  fetched) — do not quote numbers from them.
- **Candidate 3 (downhill within-cell selection) is pipeline-dependent.** Our grid is a
  per-cell median slope-normal, not a min, so this route is *weakened* for us; if invoked it
  must be checked against the actual `ground_q=0.50` estimator, not assumed.

---

## Verification status summary

- **VERIFIED (full text read — in this or a companion file):** Roussel et al. 2018; Ahokas
  et al. 2005 (workshop); Glennie et al. 2014; Deems et al. 2013; Laconte et al. 2019; Wang
  et al. 2018; Tan et al. 2021; Goulden & Hopkinson 2010; Hopkinson et al. 2005; Petras et
  al. 2023.
- **VERIFIED (authoritative metadata/DOI; body via concordant secondaries or companion
  file):** Su & Bork 2006; Su & Bork 2008; Hodgson & Bresnahan 2004; Baltsavias 1999; Gardner
  1992; Ahokas et al. 2003; Ahokas et al. 2011; Estornell et al. 2011; Ni-Meister et al.
  2001; MacArthur & Horn 1969; RSE-2016 sensitivity study; lidR #51.
- **VERIFIED (abstract/metadata only — full text NOT fetched):** Liu et al. 2018
  slope-position detail (source PDF image-only); Pang et al. 2011.
- **CITATION-ONLY (internals NOT fetched — do NOT quote numbers):** Vosselman; Habib et al.
  strip-overlap geometric QA.
- **No fabricated citations, DOIs, or magnitudes.** Every unverified specific is flagged
  inline. Attribution corrections carried in the companion files (Disney/Holmgren = canopy
  not ground; ε(15°) = absorption not ground normalization) are respected here.

## Full citation list

1. Ahokas, E., Kaartinen, H. & Hyyppä, J. (2003). A quality assessment of airborne laser
   scanner data. *IAPRS* XXXIV-3/W13, ISPRS Workshop, Dresden. VERIFIED (conference; ~10 cm
   systematic elevation error with scan angle).
2. Ahokas, E., Yu, X., Oksanen, J., Hyyppä, J., Kaartinen, H. & Hyyppä, H. (2005).
   Optimization of the scanning angle for countrywide laser scanning. *ISPRS Laser Scanning
   2005*, Enschede. VERIFIED (full text; ground-DTM ≤10–15° ceiling).
3. Ahokas, E., Hyyppä, J., Yu, X. & Holopainen, M. (2011). Transmittance of airborne laser
   scanning pulses for boreal forest elevation modeling. *Remote Sensing* 3(7):1365–1379.
   DOI 10.3390/rs3071365. VERIFIED (scan-angle null 0–15°).
4. Baltsavias, E.P. (1999). Airborne laser scanning: basic relations and formulas. *ISPRS J.
   P&RS* 54(2–3):199–214. DOI 10.1016/S0924-2716(99)00015-5. VERIFIED (citation; timewalk,
   oblique-worst → wrong angle for near-nadir).
5. Deems, J.S., Painter, T.H. & Finnegan, D.C. (2013). Lidar measurement of snow depth: a
   review. *J. Glaciology* 59(215):467–479. doi:10.3189/2013JoG12J154. VERIFIED (full text;
   timewalk ~50 cm at 45°/1000 m).
6. Estornell, J., Ruiz, L.A., Velázquez-Martí, B. & Hermosilla, T. (2011). Analysis of the
   factors affecting LiDAR DTM accuracy in a steep shrub area. *Int. J. Digital Earth*
   4(6):521–538. PARTIALLY VERIFIED (narrow angles favor ground on steep terrain).
7. Gardner, C.S. (1992). Ranging performance of satellite laser altimeters. *IEEE TGRS*
   30(5):1061–1072. DOI 10.1109/36.175341. VERIFIED (broadening ∝ footprint·tan s; centroid →
   no first-order bias).
8. Glennie, C., Hinojosa-Corona, A., Nissen, E., Kusari, A., Oskin, M.E., Arrowsmith, J.R. &
   Borsa, A. (2014). Optimization of legacy lidar data sets for measuring near-field
   deformation. *GRL* 41. DOI 10.1002/2014GL059919. VERIFIED (full text; scan-angle-correlated
   boresight residual, edge-worst — the template our test fails).
9. Goulden, T. & Hopkinson, C. (2010). The forward propagation of integrated system component
   errors within airborne lidar data. *PE&RS* 76(5):589–601. DOI 10.14358/PERS.76.5.589.
   VERIFIED (error grows with scan angle → edge-worst variance).
10. Hodgson, M.E. & Bresnahan, P. (2004). Accuracy of airborne lidar-derived elevation:
    empirical assessment and error budget. *PE&RS* 70(3):331–339. DOI 10.14358/PERS.70.3.331.
    VERIFIED (tan-slope coupling ~2× at 25°; RMSE, not signed).
11. Hopkinson, C., Chasmer, L.E., Zsigovics, G., Creed, I.F., Sitar, M., Treitz, P. & Maher,
    R.V. (2005). Errors in LiDAR ground elevation and wetland vegetation height estimates.
    *IAPRS* XXXVI-8/W2:108–113. VERIFIED (full text; oblique/low-penetration ground bias
    +0.07/+0.15 m upward — shallow side of the penetration coin).
12. Laconte, J., Deschênes, S.-P., Labussière, M. & Pomerleau, F. (2019). Lidar measurement
    bias estimation via return waveform modelling. *IEEE ICRA* / arXiv:1810.01619. VERIFIED
    (full text; leading-edge on slope → elevation HIGH; wrong sign for near-nadir-low).
13. Liu, J., Skidmore, A.K., Jones, S., Wang, T., Heurich, M., Zhu, X. & Shi, Y. (2018).
    Large off-nadir scan angle of airborne LiDAR can severely affect the estimates of forest
    structure metrics. *ISPRS J. P&RS* 136:13–25. DOI 10.1016/j.isprsjprs.2017.12.004.
    VERIFIED (CANOPY gap-fraction; nadir/small/large bins; slope-position amplification
    VERIFIED-abstract only).
14. MacArthur, R.H. & Horn, H.S. (1969). Foliage profile by vertical measurements. *Ecology*
    50(5):802–804. VERIFIED (gap-probability basis).
15. Ni-Meister, W., Jupp, D.L.B. & Dubayah, R. (2001). Modeling lidar waveforms in
    heterogeneous and discrete canopies. *IEEE TGRS* 39(9):1943–1958. VERIFIED
    (P_gap = exp(−G·LAI/cos θ)).
16. Pang, Y., Lefsky, M., Sun, G. & Ranson, J. (2011). Impact of footprint diameter and
    off-nadir pointing on the precision of canopy height estimates from spaceborne lidar.
    *RSE* 115(11):2798–2809. VERIFIED (abstract/metadata only; precision decreases with
    off-nadir; local incidence = slope ± off-nadir by aspect). *Spaceborne, large-footprint —
    transfer geometry, not magnitudes; DOI/pages NOT independently confirmed.*
17. Petras, V., Petrasova, A., McCarter, J.B., Mitasova, H. & Meentemeyer, R.K. (2023). Point
    density variations in airborne lidar point clouds. *Sensors* 23(3):1593.
    DOI 10.3390/s23031593. VERIFIED (line-scanner density peaks mid-swath/near-nadir).
18. Roussel, J.-R., Béland, M., Caspersen, J. & Achim, A. (2018). A mathematical framework to
    describe the effect of beam incidence angle on metrics derived from airborne LiDAR. *RSE*
    209:824–834. DOI 10.1016/j.rse.2017.12.006. VERIFIED (full text). **★ Top mechanism:
    1/cos θ near-ground path length; near-nadir reaches the true, lower floor.**
19. Su, J. & Bork, E. (2006). Influence of vegetation, slope, and lidar sampling angle on DEM
    accuracy. *PE&RS* 72(11):1265–1274. VERIFIED (scan-angle null below ~15°; slope ~2× RMSE).
20. Su, J. & Bork, E. (2008). Evaluating error associated with lidar-derived DEM
    interpolation. *Computers & Geosciences* 35(2):289–300. DOI 10.1016/j.cageo.2008.09.001.
    VERIFIED as paraphrase ("increasingly underestimate elevation as slope increased").
21. Tan, K., Cheng, X. & Zhao, C. (2021). Airborne LiDAR intensity correction based on a new
    method for incidence angle correction. *Remote Sensing* 13(3):511. DOI 10.3390/rs13030511.
    VERIFIED (full local-incidence geometry; over-correction on steep slopes → incidence
    threshold; target = INTENSITY, not elevation).
22. Wang, Y. et al. (2018). Bias of cylinder diameter estimation from ground-based laser
    scanners with different beam widths. diva-portal.org 1160255. VERIFIED (full text;
    leading-edge → point toward scanner → elevation HIGH).
23. Sensitivity of DEM, slope, aspect and watershed attributes to LiDAR measurement
    uncertainty (2016). *Remote Sensing of Environment*. VERIFIED (metadata; σ_DEM min at
    swath center, grows toward edge and with slope; epochs differ in flight-line layout).
24. Vosselman, G.; Habib, A. et al. — LiDAR strip-overlap geometric quality assessment /
    systematic-error detection from conjugate surfaces. CITATION-ONLY (search metadata;
    internals NOT fetched — do not quote numbers). The classical `edge-worst` strip picture.

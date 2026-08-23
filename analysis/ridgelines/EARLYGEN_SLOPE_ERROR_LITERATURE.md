# Early-generation lidar reading the ground TOO LOW on slopes — literature synthesis

Agent literature search 2026-08-22 (Andy-requested). Companion to
[GROUND_BIAS_VS_INCIDENCE_LITERATURE.md](GROUND_BIAS_VS_INCIDENCE_LITERATURE.md),
[INCIDENCE_ANGLE_LITERATURE.md](INCIDENCE_ANGLE_LITERATURE.md), and
[BARE_GROUND_SKEW_LITERATURE.md](BARE_GROUND_SKEW_LITERATURE.md). Those files treat the
**opposite-sign** vegetation-path mechanism (oblique beams stop on low near-ground
vegetation → ground reads **HIGH**). This file does the mirror-image job: it tests
whether the record documents an **early-generation, slope-dependent, DOWNWARD** ground
bias — the surface read **too low** on slopes — with the hypothesized ~25–30° onset and
tens-of-mm magnitude. It does **not** re-tread the vegetation-path material.

Citations verified per item against fetched full text or authoritative metadata; every
unverified specific is flagged. No DOIs, magnitudes, or attributions were invented.

## Our observation to explain (NOT a fact to confirm)

On stable, flow-free ridgeline reference cells, 2008 MN discrete-return lidar (gen1,
~1 pt/m², early GPS/IMU era) and 2021 3DEP (gen2, denser, modern) agree to ~+5 mm up to
~24° slope, then gen2−gen1 **jumps to +30–48 mm above ~27°** (positive = gen2 high **or
gen1 low**). Working hypothesis: the older survey reads the ground too low on steep
slopes, switching on ~25–30°. The literature's verdict on that specific claim follows.

---

## (1) Detector time-walk on the slope-broadened echo — DOWNWARD, DOCUMENTED, but large-slope/large-footprint

The one clean, published, **downward** slope mechanism. A footprint spread over an
inclined surface temporally broadens the return; on a broadened, lower-amplitude echo a
fixed **intensity-threshold / leading-edge** detector crosses its threshold **later**, so
the **recorded range increases** and the elevation is biased **low**. This is the
"time-walk" of Baltsavias.

- **Baltsavias (1999), *ISPRS J. P&RS* 54(2–3):199–214** (DOI 10.1016/S0924-2716(99)00015-5).
  VERIFIED as citation; body paywalled — quoted **verbatim through** the open-access
  review below, which I read directly.
- **Deems, Painter & Finnegan (2013), *Journal of Glaciology* 59(215):467–479**
  (doi:10.3189/2013JoG12J154). VERIFIED — primary PDF read; quote extracted directly.
  Their §3.3 / Fig. 5b, verbatim: *"This effect will spread the time distribution of the
  returned pulse, increasing the 'rise time' for the return to reach the intensity
  threshold for return signal registration and thus **increasing the recorded range
  distance**. For a 45° slope with a flight height of 1000 m, this 'time-walk' effect can
  induce a vertical error **close to 50 cm** (Baltsavias, 1999a). Smaller beam divergence
  angles produce smaller ground footprints, minimizing this effect."* Mitigation: orient
  flight lines to minimize oblique-incident shots on steep slopes. Fig. 5 caption cleanly
  **separates the two slope mechanisms**: (a) vertical error from horizontal error after
  Hodgson & Bresnahan 2004; (b) time-walk after Baltsavias 1999a.
- **Gardner (1992), *IEEE TGRS* 30(5):1061–1072** (DOI 10.1109/36.175341). VERIFIED —
  NASA full-text scan read. Supplies the **broadening physics**: mean-square pulse width
  (Eq. 18) has a slope term ∝ *z·θ_T·tan(s)* (altitude × beam-divergence half-angle ×
  tan slope) = footprint-radius × tan(slope). BUT — critical caveat — Gardner's own
  *range* biases are **centroid-based**: a symmetric footprint on a planar slope with a
  true-centroid estimator has **no first-order slope range bias**; slope enters as
  **broadening (variance)**, not bias. A signed slope range bias appears in Gardner only
  through slope × pointing-jitter coupling (his Table 2), and those are satellite-scale
  (huge-footprint) magnitudes. Take the *scaling* (∝ footprint·tan slope, growing
  faster-than-linear), not the numbers.

**Sign / magnitude / onset.** Sign is **downward** and unambiguous *for a
leading-edge/threshold detector* (Baltsavias/Deems). Magnitude is a large-slope,
large-footprint figure (~50 cm at **45°**, 1000 m AGL); the relation is **continuous**
(∝ tan s, accelerating at steep angles), **not** a step. No source states a discrete
~25–30° onset — that framing is unverified. The physical onset is where footprint span
*D·tan(s)* becomes comparable to the pulse length / detector timing resolution.

**Verdict (1): DOCUMENTED, downward, but framed at large slope / large footprint and as a
smooth tan-law, not a threshold.** Real and correctly signed; the published magnitude is
an upper-corner number, and our tens-of-mm-at-27° onset is not itself in the record.

---

## (2) Raw footprint geometry / leading-edge / first-photon — actually biases ground HIGH; downward needs an ESTIMATOR/PROCESSING route

The crux, and a **sign trap**: the naive "footprint spread → ground low" intuition is
backwards for the *raw geometry*. On a slope the **near / up-slope edge** of the footprint
is closer and returns energy **first**, so an early-triggering detector reports the range
**short → elevation HIGH**.

- **Laconte, Deschênes, Labussière & Pomerleau (2019), ICRA / arXiv:1810.01619.**
  VERIFIED — full text read. On an inclined surface the near part returns sooner, skewing
  the waveform so the detected peak is **earlier** than beam-center; a peak/leading
  detector reads range **short (toward sensor)** = elevation **high** for a
  downward-looking sensor. Bias reaches **~20 cm at high incidence** (their TLS geometry,
  range <10 m); grows monotonically with incidence.
- **Wang et al. (2018), diva-portal:1160255 (TLS beam-width bias).** VERIFIED — full text
  read. Verbatim: a low-fraction-of-maximum (leading-edge) detector on the gradual rise
  gives *"'too early' detection … the computed point is in front of the surface, biased
  towards the scanner."* Perpendicular incidence → all detector types agree; oblique →
  they diverge, leading-edge biases toward the scanner.
- **Footprint ellipse geometry** — Skyearth (2014), *"Quantitative Analysis on Geometric
  Size of LiDAR Footprint"* (cites Baltsavias as the first footprint model). VERIFIED —
  full text read. On an incline the footprint is an ellipse whose center is displaced
  from the incidence point; major axis → ∞ as incidence → grazing. The geometry gives the
  *span* (≈ *D·tan s*); it does **not** by itself pick a sign — the estimator does.

So a **DOWNWARD** ground bias on a slope is **not** the raw-geometry default; it requires
one of these **estimator / processing** routes:
1. **Fixed-threshold time-walk** on the broadened, lower-amplitude echo → crosses later →
   range long → **low** (mechanism (1); Baltsavias/Deems).
2. **Noise truncation of the weak up-slope leading tail** — the faint high-elevation
   shoulder is clipped as noise, shifting the retained centroid down-range → **low**
   (waveform terrain-slope-inversion literature, e.g. *Remote Sensing* 13(3):424, 2021;
   VERIFIED as a stated mechanism, magnitude not isolated).
3. **Min-elevation / lowest-return ground selection** picking the **down-slope (lower)
   edge** of the footprint or cell → surface **low**, magnitude up to ≈ **½·D·tan(s)**
   (see §3).

**Verdict (2): the raw geometric/leading-edge route biases ground HIGH, not low — a sign
trap to avoid.** A downward bias is real but arises from *processing* (threshold
time-walk, noise clipping, or min-return gridding), not from footprint geometry per se. No
single primary source cleanly isolates a slope→downward *waveform-ranging* bias with a
magnitude formula — a genuine gap.

---

## (3) Min-bin / lowest-return ground selection on slopes — DOWNWARD, geometrically certain, thinly published

The most likely culprit in a *gridded / classified* pipeline, and the cleanest downward
mechanism at the DTM (not ranging) stage. On a planar slope a cell spans a range of true
elevations; taking the **minimum / lowest return** in the cell selects a point **downhill
of the cell center**, pulling the surface **down**; the effect grows with slope and with
points-per-cell.

- **Su & Bork (2008), *Computers & Geosciences* 35(2):289–300** (companion to the 2006
  PE&RS paper). Lidar (~0.75 pts/m²) *"tended to increasingly underestimate terrain
  elevation as slope increased."* VERIFIED **only as a repeated paraphrase** across
  citing works — paywalled body not opened; treat as verified-attribution,
  unverified-verbatim.
- **Adams & Chandler (2002)** (lidar DTM of the Black Ven mudslide, Dorset): overall RMSE
  0.26 m, *"lidar data tended to increasingly underestimate terrain elevation as slope
  increased."* VERIFIED **only as a secondary citation**; primary not opened. The most
  direct textual sign+slope statement located.
- **Spaete et al. (2011), *Remote Sensing Letters* 2(4):317–326** (sagebrush steppe, SW
  Idaho). The cleanest **signed, slope-stratified** result: mean signed error (MSE)
  **−0.154 to +0.017 m**, **most negative (surface too LOW) on steep slopes + herbaceous
  cover**; RMSE 0.072–0.220 m. VERIFIED (journal/volume/year + MSE numbers, from
  concordant secondaries; **author list unconfirmed — 403 on primary — do not assert
  names**).
- **The geometric statement itself** is articulated cleanest as a *practitioner note*, not
  a peer-reviewed result: **lidR issue #51** (r-lidar, user jmmonnet), verbatim: *"in
  slope areas … if there are several lidar points per cell, the downhill point will be
  systematically selected, leading to a lower elevation in the case of a regular slope."*
  VERIFIED verbatim (GitHub). The tidal-marsh **minimum-bin-gridding** literature
  (Buffington et al. 2016, *RSE*) confirms min-selection's downward pull *quantitatively
  in the flat case* (it is used deliberately to cancel a positive vegetation bias) — the
  corollary being that on a **slope** that same pull becomes a slope-dependent
  **under**estimate.

**Verdict (3): DOWNWARD and geometrically certain, magnitude ≈ ½·D·tan(s) (or
½·cellsize·tan s), but under-published as a *signed* quantity.** Su & Bork 2008 and Adams
& Chandler 2002 state the sign; Spaete 2011 gives the cleanest signed slope table
(~−0.15 m worst case). The exact geometric argument survives only in a practitioner note —
a first-principles argument to make in our own writeup, empirically backed by these
papers.

> **NOTE — do NOT attribute the sign to Hodgson & Bresnahan (2004), *PE&RS*
> 70(3):331–339.** That paper supplies the **horizontal-error × tan(slope) → vertical
> error** coupling (elevation error at ~25° ≈ **2×** that at ~1.5°) and ranks horizontal
> displacement in the error budget — but it reports this as an **RMSE / magnitude** term
> (scatter), **not** a signed downward bias. tan(slope) coupling is one-directional only
> if the horizontal displacement is itself biased. Cite it for the *mechanism/magnitude
> channel*, not the sign. (VERIFIED as citation + concordant secondaries; open PDF 403.)

---

## (4) Circa-2005–2010 GPS/IMU navigation-era systematic error — SIGNED and swath-coherent (boresight/GPS), the strongest early-gen candidate

The distinction that matters for a **repeat-survey DoD**: two sub-mechanisms with
**opposite sign behavior**.

- **Random horizontal position error → vertical = δh·tan(slope): UNSIGNED scatter.** A
  mean-zero horizontal error inflates vertical *scatter* on slopes but carries no
  consistent up/down sign. This is the Hodgson & Bresnahan / Schaer channel. **Glennie
  (2007), *J. Applied Geodesy* 1(3):147–157** (DOI 10.1515/jag.2007.017; VERIFIED — full
  PDF read) is the variance-propagation basis: attitude errors dominate the **horizontal**
  budget (*"combined IMU error and boresighting error contribute from 60% to 75% of the
  overall horizontal error"*) and 25% to >50% of the vertical, growing with target range
  (altitude); it explicitly names *"incidence angle, and terrain slope"* as significant
  unmodeled error sources (future work). Treated as **scatter/covariance**, not signed
  bias.
- **Boresight / scan-mirror miscalibration and GPS trajectory drift: SIGNED,
  swath-coherent, scan-angle-correlated** — and this is the legacy-data mechanism that
  produces a systematic offset that can masquerade as change.
  **Glennie, Hinojosa-Corona, Nissen et al. (2014), *GRL* 41** (DOI 10.1002/2014GL059919).
  VERIFIED — full PDF read; quotes extracted directly. On the El Mayor–Cucapah legacy
  (2006) pre-event lidar: systematic errors were *"manifested as apparent slip"* in the
  differenced field; reprocessing the DGPS trajectory changed it by *"approximately 0.7 m
  3-D RMSE"* (vs cm-level expected); residual y-displacement *"shows a clear correlation
  with scan angle and a large 10 m difference between observations near nadir and those at
  the 30° extent of the scan,"* diagnosed as *"misalignment between the lidar system
  scanning mirror and the outgoing laser beam."* Errors *"are mainly at the edge of the
  flight line swaths (i.e., at larger scan angles)."* Full reprocessing cut them ~**25×**
  (to ~40 cm max horizontal). Verbatim on the slope coupling: *"As the angle of incidence
  increases (with increasing scan angle or as a result of topography), the estimated
  error of the lidar point increases as the intersection of the circular laser beam and
  the terrain becomes more and more elliptical."*
- **Habib et al. (strip adjustment):** boresight/lever-arm errors give strip-to-strip
  discrepancies; on **flat** terrain strip adjustment reduces to a 1-D vertical shift, so
  **sloped surfaces of varied orientation are what make the horizontal boresight error
  observable in the vertical** — i.e., slope is precisely what turns a horizontal
  calibration error into a signed vertical signal. VERIFIED thematically, not to specific
  numbers.

**Verdict (4): the strongest early-gen candidate — but its SIGN and slope-onset are
site/geometry-specific, not universal.** A **signed, slope-correlated, swath-coherent**
offset in ~2005–2010 lidar is well documented (Glennie 2014) and attributed to
boresight/scan-mirror miscalibration + GPS drift, **not** to the random-error projection
(unsigned scatter). Its diagnostic signature is correlation with **scan angle /
distance-from-flight-line / flight direction**, and concentration at **swath edges** — a
checkable prediction for our gen1 residual. It does **not** predict a fixed downward sign
at a fixed ~27° slope threshold; whether it reads *low* depends on the specific
mis-calibration and flight geometry.

---

## (5) The inter-epoch mechanism that most cleanly makes an OLD epoch differ on slopes — horizontal co-registration offset

The single most likely reason two epochs disagree **on slopes specifically**, and the one
to rule out before invoking any lidar-physics mechanism.

- **Nuth & Kääb (2011), *The Cryosphere* 5:271–290** (DOI 10.5194/tc-5-271-2011).
  VERIFIED — equation form confirmed from the paper and independently from the CNES
  `demcompare` and xDEM implementations. A horizontal misregistration between two DEMs
  produces an elevation difference

  **dh = a · cos(b − ω) · tan(α) + c**

  (α = slope, ω = aspect, a = shift magnitude, b = shift direction, c = mean vertical
  bias). The signal is **∝ tan(slope)** (zero on flat, accelerating on steep) and
  **flips sign with aspect** (positive on slopes facing the shift, negative on the
  opposite aspect). This is a **signed bias**, not scatter — and exactly the shape that
  masquerades as real change on slopes.
- **Kamp et al. (2024), *ESPL* 49** (DOI 10.1002/esp.5540). VERIFIED (open). Directly on
  point: multi-platform/multi-temporal DTM differences carry systematic vertical +
  horizontal trends *worst at steep slopes*; GPA-ICP + co-registration reduces the
  systematic DoD trend *"especially at steep slopes."*
- **Wheaton et al. (2010), *ESPL* 35(2):136–156** (DOI 10.1002/esp.1886). VERIFIED
  (open PDF). The foundational DoD-uncertainty framework — spatially variable uncertainty
  (fuzzy-inference on slope, point density, roughness) + spatial-coherence filter —
  *because* slope-correlated error is otherwise misread as geomorphic change. Frames slope
  as an **uncertainty** covariate, not a single-signed bias.

**Verdict (5): the cleanest, best-documented way an old epoch differs from a new one on
slopes is a sub-pixel HORIZONTAL co-registration offset (Nuth & Kääb), giving a signed
difference ∝ tan(slope) that flips with aspect.** This is aspect-dependent, so it is **not
a fixed downward bias** — but it is the first thing to test: if our +30–48 mm above ~27°
**flips sign between opposing aspects**, it is co-registration, not lidar physics. If it is
**one-signed across all aspects**, co-registration is excluded and a physics/processing
route (§1–3) is implicated.

---

## VERDICT — is an early-gen slope-dependent DOWNWARD ground bias documented?

**Partly, in disassembled form — and the naive version is a sign trap.** A downward
slope bias is real and documented as a *detector time-walk* effect (Baltsavias/Deems) and
as a *min-return gridding* effect (Su & Bork 2008; Spaete 2011; lidR #51), but:

1. The **raw footprint geometry / leading-edge / first-photon** route biases the ground
   **HIGH, not low** on slopes (Laconte 2019; Wang 2018). Downward requires a *processing*
   route (threshold time-walk, noise clipping, or min-return selection) — do not attribute
   a downward bias to footprint geometry per se.
2. The documented downward magnitudes are **large-slope / large-footprint** (~50 cm at
   45°, satellite-scale in Gardner) or **DTM-filter** artifacts (~−0.15 m at steep +
   herbaceous, Spaete). Our **tens-of-mm at ~27°** is smaller and not itself in the record.
3. **No source documents a discrete ~25–30° threshold/jump.** Every quantified mechanism
   is a **smooth tan-law** that accelerates at steep angles (tan 27° = 0.51, tan 45° =
   1.0). A step-onset at 27° would be a *finding of ours*, not a confirmation of the
   literature — treat the "switch-on" language as our hypothesis, and interrogate whether
   the apparent threshold is a smooth tan-curve crossing our detection floor.
4. The strongest early-gen candidate is **legacy GPS/IMU boresight/scan-mirror
   miscalibration + GPS drift** (Glennie 2014), which is **signed and swath-coherent** but
   whose **sign and slope dependence are site-specific**, diagnosable by scan-angle /
   flight-direction correlation and swath-edge concentration.
5. The **first thing to exclude** is a horizontal **co-registration offset** (Nuth & Kääb
   2011): a signed difference ∝ tan(slope) that **flips with aspect**. Test aspect
   symmetry before invoking any lidar-physics route.

**Present as our contribution (not attributed):** a *signed, slope-resolved,
tens-of-mm* gen1-low bias at moderate slope, if it (a) is one-signed across aspects
(excluding co-registration), (b) does not correlate with scan angle / swath position
(excluding boresight), and (c) turns on near ~25–30° — none of which is a documented
early-gen result at this magnitude/slope. The mechanisms above are the *candidate
explanations to test against*, not confirmations.

---

## Verification status summary

- **VERIFIED (full text read):** Deems, Painter & Finnegan 2013 (verbatim time-walk quote
  + bib); Gardner 1992 (NASA scan, broadening + centroid-bias caveat); Glennie 2007 (full
  PDF, attitude budget quotes); Glennie et al. 2014 (full PDF, "apparent slip" / 0.7 m /
  scan-angle-correlation quotes); Laconte et al. 2019 (arXiv); Wang et al. 2018
  (diva-portal); Skyearth 2014 (footprint ellipse); lidR issue #51 (verbatim).
- **VERIFIED (equation/method confirmed from paper + independent implementations):** Nuth
  & Kääb 2011 (via demcompare + xDEM).
- **VERIFIED (open text / metadata):** Kamp et al. 2024; Wheaton et al. 2010.
- **VERIFIED as citation, body paywalled (content NOT independently confirmed):**
  Baltsavias 1999 (quoted through Deems 2013); Hodgson & Bresnahan 2004 (RMSE/2×-at-25°
  via concordant secondaries; sign is scatter, not downward bias); Aguilar & Mills 2010
  (slope-dependent gridding *variance*, sign not confirmed).
- **VERIFIED only as repeated paraphrase / secondary (primary NOT opened):** Su & Bork
  2008 ("underestimate as slope increased"); Adams & Chandler 2002 (same wording); Spaete
  et al. 2011 (MSE numbers solid, **author list UNCONFIRMED — do not assert names**);
  Buffington et al. 2016 (min-bin gridding, flat-case downward pull).
- **NOT LOCATED / flagged as gaps:** a Pfeifer/Mandlburger/Hyyppä paper isolating a
  slope-dependent leading-edge low bias; a paper quantifying full-waveform's *reduction*
  of a slope low-bias; any source documenting a discrete ~25–30° threshold onset; a
  peer-reviewed (vs practitioner-note) statement of the min-return-picks-downhill-edge
  geometry.
- **No fabricated citations or magnitudes.** Every unverified specific is flagged inline.

## Full citation list

1. Adams, J.C. & Chandler, J.H. (2002). Evaluation of lidar and medium-scale
   photogrammetry for detecting soft-cliff coastal change. *Photogrammetric Record*
   17(99):405–418. Signed "underestimate elevation as slope increased." VERIFIED as
   secondary citation only; primary not opened.
2. Aguilar, F.J. & Mills, J.P. (2008). Accuracy assessment of lidar-derived digital
   elevation models. *Photogrammetric Record* 23(122):148–169. And Aguilar et al. (2010),
   Modelling vertical error in LiDAR-derived DEMs, *ISPRS J. P&RS* 65:103–110
   (DOI 10.1016/j.isprsjprs.2009.10.005). VERIFIED (citation; slope-dependent gridding
   *variance*, sign not confirmed).
3. Baltsavias, E.P. (1999). Airborne laser scanning: basic relations and formulas.
   *ISPRS J. P&RS* 54(2–3):199–214. DOI 10.1016/S0924-2716(99)00015-5. VERIFIED as
   citation; time-walk result quoted verbatim through Deems et al. 2013. Body paywalled.
4. Buffington, K.J., Dugger, B.D., Thorne, K.M. & Takekawa, J.Y. (2016). Statistical
   correction of lidar-derived DEMs (minimum-bin gridding), *RSE* 186:616–625. VERIFIED as
   secondary (min-selection downward pull in the flat case).
5. Deems, J.S., Painter, T.H. & Finnegan, D.C. (2013). Lidar measurement of snow depth: a
   review. *Journal of Glaciology* 59(215):467–479. doi:10.3189/2013JoG12J154. VERIFIED
   (full text; verbatim Baltsavias time-walk quote, ~50 cm at 45°/1000 m).
6. Gardner, C.S. (1992). Ranging performance of satellite laser altimeters. *IEEE TGRS*
   30(5):1061–1072. DOI 10.1109/36.175341. VERIFIED (full NASA scan; broadening ∝
   footprint·tan s, centroid → no first-order bias; slope×jitter coupling only).
7. Glennie, C. (2007). Rigorous 3D error analysis of kinematic scanning lidar systems.
   *J. Applied Geodesy* 1(3):147–157. DOI 10.1515/jag.2007.017. VERIFIED (full text;
   attitude 60–75% of horizontal budget; names incidence angle + terrain slope).
8. Glennie, C., Hinojosa-Corona, A., Nissen, E., Kusari, A., Oskin, M.E., Arrowsmith, J.R.
   & Borsa, A. (2014). Optimization of legacy lidar data sets for measuring near-field
   deformation. *GRL* 41. DOI 10.1002/2014GL059919. VERIFIED (full text; "apparent slip,"
   ~0.7 m trajectory RMSE, scan-angle-correlated boresight residual, ~25× reduction).
9. Hodgson, M.E. & Bresnahan, P. (2004). Accuracy of airborne lidar-derived elevation:
   empirical assessment and error budget. *PE&RS* 70(3):331–339. DOI 10.14358/PERS.70.3.331.
   VERIFIED (citation + secondaries). tan(slope) coupling, ~2× at 25°; **RMSE, NOT a
   signed bias** — do not attribute the downward sign here.
10. Kamp, N. et al. (2024). Comparability of multi-temporal DTMs derived from different
    LiDAR platforms. *Earth Surf. Process. Landf.* 49. DOI 10.1002/esp.5540. VERIFIED
    (open; systematic DoD trend worst at steep slopes, co-registration reduces it there).
11. Laconte, J., Deschênes, S.-P., Labussière, M. & Pomerleau, F. (2019). Lidar
    measurement bias estimation via return waveform modelling. *IEEE ICRA* / arXiv:1810.01619.
    VERIFIED (full text; near/up-slope edge returns first → range short → elevation HIGH;
    ~20 cm at high incidence, TLS).
12. Nuth, C. & Kääb, A. (2011). Co-registration and bias corrections of satellite
    elevation data sets. *The Cryosphere* 5:271–290. DOI 10.5194/tc-5-271-2011. VERIFIED
    (dh = a·cos(b−ω)·tan α + c; signed, ∝ tan slope, aspect-flipping).
13. Schaer, P., Skaloud, J., Landtwing, S. & Legat, K. (2007). Accuracy estimation for
    laser point cloud including scanning geometry. *5th ISPRS MMT*, Padua. VERIFIED as
    citation (per-point variance from scanning geometry; scatter, not signed bias).
14. Skyearth (2014). Quantitative analysis on geometric size of LiDAR footprint.
    skyearth.org/publication/papers/2014_qagslf.pdf. VERIFIED (full text; footprint
    ellipse center offset, major axis → ∞ at grazing).
15. Spaete, L.P. et al. (2011). Vegetation and slope effects on accuracy of a
    LiDAR-derived DEM in the sagebrush steppe. *Remote Sensing Letters* 2(4):317–326.
    DOI 10.1080/01431161.2010.515267. VERIFIED (MSE −0.154 to +0.017 m, most negative on
    steep+herbaceous). **Author list UNCONFIRMED — do not assert names.**
16. Su, J. & Bork, E. (2006). Influence of vegetation, slope, and lidar sampling angle on
    DEM accuracy. *PE&RS* 72(11):1265–1274; and Su & Bork (2008), Evaluating error
    associated with lidar-derived DEM interpolation, *Computers & Geosciences*
    35(2):289–300 (DOI 10.1016/j.cageo.2008.09.001). VERIFIED as paraphrase (RMSE ~2×
    above 10°; 2008: "increasingly underestimate terrain elevation as slope increased").
    Verbatim/exact per-slope numbers not confirmed (paywalled).
17. Wang, Y. et al. (2018). Bias of cylinder diameter estimation from ground-based laser
    scanners with different beam widths. diva-portal.org 1160255. VERIFIED (full text;
    leading-edge → "too early" → point in front of surface, toward scanner).
18. Wheaton, J.M., Brasington, J., Darby, S.E. & Sear, D.A. (2010). Accounting for
    uncertainty in DEMs from repeat topographic surveys. *ESPL* 35(2):136–156.
    DOI 10.1002/esp.1886. VERIFIED (open; spatially variable DoD uncertainty; slope as
    uncertainty covariate).
19. lidR issue #51 (r-lidar/lidR, GitHub). VERIFIED verbatim: min-in-cell selects the
    downhill point on a regular slope → lower elevation. Practitioner note, not
    peer-reviewed.

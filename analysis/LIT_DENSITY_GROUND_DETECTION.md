# Point density, order statistics, and the bias of the lidar ground estimate under vegetation

**Literature search, 2026-09-23.** Question: does the *number of returns per cell* change the
*bias* of the estimated ground elevation under vegetation — as opposed to its precision — and
is the order-statistics reason for that (the expected minimum of `n` draws falls with `n`)
anywhere in the airborne-lidar literature?

No code was written and no data processed for this document.

## Provenance tags used below

- **[FULL TEXT — this session]** — PDF retrieved and read in this session.
- **[ABSTRACT — this session]** — abstract retrieved verbatim this session (Semantic Scholar
  Graph API on the DOI, or a fetched page); body not read.
- **[METADATA ONLY]** — DOI/title/venue/volume/pages confirmed via Crossref; no abstract or
  body obtained. No numbers quoted from these.
- **[REPO RECORD]** — already recorded in this repository by an earlier session with a
  read-level tag; **not re-read in this session**. Quotes are reproduced from the repo note,
  and carry that note's own tag.
- **[NOT VERIFIED]** — could not confirm; nothing is claimed from it.

Every quantitative statement below is a quotation from a source. Nothing is paraphrased into
a number.

---

## What this repository already holds (checked before searching)

Four existing notes cover adjacent ground and were read first, so this document does not
re-find them:

| Repo note | What it settles |
|---|---|
| `analysis/ridgelines/DTM_PERCENTILE_LITERATURE.md` (1212 lines) | The cover-dependent vertical shift of the ground surface, and whether a per-cell percentile can compensate it. Its "Gaps" section already records that **no peer-reviewed study compares per-cell order statistics (min / low percentile / median / mean) for a bare-earth DTM**, and that **no study fixes cell size and tunes the percentile** — all published tuning is in *scale/window*, "conflating the quantile with the support." |
| `analysis/GROUND_RETURN_PDF_LITERATURE.md` (573 lines) | No established PDF exists for the vertical distribution of discrete returns about a ground surface; the exponential vegetation tail is derived from Beer–Lambert with named assumptions. |
| `analysis/ridgelines/RETURN_HEIGHT_DISTRIBUTION_LITERATURE.md` | Ground-return tail *form* by land cover; forest = exponential, open = power-law (the latter apparently novel). |
| `analysis/ridgelines/CSF_PARAMETER_LITERATURE.md` | Ground-filter behaviour vs. density, including the Sithole & Vosselman benchmark densities and the Štular & Lozić repeat. |
| `analysis/decimation_result.md` | Our own measurement (the thing being searched against). |

The present document adds the **density → bias** axis, which those notes do not cover.

---

## Bottom line

**1. The density-vs-ground-BIAS question is NOT settled, and is barely posed.** The
literature that decimates point clouds overwhelmingly reports *precision* — RMSE against
checkpoints, or "random variation" — and where a signed mean error is reported at all, it is
reported as a by-product rather than as the object of study. Two papers do separate bias from
scatter as a function of density (Leitold et al. 2015; Fradette et al. 2019), and they reach
*opposite* conclusions, for reasons that are exactly the order-statistics distinction we
drew: Leitold's bias comes from a ground **filter** misclassifying vegetation at low density;
Fradette's DTM, referenced against a 21 pulses/m² DTM, "did not need specific adjustment"
while the CHM — a per-cell **maximum** — did.

**2. The order-statistics framing exists in this literature, but only for the CANOPY TOP, and
it was published in 2017.** Roussel et al. (2017, *RSE*) derive the expected maximum of `n`
pulses over a discretised canopy-height histogram and recompute `hmax` "as if the density of
pulses were infinite." That is precisely the mirror image of our argument. **I found no
application of it to the ground / minimum / low quantile.** Roussel et al. explicitly take
the true maximum, not the ground, as their reference, and their DTM was vendor-supplied.

**3. The one place the field does use a low order statistic for the ground — minimum-bin
gridding in marshes, and local-minima retrieval in forest — tunes the *window size*, never
`n` at fixed window.** Cell size and sample size are therefore confounded in every published
tuning. Clark et al. (2004) come closest to naming the confound: their optimum is stated as
specific to "this tropical landscape **and lidar sampling density**."

**4. Seasonal and herbaceous-cover work exists and is consistent, but does not touch
density.** Simpson et al. (2017) and Stereńczak & Kozak (2011) are the leaf-on/leaf-off
anchors; Hopkinson et al. (2005) is the best herbaceous/wetland analogue.

**5. The vendor's delivered class-2 surface floating above true ground in proportion to
near-ground vegetation IS measured against survey control, repeatedly, in wetlands.**
Hopkinson et al. (2005) is the cleanest: +0.07 m mean for ground-classified returns over
vegetated surfaces, resolved by cover class, with the mechanism named as "the reduced
probability of LiDAR pulse penetration to true ground level."

---

## 1. Density vs. bias, as opposed to density vs. precision

### 1.1 What the decimation literature actually reports: scatter

**Hansen, E.H., Gobakken, T., & Næsset, E. (2015). Effects of Pulse Density on Digital
Terrain Models and Canopy Metrics Using Airborne Laser Scanning in a Tropical Rainforest.
*Remote Sensing*, 7(7), 8453–8468.** doi:[10.3390/rs70708453](https://doi.org/10.3390/rs70708453)
— **[ABSTRACT — this session]**

Verbatim:

> "We used a total of 612 coordinates measured with a differential dual frequency Global
> Navigation Satellite System receiver to analyze the effects on DTMs at pulse densities of
> 8, 4, 2, 1, 0.5, and 0.025 pulses·m−2."

> "**Random variation** in DTMs and canopy metrics **increased** with reduced pulse density."

This is the field's default framing, stated in a single word: the reported quantity is
*random variation*. The abstract contains no signed-bias statement for the DTM.

**Cățeanu, M., & Ciubotaru, A. (2021). The Effect of LiDAR Sampling Density on DTM Accuracy
for Areas with Heavy Forest Cover. *Forests*, 12(3), 265.**
doi:[10.3390/f12030265](https://doi.org/10.3390/f12030265) — **[ABSTRACT — this session]**
(MDPI returned 403 to every full-text fetch attempt in this session, so the body was not read.)

This is the study whose density range matches ours most closely. Verbatim:

> "To analyze its impact on the quality of ground surface modelling, the density of the
> filtered data set was artificially lowered (from 0.89 to 0.09 points/m2) by randomly
> removing point observations in 10% increments."

> "While the reduction of point density leads to a less accurate DTM in all cases (as
> expected), the exact pattern varies by algorithm. The accuracy of the LiDAR-derived DTMs is
> relatively good even when LiDAR sampling density is reduced to 0.40–0.50 points/m2 (50–60 %
> of the initial point density), as long as a suitable interpolation algorithm is used."

*(a) What the source says:* density reduction degrades accuracy; the pattern is
algorithm-dependent. *(b) What I infer:* the treatment here is the **interpolator** (IDW / NN
/ TPS) applied to points already classified as ground. That is a smoother, not an order
statistic, so this design cannot see the effect we measured even in principle. The abstract
reports no signed bias.

**Raber, G.T., Jensen, J.R., Hodgson, M.E., Tullis, J.A., Davis, B.A., & Berglund, J. (2007).
Impact of Lidar Nominal Post-spacing on DEM Accuracy and Flood Zone Delineation.
*Photogrammetric Engineering & Remote Sensing*, 73(7), 793–804.**
doi:[10.14358/PERS.73.7.793](https://doi.org/10.14358/PERS.73.7.793) — **[ABSTRACT — this session]**

The canonical systematic-decimation study. Verbatim:

> "Lidar data collected at a low post-spacing (approximately 1 to 2 m) over a piedmont study
> area in North Carolina was systematically decimated to simulate datasets with sequentially
> higher post-spacing values. Using extensive first-order ground survey information, the
> accuracy of each DEM derived from these lidar datasets was assessed and reported."

> "The results indicate that base flood elevation does not statistically change over the
> post-spacing values tested."

The abstract reports "accuracy," and the downstream test is a flood elevation. No bias
decomposition, and the site is piedmont rather than vegetated floodplain.

### 1.2 The two studies that DO separate bias — and disagree, informatively

**Leitold, V., Keller, M., Morton, D.C., Cook, B.D., & Shimabukuro, Y.E. (2015). Airborne
lidar-based estimates of tropical forest structure in complex terrain: opportunities and
trade-offs for REDD+. *Carbon Balance and Management*, 10:3.**
doi:[10.1186/s13021-015-0013-x](https://doi.org/10.1186/s13021-015-0013-x)
— **[ABSTRACT verbatim + body via fetched full text — this session]**

Abstract, verbatim:

> "The terrain model generated from full-density (~20 returns m−2) data was highly accurate
> (**mean signed error of 0.19 ± 0.97 m**), while those derived from reduced-density datasets
> (8 m−2, 4 m−2, 2 m−2 and 1 m−2) were increasingly less accurate."

From the body (quoted in the fetch of the PMC full text):

> "mean signed errors of the thinned DTM elevations increased as data density was reduced
> from D20 (0.19 ± 0.97 m) to D1 data (3.21 ± 3.12 m)"

> "DTM elevations were higher than the GNSS elevations in all cases, with increasing error
> magnitudes as data were thinned"

RMSE progression reported alongside: D20 = 0.97 m; D8 = 1.35 m; D4 = 2.30 m; D2 = 3.47 m;
D1 = 4.45 m. Thinning was random, to 8/4/2/1 returns m⁻², at 10 × 10 m resolution. Ground
classification was **progressive morphological filtering**; the fetch found no mention of a
minimum or quantile ground estimator.

*(a) What the source says:* thinning drives a large **positive** ground-elevation bias, from
+0.19 m to +3.21 m, in dense multilayered tropical montane forest.
*(b) What I infer:* this is a **filter-failure** mechanism, not an order-statistic mechanism.
The paper itself attributes the DTM errors in thinned data to "an incorrect classification of
vegetation features as ground" (as reported in the search result summary — I did not locate
that sentence in the fetched text, so treat that attribution as **[NOT VERIFIED]**, though
the direction and the filter type are verified). The magnitudes — metres — are an order of
magnitude above anything an order-statistic argument would produce on a floodplain, which is
itself evidence that a different mechanism dominates at 1–8 pts/m² under closed tropical
canopy.

**Fradette, M.-S., Leboeuf, A., Riopel, M., & Bégin, J. (2019). Method to Reduce the Bias on
Digital Terrain Model and Canopy Height Model from LiDAR Data. *Remote Sensing*, 11(7), 863.**
doi:[10.3390/rs11070863](https://doi.org/10.3390/rs11070863) — **[ABSTRACT — this session]**
(also **[REPO RECORD]** in `DTM_PERCENTILE_LITERATURE.md`, abstract-level)

This is the most directly diagnostic paper found. Verbatim:

> "Underestimation of LiDAR heights is widely known but has never been evaluated for several
> sensors and for diverse types of ecological conditions. **This underestimation is mainly
> linked to the probability of the pulse to reach the ground and the top of vegetation.**
> Main causes of this underestimation are pulse density, pattern of scan (sensors), scan
> angles, specific contract parameters (flying altitude, pulse repetition frequency) and
> characteristics of the territory (slopes, stand density and species composition)."

> "For this study, the bias of both DTM and CHM were calculated by subtracting two LiDAR
> datasets: high-density pixels with 21 pulses/m² (first return) and more (DTM or CHM
> reference value pixels) and low-density pixels (DTM or CHM value to correct). **After
> preliminary analyses, it was concluded that the DTM did not need specific adjustment. In
> contrast, the CHM needed adjustments.**"

> "Among the variables studied, three were selected for the final CHM adjustment model: the
> maximum height of the pixel (H2Corr); **the density of first returns by m2 (D_first)**; and
> the standard deviation of nine maximum heights of the neighborhood cells (H_STD9)."

*(a) What the source says:* at 1 × 1 m, with a 21 pulses/m² reference, the density-dependent
bias is in the **canopy maximum**, not in the DTM; and the correction model for the maximum
takes first-return density as an explicit predictor.
*(b) What I infer, and it is the central inference of this document:* this is our result with
the sign flipped and the estimator flipped. Their CHM is a per-cell **maximum** — an extreme
order statistic — and it moves with `n`. Their DTM is not an extreme order statistic, and it
does not. Neither of the two papers names order statistics; Fradette et al. reach the
empirical fact without the framework, and their stated mechanism is "the probability of the
pulse to reach the ground and the top of vegetation," which is a per-pulse statement, not an
`n`-dependence.

---

## 2. Is the order-statistics framing used anywhere?

### 2.1 Yes — for the canopy top, and only there

**Roussel, J.-R., Caspersen, J., Béland, M., Thomas, S., & Achim, A. (2017). Removing bias
from LiDAR-based estimates of canopy height: Accounting for the effects of pulse density and
footprint size. *Remote Sensing of Environment*, 198, 1–16.**
doi:[10.1016/j.rse.2017.05.032](https://doi.org/10.1016/j.rse.2017.05.032)
Open post-print: <https://utoronto.scholaris.ca/server/api/core/bitstreams/a4856d3b-689f-435f-a19a-1e9d54d80567/content>
— **[FULL TEXT — this session]**

Abstract, verbatim:

> "The resulting variation between and within datasets (particularly variation in pulse
> density and footprint size) can induce spurious variation in LiDAR metrics such as maximum
> height (h max) and mean height of the canopy surface model (C mean). In this study, we
> first compare two LiDAR datasets acquired with different parameters, and observe that h max
> and C mean are **56 cm and 1.0 m higher**, respectively, when calculated using the
> high-density dataset with a small footprint. Then, we present **a model that explains the
> observed bias using probability theory, and allows us to recompute the metrics as if the
> density of pulses were infinite** and the size of the two footprints were equivalent."

The conceptual statement, §2.5, verbatim:

> "Figure 5a illustrates how the bias between the observed maximum height (ĥmax) and the true
> maximum (hmax), i.e. the actual highest point of the plot, increases as pulse density
> decreases. When 21 pulses reach the canopy, the observed maximum height (ĥmax,1)
> underestimates the true maximum by the amount ∆h1 = hmax − ĥmax,1. In contrast, when only
> 11 pulses reach the canopy (i.e. after removing every second pulse), the observed maximum
> height (ĥmax,2) is even lower, and underestimates the true maximum by the amount
> hmax − ĥmax,2 > ∆h1."

The derivation, §3.3.3 "A continuous canopy," verbatim (transcribed from the PDF text layer;
subscripts flattened):

> "As shown in figure 10a), a continuous canopy can be discretized using a histogram with k
> bins, one for each height (h i) and probability (p i), where i ∈ 0, k. When sampled with n
> pulses, the expected value of h max is:
>
>     h̄(n) = P_k^n h_k + (1 − P_k^n)[ P_{k−1}^n h_{k−1} + (1 − P_{k−1}^n)(P_{k−2}^n h_{k−2} + … ) ]     (11)"

with the recursion

>     H_i^n(H) = { h_0                                  if i = 0
>                { P_i^n h_i + (1 − P_i^n) H_i^{n−1}    else                                            (12)
>     h̄(n) = H_k^n(H)                                                                                   (13)

and, earlier, the statement that

> "sampling with more than one pulse is a Binomial process."

This is the argument we made, verbatim in structure: discretise the return-height
distribution, raise the tail probability to the `n`, and the expected extreme moves
monotonically with `n`.

**It is not applied to the ground.** The paper says so directly, §2.6:

> "Since our goal is to quantify the bias between the true maximum height and the observed
> value, **our point of reference is not the ground** but the true maximum height in a given
> pixel area."

And their terrain surface was not their own product (§2.2):

> "the digital terrain model) was done by the provider and we … fied as 'ground', although we
> could not obtain further details"

(the PDF text layer is broken across a column boundary there; the readable content is that
the DTM and the ground classification came from the data provider).

*(b) What I infer:* Roussel et al. established the order-statistics correction for the upper
extreme in 2017 and it has been taken up — Fradette et al. (2019) used their 21 pulses/m²
threshold as a reference density. Nobody appears to have run the same argument downward onto
the ground. That is not surprising: the community's ground estimate is a *filter plus
interpolator*, not an order statistic, so there was no `n`-dependent extreme to correct. The
argument only becomes necessary once you *choose* a low-quantile ground estimator.

### 2.2 No — nowhere for the ground

Searched, in this session, against: `order statistics` + lidar ground/minimum elevation;
`extreme value` + lidar terrain minimum; `expected minimum of n` + ground return; `lowest
return` + number of returns per cell; `probability of ground return` / `at least one pulse
reaches ground` + per-cell pulse count. **No source was retrieved that models the ground
estimate as an order statistic with explicit `n`-dependence.** The nearest returns were
minimum-bin gridding papers (§3 below), which change `n` only by changing the window.

A search null is not a proof. The honest closure would be a citation-graph query on works
citing Roussel et al. (2017) and Schmid et al. (2011), which I could not run without
Scopus/WoS access.

### 2.3 The Beer–Lambert / penetration literature is per-pulse, not per-`n`

`GROUND_RETURN_PDF_LITERATURE.md` **[REPO RECORD]** already establishes that the gap
probability `k(N) = exp(−γ·λ·A_s/cosθ)` (He & Lyu 2023, Eq. 1) is a statement about **one**
pulse's chance of reaching the ground. Nothing in that tradition compounds it over the `n`
pulses in a cell to give the distribution of the cell minimum. Fradette et al.'s
"probability of the pulse to reach the ground" (§1.2 above) is the same per-pulse statement.
The step from a per-pulse probability to an `n`-th order statistic is the step that does not
appear to have been taken for the ground.

---

## 3. Ground-filter choice interacting with density

### 3.1 Sithole & Vosselman is the canonical filter test — and it could NOT resolve density

Verified as canonical; and the specific point is that the benchmark *declined* to answer the
density question.

**Sithole, G., & Vosselman, G. (2004). Experimental comparison of filter algorithms for
bare-Earth extraction from airborne laser scanning point clouds. *ISPRS Journal of
Photogrammetry and Remote Sensing*, 59(1–2), 85–101.**
doi:[10.1016/j.isprsjprs.2004.05.004](https://doi.org/10.1016/j.isprsjprs.2004.05.004)
— **[METADATA ONLY — this session; Semantic Scholar returned no abstract]**;
**[REPO RECORD]** in `CSF_PARAMETER_LITERATURE.md`, which records from the abstract:

> "The influence of point density could not well be determined in this experiment."

and from the 2003 ISPRS Archives proceedings version (**[REPO RECORD, FULL TEXT]**,
<https://www.isprs.org/proceedings/xxxiv/3-w13/papers/Sithole_ALSDD2003.pdf>):

> "The urban sites were at a resolution of 1-1.5m. The rural sites were at a resolution of
> 2-3.5m."

> "More tests on decreasing resolution will need to be done, as the test sites chosen have
> proved …" (§5.2, truncated in the repo record)

So: the canonical filter comparison is at roughly our density (0.08–1 pts/m² per the repo's
reading of the benchmark documentation), and it explicitly reports that it **could not**
determine the density influence.

### 3.2 The later repeat at reduced density

**Štular, B., & Lozić, E. (2020)** — recorded in `CSF_PARAMETER_LITERATURE.md`
**[REPO RECORD, FULL TEXT]** — four ALS sites spanning 1.83 to 18.79 pts/m². Verbatim from
the repo record:

> "Concerning the influence of point density on filter performance, our results are
> consistent with the results presented by Sithole and Vosselman [20], namely in that the
> lower the point density, the worse the performance of all filters. This makes processing
> low-density data much more demanding."

> "Also, the removal of vegetation on steep slopes remains problematic, but with dense data
> (AT) it is solved well by all filters, as is the discontinuity retention."

*(a) What the source says:* filter performance degrades with density, for every filter tested.
*(b) What I infer:* "performance" here is classification quality (Type I/II error), not
elevation bias. The filter-vs-density axis and the estimator-vs-density axis are different
questions, and only the former has been studied.

---

## 4. Seasonal (leaf-on / leaf-off) and low herbaceous cover

### 4.1 Leaf state

**Simpson, J.E., Smith, T.E.L., & Wooster, M.J. (2017). Assessment of Errors Caused by Forest
Vegetation Structure in Airborne LiDAR-Derived DTMs. *Remote Sensing*, 9(11), 1101.**
doi:[10.3390/rs9111101](https://doi.org/10.3390/rs9111101) — **[ABSTRACT — this session]**;
**[REPO RECORD, FULL TEXT]** in `DTM_PERCENTILE_LITERATURE.md`

Verbatim:

> "Here, we use ground survey equipment to assess digital terrain model (DTM) accuracy in a
> deciduous broadleaf forest, during both leaf-on and leaf-off conditions. … In the presence
> of leaf-on vegetation, DTM accuracy is severely reduced, **with low-stature undergrowth
> vegetation (such as ferns) causing the greatest errors (RMSE > 1 m)**. Errors are lower
> under leaf-off conditions (RMSE = 0.22 m)"

Note that the largest errors come from **low-stature undergrowth**, not canopy — which is the
regime our floodplain sits in.

**Stereńczak, K., & Kozak, J. (2011). Evaluation of digital terrain models generated in forest
conditions from airborne laser scanning data acquired in two seasons. *Scandinavian Journal of
Forest Research*, 26(4), 374–384.**
doi:[10.1080/02827581.2011.570781](https://doi.org/10.1080/02827581.2011.570781)
— **[ABSTRACT — this session]**

Verbatim:

> "For various output raster resolutions and seasons of data acquisition, **mean errors varied
> between −0.2 and 0.34 m**, and root mean square errors varied from 0.28 to 0.79 m. Errors
> increased linearly with DTM pixel size, and their variability was significantly higher in
> DTMs derived from summer data than in DTMs derived from spring data."

This is one of the few sources that reports a *signed mean error* range alongside RMSE. It
varies with **pixel size** and season, not with point density — again, scale not `n`.

### 4.2 Herbaceous / wetland cover: the closest analogue to a floodplain

**Hopkinson, C., Chasmer, L., Zsigovics, G., Creed, I.F., Sitar, M., Treitz, P., & Maher, R.V.
Errors in LiDAR ground elevation and wetland vegetation height estimates.** *International
Archives of Photogrammetry, Remote Sensing and Spatial Information Sciences*, Vol. XXXVI-8/W2,
108–113. **Conference archives, not journal peer review.**
<https://www.isprs.org/proceedings/xxxvi/8-w2/HOPKINSON.pdf> — **[FULL TEXT — this session]**
(The peer-reviewed companion is Hopkinson et al. 2005, *Can. J. Remote Sens.* 31(2), 191–206,
doi:[10.5589/m05-007](https://doi.org/10.5589/m05-007) — **[METADATA ONLY]**, no abstract
available.)

Density, verbatim:

> "The point density of ground-classified returns over the survey polygon was between 0.5 and
> 1.3 returns per m2. However, the ground class point density was highly variable, with higher
> densities over open dry ground and lower densities over areas of dense canopy and wet
> ground."

That is essentially our gen1 ground density. The instrument check, verbatim:

> "After subtracting highway RP elevations from LiDAR returns within a 0.5 m horizontal
> radius, the mean height difference was found to be 0.00 m with a standard deviation of 0.07
> m."

Table 1, transcribed from the PDF (mean height error in m, ground-classified raw LiDAR vs.
GPS reference points; second block is the raster ground DEM):

| Statistic | Hwy data | Ground (all veg.) | Aquatic | Grass/herbs | Low shrubs | Tall shrubs |
|---|--:|--:|--:|--:|--:|--:|
| Raw LiDAR mean | 0.00 | **+0.07** | **+0.15** | **+0.02** | +0.06 | +0.06 |
| Raw LiDAR St.Dev. | 0.07 | 0.16 | 0.22 | 0.10 | 0.12 | 0.03 |
| Raw LiDAR N | 95 | 127 | 35 | 45 | 43 | 4 |
| Raw LiDAR P | ND | < 0.01 | < 0.01 | ND | < 0.01 | < 0.01 |
| Raster mean | 0.00 | **+0.04** | **+0.12** | 0.00 | +0.01 | +0.11 |
| Raster St.Dev. | 0.08 | 0.14 | 0.18 | 0.10 | 0.13 | 0.16 |
| Raster N | 95 | 208 | 45 | 77 | 72 | 14 |

Verbatim text:

> "The average difference between ground classified LiDAR returns and RPs collected from
> vegetated surfaces was **+0.07 m (± 0.16 m)**. Similar results were returned for the raster
> DEM and ground control point comparison with a mean offset of **+0.04 m (± 0.14 m)**
> (Table 1). LiDAR pulses lying above RPs were also observed in vegetated areas by Töyrä et
> al. (2003) and **can be explained by the reduced probability of LiDAR pulse penetration to
> true ground level.**"

> "It is apparent that the vertical offset varied with vegetation type, with negligible
> offsets of +0.02 m (raw LiDAR) and 0.00 m (raster DEM) for the grass/herb class and
> relatively large offsets of +0.15 m (LiDAR) and +0.12 m (raster) for the aquatic vegetation
> class (Table 1)."

The paper also rules out a classification confound, verbatim:

> "for this study, the classification algorithm was kept spatially constant and any variation
> in the magnitude of vertical error can be attributed to localised ground cover conditions."

*(b) What I infer:* the grass/herb class being *statistically indistinguishable from zero*
(P = ND) at 0.5–1.3 ground pts/m² is a useful calibration on expectations — at gen1-like
density, an herbaceous offset of a few centimetres is at the edge of detectability with ~45
reference points. It does not contradict a low-tail effect; it bounds what a *mean-of-class*
comparison can resolve.

---

## 5. The delivered class-2 surface floating above true ground

Established, repeatedly, in wetlands and marshes, measured against survey control. In order
of directness:

1. **Hopkinson et al. (ISPRS Archives XXXVI-8/W2)** — §4.2 above. Ground-*classified* returns,
   +0.07 m over vegetated surfaces, resolved by cover class, mechanism named. **[FULL TEXT]**
2. **Ewald, M.J. (2013).** *Where's the Ground Surface? Elevation Bias in LIDAR-derived Digital
   Elevation Models Due to Dense Vegetation in Oregon Tidal Marshes.* **MS thesis, not
   peer-reviewed.** <https://ir.library.oregonstate.edu/downloads/1n79h8198> —
   **[REPO RECORD, FULL TEXT]**. Verbatim from the repo record:

   > "Unfortunately, our results show that LIDAR estimates of the ground surface are
   > positively biased even when the minimum-bin technique is used. This suggests that the
   > LIDAR laser pulse never reaches the ground surface within the vegetation communities we
   > studied."

   > "With 95% confidence, the DOGAMI bare-earth DEM elevation is between 2.0 cm and 3.1 cm
   > above the minimum-bin DEM elevation across the entire dataset (mean 2.5 cm, paired
   > two-sided t-test, p-value < 0.001)."

   That second quote is the only direct *delivered-product vs. low-order-statistic*
   comparison I found anywhere: the vendor bare-earth surface sits 2.0–3.1 cm above a
   minimum-bin surface built from the same cloud.
3. **Schmid, K.A., Hadley, B.C., & Wijekoon, N. (2011).** *Journal of Coastal Research*,
   27(6A), 116–132. doi:[10.2112/JCOASTRES-D-10-00188.1](https://doi.org/10.2112/JCOASTRES-D-10-00188.1)
   — **[REPO RECORD, ABSTRACT]**:

   > "By employing the minimum bin technique to the bare-earth classified LIDAR data, the
   > overall bias in the resultant surface was reduced by 12 cm, and the vertical accuracy was
   > improved by 8 cm when compared with the 'as-received' data."
4. **Hladik, C., & Alber, M. (2012).** Accuracy assessment and correction of a LIDAR-derived
   salt marsh digital elevation model. *Remote Sensing of Environment*, **121**, 224–235.
   doi:[10.1016/j.rse.2012.01.018](https://doi.org/10.1016/j.rse.2012.01.018)
   — **[METADATA ONLY — this session]**; Crossref confirms title/venue/volume/pages/date, but
   returns no abstract and ScienceDirect 403'd. Widely cited for species-class correction
   factors. **I saw specific error-reduction numbers only in a search-engine summary and am
   therefore not quoting them.**
5. **Buffington, K.J., Dugger, B.D., Thorne, K.M., & Takekawa, J.Y. (2016).** *Remote Sensing
   of Environment*, 186, 616–625. doi:[10.1016/j.rse.2016.09.020](https://doi.org/10.1016/j.rse.2016.09.020)
   — **[METADATA ONLY — this session]**; **[REPO RECORD, ABSTRACT]**:

   > "Using 17 study sites along the Pacific coast of the U.S., we achieved an average root
   > mean squared error (RMSE) of 0.072 m, with a 40–75% improvement in accuracy from the
   > lidar bare earth DEM. Results from our method compared favorably with results from three
   > other methods (minimum-bin gridding, mean error correction, and vegetation correction
   > factors)."

**What is missing from all of these:** none regresses the delivered-surface offset on a
*continuous* near-ground vegetation measure derived from the same point cloud, and none
reports how that offset varies with the number of returns in the cell. Both are things our
work does.

---

## 6. The scale/`n` confound in the low-order-statistic tradition

This is where our decimation experiment is methodologically distinct, so it is worth stating
precisely what the published designs do.

**Minimum-bin gridding (marsh tradition).** The estimator is the cell minimum. The tuning
knob is the **cell size**. Growing the cell simultaneously (i) grows `n`, which lowers the
expected minimum, and (ii) grows the spatial support, which lets the minimum be drawn from a
*different, lower place*. The literature reports the net, and names only the second effect.

Ewald (2013), **[REPO RECORD, FULL TEXT]**, verbatim:

> "DEM accuracy increased with cell size until an inflection point near 1.4 m as the influence
> of vegetation is mitigated by the minimum-bin gridding technique (Figure 2.3a). Low features
> within the landscape were captured by the gridding technique and degrade DEM performance
> after cell size enlarged beyond the optimum."

> "Even at the optimum cell size, the DEM is still positively biased when compared to known
> ground elevations. Mean LIDAR-GPS discrepancy remains positive until a cell size of 1.6 m is
> achieved. At cell sizes greater than 1.6 m, DEM are negatively biased as the minimum-bin
> method continues to capture and favor low features within the landscape."

> "The likelihood of upslope areas being assigned an elevation lower than the true ground
> elevation increases as the cell size is increased."

**Local-minima retrieval (forest tradition).** Same structure, larger scales.

Clark, M.L., Clark, D.B., & Roberts, D.A. (2004), *Remote Sensing of Environment*, 91(1),
68–89, doi:[10.1016/j.rse.2004.02.008](https://doi.org/10.1016/j.rse.2004.02.008)
— **[REPO RECORD, FULL TEXT]**, §2.5.1 verbatim:

> "the above local-minima scheme is analogous to selecting the lowest return in a square
> footprint of a specified scale (i.e., 5, 10 m, etc.)"

> "The scale with the lowest RMSE for both interpolation methods was found to be 20 m. **For
> this tropical landscape and lidar sampling density**, 20 m appears to be the near-optimum
> scale to identify ground returns with the local-minima approach"

*(a) What the source says:* the optimum window is 20 m, and that optimum is conditional on the
sampling density.
*(b) What I infer:* Clark et al. name the dependence and do not test it. A 20 m window at
their density holds some particular number of returns; nobody asks what happens if you hold
the window at 20 m and halve the returns. Our decimation does exactly that — fixed 5 m cell,
158 → 16 returns — which decouples `n` from support. The repo's own
`DTM_PERCENTILE_LITERATURE.md` gap #3 reached the same conclusion from the other direction:
"No study fixes cell size and tunes the percentile. All tuning precedent is in scale/window,
conflating the quantile with the support."

---

## 7. What appears novel, and what is well-trodden

### Well-trodden — do not claim these

- **Ground surface reads high under vegetation.** Established everywhere, forest and marsh,
  every sensor era. (`DTM_PERCENTILE_LITERATURE.md` synthesis table, 20+ studies.)
- **Low order statistics as a ground estimator.** Minimum-bin gridding (Schmid 2011; Ewald
  2013; Wang 2009) and local-minima retrieval (Clark 2004) are 20+ years old.
- **The window size of a minimum-type ground estimator must be tuned, and has a two-sided
  optimum.** Ewald 2013 is explicit and quantitative.
- **Point density degrades DTM accuracy, and ground filters degrade with it.** Hansen 2015;
  Cățeanu & Ciubotaru 2021; Raber 2007; Štular & Lozić 2020; Sithole & Vosselman 2004.
- **Density-dependent bias in the per-cell MAXIMUM, corrected by a probability model with an
  explicit `n`.** Roussel et al. 2017 owns this, completely, for the canopy top. Our argument
  is the same argument; the correct move is to cite it as the precedent and state that we
  invert it.
- **Delivered bare-earth surfaces floating above survey control in proportion to near-ground
  vegetation.** Hopkinson (ISPRS); Ewald 2013; Schmid 2011; Buffington 2016.

### Apparently novel — supported by the search nulls, which are nulls not proofs

1. **Running the order-statistics argument downward, onto the ground.** Roussel et al. (2017)
   derive `h̄(n)` for the maximum and explicitly disclaim the ground as their reference. No
   retrieved source does the minimum/low-quantile version. This is the strongest novelty
   claim available and it is a *framework transfer*, which should be stated as such rather
   than as a discovery.
2. **The estimator-dependence of the density effect, stated as a dichotomy and measured.**
   That a central statistic is `n`-invariant while a low order statistic is not is
   mathematically trivial, but I found no paper that *uses* it to explain why density matters
   for one DTM product and not another. Fradette et al. (2019) measured exactly this
   (DTM: no adjustment needed; CHM: needed) and did not frame it this way.
3. **Decoupling `n` from spatial support.** Every published tuning of a minimum-type ground
   estimator varies the window. Holding cell size fixed and decimating is, as far as this
   search reaches, not published for the ground.
4. **The paired measurement of the SAME cell's minimum and median under decimation.** Our
   "+20 mm at the lowest return, +0 mm at the median" is a within-cell contrast at fixed
   support. The closest published relative is Ewald's 2.0–3.1 cm vendor-vs-minimum-bin
   contrast, which is a *between-estimator* comparison at fixed density, not a
   *between-density* comparison at fixed estimator.

### Where we should be careful

- **Leitold et al. (2015) already report a large positive density-dependent ground bias**
  (+0.19 → +3.21 m). If we claim "density-versus-ground-bias has not been studied," that
  claim is wrong. The accurate claim is that it has been studied as a *filter-failure*
  phenomenon at metre magnitudes under closed tropical canopy, not as an *estimator/order-
  statistic* phenomenon at centimetre magnitudes.
- **A search null is not a citation-graph null.** The two closures worth buying with library
  access are (i) forward citations of Roussel et al. (2017) and (ii) forward citations of
  Schmid et al. (2011) / Ewald (2013), checked for any downward transfer of the probability
  model.

---

## 8. Not verified / not obtained in this session

| Item | Status | What NOT to claim |
|---|---|---|
| Cățeanu & Ciubotaru 2021, full text | MDPI returned 403 to WebFetch and to curl; only the abstract was obtained | Whether they report a signed mean error at all. The abstract does not. |
| Leitold et al. 2015, "incorrect classification of vegetation features as ground" | Appeared in a search-engine summary; not located in the fetched full text | Do not attribute this sentence to the paper without re-checking |
| Hladik & Alber 2012 error-reduction figures (0.10 → −0.01 m; RMSE 0.16 → 0.10 m) | Seen only in a search-engine summary | Do not quote these numbers |
| Sithole & Vosselman 2004 journal version | Abstract not retrievable this session (Semantic Scholar null); repo holds an abstract quote from an earlier session | Quote only via the repo record's tag, or re-obtain |
| Hopkinson et al. 2005 *Can. J. Remote Sens.* (the peer-reviewed version) | DOI/metadata confirmed; no abstract or body | All numbers above are from the ISPRS Archives conference version, which is **not** peer-reviewed |
| Any citation-graph sweep | Not run — no Scopus/WoS access | The order-statistics null is a search null, not a formal null |

## 9. Full reference list

Peer-reviewed journal articles unless marked otherwise. Read-level in bold.

- Cățeanu, M., & Ciubotaru, A. (2020). Accuracy of Ground Surface Interpolation from Airborne
  Laser Scanning (ALS) Data in Dense Forest Cover. *ISPRS International Journal of
  Geo-Information*, 9(4), 224. doi:[10.3390/ijgi9040224](https://doi.org/10.3390/ijgi9040224)
  — **[ABSTRACT — this session]**
- Cățeanu, M., & Ciubotaru, A. (2021). The Effect of LiDAR Sampling Density on DTM Accuracy
  for Areas with Heavy Forest Cover. *Forests*, 12(3), 265.
  doi:[10.3390/f12030265](https://doi.org/10.3390/f12030265) — **[ABSTRACT — this session]**
- Clark, M.L., Clark, D.B., & Roberts, D.A. (2004). Small-footprint lidar estimation of
  sub-canopy elevation and tree height in a tropical rain forest landscape. *Remote Sensing of
  Environment*, 91(1), 68–89. doi:[10.1016/j.rse.2004.02.008](https://doi.org/10.1016/j.rse.2004.02.008)
  — **[REPO RECORD, FULL TEXT]**
- Ewald, M.J. (2013). *Where's the Ground Surface? Elevation Bias in LIDAR-derived Digital
  Elevation Models Due to Dense Vegetation in Oregon Tidal Marshes.* MS thesis (**not
  peer-reviewed**), Oregon State University.
  <https://ir.library.oregonstate.edu/downloads/1n79h8198> — **[REPO RECORD, FULL TEXT]**
- Fradette, M.-S., Leboeuf, A., Riopel, M., & Bégin, J. (2019). Method to Reduce the Bias on
  Digital Terrain Model and Canopy Height Model from LiDAR Data. *Remote Sensing*, 11(7), 863.
  doi:[10.3390/rs11070863](https://doi.org/10.3390/rs11070863) — **[ABSTRACT — this session]**
- Hansen, E.H., Gobakken, T., & Næsset, E. (2015). Effects of Pulse Density on Digital Terrain
  Models and Canopy Metrics Using Airborne Laser Scanning in a Tropical Rainforest. *Remote
  Sensing*, 7(7), 8453–8468. doi:[10.3390/rs70708453](https://doi.org/10.3390/rs70708453)
  — **[ABSTRACT — this session]**
- Hladik, C., & Alber, M. (2012). Accuracy assessment and correction of a LIDAR-derived salt
  marsh digital elevation model. *Remote Sensing of Environment*, 121, 224–235.
  doi:[10.1016/j.rse.2012.01.018](https://doi.org/10.1016/j.rse.2012.01.018)
  — **[METADATA ONLY — this session]**
- Hopkinson, C., Chasmer, L., Zsigovics, G., Creed, I.F., Sitar, M., Treitz, P., & Maher, R.V.
  Errors in LiDAR ground elevation and wetland vegetation height estimates. *International
  Archives of Photogrammetry, Remote Sensing and Spatial Information Sciences*, XXXVI-8/W2,
  108–113. **Conference archives, not peer-reviewed.**
  <https://www.isprs.org/proceedings/xxxvi/8-w2/HOPKINSON.pdf> — **[FULL TEXT — this session]**
- Hopkinson, C., Chasmer, L., Sass, G.Z., Creed, I.F., Sitar, M., Kalbfleisch, W., & Treitz, P.
  (2005). Vegetation class dependent errors in lidar ground elevation and canopy height
  estimates in a boreal wetland environment. *Canadian Journal of Remote Sensing*, 31(2),
  191–206. doi:[10.5589/m05-007](https://doi.org/10.5589/m05-007)
  — **[METADATA ONLY — this session]**
- Leitold, V., Keller, M., Morton, D.C., Cook, B.D., & Shimabukuro, Y.E. (2015). Airborne
  lidar-based estimates of tropical forest structure in complex terrain: opportunities and
  trade-offs for REDD+. *Carbon Balance and Management*, 10:3.
  doi:[10.1186/s13021-015-0013-x](https://doi.org/10.1186/s13021-015-0013-x)
  — **[ABSTRACT + partial body — this session]**
- Raber, G.T., Jensen, J.R., Hodgson, M.E., Tullis, J.A., Davis, B.A., & Berglund, J. (2007).
  Impact of Lidar Nominal Post-spacing on DEM Accuracy and Flood Zone Delineation.
  *Photogrammetric Engineering & Remote Sensing*, 73(7), 793–804.
  doi:[10.14358/PERS.73.7.793](https://doi.org/10.14358/PERS.73.7.793)
  — **[ABSTRACT — this session]**
- Roussel, J.-R., Caspersen, J., Béland, M., Thomas, S., & Achim, A. (2017). Removing bias from
  LiDAR-based estimates of canopy height: Accounting for the effects of pulse density and
  footprint size. *Remote Sensing of Environment*, 198, 1–16.
  doi:[10.1016/j.rse.2017.05.032](https://doi.org/10.1016/j.rse.2017.05.032)
  — **[FULL TEXT — this session, open post-print]**
- Schmid, K.A., Hadley, B.C., & Wijekoon, N. (2011). Vertical Accuracy and Use of Topographic
  LIDAR Data in Coastal Marshes. *Journal of Coastal Research*, 27(6A), 116–132.
  doi:[10.2112/JCOASTRES-D-10-00188.1](https://doi.org/10.2112/JCOASTRES-D-10-00188.1)
  — **[REPO RECORD, ABSTRACT]**
- Simpson, J.E., Smith, T.E.L., & Wooster, M.J. (2017). Assessment of Errors Caused by Forest
  Vegetation Structure in Airborne LiDAR-Derived DTMs. *Remote Sensing*, 9(11), 1101.
  doi:[10.3390/rs9111101](https://doi.org/10.3390/rs9111101)
  — **[ABSTRACT — this session; REPO RECORD, FULL TEXT]**
- Sithole, G., & Vosselman, G. (2003). Comparison of Filtering Algorithms. *ISPRS Archives*
  XXXIV-3/W13, Dresden. **Conference proceedings, not peer-reviewed.**
  <https://www.isprs.org/proceedings/xxxiv/3-w13/papers/Sithole_ALSDD2003.pdf>
  — **[REPO RECORD, FULL TEXT]**
- Sithole, G., & Vosselman, G. (2004). Experimental comparison of filter algorithms for
  bare-Earth extraction from airborne laser scanning point clouds. *ISPRS Journal of
  Photogrammetry and Remote Sensing*, 59(1–2), 85–101.
  doi:[10.1016/j.isprsjprs.2004.05.004](https://doi.org/10.1016/j.isprsjprs.2004.05.004)
  — **[METADATA ONLY — this session; REPO RECORD, ABSTRACT]**
- Stereńczak, K., & Kozak, J. (2011). Evaluation of digital terrain models generated in forest
  conditions from airborne laser scanning data acquired in two seasons. *Scandinavian Journal
  of Forest Research*, 26(4), 374–384.
  doi:[10.1080/02827581.2011.570781](https://doi.org/10.1080/02827581.2011.570781)
  — **[ABSTRACT — this session]**
- Štular, B., & Lozić, E. (2020). ALS ground-filter comparison across four sites.
  — **[REPO RECORD, FULL TEXT]** in `analysis/ridgelines/CSF_PARAMETER_LITERATURE.md`; full
  citation is in that note.
- Wang, C., Menenti, M., Stoll, M.P., Feola, A., Belluco, E., & Marani, M. (2009). Separation
  of Ground and Low Vegetation Signatures in LiDAR Measurements of Salt-Marsh Environments.
  *IEEE Transactions on Geoscience and Remote Sensing*, 47(7), 2014–2023.
  doi:[10.1109/TGRS.2008.2010490](https://doi.org/10.1109/TGRS.2008.2010490)
  — **[REPO RECORD, ABSTRACT]**
- Buffington, K.J., Dugger, B.D., Thorne, K.M., & Takekawa, J.Y. (2016). Statistical correction
  of lidar-derived digital elevation models with multispectral airborne imagery in tidal
  marshes. *Remote Sensing of Environment*, 186, 616–625.
  doi:[10.1016/j.rse.2016.09.020](https://doi.org/10.1016/j.rse.2016.09.020)
  — **[METADATA ONLY — this session; REPO RECORD, ABSTRACT]**

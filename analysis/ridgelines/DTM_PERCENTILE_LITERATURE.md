# Choice and calibration of the per-cell vertical statistic in bare-earth lidar DTMs

**A literature review for the 2008 MN DNR vs. 2021 USGS 3DEP change-detection study
(forested southeastern Minnesota bluffland).**

Compiled 2026-08-26. Not committed to git.

---

## How to read this document

Every finding below carries a verification tag. Nothing is stated from a title alone,
and nothing is stated from recall.

| Tag | Meaning |
|---|---|
| **[FULL TEXT]** | I downloaded the PDF and read the passage. Quotes are transcribed from the PDF text layer. |
| **[ABSTRACT]** | I retrieved the publisher's abstract verbatim (via OpenAlex / Semantic Scholar / publisher page). Quotes are from that abstract. Body-text details are **not** verified. |
| **[SECONDARY]** | The number is quoted from a *different* paper's description of the source. The original was not opened. |
| **[NOT VERIFIED]** | Cited for completeness; I could not obtain the source text. No numbers asserted. |

Study-setting notes are given for every source, because most of the directly relevant
work is in **tidal marsh**, not forest. Marsh findings do not transfer to Minnesota
bluffland without argument (see "What this means for us").

---

## 1. Percentile / quantile choice when gridding classified ground returns

### 1.1 The short answer

**There is no published guidance on choosing a vertical percentile per cell for a
forested bare-earth DTM.** The peer-reviewed literature on "how do I turn ground
points into a raster" is almost entirely about *interpolation algorithms* (TIN, IDW,
kriging, splines, natural neighbour), not about *which order statistic of the points
inside a cell* to use. The one named percentile-like technique with a peer-reviewed
pedigree — **minimum-bin gridding** — is a coastal-marsh method, and its primary
sources are all marsh papers.

### 1.2 Minimum-bin gridding: primary sources

**Schmid, K.A., Hadley, B.C., & Wijekoon, N. (2011). Vertical Accuracy and Use of
Topographic LIDAR Data in Coastal Marshes. *Journal of Coastal Research*, 27(6A),
116–132.** doi:[10.2112/JCOASTRES-D-10-00188.1](https://doi.org/10.2112/JCOASTRES-D-10-00188.1)
— peer-reviewed journal article. **[ABSTRACT]**

This is the clearest primary definition I found. Verbatim from the abstract:

> "Custom digital elevation model (DEM) generation techniques and point classification
> processes can be used to improve estimates of ground elevations in coastal marshes.
> The simplest of these methods is minimum bin gridding, which extracts the lowest
> elevation value included within a user-specified search window and assigns that value
> to the appropriate DEM grid cell."

and the headline result:

> "By employing the minimum bin technique to the bare-earth classified LIDAR data, the
> overall bias in the resultant surface was reduced by 12 cm, and the vertical accuracy
> was improved by 8 cm when compared with the 'as-received' data."

Note the explicit cost, also verbatim:

> "Despite lowering the spatial resolution of the DEM, the application of these
> techniques significantly improves the vertical accuracy of the LIDAR-derived
> bare-earth surfaces."

*Setting:* US south-Atlantic coastal marsh (*Spartina alterniflora*, *Juncus
roemerianus*); flat; discrete-return lidar of the mid-2000s. Point density not verified.

---

**Wang, C., Menenti, M., Stoll, M.P., Feola, A., Belluco, E., & Marani, M. (2009).
Separation of Ground and Low Vegetation Signatures in LiDAR Measurements of Salt-Marsh
Environments. *IEEE Transactions on Geoscience and Remote Sensing*, 47(7), 2014–2023.**
doi:[10.1109/TGRS.2008.2010490](https://doi.org/10.1109/TGRS.2008.2010490)
— peer-reviewed journal article. **[ABSTRACT]**

The second antecedent Ewald names. Verbatim:

> "In this paper, we introduce reliable methods to remove random and systematic errors
> and to register raw data, as well as a new procedure, to determine the optimal filter
> window size to separate ground and canopy returns. A limited amount of field
> observations is used to determine the size of the filtering window which produces the
> minimally biased estimates of the digital terrain model (DTM)."

> "We apply this procedure to a study marsh within the Venice Lagoon, Italy, and obtain
> a high-accuracy DTM. The error (z_LiDAR − z_field) is 2.2 cm, with a standard
> deviation of 6.4 cm."

*Setting:* Venice Lagoon salt marsh, Italy; flat; short marsh vegetation where
"the characteristic short vegetation does not give rise to detectable differences
between first and last LiDAR returns" (verbatim). Deliberately **not** a forest method.

---

**Ewald, M.J. (2013). *Where's the Ground Surface? Elevation Bias in LIDAR-derived
Digital Elevation Models Due to Dense Vegetation in Oregon Tidal Marshes.* MS thesis,
Oregon State University.** <https://ir.library.oregonstate.edu/downloads/1n79h8198>
— **MS thesis, not peer-reviewed**, but it is the fullest published *evaluation* of the
technique. **[FULL TEXT]** (all quotes below transcribed from the PDF)

Definition as Ewald states it (p. 15):

> "Schmid et al. (2011) and Wang et al. (2009) employ a minimum-bin gridding technique
> to identify the optimum search radius for ground LIDAR returns and interpolate a
> raster DEM from the LIDAR point cloud. This gridding and interpolation technique
> selects the lowest LIDAR return within a specified search radius. As the search radius
> and cell size is increased, the probability of capturing a true ground also increases
> as more candidate LIDAR returns are considered. The minimum-bin technique is
> attractive because it is easy to implement and validate in the field."

The evaluation is the valuable part, because it documents the **two-sided failure**:

> "DEM accuracy increased with cell size until an inflection point near 1.4 m as the
> influence of vegetation is mitigated by the minimum-bin gridding technique
> (Figure 2.3a). Low features within the landscape were captured by the gridding
> technique and degrade DEM performance after cell size enlarged beyond the optimum."

> "Even at the optimum cell size, the DEM is still positively biased when compared to
> known ground elevations. Mean LIDAR-GPS discrepancy remains positive until a cell size
> of 1.6 m is achieved. At cell sizes greater than 1.6 m, DEM are negatively biased as
> the minimum-bin method continues to capture and favor low features within the
> landscape."

> "Minimum-bin LIDAR-derived DEM elevations underpredict (typically 70 cm to 90 cm
> below) the measured elevation along channels and the man-made dike that form the
> southern edge of the site adjacent to the Coquille River. These features are
> characterized by high ground, a moderate to steep slope, and low ground over a short
> horizontal distance. For example, the minimum-bin filter is likely to select a LIDAR
> return from an adjacent low riverbank rather than the surveyed wetland surface. The
> likelihood of upslope areas being assigned an elevation lower than the true ground
> elevation increases as the cell size is increased."

**This is the single most transferable result in the marsh literature for us.** The
minimum statistic reaches downhill. On a bluffland with 20–35° slopes, that failure mode
is not an edge case — it is the dominant one, and it is *exactly* the low-percentile
downhill bias already identified in this project as the source of DoD banding.

Crucially, Ewald also finds min-bin barely beats the vendor bare-earth surface:

> "Throughout the vegetation types we evaluated, the minimum-bin DEM performs slightly
> better than the DOGAMI bare-earth DEM. With 95% confidence, the DOGAMI bare-earth DEM
> elevation is between 2.0 cm and 3.1 cm above the minimum-bin DEM elevation across the
> entire dataset (mean 2.5 cm, paired two-sided t-test, p-value < 0.001)."

> "Unfortunately, our results show that LIDAR estimates of the ground surface are
> positively biased even when the minimum-bin technique is used. This suggests that the
> LIDAR laser pulse never reaches the ground surface within the vegetation communities
> we studied."

*Setting:* Ni-les'tun Unit, Bandon Marsh NWR and five other Oregon estuaries; 174 ha;
tidal marsh and diked pasture; flat except channels and dikes; > 13,000 RTK GPS points.

---

**Buffington, K.J., Dugger, B.D., Thorne, K.M., & Takekawa, J.Y. (2016). Statistical
correction of lidar-derived digital elevation models with multispectral airborne imagery
in tidal marshes. *Remote Sensing of Environment*, 186, 616–625.**
doi:[10.1016/j.rse.2016.09.020](https://doi.org/10.1016/j.rse.2016.09.020)
— peer-reviewed journal article. **[ABSTRACT]** (text obtained from the NOAA NCCOS
record page, which reproduces the publisher abstract; I could not open ScienceDirect,
which returned HTTP 403)

Verbatim from that page:

> "Using 17 study sites along the Pacific coast of the U.S., we achieved an average root
> mean squared error (RMSE) of 0.072 m, with a 40–75% improvement in accuracy from the
> lidar bare earth DEM. Results from our method compared favorably with results from
> three other methods (minimum-bin gridding, mean error correction, and vegetation
> correction factors)."

This is the most useful *benchmarking* statement available: it establishes that
minimum-bin gridding is one of the three recognised correction families, and that a
fitted statistical model beat all three. A widely-repeated figure of "118 points
necessary to calibrate a site-specific correction model" appeared in a search summary
but **I could not verify it against the source text — treat as [NOT VERIFIED]**.

*Setting:* 17 Pacific-coast tidal marshes; flat; RTK-GPS reference; NAIP imagery.

### 1.3 The principled continuous alternative: asymmetric (skew) robust interpolation

The forestry/photogrammetry tradition solved the same problem *without* a percentile,
by making the residual weight function asymmetric so the fitted surface is pulled toward
the low returns. This is the closest thing to a theory of "why not just take a low
quantile."

**Kraus, K. & Pfeifer, N. (1998). Determination of terrain models in wooded areas with
airborne laser scanner data. *ISPRS Journal of Photogrammetry and Remote Sensing*,
53(4), 193–203.** doi:[10.1016/S0924-2716(98)00009-4](https://doi.org/10.1016/S0924-2716(98)00009-4)
— peer-reviewed journal article. **[NOT VERIFIED]** — I could not obtain the full text
or a verbatim abstract. I therefore assert nothing numerical from it directly, and
instead quote the companion conference paper below, which restates the method and cites
Kraus & Pfeifer (1998) for the numbers.

**Kraus, K. & Rieger, W. (1999). Processing of laser scanning data for wooded areas.
In D. Fritsch & R. Spiller (Eds.), *Photogrammetric Week '99*, pp. 221–231. Wichmann,
Heidelberg.** <https://phowo.ifp.uni-stuttgart.de/publications/phowo99/kraus.pdf>
— **conference proceedings, not peer-reviewed journal.** **[FULL TEXT]**

The method, verbatim:

> "This algorithm estimates the skewness of the error distribution of the laser scanner
> data in forests and assigns small weights to those points that show large positive
> errors during the interpolation with filtering. The process results in a
> classification of the laser points in terrain and off-terrain (mainly vegetation)
> points."

And the accuracy law — directly relevant to our slope-dependent residual:

> "σH[cm] = ± (18 + 120·tanα)"
>
> "Equation (1) is valid for a ground penetration rate of the laser signal of at least
> 25 % and a good mixture of vegetation and ground points for the whole region."

> "The elimination of inherent systematic errors allows a significant improvement of the
> accuracy of the laser DTM particularly in flat terrain. The constant value of ± 18 cm
> in equation (1) can be reduced down to ± 10 cm (K. Kraus, N. Pfeifer, 1998)."

*Setting:* Austrian wooded areas (Danube riparian forest, Vienna Woods); mixed
deciduous/conifer; late-1990s fixed-wing scanners, ~3 m point spacing. Old sensors, but
the *form* of the accuracy law (a slope term ~10× the flat-ground term) is the published
statement closest to what this project has measured.

### 1.4 What the field actually compares (and therefore what is missing)

**Boreggio, M., Bernard, M., & Gregoretti, C. (2018). Evaluating the Differences of
Gridding Techniques for Digital Elevation Models Generation and Their Influence on the
Modeling of Stony Debris Flows Routing. *Frontiers in Earth Science*, 6, 89.**
doi:[10.3389/feart.2018.00089](https://doi.org/10.3389/feart.2018.00089)
— peer-reviewed journal article. **[ABSTRACT + body via publisher page]**

Despite the title promising "gridding techniques", the compared set is
"linear triangulation, natural neighbor, nearest neighbor, Inverse Distance to a Power,
ANUDEM, Radial Basis Functions, and ordinary kriging" (verbatim) — i.e. **interpolators
only, no per-cell order statistic.** *Setting:* Rovina di Cancia debris-flow basin,
Italian Dolomites; steep, largely unvegetated.

**Montealegre, A.L., Lamelas, M.T., & de la Riva, J. (2015). Interpolation Routines
Assessment in ALS-Derived Digital Elevation Models for Forestry Applications.
*Remote Sensing*, 7(7), 8631–8654.** doi:[10.3390/rs70708631](https://doi.org/10.3390/rs70708631)
— peer-reviewed journal article. **[ABSTRACT]**

Same story in a forest setting. Verbatim:

> "In this study, six interpolation routines were tested over a range of land cover and
> terrain roughness in order to generate a collection of DEMs with spatial resolution of
> 1 and 2 m."

> "The Triangulated Irregular Network (TIN) to raster interpolation method produced the
> best result in the validation process with the training data set while the Inverse
> Distance Weighted (IDW) routine was the best in the validation with GPS (RMSE of 2.68
> cm and RMSE of 37.10 cm, respectively)."

*Setting:* Mediterranean forest, Spain. Again: **interpolators, not statistics.** Note
the ~14× gap between self-validation and GPS validation — a caution about validating a
gridding choice against the same point cloud that produced it.

**USGS Lidar Base Specification (agency standard, not peer-reviewed).**
<https://www.usgs.gov/ngp-standards-and-specifications/lidar-base-specification-digital-elevation-model-surface>
**[partially verified]** — I read the "Digital Elevation Model Surface Treatments" page
of *Lidar Base Specification 2025 rev. A*. It specifies which points are **excluded**
(e.g. "Bare-earth lidar points (serving as mass points) that are in close proximity to
any breakline shall be classified as Ignored Ground (class 20) and shall be excluded
from the DEM generation process", verbatim). **The page I read does not state an
interpolation method and does not specify any per-cell statistic.** I did *not* verify
the separate "Data Processing and Handling Requirements" page, so I make no claim about
whether TIN is mandated elsewhere in the spec. What is safe to say: the governing
specification for 3DEP treats the DEM as an *interpolation over all valid ground points*,
not as a per-cell order statistic — so a percentile-gridded product is a departure from
the delivered-product convention and must be built identically for both epochs.

### 1.5 Does anyone use a percentile *above* the median?

**Not for ground elevation. I found no case.** I searched explicitly for 75th / 90th /
upper-quartile ground gridding and found only:

- **Canopy and crop *surface* models**, where high percentiles (Zp90, Zp95) replace the
  maximum because they are more robust to point-density variation. That is the *opposite*
  problem — the upper envelope, not the lower — and provides no argument for raising the
  ground statistic.
- **Medeiros et al. (2015)** (below), which applies a *quartile-based elevation
  adjustment*, not a quartile of the returns in a cell. The nearest published relative,
  but not the same thing.

The one *mechanism* in the literature that would justify a statistic above the minimum
is Ewald's downhill-capture failure — the minimum is biased **low** on sloping ground
adjacent to lower ground. That argues for moving *up* from the minimum, and Ewald shows
the optimum overshoots into negative bias. Nobody has taken the next step and asked which
percentile is optimal at fixed cell size.

---

## 2. Empirical vertical bias of lidar DTMs under vegetation

**Sign convention used below: positive = lidar DTM reads ABOVE surveyed ground.**
Where a source's own convention is ambiguous I say so.

### 2.1 Temperate conifer, total-station reference — the classic benchmark

**Reutebuch, S.E., McGaughey, R.J., Andersen, H.-E., & Carson, W.W. (2003). Accuracy of a
high-resolution lidar terrain model under a conifer forest canopy. *Canadian Journal of
Remote Sensing*, 29(5), 527–535.** doi:[10.5589/m03-022](https://doi.org/10.5589/m03-022)
— peer-reviewed journal article. **[FULL TEXT]**

Abstract, verbatim:

> "Conventional ground survey methods were used to collect coordinates and near-ground
> vegetation heights at 347 ground checkpoints distributed under a range of canopy
> covers. These points were used to check the DTM accuracy. The mean DTM error was
> 0.22 ± 0.24 m (mean ± SD). DTM elevation errors for four tree canopy cover classes
> were: clearcut 0.16 ± 0.23 m, heavily thinned 0.18 ± 0.14 m, lightly thinned
> 0.18 ± 0.18 m, and uncut 0.31 ± 0.29 m. These DTM errors show a slight increase with
> canopy density but the differences are strikingly small."

Table 2 (transcribed from the PDF), lidar DTM minus surveyed checkpoint, metres:

| Canopy class | Mean | SD | Min | Max | n |
|---|---|---|---|---|---|
| Clearcut | 0.16 | 0.23 | −0.48 | 0.61 | 38 |
| Heavy thinned | 0.18 | 0.14 | −0.11 | 0.41 | 21 |
| Lightly thinned | 0.18 | 0.18 | −0.63 | 0.69 | 147 |
| Uncut | 0.31 | 0.29 | −0.60 | 1.31 | 120 |

Grouped analysis (Table 3, transcribed): no near-ground vegetation, mean 0.15 m (n = 132);
any vegetation within 6 m of ground, mean 0.26 m (n = 212); slope < 18 %, mean 0.21 m
(n = 174); slope ≥ 18 %, mean 0.22 m with SD rising from 0.20 to 0.28 m (n = 173).

The **empirical offset removal** — see §3 — verbatim:

> "The observed error in the clearcut area (0.16 m) is very similar to the lidar
> manufacturer's stated accuracy of ±0.15 m (Baltsavias, 1999). If one assumes that this
> error in the open, bare-ground clearcut area is the system bias and adjusts the
> individual checkpoint errors to remove this bias, then 69% of the observed checkpoint
> errors are within ±0.22 m (the observed SD of the DTM grid error)."

*Setting:* 500 ha, western Washington State; slopes 0–45°; 70-year-old conifer with
clearcuts and thinnings; Saab TopEye helicopter lidar, spring 1999; ~4.22 raw returns/m²,
~0.58 filtered ground points/m²; 1.52 m DTM grid. **Conifer, not deciduous** — the
leaf-state axis is absent by construction.

### 2.2 Mixed land cover, deciduous vs. evergreen, GPS + total station

**Hodgson, M.E. & Bresnahan, P. (2004). Accuracy of Airborne Lidar-Derived Elevation:
Empirical Assessment and Error Budget. *Photogrammetric Engineering and Remote Sensing*,
70(3), 331–339.** doi:[10.14358/PERS.70.3.331](https://doi.org/10.14358/PERS.70.3.331)
— peer-reviewed journal article. **[ABSTRACT]**

Verbatim:

> "The variability of vertical accuracy was evaluated for six land-cover categories.
> Root-mean-squared error (RMSE) values ranged from a low of 17 to 19 cm (pavement, low
> grass, and evergreen forests) to a high of 26 cm (deciduous forests)."

> "Observed elevation error in steeper slopes (e.g., 25°) was estimated to be twice as
> large as those on low slopes (e.g., 1.5°)."

Note the method, which is unusually clean and worth emulating:

> "Rather than using an interpolation approach for gathering observed elevations at
> reference points, the x-y coordinates of lidar points were located in the field and
> these elevations were surveyed."

*Setting:* Richland County, South Carolina; Optech ALTM 1210, 1207 m AGL, 2 m nominal
posting. Mixed temperate; includes deciduous forest. **Reports RMSE by cover, not signed
bias by cover** — so it constrains scatter, not the sign of our offset.

### 2.3 Deciduous (aspen) parkland — signed bias by cover, and a nadir result

**Su, J. & Bork, E. (2006). Influence of Vegetation, Slope, and Lidar Sampling Angle on
DEM Accuracy. *Photogrammetric Engineering and Remote Sensing*, 72(11), 1265–1274.**
doi:[10.14358/PERS.72.11.1265](https://doi.org/10.14358/PERS.72.11.1265)
— peer-reviewed journal article. **[FULL TEXT]**

Abstract, verbatim:

> "Across the study area, overall signed error and RMSE were 0.02 m and 0.59 m,
> respectively. Signed errors indicated elevations were over-estimated in forest but
> under-estimated within meadow habitats. Increasing slope gradient increased vertical
> absolute errors and RMSE. In contrast, lidar sampling angle had little impact on
> measured error."

The numbers, verbatim from the Results:

> "In identifying a tendency to under- or over-estimate elevations, the mean signed
> errors in Figure 3c show that elevations within aspen forest were over-estimated
> (0.20 m) while those in lowland meadows were under-estimated (0.22 m). Examination of
> the eight detailed classes of vegetation revealed a strong tendency to over-estimate
> elevations in both closed and semi-open aspen forest (Figure 3d)."

> "Finally, RMSE values indicated the lidar-derived DEM accuracy generally decreased as
> slope gradient increased: the RMSE at slopes over 10° was twice that found when slopes
> were less than 2°. This finding was similar to Hodgson and Bresnahan (2004), who
> observed errors on slopes of 25° to be twice that found on relatively flat areas."

And — directly relevant to this project's near-nadir finding, so quoted in full including
the authors' own scepticism:

> "Signed errors and RMSEs were generally greater when lidar data were collected close to
> nadir (less than 3°) relative to those sampled in angle classes further away from the
> central flightline (Figure 3b). However, this pattern may be attributed to the presence
> of extreme errors. The mean top five signed errors near nadir were 23 times larger than
> their corresponding overall signed errors. Moreover, extreme errors were accompanied by
> high (10°) slope gradients, which may also have contributed to the observed elevation
> errors."

*Transcription note.* The above is quoted exactly as the PDF text layer renders it. Two
checks: (i) the string "(10°)" almost certainly reads "(>10°)" in print — the extracted
text of this paper contains **zero** `<` or `>` characters across all 10 pages, despite
repeated slope-threshold discussion, so those glyphs are systematically lost; (ii) "23
times" is **not** a mangled "2–3": en-dashes survive extraction elsewhere in the file
(e.g. the reference page range "3482–3486"), so "23 times larger" appears to be literal,
and it is physically plausible because the comparison is top-five extremes against an
overall mean signed error near zero (+0.02 m area-wide). **Verify both against the
printed page before quoting in the manuscript.** I initially reconstructed these as
"2–3 times" and ">10°"; the glyph audit overturned the first reconstruction.

*Setting:* Kinsella Research Station, Aspen Parkland, Alberta; 2,700 ha of knob-and-kettle
terrain with 5–10 m relief; **deciduous aspen forest**, shrubland, grassland, meadow;
last-return lidar (1998-era), IDW interpolation at 1.5 m; 256 reference plots by total
station + DGPS against 27 interconnected benchmarks. Of everything in this review, this
is the **closest analogue to our forest/open contrast** — deciduous, mixed cover, sloping,
first-generation sensor.

### 2.4 Temperate deciduous broadleaf — the strongest modern study

**Simpson, J.E., Smith, T.E.L., & Wooster, M.J. (2017). Assessment of Errors Caused by
Forest Vegetation Structure in Airborne LiDAR-Derived DTMs. *Remote Sensing*, 9(11),
1101.** doi:[10.3390/rs9111101](https://doi.org/10.3390/rs9111101)
— peer-reviewed journal article, open access. **[FULL TEXT]**

Abstract, verbatim:

> "Here, we use ground survey equipment to assess digital terrain model (DTM) accuracy in
> a deciduous broadleaf forest, during both leaf-on and leaf-off conditions. Using the
> leaf-on LiDAR dataset we quantitatively assess vertical vegetation structure, and use
> this as a categorical explanatory variable for DTM accuracy. In the presence of leaf-on
> vegetation, DTM accuracy is severely reduced, with low-stature undergrowth vegetation
> (such as ferns) causing the greatest errors (RMSE > 1 m). Errors are lower under
> leaf-off conditions (RMSE = 0.22 m)."

Conclusion, verbatim:

> "Results demonstrate that leaf-on vegetation causes greater DTM error (RMSE = 0.83 m)
> than leaf-off vegetation (RMSE = 0.22) across all vegetation categories. Furthermore,
> DTM accuracy is not affected by all vegetation structures equally; with dense
> understory vegetation such as ferns and brambles causing the greatest positive DTM
> errors. Grassland vegetation yields the most accurate DTMs."

**Sign-convention caveat, flagged because it matters.** The Methods say the reference was
subtracted from the lidar surface — "The control DTM was subtracted from the ALS-derived
DTM (DTM_ALS) to produce difference rasters" — which makes positive = lidar high. But the
caption of Figure 11 reads "Digital Terrain Model (DTM) error (DTM_TS − DTM_ALS)", the
opposite. **The paper is internally inconsistent.** The physical direction the authors
intend is nonetheless unambiguous from the Discussion:

> "The results of the present study suggest that without an adequate ground control
> scheme, such applications may be prone to significant positive biases, especially in
> areas of leaf-on, open canopy forest with large amounts of ground cover."

i.e. **lidar reads above true ground under leaf-on vegetation.**

Also verbatim, and important for us — slope was *not* a driver here:

> "The median slope within the plot was 5.7°. The relationship between vertical DTM
> residuals and slope at 1m resolution was examined using a non-parametric GAM. Slope has
> no meaningful effect on residuals, with slope explaining 0.25% of the deviance, and a
> poor goodness of fit (adjusted R2 = 0.002, Figure 10)."

That null is only informative up to ~6° — it does not speak to bluffland slopes.

And the authors' own limitation on transferability, verbatim:

> "The results presented here show how DTM accuracy is relatively affected by vegetation
> structure, and as such they cannot be applied absolutely to other forest environments
> (i.e., these are not correction factors)."

*Setting:* < 1 ha plot in UK deciduous broadleaf woodland (alder, field maple, hazel);
low relief, median slope 5.7°; NERC ARSF leaf-on 24 June 2014 and leaf-off 9 March 2009;
~3–5 returns m⁻²; DTMs by IDW at 1 m; reference DTM from total-station survey (n = 657).
Small plot, low slope, only one site — but the only study that isolates *vertical
vegetation structure* as the explanatory variable in temperate deciduous forest.

### 2.5 A compiled cross-biome table (secondary, but useful)

Simpson et al. (2017) Table 1, transcribed **[FULL TEXT — but the numbers are their
transcription of other papers, so [SECONDARY] with respect to the originals]**:

| Biome | Vertical accuracy (m) | Metric | Source as cited |
|---|---|---|---|
| Old growth tropical forest | 1.95 | RMSE | Clark et al. 2004 |
| Secondary tropical forest | 1.44 | RMSE | Clark et al. 2004 |
| Steep Mediterranean shrubland | 0.13–0.41 | RMSE | Estornell et al. 2011 |
| Temperate conifer | 0.21 | RMSE | Bao et al. 2008 |
| Temperate conifer | −0.05 / 0.12 | Mean/SD | Hyyppä et al. 2005 |
| Temperate conifer | 0.31 / 0.29 | Mean/SD | Reutebuch et al. 2003 |
| Temperate conifer | 0.59 | RMSE | Su & Bork 2006 |
| Temperate conifer | 0.24 | RMSE | Tinkham et al. 2012 |
| **Temperate deciduous and conifer** | **1.22** | RMSE | Hodgson et al. 2003 |
| Temperate grass | 0.37 | RMSE | Hodgson et al. 2003 |
| Temperate mixed | 0.38 | N/A | Wasser et al. 2013 |
| Temperate pine | 0.45 | RMSE | Hodgson et al. 2003 |
| Temperate shrub | 1.53 | RMSE | Hodgson et al. 2003 |
| Tropical forest | 1.8 | Mean | Hansen et al. 2015 |

The "1.22 m RMSE, temperate deciduous and conifer" entry is attributed to **Hodgson,
M.E., Jensen, J.R., Schmidt, L., Schill, S., & Davis, B. (2003). An evaluation of LIDAR-
and IFSAR-derived digital elevation models in leaf-on conditions with USGS Level 1 and
Level 2 DEMs. *Remote Sensing of Environment*, 84(2), 295–308.**
doi:[10.1016/S0034-4257(02)00114-1](https://doi.org/10.1016/S0034-4257(02)00114-1)
— peer-reviewed. **[NOT VERIFIED]** — the abstract is withheld by the publisher on every
aggregator I tried (Crossref, OpenAlex, Semantic Scholar all return null; ScienceDirect
403). **I could not open this paper. Do not cite the 1.22 m figure to it in the manuscript
without obtaining the original.** Its title does establish a leaf-on acquisition, and
Simpson et al. (2017) attribute the value to it, but that is a chain of two unverified
steps.

Also note the two-order-of-magnitude spread in this table. The published "vegetation bias"
in forest ranges from −0.05 m to +1.95 m depending on biome, sensor and undergrowth. There
is no canonical number.

### 2.6 A dissenting result worth keeping

**Tinkham, W.T., Smith, A.M.S., Hoffman, C., Hudak, A.T., Falkowski, M.J., Swanson, M.E.,
& Gessler, P.E. (2012). Investigating the influence of LiDAR ground surface errors on the
utility of derived forest inventories. *Canadian Journal of Forest Research*, 42(3),
413–422.** doi:[10.1139/x11-193](https://doi.org/10.1139/x11-193)
— peer-reviewed journal article. **[ABSTRACT]**

Verbatim:

> "This study combines LiDAR DEMs and 54 ground survey plots to investigate how surface
> morphology and vegetation structure influence DEM errors. The study further compared
> two LiDAR classification algorithms and found no significant difference in their
> performance. Vegetation structure was found to have no influence, whereas increased
> variability in the vertical error was observed on slopes exceeding 30°, illustrating
> that these algorithms are not limited by high-biomass western coniferous forests, but
> that slope and sensor accuracy both play important roles."

*Setting:* western US high-biomass **coniferous** forest. **Conifer, evergreen — no leaf
state to mismatch.** This is the cleanest published statement that, in conifer, slope
rather than vegetation dominates. It is consistent with, not contrary to, the deciduous
results: the leaf-state axis simply does not exist in evergreen stands.

### 2.7 Marsh numbers (for completeness; low transferability)

- **Ewald (2013) [FULL TEXT]:** "Within wetland vegetation communities, my results
  suggest that LIDAR estimates of the ground surface in tidal wetlands are typically
  10 cm to 30 cm above GPS measurements. Plant associations dominated by *Carex obnupta*
  and *Carex lyngbyei* exhibited the largest discrepancy between LIDAR and GPS
  measurements (mean discrepancies 36.6 cm and 48.8 cm respectively)." Open land cover
  showed "fundamental vertical accuracy (FVA) of the LIDAR datasets was 4.5 cm root mean
  square error (RMSE) and had no consistent positive or negative bias in open landcover."
- **Hladik & Alber (2012) [ABSTRACT]:** "We found that DEM mean vertical errors for
  different cover classes ranged from 0.03 to 0.25 m in comparison to the RTK ground
  truth data, with the larger offsets for taller vegetation."

Both are dense, uniform, sub-metre herbaceous canopies on flat ground with essentially
zero ground-return probability. Minnesota bluffland forest has an open winter canopy, high
ground-return density, and 20–35° slopes. **These magnitudes should not be carried across.**

---

## 3. Tuning a percentile or an offset against ground-truth checkpoints

Yes — this is done routinely, but **almost always on a parameter other than the percentile
itself**, and almost always in marsh.

### 3.1 Tuning the aggregation *scale* (the closest published analogue)

- **Wang et al. (2009) [ABSTRACT]:** "A limited amount of field observations is used to
  determine the size of the filtering window which produces the minimally biased estimates
  of the digital terrain model (DTM)." — the filter window is explicitly a
  calibrated parameter fitted to minimise bias.
- **Schmid et al. (2011)**, as reported by **Ewald (2013) [FULL TEXT, SECONDARY for the
  Schmid numbers]:** "Schmid et al. (2011) considered cell sizes of between 2.0 m and
  10.0 m, selecting an optimal cell size of 4.0 m in *Spartina alterniflora* and a cell
  size of 10.0 m in *Juncus roemerianus* by minimizing the Root Mean Square Error (RMSE)
  on 280 survey-grade GPS measurements in South Carolina. Wang et al. (2009) used 240
  survey-grade GPS measurements and cell sizes between 0.5 m and 6.5 m, finding that an
  optimum cell size of 3.5 m by minimizing the overall RMSE."
- **Ewald (2013) [FULL TEXT]:** "Minimizing MAE yielded (1 in Table 2.1) an optimum cell
  size of 1.4 m (1.92 m²) for this dataset. If RMSE (2 in Table 2.1) was used to define
  the optimum DEM cell size instead of MAE, the optimum cell size was 1.2 m (1.44 m²)."
  Note that the *choice of loss function moved the optimum by 0.2 m* and the bias by
  2.9 cm — a caution for anyone tuning a statistic against checkpoints.

**The structural point:** in all three, the minimum is held fixed and the *neighbourhood*
is tuned. Enlarging the neighbourhood under a minimum operator is monotonically equivalent
to lowering the effective quantile of a fixed neighbourhood. So the literature has been
tuning a quantile all along — just parameterised as cell size, which entangles it with
resolution. **Nobody has decoupled the two by fixing cell size and tuning the percentile
directly.** That decoupling is, as far as I can establish, unpublished.

### 3.2 Tuning an additive offset by cover class

- **Hladik, C. & Alber, M. (2012). Accuracy assessment and correction of a LIDAR-derived
  salt marsh digital elevation model. *Remote Sensing of Environment*, 121, 224–235.**
  doi:[10.1016/j.rse.2012.01.018](https://doi.org/10.1016/j.rse.2012.01.018)
  — peer-reviewed. **[ABSTRACT]** Verbatim: "We developed species-specific correction
  factors for ten cover classes and used these correction factors to modify the
  LIDAR-derived DEM in four areas of the study domain where vegetation boundaries were
  mapped directly in the field. Application of the derived correction factors greatly
  improved the accuracy of the LIDAR-derived DEM within these areas, reducing the overall
  mean DEM error from 0.10 ± 0.12 (SD) to − 0.01 ± 0.09 m (SD), and the Root Mean Square
  Error from 0.16 m to 0.10 m."
  *Setting:* Sapelo Island, GA salt marsh; Optech Gemini ALTM at 125 kHz; RTK GPS reference.

- **Medeiros, S., Hagen, S., Weishampel, J., & Angelo, J. (2015). Adjusting
  Lidar-Derived Digital Terrain Models in Coastal Marshes Based on Estimated Aboveground
  Biomass Density. *Remote Sensing*, 7(4), 3507–3525.**
  doi:[10.3390/rs70403507](https://doi.org/10.3390/rs70403507)
  — peer-reviewed, open access. **[ABSTRACT]** This is the **only source I found that
  treats a quantile as the tunable knob.** Verbatim: "Elevation adjustments associated
  with these classes using both median and quartile approaches were applied to adjust
  lidar-derived elevation values closer to true bare earth elevation. The performance of
  the method was tested on 229 elevation points in the lower Apalachicola River Marsh. The
  two-class quartile-based adjusted DEM produced the best results, reducing the RMS error
  in elevation from 0.65 m to 0.40 m, a 38% improvement. The raw mean errors for the lidar
  DEM and the adjusted DEM were 0.61 ± 0.24 m and 0.32 ± 0.24 m, respectively, thereby
  reducing the high bias by approximately 49%."
  *Setting:* Apalachicola River Marsh, Florida. Note the quartile here is a quantile of
  the *error distribution within a biomass class*, used to set an offset — not a quantile
  of returns within a cell.

- **Buffington et al. (2016) [ABSTRACT]:** the LEAN method, RMSE 0.072 m, "40–75%
  improvement", beating minimum-bin gridding, mean-error correction and vegetation
  correction factors (quoted in §1.2).

### 3.3 Tuning an offset in *forest* — thin, and one explicit negative

- **Reutebuch et al. (2003) [FULL TEXT]** subtract the open-clearcut mean as a "system
  bias" (quoted in §2.1). This is the forest analogue of an offset calibration against an
  open-ground control, and it is closely parallel to what this project already does with
  stable open ground. It is one sentence in a 9-page paper, presented as a sanity check
  rather than a method.
- **Kraus & Rieger (1999) [FULL TEXT]:** "The elimination of inherent systematic errors
  allows a significant improvement of the accuracy of the laser DTM particularly in flat
  terrain. The constant value of ± 18 cm in equation (1) can be reduced down to ± 10 cm."
- **Fradette, M.-S., Leboeuf, A., Riopel, M., & Bégin, J. (2019). Method to Reduce the
  Bias on Digital Terrain Model and Canopy Height Model from LiDAR Data. *Remote Sensing*,
  11(7), 863.** doi:[10.3390/rs11070863](https://doi.org/10.3390/rs11070863)
  — peer-reviewed, open access. **[ABSTRACT]** **A negative result, and worth stating
  plainly.** Verbatim: "the bias of both DTM and CHM were calculated by subtracting two
  LiDAR datasets: high-density pixels with 21 pulses/m² (first return) and more (DTM or
  CHM reference value pixels) and low-density pixels (DTM or CHM value to correct). After
  preliminary analyses, it was concluded that the DTM did not need specific adjustment.
  In contrast, the CHM needed adjustments."
  *Setting:* Quebec forests, multiple sensors, 1 m resolution. Their density contrast is
  large but both epochs are modern; they found the *canopy* model needed correction and
  the *terrain* model did not. Our situation differs (2008 vs 2021 sensors, leaf-off vs
  green-up), but this is the clearest published statement that DTM-side adjustment is not
  automatically necessary — and it should be cited as the null we are arguing against.

---

## 4. Leaf-on vs leaf-off acquisition and DTM ground elevation

This is the best-supported of the five questions, and it supports a **large** effect in
deciduous forest.

**Simpson et al. (2017) [FULL TEXT]** — the primary result, quoted in full in §2.4:
RMSE 0.83 m leaf-on vs 0.22 m leaf-off in temperate deciduous broadleaf, same site, total-
station reference, n = 1750 at 1 m. "Leaf-off conditions improved overall DTM accuracy by
61 cm (RMSE_leaf-off = 0.22 m vs. RMSE_leaf-on = 0.83 m, n = 1750) at 1 m resolution
(Figure 11), demonstrating that leaf-on vegetation induces larger positive DTM errors."
Statistically: "Leaf-on and leaf-off DTM residuals were significantly different
(F = 3086, df = 1, p < 0.001)."

Their structural breakdown, verbatim:

> "In each of the six vertical vegetation structure categories of Table 4, DTM accuracy
> (RMSE) was better in leaf-off than leaf-on conditions (Figure 12)"

and, importantly for us, the categories are not equal:

> "with dense understory vegetation such as ferns and brambles causing the greatest
> positive DTM errors."

Caveat the authors raise themselves, verbatim:

> "Errors induced by leaf litter were not quantified; DTM_TS was produced using exact
> ground points from a total station survey, however DTM_ALS will always be erroneous in
> the presence of leaf litter because measurements cannot penetrate the leaf litter. Here
> the best accuracies in leaf-off conditions were approximately 20 cm, with mean residuals
> of less than 0.04 m (1sd < 0.18 m) for structural categories with little vegetation
> < 3.5 m tall".

---

**Stereńczak, K. & Kozak, J. (2011). Evaluation of digital terrain models generated in
forest conditions from airborne laser scanning data acquired in two seasons.
*Scandinavian Journal of Forest Research*, 26(4), 374–384.**
doi:[10.1080/02827581.2011.570781](https://doi.org/10.1080/02827581.2011.570781)
— peer-reviewed journal article. **[ABSTRACT]**

Verbatim:

> "In this study, a series of DTMs were produced from ALS data, acquired twice in one year
> (spring/summer). The study was carried out in a 1000-ha forested area in Poland. Spatial
> resolutions of output DTMs, season of data acquisition, number of vegetation layers and
> tree species in the first forest floor were evaluated to assess their influence on the
> DTM errors. Surveying methods were used to collect coordinates of 95 checkpoints. For
> various output raster resolutions and seasons of data acquisition, mean errors varied
> between −0.2 and 0.34 m, and root mean square errors varied from 0.28 to 0.79 m. Errors
> increased linearly with DTM pixel size, and their variability was significantly higher
> in DTMs derived from summer data than in DTMs derived from spring data. Effects of
> seasonality were modified by both forest structure and species composition. One-layer
> stands were more sensitive to season of data acquisition than were multilayer stands, as
> were larch and alder stands in comparison to pine and oak stands."

*Setting:* 1000 ha of Polish lowland forest; pine, oak, larch, alder; spring vs summer
acquisition in the same year. This is the **key structural finding for us**: the
seasonal effect is *modified by species composition and stand structure* — i.e. it is an
interaction, not an additive season term. That is independently consistent with this
project's finding that slope and cover interact rather than add.

*(Note: the page range 374–384 is from the journal listing and is [NOT VERIFIED]; volume,
issue, year, and DOI are verified.)*

---

**Wasser, L., Day, R., Chasmer, L., & Taylor, A. (2013). Influence of Vegetation
Structure on Lidar-derived Canopy Height and Fractional Cover in Forested Riparian Buffers
During Leaf-Off and Leaf-On Conditions. *PLoS ONE*, 8(1), e54776.**
doi:[10.1371/journal.pone.0054776](https://doi.org/10.1371/journal.pone.0054776)
— peer-reviewed, open access. **[FULL TEXT via publisher, targeted read]**

**This paper is frequently cited in leaf-on/leaf-off discussions but does not answer our
question.** It compares canopy height and fractional cover, not ground elevation. A
targeted read confirmed that the document "contains no explicit discussion of ground
return density, ground-point classification differences between leaf-off and leaf-on
acquisitions, or elevation differences in the DTM itself between the two acquisition
periods" — a single leaf-off DEM is used as the reference for *both* canopy height models.
*Setting:* Spring Creek Watershed, central Pennsylvania; deciduous simple-leaved,
deciduous compound-leaved, conifer needle, and mixed; leaf-off 26–29 April 2006, leaf-on
15–18 June 2007; 1.4 m point spacing, 2 returns. Cite it for canopy metrics only.

---

## 5. Matching two lidar epochs by adjusting the ground statistic

**No published study adjusts the per-cell ground statistic to reconcile two epochs.** The
published epoch-matching methods all operate *downstream* of the DEM, by fitting a spatial
correction surface on terrain assumed stable.

### 5.1 The pseudo-geoid approach — closest to our problem

**Viedma, O. (2022). Applying a Robust Empirical Method for Comparing Repeated LiDAR Data
with Different Point Density. *Forests*, 13(3), 380.**
doi:[10.3390/f13030380](https://doi.org/10.3390/f13030380)
— peer-reviewed, open access. **[FULL TEXT]**

Verbatim from the abstract:

> "Here, we aimed to apply an improved empirical method based on DEMs of difference, that
> adjust the ground elevation of a low-density LiDAR dataset to that of a high-density
> LiDAR one for ensuring credible vegetation changes."

> "The methodology consisted of producing 'the best DEM of difference' between low- and
> high-density LiDAR data (using the classification filter, the interpolation method and
> the spatial resolution with the lowest vertical error) to generate a local 'pseudo-geoid'
> (i.e., continuous surfaces of elevation differences) that was used to correct raw
> low-density LiDAR ground points."

> "Before correction and aggregating by sites, the vertical error of DEMs ranged from 0.02
> to −2.09 m (P50), from 0.39 to 0.85 m (NMDA) and from 0.54 to 2.5 m (RMSE). The
> segmented-based filter algorithm (CSF) showed the highest error, but there were not
> significant differences among interpolation methods or spatial resolutions. After
> correction and aggregating by sites, the vertical error of DEMs dropped significantly:
> from −0.004 to −0.016 m (P50), from 0.10 to 0.06 m (NMDA) and from 0.28 to 0.46 m
> (RMSE); and the CSF filter algorithm continued showing the greatest vertical error. The
> terrain slope and the distance to the nearest geoid point were the most important
> variables for explaining vertical accuracy. After corrections, changes in vegetation
> height were decoupled from vertical errors of DEMs."

Two points matter for us. First, this is the only paper that **searches over the
classification filter × interpolation method × resolution grid to pick the combination
with lowest vertical error** — a direct precedent for treating the DEM-construction
recipe as a fitted choice. Second, **it reports CSF as the worst-performing filter**,
which deserves attention since this project uses CSF for the 2008 epoch. Third, slope is
again the dominant explanatory variable.

*Setting:* six sites, Sierra de Gredos, central Spain; Mediterranean mountain vegetation.

### 5.2 The stable-ground correction surface — the standard geomorphic practice

**DeLong, S.B., Hammer, M.N., Engle, Z.T., Richard, E.M., Breckenridge, A.J., Gran, K.B.,
Jennings, C.E., & Jalobeanu, A. (2022). Regional-Scale Landscape Response to an Extreme
Precipitation Event From Repeat Lidar and Object-Based Image Analysis. *Earth and Space
Science*, 9(12), e2022EA002420.**
doi:[10.1029/2022EA002420](https://doi.org/10.1029/2022EA002420)
— peer-reviewed, open access. **[ABSTRACT verbatim; body quotes verified via publisher
page]**

This is our closest published analogue: repeat MN DNR airborne lidar, 8,000 km², change
detection. The workflow, with verified quoted fragments:

- flightline realignment with "BayesStripAlign v2.08 software";
- inter-epoch alignment by an "iterative closest point (ICP) approach" on stable uplands;
- a correction surface built from stable areas (slope < 3°, DoD < 0.7 m, > 100 m from
  streams) interpolated with an "inverse distance weighted algorithm and a search radius
  of 400 m"; "The final correction surface had a mean value of 0.20 m and a standard
  deviation of 0.26 m";
- DEMs made by reclassifying with "industry standard TerraScan software v020", then
  ground points "gridded to 1-m resolution DEMs with an inverse weighted distance
  algorithm with an inverse distance power weight of two and a search radius of 30 m";
- a level of detection: "A level of detection of 0.30 m (minimum mean value within an
  object) was selected to classify objects as areas of landscape change";
- residual quality on stable ground: "the mean elevation difference…was 0.002 m with a
  standard deviation of 0.103 m".

**The gap this exposes is the important part.** Their 2011 epoch was flown 3 May–2 June
and their 2012 epoch 29 October–8 November — a leaf-on/leaf-off mismatch of the same kind
as ours — and **the paper does not quantify vegetation-induced vertical bias between the
epochs.** The error budget is hardware, alignment and interpolation. A targeted read
confirmed no leaf-state term. So the standard practice in the *directly adjacent*
literature is to absorb any leaf-state offset into a stable-ground correction surface and
a generous level of detection, without ever naming it.

### 5.3 The negative result

**Fradette et al. (2019) [ABSTRACT]** concluded, for two lidar densities over Quebec
forest, that "the DTM did not need specific adjustment" while the canopy height model did.
Any claim that cross-epoch DTM adjustment is necessary has to engage with this.

### 5.4 Summary of §5

| Approach | Adjusts | Source |
|---|---|---|
| Stable-ground correction surface (IDW over stable areas) | DEM, spatially | DeLong et al. 2022 |
| "Pseudo-geoid" DoD surface, applied to raw ground points | Point cloud, spatially | Viedma 2022 |
| Recipe search over filter × interpolator × resolution | DEM construction | Viedma 2022 |
| Nothing — DTM judged adequate as-is | — | Fradette et al. 2019 |
| **Per-cell ground statistic / percentile matched between epochs** | — | **Not found** |

---

## What this means for us

1. **Percentile-matching two epochs is, so far as I can establish, unpublished.** That is
   a contribution, but it also means there is no methodological precedent to lean on and
   the manuscript must carry the full justification itself.

2. **The minimum is the wrong direction on a bluffland, and the literature already says
   why.** Ewald's channel/dike result — min-bin underpredicting by 70–90 cm where "high
   ground, a moderate to steep slope, and low ground" meet over a short distance — is the
   published statement of the downhill-capture mechanism this project independently found
   as DoD banding. Cite Ewald (2013) for it; it is the strongest external corroboration
   available, even though it is a thesis and the setting is a marsh.

3. **The published tuning precedent exists, but is parameterised as cell size.** Wang
   et al. (2009), Schmid et al. (2011) and Ewald (2013) all fit an aggregation *scale* to
   minimise bias or RMSE against survey checkpoints. Under a minimum operator, enlarging
   the neighbourhood lowers the effective quantile — so the field has been tuning a
   quantile while calling it a resolution, and paying a resolution cost for it. Fixing the
   cell size and tuning the percentile directly is the same idea with the confound removed,
   and that framing gives the method a lineage rather than making it look ad hoc.

4. **The leaf-state effect is large and well-documented in deciduous forest, and it is an
   interaction, not an additive term.** Simpson et al. (2017) give the magnitude
   (RMSE 0.83 m leaf-on vs 0.22 m leaf-off, same site, total-station reference) and the
   direction (positive — lidar reads high under leaf-on vegetation). Stereńczak & Kozak
   (2011) give the structure: "Effects of seasonality were modified by both forest
   structure and species composition." Both support treating our green-up 2021 epoch as
   the biased-high one and correcting with a cover-dependent, not constant, term.

5. **Be careful which way we point the correction.** Su & Bork (2006) is the closest
   analogue — deciduous aspen, sloping, first-generation sensor — and finds
   *over*-estimation of +0.20 m in forest, with the whole-area signed error near zero
   (+0.02 m). Reutebuch et al. (2003) likewise find lidar high in conifer. The literature
   consensus is that lidar reads **high** under vegetation. Our 2008 epoch reads **low**
   on steep forested slopes. That is not the vegetation-penetration effect the literature
   describes; it is a separate, near-nadir, slope-carried effect. **The manuscript should
   not blur the two.** The near-nadir-worst behaviour we measured is genuinely at odds
   with the usual edge-worst framing, and Su & Bork's nadir paragraph — greater signed
   errors near nadir, which they attribute to extreme values co-occurring with > 10° slopes
   — is the only prior hint of it I found, and they explicitly discount it.

6. **The slope term has a published functional form worth engaging.** Kraus & Rieger
   (1999): σH[cm] = ±(18 + 120·tan α), for forested terrain with ≥ 25 % ground
   penetration. That is a *scatter* law rather than a bias law, but it establishes both the
   tan α parameterisation and the ~10:1 ratio of slope term to flat term, which is the
   right order for what we measure.

7. **Comparability discipline is the thing the literature under-serves.** DeLong et al.
   (2022) worked on the same MN DNR data family, with the same leaf-on/leaf-off epoch
   mismatch, and their error budget contains no vegetation term at all. Viedma (2022) is
   the one paper that treats the *recipe* — filter, interpolator, resolution — as something
   that must be chosen and held identical across epochs, and even she works through a
   correction surface rather than the statistic. Whatever percentile we choose, the
   defensible claim is that **both epochs were reduced by an identical rule** and the
   percentile was selected against an external criterion, not tuned to make the DoD look
   good.

8. **Do not import marsh magnitudes.** Every number in §1 and most of §3 comes from dense,
   sub-metre, near-impenetrable herbaceous canopy on flat ground. Simpson et al. state the
   general principle themselves: "these are not correction factors."

---

## Gaps / not found

Stated plainly, because these are the load-bearing absences.

1. **No peer-reviewed comparison of per-cell order statistics (min / low percentile /
   median / mean) for a bare-earth forest DTM.** Every "gridding comparison" I found —
   Boreggio et al. (2018), Montealegre et al. (2015), and the interpolation literature
   generally — compares *interpolation algorithms*, not order statistics. This is a real
   hole, not a search failure; I probed it from four different phrasings.

2. **No use of a percentile above the median for ground elevation, anywhere.** High
   percentiles (p90, p95) appear only for canopy/crop *surface* models, where they replace
   the maximum for robustness to point density — the opposite problem. Medeiros et al.
   (2015)'s quartile-based *offset* is the nearest relative and is not the same construct.

3. **No study fixes cell size and tunes the percentile.** The tuning precedent (§3.1) is
   entirely in cell size / search radius, which conflates the quantile with the resolution.

4. **Minimum-bin gridding has never been evaluated outside coastal wetland.** All primary
   sources (Schmid 2011, Wang 2009) and the only substantial evaluation (Ewald 2013) are
   marsh. I found no upland, forest, or sloping-terrain evaluation.

5. **No epoch-matching study adjusts the ground statistic.** The two methods that exist
   (DeLong et al. 2022; Viedma 2022) both fit a spatial correction surface after gridding.

6. **Nobody quantifies leaf-state bias as a term in a DoD error budget.** DeLong et al.
   (2022) had a spring/autumn epoch pair over Minnesota forest and did not include one.
   This is arguably the single clearest opening for our contribution.

7. **Very few DTM accuracy studies in temperate deciduous forest at all.** Simpson et al.
   (2017) say so themselves: "very few studies have assessed how accurately LiDAR can
   measure surface topography under forest canopies… showing there are very few DTM
   accuracy reports in temperate deciduous forest environments" and "very few studies have
   formally assessed how vertical vegetation structure can affect DTM accuracy in
   broadleaf forests." Their Table 1 has exactly two temperate-deciduous rows against six
   temperate-conifer rows.

8. **Near-nadir-worst error on slopes is not described in this literature.** Su & Bork
   (2006) observed it and discounted it as outlier-driven. I found no other treatment. This
   remains, as previously assessed in this project, an open and likely novel result.

### Sources I could not verify

| Source | Why | What NOT to claim from it |
|---|---|---|
| Hodgson et al. (2003) *RSE* 84:295–308 | Abstract withheld by publisher; null on Crossref, OpenAlex, Semantic Scholar; ScienceDirect 403 | The "1.22 m RMSE, temperate deciduous and conifer" figure — it comes only via Simpson et al.'s Table 1 |
| Kraus & Pfeifer (1998) *ISPRS JPRS* 53:193–203 | No full text or verbatim abstract obtainable | Any specific number; the ±10 cm figure is quoted here only via Kraus & Rieger (1999) |
| Buffington et al. (2016) — "118 points" calibration figure | Appeared in a search summary only; not found in any verbatim source text | The 118-point number |
| Aryal et al. (2017) *PFG* 85:243–255, "Impact of Slope, Aspect, and Habitat-Type on LiDAR-Derived DTMs in a Near Natural, Heterogeneous Temperate Forest" | Springer paywall + IdP redirect; abstract elided on all aggregators | Anything. This looks highly relevant (temperate forest, slope × aspect × habitat) and is worth obtaining through the library — flagged as the top follow-up |
| Su & Bork (2006), the "23 times" and "(10°)" figures in §2.3 | `<`/`>` glyphs systematically absent from the PDF text layer (zero occurrences in 10 pages) | Quote them only after checking the printed page; the missing `>` is near-certain, "23 times" appears literal |
| Stereńczak & Kozak (2011) page range 374–384 | From a listing, not the article | Page numbers only; volume/issue/DOI are verified |

### Suggested next acquisitions (library access)

1. Aryal, Latifi, Heurich & Hahn (2017), *PFG* 85(4):243–255, doi:10.1007/s41064-017-0023-2 — temperate forest, slope × aspect × habitat, the nearest unexamined match to our covariates.
2. Hodgson, Jensen, Schmidt, Schill & Davis (2003), *RSE* 84(2):295–308 — the leaf-on temperate deciduous+conifer benchmark everyone cites.
3. Kraus & Pfeifer (1998), *ISPRS JPRS* 53(4):193–203 — for the asymmetric weight function in its original form, which is the theoretical justification for a sub-median statistic.
4. Schmid, Hadley & Wijekoon (2011), *JCR* 27(6A):116–132 — full text, for the actual min-bin implementation details and per-species cell-size optima.

---

## Full reference list

Peer-reviewed journal articles unless marked otherwise.

- Boreggio, M., Bernard, M., & Gregoretti, C. (2018). Evaluating the Differences of Gridding Techniques for Digital Elevation Models Generation and Their Influence on the Modeling of Stony Debris Flows Routing: A Case Study From Rovina di Cancia Basin (North-Eastern Italian Alps). *Frontiers in Earth Science*, 6, 89. doi:10.3389/feart.2018.00089 — **[ABSTRACT + body]**
- Buffington, K.J., Dugger, B.D., Thorne, K.M., & Takekawa, J.Y. (2016). Statistical correction of lidar-derived digital elevation models with multispectral airborne imagery in tidal marshes. *Remote Sensing of Environment*, 186, 616–625. doi:10.1016/j.rse.2016.09.020 — **[ABSTRACT, via NOAA NCCOS record]**
- DeLong, S.B., Hammer, M.N., Engle, Z.T., Richard, E.M., Breckenridge, A.J., Gran, K.B., Jennings, C.E., & Jalobeanu, A. (2022). Regional-Scale Landscape Response to an Extreme Precipitation Event From Repeat Lidar and Object-Based Image Analysis. *Earth and Space Science*, 9(12), e2022EA002420. doi:10.1029/2022EA002420 — **[ABSTRACT + body]**
- Evans, J.S., & Hudak, A.T. (2007). A Multiscale Curvature Algorithm for Classifying Discrete Return LiDAR in Forested Environments. *IEEE Transactions on Geoscience and Remote Sensing*, 45(4), 1029–1038. doi:10.1109/TGRS.2006.890412 — **[ABSTRACT]** (context only; classifier, not gridding)
- Ewald, M.J. (2013). *Where's the Ground Surface? Elevation Bias in LIDAR-derived Digital Elevation Models Due to Dense Vegetation in Oregon Tidal Marshes.* **MS thesis (not peer-reviewed)**, Oregon State University. <https://ir.library.oregonstate.edu/downloads/1n79h8198> — **[FULL TEXT]**
- Fradette, M.-S., Leboeuf, A., Riopel, M., & Bégin, J. (2019). Method to Reduce the Bias on Digital Terrain Model and Canopy Height Model from LiDAR Data. *Remote Sensing*, 11(7), 863. doi:10.3390/rs11070863 — **[ABSTRACT]**
- Hladik, C., & Alber, M. (2012). Accuracy assessment and correction of a LIDAR-derived salt marsh digital elevation model. *Remote Sensing of Environment*, 121, 224–235. doi:10.1016/j.rse.2012.01.018 — **[ABSTRACT]**
- Hodgson, M.E., & Bresnahan, P. (2004). Accuracy of Airborne Lidar-Derived Elevation: Empirical Assessment and Error Budget. *Photogrammetric Engineering and Remote Sensing*, 70(3), 331–339. doi:10.14358/PERS.70.3.331 — **[ABSTRACT]**
- Hodgson, M.E., Jensen, J.R., Schmidt, L., Schill, S., & Davis, B. (2003). An evaluation of LIDAR- and IFSAR-derived digital elevation models in leaf-on conditions with USGS Level 1 and Level 2 DEMs. *Remote Sensing of Environment*, 84(2), 295–308. doi:10.1016/S0034-4257(02)00114-1 — **[NOT VERIFIED]**
- Kraus, K., & Pfeifer, N. (1998). Determination of terrain models in wooded areas with airborne laser scanner data. *ISPRS Journal of Photogrammetry and Remote Sensing*, 53(4), 193–203. doi:10.1016/S0924-2716(98)00009-4 — **[NOT VERIFIED]**
- Kraus, K., & Rieger, W. (1999). Processing of laser scanning data for wooded areas. In D. Fritsch & R. Spiller (Eds.), *Photogrammetric Week '99*, 221–231. Wichmann, Heidelberg. **Conference proceedings, not peer-reviewed.** <https://phowo.ifp.uni-stuttgart.de/publications/phowo99/kraus.pdf> — **[FULL TEXT]**
- Medeiros, S., Hagen, S., Weishampel, J., & Angelo, J. (2015). Adjusting Lidar-Derived Digital Terrain Models in Coastal Marshes Based on Estimated Aboveground Biomass Density. *Remote Sensing*, 7(4), 3507–3525. doi:10.3390/rs70403507 — **[ABSTRACT]**
- Montealegre, A.L., Lamelas, M.T., & de la Riva, J. (2015). Interpolation Routines Assessment in ALS-Derived Digital Elevation Models for Forestry Applications. *Remote Sensing*, 7(7), 8631–8654. doi:10.3390/rs70708631 — **[ABSTRACT]**
- Reutebuch, S.E., McGaughey, R.J., Andersen, H.-E., & Carson, W.W. (2003). Accuracy of a high-resolution lidar terrain model under a conifer forest canopy. *Canadian Journal of Remote Sensing*, 29(5), 527–535. doi:10.5589/m03-022 — **[FULL TEXT]**
- Schmid, K.A., Hadley, B.C., & Wijekoon, N. (2011). Vertical Accuracy and Use of Topographic LIDAR Data in Coastal Marshes. *Journal of Coastal Research*, 27(6A), 116–132. doi:10.2112/JCOASTRES-D-10-00188.1 — **[ABSTRACT]**
- Simpson, J.E., Smith, T.E.L., & Wooster, M.J. (2017). Assessment of Errors Caused by Forest Vegetation Structure in Airborne LiDAR-Derived DTMs. *Remote Sensing*, 9(11), 1101. doi:10.3390/rs9111101 — **[FULL TEXT]**
- Stereńczak, K., & Kozak, J. (2011). Evaluation of digital terrain models generated in forest conditions from airborne laser scanning data acquired in two seasons. *Scandinavian Journal of Forest Research*, 26(4), 374–384. doi:10.1080/02827581.2011.570781 — **[ABSTRACT]**
- Su, J., & Bork, E. (2006). Influence of Vegetation, Slope, and Lidar Sampling Angle on DEM Accuracy. *Photogrammetric Engineering and Remote Sensing*, 72(11), 1265–1274. doi:10.14358/PERS.72.11.1265 — **[FULL TEXT]**
- Tinkham, W.T., Smith, A.M.S., Hoffman, C., Hudak, A.T., Falkowski, M.J., Swanson, M.E., & Gessler, P.E. (2012). Investigating the influence of LiDAR ground surface errors on the utility of derived forest inventories. *Canadian Journal of Forest Research*, 42(3), 413–422. doi:10.1139/x11-193 — **[ABSTRACT]**
- U.S. Geological Survey. *Lidar Base Specification 2025 rev. A* — Digital Elevation Model Surface Treatments. **Agency standard, not peer-reviewed.** <https://www.usgs.gov/ngp-standards-and-specifications/lidar-base-specification-digital-elevation-model-surface> — **[partially verified]**
- Viedma, O. (2022). Applying a Robust Empirical Method for Comparing Repeated LiDAR Data with Different Point Density. *Forests*, 13(3), 380. doi:10.3390/f13030380 — **[FULL TEXT]**
- Wang, C., Menenti, M., Stoll, M.P., Feola, A., Belluco, E., & Marani, M. (2009). Separation of Ground and Low Vegetation Signatures in LiDAR Measurements of Salt-Marsh Environments. *IEEE Transactions on Geoscience and Remote Sensing*, 47(7), 2014–2023. doi:10.1109/TGRS.2008.2010490 — **[ABSTRACT]**
- Wasser, L., Day, R., Chasmer, L., & Taylor, A. (2013). Influence of Vegetation Structure on Lidar-derived Canopy Height and Fractional Cover in Forested Riparian Buffers During Leaf-Off and Leaf-On Conditions. *PLoS ONE*, 8(1), e54776. doi:10.1371/journal.pone.0054776 — **[FULL TEXT, targeted read]**

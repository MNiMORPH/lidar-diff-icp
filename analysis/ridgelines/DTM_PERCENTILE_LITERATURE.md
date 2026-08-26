# Vertical shift of the lidar ground surface with vegetation cover, and its compensation by a per-cell order statistic

**A literature review for the 2008 MN DNR vs. 2021 USGS 3DEP change-detection study
(forested southeastern Minnesota bluffland).**

Compiled 2026-08-26. Not committed to git.

---

## The question this document is organised around

> **Does the lidar-derived ground surface shift vertically as a function of forest /
> vegetation cover, and can a per-cell percentile (or equivalent order statistic) of the
> return distribution be used to compensate that shift?**

Everything below is evidence on that single question. **The sign of a reported shift is a
result *within* the phenomenon, not a criterion for admitting a study.** A paper reporting
+0.31 m under uncut conifer and a paper reporting a negative shift are both evidence that
the ground estimate *moves with cover*, which is what we are asking. Sign, magnitude,
biome, leaf state, sensor era, point density, reference method and terrain form are
**covariates to tabulate** — things that may explain why studies disagree — never reasons
to demote a study.

The same applies to biome. Marsh, shrubland, floodplain and tropical work bears directly
on the general question even where its *numbers* do not transfer to Minnesota bluffland.
Transferability is a caveat attached to a **number**, not a demotion of a **study**.

**Status of the two halves, as of this revision.** The literature supports the *first* half
(the ground surface does shift with cover) strongly and across every biome surveyed. The
*second* half (can a per-cell percentile compensate it) has since been answered **no** for
our site by our own data — see `FRAME_2026-08-26.md`, commit 5335359: the ground column is
symmetric at every cover level, so there is no skew for a shifted percentile to exploit, and
per cell the column width anti-correlates with the correction needed. Sections 2–5 below
report what the literature says regardless; "What this means for us" reconciles the two, and
notes that two published results (Clark et al. 2004; Ewald 2013) independently anticipate our
negative.

*Revision note.* An earlier draft of this file organised around the sign of the bias and
used it to rule studies in and out — demoting marsh work to a "for completeness, low
transferability" appendix, and setting aside positive-bias forest results as "a different
effect" from our negative one. That framing was wrong, and it cost real evidence: it is
why the first draft missed Clark et al. (2004) — a peer-reviewed *forest* local-minima
gridding scheme with the window scale tuned against 3859 ground-survey points — and
concluded, incorrectly, that minimum-bin-style gridding "has never been evaluated outside
coastal wetland." It also missed Cobby et al. (2001), the one study that makes the ground
statistic itself a function of vegetation class. Both are now central to §3.

---

## How to read this document

| Tag | Meaning |
|---|---|
| **[FULL TEXT]** | I downloaded the PDF and read the passage. Quotes transcribed from the PDF text layer. |
| **[ABSTRACT]** | Publisher's abstract retrieved verbatim (OpenAlex / Semantic Scholar / publisher page). Body details **not** verified. |
| **[SECONDARY]** | The number is quoted from a *different* paper's description of the source. The original was not opened. Both papers are named. |
| **[NOT VERIFIED]** | Cited for completeness; source text not obtained. No numbers asserted. |

Nothing here is stated from a title alone, and nothing from recall.

---

## 1. Synthesis table: the spread of reported ground shifts

Sign convention: **positive = lidar ground surface sits ABOVE the surveyed ground.**
Where a source's own convention is ambiguous or internally inconsistent, that is noted in
the body section and flagged here with (†).

| Study | Biome / cover type | Terrain | Leaf state | Sensor era & density | Reference method | Ground statistic / DTM construction | Reported shift (signed) | How cover was quantified | Tag |
|---|---|---|---|---|---|---|---|---|---|
| Reutebuch et al. 2003 | Temperate conifer (70-yr Douglas-fir), clearcut→uncut | W. Washington, 0–45° | Evergreen (spring) | Saab TopEye helicopter, 1999; 4.22 raw returns/m², 0.58 ground pts/m² | Total station, 347 checkpoints | Proprietary last-return filter; 1.52 m grid; bilinear | **+0.16** clearcut, **+0.18** heavy thin, **+0.18** light thin, **+0.31** uncut; all-data **+0.22 ± 0.24** | 4 managed canopy-density classes (TPH) | [FULL TEXT] |
| Su & Bork 2006 | Aspen parkland: deciduous aspen forest, shrub, grassland, meadow | Alberta, knob-and-kettle, 5–10 m relief | Deciduous, leaf state not stated | Last-return, 1998-era | Total station + DGPS, 27 benchmarks, 256 plots | IDW, 1.5 m | **+0.20** aspen forest; **−0.22** lowland meadow; **+0.02** overall | 8 field vegetation classes | [FULL TEXT] |
| Hodgson & Bresnahan 2004 | 6 land-cover classes incl. deciduous & evergreen forest | Richland Co., SC | Not stated | Optech ALTM 1210, 1207 m AGL, 2 m posting | Total station + rapid-static GPS **at the lidar point itself** | Points surveyed directly (no interpolation) | **−0.06 (±0.23)** to **+0.06 (±0.19)** across covers — **both signs**; RMSE 17–19 cm (pavement/low grass/evergreen) to 26 cm (deciduous) | 6 land-cover categories | RMSE [ABSTRACT]; signed range [SECONDARY via Hopkinson et al. 2005] |
| Hodgson et al. 2003 | Leaf-on pine/deciduous landscape | South Carolina | **Leaf-on** | Last-return | 1470 survey-grade points | TIN from ground points | RMSE 0.93 m overall; **0.33** low grass, **1.22** scrub/shrub, **1.53** deciduous; ~2 m error increase from 0–2° to 6–8° slopes in shrub/scrub | Land-use classes | [SECONDARY via Clark et al. 2004] — see §2.4 conflict |
| Clark et al. 2004 | Tropical rain forest (old-growth, secondary, agroforestry, pasture) | La Selva, Costa Rica; flat→steep | Evergreen leaf-on | FLI-MAP, 0.33 m DSM support | 3859 ground-survey points | **Local-minima cell in a grid, scale tuned 5/10/15/20/30 m**; then IDW or OK | **All DTMs positive mean-signed error**; +0.08 to +1.10 depending on scheme; RMSE 2.29 m best, 1.95 m old-growth, 0.58 m flat open canopy | Land-use / forest-age classes | [FULL TEXT] |
| Cobby et al. 2001 | Floodplain incl. deciduous forest | UK floodplain, 10–15° in forest | **Leaf-on** | Last-return DSM, 2 m support | Field survey (n small: 5 and 12) | **Local-minima in 5×5-px (10 m) windows, then algorithm tailored to short vs tall vegetation classes** | RMSE **0.17 m** short vegetation; **3.99 m** deciduous forest on steeper slopes | Short vs tall vegetation classes | [SECONDARY via Clark et al. 2004] |
| Hopkinson et al. 2005 | Boreal wetland: grass/herb, low shrub, willow, aquatic, aspen, black spruce, jack pine | N. Alberta, 75 m total relief | August (leaf-on) | Optech ALTM 2050 | GPS reference points, 127 vegetated | Ground-classified points; also rasterised | **+0.07 ± 0.16** mean over vegetated transects (+0.04 rasterised); **no significant difference** for grass/herbs → **+0.15** aquatic | Ducks Unlimited vegetation classes | [FULL TEXT proceedings] + [ABSTRACT journal] |
| Töyrä et al. 2003 | Boreal wetland: graminoid, willow scrub | N. Alberta | Not stated | — | Field survey | — | **+0.07 (±0.15)** graminoid; **+0.15 (±0.26)** willow scrub | Vegetation classes | [SECONDARY via Hopkinson et al. 2005] |
| Simpson et al. 2017 | Temperate deciduous broadleaf (alder, field maple, hazel) + grassland clearing | UK, median slope 5.7° | **Both leaf-on and leaf-off, same site** | NERC ARSF; 3–5 returns/m² | Total station + GNSS, n = 657 | lasground last-return; IDW 1 m | RMSE **0.83 m leaf-on** vs **0.22 m leaf-off**; direction stated as positive (lidar high) (†) | Pgap in 3 vertical strata → 6 structural categories | [FULL TEXT] |
| Stereńczak & Kozak 2011 | Temperate lowland forest: pine, oak, larch, alder | Poland, 1000 ha | **Spring vs summer, same year** | ALS | 95 checkpoints | Varying raster resolution | Mean errors **−0.2 to +0.34 m**; RMSE 0.28–0.79 m; summer variability > spring | Species + number of vegetation layers | [ABSTRACT] |
| Tinkham et al. 2012 | High-biomass western conifer | Western USA, slopes >30° | Evergreen | — | 54 ground survey plots | Two classification algorithms compared | **Vegetation structure: no influence**; error variability rises above 30° slope; RMSE 0.24 m | Vegetation structure metrics | [ABSTRACT] |
| Salleh et al. 2015 | Tropical: rubber plantation, mixed forest | Bentong, Malaysia | Evergreen | Riegl, "low density" | GCPs → reference DTM | ATIN ground filter | MBE by canopy density — rubber: **+0.011 / +0.008 / +0.002**; mixed forest: **+0.162 / +0.075 / +0.019** (70–80 / 81–90 / 91–100 % cover). RMSE rises with cover: 0.230→0.789 | **Canopy cover % from non-ground return fraction** | [FULL TEXT] |
| Duchan et al. 2026 | Czech forest, drainage ditches | Czech Republic | Not stated | ALS | 706 GNSS + total station | National DTM product | **+0.415 m** mean elevation error, RMSE 0.464 m; **ditch geometry, not canopy height, was the dominant predictor** | Forest height as proxy for density/closure; ground reflection density | [FULL TEXT] |
| Hladik & Alber 2012 | Salt marsh, 10 cover classes | Sapelo I., GA; flat | n/a herbaceous | Optech Gemini ALTM, 125 kHz | RTK GPS | Vendor bare-earth + per-class offsets | **+0.03 to +0.25 m** by cover class; pooled **+0.10 ± 0.12** → **−0.01 ± 0.09** after correction | Field-mapped species classes | [ABSTRACT] |
| Ewald 2013 | Oregon tidal marsh, 12 vegetation associations + pasture | Flat, with channels/dikes | n/a herbaceous | DOGAMI OLC | >13,000 RTK GPS | **Minimum-bin**, cell size tuned 0.1–6.0 m | **+0.104** to **+0.488** by association; **0** in open cover (FVA 4.5 cm RMSE); **−0.70 to −0.90 m** along channels/dikes | Percent cover by species, quadrats | [FULL TEXT] |
| Schmid et al. 2011 | Coastal marsh (*Spartina*, *Juncus*) | Flat | n/a herbaceous | — | Survey-grade GPS (280 pts) | **Minimum bin**, cell size tuned 2–10 m | Bias **reduced by 12 cm**, accuracy improved 8 cm vs as-received | Marsh species classes | [ABSTRACT] |
| Wang et al. 2009 | Salt marsh, Venice Lagoon | Flat | n/a herbaceous | — | Field observations (240 pts) | **Filter window size tuned to minimise DTM bias** | **+0.022 ± 0.064 m** residual after tuning | Marsh species | [ABSTRACT] |
| Medeiros et al. 2015 | Coastal marsh | Apalachicola, FL; flat | n/a herbaceous | — | 229 elevation points | Vendor DTM + **median vs quartile** class offsets | **+0.61 ± 0.24** raw → **+0.32 ± 0.24** adjusted; RMSE 0.65 → 0.40 m | Biomass-density classes from ASTER + IfSAR + CHM | [ABSTRACT] |
| Buffington et al. 2016 | 17 Pacific-coast tidal marshes | Flat | n/a herbaceous | — | RTK GPS | LEAN regression on NDVI; benchmarked vs min-bin | RMSE **0.072 m**, 40–75 % improvement over bare-earth DEM | NDVI from NAIP imagery | [ABSTRACT] |
| Viedma 2022 | Mediterranean mountain vegetation | Sierra de Gredos, Spain | — | Low- vs high-density lidar pair | High-density lidar as benchmark | Filter × interpolator × resolution search; **pseudo-geoid** correction | P50 **+0.02 to −2.09 m** before → **−0.004 to −0.016 m** after | Vegetation height | [FULL TEXT] |
| Fradette et al. 2019 | Quebec forest | Canada | — | High (21 pulses/m²) vs low density | High-density lidar as reference | 1 m resolution | **DTM judged to need no adjustment**; CHM did | Stand density, species composition | [ABSTRACT] |
| DeLong et al. 2022 | NE Minnesota forest | 8,000 km² | 2011 **leaf-on** (May–Jun) vs 2012 **leaf-off** (Oct–Nov) | Repeat MN DNR ALS | Stable terrain, ICP | TerraScan ground; IDW 1 m, power 2, r = 30 m | Stable-ground residual **+0.002 ± 0.103 m** after correction surface (surface itself mean 0.20 ± 0.26 m) | **Not quantified — no vegetation term** | [ABSTRACT + body] |

**What the spread shows.** Reported shifts run from **−0.22 m to +1.10 m**, and the same
land-cover label gives different signs in different studies. The variance is not noise
about a consensus; it is structure. The covariates that plausibly carry it, in rough order
of how strongly the sources implicate them:

1. **Undergrowth stature and leaf state** (Simpson: RMSE 0.83 vs 0.22 m at one site;
   Stereńczak & Kozak: season effect *modified by* species and stand layering).
2. **Terrain form at short horizontal distances** — concave/convex microtopography, not
   cover (Ewald's channels −0.70 to −0.90 m; Duchan's ditch geometry beating canopy
   height; Hodgson et al.'s ~2 m rise from 0–2° to 6–8° slopes *within* shrub/scrub).
3. **Slope** (Hodgson & Bresnahan 2× from 1.5° to 25°; Su & Bork 2× from <2° to >10°;
   Tinkham >30°; Salleh r = 0.87–0.99; Kraus's ±(18 + 120·tan α) cm law).
4. **The ground statistic and its scale** (Clark: RMSE 2.29–5.09 m across local-minima
   scales at one site with one point cloud — a *larger* range than most cover effects).
5. **Ground-return density / penetration** (Duchan: forest height reduces ground
   reflection density and increases distance to nearest ground reflection, p < 0.05).
6. **Reference method** — Hodgson & Bresnahan surveyed the lidar point itself and got the
   tightest, most symmetric numbers of anyone; GPS-under-canopy studies get the loosest.

Point 4 deserves emphasis: **at least one study finds the choice of ground statistic moves
the answer more than the cover does.** That is the strongest single argument that the
statistic is a first-class parameter rather than a processing detail.

---

## 2. Does the ground surface shift with cover?

Yes, in every study that stratified by cover. The interesting content is in *how* and
*how much*, and in the fact that the sign is not universal.

### 2.1 The shift is real, cover-dependent, and directly recommended as a target for correction

**Hopkinson, C., Chasmer, L., Sass, G.Z., Creed, I.F., Sitar, M., Kalbfleisch, W., &
Treitz, P. (2005). Vegetation class dependent errors in lidar ground elevation and canopy
height estimates in a boreal wetland environment. *Canadian Journal of Remote Sensing*,
31(2), 191–206.** doi:[10.5589/m05-007](https://doi.org/10.5589/m05-007)
— peer-reviewed journal article. **[ABSTRACT]**, plus the earlier conference version
(Hopkinson et al., ISPRS Archives XXXVI-8/W2, <https://www.isprs.org/proceedings/xxxvi/8-w2/HOPKINSON.pdf>,
**not peer-reviewed**) **[FULL TEXT]**.

The title is the finding. Verbatim from the journal abstract:

> "These data were analysed to quantify vegetation class dependent errors in lidar ground
> surface elevation and vegetation canopy surface height… Aquatic vegetation was
> associated with the largest error in lidar ground surface definition (+0.15 m, SD = 0.22,
> probability of no difference in height P < 0.01), likely a result of saturated ground
> conditions."

From the conference version's conclusions **[FULL TEXT]**:

> "ground classified LiDAR data for 127 RPs over vegetated transects, an average bias of
> +0.07 (± 0.16 m) was found (+0.04 m in rasterised LiDAR data). The observation of ground
> height errors in vegetated areas is consistent with the findings of other studies (Töyrä
> et al. 2003; Hodgson and Bresnahan, 2004). The vertical bias was found to vary with
> vegetation cover, from no significant difference for grass and herbs to +0.15 m for
> aquatic vegetation."

And then, the sentence that most directly answers our question — a peer-reviewed
recommendation that the **ground-point extraction rule itself be made cover-dependent**:

> "These observations support the rationale that ground level LiDAR point classification
> should be vegetation class dependent (e.g. Cobby et al., 2001)."

*Setting:* Utikuma Uplands, northern Alberta boreal wetland complex; 40 × 6 km transect,
total relief ~75 m; trembling aspen, jack pine, black spruce, willow, low shrub, grass/herb,
aquatic; Optech ALTM 2050, August 2002. Wetland and low-relief, so the *magnitudes* are not
ours — but the *structure of the claim* (bias varies by class; extraction should therefore
vary by class) is exactly our thesis, stated in 2005.

Also recovered from this paper, **[SECONDARY]** for Töyrä et al.:

> "known comparable statistics detailing vertical bias in ground elevation for a northern
> wetland environment range from + 0.07 m (± 0.15 m) to + 0.15 m (± 0.26 m) for graminoid
> and willow scrub, respectively (Töyra et al., 2003)."

**Töyrä, J., Pietroniro, A., & Martz, L.W. (2001). Multisensor Hydrologic Assessment of a
Freshwater Wetland. *Remote Sensing of Environment*, 75(2), 162–173.**
doi:[10.1016/S0034-4257(00)00164-4](https://doi.org/10.1016/S0034-4257(00)00164-4)
— peer-reviewed. **[NOT VERIFIED]** (abstract null on OpenAlex; note Hopkinson cites a
2003 Töyrä work, which may be a different paper — do not cite the 2001 record for the
2003 numbers without checking).

### 2.2 The sign is NOT universal — both signs occur, in the same study

This is the correction that the narrow framing had buried. From Hopkinson et al.'s
literature review **[FULL TEXT]**, describing Hodgson & Bresnahan **[SECONDARY]**:

> "Analysis carried out by Hodgson and Bresnahan (2004) in South Carolina demonstrated
> both positive and negative absolute errors, from – 0.06 m (± 0.23 m) to + 0.06 m
> (± 0.19 m), in LiDAR ground elevation for various ground covers."

So the paper most often cited for "vegetation biases lidar high" in fact reports a signed
range **straddling zero across its six land-cover classes**. The published picture is *a
cover-dependent shift of either sign*, not a one-way positive bias.

Independently, **Su & Bork (2006) [FULL TEXT]** report opposite signs in adjacent cover
types within one survey:

> "In identifying a tendency to under- or over-estimate elevations, the mean signed errors
> in Figure 3c show that elevations within aspen forest were over-estimated (0.20 m) while
> those in lowland meadows were under-estimated (0.22 m)."

and their pooled figure is near zero — "overall signed error and RMSE were 0.02 m and
0.59 m, respectively" — which is precisely what happens when two cover classes shift in
opposite directions and are averaged. **A pooled bias near zero is not evidence of no
cover effect.** That is a methodological warning for us as much as a literature finding.

And **Stereńczak & Kozak (2011) [ABSTRACT]** span both signs across resolutions and
seasons: "mean errors varied between −0.2 and 0.34 m".

Our own result — 2008 reading low on steep forested slopes — sits inside this published
range rather than outside it. The earlier draft's claim that it was "a different effect"
that "should not be blurred" was an over-reach: there is no single-signed consensus for it
to contradict.

### 2.3 Cover quantified continuously, with signed bias: the closest methodological analogue

**Salleh, M.R.M., Ismail, Z., & Abdul Rahman, M.Z. (2015). Accuracy Assessment of
Lidar-Derived Digital Terrain Model (DTM) with Different Slope and Canopy Cover in
Tropical Forest Region. *ISPRS Annals of the Photogrammetry, Remote Sensing and Spatial
Information Sciences*, II-2/W2, 183–189.**
doi:[10.5194/isprsannals-II-2-W2-183-2015](https://doi.org/10.5194/isprsannals-II-2-W2-183-2015)
— **conference annals; the paper states "This contribution has been peer-reviewed. The
double-blind peer-review was conducted on the basis of the full paper."** **[FULL TEXT]**

Method note worth copying: cover is derived **from the point cloud itself**, not from
external imagery — "the point clouds belong to non-ground are then used in determining the
relative percentage of canopy cover" (verbatim), i.e. the non-ground return fraction. That
is the same family of data-derived cover metric this project uses.

Verbatim result:

> "The results show that terrain slope has high correlation for both study area (0.993 and
> 0.870) with the RMSE of the LiDAR-derived DTM. This is similar to canopy cover where
> high value of correlation (0.989 and 0.924) obtained. This indicates that the accuracy of
> airborne LiDAR-derived DTM is significantly affected by terrain slope and canopy caver of
> study area."

Tables 3 and 4, transcribed from the PDF (MBE = mean bias error, signed):

*By slope, rubber-tree area / mixed-forest area:*

| Slope | RMSE (m) | MAE (m) | MBE (m) | | RMSE (m) | MAE (m) | MBE (m) |
|---|---|---|---|---|---|---|---|
| 0–5° | 0.613 | 0.364 | +0.020 | | 0.379 | 0.153 | +0.010 |
| 6–10° | 0.723 | 0.410 | +0.002 | | 0.589 | 0.012 | −0.024 |
| 11–15° | 0.890 | 0.619 | +0.017 | | 0.590 | 0.425 | +0.054 |

*By canopy density, rubber-tree area / mixed-forest area:*

| Canopy cover | RMSE (m) | MAE (m) | MBE (m) | | RMSE (m) | MAE (m) | MBE (m) |
|---|---|---|---|---|---|---|---|
| 70–80 % | 0.230 | 0.035 | +0.011 | | 0.333 | 0.333 | +0.162 |
| 81–90 % | 0.437 | 0.255 | +0.008 | | 0.367 | 0.179 | +0.075 |
| 91–100 % | 0.789 | 0.495 | +0.002 | | 0.576 | 0.076 | +0.019 |

Note the structure: **RMSE rises steeply with cover while the signed bias *falls*.**
Scatter and shift move in opposite directions here. Any correction fitted to RMSE rather
than to signed bias would get the sign of the adjustment wrong. (Their own summary
sentence — "Most of the slope class shows the positive bias means underestimate exists in
DTM generated at high terrain" — is confusingly worded; I quote the tables rather than
that sentence.)

*Setting:* Bentong, Pahang, Malaysia; rubber plantation and mixed tropical forest; Riegl
system, described by the authors as "low density"; ATIN ground filter.

### 2.4 Cover-stratified error in a leaf-on temperate pine/deciduous landscape

**Hodgson, M.E., Jensen, J.R., Schmidt, L., Schill, S., & Davis, B. (2003). An evaluation
of LIDAR- and IFSAR-derived digital elevation models in leaf-on conditions with USGS Level
1 and Level 2 DEMs. *Remote Sensing of Environment*, 84(2), 295–308.**
doi:[10.1016/S0034-4257(02)00114-1](https://doi.org/10.1016/S0034-4257(02)00114-1)
— peer-reviewed. **[NOT VERIFIED directly]** — abstract withheld by the publisher and null
on Crossref, OpenAlex and Semantic Scholar; ScienceDirect returns 403.

Two independent **[SECONDARY]** descriptions, which I quote because they do not agree.

From **Clark et al. (2004) [FULL TEXT]**:

> "Working with last-return lidar data flown over a leaf-on pine/deciduous forest
> landscape, Hodgson et al. (2003) identified ground points through a combination of
> proprietary software and human interpretation. A comparison of DTM elevation against
> 1470 survey-grade field measurements had an overall RMSE of 0.93 m. DTM error differed
> significantly by land use. Although RMSE was 0.33 m for low grass, it increased to 1.22
> and 1.53 m for the more structurally complex scrubs/shrub and deciduous vegetation types,
> respectively. Furthermore, these researchers found that in the dense, multi-layered
> shrub/scrub class, there was a highly significant increase in DTM error of roughly 2 m
> from lowest (0–2°) to steepest (6–8°) slopes, which the authors attributed to vertical
> inaccuracies over relatively short horizontal distances under complex canopy."

From **Simpson et al. (2017) Table 1 [FULL TEXT]**, the same source is credited with
"Temperate deciduous and conifer 1.22", "Temperate grass 0.37", "Temperate pine 0.45",
"Temperate shrub 1.53".

**These conflict.** Clark assigns 1.22 to scrub/shrub and 1.53 to deciduous; Simpson
assigns 1.22 to deciduous+conifer and 1.53 to shrub. **Do not cite either attribution
without obtaining Hodgson et al. (2003) itself.** The one thing both agree on is the
*ordering* — grass lowest by a wide margin, structurally complex vegetation ~4× higher.

The slope×cover interaction Clark reports — ~2 m of extra error from 0–2° to 6–8° *within a
single cover class* — is the closest published statement to this project's finding that
slope and cover interact rather than add.

### 2.5 Terrain form can beat cover

**Duchan, M., Mráz, V., Tichá, A., Jankovský, M., & Zlatuška, K. (2026). The Influence of
Forest Cover on the Accuracy of Aerial Laser Scanning-Derived Digital Elevation Models for
Detecting Drainage Ditches in Forests in the Czech Republic. *Forests*, 17(2), 162.**
doi:[10.3390/f17020162](https://doi.org/10.3390/f17020162)
— peer-reviewed, open access. **[FULL TEXT]**

Verbatim:

> "The results indicate a positive elevation bias, with a mean elevation error of 0.415 m
> and an RMSE of 0.464 m, 54.7% higher than the 0.3 m declared in the DTM technical report.
> Forest height, acting as a proxy for forest structural density and canopy closure, was
> significantly associated with a reduction in ground reflection density and an increase in
> the distance to the nearest ground reflection (p < 0.05)."

> "Crucially, multiple regression analysis revealed that forest height was not the primary
> driver of elevation error; instead, ditch geometry was the most significant predictor.
> Narrower ditches exhibited substantially higher errors than wider ones, regardless of the
> canopy height. Furthermore, while ground reflection density decreased in mature stands,
> this reduction did not significantly diminish DTM vertical accuracy, suggesting that some
> of the LiDAR reflections of low vegetation could be misclassified as ground reflections,
> decreasing accuracy."

Two things for us. First, **a cover-driven reduction in ground-return density did not by
itself degrade vertical accuracy** — which is a published caution against treating
penetration/return-density as a proxy for bias, and independently consistent with this
project's finding that ground-return fraction is a poor canopy proxy. Second, **concave
microtopography dominated**: narrow ditches, i.e. exactly the "low ground over a short
horizontal distance" geometry that also defeats minimum-bin gridding in Ewald.

### 2.6 The evergreen counterpoint

**Tinkham, W.T., Smith, A.M.S., Hoffman, C., Hudak, A.T., Falkowski, M.J., Swanson, M.E.,
& Gessler, P.E. (2012). Investigating the influence of LiDAR ground surface errors on the
utility of derived forest inventories. *Canadian Journal of Forest Research*, 42(3),
413–422.** doi:[10.1139/x11-193](https://doi.org/10.1139/x11-193)
— peer-reviewed. **[ABSTRACT]**

> "The study further compared two LiDAR classification algorithms and found no significant
> difference in their performance. Vegetation structure was found to have no influence,
> whereas increased variability in the vertical error was observed on slopes exceeding 30°,
> illustrating that these algorithms are not limited by high-biomass western coniferous
> forests, but that slope and sensor accuracy both play important roles."

This is a genuine null on the cover effect, and it should be reported as such. Its setting
is high-biomass **evergreen conifer** — the one forest type where no leaf-state axis exists,
and where Reutebuch also found the cover effect small ("strikingly small", verbatim). The
pattern across §2 is that the cover effect is largest where **deciduous leaf state and
low-stature undergrowth** vary, and smallest in evergreen conifer.

### 2.7 Reutebuch et al. 2003 — the conifer benchmark in full

**Reutebuch, S.E., McGaughey, R.J., Andersen, H.-E., & Carson, W.W. (2003). Accuracy of a
high-resolution lidar terrain model under a conifer forest canopy. *Canadian Journal of
Remote Sensing*, 29(5), 527–535.** doi:[10.5589/m03-022](https://doi.org/10.5589/m03-022)
— peer-reviewed. **[FULL TEXT]**

> "The mean DTM error was 0.22 ± 0.24 m (mean ± SD). DTM elevation errors for four tree
> canopy cover classes were: clearcut 0.16 ± 0.23 m, heavily thinned 0.18 ± 0.14 m, lightly
> thinned 0.18 ± 0.18 m, and uncut 0.31 ± 0.29 m. These DTM errors show a slight increase
> with canopy density but the differences are strikingly small."

Table 2 (lidar DTM minus surveyed checkpoint, m), transcribed:

| Canopy class | Mean | SD | Min | Max | n |
|---|---|---|---|---|---|
| Clearcut | 0.16 | 0.23 | −0.48 | 0.61 | 38 |
| Heavy thinned | 0.18 | 0.14 | −0.11 | 0.41 | 21 |
| Lightly thinned | 0.18 | 0.18 | −0.63 | 0.69 | 147 |
| Uncut | 0.31 | 0.29 | −0.60 | 1.31 | 120 |

Table 3 groupings, transcribed: no near-ground vegetation +0.15 m (n = 132); vegetation
within 6 m of ground +0.26 m (n = 212); slope < 18 % +0.21 m (n = 174); slope ≥ 18 %
+0.22 m with SD rising 0.20 → 0.28 (n = 173). **Note that the cover contrast (+0.15 vs
+0.26) is larger than the slope contrast (+0.21 vs +0.22) in this conifer stand** — the
opposite ordering to Tinkham. The two are reconcilable: Reutebuch's split is on
*near-ground* vegetation presence, Tinkham's on overstory structure.

And the empirical offset calibration, verbatim (see also §3.4):

> "The observed error in the clearcut area (0.16 m) is very similar to the lidar
> manufacturer's stated accuracy of ±0.15 m (Baltsavias, 1999). If one assumes that this
> error in the open, bare-ground clearcut area is the system bias and adjusts the individual
> checkpoint errors to remove this bias, then 69% of the observed checkpoint errors are
> within ±0.22 m (the observed SD of the DTM grid error)."

### 2.8 Su & Bork 2006 in full, including the near-nadir observation

**Su, J. & Bork, E. (2006). Influence of Vegetation, Slope, and Lidar Sampling Angle on DEM
Accuracy. *Photogrammetric Engineering and Remote Sensing*, 72(11), 1265–1274.**
doi:[10.14358/PERS.72.11.1265](https://doi.org/10.14358/PERS.72.11.1265)
— peer-reviewed. **[FULL TEXT]**

Beyond the signed errors quoted in §2.2:

> "Finally, RMSE values indicated the lidar-derived DEM accuracy generally decreased as
> slope gradient increased: the RMSE at slopes over 10° was twice that found when slopes
> were less than 2°. This finding was similar to Hodgson and Bresnahan (2004), who observed
> errors on slopes of 25° to be twice that found on relatively flat areas."

On sampling angle — the only prior hint of a near-nadir-worst pattern I located, quoted in
full including the authors' own scepticism:

> "Signed errors and RMSEs were generally greater when lidar data were collected close to
> nadir (less than 3°) relative to those sampled in angle classes further away from the
> central flightline (Figure 3b). However, this pattern may be attributed to the presence
> of extreme errors. The mean top five signed errors near nadir were 23 times larger than
> their corresponding overall signed errors. Moreover, extreme errors were accompanied by
> high (10°) slope gradients, which may also have contributed to the observed elevation
> errors."

*Transcription note (retained from the first draft, because the self-check matters).* That
passage is quoted exactly as the PDF text layer renders it. Two checks: (i) "(10°)" almost
certainly reads "(>10°)" in print — the extracted text of this paper contains **zero** `<`
or `>` characters across all 10 pages despite repeated slope-threshold discussion, so those
glyphs are systematically lost; (ii) "23 times" is **not** a mangled "2–3": en-dashes
survive extraction elsewhere in the same file (the reference page range "3482–3486"), so
"23 times larger" appears literal, and it is physically plausible because the comparison is
top-five extremes against an overall mean signed error near zero (+0.02 m area-wide).
**Verify both against the printed page before quoting in the manuscript.** I initially
reconstructed these as "2–3 times" and ">10°"; the glyph audit overturned the first
reconstruction.

Their abstract's own summary of the angle result — "lidar sampling angle had little impact
on measured error" — is the claim they stand behind.

### 2.9 Marsh and wetland: same phenomenon, different magnitudes

Reported here as evidence on the general question, with the transferability caveat attached
to the *numbers*.

**Ewald, M.J. (2013). *Where's the Ground Surface? Elevation Bias in LIDAR-derived Digital
Elevation Models Due to Dense Vegetation in Oregon Tidal Marshes.* MS thesis, Oregon State
University.** <https://ir.library.oregonstate.edu/downloads/1n79h8198>
— **MS thesis, not peer-reviewed.** **[FULL TEXT]**

> "The fundamental vertical accuracy (FVA) of the LIDAR datasets was 4.5 cm root mean
> square error (RMSE) and had no consistent positive or negative bias in open landcover.
> Within wetland vegetation communities, my results suggest that LIDAR estimates of the
> ground surface in tidal wetlands are typically 10 cm to 30 cm above GPS measurements.
> Plant associations dominated by *Carex obnupta* and *Carex lyngbyei* exhibited the largest
> discrepancy between LIDAR and GPS measurements (mean discrepancies 36.6 cm and 48.8 cm
> respectively). The smallest errors observed in the study were about 10 cm to 11 cm"

Per-association values, verbatim: association F (*Carex obnupta*) "36.6 cm (95% CI: 29.6 cm
to 43.6 cm) in the minimum-bin DEM and 39.4 cm (95% CI: 32.2 cm to 46.5 cm) in the DOGAMI
bare-earth DEM"; association H (*Carex lyngbyei*) "48.8 cm (95% CI: 40.3 cm to 57.3 cm) in
the minimum-bin DEM and 45.1 cm (95% CI: 36.4 cm to 53.8 cm)"; association A
(*Deschampsia*/succulents) "10.4 cm (95% CI: 5.6 cm to 15.2 cm)"; association G (*Distichlis
spicata*) "10.6 cm (95% CI: 3.9 cm to 17.3 cm)".

**The zero-bias open-cover result is as important as the vegetated numbers** — it is a
clean demonstration that the shift is *caused by cover*, with the same sensor, same
processing, same site. That is the control our own open-ground stable-terrain check plays.

**Hladik, C. & Alber, M. (2012). Accuracy assessment and correction of a LIDAR-derived salt
marsh digital elevation model. *Remote Sensing of Environment*, 121, 224–235.**
doi:[10.1016/j.rse.2012.01.018](https://doi.org/10.1016/j.rse.2012.01.018)
— peer-reviewed. **[ABSTRACT]**

> "We found that DEM mean vertical errors for different cover classes ranged from 0.03 to
> 0.25 m in comparison to the RTK ground truth data, with the larger offsets for taller
> vegetation."

Vegetation *height* as the ordering variable, cleanly stated.

**Clark, M.L., Clark, D.B., & Roberts, D.A. (2004). Small-footprint lidar estimation of
sub-canopy elevation and tree height in a tropical rain forest landscape. *Remote Sensing
of Environment*, 91(1), 68–89.**
doi:[10.1016/j.rse.2004.02.008](https://doi.org/10.1016/j.rse.2004.02.008)
— peer-reviewed. **[FULL TEXT]** (see §3.2 for its gridding method, which is the main
reason it matters here)

> "In old-growth forests, RMS error on steep slopes was 0.67 m greater than on flat slopes.
> On flatter slopes, variation in vegetation complexity associated with land use caused
> highly significant differences in DTM error distribution across the landscape. The highest
> DTM accuracy observed in this study was 0.58-m RMSE, under flat, open-canopy areas with
> relatively smooth surfaces. Lidar ground retrieval was complicated by dense, multi-layered
> evergreen canopy in old-growth forests, causing DTM overestimation that increased RMS
> error to 1.95 m."

Note the clean factorial statement: slope effect *and* cover effect, measured separately,
in the same survey.

---

## 3. Can a per-cell order statistic compensate the shift?

This is where the reorganisation pays. There are **two independent traditions** of
order-statistic ground estimation — a marsh one and a forest one — and the first draft saw
only the marsh one.

### 3.1 The marsh tradition: "minimum-bin gridding"

**Schmid, K.A., Hadley, B.C., & Wijekoon, N. (2011). Vertical Accuracy and Use of
Topographic LIDAR Data in Coastal Marshes. *Journal of Coastal Research*, 27(6A), 116–132.**
doi:[10.2112/JCOASTRES-D-10-00188.1](https://doi.org/10.2112/JCOASTRES-D-10-00188.1)
— peer-reviewed. **[ABSTRACT]**

The clearest primary definition:

> "Custom digital elevation model (DEM) generation techniques and point classification
> processes can be used to improve estimates of ground elevations in coastal marshes. The
> simplest of these methods is minimum bin gridding, which extracts the lowest elevation
> value included within a user-specified search window and assigns that value to the
> appropriate DEM grid cell."

> "By employing the minimum bin technique to the bare-earth classified LIDAR data, the
> overall bias in the resultant surface was reduced by 12 cm, and the vertical accuracy was
> improved by 8 cm when compared with the 'as-received' data."

> "Despite lowering the spatial resolution of the DEM, the application of these techniques
> significantly improves the vertical accuracy of the LIDAR-derived bare-earth surfaces."

**Wang, C., Menenti, M., Stoll, M.P., Feola, A., Belluco, E., & Marani, M. (2009).
Separation of Ground and Low Vegetation Signatures in LiDAR Measurements of Salt-Marsh
Environments. *IEEE Transactions on Geoscience and Remote Sensing*, 47(7), 2014–2023.**
doi:[10.1109/TGRS.2008.2010490](https://doi.org/10.1109/TGRS.2008.2010490)
— peer-reviewed. **[ABSTRACT]**

> "In this paper, we introduce reliable methods to remove random and systematic errors and
> to register raw data, as well as a new procedure, to determine the optimal filter window
> size to separate ground and canopy returns. A limited amount of field observations is used
> to determine the size of the filtering window which produces the minimally biased estimates
> of the digital terrain model (DTM)."

> "We apply this procedure to a study marsh within the Venice Lagoon, Italy, and obtain a
> high-accuracy DTM. The error (z_LiDAR − z_field) is 2.2 cm, with a standard deviation of
> 6.4 cm."

**Ewald (2013) [FULL TEXT]** is the fullest published evaluation, and documents the
**two-sided failure** that matters most to us:

> "DEM accuracy increased with cell size until an inflection point near 1.4 m as the
> influence of vegetation is mitigated by the minimum-bin gridding technique (Figure 2.3a).
> Low features within the landscape were captured by the gridding technique and degrade DEM
> performance after cell size enlarged beyond the optimum."

> "Even at the optimum cell size, the DEM is still positively biased when compared to known
> ground elevations. Mean LIDAR-GPS discrepancy remains positive until a cell size of 1.6 m
> is achieved. At cell sizes greater than 1.6 m, DEM are negatively biased as the minimum-bin
> method continues to capture and favor low features within the landscape."

> "Minimum-bin LIDAR-derived DEM elevations underpredict (typically 70 cm to 90 cm below)
> the measured elevation along channels and the man-made dike that form the southern edge of
> the site adjacent to the Coquille River. These features are characterized by high ground, a
> moderate to steep slope, and low ground over a short horizontal distance. For example, the
> minimum-bin filter is likely to select a LIDAR return from an adjacent low riverbank rather
> than the surveyed wetland surface. The likelihood of upslope areas being assigned an
> elevation lower than the true ground elevation increases as the cell size is increased."

And the honest limit of the method:

> "Throughout the vegetation types we evaluated, the minimum-bin DEM performs slightly
> better than the DOGAMI bare-earth DEM. With 95% confidence, the DOGAMI bare-earth DEM
> elevation is between 2.0 cm and 3.1 cm above the minimum-bin DEM elevation across the
> entire dataset (mean 2.5 cm, paired two-sided t-test, p-value < 0.001)."

> "Unfortunately, our results show that LIDAR estimates of the ground surface are positively
> biased even when the minimum-bin technique is used. This suggests that the LIDAR laser
> pulse never reaches the ground surface within the vegetation communities we studied."

**Buffington, K.J., Dugger, B.D., Thorne, K.M., & Takekawa, J.Y. (2016). Statistical
correction of lidar-derived digital elevation models with multispectral airborne imagery in
tidal marshes. *Remote Sensing of Environment*, 186, 616–625.**
doi:[10.1016/j.rse.2016.09.020](https://doi.org/10.1016/j.rse.2016.09.020)
— peer-reviewed. **[ABSTRACT via NOAA NCCOS record page; ScienceDirect returned 403]**

> "Using 17 study sites along the Pacific coast of the U.S., we achieved an average root
> mean squared error (RMSE) of 0.072 m, with a 40–75% improvement in accuracy from the lidar
> bare earth DEM. Results from our method compared favorably with results from three other
> methods (minimum-bin gridding, mean error correction, and vegetation correction factors)."

This is the benchmarking statement: minimum-bin gridding is one of three recognised
correction families, and a fitted statistical model beat all three. (A widely repeated "118
points" calibration figure appeared only in a search summary — **[NOT VERIFIED]**.)

### 3.2 The forest tradition: "local-minima" ground retrieval — and it IS tuned against ground truth

This is the material the sign-based framing caused me to miss. It is peer-reviewed, it is
in forest, and it is in *Remote Sensing of Environment*.

**Clark et al. (2004) [FULL TEXT]**, §2.5.1, verbatim:

> "The local-minima algorithm proceeded as follows: a grid of non-overlapping, square cells
> was overlaid on top of the original DSM. Within each grid cell, one local-minima DSM cell
> (0.33-m support) was selected and identified as a ground return. This procedure resulted in
> a population of ground-return cells for each of the five grid scales considered
> independently: 5, 10, 15, 20 and 30 m… the above local-minima scheme is analogous to
> selecting the lowest return in a square footprint of a specified scale (i.e., 5, 10 m,
> etc.). Ground-return cells identified at each scale were then used in separate geostatistical
> interpolation schemes (described below) that generated DTMs with a 1-m cell size. Samples
> from each DTM were compared to 3859 co-located reference points… **The overall RMS errors of
> the resulting DTMs were used as the basis for the selection of a final ground-retrieval/
> interpolation scheme.**"

The tuning result, verbatim:

> "errors ranging from 2.29 to 5.09 m, using either IDW or OK for surface interpolation (data
> not shown). The scale with the lowest RMSE for both interpolation methods was found to be
> 20 m. For this tropical landscape and lidar sampling density, 20 m appears to be the
> near-optimum scale to identify ground returns with the local-minima approach… This optimal
> scale is likely determined by the average crown dimensions and canopy gap characteristics in
> old-growth forest, which comprises 69% of the study area."

Three things follow, all of them useful to us:

1. **The choice of ground statistic and its scale moved RMSE from 2.29 m to 5.09 m** on one
   point cloud at one site — a 2.2× swing, larger than most of the *cover* effects in §2.
   The statistic is not a processing detail.
2. **The optimum scale is set by canopy geometry** — mean maximum crown diameter measured
   at 19.6 m, optimum window 20 m. That is a physical, transferable selection criterion, and
   a far better justification than "we tried some values." It also predicts that our optimum
   should track *our* crown/gap scale, not Clark's.
3. They then refined it with a **multi-scale** scheme: "the iterative-addition scheme (i.e.,
   iteratively adding local-minima from 20, 15 to 10-m scales) improved the OK-interpolation
   RMSE by 0.10 m, resulting in an overall RMSE of 2.29 m."

And, critically for the sign question:

> "All DTMs had a positive mean-signed error, and so they tended to overestimate elevation."

Mean signed errors transcribed from Table 3 — IDW (local-minima) +0.68, IDW
(iterative-addition) +0.08, OK (local-minima) +1.10, OK (iterative-addition) +0.97 m.
*Transcription caveat:* the text states the OK-vs-IDW overestimation gap is "up to 0.87 m
higher" where these figures give 0.89; treat the table values as approximate pending a
check against the printed page. **Even a pure local-minimum operator left a +0.68 to +1.10 m
positive shift under tropical old-growth** — i.e. the minimum is not a floor that
automatically removes the vegetation shift.

Also, on why they preferred kriging — directly relevant to our choice of a robust statistic:

> "OK was determined to be a superior interpolation scheme because it smoothed fine-scale
> variance created by spurious understory heights in the ground-point dataset."

> "This smoothing of the variance across space tends to minimize the influence of spurious
> understory vegetation or downed trunks that are inevitably included in the DTM
> interpolation."

Clark also situates local-minima within the general taxonomy, verbatim:

> "A relatively simple approach is to find local-minima relative to neighboring samples at a
> specified scale and/or search configuration (Cobby et al., 2001; Petzold et al., 1999).
> Resulting ground samples (i.e., local minima) must then be interpolated to form a surface."

### 3.3 The one study that makes the ground statistic a FUNCTION OF VEGETATION

**Cobby, D.M., Mason, D.C., & Davenport, I.J. (2001). Image processing of airborne scanning
laser altimetry data for improved river flood modelling. *ISPRS Journal of Photogrammetry
and Remote Sensing*, 56(2), 121–138.**
doi:[10.1016/S0924-2716(01)00039-9](https://doi.org/10.1016/S0924-2716(01)00039-9)
— peer-reviewed. **[NOT VERIFIED directly]** (abstract null on OpenAlex; full text not
obtained). Described **[SECONDARY]** by Clark et al. (2004) **[FULL TEXT]**:

> "Cobby et al. (2001) developed an automated ground-retrieval scheme for a floodplain
> environment that included deciduous forests with leaf-on conditions. An initial DTM was
> interpolated from local-minima cells retrieved from non-overlapping, 5 × 5-pixel windows
> (10-m side) overlaid on a last-return DSM (2-m support). **The final DTM was achieved by
> tailoring the ground-retrieval algorithm to short and tall vegetation classes.** While
> terrain under short vegetation could be predicted with a 0.17-m RMSE (n = 5), the RMSE was
> 3.99 m (n = 12) under deciduous forests on steeper slopes (10–15°)."

**This is the direct precedent for our approach**: a local-minimum order statistic whose
*rule is switched by vegetation class*, in leaf-on deciduous forest on slopes. It is also
independently endorsed as the right direction by **Hopkinson et al. (2005) [FULL TEXT]**:
"ground level LiDAR point classification should be vegetation class dependent (e.g. Cobby
et al., 2001)."

Two caveats, stated plainly. The validation n is tiny (5 and 12 points). And the 3.99 m
RMSE under leaf-on deciduous forest on 10–15° slopes is a warning, not an endorsement — it
is the worst forest number in this whole review, from the method closest to ours. **Getting
Cobby et al. (2001) in full is the single highest-value acquisition on this list.**

### 3.4 Tuning an offset or a statistic against ground truth

Grouped by what is tuned.

**Tuning the aggregation scale under a fixed minimum operator:**
- Clark et al. 2004 — 5/10/15/20/30 m, RMSE-minimised against 3859 points, optimum 20 m **[FULL TEXT]**
- Wang et al. 2009 — filter window size, "minimally biased estimates of the DTM" **[ABSTRACT]**
- Schmid et al. 2011 — cell size 2–10 m; per **Ewald (2013) [SECONDARY]**: "selecting an
  optimal cell size of 4.0 m in *Spartina alterniflora* and a cell size of 10.0 m in *Juncus
  roemerianus* by minimizing the Root Mean Square Error (RMSE) on 280 survey-grade GPS
  measurements in South Carolina." **Note this is a per-vegetation-class scale** — another
  instance of §3.3's idea.
- Ewald 2013 — cell size 0.1–6.0 m **[FULL TEXT]**: "Minimizing MAE yielded (1 in Table 2.1)
  an optimum cell size of 1.4 m (1.92 m²) for this dataset. If RMSE (2 in Table 2.1) was used
  to define the optimum DEM cell size instead of MAE, the optimum cell size was 1.2 m (1.44
  m²)." **The loss function moved the optimum by 0.2 m and the bias by 2.9 cm** — a caution
  for anyone tuning a statistic against checkpoints.

**Tuning an additive offset by cover class:**
- Hladik & Alber 2012 **[ABSTRACT]**: "We developed species-specific correction factors for
  ten cover classes… reducing the overall mean DEM error from 0.10 ± 0.12 (SD) to −0.01 ±
  0.09 m (SD), and the Root Mean Square Error from 0.16 m to 0.10 m."
- Medeiros, S., Hagen, S., Weishampel, J., & Angelo, J. (2015). *Remote Sensing*, 7(4),
  3507–3525. doi:[10.3390/rs70403507](https://doi.org/10.3390/rs70403507) **[ABSTRACT]** —
  the only source treating a **quantile** as the knob: "Elevation adjustments associated with
  these classes using both median and quartile approaches were applied to adjust lidar-derived
  elevation values closer to true bare earth elevation… The two-class quartile-based adjusted
  DEM produced the best results, reducing the RMS error in elevation from 0.65 m to 0.40 m, a
  38% improvement. The raw mean errors for the lidar DEM and the adjusted DEM were 0.61 ± 0.24
  m and 0.32 ± 0.24 m, respectively, thereby reducing the high bias by approximately 49%."
- Buffington et al. 2016 **[ABSTRACT]** — LEAN, NDVI as the continuous cover covariate.
- Reutebuch et al. 2003 **[FULL TEXT]** — open-clearcut mean subtracted as "the system bias"
  (§2.7). The forest analogue of an open-ground offset calibration, and the direct precedent
  for this project's use of stable open ground.
- Kraus & Rieger 1999 **[FULL TEXT]** — "The elimination of inherent systematic errors allows
  a significant improvement of the accuracy of the laser DTM particularly in flat terrain. The
  constant value of ± 18 cm in equation (1) can be reduced down to ± 10 cm (K. Kraus, N.
  Pfeifer, 1998)."

**The negative:**
- **Fradette, M.-S., Leboeuf, A., Riopel, M., & Bégin, J. (2019). Method to Reduce the Bias on
  Digital Terrain Model and Canopy Height Model from LiDAR Data. *Remote Sensing*, 11(7),
  863.** doi:[10.3390/rs11070863](https://doi.org/10.3390/rs11070863) **[ABSTRACT]**: "the bias
  of both DTM and CHM were calculated by subtracting two LiDAR datasets: high-density pixels
  with 21 pulses/m² (first return) and more… and low-density pixels… After preliminary analyses,
  it was concluded that the DTM did not need specific adjustment. In contrast, the CHM needed
  adjustments." This is the null our claim must engage with.

**The structural observation.** Under a minimum operator, enlarging the neighbourhood is
monotonically equivalent to lowering the effective quantile of a fixed neighbourhood. So the
literature has been **tuning a quantile all along** — just parameterised as cell size or window
size, which pays a resolution cost and entangles the statistic with the support. Clark's 2.29 →
5.09 m swing and Ewald's 1.2-vs-1.4 m sensitivity are both really quantile sensitivity in
disguise. Fixing the cell size and tuning the percentile directly is the same idea with the
confound removed. **That framing gives our method a lineage rather than making it look ad hoc,
and it is the single most useful thing this review produces.**

### 3.5 The continuous analogue: asymmetric (skew) robust interpolation

**Kraus, K. & Pfeifer, N. (1998). Determination of terrain models in wooded areas with airborne
laser scanner data. *ISPRS Journal of Photogrammetry and Remote Sensing*, 53(4), 193–203.**
doi:[10.1016/S0924-2716(98)00009-4](https://doi.org/10.1016/S0924-2716(98)00009-4)
— peer-reviewed. **[NOT VERIFIED]** — no full text or verbatim abstract obtained.

**Kraus, K. & Rieger, W. (1999). Processing of laser scanning data for wooded areas.** In
*Photogrammetric Week '99*, 221–231. Wichmann, Heidelberg. **Conference proceedings, not
peer-reviewed.** <https://phowo.ifp.uni-stuttgart.de/publications/phowo99/kraus.pdf> **[FULL TEXT]**

> "This algorithm estimates the skewness of the error distribution of the laser scanner data in
> forests and assigns small weights to those points that show large positive errors during the
> interpolation with filtering. The process results in a classification of the laser points in
> terrain and off-terrain (mainly vegetation) points."

This is a *soft* low-quantile: instead of taking the k-th order statistic, it down-weights high
residuals continuously, so the fitted surface settles toward the low returns without ever
committing to a single point. It is the principled alternative to a hard percentile and should
be cited as such.

The accuracy law, verbatim:

> "σH[cm] = ± (18 + 120·tanα)"
>
> "Equation (1) is valid for a ground penetration rate of the laser signal of at least 25 % and
> a good mixture of vegetation and ground points for the whole region."

Clark et al. (2004) **[FULL TEXT]** describes the same method independently: "Kraus and Pfeifer
(1998) used an automated, iterative technique that interpolated a mean surface from the lidar
cloud of xyz points and then successively removed or down-weighted points with residuals higher
than a specified threshold."

### 3.6 What the delivered-product conventions do instead

**Boreggio, M., Bernard, M., & Gregoretti, C. (2018). *Frontiers in Earth Science*, 6, 89.**
doi:[10.3389/feart.2018.00089](https://doi.org/10.3389/feart.2018.00089) **[ABSTRACT + body]** —
compares "linear triangulation, natural neighbor, nearest neighbor, Inverse Distance to a Power,
ANUDEM, Radial Basis Functions, and ordinary kriging" (verbatim): **interpolators only, no order
statistic.** *Setting:* steep, largely unvegetated Italian Dolomites debris-flow basin.

**Montealegre, A.L., Lamelas, M.T., & de la Riva, J. (2015). *Remote Sensing*, 7(7), 8631–8654.**
doi:[10.3390/rs70708631](https://doi.org/10.3390/rs70708631) **[ABSTRACT]** — same in a forest
setting: "six interpolation routines were tested"; "The Triangulated Irregular Network (TIN) to
raster interpolation method produced the best result in the validation process with the training
data set while the Inverse Distance Weighted (IDW) routine was the best in the validation with
GPS (RMSE of 2.68 cm and RMSE of 37.10 cm, respectively)." Note the ~14× gap between
self-validation and GPS validation — a caution about validating a gridding choice against the
same point cloud that produced it.

**USGS Lidar Base Specification 2025 rev. A** — agency standard, not peer-reviewed.
<https://www.usgs.gov/ngp-standards-and-specifications/lidar-base-specification-digital-elevation-model-surface>
**[partially verified]** — I read only the "Digital Elevation Model Surface Treatments" page. It
specifies which points are **excluded** ("Bare-earth lidar points (serving as mass points) that
are in close proximity to any breakline shall be classified as Ignored Ground (class 20) and
shall be excluded from the DEM generation process", verbatim). **The page I read states no
interpolation method and no per-cell statistic**, and I did not verify the separate "Data
Processing and Handling Requirements" page, so I make no claim about TIN being mandated elsewhere
in the spec. What is safe: the governing specification for 3DEP treats the DEM as an
*interpolation over all valid ground points*. A percentile-gridded product is therefore a
departure from the delivered-product convention and must be built identically for both epochs.

### 3.7 Does anyone use a percentile above the median?

**Not for ground elevation. I found no case**, across four differently-phrased searches. What
exists:

- **Canopy and crop *surface* models** use high percentiles (Zp90, Zp95) in place of the maximum,
  for robustness to point-density variation. That is the upper envelope — the opposite problem —
  and supplies no argument for raising the ground statistic.
- **Medeiros et al. (2015)** apply a *quartile* of the error distribution within a biomass class
  to set an offset. Nearest relative; not a quantile of returns within a cell.

The *mechanism* that would justify moving up from the minimum is well documented, though:
Ewald's downhill capture (§3.1), Duchan's narrow-ditch dominance (§2.5), and Clark's preference
for kriging specifically because it smooths "spurious understory heights" and "downed trunks"
(§3.2). All three say the low tail contains non-ground returns *below* true ground as well as
above it. Nobody has taken the next step of asking which percentile is optimal at fixed cell size.

**Our own answer to that unasked question is negative**, and worth recording here so the gap is
not mistaken for an opportunity: at Elba the classified-ground column is symmetric at every cover
level in both epochs, so no shifted percentile — above or below the median — has any skew to
exploit (`FRAME_2026-08-26.md`, commit 5335359). The literature's silence on percentiles above the
median is therefore probably not an oversight so much as a route that does not lead anywhere in
forest; the compensation that is actually needed is a translation, not a re-ranking.

---

## 4. Leaf-on vs leaf-off

**Simpson, J.E., Smith, T.E.L., & Wooster, M.J. (2017). Assessment of Errors Caused by Forest
Vegetation Structure in Airborne LiDAR-Derived DTMs. *Remote Sensing*, 9(11), 1101.**
doi:[10.3390/rs9111101](https://doi.org/10.3390/rs9111101) — peer-reviewed, open access.
**[FULL TEXT]**

> "In the presence of leaf-on vegetation, DTM accuracy is severely reduced, with low-stature
> undergrowth vegetation (such as ferns) causing the greatest errors (RMSE > 1 m). Errors are
> lower under leaf-off conditions (RMSE = 0.22 m)."

> "Leaf-off conditions improved overall DTM accuracy by 61 cm (RMSE_leaf-off = 0.22 m vs.
> RMSE_leaf-on = 0.83 m, n = 1750) at 1 m resolution (Figure 11), demonstrating that leaf-on
> vegetation induces larger positive DTM errors. Leaf-on and leaf-off DTM residuals were
> significantly different (F = 3086, df = 1, p < 0.001)."

> "In each of the six vertical vegetation structure categories of Table 4, DTM accuracy (RMSE)
> was better in leaf-off than leaf-on conditions (Figure 12)"

> "Results demonstrate that leaf-on vegetation causes greater DTM error (RMSE = 0.83 m) than
> leaf-off vegetation (RMSE = 0.22) across all vegetation categories. Furthermore, DTM accuracy
> is not affected by all vegetation structures equally; with dense understory vegetation such as
> ferns and brambles causing the greatest positive DTM errors. Grassland vegetation yields the
> most accurate DTMs."

**(†) Sign-convention caveat, flagged because it matters.** The Methods state the reference was
subtracted from the lidar surface — "The control DTM was subtracted from the ALS-derived DTM
(DTM_ALS) to produce difference rasters" — making positive = lidar high. But the caption of
Figure 11 reads "Digital Terrain Model (DTM) error (DTM_TS − DTM_ALS)", the opposite. **The paper
is internally inconsistent.** The intended physical direction is unambiguous from the Discussion:

> "The results of the present study suggest that without an adequate ground control scheme, such
> applications may be prone to significant positive biases, especially in areas of leaf-on, open
> canopy forest with large amounts of ground cover."

Slope was **not** a driver at this site, and the null is worth recording with its range:

> "The median slope within the plot was 5.7°. The relationship between vertical DTM residuals and
> slope at 1m resolution was examined using a non-parametric GAM. Slope has no meaningful effect
> on residuals, with slope explaining 0.25% of the deviance, and a poor goodness of fit (adjusted
> R2 = 0.002, Figure 10)."

That null is informative only up to ~6° and does not speak to bluffland slopes.

The authors' own caveats:

> "Errors induced by leaf litter were not quantified; DTM_TS was produced using exact ground
> points from a total station survey, however DTM_ALS will always be erroneous in the presence of
> leaf litter because measurements cannot penetrate the leaf litter. Here the best accuracies in
> leaf-off conditions were approximately 20 cm, with mean residuals of less than 0.04 m (1sd <
> 0.18 m) for structural categories with little vegetation < 3.5 m tall"

> "The results presented here show how DTM accuracy is relatively affected by vegetation
> structure, and as such they cannot be applied absolutely to other forest environments (i.e.,
> these are not correction factors)."

*Setting:* < 1 ha UK deciduous broadleaf plot (alder, field maple, hazel); median slope 5.7°;
NERC ARSF leaf-on 24 June 2014, leaf-off 9 March 2009; ~3–5 returns m⁻²; IDW at 1 m; total-station
reference (n = 657).

---

**Stereńczak, K. & Kozak, J. (2011). Evaluation of digital terrain models generated in forest
conditions from airborne laser scanning data acquired in two seasons. *Scandinavian Journal of
Forest Research*, 26(4), 374–384.**
doi:[10.1080/02827581.2011.570781](https://doi.org/10.1080/02827581.2011.570781) — peer-reviewed.
**[ABSTRACT]**

> "Spatial resolutions of output DTMs, season of data acquisition, number of vegetation layers and
> tree species in the first forest floor were evaluated to assess their influence on the DTM
> errors. Surveying methods were used to collect coordinates of 95 checkpoints. For various output
> raster resolutions and seasons of data acquisition, mean errors varied between −0.2 and 0.34 m,
> and root mean square errors varied from 0.28 to 0.79 m. Errors increased linearly with DTM pixel
> size, and their variability was significantly higher in DTMs derived from summer data than in
> DTMs derived from spring data. Effects of seasonality were modified by both forest structure and
> species composition. One-layer stands were more sensitive to season of data acquisition than
> were multilayer stands, as were larch and alder stands in comparison to pine and oak stands."

**The key structural finding for us**: the seasonal effect is *modified by species composition and
stand structure* — an interaction, not an additive season term. Independently consistent with this
project's slope×cover interaction result. *(Page range 374–384 is from a listing and is [NOT
VERIFIED]; volume, issue, year and DOI are verified.)*

---

**Hodgson et al. (2003)** — leaf-on temperate pine/deciduous, §2.4. **[SECONDARY]**, attributions
in conflict.

**Cobby et al. (2001)** — leaf-on deciduous floodplain forest, §3.3. **[SECONDARY]**.

**Wasser, L., Day, R., Chasmer, L., & Taylor, A. (2013). *PLoS ONE*, 8(1), e54776.**
doi:[10.1371/journal.pone.0054776](https://doi.org/10.1371/journal.pone.0054776) — peer-reviewed,
open access. **[FULL TEXT, targeted read]** — frequently cited in leaf-on/leaf-off discussions but
**does not address ground elevation**. A targeted read confirmed the document "contains no explicit
discussion of ground return density, ground-point classification differences between leaf-off and
leaf-on acquisitions, or elevation differences in the DTM itself between the two acquisition
periods"; a single leaf-off DEM is the reference for *both* canopy height models. *Setting:* Spring
Creek Watershed, central Pennsylvania; deciduous simple-leaved, deciduous compound-leaved, conifer
needle, mixed; leaf-off 26–29 April 2006, leaf-on 15–18 June 2007; 1.4 m point spacing. Cite for
canopy metrics only.

---

## 5. Matching two lidar epochs

**No published study reconciles two epochs by adjusting the per-cell ground statistic.** The
existing methods operate downstream of the DEM, fitting a spatial correction surface on terrain
assumed stable.

### 5.1 Pseudo-geoid correction, with a recipe search

**Viedma, O. (2022). Applying a Robust Empirical Method for Comparing Repeated LiDAR Data with
Different Point Density. *Forests*, 13(3), 380.**
doi:[10.3390/f13030380](https://doi.org/10.3390/f13030380) — peer-reviewed, open access.
**[FULL TEXT]**

> "Here, we aimed to apply an improved empirical method based on DEMs of difference, that adjust
> the ground elevation of a low-density LiDAR dataset to that of a high-density LiDAR one for
> ensuring credible vegetation changes."

> "The methodology consisted of producing 'the best DEM of difference' between low- and
> high-density LiDAR data (using the classification filter, the interpolation method and the
> spatial resolution with the lowest vertical error) to generate a local 'pseudo-geoid' (i.e.,
> continuous surfaces of elevation differences) that was used to correct raw low-density LiDAR
> ground points."

> "Before correction and aggregating by sites, the vertical error of DEMs ranged from 0.02 to
> −2.09 m (P50), from 0.39 to 0.85 m (NMDA) and from 0.54 to 2.5 m (RMSE). The segmented-based
> filter algorithm (CSF) showed the highest error, but there were not significant differences among
> interpolation methods or spatial resolutions. After correction and aggregating by sites, the
> vertical error of DEMs dropped significantly: from −0.004 to −0.016 m (P50), from 0.10 to 0.06 m
> (NMDA) and from 0.28 to 0.46 m (RMSE); and the CSF filter algorithm continued showing the greatest
> vertical error. The terrain slope and the distance to the nearest geoid point were the most
> important variables for explaining vertical accuracy. After corrections, changes in vegetation
> height were decoupled from vertical errors of DEMs."

Three things for us: it is the only paper that **searches the filter × interpolator × resolution
grid** and picks the lowest-error combination — the nearest precedent for treating the DEM recipe
as fitted; it reports **CSF as the worst-performing filter**, which deserves attention since this
project uses CSF for the 2008 epoch; and slope is again the dominant explanatory variable.

### 5.2 Stable-ground correction surface — and the leaf-state gap

**DeLong, S.B., Hammer, M.N., Engle, Z.T., Richard, E.M., Breckenridge, A.J., Gran, K.B., Jennings,
C.E., & Jalobeanu, A. (2022). Regional-Scale Landscape Response to an Extreme Precipitation Event
From Repeat Lidar and Object-Based Image Analysis. *Earth and Space Science*, 9(12), e2022EA002420.**
doi:[10.1029/2022EA002420](https://doi.org/10.1029/2022EA002420) — peer-reviewed, open access.
**[ABSTRACT verbatim; body quotes verified via publisher page]**

Our closest published analogue — repeat MN DNR airborne lidar, 8,000 km², change detection.
Verified quoted fragments of the workflow: flightline realignment with "BayesStripAlign v2.08
software"; inter-epoch alignment by an "iterative closest point (ICP) approach" on stable uplands;
a correction surface from stable areas (slope < 3°, DoD < 0.7 m, > 100 m from streams) interpolated
with an "inverse distance weighted algorithm and a search radius of 400 m", with "The final
correction surface had a mean value of 0.20 m and a standard deviation of 0.26 m"; ground points
reclassified with "industry standard TerraScan software v020" then "gridded to 1-m resolution DEMs
with an inverse weighted distance algorithm with an inverse distance power weight of two and a
search radius of 30 m"; "A level of detection of 0.30 m (minimum mean value within an object) was
selected to classify objects as areas of landscape change"; stable-ground residual "the mean
elevation difference…was 0.002 m with a standard deviation of 0.103 m".

**The gap is the point.** Their 2011 epoch was flown 3 May–2 June and their 2012 epoch 29 October–8
November — a leaf-on/leaf-off mismatch of the same kind as ours — and **the paper does not quantify
vegetation-induced vertical bias between the epochs.** The error budget is hardware, alignment and
interpolation. Standard practice in the directly adjacent literature is therefore to absorb any
leaf-state offset into a stable-ground correction surface and a generous level of detection, without
ever naming it.

### 5.3 Summary

| Approach | What is adjusted | Source |
|---|---|---|
| Stable-ground correction surface (IDW over stable areas) | DEM, spatially | DeLong et al. 2022 |
| "Pseudo-geoid" DoD surface applied to raw ground points | Point cloud, spatially | Viedma 2022 |
| Recipe search over filter × interpolator × resolution | DEM construction | Viedma 2022 |
| Ground-retrieval rule switched by vegetation class (single epoch) | Ground point selection | Cobby et al. 2001 |
| Nothing — DTM judged adequate as-is | — | Fradette et al. 2019 |
| **Per-cell ground statistic matched between epochs** | — | **Not found** |

---

## What this means for us

**Read this section against `FRAME_2026-08-26.md` (commit 5335359), which post-dates the
first draft of this review.** That work established, from our own data, that the
classified-ground column is **symmetric at every cover level** in both epochs (|Bowley
skew| ≤ 0.048, tail ratio 0.98–1.14), that cover changes only the column's **width**, and
that the epoch offset is therefore **a pure translation of the column, not a distortion of
it**. The consequence recorded there is that a fixed-percentile correction is closed: a
percentile shift can only deliver `location + k × width`, the need is not proportional to
width, and per cell the width *anti*-correlates with the need. I take that as
authoritative and align what follows to it.

So the second half of this review's organising question — *can a per-cell percentile
compensate the shift?* — now has a measured answer for our site: **no.** What the
literature contributes is (1) strong support for the *first* half of the question, (2) an
independent explanation of *why* the percentile route fails, which our result and the
published record agree on, and (3) the correct alternative.

### The first half of the question is well supported, and that premise survives

1. **The lidar ground surface shifts with cover — this is established, not contested.**
   Every study that stratified by cover found a shift. Hopkinson et al. (2005) state the
   corollary as a recommendation: ground point classification "should be vegetation class
   dependent." This is the premise our correction rests on, and it is uncontroversial.

2. **There is no single-signed consensus for our negative shift to contradict.** Reported
   shifts span −0.22 m to +1.10 m; Hodgson & Bresnahan's six land-cover classes span
   −0.06 to +0.06 m — both signs in one survey; Su & Bork report +0.20 m in aspen and
   −0.22 m in adjacent meadow. Our 2008-reads-low result sits inside the published range.
   The first draft's framing of it as "a different effect" that must not be blurred was an
   over-reach and should not appear in the manuscript.

3. **A pooled bias near zero is not evidence of no cover effect.** Su & Bork's whole-area
   signed error is +0.02 m *because* forest and meadow shift oppositely. Any DoD-wide
   summary we quote must be stratified or it will hide the effect we are measuring.

4. **Leaf state is large and acts as an interaction.** Simpson: RMSE 0.83 m leaf-on vs
   0.22 m leaf-off, same site, total-station reference. Stereńczak & Kozak: "Effects of
   seasonality were modified by both forest structure and species composition." Both
   support a cover-dependent, not constant, correction for our green-up 2021 epoch — and
   "cover-dependent translation" is exactly the form the frame document arrives at.

### The literature already anticipated why a low order statistic cannot fix this

This is the most useful thing the broader framing recovered, and it converges with our
negative result from a completely independent direction.

5. **A pure minimum operator does not remove the vegetation shift.** Clark et al. (2004)
   ran local-minima ground retrieval at the RMSE-optimal scale under tropical old-growth
   and still reported "All DTMs had a positive mean-signed error" — +0.68 to +1.10 m
   remaining. The lowest return in the cell was still far above the ground.

6. **Ewald (2013) says why, in one sentence.** "our results show that LIDAR estimates of
   the ground surface are positively biased even when the minimum-bin technique is used.
   This suggests that the LIDAR laser pulse never reaches the ground surface within the
   vegetation communities we studied." **If the ground is not in the return distribution,
   no order statistic of that distribution can recover it.** That is the same conclusion
   our symmetry result reaches by a different route: the column has no skew to exploit
   because the missing information is missing, not hidden in a tail.

7. **Where a low percentile *does* move the answer, it moves it for the wrong reason.**
   Ewald's min-bin underpredicts by 70–90 cm where "high ground, a moderate to steep
   slope, and low ground over a short horizontal distance" meet; Duchan et al. (2026) found
   narrow concave ditch geometry beat canopy height as the dominant error predictor; Clark
   preferred kriging specifically because it smooths "spurious understory heights" and
   "downed trunks." On a bluffland, a low percentile is a **terrain-form** operator far
   more than a **cover** operator — which is precisely the failure mode our per-cell
   anti-correlation between width and need describes quantitatively.

8. **The one method in the literature that varies the ground statistic by vegetation is a
   class switch, not a percentile shift.** Cobby et al. (2001) "tailor[ed] the
   ground-retrieval algorithm to short and tall vegetation classes" — and its worst number
   in this whole review (RMSE 3.99 m under leaf-on deciduous forest on 10–15° slopes) is a
   warning rather than an endorsement. Hopkinson's endorsement of that approach is likewise
   an endorsement of **class-dependent handling**, i.e. a translation per class, not of a
   shifted quantile.

### What to carry into the manuscript instead

9. **Frame the correction as a cover-dependent translation, and cite the lineage for
   that.** The published offset-calibration family is the right ancestry: Hladik & Alber
   (2012) per-class correction factors (mean error 0.10 ± 0.12 → −0.01 ± 0.09 m);
   Buffington et al. (2016) LEAN with NDVI as a continuous cover covariate (RMSE 0.072 m,
   40–75 % improvement); Medeiros et al. (2015) class offsets; Reutebuch et al. (2003)
   subtracting the open-clearcut mean as "the system bias" — the direct forest precedent
   for our use of stable open ground. All are translations. None is an order statistic.

10. **Report the closed percentile route as a result, not a dead end silently dropped.**
    It is a clean negative with a mechanism: symmetric columns, width-proportional
    sensitivity, and a per-cell anti-correlation between width and need. Read together with
    Clark's residual +0.68–1.10 m and Ewald's "the pulse never reaches the ground," it
    generalises into a statement worth making in print: **order-statistic gridding
    compensates vegetation shift only to the extent that ground returns exist in the tail;
    where they do not, only an external, cover-dependent translation can help.** No source
    I found states this explicitly, and several skirt it.

11. **The statistic-choice literature still matters — for comparability, not correction.**
    Clark's RMSE ranged 2.29–5.09 m across local-minima scales on one point cloud at one
    site: the construction rule can move the answer more than cover does. That is the
    argument for our discipline of reducing **both epochs by an identical rule**, which
    remains the most defensible methodological claim we have. Viedma (2022) is the only
    paper treating the recipe (filter × interpolator × resolution) as something to fix
    identically across epochs — and reports CSF as her worst-performing filter, which is
    worth engaging since we use CSF for 2008.

12. **The clearest open contribution is unchanged and is now sharper.** DeLong et al.
    (2022) worked the same MN DNR data family with the same leaf-on/leaf-off epoch mismatch
    and carried **no vegetation term at all** in their error budget. Nobody quantifies
    leaf-state bias as a term in a DoD error budget. We can — as a cover-dependent
    translation, with the percentile alternative tested and ruled out.

13. **Transferability caveats attach to numbers, not studies.** Marsh, tropical and
    boreal-wetland magnitudes should not be imported — Simpson says it best about his own
    work: "these are not correction factors." But their *structures* — cover-dependence,
    class-switched retrieval, scale sensitivity, terrain-form failure modes, and the
    "pulse never reaches the ground" limit — transfer directly, and are why those studies
    are in the body of this review rather than an appendix.

## Gaps / not found

1. **No peer-reviewed comparison of per-cell order statistics (min / low percentile / median /
   mean) for a bare-earth forest DTM.** What exists compares *interpolation algorithms* (Boreggio
   2018; Montealegre 2015) or compares *scales* of a fixed minimum operator (Clark 2004; Ewald 2013).
   The statistic itself is never the treatment. Probed from five phrasings.

2. **No use of a percentile above the median for ground elevation, anywhere.** High percentiles
   appear only for canopy/crop *surface* models. Medeiros' quartile *offset* is the nearest relative.
   *Do not read this as an opening:* we tested it and the column is symmetric at every cover level,
   so a shifted percentile has nothing to work with (commit 5335359).

3. **No study fixes cell size and tunes the percentile.** All tuning precedent is in scale/window,
   conflating the quantile with the support. Same caveat as (2) — the decoupled version is a cleaner
   experiment than anything published, and we ran it, and it came back negative. That negative, with
   its mechanism, is itself the publishable item.

4. **No epoch-matching study adjusts the ground statistic.** The two that exist (DeLong 2022;
   Viedma 2022) fit spatial correction surfaces after gridding. Cobby (2001) switches the rule by
   vegetation class but within a single epoch.

5. **Nobody quantifies leaf-state bias as a term in a DoD error budget** — including DeLong et al.
   (2022), who had a spring/autumn epoch pair over Minnesota forest. Arguably the clearest opening
   for our contribution.

6. **Very few DTM accuracy studies in temperate deciduous forest.** Simpson et al. say so
   themselves: "very few studies have formally assessed how vertical vegetation structure can affect
   DTM accuracy in broadleaf forests." Their Table 1 has two temperate-deciduous rows against six
   temperate-conifer rows. Su & Bork (2006) and Simpson et al. (2017) are the only close analogues
   in this review, and both are small.

7. **Near-nadir-worst error on slopes remains essentially undescribed.** Su & Bork observed it and
   discounted it as outlier-driven; their abstract concludes sampling angle "had little impact." I
   found no other treatment. Still likely novel.

8. **Cover is rarely quantified continuously.** Most studies use categorical classes (Reutebuch's
   four thinning treatments; Hopkinson's DU classes; Hladik's ten species classes). Continuous,
   data-derived cover appears in Salleh (non-ground return fraction), Buffington (NDVI), Simpson
   (Pgap) and Duchan (forest height + ground reflection density) — four papers. A regression of
   signed shift on continuous cover, in temperate deciduous forest, is not in the literature I found.

### Sources I could not verify

| Source | Why | What NOT to claim |
|---|---|---|
| **Cobby et al. (2001)** *ISPRS JPRS* 56(2):121–138 | Abstract null on OpenAlex; full text not obtained | Anything beyond Clark et al.'s description. **Highest-value acquisition on this list** — it is the one study that makes the ground statistic a function of vegetation class |
| Hodgson et al. (2003) *RSE* 84:295–308 | Abstract withheld; null on Crossref, OpenAlex, Semantic Scholar; ScienceDirect 403 | Either attribution of the 1.22 / 1.53 m figures — Clark and Simpson **disagree** on which cover class each belongs to |
| Kraus & Pfeifer (1998) *ISPRS JPRS* 53:193–203 | No full text or verbatim abstract | Any specific number; the ±10 cm figure comes only via Kraus & Rieger (1999) |
| Töyrä et al. (2003) | Hopkinson cites a 2003 work; OpenAlex returns a 2001 *RSE* paper by the same authors with no abstract | Do not assume the 2001 record is the source of the +0.07/+0.15 m figures |
| Buffington et al. (2016) "118 points" | Search summary only; not in any verbatim source | The 118-point number |
| Aryal et al. (2017) *PFG* 85:243–255 | Springer paywall + IdP redirect; abstract elided everywhere | Anything. Highly relevant (temperate forest, slope × aspect × habitat) |
| Su & Bork (2006) "23 times" and "(10°)" | `<`/`>` glyphs systematically absent from the PDF text layer (zero in 10 pages) | Quote only after checking the printed page; the missing `>` is near-certain, "23 times" appears literal |
| Stereńczak & Kozak (2011) pages 374–384 | From a listing, not the article | Page numbers only |

### Suggested next acquisitions (library access), in priority order

1. **Cobby, Mason & Davenport (2001)**, *ISPRS JPRS* 56(2):121–138, doi:10.1016/S0924-2716(01)00039-9 — vegetation-class-dependent ground retrieval in leaf-on deciduous forest. The direct methodological ancestor.
2. **Hodgson, Jensen, Schmidt, Schill & Davis (2003)**, *RSE* 84(2):295–308 — resolve the conflicting attributions and get the leaf-on temperate deciduous numbers first-hand.
3. **Aryal, Latifi, Heurich & Hahn (2017)**, *PFG* 85(4):243–255, doi:10.1007/s41064-017-0023-2 — temperate forest, slope × aspect × habitat.
4. **Kraus & Pfeifer (1998)**, *ISPRS JPRS* 53(4):193–203 — the asymmetric weight function in original form; the theoretical justification for a sub-median statistic.
5. **Schmid, Hadley & Wijekoon (2011)**, *JCR* 27(6A):116–132 — min-bin implementation details and per-species cell-size optima.
6. **Töyrä, Pietroniro & Hopkinson (2003)** — locate the correct record.

---

## Full reference list

Peer-reviewed journal articles unless marked otherwise.

- Boreggio, M., Bernard, M., & Gregoretti, C. (2018). Evaluating the Differences of Gridding Techniques for Digital Elevation Models Generation and Their Influence on the Modeling of Stony Debris Flows Routing: A Case Study From Rovina di Cancia Basin (North-Eastern Italian Alps). *Frontiers in Earth Science*, 6, 89. doi:10.3389/feart.2018.00089 — **[ABSTRACT + body]**
- Buffington, K.J., Dugger, B.D., Thorne, K.M., & Takekawa, J.Y. (2016). Statistical correction of lidar-derived digital elevation models with multispectral airborne imagery in tidal marshes. *Remote Sensing of Environment*, 186, 616–625. doi:10.1016/j.rse.2016.09.020 — **[ABSTRACT, via NOAA NCCOS record]**
- Clark, M.L., Clark, D.B., & Roberts, D.A. (2004). Small-footprint lidar estimation of sub-canopy elevation and tree height in a tropical rain forest landscape. *Remote Sensing of Environment*, 91(1), 68–89. doi:10.1016/j.rse.2004.02.008 — **[FULL TEXT]**
- Cobby, D.M., Mason, D.C., & Davenport, I.J. (2001). Image processing of airborne scanning laser altimetry data for improved river flood modelling. *ISPRS Journal of Photogrammetry and Remote Sensing*, 56(2), 121–138. doi:10.1016/S0924-2716(01)00039-9 — **[NOT VERIFIED; described SECONDARY via Clark et al. 2004]**
- DeLong, S.B., Hammer, M.N., Engle, Z.T., Richard, E.M., Breckenridge, A.J., Gran, K.B., Jennings, C.E., & Jalobeanu, A. (2022). Regional-Scale Landscape Response to an Extreme Precipitation Event From Repeat Lidar and Object-Based Image Analysis. *Earth and Space Science*, 9(12), e2022EA002420. doi:10.1029/2022EA002420 — **[ABSTRACT + body]**
- Duchan, M., Mráz, V., Tichá, A., Jankovský, M., & Zlatuška, K. (2026). The Influence of Forest Cover on the Accuracy of Aerial Laser Scanning-Derived Digital Elevation Models for Detecting Drainage Ditches in Forests in the Czech Republic. *Forests*, 17(2), 162. doi:10.3390/f17020162 — **[FULL TEXT]**
- Evans, J.S., & Hudak, A.T. (2007). A Multiscale Curvature Algorithm for Classifying Discrete Return LiDAR in Forested Environments. *IEEE Transactions on Geoscience and Remote Sensing*, 45(4), 1029–1038. doi:10.1109/TGRS.2006.890412 — **[ABSTRACT]** (classifier context)
- Ewald, M.J. (2013). *Where's the Ground Surface? Elevation Bias in LIDAR-derived Digital Elevation Models Due to Dense Vegetation in Oregon Tidal Marshes.* **MS thesis (not peer-reviewed)**, Oregon State University. <https://ir.library.oregonstate.edu/downloads/1n79h8198> — **[FULL TEXT]**
- Fradette, M.-S., Leboeuf, A., Riopel, M., & Bégin, J. (2019). Method to Reduce the Bias on Digital Terrain Model and Canopy Height Model from LiDAR Data. *Remote Sensing*, 11(7), 863. doi:10.3390/rs11070863 — **[ABSTRACT]**
- Hladik, C., & Alber, M. (2012). Accuracy assessment and correction of a LIDAR-derived salt marsh digital elevation model. *Remote Sensing of Environment*, 121, 224–235. doi:10.1016/j.rse.2012.01.018 — **[ABSTRACT]**
- Hodgson, M.E., & Bresnahan, P. (2004). Accuracy of Airborne Lidar-Derived Elevation: Empirical Assessment and Error Budget. *Photogrammetric Engineering and Remote Sensing*, 70(3), 331–339. doi:10.14358/PERS.70.3.331 — **[ABSTRACT]**; signed-error range **[SECONDARY via Hopkinson et al. 2005]**
- Hodgson, M.E., Jensen, J.R., Schmidt, L., Schill, S., & Davis, B. (2003). An evaluation of LIDAR- and IFSAR-derived digital elevation models in leaf-on conditions with USGS Level 1 and Level 2 DEMs. *Remote Sensing of Environment*, 84(2), 295–308. doi:10.1016/S0034-4257(02)00114-1 — **[NOT VERIFIED; conflicting SECONDARY descriptions]**
- Hopkinson, C., Chasmer, L., Sass, G.Z., Creed, I.F., Sitar, M., Kalbfleisch, W., & Treitz, P. (2005). Vegetation class dependent errors in lidar ground elevation and canopy height estimates in a boreal wetland environment. *Canadian Journal of Remote Sensing*, 31(2), 191–206. doi:10.5589/m05-007 — **[ABSTRACT]**; conference version (ISPRS Archives XXXVI-8/W2, **not peer-reviewed**) — **[FULL TEXT]**
- Kraus, K., & Pfeifer, N. (1998). Determination of terrain models in wooded areas with airborne laser scanner data. *ISPRS Journal of Photogrammetry and Remote Sensing*, 53(4), 193–203. doi:10.1016/S0924-2716(98)00009-4 — **[NOT VERIFIED]**
- Kraus, K., & Rieger, W. (1999). Processing of laser scanning data for wooded areas. In D. Fritsch & R. Spiller (Eds.), *Photogrammetric Week '99*, 221–231. Wichmann, Heidelberg. **Conference proceedings, not peer-reviewed.** <https://phowo.ifp.uni-stuttgart.de/publications/phowo99/kraus.pdf> — **[FULL TEXT]**
- Medeiros, S., Hagen, S., Weishampel, J., & Angelo, J. (2015). Adjusting Lidar-Derived Digital Terrain Models in Coastal Marshes Based on Estimated Aboveground Biomass Density. *Remote Sensing*, 7(4), 3507–3525. doi:10.3390/rs70403507 — **[ABSTRACT]**
- Montealegre, A.L., Lamelas, M.T., & de la Riva, J. (2015). Interpolation Routines Assessment in ALS-Derived Digital Elevation Models for Forestry Applications. *Remote Sensing*, 7(7), 8631–8654. doi:10.3390/rs70708631 — **[ABSTRACT]**
- Reutebuch, S.E., McGaughey, R.J., Andersen, H.-E., & Carson, W.W. (2003). Accuracy of a high-resolution lidar terrain model under a conifer forest canopy. *Canadian Journal of Remote Sensing*, 29(5), 527–535. doi:10.5589/m03-022 — **[FULL TEXT]**
- Salleh, M.R.M., Ismail, Z., & Abdul Rahman, M.Z. (2015). Accuracy Assessment of Lidar-Derived Digital Terrain Model (DTM) with Different Slope and Canopy Cover in Tropical Forest Region. *ISPRS Annals of the Photogrammetry, Remote Sensing and Spatial Information Sciences*, II-2/W2, 183–189. **Conference annals; full-paper double-blind peer-reviewed per the paper itself.** doi:10.5194/isprsannals-II-2-W2-183-2015 — **[FULL TEXT]**
- Schmid, K.A., Hadley, B.C., & Wijekoon, N. (2011). Vertical Accuracy and Use of Topographic LIDAR Data in Coastal Marshes. *Journal of Coastal Research*, 27(6A), 116–132. doi:10.2112/JCOASTRES-D-10-00188.1 — **[ABSTRACT]**
- Simpson, J.E., Smith, T.E.L., & Wooster, M.J. (2017). Assessment of Errors Caused by Forest Vegetation Structure in Airborne LiDAR-Derived DTMs. *Remote Sensing*, 9(11), 1101. doi:10.3390/rs9111101 — **[FULL TEXT]**
- Stereńczak, K., & Kozak, J. (2011). Evaluation of digital terrain models generated in forest conditions from airborne laser scanning data acquired in two seasons. *Scandinavian Journal of Forest Research*, 26(4), 374–384. doi:10.1080/02827581.2011.570781 — **[ABSTRACT]**
- Su, J., & Bork, E. (2006). Influence of Vegetation, Slope, and Lidar Sampling Angle on DEM Accuracy. *Photogrammetric Engineering and Remote Sensing*, 72(11), 1265–1274. doi:10.14358/PERS.72.11.1265 — **[FULL TEXT]**
- Tinkham, W.T., Smith, A.M.S., Hoffman, C., Hudak, A.T., Falkowski, M.J., Swanson, M.E., & Gessler, P.E. (2012). Investigating the influence of LiDAR ground surface errors on the utility of derived forest inventories. *Canadian Journal of Forest Research*, 42(3), 413–422. doi:10.1139/x11-193 — **[ABSTRACT]**
- Töyrä, J., Pietroniro, A., & Martz, L.W. (2001). Multisensor Hydrologic Assessment of a Freshwater Wetland. *Remote Sensing of Environment*, 75(2), 162–173. doi:10.1016/S0034-4257(00)00164-4 — **[NOT VERIFIED]**; the +0.07/+0.15 m figures are **[SECONDARY via Hopkinson et al. 2005]** and are attributed there to a 2003 Töyrä work
- U.S. Geological Survey. *Lidar Base Specification 2025 rev. A* — Digital Elevation Model Surface Treatments. **Agency standard, not peer-reviewed.** <https://www.usgs.gov/ngp-standards-and-specifications/lidar-base-specification-digital-elevation-model-surface> — **[partially verified]**
- Viedma, O. (2022). Applying a Robust Empirical Method for Comparing Repeated LiDAR Data with Different Point Density. *Forests*, 13(3), 380. doi:10.3390/f13030380 — **[FULL TEXT]**
- Wang, C., Menenti, M., Stoll, M.P., Feola, A., Belluco, E., & Marani, M. (2009). Separation of Ground and Low Vegetation Signatures in LiDAR Measurements of Salt-Marsh Environments. *IEEE Transactions on Geoscience and Remote Sensing*, 47(7), 2014–2023. doi:10.1109/TGRS.2008.2010490 — **[ABSTRACT]**
- Wasser, L., Day, R., Chasmer, L., & Taylor, A. (2013). Influence of Vegetation Structure on Lidar-derived Canopy Height and Fractional Cover in Forested Riparian Buffers During Leaf-Off and Leaf-On Conditions. *PLoS ONE*, 8(1), e54776. doi:10.1371/journal.pone.0054776 — **[FULL TEXT, targeted read]**

# CSF parameter selection: what the literature actually supports

**Scope.** Evidence for choosing Cloth Simulation Filter (CSF) parameters —
`class_threshold` (PDAL `threshold`), `cloth_resolution` (PDAL `resolution`),
`rigidness`, `time_step` (PDAL `step`), and slope post-processing
(`smooth` / `sloop_smooth` / "steep slope fit") — for **driftless-area bluffland,
SE Minnesota**: broad low-relief upland divides dissected by steep wooded valley
walls; ~40% of cells 0-5°, 18% 5-10°, 18% 10-20°, 19% 20-30°, 6% >30°; row-crop
uplands, deciduous + conifer valley walls (canopy cover 0 to >0.65); 2008-era
Optech ALTM Gemini / Leica ALS50-II, **leaf-off November**, ~1-2 pulses/m²,
ground-return density **~0.5 pts/m²**; gridded at 5 m with a per-cell median.

**Current pipeline setting** (re-verified against `src/lidar_diff_icp/ground.py`
at commit `a0f7321`, "Pass only the overridden CSF parameter; inherit the rest
from PDAL"): the module now passes **only `rigidness=1`** and leaves every other
CSF option unset, so the rest come from PDAL itself. Its docstring, verbatim:

> "**We override exactly one CSF parameter: ``rigidness=1``.** PDAL ships 3, the
> hard cloth documented for FLAT terrain; this is dissected bluffland. Every other
> parameter is left unset and therefore comes from PDAL itself, so we inherit its
> defaults rather than holding stale copies of them."

The **effective** values are therefore `rigidness=1` plus PDAL's defaults —
`resolution=1.0, threshold=0.5, hdiff=0.3, smooth=true, step=0.65,
iterations=500`, with `returns="last, only"`. That makes the "PDAL `filters.csf`"
row of the table in §0 the operative one for everything except rigidness.

**Headline.** The evidence is strong on three parameters and weak on the two that
matter most. `time_step = 0.65` and `class_threshold = 0.5` are well supported
and should not move; slope post-processing on is what the authors prescribe for
steep terrain. **`cloth_resolution` is the dominant parameter** (Klápště et al.
2020, explicit ranking; LiDAR RMSE spans 0.18-0.53 m across their tested grid) and
the bias-minimising direction is **coarser, not finer** — so keep 1.0 and test
2.0, not 0.5. **`rigidness = 1` is probably too soft** for a tile that is 40% flat;
Zhang's own prescription for river-bank/ditch/terrace relief is `rigidness = 2`
with post-processing on. And **CSF reads the ground high** (+0.12-0.13 m vs GNSS,
~3× SMRF's bias), with the bias growing on slopes. Full reasoning in §7; what is
missing in §8.

**Verification convention used throughout.**
- *full text read* — I downloaded the PDF/HTML and read the quoted sentence in it.
- *abstract only* — I read the publisher's abstract but not the full text.
- *NOT VERIFIED* — claim comes from a search snippet or secondary citation only.

Quoted material is verbatim. Where I reason beyond the sources, it is labelled
**[inference, not from the literature]**. Sources marked *(delegated)* were read
in full by a parallel search agent rather than by me; the load-bearing ones
(Zhang 2016, Klápště 2020, Štular & Lozić 2020, Storch 2022, Viedma 2022, Wan
2018, Cai 2019, Cai 2023) I downloaded and read myself.

---

## 0. The defaults disagree across implementations — and none of them is "the paper's"

This matters before anything else, because "we used CSF defaults" is not a
well-defined statement.

| | `rigidness` | `cloth_resolution` | `class_threshold` | `time_step` | slope post-proc | `hdiff` |
|---|---|---|---|---|---|---|
| **Zhang et al. 2016 (paper)** | **scene-dependent 1/2/3** | 0.5 (fixed) | 0.5 (fixed) | 0.65 (fixed) | **scene-dependent** | 0.3 (fixed) |
| **Upstream C++ lib** (`jianboqi/CSF`, `src/CSF.cpp`) | **3** | **1** | 0.5 | 0.65 | **true** | not exposed (hard-coded 0.3) |
| **PDAL `filters.csf`** | **3** | **1.0** | 0.5 | 0.65 | **true** | 0.3 |
| **RCSF / lidR `csf()`** | **1** | **0.5** | 0.5 | 0.65 | **FALSE** | not exposed |
| **CloudCompare `qCSF`** | **2** ("Relief") | **2.0** | 0.5 | not exposed | **false** | not exposed |

Sources, all *full text read*:

- **Zhang, W.; Qi, J.; Wan, P.; Wang, H.; Xie, D.; Wang, X.; Yan, G. (2016).**
  "An Easy-to-Use Airborne LiDAR Data Filtering Method Based on Cloth
  Simulation." *Remote Sensing* 8(6):501. DOI:
  [10.3390/rs8060501](https://doi.org/10.3390/rs8060501). Peer-reviewed, open
  access (CC-BY). PDF obtained from
  `https://mdpi-res.com/d_attachment/remotesensing/remotesensing-08-00501/article_deploy/remotesensing-08-00501.pdf`
  (www.mdpi.com 403s to both WebFetch and curl).
- **Upstream library**: `https://github.com/jianboqi/CSF`, `src/CSF.cpp` lines
  29-34 (*software source, not peer-reviewed*):
  > `params.bSloopSmooth = true;` / `params.time_step = 0.65;` /
  > `params.class_threshold = 0.5;` / `params.cloth_resolution = 1;` /
  > `params.rigidness = 3;` / `params.interations = 500;`
- **PDAL** `filters.csf` documentation (*software docs*), fetched with curl
  (pdal.io 403s WebFetch):
  > "resolution Cloth resolution. [Default: 1.0] … threshold Classification
  > threshold. [Default: 0.5] hdiff Height difference threshold. [Default: 0.3]
  > smooth Perform slope post-processing? [Default: true] step Time step.
  > [Default: 0.65] rigidness Rigidness. [Default: 3] iterations Maximum number
  > of iterations. [Default: 500]"
  >
  > Also: "returns Return types to include in output. Valid values are "first",
  > "last", "intermediate" and "only". [Default: "last, only"]"
- **RCSF** (`Jean-Romain/RCSF`, `R/CSF.R`) and **lidR** (`r-lidar/lidR`,
  `R/algorithm-gnd.R` line 196), identical signature (*software docs*):
  > `csf = function(sloop_smooth = FALSE, class_threshold = 0.5, cloth_resolution = 0.5, rigidness = 1L, iterations = 500L, time_step = 0.65)`
- **CloudCompare** `qCSF` (`plugins/core/Standard/qCSF/src/qCSF.cpp`, lines
  121-125) (*software source*):
  > `static bool PostProcessing = false;` / `static double ClothResolution = 2.0;` /
  > `static double ClassThreshold = 0.5;` / `static int Rigidness = 2;` /
  > `static int MaxIteration = 500;`
  >
  > The dialog (`ui/CSFDlg.ui`) labels the three rigidness radio buttons
  > `Steep slope` (rig1), `Relief` (rig2), `Flat` (rig3), with a separate
  > `Slope processing` checkbox.

**Two documentation errors worth knowing about.**

1. The lidR book states (*software docs*, `https://r-lidar.github.io/lidRbook/gnd.html`,
   *full text read*):
   > "The `csf()` functions use the default values proposed by Zhang et al 2016
   > and can be used without providing any arguments."

   This is **not accurate**. Zhang et al. do not propose a single default
   `rigidness`; they propose a scene-dependent 1/2/3 (their Table 2), and they
   turn slope post-processing on for two of their three scene groups. lidR's
   `rigidness = 1L, sloop_smooth = FALSE` corresponds to Zhang's group III
   **with post-processing switched off** — a combination Zhang et al. never used.
   (`cloth_resolution = 0.5`, `class_threshold = 0.5`, `time_step = 0.65` *do*
   match the paper.)

2. The often-repeated guidance that cloth resolution should equal the point
   spacing traces to the **RCSF/lidR wrapper documentation only** — see §4. It
   is not in the paper, and it is not in the CSF authors' own plugin
   documentation.

**One thing that is *not* a difference, checked because it looked like one.**
PDAL exposes `hdiff` where upstream has no such option, and `CSFilter.cpp` line
250 assigns it to a struct field called `height_threshold`, which reads alarming.
It is fine. Both implementations construct the cloth the same way — upstream
`src/CSF.cpp` lines 137-147 passes `0.3` as the sixth argument (`_smoothThreshold`)
and `9999` as the seventh (`_heightThreshold`); PDAL's vendored
`filters/private/csf/CSF.cpp` lines 150-162 passes `params.height_threshold`
(= `hdiff`, default 0.3) as the sixth and `9999` as the seventh. **PDAL simply
exposes upstream's hard-coded 0.3 as a tunable, and the default behaviour is
identical.** The field name is misleading, not the behaviour. (This 0.3 is
Zhang's h_cp, used only when slope post-processing is on.)

**Consequence for us.** Our effective `resolution=1.0` is PDAL's default, which is the
upstream C++ library's default — **not** the value Zhang et al. tested and
recommended (0.5). Our `smooth=True` is PDAL's default and *is* what Zhang et
al. used for steep scenes. Our `rigidness=1` matches Zhang's group III ("high
and steep slopes (e.g., pit, cliff)").

---

## 1. The original authors' guidance (Zhang et al. 2016)

*Peer-reviewed, open access. Full text read.*

### 1.1 What the parameters are

> "CSF mainly consists of four user-defined parameters: grid resolution (GR),
> which represents the horizontal distance between two neighboring particles;
> time step (dT), which controls the displacement of particles from gravity
> during each iteration; rigidness (RI), which controls the rigidness of the
> cloth; and an optional parameter steep slope fit factor (ST), which indicates
> whether the post-processing of handling steep slopes is required or not."

> "In addition to these user-defined parameters, two threshold parameters have
> been used in this algorithm to aid the identification of ground points. The
> first is a distance threshold (h_cc) that governs the final classification of
> the LiDAR points as BE and OBJ based on the distances to the cloth grid. This
> parameter is set as a fixed value of 0.5 m. Another threshold parameter is the
> height difference (h_cp), which is used during post-processing to determine
> whether a movable particle should be moved to the ground or not. This
> parameter is set to 0.3 m for all of the datasets."

Mapping to PDAL: `h_cc` = `threshold`, `h_cp` = `hdiff`, `GR` = `resolution`,
`dT` = `step`, `ST` = `smooth`.

### 1.2 Rigidness: what it physically is, and the recommendation

> "If RI is set to 1, the movable particle is just moved only once, and the
> displacement is half of the vertical distance (VD) between the two particles.
> If the RI is set to 2, the movable particle will be moved twice, the total
> displacement is 3/4VD. Finally, if RI is set to 3, the movable particle will
> be moved three times and the total displacement is 7/8VD. The value of 3 is
> enough to produce a very hard cloth. Thus, we constrain the rigidness to
> values of 1, 2 and 3. The larger the rigidness is, the more rigidly the cloth
> will behave."

**This is the core scene→parameter recommendation** (their §3.1):

> "If the terrain is very flat and has no steep or terraced slopes, RI is set to
> a relatively large value (RI = 3), and no post-processing is needed (ST =
> false). If steep slopes exist (e.g., river bank, ditch, and terrace), a medium
> soft cloth (RI = 2) and post-processing (ST = true) are needed. When handling
> very steep slopes, we need post-processing (ST = true) and a very soft cloth
> (RI = 1)."

Restated in their Discussion §4.2:

> "Usually, RI can be set to 1, 2 or 3 according to the features of the terrain,
> they are applied to areas with high steep slopes, terraced slopes and gentle
> slopes, respectively. ST is set to "true" or "false". "true" means that there
> exists steep slopes and post-processing is needed. "false" means
> post-processing is not needed."

Their Table 2, verbatim (with `dT = 0.65, GR = 0.5` for all groups):

| Group | Feature | Parameters | Samples |
|---|---|---|---|
| I | Flat terrain or gentle slope, no steep slopes | RI = 3, ST = false | 21, 31, 42, 51, 54 |
| II | With steep or terraced slopes (e.g., river bank, ditch, terrace) | RI = 2, ST = true | 11, 12, 22, 23, 24, 41 |
| III | High and steep slopes (e.g., pit, cliff) | RI = 1, ST = true | 52, 53, 61, 71 |

Note the explicit **error trade-off** they attach to rigidness:

> "However, if a harder cloth is used it may yield the opposite error (BE points
> around the steep slopes may be identified as OBJ). Thus, adjusting the
> parameters is necessary to balance type I and type II errors. If a large
> number of low objects exists above the ground, the cloth should be harder (RI
> should be set larger), which will guarantee that fewer object measurements are
> mistakenly classified as ground objects."

Soft cloth → fewer true-ground losses on slopes (low type I) but more low
vegetation retained as ground (high type II). Hard cloth → the reverse.

### 1.3 Slope post-processing — what it actually does

> "For steep slopes, this algorithm may yield relatively large errors because the
> simulated cloth is above the steep slopes and does not fit with the ground
> measurements very well due to the internal constraints among particles … Some
> ground measurements around steep slopes are mistakenly classified as OBJ. This
> problem can be solved by a post-processing method that smoothes the margins of
> steep slopes. This post-processing method finds an unmovable particle in the
> four adjacent neighborhoods of each movable particle and compares the height
> values of CPs. If the height difference is within a threshold (h_cp), the
> movable particle is moved to the ground and set as unmovable."

And the reason to prefer post-processing over simply softening the cloth (§4.3):

> "A direct method to mitigate this problem is to set the rigidness to a lower
> value, but some low objects may be classified as BE as a result. To balance
> these two types of errors, a post-processing method for the margin area is
> proposed in this study … Thus, we can use a relatively hard cloth and
> post-processing to remove lower objects and correctly handle steep slope
> areas."

**This is the single most decision-relevant sentence in the paper for us**: the
authors' preferred route to handling steep slopes is *post-processing plus a
relatively hard cloth*, not a soft cloth alone.

Known failure mode of post-processing, in their own results (samp11):

> "samp11 has a large group of error points (type II) that almost classifies a
> whole building as ground measurements, it occurs because the building is
> located on a slope and the roof is nearly connected to the ground, causing the
> building to be treated as ground in the post-processing step."

(Bridges are the other documented failure, §4.4. Neither is a concern for a
rural bluffland tile, though it is a warning for any structure on a hillslope.)

### 1.4 `cloth_resolution` and `time_step`: what was tested

> "The main parameters that control the results and vary with the scene type
> were RI and ST … For dT and GR, we set them as fixed values of 0.65 and 0.5,
> respectively. These two values are universally applicable to all of the
> reference datasets according to our tests."

Time step sweep (their §4.2 and Figure 13):

> "we tested all of the samples with different time steps (from 0.4 to 1.5 with
> steps of 0.05; 0.4 was chosen because the value would take too much time to
> compute when it was smaller than 0.4). … This figure indicates that the total
> error increases after an initial decline for all groups, and all of them
> achieve the lowest total error around the 0.65 time step."

Grid resolution sweep (their §4.2 and Figure 14):

> "The grid resolution (GR) parameter in the simulation process has strong
> relationship with simulation time because it determines how many cloth
> particles are created for a specific dataset. Figure 14 shows the total errors
> at different GR values. It can be seen that the accuracies of group I and
> group III are relatively stable than group II because group II usually have
> complicated terrain shape and buildings (e.g., areas with terraced slopes and
> low rise buildings). However, almost all samples get the highest accuracy
> around 0.5, which was then been used as a fixed value."

Classification threshold (their §4.2 and Figure 15):

> "h_cc governs the final classification which separates LiDAR measurements into
> BE or OBJ. Most particles will stick to ground after simulation. And OBJ
> measurements (e.g., buildings and trees) are usually taller than 0.5 m. Thus,
> we set h_cc as 0.5, which is also a fixed value. The influences of h_cc on
> total errors are illustrated in Figure 15. It shows that this value has
> limited impact on total errors."

Post-processing threshold `h_cp`:

> "We simply set this parameter to 0.3 m, which indicates the height difference
> between two adjacent ground measurements is usually less than 0.3 m on a flat
> terrain. Since this parameter is only used when post-process is enabled, and
> it also only influence the movable particles over steep slopes, the influence
> is also limited."

Iterations:

> "Usually, this user-defined value is set to 500. However in most cases, CSF
> will end according to the former criterial. … both M_HV and A_HV decreased to
> a very low value around 100 iterations. Actually, M_HV is 0.0033 when
> iteration number equal to 150."

### 1.5 ISPRS benchmark results — including the steep and vegetated samples

Reported total error (T.E.), type I (T.I, true ground lost), type II (T.II,
object retained as ground), and Kappa, from their Table 3. The **steep/vegetated
samples are the ones that matter for us**, marked below:

| Sample | Feature (their Table 1) | Group | T.I(%) | T.II(%) | T.E.(%) | Kappa(%) |
|---|---|---|---|---|---|---|
| samp11 | Mixture of vegetation and buildings on hillside ★ | II | 7.23 | 18.44 | **12.01** | 75.17 |
| samp12 | Buildings on hillside ★ | II | 1.15 | 4.9 | 2.97 | 94.04 |
| samp21 | Large buildings and bridge | I | 3.89 | 1.78 | 3.42 | 90.47 |
| samp22 | Irregularly shaped buildings | II | 1.29 | 25.9 | 8.94 | 77.72 |
| samp23 | Large, irregularly shaped buildings | II | 3.52 | 6.21 | 4.79 | 90.38 |
| samp24 | Steep slopes ★ | II | 1.03 | 7.73 | 2.87 | 92.68 |
| samp31 | Complex buildings | I | 0.96 | 2.38 | 1.61 | 96.75 |
| samp41 | Data gaps | II | 1.48 | 8.78 | 5.14 | 89.73 |
| samp42 | Railway station with trains | I | 3.28 | 0.87 | 1.58 | 96.18 |
| samp51 | Mixture of vegetation and buildings on hillside ★ | I | 2.67 | 4.57 | 3.08 | 91.13 |
| samp52 | Buildings on hillside ★ | III | 1.01 | 28.79 | 3.93 | 77.05 |
| samp53 | Large buildings and bridge | III | 3.85 | 37.08 | 5.2 | **46.86** |
| samp54 | Irregularly shaped buildings | I | 3.79 | 2.64 | 3.18 | 93.61 |
| samp61 | Large, irregularly shaped buildings | III | 0.87 | 18.94 | 1.49 | 78.1 |
| samp71 | Steep slopes ★ | III | 1.61 | 37.85 | 5.71 | 68.03 |

Their own summary of where CSF is weak:

> "For group II and group III, the total errors are relatively large (especially
> for samp11) compared to group I, which shows that CSF performs relatively
> poorly in complex regions similar to other filtering algorithms … However, in
> high relief areas with very steep slopes (e.g., pit, cliff) and low rise
> buildings, our method perform worst, because when a soft cloth fit with the
> terrain, it may also reach the rooftops of low rise buildings."

Comparative standing (their Table 4, average T.E. across all 15 samples): CSF
**4.39**, vs Pingel 2013 (SMRF) **2.97**, Hu 2014 **2.85**, Mongus 2014 **2.74**,
Chen 2013 **4.11**, Axelsson 1999 (TIN densification) **4.82**, Zhang 2013
**10.63**. So CSF is mid-pack, and *worse than SMRF* on this benchmark.

> "the mean total error (4.39) and standard deviation (2.76) of all the samples
> are relatively low compared to all of the other algorithms"

**Internal inconsistency, flagged**: the abstract of the same paper says
> "the experimental results yield an average total error of 4.58%"

while Table 3/Table 4 give a mean of 4.39% (I recomputed the mean of the Table 3
T.E. column: 4.395%). Cite 4.39% and note the abstract's figure is
irreconcilable with the tables.

### 1.6 The dense-cloud test (their §3.2) — the closest thing to a forested-hillslope test

> "we tested the performance of CSF with datasets that have more dense points
> with average point distance equal to 0.6 m–0.8 m"

| Dataset | Type | Points | Scope | Features | T.I(%) | T.II(%) | T.E.(%) |
|---|---|---|---|---|---|---|---|
| 1 | Urban | 1,559,933 | 1×1 km | Flat terrain, large and dense buildings, high vegetation coverage | 0.72 | 13.36 | 6.84 |
| 2 | Urban | 1,522,256 | 1×1 km | Flat terrain with dense bungalow areas | 5.29 | 9.29 | 7.84 |
| 3 | Rural | 2,093,506 | 2×1 km | dense vegetation coverage | **36.09** | 1.84 | 5.49 |
| 4 | Rural | 1,418,228 | 0.5×0.5 km | Large number of steep slopes | 8.57 | 22.61 | **14.09** |

> "A large T.I error also has been noted for dataset 3, this is because ground
> measurements is very sparse in this area. For dataset 4, the error mainly
> occurs around steep slopes, since this area contains large number of steep
> slopes."

> "In mountain areas, CSF performs relatively poorly, especially in dense
> vegetation areas where ground measurements are usually sparse. If the cloth is
> too soft, many object measurements may be mistakenly classified as BE.
> Otherwise, ground measurements may be classified as OBJ due to the hilly
> topography (see Figure 11). For areas with large number of steep slopes, cloth
> should be more soft and post-process is also needed."

**This is the most directly transferable result in the paper.** Dataset 3
(dense vegetation, sparse ground returns) loses **36% of true ground points**
(type I). Dataset 4 (many steep slopes) has the worst total error of the four
at 14.09%. Our tile is the *combination* of those two conditions.

### 1.7 The authors' own plugin documentation (not the paper)

CloudCompare wiki page "CSF (plugin)" (last edited 20 July 2016), text
attributed to Zhang, Qi, Wan and Wang (*software documentation, not
peer-reviewed*; `https://www.cloudcompare.org/doc/wiki/index.php/CSF_(plugin)`,
*full text read*):

> "**Scenes** — Three options are under this parameter : Steep slope, Relief, and
> Flat. This parameter help users to set scenes type of the point clouds. When
> you set up this parameter, the rigidness will be determined actually."

> "**Slope post processing for disconnected terrain** — For steep slopes, this
> algorithm may yield relatively large errors because the simulated cloth is
> above the steep slopes and does not fit with the ground measurements very well
> due to the internal constraints among particles. This problem can be solved by
> selecting this option. If there are no steep slopes in your scenes, just
> neglect it."

> "**Cloth resolution** — Cloth resolution refers to the grid size (the unit is
> same as the unit of pointclouds) of cloth which is used to cover the terrain.
> The bigger cloth resolution you have set, the coarser DTM you will get."

> "**Max iterations** — … 500 is enough for most of scenes."

> "**Classification threshold** — Classification threshold refers to a threshold
> … to classify the pointclouds into ground and non-ground parts based on the
> distances between points and the simulated terrain. 0.5 is adapted to most of
> scenes."

Note what is **absent**: no guidance tying cloth resolution to point spacing, and
no numeric recommendation for cloth resolution at all.

### 1.8 A caveat on "the CSF in your software is the CSF in the paper"

The upstream repository README states (*software docs*, `jianboqi/CSF`,
*full text read*):

> "Note: This code has been changed a lot since the publication of the
> corresponding paper. A lot of optimizations have been made. We are still
> working on it, and wish it could be better."

Neither PDAL nor CloudCompare nor RCSF documents which library revision it
vendors. Treat the paper's benchmark numbers as indicative of the *method*, not
as a guarantee about the binary you are running.

---

## 2. Parameter sensitivity studies

### 2.1 Klápště et al. 2020 — the only full CSF parameter sweep against absolute ground truth

**Klápště, P.; Fogl, M.; Barták, V.; Gdulová, K.; Urban, R.; Moudrý, V. (2020).**
"Sensitivity analysis of parameters and contrasting performance of ground
filtering algorithms with UAV photogrammetry-based and LiDAR point clouds."
*International Journal of Digital Earth* 13(12):1672-1694. DOI:
[10.1080/17538947.2020.1791267](https://doi.org/10.1080/17538947.2020.1791267).
Peer-reviewed. *Full text read — via the `r.jina.ai` text proxy, because
tandfonline.com returns 403 to both WebFetch and curl. Body text renders in full;
**the numbered tables render as captions only**, so Table 2 (parameter grids) and
Table 4 (chosen best settings) are NOT available and no per-algorithm best-setting
value can be quoted verbatim.*

**This is the single most important source for our decision**, because it is the
only study that (a) sweeps CSF's parameters systematically, (b) evaluates them
against **absolute ground truth** rather than another lidar dataset, and (c)
reports a **signed** vertical error, stratified by slope and by vegetation class.

**Setting.** Un-reclaimed brown-coal spoil heap, NW Bohemia, Czech Republic
(50°34′N, 13°34′E); 30 ha test rectangle (550 × 550 m). Terrain "remained rugged
as a result of heaping that formed a typical undulated terrain". Vegetation in
four density classes — low grass (*Calamagrostis epigejos*, *Arrhenatherum
elatius*), shrub (0-1 m), medium (1-3 m), canopy (>3 m; deciduous *Betula*,
*Salix*, *Alnus*), plus *Phragmites*/*Typha* in depressions. **ALS: Riegl
LMS-Q780, May 2017, leaf-on, 1030 m AGL, 60° FOV, 55% sidelap, discrete return,
"an average point density of 8 points per square meter".** DTM at 0.5 m by
bin-average. Reference: **"In total, 1414 checkpoints were collected"** by RTK
GNSS at 4-6 cm vertical accuracy, "chosen in a way ensuring sufficient
representation of all vegetation categories and slope types".

**Transfer limitations, stated plainly:** 8 pts/m² is ~16× our density; the site
is a spoil heap with steppe, shrub and young deciduous woodland, not mature
bluffland forest; and its steep ground carries low vegetation rather than closed
canopy. The relief is "undulated", not a dissected plateau.

**Which parameter dominates — verbatim:**
> "The most influential parameters did not differ between LiDAR and photogrammetry
> point clouds. For all algorithms, the most influential were those related to
> selection of the initial minimum elevation ground points (i.e. *Cell size* for
> ATIN, PMF, SMRF, *Step size* for PTIN and **_Cloth resolution_ for CSF**). Other
> parameters had only minor effect and were important rather for fine tuning of
> ground filtering."

**How large the effect is — verbatim:**
> "The Cloth Simulation Filter (CSF) algorithm implemented in CloudCompare uses
> four main parameters: *Slope processing* (True/False), *Cloth resolution* (m),
> *Ruggedness* (1 – 'Mountain', 2 – 'Complex', or 3 – 'Flat'), and *Classification
> threshold* (m). **The greatest effect was observed for the Cloth resolution
> parameter**, with the value of 0.1 m leading to the best performances in terms
> of random-error and bias combination. It should be noted, however, that **for
> LiDAR point cloud, the parameter Cloth resolution 2 m in combination with either
> Slope processing False, or Slope processing True and Ruggedness 2 or 3, lead to
> almost zero bias, although the random error was relatively high.** The algorithm
> generally performed slightly better for the LiDAR data than for the UAV data.
> **The parameter Classification threshold was practically irrelevant.** The RMSE
> for all tested parameters ranged between **0.18–0.53 m for LiDAR point cloud**
> and 0.23–1.26 m for photogrammetric point cloud."

**Three things to take from that paragraph, and they are the crux of this whole
review.**

1. **`cloth_resolution` is the dominant CSF parameter** — by their explicit
   ranking, and it moves LiDAR RMSE across **0.18-0.53 m**, a factor of ~3, from
   parameter choice alone.
2. **`class_threshold` is "practically irrelevant"** — an independent
   confirmation of Zhang et al.'s own "limited impact on total errors".
3. **There is a bias-versus-noise trade-off along the cloth-resolution axis, and
   it runs the way that matters for us**: fine cloth (0.1 m) minimises the
   *combination* of bias and noise; **coarse cloth (2 m) gives "almost zero bias"
   at the cost of higher random error.** So a **finer cloth carries more positive
   (high) bias**, and a coarser cloth trades that bias for scatter.

Their recommendation, verbatim (and note how heavily they hedge it):
> "In case of CSF, which is based on different principles than previous
> algorithms, the most suitable parameters are similar for both types of point
> clouds (photogrammetric and LiDAR) and **we can recommend Cloth resolution
> value from 0.1–0.2 m. Note however, that these are recommendations for initial
> testing only and that fine tuning of parameters is always necessary.**"

**Their signed vertical errors, by algorithm, verbatim** (all against the 1414
GNSS checkpoints; positive = DTM above true ground):
> "In terms of RMSE, all algorithms yielded very good results, with RMSE ranging
> from 0.13 m (SMRF algorithm with LiDAR point cloud) to 0.23 m (ATIN algorithm
> with UAV data). … regardless of the point cloud, PTIN, PMF and SMRF performed
> 0.03–0.05 m better than ATIN, ArcGIS or CSF. **All algorithms overestimated the
> terrain** … With LiDAR data, the best performing algorithm was SMRF with 0.04 ±
> 0.13 m overestimation, followed by PMF and PTIN with overestimation of 0.08 ±
> 0.12 and 0.10 ± 0.12 m, respectively. **ATIN, ArcGIS and CSF all overestimated
> the terrain by approx. 0.12–0.13 (± 0.13–0.14) m.**"

**And the slope dependence, verbatim:**
> "The terrain accuracy (in terms of both random-error and bias) decreased with
> increasing slope. **The mean bias tended to increase with slope, especially in
> the low grass vegetation class.** The same pattern was observed for example by
> Hollaus et al. (2006). However, they found a mean bias of 0.15 to 0.25 m for
> slopes steeper than approximately 30°, while **we observed the same bias for
> slopes steeper than approximately 10°.** We assume that to be due to the
> presence of vegetation close to the ground, which is more problematic to filter
> out on the relatively steep slopes. **Such vegetation is not present under
> canopy vegetation (i.e. in forests)** …"

Note that last clause carefully: they attribute the slope-growing bias to **low
vegetation on open slopes**, and explicitly say it does not apply under closed
canopy. Our steep ground *is* under closed canopy, so this particular mechanism
may not transfer — but our uplands are row-crop stubble in November, which is
exactly the "low vegetation" case, on ground that is mostly flat.

Also relevant, on the algorithm groupings:
> "the ATIN, ArcGIS, and CSF algorithms were almost indistinguishable in their
> effect on ground filtering accuracy and performed generally poorer than the
> remaining three algorithms (PTIN, PMF, SMRF). Among remaining algorithms, SMRF
> performed in most cases best when using LiDAR point cloud."

Model quality: linear mixed model with checkpoint as random intercept,
conditional R² = 0.67, random-effect SD 0.106 m, residual SD 0.087 m, with
significant interactions among algorithm, slope and vegetation.

Their own caveat on the steep end:
> "The higher random-error (i.e. the width of confidence intervals) observed on
> steep slopes can be attributed to relatively low number of validation points
> especially in case of medium vegetation class at higher slopes"

⚠ **A naming discrepancy worth recording.** They describe CloudCompare's
rigidness presets as "*Ruggedness* (1 – 'Mountain', 2 – 'Complex', or 3 –
'Flat')". The current CloudCompare source labels the same three radio buttons
`Steep slope` / `Relief` / `Flat` (§0). The mapping of integer to terrain class is
the same; only the wording changed between versions.

### 2.2 Bailey et al. 2022 — a sensitivity screen that ranks the parameters differently

**Bailey, G.; Li, Y.; McKinney, N.; Yoder, D.; Wright, W.; Herrero, H. (2022).**
"Comparison of Ground Point Filtering Algorithms for High-Density Point Clouds
Collected by Terrestrial LiDAR." *Remote Sensing* 14(19):4776. DOI:
[10.3390/rs14194776](https://doi.org/10.3390/rs14194776). Peer-reviewed, open
access. *Full text read (delegated search); I did not re-verify.*

**Setting: terrestrial LiDAR on a hillslope erosion plot, scanned before and
after vegetation removal; centimetre-scale micro-topography.** Not transferable
in magnitude — only the ranking is of interest.

> "Of these parameters, only 'cloth resolution', 'time step', and 'classification
> threshold' were found to substantially affect the classification. The parameters
> that did not affect classification were set to the recommended default values of
> 'smoothing = false', **'rigidness = 3'**, and 'iterations = 500'. For the
> parameters affecting the classification, combinations of the following value
> ranges were used: 'cloth resolution' from 0.002 m to 0.012 m; 'time step' from
> 0.1 m to 0.65 m; and 'classification threshold' from 0.005 m to 0.01 m."

Optimum: cloth resolution 0.005 m, classification threshold 0.009 m.

**This partly contradicts Klápště**: here `class_threshold` *did* matter and
`rigidness` did not. The reconciliation is scale — on a single smooth hillslope
scanned at centimetre resolution there is no relief for rigidness to respond to,
and the threshold is doing the work that cloth resolution does at airborne
scales. Treat it as evidence that **which parameter dominates is
scale-dependent**, not as a contradiction of Klápště at ALS scale.

### 2.3 Zhang et al. 2016's own sweeps

Covered in §1.4. In summary: `time_step` swept 0.4-1.5 in 0.05 steps, minimum at
0.65 for all three scene groups; `cloth_resolution` swept over an **unstated
range**, best "around 0.5" for almost all samples, with group II (complex
terrain) least stable; `class_threshold` swept, "limited impact on total errors";
`hdiff` never swept.

### 2.4 What the sensitivity literature agrees on

- **`cloth_resolution` dominates** at airborne scales (Klápště, explicitly ranked;
  Zhang, implicitly, since it is the only continuous parameter he had to fix by
  sweep). Effect size: LiDAR RMSE 0.18-0.53 m across the tested grid.
- **`class_threshold` is close to irrelevant** — agreed by Klápště ("practically
  irrelevant"), Zhang ("limited impact"), and Štular ("on average or rugged
  terrain the default value is recommended"). Bailey dissents at centimetre scale.
- **`time_step` has one published optimum (0.65) and nobody has found reason to
  move it** at airborne scales.
- **`rigidness` is the least well-characterised of the four, and nobody has swept
  it.** Klápště folds it into parameter *combinations* rather than ranking it
  separately; Bailey found no effect at centimetre scale. Wan et al. 2018 (§5.2)
  is the closest thing to a rigidness study and it contains **no sweep of error
  across rigidness 1/2/3** — no table, no figure, no type I/II split. Its
  headline 0.36 pp mean cost of automatic selection is offset by "27% of the
  samples had larger biases (>1% in Total Error)". **No one has measured
  rigidness's effect on vertical bias at all.**
- **`hdiff` has never been swept by anyone.**

## 3. Applications in steep and/or forested terrain: what values were used, and why

Ordered by transferability to a low-density, leaf-off, steep, deciduous-forested
ALS tile.

### 3.1 Štular & Lozić 2020 — the only verified ALS study spanning our density range

**Štular, B.; Lozić, E. (2020).** "Comparison of Filters for Archaeology-Specific
Ground Extraction from Airborne LiDAR Point Clouds." *Remote Sensing*
12(18):3025. DOI: [10.3390/rs12183025](https://doi.org/10.3390/rs12183025).
Peer-reviewed, open access. *Full text read (I re-downloaded and re-checked this
one myself because its rigidness result is load-bearing).*

**Setting.** Four ALS sites, 500 × 500 m each, all characterised as "Vegetation
on steep slopes" and "Sharp discontinuities", all **deciduous forest**:

| Site | Density | Cover, verbatim |
|---|---|---|
| AT (Wildon, Austria, 2009) | 18.79 pts/m² | "covered with a dense, mature deciduous forest with little undergrowth" |
| SI1 (Pivka, Slovenia, 2014) | 12.28 pts/m² (functional 6.14) | "decimated deciduous forest in the first stage of regrowth" |
| SI2 (Pivka, 2014) | 12.43 pts/m² (functional 6.14) | "meadows, each surrounded by very dense deciduous hedges" |
| **ES (Santiago de Compostela, 2015)** | **1.83 pts/m²** | "covered with open deciduous forest and with two small clearings" |

Leaf state is **not stated anywhere in the paper** (I grepped the full text for
"leaf" — no hits).

**Best CSF settings actually used, verbatim from their Table A1** ("Values used
in each of the best results"), CSF row, which I read directly from the extracted
table:

| CSF setting | Tested range | Step | AT | SI1 | SI2 | **ES (1.83 pts/m²)** |
|---|---|---|---|---|---|---|
| s (Scene / rigidness) | hard – med | / | hard | hard | hard | **medium** |
| r (Cloth resolution) | 0.3 – 1 | 0.5 (0.2 for values under 0.5) | 0.3 | 0.5 | 0.5 | **1** |
| it (Max. iterations) | 1000 – 1000 | 250 | 1000 | 1000 | 1000 | **1000** |
| th (Classification threshold) | 0.5 – 0.5 | 0.5 (0.2 under 0.5) | 0.5 | 0.5 | 0.5 | **0.5** |
| sp (Slope processing) | on – on | / | on | on | on | **on** |

Their per-parameter reasoning (Appendix A.3.2; their italics are their
experimental notes). All quotes verbatim:

> "(s) Scene: the rigidness of the cloth is determined; very soft to fit steep
> slope, medium for rugged terrain, and hard cloth for flat terrain (default:
> medium). *No discernible differences between very soft and medium were
> observed, hard was not suitable.*"

> "(r) Cloth resolution: refers to the grid size of the cloth with which the
> terrain is covered, i.e., the distance between the particles in the cloth
> (default: 2.0). *The same value as the resolution of the final DEM worked
> best, lower values (e.g., ½ the cell-size of the final DEM) introduced
> artifacts. Lowering the value increases the processing times exponentially.*"

> "(it) Max. Iterations: maximum iterations for the cloth simulation (default:
> 500 iterations). *A significant increase of the iterations (max. value 1000)
> produced modest improvements at the expense of a linear increase of the
> processing time.*"

> "(th) Classification threshold: the distance to the simulated cloth to classify
> a point cloud into ground and non-ground (default: 0.5). *On flat terrain lower
> values (e.g., 0.15) slightly reduced T2 error, but on average or rugged terrain
> the default value is recommended.*"

> "(sp) Slope processing: reduces errors during post-processing when steep slopes
> are present (default: off). *Switching on significantly improved handling of
> cliffs and low walls, but resulted in noticeable T2 errors, e.g., "tree stumps"
> on slopes and buildings on the flats.*"

**⚠ An error in this paper, which I verified directly.** The Appendix text says
"*hard was not suitable*", but Table A1 reports **hard** as the best setting at
three of the four sites, and the *tested range* for `s` is only "hard – med" —
**"very soft" was never actually run**. The appendix sentence therefore cannot
be a report of an experiment. The paper does not reconcile this. **Use Table A1
as the record of what was run; treat "no discernible differences between very
soft and medium" as unsupported by their own tested range.**

**Two further verbatim findings from this paper:**

> "Concerning the influence of point density on filter performance, our results
> are consistent with the results presented by Sithole and Vosselman [20], namely
> in that the lower the point density, the worse the performance of all filters.
> This makes processing low-density data much more demanding."

> "It is noticeable that in our assessment, filters based on newer algorithms
> (PMF, CSF, SMRF) are outperformed by some of the oldest (PTIN, SBF), which
> contradicts recently published findings [46]."

Also, on steep-slope vegetation specifically:

> "Also, the removal of vegetation on steep slopes remains problematic, but with
> dense data (AT) it is solved well by all filters, as is the discontinuity
> retention."

— i.e. the steep-slope-under-canopy problem is *density-limited*, and at our
density it does not go away.

**Note the disagreement between their two stated criteria for cloth
resolution.** Their own justification is "the same value as the resolution of
the final DEM", not point spacing. Their four chosen values also happen to track
mean point spacing (1/√density: 0.23, 0.40, 0.40, 0.74 m against r = 0.3, 0.5,
0.5, 1.0) — but since they also chose their DEM resolution from density, the two
criteria are confounded in this dataset and **cannot be separated here**. Their
DEM resolutions were 0.25-1 m; extrapolating "cloth = DEM cell size" to our 5 m
grid is far outside their tested range and is not supported.

### 3.2 Storch et al. 2022 — the cleanest flat-vs-steep decision in deciduous forest

**Storch, M.; Jarmer, T.; Adam, M.; de Lange, N. (2022).** "Systematic Approach
for Remote Sensing of Historical Conflict Landscapes with UAV-Based
Laserscanning." *Sensors* 22(1):217. DOI:
[10.3390/s22010217](https://doi.org/10.3390/s22010217). Peer-reviewed, open
access. *Full text read by delegated search; I did not re-verify.*

**Setting.** Hürtgenwald, Germany. Two contrasting sites: a **flat** war
cemetery, and the **Kall valley** — "the valley's characteristics are its steep
gorges and slopes", "on the eastern side of the valley, there is a relatively
intact deciduous forest with hardly any vegetation close to the ground".
**UAV-LiDAR, mean point density 69-412 pts/m²** (ground 33-206 pts/m²) — i.e.
~100× our density, which is the main transfer limitation. **Both leaf states
flown**: "The UAV-LiDAR data acquisition took place on 19/20 August 2020 under
leaf-on conditions and on 6 March/1 April 2021 under leaf-off conditions."

Flat site, verbatim:
> "The filter is parameterized according to the recommendations of, e.g., Zhang
> et al. [53] and Serifoglu Yilmaz et al. [54]. The maximum number of iterations
> is set to 500, time step size is set to 0.65 and rigidness is set to 3 for flat
> terrain. The cloth resolution is set to 0.2 m to account for the high
> resolution input data, and the classification threshold is set to 0.1 m to
> ensure that even small memorial stones with a height of a few decimeters
> located on the surface are not classified as ground."

Steep forested site, verbatim:
> "CSF is parameterized in a similar way to the other use cases, with the
> exception for the rigidness parameter which is set to 1 in order to take into
> account the steep terrain relief."

This is the clearest published worked example of the decision we face: same
crew, same sensor, same processing — **flat → rigidness 3, steep → rigidness 1**,
with the reason stated. It does **not**, however, address a scene containing
both. And note they used the *same* CSF parameters for leaf-on and leaf-off.

### 3.3 Other steep / forested applications (values reported, mostly unjustified)

*All from delegated search; DOIs confirmed against OpenAlex by the searcher. I
did not re-read these full texts.*

- **Micu, M. et al. (2023).** "Deciphering Complex Morphology and Structural
  Connectivity of High-Magnitude Deep-Seated Landslides via Airborne Laser
  Scanning: A Case Study in the Vrancea Seismic Region, Romanian Carpathians."
  *Remote Sensing* 15(22):5286. DOI:
  [10.3390/rs15225286](https://doi.org/10.3390/rs15225286). Peer-reviewed.
  *Full text read (delegated).* Romanian Carpathians, deep-seated landslide,
  "(coniferous and broad leaved) forest coverage", UAV-LS at "point cloud
  resolution 85–110 points/m²". Their Table 3, "CSF plugin parameters in
  CloudCompare": Scenes = **"Steep slope"**, Slope Processing = **"On"**, Cloth
  Resolution = **0.6**, Max Iterations = **600**, Classification threshold =
  **0.5**. **No justification given.**

- **Marotta, F. et al. (2021).** "Integrated Laser Scanner Techniques to Produce
  High-Resolution DTM of Vegetated Territory." *Remote Sensing* 13(13):2504.
  DOI: [10.3390/rs13132504](https://doi.org/10.3390/rs13132504). Peer-reviewed.
  *Full text read (delegated).* Italian Alps, Chiavenna valley: "both of the
  sides, the northern one with an average slope of 28° and the southern one with
  34° average slope, are covered with very dense vegetation" (deciduous:
  "chestnut, hornbeam, ash, sycamore, linden, oak, cherry, and mountain ash").
  **Long-range TLS, not ALS.** Verbatim: "After several tests, the final ground
  points were obtained selecting a 0.50 m cloth resolution and a 0.10 m
  classification threshold." Justification: "after several tests" only. Their
  transferable contribution is a *pre-processing* step: "the single + first
  echoes of each of the identified areas, each one tilted by ±30°, ±45° and ±60°
  to preserve points also in the steepest slopes. To tilt the PCs, an
  interpolating plane was first fitted to the involved area and then made
  horizontal."

- **Štroner, M. et al. (2021).** "Vegetation Filtering of a Steep Rugged Terrain:
  The Performance of Standard Algorithms and a Newly Proposed Workflow on an
  Example of a Railway Ledge." *Remote Sensing* 13(15):3050. DOI:
  [10.3390/rs13153050](https://doi.org/10.3390/rs13153050). Peer-reviewed
  (Technical Note). *Full text read (delegated).* **UAV SfM photogrammetry, not
  lidar**; vegetated rock face, Czechia, cloud subsampled to 5 cm. Exhaustive
  sweep: cloth resolution "0.1; 0.2; 0.3; 0.4; 0.5; 0.75; 1.0"; classification
  threshold "0.1; 0.2; 0.3; 0.4; 0.5; 0.75; 1.0; 2.0; 2.5"; Scene "Steep slope";
  Slope processing "Yes"; Max iterations "500". Optimum in the original tilted
  orientation: resolution **0.1**, threshold **2.0**, total error **16.0%**.
  After levelling the cloud to horizontal: resolution **0.1**, threshold
  **0.1**, total error **6.2%**.

**The one cross-cutting methodological finding.** Two independent groups
(Štroner et al. on a very steep rock face; Marotta et al. on 28-34° alpine
forest) found that **rotating the cloud so the local ground plane is horizontal
before running CSF** materially improves it — Štroner measured 16.0% → 6.2%
total error and a 20× shift in the optimal threshold. The mechanism is
structural: CSF's cloth falls under gravity *in the cloud's own frame*, so on a
uniformly steep surface the method is fighting its own assumption. Neither study
is on terrain as gentle as a dissected-plateau hillslope, and neither is ALS
under canopy, so this should be **tested, not assumed**, for us.

### 3.4 Sparse / early-generation ALS: essentially one data point

**Szabó, Z.; Tóth, C.A.; Holb, I.; Szabó, S. (2020).** "Aerial Laser Scanning
Data as a Source of Terrain Modeling in a Fluvial Environment: Biasing Factors
of Terrain Height Accuracy." *Sensors* 20(7):2063. DOI:
[10.3390/s20072063](https://doi.org/10.3390/s20072063). Peer-reviewed, open
access. *Full text read (delegated).* **ALS at 4 pts/m²**, Tisza floodplain,
NE Hungary — **flat**, floodplain mosaic (forest, grazing land); leaf state not
stated.

> "CSF parameters are dependent on the point cloud density (cloth size, i.e.,
> resolution) and the complexity of the terrain (threshold, i.e., the distance
> between points and the simulated terrain), and although suggestions exist,
> **there is no definite rule for setting these parameters. Rather, a range can
> be used to find the ideal ones.** We conducted the classification with the flat
> terrain option (as the relief was very low in the study area) and with cloth
> sizes of 2 and 5 (the larger the cloth, the coarser the DTM), according to
> [29]. A value of 2 was suggested based on our point density, but we also
> intended to test the effect of a larger value (i.e., 5) We used classification
> thresholds of 0.2, 0.5, and 1, with 500 iterations in each case."

> "the best CSF parameters were in accordance with the previous results, and the
> fewest points, the cloth size of 5, and the threshold of 0.2 resulted in the
> highest accuracy."

**Read this carefully — it cuts against the point-spacing rule.** At 4 pts/m²
(spacing 0.5 m) the *best* cloth was **5 m**, ten times coarser than the spacing
rule would give. Their terrain is flat, and a coarse cloth on flat terrain
suppresses type II error without costing anything, because there are no breaks
to follow. Combined with Štular's ES site (1.83 pts/m², **steep**, best r = 1.0)
and Cai et al. 2019 (§4.3, mixed, finer cloth = lower type I / higher type II),
a consistent picture emerges: **the right cloth resolution is set by terrain
roughness at least as much as by point density, and the two studies at ALS
densities pull in opposite directions because their terrains differ.**

### 3.5 Verified silence: driftless / loess / dissected plateau

The delegated search ran a full-text query for CSF combined with
loess / gully / badland / ravine / dissected terms across the works citing Zhang
et al. 2016 and returned **zero** hits. The one forested-gully study located
(Manić, M. et al. 2022, *Frontiers in Environmental Science* 10:897248, DOI
[10.3389/fenvs.2022.897248](https://doi.org/10.3389/fenvs.2022.897248),
peer-reviewed, *full text read (delegated)*) reports **no parameter values at
all**, saying only:

> "The set parameters depend primarily on the terrain configuration."

**There is no published CSF parameterisation for driftless-area or comparable
dissected-plateau terrain.** This is a genuine gap, not a search failure.
## 4. Sparse / early-generation point clouds, and cloth resolution vs. point spacing

*(Sub-question 4, taken here out of order because it depends only on sources I
verified directly.)*

### 4.1 The "cloth resolution ≈ point spacing" rule is a wrapper-doc claim, not a result

The rule appears **only** in the RCSF/lidR roxygen documentation
(*software docs*, `Jean-Romain/RCSF` `R/CSF.R`, and `r-lidar/lidR`
`R/algorithm-gnd.R`, *full text read*):

> "`cloth_resolution` scalar. The distance between particles in the cloth. This
> is usually set to the average distance of the points in the point cloud. The
> default value is 0.5."

It is **not** in Zhang et al. 2016 (I searched the full text for "spacing",
"point distance", "density" — the only relevant statements are §1.4 above and
the §3.2 remark about 0.6-0.8 m point distance in the dense test). It is **not**
in the authors' own CloudCompare wiki text (§1.7). No supporting study is cited
for it in either wrapper. **Treat it as folklore.**

### 4.2 What density Zhang et al. actually tuned GR = 0.5 on — and it is our density

The ISPRS WG III/3 benchmark that produced GR = 0.5 is *sparse*, roughly our
regime. Two independent statements:

- **Sithole, G.; Vosselman, G. (2003).** "Comparison of Filtering Algorithms."
  *ISPRS Archives* XXXIV-3/W13, Dresden. Stable URL:
  `https://www.isprs.org/proceedings/xxxiv/3-w13/papers/Sithole_ALSDD2003.pdf`
  (*conference proceedings, not journal peer-review; full text read*):
  > "The landscape was scanned with an Optech ALTM scanner, and the data was
  > produced by FOTONOR AS. Both first and last pulse data were recorded. Eight
  > test sites (four urban and four rural) were chosen. **The urban sites were at
  > a resolution of 1-1.5m. The rural sites were at a resolution of 2-3.5m.**"

  (The journal version of this experiment is Sithole, G.; Vosselman, G. (2004),
  "Experimental comparison of filter algorithms for bare-Earth extraction from
  airborne laser scanning point clouds," *ISPRS J. Photogramm. Remote Sens.*
  59(1-2):85-101, DOI
  [10.1016/j.isprsjprs.2004.05.004](https://doi.org/10.1016/j.isprsjprs.2004.05.004)
  — *paywalled, abstract only*. Its abstract states "The influence of point
  density could not well be determined in this experiment," and §5.2 of the 2003
  proceedings version says the same at greater length: "More tests on decreasing
  resolution will need to be done, as the test sites chosen have proved
  inadequate to obtain a conclusive picture of the effects of resolution on
  filtering.")

- **The ISPRS WG III/3 filter-test site pages themselves** (ITC / University of
  Twente, `https://www.utwente.nl/en/itc/isprs/wgIII-3/filtertest/downloadsites/`)
  — *primary benchmark documentation, not peer-reviewed; full text read*:
  > "Eight sites have been chosen for comparing the performance of filters. Four
  > represent urban landscapes (City Sites) and the other four represent rural
  > landscapes (Forest Sites). **The point density for the City and Forest sites
  > are roughly 0.67 and 0.18 points per square metre respectively.**"

  The per-site table gives "0.67 points per square metre (point spacing: 1.0 -
  1.5m)" for City and "0.18 points per square metre (point spacing: 2.0 - 3.5m)"
  for Forest, with decimated variants down to 0.01 pts/m². **Note the benchmark's
  own naming: its four rural sites are "Forest Sites".**

- **Cai, S.; Zhang, W.; Liang, X.; Wan, P.; Qi, J.; Yu, S.; Yan, G.; Shao, J.
  (2019).** "Filtering Airborne LiDAR Data Through Complementary Cloth
  Simulation and Progressive TIN Densification Filters." *Remote Sensing*
  11(9):1037. DOI: [10.3390/rs11091037](https://doi.org/10.3390/rs11091037).
  Peer-reviewed, open access. *Full text read.* Same lab as Zhang et al. 2016
  (Wuming Zhang is corresponding author):
  > "The point density is 0.4–1 points/m² for urban sites and 0.08–0.25
  > points/m² for rural sites, respectively."
  >
  > "Since the ISPRS data were collected almost two decades ago, the average
  > point density was relatively low ranging from 0.08–1 points/m²."

**Our ground-return density of ~0.5 pts/m² sits between the benchmark's Forest
sites (0.18) and its City sites (0.67).**
Zhang et al.'s GR = 0.5 m was optimised at our sparsity, on data from the same
sensor family (Optech ALTM). This is unusually good transferability, and it is
strong evidence **against** the folklore rule: at 0.08-1 pts/m² the mean point
spacing is ~1-3.5 m, so the optimal cloth was **2-7× finer than the point
spacing**, not equal to it.

**[Inference, not from the literature]** This is mechanically consistent with the
algorithm as described in Zhang et al. step 4: each cloth particle takes the
*nearest* lidar point in the horizontal plane as its "corresponding point" and
records that point's height as the floor it can fall to. A cloth finer than the
point spacing therefore does not invent detail — several particles simply share
one corresponding point — while a cloth coarser than the point spacing discards
low points that no particle happens to land nearest to. The asymmetry favours
erring fine.

### 4.3 The only direct cloth-resolution comparison at this density

Cai et al. 2019 (above) tested cloth resolution 0.5 m vs 1 m on the ISPRS
samples, i.e. at 0.08-1 pts/m². *Full text read.*

> "In CSF [44], CR is set to 0.5 m based on the principle of optimal accuracy.
> To objectively compare the proposed algorithm and CSF, the proposed algorithm
> under CR = 0.5 m was first tested. Then, the proposed algorithm under CR = 1 m
> was tested…"

> "Exploration of further details for some recognizable differences between them
> shows DTM_CR=0.5 is closer to the reference DTM than DTM_CR=1 in the terrain
> discontinuities … However, compared to DTM_CR=1, there are more non-ground
> objects in DTM_CR=0.5 … Compared to SDE_CR=1, SDE_CR=0.5 shows lower type I
> error but higher type II error. The reason can be explained by exploring the
> procedure of CS. When cloth gradually approaches terrain with more cloth
> particles, a more detailed initial terrain is constructed, thereby assisting
> PTD to correctly classify the ground point in complex terrain. However, the
> ground seed points are not 100% correctly obtained by CS. Thus more seed
> points will increase the risk of type II error."

**Caveat on Zhang's GR sweep, verified.** The paper states the time-step sweep
range numerically ("from 0.4 to 1.5 with steps of 0.05") but **never states the
range or step of the grid-resolution sweep** — only "Figure 14 shows the total
errors at different GR values … almost all samples get the highest accuracy
around 0.5". I searched the full text; the numbers exist only in the figure. So
"0.5 is optimal" is a reading of an unlabelled sweep, and how far either side of
0.5 was tested is unknown.

**The trade-off is the same shape as the rigidness trade-off:** finer cloth =
fewer true-ground losses at terrain discontinuities (lower type I), more
low vegetation retained (higher type II).

*Caveat on transfer:* this comparison is embedded in their hybrid CSF+PTD
algorithm, not in bare CSF. The mechanism they invoke ("a more detailed initial
terrain … at the terrain discontinuities") is a property of the cloth stage, but
the reported type I/type II split is for the hybrid.

### 4.4 A note on PDAL's `returns` default

PDAL's `filters.csf` defaults to `returns: "last, only"` — it feeds CSF only
last and single returns. That is the correct choice for bare-earth on
leaf-off data and is worth stating explicitly, because it is a PDAL-specific
behaviour with no counterpart in the paper, RCSF, or CloudCompare, all of which
filter whatever cloud you hand them. **This is an implementation difference
that changes the effective point density seen by the cloth.**

---

## 5. Heterogeneous terrain in one tile: a single rigidness cannot be right for both

This is the crux of our problem, and the literature **recognises it explicitly**
but the published fixes are all *new algorithms*, not parameter advice for stock
CSF.

### 5.1 The problem is stated in the literature

**Cai, S.; Yu, S.; Hui, Z.; Tang, Z. (2023).** "ICSF: An Improved Cloth
Simulation Filtering Algorithm for Airborne LiDAR Data Based on Morphological
Operations." *Forests* 14(8):1520. DOI:
[10.3390/f14081520](https://doi.org/10.3390/f14081520). Peer-reviewed, open
access. *Full text read.*

> "However, CSF also has limitations. First, it is difficult for cloth to cover
> steep slopes, where accurate reference terrain cannot be obtained to
> distinguish ground points from non-ground points [26]. … **Second, a single
> cloth rigidness is unreasonable, because multiple terrain features (i.e., flat
> terrain, slopes and raised terrain) are usually contained in a landscape**
> (Figure 1b) [50]. Third, the points on rugged terrain are often misclassified
> when using a fixed height difference threshold, since terrain slopes are not
> considered [26,50]."

And its summary of the two prior adaptive approaches:

> "To solve these problems, Yang et al. [50] partitioned a point cloud into
> multiple regions, where the richness of terrain features became lower relative
> to in the entire scene, and thus the negative impact of cloth rigidness was
> reduced. In addition, the reference terrain of steep slopes was constructed
> using a bidirectional cloth simulation method. Finally, ground points were
> extracted based on adaptive height difference thresholds, which were
> calculated using a weighted sum of the height differences of unclassified
> points and their neighboring points. Wan et al. [51] proposed a terrain relief
> index to automatically estimate the cloth rigidness applicable to a scene.
> **These improved algorithms improve the filtering accuracy of CSF, but
> sacrifice its ease of use, due to the introduction of many additional
> parameters.**"

References [50] and [51] are, verbatim from the ICSF reference list:

- Yang, A.; Wu, Z.; Yang, F.; Su, D.; Ma, Y.; Zhao, D.; Qi, C. "Filtering of
  airborne LiDAR bathymetry based on bidirectional cloth simulation." *ISPRS J.
  Photogramm. Remote Sens.* **2020**, 163, 49–61. (*Bathymetry — low
  transferability to our case. NOT VERIFIED beyond this citation.*)
- Wan, P.; Zhang, W.; Skidmore, A.K.; Qi, J.; Jin, X.; Yan, G.; Wang, T. "A
  simple terrain relief index for tuning slope-related parameters of LiDAR
  ground filtering algorithms." *ISPRS J. Photogramm. Remote Sens.* **2018**,
  143, 181–190. — see §5.2.

### 5.2 Automatic scene-level rigidness selection, and how much it costs

**Wan, P.; Zhang, W.; Skidmore, A.K.; Qi, J.; Jin, X.; Yan, G.; Wang, T.
(2018).** "A simple terrain relief index for tuning slope-related parameters of
LiDAR ground filtering algorithms." *ISPRS J. Photogramm. Remote Sens.*
143:181–190. DOI:
[10.1016/j.isprsjprs.2018.03.020](https://doi.org/10.1016/j.isprsjprs.2018.03.020).
Peer-reviewed. *Full text read* (author copy at
`https://ris.utwente.nl/ws/files/72986973/Wan2018simple.pdf`; the
research.utwente.nl mirror and ScienceDirect are Cloudflare/paywall blocked).
Two authors (Wan, Zhang) are CSF co-authors, so this is quasi-authoritative.

Their Table 1, verbatim — the continuous quantity behind the integer rigidness:

| Parameter | Value | | |
|---|---|---|---|
| Terrain relief amplitude | sharply undulating | gently undulating | flat |
| Cloth rigidness | 1 | 2 | 3 |
| Spring-back ratio | 30% | 51% | 65.70% |

> "The cloth rigidness is a key parameter for CSF, which determines the
> spring-back ratio of the cloth in each iteration. **The cloth rigidness is set
> by users based on the terrain relief amplitude.** The spring-back ratio has a
> one-to-one correspondence to the cloth rigidness. This correspondence was
> predetermined by the developers, and can be found in Table 1."

> "With the increase of terrain relief amplitude, the terrain relief index will
> be increased, and the spring-back value should be accordingly decreased to
> make the cloth softer."

Result (their §4.4.1):

> "The CSF with optimal tuned parameters achieved the lowest mean Total Error
> (4.394%) and highest mean Kappa coefficient (83.86%) with all samples,
> followed by the CSF with TRIp(10), for which the Total Error was 4.758%, and
> the Kappa coefficient was 82.73%. The CSF with TRIc had a higher Total Error
> (4.762%) and lower Kappa coefficient (82.23%) on average."

**Three things to take from this — and the first one needs care.**

(a) The *mean* cost of getting rigidness "auto" rather than hand-tuned is **0.36
percentage points of total error** (4.758 vs 4.394). **Do not read that mean as
"rigidness is low-leverage."** The same paper gives the spread, verbatim:
> "In the test of the TRI-derived cloth rigidness for CSF (see Table 4),
> one-third of the samples achieved better filtering accuracies than the CSF with
> optimized parameters, and **27% of the samples had larger biases (>1% in Total
> Error).**"

So roughly a quarter of scenes moved by more than a full percentage point of
total error on a rigidness change. The small mean hides a real per-scene tail.

(b) Their index is also **less reliable in forest**, verbatim:
> "the discrepancy tends to be large in forest samples. The average absolute
> difference between TRIc and TRIp(10) in forest samples was 0.77%, while the
> value in city samples was 0.09%."

(That 0.77% vs 0.09% is the disagreement between their two *versions of the
index*, not the rigidness error itself — but it says the index is 8× noisier in
the land cover we have.)

(c) Their index is **scene-level**: it still assigns *one* rigidness to a whole
tile. It does not solve within-tile heterogeneity; it automates the same
one-value choice.

Also note their abstract statement (*full text read*):
> "the correlations were higher in city areas (r = 0.911 to 0.993) than in
> forest areas (r = 0.848 to 0.905)"
— the index is somewhat weaker in forest, which is our case.

### 5.3 The per-scene-region and per-particle approaches

**ICSF** (Cai et al. 2023, cited above) estimates rigidness *within* the scene
rather than taking one value, and uses terrain-adaptive (slope-aware) height
thresholds instead of a fixed `class_threshold`. It is the closest published
answer to "one rigidness cannot be right for both."

Its forested-sample test is the most transferable evaluation I found. Their
Table 2, verbatim (Open Topography, mountainous forested):

| Sample | Feature | Vegetation Cover (%) | Ref. ground pts | Ref. non-ground pts |
|---|---|---|---|---|
| S1 | Gentle slopes and dense vegetation | 80.08 | 473,538 | 1,493,507 |
| S2 | Undulating terrain and dense vegetation | 84.67 | 261,242 | 518,792 |
| S3 | Steep slopes and dense vegetation | 74.45 | 514,925 | 819,264 |
| S4 | Steep slopes and discontinuities | 8.74 | 298,955 | 44,499 |

Results (their Table 6 and §4.2):

| Sample | ICSF T.I(%) | ICSF T.II(%) | ICSF T.E.(%) |
|---|---|---|---|
| S1 | 0.12 | 2.75 | 2.12 |
| S2 | 0.41 | 2.94 | 2.09 |
| S3 | 5.14 | 8.44 | 7.17 |
| S4 | 7.30 | 5.95 | 7.12 |
| Average | 3.24 | 5.02 | 4.62 |

> "All algorithms worked well in the samples with gentle slopes (i.e., S1 and
> S2), with a total error of less than 5%. ICSF obtained significantly lower
> total errors in the samples with steep slopes and discontinuities (i.e., S3
> and S4). Overall, ICSF outperformed other filtering algorithms. More
> specifically, compared with MSF, PMF and CSF, ICSF reduced the average total
> error by 55.11%, 30.41% and 34.47%, respectively."

**The pattern to take away**: error roughly **triples** going from gentle-slope
dense forest (~2.1%) to steep-slope dense forest (~7.2%), even for a
slope-adaptive filter. Slope, not canopy, is the dominant error driver here —
S1/S2 have *higher* vegetation cover (80-85%) than S3 (74%) and far higher than
S4 (8.7%), yet much lower error.

**Important caveat, verified negative**: ICSF does **not** state the CSF
parameter values used for its baseline. I searched the full text — "rigidness"
occurs 14 times, all in the introduction/method exposition, never in an
experimental-settings statement, and "cloth resolution"/"grid size" never appear
as a setting. Its baseline CSF average total error on the ISPRS samples is
6.62% (their Table 4), notably worse than Zhang et al.'s own 4.39% on the same
samples — consistent with a single un-tuned parameter set being used for the
baseline. **Do not read ICSF's CSF-vs-ICSF margin as the margin against
well-tuned CSF.**

### 5.4 What is NOT in the literature

I found **no** published guidance on how to choose a single `rigidness` for a
scene that is genuinely bimodal in slope (broad flat uplands *plus* steep
valley walls) using stock CSF. Every treatment of the problem I located
(§5.1-5.3) responds by *modifying the algorithm*. There is also no published
recommendation to *tile a scene by slope and run CSF twice with different
rigidness* — the closest is Yang et al. 2020's region partitioning, which is
inside a new bathymetric algorithm and which I did not verify.

---

## 6. CSF against other ground filters, and what is known about vertical bias

### 6.0 The one study that measures a CSF vertical bias against absolute truth

**Klápště et al. 2020** is the only source I found that reports a *signed*
vertical error for CSF alongside competing filters, against surveyed ground
control, with each algorithm's parameters tuned separately. It is written up in
full in **§2.1** because it is also the only systematic CSF parameter sweep. The
headline for this section:

> "All algorithms overestimated the terrain … With LiDAR data, the best
> performing algorithm was SMRF with 0.04 ± 0.13 m overestimation, followed by
> PMF and PTIN with overestimation of 0.08 ± 0.12 and 0.10 ± 0.12 m,
> respectively. **ATIN, ArcGIS and CSF all overestimated the terrain by approx.
> 0.12–0.13 (± 0.13–0.14) m.**"

> "the ATIN, ArcGIS, and CSF algorithms were almost indistinguishable in their
> effect on ground filtering accuracy and **performed generally poorer than the
> remaining three algorithms (PTIN, PMF, SMRF)**. Among remaining algorithms,
> SMRF performed in most cases best when using LiDAR point cloud."

> "**The mean bias tended to increase with slope**, especially in the low grass
> vegetation class."

**CSF reads the ground HIGH, by roughly 3× SMRF's bias, and the bias grows with
slope.** Setting caveats — 8 pts/m², a vegetated spoil heap, "undulated" rather
than dissected relief, leaf-on ALS — are in §2.1.

**The single most directly useful number for us**: on the same cloud,
**CSF − SMRF ≈ +0.08-0.09 m** and **CSF − PTIN (`lasground`) ≈ +0.02-0.03 m**.
Since gen2 here is vendor class-2 from the Axelsson/TerraScan lineage, that
second figure is the published order of magnitude for the *filter-family
mismatch* baked into our DoD — and its sign says gen1 (CSF) should read
**higher**, i.e. it pushes gen2 − gen1 *negative*.

### 6.1 The best-matched study by design — and it is already in this repo

**Viedma, O. (2022).** "Applying a Robust Empirical Method for Comparing Repeated
LiDAR Data with Different Point Density." *Forests* 13(3):380. DOI:
[10.3390/f13030380](https://doi.org/10.3390/f13030380). Peer-reviewed, open
access. *Full text read.*

**This paper is already cited in `DTM_PERCENTILE_LITERATURE.md` §5.1**, where it
is flagged as reporting "CSF as the worst-performing filter, which deserves
attention since this project uses CSF for the 2008 epoch". I re-read it here for
the parameter question, and it turns out to be the closest published analogue of
our whole experimental design.

**Setting.** Six sites, Sierra de Gredos, Central Spain, mean 2.77 ha each.
**Low-density lidar = Spanish PNOA at "0.5 points/m²"** (stated twice in the
text; note their Table 2 says "1.5" for the same dataset — an internal
inconsistency in the paper), differenced against a **high-density benchmark at
"around 300 points/m²"**. Mean slope by site: 25.2°, 8.6°, 14.2°, 24.7°, 16.6°,
15.9° — our slope range. Vegetation: Mediterranean *Quercus ilex / suber /
pyrenaica* and *Pinus pinaster* with heather/rockrose/broom shrubland; five of
six sites burned between 1985 and 2019, so vegetation height is <4 m at four
sites and >15 m at two. **That is the main transfer limitation** — it is not
temperate deciduous forest, and leaf state is not stated.

**Sign convention, verbatim:**
> "Positive values in elevation differences indicated that the surface derived
> from the low-density LiDAR was higher than the benchmark elevation and then,
> over-predicted it, whereas negative values indicated the opposite."

**Filters compared:** progressive TIN densification via LAStools `lasground`
(five parameterisations: DEF, SF, SW2, WILD), PMF via lidR, and CSF via lidR.

**Results, uncorrected (their Table 5), by algorithm.** Letters are Dunn post-hoc
groups at p < 0.05; same letter = not significantly different.

| Algorithm | P50 (m) | NMAD (m) | RMSE (m) |
|---|---|---|---|
| **CSF** | **−1.19 ± 1.10** a | **0.90 ± 0.19** a | **1.84 ± 0.76** a |
| DEF (lasground default TIN) | −0.57 ± 0.84 a | 0.49 ± 0.16 b | 1.06 ± 0.53 b |
| PMF | −0.74 ± 0.88 a | 0.55 ± 0.20 b | 1.26 ± 0.69 b |
| SF | −0.65 ± 0.85 a | 0.53 ± 0.17 b | 1.14 ± 0.55 b |
| SW2 | −0.66 ± 0.85 a | 0.51 ± 0.19 b | 1.15 ± 0.65 b |
| WILD | −0.70 ± 0.87 a | 0.61 ± 0.17 b | 1.31 ± 0.50 ab |

**Read the letters carefully — this is where it is easy to overclaim.** CSF's
median (P50) is the most negative by ~0.5 m, but **all six algorithms share
letter "a" on P50, so the median differences between filters are NOT
statistically significant** in this design. What *is* significant is **scatter**:
CSF is in its own group on both NMAD (0.90 vs 0.49-0.61) and RMSE (1.84 vs
1.06-1.31).

**And the CSF penalty survives their correction surface (their Table 6):**

| Algorithm | P50 (m) | NMAD (m) | RMSE (m) |
|---|---|---|---|
| **CSF** | −0.019 ± 0.034 a | **0.15 ± 0.24** a | **0.62 ± 0.56** a |
| DEF | −0.011 ± 0.005 a | 0.07 ± 0.02 ab | 0.25 ± 0.09 b |
| PMF | −0.004 ± 0.006 a | 0.07 ± 0.02 ab | 0.30 ± 0.14 bc |
| SF | −0.008 ± 0.005 a | 0.06 ± 0.02 b | 0.26 ± 0.10 b |
| SW2 | −0.007 ± 0.004 a | 0.07 ± 0.02 ab | 0.26 ± 0.09 b |
| WILD | 0.001 ± 0.012 a | 0.07 ± 0.02 ab | 0.38 ± 0.16 ac |

After correction, the systematic offset is gone for every filter (P50 ≈ 0), but
**CSF still carries roughly twice the residual NMAD and RMSE of the TIN
variants** — i.e. its penalty is in the part a correction surface cannot remove.

Their qualitative diagnosis, verbatim:
> "the TIN densification algorithms using default parameters (DEF) or low spikes
> (0.5) (SF, SW2) were smoother (Figures 3 and S3) showing lower vertical errors
> (Table 5) whereas the CSF and the TIN densification algorithms using large
> spike (3) (i.e., WILD), followed by the morphological filter (PMF), showed a
> lot of bumps and spikes (Figures 3 and S3); and also, high elevation
> differences (Table 5)."

And on slope, verbatim (**note: this is the global GAM across all filters, not a
CSF-specific result**):
> "The slope and the distance to the nearest geoid point were the most important
> explanatory variables. Overall, as the slope increased, the elevation residuals
> were more negative (underestimation of benchmark elevations) whereas as the
> distance to the nearest geoid point increased, the elevation residuals were more
> positive (overestimation of benchmark elevations)."

> "all these GAMs showed low fitting ability (from 13% to 57% of deviance
> explained)"

**⚠ Verified negative, and it is a serious limitation.** Viedma **does not state
the CSF parameters used**. For PMF the paper gives them explicitly ("we used the
settings proposed by Zhang et al. [39]: (a) ws … seq (3, 12, 4); and (b) th …
seq (0.1, 1.5, length.out = length (ws))") and for `lasground` it gives step
size, spike thresholds, offset and bulge — but for CSF it says only that it "was
implemented in the lidR package". Presumably that means lidR's defaults
(`rigidness = 1`, `cloth_resolution = 0.5`, `sloop_smooth = FALSE`), which would
make the worst-performing filter in this comparison a *near-neighbour of our
current configuration* with slope post-processing off — but **that is an
inference, not something the paper states**. The comparison is therefore between
a carefully parameterised `lasground` and an un-documented CSF, which limits how
much of the gap can be attributed to the algorithm rather than to its settings.

### 6.2 CSF ranked against eight other filters in steep deciduous forest

Štular & Lozić 2020 (full citation in §3.1; *full text read*) rank nine filters
on the same steep deciduous-forested ALS data at 1.83-18.79 pts/m². Verbatim:

> "The best score in total error was obtained by SMRF, followed by a close group
> of BMHF, PTIN, SBF, and CSF. At the bottom of the group were PMF, SegBF, and,
> at some distance, WLS. The most robust performers were PTIN and BMHF (twice and
> once second best, respectively; never among the two worst), followed by SBF and
> CSF (never among the worst). **The T1 error results are almost identical to
> total error.**"

> "**The T2 error is very different, though. MCC and WLS are decidedly the best,
> followed by a narrow group of PTIN, BMHF, and PMF. SBF, CSF, SegBF, and at some
> distance SMRF are the worst.**"

And on low-density data specifically (qualitative assessment):
> "PTIN and BMHF performed best on medium density data, followed by CSF and SMRF.
> For low density data, PMF, PTIN, and WLS were the best."

> "It is noticeable that in our assessment, filters based on newer algorithms
> (PMF, CSF, SMRF) are outperformed by some of the oldest (PTIN, SBF), which
> contradicts recently published findings [46]."

**So CSF's signature in steep deciduous forest is: competitive on total and type
I error, among the worst on type II — it keeps non-ground.** This paper reports
*classification* error, not vertical bias, so the elevation consequence is not
measured; but the mechanism (retained non-ground sits above the surface) is
unambiguous in direction.

### 6.3 CSF against MCC on steep mountain forest (classification only)

**Fan, W.; Liu, X.; Zhang, Y.; Yue, D.; Wang, S.; Zhong, J. (2024).** *ISPRS
Annals* X-2-2024:73-79. DOI:
[10.5194/isprs-annals-X-2-2024-73-2024](https://doi.org/10.5194/isprs-annals-X-2-2024-73-2024).
Peer-reviewed. *Full text read (delegated search); I did not re-verify.* OpenGF
benchmark, 500 × 500 m scenes, including **S8 = steep mountain + sparse
vegetation and S9 = steep mountain + dense vegetation**. Reports OA / Kappa /
IoU only. **CSF (CloudCompare) loses to MCC in both steep-mountain scenes: OA
77.89 vs 83.11 (S8) and 88.93 vs 92.53 (S9).** Verbatim:
> "For steep mountainous forest terrain (Area S8, S9), SGSF shows a significant
> advantage over SGF and CSF algorithms, but its overall performance is not as
> effective as the MCC filtering algorithm designed specifically for forest
> environments."

This chimes with **Zhao, X.; Su, Y.; Li, W.; Hu, T.; Liu, J.; Guo, Q. (2018)**,
"A Comparison of LiDAR Filtering Algorithms in Vegetated Mountain Areas,"
*Canadian Journal of Remote Sensing* 44(4):287-298, DOI
[10.1080/07038992.2018.1481738](https://doi.org/10.1080/07038992.2018.1481738)
(*abstract only, paywalled; CSF not among the filters tested*), whose abstract
states: "The MCC works well in steep and dense forests; IBF and MCC outperform
the rest of filtering algorithms in areas with steep terrain but low vegetation
coverage; and **PTDF is more reliable for low-density LiDAR data.**"

**MCC (multiscale curvature classification) is available in lidR as `mcc()` but
is not in PDAL.** PDAL does ship `filters.smrf`, which is the filter Klápště
found best on ALS.

### 6.4 On the ISPRS benchmark, CSF loses to SMRF and to TIN densification

From Zhang et al.'s own Table 4 (§1.5), mean total error over the 15 samples:
**CSF 4.39** vs Pingel 2013 (SMRF) **2.97**, Hu 2014 **2.85**, Mongus 2014
**2.74**, Chen 2013 **4.11**, Axelsson 1999 (TIN densification, = TerraScan)
**4.82**. CSF is mid-pack and is beaten by SMRF, which PDAL also ships as
`filters.smrf`.

### 6.5 Software-maintainer opinion (not a benchmark)

The lidR book, current text (*software docs*, *full text read*):
> "Progressive TIN Densification (PTD) was proposed by Axelsson (2000). It the
> most widely used algorithms in industry and in closed-source software such as
> LAStools or TerraScan (with potentially undocumented variations). It generally
> outperforms PMF, CSF, and MCC, especially in complex terrain such as mountains."
> … "In early 2026, we added PTD to lidR and we now recommend using it
> exclusively."

This is a maintainer's claim, not a study, but it points the same way as Štular
& Lozić and Viedma: **progressive TIN densification beats CSF in complex
terrain.** Note that DeLong et al. 2022 — our closest published analogue for
repeat MN DNR lidar, already recorded in `DTM_PERCENTILE_LITERATURE.md` §5.2 —
used "industry standard TerraScan software v020", i.e. TIN densification, not
CSF.

### 6.6 What is established, and what is not

**Established.**
- CSF's characteristic classification failure in steep forest is **type II** —
  retaining non-ground (Štular & Lozić, verbatim).
- CSF-derived DEMs from ~0.5 pts/m² lidar in 9-26° terrain carry **significantly
  larger scatter (NMAD, RMSE) than progressive-TIN alternatives, before and
  after a correction surface** (Viedma, Dunn groups).
- CSF is beaten by SMRF and by several other filters on the ISPRS benchmark
  (Zhang et al.'s own Table 4).

- **CSF-derived DTMs read HIGH against surveyed ground control**: +0.12-0.13 m
  vs SMRF's +0.04 m on the same ALS cloud with per-algorithm tuning (Klápště,
  1414 GNSS checkpoints), and +0.09 to +0.16 m on flat vegetated floodplain
  (Szabó 2020, 604 RTK points, CSF only — sign convention verbatim: "We
  extracted the values of all models where ground control measurements were
  available (604 RTK points) and subtracted them from the measured values of the
  RTK points", so negative mean = model above ground).
- **CSF sits in the worse-performing group** — statistically indistinguishable
  from ATIN and ArcGIS, and 0.03-0.05 m worse in RMSE than PTIN, PMF and SMRF
  (Klápště).

**Not established.**
- **No signed CSF vertical bias has been measured in steep mature forest.**
  Klápště's site is a spoil heap with steppe, shrub and young deciduous woodland
  at 8 pts/m², and he attributes the slope-growing bias to *low vegetation on
  open slopes*, explicitly noting "Such vegetation is not present under canopy
  vegetation (i.e. in forests)". The one study that *does* run CSF at 28° and 37°
  under dense forest — **Zhang, S. et al. (2025), "Integration of Physical
  Features and Machine Learning: CSF-RF Framework…", *Sensors* 25(19):5950, DOI
  [10.3390/s25195950](https://doi.org/10.3390/s25195950)**, Tahoe National
  Forest — reports **classification errors only and never converts them to
  elevation** (*read via delegated web fetch; not re-verified by me*).
- **No study reports CSF vertical bias as a function of canopy cover.**
  Canopy-stratified DTM bias exists for other filters, not for CSF.
- **No study differences a CSF-derived DTM against an Axelsson/TerraScan- or
  `lasground`-derived DTM from the same cloud and reports the elevation offset.**
  That is exactly the artefact inherited when a CSF-classified epoch is compared
  against a class-2-classified epoch. Klápště's ME table is the only proxy.
- Viedma, the one head-to-head *vertical-error* comparison at our density,
  **does not report the CSF parameters it used**, so its CSF-vs-`lasground` gap
  cannot be cleanly attributed to the algorithm rather than to its settings.
## 7. Recommended parameter set for Elba / Whitewater bluffland

### 7.1 First, the framing the literature does not give you

Every study above optimises a **classification** metric — type I error, type II
error, total error, Kappa. Our objective is different: **an unbiased per-cell
median ground elevation on a 5 m grid, stable across slope and canopy classes**,
because the product is a difference of two epochs. Total error is a poor proxy
for that, and the two error types are not symmetric for us.

**[Inference, not from the literature — reasoned from the algorithm as described
in Zhang et al. §2.3-2.4, and flagged as such.]**

- **Type II** (non-ground retained as ground) adds points **above** the true
  surface. At 0.5 pts/m² a 5 m cell holds ~12 ground returns, so the median
  tolerates a minority of contaminants but shifts steadily as the contaminated
  fraction grows. On leaf-off November data the contaminants are crop residue
  and stubble on the uplands and understory/deadfall in the forest — both bias
  the median **high**, and both are *land-cover-correlated*, which is exactly the
  kind of bias a DoD cannot absorb.
- **Type I** (true ground discarded) biases the median only if the loss is
  systematic *in elevation within a cell*. On a planar steep face the cloth-to-
  ground offset is roughly constant across a 5 m cell, so type I loss there tends
  to be spatially coherent — whole-cell dropouts and thinning, not a within-cell
  bias. Near **convex breaks** (ridge crests, bluff shoulders), where the cloth
  bridges and its offset varies sharply within one cell, the retained ground
  should be biased **high**.

So the two error modes push the *same* direction — high — in the places where
they are worst.

**The measured evidence supports that direction, though two studies look at first
glance like they disagree.** The disagreement is only apparent, and resolving it
matters:

- **Klápště et al. 2020 (§2.1)** measured signed error against **1414 RTK-GNSS
  checkpoints** with each algorithm's parameters tuned separately. Verbatim:
  "All algorithms overestimated the terrain … ATIN, ArcGIS and CSF all
  overestimated the terrain by approx. 0.12–0.13 (± 0.13–0.14) m", against SMRF's
  0.04 m — and "the mean bias tended to increase with slope". **CSF reads HIGH,
  by roughly 3× the best available filter, and more so on slopes.**
- **Viedma 2022 (§6.1)** differenced 0.5 pts/m² lidar against a 300 pts/m² lidar
  benchmark over 9-26° slopes and found the CSF DEM the *lowest* of six filters
  (P50 −1.19 m). But **every** filter read low there, the filter-to-filter medians
  were **not statistically significant** (all Dunn group "a"), and the reference
  was another lidar dataset rather than surveyed ground.

**Klápště is the better test of the filter's own contribution** — absolute
reference, per-algorithm tuning, 1414 checkpoints, mixed model with conditional
R² = 0.67. Viedma's negative offset is a property of his low-versus-high-density
comparison and his geoid-distance term, shared by all six filters. What Viedma
*does* establish about CSF specifically is **scatter**: roughly twice the NMAD and
RMSE of the TIN variants, before and after a correction surface.

**So: CSF's own vertical signature is high-biased and noisy, and the bias grows
with slope.** Caveat the transfer honestly — Klápště's site is a vegetated spoil
heap at 8 pts/m², and he attributes the slope-growing bias to *low vegetation on
open slopes*, adding that "Such vegetation is not present under canopy vegetation
(i.e. in forests)". Our steep ground is under canopy; our flat ground is November
crop stubble.

**A structural caution specific to this pipeline.** CSF is applied to the gen1
cloud only (gen2 uses ASPRS class 2). Any change to these parameters therefore
moves **one epoch only**, and propagates directly into the DoD. That is the
comparability trap: a parameter change is not neutral here even if it improves
classification, because it changes the difference. Any change must be evaluated
on the DoD stratified by slope and cover, not on a single-epoch classification
score.

### 7.2 The recommendations

| Parameter (PDAL name) | Recommended | Currently | Confidence | Primary evidence |
|---|---|---|---|---|
| `step` (time_step) | **0.65** | 0.65 | **High** | Zhang swept 0.4-1.5 in 0.05 steps; minimum total error at 0.65 for all three scene groups. No applied study varies it. |
| `threshold` (class_threshold) | **0.5** | 0.5 | **High** | Zhang: justified ("OBJ measurements … are usually taller than 0.5 m"), swept, "limited impact on total errors". Štular tested lower values and kept 0.5 "on average or rugged terrain". |
| `smooth` (slope post-processing) | **true** | true | **Medium-high** | Zhang turns it on for *both* groups with steep slopes (II and III); it is his preferred alternative to softening the cloth. Štular: on at all four steep sites. Micu: "On". |
| `hdiff` | **0.3** | 0.3 | **Medium** (untested, but low-stakes) | Zhang's fixed value; "only influence the movable particles over steep slopes, the influence is also limited". Nobody has tested it. |
| `iterations` | **500**, or 1000 | 500 | **Medium** | Zhang shows M_HV converged by ~150 iterations. Štular used 1000 at all four steep sites: "modest improvements at the expense of a linear increase of the processing time". Low risk either way. |
| `resolution` (cloth_resolution) | **keep 1.0; test 2.0, not 0.5** | 1.0 (inherited) | **Medium** | This is the dominant parameter (Klápště, explicit ranking; LiDAR RMSE 0.18-0.53 m across the grid). Finer cloth = more positive bias; coarser cloth = "almost zero bias" with more noise. Our metric is bias. See §7.3. |
| `rigidness` | **2**, and test {1, 2, 3} | 1 | **Low-medium** | See below. This is the one parameter the literature cannot settle for us, and nobody has swept it. |

### 7.3 The two that need argument

**`resolution`: it is the dominant parameter, and the bias-minimising direction
is *coarser*, not finer.**

*I changed my mind on this while writing, and the reason is worth recording.* My
first reading favoured moving from 1.0 to Zhang's 0.5, because 0.5 is the only
value chosen by an explicit sweep at our density (0.18-0.67 pts/m², Optech ALTM
— §4.2) rather than inherited as a library default. That argument is about
**classification total error**, which is what Zhang optimised.

Klápště et al. optimised something else — **signed vertical error against 1414
GNSS checkpoints** — and it points the other way:

> "The greatest effect was observed for the Cloth resolution parameter, with the
> value of 0.1 m leading to the best performances in terms of random-error and
> bias combination. It should be noted, however, that for LiDAR point cloud, the
> parameter Cloth resolution 2 m … lead to almost zero bias, although the random
> error was relatively high."

So along the cloth-resolution axis there is a **bias-versus-noise trade**: fine
cloth conforms closely, drapes onto low vegetation, and carries more **positive**
bias; coarse cloth bridges over it and approaches zero bias at the cost of
scatter. Cai et al. 2019 report the same direction in classification terms on the
ISPRS data (finer cloth → lower type I, higher type II), and type II is the error
that lifts a ground surface.

**We are the case that should take the low-bias end of that trade.** We median
~12 ground returns per 5 m cell and then difference two epochs; per-point scatter
is largely absorbed by the median, while a bias is not — and worse, our CSF stage
touches gen1 only, so any bias lands directly in the DoD (§7.1).

**Therefore: keep `resolution = 1.0`, and if you test anything, test 2.0 rather
than 0.5.**

**Two limits on that, both real.**

1. **Klápště's "2 m → almost zero bias" was measured at 8 pts/m².** There, a 2 m
   cloth cell contains ~32 points. At our 0.5 pts/m² ground density a 2 m cloth
   cell contains ~2 ground returns. **[Inference]** Coarsening has a floor set by
   density: past some point the cloth stops finding ground at all, type I rises,
   cells thin and then drop out — and on steep slopes, which are where we can
   least afford it. The published low-density steep-forest data point (Štular's
   ES site, 1.83 pts/m²) chose **1.0**, not 2.0. That is the strongest argument
   for staying where we are.
2. **The noise coarsening buys is not free at 12 points per cell.** Losing ground
   returns thins the very cells whose median we depend on.

So: `1.0` is defensible on the evidence, `2.0` is the experiment worth running,
and `0.5` is the direction to *avoid* if bias is the objective.

**`rigidness`: 1 is probably too soft, and no single value is right.**

The case against the current `rigidness=1`:

1. **Zhang reserves 1 for terrain we do not have.** His group III is "High and
   steep slopes (e.g., pit, cliff)". A driftless bluffland is his **group II** —
   "With steep or terraced slopes (e.g., river bank, ditch, terrace)" — for which
   he prescribes **RI = 2, ST = true**.
2. **Zhang explicitly prefers the opposite remedy.** "A direct method to mitigate
   this problem is to set the rigidness to a lower value, but some low objects
   may be classified as BE as a result. … Thus, we can use a relatively hard
   cloth and post-processing to remove lower objects and correctly handle steep
   slope areas."
3. **The only ALS study spanning our density range never found soft best.**
   Štular & Lozić's Table A1 gives **hard** at three of four steep deciduous
   forested sites and **medium** at the sparsest; soft was not even in their
   tested range.
4. **Soft cloth maximises exactly the error that biases a median.** Rigidness 1
   gives the highest type II, i.e. the most low vegetation retained as ground —
   and that is the error mode that moves the per-cell median, in a
   land-cover-correlated way.
5. **CSF's characteristic failure mode is exactly the one that biases a median.**
   Štular & Lozić rank nine filters on the same steep deciduous-forested ALS
   data and find, verbatim: "The best score in total error was obtained by SMRF,
   followed by a close group of BMHF, PTIN, SBF, and CSF. … The T1 error results
   are almost identical to total error. … The T2 error is very different, though.
   MCC and WLS are decidedly the best, followed by a narrow group of PTIN, BMHF,
   and PMF. **SBF, CSF, SegBF, and at some distance SMRF are the worst.**" CSF is
   near the top on total/type I and near the bottom on type II — it keeps
   non-ground. Softening the cloth pushes further in that direction.
6. **40% of our cells are 0-5°**, where every source says hard.
7. Wan et al.'s mapping puts rigidness 2 at "gently undulating" and 1 at
   "sharply undulating"; a tile that is 40% under 5° is not, on a scene average,
   sharply undulating.

The case for 1: 25% of our cells exceed 20°, and those are the cells whose
ground we can least afford to lose.

**Recommendation: `rigidness = 2` with `smooth = true`** as the single-value
default. That is Zhang's own group-II prescription, it is CloudCompare's default
and its "Relief" scene, and it is the middle of the range Štular found best on
steep forested ALS.

**But treat this as a hypothesis to test, not a settled answer**, because the
literature is explicit that a single rigidness cannot be right for a scene like
ours — Cai et al. 2023, verbatim: "a single cloth rigidness is unreasonable,
because multiple terrain features (i.e., flat terrain, slopes and raised terrain)
are usually contained in a landscape". Run `rigidness` ∈ {1, 2, 3} and compare on
the DoD, stratified by slope band and canopy class.

**And do not be reassured by the one number that looks reassuring.** Wan et al.'s
headline — automatic scene-level rigidness selection costs 0.36 percentage points
of mean total error against hand-tuning (4.758% vs 4.394%) — is a *mean* over 15
samples, and the same paper reports that "27% of the samples had larger biases
(>1% in Total Error)". A quarter of scenes moved by more than a full percentage
point. **Nobody has swept rigidness against error at all** (Wan has no such
table or figure), and **nobody has measured its effect on vertical bias.** The
leverage of this parameter for our purpose is simply unknown.

### 7.4 One thing worth testing that is not a parameter

Two independent groups found that **levelling the cloud before running CSF on
steep terrain** materially improves it (Štroner et al. 2021: total error 16.0% →
6.2%; Marotta et al. 2021: tilting by ±30/45/60° "to preserve points also in the
steepest slopes"). The mechanism is structural — the cloth falls under gravity in
the cloud's own frame, so on a steep face the algorithm works against its own
assumption. Neither study is ALS-under-canopy and neither terrain is as gentle as
a bluffland valley wall, so this is a **candidate to test, not a recommendation**.
It is also awkward here: levelling is a per-region operation, and doing it
per-region reintroduces exactly the tiling problem of §5.

### 7.5 What not to change

- **Do not raise `threshold` above 0.5.** The module docstring records that an
  earlier version ran `threshold=1.5` and that it "retained sub-ground returns
  and biased the gen1 ground low under canopy". That is precisely the direction
  the literature predicts for a wider threshold, and Zhang's justification for
  0.5 (objects are usually taller than 0.5 m) does not survive widening.
- **Do not turn `smooth` off** — with one caveat now on the record. It is the
  authors' designated remedy for the steep-slope failure mode, and every
  steep-terrain study that reports it has it on. Its known cost — Štular,
  verbatim: "noticeable T2 errors, e.g., 'tree stumps' on slopes" — is real but is
  a type II cost, which a median over 12 points is comparatively well placed to
  absorb. **The caveat:** Klápště's near-zero-bias combination was "Cloth
  resolution 2 m in combination with **either Slope processing False, or Slope
  processing True and Ruggedness 2 or 3**" — so at coarse cloth, slope
  post-processing off *also* reached low bias. If you test `resolution = 2.0`,
  test `smooth` both ways alongside it rather than holding it fixed.
- **Keep PDAL's `returns` default ("last, only").** It is a PDAL-specific
  behaviour with no counterpart in the paper or the other wrappers, and it is the
  right one for bare earth on leaf-off data.
## 8. Gaps: what the literature does not say

Listed roughly in order of how much they cost us.

1. **No CSF parameterisation exists for driftless / loess / dissected-plateau
   terrain.** A full-text query for CSF combined with loess, gully, badland,
   ravine and dissected terms across the ~1,560 works citing Zhang et al. 2016
   returned zero hits. The one forested-gully study located (Manić et al. 2022)
   reports no parameter values at all. This is verified silence, not a search
   failure.

2. **Nobody has published CSF parameters for ~0.5 pts/m² ALS in steep forest.**
   The closest verified points are Štular & Lozić's ES site (1.83 pts/m², steep
   deciduous, best `r = 1.0`, scene = medium) and Viedma's PNOA data (0.5
   pts/m², 9-26° slopes) — and Viedma **does not report his CSF settings**.
   Zhang's own tuning was done at 0.08-1 pts/m², which is the right density, but
   on European urban/rural benchmark scenes, not on a forested dissected
   plateau.

3. **No signed CSF vertical bias exists for steep mature forest, and none exists
   as a function of canopy cover.** Klápště et al. 2020 give the only signed,
   absolutely-referenced CSF bias (+0.12-0.13 m, growing with slope) — on a
   vegetated spoil heap at 8 pts/m² whose slope effect they attribute to *low
   vegetation on open slopes*, a mechanism they say does not operate under
   canopy. The one study that runs CSF at 28° and 37° under dense forest (Zhang
   et al. 2025, Tahoe National Forest) reports classification errors and never
   converts them to elevation. **The measurement we need has not been made.**

4. **No guidance for a scene that is bimodal in slope.** The literature states
   the problem plainly — Cai et al. 2023: "a single cloth rigidness is
   unreasonable, because multiple terrain features … are usually contained in a
   landscape" — but every published response is a *new algorithm* (ICSF; Yang et
   al.'s region partitioning; Wan et al.'s scene-level relief index; the
   multi-resolution hierarchical variants). **There is no published
   recommendation for how to set stock CSF's single `rigidness` on such a scene,
   and no published test of tiling a scene by slope and running CSF twice.**

5. **`rigidness` in steep forest is genuinely contested, not merely unstudied.**
   Zhang's Table 2 says 1 for "high and steep slopes"; Zhang's §4.3 says prefer a
   hard cloth plus post-processing; Štular & Lozić found hard best at three of
   four steep forested sites (and never tested soft); Storch et al. used 1 with a
   stated reason; Bailey et al. found rigidness had no effect at all. No study
   isolates rigidness on low-density ALS in steep forest.

6. **`hdiff` (Zhang's h_cp) has never been tested by anyone**, including the
   original authors — "We simply set this parameter to 0.3 m". It only acts when
   slope post-processing is on, which for us is always.

7. **Zhang's cloth-resolution sweep range is not reported numerically.** The
   time-step sweep is given as "from 0.4 to 1.5 with steps of 0.05"; the grid
   resolution sweep is only "different GR values" with the result read off an
   unlabelled figure. How far either side of 0.5 was tested is unknown.

8. **Leaf state is almost never reported.** Of the studies with usable
   parameters, only Storch et al. state leaf condition — and they used the *same*
   CSF parameters for leaf-on and leaf-off. Štular & Lozić's four deciduous ALS
   sites never mention it (verified: no occurrence of "leaf" in the full text).
   **There is no published evidence that CSF parameters should differ between
   leaf-on and leaf-off acquisitions.**

9. **No study evaluates CSF against the metric we care about** — the stability of
   a per-cell median ground elevation across slope and cover strata. All optimise
   classification error, or a DEM-wide RMSE.

10. **Implementation provenance is undocumented.** The upstream repository warns
    that "This code has been changed a lot since the publication of the
    corresponding paper", and none of PDAL, CloudCompare or RCSF records which
    revision it vendors. Benchmark numbers from 2016 may not describe the binary
    in the conda env.

### Sources I could not obtain

- **Klápště, P.; Fogl, M.; Barták, V.; Gdulová, K.; Urban, R.; Moudrý, V.
  (2020).** "Sensitivity analysis of parameters and contrasting performance of
  ground filtering algorithms with UAV photogrammetry-based and LiDAR point
  clouds." *International Journal of Digital Earth* 13(12):1672-1694. DOI:
  [10.1080/17538947.2020.1791267](https://doi.org/10.1080/17538947.2020.1791267).
  Peer-reviewed. **Abstract only** — tandfonline returns HTTP 403 to WebFetch and
  to curl with a browser user-agent; DOAJ lists it as OA but its only fulltext
  link is the blocked DOI. On paper this is the single most relevant sensitivity
  study: it compares six filters (CloudCompare's CSF among them) on **leaf-off**
  UAV photogrammetry and **leaf-on** airborne lidar of the same area, and
  explicitly evaluates "the effect of vegetation density and terrain slope on
  filtering accuracy". Verbatim from the abstract:
  > "Our results show that the performance of filtering algorithms was affected
  > by the point cloud type, terrain slope and vegetation cover. The results were
  > generally better for LiDAR (RMSE 0.13–0.19 m) than for photogrammetric (RMSE
  > 0.19–0.23 m) point clouds. … **Parameters related to the selection of the
  > initial minimum elevation ground points were the most influential in all
  > algorithms and point clouds.**"

  **No CSF parameter values retrieved.** Worth requesting through a library.

- **Sithole, G.; Vosselman, G. (2004).** "Experimental comparison of filter
  algorithms for bare-Earth extraction from airborne laser scanning point
  clouds." *ISPRS J. Photogramm. Remote Sens.* 59(1-2):85-101. DOI:
  [10.1016/j.isprsjprs.2004.05.004](https://doi.org/10.1016/j.isprsjprs.2004.05.004).
  **Abstract only** (Elsevier paywall). The 2003 ISPRS proceedings version of the
  same experiment is open and was read in full; the point-spacing figures quoted
  in §4.2 come from it.

- **Yang, A. et al. (2020).** "Filtering of airborne LiDAR bathymetry based on
  bidirectional cloth simulation." *ISPRS J. Photogramm. Remote Sens.* 163:49-61.
  **NOT VERIFIED** — known only from Cai et al. 2023's reference list and their
  one-sentence description of it. It is bathymetry, so transferability is low
  regardless, but it is the origin of the region-partitioning idea in §5.

- **Yilmaz, V. (2021)** — a metaheuristic optimisation of CSF parameters (Grey
  Wolf Optimizer / Jaya). **NOT VERIFIED — abstract only.** Every open route was
  exhausted: ScienceDirect (Cloudflare captcha), OpenAlex / Unpaywall / Semantic
  Scholar (all report CLOSED, zero OA locations), CORE, BASE,
  scholar.archive.org, fatcat, ADS gateways, and the author's institutional
  repositories at Karadeniz Technical University and Artvin (the KTU AVESİS
  record is metadata-only with no file attached). Its search bounds, objective
  function and per-site error numbers are unobtainable without interlibrary loan
  or contacting the author. **This is the paper that would tell us what parameter
  *ranges* practitioners consider plausible**, and it is the one real hole left in
  this survey.

- **Grube, G.; Talbot, B.; Grigolato, S. (2026).** "How much data is enough?
  Sensor choice, scale, and point-density effects on terrain metrics in steep
  mountain forests." *Scandinavian Journal of Forest Research*. DOI:
  [10.1080/02827581.2026.2715106](https://doi.org/10.1080/02827581.2026.2715106).
  **Abstract only** (tandfonline 403). On-topic for density effects in steep
  closed-canopy forest; whether it reports CSF parameters is unknown.

### Method notes for anyone repeating this search

- `www.mdpi.com` returns 403 to both WebFetch and curl. The article PDFs are
  served without challenge from
  `https://mdpi-res.com/d_attachment/<journal-slug>/<slug>-<vol4>-<art5>/article_deploy/<slug>-<vol4>-<art5>.pdf`
  (e.g. `remotesensing-08-00501`, `forests-14-01520`, `sensors-22-00217`).
- `pdal.io` and `research.utwente.nl` 403 WebFetch; the former yields to curl
  with a browser user-agent, the latter does not (Cloudflare). The Twente
  author copy *is* reachable at `https://ris.utwente.nl/ws/files/<id>/<name>.pdf`.
- `tandfonline.com` and `sciencedirect.com` defeated every route tried.
- `https://api.openalex.org/works/doi:<DOI>` and
  `https://api.semanticscholar.org/graph/v1/paper/DOI:<DOI>?fields=...` are the
  fastest way to confirm a DOI resolves to the title claimed, and to find
  repository copies (`locations[].pdf_url`).

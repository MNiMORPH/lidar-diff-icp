# gen2's absolute datum: the 2021 project publishes per-point residuals

**Date:** 2026-08-27
**Driver:** `analysis/gen2_absolute_datum.py` – one run prints every table below
**Parser:** `analysis/groundtruth/parse_gen2_control.py` – rebuilds the bundled CSV
**Data:** `src/lidar_diff_icp/groundtruth/data/mn_se_driftless_2021_control.csv`
(**534** marks; **390** of them carry a residual)
**Ledger:** `.trust/runs/20260827T072813-1731650.json`
**Bears on:** `analysis/CONTROL_RESIDUAL_FIELD.md` (gen1's side of the same measurement),
`analysis/groundtruth/gen2_checkpoint_tie.py` (the six-mark estimate this supersedes)

Sign convention throughout, unchanged from `groundtruth.tie`: **positive = the surface
reads LOW**, i.e. the number is the constant to ADD. Section 1 re-derives it from the
data rather than trusting a column name.

Every number below is pasted from the driver's output. Nothing was downloaded but
documents; no 3DEP tile was fetched.

---

## The answer to the question that was asked

**Yes. The 2021 project publishes a per-point residual for every held-out checkpoint,
and the file had not been opened.** It is not in the survey report, not in the
contractor's checkpoint shapefile and not in the VA text files – all three were checked
before, and all three genuinely lack it. It is in the **USGS** side of the same
directory, in the NGTOC "VATool" output shapefiles:

```
https://rockyweb.usgs.gov/vdelivery/Datasets/Staged/Elevation/metadata/
    MN_SE_Driftless_2021_B21/Vertical_Accuracy/USGS/
    USGS_MN_SE_Driftless_2021_B21_QL{0,1}.{shp,shx,dbf,prj}
```

whose `.dbf` carries, one row per checkpoint:

```
proj_name, srcChkptId, X, Y, Z, Lndcover, DEMz, zdiff, zdiffSq, LAZz, LAZzdiff, LAZzdiffSq
```

`Z` is the surveyed orthometric height, **`DEMz` is the delivered OPR DEM read at the
mark** and **`LAZz` the delivered classified point cloud read at the mark**. So gen2's
field comes for free, with **two** delivered surfaces where gen1's validation reports
give one, and on exactly the marks that are *not* circular: the VATool tested only the
NVA/VVA checkpoints, never the LCPs that calibrated the data.

**390 unique checkpoints, 395 rows** (five marks were tested against both quality-level
blocks). Item 1 makes items 3–4 answerable without touching a point cloud, exactly as
hoped.

---

## Bottom line

1. **gen2's level on open ground is indistinguishable from zero, and it is now measured
   on 227 marks instead of four.** Against the delivered point cloud, project-wide:

   | class | n | mean (mm) | median (mm) | sd (mm) | SE of the mean (mm) |
   |---|---|---|---|---|---|
   | **NVA (non-vegetated)** | 227 | **−2.37** | −0.32 | 35.66 | **2.37** |
   | VVA (vegetated) | 163 | −74.60 | −58.70 | 85.45 | 6.69 |

   The QL1 block, the one containing Elba, gives NVA **−6.83 ± 2.96 mm** and VVA
   **−99.14 ± 9.48 mm**. This is what the vendor bias adjustment was for: gen2's open
   ground was tuned onto its control and it sits there.

2. **The canopy discriminator is real, large, and it points the way a leaf-on flight
   should.** NVA minus VVA, on the point cloud:

   | block | n NVA | n VVA | difference (mm) | SE (mm) | Welch t | df |
   |---|---|---|---|---|---|---|
   | all blocks | 227 | 163 | **+72.23** | 7.10 | **+10.174** | 202.8 |
   | QL1 | 139 | 99 | **+92.31** | 9.93 | +9.293 | 117.2 |
   | QL0 | 88 | 64 | +41.30 | 7.35 | +5.619 | 108.0 |

   The delivered 2021 ground reads **72 mm higher under vegetation than on open ground**,
   and **92 mm** in the block Elba sits in. The six-mark estimate this supersedes –
   `gen2_checkpoint_tie.py`, re-run this session, NVA median **+17.2 mm** on four marks
   against VVA median **−40.6 mm** on two – had the direction right and put the contrast
   at 57.8 mm, on a sample far too small to call it significant.

3. **Our own estimator and USGS's disagree, at the same six marks, on the same delivered
   cloud, by more than the level being measured.** `estimate_tie` minus `LAZzdiff`:
   mean **+5.3 mm**, sd **55.4 mm**, SE **22.6 mm**, r = **+0.453**, n = 6. That is the
   honest size of the estimator-plus-siting term, and it is the reason to prefer the
   390-mark published field over six locally measured ties.

4. **The point cloud on disk reaches almost none of these marks, and does not need to.**
   Of the 390 residual-bearing checkpoints, **exactly six** have a gen2 box on disk and
   **384 do not**. The only control point inside `data/after/elbaext_*` is LCP 1079,
   which is a **calibration** point and therefore circular. Fetching the rest would cost
   0.7–2.3 GB per tile to re-measure a quantity USGS has already published.

5. **The epoch difference is now bounded on the gen2 side, and the gen1 side is what is
   loose.** Like-for-like, delivered vendor DEM against delivered vendor DEM,
   project-wide means:

   | stratum | gen1 (mm) | gen2 (mm) | gen1 − gen2 (mm) | n |
   |---|---|---|---|---|
   | open / NVA | +13.79 | −3.38 | **+17.16 ± 7.50** | 230 / 227 |
   | vegetated / VVA | −86.51 | −74.15 | **−12.36 ± 8.90** | 534 / 163 |

   gen2's contribution to that difference now carries an SE of 2.3 mm; gen1's carries
   7.1 mm on open cover and its *cover treatment* is still worth 71 mm
   (`CONTROL_RESIDUAL_FIELD.md`). **The weaker side is no longer gen2.**

6. **Per-line structure at Elba is not measurable from the data on disk, and that is the
   result.** At five of the six boxed marks only **one** flight line puts returns inside
   the fitting radius; the sixth's second line contributes 95 returns and returns
   +2131 mm with a σ of 720 mm, which is a swath-edge artefact and not a line offset.

---

## 1. The sign, re-derived on every row

The driver tests both subtractions on all 395 rows of both blocks, at a caller-supplied
tolerance of 1e-9 m:

```
  dem QL0: surveyed-surface 157/157 rows (max |resid| 4.72e-16 m)   surface-surveyed 0/157 (max 3.86e-01 m)
  dem QL1: surveyed-surface 238/238 rows (max |resid| 4.72e-16 m)   surface-surveyed 0/238 (max 1.06e+00 m)
  laz QL0: surveyed-surface 157/157 rows (max |resid| 5.00e-16 m)   surface-surveyed 0/157 (max 4.00e-01 m)
  laz QL1: surveyed-surface 238/238 rows (max |resid| 5.00e-16 m)   surface-surveyed 0/238 (max 1.10e+00 m)
```

`zdiff` and `LAZzdiff` are **surveyed minus surface**. That is the same sign family as
`tie = surveyed − z_lidar` and as gen1's `dnr_error_m`, so the two epochs' residuals can
be differenced with no flip.

**The geoid.** All 534 marks read NAVD88 / GEOID18, asserted **per mark** from the
contractor shapefile's `geoid` attribute on the 390 that appear in one, and from the
report's own per-table header for the 144 that do not. gen2 is delivered on GEOID18, so
the geoid term is exactly zero and the driver raises rather than converting.

---

## 2. The control set, recovered in full

**534 points: 143 LCP, 227 NVA, 164 VVA.** That matches the survey report's §1.3 text
exactly, which no previous count did. Three things had to be got right:

1. The report's coordinate tables hold **533**, not 532. One NVA id carries a letter
   suffix, `2198A_2022_MN`, and the regex in `ADDITIONAL_GROUND_CONTROL.md` §1.2 has no
   `[A-Z]?` in it, so it drops that mark silently.
2. The 164th VVA, **`3000_2021_MN`**, is missing from the report's tables altogether and
   is recovered from the USGS shapefile. Report plus shapefile together are complete;
   neither is complete alone.
3. One report-table VVA, `3021_2021_MN`, was never tested by the VATool, so it has
   coordinates and no residual.

**The two sources agree exactly.** On the 389 marks they share, max |dE| = 0.0000 m,
max |dN| = 0.0000 m, max |dZ| = 0.0000 m, and the land-cover code matches on every one.

**The parse validates itself against all eight published aggregates:**

| block | surface | n NVA | RMSEz cm | published | n VVA | 95th pct cm | published |
|---|---|---|---|---|---|---|---|
| QL1 | DEM | 139 | 3.5094 | 3.51 | 99 | 25.4797 | 25.48 |
| QL1 | LAZ | 139 | 3.5395 | 3.54 | 99 | 27.1359 | 27.14 |
| QL0 | DEM | 91 | 3.5093 | 3.51 | 66 | 13.3262 | 13.33 |
| QL0 | LAZ | 91 | 3.5535 | 3.55 | 66 | 12.6388 | 12.64 |

The VVA figure is the 95th percentile of |residual| and numpy's default `linear`
interpolation is the one that reproduces it; the alternatives are printed beside it
(QL1 LAZ: lower 26.91, higher 29.19, nearest 26.91, midpoint 28.05) rather than one
being chosen quietly.

**`role` is the column that keeps the check honest.** The vendor FGDC metadata says the
LCPs calibrated the lidar and the NVA/VVA checkpoints "were not used to calibrate or
post process the data". The 143 LCPs therefore carry `role=calibration`, carry no
residual, and enter **no** statement below except the coverage count.

---

## 3. gen2's level

Against the delivered classified point cloud (`laz`), one row per mark, marks tested in
both blocks taken from QL1:

```
  surface       block  class    n  mean_mm  median_mm  sd_mm  se_mm  rmse_mm
      laz  all blocks    all  390   -32.56     -19.92  71.08   3.60    78.10
      laz  all blocks    NVA  227    -2.37      -0.32  35.66   2.37    35.66
      laz  all blocks    VVA  163   -74.60     -58.70  85.45   6.69   113.24
      laz         QL1    all  238   -45.23     -27.60  80.43   5.21    92.12
      laz         QL1    NVA  139    -6.83      -7.89  34.86   2.96    35.39
      laz         QL1    VVA   99   -99.14     -90.63  94.36   9.48   136.54
      laz         QL0    all  152   -12.72      -8.73  47.05   3.82    48.59
      laz         QL0    NVA   88    +4.66      +3.51  35.98   3.84    36.08
      laz         QL0    VVA   64   -36.63     -32.39  50.16   6.27    61.79
```

Against the delivered OPR DEM the same table reads NVA **−3.38 ± 2.34**, VVA
**−74.15 ± 6.72** over all blocks, and NVA **−7.73 ± 2.91**, VVA **−99.34 ± 9.37** in
QL1. **The two surfaces agree to within 1.5 mm in every stratum**, which says the residual
is a property of the classified ground and not of the gridding.

**Read `se_mm` for what it is.** It is `sd/sqrt(n)` of the mean over *these* marks. It
is not the uncertainty of gen2's level at Elba; §5 is that.

**A pooled number over all 390 marks is meaningless here** and is printed only to show
why: at −32.56 mm it is a weighted average of two populations 72 mm apart, and its value
tracks the NVA/VVA ratio of whatever mark set is pooled. The QL1-vs-QL0 gap in the `all`
row (−45.23 against −12.72) is mostly that ratio, not a block offset: the NVA rows differ
by only 11.5 mm.

---

## 4. NVA versus VVA – the canopy discriminator

This is the cleanest cover contrast available in either epoch, because the surveyors
assigned the classes in the field, before any lidar was processed.

```
  surface       block    n    n  diff_mm  se_diff_mm        t     df
      laz  all blocks  227  163   +72.23        7.10  +10.174  202.8
      laz         QL1  139   99   +92.31        9.93   +9.293  117.2
      laz         QL0   88   64   +41.30        7.35   +5.619  108.0
      dem  all blocks  227  163   +70.77        7.11   +9.951  201.5
      dem         QL1  139   99   +91.61        9.81   +9.337  117.1
      dem         QL0   88   64   +38.68        7.68   +5.035  102.5
```

**Measurement, stated bare:** the delivered 2021 ground surface sits 72 mm higher under
vegetation than on open ground project-wide, and 92 mm higher in the QL1 block.

*Interpretation, marked as such and not tested here:* this is the direction the leaf-on
acquisition would produce – gen2 was flown 2021-05-01 at green-up, NDVI 0.49,
contradicting its own vendor "leaf-off" spec (that flight-date and NDVI finding is
quoted, not re-checked here) – and it carries the same sign as the forest-versus-open
differential the return-structure work reports. **It is not proof of that mechanism.** VVA marks also sit on rougher, softer
ground than NVA marks, and this contrast cannot separate the two.

**The QL1/QL0 gap is worth noticing and I have not chased it.** The same contrast is
92 mm in one block and 41 mm in the other, and the two blocks were flown by different
sensors in different seasons. That is a testable question the collection dates in the
bundled CSV would answer.

---

## 5. The residual as a spatial field at Elba

Fitted and kriged with the **same** `groundtruth.residual_field` machinery and the
**same** swept nuisance grid that produced gen1's site value in
`CONTROL_RESIDUAL_FIELD.md` – `max_lag` 20/40/60/100/184 km, `n_lags` 15 and 30,
`n_pairs` 200 000 and 800 000, estimators dowd and matheron, seeds 0/1/2. Ten rows per
variant; the median and the full min–max of the sweep:

```
surf  variant                  nrow  pred med      pred min..max  sdfield med  sdmark med
laz   NVA only                   10      -2.3       -2.4..    +1.1         25.0        35.5
laz   VVA only                   10     -71.5      -80.4..   -44.6         39.8        52.5
laz   class covariate -> NVA     10      -1.8       -4.9..    +2.9         22.1        52.6
laz   class covariate -> VVA     10     -72.5      -76.2..   -67.7         22.2        52.8
dem   NVA only                   10      -3.3       -3.4..    +2.5         24.1        34.7
dem   VVA only                   10     -67.5      -75.2..   -49.7         40.2        56.5
dem   class covariate -> NVA     10      +1.6       -4.4..    +6.2         23.2        52.8
dem   class covariate -> VVA     10     -67.7      -77.0..   -63.0         23.4        52.9
```

Three things follow, and one of them is the opposite of gen1's story.

1. **The kriging barely moves the answer.** gen2's NVA prediction at Elba, −2.3 mm, is
   within 0.1 mm of the plain NVA mean over all 227 marks. On gen1 the same exercise
   moved the site value by tens of millimetres. gen2's open-ground residual has almost
   no spatial structure to exploit – which is what a well-adjusted delivery should look
   like.
2. **The cover treatment is still the whole estimate.** NVA-only and VVA-only differ by
   69 mm at Elba, dwarfing every `sd_field` in the table. gen1's three treatments spanned
   71 mm. **Both epochs' site values are set by the same unresolved choice, and it is not
   a modelling choice – it is the question of which cover stratum defines the datum.**
3. **`sd_field` and `sd_mark` are different quantities and both are printed.** The
   ~22–25 mm `sd_field` is the error in predicting the correlated part of the field at
   Elba; the ~35 mm (NVA) to ~53 mm (pooled) `sd_mark` is what a single new mark placed
   there would scatter by. Do not quote the first as if it were the second.

The fitted variogram parameters are in the driver's output. **The NVA field's range is
not determined by this control set**: over the sweep it runs 3 892 m to 80 355 m with
nuggets from 17 to 1 192 mm², i.e. it tracks `max_lag` rather than the data, which is the
same pathology gen1's field showed and is reported for the same reason.

---

## 6. Our estimator against USGS's, at the six marks we can reach

Both read the **same** delivered class-2 cloud at the **same** mark, so any difference is
estimator plus siting, with no epoch and no datum in it.

```
      point_id  class  usgs_mm  ours_mm  delta_mm  sigma_mm
  2210_2021_MN    NVA    +23.5     -0.7     -24.2      11.5
  3056_2021_MN    VVA    -18.7    +27.9     +46.6      20.7
  2024_2021_MN    NVA    -22.0    +35.1     +57.1      39.4
  2036_2021_MN    NVA    +15.4    +59.4     +44.0      27.4
  2099_2021_MN    NVA    -38.8    -43.7      -4.9       8.7
  3089_2021_MN    VVA    -22.4   -109.2     -86.7      36.9
  n=6  mean(ours-USGS) +5.3 mm  median +19.6 mm  sd 55.4 mm  RMS 50.8 mm
  correlation r = +0.453
```

**This is a negative result and it matters.** Two competent estimators, on identical
data, at identical points, scatter against each other with a 55 mm sd – larger than
gen2's entire open-ground level. The six-mark inverse-variance weighting was flagged as invalid in the brief that set
this task (χ² = 31.9 on 5 dof, p = 6.30e-06, inflating to −5 ± 24 mm; **UNVERIFIED
here** – those three numbers come from the brief and I did not recompute them). It was
reading this scatter. **USGS's own residual at those same six marks is −10.50 ± 9.95 mm**
(mean ± SE; median −20.35, sd 24.37), and over all 227 NVA marks it is −2.37 ± 2.37 mm.

`estimate_tie` is not thereby wrong. It answers a different question – it fits a local
quadratic and reads it at the mark, deliberately, because a plane walks by a metre with
radius on these marks – whereas the VATool reads the delivered surface. But **carrying
gen2's published level into a comparison with a number our estimator produced needs the
+5.3 ± 22.6 mm bridge above**, and 22.6 mm on n = 6 is not a bridge worth much weight.

**The 143 LCPs, which is where the temptation lies.** LCP 1079 sits inside
`elbaext_3dep_fd_class2.laz`, is 2.18 km from the reference point, and has both epochs on
disk. It is a **calibration** point. `ADDITIONAL_GROUND_CONTROL.md` §1b already measured
gen2 there at −11 to −26 mm and read it as evidence about siting; as evidence about
**accuracy** it is circular, and it is excluded here.

---

## 7. Per-line structure

gen2's lines over Elba run **east–west** (bearing 90.0–90.8°, from the vendor swath
polygons – `LIDAR_DOCUMENTATION_MINE.md` §2.4), perpendicular to gen1's north–south
lines. Our pipeline reads gen2's delivered class 2 and performs **no** swath alignment on
it, so anything here is descriptive.

```
      point_id  class  line   n_line  ours_mm  sigma_mm     n
  2210_2021_MN    NVA  3041    81281     +nan       nan     0
  2210_2021_MN    NVA  3042   887569     -0.7      11.5  1293
  3056_2021_MN    VVA  3041   231726  +2131.2     720.3    95
  3056_2021_MN    VVA  3042  1155870    +20.6      17.9  1474
  2024_2021_MN    NVA  3038  1032036    +35.1      39.4  1334
  2024_2021_MN    NVA  3039    63518     +nan       nan     0
  2036_2021_MN    NVA  3041   118146     +nan       nan     0
  2036_2021_MN    NVA  3042  1054690    +59.4      27.4  1157
  2099_2021_MN    NVA  3032   820208    -43.7       8.7   927
  3089_2021_MN    VVA  3030    12694     +nan       nan     0
  3089_2021_MN    VVA  3031   530208   -109.2      36.9   803
```

**The measurement cannot be made at these marks.** At five of the six, the second line in
the 400 m box puts **zero** returns inside the fitting radius – the marks sit under one
line's swath, not in an overlap. The sixth, 3056, has 95 second-line returns at the mark
and returns +2131 mm with σ 720 mm; that is a swath edge, not a line offset, and it must
not be quoted as one.

A per-line statement about gen2 needs marks chosen for sitting in overlap, and none of
the six was.

---

## 8. Coverage, and what closing the gap would cost

```
  gen2 point cloud on disk covers Elba/elbaext only. Marks inside it:
    1079_2021_MN    LCP  role=calibration  in data/after/elbaext_3dep_fd_class2.laz
  plus 6 marks with a separately fetched 400 m box under data/after/checkpoints/
  390 marks carry a USGS per-point residual; 384 of them are NOT reachable from any point cloud on disk.
```

The 11 GB under `data/after/` covers Elba and elbaext and nothing else. Reaching the
other 384 marks means a 3DEP tile at 0.7–2.3 GB for each one, i.e. **hundreds of
gigabytes and hundreds of fetches even where marks share a tile**, to produce a second
opinion on a quantity USGS has already
measured and published at every one of them – with, per §6, a 55 mm sd between the two
opinions. **Recommendation: do not fetch.** The published field is the better instrument
here, and it was free.

---

## 9. What this does and does not settle for the DoD

**Settled.** gen2's level is no longer the weak side. On open ground, project-wide,
gen2 sits **−2.37 ± 2.37 mm** against its own held-out control, on 227 marks, with the
geoid cancelling exactly and no chain, no lateral shift and no conversion in the way.

**Not settled, in order of size:**

1. **Which cover stratum defines the datum.** 72 mm on the gen2 side, 71 mm on the gen1
   side, and it is the same question in both. Neither epoch's control answers it; it is a
   decision about what "the ground" means, and it is Andy's.
2. **The bridge from the published residual to our own surface.** gen1 has one, measured:
   **−7.2 ± 10.8 mm** over 18 marks (`CONTROL_RESIDUAL_FIELD.md` §5). gen2's is
   **+5.3 ± 22.6 mm** over 6 (§6 above), which is too thin. The two epochs' bridges are
   the least comparable part of the whole chain.
3. **Both epochs carry an unpublished vendor bias adjustment.** gen1's is in a paper-only
   accuracy report; gen2's is stated in the Lidar Mapping Report p. 15 without a value.
   Neither is recoverable from the data. Every level above is the level *after*
   adjustment, which is the right thing for a DoD but leaves the epoch difference
   resting on two constants nobody has published.
4. **The QL1/QL0 NVA−VVA gap** (92 mm against 41 mm), which the collection dates in the
   bundled CSV could probably explain.

**The arithmetic the DoD wants, stated plainly and with its caveat:**

```
    open / NVA         gen1  +13.79 - gen2   -3.38 =  +17.16 +/- 7.50 mm  (n 230/227)
    vegetated / VVA    gen1  -86.51 - gen2  -74.15 =  -12.36 +/- 8.90 mm  (n 534/163)
```

Delivered vendor DEM against delivered vendor DEM, which is the only like-for-like
pairing the two epochs offer – gen1 publishes no point-cloud residual. **These are
project-wide means over two differently distributed mark sets, not the value at Elba**,
and the gen1 open/vegetated grouping (`L1O` open; `L2T`+`L3B`+`L4F` vegetated) is a
choice I made to put the two taxonomies on one contrast. It is a proposal, and it
changes the vegetated row.

---

## 10. Reproducing

```
python analysis/groundtruth/parse_gen2_control.py \
  --report-pdf <MN_SE_Driftless_2021_B21_Ground_Control_Survey_Report.pdf> \
  --usgs-dir <dir with USGS_MN_SE_Driftless_2021_B21_QL{0,1}.dbf> \
  --contractor-dir <dir with MN_Driftless_NVA_VVA_UTM15_QL{0,1}.dbf> \
  --check --tol-m 1e-6 \
  --out src/lidar_diff_icp/groundtruth/data/mn_se_driftless_2021_control.csv

TRUST_GIT_REV=$(git rev-parse --short HEAD) env -u PROJ_DATA -u GDAL_DATA \
  ./lidar-icp/bin/python analysis/gen2_absolute_datum.py \
  --site-name Elba --site-easting 578762.8 --site-northing 4884487.6 \
  --sign-tol-m 1e-9 --surfaces laz,dem --block-preference QL1,QL0 \
  --band-radii-km 5,10,15,20,30,50,100,200 \
  --max-lag-m 20000,40000,60000,100000,184000 --n-lags 15,30 --n-pairs 200000,800000 \
  --estimators dowd,matheron --seeds 0,1,2 \
  --boxes-dir data/after/checkpoints --line-half-width-m 200 --res-m 5.0 \
  --gen1-csv src/lidar_diff_icp/groundtruth/data/mn_dnr_2008_control_semn.csv \
  --gen1-open-classes L1O --gen1-veg-classes L2T,L3B,L4F \
  --tie-json data/derived/groundtruth/gen2_checkpoint_tie.json
```

The `--tie-json` input is the output of
`analysis/groundtruth/gen2_checkpoint_tie.py --json`, run unchanged.

---

## 11. Verification status of every claim

**Verified this session, by a command whose output is pasted above:**

| claim | how |
|---|---|
| the USGS VATool shapefiles carry per-point `Z`, `DEMz`, `zdiff`, `LAZz`, `LAZzdiff` | both `.dbf` files downloaded and their schemas dumped |
| the residual is `surveyed − surface` | both subtractions tested on all 395 rows, 1e-9 m tolerance |
| the parse reproduces all eight published aggregates | driver §2, against the VATool's own `_VA.txt` and `_las_checkpoint_report.txt` |
| report tables hold 533 points, and `2198A_2022_MN` is the one a naive regex drops | id tokens counted in `pdftotext -layout` output |
| `3000_2021_MN` is in the shapefile and not in the report | set difference, printed by the parser |
| report and shapefile agree to 0.0000 m on 389 shared marks | parser `--check` |
| the geoid is GEOID18 on all 395 shapefile rows | `geoid` attribute, all rows; the driver raises otherwise |
| every level, contrast, band and kriged value in §3–§5, §7 | driver output, ledger `.trust/runs/20260827T072813-1731650.json` |
| six marks reachable on disk, 384 not | box bounds read from the LAZ headers |

**Quoted from a document, not re-derived here:**

| claim | source |
|---|---|
| the 143 LCPs calibrated gen2; NVA/VVA were held out | vendor FGDC `MN_SEDriftless_2_2021_Classified_Point_Cloud_Metadata.xml`, Ground Conditions, verbatim in `LIDAR_DOCUMENTATION_MINE.md` §1.2 |
| gen2 carries an unpublished vendor bias adjustment | Lidar Mapping Report p. 15, verbatim, same source |
| gen2's lines over Elba run east–west at 90.0–90.8° | vendor swath polygons, `LIDAR_DOCUMENTATION_MINE.md` §2.4 |
| gen2 was flown 2021-05-01 at green-up, NDVI 0.49 | the leaf-state work; not re-checked here |
| gen1's cover treatment is worth 71 mm at Elba | `CONTROL_RESIDUAL_FIELD.md`, its headline table |
| gen1's surface-to-surface bridge, −7.2 ± 10.8 mm | `CONTROL_RESIDUAL_FIELD.md` §5 |

**Explicitly a choice of mine, not a result:**

| choice | what it decides |
|---|---|
| gen1 open = `L1O`; vegetated = `L2T`+`L3B`+`L4F` (§9) | which gen1 marks are called vegetated, hence the vegetated row of the epoch difference. `L5U` urban is in neither group. |
| marks tested in both blocks taken from QL1 | 5 of 390 marks; the two blocks' DEMs differ at them by up to 3 cm |
| the 200 m half-width and `res` 5 m in §7 | how much of each box the per-line read sees |

**Not attempted:**

- The **QL0** block's own site value. Elba is in QL1; QL0 enters only as a contrast.
- Any use of the 143 LCPs as an accuracy check. They calibrated the data.
- Any 3DEP tile fetch. §8 says why.

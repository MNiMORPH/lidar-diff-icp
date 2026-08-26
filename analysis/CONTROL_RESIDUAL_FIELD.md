# The 2008 control residual as a spatial FIELD, and its value at Elba

**Date:** 2026-08-26
**Module:** `src/lidar_diff_icp/groundtruth/residual_field.py`
**Tests:** `tests/test_groundtruth_residual_field.py` (21)
**Driver:** `analysis/control_residual_field.py` — one run prints every table below
**Data:** `src/lidar_diff_icp/groundtruth/data/mn_dnr_2008_control_semn.csv`
(1 004 rows → **963 marks**), plus two committed per-mark tables read but not recomputed:
`analysis/GEN1_DATUM_MORE_MARKS.md` §5 (18 marks) and `analysis/GEN1_OWN_CONTROL_TIE.md`
§3 (16 marks)
**Figure:** `figures/control_residual_field.png`
**Bears on:** `analysis/GEN1_OWN_CONTROL_TIE.md` §6 (the plane fits this replaces),
`analysis/GEN1_DATUM_MODULE.md` §1 (the `sd/sqrt(n)` estimates this is an alternative to),
`analysis/FRAME_2026-08-26-PM.md` (the **+53.6 ± 13.0 mm** evening figure)

**No point cloud was read and nothing was downloaded.** This is a 963-point tabular
problem end to end.

Sign convention throughout, unchanged from `groundtruth.tie`: **positive = the surface
reads LOW**, i.e. the number is the constant to ADD. The vendor's `dnr_error_m` is
`Control Z − Surface Z`, re-derived below as exact on 1 004 of 1 004 rows, so it is the
same sign family as our `tie = surveyed − z_lidar` with no flip.

Every number below is pasted from the driver's output. The two places where I could not
reproduce a figure from the brief are marked `UNVERIFIED` in §12 with the source named.

---

## The number asked for

Kriged value of the 2008 control residual **at Elba's tile centre (E 578 762.8,
N 4 884 487.6, EPSG:26915)**, one row per cover treatment, summarised over the whole
10-row variogram sweep for that treatment so that no single `max_lag` or estimator is
privileged. `carried18` adds the **−7.2 ± 10.8 mm** surface-to-surface offset of §9 to
put the number on *our* reconstructed gen1 surface rather than the vendor's delivered
one. Positive = the surface reads LOW, i.e. the constant to ADD.

```
variant          n rows  pred med    pred min..max  sdfield med  carried18 med   carried18 min..max  sd18 med
open                 10     -13.1    -32.3..   -1.8         33.7          -20.4        -39.5..   -9.1      35.4
open+urban           10     +19.4     -0.6..  +28.4         32.5          +12.1         -7.8..  +21.1      34.3
cover-covariate      10     +58.0    +35.5..  +62.0         32.2          +50.7        +28.2..  +54.8      34.0
```

*(This block is the only table in the document that the driver does not print itself: it
is the per-variant median and min–max of the ten rows of §11, which the driver does print
in full. Every one of its numbers can be checked against §11 by eye.)*

**`sdfield med` is the uncertainty of one thing and one thing only:** the error in
predicting *the spatially correlated component of the vendor's control residual at that
coordinate*, given 963 marks and the fitted spherical variogram, with the nugget treated
as uncorrelated noise and filtered out. It is not an SE of a mean over marks, and it is
not the spread a single new control mark at Elba would show (that is `sd_mark_mm`,
47.5–104.9 mm).

**The choice among the three rows is worth 71 mm in the median and I have not made it.** §7 and §8
measure which treatment predicts an open mark better, on a common set of 230 open marks,
and the cover-covariate treatment wins at every block size — that is evidence for the
bottom row, and it is Andy's call whether to take it.

---

## Bottom line

1. **The answer at Elba is set by the cover treatment, not by the kriging.** Kriged
   vendor residual at Elba, across the whole variogram sweep:

   | treatment | marks | prediction range over the sweep | kriging sd of the field |
   |---|---|---|---|
   | open (`L1O`) | 230 | **−1.8 to −32.3 mm** | 29.3–71.8 mm |
   | open+urban (`L1O`+`L5U`) | 397 | **−0.6 to +28.4 mm** | 28.8–53.0 mm |
   | cover as a covariate (all 963) | 963 | **+35.5 to +62.0 mm** | 23.6–37.0 mm |

   The **94 mm** span between the treatments is larger than any single treatment's
   prediction sd. Choosing among them is the whole estimate; the kriging is the small
   part. I have not chosen one — §5 reports all three side by side and §7–§8 measure
   which of them predicts an open mark better.

2. **The field model does beat a single global constant, and the margin is largest for
   the cover-covariate treatment.** Scored on the *same* 230 open marks in every case
   (RMSE over three different mark sets is not a comparison), leave-one-out skill
   `1 − rmse_krige/rmse_null` is **0.110–0.276 (open)**, **0.167–0.320 (open+urban)**,
   **0.369–0.400 (cover-covariate)**.

3. **Under spatially blocked cross-validation that skill decays with block size and is
   effectively gone at 80 km.** On the same 230 marks, cover-covariate skill runs
   **0.261–0.304** at 15 km blocks, **0.184–0.247** at 25 km, **0.138–0.260** at 50 km
   and **0.013–0.095** at 80 km; the open treatment reaches **0.026–0.207** at 15 km and
   **0.004–0.066** at 80 km, and open+urban goes slightly **negative** (−0.014) in one
   80 km cell. Local prediction is worth something at tens of kilometres and worth
   nothing at the scale of the whole acquisition.

4. **The variogram never reaches a sill inside the data, so "the range" is not a
   number this control set determines.** The fitted range sits *at* the upper bound of
   its own fitting window at `max_lag` = 40 km (39 000 m) and 100 km (97 500 m), and the
   empirical semivariance is still climbing at 100 km before the far bins turn over on
   ~10³-pair counts. The whole sweep is reported and the empirical points are shown, in
   §5b and in the figure, rather than one fitted range.

5. **Carried to our own reconstructed surface, the offset is small and it is not zero.**
   Our tie minus the vendor residual, at the marks where both have been read:
   **−7.2 ± 10.8 mm** over the 18-mark set (r = +0.807) and **−16.4 ± 23.9 mm** over the
   16-mark set (r = +0.701). Negative means our reconstruction sits *above* the delivered
   2008 surface. Carrying the 18-mark offset, the cover-covariate prediction at Elba
   becomes **+28.2 to +54.8 mm**, the open+urban prediction **−7.8 to +21.1 mm**, and the
   open prediction **−39.5 to −9.1 mm**.

6. **The marks the recent estimates rest on sit well above their own local field, and how
   far above depends on which field.** The 18 open/urban marks on lines 133–138 average
   **+54.6 mm** of vendor residual against a population mean of **−18.7 mm** within 10 km
   (excess **+73.2 mm**) and **−69.5 mm** within 30 km (excess **+124.1 mm**). Against a
   leave-one-out kriged field evaluated at those same locations, the excess is
   **+31.7 mm** under open+urban and **+13.1 mm** under cover-covariate — both over all
   18 marks. The open treatment can only hold the **8** of them that are `L1O`, and over
   those 8 the excess is **+55.4 mm**; that is a different mark set and is not a
   like-for-like third number.

7. **County structure is real and is not a parameter.** On open cover, one-way ANOVA over
   8 counties gives **F = 17.13, p = 4.704e-18**, sd of the county means **59.0 mm** —
   against a pooled within-county sd of **88.5 mm**, which is larger. County is a
   diagnostic that a field exists, not a level to fit.

8. **What the uncertainty is the uncertainty of.** `sd_field_mm` is the standard
   deviation of the error in predicting *the spatially correlated component of the vendor
   residual at Elba's tile centre*, from the 963 marks under the stated variogram, with
   the nugget treated as uncorrelated noise and filtered out. It is **not** an SE of a
   mean over marks, and it is **not** the sd of what a single new control mark at Elba
   would read — that second quantity is reported separately as `sd_mark_mm` and is
   47.5–104.9 mm, larger by exactly the nugget.

---

## 1. What the file is, re-derived

```
                                                   what      value
  -----------------------------------------------------  ---------
                                            rows in CSV       1004
              unique marks after (E,N,Z) de-duplication        963
                                rows dropped as repeats         41
                         groups those repeats came from         39
  rows where elevation - dnr_surface_z_m == dnr_error_m  1004/1004
                     max |residual| of that identity, m  5.411e-14
  rows where dnr_surface_z_m - elevation == dnr_error_m     0/1004
                 max |residual| of the reverse order, m   1.080000
```

The de-duplication is on the exact `(easting, northing, elevation)` triple: a mark on a
county line is printed in both counties' validation reports. No tolerance is used,
because a tolerance would be a parameter nobody asked for.

## 2. The sample mean moves further than its own SE

```
  radius_km  n_all  mean_all_mm  se_all_mm  n_ou  mean_ou_mm  se_ou_mm
  ---------  -----  -----------  ---------  ----  ----------  --------
          5     11        -15.4       27.5     3        55.0      46.3
         10     34        -18.7       12.9     9        -3.1      26.0
         15     66        -31.2       11.5    19        -3.0      20.3
         20     99        -54.2       10.4    34       -37.2      16.2
         30    207        -69.5        7.6    66       -36.9      11.5
         50    431        -85.7        5.7   167       -47.3       7.3
        100    912        -48.0        4.4   372         6.0       5.7
        200    963        -43.4        4.3   397        10.3       5.5
```

`se_all_mm` and `se_ou_mm` are `sd/sqrt(n)` of the sample mean over the marks inside the
radius. They are the quantity this document exists to replace: between 10 and 50 km the
all-cover mean moves from **−18.7 to −85.7 mm** while its SE only ever reaches 12.9 mm.

**One discrepancy against the brief, and its cause.** The brief's open+urban column
(n = 18 at 15 km, mean +1.8, SE 20.8; n = 32 at 20 km, −36.6 ± 17.1; n = 60 at 30 km,
−37.2 ± 12.6) is reproduced **exactly** by classifying cover from the leading three
characters of `point_id` instead of from the `point_type` column. Twenty-two marks are
printed in the validation reports as `L10…` — digit zero — rather than `L1O`, so a prefix
match silently drops them from open terrain. `point_type` already carries them as `L1O`
(241 rows / 230 marks vs 219 rows / 209 marks by prefix), and that is the column used
everywhere here. **The same bug explains the brief's county ANOVA** — see §12, where its
`F = 15.71, p = 2.202e-16` reproduces exactly on the 209 prefix-matched marks.

## 3. By cover class, over all 963 marks

```
  cover                  meaning    n  mean_mm  median_mm  sd_mm  se_mm
  -----  -----------------------  ---  -------  ---------  -----  -----
    L1O             open terrain  230     13.8       26.0  108.1    7.1
    L2T     tall weeds and crops  214    -94.8      -79.5  134.5    9.2
    L3B      brush and low trees  164   -115.2     -112.0  130.1   10.2
    L4F                 forested  156    -45.0      -22.0  131.2   10.5
    L5U                    urban  167      5.5       18.0  110.7    8.6
  other  unclassed (FSA targets)   32      9.7        1.0   87.5   15.5
```

Open terrain is the only class whose mean is positive, and brush/low trees the most
negative. That ordering is what `Control − Surface` requires if vegetation left in the
ground class lifts the surface, and it is the reason a pooled field mixes a sensor effect
with a spatial one.

## 4. County as a diagnostic

```
          variant  n_counties    n      F          p  sd_group_means_mm  sd_within_mm
  ---------------  ----------  ---  -----  ---------  -----------------  ------------
             open           8  230  17.13  4.704e-18               59.0          88.5
       open+urban           8  397  27.26  2.008e-30               60.3          90.2
  cover-covariate           8  963  26.97  6.713e-34               53.4         121.5
```

## 5. The variogram sweep and the site prediction

`max_lag_m` is a row axis because the fitted range tracks it. `n_lags` (15, 30),
`n_pairs` (200 000, 800 000) and `seed` (0, 1, 2) are a 12-point nuisance grid: the
fitted parameters and the prediction are reported as the **median** over that grid, with
`pred_lo_mm`/`pred_hi_mm` giving its full min–max, so no single binning is privileged.

```
          variant  estimator  max_lag_m  nugget_mm2  sill_mm2  range_m  sqrt_tot_mm  pred_mm  pred_lo_mm  pred_hi_mm  sd_field_mm  sd_mark_mm
  ---------------  ---------  ---------  ----------  --------  -------  -----------  -------  ----------  ----------  -----------  ----------
             open       dowd      20000        1404      3897     8157         72.8     -1.8        -9.2         1.4         50.5        65.0
             open       dowd      40000        2012      6240    39000         90.8    -10.0       -10.7        -8.9         34.5        56.6
             open       dowd      60000        2252      7672    58000         99.6    -11.2       -12.0        -9.8         32.9        57.8
             open       dowd     100000        3374      8824    97500        110.4    -21.6       -22.9       -20.8         31.1        65.9
             open       dowd     184000        1397     18791   160427        142.1    -10.4       -12.5        -9.8         29.3        47.5
             open   matheron      20000         163      6921     5239         84.2     -2.1        -2.5         1.6         71.8        73.6
             open   matheron      40000        4423      5438    39000         99.3    -15.1       -15.7       -14.5         37.8        76.5
             open   matheron      60000        4937      5404    50552        101.7    -18.6       -21.8       -16.7         35.4        78.8
             open   matheron     100000        5758      6834    97500        112.2    -32.3       -32.8       -31.7         31.5        82.2
             open   matheron     184000        5289     10209   135078        124.5    -28.3       -28.8       -27.7         31.8        79.4
       open+urban       dowd      20000        2224      3468     5544         75.4     19.4        16.9        23.8         53.0        70.9
       open+urban       dowd      40000        2516      6492    39000         94.9     26.8        25.1        27.2         33.7        60.5
       open+urban       dowd      60000        3441      6626    57125        100.3     18.4        15.6        25.6         31.4        66.6
       open+urban       dowd     100000        3900      8917    97500        113.2     11.8        10.7        13.1         29.7        69.2
       open+urban       dowd     184000        2634     18532   179400        145.5     19.5        19.1        20.6         28.8        58.9
       open+urban   matheron      20000        4093      3844    19333         89.1     28.4        25.5        31.0         39.5        75.1
       open+urban   matheron      40000        4306      5876    39000        100.9     21.2        20.4        21.6         35.5        74.7
       open+urban   matheron      60000        4555      5651    41008        101.0     20.1        19.5        20.5         34.8        75.9
       open+urban   matheron     100000        5968      6684    97500        112.5     -0.6        -1.4         0.2         29.1        82.6
       open+urban   matheron     184000        5224     11805   155918        130.5      4.7         4.2         5.2         29.4        78.0
  cover-covariate       dowd      20000        5982      4724    19333        103.5     62.0        58.3        64.8         35.4        85.1
  cover-covariate       dowd      40000        6226      7174    39000        115.8     58.6        57.7        59.5         32.3        85.2
  cover-covariate       dowd      60000        6338      6904    39286        115.1     58.4        57.6        59.8         32.1        85.9
  cover-covariate       dowd     100000        8344      6005    66169        119.8     45.3        39.9        58.0         27.0        95.3
  cover-covariate       dowd     184000        8386     12151   179400        143.3     43.7        42.3        44.2         24.8        94.9
  cover-covariate   matheron      20000        7018      4982    19500        109.5     61.7        60.1        63.5         37.0        91.6
  cover-covariate   matheron      40000        7783      6358    36843        118.9     58.1        55.5        58.5         32.9        94.2
  cover-covariate   matheron      60000        7676      6218    34050        117.9     57.8        56.8        58.5         33.4        93.8
  cover-covariate   matheron     100000        8047      6408    40926        120.2     56.0        55.4        56.6         32.0        95.2
  cover-covariate   matheron     184000       10450      8215   164683        136.6     35.5        34.3        36.8         23.6       104.9
```

Three things to read off this.

* **The nuisance grid is not what moves the answer.** `pred_hi − pred_lo` is at most
  **18.1 mm** (cover-covariate, dowd, 100 km), with a median of **2.1 mm** over the 30
  rows and 22 of 30 rows under 4 mm. The estimator
  choice moves it more (open at 100 km: −21.6 dowd vs −32.3 matheron), and `max_lag`
  more again.
* **The range is pinned at its own fitting bound in nine of the twelve `max_lag` = 40 km
  and 100 km rows** — 39 000 m and 97 500 m are the largest lag centre,
  `max_lag × (1 − 0.5/n_lags)`, for the two `n_lags` values. The three that are not
  pinned are all `cover-covariate` (dowd at 100 km, 66 169 m; matheron at 40 km,
  36 843 m, and at 100 km, 40 926 m), which is the only variant whose variogram is
  computed on a detrended residual. A spherical model fitted to a variogram that has not
  reached a sill returns the window, not the field.
* **The nugget is the honest floor, and it is not small.** `sqrt(nugget)` runs from
  **12.8 mm** (open, matheron, 20 km, fitted nugget 163 mm²) to **102.2 mm**
  (cover-covariate, matheron, 184 km, 10 450 mm²), with the other 28 rows spanning
  37.4 to 91.6 mm. No amount of control at this density predicts what a *single* mark will
  read better than that.

### 5b. The empirical variogram itself

The full sweep — every variant, estimator, `max_lag`, `n_lags`, `n_pairs` and seed, one
row per lag bin — is written to
`data/derived/control_residual_field/empirical_variogram_sweep.csv` (gitignored,
regenerable). The figure draws all of it. Printed below is one slice of it, the open
variant at the widest swept `max_lag`, so the shape is on the page and not only in the
model:

```
          variant  max_lag_m   lag_m  gamma_mm2  pairs
  ---------------  ---------  ------  ---------  -----
             open     184000    6133       4362   8370
             open     184000   18400       5082  18888
             open     184000   30667       7940  22798
             open     184000   42933       9505  21870
             open     184000   55200      10128  22356
             open     184000   67467      10340  21196
             open     184000   79733      11659  20673
             open     184000   92000      12348  17993
             open     184000  104267      15044  14745
             open     184000  116533      29920  10549
             open     184000  128800      32513   8004
             open     184000  141067      16898   5477
             open     184000  153333       9302   4096
             open     184000  165600       7390   1786
             open     184000  177867       3571    289
```

Semivariance climbs monotonically across all eleven bins out to 128.8 km and then falls
away over the last four. Those four falling bins carry **289 to 5 477** pairs, against
**17 993 to 22 798** in the middle of the range. *(Reading,
offered as hypothesis and not tested here: at those pair counts the far-field bins are the
sparse corners of the point set rather than a resolved hole effect. The `--n-pairs` sweep
does not settle it, because the shortage is of mark pairs, not of sampled pairs.)*

### 5c. Under the cover-covariate treatment, every class predicted at Elba

```
  cover                  meaning  max_lag_m  pred_mm  pred_lo_mm  pred_hi_mm  sd_field_mm
  -----  -----------------------  ---------  -------  ----------  ----------  -----------
    L1O             open terrain      20000     62.0        58.3        64.8         35.4
    L2T     tall weeds and crops      20000    -34.3       -37.9       -31.4         35.2
    L3B      brush and low trees      20000    -54.5       -58.5       -51.9         35.4
    L4F                 forested      20000     19.9        16.6        23.3         35.6
    L5U                    urban      20000     78.0        74.0        81.2         35.6
  other  unclassed (FSA targets)      20000     72.1        67.5        75.5         38.4
    L1O             open terrain     184000     43.7        42.3        44.2         24.8
    L2T     tall weeds and crops     184000    -48.1       -49.3       -47.6         24.5
    L3B      brush and low trees     184000    -74.0       -75.5       -73.3         24.9
    L4F                 forested     184000      1.1        -0.4         1.7         25.1
    L5U                    urban     184000     56.4        54.8        57.1         25.2
  other  unclassed (FSA targets)     184000     48.6        46.9        49.3         29.6
```

The drift term is a **per-class constant estimated over the whole project**; the spatial
part is then common to every class, which is why `sd_field_mm` is the same down each
block. The open-class prediction here (+62.0 at 20 km, +43.7 at 184 km) is what the whole
963-mark field says about open ground at Elba; the open-only treatment, using the 230
open marks alone, says −1.8 and −10.4 at the same two settings. §7 and §8 score both out
of sample. *(Hypothesis, not tested here: a per-class offset that varies regionally would
break this treatment in exactly the direction of that disagreement. A county × cover
interaction test would separate it.)*

## 6. The cross-validation arithmetic, checked

```
          variant  n  max_abs_err_mm  max_abs_var_mm2
  ---------------  -  --------------  ---------------
             open  6        1.31e-12         4.55e-12
       open+urban  6        8.29e-13         7.28e-12
  cover-covariate  6        2.75e-12         2.00e-11
```

Leave-one-out below is computed from a single inverse of the augmented kriging system
rather than 963 refits. It is checked against a genuine refit at six marks per variant,
and by `tests/test_groundtruth_residual_field.py` at every mark of a synthetic set, for
both ordinary kriging and a cover drift. Both tests were shown to FAIL when the identity's
sign is flipped.

## 7. Leave-one-out against the null of one global constant

```
          variant  estimator  max_lag_m    n  rmse_krige_mm  rmse_null_mm  mae_krige_mm  mae_null_mm  skill  n_L1O  rmse_L1O_mm  null_L1O_mm  skill_L1O
  ---------------  ---------  ---------  ---  -------------  ------------  ------------  -----------  -----  -----  -----------  -----------  ---------
             open       dowd      20000  230           91.4         108.4          68.5         82.9  0.156    230         91.4        108.4      0.156
             open       dowd      40000  230           78.5         108.4          54.6         82.9  0.276    230         78.5        108.4      0.276
             open       dowd      60000  230           79.4         108.4          55.6         82.9  0.267    230         79.4        108.4      0.267
             open       dowd     100000  230           79.8         108.4          56.0         82.9  0.263    230         79.8        108.4      0.263
             open       dowd     184000  230           78.9         108.4          55.0         82.9  0.272    230         78.9        108.4      0.272
             open   matheron      20000  230           96.5         108.4          72.7         82.9  0.110    230         96.5        108.4      0.110
             open   matheron      40000  230           79.0         108.4          55.2         82.9  0.271    230         79.0        108.4      0.271
             open   matheron      60000  230           79.8         108.4          56.1         82.9  0.264    230         79.8        108.4      0.264
             open   matheron     100000  230           81.2         108.4          57.6         82.9  0.251    230         81.2        108.4      0.251
             open   matheron     184000  230           80.9         108.4          57.3         82.9  0.253    230         80.9        108.4      0.253
       open+urban       dowd      20000  397           86.6         109.3          65.3         84.4  0.207    230         90.2        108.2      0.167
       open+urban       dowd      40000  397           71.8         109.3          53.9         84.4  0.343    230         73.6        108.2      0.320
       open+urban       dowd      60000  397           73.3         109.3          55.2         84.4  0.329    230         75.0        108.2      0.307
       open+urban       dowd     100000  397           74.1         109.3          55.6         84.4  0.322    230         75.7        108.2      0.300
       open+urban       dowd     184000  397           72.8         109.3          54.7         84.4  0.334    230         74.8        108.2      0.309
       open+urban   matheron      20000  397           72.9         109.3          54.9         84.4  0.333    230         75.0        108.2      0.307
       open+urban   matheron      40000  397           73.1         109.3          55.0         84.4  0.332    230         74.4        108.2      0.313
       open+urban   matheron      60000  397           73.4         109.3          55.3         84.4  0.329    230         74.6        108.2      0.311
       open+urban   matheron     100000  397           76.6         109.3          57.5         84.4  0.299    230         77.6        108.2      0.283
       open+urban   matheron     184000  397           75.9         109.3          57.1         84.4  0.305    230         77.2        108.2      0.287
  cover-covariate       dowd      20000  963           92.8         132.6          70.4        104.7  0.300    230         74.5        122.3      0.391
  cover-covariate       dowd      40000  963           93.4         132.6          70.8        104.7  0.296    230         73.4        122.3      0.400
  cover-covariate       dowd      60000  963           93.4         132.6          70.8        104.7  0.295    230         73.4        122.3      0.399
  cover-covariate       dowd     100000  963           95.1         132.6          72.0        104.7  0.282    230         75.1        122.3      0.386
  cover-covariate       dowd     184000  963           95.8         132.6          72.6        104.7  0.277    230         75.6        122.3      0.382
  cover-covariate   matheron      20000  963           92.8         132.6          70.4        104.7  0.300    230         74.4        122.3      0.391
  cover-covariate   matheron      40000  963           93.8         132.6          71.0        104.7  0.293    230         73.7        122.3      0.397
  cover-covariate   matheron      60000  963           93.7         132.6          70.9        104.7  0.293    230         74.0        122.3      0.395
  cover-covariate   matheron     100000  963           93.9         132.6          71.2        104.7  0.292    230         73.6        122.3      0.398
  cover-covariate   matheron     184000  963           97.6         132.6          74.3        104.7  0.263    230         77.1        122.3      0.369

  The last four columns score all three variants on THE SAME marks (every
  L1O mark each variant contains), because rmse over different mark sets is not
  a comparison. The open variant's restricted and unrestricted columns are the
  same numbers by construction.
```

## 8. Spatially blocked cross-validation

The variogram is re-estimated inside every training fold, which is the honest version;
the block side is swept.

```
          variant  block_m  n_blocks  max_lag_m    n  rmse_krige_mm  rmse_null_mm  mae_krige_mm  mae_null_mm   skill  n_L1O  rmse_L1O_mm  null_L1O_mm  skill_L1O
  ---------------  -------  --------  ---------  ---  -------------  ------------  ------------  -----------  ------  -----  -----------  -----------  ---------
             open    15000        61      20000  230          106.8         109.7          81.7         82.9   0.026    230        106.8        109.7      0.026
             open    15000        61      40000  230           89.4         109.7          64.9         82.9   0.185    230         89.4        109.7      0.185
             open    15000        61      60000  230           87.7         109.7          64.1         82.9   0.200    230         87.7        109.7      0.200
             open    15000        61     100000  230           86.9         109.7          63.3         82.9   0.207    230         86.9        109.7      0.207
             open    15000        61     184000  230           88.3         109.7          64.6         82.9   0.195    230         88.3        109.7      0.195
             open    25000        24      20000  230          108.0         111.5          83.0         83.9   0.031    230        108.0        111.5      0.031
             open    25000        24      40000  230           94.0         111.5          69.6         83.9   0.157    230         94.0        111.5      0.157
             open    25000        24      60000  230           93.5         111.5          70.4         83.9   0.161    230         93.5        111.5      0.161
             open    25000        24     100000  230           94.4         111.5          71.2         83.9   0.154    230         94.4        111.5      0.154
             open    25000        24     184000  230           94.6         111.5          71.2         83.9   0.151    230         94.6        111.5      0.151
             open    50000         9      20000  230          114.3         114.5          88.5         86.6   0.002    230        114.3        114.5      0.002
             open    50000         9      40000  230          108.8         114.5          82.5         86.6   0.050    230        108.8        114.5      0.050
             open    50000         9      60000  230          108.0         114.5          82.1         86.6   0.057    230        108.0        114.5      0.057
             open    50000         9     100000  230           98.0         114.5          72.3         86.6   0.144    230         98.0        114.5      0.144
             open    50000         9     184000  230           96.0         114.5          69.4         86.6   0.162    230         96.0        114.5      0.162
             open    80000         5      20000  230          127.8         132.5         102.7        103.8   0.035    230        127.8        132.5      0.035
             open    80000         5      40000  230          132.0         132.5         107.7        103.8   0.004    230        132.0        132.5      0.004
             open    80000         5      60000  230          131.1         132.5         107.4        103.8   0.011    230        131.1        132.5      0.011
             open    80000         5     100000  230          128.2         132.5         106.7        103.8   0.032    230        128.2        132.5      0.032
             open    80000         5     184000  230          123.7         132.5         102.4        103.8   0.066    230        123.7        132.5      0.066
       open+urban    15000        61      20000  397          105.6         110.8          82.2         84.9   0.047    230        105.5        109.6      0.038
       open+urban    15000        61      40000  397           86.2         110.8          65.5         84.9   0.222    230         85.6        109.6      0.219
       open+urban    15000        61      60000  397           85.5         110.8          64.5         84.9   0.228    230         84.9        109.6      0.226
       open+urban    15000        61     100000  397           85.3         110.8          63.9         84.9   0.230    230         84.3        109.6      0.231
       open+urban    15000        61     184000  397           85.6         110.8          64.5         84.9   0.227    230         85.0        109.6      0.224
       open+urban    25000        24      20000  397          109.0         112.4          84.7         86.6   0.030    230        107.6        111.4      0.034
       open+urban    25000        24      40000  397           95.4         112.4          72.7         86.6   0.152    230         93.4        111.4      0.162
       open+urban    25000        24      60000  397           93.2         112.4          71.0         86.6   0.171    230         91.7        111.4      0.177
       open+urban    25000        24     100000  397           93.7         112.4          71.6         86.6   0.166    230         92.8        111.4      0.167
       open+urban    25000        24     184000  397           93.9         112.4          71.7         86.6   0.165    230         93.3        111.4      0.162
       open+urban    50000         9      20000  397          113.7         114.9          88.8         89.2   0.010    230        113.0        114.1      0.010
       open+urban    50000         9      40000  397          104.3         114.9          80.2         89.2   0.093    230        101.8        114.1      0.108
       open+urban    50000         9      60000  397          102.4         114.9          78.4         89.2   0.109    230         99.8        114.1      0.125
       open+urban    50000         9     100000  397           95.1         114.9          70.4         89.2   0.173    230         93.0        114.1      0.184
       open+urban    50000         9     184000  397           95.0         114.9          71.0         89.2   0.173    230         92.9        114.1      0.186
       open+urban    80000         5      20000  397          127.9         132.0         102.6        103.6   0.031    230        127.3        131.7      0.033
       open+urban    80000         5      40000  397          134.1         132.0         109.1        103.6  -0.016    230        133.5        131.7     -0.014
       open+urban    80000         5      60000  397          130.4         132.0         107.1        103.6   0.012    230        131.2        131.7      0.004
       open+urban    80000         5     100000  397          129.7         132.0         106.1        103.6   0.017    230        130.3        131.7      0.010
       open+urban    80000         5     184000  397          126.6         132.0         103.9        103.6   0.041    230        127.2        131.7      0.034
  cover-covariate    15000        65      20000  963          108.1         133.7          83.7        105.6   0.192    230         91.4        123.6      0.261
  cover-covariate    15000        65      40000  963          107.1         133.7          81.4        105.6   0.199    230         86.5        123.6      0.300
  cover-covariate    15000        65      60000  963          107.0         133.7          81.5        105.6   0.199    230         86.0        123.6      0.304
  cover-covariate    15000        65     100000  963          108.7         133.7          82.9        105.6   0.186    230         87.0        123.6      0.296
  cover-covariate    15000        65     184000  963          107.4         133.7          82.1        105.6   0.196    230         86.0        123.6      0.304
  cover-covariate    25000        25      20000  963          118.7         134.7          92.0        106.4   0.119    230        102.2        125.2      0.184
  cover-covariate    25000        25      40000  963          115.4         134.7          89.4        106.4   0.144    230         94.5        125.2      0.245
  cover-covariate    25000        25      60000  963          116.8         134.7          90.8        106.4   0.133    230         96.0        125.2      0.233
  cover-covariate    25000        25     100000  963          118.4         134.7          92.1        106.4   0.121    230         97.4        125.2      0.222
  cover-covariate    25000        25     184000  963          115.0         134.7          89.4        106.4   0.147    230         94.2        125.2      0.247
  cover-covariate    50000        10      20000  963          122.1         136.2          96.2        108.3   0.104    230        109.6        127.1      0.138
  cover-covariate    50000        10      40000  963          120.2         136.2          93.8        108.3   0.118    230        100.3        127.1      0.211
  cover-covariate    50000        10      60000  963          122.7         136.2          95.7        108.3   0.099    230        102.7        127.1      0.192
  cover-covariate    50000        10     100000  963          122.8         136.2          95.1        108.3   0.099    230        100.6        127.1      0.208
  cover-covariate    50000        10     184000  963          117.7         136.2          90.0        108.3   0.136    230         94.1        127.1      0.260
  cover-covariate    80000         5      20000  963          132.5         145.4         105.8        117.3   0.089    230        126.8        140.2      0.095
  cover-covariate    80000         5      40000  963          138.6         145.4         112.8        117.3   0.046    230        133.2        140.2      0.049
  cover-covariate    80000         5      60000  963          138.6         145.4         112.7        117.3   0.046    230        134.1        140.2      0.043
  cover-covariate    80000         5     100000  963          142.2         145.4         115.6        117.3   0.021    230        138.3        140.2      0.013
  cover-covariate    80000         5     184000  963          140.6         145.4         114.8        117.3   0.033    230        136.7        140.2      0.024
```

**This is the result that limits the whole exercise.** Every treatment's skill falls as
the held-out block grows, and at 80 km blocks — five folds, i.e. asking the field to
extrapolate across a county-sized gap — nothing beats the constant by more than 0.10, and
one cell is negative. The field is real at tens of kilometres; it is not a law that can
be carried across the acquisition.

## 9. Our surface against the delivered one

```
                             set   n  dmean_mm  dmedian_mm  dsd_mm  dse_mm      r
  ------------------------------  --  --------  ----------  ------  ------  -----
  18mk-lines133-138-class2-nolat  18      -7.2        -2.8    45.8    10.8  0.807
      16mk-tiles-on-disk-CSF-lat  16     -16.4       -13.7    95.7    23.9  0.701

  (our tie - vendor residual) = (surveyed - z_ours) - (surveyed - z_vendor)
  = z_vendor - z_ours.  POSITIVE means our reconstructed surface sits BELOW the
  delivered 2008 surface, so this is the constant to ADD to a field prediction
  made on the vendor residual before it applies to our surface.
```

`(our tie − vendor residual) = (surveyed − z_ours) − (surveyed − z_vendor) = z_vendor −
z_ours`. The two sets are two different reconstructions — the 18-mark set uses the
vendor's class-2 ground with no lateral shift, the 16-mark set our CSF ground with the
elbaext Nuth & Kääb shift — and they are reported separately for that reason, not
pooled. The 18-mark set is both larger and tighter (sd 45.8 vs 95.7 mm).

## 10. Where the selected marks sit relative to the field

```
          variant                            set2   n  sel_mean_mm  field_at_sel_mm  excess_mm
  ---------------  ------------------------------  --  -----------  ---------------  ---------
             open  18mk-lines133-138-class2-nolat   8         68.2             12.8       55.4
             open      16mk-tiles-on-disk-CSF-lat   4        -19.5              1.4      -20.9
       open+urban  18mk-lines133-138-class2-nolat  18         54.6             22.9       31.7
       open+urban      16mk-tiles-on-disk-CSF-lat   6         30.5             20.9        9.6
  cover-covariate  18mk-lines133-138-class2-nolat  18         54.6             41.5       13.1
  cover-covariate      16mk-tiles-on-disk-CSF-lat  16        -43.2            -48.6        5.4

  The field prediction at each selected mark is LEAVE-ONE-OUT, so the mark
  itself does not enter the field it is compared against.
```

**Read the `n` column before comparing rows.** A variant can only evaluate the marks it
contains: the open treatment holds 8 of the 18 and 4 of the 16, open+urban 18 and 6, the
cover-covariate treatment all of both. Only the open+urban and cover-covariate rows of
the 18-mark set are the same marks, and it is those two — **+31.7** and **+13.1 mm** —
that can be set against each other.

```
                             set   n  radius_km  sel_mean_mm  pop_mean_mm  excess_mm
  ------------------------------  --  ---------  -----------  -----------  ---------
  18mk-lines133-138-class2-nolat  18          5         54.6        -15.4       69.9
  18mk-lines133-138-class2-nolat  18         10         54.6        -18.7       73.2
  18mk-lines133-138-class2-nolat  18         15         54.6        -31.2       85.8
  18mk-lines133-138-class2-nolat  18         20         54.6        -54.2      108.7
  18mk-lines133-138-class2-nolat  18         30         54.6        -69.5      124.1
  18mk-lines133-138-class2-nolat  18         50         54.6        -85.7      140.2
  18mk-lines133-138-class2-nolat  18        100         54.6        -48.0      102.6
  18mk-lines133-138-class2-nolat  18        200         54.6        -43.4       98.0
      16mk-tiles-on-disk-CSF-lat  16          5        -43.2        -15.4      -27.8
      16mk-tiles-on-disk-CSF-lat  16         10        -43.2        -18.7      -24.5
      16mk-tiles-on-disk-CSF-lat  16         15        -43.2        -31.2      -12.0
      16mk-tiles-on-disk-CSF-lat  16         20        -43.2        -54.2       11.0
      16mk-tiles-on-disk-CSF-lat  16         30        -43.2        -69.5       26.4
      16mk-tiles-on-disk-CSF-lat  16         50        -43.2        -85.7       42.5
      16mk-tiles-on-disk-CSF-lat  16        100        -43.2        -48.0        4.8
      16mk-tiles-on-disk-CSF-lat  16        200        -43.2        -43.4        0.2
```

The point of the two tables together: the +54.6 mm mean of the 18 marks is **+73.2 mm**
above the plain population mean within 10 km, but only **+13.1 mm** above a
cover-covariate field evaluated at those same 18 locations with each mark left out.
*(Reading, offered as hypothesis: most of the apparent selection excess is the cover
composition of the neighbourhood rather than the siting of the marks, since the treatment
that carries a cover term removes most of it. What would test it is a field fitted on
open marks alone at a density that can resolve Elba's neighbourhood, which this control
set does not supply — the open treatment leaves **+55.4 mm** standing over its 8
available marks.)*

## 11. The prediction carried to our own reconstructed surface

```
          variant  estimator  max_lag_m  pred_mm  sd_field_mm  carried18_mm  sd18_mm  carried16_mm  sd16_mm
  ---------------  ---------  ---------  -------  -----------  ------------  -------  ------------  -------
             open       dowd      20000     -1.8         50.5          -9.1     51.6         -18.2     55.9
             open       dowd      40000    -10.0         34.5         -17.2     36.1         -26.4     42.0
             open       dowd      60000    -11.2         32.9         -18.4     34.6         -27.6     40.7
             open       dowd     100000    -21.6         31.1         -28.9     32.9         -38.0     39.2
             open       dowd     184000    -10.4         29.3         -17.6     31.2         -26.8     37.8
             open   matheron      20000     -2.1         71.8          -9.3     72.6         -18.5     75.7
             open   matheron      40000    -15.1         37.8         -22.4     39.3         -31.5     44.7
             open   matheron      60000    -18.6         35.4         -25.9     37.0         -35.0     42.7
             open   matheron     100000    -32.3         31.5         -39.5     33.3         -48.7     39.6
             open   matheron     184000    -28.3         31.8         -35.5     33.6         -44.7     39.8
       open+urban       dowd      20000     19.4         53.0          12.1     54.1           3.0     58.1
       open+urban       dowd      40000     26.8         33.7          19.6     35.4          10.4     41.3
       open+urban       dowd      60000     18.4         31.4          11.1     33.2           2.0     39.5
       open+urban       dowd     100000     11.8         29.7           4.6     31.6          -4.6     38.1
       open+urban       dowd     184000     19.5         28.8          12.2     30.8           3.1     37.4
       open+urban   matheron      20000     28.4         39.5          21.1     40.9          12.0     46.2
       open+urban   matheron      40000     21.2         35.5          13.9     37.1           4.8     42.8
       open+urban   matheron      60000     20.1         34.8          12.9     36.4           3.7     42.2
       open+urban   matheron     100000     -0.6         29.1          -7.8     31.0         -17.0     37.7
       open+urban   matheron     184000      4.7         29.4          -2.5     31.3         -11.7     37.9
  cover-covariate       dowd      20000     62.0         35.4          54.8     37.0          45.6     42.7
  cover-covariate       dowd      40000     58.6         32.3          51.4     34.1          42.2     40.2
  cover-covariate       dowd      60000     58.4         32.1          51.1     33.9          42.0     40.0
  cover-covariate       dowd     100000     45.3         27.0          38.0     29.1          28.9     36.1
  cover-covariate       dowd     184000     43.7         24.8          36.5     27.0          27.3     34.5
  cover-covariate   matheron      20000     61.7         37.0          54.5     38.5          45.3     44.1
  cover-covariate   matheron      40000     58.1         32.9          50.9     34.6          41.7     40.7
  cover-covariate   matheron      60000     57.8         33.4          50.5     35.1          41.4     41.1
  cover-covariate   matheron     100000     56.0         32.0          48.8     33.8          39.6     39.9
  cover-covariate   matheron     184000     35.5         23.6          28.2     25.9          19.1     33.6

  offsets carried: 18-mark set -7.2 +/- 10.8 mm, 16-mark set -16.4 +/- 23.9 mm (mean +/- sd/sqrt(n) of
  our tie minus the vendor residual at the marks where both have been read).
  The offset own uncertainty is an SE OF A MEAN over marks; only the field
  term is a prediction sd at the site, and the two are added in quadrature.
```

`sd18_mm` and `sd16_mm` add a prediction sd and an SE of a mean in quadrature. They are
different kinds of uncertainty and the sum is only as meaningful as that mixture; the two
terms are printed separately in §5 and §9 so the mixture can be undone.

## 12. What I could not reproduce, and what I chose

**Every number the brief supplied, checked.** Searched across `analysis/*.md` and
re-derived from the CSV where a construction could be guessed:

| brief's number | status |
|---|---|
| gen1 datum **+50.6 ± 30.5 mm** from 14 marks | `UNVERIFIED` — not found. `analysis/GEN1_DATUM_MODULE.md` §1 gives **+44.4 mm, SE 30.5, 14 marks, 6 lines** for `common_datum` in the elbaext frame; the SE matches and the value does not. `analysis/FRAME_2026-08-26-PM.md` gives **+53.6 ± 13.0 mm** over 18 marks |
| those 14 marks' vendor residuals averaging **+42.1 mm** | `UNVERIFIED` — the 14-mark set is not identified in any committed table. The 18-mark set of `GEN1_DATUM_MORE_MARKS.md` §5 averages **+54.6 mm**, computed here |
| ours − vendor at 14 marks: mean **−4.4**, median **−11.5**, sd **49.2**, r **+0.849** | `UNVERIFIED` — same reason. The two committed sets give **−7.2 / −2.8 / 45.8 / +0.807** (18 marks) and **−16.4 / −13.7 / 95.7 / +0.701** (16 marks) |
| ANOVA on open cover **F = 15.71, p = 2.202e-16**, sd of county means **58.4**, within **91.5** | **reproduced exactly** — on the de-duplicated 963 marks with open terrain taken from the `point_id` prefix (**209** marks). It is the same `L10`/`L1O` prefix bug as §2. Using the `point_type` column (**230** marks) the same statistics are **F = 17.13, p = 4.704e-18, 59.0, 88.5** |

**Parameters I chose, none of which drops an observation.** All are printed in the run
banner tagged `MINE`, and all are swept:

| what | value | what it does |
|---|---|---|
| `band_radii_km` | 5, 10, 15, 20, 30, 50, 100, 200 | rows of the §2 table only; selects nothing downstream |
| `max_lag_m` | 20, 40, 60, 100, 184 km | a row axis, because the fitted range tracks it |
| `n_lags` / `n_pairs` / `seeds` | (15, 30) / (200 k, 800 k) / (0, 1, 2) | nuisance grid, reported as median with full min–max |
| `block_m` | 15, 25, 50, 80 km | swept block sides for §8; no mark is dropped at any of them |
| `sign_tol_m` | 1e-9 | when the `Control − Surface` identity is called exact; the max residual is printed beside it |
| `loo_verify_idx` | 0, 1, 7, 42, 120, 229 | marks at which the LOO shortcut is checked against a refit |

There is **no** minimum `n`, no maximum slope, no spread or siting screen, no distance
cut and no search neighbourhood anywhere in the module or the driver.
`tests/test_groundtruth_residual_field.py::test_every_tunable_is_required` fails the
moment one of them acquires a default.

## 13. What this does not settle

* **Which cover treatment is right.** §7 and §8 say the cover-covariate field predicts an
  open mark better than the open-only field does, on the common set of 230 open marks and
  at every block size. That is evidence, not a decision, and the decision moves the answer
  at Elba by 71 mm in the median of the sweep. *(Hypothesis, untested here: the cover-covariate field wins because
  963 marks constrain the neighbourhood of Elba where 230 do not, and it would lose if the
  per-class offset varied regionally. A county × cover interaction test would separate
  those; it has not been run.)*
* **Whether the vendor's delivered surface and our reconstruction differ by a constant.**
  §9 estimates one offset per reconstruction, each an SE of a mean over 16–18 marks; the
  16-mark sd of 95.7 mm says the per-mark agreement is not good, as
  `GEN1_OWN_CONTROL_TIE.md` §5 already reported.
* **Anything at all about gen2.** No cross-epoch term enters here.

## 14. Reproducing

```
TRUST_GIT_REV=$(git rev-parse --short HEAD) ./lidar-icp/bin/python analysis/control_residual_field.py \
  --site-name Elba --site-easting 578762.8 --site-northing 4884487.6 \
  --sign-tol-m 1e-9 --band-radii-km 5,10,15,20,30,50,100,200 \
  --max-lag-m 20000,40000,60000,100000,184000 --n-lags 15,30 --n-pairs 200000,800000 \
  --estimators dowd,matheron --seeds 0,1,2 \
  --block-m 15000,25000,50000,80000 --loo-verify-idx 0,1,7,42,120,229 \
  --vgm-csv data/derived/control_residual_field/empirical_variogram_sweep.csv \
  --fig figures/control_residual_field.png \
  --tie-table-18 analysis/GEN1_DATUM_MORE_MARKS.md \
  --tie-table-16 analysis/GEN1_OWN_CONTROL_TIE.md
```

Runtime about twelve minutes, all of it in `numpy`'s dense solves; peak memory is a
969 × 969 matrix.

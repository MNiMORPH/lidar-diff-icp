# Enlarging the 2008 control on gen1's own flight lines 133-138

**Date:** 2026-08-26
**Scripts:** `analysis/groundtruth/gen1_line_tracks.py` (flight-line ground tracks),
`analysis/groundtruth/gen1_swath_seam.py` (where the vendor's ground class is cut),
`analysis/groundtruth/gen1_more_marks_tie.py` (every tie below),
`analysis/groundtruth/gen1_more_marks_report.py` (every aggregate below)
**Data:** `src/lidar_diff_icp/groundtruth/data/mn_dnr_2008_control_semn.csv`
(1 004 rows, 986 distinct `point_id`), plus 46 gen1 tiles on disk
**Follows:** `analysis/GEN1_OWN_CONTROL_TIE.md` (sign conventions, the 16-mark run),
`analysis/ADDITIONAL_GROUND_CONTROL.md` §7.1 (the siting screen)

Sign convention, unchanged: `tie = surveyed - z_lidar`, **positive = gen1 reads BELOW
the mark**. `geoid_shift_m = 0` is asserted, not approximated: the 2008 control is
NAVD88(GEOID03) and so is the raw gen1 cloud, so this is one frame and the geoid has
nothing to convert. `swath_shift_m = (0, 0, 0)`: no gen2-derived lateral or vertical
term enters a gen1-against-its-own-control comparison. Ground source is **vendor class
2**, not CSF.

---

## Bottom line

1. **The catchment cannot be widened. There is a hard physical boundary at 481 m and the
   old 500 m rule was already sitting on it.** The delivered 2008 tiles carry an explicit
   **overlap class (class 12)**: where two swaths see the same ground, one line's returns
   stay in the ground class and the other line's are moved out of it. Bare-earth ground is
   therefore cut at a **seam at half the line spacing**, not at the edge of the swath.
   Measured: swath half-width **718-792 m**, class-2 ground half-width **473-577 m**,
   line spacing **942/987/932/988/956 m, mean 961, half-spacing 481 m**.

2. **Going from 500 m to 2 km added 30 candidates and 0 marks.** Every one of the 12
   extra candidates carries ground from lines **131, 132 or 140** and none from 133-138.
   The furthest mark that does carry a 133-138 line sits at **479.3 m**, i.e. on the seam.
   Five tiles (107 008 578 bytes) were fetched to establish this and they bought nothing.

3. **The nearest-centreline proxy got 0 of the 18 wrong, and that is not luck.** The
   vendor's seam is the perpendicular bisector between adjacent flight lines, which is
   exactly the partition "nearest centreline" computes. Read off `point_source_id`,
   **every one of the 18 marks carries exactly one line** — so there are no mixed marks,
   no line-pairs at a mark, and the per-line ties are **bit-identical** to the pooled ties
   (max |difference| **0.000000 mm** over 18 marks).

4. **The estimate is therefore unchanged: +55.0 mm, SE over the six lines 16.6 mm**
   (marks as independent, **+47.3 ± 14.1 mm**, n = 18; median **+46.4 mm**). More marks
   did not tighten it because there are no more marks to have.

5. **The line structure claimed from the 43-mark screen does not survive inside
   133-138.** Over the six lines, ANOVA gives **F = 2.12, p = 0.133, df = (5, 12)** — not
   significant. The **F = 8.63** figure is reproduced exactly, and it is an ANOVA over
   **lines 128-149** with the singleton lines dropped (10 groups, 38 marks,
   p = 4.61e-06). Whatever is organised by line is organised over the 22-line span, not
   within these six.

6. **The larger sample gives the pipeline-swath test no power, and the sign of the
   previously reported correlation does not reproduce.** With the natural sign, **every**
   reference line gives a *positive* correlation (r = +0.32 to +0.70, none significant);
   the headline reference-free version over all six lines is **r = +0.505, p = 0.306**.
   The reported **r = -0.605, p = 0.279** is reproduced exactly by taking line **136** as
   the reference and negating the mark side — not line 135 as labelled. §8.

7. **A new caution, and it is the most interesting thing here.** The three *adjacent*
   lines picked up by the widened catchment — 131, 132, 140, same acquisition, same
   estimator, same counties — give **-8.0 mm, SE over lines 16.8**, against +55.0 for
   133-138. About half that gap is the known county-level pattern; the county-centred
   difference is **+26.4 mm and not significant (p = 0.136)**. But it means the +55 mm is
   a statement about six flight lines, not about gen1. §9.

---

## 1. What was fetched, and what it cost

41 gen1 tiles were already on disk. **Five were downloaded**, sequentially, ~3 s apart,
one at a time: `4342-20-62` (wabasha), `5142-04-62` (winona), `5142-11-62` (fillmore),
`5158-04-01` (winona), `5158-10-01` (fillmore) — **107 008 578 bytes**, 2.0-2.5 s each.
`county_for_lonlat` was not used; the county was taken from the `county` column of the
control mark that falls in each tile, with a candidate-county fallback list that was never
needed. Peak memory during the run stayed at 15 GiB used of 46 GiB with swap at 7 MiB.

**Those five tiles bought zero additional marks** (§3). That is the finding, not a
failure of the fetch.

---

## 2. The line assignment, taken from the data

The first pass assigned a mark to the flight line whose fitted centreline was nearest.
That is a proxy, and at a 961 m line spacing it was expected to mislabel as the radius
widened. It does not, and the reason is worth stating because it makes the rest of the
document simple.

**The tracks.** `gen1_line_tracks.py` fits one straight ground track per
`point_source_id` from **near-nadir returns only** (`|scan_angle_rank| <= 2`). A first cut
took the mean of *all* returns in each gps_time bin and got "track residuals" of
100-450 m; those were **tile edges**, not flight — where a tile boundary clips the swath
asymmetrically the bin mean is displaced by up to a swath half-width. Near-nadir returns
are either in the tile or not.

`point_source_id` **is reused across missions** in this acquisition: psid 151 spans
**337 385.7 s of gps_time (3.9 days)** over four tiles. Each psid is therefore split into
passes at gps_time gaps. **Lines 133-138 are not affected** — each is a single continuous
pass:

| psid | near-nadir bins | span (km) | heading | resid median (m) | resid p95 (m) | resid max (m) |
|---|---|---|---|---|---|---|
| 133 | 391 | 97.1 | 179.3 | 15.8 | 65.2 | 139.4 |
| 134 | 388 | 97.1 | 359.3 | 16.9 | 71.3 | 121.9 |
| 135 | 388 | 97.2 | 179.3 | 13.9 | 56.8 | 151.1 |
| 136 | 272 | 93.7 | 359.3 | 16.9 | 57.2 | 137.8 |
| 137 | 260 | 93.7 | 179.4 | 16.7 | 59.8 | 109.8 |
| 138 | 360 | 93.7 | 359.3 | 34.8 | 92.0 | 105.7 |

These are 94-97 km straight lines held to 14-35 m in the median. The track model is used
for **targeting only** — which tiles to fetch. Membership is decided by
`point_source_id`.

**The assignment, corrected.** Over the 18 marks that carry 133-138 ground:

```
line assignment: 18 marks carry 133-138 ground; the nearest-centreline proxy disagrees
                 with point_source_id on 0 of them
mixed-line tie vs single-line tie, marks with one line: max |difference| 0.000000 mm over 18 marks
marks carrying MORE THAN ONE of lines 133-138 in the report radius: 0
```

**How many of the existing 18 the proxy got wrong: none.** And no mark is mixed, so the
per-line ties equal the pooled ties exactly. §3 says why.

*(This differs from `GEN1_OWN_CONTROL_TIE.md` §3, which found seven marks with two flight
lines inside the report radius. That run used **CSF** ground, and CSF reclassifies the
class-12 overlap returns as ground. On the vendor's own bare-earth class the overlap is
simply not there.)*

---

## 3. The seam: why the catchment cannot be widened

`gen1_swath_seam.py`, cross-track distance from the whole-line fit, over eight tiles:

```
 line   subset   p0.05_m  p99.95_m   half_m             n
  133      all    -179.0     718.1    718.1     8,333,069
  133  class12     434.1     731.9    731.9     1,953,510
  133   class2    -179.0     473.2    473.2     5,418,112
  134      all    -766.8     739.7    766.8    13,592,742
  134  class12    -779.4     751.9    779.4     3,900,492
  134   class2    -499.2     490.5    499.2     8,070,444
  135      all    -771.3     728.0    771.3    12,660,200
  135  class12    -782.7     729.1    782.7     3,119,142
  135   class2    -514.1     577.0    577.0     7,865,751
  136      all    -791.6     767.7    791.6     9,391,754
  136  class12    -791.6     785.5    791.6     3,129,395
  136   class2    -749.9     504.1    749.9     5,160,328
  137      all    -755.8     733.8    755.8     8,534,547
  137  class12    -762.8     748.0    762.8     2,502,602
  137   class2    -510.6     475.6    510.6     4,918,106
  138      all    -758.8      25.0    758.8     4,670,621
  138  class12    -762.9    -451.7    762.9     1,330,771
  138   class2    -504.1      25.2    504.1     2,627,301

line spacing at N = 4883678: 942, 987, 932, 988, 956 m; mean 961 m, half-spacing 481 m
```

Read the three rows of any line together: the returns reach **~760 m** either side, the
**overlap class occupies the outer band** from ~450 m out, and **ground stops at ~500 m**.
The single-sided rows (133 at -179, 138 at +25) are tiles that do not extend that far;
136 reaching -749.9 in ground is where its neighbour is absent so nothing was cut.

**The vendor's seam is the perpendicular bisector between adjacent flight lines.** That is
why the "nearest centreline" proxy is not a proxy at all — it computes the same partition.
And it is why no radius can enlarge the set: past 481 m the ground belongs to the next
line, by construction, in the delivered file.

---

## 4. The marks-versus-tiles curve

Candidates are open (`L1O`) and urban (`L5U`) 2008 control marks within `R_m` of the
nearest of the six centrelines. `n_new_tiles` counts tiles not on disk before this run.

```
   R_m  n_cand  n_on_lines  n_tiles  n_new_tiles
  ----  ------  ----------  -------  -----------
   250       7           7        7            0
   500      18          18       13            0
   700      19          18       14            1
  1000      21          18       15            1
  1500      24          18       17            3
  2000      30          18       21            5
```

**`n_on_lines` is flat at 18 from 500 m outward.** The `n_cand` column keeps growing and
buys nothing.

**A number in the brief I could not reproduce.** The brief states that widening to 1 km
takes the set "from 18 to ~35 open/urban marks." Measured here it is **21 candidates, of
which 18 on lines**. The accompanying claim that 1 km "needs only ONE new tile" **is**
reproduced: `5142-04-62`. I do not know where ~35 came from; my count uses the 986
distinct `point_id` in the bundled CSV, `L1O`+`L5U` only, and the whole-line track fits
above.

**The 18 are the complete set.** Relaxing the requirement that a mark lie inside each
line's *observed* northing span changes nothing (the six lines span N 4 820 272 -
4 917 429 against a control-mark range of N 4 817 327 - 4 922 574, and **0** open/urban
marks within 500 m of a centreline fall outside it). Every tile needed at R = 500 m was
already on disk. There is no unread open/urban 2008 control on lines 133-138 anywhere in
the eight counties.

---

## 5. The screen table

Vendor class-2 ground, 300 m crop, `res = 5 m`, report radius 7.5 m, siting radius 5 m.
`line_proxy` is the old nearest-centreline rule, `line_psid` the `point_source_id` at the
mark; they agree everywhere.

```
                  point  cls     km  d_near_m  line_proxy  line_psid  n_report  tie_mm  sigma_mm  ladder_mm  spread_mm  ctrl_pct  dnr_mm
  ---------------------  ---  -----  --------  ----------  ---------  --------  ------  --------  ---------  ---------  --------  ------
                 L5U135  L5U  16.53     314.0         133      133.0        90    12.9      16.2       32.4      182.0      56.0   -18.0
                  L1O43  L1O  18.21     448.3         133      133.0        92   -32.7      42.6       85.2      380.0      32.0    16.0
  L1O-2124 Fillmore VRS  L1O  62.89      63.5         133      133.0        96   121.0       4.5        9.0      148.0     100.0   129.0
                 L5U172  L5U   8.18     148.0         134      134.0        72    67.9      13.7       27.3      385.0      61.0    97.0
   L1O-6123 Wabasha RTK  L1O  14.42     321.7         134      134.0        91   154.5      17.7       35.5      296.0      95.0   173.0
   L1O-6196 Wabasha RTK  L1O  24.48     170.4         134      134.0        93   118.2      15.4       30.8      189.0     100.0   104.0
   L1O-6189 Wabasha RTK  L1O  27.46     175.8         135      135.0        92    89.4      17.7       35.4      471.0      60.0    97.0
   L5U-6187 Wabasha RTK  L5U  30.61     208.2         135      135.0        83   -16.0      16.8       33.6      190.0      53.0    32.0
   L5U-6188 Wabasha RTK  L5U  29.35     117.5         136      136.0        87    82.0      30.8       61.6      360.0      90.0    80.0
                 L1O101  L1O   1.79     265.3         137      137.0        85   -11.3       9.2       18.3      110.0      63.0   -14.0
                 L5U171  L5U   2.32     479.3         137      137.0        62    47.1      29.5       59.1      345.0      77.0   143.0
   L1O-6182 Wabasha RTK  L1O   27.4     251.1         137      137.0        93     7.7      11.6       23.2      131.0      68.0     5.0
  L5U-2122 Fillmore VRS  L5U  60.63     123.9         137      137.0        90     4.9      35.7       71.3      295.0      60.0   -64.0
  L5U-2121 Fillmore VRS  L5U  62.81     365.1         137      137.0        63   -57.0      69.6      139.1      120.0      12.0   -19.0
                 L1O173  L1O   4.45     409.2         138      138.0        94    45.7       4.8        9.5      260.0      76.0    36.0
   L5U-6106 Wabasha RTK  L5U  22.07     281.6         138      138.0        80    71.4       8.1       16.3      195.0      75.0    54.0
  L5U-2119 Fillmore VRS  L5U   61.8     333.8         138      138.0        93    26.5      37.8       75.6      241.0      75.0   -63.0
  L5U-2120 Fillmore VRS  L5U  61.99     430.0         138      138.0        96   119.3      43.2       86.5      296.0      93.0   194.0
```

Two marks sit at `ctrl_pct = 100` (`L1O-2124`, `L1O-6196`) — the surveyed height is above
every ground return within 5 m — and one at `ctrl_pct = 12` (`L5U-2121`). These are exactly
the §7.1 siting failure mode. **They are not removed**; §7 measures what removing them
would do.

Column definitions are printed by the script, not typed here: `tie_mm` = surveyed minus
lidar at the report radius from that line's returns only; `sigma_mm` = half the tie spread
across the pipeline-scale radii 2.5-10 m; `ladder_mm` = the full spread; `spread_mm` =
local p05-p95 of ground z within 5 m; `ctrl_pct` = percentile of the surveyed height in
that local distribution; `dnr_mm` = MnDNR's own `Control Z - Surface Z` for the mark.

---

## 6. Per line

```
  line  n_marks  mean_mm  sd_mm  SE_mm  median_mm  km_min  km_max
  ----  -------  -------  -----  -----  ---------  ------  ------
   133        3     33.7   78.9   45.6       12.9    16.5    62.9
   134        3    113.5   43.5   25.1      118.2     8.2    24.5
   135        2     36.7   74.5   52.7       36.7    27.5    30.6
   136        1     82.0    nan    nan       82.0    29.4    29.4
   137        5     -1.7   37.6   16.8        4.9     1.8    62.8
   138        4     65.7   40.2   20.1       58.6     4.4    62.0
```

`SE_mm` is the standard error of the **mean of the per-mark ties on that line**. Line 136
rests on a single mark and has no SE. The five lines that have one give a median SE of
**25.1 mm**, which is the resolution these marks bring to a per-line question.

---

## 7. Combined, both ways, and whether screening helps

```
                             quantity     value                                           note
  -----------------------------------  --------  ---------------------------------------------
           mean of the six line means  +55.0 mm                  SE over the six lines 16.6 mm
   mean over marks, marks independent  +47.3 mm                             SE 14.1 mm, n = 18
                    median over marks  +46.4 mm                                         n = 18
             ANOVA over the six lines  F = 2.12                       p = 0.1330, df = (5, 12)
  sigma_site (within-line, all marks)   51.9 mm  pooled sd of a mark about its own line's mean
         sd over marks, ignoring line   59.8 mm                                         n = 18
```

The mean of the six line means weights line 136 (one mark) equally with line 137 (five).
That definition is carried unchanged from the previous run so the two numbers are
comparable; the marks-as-independent figure is the one that does not.

**The ANOVA does not support line structure inside 133-138** (F = 2.12, p = 0.133). The
**F = 8.63** in the brief reproduces exactly, but it is a different test on a different
set: `scipy.stats.f_oneway` over the 43-mark screen, **lines 128-149**, singleton lines
dropped — 10 groups, 38 marks, **F = 8.629, p = 4.61e-06**. Keeping the singletons and
computing the ANOVA by hand over all 15 groups gives **F = 5.94, p = 3.23e-05,
df = (14, 28)**. Both say the same thing about the 22-line span and neither says anything
about these six.

**Does screening on radius spread reduce site scatter? No.** Sweeping a cut on
`ladder_mm` over the 18 (this is a **measurement of what cutting does**; the headline
numbers above use every mark):

```
  cut_mm  n_kept  n_lines_kept  sigma_site_mm  sd_all_mm  mean_mm_kept
  ------  ------  ------------  -------------  ---------  ------------
      15       2             2            nan        nan           nan
      20       4             3           18.2       55.1          56.7
      25       5             3           16.0       52.5          46.9
      30       6             4           16.0       47.7          50.4
      40      11             5           51.1       57.5          60.1
      50      11             5           51.1       57.5          60.1
      75      14             6           46.6       53.1          56.8
     100      17             6           50.9       55.5          53.4
  no cut      18             6           51.9       59.8          47.3
```

`sd_all_mm` moves between **47.7 and 59.8 mm** across every cut from 20 mm to none, and
the mean moves between **+46.9 and +60.1 mm** — inside its own 14 mm SE. **The screen does
not buy precision**; what it buys, at cuts below 30 mm, is fewer than seven marks on four
or fewer lines.

The same sweep on the 43-mark screen gives `sd_all_mm` of **83.2, 86.2, 95.9, 97.7, 99.9,
98.1, 94.4, 91.9, 94.1 mm** at cuts of 15, 20, 25, 30, 40, 50, 75, 100 and none — which
is where the brief's "sigma_site 83-97 mm at every cut" came from. **That quantity is the
sd of the ties ignoring line, not the within-line pooled sd**; the within-line pooled sd
over the same 43 is **57.9 mm**. Both are reported here under names that say which is
which.

---

## 8. Against the pipeline's own swath constants

`data/derived/elbaext/corrections_geoid.json`,
`per_swath_internal_alignment_dxdydz_m`, vertical component, in mm — beside the
mark-derived per-line means:

```
  line  dz_pipeline_mm  mark_mean_mm
  ----  --------------  ------------
   133             0.0          33.7
   134            22.0         113.5
   135             6.2          36.7
   136            -9.8          82.0
   137           -18.4          -1.7
   138           -22.6          65.7
```

Note **line 133 is the pipeline's reference** (dz = 0), not 135.

```
                          quantity     value                                                              note
  --------------------------------  --------  ----------------------------------------------------------------
     Pearson r, all six lines, raw    +0.505              p = 0.306, n = 6; correlation is invariant to the reference line
   Pearson r, relative to line 135    +0.611                                                  p = 0.273, n = 5
   Pearson r, relative to line 133    +0.556    p = 0.331, n = 5; 133 is the pipeline's own reference (dz = 0)
     RMS residual, relative to 135   47.0 mm                                        mark minus pipeline, n = 5
   RMS residual, both mean-removed   32.2 mm                                                             n = 6
         spread of the pipeline dz   44.6 mm                                          max - min over six lines
     spread of the mark line means  115.3 mm                                          max - min over six lines
  SE of a single line mean, median   25.1 mm  over the five lines with more than one mark; line 136 has one mark and no SE
```

**Plainly: no, the larger sample gives it no power — and there is no larger sample.** The
set is the same 18 marks. The reference-free statement is **r = +0.505, p = 0.306, n = 6**.
A signal of 44.6 mm peak-to-peak is being probed with per-line means whose median SE is
25.1 mm and whose own spread is 115.3 mm. The test cannot decide this, and nothing on
lines 133-138 will make it able to.

**The previously reported r = -0.605 does not reproduce under the convention it was
labelled with.** Every reference line and both signs:

```
  ref  sign   n        r        p   rms_mm
  133    +1   5   +0.556    0.331     44.7
  133    -1   5   -0.556    0.331     54.6
  134    +1   5   -0.064    0.919     50.4
  134    -1   5   +0.064    0.919    105.6
  135    +1   5   +0.611    0.273     47.0
  135    -1   5   -0.611    0.273     52.0
  136    +1   5   +0.605    0.279     50.6
  136    -1   5   -0.605    0.279     56.1
  137    +1   5   +0.320    0.599     58.1
  137    -1   5   -0.320    0.599     93.4
  138    +1   5   +0.700    0.188     47.9
  138    -1   5   -0.700    0.188     51.9
 None    +1   6   +0.505    0.306     32.2
 None    -1   6   -0.505    0.306     46.7
```

**`r = -0.605, p = 0.279` is the `ref = 136, sign = -1` row**, not `ref = 135`. The
brief's accompanying `RMS residual 50.3 mm` matches no row (the nearest is 50.4 and 50.6).
With the natural sign — a line that reads low needs a positive `dz` added and produces a
positive tie — **every reference gives a positive r**, from +0.064 to +0.700, none of
them significant. The earlier statement "no power" was right; its sign was an artefact of
the reference and sign convention it was computed under.

---

## 9. What the adjacent lines say — the reason not to believe +55 mm is gen1's

The widened catchment did not add marks to 133-138, but it did *measure* 12 open/urban
marks on the immediately adjacent lines. Same acquisition, same estimator, same crop, same
counties, one flight line over:

```
      n  mean_mm  sd_mm
line
131   2    -29.1   29.7
132   6     25.3   51.6
140   4    -20.3   36.4

mean of the three line means -8.0 mm, SE over lines 16.8
marks as independent +1.0 mm, SE 13.9, n=12
```

Against **+55.0 ± 16.6** and **+47.3 ± 14.1** on 133-138. Welch on the marks:
**t = 2.341, p = 0.0269, df = 26.9** (A: n = 18, +47.3, sd 59.8; B: n = 12, +1.0,
sd 48.0).

**Most of that is not the flight lines.** The 30 marks split by county:

```
           n  mean_mm  sd_mm
fillmore  10     20.6   65.1
wabasha    8     76.4   56.2
winona    12      3.9   37.1
```

and Wabasha supplies **7 of the 18** on 133-138 against **1 of the 12** on the adjacent
lines. Matched within county:

```
   county  n_A    meanA  n_B    meanB      A-B
 fillmore    5    +42.9    5     -1.8    +44.8
  wabasha    7    +72.5    1   +103.9    -31.5
   winona    6    +21.6    6    -13.7    +35.3
```

County-centred, the difference falls from **+46.3 to +26.4 mm** and is **not significant
(Welch t = 1.535, p = 0.136, df = 28.0)**. The county-scale pattern is the one
`GEN1_OWN_CONTROL_TIE.md` §6 already found in MnDNR's own residuals (Winona -101 mm,
Dodge +89 mm on open cover).

**Measured, not interpreted:** two adjacent groups of flight lines from the same 2008
acquisition give +47.3 ± 14.1 and +1.0 ± 13.9 mm against the same control set, and about
half of that separation is accounted for by which counties the marks are in.

Pooling all nine lines: **30 marks, mean of nine line means +34.0 mm, SE over lines
15.8; marks as independent +28.8 ± 10.8 mm**. ANOVA over the nine lines (singletons
dropped, as in §7): **F = 2.861, p = 0.0291**.

---

## 10. The one enlargement axis that is left, and what it costs

Radius is exhausted. The only other way to add marks is to drop the open/urban
restriction. Measured, at R = 500 m, all cover classes: **45 candidates, of which 31 are
readable and 14 need tiles not on disk** (14 more tiles, ~350 MB, all for vegetated
marks). The 31 readable ones, each again carrying exactly one line:

```
        n  mean_mm  median_mm  sd_mm  SE_mm
L1O     8     61.5       67.5   69.2   24.5
L2T     5     62.5       50.6   63.8   28.5
L3B     4    -20.6      -68.7  166.3   83.1
L4F     3     -1.2       -3.5   15.2    8.8
L5U    10     35.9       36.8   52.1   16.5
other   1    147.8      147.8    NaN    NaN
```

Per line, all covers against open/urban only:

```
     all_covers        open_urban
           size   mean       size   mean
line
133           5   17.1          3   33.7
134           6  110.6          3  113.5
135           4   48.7          2   36.7
136           3   39.3          1   82.0
137           6    1.1          5   -1.7
138           7   22.5          4   65.7

mean of line means: all covers +39.9 (SE 15.7), open/urban +55.0 (SE 16.6)
marks as independent: all covers +39.5 (SE 14.5, n=31)
```

**Thirteen more marks move the SE from 16.6 to 15.7 mm and the estimate by 15 mm.** The
`L3B` stratum alone has an sd of 166 mm. This is the trade in numbers; the decision — and
whether to spend 14 tiles on the rest of it — is Andy's. It is reported here, not adopted:
**the headline numbers in this document are open/urban only.**

---

## 11. Parameters, and the ones that are mine

| parameter | value | source | effect |
|---|---|---|---|
| `lines` | 133-138 | **andy** | the scope of the task |
| `cover_classes` | `L1O`, `L5U` | repo | `GEN1_OWN_CONTROL_TIE.md` §4; the alternative is measured in §10 and not adopted |
| `ground_source` | `vendor_class2` | **andy** | CSF is ~460 s/tile. Note it is *also* what makes each mark single-line (§3); CSF would reclassify the class-12 overlap as ground and reintroduce mixing |
| `geoid_shift_m` | 0.0 | repo | asserted, not approximated — both sides NAVD88(GEOID03) |
| `swath_shift_m` | (0, 0, 0) | repo | no gen2-derived term |
| `res_m` / report radius / siting radius | 5.0 / 7.5 / 5.0 m | repo | `corrections_geoid.json`, `tie.py`, §7.1 |
| `csf_crop_halfwidth_m` | 300.0 | repo | carried unchanged from `elba_absolute_tie.py` so the runs are one method |
| **`catchment_radius_m`** | **2000.0** | **MINE** | widened from 500 m *to measure where the catchment ends*. It selects candidates to READ; membership is decided by `point_source_id`. Effect: 30 candidates instead of 18, **18 on lines either way**, 5 tiles fetched, 107 MB, 0 marks gained. Flagged in the run banner. |
| **`NADIR_DEG`** | **2.0** | **MINE** | the `|scan_angle_rank|` band used for the nadir track. Targeting only; it never touches a tie |
| **`GAP_S`** | **120.0** | **MINE** | splits a psid into passes at gps_time gaps. Targeting only. Lines 133-138 are single passes and are unaffected |
| **`STRIDE`** | **7** | **MINE** | subsampling for the track fit. Targeting only |
| **`ladder_mm` cut sweep** | 15-100 mm | **MINE** | a **measurement** (§7), never applied to a headline number |

The run banner prints the `MINE` warning; the ledger is in `.trust/runs/`.

---

## 12. What is still open

* **The +55 mm is a six-flight-line statement, not gen1's datum** (§9). The adjacent lines
  give +1.0 ± 13.9 mm on 12 marks. Whether that is line-to-line offset, a county-scale
  pattern in the delivered surface, or site scatter is **not resolved** by 30 marks, and
  the county-matched test (p = 0.136) does not decide it.
* **The all-cover enlargement** (§10) is measured for the 31 marks whose tiles are on
  disk; 14 more tiles would complete it. Not fetched.
* **The pipeline-swath test is not merely underpowered, it is closed on these lines** —
  there is no further control to add, so the 44.6 mm swath signal cannot be checked
  against the 2008 marks at all. If it needs checking, it needs a different observable.
* **Everything `GEN1_OWN_CONTROL_TIE.md` §10 left open is still open**, in particular the
  NAD83 realization term and the unpublished 2008 bias adjustment.

---

## Appendix — verification status

| claim | how |
|---|---|
| the seam, class-12 overlap, and the 481 m half-spacing | `gen1_swath_seam.py`, pasted verbatim in §3 |
| lines 133-138 are single 94-97 km passes | `gen1_line_tracks.py`, per-pass and whole-line fits |
| psid reuse across missions (psid 151, 337 385.7 s) | `line_tracks.json` `tile_lines`, min t0 to max t1 over the four tiles holding it |
| 0 of 18 marks mislabelled by the proxy; 0 mixed marks | `gen1_more_marks_report.py`, printed |
| per-line ties identical to pooled ties | max \|difference\| **0.000000 mm**, printed by the same script |
| every tie, screen statistic and per-line mean | `gen1_more_marks_tie.py --radius-m 2000`, vendor ground, one run |
| the marks-vs-tiles curve and `n_new_tiles` | the same run, against the 41-tile list held before the five downloads |
| **the brief's "~35 marks at 1 km"** | **NOT reproduced** — measured 21 candidates / 18 on lines. Source of ~35 unknown |
| **the brief's r = -0.605 "relative to line 135"** | **reproduced only as `ref = 136, sign = -1`**; §8 gives all fourteen variants |
| **the brief's RMS residual 50.3 mm** | **NOT reproduced** — no reference/sign variant gives it; nearest are 50.4 and 50.6 |
| the brief's F = 8.63 | **reproduced exactly** as `f_oneway` over lines 128-149 with singletons dropped (10 groups, 38 marks, p = 4.61e-06) |
| the brief's "sigma_site 83-97 mm" | **reproduced** as the sd of the 43 ties *ignoring line* (83.2-99.9 across the cut sweep). The within-line pooled sd over the same 43 is 57.9 mm |
| the county pattern the §9 confound rests on | measured here on 30 marks; consistent with `GEN1_OWN_CONTROL_TIE.md` §6, which was not re-run |
| the 2008 control's datum | still a **dataset-level assertion** (`lidar_semn2008.html`); the validation reports state no datum and no geoid |

# gen1 against its OWN 2008 ground control

**Date:** 2026-08-26
**Scripts:** `analysis/groundtruth/parse_mndnr_2008_control.py` (the reports → a bundled
checkpoint CSV, with the sign test), `analysis/groundtruth/gen1_own_control_tie.py`
(every measurement below)
**Data added:** `src/lidar_diff_icp/groundtruth/data/mn_dnr_2008_control_semn.csv`,
1 004 surveyed control points across gen1's eight-county footprint
**Bears on:** `analysis/ABSOLUTE_BASIS_ELBA.md` (the +22.7 mm constant),
`data/derived/elba_fulldensity/z_before_absolute.{npy,json}`,
`analysis/ADDITIONAL_GROUND_CONTROL.md` §5 (which flagged the sign tension as a task)

Every table below is printed by the script that computed it, through
`trust/provenance.py`. Nothing in `pipeline.py` or `coreg.py` was touched. The one
parameter I chose is flagged `MINE` in the run banner and its effect is measured (§8).

---

## Bottom line

1. **The sign convention is settled, by arithmetic.** The MnGeo validation tables'
   `Error` column is **`Control Z - Surface Z`** – exact on **1 022 of 1 022** rows of
   the eight SE-MN reports; the other order misses by up to **1.08 m**. That is the *same* sign
   family as `groundtruth.tie`'s `tie = surveyed - z_lidar`. **No sign flip separates the
   two numbers.** They really do point opposite ways.

2. **Our own point-cloud read at gen1's own control, screened for siting, is negative.**
   At the two best-sited open marks near Elba the tie is **-19.4 ± 7.9 mm** (`L1O101`,
   1.34 km, local spread 110 mm) and **-37.8 ± 4.8 mm** (`L1O59`, 6.91 km, spread
   140 mm) – gen1 reads **above** its own control, no chain, no geoid conversion, no
   lateral extrapolation.

3. **The +22.7 mm anchor is not refuted, it is unresolved.** Anchor minus the median of
   the four open 2008 marks is **+51.3 mm, 1.29x the anchor's own sigma_total of
   39.7 mm**. The two control sets do not contradict each other at any useful confidence;
   what they do is show that a constant whose uncertainty already exceeds its value can
   change sign under an equally defensible, better-sited control set.

4. **The two anchor marks fail the §7.1 siting screen, read with this same estimator.**
   Local p05-p95 at R = 5 m: **442 mm at 2210** with the surveyed height at **p93**, and
   **865 mm at 2036** at **p90**. The four open 2008 marks sit at p18-p88 with spreads
   110-650 mm, and the two used above at p45 and p18.

5. **The vendor's `Surface Z` is not the problem.** Our reconstruction tracks the MnDNR
   DEM at the same 16 marks: median difference **-13.7 mm**, Pearson **r = +0.701**. The
   `Error` column is a fair proxy for our residual to about a centimetre in the median,
   which is why the whole 1 004-point set can be read at all.

6. **`z_before_absolute.npy` does not need its arithmetic corrected – it needs its claim
   weakened.** The constant, the sign convention and the code are right. What is wrong is
   reading `+22.7` as a determination; §7.

7. **The bias adjustment is real, applied, and its value is published nowhere.** Metadata
   process step 8 states it; none of the eight validation reports, the two data READMEs or
   the dataset metadata gives a number. **gen2 carries the same unpublished constant**
   (`analysis/LIDAR_DOCUMENTATION_MINE.md` §1.2). The request is drafted in §9.

---

## 1. What was fetched, and what it cost

Six county validation reports, one request each, spaced 3 s: dodge, fillmore, freeborn,
houston, mower, steele – **3 557 501 bytes** total (595 057 + 613 608 + 615 853 +
564 743 + 580 050 + 588 190). Plus six directory listings (~21 kB), the dataset metadata
page `lidar_semn2008.html` (25 563 bytes), one MnGeo directory listing, and Winona's two
data READMEs (54 265 + 1 037 759 bytes). Winona, Wabasha and Olmsted were already in the
scratchpad from `ADDITIONAL_GROUND_CONTROL.md`. **No lidar was downloaded**; every point
cloud read was already on disk.

**Freeborn is excluded.** Its report exists in the same tree, but Freeborn is not among
the eight counties the dataset metadata lists a per-county RMSE for, so it is not this
acquisition. The bundled CSV covers Dodge, Fillmore, Houston, Mower, Olmsted, Steele,
Wabasha and Winona.

---

## 2. The sign convention, pinned three ways

**(a) By arithmetic, which is decisive.** Over every table of all eight SE-MN reports:

```
Control - Surface : 1022/1022 rows agree, max |resid| 0.000000 m
Surface - Control : 0/1022 rows agree, max |resid| 1.080000 m
```

**(b) By reproducing what the reports print.** Parsing each report's *overall* table and
squaring the parsed `Error` column reproduces the RMSE printed at the foot of that same
table:

```
county         n  n_pub    RMSE  in-report  published
dodge        120    121  0.1295      0.129      0.129
fillmore     127    128  0.1546      0.155      0.155
houston      109    134  0.1344      0.134      0.110
mower        114    115  0.1620      0.161      0.161
olmsted      125    125  0.1167      0.117      0.117
steele       137    137  0.1253      0.125      0.125
wabasha       97     97  0.1057      0.067      0.106
winona       175    176  0.1601      0.161      0.161
```

`n_pub` and `published` are from the dataset metadata's "MnDNR's Tests" paragraph.
Dodge, Fillmore, Mower and Winona come up one row short – rows split across a page
break, the same shortfall `ADDITIONAL_GROUND_CONTROL.md` §1.1 reports for the 2021 survey
report. `wabasha`'s `in-report` 0.067 is the first RMSE its file prints, which belongs to
a per-class table, not the overall one; the parsed 0.1057 matches its published 0.106.
**Houston is a genuine discrepancy in the published record**: the report's own overall
table holds 109 rows and prints 0.134, while the metadata says 0.110 m over 134 points.
That goes in the email.

**(c) By physics, which is the check that the column *labels* are not swapped.** If
`Error = Control - Surface`, then vegetation left in the ground class must push the
*surface* up and the *Error* negative. Over the whole 1 004-point set:

```
  cover              what    n  mean_mm  median_mm  SE_mm
  -----  ----------------  ---  -------  ---------  -----
    L1O      open terrain  241    +17.1      +28.0    7.0
    L2T  tall weeds/crops  223    -92.2      -77.0    8.9
    L3B   brush/low trees  168   -114.9     -112.0    9.9
    L4F          forested  165    -45.9      -21.0   10.2
    L5U             urban  175     +5.1      +15.0    8.6
  other         unclassed   32     +9.7       +1.0   15.5
```

Brush and low trees is the most negative class (-114.9 ± 9.9 mm), tall weeds and crops
next (-92.2 ± 8.9), and open terrain is the only class that is positive
(+17.1 ± 7.0). That is the correct physical ordering for `Control - Surface`, and it
would be inverted if the labels were the other way round. **So: negative `Error` means
the delivered 2008 surface reads ABOVE the mark, and `Error` is directly comparable to
`tie_mm` without a sign change.**

**Where the previous reading went.** `ADDITIONAL_GROUND_CONTROL.md` §5 reported a
**-22.4 mm** intercept within 10 km, pooled over all cover classes, and read it as "the
2008 DEM sits ~22 mm above the surveyed control near Elba". The sign of that statement is
**right**. §5 was also right to refuse to call it a disagreement with the anchor.

---

## 3. The estimator, and what it was pointed at

Our own residual is read with `groundtruth.tie.estimate_tie` – the same
curvature-unbiased order-2 local surface, the same radius ladder, the same
`ground_source = csf`, the same 300 m CSF crop, and the same 5 m grid as the run that
produced the anchor. Two things differ from the anchor run, and both *remove* a term:

* **`geoid_shift_m = 0`.** The 2008 control is NAVD88(GEOID03) and so is the raw gen1
  cloud, so no conversion enters. This is also why the comparison is fair: the anchor's
  tie is `H18_surveyed - (H03_gen1 + (N03 - N18))`, which is `H03_surveyed - H03_gen1`.
  **The geoid cancels out of the anchor's tie**, so `tie21` and `tie08` are the same
  quantity, expressed in the same frame, and the +67 mm GEOID03/GEOID18 step is not
  available to explain the difference between them.
* **No chain.** Every mark below lies under a gen1 flight line inside a tile already on
  disk. That is the whole value of this control set.

**16 of the 1 004 control points fall inside the 13 gen1 tiles on disk.** Each is
reported separately, before any aggregation.

```
                 point  cover     km        tile    n  lines  tie_mm  sigma_mm  tie_nolat_mm  dnr_mm  spread_mm  ctrl_pct
  --------------------  -----  -----  ----------  ---  -----  ------  --------  ------------  ------  ---------  --------
                L1O101    L1O   1.34  4342-30-64  166      2   -19.4       7.9         -17.3   -14.0        110        45
                 L3B99    L3B   1.70  4358-29-01  109      1  -120.2      80.0        -117.3   +15.0        468        29
                L5U171    L5U   1.88  4342-30-64  141      2   +73.4      42.5         +68.3  +143.0        425        80
                 L2T51    L2T   3.20  4342-29-63  191      2   +72.6      10.2         +45.8   +45.0        386        61
                 L2T58    L2T   3.69  4358-29-01  187      2  -140.4      11.6        -136.5   -98.0        349         9
                 L2T54    L2T   4.62  4342-28-64  108      1   +24.4      13.2         +23.7    -5.0        403        52
                L3B143    L3B   4.77  4342-30-63  101      1  -157.6      42.9        -174.1  -184.0        533        21
                 L1O59    L1O   6.91  4358-29-03   95      1   -37.8       4.8         -45.6   -39.0        140        18
                L1O144    L1O   7.09  4358-29-03  101      1   -75.0       4.6         -51.6   -52.0        337        34
  L4F-6127 Wabasha RTK    L4F   7.65  4342-29-61   98      1  -247.8      41.5        -169.3   -48.0        909        36
  L2T-6126 Wabasha RTK    L2T   7.65  4342-29-61  100      1  -274.6      55.4        -203.0  -125.0        875        24
                 L3B94    L3B  11.62  4358-26-03  101      1  -163.0      36.0        -124.6  -111.0        432         9
  L1O-6104 Wabasha RTK    L1O  13.46  4358-26-03  203      2  +196.3      86.7        +194.4   +27.0        650        88
  L2T-6102 Wabasha RTK    L2T  13.64  4358-26-03  201      2  -171.6      20.4        -177.0  -302.0        358        10
  L3B-6103 Wabasha RTK    L3B  13.70  4358-26-03  199      2    -9.5      18.9         -89.4   -61.0       1899        51
  L5U-6101 Wabasha RTK    L5U  13.93  4358-26-03  100      1   +96.1      42.6         +64.0  +118.0        478        63
```

`tie_nolat_mm` is the same read with the elbaext Nuth & Kaaeb shift withheld; the lateral
term moves these marks by -78.5 to +79.9 mm, which is larger than its -6.8/+12.3 mm effect at
the anchors and is a warning about extrapolating one shift vector.

The full radius ladder, which is the honest uncertainty on this data:

```
                 point  cover  R=2.5   R=5  R=7.5  R=10  R=15  R=20  R=25
  --------------------  -----  -----  ----  -----  ----  ----  ----  ----
                L1O101    L1O     -5    -4    -19   -19    +2   +10   +21
                 L3B99    L3B    +27   -18   -120  -133   -59   -73   -86
                L5U171    L5U   +151  +127    +73   +66   +74   +78   +74
                 L2T51    L2T    +73   +61    +73   +81  +110   +97   +91
                 L2T58    L2T   -121  -138   -140  -145  -152  -172  -201
                 L2T54    L2T     +3   +30    +24   +12   -29    -5   +32
                L3B143    L3B   -206  -215   -158  -129   -28   +80  +156
                 L1O59    L1O    -40   -40    -38   -31   +43  +123   +58
                L1O144    L1O    -80   -83    -75   -74   -95   -57    +6
  L4F-6127 Wabasha RTK    L4F   -170  -253   -248  -244  -194  -172  -210
  L2T-6126 Wabasha RTK    L2T   -182  -259   -275  -293  -299  -294  -289
                 L3B94    L3B   -163  -120   -163  -192  -235  -297  -357
  L1O-6104 Wabasha RTK    L1O    +61   +89   +196  +235  +218   +94  +216
  L2T-6102 Wabasha RTK    L2T   -181  -165   -172  -206  -164   -77   -22
  L3B-6103 Wabasha RTK    L3B    +25   +28     -9    -1   +35   +53   -66
  L5U-6101 Wabasha RTK    L5U   +162  +117    +96   +77   +90  +106   +67
```

`L1O101` is flat to 15 mm from R = 2.5 to 10 m (-5, -4, -19, -19) and `L1O59` to 9 mm
(-40, -40, -38, -31). The vegetated marks are not: `L3B143` walks 371 mm across the
ladder and `L3B94` 237 mm.

**Seven marks have two flight lines inside the report radius.** Read one line at a time:

```
                 point  cover  line    n  tie_mm  sigma_mm  mixed_tie_mm
  --------------------  -----  ----  ---  ------  --------  ------------
                L1O101    L1O   136   78   -20.1       9.1         -19.4
                L1O101    L1O   137   88   -19.1      10.5         -19.4
                L5U171    L5U   136   64   +75.8      17.8         +73.4
                L5U171    L5U   137   77   +80.2     157.3         +73.4
                 L2T51    L2T   134  104   +78.8      15.5         +72.6
                 L2T51    L2T   135   87   +60.7       6.5         +72.6
                 L2T58    L2T   140   94  -169.6      52.9        -140.4
                 L2T58    L2T   141   93  -112.0      11.9        -140.4
  L1O-6104 Wabasha RTK    L1O   145  114  +163.9     109.0        +196.3
  L1O-6104 Wabasha RTK    L1O   146   89  +214.7      92.9        +196.3
  L2T-6102 Wabasha RTK    L2T   143  112  -201.2      15.9        -171.6
  L2T-6102 Wabasha RTK    L2T   144   89  -129.5      42.5        -171.6
  L3B-6103 Wabasha RTK    L3B   143  107   -37.4      39.1          -9.5
  L3B-6103 Wabasha RTK    L3B   144   92   +48.0      17.7          -9.5
```

At `L1O101` the two lines agree to **1.0 mm** (-20.1 on line 136, -19.1 on 137), so the
headline mark carries no line-mixing term at all. `L2T58` and `L3B-6103` disagree by 58
and 85 mm between lines and should not be pooled.

---

## 4. By land cover – and which classes are usable

```
  cover              what  n  mean_tie_mm  median_tie_mm  SE_mm  mean_spread_mm  mean_dnr_mm
  -----  ----------------  -  -----------  -------------  -----  --------------  -----------
    L1O      open terrain  4        +16.0          -28.6   61.2             309        -19.5
    L2T  tall weeds/crops  5        -97.9         -140.4   64.2             474        -97.0
    L3B   brush/low trees  4       -112.6         -138.9   35.7             833        -85.2
    L4F          forested  1       -247.8         -247.8     --             909        -48.0
    L5U             urban  2        +84.7          +84.7   11.4             452       +130.5
```

**Only `L1O` (open terrain) is usable, and only after the §7.1 screen.** The mean over
four open marks is +16.0 mm and the median -28.6 mm, which is not a contradiction so much
as a statement that one mark (`L1O-6104`, +196.3 ± 86.7 mm, spread 650 mm, surveyed
value at p88) dominates a four-point mean. The classes with vegetation in them are worse
than useless as a datum: `L4F` gives -247.8 mm, `L2T` a mean of -97.9, `L3B` -112.6, and
their mean local spreads are 909, 474 and 833 mm.

This reproduces, at 16 marks, exactly what the whole 1 004-point set says by class in §2:
vegetation lifts the ground surface, and any pooled statistic over these classes measures
canopy, not datum.

**The siting screen (`ADDITIONAL_GROUND_CONTROL.md` §7.1) applied to the anchor marks,
with this same estimator and the same 5 m radius:**

```
          mark  type        tile   n  spread_mm  ctrl_pct         sigma_mm
  ------------  ----  ----------  --  ---------  --------  ---------------
  2210_2021_MN   NVA  4342-29-61  56        442        93              5.3
  3056_2021_MN   VVA  4342-29-61  55        862        65             15.9
  2024_2021_MN   NVA  4342-28-61  84        620        95             47.7
  2036_2021_MN   NVA  4358-29-03  48        865        90             41.6
  2099_2021_MN   NVA          --  --         --        --  no tile on disk
  3089_2021_MN   VVA          --  --         --        --  no tile on disk
```

Both anchors sit in the top decile of their own local return distribution (p93, p90) with
local spreads of 442 and 865 mm. The two 2008 marks the reconciliation leans on sit at
p45 and p18 with spreads of 110 and 140 mm.

A caution about reading too much into that. Over the 16 control marks, `ctrl_pct` and
`tie_mm` are strongly related:

```
                                        stat           value
  ------------------------------------------  --------------
                                     n marks              16
              Pearson r (ctrl_pct vs tie_mm)          +0.836
                    slope, mm per percentile  +4.44 +/- 0.78
                       fitted tie at p50, mm           -12.5
  fitted tie at p93 (mark 2210's siting), mm          +178.6
  fitted tie at p90 (mark 2036's siting), mm          +165.3
```

but this relation is **partly tautological** – `tie = control - z_lidar` and `ctrl_pct`
is the rank of `control` in the local `z` distribution, so both move with `control` by
construction. The fitted slope of +4.44 mm per percentile is the size the arithmetic
alone predicts: the mean p05-p95 spread over these 16 marks is 547 mm, which is
6.08 mm per percentile. It is **not** independent evidence that the anchors' ties are
inflated, and the extrapolated +178.6 mm at p93 is 8.4 times the anchor's measured
+21.3 mm, so the relation plainly does not transfer between control sets. What
survives is the plain observation: the marks the anchor rests on are sited where the
estimator is least stable, and the 2008 marks used here are not.

---

## 5. Is the vendor's `Surface Z` a different surface from our reconstruction?

```
                       stat   value
  -------------------------  ------
                    n marks      16
  median(ours - theirs), mm   -13.7
    mean(ours - theirs), mm   -16.4
     RMS(ours - theirs), mm    94.1
                  Pearson r  +0.701
```

Our tie and the MnDNR table's own `Error` agree in the median to **-13.7 mm** and
correlate at **r = +0.701** over the 16 marks, with a 94.1 mm RMS difference driven by
the vegetated marks where the two surfaces genuinely differ. **The vendor column is a
usable proxy for our residual at the ~centimetre level in the median, and not at the
per-mark level.** That is enough to read the 1 004-point set for spatial structure, and
not enough to use any single row of it as a tie.

---

## 6. What the whole 1 004-point set says

The MnDNR residual (`Control Z - Surface Z`, positive when the delivered DEM reads below
the mark) by distance from the Elba reference point, with the open-cover subset broken
out:

```
      km    n  mean_mm  median_mm  SE_mm  L1O_mean_mm  n_L1O
  ------  ---  -------  ---------  -----  -----------  -----
     0-5   11    -20.4       -6.0   27.0        -14.0      1
    5-10   25    -22.6      -35.0   14.1        -37.6      5
   10-20   63    -74.2      -61.0   14.2        -69.7     16
   20-40  238   -102.1      -97.0    7.8        -59.3     43
  40-200  667    -18.5       -1.0    5.2        +45.3    176
```

Plane fits to the same residual about the Elba reference point. The pooled fit is the one
that produced the "2008 surface sits ~22 mm above its control" reading; the open-cover
fit is the stratum an NVA mark belongs to:

```
        stratum  radius_km    n    intercept_mm    dE_mm_per_km     dN_mm_per_km  resid_rms_mm
  -------------  ---------  ---  --------------  --------------  ---------------  ------------
     all covers         10   36  -22.4 +/- 13.0  -1.32 +/- 3.06   -0.40 +/- 2.67            74
     all covers         20   99   -45.3 +/- 9.6  -5.35 +/- 1.12   +1.81 +/- 1.01            92
     all covers         30  208   -69.3 +/- 6.6  -3.70 +/- 0.48   +1.67 +/- 0.44            92
     all covers         50  445   -81.9 +/- 5.8  -1.40 +/- 0.24   +0.97 +/- 0.26           108
  L1O open only         10    6  -27.3 +/- 12.0  -2.45 +/- 2.35  +10.51 +/- 3.81            19
  L1O open only         20   22  -34.1 +/- 17.2  -6.77 +/- 1.77   +2.38 +/- 1.72            68
  L1O open only         30   42  -51.9 +/- 11.0  -4.12 +/- 0.81   +1.89 +/- 0.74            67
  L1O open only         50   90   -47.5 +/- 9.0  -2.00 +/- 0.37   +0.67 +/- 0.40            77
```

Two things to take from this, and one not to.

* **The pooled -22.4 mm intercept within 10 km is not a pooling artefact.** Restricted to
  open cover it is **-27.3 ± 12.0 mm** (n = 6), with a residual RMS of 19 mm against the
  pooled fit's 74 mm. The open stratum is *cleaner* and says the same thing slightly more
  strongly.
* **The 2008 surface does not sit at one level across the project**, so this is a local
  statement, not a project constant. Per-county mean of the MnDNR residual on open cover:

```
    county   n  mean_mm  median_mm  SE_mm  km
  --------  --  -------  ---------  -----  --
     dodge  35    +89.0      +87.0   18.1  72
  fillmore  23    -22.2      +23.0   22.1  49
   houston  23     -4.1       -1.0   20.7  66
     mower  25    +26.0      +74.0   23.8  76
   olmsted  23     -5.6       -6.0   12.2  35
    steele  55    +78.8      +68.0   10.3  98
   wabasha  19    +30.4      +24.0   14.2  33
    winona  38   -100.9      -88.5   12.7  25
```

  Winona, the county Elba is in, is **-100.9 ± 12.7 mm** on open cover, while Dodge is
  **+89.0 ± 18.1** and Steele **+78.8 ± 10.3**. A spread of ~190 mm between counties
  in the delivered surface's residual against its own control is the single most
  interesting number the whole 1 004-point set contains, and it is the reason the email in
  §9 asks whether the bias adjustment was applied once or per block.
* **What not to take from it:** the near-Elba open-cover fit has n = 6 and its dN gradient
  (+10.51 ± 3.81 mm/km) is not stable against the 20 km fit (+2.38 ± 1.72). Six marks
  do not constrain a plane. The intercept is the only term worth quoting.

---

## 7. The reconciliation, and the verdict on +22.7 mm

```
                                             quantity               mm                                                            note
  ---------------------------------------------------  ---------------  --------------------------------------------------------------
                 2008 control, L1O101 (1.34 km, open)    -19.4 +/- 7.9                                            spread 110 mm at p45
                  2008 control, L1O59 (6.91 km, open)    -37.8 +/- 4.8                                            spread 140 mm at p18
                 2008 control, L1O144 (7.09 km, open)    -75.0 +/- 4.6                                            spread 337 mm at p34
  2008 control, L1O-6104 Wabasha RTK (13.46 km, open)  +196.3 +/- 86.7                                            spread 650 mm at p88
  2008 control, open marks, inverse-variance weighted            -51.1  weights are each mark's own radius sigma, as combine_ties does
                 2008 control, open marks, plain mean            +16.0                             SE of the mean over 4 marks 61.2 mm
                     2008 control, open marks, median            -28.6                                                         4 marks
  2021 anchor, 2210_2021_MN (chain + geoid + lateral)   +21.3 +/- 12.4                  the constant z_before_absolute.npy is built on
  2021 anchor, 2036_2021_MN (chain + geoid + lateral)   +28.9 +/- 27.0                  the constant z_before_absolute.npy is built on
       2021 anchor, combined (z_before_absolute.json)   +22.7 +/- 39.7                         42.6 mm unmodelled bound held beside it
    anchor - (2008 control, median of the open marks)            +51.3                             1.29 x the anchor's own sigma_total
```

**Same sign convention on both sides, the same estimator, the same ground source, the
same 5 m grid, and the geoid cancelling out of both.** The two numbers are directly
comparable, and they differ by **+51.3 mm, which is 1.29x the anchor's own sigma_total**.
Repeating the entire run on vendor class-2 ground instead of CSF gives **+46.1 mm,
1.16 sigma** (§8), so nothing below turns on the ground source.

So, taking the four possibilities in the order they were asked:

* **Wrong in sign?** No. The sign convention is identical on both sides and verified two
  ways. The `Error` column being `Control - Surface` does *not* rescue the anchor.
* **Wrong in magnitude?** Not demonstrably. +22.7 ± 39.7 mm and a 2008-control median of
  -28.6 mm over four open marks are 1.29 sigma apart on the anchor's own budget alone,
  before the 2008 control's own uncertainty is added. **Neither set resolves gen1's
  absolute level from zero.**
* **An artefact of the two marks' poor siting?** Partly, and this is where the weight of
  the evidence sits. Both anchor marks fail the §7.1 screen (p93/p90, spreads 442/865 mm)
  where the 2008 marks used do not (p45/p18, spreads 110/140 mm). But the tautology in §4
  means the siting statistic cannot be turned into a correction, and gen2's own read at
  2210 – no chain, no geoid, no lateral term – is **-0.7 ± 11.5 mm**
  (`z_before_absolute.json`, `gen2_absolute_offsets`), which is *not* what a mark that
  systematically manufactures a large positive tie looks like.
* **Genuinely different from what gen1's own control says?** Yes, in the point estimate,
  and the most likely reason is a term neither set can see. The 2008 control's heights
  are GPS-derived in the NAD83 realization of 2008 and the 2021 marks in NAD83(2011); the
  anchor's tie is `h_2011(mark) - h_2008(lidar)` after the geoid cancels, so it carries
  the **ellipsoid-height realization difference** in full while the 2008 comparison
  cancels it. I tried to measure that term with `pyproj` 3.6.1 / PROJ 9.4.0 and
  **could not**: NAD83(CORS96)->NAD83(2011) resolves to a Ballpark (null) operation and
  NSRS2007->2011 pivots through WGS 84 for a null result, so PROJ here has no grid for it.
  **This is a named, unmeasured candidate, not a measured explanation**, and the existing
  budget already says so: its "horizontal datum" row records that the NAD83 realization
  difference "is not modelled by PROJ here and is therefore not in this number."

### Does `z_before_absolute.npy` need correcting?

**Its arithmetic does not. Its interpretation does.** Concretely:

* The array, the constant `datum_constant_mm = 22.66482504973182`, the sign convention and
  the budget are all internally correct and reproduce from the code. **No recomputation is
  required and none should be done on the strength of this document.**
* What should change is what the product is allowed to claim. It is currently described as
  "gen1 ... on the surveyed NAVD88(GEOID18) datum". A second, independent, better-sited
  control set on gen1's own geoid puts the same quantity at **-28.6 mm**. The honest
  description is that **gen1's absolute level at Elba is unresolved within about
  ±50 mm**, and `+22.7` is one of two defensible estimates of opposite sign.
* The concrete consequence for the pending decisions: **`z_before_absolute` must not be
  the thing that settles the re-gauge**. `FRAME_2026-08-26-PM.md` §"Do NOT bundle the
  re-gauge" notes that re-gauging on swath 137 raises gen1 by 32.5 mm – which is smaller
  than the disagreement between the two control sets. A LEVEL decision cannot be taken
  against a level that is not determined.
* Nothing here touches the **DoD**. `dod.npy` is unchanged, `implied_dod_shift_mm` is a
  derived note, and the divides still pin the epoch difference 20x more precisely than any
  mark. The finding is about the absolute level only.

**Recommended edit, not made here** (it changes a shipped product's sidecar and belongs to
Andy): add a `caveats` entry to `z_before_absolute.json` naming this document, the
-28.6 mm figure, and the 1.29 sigma separation.

---

## 8. Parameters, and the one that is mine

| parameter | value | source | effect, measured |
|---|---|---|---|
| `ground_source` | `csf` | `corrections_geoid.json` | the whole run repeated on vendor class 2: median absolute change **6.5 mm** over the 16 marks, max 36.2 mm. On the open marks: `L1O101` -19.4 → -12.3, `L1O59` -37.8 → -34.6, `L1O144` -75.0 → -75.0. The open-mark median moves -28.6 → **-23.4 mm** and the anchor gap +51.3 → **+46.1 mm** (1.29 → 1.16 sigma). **The verdict does not depend on the ground source.** |
| `res_m` | 5.0 | `corrections_geoid.json` | sets the radius ladder; the whole ladder is reported |
| `siting_radius_m` | 5.0 | `ADDITIONAL_GROUND_CONTROL.md` §7.1 | kept at 5 m so these 16 measurements are comparable to the 17 there |
| `geoid_shift_m` | 0.0 | `lidar_semn2008.html` | both sides are GEOID03; §3 shows the geoid cancels from the anchor too |
| `lateral_shift_m` | (-0.7498, -0.1893) | `corrections_geoid.json` | reported applied *and* withheld in every row: -78.5 to +79.9 mm across the 16 marks |
| **`csf_crop_halfwidth_m`** | **300.0** | **MINE** | copied unchanged from `elba_absolute_tie.py` so the two runs are one method. Its effect is bounded by the vendor-ground re-read above, which uses no crop at all: +7.1 mm at `L1O101` and +3.2 mm at `L1O59`, the two marks the conclusion rests on |

The run banner prints the `MINE` warning; the ledger is in `.trust/runs/`.

---

## 9. The bias adjustment

**What the record says, verbatim** (`lidar_semn2008.html`, process step 8):

> AeroMetric provided Quality Assurance and Quality Control (QA/QC) data for this project.
> AeroMetric captured 127 QA/QC points in multiple land cover categories that were used to
> test the accuracy of the lidar ground surface. TerraScan's Output Control Report (OCR)
> was used to compare the QA/QC data to the lidar data. This routine searches the lidar
> dataset by X and Y coordinate, finds the closest lidar point and compares the vertical
> (Z) values to the known data collected in the field. **Based on the QA/QC data, a bias
> adjustment was determined, and the results were applied to the lidar data.** A final OCR
> was performed with a resulting RMSE of 0.109 meters.

**So: one was applied, and its value is stated nowhere I can find.** Checked and found
silent on it: all eight SE-MN county validation reports (they are pure tables and charts – no
prose, no datum, no geoid, no methods section), `lidar_semn2008.html` itself beyond the
paragraph above, Winona's `raw_LiDAR_Data_README.rtf` and
`county_mosaic_LiDAR_Data_README.rtf`, and the MnGeo documentation directory listing (no
project or accuracy report is published there). The metadata's Final Deliverables lists a
**"Lidar Accuracy Assessment Report"** as *"one paper copy"*, which is the document most
likely to hold the number.

**This is the one quantity in this whole investigation that no amount of our own
processing can recover**, because it was applied to the delivered points before we ever
saw them.

**And gen2 has the same unknown, found independently while this was being written.**
`analysis/LIDAR_DOCUMENTATION_MINE.md` §1.2 quotes the 2021 Lidar Mapping Report p. 15:
*"Based on the statistical analysis, the lidar data was then adjusted to reduce the
vertical bias when compared to the survey ground control of higher accuracy."* Its value
is not published either. Two unpublished constants, one on each side of the DoD, is a
sharper statement of §7's verdict than anything measured here: the absolute level of
*neither* epoch is recoverable from the delivered data. That document also establishes
that the 143 LCPs **calibrated** gen2 while the NVA/VVA checkpoints were held out, so an
LCP tie is not independent of gen2's own calibration -- which is a further reason the
2008 control, held out from nothing we use, is the cleaner reference.

### Draft request, ready to send


> **Subject:** Request – 2008 SE Minnesota lidar: accuracy assessment report, QA/QC checkpoints, and the applied bias adjustment
>
> Dear Sean and MnGeo staff,
>
> I am using the 2008 Southeast Minnesota lidar (MN DNR / AeroMetric) together with the 2021 3DEP coverage to measure landscape change near Elba, in Winona County. Separating real change from a datum difference between the two epochs requires the original ground control, and I have got as far as the public record goes.
>
> The per-county validation reports on the MnGeo resources site have been extremely useful. I have all nine of them, and parsing the eight that belong to this project gives 1,004 control points whose per-county RMSE reproduces the figure printed in each report. Four things are not public, and I would like to request them:
>
> 1. **The Lidar Accuracy Assessment Report**, listed in the project metadata under Final Deliverables as "one paper copy." A scan is fine.
> 2. **The 127 AeroMetric QA/QC checkpoints**, as a shapefile, CSV, or table of coordinates and elevations. If they exist as a shapefile, the 2007 Pine County lidar checkpoint dataset you already publish is exactly the format I need.
> 3. **The value of the bias adjustment** described in process step 8 of the metadata: "Based on the QA/QC data, a bias adjustment was determined, and the results were applied to the lidar data." The magnitude is not published anywhere I can find, and no amount of processing on my end can recover it. I would also like to know **whether it was a single value for the whole project or determined separately per block, lift, or county** – the reason I ask is below.
> 4. **The vertical datum and geoid model of the MnDNR control points in the county validation reports.** The reports themselves state neither. The dataset metadata says the lidar is NAVD88 (Geoid03); I have assumed the control is on the same model, but I would rather not assume it, and I would also like the NAD83 realization (CORS96, NSRS2007, or HARN) the GPS positions were computed in.
>
> Two things you may want to know, both from parsing your own published tables:
>
> – **The per-county checkpoint counts in the metadata sum to 1,033, while the text says 1,009.** Also, the metadata gives Houston County as "0.110 m, 134 points", but the Houston validation report's own overall table has 109 points and prints an RMSE of 0.134 – the two numbers look transposed somewhere.
> – **The mean residual on open-terrain control varies a lot between counties** – Winona -101 mm and Olmsted -6 mm against Dodge +89 mm and Steele +79 mm, on the reports' own Control-minus-Surface convention. That county-to-county pattern is what makes me want to know whether the bias adjustment was applied once or per block.
>
> I am happy to cover copying or scanning costs, and to send you what I find – the 2008 control has turned out to be the best absolute reference available for this work.
>
> Thank you,
>
> Andy Wickert
> University of Minnesota

**Who to ask.** MN DNR lidar data steward **Sean Vaughn**, sean.vaughn@state.mn.us,
763-284-7223; MnGeo, gisinfo.mngeo@state.mn.us, 651-201-2499. *(Contacts carried from
`ADDITIONAL_GROUND_CONTROL.md` §5, which took them from a sub-agent fetch; I did not
re-verify them and they should be checked before sending.)*

---

## 10. What is still open

* **The NAD83 realization term.** Named in §7, unmeasured – PROJ 9.4.0 here has no
  operation for NAD83(CORS96)->NAD83(2011) ellipsoid heights. NGS's NADCON 5.0 does. It is
  the only candidate identified that would move the anchor and not the 2008 comparison,
  and it is worth one hour with the NGS tool.
* **The 2008 control's datum is still a dataset-level assertion.** The validation reports
  state no datum and no geoid; the GEOID03 linkage comes from `lidar_semn2008.html` alone.
  Request 4 of the email asks for it per-mark. If it is wrong, §7 changes completely.
* **Only four open marks are locally measurable.** 241 open-cover control points exist in
  the eight counties; 36 control points of all classes lie within 10 km of Elba. Screening
  the rest by §7.1 needs their tiles, one at a time. `ADDITIONAL_GROUND_CONTROL.md` §7.5
  step 5 already has this in the right order.
* **The per-county 190 mm spread in the delivered surface's residual** (§6) is
  unexplained and is a statewide-workflow question, not an Elba one.

---

## Appendix – verification status

| claim | how |
|---|---|
| `Error = Control Z - Surface Z` | arithmetic, 1 022/1 022 rows, max residual 0.000000 m; plus the physical class ordering in §2c |
| parsed tables reproduce the published RMSE | §2b, eight counties, each against the figure printed in its own report |
| 1 004 control points, eight counties | `parse_mndnr_2008_control.py`, the bundled CSV |
| 16 of them inside tiles on disk | tile header bounds vs control coordinates, printed by the run |
| every tie, ladder, spread and percentile in §3-§5 | `gen1_own_control_tie.py --ground csf`, one provenance run, ledger in `.trust/runs/` |
| the ground-source sensitivity | the same script re-run `--ground vendor` |
| gen2's -0.7 ± 11.5 mm at mark 2210 | read from `z_before_absolute.json`, `gen2_absolute_offsets`; produced by `analysis/groundtruth/gen2_checkpoint_tie.py`, not re-run here |
| the anchor's +22.7 ± 39.7 mm and its budget | read from `z_before_absolute.json`; not recomputed here |
| the NAD83 realization term | **attempted and failed** – PROJ returns a Ballpark/null operation. Named, not measured. |
| the bias adjustment value | **searched and not found** in any public source listed in §9 |
| MnGeo/DNR contact details | **not verified by me**; carried from `ADDITIONAL_GROUND_CONTROL.md` §5 |
| the ctrl_pct/tie_mm relation | measured, and flagged in §4 as partly tautological – it is not used as evidence |

# Null check of the half-sample mode on open control marks -- and the tie defect it caught

**Date:** 2026-10-16. **Script:** `analysis/control_mode_null_check.py`. **Output:**
`data/derived/control_mode_null_check.csv` (one row per mark, epoch, radius). Every number
below is pasted from the script's stdout or from the synthetic test recorded in
`analysis/hsm.py`.

## Why this exists

Andy, 2026-10-16: the control marks are not on the floodplain, so they cannot validate the
mode for vegetation. What they can do is **falsify the estimator**: on open, well-observed
ground the half-sample mode and the median target the same population and must agree. A
disagreement there is a defect in the estimator, and it would propagate into every result
that leans on the mode -- the coarsened-gen1 floodplain result and the gen2 band table in
`FRAME_2026-09-23.md` included. This is a null check, not a validation.

## Data paths

* **gen2**: the stored near-ground histograms at every control mark (`control_structure`,
  bins -1.0..2.0 m at 20 mm, radius 7.5 m), 519 marks of the 2021 set and 269 of the 2008
  set. The boxes at the 2008 marks are gen2 clouds too (LAS point format 7, GPS time
  3.03e8 = 2021), so there the surveyed level is GEOID03 against a GEOID18 cloud and only
  mode-minus-median is meaningful. At the 2021 marks both are GEOID18.
* **gen1**: the 88 marks of the 2008 set whose gen1 tile is on disk under `data/before`:
  29 L1O (open) plus 59 of other covers. Returns read from the tile, an order-2 local
  surface through vendor class-2 within R, slope-normal residuals of all returns but class
  7, histogrammed at the same 20 mm bins the floodplain work used. gen1 and its control
  share NAVD88(GEOID03). Two radii: 7.5 m (the gen2 box) and 2.5 m (one 5 m cell, where
  gen1 holds ~15 returns).

Columns: `hsm_med` = mode minus median; `med_surv`, `hsm_surv` = each minus surveyed ground;
`hsmraw_med` = the mode on the UNBINNED returns minus the median (gen1 only). mm; median
over marks with NMAD in parentheses.

## First run: the null FAILED, and the failure was the estimator

    gen2 at the 2021 marks        n marks   hsm_med        med_surv       hsm_surv
      NVA (open)                      227   -20.0 (29.7)    0.9 (47.7)    -24.3 (45.9)
    gen2 at the 2008 marks
      L1O (open)                       60   -20.0 ( 0.0)
      L2T                              88   -20.0 ( 0.0)
    gen1 at the 2008 marks, R 7.5 m        hsm_med        hsmraw_med
      L1O (open)                       29   -26.3 (14.5)   -6.2 (18.1)
    gen1 at the 2008 marks, R 2.5 m
      L1O (open)                       29    -8.9 ( 8.1)    0.0 ( 3.6)

The histogram mode sat exactly one 20 mm bin below the median at more than half the open
marks -- NMAD 0.0 -- while the mode on the unbinned returns agreed with the median. That
is not a property of the ground; it is the estimator.

**The cause.** The histogram form recurses on a sample of REPEATED bin centres. Many
half-windows then tie at the minimum width, and `np.argmin` returns the FIRST of them --
the lowest-starting window -- at every level of the recursion. Synthetic test (normal
samples, 20 mm bins, 400 draws each; mean of each estimator, mm):

    sigma   n     median   mode on raw   histogram mode, first-of-ties   middle-of-ties
      60   150      0.0        -0.2               -29.5                      +1.1
      60    16      1.3        +0.8               -15.7                      +7.5
     150   150     -0.6        -1.4               -58.1                      +5.8
      30  3000      0.0        -0.3               -18.5                      -0.5

The bias is downward and **grows with the spread of the distribution**. That is the very
signature the floodplain work was looking for in the data: a gen1 ground estimate that
drops further in the cells with the widest return spread. So the two results that used the
first-of-ties form on binned values had to be re-run.

**The fix.** `analysis/hsm.py`: break ties at the MIDDLE tied window. One implementation,
imported by the null check, `gen1_coarsened_ground.py` and `gen2_density_experiment.py`.
A residual of order +7 mm remains at n = 16 on binned values; the mode on raw values does
not have it, so prefer raw values when they are available.

## Second run: the null PASSES

    gen2 at the 2021 marks        n marks   n ret   hsm_med        med_surv       hsm_surv
      NVA (open)                      227    3179    0.0 ( 0.0)     0.9 (47.7)    -2.4 (48.8)
      VVA (vegetated)                 162    3232    0.0 (29.7)   113.1 (139.1)  103.0 (133.5)
      LCP                             130    4971  -20.0 (29.7)   -23.0 (48.7)   -38.6 (56.2)
    gen2 at the 2008 marks (level is gauge)
      L1O (open)                       60    1692    0.0 ( 0.0)
      L2T                              88    1712    0.0 ( 0.0)
      L3B                              44    1641    0.0 ( 0.0)
      L4F (forest)                     40    1604  -20.0 (29.7)
      L5U (urban)                      37    1650    0.0 (29.7)
    gen1 at the 2008 marks, R 7.5 m   n ret   hsm_med        hsmraw_med     med_surv       hsm_surv
      L1O (open)                 29     119    0.1 (12.4)    -6.2 (18.1)    18.5 (106.3)    9.3 (88.5)
      L2T                        16     144    1.2 ( 8.7)    -1.5 (11.6)    50.9 (115.5)   53.5 (124.2)
      L3B                        10     104    9.3 (34.9)     2.6 (36.1)   119.7 (83.4)   123.5 (104.2)
      L4F (forest)                6     104   -6.3 (16.2)   -16.6 (15.0)    30.4 (61.6)    24.0 (78.7)
      L5U (urban)                26     171    5.5 (15.8)    -5.9 (10.6)    -7.3 (76.3)    15.2 (73.2)
    gen1 at the 2008 marks, R 2.5 m
      L1O (open)                 29      15   -0.5 (11.8)     0.0 ( 3.6)     9.4 (69.6)    14.1 (88.8)
      L2T                        16      16    0.3 (10.3)    -0.3 ( 8.9)    71.4 (106.3)   58.7 (100.9)
      L3B                        10      12    4.1 ( 9.5)    -3.1 (15.3)    85.9 (104.2)   90.2 (106.2)
      L4F (forest)                6      12    9.0 (19.2)     2.8 (22.6)    25.8 (53.9)    46.5 (89.9)
      L5U (urban)                26      19    0.5 (13.9)     1.4 (10.8)    24.2 (66.0)    33.9 (58.2)

On open ground the mode equals the median (gen2 NVA: 0.0, NMAD 0.0; gen1 L1O: 0.1 at 7.5 m,
-0.5 at 2.5 m), and each sits the same distance from the surveyed elevation within the
mark-to-mark scatter. At one 5 m cell with ~15 gen1 returns the unbinned mode agrees with
the median to 0.0 mm (NMAD 3.6): **the small-n bias the resampling experiment measured in
heavy vegetation does not appear on open ground**, so it is a vegetation-and-count effect,
not a small-n property of the estimator.

What the marks do NOT show: at the vegetated marks (VVA, L4F) the mode sheds at most one bin
relative to the median. These marks are open enough to survey, so their vegetation tails
are thin; that says nothing about the floodplain sedge, which is the question the marks
cannot answer.

## Consequences: the two recorded results, re-run with the fixed estimator

**Coarsened gen1 on the flat floodplain** (`gen1_coarsened_ground.py`, 53,436 cells,
17 -> 144 returns per cell by 3x3 pooling). DoD by gen1-spread quartile, mm, and the
DoD-vs-spread regression coefficient, mm/mm:

    gen1 ground estimator               Q1      Q2      Q3       Q4    Q1-Q4   coefficient
    current: median of own cell       10.2    -6.0   -45.8   -119.2    129.4     -0.360
    pooled 3x3: median                 5.7   -13.3   -49.3   -101.1    106.8     -0.277
    pooled 3x3: mode, first-of-ties   34.8    32.0    12.4    -10.5     45.3     -0.187   (FRAME 2026-09-23; WRONG)
    pooled 3x3: mode, fixed            7.1    -6.4   -30.2    -58.5     65.6     -0.226

The mode still does real work beyond the pooled median -- Q4 rises from -101.1 to -58.5 and
the quartile range falls from 106.8 to 65.6 mm -- but the artifact's structure collapses by
**49%**, not the 65% the frame recorded, and the coefficient is -0.226, not -0.187. The old
+21 mm level shift was the tie bias; the fixed mode's level is -14.7 against the median's
-23.5. The remaining 65.6 mm of structure is the measured fixable-versus-maskable boundary.

**gen2 full-density estimators across near-ground vegetation bands**
(`gen2_density_experiment.py --population highveg`, flat floodplain, >= 64 returns):

    lowveg band       cells     q50     hsm (old)   hsm (fixed)   q50-hsm (old)   q50-hsm (fixed)
    0.0000-0.0000     6,122     6.9       -1.9          7.7            8.8             -0.8
    0.0000-0.0237     2,334     9.7       -5.6          8.7           15.3              1.0
    0.0237-0.1019     4,219    15.9       -6.0         11.6           21.8              4.3
    0.1019-9.0000     4,225    39.2       -6.0         16.5           45.2             22.7

This is a cleaner result than the old one, not a weaker one: with no near-ground vegetation
the mode now equals the median (-0.8 mm), which is the correct null, and the separation is
entirely vegetation-driven, 0 -> 22.7 mm. The old table showed 8.8 mm of "separation" with
zero vegetation -- that was the tie bias. The mode moves 8.8 mm across the bands where the
median moves 32.3 (the frame said 4.1). Resampling in the high-vegetation band: the mode's
bias against full density is +11.5 mm at n = 16 (was +19.7) and +3.4 at n = 64 (was +7.0);
per-cell scatter 56.8 and 29.1 mm (unchanged).

## Reproduce

    ./lidar-icp/bin/python analysis/control_mode_null_check.py
    ./lidar-icp/bin/python analysis/gen1_coarsened_ground.py
    ./lidar-icp/bin/python analysis/gen2_density_experiment.py --population highveg

## Open

* `gen1_coarsened_ground.py` pools histograms with `uniform_filter * 9` in float32 and the
  mode floors the counts to int; the median does not. Both forms were used before and after
  the fix, so the before/after comparison is like for like, but the floor should go.
* At the 29 open 2008 marks gen1's median sits +18.5 mm (NMAD 106) ABOVE the surveyed
  ground at 7.5 m; the adopted delivered-surface offset at Elba's open marks is +62.7 mm
  reads-LOW from 8 marks assigned by flight line. Different population, different
  estimator, both within their scatter; noted, not chased.

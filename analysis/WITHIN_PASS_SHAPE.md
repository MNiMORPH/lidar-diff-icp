# Shape of the within-pass horizontal offset, gen1 at elbaext

Read-only measurement. No pipeline code was changed, no tile was rebuilt, and no point
cloud was opened -- everything below comes from
`data/derived/elbaext_delong/beam_offset_table.parquet` (five columns) and the grids
named in "Inputs", and nothing else.

Scripts are throwaway and live in the session scratchpad
(`prep.py`, `common.py`, `t0_baseline.py`, `t0b_checks.py`, `t1_shape.py`,
`t1_shape_boot.py`, `t2_boot_sens.py`, `t4_terrain.py`, `t4b_space.py`,
`final_slopes.py`, plus two ad-hoc one-liners for the semivariogram in Section 1 and the
correction-surface projection in Section 5); every number here is pasted from their
stdout, and every number inside a code fence was checked back against the captured
stdout programmatically.

---

## Inputs, and exactly what each supplies

| artifact | role |
|---|---|
| `beam_offset_table.parquet`, columns `cell`, `point_source_id`, `gps_time`, `d_mm_corr`, `in_grid` | the response and its time/space/swath labels. 15,581,591 rows read; 13,928,563 are `in_grid` with finite `d_mm_corr` |
| `stable.npy` | the cell mask. 451,869 of 720,900 cells |
| `z_after_differenced.npy` | the gen2 reference whose **gradient** the fit regresses on. Nearest-filled, then Gaussian-smoothed at sigma 1.2 cells, then `np.gradient(..., 5 m)` -> `(dz/dnorth, dz/deast)` |
| `z_after.npy` | **only** for `nnorm = sqrt(gx^2+gy^2+1)`, to turn the table's slope-normal mm into vertical mm. This is the grid the beam table's own producer used for `nnorm` (`gen1_save_angles_slope.py` line 56-59), so the conversion is exactly the inverse of the one that made the column |
| `applied_correction_surface.npy` | used **only** as a control (Section 5), never as an input to the primary fit |

Grid read from `corrections.json`: NY=810, NX=890, res 5.0 m, origin (575600.0, 4882200.0),
`cell = iy*NX + ix`, rows increasing northward.

### One thing worth flagging about the response column

`d_mm_corr` is **slope-normal** millimetres -- `gen1_save_angles_slope.py` divides by
`nnorm`. The gradient regression is a statement about **vertical** differences, so every
value is multiplied back by `nnorm` at its cell before fitting. Skipping that would
shrink the fitted shifts by `cos(slope)` in a slope-dependent way, i.e. it would
manufacture exactly the terrain dependence T4 is trying to test for.

`d_mm_corr` carries four registration terms (`geoid`, `lateral`, `swath`, `drift`). On
this product `along_track_drift = False` and `geoid_applied = False`, so `dz_drift_mm`
is identically zero and there is no geoid term; the lateral cross-epoch shift
(-0.7519, -0.1896 m) and the per-swath internal alignment are in. It does **not** carry
the DeLong correction surface -- that is applied to the gridded product downstream, not
to the per-return table. Section 5 adds it back as a control.

### The estimator

For a set of (cell, swath) units,

    d_vertical_mm  =  dx * dz/deast  +  dy * dz/dnorth  +  c

by ordinary least squares; `dx`, `dy` come out in millimetres directly (the response is
in mm and the gradient is dimensionless -- unlike `starframe.fit_star_shifts`, which
fits metres and scales by 1000). Sign convention is `starframe`'s: the fitted
coefficient is MINUS the applied shift, i.e. it is the shift **ready to add** to the
gen1 coordinates.

### Choices I made, stated because they are choices

1. **`>= 3` returns per (cell, swath)** -- handed to me in the brief. Sensitivity in
   Section 8: it moves nothing.
2. **8 blocks per pass**, equal-COUNT in `gps_time` (not equal-width), so every block
   carries the same number of cell-swath units. Sensitivity in Section 8.
3. **Terrain strata for T4**: 10 slope quantile bins (computed per pass) x 8 aspect
   sectors of 45 deg, with aspect-undefined (flat) cells as their own stratum. I picked
   10 x 8 to get strata that are still populated in a single time block; I did not tune
   it. If you want a different stratification, that is a one-line change.
4. **Superblock size for the spatial bootstrap**: 50 cells = 250 m, with 100 cells =
   500 m run alongside as a convergence check. Both are reported because neither is
   clearly right -- see the caveat in Section 3.

No minimum cell count, no slope cut, no aspect cut, no outlier rejection was applied
anywhere. Where a block is thin, its bootstrap SE says so.

---

## 1. Baseline: the whole-pass fit, and what the noise floor actually is

663,207 (cell, swath) units survive `stable & n>=3`.

```
 swath    cells          t0          t1   dur_s    dx_mm    SE    dy_mm    SE     c_mm    SE     R2
   133   33,913    249658.5    249708.3    49.7    590.6   4.1    116.4   3.2    -88.2   0.5  0.392
   134  155,049    248598.6    248650.1    51.5    437.6   1.8    -13.5   1.6    -96.7   0.2  0.281
   135  171,660    246671.4    246723.1    51.7    175.5   1.6    -55.9   1.3    -95.4   0.2  0.072
   136  153,551    245641.6    245694.6    53.0    -49.9   1.5     28.4   1.4    -78.0   0.3  0.011
   137  113,572    243650.7    243701.9    51.2   -333.5   1.7     48.1   1.4    -82.1   0.3  0.297
   138   35,462    237352.5    237404.6    52.1   -478.3   2.6    -29.4   2.6    -70.1   0.7  0.510
```

Each pass's `[t0, t1]` window brackets the corresponding mid gps_time given in the
brief. Two small departures from the brief's framing, noted rather than smoothed over:
the passes are **49.7 to 53.0 s** long, not ~55 s, and the grid is **4050 m** tall
(bounds 4882200-4886250 in `corrections.json`), not 4350 m. Flight direction,
**measured** rather than assumed, from `d(north)/d(gps_time)` per pass:

```
  swath 133 (S): -80.94 m/s     swath 136 (N): +78.67 m/s
  swath 134 (N): +81.25 m/s     swath 137 (S): -81.59 m/s
  swath 135 (S): -80.81 m/s     swath 138 (N): +77.80 m/s
```

### The handed noise floor is a random-split floor, and it is not the right scale

Reproducing the handed floor the obvious way -- random half-split of the cell-swath
units, 20 draws x 6 swaths = 120 splits:

```
half-split dx: RMS(halfA-halfB) =   8.50 mm  -> per-half RMS  6.01 mm
half-split dy: RMS(halfA-halfB) =   6.70 mm  -> per-half RMS  4.74 mm
half-split dz: RMS(halfA-halfB) =   0.82 mm  -> per-half RMS  0.58 mm
```

Close to the handed dx 4.7 / dy 3.9 / dz 0.6 (mine are a little larger; the two
conventions -- half the difference, or the difference over sqrt(2) -- straddle the
handed values, so I take this as agreement rather than a discrepancy).

**But a random half-split is the wrong floor for a BLOCK.** It splits interleaved cells,
so both halves sample the same terrain patches and the split measures only within-patch
return noise. The residual after the three-parameter fit is spatially correlated --
measured here, not assumed. Empirical semivariogram of the post-fit residual over
4,000,000 random cell-swath pairs per swath (mm^2):

```
 swath  sill(var)      25m      75m     150m     300m     600m    1200m    2400m
   134       7917     4723     5472     6398     6946     6974     7920    11901
   135       8453     5564     6383     7167     7565     7695     8525    10709
   136      11052     8209     8421     9865    10777    10517    10837    13066
   137      12175     9014     9295    10236    10814    11383    12790    15476
```

gamma is still below the sill at 300 m and reaches it between 600 and 1200 m (beyond
that it overshoots, which is the tile-scale trend, not correlation). So the practical
correlation range is several hundred metres. A time block here is ~500 m of northing by
roughly one swath width -- **of order one correlation length**. Section 3 measures the
consequence: the honest per-block uncertainty is 20-100 mm, not 5 mm.

This is the first thing in the brief that the measurement contradicts. "Anything much
above [4-5 mm] is real" is true of within-patch noise and false of a block-to-block
difference.

### Reproducing the handed 4-blocks-by-northing slopes

Same estimator, 4 equal-count blocks by northing, weighted linear slope in mm per block:

```
 swath  dx_mm/block     SE  dy_mm/block     SE        handed dx    handed dy
   133         16.0    4.0        -11.8    3.0            +15.0         -0.7
   134          2.6    1.6        -34.3    1.4            +13.2        -32.1
   135          6.5    1.5         22.8    1.2             +8.3        +30.1
   136          9.4    1.6        -16.4    1.2            +12.0        -10.5
   137        -12.0    1.6         -0.9    1.2             -9.1         +6.7
   138         -8.6    2.8        -38.0    2.4             -5.0        -20.7
```

Same sign on **6 of 6** for dx and **5 of 6** for dy. The one dy disagreement is swath
137, where the handed value is +6.7 and mine -0.9 -- both within 7 mm of zero. I treat
the handed numbers as reproduced.

---

## 2. T1 -- within-pass shape: the block curves

8 equal-count blocks by `gps_time` per pass. `t-tmid` is seconds from the pass midpoint,
`north` the mean northing of the block's cells (it runs the opposite way in the two
headings, which is the whole point).

```
swath 133 (S)  [4239 units per block]        swath 134 (N)  [19381 per block]
 blk  t-tmid   north      dx     dy           blk  t-tmid   north      dx     dy
   0   -21.7 4885926   720.2   89.2             0   -20.3 4882445   407.0  -37.4
   1   -14.5 4885307   581.5  111.0             1   -13.7 4882987   405.4  156.8
   2    -8.4 4884825   659.8  261.2             2    -8.4 4883437   562.7   52.0
   3    -3.0 4884400   455.5  -30.9             3    -2.7 4883883   546.2  -29.3
   4     2.5 4883969   651.5   68.4             4     2.5 4884293   448.4  -83.8
   5     9.2 4883418   533.3   89.3             5     7.6 4884703   485.3   51.3
   6    15.4 4882924   543.2  156.0             6    13.3 4885188   506.2   -4.9
   7    20.6 4882470   735.8  139.0             7    21.7 4885873   439.2  -57.8

swath 135 (S)  [~21457 per block]            swath 136 (N)  [19194 per block]
 blk  t-tmid   north      dx     dy           blk  t-tmid   north      dx     dy
   0   -22.7 4885972   148.5   -5.3             0   -21.8 4882551   -75.2  160.0
   1   -15.0 4885347   255.1  -37.1             1   -15.6 4883039  -143.7   24.2
   2    -7.6 4884764   191.1  -49.8             2    -8.9 4883551   -56.6   43.4
   3    -2.0 4884311   169.5    0.1             3    -3.2 4884005    61.6  -30.9
   4     3.6 4883847   190.6 -122.7             4     2.8 4884482   -55.3    8.2
   5     8.9 4883417   198.6 -308.7             5     8.7 4884943   -96.1   43.6
   6    14.3 4882996    33.0 -173.8             6    15.6 4885504   -17.2   -9.6
   7    20.5 4882480   207.4  -35.8             7    22.3 4886010    23.5  -56.7

swath 137 (S)  [~14196 per block]            swath 138 (N)  [~4433 per block]
 blk  t-tmid   north      dx     dy           blk  t-tmid   north      dx     dy
   0   -17.3 4885998  -344.0    5.3             0   -20.8 4882643  -499.9   71.1
   1   -11.8 4885527  -404.5   24.9             1   -16.2 4883019  -498.9  -61.3
   2    -7.3 4885151  -297.4   66.6             2    -6.5 4883754  -470.1  -91.7
   3    -3.4 4884818  -343.0   96.5             3     0.4 4884278  -558.4 -114.3
   4     0.5 4884517  -489.3   35.1             4     4.7 4884624  -573.9  -42.2
   5     5.5 4884115  -289.0   15.5             5     8.4 4884911  -487.5   49.9
   6    11.6 4883631  -239.7  102.1             6    12.3 4885218  -475.3 -131.9
   7    22.1 4882763  -372.6   42.5             7    17.8 4885668  -576.0 -120.7
```

Block-to-block peak-to-peak, mm: dx 105.9 (138) to 280.3 (133); dy 96.9 (137) to
308.8 (135). So the offset is emphatically **not** constant along a pass -- that part of
the brief is confirmed. The question is what the variation is a function of.

---

## 3. The uncertainty that decides everything: OLS SE vs spatial block bootstrap

The OLS standard error of a block fit treats each 5 m cell as independent. It is not:
after the three-parameter fit the residual is a spatially correlated field. Resampling
250 m superblocks (50x50 cells) with replacement inside each time block, 200 draws:

```
 swath  superblocks/blk  OLS SE dx  boot SE dx  ratio  OLS SE dy  boot SE dy  ratio
   133              6.1      14.65      101.63    6.9      13.37       76.31    5.7
   134             16.5       5.82       51.36    8.8       5.47       55.87   10.2
   135             18.6       5.29       49.24    9.3       4.63       41.53    9.0
   136             19.4       4.68       53.76   11.5       4.30       28.81    6.7
   137             19.4       5.18       38.95    7.5       4.44       29.05    6.5
   138              9.4      10.68       34.34    3.2      10.77       37.67    3.5
```

**The OLS SE understates the block uncertainty by 3.2x to 11.5x.** Honest per-block
uncertainties are 20-100 mm, against block-to-block spreads of 100-300 mm. Every
significance statement below uses the bootstrap SE.

Stated as a limitation, not hidden: 250 m is SHORTER than the correlation range just
measured (gamma still climbing at 300 m),
so even this bootstrap is a **lower bound** on the uncertainty. I tried 500 m
superblocks and it cannot be done -- a time block then contains 1 to 3 superblocks
(swath 133: `[3, 2, 2, 2, 2, 3, 2, 1]`), the bootstrap degenerates (one block returned
SE exactly 0.0, another 101470 mm) and the weighted fit fails outright. That failure is
itself the finding: **at the correlation length of the residual field, a 6.5 s time
block contains only one or two independent terrain samples.** There is no way to make a
within-pass block estimate of a lateral shift here that is better-determined than a few
tens of millimetres, whatever the return count.

---

## 4. T1 -- the shape verdict: line, quadratic, or neither

For each pass and each component I fit a constant, a line and a quadratic in within-pass
time, weighting the 8 blocks by their bootstrap SEs, and report chi-square with its
degrees of freedom and p-value. The classification rule -- **my choice, stated so you
can overrule it**:

* **no resolved shape** -- a CONSTANT already fits (chi-square p > 0.05): the block
  scatter is within the honest uncertainty, nothing to explain;
* **linear** -- the constant is rejected, the LINE fits, and the quadratic term adds
  nothing;
* **curved** -- the line is rejected and the quadratic fits with a significant
  second-order term;
* **erratic** -- the constant is rejected and NEITHER the line nor the quadratic fits.

```
swath 133 (S)  dx: ptp  280.3 mm | const chi2  10.95/7 p=0.141 | linear slope  -3.236+-1.452 mm/s (t=-2.23), chi2   5.98/6 p=0.425 | quad c2 +230.887+-199.214 mm/s^2 x1e-3 (t=+1.16), chi2   4.64/5 p=0.462
swath 133 (S)  dy: ptp  292.1 mm | const chi2   7.29/7 p=0.399 | linear slope  +0.890+-0.936 mm/s (t=+0.95), chi2   6.39/6 p=0.381 | quad c2 +112.263+-132.493 mm/s^2 x1e-3 (t=+0.85), chi2   5.67/5 p=0.340
swath 134 (N)  dx: ptp  157.3 mm | const chi2  10.82/7 p=0.147 | linear slope  +0.954+-1.249 mm/s (t=+0.76), chi2  10.24/6 p=0.115 | quad c2 -190.124+-96.253 mm/s^2 x1e-3 (t=-1.98), chi2   6.34/5 p=0.275
swath 134 (N)  dy: ptp  240.6 mm | const chi2  16.99/7 p=0.017 | linear slope  -1.870+-1.145 mm/s (t=-1.63), chi2  14.33/6 p=0.026 | quad c2 -107.161+-97.224 mm/s^2 x1e-3 (t=-1.10), chi2  13.11/5 p=0.022
swath 135 (S)  dx: ptp  222.1 mm | const chi2  10.42/7 p=0.166 | linear slope  -0.221+-1.052 mm/s (t=-0.21), chi2  10.37/6 p=0.110 | quad c2   -8.909+-79.550 mm/s^2 x1e-3 (t=-0.11), chi2  10.36/5 p=0.066
swath 135 (S)  dy: ptp  308.8 mm | const chi2  28.32/7 p=0.000 | linear slope  -2.875+-1.096 mm/s (t=-2.62), chi2  21.44/6 p=0.002 | quad c2   +0.617+-82.732 mm/s^2 x1e-3 (t=+0.01), chi2  21.44/5 p=0.001
swath 136 (N)  dx: ptp  205.4 mm | const chi2  14.17/7 p=0.048 | linear slope  +2.424+-1.087 mm/s (t=+2.23), chi2   9.20/6 p=0.162 | quad c2  +17.570+-88.906 mm/s^2 x1e-3 (t=+0.20), chi2   9.16/5 p=0.103
swath 136 (N)  dy: ptp  216.8 mm | const chi2  24.48/7 p=0.001 | linear slope  -2.201+-0.596 mm/s (t=-3.69), chi2  10.84/6 p=0.093 | quad c2  -17.381+-57.085 mm/s^2 x1e-3 (t=-0.30), chi2  10.75/5 p=0.057
swath 137 (S)  dx: ptp  249.6 mm | const chi2  27.15/7 p=0.000 | linear slope  +1.022+-1.220 mm/s (t=+0.84), chi2  26.45/6 p=0.000 | quad c2  -58.093+-99.984 mm/s^2 x1e-3 (t=-0.58), chi2  26.11/5 p=0.000
swath 137 (S)  dy: ptp   96.9 mm | const chi2  11.31/7 p=0.126 | linear slope  +0.799+-0.873 mm/s (t=+0.91), chi2  10.48/6 p=0.106 | quad c2 -102.144+-67.478 mm/s^2 x1e-3 (t=-1.51), chi2   8.19/5 p=0.146
swath 138 (N)  dx: ptp  105.9 mm | const chi2  18.77/7 p=0.009 | linear slope  +0.373+-0.868 mm/s (t=+0.43), chi2  18.59/6 p=0.005 | quad c2 +129.210+-93.793 mm/s^2 x1e-3 (t=+1.38), chi2  16.69/5 p=0.005
swath 138 (N)  dy: ptp  203.0 mm | const chi2  56.73/7 p=0.000 | linear slope  -3.230+-0.742 mm/s (t=-4.35), chi2  37.79/6 p=0.000 | quad c2  +59.532+-91.281 mm/s^2 x1e-3 (t=+0.65), chi2  37.37/5 p=0.000
```
(Pasted verbatim from `t1_shape_boot.py`, reordered so the two components of a pass sit
together. `ptp` is the block peak-to-peak; the quadratic coefficient is in mm/s^2 x 1e-3.)

### Verdict per pass

| pass | heading | **dy (along-track)** | **dx (cross-track)** |
|---|---|---|---|
| 138 | N | **erratic**, with a strong linear part (-3.23 +- 0.74 mm/s) that does not exhaust it | **erratic** |
| 137 | S | **no resolved shape** (a constant fits) | **erratic** |
| 136 | N | **linear**, -2.20 +- 0.60 mm/s | **linear**, +2.42 +- 1.09 mm/s (constant only marginally rejected) |
| 135 | S | **erratic**, with a linear part (-2.88 +- 1.10 mm/s); driven by one block at -308.7 mm | **no resolved shape** |
| 134 | N | **erratic** (constant rejected p=0.017, neither line nor quadratic fits) | **no resolved shape** |
| 133 | S | **no resolved shape** (a constant fits) | **no resolved shape** -- a constant already fits (p=0.141), so the t=-2.23 on the line is a term the data do not need |

**No pass anywhere shows a supported quadratic.** The largest second-order t-statistic
in dy is -1.51 (137) and in dx is -1.98 (134); every other one is below 1.4. The
"oscillating error on a minutes timescale" hypothesis is **not supported by any pass**.

One clean linear pass (136 dy), two passes with no resolved within-pass shape at all
(133, 137 dy), three erratic (134, 135, 138 dy). That is the honest answer to T1: the
within-pass shape is mostly **erratic**, and where a linear part is present it does not
account for the block scatter.

---

## 5. T4 -- terrain-matched blocks

### How the match was made

Within each pass, cells are stratified by **10 slope quantile bins** (quantiles of that
pass's own stable cells, so the bins are populated) **x 8 aspect sectors of 45 degrees**,
with aspect-undefined (flat) cells forming their own stratum -- 90 strata. The target
histogram is the whole pass's own joint (slope, aspect) distribution. A unit in time
block k gets weight `target_frac(stratum) / block_k_frac(stratum)`, renormalised to mean
1 within the block, and the fit is weighted least squares with those weights. Strata
absent from a block cannot be filled; the `cov%` column is the share of the target
histogram that the block can reach, so the quality of the match is visible.

I also ran a second control in parallel: adding back the **DeLong correction surface**
(`applied_correction_surface.npy`), which the pipeline adds to the gen1 elevations but
which is NOT in `d_mm_corr`. Regressed on the same gradient basis over stable cells it
is worth only `dx = -5.5 +- 0.4 mm`, `dy = +5.0 +- 0.3 mm` tile-wide (R2 = 0.0011,
n = 448,550; the surface itself is median 93.5 mm, IQR 50.0 mm, peak-to-peak 1378.3 mm),
so it is a level field, not a disguised lateral shift -- but it is spatially structured
and therefore a candidate for manufacturing a within-pass trend.

### The terrain statistics per block, so the match is checkable

`dx_raw`/`dy_raw` are the unweighted fits from Section 2; `dx_mat`/`dy_mat` the
terrain-matched ones; `dx_+cs`/`dy_+cs` the ones with the correction surface added back.

```
swath 133 (S)  blk       n mean_slope  mean|gE|  mean|gN|  mean gE  mean gN   cov%   dx_raw   dx_mat   dy_raw   dy_mat   dx_+cs   dy_+cs
                 0   4,239      11.79     0.121     0.162   -0.032   -0.024  100.0    720.2    720.1     89.2    110.9    791.3    108.1
                 1   4,239       9.34     0.112     0.086   -0.025    0.009  100.0    581.5    570.4    111.0     66.4    628.4    119.6
                 2   4,239       4.62     0.049     0.052   -0.019   -0.035   70.3    659.8    645.8    261.2    273.0    674.3    301.6
                 3   4,239       5.15     0.054     0.052   -0.010    0.032   82.7    455.5    447.5    -30.9    120.0    250.7     91.3
                 4   4,239       8.88     0.112     0.082    0.012    0.025  100.0    651.5    647.5     68.4    122.3    655.7     37.5
                 5   4,239      10.94     0.109     0.130    0.057   -0.032   87.9    533.3    564.7     89.3    118.1    534.0     89.3
                 6   4,239      13.52     0.107     0.197   -0.057    0.176   92.2    543.2    543.4    156.0    165.4    534.8    140.6
                 7   4,240       4.60     0.029     0.066    0.002   -0.030   96.3    735.8    797.5    139.0    234.8    685.3    193.9

swath 134 (N)  blk       n mean_slope  mean|gE|  mean|gN|  mean gE  mean gN   cov%   dx_raw   dx_mat   dy_raw   dy_mat   dx_+cs   dy_+cs
                 0  19,381       6.39     0.061     0.084    0.021    0.052  100.0    407.0    450.3    -37.4     15.6    417.4    -44.1
                 1  19,381       9.91     0.092     0.131   -0.005    0.103   98.0    405.4    428.4    156.8    100.8    394.6    234.7
                 2  19,381       6.25     0.056     0.076    0.020   -0.018  100.0    562.7    566.5     52.0    134.8    630.7     32.1
                 3  19,381       4.94     0.056     0.049   -0.000    0.015  100.0    546.2    588.7    -29.3     45.6    535.3    -44.1
                 4  19,381       4.42     0.046     0.047   -0.013    0.018   91.2    448.4    533.6    -83.8     -6.5    316.5   -183.3
                 5  19,381       5.80     0.062     0.061   -0.019   -0.021  100.0    485.3    586.2     51.3    -88.2    517.4    104.0
                 6  19,381       8.60     0.087     0.104   -0.014   -0.049  100.0    506.2    520.4     -4.9      7.4    597.9     52.4
                 7  19,382      15.18     0.172     0.180   -0.037    0.033  100.0    439.2    467.3    -57.8    -64.7    463.6    -75.7

swath 135 (S)  blk       n mean_slope  mean|gE|  mean|gN|  mean gE  mean gN   cov%   dx_raw   dx_mat   dy_raw   dy_mat   dx_+cs   dy_+cs
                 0  21,458      12.34     0.150     0.136   -0.010    0.080  100.0    148.5    134.9     -5.3    -10.9    175.5    -98.8
                 1  21,457      12.06     0.097     0.170    0.005   -0.063  100.0    255.1    281.5    -37.1    -11.7    297.1    -25.5
                 2  21,455       8.42     0.090     0.099    0.001   -0.022  100.0    191.1    173.0    -49.8    -88.0    233.7    -39.0
                 3  21,460       6.71     0.075     0.072   -0.022    0.014  100.0    169.5    102.5      0.1     25.6    117.9    -24.3
                 4  21,457       6.37     0.055     0.084   -0.007   -0.017  100.0    190.6    177.8   -122.7    -70.0    281.7   -147.4
                 5  21,458       4.61     0.043     0.054   -0.006   -0.008   89.0    198.6    267.2   -308.7   -218.5    219.6   -364.9
                 6  21,451       6.40     0.065     0.073   -0.006    0.005  100.0     33.0    122.0   -173.8   -121.2     29.3   -241.1
                 7  21,464      12.95     0.110     0.185    0.017    0.176   98.8    207.4    216.8    -35.8    -92.9    216.8    -42.1

swath 136 (N)  blk       n mean_slope  mean|gE|  mean|gN|  mean gE  mean gN   cov%   dx_raw   dx_mat   dy_raw   dy_mat   dx_+cs   dy_+cs
                 0  19,194      14.48     0.156     0.177   -0.066    0.152   97.6    -75.2    -25.9    160.0    202.6    -12.1    101.8
                 1  19,194      10.32     0.091     0.137   -0.046   -0.025  100.0   -143.7   -137.2     24.2     15.5   -120.4     52.9
                 2  19,194       9.79     0.102     0.119   -0.016    0.024   99.0    -56.6    -40.1     43.4     43.8    -85.7    100.4
                 3  19,193      10.35     0.106     0.124   -0.029    0.014  100.0     61.6     69.9    -30.9    -36.4     -4.1      5.9
                 4  19,194       8.58     0.106     0.095    0.007    0.006   99.3    -55.3    -31.0      8.2     30.6    -76.4    -13.8
                 5  19,194      12.71     0.169     0.122    0.017   -0.098  100.0    -96.1    -55.5     43.6     37.3   -121.3    -64.5
                 6  19,194      11.37     0.094     0.160   -0.014   -0.026  100.0    -17.2    -36.4     -9.6    -37.7    -58.9    -41.8
                 7  19,194      12.95     0.120     0.178   -0.006    0.104  100.0     23.5    -27.4    -56.7    -46.4     -8.5     36.6

swath 137 (S)  blk       n mean_slope  mean|gE|  mean|gN|  mean gE  mean gN   cov%   dx_raw   dx_mat   dy_raw   dy_mat   dx_+cs   dy_+cs
                 0  14,197      15.77     0.123     0.243   -0.064    0.051   99.2   -344.0   -373.0      5.3      1.8   -436.6     52.3
                 1  14,195      13.62     0.107     0.206    0.034   -0.078  100.0   -404.5   -423.5     24.9     -2.4   -425.1     24.0
                 2  14,198      14.06     0.162     0.156   -0.022   -0.133   97.5   -297.4   -313.4     66.6     43.6   -299.2    -14.7
                 3  14,196       8.72     0.100     0.092   -0.045    0.007  100.0   -343.0   -344.5     96.5     53.2   -333.2     16.7
                 4  14,196       8.53     0.089     0.105   -0.025    0.042  100.0   -489.3   -509.1     35.1     42.2   -499.9     84.8
                 5  14,197      16.37     0.161     0.222   -0.042    0.010  100.0   -289.0   -308.9     15.5      9.2   -266.1     16.4
                 6  14,196      16.05     0.202     0.190   -0.063    0.103  100.0   -239.7   -241.9    102.1     81.3   -257.5    142.7
                 7  14,197      20.70     0.252     0.236   -0.084    0.019  100.0   -372.6   -373.5     42.5     44.8   -302.7     92.2

swath 138 (N)  blk       n mean_slope  mean|gE|  mean|gN|  mean gE  mean gN   cov%   dx_raw   dx_mat   dy_raw   dy_mat   dx_+cs   dy_+cs
                 0   4,433      20.98     0.200     0.292    0.120    0.197  100.0   -499.9   -459.5     71.1     69.4   -464.3    100.7
                 1   4,433      19.33     0.287     0.168    0.248   -0.121   63.3   -498.9   -434.7    -61.3     35.5   -481.1    -89.5
                 2   4,432      20.83     0.256     0.277   -0.143    0.042   87.6   -470.1   -517.8    -91.7   -117.3   -460.6    -70.9
                 3   4,433      20.31     0.255     0.251    0.018   -0.118   96.2   -558.4   -654.4   -114.3   -212.2   -537.4   -118.4
                 4   4,432      15.09     0.176     0.169   -0.136    0.078   98.2   -573.9   -574.8    -42.2    -54.2   -589.3    -32.8
                 5   4,433      13.64     0.199     0.089   -0.191    0.061   83.8   -487.5   -602.3     49.9    115.4   -497.3     50.2
                 6   4,433      13.42     0.175     0.144   -0.103   -0.057   81.7   -475.3   -487.8   -131.9   -122.7   -479.1   -131.2
                 7   4,433      22.06     0.189     0.333   -0.024   -0.210   89.5   -576.0   -584.6   -120.7   -151.3   -635.9   -125.2
```

The terrain really does differ block to block: mean slope ranges 4.42-15.18 deg within
swath 134 alone, and `mean gN` swings from -0.210 to +0.197 within swath 138. That is
exactly the confound worth removing. Coverage is exactly 100% in 27 of the 48 blocks,
falls below 90% in nine, and below 80% in only two: 138 block 1 at 63.3% and 133 block 2
at 70.3% -- both in the two thinnest swaths.

---

## 6. A second terrain test the geometry hands us for free: space-locked or time-locked?

The six passes overlap and the heading alternates, so the **same northing band is early
in a north-bound pass and late in a south-bound one**. Fitting the shift in COMMON
northing bands (8 equal-width bands of 506 m spanning 4882200-4886245 m N) therefore
separates the two hypotheses directly:

* an error that is a function of **place** -- terrain, or any unmodelled field of
  position -- gives the same band profile in every pass, whichever way the plane flew;
* an error that is a function of **within-pass time** gives band profiles that reverse
  between headings.

Pairwise Pearson r of the band profile (pass mean removed), 15 pairs:

```
dy:  SAME-heading pairs (n=6): mean r = +0.312   median +0.375
     OPPOSITE-heading pairs (n=9): mean r = +0.092   median +0.121
dx:  SAME-heading pairs (n=6): mean r = +0.061   median +0.157
     OPPOSITE-heading pairs (n=9): mean r = +0.142   median +0.315
```

Individual pairs, across both components, run from -0.734 to +0.893. With 8 bands, the sampling SE of a
correlation under the null is about 1/sqrt(5) = 0.45, which is larger than every mean
above. **This test does not discriminate.** It is reported because it was run, and
because "same-heading correlates better than opposite-heading, for dy" is the kind of
0.3-vs-0.1 difference that would be easy to sell as evidence and is not.

### T4 result: the along-track trend survives terrain matching

Per-pass linear trend of the block values against within-pass time, blocks weighted by
their 250 m spatial-block-bootstrap SEs (150 draws), averaged over three bootstrap
seeds; `seed sd` is the spread of the trend across those three seeds, i.e. how much of
the number is bootstrap noise.

```
--- raw ---
 swath dir   t(h)    dy/dt     SE  seed sd      t    dx/dt     SE  seed sd      t
   138   N  0.000   -2.971  0.697    0.310  -4.26    0.231  0.825    0.197   0.28
   137   S  1.749    0.852  0.806    0.088   1.06    1.248  1.174    0.235   1.06
   136   N  2.303   -2.257  0.593    0.089  -3.81    2.441  1.115    0.050   2.19
   135   S  2.589   -2.863  1.060    0.128  -2.70   -0.307  1.061    0.184  -0.29
   134   N  3.124   -2.149  1.130    0.248  -1.90    0.828  1.240    0.079   0.67
   133   S  3.418    1.024  0.859    0.053   1.19   -3.119  1.366    0.540  -2.28

--- terrain-matched ---
   138   N  0.000   -4.459  1.632    0.042  -2.73   -2.430  1.371    0.361  -1.77
   137   S  1.749    1.294  0.859    0.025   1.51    0.998  1.164    0.149   0.86
   136   N  2.303   -1.942  0.716    0.074  -2.71    1.044  1.464    0.135   0.71
   135   S  2.589   -2.142  1.003    0.058  -2.13    1.009  1.196    0.053   0.84
   134   N  3.124   -3.259  1.176    0.090  -2.77    0.314  1.403    0.281   0.22
   133   S  3.418    2.176  1.006    0.322   2.16   -2.301  1.570    0.392  -1.47

--- with the DeLong correction surface added back ---
   138   N  0.000   -1.955  0.982    0.403  -1.99   -1.042  0.940    0.195  -1.11
   137   S  1.749    2.385  0.724    0.048   3.29    4.018  1.310    0.144   3.07
   136   N  2.303   -1.727  0.498    0.013  -3.47    0.537  1.008    0.188   0.53
   135   S  2.589   -3.751  0.973    0.272  -3.85    0.365  1.011    0.118   0.36
   134   N  3.124   -1.464  0.755    0.212  -1.94    1.142  1.185    0.518   0.96
   133   S  3.418    2.185  0.580    0.057   3.77   -2.393  1.611    0.136  -1.49
```

**Measured.** For **dy**, all six passes keep the sign of their trend under both
controls -- terrain matching and the correction surface. Magnitudes move by up to
1.5 mm/s (the two largest moves: 138 from -2.971 to -4.459 under matching, 137 from
+0.852 to +2.385 with the surface), which is about one standard error, but nothing flips
and nothing collapses to zero. Group means:
raw N `-2.459 (sd 0.447)`, matched N `-3.220 (sd 1.259)`, surface-controlled N
`-1.715 (sd 0.246)`.

For **dx** the picture is the opposite: signs do NOT survive. 138 goes `+0.231 ->
-2.430`, 135 `-0.307 -> +1.009`, 136 `+2.441 -> +1.044`. The cross-track trend is a
function of which terrain the block happened to contain.

**Inferred.** The along-track (dy) within-pass trend is not an artifact of block terrain
composition, and it is not the missing DeLong correction surface. The cross-track (dx)
within-pass trend is not established as anything beyond terrain and block noise.

---

## 7. T2 -- does |slope| track pass time within a direction group?

Using the raw variant above.

```
N-bound: pass times [0.000 2.303 3.124]  |dy/dt| [2.971 2.257 2.149]  SEs [0.697 0.593 1.130]
   weighted slope of |dy/dt| vs pass time = -0.289 +- 0.343 mm/s per hour  (t=-0.84, n=3)
   first vs last pass: 2.971 - 2.149 = +0.822 +- 1.328  (t=+0.62)
S-bound: pass times [1.749 2.589 3.418]  |dy/dt| [0.852 2.863 1.024]  SEs [0.806 1.060 0.859]
   weighted slope of |dy/dt| vs pass time = +0.141 +- 0.705 mm/s per hour  (t=+0.20, n=3)
   first vs last pass: 0.852 - 1.024 = -0.172 +- 1.178  (t=-0.15)
```

**Measured.** The north-bound group's |slope| decreases monotonically with pass time
(2.971 -> 2.257 -> 2.149 mm/s over 3.1 hours), but the decrease is -0.289 +- 0.343
mm/s per hour: **not resolved**. The south-bound group is non-monotonic (0.852 ->
2.863 -> 1.024) and its trend is +0.141 +- 0.705, also not resolved, and with the
opposite sign to the north group's.

**This is n=3 per group and it is underpowered.** I am not computing a p-value on it
and I am not leaning on the monotonicity of the north group: three points in order
happens one time in six by chance, and the one group that could contradict it does.

### The heading contrast, and why a per-pass slope is not a discrimination

```
dy/dt: N mean -2.459 (sd 0.447, n=3) vs S mean -0.329 (sd 2.196, n=3)
  difference -2.130 +- 1.294 (Welch t=-1.65)
  S within-group sd 2.196  vs  |difference| 2.130
```

The **south group's internal scatter is as large as the north-south difference**. In
the space frame the same thing holds: `dy/dNorth` N mean -0.0311 (sd 0.0062) vs S mean
+0.0041 (sd 0.0271), difference 0.0352 against an S scatter of 0.0271. So neither
"heading-dependent" nor "space-locked" is established, and a single pass's slope carries
no discrimination at all -- it is one draw from a distribution whose width is the size
of the effect being tested.

---

## 8. Sensitivity to the two choices

Block count K (OLS-weighted trend, mm/s, so that the numbers are comparable across K
without bootstrap noise):

```
             K=4       K=6       K=8      K=10      K=12
   133     0.926     0.781     0.749    -0.126    -0.176
   134    -3.118    -2.643    -2.254    -1.887    -1.743
   135    -2.053    -2.228    -2.348    -1.775    -2.730
   136    -1.512    -3.071    -3.035    -2.702    -2.578
   137     0.185     0.516     1.190     0.893     0.758
   138    -3.330    -3.492    -3.553    -3.089    -3.163
```

The three north-bound passes stay negative at every K. 133 changes sign between K=8 and
K=10 and 137 varies by a factor of six -- the two passes T1 already classified as having
no resolved within-pass shape. The K-spread on 134 (-3.12 to -1.74) is comparable to its
bootstrap SE (1.13), so K is a real but not dominant source of variation.

Minimum returns per (cell, swath) unit:

```
min_n= 1 (668,275 units): 133: +0.690 134: -2.347 135: -2.423 136: -3.086 137: +1.238 138: -3.519
min_n= 2 (665,891 units): 133: +0.747 134: -2.301 135: -2.394 136: -3.066 137: +1.210 138: -3.537
min_n= 3 (663,207 units): 133: +0.749 134: -2.254 135: -2.348 136: -3.035 137: +1.190 138: -3.553
min_n= 5 (657,515 units): 133: +0.679 134: -2.229 135: -2.326 136: -2.937 137: +1.121 138: -3.593
min_n=10 (630,329 units): 133: -0.040 134: -2.074 135: -1.974 136: -2.821 137: +0.786 138: -3.410
```

This choice is immaterial: 1 through 10 returns changes nothing except 133, which is
consistent with zero throughout.

---

## 9. What I measured, and what I infer

### Measured

1. The within-pass horizontal offset is **not constant**: block-to-block peak-to-peak is
   97-309 mm in dy and 106-280 mm in dx.
2. The honest per-block uncertainty is **20-100 mm**, 3.2x to 11.5x the OLS standard
   error, because the fit residual is spatially correlated over several hundred metres
   (variogram, Section 1) and a time block holds only one or two independent terrain
   samples. The handed 4.7/3.9/0.6 mm floor is
   a random-half-split floor (I reproduce it: 6.0/4.7/0.6 mm per half) and it is not the
   right scale against which to call a block difference real.
3. Against that scale: **no pass shows a supported quadratic** in either component
   (largest |t| on the second-order term is 1.98). The "oscillating on a minutes
   timescale" shape is absent.
4. dy shape by pass: 136 **linear**; 134, 135, 138 **erratic**; 133, 137 **no resolved
   shape**. dx shape by pass: 136 linear (marginally); 137, 138 erratic; 133, 134, 135
   no resolved shape.
5. The dy trend **survives terrain matching** (10 slope quantiles x 8 aspect sectors,
   post-stratification reweighting) and **survives** adding the DeLong correction
   surface back: all six passes keep their sign under both. The dx trend does not
   survive either control.
6. |dy/dt| declines monotonically across the three north-bound passes (2.971 -> 2.257 ->
   2.149 mm/s over 3.1 h) but the decline is -0.289 +- 0.343 mm/s per hour, not
   resolved; the south-bound group is non-monotonic with the opposite-signed,
   equally-unresolved trend.
7. The north-south contrast in dy/dt is -2.130 +- 1.294 mm/s, and the south group's own
   internal scatter (sd 2.196) is as large as that contrast.
8. The common-northing-band profile test does not discriminate space-locked from
   time-locked: same-heading mean r = +0.312 vs opposite-heading +0.092 for dy, against
   a null sampling SE of about 0.45.

### Inferred (hypotheses, not results)

* The along-track component carries something real that terrain composition does not
  explain, of order 2-3 mm/s within a pass (100-150 mm over a 52 s pass). Its shape is
  mostly **not** a clean drift: in three of six passes a straight line leaves
  block-to-block structure that neither it nor a quadratic accounts for.
* The residual that defeats the line is itself consistent with an unmodelled field --
  terrain-correlated error that my slope-and-aspect stratification does not capture
  (curvature, cover, incidence angle are all uncontrolled here), or a real
  spatially-correlated difference between the epochs. This analysis cannot separate
  those.
* Nothing here supports a per-pass, heading-dependent instrument error over a
  time-varying GNSS error, or the reverse. Both groups are n=3 and the south group's
  scatter swallows the contrast.

### What would actually decide it

The limit is not return count, it is **independent terrain samples per unit time**: at
the residual correlation range measured in Section 1 a 6.5 s block contains one or two. More passes
over the same ground, or a longer pass, buys discrimination; more returns does not. If
you want the heading-vs-time question settled, the measurement to make is on a tile
where the same flight line is flown twice in opposite directions, or across enough
tiles that the per-heading group has n much greater than 3.

### Contradictions of the numbers I was handed

1. **The noise floor.** dx 4.7 / dy 3.9 / dz 0.6 mm is reproduced as a random-half-split
   floor, but it is 3-11x too small for a block estimate. "Anything much above that is
   real" does not hold: the honest floor is 20-100 mm per block.
2. **Swath 137's along-track sign.** The handed 4-block-by-northing slope is +6.7
   mm/block; I get -0.9 +- 1.2 mm/block with the same blocking. Both are consistent with
   zero, and T1 classifies 137 dy as having no resolved within-pass shape -- so the sign
   of this one should not be used.
3. **Swath 134's cross-track slope.** Handed +13.2 mm/block by northing; I get
   +2.6 +- 1.6. Same sign, one fifth the size.
4. **The dy slopes generally are larger than the handed 4-block ones** when measured in
   time with 8 blocks -- e.g. 138's K=8 OLS-weighted slope of -3.553 mm/s (Section 8) is
   -185.2 mm over its 52.1 s pass, against the handed -20.7 mm/block, which over the
   three block-steps of a 4-block fit is -62 mm across the tile. Different blocking, and
   part of it is that a 4-block fit in northing averages over structure that an 8-block
   fit in time resolves. Not a contradiction, but the two sets of numbers are not
   interchangeable.
5. **Two framing details.** The passes are 49.7-53.0 s long, not ~55 s, and the grid is
   4050 m tall, not 4350 m. Neither changes anything, but both are in the brief.

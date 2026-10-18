# Pre-test: the statistical structure of gen1's floodplain returns

**Date:** 2026-10-16. **Script:** `analysis/floodplain_return_pdf.py`. **Figures:**
`figures/floodplain_return_pdf_{A_scales,B_csf_vs_all,C_geometry}.png`. Every number is
pasted from the script's stdout.

**Reproduce:**

    ./lidar-icp/bin/python analysis/floodplain_return_pdf.py
    ./lidar-icp/bin/python analysis/floodplain_return_pdf.py --all-table data/derived/elbaext_star_allreturns

The all-returns table is built by the product's own two producers on the raw merged cloud
in a scratch tile of symlinks (`data/derived/elbaext_star_allreturns`, every input a
symlink to `elbaext_delong_star_nosurf`):

    ./lidar-icp/bin/python analysis/modules/vegetation_correction/gen1_save_angles_slope.py elbaext_star_allreturns data/before/elbaext_gen1_merged.laz
    ./lidar-icp/bin/python analysis/modules/vegetation_correction/beam_offset_table.py data/derived/elbaext_star_allreturns data/before/elbaext_gen1_merged.laz

## Why (Andy, 2026-10-16)

After 3x3 pooling with the fixed half-sample mode, ~66 mm of DoD structure across
gen1-spread quartiles remains on the flat floodplain; it is neither the estimator's count
bias nor the absence of returns at the gen2 surface. Before building a lowest-cluster
estimator or settling for the mask, look at the PDF itself: does a ground peak exist below
the mat at all, does the cloth's point selection hide one, and does beam geometry sort the
returns?

Population: 53,590 flat floodplain cells (floodplain mask AND terrain slope <= 2 deg,
>= 3 CSF returns), gen1-spread quartile edges 140.2 / 209.2 / 306.6 mm; 1,060,995 CSF
returns, 17 per cell. The residual is `d_mm_corr`: each return's slope-normal height above
the gen2 surface, registration applied (positive = gen1 above gen2).

Parameters that are mine, swept: peak-finding smoothing 20 and 40 mm; second-peak
prominence 10% and 20% of the main peak; |scan angle| bins of 5 deg.

## A. PDF over growing areas -- UNIMODAL at every width

Per-cell statistics, median over cells (mm; "below" = share of returns at or below the gen2
surface; "2nd" = share of cells with a qualifying local maximum BELOW the main peak, at
40 mm smoothing and 10% prominence, the gentler of the four settings):

    quartile     scale   n/cell   mode    med   below   2nd peak
    Q1 <140        1x1      15    -10    -10    0.63      1.1%
    Q1            15x15   3417    -10    -10    0.54      1.9%
    Q2 140-209     1x1      22     10    -10    0.50      9.2%
    Q2            15x15   4837    -10     10    0.47      1.1%
    Q3 209-307     1x1      19     30     30    0.33     29.0%
    Q3            15x15   4347     10     30    0.40      1.0%
    Q4 >307        1x1      18     90    110    0.21     51.3%
    Q4             3x3     147     50     90    0.25     19.8%
    Q4             5x5     404     30     90    0.27     10.9%
    Q4             9x9    1336     30     70    0.29      6.0%
    Q4            15x15   3854     30     70    0.30      3.5%

The full table with all four peak-finding settings is in the script output. The pattern
is the same at every setting: **the share of cells with a second peak FALLS with pooling
width**, at every quartile, toward 1-4%. Secondary peaks are small-sample wiggles. No
ground peak emerges beneath the mat as the count grows; the quartile-average PDF
(`figures/floodplain_return_pdf_A_scales.png`) is a single peak at every width.

What changes across quartiles is (1) the single peak shifts from about -10 mm to about
+30 mm, (2) an upper tail grows that carries the median to +70-110 mm, and (3) the share
of returns at or below the gen2 surface falls from 0.54 to 0.30. From 3x3 outward the mode
is stable (Q4: 50 / 30 / 30 / 30), so pooling past 3x3 is not finding anything new.

## B. CSF-selected vs all returns -- IDENTICAL on the flat floodplain

Identity check: 100.000% of CSF returns matched in the all-returns table; max
|d_mm_corr difference| 0.000e+00 mm. The all-returns table is the product's own producer on
the superset cloud.

    quartile     scale   table   n/cell   mode    med   below   2nd (40 mm, 10%)
    Q1 <140        1x1    CSF       15    -10    -10    0.63      1.1%
    Q1             1x1    all       15    -10    -10    0.63      1.1%
    Q4 >307        1x1    CSF       18     90    110    0.21     51.3%
    Q4             1x1    all       18     90    110    0.21     51.3%
    Q4             3x3    CSF      147     50     90    0.25     19.8%
    Q4             3x3    all      147     50     90    0.25     19.8%

Every row of the table is identical to the digit (Q2, Q3 included). Within the +/-1.5 m
window the cloth discards NOTHING on the flat floodplain: 17 returns per cell both ways.
The 1.77 M points CSF removes from this tile are elsewhere (canopy, buildings). So the
floodplain's return population IS the whole population, and there is nothing lower for a
different point selection to recover.

## C. Beam geometry -- scan angle: nothing; intensity: the tail only

Population PDF of CSF returns per quartile (mm):

    Q4 >307      bin             n        mode    med   below
                 |scan| 0-5    115,543     50    110    0.28
                 |scan| 5-10    70,782     50    130    0.25
                 |scan| 10-15   79,036     50    130    0.25
                 |scan| 15-90    9,229     70    130    0.27
                 intensity low  77,054     30     90    0.33
                 intensity mid  77,850     30     90    0.28
                 intensity high 119,686    70    170    0.21

Scan angle does not sort the returns: the four bins lie on top of each other in every
quartile (`figures/floodplain_return_pdf_C_geometry.png`). Near-nadir beams do not fall
through the sedge. Intensity separates only the UPPER TAIL in Q3 and Q4: the high-intensity
tercile's median sits 80 mm above the low tercile's in Q4 while all three terciles share the
same peak. Intensity is a weak tail classifier, not a ground classifier.

## What this decides -- measurement, then the hypothesis separately

**Measured.** The PDF of gen1's flat-floodplain returns is unimodal at every pooling width
from 5 to 75 m; the cloth's point selection and the whole-cloud population are the same
population here; beam geometry does not sort the returns. Across gen1-spread quartiles the
single peak shifts ~40 mm upward and an upper tail grows; the mode captures the peak and
is stable from 3x3 pooling outward.

**Consequence for the estimator question.** A lowest-cluster estimator has no target: there
is no lower cluster. The per-epoch mode is the right estimator for what gen1 can see, and
what it cannot see -- the ~40 mm peak shift plus the tail in the most-vegetated cells -- is
not recoverable from gen1's returns by any point selection, pooling width or beam filter
tested. That leaves the mask (the span LoD, gen1's own p90-p10, no fitted number) for those
cells.

**Hypothesis, not measured.** The peak that gen1 sees in Q4 sits ~40 mm above the one it
sees in Q1 relative to the gen2 surface. That is consistent with gen1 returning from a
dense litter or sedge mat flush with the ground (a surface, not a tail), and partly with
gen2's own reference surface reading high in those same cells (its median sits 23 mm above
its mode in the most-vegetated band). Neither is separable from this data alone.

## D. Shape vs pooling width: skew, half-widths, tail form, and where gen2 sits (2026-10-17)

Andy: "How skewed is the PDF at each pooling width? Do we see differences in its shape? Can
we characterize its tail? At what percentile would the gen2 surface land within the gen1
cloud at each size of cell?" Script: `analysis/floodplain_return_shape.py`; figure
`figures/floodplain_return_shape_logpdf.png` (log axis). Pasted from stdout.

**Per cell (median over cells).** Bowley skewness is quartile-based and holds at 17 returns;
moment skewness needs the tails, which 17 returns cannot sample, so it climbs with pooling
for count reasons alone (Q1: 0.07 -> 1.37) and is not a shape statement at small n.

    quartile    scale   n/cell   bowley   pct0   p10   p50   p90    sd
    Q1 <140       1x1      15     0.00     63    -70   -10    50    44
    Q1          15x15    3417     0.00     54   -110   -10   130   125
    Q2 140-209    1x1      22     0.00     50    -90   -10   110    73
    Q2          15x15    4837     0.11     47   -110    10   210   151
    Q3 209-307    1x1      19     0.11     33    -70    30   190   104
    Q3          15x15    4347     0.11     40   -110    30   270   170
    Q4 >307       1x1      18     0.14     21    -50   110   390   166
    Q4            3x3     147     0.12     25    -70    90   370   178
    Q4            5x5     404     0.13     27    -70    90   350   180
    Q4          15x15    3854     0.14     30    -90    70   350   191

**Quartile-average PDF** (cells weighted equally; half-widths at 40 mm smoothing; upper
tail fitted over mode+60..mode+1000 mm as exponential, e-folding length L, and as power
law, exponent a; R^2 in log-density space):

    quartile    scale   mode   bowley   pct0   hw-   hw+   L exp   R2     a pow   R2
    Q1 <140       3x3    -10    0.00     56    100   100    139   0.953   -2.86  0.980
    Q1          15x15    -10    0.33     52    100   100    171   0.955   -2.31  0.978
    Q2 140-209    3x3    -10    0.14     46    120   140    133   0.987   -2.85  0.932
    Q2          15x15    -10    0.25     44    120   160    161   0.990   -2.36  0.931
    Q3 209-307    3x3     10    0.00     37    140   200    140   0.996   -2.63  0.878
    Q3          15x15     10    0.20     39    140   180    165   0.996   -2.24  0.895
    Q4 >307       1x1     50    0.20     27    180   300    196   0.983   -1.82  0.821
    Q4            3x3     30    0.14     30    160   240    189   0.991   -1.92  0.851
    Q4            5x5     30    0.08     31    160   220    186   0.993   -1.95  0.861
    Q4          15x15     10    0.17     34    140   220    183   0.995   -1.99  0.866

**Skew.** Modest and set by vegetation, not by width: Bowley 0.00 in Q1 and Q2, 0.11-0.14 in
Q3 and Q4, the same at every pooling width per cell. The distribution is right-skewed only
through its tail; its central quartiles are nearly symmetric.

**Shape vs width.** The core does not change with pooling: mode, Bowley skew and both
half-widths are the same from 3x3 to 15x15 in every quartile (Q4 half-widths 160/240 ->
140/220 mm). What pooling does is fill in the tails that 17 returns cannot sample -- on the
log axis the 1x1 curve drops away early on both sides and each wider pool lifts both limbs.
In Q4 even that is small: its curves nearly coincide at all widths, because 18 returns per
cell already reach a tail that extends 400 mm. The shape difference that matters is ACROSS
quartiles: the peak widens on both sides with vegetation (hw- 100 -> 160, hw+ 100 -> 240 mm
at 3x3) and shifts up ~40 mm; the lower side widening says the sub-surface part of gen1's
cloud is broader in vegetated cells too.

**Tail.** EXPONENTIAL, in every quartile, at every width from 3x3 up: the upper limb is a
straight line on the log axis, R^2 0.95-0.996 against 0.85-0.98 for a power law (the power
law competes only in Q1, whose tail is short). The e-folding length is ~140 mm in Q1-Q3 and
~185 mm in Q4, and in Q4 it is stable across widths (196 / 189 / 186 / 183). This is the
canopy-type tail the forest work found (an exponential, Beer-Lambert-like attenuation
through scatterers), not the power-law tail of bare and cropped ground: standing sedge
behaves like a low canopy.

**Where gen2 sits in the gen1 cloud.** Per cell, the gen2 surface is at the 63rd percentile
of gen1's returns in the least-vegetated quartile and the 21st in the most-vegetated at
5 m, converging to the 54th and 30th as pooling widens (the small-n percentiles are biased
toward the extremes by unsampled tails). In Q1 gen1's median sits 10 mm BELOW gen2; in Q4
seven of ten gen1 returns sit ABOVE the gen2 surface.

## E. The transfer function pct0(span), by pooling window (2026-10-18)

Andy: "What I really want is a continuous function representing that transfer. Please
generate that relationship at a range of different spatial window sizes and see which one
is the best." Scripts: `analysis/floodplain_gen2_percentile.py` (binned curve, figure
`figures/floodplain_gen2_percentile_vs_span.png`) and
`analysis/floodplain_transfer_function.py` (fit; `figures/floodplain_transfer_function.png`;
parameters in `data/derived/floodplain_transfer_function.json`).

**Form.** pct0(span) = p_lo + (p_hi - p_lo) / (1 + exp((span - s0)/w)): p_hi the percentile
at zero vegetation, p_lo under the densest sedge, s0 the span at the midpoint of the fall,
w the transition width. Robust least squares on all 53,590 cells; span and pct0 both read
from the same k x k pooled cloud.

**Criterion for "best" (mine, a proposal):** out-of-sample skill under 5-fold CV with folds
by 500 m block, so no held-out cell shares pooled returns with a training cell.

    window   returns   CV R2   NMAD (pct pts)   p_hi   p_lo    s0     w
    1x1         17     0.130       31.5         78.7   21.1   182    74
    3x3        143     0.147       24.5         63.0   24.8   267    49
    5x5        399     0.124       22.8         60.0   27.7   284    43
    7x7        788     0.108       21.9         58.4   29.8   290    38
    9x9       1331     0.098       21.3         57.2   31.4   294    34
    15x15     3916     0.079       20.7         54.5   34.2   308    23

(Own-cell span as x gives the same ranking; s0 then sits at ~175 mm at every window because
the x axis no longer pools.)

**Best window: 3x3**, by the criterion -- the only window whose CV R^2 beats 1x1, with the
per-cell sampling noise of 17 returns removed (NMAD 31.5 -> 24.5) before spatial smoothing
starts blurring the response (R^2 falls from 5x5 on as pooling drives every cell's pct0
toward the neighbourhood mean). Its function: 63.0 at zero span, 24.8 at high span, midpoint
267 mm, width 49 mm.

**The honest number is the R^2.** The binned median curve is clean, but span explains only
13-15% of cell-to-cell variance in where gen2 sits; the 25-75% band is +/-25 percentile
points at every span. At the same vegetation span, cells differ widely in where the gen2
surface falls in gen1's cloud. The function is a statement about the MEDIAN cell, not a
per-cell predictor, and any per-cell use of it carries that scatter.

**What it is not (stated once, in the script too).** Applied to choose gen1's ground
quantile per cell, it would set gen1's ground to gen2's surface at every span, and the
floodplain DoD would be zero by construction -- the cross-epoch calibration that
STRATEGY_2026-09-23.md argues against, and the forest-vs-open degeneracy again. Its use is a
decision for Andy.

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

## F. A more appropriate form? Tried, and the logistic stays (2026-10-18)

Andy asked for a more appropriate, perhaps elegant, fitting function. The physically
motivated candidate: since the PDF is a symmetric core plus an exponential tail (section
D), let span and pct0 both be readouts of one latent tail fraction f in the mixture
(1-f) Normal(mu, sigma) + f ExpNormal(mu, sigma, L) -- three physical parameters against
the logistic's four. Script `analysis/floodplain_transfer_mixture.py`, same block CV.

    window   mixture CV R2   logistic CV R2
    1x1          0.114           0.130
    3x3          0.105           0.147
    5x5          0.063           0.124

It loses, and the reason is the useful result: a fixed core width cannot reach the data.
The fitted core (sigma ~70 mm) alone spans 180 mm, but bare cells span 100; the mixture
saturates at span 436, but vegetated cells reach 660. The core widens with vegetation --
per cell at 3x3, p50-p10 runs 60 -> 180 mm from the least- to the most-vegetated quartile
-- so span is set by two latent quantities, core width and tail fraction. No one-parameter
family through (span, pct0) can follow that, which is also why span explains only ~14% of
per-cell variance under any form. The logistic's four parameters remain the appropriate
description of the median tendency.

## G. Is the plateau real? Three tests; no global-parameter, no-plateau form survives (2026-10-18)

Andy doubted the logistic's floor: denser vegetation should keep blocking the beam, so pct0
should keep falling, and a form with global parameters only should follow. Script
`analysis/floodplain_transfer_plateau.py` (three sections). All 3x3, block CV as before.

**1. Not the window.** Widening the histogram window from +/-1.5 to +/-5 m recovers 0.74% of
returns and changes nothing; fewer than 0.3% of cells below the top bin have p90 > 1.5 m.

**2. No-plateau forms in span lose.**

    form                          CV R2   NMAD   behaviour at span 600 / 900
    logistic (4)                  0.153   24.5   24.0 / 23.9
    exp-to-floor (4)              0.153   24.5   23.2 / 22.3
    power law (3)                 0.135   24.8   21.3 / 13.7
    stretched exponential (3)     0.126   24.8   19.8 / 10.8
    threshold exponential (3)     0.115   24.4   16.2 /  6.6

Between spans of 440 and 580 mm the binned pct0 is flat at 20-23; above 600 mm (4.3% of
cells) it RISES to 30-42, a different population, not denser sedge. Span cannot register
density past ~570 mm for an algebraic reason: for a core plus exponential tail, p90 - p10
saturates near core + L ln(10) as the tail fraction approaches 1. The proxy plateaus; the
physics need not.

**3. A proxy that keeps registering, and the null that kills it.** Beer-Lambert
p_hi exp(-x/h0) in the p90 HEIGHT scores CV R2 0.656 with two global parameters and keeps
falling at the top (12.9 -> 8.1). But under the algebraic null, every cell given the same
bare-floodplain shape and shifted only to reproduce its own pct0, p90 height gives R2 1.000
and the tail integral from gen2 gives 0.99. Both are measured relative to the gen2 surface,
so they read gen2's LEVEL in gen1's cloud: they are pct0 restated. The translation-invariant
candidates (p90 - mode; tail length as mean excess above the mode; tail mass above the mode;
their product) score 0.05-0.13, below span, and none keeps falling in its top decile.

**What this means.** gen1's own cloud contains no vegetation-density measure that is both
independent of the gen2 level and able to register past span's saturation. Beer-Lambert
P = exp(-k tau) is the right physics, but tau must come from outside gen1's returns, and
the cloud's own P is what pct0 already measures (pct0 ~ 0.63 P, from the bare-ground 0.63).
The logistic in span stays, read correctly: its floor is the floor of the observed range,
not a statement that denser vegetation would stop blocking the beam.

## H. Is the plateau a mat height? The core does not decompose (2026-10-18)

Andy: is the floor telling us the height of a low vegetation component, so that gen1 is
hitting a mat rather than bare earth? Script `analysis/floodplain_mat_decomposition.py`.

Test: fix a ground component's shape from the barest span decile (Gaussian, centre -17.6,
width 49.5 mm against gen2) and, per span decile, fit ground (free weight) + a free mat
Gaussian + an exponential tail above the mat.

    span decile   span   obs pct0   w_ground   core width (1-comp)   R2 2-comp   R2 1-comp
    0 (barest)     120     62.5       0.72          50                1.000       1.000
    2              200     56.9       0.11          82                1.000       0.999
    5              300     38.8       0.07         105                1.000       0.999
    8              460     21.8       0.07         134                0.997       0.994

Not identifiable. From the second decile up the ground weight falls to 0.07-0.09 and is
unconstrained, and a single broad core plus tail fits as well. The dense-cell core is ONE
hump, ~130 mm wide against 50 on bare ground, with no resolvable narrow ground peak under
it. A sharp ground surface carrying 20% of returns beneath a mat would show as a shoulder;
none does.

**Consistent with the hypothesis, not a confirmation of it (hypothesis, not measured):** the
surface gen1 sees in dense sedge is a single broad one whose peak sits ~30-50 mm above the
gen2 surface (section A), and gen2's own surface there sits ~20 mm above its mode (gen2
density experiment), so that surface is of order 50-70 mm above where gen2's mode puts
the ground, with a 130 mm spread -- a mat or litter layer with no sharp base, inside which
the gen2 surface lies. The floor of ~25% is then the lower flank of that hump below gen2's
level, not a count of bare-earth hits. The data cannot separate "few ground hits under a
mat" from "a mat with no base", and no independent ground truth exists on this floodplain.

## I. The correction, applied (2026-10-18)

Andy accepted the mat hypothesis as likely correct and asked for the 3x3 logistic to be
applied: erosion side only, tapered, with success judged by whether the vegetation
"erosion" patches fade into the DoD's background. Script
`analysis/floodplain_vegcorr_apply.py`; outputs next to the product:
`dod_floodplain_vegcorr.npy`, `floodplain_vegcorr_delta.npy`, `floodplain_vegcorr.json`;
figures `floodplain_vegcorr_maps.png`, `floodplain_vegcorr_vs_span.png`.

Stated once: this is a cross-epoch calibration. The mean floodplain DoD at a given span goes
to ~0 by construction and real erosion in vegetated cells would be absorbed. Numbers are
site-fitted; the method travels, the numbers do not.

**Mechanics.** Flat floodplain (mask AND slope <= 2 deg), 54,372 cells. gen1's ground becomes
the logistic-predicted percentile (refitted on these histograms: 63.1 / 24.5 / 268 / 50) of
the 3x3 pooled cloud; delta = own-cell median minus that percentile. Only delta > 0 is
applied, weighted by clip(-DoD_raw / T, 0, 1) with T = 55.5 mm, the NMAD of the DoD on
stable ground (measured). 42% of flat cells receive a correction; median 60 mm, p90 220.

**Before / after, flat floodplain (mm):**

    DoD by pooled-span quartile       Q1      Q2      Q3       Q4    coeff   NMAD   median
    raw                             14.3    -0.6   -45.1   -111.6   -0.353   89.3    -18.8
    corrected (tapered)             14.4     2.5    -5.5    -13.4   -0.134   55.7      0.7
    full correction, no taper       16.5    10.5     0.5    -11.7   -0.108   64.9      6.8
    stable ground (background)                                               55.5      5.0

    cells beyond the LoD            below -LoD          largest patch   above +LoD
    stable ground (background)       2.8%                   359 cells     4.5%
    floodplain raw                  28.6%                 3,405           7.3%
    floodplain corrected            13.1%                 1,027           8.0%

The span signature is gone: quartile medians flatten to within +/-14 mm. The tapered
floodplain NMAD (55.7) matching the background (55.5) is PARTLY mechanical -- T sets the
scale of what stays negative -- so the taper-free row (64.9 vs 55.5) is the honest
measure: close to background, not at it. Erosion patches halve in count and in largest
size but remain ~5x the background rate (13.1% vs 2.8%; largest 1,027 vs 359 cells).
Partial success, measured.

**Two defects, not fixed, for Andy:**
1. **Overshoot.** The taper is not capped: 382 cells (0.7%) are lifted from erosion to above
   +LoD, so "above +LoD" rose 3,998 -> 4,380. They sit where the pooled cloud is
   multi-surface (823 cells get delta > 300 mm; their own-cell span median is 567 mm against
   209): the logistic's percentile jumps across a bimodal cloud. A cap at zero, or at the
   LoD, is one line -- but it is a parameter.
2. **The slope cut is a hard cutoff**, which Andy asked to avoid. The floodplain mask holds
   96,716 cells; only the 61% at <= 2 deg were corrected. The 36,254 steeper cells (valley
   walls and banks inside the elevation cut) are untouched and are most of the dark red in
   the corrected map: 19.9% of them sit below -LoD. Tapering the correction in slope, or
   extending the population, is a decision.

## J. The patches that survive the correction: the NE patch is not vegetation (2026-10-18)

Andy: "The correction is in the right direction. There are still significant 'erosional'
patches, especially in the NE corner." Script `analysis/floodplain_ne_patch.py` (every number
below is its stdout); figures `floodplain_ne_patches.png`, `floodplain_ne_intensity.png`.

**The remaining patches are flat, so the slope cut is not the reason.** The eight largest
are 64-99% flat cells (two exceptions are the channel itself, at -0.9 to -1.7 m). What they
share is penetration failure: 44-87% of their cells have NO gen1 return at the gen2 surface
(own cell), so the logistic's 25th percentile still sits above gen2 and the correction is
bounded by what the cloud holds (delta 0-80 mm against raw DoD -130 to -260).

**The largest, patch 1 (NE corner; 2,355 cells = 5.9 ha; E 579925 N 4885765; 72% flat,
median slope 0.9 deg; raw DoD -246, corrected -179 mm):**

    gen1 return heights above gen2 (mm)     p5    p25    p50    p75    p95   share <= 0   intensity
      patch 1                                 8    136    245    368    568        4%          27
      rest of flat floodplain               -156    -49     20    121    390       43%          19

gen1's ENTIRE cloud sits ~25 cm above gen2's surface -- not a tail above a core but the whole
distribution shifted, with the same 99% single-return, 81% vendor class-2 character as the
floodplain. Same flight line (138), same pass, same near-nadir scan angle as the ring of
cells around it: not a flight artifact. gen1's intensity is 27 inside against 19-21 outside,
and the boundary is SHARP: across 20 m the raw DoD goes -216 -> -56 and intensity 27 -> 19.

gen2 (May 2021) sees clean, open ground there: near-ground vegetation share 0.030 against a
floodplain median of 0.18, class-2 share 0.81, p10/p50/p90 of -50/+10/+50 mm. No water
class in either epoch. NAIP 2010 NDVI 0.23 (sparse) against 0.19 for the rest of the
floodplain; the patch sits in the LOWEST zone of the floodplain (gen2 ground 217.9 m against
~223) at the tile's east edge beside the channel.

**Reading, as hypothesis.** This is not the sedge-mat signature: the vegetated-floodplain
signature is a tail above a core with gen2 also vegetated, and this is a whole-cloud offset
with a sharp edge, brighter in 2008, and open ground in 2021. In November 2008 gen1 returned
from a surface 25 cm above the 2021 ground over a sharply bounded, low-lying, bright area:
a frozen or snow-covered backwater/pond, or a bright sand sheet (the Whitewater's August
2007 flood deposited sand on these floodplains) that was gone by 2021. Which of these cannot
be settled from this data. Either way it is a SURFACE CHANGE of a different kind from the
vegetation artifact, and the vegetation correction is right to leave it: these patches are
the thing to FLAG for interpretation, not to correct away.

The deposition patch 1 (870 cells, +158 mm, E 579478 N 4885998) lies just NW of it along
the channel; a channel-migration pairing is possible and untested.

## K. The NE patch, resolved as far as the data allow: a closed grass canopy gen1 could not penetrate (2026-10-18)

Andy rejected both candidates in J (a sand sheet does not vanish; snow or a frozen
backwater is against spec) and asked for imagery and other ideas. Scripts
`analysis/floodplain_ne_patch.py` (geometry tests, appended to its stdout below) and
`analysis/floodplain_ne_imagery.py` (MnGeo WMS, figure `floodplain_ne_imagery.png`).

**Water is excluded by geometry.** A water surface is a horizontal plane. Inside the patch
gen1's absolute elevation is NOT constant: sd 309 mm, correlation with gen2's bed +0.83,
and corr(d, z2) = -0.02 -- the 25 cm offset is UNIFORM and follows the ground, which a pond
cannot do. The patch level (218.10 m) is 63 cm above the river's 2008 water surface from
gen1's class-9 returns in the channel (217.47 m), so it is not a connected backwater either.
The patch edge does not follow a gen2 contour (43% of nearby cells below the patch's 95th
percentile elevation lie inside it).

**One inference of mine is withdrawn.** The patch's returns are 87% within 2 deg of nadir
against 61% for the ring around it, which I read as specular drop-out. Scan angle along a
line is across-track position (corr(scan angle, easting) = +0.995 on line 138 here); the
patch simply lies under the nadir track. Geometry, not surface.

**Imagery, 2008-2023 (MnGeo; 2008 CIR, NAIP 2008/09/10, FSA 2013-2023).** The patch is an
open herbaceous meadow in every year: no water body, no bare sand, no crop rows, a farmstead
and road to its west, scattered shrubs increasing over time, tan (cut hay) in 2015, 2019
and 2021 and green in other years. gen2's May 2021 cloud there is clean open ground
(near-ground vegetation share 0.03, class-2 share 0.81).

**What it is, then (hypothesis, consistent with everything measured).** A dense, closed,
short grass canopy -- a floodplain hay meadow, reed canary grass being the likely species --
that in November 2008 stood about 25 cm tall with no gaps: gen1, at ~1 return/m^2 and 99%
single-return, read its top as the ground over the whole patch (uniform offset that follows
the terrain; dry grass is bright at 1064 nm, intensity 27 against 19-21). By May 2021 the
sward was short or flattened and gen2 at ~200 returns per cell read the ground. So it IS
vegetation -- my "not vegetation" in J was wrong in substance -- but it is a different
signature from the sedge: a closed canopy with ZERO penetration, where the sedge lets ~25%
of beams through. That is why the logistic correction cannot reach it: 70% of its cells
already carry the product's own penetration-failure flag (no gen1 return at the gen2
surface) against 9% elsewhere on the flat floodplain.

**Consequence.** The surviving patches are the UNOBSERVED-GROUND case the strategy note
anticipated. The right treatment is the flag the product already has, `penetration_failure`,
not a correction: there is nothing in gen1's cloud for any estimator to reach. The
deposition patch NW of it (870 cells, +158 mm) is in the forest and remains untested.

**Attribution (Andy, 2026-10-18, local knowledge):** "There is commonly reed canary grass.
These fields were since planted and used for hay. This might be our culprit." That matches
the measurements: a closed reed canary grass sward in November 2008 (uniform 25 cm
terrain-following offset, bright, zero penetration at gen1's density) and a planted,
hayed field by 2021 (cut-hay tan in the 2015/2019/2021 photos; clean ground in gen2's May
2021 cloud). The land-use change between the epochs is why the two surveys read different
surfaces here, and why no estimator on gen1's returns can recover the 2008 ground.

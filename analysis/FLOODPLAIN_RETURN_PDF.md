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

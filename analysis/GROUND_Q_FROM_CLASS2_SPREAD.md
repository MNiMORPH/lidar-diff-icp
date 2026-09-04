# The ground percentile from the class-2 spread — measured, and NOT adopted

Adopted 2026-09-03 for gen2 at Elba. **Withdrawn as the default 2026-09-04**, on measurement.
This document was rewritten that day: the earlier version reported the pooled result as a
success, and the code has overruled it. What follows leads with the finding that decides the
question.

## The finding

**On open ground the class-2 median is already right, and the spread-dependent curve makes it
worse.** Held out on 5 folds of 10 km spatially blocked marks, over the 227 NVA
(non-vegetated) control marks:

                          ground estimator  median err  |median|      RMS  p90|err|
               q = 0.50 (pipeline default)        -3.5       3.5     49.1      73.1
           q = 0.527 constant (calibrated)         0.1       0.1     48.7      73.9
    q = isotonic(log class-2 SD), held out        -5.8       5.8     52.5      76.5

`ground_q = 0.50` lands 3.5 mm from surveyed truth. The curve costs 3.4 mm of RMS to fix a
3.5 mm bias, and does not fix it.

So `difference_dem` takes `ground_q = 0.50` by default, and `ground_q="calibrated"` requires
`gen2_curve=` — a curve must be asked for by name, never arrived at by default.

## Why the earlier result looked good, and what was wrong with it

The version of this document written on 2026-09-03 reported a real gain:

                          ground estimator  median err  |median|      RMS  p90|err|
               q = 0.50 (pipeline default)         8.1       8.1    124.5     202.7
           q = 0.453 constant (calibrated)         0.0       0.0    117.7     189.2
    q = isotonic(log class-2 SD), held out        -2.4       2.4    104.6     171.0

Those numbers are correct and they are still reproducible. They were measured on **all 519
marks pooled**, and the control set is three different populations:

    point_type   n fitted   class-2 median minus surveyed truth
    NVA               227          -3.5 mm      non-vegetated, open ground
    VVA               162        +103.3 mm      sited UNDER VEGETATION, by design
    LCP               130         -23.1 mm      the acquisition's own calibration points

VVA checkpoints exist to measure how badly lidar performs under canopy. Pooling them into a
calibration of "where the ground sits in the return column" imports that canopy response into
every cell of every tile. The pooled curve's entire falling limb is the VVA population:

      SD    NVA only    all three pooled
      20 mm    0.494       0.506
      60 mm    0.476       0.445
     120 mm    0.476       0.395
     200 mm    0.468       0.237
     400 mm    0.392       0.101

On open ground the curve is nearly flat. The +8.1 mm bias that motivated the whole correction
is a pooling artifact: on NVA marks alone the bias is -3.5 mm.

## Why the curve is REGIONAL, and why that is a requirement

Not a compromise forced by thin data -- a condition for the measurement to mean anything.
The offset being measured is a ground-cover effect of a few tens of mm sitting on a per-mark
scatter of 49.1 mm (NVA, n=227). Resolving it needs many marks spanning the whole cover
range, and no single flight line or site supplies that:

    519 marks span 397 distinct flight lines over a 228 x 127 km survey
    marks per line: median 2, mean 2.1, min 1, max 8
    40% of lines carry exactly ONE mark; 89% carry three or fewer

A per-line fit would be two points. The REGIONAL COLLECTION IS REQUIRED for a robust
measurement of cover-based offsets (Andy, 2026-09-04), and the curve is therefore correctly
one curve per epoch, applied to every tile.

The cost is stated rather than hidden: a per-line level error cannot be separated from the
cover signal, because no line carries enough marks to estimate its own level. That is a
limit on attributing a residual to a line, not on the cover relation itself.

## The rule that already existed, and that this broke

`ground_control/run_bridge_gen2.py` makes `--point-types` a **required** argument —
"`--point-types NVA` is the consequence, not a preference". `README.md` lists *open ground
only* among the three rules governing the datum. `analysis/control_lowveg_offset.py` excludes
LCPs and says so at the top of the file. The original calibration here pooled all three types
and recorded the choice nowhere, because it did not know it was making one.

`analysis/calibrate_ground_q.py --point-types` is now required with no default, the types go
into the curve's filename and its provenance, and `groundq.load_curve` refuses a curve that
records none.

## How the defect surfaced

Not by inspection — by re-measuring the curve on the flight lines of a second site, which is
what should have happened before it was ever applied there. Whitewater's gen2 is five lines
(`point_source_id` 3029–3033). Nine control marks touch them, and their ranks are bimodal with
nothing in between:

      point_id     psid    rank   class-2 median - truth
    3149_2021_MN   3030  0.0006          +129.0 mm
    3089_2021_MN   3031  0.0235          +104.0
    2099_2021_MN   3032  0.0269           +43.7
    3085_2021_MN   3028  0.0335          +195.9
    1077_2021_MN   3030  0.7731          -133.1
    2080_2021_MN   3029  0.7813          -107.3
    1126_2021_MN   3032  0.8414           -68.1
    1019_2021_MN   3030  0.8537          -173.3
    2045_2021_MN   3032  0.8607           -40.2

    middle band pooled : 65.3% of 519 marks
    middle band on Whitewater's lines : 0 of 9
    binomial p = 7.260162231249149e-05

The split is *within* lines — psid 3030 carries 0.001, 0.773 and 0.854 — so it is not a
flight-line property. It is `point_type`, and it is systematic across all 519:

    point_type   rank<0.05   middle   rank>0.75
    LCP                  1       82          47
    NVA                 16      162          49
    VVA                 63       95           4

(The 0.05 / 0.75 band edges are chosen for this table and the binomial test only. They exclude
no marks; every mark appears in `figures/sites/whitewater_lines_rank_vs_spread.png`.)

## What applying it actually did

Measured on the cells both runs produced, so the comparison is like for like:

    ELBA        340,020 cells    gen2 at its median    NMAD 74.77 mm   median  +6.10 mm
                                 gen2 at calibrated q  NMAD 79.09      median  -3.33

    WHITEWATER  241,747 cells    without correction    NMAD 85.01 mm   median -11.75 mm
                                 with correction       NMAD 92.35      median -22.29

The correction moves the level and widens the scatter, at both sites. Elba had only ever been
checked on the median, which is why this went unnoticed for a day.

## The measurement that stands, restated honestly

The relation between the class-2 spread and where the ground sits **is real** — it is simply a
statement about vegetated ground, because that is where the marks that carry it were sited.
Pooled over the 519 marks:

    RANK OF SURVEYED GROUND WITHIN THE CLASS-2 RETURNS  (pooled -- see above)
      median 0.4532   mean 0.4354   p16 0.0550   p84 0.7798
      at exactly 0 (truth below every ground return): 16   at 1 (above all): 1

    SD bin mm     n  median rank      SE           95% CI
         0-30    58        0.344   0.137 [0.200,0.663]
        30-45    63        0.610   0.065 [0.507,0.728]
        45-60    86        0.587   0.042 [0.523,0.677]
        60-80    89        0.431   0.052 [0.299,0.547]
       80-110    94        0.382   0.041 [0.307,0.458]
      110-160    82        0.413   0.059 [0.255,0.496]
      160-250    44        0.211   0.111 [0.066,0.435]

    within <60 mm:   rho +0.006  p 0.935
    within >=60 mm:  rho -0.183  p 1.166e-03
    below 60: median rank 0.571 | 60 and above: 0.390 | Mann-Whitney p 3.461e-06

A curve fitted on VVA marks may well be the right tool for correcting ground *under canopy*,
which is a different and narrower claim than the one this document used to make. It has not
been tested that way.

## Reproduce all of it

    ./lidar-icp/bin/python analysis/calibrate_ground_q.py --set gen2_2021_control \
        --point-types NVA --diagnostics
    ./lidar-icp/bin/python analysis/calibrate_ground_q.py --set gen2_2021_control \
        --point-types NVA VVA LCP --diagnostics
    ./lidar-icp/bin/python analysis/marks_by_flight_line.py

`marks_by_flight_line.py` writes the per-mark table with each mark's flight line, which is what
makes the per-line check above a one-liner rather than a project.

## Code

    analysis/calibrate_ground_q.py --set <epoch> --point-types <TYPES...>
        -> data/derived/ground_q_vs_class2sd_<epoch>_<TYPES>.npz

    src/lidar_diff_icp/groundq.py
        calibrate : mark_statistics . spatial_folds . fit_curve . save_curve
        apply     : surface_from_grid . reference_surface . column_histogram .
                    spread_from_histogram . ground_at_q . ground_at_median . correct_gen2

    analysis/ridgelines/dod_cover_corrected.py --q-from-class2-spread <curve.npz>
    pipeline.difference_dem(ground_q="calibrated", gen2_curve=<curve.npz>)

Both callers use `groundq.correct_gen2`; there is no second implementation. The application is
intrinsically two-pass — the reference surface IS the q = 0.50 grid — and in the pipeline the
correction is applied as a post-registration delta so the gen1 tie cannot absorb it.

## Limits, all measured rather than assumed

- **Population, not cell.** Even pooled, held-out RMS is 104.6 mm and p90 error 171.0 mm.
- **16 marks of 519 have surveyed ground BELOW every class-2 return.** No percentile reaches
  the ground there.
- **Support mismatch.** Calibrated on 7.5 m discs, applied to 5 m cells.
- **Cells with < 20 class-2 returns are declined**, and a declined cell is dropped, never
  given the default. `ground_at_q` used to read a NaN percentile as bin 0 and return the FLOOR
  of the column: 1,131 Elba cells carried a median DoD of -1.091 m of fabricated erosion.
  Fixed 2026-09-04 with a regression test.
- **Per epoch, and now per point type.** Both are refusals in `load_curve`, not conventions.

## gen1: calibrated, and it does not work — leave it

    gen1, 269 marks (of 954; the rest lack 20 class-2 returns)
           q = 0.50 (pipeline default)        89.8      89.8    146.2     216.3
    q = isotonic(log class-2 SD), held out    69.0      69.0    138.1     197.0
    curve: q = 0.219 at EVERY spread, 20 mm to 400 mm

The isotonic fit is flat: in gen1 the class-2 spread carries no information about where the
ground sits. The +89.8 mm is NOT interpretable as vegetation — the 2008 control's own vertical
reference is uncertain, and the bundled CSV says so in its `verified` column. **This is a
statement about the 2008 control, not about the 2008 lidar.** Not pursued further
(Andy, 2026-09-03). Note this fit predates `--point-types` and pooled whatever the 2008 set
contains.

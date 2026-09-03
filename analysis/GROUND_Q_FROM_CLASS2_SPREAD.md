# The ground percentile from the class-2 spread — gen2, calibrated on surveyed control

Adopted 2026-09-03. Replaces the cover-relation route for gen2's ground estimate at Elba.

## What it replaces, and why

The pipeline takes the per-cell **median** of class-2 returns as the ground, `ground_q = 0.50`.
Every previous attempt to improve on that fitted a percentile against a COVER metric —
`q2 = 0.5 − 0.1922·cover`, then `q_ctrl = 0.4558 − 0.9264·lowveg` — and each ran into the same
three problems: the cover metric was built from chosen windows, the relation was fitted as a
straight line through something that is not straight, and applying it drove `q` outside [0, 1]
on 8–11% of tile cells, which then had to be clipped.

This route uses no cover metric at all. The cell's own class-2 returns supply both the
covariate and the estimate.

## The measurement

At each of 519 gen2 SE-Driftless control marks: take the class-2 returns within 7.5 m, express
them as slope-normal heights above the mark's order-2 surface, and ask **what rank the surveyed
elevation occupies among them**. That is the percentile the pipeline should be using at that
mark, measured rather than assumed.

    RANK OF SURVEYED GROUND WITHIN THE CLASS-2 RETURNS
      median 0.4532   mean 0.4354   p16 0.0550   p84 0.7798
      at exactly 0 (truth below every ground return): 16   at 1 (above all): 1
    class-2 MEDIAN minus truth: median +8.1 mm   mean +36.3   sd 119.2

## The covariate: the class-2 spread, which needs no windows

`lowveg` is a ratio of two chosen bands, `(0.15, 2.00]` over `(−1.00, +2.00]`, and its own
document records that the value moves ~50× as the lower edge is swept while the rank
correlation moves 0.07. The spread of the class-2 returns has no bands in it. Measured against
the same response, on the same marks, they are indistinguishable:

         covariate   windows?      rho           p
              sd_g       none   -0.225    2.11e-07
            nmad_g       none   -0.210    1.36e-06
            skew_g       none   +0.242    2.49e-08
            lowveg        TWO   -0.224    2.42e-07

Equal performance, one fewer arbitrary choice, and an axis in millimetres that means
something. `skew_g` is marginally the strongest single measure, and its sign says the ground
class widens **upward** — but it was not adopted, because the spread is the quantity with a
physical scale to compare against.

## The shape: flat, then falling (Andy's reading, tested)

    SD bin mm     n  median rank      SE           95% CI
         0-30    58        0.344   0.137 [0.200,0.663]
        30-45    63        0.610   0.065 [0.507,0.728]
        45-60    86        0.587   0.042 [0.523,0.677]
        60-80    89        0.431   0.052 [0.299,0.547]
       80-110    94        0.382   0.041 [0.307,0.458]
      110-160    82        0.413   0.059 [0.255,0.496]
      160-250    44        0.211   0.111 [0.066,0.435]

    within <60 mm:   rho +0.006  p 0.935        (no trend at all)
    within >=60 mm:  rho -0.183  p 1.166e-03
    below 60: median rank 0.571 | 60 and above: 0.390 | Mann-Whitney p 3.461e-06

The break sits near **60 mm**, which is close to the 59.3 mm bare-ground class-2 NMAD recorded
in `CONTROL_LOWVEG_OFFSET.md` as the surface's own noise. The rank holds while the ground class
is no wider than bare-ground noise, and falls once it is wider — i.e. once something that is
not ground is in the class.

A dip I initially read at 20 mm was noise: that bin's CI is [0.200, 0.663], wide enough to
contain everything around it.

## The calibration: isotonic, so no break is imposed

`analysis/calibrate_ground_q.py` fits a **monotone non-increasing** regression of rank on
log(class-2 SD). That reproduces flat-then-falling without a threshold to defend, without a
functional form, and without any cutoff. The only constraint is physical: more contamination
cannot mean a higher ground rank.

      SD   20 mm -> q = 0.506        SD  120 mm -> q = 0.395
      SD   40 mm -> q = 0.506        SD  200 mm -> q = 0.237
      SD   60 mm -> q = 0.445        SD  400 mm -> q = 0.101
      SD   80 mm -> q = 0.417

Note the plateau is at **0.506**, not the 0.571 quoted from the pooled sub-60 mm median, which
included the shoulder where the decline has begun. So `ground_q = 0.50` is RIGHT on clean
ground, and the whole gain comes from lowering `q` where the spread says the class is dirty.

## The validation: held out, spatially blocked

5-fold CV on 10 km blocks, so the curve never sees the marks it is scored on.

                      ground estimator  median err  |median|      RMS  p90|err|
           q = 0.50 (pipeline default)         8.1       8.1    124.5     202.7
       q = 0.453 constant (calibrated)         0.0       0.0    117.7     189.2
    q = isotonic(log class-2 SD), held out    -2.4       2.4    104.6     171.0

Against `q = 0.50`: the +8.1 mm median bias goes, RMS falls 16% and p90 error falls 16%. The
constant-`q` row is the control that matters — it gets about half the gain, so the spread
dependence is doing real work rather than dressing up a level shift.

## Applied to Elba

    class-2 spread: 353,296 cells with >=20 returns; median 45.0 mm  p90 141.5
    q from the curve: median 0.506  p10 0.329  min 0.101  (118,064 cells below 0.45)
    median DoD (mm): gen2 median +6.1  ->  cover-corrected -3.6
    existing dod.npy on the same cells: +7.6 mm  (n=340,961)

The median cell is left alone. A third of the tile — 118,064 cells — gets a lowered percentile,
the worst dropping to 0.101. Elba's median DoD moves +7.6 → −3.6 mm, a far gentler correction
than the cover route's +40.4 mm, and it acts only where the ground class is measurably dirty.

## Limits, all measured rather than assumed

- **Population, not cell.** Held-out RMS is 104.6 mm and p90 error 171.0 mm. This corrects a
  distribution; it does not fix an individual cell.
- **16 marks of 519 have surveyed ground BELOW every class-2 return.** No percentile reaches
  the ground there. Those cells are unrecoverable by any choice of `q`.
- **Support mismatch.** Calibrated on 7.5 m discs, applied to 5 m cells.
- **Cells with < 20 class-2 returns get NaN**, not the default — the curve says nothing about
  them, and 20 is the same minimum the calibration required of a mark.
- **gen2 only.** See below.

## gen1: calibrated, and it does not work — leave it

    gen1, 269 marks (of 954; the rest lack 20 class-2 returns)
           q = 0.50 (pipeline default)        89.8      89.8    146.2     216.3
    q = isotonic(log class-2 SD), held out    69.0      69.0    138.1     197.0
    curve: q = 0.219 at EVERY spread, 20 mm to 400 mm

The isotonic fit is flat: in gen1 the class-2 spread carries no information about where the
ground sits. And `q = 0.50` is biased by +89.8 mm, an order of magnitude more than gen2's
+8.1 mm — which is NOT interpretable as vegetation, because the 2008 control's own vertical
reference is uncertain. The bundled control CSV says so in its `verified` column: the datum is
a dataset-level assertion from `lidar_semn2008.html`, and "the validation reports themselves
state no datum and no geoid". The values are transcribed from tables labelled `Control Z` —
the control used in the adjustment, not independent checkpoints. A geoid-model ambiguity is
the right size to explain the offset, and independent gen1 datum work puts the surface at
+55.0 ± 16.6 mm and +51.3 mm by other routes.

**So this is a statement about the 2008 control, not about the 2008 lidar.** gen1's absolute
level should keep coming from the swath ties and the DoD closure. Not pursued further
(Andy, 2026-09-03).

## Code

    analysis/calibrate_ground_q.py --set gen2_2021_control
        -> data/derived/ground_q_vs_class2sd_gen2_2021_control.npz  (curve + its provenance)

    analysis/ridgelines/dod_cover_corrected.py --tile data/derived/elba \
        --gen2 data/after/3dep2021_fulldensity.laz \
        --q-from-class2-spread data/derived/ground_q_vs_class2sd_gen2_2021_control.npz
        -> dod_cover_q2.npy, gen2_q2_used.npy, class2_sd_mm.npy

The curve carries its own provenance — population, response, covariate, shape, CV scheme and
known limits — and the application echoes them on every run. `--q-from-class2-spread` refuses
to be combined with `--relation` or `--slope`, since all three set the same quantity.

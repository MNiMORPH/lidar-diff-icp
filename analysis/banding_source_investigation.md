# Source of the pervasive undulating banding in the DoD — and the fix

Investigation prompted by the ~40-m (red-to-blue), terrain-draped, feathered
banding visible across the DoD (`figures/final_dod.png`). Andy's framing: a
repeating, pervasive pattern must have a consistent source. It does — and it turns
out to be **largely correctable**. This note records what the banding is, the long
chain of things it is *not*, and the estimator change that removes most of it.
Numbers are per-tile on the Elba pilot (`4342-29-64`), 5 m grid, stable ground.

## Bottom line

The banding is a **processing artifact of horizontal low-percentile ground
estimation on sloped terrain.** Taking the lowest ~10% of points in a horizontal
x,y cell on a slope necessarily selects the **downhill** points, biasing the
"ground" downward by ~(offset·slope). The sparse 2008 cloud (19 pts per 5 m cell)
realizes that bias differently than the 12x-denser 2021 (232 pts), so the
mismatch shows up as coherent, slope-keyed hillslope bands. It is *not* real
topography (it does not correlate with actual landforms), and it is *not* a
navigation error.

**The fix (Andy's idea): take the low percentile normal to the local surface, for
both epochs against a common regional plane.** Detrending each point by the
regional slope before the low-pick removes the downhill bias at its root, and
doing it to *both* epochs against a *shared* plane means the tilt cancels and
neither side carries a residual bias to mismatch. Measured:

Validated in the **full pipeline** (registration applied), Elba tile:

| ground | stable σ | band-pass | convex false-dep (>LoD) |
|---|--:|--:|--:|
| horizontal low-10% (current) | 0.091 m | 0.035 m | 3.4% |
| **slope-normal (both epochs)** | **0.059 m (−35%)** | **0.028 m (−22%)** | **4.5%** |

The stable scatter drops a third and the bands a fifth, while the convex-hillslope
false-deposition metric that motivated low-10% (median gives 16-32%) is preserved
(3.4 -> 4.5%). A synthetic sloped surface with a known +0.80 m change recovers it
**identically** (0.726 m both methods) — so it is a targeted debias, not smoothing.
The one residue is a +14 mm shift on convex hillslopes (well below the LoD), most
likely slope-normal *removing* a density-driven false-erosion bias: the denser 2021
finds the true low point more reliably, so its ground sits lower and `z21 - z08`
skews erosional. A small real **roll ripple** (~5%) remains and is separately
correctable (`Δz ~ k(gps_time)·scan_angle`).

(An earlier scratchpad comparison reported a convex *regression* to 7.8% and a +8 cm
bias; that was an artifact of differencing without the pipeline's alignment/tie/drift
registration. Band-pass rstd is high-pass and immune to it; the convex metric is not,
so it must be measured inside the full pipeline, as above.)

## Mechanism (the confirmed cause)

On a planar slope, the true ground at cell center is the mean; the low percentile
picks points that are physically downhill and therefore lower. Raw magnitude of
this bias here: **median 0.33 m.** Because the 2008 low-pick is defined by only
~2 points per cell, *where* those points sit (hence the bias) varies cell-to-cell
with the scan pattern — deterministic given the points, coherent in space, and
different between epochs. It concentrates on hillslopes because that is where a
horizontal cell spans a large vertical range (slope × cell size). The band
*amplitude* tracks the 2008 within-cell elevation range (corr 0.25), not raw
density (0.03); the band *sign* is set by the sub-cell point arrangement, which is
why it is uncorrelated with any smooth terrain derivative.

## What the banding is NOT (each tested and refuted)

| hypothesis | test | verdict |
|---|---|---|
| real topographic variability | corr(DoD bands, real 2021 detail) | 0.01 — no (bands as strong on smooth-but-sloped ground) |
| DEM elevation quantization | DoD vs elevation-phase, all Δz | corr 0.026 — no |
| swath seams / overlap | 3 seams vs many bands; σ same 1 vs 2 swaths | no |
| along-track vertical drift f(gps_time) | overlap variance decomposition | only 0.02 of 0.096 m along-track — no |
| scan angle g(θ) | overlaps + per-swath vs 2021 | not identifiable / flat — no |
| point density (linear covariate) | corr(DoD, count) | 0.06 — no |
| intensity | corr(DoD, intensity) | −0.03 — no |
| horizontal misregistration | DoD ~ a·∂z/∂x + b·∂z/∂y, global / local-CV / all lags | R²≈0, CV −0.04 — no |
| single periodic carrier | 2-D autocorrelation | monotonic decay, no peak — no |
| navigation roll ripple | time-resolved cross-track tilt k(gps_time) | real but only ~5% of the bands |

## The removal-attempt arc — why everything failed until slope-normal

The decisive clue was that **one-sided** fixes made it *worse*: the horizontal
low-pick already partially cancels in `z21 − z08` (both biased downhill in the
same direction), so removing the bias from only one epoch breaks the cancellation.

| method | Δ bands | why |
|---|--:|---|
| directional band-pass filtering | −4% | no consistent azimuth (bands are reticulate, not striped) |
| per-cell plane detrend | +214% | overfits ~2 low points on sparse data |
| wider-support low-percentile | +445…832% | pooling grows the downhill bias with radius |
| covariate / density-bias regression | +6% | removes only the mean, not the per-cell realization |
| matched density (decimate 2021) | +11% | independent realizations don't cancel |
| centroid de-bias (1-sided) | +196% | breaks the shared-bias cancellation |
| regional detrend (1-sided) | +166% | same |
| **slope-normal low-10%, BOTH epochs, common plane** | **−22%** | removes the bias from both; nothing to mismatch |

## The fix, in detail

For each point, subtract the regional plane (from a lightly smoothed 2021 surface,
so the slope is robust and not a per-cell overfit): `resid = (z − [Z_reg(cell) +
Δx·∂z/∂east + Δy·∂z/∂north]) · cos(slope)`. Grid the 10th percentile of `resid`
per cell, per epoch. The shared plane cancels in the difference, so
`DoD = low10(resid_2021) − low10(resid_2008)` preserves real change while dropping
the slope-driven bias. Reproducible driver: `analysis/slope_normal_ground.py`.
Figure: `figures/slope_normal_dod.png`.

## Status and remaining work

- **Integrated ✓** behind `pipeline.difference_dem(..., ground="slope_normal")`
  (default `"low_q"` is byte-identical; the 4 pipeline tests pass). Validated in the
  full pipeline (table above).
- **Slope scale ✓** swept: tighter is slightly better (3 m −24%, 6 m −22%, 15 m
  −19%); a smoothed *varying* slope (~6 m) is used — not a constant, and not a
  per-cell fit (which is too textured: only −3%).
- **Curvature ✓** ruled out as the lever: a curvature-aware (bilinear smoothed-
  surface) reference gives identical band rstd and the same convex metric as the
  linear plane, so the simple linear plane suffices.
- **Memory ✓** `stream=True` grids the after cloud in chunks via a per-cell
  histogram (O(cells) RAM, not O(points)); validated to reproduce the exact product
  to <4 mm on stable ground with identical stable σ. Sparse cells need a min-count
  mask. In-memory default (`stream=False`) unchanged.
- **Roll ripple** (~5%, `k(gps_time)·scan_angle`) is orthogonal and separately worth
  adding to the statewide stack.
- **Default:** `low_q` remains the default pending broader validation; `slope_normal`
  is opt-in.

## One-line answer

Horizontal low-percentile ground estimation is downhill-biased on slopes; the
sparse-2008 / dense-2021 mismatch turns that into coherent hillslope bands.
Measuring the low percentile **normal to the surface, for both epochs against a
common plane**, removes ~a third of it while preserving real change — a
deterministic, generalizable correction for the statewide workflow.

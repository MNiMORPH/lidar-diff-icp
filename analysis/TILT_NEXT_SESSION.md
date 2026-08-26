# The tilt — next session's job

Measured 2026-08-26, in-session, not yet written by a committed script. **Re-run and confirm
before building on it.**

## The measurement

`dod_cover_q2.npy` (gen2 − gen1, mm, + = elevation rose) regressed on easting and northing
over stable reference cells: `reference_cells()` defaults (ridge, |curv| <= 0.015, slope < 12
deg, no buildings, not clear-cut, not blufftop margin, |DoD| <= 500 mm), further restricted to
`canopy_cover_pfs <= 0.02`, `floodplain_mask` removed. n = 24,287 cells. Errors are block
bootstrap over 50 m blocks.

| term | estimate | block-boot SE | t |
|---|---|---|---|
| intercept | −10.30 mm | 3.09 | −3.33 |
| dE | **−14.19 mm/km** | 5.15 | −2.75 |
| dN | **−16.70 mm/km** | 3.65 | −4.57 |

Across the tile: **~58 mm of ramp in northing** (3.50 km) and ~36 mm in easting (2.54 km).
That is larger than the canopy term (5–60 mm), larger than the epoch-difference argument
(5.5 mm), and it means the −4.62 mm tile-median offset on divides is the average of a plane,
not a constant.

## Independent, weaker corroboration

gen2's own six surveyed marks, fitted as a plane (weighted, 4 NVA and all 6):

| fit | dE mm/km | dN mm/km | resid RMS | dof |
|---|---|---|---|---|
| NVA only (4) | +3.8 ± 2.4 (t 1.55) | **−5.3 ± 1.8 (t −2.86)** | 24.6 mm | 1 |
| all six | +1.4 ± 2.7 (t 0.51) | **−6.4 ± 2.1 (t −2.99)** | 36.7 mm | 3 |

Same sign in northing, about 2.6× shallower, on 4–6 points fitting 3 parameters. Suggestive,
not established. East–west not significant in the marks.

For reference, the tilt term our registration DOES apply (`cross_epoch_datum`) is
+0.78 and −0.57 mm/km — twenty times smaller than what is left over.

## The leading hypothesis, and why it is testable

gen1's flight lines run NORTH–SOUTH, spaced ~960 m. The per-swath constants at Elba are
0, −23.9, −32.5, −43.7 mm running west to east across ~2 km — **about −22 mm/km in EASTING**,
close to the measured dE. And the along-track drift, being along-track, runs NORTH–SOUTH,
the same axis as dN. So both components of the "tilt" may be residuals of the two registration
terms already in the model, projected onto the map axes.

**The separating test:** fit the plane WITHIN a single swath, where the across-swath constant
is fixed by construction. If dE collapses within a swath and dN does not, the two mechanisms
are cleanly distinguished. Agent `a8831934b5349b50e` was running this at session end —
check `analysis/STABLE_POINT_TILT_AUDIT.md`.

## What must be checked before believing the tilt is real

1. **Where the stable cells are.** They are ridge-top, low-curvature, low-slope, treeless —
   in a dissected landscape a spatially non-uniform, largely agricultural-upland population.
   A plane fitted to a clustered or elevation-biased subset is not a plane over the tile.
2. Whether the tilt is stable across halves/quadrants, or carried by one region.
3. Whether it is really a dependence on ELEVATION or distance from the valley.
4. The within-swath test above.

## Consequence if it survives

A two-parameter correction we currently do not make, worth ~58 mm across the site. It would
also mean much of what we have been calling a diffuse "residual field" (block medians −107 to
+50 mm over 500 m blocks, real between-block sd 31.7 mm after removing 9.6 mm sampling noise)
is one plane rather than noise.

## Reproduce

The regression above was run inline and is NOT yet a committed script. Re-derive it, as a
script using `trust/provenance.py`, before using the numbers.

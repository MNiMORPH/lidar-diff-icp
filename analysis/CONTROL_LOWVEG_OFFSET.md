# gen2's ground offset scales with near-ground vegetation, measured on surveyed control

**Reproduce:** `./lidar-icp/bin/python analysis/control_lowveg_offset.py --sweep`
(local data only; no network). Acquisition: `analysis/cover_at_control_marks.py`.

## The design

Two quantities at the same 389 points, independent by construction:

| | source | measures |
|---|---|---|
| **density** | the cloud's own shape, self-referenced | how much low vegetation is present |
| **offset** | surveyed elevation vs the vendor's delivered surface | how wrong the lidar ground is there |

The density metric is the fraction of returns in (0.15, 2.0] m above an order-2 least-squares
surface fitted to class-2 returns within 7.5 m. Because that surface comes from the box's own
returns, the histogram is **invariant to any vertical shift of the cloud** — the metric
structurally cannot carry offset information. The offset is USGS's published
`surveyed_Z − delivered_LAZ_Z`, from an estimator we had no hand in. Single epoch: **erosion
cannot enter.**

LCPs are excluded — the 143 LiDAR Control Points calibrated the acquisition, NVA/VVA were held
out. The script asserts no LCP carries a residual.

## Result

```
          bin    n   median     mean      SE
  0.00-0.06  228     -8.1    -13.0     2.9
  0.06-0.12   61    -26.2    -35.8    11.2
  0.12-0.18   40    -21.6    -39.0    10.8
  0.18-0.24   20    -84.5    -78.6    17.1
  0.24-0.30   19    -99.3   -118.6    28.0
  0.30-0.36   10    -95.4   -110.0    23.9
  0.36-0.42    6    -60.0    -61.9    19.5
  0.42-0.48    4   -132.1   -131.3    16.9
  0.48-0.54    1   -234.8   -234.8     nan

  DESIGN-weighted  1/SE^2 : intercept    -5.7 +/- 4.3   slope   -260.7 +/- 37.1 mm per unit
  ABUNDANCE-weighted by n : intercept    -4.7 +/- 6.4   slope   -305.1 +/- 49.0 mm per unit
  per-mark slope:  naive -293.1 +/- 31.8 (p 2.2e-18)
                   block bootstrap, 175 blocks of 10 km: -293.0 +/- 37.8  (SE x1.19)
```

**~ −29 mm per 0.1 of low-vegetation fraction. Intercept −5.7 ± 4.3 mm: on bare ground the
delivered surface is unbiased against survey.** The bias is vegetation, essentially all of it.

## Checks

**Not a block or phenology artifact.** All five EPT blocks give a negative slope, four
significant, despite spanning 2021-04-16 to 2022-06-05:
`_1_ −468.0 ± 107.6 · _2_ −242.8 ± 50.1 · _3_ −202.2 ± 53.6 · _4_ −144.3 ± 76.8 · _5_ −188.3 ± 51.3`.

**Weighting does not decide it.** Design vs abundance differ by 44 mm per unit, inside
overlapping errors.

**Spatial correlation costs little.** SE inflation 1.19x, against the 1.38–1.42x that
flight-line clustering produced for the datum work. Control marks are dispersed.

## Two qualifications that must travel with the coefficient

**1. The metric is ORDINAL.** Moving the band's lower edge over ±0.10 m changes the metric's
value by a factor of ~50 while `rho` stays in −0.30 … −0.37:

```
 lower edge  median lowveg   rho vs offset          p
       0.05 m         0.3056          -0.346   2.39e-12
       0.15 m         0.0405          -0.366   9.69e-14
       0.25 m         0.0059          -0.302   1.15e-09
```

The ordering of marks is stable; the scale is not. **A slope in mm-per-unit is meaningless
without the band definition attached.**

**2. Canopy cover is the WRONG covariate here, and this is why.** `cover_r7.5` at
`min_height = 2.0` gives `rho −0.053, p 2.94e-01` on the same marks. 131 of 162 surveyor-classed
VVA (vegetated) marks show canopy cover below 0.05, while carrying 5x the near-ground returns
(median 0.1117 vs NVA's 0.0218). The vegetation is **below 2 m** — and lowering `min_height`
does not fix it: `mh=0.15` and `mh=0.5` are identical row-for-row, because `voxel_z = 1.0`
puts everything from bare ground to 1 m in one voxel. The voxel height, not the threshold, is
the binding constraint.

## Sign convention — read this before any number above

`offset = surveyed_Z − delivered_lidar_Z`, mm.

* **positive** = surveyed ground is ABOVE the lidar surface → the lidar reads **low**
* **negative** = surveyed ground is BELOW the lidar surface → the lidar reads **high**

A **negative correlation with `lowveg`** therefore says: more low vegetation → true ground
sits further beneath the lidar surface → **the lidar's ground floats up**. Vegetation returns
enter the ground class and drag the fitted surface above the soil. About −8 mm on bare ground
(unbiased) to −235 mm at the densest mark. For the DoD this is the sign that makes vegetated
ground look like it gained elevation.

## Scatter is a symptom, not a second driver

`--scatter`. Class-2 IQR runs 60 mm on bare ground to 630 mm in the densest bin, and tracks
vegetation almost perfectly (Spearman +0.83 to +0.89). Against the offset it does no better
than `lowveg` itself (−0.36 to −0.40 vs −0.366), and the partials settle it:

```
  class-2 NMAD | controlling for lowveg : r +0.035  p 4.85e-01
  lowveg       | controlling for NMAD   : r -0.253  p 4.10e-07
```

Spread carries NOTHING about the offset once vegetation is accounted for. It is the tempting
covariate — computable anywhere, no classification decisions — and it is the worse one.
(IQR values are quantised at the 20 mm bin width; treat the ladder as ordinal.)

## Slope: removed on one side, negligible on the other, and not the cause

`--slope-check`. The scatter is slope-normal by construction — an order-2 surface removes
slope AND curvature as a trend, and the `|n|` division converts vertical to perpendicular.
The OFFSET is a vertical difference and is not converted, but median slope at these marks is
2.59°, so the conversion is 0.10% typical / 6.7% worst and `rho` is `-0.366` either way.

Slope does covary with vegetation (`+0.394`) — steeper ground is less farmed — so it is a
real candidate confound. It survives: `lowveg vs offset controlling for slope, r -0.349,
p 1.31e-12`.

## It is the lowest layer, not the canopy

`--strata`. Extending the upper edge from 0.5 m to 45 m makes the relation WEAKER, and the
0.25–0.5 m band alone is the strongest of the sweep:

```
  0.25-  0.5 m       0.0016          -0.316    1.88e-10
  0.25-  2.0 m       0.0050          -0.298    2.04e-09
  0.25- 45.0 m       0.0087          -0.265    1.07e-07

   0.25-  2.0 m       0.0050          -0.298    2.04e-09
   2.00-  5.0 m       0.0000          -0.109    3.08e-02
   5.00- 10.0 m       0.0000          -0.076    1.34e-01
  10.00- 45.0 m       0.0000          -0.069    1.75e-01
```

Ankle-to-knee height does the work: stubble, grass, low brush. **Limit:** the tall strata have
a median fraction of exactly 0.0000 — surveyors do not place marks under closed canopy — so
this cannot test tall canopy. It shows only that the near-ground layer alone reproduces the
full signal.

## Open

- **Within-block slopes span −144 to −468, a factor of three.** Phenology is the obvious
  suspect and is now testable: every mark carries its own `gps_utc_min`.
- **The marks under-reach the tile.** At Elba, 5.4% of reference cells exceed any mark's
  vegetation density, and those cells carry 17–24% of the total correction depending on
  whether the relation is held flat or extrapolated beyond the fitted range. State the range
  of validity; do not apply blind to the dense tail.
- Use the **published** residual, never ours: our least-squares surface diverges from the
  vendor's by ~62 mm in vegetation and ~2 mm in the open — i.e. as a function of the very
  covariate being fitted.

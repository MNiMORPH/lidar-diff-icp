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

# GEN2 incidence-angle test: does the beam-to-surface incidence control apparent floor elevation?

**Question.** For gen1 (2008 MN lidar) the beam-to-local-surface **incidence angle** controls the
apparent forest-floor elevation: near-slope-perpendicular beams (low incidence) read the ground
DEEPER, oblique-to-surface beams (high incidence) read it HIGHER (~+3.5 mm/deg; ~35 mm gap across
16-24 deg slopes). **Does gen2 (2021 USGS 3DEP, leaf-on, ~24x denser) show the same dependence?**

**Answer: NO.** Gen2 shows **no meaningful incidence dependence** of the forest-floor elevation.
At the controlled fixed-slope band (core forest, 16-24 deg slope) gen2 gives **+0.05 mm/deg**
(perpendicular-to-oblique gap **3.6 mm**, non-monotonic), against gen1's **+2.97 mm/deg**
(**+90 mm gap**, monotonic) computed the identical way. Relative to the ~20 mm signal budget, the
gen2 incidence effect is **below threshold** (~3-4 mm, comparable to noise), whereas the gen1 effect
(~90 mm at fixed slope; ~35 mm quoted over the slope range) is a **first-order** contaminant. The
2008 incidence artifact is a gen1-specific problem and does **not** reappear in the modern survey.

## Data sources (gen2), and a data-provenance caveat that shaped the method

- **Per-return angle/time/flight-line fields survive ONLY in the raw EPT tiles**
  `data/after/_ept_tiles/*.laz` (3DEP EPT pull, **EPSG:3857** web-mercator, LAS point format 1,
  6855 tiles; 4184 overlap the study grid). The merged products `data/after/3dep2021_*.laz` have
  **`scan_angle_rank` and `gps_time` zeroed and `point_source_id` collapsed** by the entwine/PDAL
  merge that produced them (verified: every chunk reads exactly 0.0 for both fields). They are
  **unusable** for incidence reconstruction; this test therefore streams the EPT tiles directly and
  reprojects **3857 -> 26915 (UTM 15N)** onto the grid. (pyproj is broken in the venv; this must be
  run in the conda `lidar-icp` env with `PROJ_DATA` set — see the script header.)
- **`scan_angle_rank`**: point format 1, 1 deg/unit, signed; observed range **-20 .. 20 deg** on the
  in-grid ground returns (wider/coarser-quantized than gen1's 0.006-deg `scan_angle`; range is
  sensible). **Class 2 = ground** (n = **51,419,059** in-grid returns); **class 7 = noise, excluded**
  (20,282 seen).
- **Reference bare earth** = `data/derived/elba_fulldensity/z_after.npy` (this IS gen2's own bare
  earth, so `d_mm` centers near 0 by construction; the test is the **relative** dependence of `d_mm`
  on incidence, ideally at fixed slope). Slope from `np.gradient(z_after)`.
- **Leaf state**: gen2 is 2021-05-01 green-up / **leaf-on** (per project memory: NDVI ~0.49);
  gen1 is Nov-2008 dormant. So this test asks the incidence question in gen2's leaf-on regime.

## Method (reconstructed the SAME validated way as gen1)

Per flight line (`point_source_id`): heading H from fitting (x,y) vs `gps_time`; cross-track unit
`c=(-sinH,cosH)`; side of +scan_angle from `sign(corr(cross-track pos, scan_angle))`. Beam
horizontal unit (ground->sensor) `= -sign(scan_angle)*sgn*c`; beam `b = sin|θ|·ĥ + cos|θ|·ẑ`.
Surface normal `n=(-gx,-gy,1)/|.|`; `incidence = arccos(b·n)`. `d_mm` = slope-normal distance of
each ground return to the `z_after` plane, ×1000. (`analysis/ridgelines/incidence_angle.py`,
`gen1_save_angles_slope.py`.) The in-grid returns belong to flight lines 3039-3043;
**corr(cross-track, scan_angle) = -1.00** for the four large lines (-0.83 for the smallest),
confirming a clean side-sign reconstruction.

### Flat-ground validation (incidence must reduce to |scan angle|)

On flat OPEN farmland (slope < 2 deg), median reconstructed incidence tracks |scan angle|:

| \|scan\| band (deg) | median incidence (deg) | expected | n |
|---|---|---|---|
| 0-2   | 1.2  | ~1  | 554,588 |
| 4-6   | 4.5  | ~5  | 663,403 |
| 8-10  | 8.5  | ~9  | 543,224 |
| 12-16 | 12.9 | ~14 | 260,911 |
| 16-24 | 18.2 | ~20 | 215,840 |

(The small low bias at high bands is expected: on a near-flat cell incidence ≈ |scan|, and the
median |scan| within a band sits below the band's upper-center; the 1-deg quantization of
`scan_angle_rank` also broadens the mapping.) On steep forest (slope > 20 deg) median |scan| = 11.0
but median incidence = 28.2 deg — incidence correctly departs from scan angle on slopes. **Validation
passes.**

## Results — median gen2 `d_mm` per incidence band

**Stratification 1 — CORE farmland + CORE forest combined** (`core_open` | `core_forest`):

| inc band (deg) | median d (mm) | n |
|---|---|---|
| 0-5   | -0.6 | 2,819,959 |
| 5-10  | -0.6 | 3,446,503 |
| 10-15 | -0.2 | 1,753,344 |
| 15-25 | -0.2 | 1,493,641 |
| 25-40 | -2.8 |   845,329 |

**gen2 combined: -0.07 mm/deg** (flat; dominated by the abundant flat farmland).

**Stratification 2 — CORE forest alone** (`core_forest`):

| inc band (deg) | median d (mm) | n |
|---|---|---|
| 0-5   | -5.0 |  77,649 |
| 5-10  | -2.6 | 228,911 |
| 10-15 | -1.0 | 237,640 |
| 15-25 | -0.9 | 679,573 |
| 25-40 | -2.8 | 841,174 |

**gen2 core forest: +0.06 mm/deg** (all-slopes; band medians span only ~4 mm). Adding the flat
farmland (Stratification 1) does not change the conclusion — both are ~0 mm/deg — but it does confirm
the farmland is flat and near-zero, so the forest-only number is not being propped up by open ground.

### Controlled test — CORE forest at FIXED slope 16-24 deg (removes incidence-slope collinearity)

Incidence and slope are collinear (oblique incidence concentrates on steep cells), so the honest test
holds slope fixed. This is the gen1 headline band.

| inc band (deg) | gen2 median d (mm), n | gen1 median d (mm), n |
|---|---|---|
| 0-5   | **-8.8**  (13,044)  | **-132.6** (284)    |
| 5-10  | -3.1  (83,986)      | -110.6 (3,728)      |
| 10-15 | +0.3  (101,943)     | -112.5 (9,154)      |
| 15-25 | -3.1  (130,841)     | -79.5  (51,535)     |
| 25-40 | -5.1  (248,054)     | -42.2  (19,387)     |
| **fit** | **+0.05 mm/deg** | **+2.97 mm/deg** |
| **perp→oblique gap** | **+3.6 mm** (non-monotonic) | **+90 mm** (monotonic) |

Gen1's +2.97 mm/deg fixed-slope figure reproduces the previously reported ~+3.5 mm/deg. Gen2 is flat
and non-monotonic at ±5 mm — i.e. **null**.

## Interpretation, in mm relative to the ~20 mm signal budget

- **gen1**: the incidence artifact is **~90 mm** at fixed 16-24 deg slope (~35 mm over the slope
  range as previously quoted) — several times the ~20 mm signal budget, a **first-order** contaminant
  that must be corrected.
- **gen2**: **~3-4 mm** span and no consistent sign (+0.05 mm/deg fixed-slope; +0.06 mm/deg
  all-slope forest; -0.07 mm/deg combined) — **within the noise / below the ~20 mm budget**. Gen2 does
  **not** exhibit the perpendicular-deeper / oblique-higher behavior.

**Why plausibly different.** The gen1 effect was tied to 2008 pulse/waveform behavior (broad returns
penetrating deeper near-perpendicular). Gen2 is a modern, ~24x-denser 3DEP survey with a narrower
effective footprint and dense multi-return ground; the leaf-on canopy adds a separate, incidence-
**independent** near-zero-to-slightly-negative floor offset rather than the strong incidence ramp.
Whatever residual structure gen2's forest floor carries (median -1.6 mm, MAD 62 mm over 2.26 M core-
forest returns), **it is not organized by beam-to-surface incidence.** The 2008 incidence correction
is a gen1-only fix.

## Reproduce

```
PROJ_DATA=/home/awickert/anaconda3/envs/lidar-icp/share/proj \
  /home/awickert/anaconda3/envs/lidar-icp/bin/python analysis/ridgelines/gen2_incidence_test.py
```

Outputs: `data/derived/elba_fulldensity/gen2_csf_angles.npz` (incidence, scan_angle, slope, d_mm,
cell, point_source_id, stratum, core_forest, core_open) and
`figures/refdatum/gen2_ground_pdf_vs_incidence.png`.

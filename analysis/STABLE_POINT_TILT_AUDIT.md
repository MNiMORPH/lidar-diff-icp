# Where do the "stable" points come from, and is the DoD tilt over them real?

**Date:** 2026-08-26
**Script:** `analysis/stable_point_tilt_audit.py`
**Run:** `env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python analysis/stable_point_tilt_audit.py`
**Inputs:** `dod_cover_q2.npy`, `z_after.npy`, `canopy_cover_pfs.npy`, `floodplain_mask.npy`,
`ag_region.npy`, `corrections.json` and `beam_offset_table.parquet` for `elba_fulldensity`.
Every input, parameter, mask and column is declared through `trust/provenance.py`; the full
banner prints with the run. `pipeline.py` and `coreg.py` are imported, not modified. No
minimum count is imposed anywhere; the three parameters the assistant chose unasked
(`block_m` ladder, `n_boot`, `n_elev_bins`) are printed as such in the banner and none of
them removes an observation.

## The claim under audit

A plane fitted to `dod_cover_q2.npy` (gen2 − gen1, mm, + = elevation rose) over the stable
open reference cells returns

    intercept  -10.30 mm      dE  -14.19 mm/km      dN  -16.70 mm/km

**Reproduced exactly**: 24,287 cells, `intercept -10.29884, dE -14.18607, dN -16.70099`.
The population is `refcells.reference_cells()` at its repo defaults (divide cells,
|curv| ≤ 0.015, slope < 12°, no building returns either epoch, not clear-cut, not in the
blufftop margin, |DoD| ≤ 500 mm), then `canopy_cover_pfs ≤ 0.02` and `floodplain_mask`
removed. Cluster-robust standard errors at 50 m blocks are ±1.77 / ±3.24 / ±2.10 and a
500-replicate block bootstrap over the same blocks gives ±1.81 / ±3.37 / ±2.21 – the two
agree, and both are *smaller* than the ±3.09 / ±5.15 / ±3.65 originally quoted, which sit
between the 50 m and 100 m rungs of the ladder in §2.

## Answer in one paragraph

The easting term is not a property of the tile. It is carried entirely by the 19% of the
population that sits on the valley floor rather than the upland: on the upland limb
`dE = +3.5 ± 6.6 mm/km` (nothing), on the valley limb `dE = -84.8 ± 10.8`. It also
disappears inside the two flight lines that have any easting leverage (`-4.0 ± 4.4` and
`+0.7 ± 7.2`), and it is not significant at any block size ≥ 250 m even before those
controls. The northing term survives every control – elevation, distance from the valley,
land use, L1 in place of L2, checkerboard halves, block size out to 1 km, and refitting
inside each of the four gen1 flight lines – at about `-16 to -18 mm/km`. It is **not**
residual along-track drift: consecutive gen1 lines fly opposite directions, so the sign of
a drift residual in northing must alternate line to line, and it does not (all four lines
give a negative `dN`, while their slopes in `gps_time` alternate `+1.09, -1.11, +0.96,
-2.21 mm/s`). What the northing term actually is, however, is a straight line drawn through
a **patchy field, not a ramp**: 500 m northing bands on the upland run `+12, +34, +6, -30,
-7, (+1), -22 mm`, and the residual is still correlated at 400 m and anticorrelated beyond
1.6 km. A linear tilt is a poor description of it.

---

## 1. Where the cells are

24,287 cells = **6.83%** of the 355,600-cell tile (700 × 508 at 5 m). They occupy 122 of the
154 250 m blocks (79%), but very unevenly: of the 2,500 cells a full block could hold, the
occupied blocks carry a median of **111**, quartiles 34–321, maximum 821. Fourteen of the 122
occupied blocks hold fewer than ten cells, and 44 hold fewer than fifty.

    cells per 250 m block (north at top, west at left)
       233    1  262  308  290   37    .    .  144  330    .
       425    2    1   84  233   80    2   46  258  455    6
         .    .   61  101  254   41  219   34    6   38    2
         .    .    .    .    .   89    .   93    .   32   11
         .  203    .  261   40    .    4  821  103   90    .
       355  511    .  380  403  381  462  312  109  127    3
       278  228    .  686  473  647  557   41   24    .    .
       325   16  485  661   82   37   35    .    .    3   20
        34  150  592   31    .    .    .    2    .  127    2
       422  656  416   58  148  190  293   45    .   71    3
       543  494  695  176    2  406  535   24  109   87    .
       430  272   10   55    .   32   25  217  114    .    .
       337   78   87   44   31  151  514  156   43   34   31
       292  401  682  115   20  219  181   10    .   26   28

The cells span the tile's full easting (577,495–580,030 against bounds 577,493–580,033) and
northing (4,882,740–4,886,235 against 4,882,738–4,886,238), so the plane is not being
extrapolated. **Elevation is where the population is unrepresentative:**

| population | n | p0 | p5 | p25 | p50 | p75 | p95 | p100 |
|---|---|---|---|---|---|---|---|---|
| stable cells under audit | 24,287 | 216.5 | 219.5 | **332.0** | 338.1 | 340.8 | 344.3 | 352.2 |
| reference cells (all cover) | 57,056 | 216.4 | 218.6 | 224.1 | 300.4 | 337.3 | 342.7 | 352.2 |
| whole tile | 354,717 | 214.1 | 219.8 | 242.6 | **293.8** | 331.3 | 341.4 | 353.4 |

Three quarters of the audited cells sit above 332 m while the tile's own median is 293.8 m.
The distribution is **bimodal** – modes at 220–225 m and 335–340 m with the antimode at
270 m – so it is two populations: an **upland limb of 19,683 cells (81%)** and a **valley
limb of 4,604 cells (19%)** on the terraces below, with almost nothing on the hillslopes
between (they are steeper than the 12° cut, or forested). And in this dissected valley the
uplands lie west: **corr(z, easting) = −0.623**, against corr(z, northing) = −0.084. **Any
easting term fitted on this population is confounded with elevation by construction.**

## 2. The plane, and what its uncertainty depends on

### 2a. The standard error depends strongly on the block size

| block_m | n | blocks | dE | dE_se | t_dE | dN | dN_se | t_dN |
|---|---|---|---|---|---|---|---|---|
| 50 | 24,287 | 1,152 | −14.19 | 3.24 | −4.38 | −16.70 | 2.10 | −7.94 |
| 100 | 24,287 | 438 | −14.19 | 4.56 | −3.11 | −16.70 | 2.97 | −5.61 |
| 250 | 24,287 | 122 | −14.19 | 8.39 | **−1.69** | −16.70 | 5.06 | −3.30 |
| 500 | 24,287 | 40 | −14.19 | 10.34 | **−1.37** | −16.70 | 6.08 | −2.75 |
| 1000 | 24,287 | 12 | −14.19 | 11.33 | **−1.25** | −16.70 | 6.50 | −2.57 |

r² of the plane is 0.0880; mean DoD −10.30 mm, median −4.69, NMAD 51.5, sd 69.7.

### 2b. The correlation length says which rung to read

Semivariogram of the plane residual (400,000 random cell pairs, residual sd 66.6 mm):

| lag (m) | 0–25 | 25–50 | 50–75 | 100–150 | 200–300 | 300–400 | 600–800 | 800–1200 | 1600–2400 | 2400–3600 |
|---|---|---|---|---|---|---|---|---|---|---|
| ρ | +0.64 | +0.61 | +0.52 | +0.37 | +0.25 | **+0.26** | +0.22 | +0.11 | −0.07 | −0.61 |

The field is still correlated at +0.25 past 400 m and only crosses zero near 1.6 km, going
strongly negative beyond – the signature of a large-scale undulation, which is exactly what
a linear term picks up. **A 50 m block therefore treats correlated cells as independent, and
the 50 m row above is the optimistic end.** Read the 250–500 m rows: `dE` is then not
distinguishable from zero (t = −1.7, −1.4) and `dN` is 2.8–3.3σ.

### 2c. L2 against L1, and one vote per block

|  | intercept | dE | dN |
|---|---|---|---|
| OLS, per cell | −10.30 | **−14.19** | −16.70 |
| LAD, per cell | −4.04 | **−1.38** | −15.38 |
| OLS on the 122 250 m block medians, one vote per block | −6.88 ± 5.25 | −18.37 ± 6.31 | −11.77 ± 4.74 |

The repo's standing preference is L1 (right-skewed residuals give L2 a ~10 mm median bias,
`FRAME_2026-08-26.md`). Under L1 **the easting term collapses from −14.19 to −1.38** while
the northing term is unchanged (−15.38). The mean-minus-median gap (−10.30 against −4.69)
is a left-heavy tail: p1 = −269 mm against p99 = +131 mm.

## 3. Is the tilt a property of the tile? No – the coefficients are not stable across it

| subset | n | blocks | mean_mm | dE | dE_se | dN | dN_se |
|---|---|---|---|---|---|---|---|
| ALL | 24,287 | 1,152 | −10.30 | −14.19 | 3.24 | −16.70 | 2.10 |
| west half | 13,405 | 575 | −4.37 | −11.39 | 5.68 | −12.48 | 1.95 |
| east half | 10,882 | 593 | −17.60 | **−55.64** | 8.78 | −20.87 | 3.52 |
| south half | 11,666 | 593 | +7.43 | −1.02 | 3.50 | −9.60 | 5.50 |
| north half | 12,621 | 570 | −26.69 | −24.77 | 4.87 | −12.15 | 4.65 |
| quadrant SW | 7,801 | 342 | +5.77 | **+6.41** | 11.21 | −27.17 | 6.45 |
| quadrant SE | 3,865 | 252 | +10.78 | **−45.95** | 10.73 | **+54.75** | 10.63 |
| quadrant NW | 5,604 | 244 | −18.49 | −13.85 | 6.15 | **+7.19** | 3.42 |
| quadrant NE | 7,017 | 341 | −33.23 | **−60.59** | 10.36 | −23.28 | 7.40 |
| checkerboard A (50 m) | 12,196 | 571 | −10.72 | −14.76 | 4.96 | −20.54 | 3.17 |
| checkerboard B (50 m) | 12,091 | 581 | −9.87 | −13.45 | 4.03 | −12.90 | 2.69 |

Twenty random 50%-of-blocks subsamples give dE −13.54 ± 3.55, dN −16.52 ± 2.05.

**The two checkerboard halves and the random subsamples agree with the whole, and the four
quadrants do not.** That distinction is the finding: interleaved subsamples reproduce the
fit because they sample the same field, while *spatially disjoint* subregions disagree
wildly – dE from +6.4 to −60.6, dN from +54.8 to −27.2, several of them many standard errors
apart. A single plane is not describing this surface; it is averaging regions that differ.

## 4. Is it landform? For the easting term, entirely

Nested models, cluster-robust at 50 m blocks (the optimistic rung; the ordering is what
matters here, not the significance):

| model | dE | dN | z (mm/100 m) | dist (mm/km) | ag (mm) |
|---|---|---|---|---|---|
| E, N | −14.19 ± 3.24 | −16.70 ± 2.10 | | | |
| z only | | | +27.59 ± 6.84 | | |
| dist only | | | | −123.14 ± 47.34 | |
| E, N, z | **−6.73 ± 3.03** | −17.35 ± 2.20 | +17.98 ± 6.48 | | |
| E, N, dist | −14.31 ± 3.25 | −16.22 ± 2.16 | | −37.64 ± 43.52 | |
| E, N, ag | −13.84 ± 3.26 | −16.73 ± 2.12 | | | +1.42 ± 2.69 |
| E, N, z, dist | **−3.83 ± 3.17** | −15.96 ± 2.13 | +26.01 ± 7.55 | −131.46 ± 51.41 | |
| E, N, z, dist, ag | **−3.99 ± 3.22** | −15.95 ± 2.13 | +26.36 ± 7.61 | −130.94 ± 51.40 | −1.24 ± 2.62 |

`dist` is the Euclidean distance from each cell to the nearest `floodplain_mask` cell
(median 92 m, range 5–233 m). `ag_region` covers 43% of the population and carries **no**
signal (+1.4 ± 2.7 mm); note also that `ag_region` is thresholded on `penetration`, which
this project has already shown is dominated by scan angle and flight-line overlap rather
than canopy, so it is a weak land-use covariate here in any case.

**Adding elevation removes half the easting term and adding valley distance removes three
quarters of it. Neither touches the northing term.**

Inside equal-count elevation bins (`binstats.quantile_edges`, spanning the full range, no
truncation):

| elevation | n | mean_mm | dE | dE_se | dN | dN_se |
|---|---|---|---|---|---|---|
| 217–223 m | 3,036 | −54.98 | −72.30 | 33.54 | −28.12 | 8.63 |
| 223–332 m | 3,036 | +11.76 | −27.67 | 7.50 | **+10.40** | 3.76 |
| 332–336 m | 3,036 | −1.79 | +1.66 | 4.43 | −21.38 | 2.58 |
| 336–338 m | 3,035 | −9.09 | **+18.26** | 5.70 | −22.68 | 2.88 |
| 338–339 m | 3,036 | −4.44 | **+16.53** | 4.06 | −21.37 | 2.92 |
| 339–341 m | 3,036 | −7.61 | −6.44 | 4.49 | −11.70 | 2.77 |
| 341–343 m | 3,036 | −4.81 | −17.39 | 5.22 | −6.34 | 3.94 |
| 343–352 m | 3,036 | −11.43 | −33.90 | 8.69 | −4.61 | 5.52 |

dE **changes sign** across the bands (−72 → +18 → −34). dN is negative in seven of eight.

Split at the antimode of the population's own elevation histogram:

| limb | n | mean_mm | dE | dE_se | dN | dN_se |
|---|---|---|---|---|---|---|
| upland z ≥ 270 m | 19,683 (81%) | −4.74 ± 1.44 | **+3.51 ± 2.61** | | **−17.90 ± 1.76** | |
| valley z < 270 m | 4,604 (19%) | −34.08 ± 5.94 | **−84.82 ± 10.84** | | −11.97 ± 4.34 | |

**The whole easting term lives in the 19% of cells on the valley floor**, where the DoD also
averages −34 mm rather than −5. Those are terrace and valley-margin cells that pass the
divide/curvature/slope test and lie outside `floodplain_mask`; they are a different landform
population, and they are concentrated at the tile's eastern side.

## 5. Is it the registration? The decisive within-swath test

`beam_offset_table.parquet` carries `point_source_id` and `gps_time` per return, so the
pipeline's own ground estimator (median of `d_mm_corr`, `ground_q = 0.50`) can be evaluated
one flight line at a time: **35,422 (cell, flight-line) estimates over the 24,287 cells**.
Inside one line the across-swath constant is fixed by construction, so any surviving `dE` is
not residual per-swath structure.

| fit | n | mean_mm | dE | dE_se | dN | dN_se |
|---|---|---|---|---|---|---|
| (cell,line) rows, no swath dummies | 35,422 | −11.06 | −11.47 | 2.90 | −15.01 | 1.87 |
| (cell,line) rows, **with swath dummies** | 35,422 | −8.72 | −10.24 | 3.66 | **−14.82** | 1.83 |

Swath fixed effects relative to line 135: 136 −4.08 ± 2.78, 137 +2.85 ± 4.77,
138 −13.69 ± 6.82 mm.

### 5a. Per flight line

| swath | rows | E_span (km) | h (m) | dir | dE | dE_se | dN | dN_se | mm/s of gps_time | c_tan |
|---|---|---|---|---|---|---|---|---|---|---|
| 135 | 4,466 | **0.33** | −2367 | south | +200.23 | 31.45 | −9.17 | 2.93 | **+1.087 ± 0.209** | −474 |
| 136 | 15,037 | 1.32 | +2595 | north | **−4.02** | 4.40 | −13.57 | 1.87 | **−1.107 ± 0.144** | −10 |
| 137 | 11,873 | 1.52 | −2602 | south | **+0.69** | 7.20 | −11.56 | 3.19 | **+0.958 ± 0.253** | −2 |
| 138 | 4,046 | **0.78** | +2525 | north | −159.87 | 24.68 | −21.48 | 6.28 | **−2.211 ± 0.569** | −404 |

`h` is fitted within each line from `x = a + b·y + h·tan(scan)` over all its in-grid returns
(r² 0.95–0.99). Its sign alternates because `scan_angle` is body-fixed and the aircraft turns
around; the magnitudes (2367–2602 m) reproduce the independent fit in
`SWATH_ACROSS_TRACK_TEST.md` §0 (−2424, +2564, −2626, +2473) to within 60 m. `dir` is
`corr(gps_time, northing)` within the line: −1.00, +1.00, −1.00, +1.00 – **consecutive lines
fly opposite directions**, which is what makes the drift test below work.

**dE.** The two lines that have any easting leverage – 136 and 137, together 79% of the rows,
with 1.32 and 1.52 km of easting span – give `−4.02 ± 4.40` and `+0.69 ± 7.20`. **Within a
swath, on the ground the swath actually covers, there is no easting tilt.** The two large
values belong to the lines the tile clips: 135 keeps only a 0.33 km sliver and 138 only
0.78 km, over which +200 and −160 mm/km amount to +66 and −125 mm total – the amplitude of
the patchy field itself, fitted over a third of its correlation length. The pooled within-swath
estimate is theirs, not the interior lines': the within estimator weights each line by
`n·Var(E)`, which is 0.008 / 0.484 / 0.447 / 0.061 for lines 135–138, and the weighted
average of the four per-line slopes is −9.70 against the −10.24 the swath-dummy model
returns – of which **line 138 alone supplies −9.70 of that, on 6% of the weight**, because its slope
is forty times the interior lines'. Re-expressed in each line's own body
frame (`c_tan = dE·h`, mm per unit tan θ) the four lines give −474, −10, −2, −404, which is
not one across-track law and does not match the +34 to +193 measured from the swath overlaps
in `SWATH_ACROSS_TRACK_TEST.md` – as expected, because in the overlap difference the spatial
field cancels exactly and in a DoD it does not.

**dN.** All four lines give a negative northing gradient, −9.2 to −21.5 mm/km, and the
pooled within-swath value (−14.82 ± 1.83) is within one standard error of the tile-wide one.
**Across-swath structure explains none of it.**

### 5b. Along-track drift is excluded by the flight direction

A residual along-track drift is a function of `gps_time` within a line. Because consecutive
lines fly opposite directions, such a residual **must alternate its sign in northing** line
to line. What is observed is the opposite pattern: the slope in `gps_time` alternates
(+1.09, −1.11, +0.96, −2.21 mm/s, each many standard errors from zero) while the slope in
**northing does not** (−9.2, −13.6, −11.6, −21.5, all negative). The gradient is fixed in
ground coordinates, not in mission time. That is a drift-shaped signature only in the sense
that a ground-fixed north-south field *projects onto* `gps_time` with the sign of whichever
way the aircraft was flying.

### 5c. What each registration term is worth to the tilt (250 m blocks)

| term subtracted back out | mean_mm | dE | dE_se | dN | dN_se |
|---|---|---|---|---|---|
| nothing (as registered) | −10.30 | −14.19 | 8.39 | −16.70 | 5.06 |
| along-track drift | +5.35 | −14.23 | 10.28 | **−14.64** | 5.99 |
| per-swath constants | −36.16 | **−24.25** | 8.82 | −15.11 | 5.23 |

The per-swath constants are **doing their job**: without them the easting gradient is −24.25
mm/km, so they remove about 10 mm/km of it, and the gradient implied by the constants
themselves – 0, −23.9, −32.5, −43.7 mm over the 2.16 km spanned by the four lines' mean
eastings, i.e. about −20 mm/km – is of the right size and sign. The
drift curves move the northing gradient by only 2 mm/km, in either direction – consistent
with §5b: they are not making `dN`, and turning them off does not remove it.

### 5d. Upland limb, per line (landform and swath controlled together)

| fit | n | mean_mm | dE | dE_se | dN | dN_se |
|---|---|---|---|---|---|---|
| upland, all lines | 29,158 | −7.16 | **+3.22 ± 2.43** | | **−18.01 ± 1.62** | |
| upland, swath 135 | 4,466 | +4.25 | +200.23 | 31.45 | −9.17 | 2.93 |
| upland, swath 136 | 14,556 | −7.30 | −6.76 | 4.37 | −15.93 | 1.89 |
| upland, swath 137 | 8,746 | −14.02 | +24.67 | 5.86 | −27.02 | 3.41 |
| upland, swath 138 | 1,390 | +0.97 | −85.98 | 43.28 | −19.53 | 11.89 |

## 6. What is left

**dE: nothing.** On the upland limb, at every block size:

| block_m | dE | dE_se | t | dN | dN_se | t |
|---|---|---|---|---|---|---|
| 50 | +3.51 | 2.61 | +1.34 | −17.90 | 1.76 | −10.19 |
| 100 | +3.51 | 3.93 | +0.89 | −17.90 | 2.50 | −7.16 |
| 250 | +3.51 | 6.55 | +0.54 | −17.90 | 3.90 | −4.59 |
| 500 | +3.51 | 7.32 | +0.48 | −17.90 | 4.33 | −4.13 |
| 1000 | +3.51 | 9.18 | +0.38 | −17.90 | 4.74 | −3.78 |

**dN: about −18 mm/km, but it is not a ramp.** Median DoD by 500 m northing band on the
upland limb:

| band (m north of the tile's southern edge) | n | median (mm) | mean (mm) |
|---|---|---|---|
| 0–500 | 2,159 | +11.83 | +11.79 |
| 500–1000 | 3,650 | **+33.71** | +30.65 |
| 1000–1500 | 3,037 | +6.28 | +1.12 |
| 1500–2000 | 4,551 | **−30.45** | −33.77 |
| 2000–2500 | 4,237 | −6.70 | −7.82 |
| 2500–3000 | 93 | +1.35 | −2.24 |
| 3000–3500 | 1,956 | −21.57 | −23.99 |

Two patches of opposite sign at ~1 km scale and ±30 mm amplitude, with the fitted line drawn
through them. That is the same object the semivariogram describes: a field correlated to
400 m and anticorrelated past 1.6 km. **Quoting `dN = −16.7 mm/km` as a tilt over-describes
it; the honest statement is a patchy residual of about ±30 mm at 0.5–1.5 km scale whose
linear projection in northing is −18 ± 4 mm/km (250 m blocks, upland limb).**

Two things it is *not*, both tested here rather than assumed: it is not the geoid datum
(`registration.geoid_term` applies the tilt as well as the constant, and the northing
component of that tilt, `tilt_c = −0.57 mm/km` in `corrections.json`, is thirty times too
small in any case), and it is not the gen1 along-track drift (§5b). Whether the remaining field is
gen2-side structure, gen1-side structure the four terms do not model, or real ground change
on the upland is **not decided by this run** and would need the same test on `elbaext`,
whose flight lines cover different ground.

## What this changes

1. **Do not carry `dE = −14.19 mm/km` forward.** It is a landform confound, and it is not
   present within a flight line or on the upland limb. Anything built on it – a fitted
   across-tile ramp, a "residual swath structure" term – would be fitting the valley-floor
   terraces.
2. **The per-swath constants and the drift curves are exonerated by this test.** The
   constants remove roughly 10 mm/km of real easting gradient; the drift curves account for
   about 2 mm/km of the northing one and, decisively, the remaining northing gradient is
   fixed in ground coordinates rather than in mission time.
3. **Cluster-robust errors at 50 m are not conservative enough for tile-scale coefficients.**
   The correlation length of this residual field is at least 400 m. For anything fitted over
   the whole tile, quote the 250–500 m rung.
4. **The audited population is 81% upland plateau and 19% valley terrace, with the hillslopes
   between missing.** Every tile-scale statement made on it inherits that split; report the
   two limbs separately, or state which one the number belongs to.

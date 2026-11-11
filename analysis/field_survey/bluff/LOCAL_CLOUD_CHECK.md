STATUS: COMPLETE

# Bluff RTK points against the POINT CLOUDS, not the 5 m grids

Script: `analysis/field_survey/bluff/local_cloud_check_rtk.py`
Outputs: `data/derived/4342-29-64/bluff_check/local/`
Inputs: RTK `data/field/2026-10-25/rtk_points_all.csv` (status FIX*, group hillslope*, 39 rows);
gen2 class-2 `data/after/valley/4342-29-64_3dep_class2.laz`; gen1 raw `data/before/4342-29-64.laz`;
gen1 CSF ground `data/csf_cache/4342-29-64.las`; shifts `data/derived/4342-29-64/corrections.json`.

## Sign convention and what the product actually applied (verified in code)

- `chain.StarFrame.apply`: `ctx.x[m] += dx; ctx.y[m] += dy; ctx.z[m] += dz` per point_source_id,
  with `per_swath_star_frame_dxdydz_m`.
- `chain.LateralShift.apply`: `ctx.x += a0; ctx.y += b0` with `cross_epoch_datum.horizontal_shift_m
  = [-0.0373, -0.0338]` (order 0, one constant).
- Under `swath_frame = "star"` (this product) `per_swath_internal_alignment_dxdydz_m` is SOLVED
  BUT NOT APPLIED (`pipeline.py` line 1237: `apply_swath_corr=(swath_frame == "chain")`;
  `registration.read_swath_alignment` refuses that key for a star product). So the registered
  gen1 position used here is raw + star + lateral. The "internal then star" composition named in
  the task brief is computed too, as a labeled variant, because it is NOT what the product applied.
- The CSF cache file has the same x,y,z minima as the raw file, so it holds raw (unshifted)
  positions; it is shifted here with the same per-swath star + lateral terms.

## Method (step 2)

For each RTK point and radius r (3 m = Andy's proposal; 5 m also run), the cloud points within r of
the tied E,N get a least-squares plane z = a + b(x-x0) + c(y-y0); the residual is navd88_tied - a.
Plane slope = atan(hypot(b,c)); downslope aspect = compass azimuth of -(b,c). A fit needs >= 3 points
(the brief's CSF minimum, applied to every population: 3 points define a plane). gen2 z and
navd88_tied_m are compared as they stand, with NO datum constant applied, exactly as the grid columns
tied_minus_gen2_m were made, so the two are directly comparable. Populations: gen2 (2021 class 2);
gen1_raw (2008 last returns, raw positions); gen1_reg (same at raw + star + lateral, what the product
applied); gen1_is (VARIANT raw + internal + star + lateral, the composition named in the brief, NOT
applied by the product); csf_reg (CSF-selected ground, registered).

PROPOSAL added after the first run (labeled in the tables): a 3-parameter plane on 3-5 CSF points
that lie nearly along one line, or on 12-21 gen1 last returns of which some are canopy single
returns, extrapolates to nonsense at the RTK x,y (-267 m, -11 m in the per-point table). So the
per-point gen1 offset is ALSO measured with the 2021 plane's slope borrowed: height of each gen1
return above the gen2 plane at the return's own x,y, median over the disk ("h>gen2 median"), and
residual = navd88_tied - (gen2 plane + that median). That is the same quantity step 3 reports for
the whole population. The plane fits the brief asked for are kept and reported.

## Script output (pasted verbatim from data/derived/4342-29-64/bluff_check/local/results.md)

```
RTK rows: 39 (status FIX*, group hillslope*)
RTK bbox E 579743.1-579843.1, N 4885355.5-4885648.7; read pad 14.0 m
swath_frame=star; star shifts (ADDED to x,y,z) {'135': [-0.1695, 0.0798, 0.0806], '136': [-0.1281, -0.0733, 0.0727], '137': [-0.1106, -0.1265, 0.066], '138': [0.1579, 0.084, 0.0635]}
internal alignment (NOT applied under star) {'135': [0.0, 0.0, 0.0], '136': [0.3176, -0.0821, -0.0179], '137': [0.5906, -0.1707, -0.026], '138': [1.0359, 0.0964, -0.0364]}
lateral shift (ADDED to x,y) [-0.0373, -0.0338]
gen2 points in bbox: 242,852; classification values [2]
gen2 class-2 points in bbox: 242,852
gen1 raw returns in bbox: 27,379 of 7,728,747; point_source_ids [138]; classification values [0, 2, 5, 6, 8]
gen1 last returns in bbox: 24,160 (88.2%); single returns (nr==1): 74.6% of all returns
gen1 CSF ground points in bbox: 3,275; point_source_ids [138]
CSF membership of gen1 last returns in bbox: 3,275 of 24,160 (13.6%) matched by integer X,Y,Z

Wrote /home/awickert/projects/lidar-diff-icp/data/derived/4342-29-64/bluff_check/local/per_point_planes.csv (78 rows = 39 points x 2 radii)

## Step 2. Residual = navd88_tied - local LS plane, mm (median, NMAD, mean, SD; n = points with a fit)

### r = 3 m
| population / group | n | median | NMAD | mean | SD |
|---|---|---|---|---|---|
| gen2, ALL | 39 | 35.1 | 135.1 | -18.9 | 286.1 |
| gen2, hillslope_RS3 | 18 | 6.2 | 180.7 | -115.8 | 382.8 |
| gen2, hillslope_T4a | 21 | 43.7 | 135.4 | 64.2 | 123.5 |
| gen1_raw, ALL | 39 | 76.0 | 488.6 | -604.1 | 2124.0 |
| gen1_raw, hillslope_RS3 | 18 | 48.5 | 375.5 | -127.1 | 635.6 |
| gen1_raw, hillslope_T4a | 21 | 82.1 | 516.6 | -1012.9 | 2801.5 |
| gen1_reg, ALL | 39 | 8.1 | 357.8 | -690.7 | 2006.2 |
| gen1_reg, hillslope_RS3 | 18 | -52.9 | 315.3 | -267.4 | 667.2 |
| gen1_reg, hillslope_T4a | 21 | 49.4 | 362.7 | -1053.6 | 2639.9 |
| gen1_is, ALL | 39 | -462.0 | 548.9 | -997.9 | 1512.2 |
| gen1_is, hillslope_RS3 | 18 | -464.9 | 304.4 | -740.8 | 873.4 |
| gen1_is, hillslope_T4a | 21 | -462.0 | 632.3 | -1218.2 | 1893.7 |
| csf_reg, ALL | 6 | 4.3 | 64855.3 | -41506.5 | 126398.3 |
| csf_reg, hillslope_RS3 | 2 | -177342.5 | 133336.3 | -177342.5 | 127186.1 |
| csf_reg, hillslope_T4a | 4 | 44.1 | 59.0 | 26411.5 | 52763.9 |

Point counts within r=3 m (median [min-max] over the 39 points): gen2 185 [81-239]; gen1_raw 17 [12-21]; gen1_reg 17 [12-21]; gen1_is 16 [12-20]; csf_reg 0 [0-16]
Plane-fit RMS about the plane, m (median over points): gen2 0.049; gen1_raw 0.118; gen1_reg 0.115; gen1_is 0.123; csf_reg 0.014
gen2 plane slope deg (median [min-max]): 31.2 [7.0-36.3]; downslope aspect deg (median): 99; grid slope_deg column (median): 32.2
gen1_reg plane slope deg (median [min-max]): 32.9 [6.0-79.1]; downslope aspect deg (median): 94
gen2 cloud residual MINUS grid-based tied_minus_gen2_m: median -7.9 mm, NMAD 78.7 mm, n 39; grid column itself: median 34.7 mm, NMAD 138.8 mm
Points where the gen1 GRID has a value (n=6): grid tied_minus_gen1 median 149.0 mm; cloud csf_reg residual at those points median 44.1 mm (n fit 4); gen1_reg residual median -215.0 mm
gen1_raw plane minus gen2 plane at the RTK x,y: median -197.2 mm, NMAD 145.7 mm, n 39
gen1_reg plane minus gen2 plane at the RTK x,y: median -97.5 mm, NMAD 76.9 mm, n 39
gen1_is plane minus gen2 plane at the RTK x,y: median 504.0 mm, NMAD 533.0 mm, n 39
csf_reg plane minus gen2 plane at the RTK x,y: median -66.0 mm, NMAD 64757.5 mm, n 6
PROPOSAL (gen2-slope offset): per-point MEDIAN height of gen1 returns above the gen2 plane, then residual = navd88_tied - (gen2 plane + that median)
  gen1_reg: median-height-above-gen2 over points: median -113.3 mm, NMAD 50.6 mm, n 39; p25-height median over points -173.7 mm; residual RTK minus (gen2 plane + median height): median 129.5 mm, NMAD 175.0 mm
  csf_reg: median-height-above-gen2 over points: median -200.3 mm, NMAD 50.5 mm, n 8; p25-height median over points -203.5 mm; residual RTK minus (gen2 plane + median height): median 93.3 mm, NMAD 119.2 mm

### r = 5 m
| population / group | n | median | NMAD | mean | SD |
|---|---|---|---|---|---|
| gen2, ALL | 39 | 62.9 | 136.9 | -5.6 | 297.4 |
| gen2, hillslope_RS3 | 18 | 14.0 | 228.5 | -108.8 | 382.4 |
| gen2, hillslope_T4a | 21 | 88.2 | 110.4 | 82.9 | 161.0 |
| gen1_raw, ALL | 39 | -54.0 | 730.0 | -707.9 | 1788.7 |
| gen1_raw, hillslope_RS3 | 18 | -94.5 | 668.8 | -368.3 | 696.4 |
| gen1_raw, hillslope_T4a | 21 | 59.3 | 742.1 | -999.0 | 2339.6 |
| gen1_reg, ALL | 39 | -152.9 | 693.0 | -831.0 | 1741.6 |
| gen1_reg, hillslope_RS3 | 18 | -210.8 | 795.6 | -506.6 | 693.3 |
| gen1_reg, hillslope_T4a | 21 | -50.3 | 826.6 | -1109.0 | 2275.6 |
| gen1_is, ALL | 39 | -816.9 | 921.8 | -1253.7 | 1615.9 |
| gen1_is, hillslope_RS3 | 18 | -761.6 | 720.0 | -920.2 | 828.5 |
| gen1_is, hillslope_T4a | 21 | -945.4 | 1112.3 | -1539.7 | 2047.4 |
| csf_reg, ALL | 7 | -27.9 | 4017.0 | 681.7 | 51657.6 |
| csf_reg, hillslope_RS3 | 3 | -43940.0 | 15157.4 | -33613.6 | 27223.8 |
| csf_reg, hillslope_T4a | 4 | 41.7 | 61.7 | 26403.2 | 52769.5 |

Point counts within r=5 m (median [min-max] over the 39 points): gen2 492 [257-675]; gen1_raw 47 [42-52]; gen1_reg 48 [43-56]; gen1_is 46 [41-53]; csf_reg 0 [0-49]
Plane-fit RMS about the plane, m (median over points): gen2 0.088; gen1_raw 1.495; gen1_reg 1.444; gen1_is 1.672; csf_reg 0.037
gen2 plane slope deg (median [min-max]): 31.9 [5.5-36.9]; downslope aspect deg (median): 100; grid slope_deg column (median): 32.2
gen1_reg plane slope deg (median [min-max]): 34.1 [8.3-68.3]; downslope aspect deg (median): 91
gen2 cloud residual MINUS grid-based tied_minus_gen2_m: median 10.6 mm, NMAD 60.4 mm, n 39; grid column itself: median 34.7 mm, NMAD 138.8 mm
Points where the gen1 GRID has a value (n=6): grid tied_minus_gen1 median 149.0 mm; cloud csf_reg residual at those points median 28.1 mm (n fit 5); gen1_reg residual median -1072.1 mm
gen1_raw plane minus gen2 plane at the RTK x,y: median 30.4 mm, NMAD 471.6 mm, n 39
gen1_reg plane minus gen2 plane at the RTK x,y: median 134.6 mm, NMAD 423.2 mm, n 39
gen1_is plane minus gen2 plane at the RTK x,y: median 695.0 mm, NMAD 787.4 mm, n 39
csf_reg plane minus gen2 plane at the RTK x,y: median -56.0 mm, NMAD 4259.8 mm, n 7
PROPOSAL (gen2-slope offset): per-point MEDIAN height of gen1 returns above the gen2 plane, then residual = navd88_tied - (gen2 plane + that median)
  gen1_reg: median-height-above-gen2 over points: median -122.2 mm, NMAD 50.3 mm, n 39; p25-height median over points -233.2 mm; residual RTK minus (gen2 plane + median height): median 116.5 mm, NMAD 170.7 mm
  csf_reg: median-height-above-gen2 over points: median -152.2 mm, NMAD 131.3 mm, n 15; p25-height median over points -162.0 mm; residual RTK minus (gen2 plane + median height): median 132.7 mm, NMAD 130.6 mm

### Per-point table, r = 3 m (mm; n = points in the fit)
| group | pt | gen2 n | gen2 resid | grid resid | gen2 slope | gen1_reg n | gen1_reg plane resid | gen1_reg-gen2 plane | gen1_reg h>gen2 median | gen1_reg h>gen2 p25 | csf n | csf plane resid | csf h>gen2 median |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| RS3 | 7 | 142 | 123 | 172 | 27.9 | 20 | -1110 | 1234 | 31 | -15 | 0 | - | - |
| RS3 | 8 | 173 | -53 | 3 | 27.7 | 17 | 96 | -149 | -183 | -295 | 1 | - | -230 |
| RS3 | 9 | 198 | 115 | 181 | 30.9 | 18 | -1757 | 1872 | 1 | -92 | 0 | - | - |
| RS3 | 10 | 154 | -1463 | -1478 | 29.6 | 15 | -1697 | 234 | -89 | -130 | 0 | - | - |
| RS3 | 11 | 138 | 76 | 23 | 28.0 | 17 | -430 | 507 | -92 | -146 | 0 | - | - |
| RS3 | 12 | 164 | 133 | 28 | 29.2 | 17 | 282 | -149 | -142 | -204 | 0 | - | - |
| RS3 | 13 | 156 | -5 | -12 | 34.9 | 16 | 96 | -100 | -85 | -157 | 0 | - | - |
| RS3 | 14 | 165 | 73 | 14 | 35.7 | 18 | 170 | -98 | -72 | -158 | 0 | - | - |
| RS3 | 15 | 180 | 33 | 70 | 34.2 | 17 | 141 | -108 | -133 | -163 | 0 | - | - |
| RS3 | 16 | 198 | 17 | 70 | 31.2 | 14 | 149 | -132 | -123 | -194 | 0 | - | - |
| RS3 | 17 | 176 | 35 | 43 | 34.6 | 18 | 249 | -214 | -191 | -266 | 0 | - | - |
| RS3 | 18 | 167 | -226 | -260 | 30.9 | 14 | -77 | -149 | -155 | -189 | 0 | - | - |
| RS3 | 19 | 238 | -196 | -151 | 32.0 | 15 | -68 | -128 | -122 | -182 | 0 | - | - |
| RS3 | 20 | 220 | -434 | -434 | 31.1 | 16 | -996 | 561 | -119 | -155 | 0 | - | - |
| RS3 | 21 | 227 | -209 | -297 | 32.9 | 18 | -74 | -136 | -151 | -193 | 0 | - | - |
| RS3 | 22 | 222 | -274 | -235 | 32.8 | 19 | -191 | -84 | -100 | -139 | 0 | - | - |
| RS3 | 23 | 225 | -158 | 15 | 33.3 | 19 | -38 | -120 | -145 | -268 | 5 | -87408 | -227 |
| RS3 | 24 | 234 | 329 | 556 | 31.9 | 19 | 441 | -111 | -136 | -232 | 5 | -267277 | -238 |
| T4a | 9 | 236 | -26 | -5 | 18.9 | 17 | 49 | -76 | -70 | -104 | 7 | 80 | -70 |
| T4a | 10 | 94 | 204 | 199 | 28.9 | 15 | -3017 | 3221 | 3193 | -49 | 3 | 105557 | -177 |
| T4a | 11 | 157 | -106 | -44 | 27.5 | 17 | 86 | -192 | -224 | -288 | 1 | - | -224 |
| T4a | 12 | 154 | 144 | 124 | 31.2 | 18 | -1527 | 1671 | -118 | -222 | 0 | - | - |
| T4a | 13 | 159 | 42 | 73 | 34.0 | 14 | 120 | -78 | -82 | -160 | 0 | - | - |
| T4a | 14 | 116 | 24 | 149 | 28.8 | 15 | -4480 | 4504 | 61 | -55 | 0 | - | - |
| T4a | 15 | 179 | 11 | 79 | 25.2 | 14 | -1397 | 1408 | -118 | -124 | 0 | - | - |
| T4a | 16 | 204 | 135 | 84 | 32.5 | 19 | -108 | 244 | -112 | -290 | 0 | - | - |
| T4a | 17 | 168 | 65 | 42 | 32.2 | 17 | -1710 | 1774 | -77 | -159 | 0 | - | - |
| T4a | 18 | 185 | 217 | 237 | 34.7 | 21 | 298 | -81 | -40 | -174 | 0 | - | - |
| T4a | 19 | 203 | 97 | 296 | 32.7 | 18 | 265 | -168 | -159 | -259 | 0 | - | - |
| T4a | 20 | 220 | 40 | 134 | 33.9 | 16 | 139 | -99 | -101 | -200 | 0 | - | - |
| T4a | 21 | 212 | 135 | 86 | 30.6 | 17 | 305 | -170 | -147 | -243 | 0 | - | - |
| T4a | 22 | 223 | 119 | 18 | 35.4 | 15 | 294 | -175 | -157 | -223 | 0 | - | - |
| T4a | 23 | 222 | 414 | 585 | 30.0 | 16 | 558 | -144 | -113 | -265 | 0 | - | - |
| T4a | 24 | 239 | -67 | -73 | 13.5 | 15 | 8 | -75 | -79 | -144 | 15 | 8 | -79 |
| T4a | 25 | 214 | -62 | -156 | 35.5 | 15 | 119 | -181 | -160 | -248 | 0 | - | - |
| T4a | 26 | 216 | 44 | 196 | 36.3 | 13 | -1220 | 1264 | -97 | -170 | 0 | - | - |
| T4a | 27 | 209 | 80 | 35 | 32.6 | 15 | 229 | -149 | -126 | -177 | 0 | - | - |
| T4a | 28 | 127 | -105 | -199 | 28.1 | 12 | -11138 | 11034 | 9822 | -43 | 0 | - | - |
| T4a | 29 | 81 | -56 | -59 | 7.0 | 16 | 1 | -57 | -71 | -107 | 16 | 1 | -71 |

## Step 3. 2008 last returns within 5 m of any RTK point: height above the LOCAL gen2 plane (r=3 m fit at the nearest RTK point)
| population | n | p5 | p25 | p50 | p75 | p95 |  (mm above gen2 plane)
|---|---|---|---|---|---|---|
| gen1 last returns, REGISTERED positions | 1672 | -431 | -224 | -119 | 5 | 7978 |
|   ...of which CSF-selected (matched) | 123 | -275 | -165 | -101 | -2 | 123 |
|   ...NOT CSF-selected | 1549 | -443 | -228 | -121 | 5 | 8253 |
| gen1 last returns, RAW positions | 1665 | -554 | -349 | -236 | -117 | 7842 |
| gen1 last returns, internal+star VARIANT (not applied) | 1625 | -188 | 97 | 325 | 518 | 8565 |
| CSF ground cloud, registered (direct) | 123 | -275 | -165 | -101 | -2 | 123 |
| gen2 class-2 points about the same planes (scale) | 17169 | -230 | -68 | -10 | 40 | 168 |
Fractions of registered gen1 last returns (n=1672): below gen2 plane by >0.1 m 55.3%, within +/-0.1 m 28.9%, above by >0.1 m 15.8%, above by >1 m 10.2%
CSF-selected share of those returns: 7.4%; CSF-selected median -101 mm vs not-selected median -121 mm
All 2008 returns within 5 m (raw positions): n 1,885; number_of_returns histogram {1: 1404, 2: 472, 3: 9}; last-return share 88.3%; last returns per RTK 5 m disk: median 48
Wrote /home/awickert/projects/lidar-diff-icp/data/derived/4342-29-64/bluff_check/local/gen1_last_returns_near_rtk.csv

## Step 4. gen2 cloud residual (r=3 m) under horizontal shifts of the RTK positions, mm
| shift | dE: median | dE: NMAD | dN: median | dN: NMAD |
|---|---|---|---|---|
| -0.5 | -187.5 | 142.4 | 13.3 | 262.7 |
| -0.3 | -119.6 | 114.8 | 7.0 | 191.4 |
| -0.2 | -71.1 | 125.0 | -2.2 | 167.3 |
| -0.1 | -24.0 | 133.2 | 21.1 | 152.8 |
| +0.0 | 35.1 | 135.1 | 35.1 | 135.1 |
| +0.1 | 96.2 | 137.1 | 26.7 | 122.2 |
| +0.2 | 137.3 | 147.4 | 37.4 | 102.7 |
| +0.3 | 193.4 | 160.6 | 37.3 | 117.8 |
| +0.5 | 277.9 | 212.7 | 49.5 | 177.7 |
Sensitivity of the median residual to dE: 488 mm per m (LS slope over the 9 shifts)
Sensitivity of the median residual to dN: 46 mm per m (LS slope over the 9 shifts)
Expected from the median gen2 plane (slope 31.2 deg, downslope aspect 99 deg): d(resid)/dE = +tan(slope)*sin(aspect) = 598 mm/m, d(resid)/dN = +tan(slope)*cos(aspect) = -94 mm/m (moving the RTK point downslope lowers the plane under it, so the residual rises)
```

Files: `per_point_planes.csv` (39 points x 2 radii, every column above plus plane heights, slopes,
aspects, RMS), `gen1_last_returns_near_rtk.csv` (the 1,672 registered 2008 last returns within 5 m
of an RTK point with their height above the local 2021 plane and CSF membership),
`shift_sensitivity.csv`, `run.log`, `results.md`.

## Step 5. Reading of the numbers

**Does the 2021 cloud match the surveyed ground?** In the median, yes: RTK minus the local 2021
plane is +35.1 mm (r = 3 m, n = 39) and +62.9 mm (r = 5 m), against +34.7 mm from the 5 m grid. The
cloud and the grid give the same number per point: cloud minus grid is -7.9 mm median, 78.7 mm NMAD
(r = 3 m) and +10.6 mm, 60.4 mm NMAD (r = 5 m). Gridding, the slope-normal median and bilinear
sampling are therefore not where the scatter comes from: the per-point scatter is 135.1 mm NMAD
against the cloud and 138.8 mm against the grid, and the two groups persist in the cloud
(hillslope_RS3 +6.2 mm, hillslope_T4a +43.7 mm at r = 3 m; +14.0 and +88.2 mm at r = 5 m). The
biggest residuals are the same points in both (RS3 point 10: -1463 cloud / -1478 grid; RS3 24: +329
/ +556; T4a 23: +414 / +585), so they are properties of the RTK position or the ground there, not
of the product. The 2021 cloud itself is smooth at this scale: its own points scatter 0.049 m RMS
about the 3 m planes (0.088 m at 5 m), with a plane slope of 31.2 deg [7.0-36.3] and downslope
aspect 99 deg, matching the grid slope column (32.2 deg) and the ~99 deg aspect the brief quotes from the grid.

**With what uncertainty.** The per-point agreement is 135 mm NMAD, of which the 2021 surface's own
roughness accounts for ~50 mm RMS. The rest is the horizontal degeneracy of step 4: the median
residual moves 488 mm per meter of east shift (LS slope over the nine shifts; the median plane
predicts 598 mm/m = tan 31.2 deg x sin 99 deg) and 46 mm per meter of north shift (plane predicts
-94 mm/m, the measured north response is within its own NMAD and not resolved). A 135 mm vertical
NMAD is therefore equivalent to ~0.23-0.28 m of across-slope (east-west) position scatter. The RTK
rms_e column is ~0.015 m, so either the real horizontal error under canopy on a 32 deg slope is an
order of magnitude larger than the receiver's figure, or ~0.1 m vertical disagreement is real
microtopography that neither a 3 m plane nor the lidar ground resolves. This check cannot separate
the two. What it does fix: the median level, +35 mm (r = 3 m), moves by only +/-25 mm for a 0.1 m
east error (-24.0 / +96.2 mm at -0.1 / +0.1 m), so the group medians are robust at the few-cm
level and the per-point values are not.

**What the 2008 returns are doing there.** At registered positions, the 1,672 last returns within 5
m of an RTK point sit BELOW the 2021 surface as a body: p5/p25/p50/p75/p95 = -431/-224/-119/+5/+7978
mm. 55.3% are more than 0.1 m below the 2021 plane, 28.9% within +/-0.1 m, 15.8% more than 0.1 m
above and 10.2% more than 1 m above (canopy: the bluff's 2008 returns are 74.6% single returns and
88.2% last returns, 1404 of 1885 returns within 5 m have number_of_returns = 1). For scale, the
2021 class-2 points about the same planes are -230/-68/-10/+40/+168 mm. So the 2008 last-return
population is centered ~0.12 m below the 2021 ground with about twice the spread, and the CSF
subset is NOT its low tail: the 123 CSF-selected returns (7.4% of the population) are
-275/-165/-101/-2/+123 mm, a median 20 mm HIGHER than the non-selected returns (-121 mm) with the
canopy tail removed. The per-point version agrees: the median 2008 height above the 2021 plane is
-113.3 mm (NMAD 50.6 mm over the 39 points, r = 3 m) for all last returns and -200.3 mm for the 8
points that have CSF returns. The plane-to-plane difference at the RTK x,y says the same thing
where the 2008 plane is clean: gen1_reg plane minus gen2 plane -97.5 mm median, 76.9 mm NMAD (r =
3 m). RTK minus the 2008 level (gen2 plane + 2008 median height) is +129.5 mm median, 175.0 mm NMAD,
against +149.0 mm from the six grid cells; the brief's "0.4 m below at 30 deg" is T4a points 10
and 15 (grid +456 / +426 mm), where the cloud gives 2008 medians of +3193 mm (point 10: 15 returns,
most in canopy, p25 -49 mm) and -118 mm (point 15), so those two grid values are a 3-return CSF
cell and a canopy-dominated disk, not a 0.4 m ground offset.

**The registration check that fell out.** The composition named in the brief (internal alignment,
then star) puts the 2008 returns +504.0 mm ABOVE the 2021 plane (median, NMAD 533.0 mm) and the
2008 population at +325 mm p50, because the internal dx for swath 138 (+1.036 m east) moves every
return 1 m downslope on a 32 deg east-facing slope. The product did not apply it (verified in
`pipeline.py` and `registration.py`), and the star + lateral composition is the tightest of the
three (gen1_reg minus gen2 plane NMAD 76.9 mm against 145.7 mm raw and 533.0 mm with the internal
term). All 2008 returns here are from one swath, point_source_id 138; the 2008 class codes present
are 0, 2, 5, 6, 8.

**What this check could not determine.** (1) Whether the 135 mm per-point scatter is RTK position
error under canopy or real sub-3-m ground form: both produce the same signature and the check has
no independent position. (2) The absolute sign of the +35 mm: the floodplain survey's ~0.1 m datum
constant between navd88_tied_m and gen2 z is not applied here, so the median is in the same frame
as the grid columns and no more absolute than they are. (3) Why hillslope_RS3 and hillslope_T4a
differ by ~40-70 mm in the median; the cloud reproduces the split, so it is not a gridding
artifact; whether the two surveys occupied any common points, which would test a per-survey tie, was not checked here. (4) Whether the 2008
population's -0.12 m is a 2008 ranging/penetration bias on the slope or 13 years of real change:
the 2021-2008 difference at the bluff is the same sign and size as the canonical gen1-reads-low-on-
slopes signal, and this survey cannot tell them apart because it has no 2008-contemporaneous ground.
(5) The north-shift sensitivity: the measured 46 mm/m sits inside the NMAD of the shifted residuals
and disagrees in sign with the plane's -94 mm/m, so only the east response is resolved.

Offers, not executed: (a) regress the per-point gen2 residual on the RTK rms_e, rms_n, tilt_deg and
antenna_m columns to see whether the receiver's own quality flags predict the scatter; (b) repeat
step 4 along the local downslope direction per point instead of fixed E/N axes.

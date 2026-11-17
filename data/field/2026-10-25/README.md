# RTK-GNSS ground survey, Whitewater WMA north of Elba, 25 October 2026

ESCI 4701/8701 Geomorphology field trip. Emlid Reach rovers; one on MnDOT CORS (Andy
Wickert), the rest on a class base station that was set twice. Compiled 2026-11-11 by
`analysis/field_survey/compile_rtk_2026-10-25.py`; every number here is pasted from its
log (`compile_log.txt`).

## The table: `rtk_points_all.csv` (220 rows); the fitted table `rtk_points_all_fitted.csv`

One row per logged point, nothing filtered. Columns: `group`, `people`, `device`, `serial`,
`point`, `t_start` (CDT), `base` (B1, B or CORS), `baseline_m`, `status` (FIX / FLOAT /
SINGLE as exported), `lon`, `lat` (NAD83(2011)), `easting`, `northing` (NAD83(2011) / UTM 15N, EPSG:6344, projection only; see Frame),
`h_ell_m` (ellipsoid height as exported), **`h_ell_tied_m`** (after the base tie below),
`tie_source`, `easting_tied`, `northing_tied` (B shifted onto the monument, see the tie; the lidar is read here), `navd88_m`, **`navd88_tied_m`** (GEOID18 through pyproj, checked at DR9925:
218.748 vs published 218.75), `elev_app_m` (the Flow app's own NAVD88 where exported),
`antenna_m` (as exported), `antenna_corr_m` (the pole-height correction, see the check), RMS columns, `tile`, `gen2_z_m` (2021 lidar ground, the horizontal per-cell median grid, bilinear from 5 m;
the point-cloud check reproduces it at the survey points to -8 mm), `gen2_zref_m` (the
slope-normal 2021 grid the difference was built on), `gen1_z_m` (2008 lidar ground =
`gen2_zref_m` minus the tile's difference; until 2026-11-11 it was taken from `gen2_z_m`, which
is wrong by up to 0.13 m at single bluff points and by 8 mm in the group medians), `dod_2021_2008_m`,
`rtk_minus_gen2_m`, `rtk_minus_gen1_m` (untied), **`tied_minus_gen2_m`,
`tied_minus_gen1_m`** (RTK minus lidar), **`gen2_minus_rtk_m`, `gen1_minus_rtk_m`** (their negatives:
the LIDAR misfit re-referenced to the surveyed ground, positive = lidar above the ground; the
figures use these), `slope_deg`, `lowveg_gen2`, `clubb_floodplain`,
`penetration_failure_2008`, `source`.

Use FIX rows and the `tied_` columns. Groups:

`rtk_points_all_fitted.csv` (same columns, 2026-11-11) applies the fitted offset from the
weighted inversion, dE -0.05, dN -0.05, dZ +0.13 m (an EMPIRICAL offset between the May 2021
ground class and the pole tips on vegetated covers, not a datum correction: see "Is the fitted
offset vegetation?" below), to every point (B and CORS alike), with every
lidar layer re-read at the shifted position by the same compile code; its `gen2_minus_rtk_m` is the
residual after that shift, and `tie_source` names the shift. The main table is unchanged.

| group | people | device | n | where |
|---|---|---|---|---|
| hillslope_RS3 | T3a, T3b, T3c | Reach RS3 | 18 (24 logged; 6 under B1 dropped) | bluff west of the road, up and north |
| hillslope_T4a | T4a, T4b, T3a, T3b, T3c | Reach RS3 (a second unit) | 21 (29 logged; 8 under B1 dropped) | same bluff, parallel traverse (Flow export from T4b's Canvas submission, received 2026-11-12; it reproduces the F07 .docx table to the millimeter) |
| floodplain_RSplus_A | T1a, T1b | Reach RS+ | 62 (63 logged; 1 under B1 dropped) | the field east of the road (R1-v140), loop to the north; 30 FLOAT/SINGLE where the link dropped |
| floodplain_RSplus_B | T2a (+ group) | Reach RS+ | 100 | the same field, a tighter loop |
| floodplain_CORS_Wickert | Wickert | Reach RS2+, MnDOT CORS | 17 | SW corner of the same field |
| benchmark_DR9925 | T5a, T5b, T5c | Reach RS2 | 2 | monument DR9925 "8508 Y", two 60-s occupations under B |

## Coordinate system and datums (what every column is in)

Horizontal: NAD83(2011) / UTM zone 15N, meters (EPSG:6344), for the RTK points and both
lidar grids (3DEP 2021 is NAD83(2011); the 2008 MN DNR tiles are NAD83(NSRS2007), the same
to centimeters; the pipeline labels the grids "EPSG:26915", the same projection). B points
are shifted onto the monument (see the tie), CORS points are as delivered. Vertical: NAVD 88
orthometric heights on GEOID18. The CORS rover delivers NAD83(2011) ellipsoid heights
directly; every B point is brought to the same ellipsoid through DR9925's published
NAD83(2011) height; both are then converted with GEOID18 through pyproj at each point's own
position (`navd88_tied_m`; the published 218.75 at DR9925 is reproduced to 2 mm). The 2021
lidar is NAVD 88 on GEOID18 as delivered (`gen2_z_m`, `gen2_zref_m`). The 2008 surface
(`gen1_z_m`) is the 2021 slope-normal grid minus the tile's difference, so it is in the
2021 datum BY CONSTRUCTION: the pipeline registered each 2008 swath onto the 2021 surface
(star frame, vertical shifts +0.064 to +0.081 m per swath, `corrections.json`) and applied
no geoid grid (`geoid_applied: false`); the GEOID03-to-GEOID18 difference went into those
shifts, and the 2008 column's absolute level is the registration's, not a measurement.

## Bases and the tie

The student exports show two base positions. **B1 was not the class base**: a student's
receiver was accidentally configured as a base and confused the radio link (Andy,
2026-11-11), so the 15 points corrected from it between 13:20 and 13:30 CDT, T3a 1-6,
T4a 1-8 and T1a 1, are dropped from the table; the compile log lists them. **B**
(the field east of the road, lon -91.998351, lat 44.117216, h 191.032, from 14:40 on) is
the class base; its position was the base's own average, so it carries an unknown offset
that the ties below remove. Ten points whose export lost the base columns were assigned to
B by time (all after 14:49; listed in the log).

**B is tied to monument DR9925 in both axes** (Andy, 2026-11-11: "Let's trust the reference
point and not on my comparison with CORS"). T5a's two 60-second occupations of DR9925
(unrounded Flow export, `raw/benchmark_DR9925_occupation_T5a_flow_export.csv`; Andy: "The
reference point was made on the metal disk, properly") read 189.202 and 189.210 m ellipsoid
against the published 188.669 m (NGS datasheet, `raw/DR9925_NGS_datasheet_2026-11-11.txt`,
network accuracy 0.61 cm horizontal, 1.49 cm height): B high by **+0.537 m**, removed from
every B point in `h_ell_tied_m`. Horizontally they sit at E 579636.167 / 579636.176,
N 4885873.225 / 4885873.203 against the published E 579636.942, N 4885872.102, so every B
point is shifted by **dE +0.771 m, dN -1.112 m** (`easting_tied`, `northing_tied`), and the
lidar is read at the tied position. After the tie the two benchmark rows land within 1 cm of
the monument and read 218.744 and 218.752 m NAVD 88 against the published 218.75.

**The check: the CORS-11 pair.** T1a point 63 (15:53:33) and Andy's CORS point 11
(15:53:52), the same spot one after the other, give an independent tie: CORS minus T1a
dE +0.802, dN -1.123 m, within 0.03 m of the benchmark in both axes; vertically 188.569 minus
188.096 = +0.473 m, 64 mm short of the benchmark. The pair was the applied tie until the
unrounded export arrived. Where the 64 mm comes from, read off the exports: T1a's Reach
RS+ was logged with antenna height 2.134 m, T2a's Reach RS+ with 2.065 m, and every
RS2 / RS2+ / RS3 with 2.134 m. Emlid's offset from the bottom of the receiver to the phase
center is 65 mm for the RS/RS+ and 134.9 mm for the RS2/RS2+/RS3, so the two RS+ entries
differ by 69 mm with T1a's reading lower, and a T1a 63 that is 69 mm low makes the
pair tie 69 mm short, against the 64 mm observed. Under the benchmark tie T2a's 98 points
sit -129 mm from the 2021 lidar and Andy's 17 CORS points -122 mm, 7 mm apart, while
T1a's 32 FIX points sit -228 mm, 99 mm below T2a on the same field. APPLIED (Andy,
2026-11-11: "after correcting the pole heights"): column `antenna_corr_m` raises each point by its
rover's implied-pole excess over the others' 2.000 m, +0.069 m on T1a's 62 B points and 0
everywhere else, and enters `h_ell_tied_m` and the `tied_` columns; `h_ell_m`, `navd88_m` and
`antenna_m` stay as exported. T1a's group moves to -159 mm, and the pair tie, recomputed
with the corrected T1a 63, becomes +0.542 m, 5 mm from the benchmark. The compile now prints this check on every run (Andy:
"checking antenna heights seems a necessary thing in all cases"): exported antenna height
minus the model's offset, one row per distinct value, flagged where it differs from the
others' 2.000 m pole; T1a's one B1-era point carried 2.086 m (+21 mm), the 62 B points
2.134 m (+69 mm); T4a's export (received 2026-11-12) carries 2.134 m on a Reach RS3, the same
entry as T3a's RS3, so no correction applies to either hillslope rover.

**Frame (fixed 2026-11-11).** The lidar grids are NAD83(2011) UTM 15N (3DEP; the pipeline's
"EPSG:26915" labels the same meters). The compile had projected lon/lat with EPSG:26915,
which made pyproj insert a NAD83(2011) to NAD83 datum shift of +0.012 m E, -0.242 m N here,
so every surveyed point sat 0.24 m south of the lidar. It now projects with EPSG:6344
(projection only). Floodplain medians moved by at most 3 mm, the RS3 bluff median by 34 mm,
single bluff points by up to 167 mm.

**B1** is dropped (above); no point in the table depends on it.

**Earlier horizontal test** (`analysis/field_survey/test_b2_horizontal_tie.py`, now in the
EPSG:6344 frame): the hillslope groups' NMAD against the 2021 lidar is 604 mm with no shift
and 136 mm with the pair shift, and a free search over +/-2.5 m on a 0.25 m grid puts the
minimum (116 mm) at dE +0.50, dN -1.25, within 0.3 m of the ties; the CORS-10 alternative
gives 193 mm.

## Result: RTK (tied) minus the lidar surfaces, FIX points, mm [median, NMAD, n]

(The figures show the negative, lidar minus RTK, so that the misfit reads as the lidar's.)

| group | base | vs 2021 lidar (gen2) | vs 2008 lidar (gen1) |
|---|---|---|---|
| benchmark_DR9925 | B | -96, 6, 2 | -54, 6, 2 |
| floodplain_CORS_Wickert | CORS | **-122, 33, 17** | -145, 57, 6 |
| floodplain_RSplus_B | B | **-129, 31, 98** | **-357, 81, 98** |
| floodplain_RSplus_A | B | -159, 50, 32 (pole-corrected, +0.069 m) | -232, 102, 22 |
| hillslope_T4a | B (shifted) | -2, 115, 21 | +50, 180, 5 |
| hillslope_RS3 | B (shifted) | -98, 162, 18 | +103, -, 1 |

Bare: in the sward field the surveyed ground lies 0.12-0.16 m below the May-2021 lidar
surface (CORS rover, which needs no tie, 0.12; T2a's RS+ 0.13; T1a's RS+ after the pole
correction 0.16), and 0.15-0.36 m below the November-2008 one.
On the bluff, with B shifted onto the monument, the two traverses sit -2 and -98 mm from the
2021 lidar with 115-162 mm scatter (the gen1 comparison there has 1-5 points: the 2008
difference is NaN under canopy on most of the traverse). At the monument itself, on the
mown roadside, the 2021 lidar reads 96 mm above the disk and the 2008 lidar 54 mm above it.

**Why the 2008 grid is blank on the bluff (run down 2026-11-11).** The blank is not the
no-2008-ground-return mask (`penetration_failure` is 0 at every blank point) and the 2008 raw
cloud is not empty there: in the 67 x 29-cell window around the traverse the blank cells hold a
median of 17 raw 2008 returns per 5 m cell (15 last returns), the same as the 16 (15) in the
cells that have a value. What is empty is the CSF ground selection: the cached CSF ground
(`data/csf_cache/4342-29-64.laz`, rigidness 1, PDAL defaults otherwise) holds a median of 0
points in the blank cells against 11-16 in the finite ones, and the blank fraction of the tile
climbs with the 2021 slope: 4 % at 0-5 deg, 15 % at 15-20, 30 % at 20-25, 47 % at 25-30, 64 %
at 30-35 and 79 % above 35 deg. The cloth does not drape the bluff. What it keeps is NOT the low tail:
against local 2021 planes at the 39 points (`analysis/field_survey/bluff/LOCAL_CLOUD_CHECK.md`),
the registered 2008 last returns within 5 m sit at a median of -119 mm (p25 -224, p75 +5) and the
CSF-selected 7.4 % of them at -101 mm, the central body with the canopy tail removed. The two
grid cells that read the 2008 surface 0.4 m below the ground at 29-30 deg are a 3-return CSF cell
and a disk whose 2008 returns are mostly canopy; the 2008 population at the bluff is ~0.12 m below
the 2021 surface, not 0.4.

**Is the bluff offset horizontal misregistration? No (`analysis/field_survey/bluff/NUTH_KAAB_RAW.md`).**
A 0.40 m vertical offset on this bluff (aspect 116 deg, 30 deg) needs 0.69 m of shift along the
downslope direction. The Nuth & Kaab (2011) fit of the registered 2008 raw last returns against the
slope-normal 2021 grid on all non-floodplain cells steeper than 15 deg leaves a residual amplitude
of 0.0194 m pooled (bootstrap 95 % 0.0082-0.0378) and 0.0695 m for swath 138, the only swath on
the traverse; projected onto the bluff those residuals predict +0.0004 +/- 0.0059 m (pooled) and
-0.0055 +/- 0.0067 m (swath 138) of gen1 minus gen2. The pipeline's applied per-swath shifts agree
with the independent raw-cloud fit to 2-6 cm. In the traverse window the registered raw last
returns sit -0.0611 m from the 2021 surface on cells under 10 deg and -0.1015 m on cells over 25
deg, and the pipeline's own difference grid on its 290 finite steep cells reads +0.2473 m against
+0.1239 m from a raw-last-return proxy on the same cells: the extra 0.12 m, and the blank cells,
are what the CSF selection calls ground there, not where the swath sits. Tile-wide, after
registration, the raw last returns show no reads-low-on-steep-slopes trend: slope-band medians
relative to 0-5 deg are within +/-0.011 m up to 35 deg and +0.029 m at 35-45 (agent's cuts at
15, 10 and 25 deg are proposals). One detail not
closed: tile-wide, 30,475 blank cells hold 1-9 CSF points and no finite cell holds fewer than
11, so a floor near ten acts somewhere after the CSF cache; `pipeline.difference` has no such
floor (`estimate_ground` blanks only empty cells).

## Not included

T3a's other CSV (9 October 2026, UTM 481 km E) is a practice session in St Paul.
Clast counts and OSL sample locations are in the submissions but are not elevation data.

Figure: `figures/rtk_2026-10-25_pdf_polecorr.png` (and `.pdf`, from
`analysis/field_survey/plot_rtk_pdf_2026-10-25.py _polecorr`; the unsuffixed file is the plot
before the pole correction, kept) shows these as probability densities,
floodplain and hillslope rows, minus 2021 and minus 2008 columns, one Gaussian KDE per group
over a rug of the points, nothing cut, nothing annotated.

Map: `figures/rtk_2026-10-25_diffmap_polecorr.png` (from `analysis/field_survey/plot_rtk_diffmap_2026-10-25.py _polecorr`;
also `_profiles_polecorr.png`, the elevation transects, from `plot_rtk_2026-10-25.py`; unsuffixed files
predate the pole correction; the one-panel map that script once drew is dropped). The figures carry no
titles, medians or provenance text (Andy, 2026-11-11): the tie and the corrections are stated here and in
the compile log.
puts a colored marker at every FIX point on 2025 FSA imagery, minus 2021 on the left and minus
2008 on the right; points with no 2008 difference (most of the steep bluff) are open gray.

Misfit transects: `figures/rtk_2026-10-25_profiles_misfit_polecorr.png` (from
`analysis/field_survey/plot_rtk_profiles_misfit_2026-10-25.py _polecorr`) is the elevation transect
figure with tied RTK minus each lidar surface on the y-axis, points only, nothing cut (one RS3
point at -1.6 m sets the hillslope panel's range).

**Why T3a's RS3 line reads farther from the lidar than T4a's (asked 2026-11-11).** Nothing
measurable separates the two traverses: slope at the points 31.2 vs 32.2 deg (medians), time
14:43-15:12 vs 14:44-15:34, rover vertical RMS 0.020 vs 0.018 m, 2021 class-2 points within 3 m
178 vs 203, local 2021 plane rms 0.052 vs 0.049 m. Pairing each RS3 point with its nearest T4a
point (19 m away on median) the RS3 point reads +100 mm more lidar-high, steadily so on the clean
east-facing part of the face. The RS3 export records pole tilt (median 2.4 deg) and the misfit does
not track it (correlation +0.05). Vegetation: the RS3 line sits in more 2021 low vegetation
(lowveg_gen2 0.215 vs 0.163), but within each group the misfit falls with vegetation (correlation
-0.41 and -0.35) and the group term with vegetation held fixed is +87 mm, so the 2021 proxy does
not explain it. Andy's notes say both groups used the same pole and the same model; both exports
carry antenna height 2.134 m (a 2.000 m pole entry) on a Reach RS3 (T4a's export arrived
2026-11-12; see "T4a's export" below).

**One or more translations (Andy, 2026-11-11), `analysis/field_survey/fit_hillslope_offsets_by_group.py`.**
For each group and for both pooled, the lateral shift of the TIED points that minimizes the NMAD of
lidar minus RTK against the 2021 grid, with the vertical offset read at that shift (grid +/-1.5 m,
0.05 m steps; bootstrap 500 draws; nothing cut). `figures/rtk_2026-10-25_hillslope_offsets_by_group.png`;
grids in `data/derived/4342-29-64/bluff_check/offsets_by_group.csv`.

| fit | n | shift that minimizes the spread (dE, dN) | NMAD before -> after | lidar minus RTK at that shift |
|---|---|---|---|---|
| RS3 (T3a) | 18 | -0.35, -0.10 m (bootstrap dE -0.45 to -0.05) | 162 -> 92 mm | +268 mm (p16-p84 +88 to +324) |
| T4a | 21 | -0.35, -0.10 m (bootstrap dE -0.55 to -0.25) | 115 -> 53 mm | +156 mm (+99 to +201) |
| both | 39 | -0.35, -0.15 m (bootstrap dE -0.35 to -0.05) | 131 -> 99 mm | +199 mm (+85 to +215) |

**The combined shift (Andy: "Is there a combined horizontal and vertical shift that best fits all
these points?").** Minimizing the mean absolute residual with the vertical closed-form (median), the
same script:

| points | n | dE, dN | lidar minus RTK at that shift | mean abs residual before -> after |
|---|---|---|---|---|
| hillslope, both groups | 39 | -0.20, -0.05 m | +122 mm | 163 -> 158 mm |
| every FIX point | 188 | -0.20, -0.05 m | +127 mm | 69 -> 64 mm |
| RS3 alone | 18 | -0.05, 0.00 m | +128 mm | 204 -> 204 mm |
| T4a alone | 21 | -0.30, -0.05 m | +139 mm | 114 -> 95 mm |

So one rigid shift, 0.20 m west and 0.05 m south of the tied positions with the lidar 0.12-0.13 m
above the surveyed ground, is the best fit to the whole survey and to the hillslope alone; the
floodplain's 149 flat points fix the vertical (their +0.125 m hardly moves with any shift) and
the 39 hillslope points fix the horizontal. The gain over no shift is small (mean absolute
residual 69 -> 64 mm overall) because the hillslope's scatter is wide; the minimum-NMAD criterion
above prefers a larger westward shift (0.35 m) that the mean-absolute criterion does not.

**All three axes free (Andy: "I am asking E, N, and z"), `analysis/field_survey/fit_shift_3d.py`.**
A full grid search over dE, dN (+/-1.5 m by 0.05) and dZ (-0.2 to +0.6 m by 0.01) added to the
tied RTK points, residual = lidar minus shifted RTK, nothing cut; `figures/rtk_2026-10-25_shift_3d.png`
shows RMS slices through the minimum.

| points | n | RMS minimum (dE, dN, dZ) | RMS before -> after | mean-abs minimum (dE, dN, dZ) | within 5 % of the RMS minimum |
|---|---|---|---|---|---|
| every FIX point | 188 | -0.10, -0.15, +0.13 m | 194 -> 150 mm | -0.20, -0.05, +0.13 m | dE -0.30 to +0.15, dZ +0.08 to +0.18 |
| hillslope, both | 39 | +0.05, -0.10, +0.05 m | 318 -> 309 mm | -0.20, -0.05, +0.12 m | dE -0.35 to +0.45, dZ -0.13 to +0.23 |
| RS3 alone | 18 | +0.25, -0.15, +0.06 m | 433 -> 395 mm | -0.05, -0.05, +0.12 m | dE -0.20 to +0.75 |
| T4a alone | 21 | -0.40, -0.15, +0.15 m | 163 -> 152 mm | -0.30, -0.05, +0.14 m | dE -0.60 to -0.20 |

With the vertical free the answer for the whole survey is dZ +0.13 m (the lidar 0.13 m above the
surveyed ground) with a small westward and southward lateral shift of one to two decimeters; the
vertical is tight (0.08 to 0.18 m at 5 %) and the lateral is loose. The hillslope on its own cannot
fix the shift by RMS: the -1.6 m RS3 point dominates the squares (RMS 318 mm against a mean absolute
residual of 166), and the 5 % region spans 0.8 m in dE and 0.36 m in dZ along the dE-dZ ridge; the
mean-absolute criterion, which that point cannot dominate, puts the hillslope at the same
(-0.20, -0.05, +0.12) as the whole survey.

**Weighted inversion, every logged point (Andy: "Instead of classifying FIX points, perform the
inversion taking all uncertainties into account"), `analysis/field_survey/invert_shift_weighted.py`.**
Residual = lidar minus (tied RTK + dZ) at the shifted position; sigma_i^2 = (reported vertical RMS)^2
+ (tan(slope) x reported horizontal RMS)^2 + (rms of the 2021 class-2 returns within 3 m about their
plane)^2; chi2 = sum (r/sigma)^2 over the full dE, dN, dZ grid; 1-sigma where the profile chi2 rises by
1 (the 3 m disk, the grid and the error model are mine). 203 of 220 points enter: CORS 1-3 and T1a
27-40 sit where the 2021 ground class is empty, 12-18 m to the nearest class-2 return, so the grid
there is filled from neighbors and they have no lidar datum (listed in the log, not a cut).
`figures/rtk_2026-10-25_shift_weighted.png` (chi2 slices through the calibrated minimum).

| solution | dE, dN, dZ | 1-sigma | reduced chi2 | note |
|---|---|---|---|---|
| sigmas as reported | 0.00, 0.00, +0.10 m | dE 0.00 to +0.05, dN -0.05 to 0.00, dZ +0.10 | 620 | T1a's 14 FLOAT points sit 4-24 m off the lidar with reported vertical RMS 0.05-0.24 m and supply 98.6 % of chi2: the reported FLOAT RMS is not an uncertainty |
| sigmas calibrated by status (reported x robust spread of r/sigma: FIX x2.3, FLOAT x79.5, SINGLE x0.8) | -0.05, -0.05, +0.13 m | dE -0.05 to 0.00, dN -0.05 to 0.00, dZ +0.13 | 1.63 | 2 of 203 points beyond 3 sigma; floodplain carries 0.95 of the weight |

The same inversion on the HILLSLOPE ALONE (`invert_shift_weighted.py hillslope`; Andy: "I meant the
hillslope lidar alone"): all 39 points enter (every one has class-2 returns within 3 m); reported
sigmas give dE -0.05, dN -0.05, dZ +0.14 m with reduced chi2 31.6; calibrated (FIX x3.3) the
minimum stays at -0.05, -0.05, +0.14 m with reduced chi2 2.93, 1-sigma dE -0.20 to +0.05, dN
-0.15 to 0.00, dZ +0.09 to +0.21 m, 1 of 39 points beyond 3 sigma
(`figures/rtk_2026-10-25_shift_weighted_hillslope.png`). The weighting is what the unweighted
hillslope fits lacked: each point's lidar roughness and reported RMS set its weight, the -1.6 m RS3
point and the rough north end count for little, and the hillslope then gives the same vertical as
the floodplain, +0.14 against +0.13 m, with the lateral residual within a grid step of zero.

The weighted answer is the 2021 lidar 0.13 m above the surveyed ground with no lateral shift
beyond the grid step, formally tight because 185 FIX points at ~0.05 m each enter independently;
the formal ranges do not include whatever is systematic between a pole tip and a May lidar ground
(the FIX factor of 2.3 is that excess over the reported RMS). The reported RMS of a FLOAT solution
understates its error by a factor near 80 on this survey.

**Why the same +0.13 m on a grassy floodplain and a forested bluff (Andy's question, 2026-11-11).**
First, the agreement is at the hillslope's 1-sigma, +0.09 to +0.21 m, so it is a match at +/-6 cm, not
to the centimeter. Second, what it is not: (1) not the geoid conversion, since pyproj's GEOID18 and
Flow's own NAVD 88 agree on Andy's CORS rover to 0.1 mm (max 0.8 mm, 17 points); (2) not the B tie,
since the CORS rover, which needs no tie, reads +122 mm on the same field as T2a's +129; (3) not a
common pole or antenna constant on the survey side, because such a constant c would enter the tied
B points as zero (it cancels through T5a's occupation) and the CORS rover as c, so T2a and CORS
would differ by c, and they differ by 7 mm. The +0.13 m is in the lidar relative to the pole tips, and
the flush monument disk in roadside grass reads +96 mm the same way. The 2021 low-vegetation proxy
cannot test the floodplain (it counts returns 0.15-2 m up and reads 0.00-0.01 over the sward) and is
not monotonic on the hillslope (quartile medians +31, +102, +110, -28 mm). What remains is that the
May 2021 ground surface sits about a decimeter above where a shod pole tip rests on both covers,
spring grass and forest floor, which is either one property of that acquisition's ground or two
causes of similar size; the project's bare control marks (-3.5 mm) say it is not there on bare
gravel. The measurement that separates the two is a bare surface occupied with the same rovers; this
survey has none (the gravel-lot points were under B1 and are dropped).

**Does the hillslope's lateral shift manufacture its vertical (Andy, 2026-11-11)?** The hillslope's
weighted chi2, profiled along dE with dN and dZ free: the ridge runs -400 mm of dZ per m of dE, and
within delta chi2 < 4 (2 sigma) it spans dE -0.40 to +0.20 m and dZ +0.28 to +0.04 m. With no lateral
shift allowed at all the hillslope's best dZ is +0.11 m (delta chi2 0.4 above the minimum); at the
floodplain's +0.13 the hillslope wants dE -0.05 (delta chi2 0.1); for dZ = 0 it needs dE +0.25 m
(downslope) at delta chi2 4.9, just outside 2 sigma. So the agreement is not made by the lateral
shift, which is 0.05 m at the minimum and 0.00 m gives +0.11, but the hillslope's vertical is weakly
held: on its own it says +0.14 with a 2-sigma range of +0.04 to +0.24 m, and the floodplain's +0.13
sits at its center rather than being tested by it.

**Is the fitted offset a coordinate-system fix or vegetation (Andy, 2026-11-11)? Vegetation, by the
acquisition's own checkpoints.** The 3DEP MN_SEDriftless 2021 acquisition's 519 surveyed checkpoints
(`analysis/GROUND_Q_FROM_CLASS2_SPREAD.md`, class-2 median minus surveyed truth): NVA, non-vegetated
open ground, n 227, **-3.5 mm**; VVA, sited under vegetation by design, n 162, **+103.3 mm**; LCP, the
acquisition's calibration points, n 130, -23.1 mm. On open ground the 2021 lidar is on its datum to
a few millimeters; under vegetation it reads about a decimeter high. This survey has no open-ground
point (the gravel-lot points were under B1 and are dropped), its floodplain is sward and its bluff is
forest floor, and its +0.13 m is the VVA number. The vertical part of the fitted offset is therefore
the lidar's vegetated-surface response, not a datum error, and the fitted table should be read as
"residual after removing the vegetated-surface offset", not as a corrected coordinate system. The
geoid, the tie and the rovers were excluded as its source above; the measurement that would show it
directly, the same rovers on a bare surface, is still missing from this survey.

**After the settled shift (fitted table), 2021 lidar minus RTK, FIX, mm [median, NMAD, n]:**
CORS -8, 34, 17; RS+ B -2, 31, 98; RS+ A +28, 50, 32; benchmark -35, 6, 2; T4a -108, 115, 21;
RS3 -7, 143, 18. Against 2008: CORS +15, RS+ A +101, RS+ B +224, benchmark -78, T4a -130 (4),
RS3 -216 (1). Figures from this table carry the suffix `_fitted` (`rtk_2026-10-25_{profiles_misfit,
pdf,diffmap}_fitted.png`); the `_polecorr` figures are the same plots before the shift.

**Does the translation close the gap between the two hillslope teams (Andy, 2026-11-11)? No.** T4a
-108 against RS3 -7 after the shift, a 101 mm difference, against -2 and +98 (100 mm) before. A
single rigid shift moves both traverses together; its 5 cm lateral part changes each by a few mm
through the re-read lidar and cannot change what separates them. The team difference stays a
vertical constant on one rover (above).

Which rover, and what kind of error (Andy's hypothesis, 2026-11-11: T4a "set up for a RS+ on a 2-m
pole with adapter but instead used a RS2"). Sign first: T4a's RTK heights sit about 0.10 m ABOVE
T3a's for the same ground, so whichever rover is wrong, T4a's entered antenna height would be
~0.10 m too SHORT, or T3a's ~0.10 m too LONG. A receiver-model mix-up in Flow is not available as a
setting: Emlid's documentation (read 2026-11-11) says the app applies the phase-center offset of the
connected receiver and the user enters only the pole height; if an RS+ offset (65 mm) had somehow been
applied to an RS2 (134.9 mm) the ground would come out 70 mm too high, T4a's sign but short of the
100 mm. A pole adapter of ~0.1 m under the receiver that was not added to the entered pole height fits
T4a's sign and size exactly; a 1.9 m pole entered as 2.0 on the RS3 fits T3a's. One circumstantial
pointer: the acquisition's vegetated checkpoints read +103 mm, and before the fitted offset T3a's
forest-floor traverse read +98 and T4a's +2, so T3a's is the one that looks like the rest of the
survey. The Flow export from T4a's unit (device, serial, antenna height) decides it.

**T4a's export (received 2026-11-12, `raw/rover_RS3_82434D3F6234D057_hillslope_T4a_T4b.csv`).**
It excludes every entry-side explanation above. T4a's team ran a second Reach RS3 (serial
82434D3F6234D057 against T3a's 8243D2AFB94BDAE3) with antenna height 2.134 m, the same 2.000 m
pole entry as T3a's, so neither rover's entered height differs from the other's and the
compile's antenna check now passes both. Both exports carry base B at the same coordinates
(E 580148.161 N 4885379.597 h 191.032) and the same coordinate-system label, NAD83(HARN) / UTM 15N
(the base's entered frame; T3a's adds "+ NAVD88 height" for its Elevation column). An RTK rover
reports base position plus baseline, so the two teams share whatever datum error the base carries,
and the benchmark tie removes it from both at once: a datum error cannot separate the teams. The
export also confirms the B1 assignment of T4a 1-8 from the base columns rather than from time,
and shows 26-sample averages on every point but the first (1 sample). Three points (24, 25, 28)
have no base columns, and 25 re-occupied 22 nine minutes later within 6 mm in height. Tilt is
logged: median 4.6 deg, maximum 14.1 over the 29 points, against T3a's 2.0 and 5.8 over 24; with the RS3's tilt
compensation the exported position is the pole tip, and without it 4.6 deg on a 2 m pole is 6 mm
of height, so tilt does not carry 100 mm either way. What remains, as hypotheses: a physical
difference at the pole tip (a pointed tip through ~0.1 m of fresh leaf litter on one traverse
against a topo shoe resting on it on the other, which gives T3a's sign: T3a's ground reads
0.1 m below the lidar surface and T4a's on it), a pole not fully extended on one unit, or a
surface difference between the two traverses that the 2021 lowveg proxy does not see. The first
two can be settled by asking the teams; the third by re-occupying both lines with one rover.

**Were the teams at different heights on the hill (Andy, 2026-11-12)? Yes, and it confounds the team
difference.** The lines are not parallel traverses of the same ground: T4a's 21 FIX points span
223-303 m NAVD88 with 14 above 275 m (upper slope and crest), 3 in 245-275 m and 4 at the toe below
245 m; T3a's 18 span 244-273 m with 17 in 245-275 m (mid-slope). 2021 lidar minus RTK by
position (median, NMAD, n, mm): T4a upper -10, 82, 14; T4a mid +12, 75, 3; T4a toe +101, 136,
4; RS3 mid +102, 157, 17; RS3 toe -28 (1 point). T4a's near-zero median is an upper-slope number
and T3a's +100 a mid-slope one; T4a's toe reads like the floodplain. In the shared band T4a
has three points and RS3 minus T4a is +87 mm, bootstrap p16-p84 -23 to +133, 95 % -91 to +165 mm,
a quarter of the draws at or below zero: three points cannot separate "team" from "position". Within
each team the residual does not trend with elevation (Theil-Sen -1.6 mm/m, 95 % -4.7 to +2.6, and
-0.5, -22.4 to +12.8), so a position effect would be a step between mid-slope and upper slope, not a
gradient; slope (31-33 deg) and lowveg_gen2 (0.20-0.21) are alike across the classes. This corrects
the paragraph above that called the lines parallel with nothing measurable between them: slope was
alike, elevation was not. The one-rover re-occupation should cover the mid-slope and the upper slope
with the same pole tip.

Both groups minimize at the same lateral shift, 0.35 m west (upslope on this east-facing bluff) of
the tied positions, and the shift halves T4a's spread; at that shift the two groups differ by
about 110 mm in level, as they do with no shift, so the RS3-T4a difference is vertical and not
positional. The minimum-clipped-rms criterion lands 0.1-0.2 m from the minimum-NMAD one (the dE-dZ
trade-off, 486 mm of level per m of easting, is still there: a 0.35 m shift moves the level by
+0.17 m, which is why the offsets at the minimum are larger than with no shift).

Translation fit (PARKED, Andy 2026-11-11: "more exploration than believing it to be correct"; an idea to
come back to): `analysis/field_survey/fit_hillslope_translation.py` fits a rigid (dE, dN, dZ)
of the UNTIED hillslope points onto the 2021 grid; `figures/rtk_2026-10-25_hillslope_translation.png`
(rerun 2026-11-11 in the EPSG:6344 frame). L1 minimum dE +0.55, dN -1.15, dZ -0.411 m (bootstrap
p16-p84: dE 0.55-0.85, dN -1.15 to -1.05, dZ -0.561 to -0.400; corr(dE, dZ) -0.96, -464 mm of
dZ per m of dE). At the benchmark tie's horizontal (+0.771, -1.112) the hillslope's own best dZ
is -0.503 against the tie's -0.537, i.e. the tied ground 34 mm below the 2021 surface; at the
pole pair's (+0.802, -1.123) it is -0.515 against the pair's -0.473, 42 mm above. The fit's
dE-dZ degeneracy is why the bluff cannot choose between the two ties on its own.

**Were the two RS3 rovers ever on the same ground (Andy, 2026-11-17)? Once, under B1, and the one pair is
confounded by tilt.** No point in the table is shared: the closest T4a-T3a pair on the bluff is 2.70 m
apart, and neither RS3 came within 3 m of any other rover's point. But between 13:24 and 13:28 CDT both teams
walked the toe of the hill side by side under B1 (T4a 4-8, T3a 1-5, logged within 2-90 s of each other),
and B1's position error is common to both rovers, so rover minus rover is valid there even though the points
are dropped from the table (`analysis/field_survey/b1_rover_pair_2026-10-25.py`, log `b1_rover_pair_log.txt`;
the 2021 lidar is read at the raw position and, as a sensitivity, at the position shifted by B's tie, since
B1's own error is unknown). T4a 4 and T3a 1: 2.93 m apart, 2 s apart, both FIX with 26 samples, on
1 deg ground; T4a's height minus T3a's is -46 mm, the lidar says T4a's spot is 21 mm higher, so
T4a's rover reads **67 mm LOW** against T3a's on the same ground (rover vertical RMS 0.026 and 0.030
m), the OPPOSITE sign to the bluff, where T4a reads ~100 mm high. A constant pole or antenna difference
would carry one sign everywhere, so taken at face value the pair argues against the pole and for a
ground-dependent effect. It cannot be taken at face value: T4a 4 carries the export's largest tilt,
14.1 deg, and if the RS3's tilt compensation was not active a 14.1 deg lean on 2.134 m puts the computed
ground 2.134 x (1 - cos 14.1 deg) = 64 mm below the tip, which is the 67 mm. Both RS3 exports carry a
populated Tilt angle column and the RS+/RS2/RS2+ exports carry none; Emlid staff on the community forum say
the column is written "for points collected with tilt" and that such coordinates "are already
tilt-compensated", which suggests the IMU was engaged but is not a statement about this unit's setting
(Emlid's docs pages on the export format were not reachable 2026-11-17). The next pairs are too far apart to
help: T3a 2 / T4a 5 (8.2 m, 8-11 deg, T4a 5 tilt 4.3 deg) gives T4a -56 mm raw and -134 mm with
the shift, T3a 3 / T4a 6 (7.3 m, 23 deg) +32 / +42, T3a 4 / T4a 8 (3.5 m but 1.3 m of relief at
20 deg) +224 / +185; the raw-to-shifted spread is the terrain correction's uncertainty. Lidar minus RTK at
these B1 points is +1.8 to +2.7 m, B1's own height error, the same for both rovers. What would settle the pair:
whether T4a's unit had tilt compensation on (ask T4b or T4a; the Flow project on that unit records it).

**Test: what if the RS3 exports were NOT tilt-compensated (Andy, 2026-11-17)?**
`analysis/field_survey/tilt_uncompensated_test_2026-10-25.py`, log `tilt_uncompensated_test_log.txt`. The
export is treated as (phase center) minus a vertical 2.134 m, so the tip would be 2.134 x (1 - cos tilt) higher
than the exported ground and 2.134 x sin(tilt) away horizontally in a lean direction the export does not
record. Vertical part, 2021 lidar minus RTK [median, NMAD, n, mm]: T4a +2, 116, 21 as compiled and -21, 100
under the hypothesis (tilt median 4.8 deg, correction median 7 mm, max 46); T3a +98, 162, 18 and +97, 164
(tilt median 2.4 deg, correction median 2 mm, max 11); the team gap goes from +97 to +117 mm, i.e. the
hypothesis widens it by 20 mm and cannot be the bluff difference. The B1 pair closes under it: T4a 4's 64 mm
and T3a 1's 1 mm turn -67 mm into **-4 mm** on the same ground. The bluff data cannot tell the two
hypotheses apart: regressing lidar minus RTK on the predicted correction gives slopes +0.66 (se 2.58) for T4a
and +1.60 (se 37.72) for T3a against the expected 1 (off) or 0 (on), because the corrections span only
1-46 mm under 100-160 mm of scatter. The horizontal part is the larger unknown: on the 30 deg bluff T4a's
leans of 0.07-0.44 m could move the lidar reading by up to 124 mm (median of the per-point bounds, max 278)
in an unknown direction if uncompensated, T3a's by 54 (max 129); the scatter does not track that bound
(Spearman rho of |residual from the team median| on the bound: T4a -0.43, p 0.05, n 21; T3a -0.37, p 0.13, n 18, the WRONG sign for uncompensated leans, which would make the scatter grow with the bound), and T3a's line, with the smaller leans, has the
larger NMAD. Reading: the hypothesis would resolve the B1 pair exactly and leave the bluff gap untouched, so
the sign flip between the toe and the bluff is NOT evidence against the pole unless T4a's unit had
compensation on; the bluff points themselves cannot say, and the unit's Flow project can. Pole handling, from the exports (2026-11-17): T3a's 24 points median 2.0 deg, p25-p75 1.4-2.6, max
5.8, one point over 5 deg; T4a's 29 points median 4.6, p25-p75 3.4-6.9, max 14.1, thirteen over 5 deg and five
over 10. T3a's team held the pole vertical; T4a's let it lean, as one does when the receiver is expected to
correct it, and (Andy, 2026-11-17) T4a's line was the high, steep part of the hillslope where a vertical pole is
hard to hold. The behavior reads as compensation on for T4a's unit; the setting itself is still to be asked for.

**T4a's answer (2026-11-17): "We used one of the blue 2 m poles, and we turned on tilt compensation."** So
the exported RS3 coordinates are the pole tip, the tilt-off hypothesis above is excluded, and the B1 pair stands:
on the flat toe, same base and same two seconds, T4a's rover read **67 mm LOW** against T3a's (rover
vertical RMS 0.026 and 0.030 m, plus whatever microtopography lies in the 2.93 m between the two tips, which the
5 m lidar grid cannot resolve; the grid itself puts T4a's spot 21 mm higher). On the bluff T4a reads ~100 mm
HIGH against T3a. Both teams carried the same 2 m pole and the same receiver model with the same 2.134 m
entry, and no constant, a pole length, an antenna offset or a base datum, can read low at the toe and high on
the bluff. The one way a constant survives is a pole re-extended between 13:28 and 14:44, when the base was
replaced; the teams have not been asked. What the gap can still be is in the ground: the pointed tip through
fresh litter on one line against a tip, or a shoe, resting on it on the other, or a true surface difference
between T3a's mid-slope and T4a's upper slope (the two lines share only three points in one elevation
band, so team and position are still confounded). The hypotheses standing after this, in the order of the
evidence: (1) a surface or litter difference between mid-slope and upper slope, the tip mechanism included;
(2) a pole change between the B1 and B sessions. The one-rover re-occupation of both lines, mid-slope and
upper slope with the same tip, decides (1); a question to both teams decides (2).

**How close in space and time was the B1 pair, and could B1 have been moving (Andy, 2026-11-17)?** T4a 4 and
T3a 1 are 2.93 m apart (dE -1.052, dN +2.733 m) on 1.1-1.4 deg ground; each is a 5.0 s average of 26
samples, the windows start 1.8 s apart and overlap for 3.2 s; both FIX, PDOP 1.5 and 1.3, the same 7 GPS,
6 GLONASS and 6 Galileo satellites, baselines 93.8 and 92.6 m. B1 was a student's rover set as a base, so it
may have been carried while broadcasting. Its broadcast position is the same to full export precision on every
one of the 14 B1 points in both files (lon -92.00038215, lat 44.11685493, h 187.878): the coordinates it sent
did not update, so a carried base puts the antenna's displacement from that fixed point into every rover
position as a common error that cancels in rover minus rover, up to what the base moved in the 1.8 s between the
two windows. The flat-ground points of the B1 session show no drift the pair could be blamed on: lidar minus
RTK on the five points under 2 deg of slope is 2.622, 2.611, 2.650 (the lot, 13:20:48-13:21:06), 2.672 (T4a 4,
13:24:47) and 2.605 m (T3a 1, 13:24:45), a 67 mm range over four minutes that includes the pair itself.
The same five points show the pair's real uncertainty, which the paragraph above understates: T4a 2 and T4a
3, 0.08 m apart and 8 s apart on the same lot gravel, differ by 39 mm in height with the same lidar under both;
the rovers' own RMS are 26 and 30 mm; and the 2.93 m between the two tips can hold microtopography the 5 m grid
does not see. So the pair separates "T4a reads 100 mm high" from "T4a reads 67 mm low" by 167 mm against
something like 70 mm of combined noise and terrain, roughly 2 sigma: a constant of +100 mm is disfavored, not
excluded, and a constant of a few tens of millimeters would be invisible here. The sentence above that no
constant can do both is corrected to that.

**Andy's reading (2026-11-17).** The lot points (T4a 1-3, 13:20-13:21) were the students' separate tests on flat
ground by the vehicles, nearby but not replicating one another; the teams then went onto the hill. "These
measurements are likely all to be trusted": the hillslope points stand as measured, and the ~100 mm between the
two lines is taken as a property of the ground they occupied (litter under the tip, or a true mid-slope versus
upper-slope difference), not of the instruments. The one-rover re-occupation of both lines, mid-slope and upper
slope with the same tip, is the measurement that resolves which.

**What the surfaces under the two lines look like from above (Andy, 2026-11-17: "look on satellite imagery to
characterize the surfaces").** Scripts `analysis/field_survey/surface_from_imagery_2026-10-25.py` (Sentinel-2
L2A, 10 m, Planetary Computer, the access path of mn-3dep-audit; NDVI at every hillslope FIX point for every clear
scene in three windows), `surface_from_imagery_extra_2026-10-25.py` (per-point test) and `naip_chips_2026-10-25.py`
(NAIP 0.6-1 m aerial chips with the tied points drawn, and NDVI within 3 m of each point); outputs and logs in
`imagery/`. Positions as above: upper > 275 m NAVD88 (T4a, 14 points), mid 245-275 (T3a, 17), toe < 245.

Sentinel-2 NDVI, median over the points, T4a upper vs T3a mid, with the per-point Mann-Whitney p: spring
leaf-off 2026-04-05 0.238 vs 0.233 (p 0.17), 04-10 0.218 vs 0.218 (0.51), 04-20 0.251 vs 0.241 (0.02); summer
2026-07-14 0.608 vs 0.622 (0.10), 07-29 0.567 vs 0.581 (0.13), 08-03 0.574 vs 0.585 (0.12), 08-08 0.609 vs 0.634
(0.08); around the survey, canopy mostly down, 2026-10-17 0.372 vs 0.303 (0.49), 10-22 0.298 vs 0.255 (0.57),
11-01 0.273 vs 0.223 (**0.004**), 11-08 0.250 vs 0.222 (0.06). In summer and in April the two lines are the same
to 0.01-0.02 (T3a's a shade greener in summer); after leaf fall T4a's upper slope keeps 0.03-0.07 more
NDVI than T3a's mid-slope. T4a's toe (4 points, the grass by the road) is the greenest class in autumn,
0.35-0.41.

NAIP NDVI within 3 m of each point (uncalibrated digital numbers, so relative only), T4a upper vs T3a mid:
2015-11-03 (leaf-off) 0.068 vs 0.095 (p 0.44); 2017-10-16 0.153 vs 0.112 (0.68); **2021-08-15, the lidar year,
0.209 vs 0.315 (p 0.002)**; 2025-09-29 0.270 vs 0.315 (0.80). At the scale of the lines, T3a's mid-slope sat
under a greener canopy in the August of the lidar year; the two autumn flights see no difference between the
lines, and both see the bare-canopy forest floor as gray-brown.

The chips (`imagery/naip_*_rgb.png`, `_cir.png`): both lines lie inside one closed deciduous stand on the
east-facing bluff, 100 m west of the road, 20-40 m apart; neither crosses an opening, a conifer patch or a
distinct stand. The conifers of the area are a band along the crest road to the north and scattered dark crowns
in the northwest, outside both lines. In the November 2015 image the floor between the bare crowns is gray-brown
litter under both lines, with light-gray patches west of and along T4a's upper line (bare ground, rock or
light trunks; not resolvable at 1 m) that do not appear along T3a's.

Measurement first: from above, the two lines are the same forest; the two differences that reach significance are
T3a's greener canopy at 3 m in August 2021 and T4a's greener floor at 10 m after leaf fall in 2026. As
hypotheses, both point the same way as the residuals: more canopy over T3a's line in the lidar year is the
side the 2021 lidar reads high on (its vegetated checkpoints are +103 mm), and a browner floor under T3a's
line in late October 2026 is a floor of fresh litter, the surface a pole tip goes into. Neither image can see
litter depth, and a 10 m pixel after leaf fall mixes shrubs, late-holding oaks and the floor; the re-occupation is
still the test.

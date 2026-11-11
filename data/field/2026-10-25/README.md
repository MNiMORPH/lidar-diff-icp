# RTK-GNSS ground survey, Whitewater WMA north of Elba, 25 October 2026

ESCI 4701/8701 Geomorphology field trip. Emlid Reach rovers; one on MnDOT CORS (Andy
Wickert), the rest on a class base station that was set twice. Compiled 2026-11-11 by
`analysis/field_survey/compile_rtk_2026-10-25.py`; every number here is pasted from its
log (`compile_log.txt`).

## The table: `rtk_points_all.csv` (220 rows)

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
`tied_minus_gen1_m`** (the numbers to use), `slope_deg`, `lowveg_gen2`, `clubb_floodplain`,
`penetration_failure_2008`, `source`.

Use FIX rows and the `tied_` columns. Groups:

| group | people | device | n | where |
|---|---|---|---|---|
| hillslope_RS3 | T3a, T3b, T3c | Reach RS3 | 18 (24 logged; 6 under B1 dropped) | bluff west of the road, up and north |
| hillslope_T4a | T4a, T4b | Reach | 21 (29 logged; 8 under B1 dropped) | same bluff, parallel traverse (table from the F07 .docx; raw export not uploaded) |
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
2.134 m (+69 mm); T4a's table exports no antenna height or model.

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
(`data/csf_cache/4342-29-64.las`, rigidness 1, PDAL defaults otherwise) holds a median of 0
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
over a rug of the points, nothing cut.

Map: `figures/rtk_2026-10-25_diffmap_polecorr.png` (from `analysis/field_survey/plot_rtk_diffmap_2026-10-25.py _polecorr`;
also `_map_polecorr.png` and `_profiles_polecorr.png` from `plot_rtk_2026-10-25.py`; unsuffixed files predate the pole correction)
puts a colored marker at every FIX point on 2025 FSA imagery, minus 2021 on the left and minus
2008 on the right; points with no 2008 difference (most of the steep bluff) are open gray.

Misfit transects: `figures/rtk_2026-10-25_profiles_misfit_polecorr.png` (from
`analysis/field_survey/plot_rtk_profiles_misfit_2026-10-25.py _polecorr`) is the elevation transect
figure with tied RTK minus each lidar surface on the y-axis, group medians dashed, nothing cut (one
RS3 point at -1.6 m sets the hillslope panel's range).

Translation fit: `analysis/field_survey/fit_hillslope_translation.py` fits a rigid (dE, dN, dZ)
of the UNTIED hillslope points onto the 2021 grid; `figures/rtk_2026-10-25_hillslope_translation.png`
(rerun 2026-11-11 in the EPSG:6344 frame). L1 minimum dE +0.55, dN -1.15, dZ -0.411 m (bootstrap
p16-p84: dE 0.55-0.85, dN -1.15 to -1.05, dZ -0.561 to -0.400; corr(dE, dZ) -0.96, -464 mm of
dZ per m of dE). At the benchmark tie's horizontal (+0.771, -1.112) the hillslope's own best dZ
is -0.503 against the tie's -0.537, i.e. the tied ground 34 mm below the 2021 surface; at the
pole pair's (+0.802, -1.123) it is -0.515 against the pair's -0.473, 42 mm above. The fit's
dE-dZ degeneracy is why the bluff cannot choose between the two ties on its own.

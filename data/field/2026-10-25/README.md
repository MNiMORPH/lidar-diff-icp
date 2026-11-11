# RTK-GNSS ground survey, Whitewater WMA north of Elba, 25 October 2026

ESCI 4701/8701 Geomorphology field trip. Emlid Reach rovers; one on MnDOT CORS (Andy
Wickert), the rest on a class base station that was set twice. Compiled 2026-11-11 by
`analysis/field_survey/compile_rtk_2026-10-25.py`; every number here is pasted from its
log (`compile_log.txt`).

## The table: `rtk_points_all.csv` (220 rows)

One row per logged point, nothing filtered. Columns: `group`, `people`, `device`, `serial`,
`point`, `t_start` (CDT), `base` (B1, B2 or CORS), `baseline_m`, `status` (FIX / FLOAT /
SINGLE as exported), `lon`, `lat` (NAD83(2011)), `easting`, `northing` (UTM 15N),
`h_ell_m` (ellipsoid height as exported), **`h_ell_tied_m`** (after the base tie below),
`tie_source`, `easting_tied`, `northing_tied` (B2 shifted, see Horizontal), `navd88_m`, **`navd88_tied_m`** (GEOID18 through pyproj, checked at DR9925:
218.748 vs published 218.75), `elev_app_m` (the Flow app's own NAVD88 where exported),
`antenna_m`, RMS columns, `tile`, `gen2_z_m` (2021 lidar ground, the horizontal per-cell median grid, bilinear from 5 m;
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
| benchmark_DR9925 | T5a, T5b, T5c | Reach RS2 | 2 | monument DR9925 "8508 Y", two 60-s occupations under B2 |

## Bases and the tie

The student exports show two base positions. **B1 was not the class base**: a student's
receiver was accidentally configured as a base and confused the radio link (Andy,
2026-11-11), so the 15 points corrected from it between 13:20 and 13:30 CDT, T3a 1-6,
T4a 1-8 and T1a 1, are dropped from the table; the compile log lists them. **B2**
(the field east of the road, lon -91.998351, lat 44.117216, h 191.032, from 14:40 on) is
the class base; its position was the base's own average, so it carries an unknown offset
that the ties below remove. Ten points whose export lost the base columns were assigned to
B2 by time (all after 14:49; listed in the log).

**B2 vertical: two ties, 64 mm apart; the CORS-11 pair is applied (Andy's call).** (1) T5a
T5a's two 60-second occupations of DR9925 read 189.202 and 189.210 m ellipsoid against the
published 188.669 m: B2 high by +0.537 m. (2) T1a point 63 and Andy's CORS point 11, the
same spot occupied one after the other: 188.569 minus 188.096, B2 high by **+0.473 m**. Andy
(2026-11-11): the two poles "at the exact same point on the ground, accurate to within a
centimeter", so the pair is a tie of benchmark quality in both axes, and +0.473 m is what
`h_ell_tied_m` removes from every B2 point; `tie_source` carries both numbers. What each tie implies: under the benchmark tie
T2a's 98 floodplain points sat -124 mm from the 2021 lidar and Andy's 17 CORS points
-121 mm, 3 mm apart; under the pair tie T2a sits -65 mm and the benchmark itself reads
64 mm above its published height. With the pair fixed to a centimetre, the 64 mm is a
question about the monument occupation, not the pair: DR9925 is a rod inside a PVC access
cover, and a pole tip resting on the cover rim or on the sleeve instead of the rod's disk
would read high by exactly that kind of amount; an antenna-height slip would do the same.
T5a's raw export gives the occupation's position unrounded and may show which. Switching
the table to the benchmark tie is one constant if that is ever warranted.

**B1 is tied PROVISIONALLY.** No benchmark was occupied under B1. The four FIX points on the
gravel of lot 204 read -2.636 m against the 2021 lidar, and that constant is removed from
every B1 point (`tie_source` says so). It rests on gen2 reading bare gravel correctly (the
control marks say -3.5 mm) and on those four points; a raw Flow export from any B1-era
rover that later re-occupied a B2 point would replace it.

**Horizontal (provisional).** The benchmark PDF prints the measured position rounded (0.1 m
E, 1 m N), so it cannot fix B2 horizontally. Andy's CORS point 11 (15:53:52) and T1a
point 63 (15:53:33) were the same spot occupied one after the other (Andy: "I am guessing
our tie point is 11"); the recorded difference, CORS minus T1a, is **dE +0.802 m,
dN -1.123 m**, and it is applied to every B2 point as `easting_tied`, `northing_tied`; the
lidar is read at the tied position. Test (`analysis/field_survey/test_b2_horizontal_tie.py`):
the hillslope groups' NMAD against the 2021 lidar falls from 540 mm with no shift to 139 mm
with this one, and a free search over +/-2.5 m on a 0.25 m grid puts the minimum (127 mm)
at dE +0.75, dN -1.00, within 0.1 m of the CORS-11 shift; the CORS-10 alternative gives
200 mm. T5a's raw Flow export of the benchmark occupation (requested) will replace this
with a direct tie to the monument.

## Result: RTK (tied) minus the lidar surfaces, FIX points, mm [median, NMAD, n]

| group | base | vs 2021 lidar (gen2) | vs 2008 lidar (gen1) |
|---|---|---|---|
| benchmark_DR9925 | B2 | -45, 6, 2 | -3, 6, 2 |
| floodplain_CORS_Wickert | CORS | **-121, 31, 17** | -146, 60, 6 |
| floodplain_RSplus_B | B2 | **-65, 32, 98** (-129 under the benchmark tie) | **-290, 92, 98** (-354) |
| floodplain_RSplus_A | B2 | **-163, 51, 32** (-227) | -233, 95, 22 (-297) |
| hillslope_T4a | B2 (shifted) | +79, 124, 21 | +125, 198, 5 |
| hillslope_RS3 | B2 (shifted) | +15, 158, 18 | +233, -, 1 |

Bare: in the sward field the surveyed ground lies 0.07-0.16 m (pair tie) or 0.12-0.23 m
(benchmark tie) below the May-2021 lidar surface, and 0.23-0.29 m (pair) or 0.30-0.35 m
(benchmark) below the November-2008 one; the CORS rover, which needs no tie, says 0.12 m
below 2021 and 0.14 m below 2008 where it has 2008 data. On the bluff, with B2 shifted
onto CORS 11, the two traverses sit +15 and -49 mm from the 2021 lidar with 124-158 mm
scatter (the gen1 comparison there has 1-5 points: the 2008 difference is NaN under
canopy on most of the traverse).

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

Figure: `figures/rtk_2026-10-25_pdf.png` (and `.pdf`, from
`analysis/field_survey/plot_rtk_pdf_2026-10-25.py`) shows these as probability densities,
floodplain and hillslope rows, minus 2021 and minus 2008 columns, one Gaussian KDE per group
over a rug of the points, nothing cut.

Map: `figures/rtk_2026-10-25_diffmap.png` (from `analysis/field_survey/plot_rtk_diffmap_2026-10-25.py`)
puts a colored marker at every FIX point on 2025 FSA imagery, minus 2021 on the left and minus
2008 on the right; points with no 2008 difference (most of the steep bluff) are open gray.

Translation fit: `analysis/field_survey/fit_hillslope_translation.py` fits a rigid (dE, dN, dZ)
of the UNTIED hillslope points onto the 2021 grid; `figures/rtk_2026-10-25_hillslope_translation.png`.
L1 minimum dE +0.65, dN -1.05, dZ -0.436 m (bootstrap p16-p84: dE 0.55-0.75, dN -1.15 to -0.95,
dZ -0.506 to -0.405; corr(dE, dZ) -0.93, -466 mm of dZ per m of dE); the pole-pair tie
(+0.802, -1.123, -0.473) sits at the edge of that band, and at its horizontal the hillslope's own
best dZ is -0.508, i.e. the ground 35 mm above the 2021 surface, the same +35 mm the cloud gives.

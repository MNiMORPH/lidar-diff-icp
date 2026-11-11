# RTK-GNSS ground survey, Whitewater WMA north of Elba, 25 October 2026

ESCI 4701/8701 Geomorphology field trip. Emlid Reach rovers; one on MnDOT CORS (Andy
Wickert), the rest on a class base station that was set twice. Compiled 2026-11-11 by
`analysis/field_survey/compile_rtk_2026-10-25.py`; every number here is pasted from its
log (`compile_log.txt`).

## The table: `rtk_points_all.csv` (235 rows)

One row per logged point, nothing filtered. Columns: `group`, `people`, `device`, `serial`,
`point`, `t_start` (CDT), `base` (B1, B2 or CORS), `baseline_m`, `status` (FIX / FLOAT /
SINGLE as exported), `lon`, `lat` (NAD83(2011)), `easting`, `northing` (UTM 15N),
`h_ell_m` (ellipsoid height as exported), **`h_ell_tied_m`** (after the base tie below),
`tie_source`, `easting_tied`, `northing_tied` (B2 shifted, see Horizontal), `navd88_m`, **`navd88_tied_m`** (GEOID18 through pyproj, checked at DR9925:
218.748 vs published 218.75), `elev_app_m` (the Flow app's own NAVD88 where exported),
`antenna_m`, RMS columns, `tile`, `gen2_z_m` (2021 lidar ground, bilinear from the 5 m
grid), `gen1_z_m` (2008 lidar ground = gen2 minus the tile's difference), `dod_2021_2008_m`,
`rtk_minus_gen2_m`, `rtk_minus_gen1_m` (untied), **`tied_minus_gen2_m`,
`tied_minus_gen1_m`** (the numbers to use), `slope_deg`, `lowveg_gen2`, `clubb_floodplain`,
`penetration_failure_2008`, `source`.

Use FIX rows and the `tied_` columns. Groups:

| group | people | device | n | where |
|---|---|---|---|---|
| hillslope_RS3 | T3a, T3b, T3c | Reach RS3 | 24 | bluff west of the road, from lot 204 up and north |
| hillslope_T4a | T4a, T4b | Reach | 29 | same bluff, parallel traverse (table from the F07 .docx; raw export not uploaded) |
| floodplain_RSplus_A | T1a, T1b | Reach RS+ | 63 | the field east of the road (R1-v140), loop to the north; 30 FLOAT/SINGLE where the link dropped |
| floodplain_RSplus_B | T2a (+ group) | Reach RS+ | 100 | the same field, a tighter loop |
| floodplain_CORS_Wickert | Wickert | Reach RS2+, MnDOT CORS | 17 | SW corner of the same field |
| benchmark_DR9925 | T5a, T5b, T5c | Reach RS2 | 2 | monument DR9925 "8508 Y", two 60-s occupations under B2 |

## Bases and the tie

B1 (lot 204, lon -92.000382, lat 44.116855, h 187.878) served ~13:20-13:30 CDT; B2 (the
field east of the road, lon -91.998351, lat 44.117216, h 191.032) from 14:40 on. Both
positions were the base's own average, so each carries an unknown absolute error and the
two differ by an unknown constant. Ten points whose export lost the base columns were
assigned by time (listed in the log).

**B2 is tied to the benchmark.** T5a's two occupations of DR9925 read 189.202
and 189.210 m ellipsoid against the published 188.669 m: B2 is high by **+0.537 m**, removed
from every B2 point. Independent check: five of Andy's CORS points lie within 1.4-2.0 m of
T1a points 62-63 on 1 deg ground; raw B2 minus CORS there is **+0.508 m** (NMAD 6 mm), 29 mm
from the benchmark value on ground that varies over that separation. The benchmark after
the tie sits 93 mm below the 2021 lidar cell it falls in (the monument is 11 m off the
highway shoulder).

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
| benchmark_DR9925 | B2 | -93, 6, 2 | -34, 6, 2 |
| floodplain_CORS_Wickert | CORS | **-121, 31, 17** | -144, 56, 6 |
| floodplain_RSplus_B | B2 | **-124, 34, 98** | **-347, 94, 98** |
| floodplain_RSplus_A | B2 | **-227, 48, 32** | -309, 125, 22 |
| hillslope_T4a | B1 (provisional) | +86, 165, 8 | +163, 209, 8 |
| hillslope_RS3 | B1 (provisional) | +285, 256, 6 | +343, 210, 5 |
| hillslope_T4a | B2 (shifted) | +15, 124, 21 | +69, 187, 5 |
| hillslope_RS3 | B2 (shifted) | -49, 158, 18 | +101, -, 1 |

Bare: in the sward field the surveyed ground lies 0.12 m (CORS, T2a) to 0.23 m
(T1a) below the May-2021 lidar surface and 0.30-0.35 m below the November-2008 one, on
two independent ties (benchmark and CORS) and three rovers. On the bluff, with B2 shifted
onto CORS 11, the two traverses sit +15 and -49 mm from the 2021 lidar with 124-158 mm
scatter (the gen1 comparison there has 1-5 points: the 2008 difference is NaN under
canopy on most of the traverse). The B1 rows (+86, +285) rest on the provisional B1 tie.

## Not included

T3a's other CSV (9 October 2026, UTM 481 km E) is a practice session in St Paul.
Clast counts and OSL sample locations are in the submissions but are not elevation data.

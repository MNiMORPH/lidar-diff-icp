# Additional absolute ground control for the Elba gen1/gen2 datum

**Why.** The gen1 datum is currently pinned by two surveyed 3DEP checkpoints –
2210 (+21.3 ± 12.4 mm, west chain) and 2036 (+28.9 ± 27.0 mm, east chain) –
combining to **+22.7 mm with σ_total 39.7 mm**. The dominant term is
**estimator + mark scatter, 40.8 mm RMS**, measured on gen2 against four NVA
checkpoints. Two of those four fail 3DEP's own 35 mm RMSEz on radius spread
alone. Separately, the DoD over stable ground shows a **tilt** of
dE −14.19 ± 5.15 mm/km and dN −16.70 ± 3.65 mm/km on 24 287 cells; four marks
cannot constrain a three-parameter plane.

This document extends `analysis/ridgelines/ABSOLUTE_ELEVATION_REFS.md`
(compiled 2026-08-22) and **corrects one of its conclusions**. Its airport search
and its assessment of the six in-AOI NGS marks stand and are not repeated. Its
finding that the 2008 gen1 checkpoints are not publicly downloadable **does not
stand** – see §5.

Compiled 2026-08-26. Every number is tagged with how it was checked.
Nothing committed to git; downloads went to the session scratchpad.

---

## Bottom line

1. **Two source documents had not been fully read, and both hold what was being
   searched for.**
   - The 2021 3DEP project's **Ground Control Survey Report PDF** holds **143
     LiDAR Control Points that appear in no shapefile**. The previous inventory
     read one shapefile (238 points) and concluded "zero checkpoints fall inside
     the AOI box." The project actually holds **534 surveyed points**, and one
     of them sits **inside the elbaext grid**.
   - The **MnGeo county validation reports** hold the **2008 gen1 checkpoint
     coordinates in full** – `Name, Control X, Control Y, Control Z, Surface Z,
     Error`, per point. `ABSOLUTE_ELEVATION_REFS.md` §1b called this a "dead end
     for public raw coordinates." **That is wrong**, and §5 shows the parse
     reproducing the published per-county RMSE exactly (Winona 0.160 vs 0.161 m).

2. **There is surveyed control 1.34 km from the reference point, on gen1's own
   geoid.** `L1O101`, an open-cover 2008 DNR control point at
   E 578 790.555, N 4 882 696.274, Control Z 223.291 m. Our local gen1 cloud
   gives a class-2 median of **223.290 m – a 1 mm difference at R = 5 m**, with a
   p05–p95 spread of only **0.110 m**. Same datum, same geoid (GEOID03), no
   chain, no conversion, **no download**. Three more 2008 control points and one
   2021 LCP also fall inside data already on disk.

3. **Counts.** **321 of 532** 2021 points fall inside gen1's eight-county
   footprint (137 NVA, 98 VVA, 86 LCP); 48 within 25 km of Elba, of which **16
   NVA/LCP lie under a flight line of the Elba network**, spanning −23.7 to
   +20.3 km in northing. **398 of ~1 033** 2008 gen1 control points are already
   parsed from three counties; 11 within 5 km, 36 within 10 km.

4. **The 40.8 mm floor is a siting artefact, and siting is now measurable.**
   Local p05–p95 spread at R = 5 m: **0.110 m** on open ground (L1O101),
   0.41–0.67 m on roads (the LCPs), metres at mark 2210 – where the surveyed
   elevation sits at **p95**. §7.1 turns this into a screening rule computable
   before any tie is attempted. **Note a prediction of mine that failed:** I
   argued from site photographs that LCP 1080's asphalt would be best-sited; the
   point cloud says it is on a raised embankment and sits at p82 (§1c). A
   photograph shows the surface, not the local relief.

5. **The tilt.** Two independent lines now bear on it. Sixteen on-line 2021
   marks would give σ(dE) 2.24, σ(dN) 0.84 mm/km against the ≤ 4.7 needed (§1.7).
   And the 2008 control *already* shows **no significant tilt within 10 km**
   (−1.32 ± 3.06 and −0.40 ± 2.67 mm/km, n = 36), with a real but smaller
   −4 to −5 mm/km east gradient appearing only at 20–30 km (§5). Together with
   the published 161 mm project RMSE, a coherent regional tilt of −14 mm/km is
   hard to sustain – **expect these marks to localise the tilt, not confirm it**.

6. **Binding constraint: siting, then spread. Count was never the problem** – and
   §5 raises a sign tension between the 2008 control residual and the +22.7 mm
   anchor that must be resolved before any of this is believed. See §7.5.

---

## 1. MN_SE_Driftless_2021_B21 – the full control set (Priority 1)

### 1.1 What exists, and where

**Datum, verified from the Ground Control Survey Report §1.8.4 (Woolpert,
January 2023, contract G16PC00022, task order 140G0221F0253) and repeated in
every coordinate-table header (§2.1–2.4):**

> "The spatial reference system for the project is NAD83 2011 (2010.00 epoch).
> Orthometric heights are based on NAVD88 vertical datum. **Geoid18** was used to
> determine the orthometric heights from the ellipsoid heights. The projected
> coordinates are displayed in Universal Transverse Mercator, Zone 15 North
> (UTM15N). Units for both the horizontal and vertical datums are expressed in
> meters to three (3) decimal places."

The shapefile `geoid` attribute independently reads `Geoid 18` on all 241
in-coverage rows. **Geoid model verified two ways.** (Note the shapefile's
`source_v_1` field still says "US Feet"; that mislabel was already caught and
disproved in `ABSOLUTE_ELEVATION_REFS.md` – the values are metres.)

**Point counts, verified from the report §1.3 and by parsing its tables:**

| class | report §1.3 | parsed from report tables | in a shapefile? | method |
|---|---|---|---|---|
| LiDAR Control Points (LCP) | 143 | 143 | **no – PDF only** | static GPS, ≥30 min, multi-MNCORS baselines (§1.8.1) |
| Non-vegetated check (NVA) | 227 | 226 | yes | RTK, 180 s, occupied twice (§1.8.2) |
| Vegetated check (VVA) | 164 | 163 | yes | RTK, 180 s, occupied twice (§1.8.2) |
| **total** | **534** | **532** | | |

(The two-point shortfall is my regex missing rows split across a page break; it
does not affect any candidate below. The shapefiles hold 157 QL0 + 238 QL1 =
395 rows, with a handful of points duplicated across both blocks.)

### 1.2 Download paths (all verified by fetching)

Root: `https://rockyweb.usgs.gov/vdelivery/Datasets/Staged/Elevation/metadata/MN_SE_Driftless_2021_B21/Vertical_Accuracy/`

- `contractor_provided/MN_Driftless_NVA_VVA_UTM15_QL1.{shp,shx,dbf,prj,gpkg}` – 238 pts *(previously used)*
- `contractor_provided/MN_Driftless_NVA_VVA_UTM15_QL0.{shp,shx,dbf,prj,gpkg}` – **157 pts, not previously used**
- `contractor_provided/Survey_Report/MN_SE_Driftless_2021_B21_Ground_Control_Survey_Report.pdf` – **1.67 MB; the only source of the 143 LCPs**, plus MNCORS station list and control diagram
- `contractor_provided/Survey_Report/MN_SE_Driftless_2021_B21_Ground_Control_Photos/{LCP,QL-0_Check_Points,QL-1_Check_Points}.zip` – **98.7 / 95.0 / 233.0 MB**; five photographs per point (close-up + N/E/S/W views). **This is the siting oracle.** A ZIP central directory can be range-read from the tail (~200 kB) and individual photos pulled by byte range – I did exactly that for three points, at ~1.5 MB total.
- `USGS/USGS_..._QL0_VA.txt` – QL0 block: NVA RMSEz **3.51 cm**, n = 91; VVA 95th pct **13.33 cm**, n = 66. (QL1 block figures were already recorded: NVA RMSEz 3.54 cm, n = 139.) Aggregate only; **still no per-point residuals**.

Recipe to regenerate the 532-point table:
`pdftotext -layout <report>.pdf` then match
`^\s*(\d{4}_20\d\d_MN)\s+(\d{7}\.\d+)\s+(\d{6}\.\d+)\s+(\d+\.\d+)\s+(LCP|NVA|VVA)\s*$`
– northing, easting, orthometric height, code, in that column order.

### 1.3 Which points gen1 actually overflew

gen1 coverage was taken as the union of MnGeo first-generation lidar tiles for
the eight counties the 2008 SE-MN project covers (Dodge, Fillmore, Houston,
Mower, Olmsted, Steele, Wabasha, Winona – the counties for which InPort 68818
publishes per-county RMSE). Tile lists fetched from
`resources.gisdata.mn.gov/pub/data/elevation/lidar/county/<county>/<county>_tile_list.txt`
(8 files), joined to the cached statewide centroid index
`data/mn_tile_centroids.csv`, giving **1 529 tiles**; containment tested against
the measured tile size (2 540 × 3 500 m).

**321 of 532 points fall inside gen1 coverage: 137 NVA, 98 VVA, 86 LCP.**
Within 15 km of Elba: 19 (9 NVA, 5 VVA, 5 LCP). Within 25 km: 48. Within 40 km: 132.

### 1.4 Flight-line reachability

Cross-track distance is computed against the nadir tracks fitted in
`analysis/ELBAEXT2_SCOPE.md` §2 (easting where each line crosses N 4 884 126),
with the measured swath half-width of ~730 m as the criterion. **Caveat, stated
plainly:** those tracks were fitted near Elba's latitude and the headings vary
179–183° / 357–359°, so a 1–3° heading error displaces a track by 350–1 050 m
over 20 km. **Line assignment is reliable within a few km of Elba and is an
extrapolation beyond that** – for the distant marks it must be confirmed by
probing the tile, exactly as `ELBAEXT2_SCOPE.md` §8 already flags for 2099.

`on-line? = YES` means |cross-track| < 730 m – the mark is under that line's
swath, needing **no chain link at all** if the line is one Elba already uses,
or a short chain otherwise.

All 48 surveyed points within 25 km of the Elba reference point (E 579 705.72,
N 4 883 677.71) that lie inside gen1 coverage. All NAVD88 / **GEOID18**, NAD83(2011)
epoch 2010.00, UTM 15N, metres. `dN` is northing offset from Elba (the tilt lever arm).

| id | code | UTM15 E | UTM15 N | lat | lon | NAVD88 (m) | dist km | dN km | nearest line (x-track) | gen1 tile | on-line? |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1079 | LCP | 579659.8 | 4885857.1 | 44.12157 | -92.00438 | 219.366 | 2.18 | +2.18 | **138** (-329 m) | `4342-29-64` | **YES** |
| 1080 | LCP | 577786.5 | 4881179.0 | 44.07966 | -92.02848 | 231.671 | 3.15 | -2.50 | **136** (-267 m) | `4342-30-64` | **YES** |
| 1078 | LCP | 579942.2 | 4888142.0 | 44.14211 | -92.00051 | 216.438 | 4.47 | +4.46 | **138** (-47 m) | `4342-28-64` | **YES** |
| 2036 | NVA | 585982.1 | 4884249.8 | 44.10638 | -91.92564 | 353.119 | 6.30 | +0.57 | **144** (+226 m) | `4358-29-03` | **YES** |
| 2210 | NVA | 570492.1 | 4884126.1 | 44.10693 | -92.11918 | 349.288 | 9.22 | +0.45 | **128** (+147 m) | `4342-29-61` | **YES** |
| 3056 | VVA | 570473.8 | 4884246.4 | 44.10801 | -92.11939 | 353.259 | 9.25 | +0.57 | **128** (+129 m) | `4342-29-61` | **YES** |
| 2024 | NVA | 571243.8 | 4887693.8 | 44.13897 | -92.10931 | 344.735 | 9.37 | +4.02 | **129** (-67 m) | `4342-28-61` | **YES** |
| 3088 | VVA | 590265.2 | 4882485.1 | 44.08998 | -91.87243 | 252.783 | 10.63 | -1.19 | 145 (+3501 m) | `4358-30-05` | no |
| 2073 | NVA | 571650.2 | 4876582.6 | 44.03890 | -92.10573 | 325.972 | 10.73 | -7.10 | **129** (+339 m) | `4342-31-61` | **YES** |
| 3089 | VVA | 582881.0 | 4894754.2 | 44.20130 | -91.96273 | 205.737 | 11.52 | +11.08 | 144 (-2875 m) | `4358-26-02` | no |
| 3050 | VVA | 570782.5 | 4876109.4 | 44.03473 | -92.11663 | 346.221 | 11.70 | -7.57 | **128** (+437 m) | `4342-31-61` | **YES** |
| 1081 | LCP | 574596.5 | 4872136.4 | 43.99859 | -92.06959 | 357.652 | 12.62 | -11.54 | **132** (+375 m) | `5142-01-62` | **YES** |
| 2099 | NVA | 571954.9 | 4893677.0 | 44.19276 | -92.09960 | 355.053 | 12.65 | +10.00 | **130** (-320 m) | `4342-26-61` | **YES** |
| 2035 | NVA | 588253.4 | 4873487.3 | 44.00923 | -91.89906 | 373.207 | 13.30 | -10.19 | 145 (+1489 m) | `4358-32-04` | no |
| 2160 | NVA | 569602.3 | 4874508.5 | 44.02043 | -92.13156 | 341.953 | 13.64 | -9.17 | 128 (-743 m) | `4342-32-60` | no |
| 3085 | VVA | 577464.7 | 4897146.0 | 44.22343 | -92.03015 | 321.617 | 13.65 | +13.47 | **135** (+346 m) | `4342-25-64` | **YES** |
| 2046 | NVA | 566833.7 | 4878515.0 | 44.05676 | -92.16560 | 343.499 | 13.87 | -5.16 | 128 (-3511 m) | `4342-31-59` | no |
| 2158 | NVA | 581228.6 | 4897478.9 | 44.22601 | -91.98298 | 234.981 | 13.89 | +13.80 | 138 (+1240 m) | `4358-25-01` | no |
| 1082 | LCP | 583677.3 | 4869939.6 | 43.97783 | -91.95670 | 361.280 | 14.30 | -13.74 | 144 (-2079 m) | `5158-01-02` | no |
| 2045 | NVA | 591601.0 | 4893330.1 | 44.18744 | -91.85386 | 210.338 | 15.32 | +9.65 | 145 (+4837 m) | `4358-26-05` | no |
| 3072 | VVA | 568622.7 | 4872893.8 | 44.00599 | -92.14399 | 349.882 | 15.46 | -10.78 | 128 (-1722 m) | `4342-32-60` | no |
| 2010 | NVA | 570565.3 | 4870504.0 | 43.98429 | -92.12008 | 344.775 | 16.03 | -13.17 | **128** (+220 m) | `5142-01-61` | **YES** |
| 3122 | VVA | 589895.7 | 4870855.1 | 43.98534 | -91.87902 | 368.187 | 16.38 | -12.82 | 145 (+3132 m) | `5158-01-04` | no |
| 2071 | NVA | 581687.5 | 4866623.5 | 43.94820 | -91.98201 | 363.150 | 17.17 | -17.05 | 138 (+1698 m) | `5158-02-01` | no |
| 2080 | NVA | 567913.3 | 4897041.7 | 44.22344 | -92.14973 | 344.449 | 17.82 | +13.36 | 128 (-2432 m) | `4342-25-60` | no |
| 3036 | VVA | 572706.0 | 4867278.1 | 43.95504 | -92.09383 | 380.723 | 17.83 | -16.40 | **130** (+431 m) | `5142-02-61` | **YES** |
| 2065 | NVA | 576968.6 | 4865962.2 | 43.94276 | -92.04091 | 408.076 | 17.93 | -17.72 | **135** (-150 m) | `5142-02-63` | **YES** |
| 3055 | VVA | 564169.9 | 4874578.4 | 44.02155 | -92.19933 | 363.306 | 18.00 | -9.10 | 128 (-6175 m) | `4342-32-58` | no |
| 2209 | NVA | 564185.3 | 4874395.9 | 44.01991 | -92.19916 | 360.815 | 18.08 | -9.28 | 128 (-6160 m) | `4342-32-58` | no |
| 2197 | NVA | 560732.0 | 4884572.2 | 44.11182 | -92.24107 | 359.272 | 18.99 | +0.89 | 128 (-9613 m) | `4342-29-57` | no |
| 2059 | NVA | 578714.5 | 4864537.0 | 43.92974 | -92.01936 | 379.543 | 19.17 | -19.14 | **137** (-331 m) | `5142-03-64` | **YES** |
| 3037 | VVA | 560056.7 | 4884501.9 | 44.11124 | -92.24952 | 331.358 | 19.67 | +0.82 | 128 (-10288 m) | `4342-29-57` | no |
| 3149 | VVA | 563266.9 | 4895290.1 | 44.20809 | -92.20811 | 335.818 | 20.13 | +11.61 | 128 (-7078 m) | `4342-26-58` | no |
| 2025 | NVA | 578387.7 | 4903992.4 | 44.28496 | -92.01757 | 316.371 | 20.36 | +20.31 | **136** (+334 m) | `4342-23-64` | **YES** |
| 2090 | NVA | 570320.8 | 4864879.7 | 43.93368 | -92.12388 | 401.174 | 21.01 | -18.80 | **128** (-24 m) | `5142-03-61` | **YES** |
| 2004 | NVA | 585615.8 | 4903994.3 | 44.28416 | -91.92699 | 203.049 | 21.16 | +20.32 | **144** (-140 m) | `4358-23-03` | **YES** |
| 2200 | NVA | 558225.5 | 4879269.2 | 44.06428 | -92.27298 | 348.776 | 21.93 | -4.41 | 128 (-12119 m) | `4342-30-56` | no |
| 3044 | VVA | 558519.3 | 4877713.9 | 44.05025 | -92.26948 | 356.646 | 22.01 | -5.96 | 128 (-11826 m) | `4342-31-56` | no |
| 2221 | NVA | 561791.6 | 4870885.3 | 43.98851 | -92.22943 | 379.392 | 22.01 | -12.79 | 128 (-8553 m) | `5142-01-57` | no |
| 2201 | NVA | 559674.4 | 4874502.1 | 44.02125 | -92.25542 | 400.273 | 22.03 | -9.18 | 128 (-10671 m) | `4342-32-56` | no |
| 2109 | NVA | 588868.7 | 4863311.7 | 43.91756 | -91.89309 | 358.352 | 22.33 | -20.37 | 145 (+2105 m) | `5158-03-04` | no |
| 3148 | VVA | 584661.8 | 4861804.9 | 43.90449 | -91.94572 | 347.800 | 22.43 | -21.87 | 144 (-1094 m) | `5158-04-02` | no |
| 2105 | NVA | 602120.9 | 4877279.6 | 44.04157 | -91.72537 | 372.329 | 23.31 | -6.40 | 145 (+15357 m) | `4358-31-09` | no |
| 1094 | LCP | 594089.0 | 4865197.3 | 43.93388 | -91.82775 | 370.205 | 23.42 | -18.48 | 145 (+7325 m) | `5158-03-06` | no |
| 2122 | NVA | 564290.5 | 4865339.3 | 43.93837 | -92.19894 | 396.582 | 23.96 | -18.34 | 128 (-6055 m) | `5142-02-58` | no |
| 3059 | VVA | 571570.6 | 4860912.1 | 43.89784 | -92.10884 | 346.999 | 24.18 | -22.77 | **129** (+260 m) | `5142-04-61` | **YES** |
| 1083 | LCP | 574738.3 | 4860012.0 | 43.88942 | -92.06953 | 325.632 | 24.18 | -23.67 | **133** (-447 m) | `5142-04-62` | **YES** |
| 3051 | VVA | 565897.9 | 4863246.9 | 43.91939 | -92.17918 | 401.237 | 24.66 | -20.43 | 128 (-4447 m) | `5142-03-59` | no |

### 1.5 The three near LCPs – siting from photographs (**superseded by §1c**)

I range-read the close-up and north-view photographs for 1078, 1079 and 1080 out
of `LCP.zip` (six images, ~1.5 MB of a 98.7 MB archive) and looked at them.

| id | UTM15 E, N | NAVD88 (m) | surface, from the photograph | judgement |
|---|---|---|---|---|
| **1080** | 577 786.5, 4 881 179.0 | 231.671 | **asphalt shoulder of a state highway**, tripod on the pavement just inside the fog line; grass ditch begins ~1 m west; wide flat carriageway east | *(my judgement at this stage: "best siting of the three" – paved, extended, dimensionally stable. **§1c shows this was wrong.**)* |
| **1079** | 579 659.8, 4 885 857.1 | 219.366 | open **gravel road**, ~5–6 m wide, mild crown, grass ditch and trees on the west side | good – flat and open, but a maintained gravel surface |
| **1078** | 579 942.2, 4 888 142.0 | 216.438 | open **gravel road**, similar, woods on the west side | good – same caveat |

**Caveat I will not soften: gravel roads get regraded.** 1078 and 1079 were
surveyed in September 2021; gen1 flew in November 2008. A regrade or resurfacing
in thirteen years moves the surface by centimetres to decimetres, and that motion
is real ground change, not datum error. 1080's asphalt does not have that
problem, and that argument still holds – but it is an argument about *stability*,
not about *siting*. **§1c measures the siting and overturns the ranking:** 1080
sits at p82 of gen1's local returns because its carriageway is a raised
embankment. The table and judgement above are left as written, with the
correction attached, rather than quietly revised.

### 1.6 Verification at 1079 against the local gen1 tile

Run on `data/las_local/4342-29-64.las` (already on disk – **no download**),
class 2, `laspy` chunked read:

| radius | n class-2 | flight lines present | median &#124;scan_angle_rank&#124; | gen1 z p05/p50/p95 (m) | surveyed sits at |
|---|---|---|---|---|---|
| 2 m | 8 | **138 only** | 6° | 219.14 / 219.26 / 219.42 | p88 |
| 3 m | 15 | **138 only** | 6° | 219.07 / 219.24 / 219.47 | p73 |
| 5 m | 42 | **138 only** | 6° | 218.82 / 219.27 / 219.49 | p62 |
| 10 m | 164 | **138 only** | 6° | 218.31 / 219.23 / 219.53 | p60 |

Three things this establishes. (1) The mark is genuinely overflown, by a single
unmixed line, at near-nadir geometry. (2) **Line 138 is already in the Elba swath
network** and already aligned in `data/derived/elbaext/corrections.json` – so
this tie needs **no chain at all**, unlike the five-link west chain and six-link
east chain. (3) The siting is materially better than mark 2210: surveyed
elevation at p62 (R = 5 m) rather than p95, and a p05–p95 spread of 0.67 m at
R = 5 m on a road corridor rather than a local high.

**No offset is quoted here, deliberately.** Turning this into a number requires
the project's slope-normal ground estimator at a justified radius, the
per-point GEOID03→GEOID18 conversion via `references.geoid_difference`, and the
swath correction – the same three requirements `ELBAEXT2_SCOPE.md` §4 sets out.
The raw percentiles above are sampling diagnostics, not a tie.

### 1.7 What the spatial spread buys – the tilt

Least-squares plane through the marks, 40.8 mm per mark, σ on the fitted
coefficients (this is geometry only; it assumes each mark contributes an
independent 40.8 mm and that a plane is the right model):

| mark set | n | σ(mean) mm | σ(dE) mm/km | σ(dN) mm/km |
|---|---|---|---|---|
| current: 2210 + 2036 | 2 | – | **cannot fit a plane** | – |
| + LCP 1078, 1079, 1080 | 5 | 20.6 | 3.73 | 8.11 |
| + the five on-line N/S marks at 18–21 km (2025, 2004, 2059, 2065, 2090) | 10 | 13.4 | 2.99 | **1.08** |
| all 16 on-line NVA + LCP within 25 km | 16 | **12.6** | **2.24** | **0.84** |

Target for a 3σ detection of the observed −14.19 mm/km east gradient:
σ(dE) ≤ 4.7 mm/km. **Every option except the current pair clears it.** The
near-LCP trio alone clears the east gradient (3.73) but is weak north–south
(8.11) because all three sit within ±4.5 km of Elba in northing; the ±20 km
on-line marks are what collapse σ(dN) to ~1 mm/km.

**A reasoning check on the tilt, not a measurement.** −16.70 mm/km sustained
over the 40 km northing span of the on-line set would be a 670 mm gen1–gen2
difference end to end. The 2008 project's own Winona-county vertical RMSE is
161 mm. So the tilt **cannot be a coherent regional plane**; it is local to the
tile or to the swath network's registration residual. That is precisely what the
distant on-line marks would settle – they bound the regional term hard, and
whatever survives is local. Treat this paragraph as an argument from published
accuracy, not as a result. **§5 turns it into a measurement:** gen1's own 2008
control shows no significant gradient within 10 km.

---

## 1b. gen2 at LCP 1079 – the siting argument, measured

gen2 is delivered on NAVD88 / GEOID18, the **same** datum and geoid as the
survey, so gen2-minus-surveyed needs no conversion and is a direct absolute
check. Run on the local `data/after/elbaext_3dep_fd_class2.laz` (**no download**):

| radius | gen2 n | gen2 density | gen2 median − surveyed | surveyed sits at | gen1 n | gen1 median − surveyed (**raw**) |
|---|---|---|---|---|---|---|
| 2 m | 94 | 7.5 pts/m² | **−11 mm** | p51 | 8 | −106 mm |
| 3 m | 220 | 7.8 | **−26 mm** | p54 | 15 | −126 mm |
| 5 m | 592 | 7.5 | **−26 mm** | p53 | 42 | −101 mm |
| 10 m | 2 376 | 7.6 | **−16 mm** | p51 | 164 | −136 mm |

**The gen1 column is raw** – no GEOID03→GEOID18, no swath correction, no
slope-normal estimator, plain median of vendor class 2. **It is not a tie and
must not be quoted as one.** It is here only to show the sampling is ample.

**The gen2 column is the finding.** Two things in it:

1. **Radius stability.** gen2's median moves by **15 mm across R = 2–10 m**. At
   mark 2210 the equivalent exercise moved by **hundreds of millimetres**
   (`ELBAEXT2_SCOPE.md` §4: −200 / −589 / −1 169 mm at R = 5/10/20 m, with a
   543 mm plane residual). The surveyed elevation sits at gen2's **median**
   here, versus 2210's **p95**. This is the "marks on local highs" diagnosis
   confirmed from the other side: on a flat road surface the radius sensitivity
   that generates the 40.8 mm floor largely disappears.
2. **Caveat on comparability.** A plain median within a radius is *not* the
   estimator that produced the 40.8 mm RMS figure, so −11/−26 mm is not a
   like-for-like replacement for that number. The **radius sensitivity** is the
   part that transfers, and it is what the floor is made of. Re-running the
   actual estimator at 1079 is the next step, and it costs nothing to run.

---

## 1c. Measured siting – the photograph argument does not survive the data

I predicted from the photographs that **1080's asphalt would be the best-sited
mark**. It is not. Both remaining LCP tiles turned out to be **already on disk**
(`data/before/4342-30-64.laz` and `4342-28-64.laz`), so I measured instead of
guessing:

| point | surface | R | n | lines | gen1 median − surveyed (**raw**) | p05–p95 | surveyed sits at |
|---|---|---|---|---|---|---|---|
| **1079** | gravel road, valley floor | 5 m | 42 | 138 only | −101 mm | 0.670 m | **p62** |
| 1080 | asphalt shoulder, TH 74 | 5 m | 39 | 136 only | −111 mm | 0.674 m | **p82** |
| 1078 | gravel road | 5 m | 43 | 138 only | −98 mm | 0.405 m | **p84** |

**The asphalt is on a raised road embankment.** The photograph shows the
carriageway sitting above its grass ditch – the surveyed point is on the
pavement, high relative to the ground within 5 m, which is the same p95 pathology
as mark 2210, just milder. **A photograph shows the surface; it does not show the
local relief.** Only the point cloud does. Correct §1.5 accordingly: 1080 remains
valuable for *stability* across thirteen years, but its siting is **worse** than
1079's, not better.

**All three LCP gen1 tiles are on disk.** `4342-29-64` (1079), `4342-30-64`
(1080), `4342-28-64` (1078). Tier 1 of the original plan costs **nothing**.

---

## 2. Adjacent 3DEP projects

The cached EPT boundary index (`data/ept_boundaries.json`) plus the USGS
metadata listing `.../Staged/Elevation/metadata/` give the complete set of
lidar projects touching southeast Minnesota. Verified by fetching the metadata
index and each candidate's directory.

| project | overlaps gen1 8-county footprint? | vertical accuracy dir? | verdict |
|---|---|---|---|
| `MN_SE_Driftless_2021_B21` | yes – it *is* gen2 | `Vertical_Accuracy/` | **§1, the whole story** |
| `MN_GoodhueCounty_2020_A20` | **yes, marginally** – its checkpoint line along the Goodhue county boundary falls into the northernmost gen1 tiles of Wabasha/Olmsted/Dodge/Steele | `vertical_accuracy/` | **secondary, see below** |
| `MN_UpperMissRiver_B22` | no – lat 45.75–47.92, central/northern MN | `vertical_accuracy/` (per-block `SurveyControl/`) | out of area |
| `MN_CentralMissRiver_B22`, `MN_RiverEast_B23`, `MN_RiverWest_B23` | no | – | out of area |
| `IA_EasternIA_2019_B19` | no – south of the MN state line, which is a hard edge of gen1 coverage | `vertical_accuracy/` | unusable for gen1 |
| `WI_12County_B22`, `WI_2County_B23`, `WI_8County_2020_A20`, `WI_Statewide_*` | no – east of the Mississippi; the MN-side extent in their EPT bounding boxes is convex-hull artefact, not coverage | – | unusable for gen1 |
| `MN_SEDriftless_1/2/3/4_2021` | these are **EPT delivery blocks of `MN_SE_Driftless_2021_B21`**, not separate projects | – | already counted |

**`MN_GoodhueCounty_2020_A20`, downloaded and checked:**
`.../metadata/MN_GoodhueCounty_2020_A20/vertical_accuracy/contractor_provided/`
holds `GoodhueCo_MN_NVA_VVA.shp` (**133 check points**) and
`GoodhueCo_MN__Ground_Control.shp` (**43 ground control points**), plus
`Ground Control Report ..._Ayres.pdf` and `Vertical Accuracy Report ..._Ayres.pdf`.

- **Coordinates are in `NAD_1983_HARN_Adj_MN_Goodhue_Feet`** – a county Lambert
  system on a modified spheroid, US survey feet (verified from the `.prj`).
  Converted to EPSG:26915 with pyproj.
- **Geoid model: GEOID12B**, verified from the project report
  (`USGS_MN_GoodhueCounty_2020_A20_ProjectReport.pdf`, "Geoid Model: GEOID12B").
  **Not GEOID18** – so these need a GEOID12B→GEOID18 conversion before they can
  be mixed with §1, computable per point from PROJ grids the same way
  `references.geoid_difference` already does GEOID03→GEOID18. Do not mix them
  without it.
- **38 check points and 8 ground control points fall inside gen1 coverage**,
  strung along the Goodhue county line from E 496 000 to E 560 000 at
  N ≈ 4 893 700, plus a few in northern Wabasha.
- **Distance from Elba: 20–90 km, all west or northwest.** No lever arm on the
  Elba tile that §1 does not already provide better, on a different geoid, with
  a coordinate conversion in the way. **Low priority – use only if a
  project-wide gen1 datum map is wanted.**

**Nothing else adjacent is usable.** No MN project other than these two overlaps
gen1's footprint, and the state line and the Mississippi cut off Iowa and
Wisconsin respectively.

---

## 3. NGS marks beyond the AOI box – widened, and filtered explicitly

Searched the whole eight-county footprint (lon −93.45…−91.15, lat 43.45…44.50)
via the NGS bounds API `geodesy.noaa.gov/api/nde/bounds` (12 tiles, plus 3 probe
queries; result cap disproved by subdividing the largest tile – 282 + 202 = 484,
exactly the parent count), then bulk datasheets through
`cgi-bin/ds_county.prl` (`MarkSelected` accepts many PIDs per POST, so 126 full
datasheets cost 8 requests, not 126). 3 s between every request.
*This section is the sub-agent's work; I have not independently re-run it.*

**Filter cascade, with what each step removed:**

| step | left | removed |
|---|---|---|
| bounds harvest | 4 362 | – |
| restrict to the eight counties | 3 012 | 1 350 (Goodhue, Rice, Freeborn, WI, IA…) |
| **(a) LEVELED height** (`ADJUSTED`) | 2 555 | 457: 251 GPS OBS, 100 VERTCON3, 85 blank, 16 SCALED, 5 other |
| drop "MARK NOT FOUND" / "POOR" | 2 450 | 105 |
| **(b) setting code could be flat** (31/35/36/40/66) | **126** | 2 324: 1 891 rods, 1 075 concrete-monument tops, 377 posts, 316 bridge abutments/piers, 160 footings, 79 retaining walls |
| drop narrow/vertical settings | 48 | 78 (49 bridge railing, 6 parapet, 2 curb, wingwalls, steps, pedestals, tank bases, doorsill, culvert) |
| read the description | **18** | 30: all 16 rock-outcrop marks are set *vertically* in highway-cut faces or are 0.6–3 m² ledges; 6 "BRIDGE" are headwalls or abutment corners; 3 "BUILDING" are set vertically in brick; 6 whose text contradicts the setting code |

A completeness check found 573 marks in the county lists but absent from the
bounds API – all condition X (destroyed/not found) except 16, of which 5 are
rods/posts and 11 are NGS non-published. **No usable mark was missed.**

All 18 survivors are NAVD 88 "determined by differential leveling", vertical
order SECOND CLASS I (OO0460 is FIRST CLASS II). GEOID18 appears on each
datasheet but applies only to the *geoid height* line, not to the leveled
orthometric height – **so these are geoid-model-independent**, which is what we
want.

### 3.1 Tier 1 – genuinely extended flat horizontal surfaces (6)

| PID | desig | UTM15 E | UTM15 N | NAVD88 (m) | levelled | county | recovery | from Elba | horiz. pos | surface |
|---|---|---|---|---|---|---|---|---|---|---|
| DK7698 | 2307 B 1 | 599 481 | 4 852 912 | 243.615 | 2009-05 | Fillmore | GOOD 2024 | 37 km, 148° | ±3 m | concrete radio-antenna pad, MnDOT truck shop, Rushford; a few m across |
| CN9146 | RST B | 541 022 | 4 861 847 | 392.411 | 2011-03 | Olmsted | GOOD 2010 | 44 km, 241° | **cm** | drain top at the centre of the **Rochester Intl airport apron**; apron vast, but the drain may sit below apron level |
| OO0025 | A 243 | 540 771 | 4 862 187 | 393.550 | 1991-06 | Olmsted | GOOD 2010 | 44 km, 242° | **±180 m** | disk flush with concrete, Rochester airport terminal |
| OO0028 | P 24 RESET | 541 348 | 4 855 557 | 378.560 | 2017-07 | Olmsted | GOOD 2013 | 48 km, 234° | **±180 m** | flush in a 4 × 7 ft (2.6 m²) slab, Stewartville elevator yard |
| ON1221 | 2309 L | 586 711 | 4 820 370 | 408.627 | 2024-02 | Fillmore | GOOD 2023 | 64 km, 174° | ±3 m | concrete **roof slab** of a boiler room, Canton schoolhouse |
| DM4186 | 5080 X 2 | 522 734 | 4 834 057 | 407.000 | 2011-01 | Mower | MONUMENTED 2007 | 76 km, 230° | **±180 m** | **bridge deck**, 210th St over Schwerin Creek – the best surface type in the set, never recovered since 2007 |

### 3.2 Tier 2 – bridge sidewalks (11), and the verdict

Eleven more are on bridge walkways: flat and horizontal, but 1.5–2 m strips with
railing along one edge, so returns are contaminated at both edges and each is
worth a handful of clean points with manual masking. They run 37–100 km from
Elba: DN5956, DN5960, DK7670, DN5963, DO9082, DP7526, DP7522, DP7512, OO0460,
DL4702, PP1657. All are hand-held-GPS positions (±3 m).

**Verdict – the widened NGS search does not help at Elba.**

- **Winona County, the county containing Elba, has exactly one usable mark:
  DG8385**, the flagpole-base slab already inventoried. The next nearest
  anything is **37 km** away.
- **Badly clustered:** 10 of the 18 are in Olmsted County, most in a ~12 km
  corridor along TH 63 through Rochester. The envelope reaches 128 km, so a
  regional tilt is constrainable *in principle*, but the design matrix would rest
  almost entirely on three or four singleton outliers.
- **Three of the six Tier-1 marks have SCALED horizontal positions (±180 m)** –
  they cannot be located in a point cloud from their coordinates at all, only
  from descriptive text, and DM4186 has not been recovered since 2007.
- **Surface extent is inferred, never measured** – NGS descriptions almost never
  give dimensions. Every "how large" judgement above is from the setting code and
  wording. Given §1c, that inference should be distrusted until checked against
  the point cloud.

---

## 4. MnDOT / MN DNR / county survey control

**MnDOT is the state's sole steward.** MnGeo's MSDI framework page names MnDOT as
the only data steward for geodetic control; **MnGeo and the DNR maintain no
separate control-point dataset**. Of the eight counties, none publishes a survey
or geodetic control layer on the Geospatial Commons.

**The Geospatial Commons package `loc_geodetic` is metadata only** – 12 kB of
XML, no points. The real data is on MnDOT's own server:

| route | URL | what |
|---|---|---|
| per-county CSV | `https://www.olmweb.dot.state.mn.us/geod/CSV%20Coords/{DODG1,FILL1,HOUS1,MOWE1,OLMS1,STEE1,WABA1,WINO1}.csv` | NAME, NAD83 lat/lon at five epochs, **NAVD88**, NGVD29, county zone, GSID |
| per-county datasheet PDF | `https://www.olmweb.dot.state.mn.us/geod/pdf/WINO1.pdf` | full sheets (Winona is 1 228 pp) |
| per-station datasheet | `.../Geod/PDF%20Metric%20Sheet/M_GSID_<gsid>.pdf` | one station |
| ArcGIS REST | `https://dotapp9.dot.state.mn.us/egis12/rest/services/OLM/mndot_geodetic_agol_operational_lyrs1/MapServer` | locations + order, **no heights** |
| index | `https://www.dot.state.mn.us/surveying/geodetics/geoindex.html` | |

**Verified by fetching `WINO1.csv`:** 1 017 rows for Winona County, **621
carrying a NAVD88 elevation**, in US survey feet. Across the eight counties the
agent's count is **5 867 stations, 3 411 with NAVD88** (that aggregate is from
the sub-agent and I re-verified only the Winona file).

**The critical limitation for our purpose.** The CSV carries **no geoid model and
no leveled-versus-GPS flag**. Those live only in the per-station datasheet, which
prints separate *"Leveling-Derived Orthometric Heights"* and *"Non
Leveling-Derived"* blocks with a `Determination Method` (VERTICAL CONTROL SURVEY
/ VERTICAL ADJUSTMENT / NON-RECIPROCAL VERTICAL ANGLE / GPS–STATIC / GPS–RTRN),
and geoid separations for **GEOID18, 12B, 09 and 03** at each station. So the
geoid model *is* recoverable per point – but at one PDF fetch per mark.

**And they are the same class of object as the NGS marks:** rods, disks, and
monument tops, set for surveyors, not flat surfaces a lidar can measure. Near
Elba the closest MnDOT stations are the same physical marks already inventoried
– **8508 G, GSID 41203, is NGS PID DG8387**, a disk on a 64-ft aluminium rod
(2nd order Class 1, ±0.016 ft, leveled, GOOD 2024) – already rejected as a rod
cap in `ABSOLUTE_ELEVATION_REFS.md`.

**Verdict: a large network, poorly suited to lidar.** Use it as a *lookup* to
check the geoid model and height source of a mark found some other way, not as a
source of new lidar-measurable targets. Contact: Geoff Bitner,
geoffrey.bitner@state.mn.us, 612-749-2113; olm.geodetic.support.dot@state.mn.us.

---

## 5. The 2008 gen1 checkpoints – **the premise was out of date; they are public**

`ABSOLUTE_ELEVATION_REFS.md` §1b concluded "not publicly downloadable as
coordinates … dead end for public raw coordinates." **That is wrong, and I have
verified it is wrong.** The per-county *validation reports* sit in the open
MnGeo lidar tree and contain the **full checkpoint tables**, not aggregate
statistics:

```
https://resources.gisdata.mn.gov/pub/data/elevation/lidar/county/<county>/<County>_county_validation_report.pdf
```

Capitalisation varies (`Winona_county_validation_report.pdf`,
`Wabasha_county_validation_report.pdf`,
`Olmsted_County_Validation_report.pdf`, lowercase for dodge/fillmore). Reports
exist for all nine counties: dodge, fillmore, freeborn, houston, mower, olmsted,
steele, wabasha, winona.

**Verified by downloading and parsing three of them.** Each row gives
`Name, Control X, Control Y, Control Z, Surface Z, Error`, UTM 15N metres,
Z in metres, with `Error = Control Z − Surface Z` (checked to 0.0000 m on every
row). Land cover is encoded in the point name: **L1O** open, **L2T** tall
weeds/crops, **L3B** brush, **L4F** forest, **L5U** urban.

| county | unique points | mean error | RMSE | published RMSE (InPort 68818) |
|---|---|---|---|---|
| Winona | **176** | −0.115 m | **0.160 m** | **0.161 m** ✓ |
| Wabasha | 97 | −0.027 m | 0.106 m | 0.106 m ✓ |
| Olmsted | 125 | −0.056 m | 0.117 m | 0.117 m ✓ |

The parsed RMSE reproduces the published per-county figure in all three cases –
that is the parse validating itself.

**Eleven control points lie within 5 km of Elba; four are inside the elbaext
grid** (E 575 450–580 200, N 4 882 050–4 886 400), i.e. inside gen1 data
**already on disk** as `data/before/elbaext_gen1_merged.laz`:

| name | cover | UTM15 E | UTM15 N | Control Z (m) | DNR error (mm) | km from ref pt |
|---|---|---|---|---|---|---|
| L1O101 | open | 578 790.555 | 4 882 696.274 | 223.291 | −14 | 1.34 |
| L3B99 | brush | 580 050.679 | 4 885 346.715 | 218.417 | +15 | 1.70 |
| L5U171 | urban | 578 582.791 | 4 882 169.886 | 225.482 | +143 | 1.88 |
| L2T51 | crops | 576 507.634 | 4 883 540.432 | 344.129 | +45 | 3.20 |

**Datum – and this is the one thing that is *not* verified per point.** The
validation reports state **no datum and no geoid**. The linkage comes only from
the dataset-level metadata,
`resources.gisdata.mn.gov/pub/data/elevation/lidar/documentation/lidar_semn2008.html`:
*"Vertical datum: NAVD88 (Geoid03)"*, and process step 5, *"Geoid Model used to
reduce satellite derived elevations to orthometric heights – NGS Geoid03."*
That is a dataset-level assertion, not a per-mark one. The Wabasha point names
carry the string `RTK` and the Olmsted names `VRS`, so these are **GPS-derived,
not leveled** – they inherit GEOID03 fully. **Any use against our GEOID18 gen1
surface must apply `references.geoid_difference` per point.**

**What they measure, and a sign warning.** `Surface Z` is the **delivered 2008
DNR DEM**, not our reconstruction from the point cloud, so the DNR `Error`
column is not our residual. `Control Z` is the useful column: 398 surveyed
elevations we can evaluate our own gen1 surface against.

Statistics of the DNR residual (`Control − Surface`; **negative means the 2008
surface reads high**), three counties pooled:

| band from Elba | n | mean (mm) | median (mm) | RMSE (mm) |
|---|---|---|---|---|
| 0–5 km | 11 | −20.4 | −6.0 | 87.8 |
| 5–10 km | 25 | −22.6 | −35.0 | 72.5 |
| 10–20 km | 64 | −74.4 | −64.5 | 133.7 |
| 20–40 km | 195 | −87.2 | −76.0 | 145.9 |

Plane fits to the same residual (mm and mm/km, from the Elba reference point):

| radius | n | intercept | dE | dN | resid RMS |
|---|---|---|---|---|---|
| 10 km | 36 | **−22.4 ± 13.0 mm** | −1.32 ± 3.06 | −0.40 ± 2.67 | 74 mm |
| 20 km | 100 | −45.3 ± 9.5 mm | **−5.35 ± 1.11** | +1.80 ± 0.99 | 91 mm |
| 30 km | 208 | −69.4 ± 6.6 mm | −3.70 ± 0.48 | +1.69 ± 0.45 | 92 mm |

**Read this carefully, and do not over-read it.** (1) Within 10 km there is **no
significant tilt** in gen1's own control residual – both gradients are
consistent with zero, which is independent support for the §1.7 argument that a
−14 to −17 mm/km tilt cannot be regional. (2) At 20–30 km a real east gradient
of −4 to −5 mm/km appears, about a third of the DoD's −14.19 mm/km. (3) The
intercept sign says the **2008 DEM sits ~22 mm above the surveyed control near
Elba, on GEOID03**, whereas the 3DEP anchors give **+22.7 mm** – **I am not
claiming these disagree.** They are different surfaces (vendor DEM vs our
slope-normal reconstruction), different geoids (GEOID03 vs GEOID18, a +67 mm
difference at Elba), and I have not checked our anchor's sign convention against
this one. **Reconciling them is a task, not a result.** It is, however, the most
interesting thing this search turned up.

**Still not public, and worth requesting:**
1. The **"Lidar Accuracy Assessment Report"** – listed in the metadata's Final
   Deliverables as *"One paper copy"*, so plausibly never digitised.
2. The **127 AeroMetric QA/QC checkpoint coordinates**.
3. **The value of the bias adjustment.** Metadata process step 8: AeroMetric's
   127 points were run through TerraScan's Output Control Report, *"a bias
   adjustment was determined, and the results were applied to the lidar data"*,
   final OCR RMSE 0.109 m. **The magnitude of that applied bias is not
   published.** For a project about gen1's floating datum, this is the single
   most valuable unpublished number in the whole search.

Minor discrepancy worth raising in the same email: the metadata's per-county
counts sum to **1 033**, while its text says **1 009**.

**Precedent to cite when asking:** the DNR *has* published checkpoints as a
shapefile before – `lidar_checkpts_pine2007`, 100 survey-grade GPS points,
documented at
`resources.gisdata.mn.gov/pub/data/elevation/lidar/documentation/lidar_checkpts_pine2007.html`.

**Who to ask.** MN DNR lidar data steward **Sean Vaughn**,
sean.vaughn@state.mn.us, 763-284-7223. MnGeo: gisinfo.mngeo@state.mn.us,
651-201-2499; director **Alison Slaats**, alison.slaats@state.mn.us.
*(Contacts are from the sub-agent's fetch of the 2025 statewide lidar DEM
metadata and MnGeo contact pages; I did not independently re-verify them.)*

**Data Practices Act?** Not as the first move – this is public government data
under Minn. Stat. ch. 13 and a plain email normally suffices; a paper-only 2008
deliverable is a retrieval problem more than a disclosure one. Escalate only if
the informal request stalls.

### 5.1 Measured: our own gen1 cloud against the 2008 control, at Elba

`data/before/elbaext_gen1_merged.laz` is a **pure crop-and-concatenate** of the
downloaded tiles (`scripts/merge_gen1_tiles.py` copies dimensions and touches no
elevation), so it carries gen1 exactly as delivered: NAVD88 on **GEOID03** – the
same datum and geoid as the control. **No conversion is needed for this
comparison, which is why it can be made today.** Plain median of vendor class 2:

| point | cover | R | n | gen1 median | Control Z | median − control | p05–p95 | control at |
|---|---|---|---|---|---|---|---|---|
| **L1O101** | **open** | 3 m | 14 | 223.285 | 223.291 | **−6 mm** | 0.091 m | p71 |
| **L1O101** | **open** | 5 m | 38 | 223.290 | 223.291 | **−1 mm** | **0.110 m** | p63 |
| L5U171 | urban | 5 m | 31 | 225.450 | 225.482 | −32 mm | 0.345 m | p77 |
| L2T51 | crops | 5 m | 41 | 344.080 | 344.129 | −49 mm | 0.370 m | p61 |
| L3B99 | brush | – | – | – | 218.417 | – | – | **no data**: 7 m east of the merged file's east edge (max E 580 043.5) |

**L1O101 is the best-conditioned absolute target found anywhere in this search:**
open cover, 1.34 km from the reference point, a p05–p95 spread of **110 mm** at
R = 5 m against 670 mm at LCP 1079 and metres at mark 2210 – and gen1 lands on
it to **1 mm**.

**Now interrogate that, because a 1 mm agreement is exactly the kind of pleasing
number that is usually an artefact:**

- **n = 38 at R = 5 m.** One point, no averaging. The single-point scatter is the
  110 mm spread, not 1 mm; 1 mm is where the median happened to fall.
- **A plain median is not the project's slope-normal estimator.** Same caveat as
  §1b.
- **The control itself is RTK/VRS**, not leveled – its own uncertainty is
  plausibly a few centimetres, and it inherits GEOID03 in full.
- **Independence.** The DNR's 1 009-point set is independent of the AeroMetric
  127 points that produced the applied bias adjustment, but it is not
  independent of the delivered data's *acceptance*: this control is what
  certified gen1. Agreement is partly expected by construction.
- **Cross-check that does hold:** the DNR's own reported error for L1O101 is
  −14 mm against their DEM; we get +1 mm against the point cloud. Agreement to
  15 mm between two different surfaces built from the same data is a sanity
  check on both our read and their table.

**What this does not do:** it does not settle the +22.7 mm anchor, because that
lives on GEOID18 after swath correction and this does not. It **does** say that
gen1's raw, uncorrected surface near Elba is already within a few centimetres of
its own 2008 control on its own geoid – which is a strong constraint to carry
into the reconciliation.

### 5.2 Draft request email (send as-is)

> **Subject:** Request – 2008 SE Minnesota lidar accuracy assessment report and checkpoint coordinates
>
> Dear Sean and MnGeo staff,
>
> I am using the 2008 Southeast Minnesota lidar (MN DNR / AeroMetric) together with the 2021 3DEP coverage to measure landscape change near Elba, in Winona County. Resolving a systematic vertical offset between the two epochs requires the original ground control.
>
> The county validation reports on the MnGeo resources site have been very useful – they give the full DNR checkpoint tables, and I have those. Three things are not public, and I would like to request them:
>
> 1. The Lidar Accuracy Assessment Report, listed in the project metadata under Final Deliverables as "one paper copy." A scan is fine.
> 2. The 127 AeroMetric QA/QC checkpoints, as a shapefile, CSV, or table of coordinates and elevations.
> 3. The value of the bias adjustment described in process step 8 of the metadata, which was determined from those points and applied to the delivered data.
>
> If the checkpoints exist as a shapefile, the 2007 Pine County lidar checkpoint dataset you already publish is exactly the format I need.
>
> One small thing you may want to know: the per-county checkpoint counts in the metadata sum to 1,033, while the text says 1,009.
>
> I am happy to cover copying or scanning costs, and to tell you what I find.
>
> Thank you,
>
> Andy Wickert
> University of Minnesota

---

## 6. Other absolute vertical references

*(Airports were checked in `ABSOLUTE_ELEVATION_REFS.md` §3 – none within 30 km –
and are not revisited.)*

**CORS / continuous GNSS – MnCORS.** Five stations bear on the area
(from `https://www.dot.state.mn.us/surveying/cors/mncors_site_info.pdf` and the
MnDOT metric datasheets, fetched by the sub-agent; **I did not re-verify these
numbers**):

| site | NGS ID | lat, lon | NAVD88 (m) | mount |
|---|---|---|---|---|
| Wabasha (~31 km N) | MNWB | 44.379541, −92.030153 | 228.89 | **roof of Wabasha County Courthouse** |
| Winona | MNWN | 44.065002, −91.712633 | 210.12 | attached to the **side** of a building |
| Eyota | MNEY | 43.955624, −92.205389 | 401.79 | RWIS mast |
| Stewartville | MNSV | 43.902235, −92.482481 | 392.16 | RWIS mast |
| Rushford | MNRF | 43.822613, −91.762760 | 252.56 | truck-station building |

MNWB's datasheet (GSID 62781, NGS PID DS3253) gives NAVD88 228.896 m from a
**non-reciprocal vertical angle, 3rd order, 2006**, and prints geoid separations
for GEOID18/12B/09/03. **All CORS heights are the antenna reference point**,
metres above the roof; the lidar measures the roof, and **the ARP-to-roof offset
is not published**. Third-order vertical on a rooftop, with an unpublished
offset, is worse than what we already have. **Not usable.**

**Locks and dams – not usable without an independent tie.** The USACE water
control manual for **Lock and Dam 5** (44°09′42″N 91°48′42″W, ~17 km ENE of
Elba) publishes structure elevations – top of lock walls 665.0 ft, top of bridge
deck 688.0 ft, roller gate sill 640.0 ft, pool 660.0 ft – but the datum is
**"MSL – 1912 adjustment"** throughout, and the strings "NAVD" and "NGVD" appear
**zero times** in the manual. These are also design elevations, not as-built
survey. The 1912 river datum is district-wide, so L&D 4, 5A and 6 carry the same
problem (*unverified individually*). **Unknown offset to NAVD88 – reject.**

**Lidar calibration range in Minnesota: none found.** MnROAD (Albertville,
~200 km NW) is a pavement research facility and there is no statement that it
serves as a lidar vertical calibration range. Reported as a negative result.

---

## 7. Ranking, and what to acquire

### 7.1 Measured siting quality – the screening criterion

Everything above converges on one point: **the marks differ far more in *siting*
than in *survey quality*, and siting is measurable from the point cloud.** For
every mark reachable in local gen1 data I measured, at R = 5 m on vendor class 2,
the p05–p95 spread and where the surveyed elevation falls in the local return
distribution:

| mark | source | cover / surface | p05–p95 at R = 5 m | surveyed at |
|---|---|---|---|---|
| **L1O101** | 2008 DNR control | **open ground** | **0.110 m** | **p63** |
| L1O59 | 2008 DNR control | open ground | 0.139 m | p19 |
| L1O144 | 2008 DNR control | open ground | 0.339 m | p36 |
| 1079 | 2021 LCP | gravel road, valley floor | 0.670 m | p62 |
| 1078 | 2021 LCP | gravel road | 0.405 m | p84 |
| 1080 | 2021 LCP | asphalt shoulder (raised) | 0.674 m | p82 |
| L2T/L3B/L4F points | 2008 DNR control | crops / brush / forest | 0.29–1.83 m | p4–p49 |
| 2210 | 2021 NVA | road shoulder | metres (plane resid 543 mm) | **p95** |

**A screening rule falls out of this, and it is the most useful thing in this
document:** accept a mark only if (i) its local p05–p95 spread at R = 5 m is
small, and (ii) the surveyed elevation falls near the middle of the local
distribution rather than in a tail. Both are computable *before* any tie is
attempted, from data we have. The threshold itself should be chosen from the
distribution over the full candidate set, not guessed – I have 17 measurements
here, which is enough to show the criterion works and **not** enough to set the
cut. Set it properly once all candidates are measured.

**Measured against 17 local 2008 control points**, gen1 median − control is
mean **+57 mm, sd 107 mm, SE 26 mm**; on the four open (L1O) points alone,
mean −24 mm, SE 49 mm. **The per-point scatter of the 2008 control (107 mm) is
worse than the 3DEP marks' 40.8 mm** – it is dominated by the vegetated classes
(L2T crops median +154 mm, exactly the November-stubble contamination already
documented in the bare-ground-skew work). Open-cover points are the usable
subset, and there are far fewer of them.

### 7.2 What each source actually buys

| source | n near Elba | per-point quality | geoid | chain | cost |
|---|---|---|---|---|---|
| **2008 DNR control, open (L1O)** | 4 within 14 km measurable locally; 76 L1O in Winona alone | spread 0.11–0.34 m; **same geoid as gen1 (GEOID03) – no conversion** | GEOID03 (dataset-level only) | **none** | **free – PDFs** |
| **2021 LCPs** | 3 within 4.5 km, all on Elba lines, **all three tiles on disk** | spread 0.41–0.67 m | GEOID18, verified twice | **none** | **free** |
| 2021 NVA/VVA on Elba lines | 13 more within 25 km, ±24 km northing | 40.8 mm RMS floor, siting varies | GEOID18 | none, if the line assignment holds | ~25 MB per tile |
| Goodhue 2020 | 46, all 20–90 km W | unknown | **GEOID12B – conversion required** | n/a | small |
| NGS leveled + flat | **1 in Winona (DG8385)**; 6 Tier-1 at 37–76 km | geoid-independent, but 3 of 6 are ±180 m horizontally | none needed | n/a | small |
| MnDOT geodetic | 621 with NAVD88 in Winona | rods and monument tops – not lidar-measurable | per-station datasheet | n/a | 1 PDF per mark |
| CORS, locks and dams | – | ARP offsets unpublished; 1912 river datum | – | – | **reject** |

### 7.3 The single highest-value acquisition

**The remaining six 2008 county validation reports** – dodge, fillmore,
freeborn, houston, mower, steele – roughly **4 MB of PDF**, completing a set of
about **1 033 surveyed control points across gen1's whole footprint**, of which
I have already parsed 398.

Why this beats every lidar tile:

1. **It is gen1's own control, on gen1's own geoid.** No GEOID03→GEOID18
   conversion enters the comparison, no swath chain, no flight-line
   extrapolation. Every other candidate needs at least one of those.
2. **It is the only source that supplies both things at once** – density near
   Elba (11 points within 5 km, 36 within 10 km) *and* the tens-of-km spread the
   tilt question needs. The 3DEP marks give spread only along flight lines; the
   NGS marks give spread only 37 km away and clustered.
3. **It already produced a result** (§5): within 10 km, gen1's control residual
   shows **no significant tilt** (−1.32 ± 3.06 and −0.40 ± 2.67 mm/km), which is
   independent evidence against the DoD's −14 mm/km being regional.
4. **It costs almost nothing and touches no point cloud.**

**Second, and it is genuinely a close second: the email in §5.2** – specifically
request 3, the **value of the bias adjustment applied to gen1 in 2008**. That is
one number, it is not recoverable by any amount of our own processing, and it
bears directly on why gen1's datum floats. Nothing else on this list is
irreplaceable in that way.

**Third, if a lidar acquisition is wanted:** the ~16 gen1 tiles under the 36
control points within 10 km (~400 MB; 8 of the 24 tiles involved are already on
disk). Fetch them tile-by-tile, one at a time, and screen by §7.1 before running
any estimator.

### 7.4 Count, siting, or spread?

**None of the three was binding. What was binding is that two of the three source
documents had not been fully read** – the 3DEP survey report (which holds 143
control points absent from every shapefile) and the MnGeo county validation
reports (which hold gen1's own checkpoint coordinates, contradicting the earlier
"dead end" finding).

Having read them, in order of what now limits the answer:

1. **Siting is binding, and it is now measurable.** The 40.8 mm floor is a
   siting artefact – surveyed points on road crowns, embankments and local highs,
   where the estimator's answer depends on the radius. On genuinely flat open
   ground the local spread drops by a factor of five (0.110 m at L1O101 versus
   metres at 2210). **Screening by measured local spread will buy more than
   quadrupling the mark count.**
2. **Spread is binding for the tilt, and it is now solved on paper.** The 2008
   control gives ±30 km coverage for free; the on-line 3DEP marks give ±24 km at
   one tile each. Both clear the ≤ 4.7 mm/km needed.
3. **Count was never binding.** 321 3DEP points and ~1 033 gen1 control points
   sit inside gen1's footprint. Count is now so abundant that the risk has
   inverted: the danger is averaging in badly-sited marks and manufacturing
   precision. **Screen first, then count.**

### 7.5 What to do next, in order

1. **Today, no downloads.** Run the project's slope-normal estimator plus the
   proper geoid handling at: **L1O101** (1.34 km, open, spread 0.110 m),
   **LCP 1079** (2.18 km, both epochs on disk), and the other two in-grid
   control points. Compare against the existing +22.7 ± 39.7 mm and reconcile
   the sign conventions – §5 flags a possible sign tension that must be resolved
   before anything else is believed.
2. **Fetch the six remaining validation reports** (~4 MB) and parse them with
   the §5 recipe.
3. **Measure the §7.1 screening statistics** for every candidate that falls in a
   local tile, and choose the threshold from that distribution.
4. **Send the §5.2 email.**
5. **Only then** fetch gen1 tiles, one at a time, for the screened survivors.

---

## Appendix A – what was fetched (be gentle with the server)

All requests spaced 2–3 s apart. **No bulk lidar was downloaded**; every point
cloud touched was already on disk.

| what | host | bytes |
|---|---|---|
| `Vertical_Accuracy/` directory listings (5 dirs) | rockyweb.usgs.gov | ~40 kB |
| `MN_Driftless_NVA_VVA_UTM15_QL0` + `QL1` (shp/shx/dbf/prj) | rockyweb | 274 kB |
| `USGS_..._QL0_VA.txt`, `..._QL1_VA.txt` | rockyweb | 84 kB |
| **`MN_SE_Driftless_2021_B21_Ground_Control_Survey_Report.pdf`** | rockyweb | 1.67 MB |
| `LCP.zip` – **tail 200 kB only** (central directory), then 6 photos by byte range | rockyweb | 1.74 MB of a 98.7 MB archive |
| USGS metadata project index + 5 adjacent-project dir listings | rockyweb | ~60 kB |
| `GoodhueCo_MN_NVA_VVA` + `GoodhueCo_MN__Ground_Control` (shp/shx/dbf/prj) | rockyweb | 92 kB |
| `USGS_MN_GoodhueCounty_2020_A20_ProjectReport.pdf` | rockyweb | 188 kB |
| 8 county tile lists | resources.gisdata.mn.gov | ~35 kB |
| **Winona / Wabasha / Olmsted county validation reports** (PDF) | resources.gisdata.mn.gov | 1.90 MB |
| `WINO1.csv` (MnDOT geodetic, Winona) | olmweb.dot.state.mn.us | 76 kB |
| NGS bounds API, 15 queries; 8 bulk datasheet POSTs (126 sheets); 13 single-mark probes | geodesy.noaa.gov | *(sub-agent)* |
| MnDOT / MnGeo / DNR / USACE pages, L&D 5 water control manual | various | *(sub-agent; the manual is 6 MB)* |
| tile centroid index, EPT project boundaries | **cached locally** | 0 requests |
| gen1 at 1078/1079/1080 and at 17 control points; gen2 at 1079 | **local disk** | 0 requests |

**Fetched by me directly: about 6.2 MB.** Sub-agent fetches are marked as such.

Scratchpad artefacts (session-scoped, not committed): `va/` (both checkpoint
shapefiles + VA reports), `survey_report.pdf`/`sr.txt`, `sr_points.pkl` (532
parsed 2021 points), `winona_val.pdf`/`wabasha_val.pdf`/`olmsted_val.pdf` and
`gen1_control.pkl` (398 parsed 2008 control points), `gen1_tiles.pkl` (1 529
tiles), `gen1_cps.pkl` (321 in-coverage 2021 points), `local_control_meas.pkl`
(17 measured), `gh/` (Goodhue), `lcpimg/` (6 site photographs), `tilelists/`,
`WINO1.csv`.

## Appendix B – verification status of every claim

**Verified by me, this session:**

| claim | how |
|---|---|
| 534 control points in the 2021 project; 143 are LCPs, absent from every shapefile | survey report §1.3; 532 parsed from its tables; both shapefiles read (157 + 238 rows, no LCP codes) |
| all §1 points are NAVD88 / **GEOID18** / NAD83(2011) 2010.00 | report §1.8.4 **and** the shapefile `geoid` attribute – two independent sources |
| LCPs are static GPS ≥30 min; NVA/VVA are 180 s RTK | report §1.8.1, §1.8.2 |
| 321 of 532 inside gen1 coverage | 8 county tile lists × cached centroid index, measured tile size |
| 1079 inside the elbaext grid, line 138 only, 42 class-2 returns at R = 5 m | chunked read, local gen1 tile |
| gen2 median − surveyed at 1079 = −11 to −26 mm, stable over R = 2–10 m | chunked read, local gen2 file |
| 1080 on asphalt, 1078/1079 on gravel – **and 1080 sits at p82, not best-sited** | photographs viewed **and** measured against the local tiles (§1c) |
| Goodhue 2020 is on **GEOID12B**, county-foot coordinates | project report; `.prj` |
| **2008 gen1 checkpoint coordinates are public** in the county validation reports | 3 PDFs downloaded and parsed; parsed RMSE reproduces published per-county RMSE in all three |
| `Error = Control Z − Surface Z` | checked to 0.0000 m on every row |
| 4 control points inside the elbaext grid; L1O101 gen1 median − control = −1 mm at R = 5 m | chunked read of `elbaext_gen1_merged.laz` |
| `elbaext_gen1_merged.laz` is unmodified gen1 (no elevation touched) | read `scripts/merge_gen1_tiles.py` |
| 17 local control points measured: mean +57 mm, sd 107 mm | chunked reads of 13 local tiles |
| MnDOT `WINO1.csv` exists: 1 017 rows, 621 with NAVD88, **no geoid or height-source field** | file downloaded and parsed |
| 13 gen1 tiles already on disk, including all three LCP tiles | `ls data/before/` |

**From the sub-agents, not independently re-verified by me:** the NGS filter
cascade and its 18 survivors (§3); MnDOT's 8-county aggregate counts and the
datasheet field structure (§4); the MnGeo/DNR contact names (§5); the CORS table
and the Lock & Dam 5 datum finding (§6).

**Explicitly unverified or reasoning-only:**

| claim | status |
|---|---|
| line assignment for marks 10–25 km from Elba | **extrapolation.** Heading 179–183° displaces a track 350–1 050 m over 20 km. Must be probed per tile. |
| the 2008 control's datum is NAVD88/GEOID03 | **dataset-level metadata only.** The validation reports state no datum. Names carry `RTK`/`VRS`, so GPS-derived, not leveled. |
| σ values in §1.7 | **design calculation**, assuming 40.8 mm independent per mark and a planar model – not a measurement |
| "the tilt cannot be regional" | **reasoning** from the published 161 mm Winona RMSE, plus §5's measured near-zero gradient within 10 km |
| gen1 raw medians at 1078/1079/1080 vs the 2021 survey | **raw** – no geoid conversion, no swath correction, no slope-normal estimator. **Not ties.** |
| the sign relation between §5's −22 mm intercept and the +22.7 mm anchor | **unresolved.** Different surfaces, different geoids, and I did not check our anchor's sign convention. **A task, not a result.** |
| NGS surface extents in §3 | inferred from setting codes and wording, never measured – and §1c is a warning about exactly that kind of inference |

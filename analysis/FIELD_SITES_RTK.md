# Candidate RTK-GNSS survey sites, Whitewater WMA, Elba tile (2026-10-18)

For Andy's geomorphology class. Script `analysis/field_sites_rtk.py` (public-land layers
from `analysis/fetch_public_land.py`); table `data/derived/field_sites_rtk.csv`; map
`figures/field_sites_rtk.png`. Every number below is the script's stdout. Public land only:
the two Whitewater WMA units that intersect the tile; no state park statutory boundary
intersects it. DNR says the WMA boundaries are for general reference; verify at the gate.

## What the survey measures, by site class

* **F1 sedge-correction test** (flat floodplain, most-vegetated pooled-span quartile, read
  erosion, was corrected). At each point: RTK on the GROUND (probe through the mat) and on
  the MAT TOP, late autumn. Ground vs gen2 tests gen2's floodplain surface (expected ~+20
  mm high); mat top minus ground is the sedge offset the logistic encodes.
* **F2 zero penetration** (no 2008 ground return; the label was "closed canopy" until Andy, 2026-10-20: "Zero penetration works. I was getting confused and thought you were discussing forests") (the reed-canary-grass-type cells). Same two
  measurements: this is the species-and-season offset that the marsh literature calibrates
  by survey and that cannot come from the lidar. F2-1 and F2-3 are the NE patch's public
  part.
* **F3 bare controls** (least-vegetated, |DoD| < 30 mm): the survey's zero point; RTK must
  agree with gen2 to within its own precision here.
* **H2 steep slopes, gen1 reads low** (> 25 deg, DoD < -LoD): the near-nadir steep-slope
  effect; RTK ground decides which epoch is right. Forested: leaf-off, long occupations.
* **H3 open hillslope controls** (10-35 deg, open cover, |DoD| < LoD/2): hillslope zero.
* **H1 forest hillslope, gen2 reads high**: at the stated criterion (DoD > +LoD, >= 40
  cells) there are NO clusters on WMA land; forest hillslope cells there number 2,840 with
  median DoD +63 mm. Relaxed to DoD > +50 mm and >= 40 cells there are 5 (largest 92 cells
  at E 576537 N 4882490, median +120 mm, slope 29 deg). Whether to relax is Andy's call.

## Sites (UTM 15N, EPSG:26915; lat/lon WGS84; DoD in mm, raw = star-only product)

    site  cells   ha   easting   northing   lat        lon         rawDoD  corrDoD  span  slope  penfail  parking m
    F1-1    911  2.28  579481.1  4883152.8  44.097243  -92.007023   -179     -61    423   0.6     0%     521
    F1-2    300  0.75  579639.6  4883352.4  44.099023  -92.005012   -179     -49    412   0.6     0%     316
    F1-3    258  0.64  579716.4  4883513.1  44.100461  -92.004028   -166     -40    418   0.7     0%     199
    F1-4     88  0.22  579892.4  4883868.8  44.103643  -92.001776   -197     -33    508   1.0     0%     369
    F2-1    680  1.70  579959.4  4885620.6  44.119406  -92.000674   -271    -177    282   0.6   100%     281
    F2-2    608  1.52  579423.1  4882397.4  44.090449  -92.007860   -225    -192    211   0.5   100%     230
    F2-3    463  1.16  579812.5  4886050.8  44.123295  -92.002444   -192    -141    226   0.3   100%     547
    F2-4    324  0.81  579254.0  4882549.1  44.091833  -92.009950   -254    -153    279   0.5   100%      93
    F3-1    725  1.81  579940.7  4884134.1  44.106027  -92.001133     -2      -2     91   0.4     0%     534
    F3-2    249  0.62  579816.0  4884004.2  44.104871  -92.002711     -8      -7     95   0.5     0%     411
    F3-3     97  0.24  580001.1  4884367.6  44.108122  -92.000343     15      15     86   0.5     0%     297
    F3-4     64  0.16  579266.8  4883155.9  44.097295  -92.009699     16      16     93   0.8     0%     599
    H2-1    163  0.41  575876.1  4885946.5  44.122775  -92.051652   -154      --    266  31.8    34%     623
    H2-2    146  0.36  579217.7  4886175.4  44.124481  -92.009858   -156      --    227  28.3    43%     455
    H2-3    135  0.34  577786.6  4886159.1  44.124488  -92.027746   -150      --    258  30.5    25%     700
    H2-4     75  0.19  579313.0  4884479.3  44.109202  -92.008922   -156      --    238  27.2    35%     711
    H3-1    100  0.25  579959.6  4882781.1  44.093845  -92.001102     -6      --    140  21.5     3%     816
    H3-2     98  0.24  578713.4  4886124.7  44.124079  -92.016168      8      --    191  16.5     6%     750
    H3-3     55  0.14  579523.4  4884563.2  44.109935  -92.006281    -43      --    165  11.0     5%     487
    H3-4     55  0.14  579774.3  4884132.4  44.106030  -92.003211     10      --    350  14.8     4%     505

A compact day: F1-3, F1-2, F3-2, F3-1, F1-4, F3-3, H3-4, H3-3 all lie within ~600 m of the
WMA parking lot at E 579999 N 4884665 along the east-side road; F2-4 and F2-2 are 90-230 m
from the lot at E 579205 N 4882471 near Elba; F2-1 and F2-3 (the NE patch) are ~300-550 m
from the lot at E 579979 N 4885340.

## The datum tie

* **1079_2021_MN** (2021 project LCP calibration mark): E 579659.824 N 4885857.071,
  lat 44.121567 lon -92.004381, 219.366 m NAVD88 (GEOID18) -- ON WMA LAND, 250 m from F2-1.
  Occupy it: it ties the RTK directly to gen2's datum.
* **L1O101** (2008 open control mark): E 578790.555 N 4882696.274, 223.291 m NAVD88
  (GEOID03) -- NOT on WMA land; ask before occupying.
* RTK via MnCORS network RTK gives NAD83(2011) ellipsoid heights; convert with GEOID18 to
  compare with gen2 directly. gen1 is GEOID03; the pipeline handles that difference
  (`cross_epoch_datum` in corrections.json), so compare RTK to gen2 and let the DoD carry
  gen1.
* Expect RTK vertical precision of 20-30 mm per point; the floodplain offsets to be
  measured are 100-270 mm, so a few dozen points per site suffice; the hillslope effects
  (~150 mm at H2) likewise; a ~12 mm ridge-crest effect is NOT resolvable and is not here.

## Caveats, stated once

* Public-land status is from the DNR layer dated 2026-10-18; the layer is "for general
  reference only". Verify on site.
* F2 sites have hay use in the imagery: the autumn thatch may be absent in a cut year, in
  which case the survey still gives the ground (for gen2) but not the offset.
* H sites are under forest: leaf-off season, multi-constellation RTK, longer occupations,
  and expect some points to fail to fix.

## Valley-wide extension (2026-10-20)

The same selection for the whole Whitewater valley, from imagery only (no lidar off the
Elba tile yet), is `analysis/PLAN_VALLEY_FIELD_SITES.md` section 8: 336 clusters on public
land in five classes, `data/derived/valley/valley_sites.csv`, map `figures/valley_sites.png`.
The Elba sites above keep the lidar criteria; the valley sites carry the per-year NE-like
class, spring vigour, and distance to WMA parking and roads.

## One field day (Andy, 2026-10-20: "We have one day in the field just now")

Produced by `analysis/rcg/field_day_plan.py` (log `data/derived/valley/field_day_plan.log`,
table `data/derived/valley/field_day_sites.csv`, map `figures/field_day.png`). Every WMA
parking lot (53) was scored by the cover types reachable within a straight-line walking
radius (600 and 1000 m; the radius is mine), using the valley imagery clusters (post-hoc
view, D4), the 20 Elba lidar sites above, and the five 2021 control marks on public land.

**The minimal set: two lots 1.1 km apart on the WMA road north of Elba, both on the Elba
tile.** Lot 502 (E 579942.8 N 4884848.7; 44.112459 N 92.000998 W, pyproj EPSG:26915 -> 4326)
reaches nine types within 1 km: C, F1, F2, F3, H2, H3, R1, R2, S, including the
both-epoch sward. Lot 501 (E 579359.8 N 4885743.3; 44.120576 N 92.008147 W) reaches eight (adds E, lacks
F1, F3) and has mark **1079_2021_MN** 321 m away (E 579659.8 N 4885857.1; 44.121567 N
92.004381 W; LCP) for the datum tie. Two
further lots sit on the same road between them.

Nearest two sites of each type per lot (distance from the lot; NE-like % 2008/2021/2025/2026;
spring and October 2026 NDVI; largest 2026 NDVI drop and its day-of-year; Elba lidar: raw
DoD mm, gen1 span mm, slope deg):

| lot | type | site | m | ha | E | N | NE08/21/25/26 | spr/Oct NDVI | drop (DOY) | DoD / span / slope |
|---|---|---|---|---|---|---|---|---|---|---|
| 502 | C | C-v29 | 418 | 1.32 | 580125 | 4884472 | 18/5/0/0 | 0.18/0.16 | 0.25 (240) | |
| 502 | C | C-v61 | 722 | 0.70 | 579774 | 4884146 | 10/17/0/0 | 0.15/0.20 | 0.11 (240) | |
| 502 | F1 | F1-4 | 981 | 0.22 | 579892 | 4883869 | | | | -197 / 508 / 1.0 |
| 502 | F2 | F2-1 | 772 | 1.70 | 579959 | 4885621 | | | | -271 / 282 / 0.6 |
| 502 | F3 | F3-3 | 485 | 0.24 | 580001 | 4884368 | | | | 15 / 86 / 0.5 |
| 502 | F3 | F3-1 | 715 | 1.81 | 579941 | 4884134 | | | | -2 / 91 / 0.4 |
| 502 | H2 | H2-4 | 730 | 0.19 | 579313 | 4884479 | | | | -156 / 238 / 27.2 |
| 502 | H3 | H3-3 | 507 | 0.14 | 579523 | 4884563 | | | | -43 / 165 / 11.0 |
| 502 | H3 | H3-4 | 736 | 0.14 | 579774 | 4884132 | | | | 10 / 350 / 14.8 |
| 502 | R1 | R1-v149 | 442 | 0.74 | 580056 | 4885276 | 97/1/100/100 | 0.23/0.49 | 0.08 (200) | |
| 502 | R1 | R1-v140 (both epochs; the NE training field) | 899 | 11.50 | 579981 | 4885747 | 97/91/100/100 | 0.27/0.57 | 0.08 (240) | |
| 502 | R2 | R2-v218 | 764 | 0.69 | 580675 | 4885067 | 38/77/3/97 | 0.22/0.28 | 0.09 (240) | |
| 502 | S | S-v324 | 957 | 0.51 | 579465 | 4885678 | 49/37/0/0 | 0.25/0.36 | 0.12 (275) | |
| 501 | C | C-v52 | 129 | 1.09 | 579489 | 4885742 | 26/6/0/0 | 0.15/0.21 | 0.23 (275) | |
| 501 | E | E-v131 | 648 | 0.61 | 579367 | 4886391 | 90/7/0/0 | 0.19/0.49 | 0.07 (240) | |
| 501 | F2 | F2-3 | 547 | 1.16 | 579813 | 4886051 | | | | -192 / 226 / 0.3 |
| 501 | H2 | H2-2 | 455 | 0.36 | 579218 | 4886175 | | | | -156 / 227 / 28.3 |
| 501 | H3 | H3-2 | 751 | 0.24 | 578713 | 4886125 | | | | 8 / 191 / 16.5 |
| 501 | R2 | R2-v237 (both epochs) | 394 | 0.54 | 579681 | 4885972 | 56/93/100/0 | 0.16/0.37 | 0.11 (240) | |
| 501 | R2 | R2-v196 | 657 | 1.48 | 579708 | 4886300 | 6/100/100/0 | 0.17/0.36 | 0.11 (275) | |
| 501 | S | S-v259 | 744 | 2.85 | 579273 | 4886483 | 90/39/0/0 | 0.23/0.39 | 0.10 (275) | |

What the day gives: the ground and the sward top under thick RCG-type sward in a field that
was NE-like in 2008, 2021, 2025 and 2026 (R1-v140, 11.5 ha, uncut in 2026: largest NDVI
drop 0.08), a sward that stood in 2008 and 2021 but was cut in 2026 (R2-v237, drop 0.11 is
small too), a late-green non-NE-like meadow (S), bare and cropped controls (C, F3), the
Elba sedge and zero-penetration lidar sites (F1, F2), a steep gen1-low slope and an open
hillslope (H2, H3), a forest edge (E), and the mark. Species, stand height and thatch
thickness at each point are the quantities nobody has (sections R, S, T, W).

**Road-first (Andy, 2026-10-20: "Walks over open land can take a long time and include
noxious plants ... Close places to the road are the best").** Exact shapely distances from
each site's nearest pixel to the WMA road (the lot-501 list, recomputed 2026-10-20; the
table's `dist_parking_m` was a raster artefact and is being fixed in `valley_rank.py`):

| site | type | m from WMA road | nearest lot (m from site) | note |
|---|---|---|---|---|
| C-v52 | C bare / cropped | 1 | 501 (40) | |
| S-v324 | S late-green meadow | 43 | 501 (52) | |
| mark 1079_2021_MN | datum tie | 41 | 501 (321) | |
| R2-v237 | R2, sward 2008 + 2021, not 2026 | 62 | 501 (355) | both epochs |
| R1-v140 | R1, sward 2008/2021/2025/2026, uncut 2026, 11.5 ha | 154 | 204 (69 from the field edge) | both epochs; the NE training field |
| F2-3 | F2 zero penetration | 246 | 501 (547) | Elba lidar |
| H3-2 | H3 open hillslope control | 303 | 501 (751) | Elba lidar |
| H2-2 | H2 steep gen1-low slope | 456 | 501 (455) | Elba lidar |
| E-v131 | E forest edge | 501 | 501 (648) | the only type > 400 m off the road |

So the day is one road: park at lot 501 for C, S, the mark and R2-v237 (all within ~60 m of
the road), drive 500 m south to lot 204 for the big R1 field (its edge is 69 m from that
lot), and take F2-3 and H3-2 as 250-300 m walks off the road. Everything else is optional.

**Hillslope gen2-high patches by lot 204 (Andy, 2026-10-20: "where we saw issues on
hillslopes with low veg on the [gen2] lidar biasing those returns high").** gen2 is the 3DEP
flight of 2021-05-01 (green-up), not 2024. On the Elba tile's public hillslopes (10-35 deg,
not floodplain, 169,697 cells) the DoD rises with PyForestScan canopy cover: median +1.8,
+6.8, +10.8, +19.6 mm by cover quartile (edges 0.22 / 0.28 / 0.35), share above +LoD 6.9,
7.2, 7.5, 9.6%. Clusters of DoD > +LoD (>= 20 cells) within 300 m of lot 204, nearest cell:

| ha | m from lot 204 | m from road | E | N | lat | lon | DoD (mm) | slope | canopy cover |
|---|---|---|---|---|---|---|---|---|---|
| 0.06 | 124 | 446 | 579862 | 4885382 | 44.117273 | -92.001920 | +146 | 20 | 0.21 |
| 0.12 | 137 | 481 | 579842 | 4885338 | 44.116870 | -92.002177 | +170 | 23 | 0.29 |
| 0.61 | 166 | 239 | 579868 | 4885462 | 44.117993 | -92.001846 | +147 | 24 | 0.27 |
| 0.18 | 192 | 526 | 579798 | 4885278 | 44.116335 | -92.002748 | +164 | 21 | 0.26 |
| 0.58 | 257 | 383 | 579738 | 4885252 | 44.116117 | -92.003502 | +147 | 27 | 0.34 |

and the hillslope control beside them: 12.46 ha of |DoD| < LoD/2 (median +19 mm, slope
22 deg, canopy cover 0.27) 271 m from lot 204 and 102 m from the road (nearest cell
E 579738 N 4885462; 44.118007 N 92.003470 W). The gen2-high patches and the control have
the SAME canopy cover (0.21-0.34 vs 0.27): cover does not separate them, which is what
`gen2-vegetation-offset` found at the marks (canopy cover rho -0.053; the near-ground
composition `lowveg` is what predicts the offset, -268 mm per unit). There is NO lowveg grid
for the Elba tile yet (the only one, `data/derived/whitewater/lowveg_grid.npy`, is the
Weaver-end tile 10 km NE); `scripts/make_lowveg.py --gen2` on the elbaext class-2 cloud
makes it, a streaming pass, queued behind the machine pause.

**gen2 lowveg on the Elba hillslopes (built 2026-10-20, Andy: "Go ahead and find those
patches").** `scripts/make_lowveg.py --tile elbaext_independent --gen2
data/after/elbaext_3dep_fulldensity.laz` wrote `data/derived/elbaext_independent/lowveg_grid.npy`
(229,022,726 returns placed; 719,424 of 720,900 cells; lowveg p10 0.000, median 0.050, p90
0.413). On the 169,697 public hillslope cells (10-35 deg, not floodplain), lowveg quartiles
0.068 / 0.161 / 0.339 (p90 0.571), and the DoD by lowveg bin (log
`data/derived/valley/hillslope_lowveg_patches.log`):

| lowveg bin | median DoD (mm) | share > +LoD | n |
|---|---|---|---|
| 0.000-0.068 | 3.7 | 5.2% | 42,412 |
| 0.068-0.161 | 11.0 | 6.7% | 42,416 |
| 0.161-0.339 | 15.0 | 8.2% | 42,444 |
| 0.339-0.571 | 13.6 | 11.6% | 25,453 |
| 0.571-1.000 | 3.9 | 10.4% | 16,972 |

Linear DoD = 15.7 - 40.9 * lowveg mm, r -0.066; canopy cover vs lowveg r 0.035. THE RESULT:
the hillslope DoD on this tile does not carry the -268 mm per unit lowveg offset that the
389 surveyed marks showed for gen2 alone (`gen2-vegetation-offset`). DoD = gen2 - gen1,
so a gen2 offset shows in the DoD only if gen1 does not share it; a gen1 that reads high in
the same vegetation (`gen1-cannot-see-its-own-vegetation`, the floodplain sward) cancels it.
That is a hypothesis; the RTK transect is what tests it, because the survey compares each
epoch to the ground separately. Patches of high lowveg (>= p75) with DoD > +LoD within
800 m of the lots exist but are small (0.04-0.05 ha; e.g. E 579822 N 4885328, 157 m from
lot 204, lowveg 0.47, DoD +166); the full list with low-lowveg controls is
`data/derived/valley/hillslope_lowveg_patches.csv` (cut-offs mine, for locating only).

**The class transect (Andy: "a hillslope transect that the class could walk to look for
spatial variability in land-surface elevation to compare against the lidar data").** Of 47
walkable 300 m lines near the two lots (all on public land, max slope <= 25 deg, mean >= 8,
<= 20% floodplain; cut-offs mine), the one chosen: start E 579259.8 N 4885643.3 (44.119686 N
92.009412 W), 141 m west of lot 501, bearing 300 deg, 300 m, never more than 115 m from the
WMA road; rises 34 m; median slope 7 deg, max 25. It crosses a vegetated toe-slope (25-50 m:
lowveg 0.54-0.79, DoD -274 then +110 mm), a bare bench (125-200 m: lowveg 0.000, DoD +26 to
+70 mm) and a vegetated upper slope (225-300 m: lowveg 0.17-0.42, DoD +59 to +108 mm).
Waypoints every 25 m with lat/lon: `data/derived/valley/hillslope_transect_lot501.csv`;
profile `figures/hillslope_transect_lot501.png`; producer `analysis/hillslope_transect.py
--start 579259.8 4885643.3 --az 300 --length 300`. RTK every 10 m on the ground, plus the
vegetation top where there is one, gives both epochs' offsets as a function of position.

**Site maps (Andy, 2026-10-20: "Make site maps of where we should park and go").**
`analysis/rcg/field_day_maps.py` -> `figures/field_day_lot204.png` (stop 1: the RCG field
R1-v140 and R1-v149, F2-1, the bluff lowveg stand coloured by slope, its least-steep entry
route, the gen2-high patches, the hillslope control) and `figures/field_day_lot501.png`
(stop 2: C-v52, S-v324, mark 1079, R2-v237, F2-3, the class transect with 25 m dots), both
on 2025 FSA imagery with the WMA road and a 200 m bar; and
`data/derived/valley/field_day_waypoints.gpx` (27 waypoints, the transect and the bluff
route as tracks) for a handheld GPS. Seen on the imagery: the bluff lowveg stand lies
UNDER tree canopy although PyForestScan's canopy cover reads 0.23 there (the voxel
metric under-reads, `gen2-vegetation-offset`), so expect a wooded slope, not shrubland.

**Navigating to the regions with Emlid RTK and the free Emlid Flow (Andy, 2026-10-20).**
Emlid's documentation (read 2026-10-20: docs.emlid.com "Create or import project"; the
Emlid Flow plan comparison as summarised by two dealers) says the free plan imports and
exports CSV, DXF, Shapefile and KML, collects and stakes out POINTS, and shows them on the
app's default street map; lines and polygons, and switching the background to satellite or
a WMS layer, need the paid Survey plan. So:

1. `analysis/rcg/field_day_emlid_export.py` turns every region outline into a ring of named
   points every ~20 m plus a CENTRE point (`R1v140-01..`, `R1v140-CENTRE`, ...), with the
   lots, the mark and the transect's 25 m points: `data/derived/valley/field_day_points_emlid.csv`
   (Name, Longitude, Latitude, WGS84), the same as `.kml` and as a point Shapefile. Import
   one of these into the Flow project; the rover's position then sits among the region's
   outline points on the map screen, and "Stake out" to a CENTRE or an outline point gives
   distance and bearing. Being inside the ring of points is being inside the region.
2. A phone beside the rover shows the true polygons: `field_day_regions_phone.kml` (polygons,
   centres, waypoints) or `field_day_waypoints.gpx` (outlines as tracks) in any free map app
   that reads KML/GPX over imagery (Google Earth, OsmAnd, Avenza, Gaia); phone GPS at a few
   metres is enough to walk into a 0.5-12 ha region, and the Emlid does the measuring.
3. Printed `figures/field_day_lot204.png` and `figures/field_day_lot501.png` as the
   fallback; each has the WMA road, the lots and a 200 m bar.
Project CRS in Flow: NAD83(2011) / UTM 15N with GEOID18 heights if the MnCORS network is
used (the plan's datum paragraph); the navigation points are WGS84 lat/lon and import into
any project CRS.

**The day (Andy, 2026-10-20):** a meander bend in the morning, the GNSS mapping in the
afternoon, then home. Afternoon order, by what the lidar question needs most: lot 204
first (R1-v140 ground + sward top, F2-1 inside it, R1-v149 by the lot; the bluff stand
only if the slope looks right from the field edge), then lot 501 (mark 1079 first for the
tie, then C-v52 and S-v324 beside the lot, R2-v237, and the transect as the closing walk).

**Monuments for a base station (Andy, 2026-10-20).** Sources: NGS datasheet API
(`geodesy.noaa.gov/api/nde/bounds`, 45 marks in 44.05-44.20 N / 92.10-91.90 W, saved
`data/derived/valley/monuments/ngs_marks.csv`) and the MnDOT Geodetic Control REST service
(`dotapp9.dot.state.mn.us/egis12/rest/services/OLM/mndot_geodetic_agol_operational_lyrs1`,
54 active stations in the box, `mndot_geodetic_box.csv`); datasheets DR9925 and DG8385 read.

| mark | what | position | height | from lot 501 / 204 |
|---|---|---|---|---|
| DR9925 "8508 Y" (MnDOT, set 2023) | flange-encased aluminium rod, 18.9 m deep, removable disk, 37.9 ft SW of TH 74 | NAD 83(2011) ADJUSTED, network accuracy 0.61 cm horizontal, 1.49 cm ellipsoid height (95%) | ellipsoid 188.669 m; NAVD 88 218.75 m (GPS-derived); stability B | 0.31 / 0.63 km |
| DG8385 "9 DRL" (USGS 1972) | survey disk in a flagpole base / concrete footing, by TH 74 | HD_HELD1: hand-held GPS, +/- 3 m | NAVD 88 223.352 m, leveled 2nd order class I; stability C | 0.37 / 0.45 km |
| DR9926 "8508 F RESET" | rod, vertical only published | hand-held horizontal | NAVD 88 221.492 m, 2nd order | 1.28 / 0.70 km |

DR9925 is the base-station monument: a published NAD 83(2011) position and ellipsoid height
at the centimetre, "suitable for satellite observations", 310 m from lot 501 and within
1.3 km of every site, in the public TH 74 right-of-way. Base over it with the published
coordinates entered as the base position in Flow gives rover points in NAD 83(2011)
directly, NAVD 88 by GEOID18 -- the same frame as gen2. DG8385 is the leveled benchmark
of `leveled-benchmark-cross-epoch` (gen2 - gen1 = +73 mm there): not a base (3 m
horizontal) but the best vertical CHECK shot of the day. 1079_2021_MN is the lidar
vendor's control, a check point, not a monument. Alternative that needs no base: MnCORS
network RTK over NTRIP, if there is cell data in the valley (UNVERIFIED); then occupy
DR9925 as the check instead.
Andy, 2026-10-20: the buildings at the junction are the DNR (Whitewater WMA) headquarters, so
a base left on DR9925 beside them "is probably safe"; the base goes there for the afternoon.

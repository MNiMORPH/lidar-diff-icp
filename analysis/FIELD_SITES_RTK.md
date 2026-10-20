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
* **F2 closed canopy, zero penetration** (the reed-canary-grass-type cells). Same two
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
Elba sedge and closed-canopy lidar sites (F1, F2), a steep gen1-low slope and an open
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
| F2-3 | F2 closed canopy, zero penetration | 246 | 501 (547) | Elba lidar |
| H3-2 | H3 open hillslope control | 303 | 501 (751) | Elba lidar |
| H2-2 | H2 steep gen1-low slope | 456 | 501 (455) | Elba lidar |
| E-v131 | E forest edge | 501 | 501 (648) | the only type > 400 m off the road |

So the day is one road: park at lot 501 for C, S, the mark and R2-v237 (all within ~60 m of
the road), drive 500 m south to lot 204 for the big R1 field (its edge is 69 m from that
lot), and take F2-3 and H3-2 as 250-300 m walks off the road. Everything else is optional.

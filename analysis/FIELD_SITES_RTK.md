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

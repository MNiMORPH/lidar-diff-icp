# RCG phenology: the record (Phases 0-2 and 4 of PLAN_RCG_REMOTE_SENSING.md)

**Date:** 2026-10-19. Scripts in `analysis/rcg/`; layers and tables in
`data/derived/elbaext_delong_star_nosurf/rcg/`; figures `figures/rcg_phenology_2024.png`,
`figures/rcg_phenology_2008.png`. Every number is pasted from the scripts' stdout.

## Phase 0 -- scenes (done)
23 scenes as masked NDVI GeoTIFFs in EPSG:26915 (`fetch_scenes.py`, index `scenes.json`).
Landsat-5 2008-10-29: clear 1.000 over tile and floodplain, snow 0.0000, gap-free. Sentinel-2
2024: 04-05, 04-15, 05-05, 05-10, 08-03, 08-08, 09-27, 10-02, 10-17, 10-22 fully clear;
10-27 and 11-01 partly clouded (0.71, 0.29 clear). 2023-10-23 clear.

## Phase 1 -- present-day signature (Sentinel-2 2024), null check first

Metric definitions: apr = mean(04-05, 04-15); sum = max(08-03, 08-08); oct = mean(10-02,
10-17); retained = oct/sum; greenup = apr/sum; senesc_doy = first autumn date with NDVI
below half the summer maximum (999 = never by 11-01). Medians over Sentinel-2 pixels (NMAD).

    population                              px    apr    sum    oct   retained  greenup  senesc  oct2023
    KNOWN tree cover (WorldCover 10)   154,261   0.22   0.56   0.38    0.69      0.39     296     0.28
    KNOWN cropland (40)                 64,947   0.12   0.63   0.27    0.47      0.20     291     0.14
    KNOWN grassland (30)                19,305   0.22   0.52   0.39    0.77      0.46     999     0.35
    KNOWN forest hillslope (product)     1,388   0.30   0.50   0.37    0.76      0.62     999     0.32
    KNOWN F3 bare floodplain control     1,597   0.13   0.60   0.21    0.46      0.26     276     0.17
    all flat floodplain                 13,609   0.22   0.57   0.35    0.64      0.39     296     0.27
    F1 sedge cells (corrected)           1,208   0.22   0.59   0.40    0.68      0.38     999     0.29
    F2 closed canopy (penfail)           1,634   0.22   0.60   0.37    0.65      0.39     301     0.29
    NE patch                               601   0.36   0.61   0.18    0.30      0.58     271     0.41

The null check behaves: cropland is bare by October (retained 0.47, NMAD 0.26 from the corn/
soy split), grassland and the forest hillslope never fall below half their summer maximum,
and the F3 "bare" floodplain controls behave exactly like cropland -- they ARE cropped
fields, which is why gen1 saw bare ground there.

The NE patch in 2024: the EARLIEST green-up of any population (April 0.36 against 0.22;
greenup 0.58) and then BROWN by late September (October 0.18, senescence day 271): a hay
field cut in late summer. In October 2023 the same patch was the GREENEST population on
the tile (0.41). Management decides its autumn signal year to year -- the risk the plan
stated. Within F2, the 34 patches with >= 5 pixels split into two behaviours in 2024:
never-senescing, early-and-late green (e.g. 87 px at E 579823 N 4886068: retained 0.90,
greenup 0.74; 42 px: 0.94/0.69; 31 px: 0.89/0.66 -- the reed-canary-grass-type curve) and
early-green-then-cut (the NE main patch, 210 px: 0.24/0.55; 146 px at E 579422 N 4882394:
0.23/0.16). F1 (sedge) sits with grassland: never senesces, retained 0.68.

## Phase 2 -- the 2008 state (Landsat 5/7, 30 m)

Metrics: spr = mean(04-04, 05-06, 05-14, 05-22); sum = max(08-02, 08-10, 08-18); sep =
09-19; oct = 10-29; retained = oct/sum; greenup = mean(04-04, 05-06)/sum. A lidar-defined
class enters a 30 m pixel only if it covers > 50% of it. Levels are not comparable to
Sentinel-2; ranks are.

    population                            px    spr    sum    sep    oct   retained  greenup
    KNOWN tree cover                  17,098   0.49   0.91   0.80   0.41    0.46      0.37
    KNOWN cropland                     7,217   0.23   0.88   0.70   0.26    0.29      0.23
    KNOWN grassland                    2,134   0.43   0.80   0.68   0.40    0.50      0.39
    KNOWN forest hillslope                44   0.55   0.85   0.78   0.60    0.73      0.54
    KNOWN F3 bare control                 84   0.32   0.90   0.77   0.27    0.30      0.30
    all flat floodplain                1,466   0.50   0.88   0.75   0.40    0.46      0.39
    F1 sedge cells                        34   0.51   0.85   0.75   0.41    0.48      0.39
    F2 closed canopy                      98   0.45   0.85   0.71   0.39    0.47      0.37
    NE patch                              61   0.45   0.82   0.68   0.41    0.50      0.35

    percentile of the class's median retained greenness within the flat floodplain, Oct 2008:
      F3 bare 12th   F2 closed canopy 54th   F1 sedge 60th   NE patch 69th

The null check holds again (cropland 0.29 retained, 12th percentile for the cropped
controls). THE ANSWER TO THE QUESTION: on 29 October 2008, two weeks before the gen1
flight, the NE patch and the closed-canopy class were ORDINARY flat floodplain --
NDVI 0.41 against 0.40, 69th and 54th percentile of retained greenness -- and
indistinguishable from the sedge cells (0.48 vs 0.47). They were still moderately green
(a live standing canopy, not a harvested field), but phenology at 30 m does not single
them out. Per patch: the six F2 patches with >= 5 pixels have retained 0.34-0.60, within
the floodplain's range.

## Phase 4 -- the DoD against October-2008 greenness: NEGATIVE

`dod_vs_greenness.py`: 1,466 pixels of 30 m, > 50% flat floodplain; y aggregated from
5 m; x from the imagery, independent of the lidar.

    x (2008)      y = raw DoD   r       slope        DoD by x-quartile Q1..Q4 (mm)
    ndvi_oct08    dod          -0.013   -16 mm/unit   -12   -82   -60    -2
    retained08    dod           0.023    20           -15   -65   -62    -3
    greenup08     dod           0.100   117           -21   -47   -55    -2
    ndvi_sep08    dod          -0.098  -109             0  -101   -34   -13
    ndvi_oct08    gen1 span     0.211   207 mm/unit   160   258   276   237
    (2024 Sentinel-2 October, wrong epoch, for contrast: r -0.132 with raw DoD)

October-2008 greenness does not predict the floodplain DoD: r within +/-0.1 for every 2008
metric, quartile medians non-monotonic. It barely predicts gen1's own span (r 0.21).
The premise of Phase 4(c) -- NDVI as the Beer-Lambert tau for the gen1 vegetation bias
-- FAILS on this floodplain at 30 m.

**Why, as hypothesis.** NDVI measures live green leaf area. The lidar bias comes from
standing STRUCTURE -- senesced sedge and grass, brown but tall and dense -- which an
October NDVI does not see; a cut hay field and a brown sedge mat can share an NDVI and
differ by 20 cm of canopy. The quantity needed is autumn vegetation HEIGHT x density, which
optical greenness is not. (That is also why the lidar itself, which does see structure,
remains the better vegetation proxy for the lidar's own bias -- the circularity of section
I of FLOODPLAIN_RETURN_PDF.md is the price of the only instrument that sees the mat.)

## What stands after Phases 1, 2, 4

* Phenology from imagery CAN separate cropland, grassland and forest here (null checks
  pass in both epochs) and shows early-green-up, never-senescing patches on the floodplain
  in 2024 that look like the reed-canary-grass curve -- but the autumn signal flips with
  cutting from year to year, and in 2008 the lidar's closed-canopy cells were not
  distinctive.
* An October-2008 NDVI is NOT a usable covariate for the gen1 vegetation bias. Phase 4(c)
  is closed; 4(a), a mask class, has no imagery support either, since the closed-canopy
  cells are not separable in 2008 imagery.
* What remains live from the plan: Phase 3 as a SPECIES question (truth from outside the
  lidar; the field day), and Phase 4b, finding reed-canary-grass stands to survey from the
  2024 Sentinel-2 signature across the whole valley -- that uses the early-green-up +
  never-senescing curve, which the tile does show, and needs no lidar.

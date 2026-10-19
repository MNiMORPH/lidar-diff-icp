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

## N. The approach WAS wrong: what the Wisconsin method actually did, and the test redone (2026-10-19)

Andy: "Could it be that your approach was not right? Do we need not just greenness, but the
type of green? How do prior studies differentiate reed canary grass?" Both reports read in
full this time (Bernthal & Willis 2004, pilot; Hatch & Bernthal 2008, statewide):

* **No NDVI, no multi-date.** "They concluded that multiple dates and use of NDVI are not
  needed." One scene from Oct 9-31, <= 10% cloud, ALL reflective bands (Landsat 1,2,3,4,5,7;
  thermal dropped).
* **Wetland mask first** (Wisconsin Wetlands Inventory removes uplands), then ISODATA with 75
  clusters, each cluster assigned to >= 50% / < 50% reed canary grass by field truth, then
  maximum likelihood. 269 sites, 1,388 polygons; overall 61-83%, user's 72-92%.
* **"Band 5, the midIR band (1.55-1.75 um), showed the greatest spectral separation."**
* Confusions: heavily grazed pasture classified as RCG (short-cropped, geometric fields);
  cranberry bogs; a single shrub or standing water in a pixel breaks it (NIR-driven); only
  open-canopy wetlands; minimum mapping unit 0.5 acre (two pixels).

So "the type of green" is the right phrase: the signature is the whole reflectance vector,
and NDVI -- a ratio of two of six bands -- throws away the level and the SWIR. Redone
(`fetch_bands.py`, `band_separability.py`, `dod_vs_bands.py`; figure
`figures/rcg_spectra_october.png`):

**Landsat-5 2008-10-29, median surface reflectance (NMAD):**

    population               blue     green     red      nir08    swir16   swir22   NDVI
    all flat floodplain      0.039    0.064    0.067    0.161    0.203    0.130    0.398
    F1 sedge                 0.042    0.064    0.080    0.193    0.227    0.133    0.411
    NE patch                 0.052    0.090    0.108    0.254    0.271    0.159    0.405
    KNOWN grassland          0.045    0.071    0.082    0.202    0.239    0.147    0.396
    KNOWN cropland           0.064    0.096    0.123    0.225    0.275    0.194    0.256

**Separability d' (difference of medians over pooled NMAD), with KNOWN pairs as the null:**

    pair                                  blue  green   red   nir08  swir16  swir22  NDVI  NDMI
    NULL cropland vs grassland            1.06   1.06  1.14   0.41    0.74    1.00   1.39  0.11
    NULL tree vs grassland                1.35   1.02  1.32   1.22    1.23    1.14   0.16  0.39
    NE patch vs all flat floodplain       1.70   1.71  2.12   2.62    2.25    1.52   0.15  1.61
    NE patch vs F1 sedge                  2.04   2.70  2.10   2.37    1.98    2.39   0.26  2.21
    F2 closed canopy vs flat floodplain   0.67   0.38  1.31   0.67    0.78    0.45   0.20  0.63

In October 2008 the NE patch IS spectrally distinct: brighter in every band at the same
NDVI, d' 2.1-2.6 in red, NIR and SWIR1 against the floodplain -- more separable than
cropland is from grassland in the same scene -- and invisible in NDVI (0.15). Wisconsin-
style k-means (k = 30, six bands, floodplain mask only) puts 87% of the patch in two
clusters of purity 40% and 22% against a 2.2% base rate. Sentinel-2 confirms it both years
(2023: d' 2.6-3.0 in NIR / red-edge, late-green; 2024: d' 2.2-3.7 in visible / SWIR, cut).
The wider F2 class is heterogeneous and only weakly distinct (d' 0.4-1.3).

**What it is, spectrally (hypothesis).** High NIR AND high SWIR at an ordinary NDVI: a dense,
bright, dryish herbaceous canopy -- the standing-structure signal that NDVI cancels. That is
the quantity the lidar bias follows, and it was in the bands all along.

**But the floodplain-wide DoD is still not predicted** (`dod_vs_bands.py`, 1,466 pixels):

    x (2008-10-29)   r with raw DoD   DoD by x-quartile
    swir16              -0.147         -3  -37  -74  -20
    NBR                 -0.154         -3  -24  -60  -60
    nir08               -0.113         -3  -45  -48  -29
    NDVI                -0.013        -12  -82  -60   -2
    six bands, OLS, CV by 500 m block:  CV R2 0.090 (raw DoD), 0.066 (corrected), 0.088 (gen1 span)

The bands do better than NDVI (0.09 against 0.00) but the NE patch is one spectral object,
and across the floodplain the sedge-type bias is not tracked by October reflectance at
30 m. Phase 4(c) stays closed. What reopens: Phase 3 and 4b -- a Wisconsin-style
classification (wetland mask, all bands, clusters, outside truth) is feasible on this
tile and the valley, and it finds the NE-type stands.

## O. 2025 vs 2008: the grass DID change (2026-10-19)

Andy: "Can you check 2025 vs. 2008 reflectance and if the grass changed?" Script
`analysis/rcg/compare_2008_2025.py`; figure `figures/rcg_relative_spectra_2008_2025.png`;
table `rcg/relative_spectra_2008_2025.csv`. Scenes added: Landsat-8/9 2025-09-26 (clear
1.000), 2025-10-04 (0.992), 2025-11-05 (1.000); Sentinel-2 2025-10-02 and 10-27 (clear).

Gauge-free comparison: each population's median reflectance divided by the flat
floodplain's median in the same scene and band (sensors and years are not comparable in
level). The known classes are the null: their relative spectra barely move between years.

    NE patch / flat floodplain     blue   green    red   nir08  swir16  swir22   NDVI (patch / floodplain)
    Landsat-5  2008-10-29          1.32   1.40    1.60   1.58    1.34    1.22    0.41 / 0.40
    Landsat-8  2025-09-26          1.12   1.56    1.22   1.84    1.20    1.07    0.82 / 0.73
    Landsat-9  2025-10-04          1.04   1.40    1.01   2.05    1.15    0.99    0.81 / 0.66
    Landsat-9  2025-11-05          1.12   1.36    1.18   2.03    1.10    0.94    0.69 / 0.56
    null: grassland, 2008 / 2025   1.16-1.25 in every band, both years (stable)

Per pixel, log(2025/2008) minus the floodplain's own log-ratio (0 = changed like its
surroundings), 2008-10-29 -> 2025-10-04, median (NMAD):

    population        blue          green         red           nir08         swir16        swir22
    KNOWN grassland   -0.04(0.25)   -0.04(0.17)   -0.03(0.30)   -0.05(0.21)   -0.02(0.14)   -0.04(0.22)
    KNOWN tree        -0.05(0.20)   +0.01(0.16)   -0.07(0.23)   +0.06(0.19)   -0.04(0.18)   -0.10(0.22)
    F1 sedge          -0.11(0.16)   -0.16(0.15)   -0.28(0.21)   -0.10(0.09)   -0.14(0.07)   -0.11(0.11)
    NE patch          -0.25(0.17)   -0.05(0.09)   -0.46(0.23)   +0.21(0.14)   -0.13(0.11)   -0.23(0.13)

The null classes sit at 0 within +/-0.1. The NE patch does not: red fell by 0.46 in the
log (a factor of 0.63), SWIR2 by 0.23, while NIR ROSE by 0.21 -- the same sign and size on
all three 2025 dates. Sentinel-2 agrees: in October 2023 and 2025 the patch's NIR is 1.5-1.8x
the floodplain's with NDVI 0.42-0.58 against 0.27-0.36; only 2024 (the cut year) differs.

**What changed (measurement, then reading).** In October 2008 the patch was brighter than
its surroundings in EVERY band with the floodplain's own NDVI: the spectrum of a dense,
senesced, high-albedo canopy -- standing dead material and litter. In October 2025 it is
dark in red and bright in NIR with an NDVI well above the floodplain's: a live, green,
dense sward in autumn. The grass changed from a standing-dead sward to a green one, which
is what a field "since planted and used for hay" looks like from above: cut in summer,
regrowing green into October. The sedge cells (F1) moved a little in the same direction
(red -0.28) but far less and from a relative spectrum near 1.0, so they remain ordinary
floodplain in both years.

For the lidar this is the mechanism in one picture: the 2008 surface gen1 read 25 cm up
was a dense standing-dead canopy (bright, dry, closed); the 2021 surface gen2 read as
ground was under short green regrowth.

## P. Mapping the NE-like class in 2008 and 2025, and the change (2026-10-19)

Andy: identify the regions with reed canary grass at each time, and the difference. I had
NOT mapped it before this; what existed was the NE patch's separability. Script
`analysis/rcg/map_ne_like.py`; layers `rcg/nelike_p_{ls_20081029,ls_20251004,s2_20251002}.tif`,
`rcg/nelike_change_2008_2025_30m.tif`; table `rcg/nelike_change_patches.csv`; figure
`figures/rcg_nelike_change.png`.

**What is mapped, stated once.** "NE-like" = spectrally like the NE patch in THAT year's
October scene: two-class maximum likelihood on all bands (Wisconsin's classifier shape),
class 1 = the patch's pixels, class 0 = the rest of the flat floodplain, equal priors, p >
0.5. In 2008 that class is a bright senesced dense sward; in 2025 a green high-NIR late
sward (section O). Both are states reed canary grass takes; neither is a species label,
and there is no field truth. Domain: the floodplain mask.

**Confusion check** (share of each WorldCover class over the whole tile classified
NE-like): 2008 -- tree 1.6%, grassland 9.6%, cropland 4.9%; 2025 Landsat -- tree 0.9%,
grassland 21.1%, cropland 10.3%; 2025 Sentinel-2 -- grassland 26.0%, cropland 10.7%. The
2025 class is a quarter of all grassland: a green dense sward is not rare in October, which
is the Wisconsin confusion (pasture, hay). The 2008 class is more specific.

**Areas on the floodplain mask (241.9 ha), 30 m:**

    NE-like 2008            17.8 ha  (198 px)
    NE-like 2025            20.6 ha  (229 px)     [Sentinel-2 10 m: 22.2 ha]
    kept (both)              8.9 ha  -- essentially the NE patch itself (8.3 ha, p 0.98 / 1.00)
    gained (2025 not 2008)  11.7 ha
    lost (2008 not 2025)     8.9 ha

**Largest gained** (NE-like now, not in 2008): 2.7 ha at E 578812 N 4882388 near Elba
(private, grassland); 2.2 ha at E 579764 N 4885865 and 1.8 ha at E 579866 N 4886054, both
WMA, immediately north of the NE patch; 1.5 ha at E 578824 N 4882649 (private).
**Largest lost** (NE-like in 2008, not now): 1.8 ha at E 579548 N 4883110 and 1.2 ha at
E 579730 N 4883491, both WMA, both WorldCover "tree" in 2021 -- open bright swards in 2008
that are woody now; 0.6 ha at E 579649 N 4883329 (cropland now).

**Reading, as hypothesis.** The 2008 class outside the NE patch sat on the WMA floodplain
south of it and has since gone to trees and shrubs (succession or planting); the 2025
class has spread north of the NE patch along the river on WMA land and appears near Elba.
The NE patch is "kept" by construction -- it is the training reference in both years --
and is the only place where the two different states are known to be the same ground.
The field day can turn these from spectral classes into species: the gained WMA patches
north of the NE patch are 250-450 m from the parking lot at E 579979 N 4885340.

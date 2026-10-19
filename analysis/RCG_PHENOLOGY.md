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

## Q. The lidar fingerprint of the class, per survey (2026-10-19)

Andy: find the patches with and without in each time, then the characteristic lidar
signature from each time -- "a fingerprint of reed canary grass in each LiDAR survey
return set", valuable on its own, then possibly a coarse correction. Scripts
`analysis/rcg/map_ne_like_2021.py`, `analysis/rcg/lidar_fingerprint.py`; table
`rcg/lidar_fingerprint.csv`; gen2 per-cell statistics cached in `rcg/gen2_cellstats.npz`.

**Matching each survey to its own year (Andy's correction, caught in passing):** gen1 (Nov
2008) to the NE-like class on Landsat-5 2008-10-29; gen2 (May 2021) to the class on
Sentinel-2 2021-10-18 (0.1% cloud). In October 2021 the NE patch was in the SAME bright
senesced state as 2008 (1.5-2.0x the floodplain in every band), not the green 2025 state.
NE-like on the floodplain mask: 17.1 ha (2021 S2), 19.2 ha (2021 L8).

**Design.** The class is defined by imagery alone; the lidar is read inside and outside it.
Two passes: (A) all flat floodplain with the NE patch included -- inflated, because the
patch was found by the lidar and trains the class; (B) the independent test: OPEN flat
floodplain only (no gen2 return above 2 m, so sward against sward, not sward against
trees), NE patch EXCLUDED, 12,215 cells, 937 class members for 2008 and 854 for 2021.
Null for scale: F3 cropped cells against the rest. d' = difference of medians over pooled
NMAD. Then a lidar-only LDA predicting the imagery class, 5-fold CV by 500 m block.

**gen1 (Nov 2008), pass B -- the fingerprint of the standing-dead sward:**

    feature (per 5 m cell)              in class        out of class     d'    null d' (F3)
    CSF ground median above gen2 (mm)   +106 (124)      +6 (60)         1.03     0.46
    share of returns at/below gen2      0.13 (0.20)     0.46 (0.39)     1.06     0.53
    p90 of all returns above gen2 (mm)  220 (199)       24 (83)         1.29     0.76
    CSF p90 (mm)                        267 (202)       72 (82)         1.27     0.74
    return intensity (median)           26 (4)          21 (6)          0.95     0.36
    vendor class-2 share                0.64 (0.32)     0.86 (0.21)     0.79     0.60
    gen1 p90-p10 span (mm)              251 (159)       137 (67)        0.93     1.22
    lidar-only LDA, block CV:  AUC 0.64   (pass A, patch included: 0.77)

The sward lifts gen1's ground reading by ~100 mm, lets only 13% of returns reach the gen2
surface, returns brighter, and is classified ground less often by the vendor. Every one of
those separates the class better than cropped fields separate from the rest, except span.
The lidar alone recovers the imagery class only weakly (AUC 0.64): at 16 returns per cell
the per-cell statistics are noisy, which is the same sparsity that limits the correction.

**gen2 (May 2021, full density), pass B -- the fingerprint of the cut sward:**

    feature (per 5 m cell)              in class        out of class     d'    null d' (F3)
    vendor class-2 share                0.89 (0.06)     0.52 (0.11)     4.10     0.10
    height median above gen2 (mm)       0 (0)           10 (15)         0.95     0.00
    share at/below gen2                 0.47 (0.09)     0.40 (0.13)     0.70     0.13
    returns per cell                    227 (7)         218 (19)        0.62     0.35
    near-ground 0.15-2 m share          0.00            0.00            0.00     0.00
    lidar-only LDA, block CV:  AUC 0.68

In May 2021 the class is an exceptionally CLEAN surface: 89% of returns classified ground
against 52% on the rest of the open floodplain and 54% on cropped fields, no near-ground
vegetation at all, every return within a few cm of the surface. That is a short, cut or
flattened sward -- a hay field in spring -- not standing grass. gen2's fingerprint of this
class is therefore "too clean", and gen2's own vegetation offset there is ~0.

**Offsets by class, flat floodplain, NE patch excluded (mm, median (NMAD)):**

    class                      n       raw DoD        corrected DoD   gen1 span   penfail
    NE-like 2008            2,828   -120 (124)        -51 (64)         298        25%
    NE-like 2021            2,604    -81 (114)        -34 (60)         244        22%
    both years                570   -104 (99)         -47 (42)         299        21%
    neither, open          10,703     -3 (57)          +1 (48)         133        12%

Outside the class the open floodplain reads zero. Inside, gen1 reads ~100-120 mm high in
2008, and the logistic correction removes about half of that; the NE patch itself, where
the sward was closed, reads -246 raw and is mostly unreachable (penetration failure).

**What this gives, and what it does not.** There IS a fingerprint in each survey, and
they are different because the sward was in different states: in gen1 a standing-dead
canopy ~100 mm up that few beams penetrate; in gen2 a flattened surface so clean that
near-ground vegetation vanishes. A coarse correction built on the imagery class is now
possible in principle -- the class-level offset (-120 raw, -51 after the logistic) is a
measured number with a spread -- but it would be a class-mean correction with NMAD
~120 mm per cell, and the class itself is a spectral state that a quarter of grassland
shares in some years. The field-day ground truth decides whether the class is reed canary
grass and what the true ground under it is; the fingerprint says what to expect in the
returns when it is.

## R. Mode shift or fixed percentile? Head to head on the class cells (2026-10-19)

Andy: does reed canary grass create a MODE ~100 mm above the ground, or is there a FIXED
PERCENTILE between the returns and the real ground, with real ground returns present --
and which recovers a bare-ground surface? Script `analysis/rcg/mode_vs_percentile.py`.
Returns pooled 3x3 (~150 per cell); residual = estimate minus the gen2 surface (mm), gen2
being the ground reference (it reads a clean cut sward in this class, section Q). Model
constants (c for the mode shift, q for the percentile) are fitted on training blocks and
scored on held-out blocks, 5-fold by 500 m block.

**Diagnostic first.** A rigid shift keeps the distribution's widths; a mat over ground
widens the lower side.

    population                              cells   median  mode   p10   med-p10  p90-med  share<=0
    open floodplain outside both classes   10,651     10     10   -70      80       80      0.48
    NE-like 2008, open, no NE patch           932    130     70   -30     160      220      0.16
    both years, open, no NE patch             251    150     90   -10     160      240      0.14
    NE patch itself (closed canopy)         1,690    250    230    70     180      200      0.02

The lower width DOUBLES in the class (160 vs 80 mm): the distribution does not translate,
it stretches, with the mode at +70 and a tail that still reaches the surface (16% of
returns at or below it). In the closed NE patch nothing reaches the surface (2%).

**Head to head, held-out (bias = median residual; NMAD = per-cell error; slope = residual
on gen1 span, mm/mm -- the vegetation dependence left behind):**

    NE-like 2008, open, no NE patch (932 cells)     bias   NMAD   slope
      raw median (no correction)                     130    119   0.344
      raw mode                                        70     89   0.348
      MODEL 1  median - c  (c ~ 90-130 fitted)        20    119   0.368
      MODEL 1  mode - c    (c ~ 70 fitted)             0     89   0.369
      MODEL 2  fixed p10   (no fit)                  -30     59   0.072
      MODEL 2  fixed p15   (no fit)                  -10     59   0.112
      MODEL 2  q fitted per fold (0.11-0.14)         -10     59   0.109
      minimum                                       -130     59  -0.043
    both years, open, no NE patch (251 cells)
      MODEL 1  mode - c                                0     89   0.119
      MODEL 2  fixed p15 / p20                     10/30     59  -0.07/0.01
    NE patch itself (closed canopy)
      MODEL 1  mode - c (c ~ 250)                      0     89   0.180
      MODEL 2  fixed p10                              70     59   0.024
      MODEL 2  q fitted per fold                     -10     59   -- but q = 0.02-0.03

**Answer.** The FIXED PERCENTILE is the effective one. A low percentile of the pooled
returns, about the 10th to 15th, recovers the gen2 surface with held-out bias within
+/-10 mm, HALF the per-cell error of the mode shift (59 against 89-119 mm), and -- the
decisive part -- it removes the dependence on sward density (slope 0.07-0.11 against
0.37): denser sward lifts the mode further, but the ground returns stay where the ground
is. The mode shift corrects only the class MEAN; every cell keeps its own error and the
densest cells stay the most wrong. Both agree the mode sits ~70-90 mm above the ground in
open sward, so the "mode ~100 mm up" picture is right as a description and wrong as a
correction.

Two limits. (1) The fixed percentile works because ground returns EXIST: in the closed NE
patch the fitted q collapses to 0.02-0.03, which is the minimum by another name, has no
population target, and reads +70 at p10 -- there the answer stays the mask. (2) The
percentile is CLASS-conditional: ~0.12 inside the class against ~0.48 on the open
floodplain outside it (the share at or below the surface), so applying it needs the
imagery class or an equivalent lidar fingerprint to say where. The 59 mm NMAD is a floor
shared by every low-quantile estimator here (20 mm histogram bins, gen2's own noise).

## S. Is there a characteristic return SHAPE -- detectable without imagery? (2026-10-19)

Andy: is there a characteristic shape in the spreading that corresponds to presence /
absence, or to cover intensity, so the grass could be detected from return structure alone
in 2008 and/or 2021, with the imagery maps as the independent truth for learning it?
Scripts `analysis/rcg/return_shape.py`, `analysis/rcg/return_shape_fine.py`; figure
`figures/rcg_return_shape.png`; table `rcg/return_shape.csv`.

**Design.** Shape features measured about the cell's OWN mode, so nothing cross-epoch
enters: lower width (median - p10), upper width (p90 - median), mode - median, Bowley
skew, peak share within +/-50 mm of the mode, upper-tail e-folding length and mass beyond
+100 mm, lower mass beyond -100 mm. gen1: 3x3-pooled CSF returns, 20 mm bins. gen2: the
full-density cloud streamed into per-cell histograms, first at 20 mm then at 5 mm bins
(windows scaled to +/-20 and 50 mm). Truth: the imagery class of the survey's own year
(presence) and two intensity proxies (logit of the class probability; NIR reflectance
relative to the floodplain). Domain: open flat floodplain, NE patch excluded (and
included, for reference). Scoring by 5-fold CV on 500 m blocks.

**gen1 (Nov 2008): YES for presence, as a class; NO for a cover gradient.**

    shape feature (about own mode)     in class        out of class    d'
    upper width p90 - median           220 (89)        80 (30)        2.11
    peak share within +/-50 mm         0.30 (0.12)     0.56 (0.20)    1.52
    upper-tail mass beyond +100 mm     0.26 (0.18)     0.07 (0.09)    1.38
    lower width median - p10           140 (59)        80 (30)        1.28
    upper-tail e-folding length (mm)   126 (75)        53 (49)        1.14
    mode - median                      -20 (30)        0 (0)          0.95
    lower mass beyond -100 mm          0.08 (0.09)     0.03 (0.04)    0.69
    Bowley skew                        0.09            0.00           0.38
    lidar-only LDA, shape features alone:  AUC 0.71 (patch excluded), 0.82 (patch included)

The characteristic shape: the peak collapses to half its share, the distribution widens on
BOTH sides, an exponential upper tail of ~125 mm e-folding carries a quarter of the
returns, and the mode slips 20 mm below the median. Adding intensity and return count
changes the AUC by < 0.01: the shape IS the signal. The mode-centred mean PDFs
(figure) show it as one coherent form, and the in-class PDF deepens monotonically from
p 0.5-0.95 to p > 0.95.

Cover INTENSITY is not recoverable from the shape: block-CV R^2 of logit(p) on the shape
features is 0.03 (patch excluded; 0.25 with it, which is presence again); R^2 of the NIR
index is negative; within the 932 class cells the rank correlations of every shape feature
with NIR, SWIR1 or logit(p) are |rho| <= 0.18. The sward reads as present or absent, not
graded -- or the imagery's gradient is not a cover gradient either.

**gen2 (May 2021): NO shape signal at all.** At 20 mm and at 5 mm bins every shape
feature is identical inside and outside the class (widths 40 mm, peak share 0.54 vs 0.55,
tail mass 0.03 vs 0.03; d' <= 0.35; shape-only AUC 0.18-0.20, i.e. nothing). The
May-2021 return GEOMETRY of a cut sward is the same as the rest of the open floodplain.
What differs is the vendor's classification: class-2 share 0.89 vs 0.52 (d' 4.1), which
is a product of their classifier, not of the return heights, and carries the whole of
gen2's AUC (0.67 with it).

**What this means.** In the leaf-off 2008 survey the standing sward leaves a coherent,
translation-invariant shape in the returns that detects the class without imagery, with
an AUC of 0.71-0.82 on the imagery's own definition of it. In the May-2021 survey the
same ground was a cut sward with no return-structure signature; only the vendor's
classifier noticed something. So a lidar-only detector is a 2008 (standing, senesced
sward) capability, not a 2021 one -- and it gives presence, not cover.

## T. The conjugate test: class in 2021 but not 2008 (2026-10-19)

Andy: look for areas with the class in 2021 but not in 2008 and try the conjugate test.
Script `analysis/rcg/conjugate_test.py` (`CONJ_2021=l8` for the 30 m 2021 class);
figure `figures/rcg_conjugate_groups.png`.

**The prediction.** If gen1's fingerprint and offset belong to the sward's 2008 STATE, the
cells flagged only in 2021 must look like ordinary floodplain in gen1 (Nov 2008), and the
cells flagged only in 2008 must carry the fingerprint anyway.

**The result: the prediction FAILS on the 2021-only side.** Open flat floodplain, NE patch
excluded, 2021 class from Sentinel-2 (10 m) and, like-with-like, from Landsat-8 (30 m):

    group (30 m both years)     n      peak   hi_w   tail_mass   gen1 median vs gen2   share<=0   raw DoD    gen2 class-2
    neither                 10,521    0.57    80     0.06            +10 mm             0.48       -2         0.52
    2008 only                  680    0.34   180     0.23            +90                0.18      -95         0.50
    2021 only                  705    0.35   160     0.20           +130                0.12     -124         0.88
    both years                 252    0.27   260     0.35           +150                0.13     -134         0.87
    d' vs neither, 2021-only:        1.43   1.71     1.26            1.59               1.54      1.28        3.64

The 2021-only cells carry the FULL gen1 fingerprint in November 2008 -- as strong as the
2008-only cells, nearly as strong as 'both' -- and the full offset (-124 mm). The gen1
shape-only LDA trained on 2008-class-vs-neither scores them apart from 'neither' at AUC
0.80 (Sentinel-2 version 0.81). The 2008-only side behaves as predicted (fingerprint
present, AUC 0.78-0.79), so the 2008 class is not spurious; the 2021-only side says it is
INCOMPLETE.

**Why, checked.** Not resolution: the result is the same with both classes at 30 m. Not a
boundary effect: the 2021-only cells' 2008 class probability has median 0.002, 75th
percentile 0.03, and their October-2008 reflectance is ordinary floodplain in every band
(relative NIR 1.00, SWIR1 1.04, red 0.89). In October 2008 they were spectrally NOT the
bright senesced sward -- yet in November 2008 gen1 read a standing canopy ~130 mm up
with one return in eight reaching the surface. And in May 2021 they were as clean-cut as
the NE patch (gen2 class-2 share 0.88).

**Reading (hypothesis).** gen1's fingerprint marks STANDING HERBACEOUS STRUCTURE in
November; the 2008 imagery class marks one spectral state of it -- the bright, senesced,
high-albedo sward -- and misses stands that were still green or otherwise unremarkable on
29 October (the Wisconsin reports' own caveat: "errors of omission were not discovered.
They would be a lot more difficult to identify"). The 2021 imagery flags these cells
because by October 2021 they had been brought into the same hay management as the NE
patch and were in its bright state. So the imagery class is phenology-state-dependent and
under-detects; the lidar fingerprint is more inclusive, and it is a signature of standing
canopy rather than of a species.

**Consequences.** (1) For the lidar: the 2008 fingerprint-and-percentile route (sections
R, S) reaches more sward than the 2008 imagery class does -- the lidar-only detector finds
705 cells (0.4% of the tile, 1.8 ha) the 2008 imagery missed, with the same offset. (2) For
the imagery: a single late-October scene under-detects; the 2021 and 2025 maps show that
the same ground flips state with management. (3) For the field day: whether these cells
are reed canary grass or another standing herb (sedge meadow is the alternative on this
floodplain) is exactly the species question only the ground can answer.

**Scale check (Andy, 2026-10-19; inline run, numbers pasted).** Could the 2021-only result be
the 30 m pixel against 5 m cells? Partly, and measurably. (1) The 2021-only cells are not
mainly fringe: median distance to the nearest 2008-class pixel 55 m; 29% within one pixel
(against 16% for 'neither' cells). (2) The 2008 classifier DOES respond to sub-pixel cover:
among 258 pixels at least half on open floodplain, those with >= 80% of their 5 m cells
carrying the gen1 fingerprint are classified 2008 in 100% of cases, 60-80% in 50%, 20-40% in
22%, < 20% in 4%. The 2021-only cells sit in pixels whose fingerprint fraction has median
0.18 against 0.00 elsewhere: sub-pixel clumps of sward inside pixels that were mostly not
sward in 2008, below what a 30 m mixed pixel can register. So the imagery's omission is
largely a resolution effect -- the class needs most of a pixel to be sward -- and the lidar
at 5 m sees the clumps. Whether those pixels were FLAGGED in 2021 because the sward had
grown to fill them or because the cut-hay state is uniform across a field is not separable
here; both are consistent with the 2025 change map.

# PLAN: reed canary grass from phenology, alongside the lidar work

**Date:** 2026-10-19. **Status:** the focus for now (Andy). **Companion to**
`FRAME_2026-09-23.md` and `FLOODPLAIN_RETURN_PDF.md` (sections J-M), which say why.

## Why this exists, in one paragraph

The floodplain DoD's surviving "erosion" patches are closed grass canopies in November 2008
that gen1 could not penetrate; the species is presumed (reed canary grass, RCG) and the
2008 ground under them is unobservable from the lidar. Two things follow. (1) The
presumption is testable from imagery, because RCG has a phenological signature -- it greens
up earlier and senesces much later than native floodplain vegetation -- which the Wisconsin
DNR used to map it statewide from mid-October Landsat. (2) More important for the lidar
work: a 2008 autumn vegetation map is an INDEPENDENT vegetation measure for the gen1
epoch. Everything in sections G and L failed for want of exactly that -- gen1's own cloud
cannot measure its vegetation without reading the answer. A Landsat index from 29 October
2008 is outside the lidar entirely. It is the tau that Beer-Lambert needs.

## What exists (verified 2026-10-19 by STAC query, Planetary Computer)

    2008, Landsat Collection 2 L2, path 26 row 29, over Elba (cloud cover as catalogued):
      2008-09-19  Landsat 7   2%     2008-10-21  Landsat 7   1%  (SLC-off gaps)
      2008-10-29  Landsat 5   0%     <- two weeks before the gen1 flight; THE scene
      2008-04-04  Landsat 5  10%     2008-05-06  Landsat 5  15%
      2008-05-14  Landsat 7   0%     2008-05-22  Landsat 5   9%
    PHASE 0 DONE 2026-10-19: 23 scenes fetched as NDVI GeoTIFFs + clear masks into
    data/derived/elbaext_delong_star_nosurf/rcg/ (index scenes.json). 2008-10-29 L5: clear
    1.000 over tile and floodplain, snow 0.0000. Summer added: S2 2024-08-03 (0%), 2024-08-08;
    Landsat 2008-08-02/08-10/08-18. 2024-10-27 and 11-01 are partly clouded (clear 0.71, 0.29).
    Sentinel-2 L2A, tile 15TWJ, present day:
      2024-09-27 0%, 2024-10-02 0%, 2024-10-17 1%, 2024-10-22 17%, 2024-10-27 19%, 2024-11-01 23%
      2024-04-05 3%, 2024-04-15 7%, 2024-05-05 1%, 2024-05-10 6%
      2023-10-23 1% (the only clean October 2023 date)

Tooling: `~/projects/mn-3dep-audit/mn_3dep_audit/ndvi.py` already does Sentinel-2 NDVI for
this tile via pystac-client + odc-stac + planetary-computer, with SCL cloud masking and an
ESA WorldCover overlay; its docstring notes the same pattern works on `landsat-c2-l2`.
Ground layers in this repo: the floodplain mask, slope, the penetration-failure flag, gen1
span quartiles, the vegetation-corrected DoD, WMA boundaries, field-site table.

## Principles carried over from the lidar work (each has bitten once)

* A method travels, a number does not: thresholds are proposals, swept, and stated.
* Null-check every metric on populations whose answer is KNOWN before reading it on the
  population in question (forest senesces early October; harvested corn is bare by then).
* Level is gauge: compare Landsat to Sentinel-2 by relative position among classes within a
  scene, never by absolute NDVI across sensors.
* Like with like: aggregate the 5 m DoD to the 30 m Landsat grid before any comparison;
  never regress a 5 m field on a 30 m one.
* The producer labels the number: every map and table prints its scenes, masks, and metric
  definitions.
* Do not train a classifier on the lidar's own patches and then "test" the lidar with it.
  Training truth comes from outside the lidar (Phase 3).

## STATUS 2026-10-19 (see RCG_PHENOLOGY.md): Phases 0, 1, 2, 4 RUN. Phase 4 is NEGATIVE --
October-2008 Landsat greenness does not predict the floodplain DoD (r within +/-0.1), and in
2008 the closed-canopy cells were ordinary floodplain by phenology. 4(a) and 4(c) are closed.
Live: Phase 3 as a species question with truth from outside the lidar, and Phase 4b (valley-
wide survey sites from the 2024 early-green + never-senescing signature). Phase 5 is moot
for the DoD unless a structure-sensitive (not greenness) 2008 layer is found.

## CORRECTION 2026-10-19 (Andy's challenge; RCG_PHENOLOGY.md section N): the Wisconsin method
uses ALL bands on one mid-October scene, SWIR first, no NDVI. Redone that way, the NE patch IS
spectrally distinct in Oct 2008 (d' 2.1-2.6 in red/NIR/SWIR1 vs the floodplain; NDVI 0.15),
brighter in every band at the same NDVI, and k-means isolates it at 10-20x the base rate.
The floodplain-wide DoD is still not predicted (six-band CV R2 0.09): 4(c) stays closed;
Phase 3 (classification with outside truth) and 4b (valley sites) are back on, done the
Wisconsin way: wetland/floodplain mask, all bands, clustering, field truth.

## Phases, in order, each with the result that would stop the next

### Phase 0 -- environment and scene access (half a session)
Install or reuse the STAC stack (the audit repo's venv or this one); load one Sentinel-2
date and one Landsat date for the tile, write both as 10 m / 30 m NDVI GeoTIFFs in
EPSG:26915 with their SCL / QA masks. Check the 2008-10-29 scene's QA for snow and for
the SLC-off gaps on the L7 dates. STOP if the Landsat-5 scene is snow-covered or missing.

### Phase 1 -- the present-day signature, and its null check (one session)
Per 10 m cell, from 2024 (replicate on 2023 where clean): NDVI on each clean date;
metrics -- late-October NDVI, the ratio October/summer (retained greenness), early-April
NDVI (early green-up), and the senescence date (first autumn date below half the summer
maximum). Populations whose state is known, as the null check: forest, cropland (ESA
WorldCover + the 2023 photo), the F3 bare floodplain controls. THEN the populations in
question: the penetration-failure (F2) patches against the sedge-quartile (F1) cells. The
result is a table: metric medians per population with NMAD. If F2 does not separate from
F1 on any metric, the presumption is not supported by present-day phenology and Phase 3's
training must come entirely from outside; Phase 2 still runs because 2008 is the epoch
that matters.

### Phase 2 -- the 2008 state (one session)
The same metrics at 30 m from Landsat 5/7 2008: 29 October against 19 September and the
May scenes. Same null check (forest, crop), then the F2 patches (11-19 pixels each; report
per-patch pixel counts and treat < 10 as unresolved). Compare the F2 patches' October-2008
retained greenness with the sedge cells' and the bare controls'. This is the direct test:
were the surviving patches late-green in the autumn the lidar flew? Also produce the
tile-wide 2008 October-greenness map as a layer in the product directory, labelled.

### Phase 3 -- a classifier with independent truth (one to two sessions)
Truth from outside the lidar: EDDMapS / MN DNR invasive-species occurrence points for
Phalaris in Winona County; Whitewater WMA management documents if they map RCG stands; the
Wisconsin DNR thresholds as a published starting point; and the class field day, which can
GPS stand boundaries while surveying (double duty with the RTK plan). Method: a simple,
explainable rule on the phenology metrics first (as Wisconsin did), a random forest only if
the rule fails, and NEVER more complexity than the truth can score. Score by held-out truth,
report confusion with the F2 patches as a separate, untouched test set. Output: an RCG
likelihood map at 10 m (present) and 30 m (2008).

### Phase 4 -- feed it back to the lidar (one session)
Aggregate the star-only DoD and the gen1 span to the 30 m Landsat grid. Then, on flat
floodplain pixels: DoD against October-2008 greenness -- a regression whose x is
INDEPENDENT of the lidar, the first such on this floodplain. Three uses, in order of
safety: (a) a 2008 RCG-likelihood MASK class for the penetration-failure cells, replacing
"unobserved" with "unobserved, closed canopy, RCG-likely"; (b) a per-class LoD; (c) only if
the regression is clean and Andy chooses it, a per-pixel offset model with an independent
x. The circularity warning of section I does not apply to (c): x is not the DoD.

### Phase 4b -- survey sites from imagery alone, throughout the valley (additional goal, Andy 2026-10-19)
Once the present-day signature is validated (Phase 1) and, if it survives, classified
(Phase 3), run it over the WHOLE Whitewater valley at 10 m -- not only the tile -- and
intersect with the WMA / state-park boundaries and the WMA parking and road layers, to
rank reed-canary-grass stands by size, access and public status as places to survey
through the grass. This needs no lidar at all; it is the imagery product the field class
uses. Output: `rcg_survey_sites_valley.csv` and a map, same form as FIELD_SITES_RTK.md.

### Phase 5 -- what travels (later)
The Wisconsin method ran statewide. A 2008-season Landsat retained-greenness layer over
the whole SE MN 2008 lidar block (same path/row family, same autumn) would flag closed
canopies for every tile before any DoD is built. That is the Minnesota-scale product; it is
not started until Phases 1-4 have produced a result on this tile.

## Deliverables, as named artifacts
* `analysis/rcg/` scripts: `fetch_scenes.py`, `phenology_metrics.py`, `null_check.py`,
  `landsat_2008.py`, `classify.py`, `dod_vs_greenness.py` -- each prints its inputs, scenes,
  masks and metric definitions above its result.
* Layers in `data/derived/elbaext_delong_star_nosurf/`: `s2_oct_ndvi_2024.tif`,
  `s2_retained_greenness_2024.tif`, `landsat_oct2008_ndvi_30m.tif`, `rcg_likelihood_*.tif`.
* `analysis/RCG_PHENOLOGY.md`: the record, pasted numbers only, negatives included.
* Field day addendum to `FIELD_SITES_RTK.md`: stand-boundary GPS protocol.

## Decisions that are Andy's (proposals, not defaults)
1. Metric definitions: "late October" as the single clean date vs a window median; the
   senescence threshold at half the summer maximum.
2. Whether to pool 2023 and 2024 for the present-day signature.
3. The minimum pixel count for a 2008 patch to count as resolved (proposed 10).
4. Whether the classifier's truth includes the field-day GPS boundaries (yes, if the day
   happens before Phase 3 closes).
5. Which of Phase 4 (a)-(c) to adopt; (c) changes the product and is a separate decision.

## Risks, stated once
* Hay cutting: a field cut in late summer regrows and is green in October regardless of
  species; the October signal alone cannot separate "RCG" from "any cut hay field". The
  April green-up metric and the Wisconsin thresholds are the discriminators; if they fail,
  the honest output is "late-green grassland", which is still the class the lidar needs.
* 30 m against 1-2 ha patches: marginal; the per-patch pixel count is reported with every
  number and sub-10-pixel patches are not scored.
* Landsat 7 SLC-off gaps may fall across the patches; the Landsat-5 2008-10-29 scene is
  gap-free and is the primary.
* Cross-sensor NDVI (TM vs MSI) is not comparable in level; only within-scene ranks are.

## Sequence with the lidar work
The lidar side pauses at: the uncapped taper, the slope cut, the terrain-slope test of the
mode, the hillslope application, and the four untracked September files. None blocks this
plan, and Phase 4 is where the two rejoin. The field day (RTK + stand boundaries) serves
both and is scheduled by Andy.

# PLAN: RTK field-verification sites up and down the Whitewater valley (2026-10-20)

A plan, not a run. Nothing here was downloaded, processed or committed; every number is
pasted from a read-only check made while writing (listed under "Checks run"). Andy's intent,
verbatim: "identify a number of locations up and down the Whitewater valley for field
verification of land-surface elevation ... primarily in thick vegetation (RCG = reed canary
grass), but also in a range of veg conditions ... 'Erosion' will help us to highlight these
regions. Likewise, obtain imagery (Sentinel-2 for resolution) to evaluate veg cover."

What this extends: one-tile site selection is done (`analysis/FIELD_SITES_RTK.md`,
`analysis/field_sites_rtk.py`, 20 sites in six classes on WMA land in the Elba tile);
public-land layers from `analysis/fetch_public_land.py`; the imagery program and its limits
in `analysis/RCG_PHENOLOGY.md` (sections N, S, T, U, V, W, X, Y, Z) and `analysis/rcg/`; the
frame `analysis/FRAME_2026-09-23.md` (2026-10-18/19/20 sections); Phase 4b of
`analysis/PLAN_RCG_REMOTE_SENSING.md`, which named the valley-wide imagery product. The
governing constraint is `method-development-not-rollout.md`: one tile is the instrument; a
valley-wide lidar run is production and is chosen deliberately, with its cost stated.

## 1. The AOI

**Definition.** The valley floor of the Whitewater River: Main stem (DNR Kittle M-031), North
Branch (M-031-018), Middle Branch (M-031-019) and South Branch (M-031-017), from the
headwaters near Dover / St. Charles / Elgin to the Mississippi at Weaver. Survey access is
restricted to public land: Whitewater WMA (three units), Whitewater State Park, Carley State
Park (on the North Branch), and any other DNR unit the layers return.

**Boundaries, from the services already used** (`analysis/fetch_public_land.py`, Geospatial
Commons ArcGIS REST, EPSG:26915), queried 2026-10-20 over the box E 552000-592000,
N 4862000-4900000:

* WMA polygons: `us_mn_state_dnr/bdry_dnr_wildlife_mgmt_areas_pub/FeatureServer/0`. Returns
  six units; the three Whitewater units are Main Branch 8767 ha, North Branch 1508 ha,
  South Branch 834 ha (shapely area of the returned geometry). McCarthy Lake WMA (1320 ha)
  is Mississippi floodplain north of the mouth, outside the Whitewater; Haverhill and
  Vermilya are upland, outside.
* State parks: `bdry_dnr_lrs_prk/FeatureServer/2` ("Dnr Management Units Prk", POLYGONS).
  Layer 0 is reference points and returns zero-area features; `fetch_public_land.py` queries
  layer 0 and must be changed to layer 2. Whitewater State Park 826 ha, Carley State Park
  84 ha, John A. Latsch 154 ha (Mississippi bluffs, outside).
* WMA parking and roads: `struc_dnr_wld_mgmt_areas_pub` layers 0 and 1. 53 WMA parking lots
  in the box: Main Branch 40, North Branch 7, South Branch 3.
* River centrelines: `water_dnr_hydrography/FeatureServer/0` ("Dnr Rivers And Streams"),
  field `kittle_name LIKE '%WHITEWATER%'`: 7 features, 177.3 km, bounds E 557491-587128,
  N 4868983-4895859.

**Valley-floor polygon.** For COUNTING, a centreline buffer is the proxy: half-width 300 m
gives 88.8 km2, 500 m 138.7 km2, 800 m 208.5 km2, 1000 m 253.3 km2. For the PRODUCT the valley
floor must come from where it comes from at Elba: the pipeline's `floodplain_mask.npy`
(convexity step, `src/lidar_diff_icp/workflow.py`; valley top from `refcells.py`; the
`whitewater` Site uses the "histogram" mode). One elevation cut does not travel up a valley
whose floor rises from Weaver to Dover, so the per-tile histogram or `landscape_of` is the
mechanism, with a HAND-type height-above-channel mask from the gen2 DTM as the fallback
where the histogram has no minimum (the Cook case in `sites.py`). Andy's call (D1). Public
land inside the 800 m proxy: 6562 ha (Main Branch 3961, North Branch 1330, South Branch 633,
Whitewater SP 675, Carley 84).

**Tiles.** 2008 DNR tiles are ~2497 x 3471 m (`lidar-2008-data-gotchas.md`), named
AAAA-BB-CC. Programmatic listing: `tiles.centroid_index()` (`src/lidar_diff_icp/tiles.py`,
cached `data/mn_tile_centroids.csv`, 55,296 tiles) gives every tile centroid; keep those
within a half-diagonal (2138 m) of the valley polygon. Counts against the proxy: 66 tiles at
300 m, **71 at 500 m** (listed in "Checks run"), 79 at 800 m, 82 at 1000 m. Tiles whose cell
intersects PUBLIC land inside the proxy: **29** (4342-27-64, 28-64, 29-59..64, 30-60..64,
31-62, 31-63; 4358-26-01..03, 27-01, 27-02, 28-01, 29-01, 29-02, 30-01, 30-02, 31-01, 31-02,
32-01, 32-02). Already in `data/before/`: 19 of the 71, including 14 of the 29.

**3DEP (gen2).** `threedep.find_projects` at five points: the mouth, Elba, the South Branch at
St. Charles and the North Branch at Elgin resolve to **MN_SEDriftless_2_2021**; the Middle
Branch headwaters at Dover resolve to **MN_SEDriftless_4_2021**. From
`data/3dep_boundaries.geojson`, project 4's part of the box lies west of E ~574686, so the
upper North and Middle branches need a second EPT project, and the two must be treated as
two gen2 acquisitions (vintage check per `threedep.resolve_reference`, which refuses to
downgrade silently).

## 2. Data acquisition

**gen1 (2008 SE MN DNR).** Source verified 2026-10-20:
`https://resources.gisdata.mn.gov/pub/data/elevation/lidar/county/<county>/laz/<AAAA-BB-CC>.laz`,
with `<county>_tile_list.txt` beside it; HEAD on `winona/laz/4342-28-62.laz` returned 200,
23,379,726 bytes, `accept-ranges: bytes`. The project metadata page
(`winona/lidar_semn2008.html`) says nine counties including Olmsted, Wabasha and Winona,
acquired 18-24 November 2008, AeroMetric, 2400 m AMT, NPS 1.0 m. Of the 71 proxy tiles, the
county listings place 27 in winona only, 23 in olmsted only, 4 in wabasha only, and 17 in
two or three listings (a tile on a county line is served from more than one directory;
none is missing from all three). Per tile 7.7-33 MB on disk (`data/before/`, median ~22 MB):
71 tiles ~1.6 GB, 29 tiles ~0.65 GB, 52 to download, via `scripts/fetch_tile.py --tile NAME
--county C` (retry-aware). Gotchas (`lidar-2008-data-gotchas.md`): no CRS in the files
(assign EPSG:26915); PDAL cannot read the old LAZ (laspy transcodes); GEOID03, one survey
(lidar_semn2008) for the whole valley. Merge multi-tile AOIs with
`scripts/merge_gen1_tiles.py`, which keeps every dimension.

**gen2 (3DEP).** Two routes. (a) `scripts/fetch_3dep_curl.py --auto --bounds ...` walks the
EPT hierarchy and downloads every node at every depth (no depth cap since the whitewater
seam incident), 16 workers, streaming to one LAS. Measured volumes: whitewater tile
148,050,625 points, 814 MB full-density LAZ (`data/derived/whitewater/data_completeness.json`,
ratio 0.971); elbaext (18 km2) 415 M points, 2.3 GB full and 665 MB class-2 only; node cache
`data/after/_ept_tiles/` 13,441 files, 3.3 GB. Per tile: ~150 M points, ~0.8 GB full,
~0.25 GB class 2, ~1.5 GB of deletable node cache. Wall time per tile is recorded nowhere I
read (UNVERIFIED; time the first tile). (b) Indexed box reads (`gen2-ept-indexed-reads.md`):
~1 s per box, flat in box size to 100 m, classification intact; for point checks at sites
and marks, not for a DoD. The elbaext base product was built against the CLASS-2 gen2 cloud
(`data/derived/elbaext/meta.json`: `after = elbaext_3dep_fd_class2.laz`); full density served
only PyForestScan, canopy_struct and the gen2 histograms, none needed here (section 3). So
fetch full, write class 2, keep full only if Andy wants gen2 histograms later (D2).

**Budgets on the shared laptop** (`machine-memory-gentle-processing.md`). At the time of
writing `free -h` showed 35 of 46 GB RAM in use and swap FULL (8.0 GB used, 48 KB free):
nothing heavy should start until that clears. Disk: 191 GB free. Per tile: gen1 22 MB, gen2
class-2 0.25 GB (full 0.8), CSF cache 0.21-0.25 GB (`data/csf_cache/whitewater.las` 211 MB),
derived products 0.37 GB with the beam table (`data/derived/whitewater`) and 0.67 GB at
elbaext of which `beam_offset_table.parquet` is 541 MB. Whole valley, 71 tiles, keeping
class-2 gen2 and the beam tables: ~1.6 + 18 + 16 + 26 = ~60 GB; with full-density gen2 kept,
~100 GB. The 29-tile public-land sample: ~25 GB (~40 GB with full density). Cache: the CSF
ground per tile (`cache-expensive-pipeline-stages.md`: ~460 s, deterministic, never re-run),
the EPT boundary index, the tile centroid CSV (already cached), the class-2 gen2 per tile.
Delete after use: the EPT node cache, the full-density gen2 if D2 says class 2 suffices.

## 3. Processing, the minimum per tile

From the step graph in `src/lidar_diff_icp/workflow.py` (`lidar-diff-workflow --tile ...
--plan` prints the order) and the Elba build (`scripts/build_elbaext_routes.py --route
delong --swath-frame star --no-correction-surface`, the adopted product per
`star-only-adopted-floodplain-is-vegetation.md`):

1. `completeness` (EPT hierarchy count, no point download) BEFORE anything else: the
   whitewater tile shipped for weeks at 5.52 returns/m2 east of a seam against 15.45 west.
2. `base` = `pipeline.difference_dem` with the Elba `SHARED` kwargs (res 5 m, ground_q 0.50,
   slope_normal, CSF, robust_stable, swath_tie intercept), `swath_frame="star"`,
   `correction_surface=False`, valley top per D1. Produces `dod.npy`, `lod.npy`,
   `z_after.npy`, `stable.npy`, `corrections.json`. Cost: CSF ~460 s per tile (recorded);
   the rest of difference_dem is not timed in anything I read (UNVERIFIED; the first tile
   sets it).
3. `slope`, `ridge_mask`, `convexity` (produces `floodplain_mask.npy`), `curvature`
   (`curv_laplacian.npy`, required by the beam table). Minutes each.
4. `gen1_angles` and `beam_table` (per-return gen1 offsets against the gen2 surface with the
   star registration applied). Then, from the beam table: `vegetation_lod.gen1_spread`
   (p90-p10 per cell) and `vegetation_lod.penetration_failure` (NO gen1 return at or below
   gen2's surface, min 10 returns; a yes/no test, no threshold), as
   `scripts/figures_with_vegetation_lod.py` writes them; and the translation-invariant
   shape features of `analysis/rcg/return_shape.py::shape_feats` on the 3x3-pooled gen1
   histograms (peak share, upper width, upper-tail mass and e-folding length, mode minus
   median). The LDA fingerprint detector trained at Elba (AUC 0.71-0.82, section S) is
   applied, not refitted.
5. OPTIONAL, cheap, Andy's call (D3): `analysis/floodplain_vegcorr_apply.py`, the Elba
   logistic pct0(span) correction, to fill the "corrected DoD" column of the site table. Its
   four parameters were fitted at Elba; applying them elsewhere is a transfer, and the table
   should say so.

**Skipped, with the reason:** DeLong surface (imprints); `pfs_cover` (opt-in, needs an
untwine COPC and ~14 GB scratch per 1 GB cloud); `nearground*`, `class2_spread`, `q2_fit`,
`dod_cover`, `lod_cover`, `canopy_struct` (superseded cover route, gen2-side); boresight
(retracted); gen2 return histograms (section S: no shape signal); the half-sample-mode
rebuild (analysis scripts only, not in the pipeline).

**Gentle scheduling.** One tile at a time, no overlapping point-cloud jobs; `free -h` before
each tile, stop if swap fills; an overnight serial loop with a per-tile log that skips
tiles whose products exist (the workflow's STALE/OK state); CSF to
`data/csf_cache/<tile>.las` so a crash resumes without re-classifying; delete the EPT node
cache per tile once the LAS is written; the 29 public-land tiles first (section 6).

## 4. Screening for candidate sites

**Lidar layers, each a 5 m grid per tile** (all exist at Elba; `field_sites_rtk.py` is the
template): raw star-only DoD; `floodplain_mask` AND terrain slope <= 2 deg (flat
floodplain); `gen1_spread_mm` (pooled 3x3 as the correction used it); `penetration_failure`;
the gen1 return-shape fingerprint score; LoD. "Erosion" patches (DoD < -LoD, flat
floodplain) are the highlighter Andy asked for; at Elba they are vegetation by spread
quartile (+10.5 / -6.1 / -46.3 / -120.1 mm) and water and sand are excluded first.

**Sentinel-2 program** (Planetary Computer STAC, `analysis/rcg/fetch_scenes.py` generalised
to the valley bbox). The valley straddles MGRS tiles 15TWJ and 15TWK, so every date needs
both. STAC count 2026-10-20, cloud < 40%: late Oct-Nov 2025 twelve scenes, with 2025-10-02,
2025-10-27 and 2025-11-16 at 0% on both tiles; spring 2026 sixteen, with 2026-04-10,
2026-04-15 and 2026-05-15 clear; autumn 2026 seven, with 2026-10-02, 10-07 and 10-17 clear
on both. The three uses, each a different season (sections N, W, X):

* PRESENCE: one mid/late-October scene per year, all ten bands, the Wisconsin way (no NDVI;
  SWIR decisive). 2026-10-17 is the current state; 2025-10-27, 2021-10-18 and Landsat-5
  2008-10-29 (whole-valley coverage of WRS-2 026/029 UNVERIFIED) give the history.
* DENSITY: April-May green-up vigour (section W, rho ~0.5 with gen1's November canopy
  signal, untested for block structure): 2026-04-10/04-15/05-15.
* CUT HISTORY: within-season NDVI drops (section V: a cut shows the next clear day, NDVI
  halves, regrowth in a month); a per-field Jul-Oct 2026 timeline says which stands are
  hayed, hence whether thatch will be present on the survey day.

**Classifier.** `analysis/rcg/map_ne_like.py`: two-class maximum likelihood on all bands,
trained on the NE patch's pixels of the SAME year, domain the floodplain mask, WorldCover
confusion check printed. Confusions kept in the table, not hidden: the bright October state
includes hayed fields regrown since August (V); grazed pasture and cut hay (N); a shrub or
standing water breaks a pixel; the 2008 class under-detects sub-pixel sward at 30 m (T).
The class is "NE-like", not a species; the field day supplies the species.

**Combining into a RANGE of conditions.** Rank within class, not across, so thick RCG does
not crowd out the controls. Proposed classes (MINE; all cut-offs are proposals, D4):

| class | lidar criterion | imagery criterion |
|---|---|---|
| R1 closed RCG sward | flat, penetration_failure, DoD < -LoD | NE-like in 2026-10 and in 2025-10; bright or green state |
| R2 open sward | flat, fingerprint score high, NOT penetration_failure, DoD < -LoD | NE-like in at least one of 2025/2026 |
| S sedge meadow | flat, top spread quartile, fingerprint low, DoD < -LoD | NOT NE-like; WorldCover herbaceous wetland or NWI emergent |
| C cropped / bare control | flat, bottom spread quartile, abs(DoD) < 30 mm | low October NDVI; bare in July-Aug 2026 |
| U unfarmed bare | gravel lots and roads (section Z: DoD +2.6 / +6.8 mm) | n/a |
| E forest edge / shrub | flat, within 25 m of forest_pfs or gen2 above-2 m returns | tree or shrub in WorldCover |
| H2 / H3 hillslopes | as in `field_sites_rtk.py` (optional, Andy) | n/a |

Each candidate carries the per-year class membership 2008 / 2021 / 2025 / 2026 and the 2026
cut date, so a site can be chosen as "sward in both epochs" (the cross-epoch calibration
population, section Y) or "sward in one epoch only" (the conjugate population, section T).

## 5. Site deliverable

Same table as `analysis/FIELD_SITES_RTK.md`, one row per cluster (>= 40 cells = 0.1 ha,
MINE), columns: site, class, cells, ha, easting, northing (EPSG:26915), lat, lon (WGS84),
raw DoD (median, mm), corrected DoD (if D3), gen1 span, slope, canopy cover, penetration
failure %, fingerprint score, NE-like 2008/2021/2025/2026, 2026 cut date, distance to the
nearest WMA parking lot and road, public unit name, tile. Written as
`data/derived/field_sites_valley.csv` with a per-unit map; the producer prints every input
layer, criterion and count above its result (`trust/provenance.py` rule).

**RTK protocol, as in FIELD_SITES_RTK.md:** at each floodplain point, RTK on the GROUND
(probe through the thatch) AND on the thatch / sward top; late autumn; stand height and
stem density at the point (section W's field quantity); a photo and species call per site.
Datum: MnCORS network RTK gives NAD83(2011) ellipsoid heights; convert with GEOID18 and
compare to gen2; the DoD carries gen1's GEOID03. Tie: occupy a 2021 project mark. Of 519 in
`data/derived/control_marks_gen2_2021_control_LCP-NVA-VVA_marks.csv`, 47 lie in the valley
box and 5 on Whitewater public land: 1078_2021_MN (E 579942 N 4888142, LCP), 1079_2021_MN
(E 579660 N 4885857, LCP, 250 m from Elba's F2-1), 1080_2021_MN (E 577787 N 4881179, LCP),
3056_2021_MN (E 570474 N 4884246, VVA, North Branch), 3089_2021_MN (E 582881 N 4894754, VVA,
lower main stem). Precision 20-30 mm per point against 100-270 mm effects: a few dozen
points per site.

**Additions for the valley (proposals, D5):** 3-4 sites per class per reach (upper, middle,
lower main stem; one each on the North and South branches), spaced > 1 km so no two share a
2008 flight line's offset (`gen1-level-is-per-line.md`: +/- 18 mm per line); one day per
parking cluster, each with a C or U control and a mark occupation; the Elba day is the
template.

## 6. Sequencing, with costs and what stops the next phase

* **Phase A, imagery, whole valley (no lidar; one session).** Extend `fetch_scenes.py`
  and `map_ne_like.py` to the valley bbox and both MGRS tiles; per-year NE-like maps 2008
  (Landsat), 2021, 2025, 2026 at 10 m; spring vigour; 2026 cut timeline; intersect with
  public land and parking. Produces the imagery-only ranking Phase 4b asked for. Stops
  Phase B if the class is not specific valley-wide (the WorldCover confusion check reports
  cropland or grassland mostly NE-like).
* **Phase B, lidar on the 29 public-land tiles (overnight runs).** 14 gen1 tiles on disk,
  15 to fetch (~0.35 GB); 29 gen2 fetches (~25 GB, class 2 kept); 29 CSF runs (~3.7 h of
  CSF alone); then the base and the four cheap steps. Produces the section-5 table for every
  public parcel: the field deliverable. Stops Phase C if the star frame fails on the wide
  lower valley (`robust-stable-wide-floodplain.md`: 37% stable-mask contamination on a
  >600 m valley; the Weaver reach is wider) or the two 3DEP projects disagree at their seam.
* **Phase C, the remaining 42 tiles (D6).** ~35 GB, ~5.5 h of CSF plus the untimed base
  runs; gives private-land context and the valley-long erosion map but no new survey sites.
  Under `method-development-not-rollout.md` this is production, not the next action.
* **Phase D, field days**, scheduled by Andy; the data return to sections R/S/T as the
  species truth and the thatch thickness nobody measured.

**Decisions that are Andy's (none defaulted):** D1 the valley-floor mechanism (per-tile
histogram valley top, landscape, or HAND) and the slope <= 2 deg cut; D2 keep full-density
gen2 or class 2 only; D3 apply the Elba logistic correction off-tile; D4 every class cut-off
in the section-4 table (spread quartile edges 140/209/307 mm are Elba's; the 30 mm control
window; the 40-cell minimum; NE-like p > 0.5); D5 sites per class, spacing, days; D6 Phase
C at all, and whether the 500 m or 800 m proxy defines "the valley" for tile selection
(71 vs 79 tiles); D7 whether to re-time the pipeline on the first tile before committing
to overnight loops.

## 7. Risks

* **Private land.** 6562 ha of the proxy's ~21,000 ha is public; most R-class stands will be
  private. The table carries public status; survey only inside the DNR polygons, which DNR
  labels "for general reference only".
* **Water and wetland cells.** Penetration failure also fires on open water and saturated
  ground; exclude water (NWI, `water_dnr_hydrography` layer 1, gen2 class 9) before ranking,
  as Elba did. The Weaver reach is Mississippi backwater.
* **Tile edges and seams.** Flag sites within 100 m of a tile edge; record the completeness
  ratio per tile, never judge it; the SEDriftless 2 / 4 boundary crosses the upper North
  and Middle branches, so project 4's vintage and seam are checked before any tile west of
  E ~574,700 is differenced.
* **30 m vs 10 m (T).** The 2008 class needs most of a 30 m pixel to be sward; Sentinel-2
  at 10 m has no 2008 counterpart. Sites chosen on the 2008 class favour large fields; the
  5 m fingerprint layer finds the clumps.
* **Management change 2008 -> 2021 -> now (U, V).** A 2008 sward may be hay or crop today;
  the survey measures the present-day analogue and the climatology (X) places the day on
  the cycle. A cut year yields the ground but not the thatch offset.
* **The shared laptop.** Swap was full at the time of writing; the loops are serial,
  per-tile, resumable, and never start while swap is in use.

## Checks run while writing (read-only; nothing processed)

The three county tile lists (winona 223, olmsted 229, wabasha 194 lines); HEAD of one tile
URL; the semn2008 metadata page; the Geospatial Commons queries of section 1;
`threedep.find_projects` at five points; `data/3dep_boundaries.geojson` against the box;
tile-centroid buffer counts; a STAC count on sentinel-2-l2a for three windows; `ls`/`du` of
`data/before/`, `data/after/`, `data/csf_cache/`, `data/derived/whitewater`,
`data/derived/elbaext_delong_star_nosurf`; `df -h`, `free -h`; the marks CSV. The 71 tiles
at 500 m: 4342-27-64, 28-56..59, 28-64, 29-55..64, 30-55, 30-56, 30-60..64, 31-55..57,
31-61..64, 32-56, 32-57, 32-59..64; 4358-26-01..04, 27-01, 27-02, 28-01..03, 29-01..03,
30-01, 30-02, 31-01, 31-02, 32-01, 32-02; 5142-01-56..64, 02-58..62; 5158-01-01.
UNVERIFIED in the text: the gen2 fetch wall time, the base-run wall time beyond CSF, and
the Landsat 2008-10-29 scene's coverage of the whole valley.

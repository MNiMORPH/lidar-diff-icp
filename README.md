# lidar-diff-icp

A reusable workflow to build **bare-earth DEMs of Difference** between Minnesota's
**First-Generation** statewide lidar (`gen1`, flown 2008-2012) and the modern
**Second-Generation** USGS **3DEP** survey (`gen2`, 2020s), correcting the
first-generation survey's navigation error so that what remains is real geomorphic
change. ("First/Second Generation" is MnGeo's own terminology for the two eras.)

The Elba pilot (below) is the test case; the method is written to generalize to
any first-generation MN tile with 3DEP coverage.

## The problem

The first-generation lidar carries navigation (GNSS/IMU trajectory) error that
differs from one flight line to the next — adjacent swaths over a single tile were
flown up to ~1.75 h apart, so each pass has its own, largely independent error state. A single
rigid alignment of the whole cloud is therefore wrong: the error is piecewise per
swath, and its dominant remaining component is a smooth **along-track drift** that
follows the flight path. Two facts make this tractable:

1. **The overlaps are self-calibrating.** Adjacent swaths were flown minutes to
   hours apart, so there is *no real land-surface change* between them — every bit
   of elevation discrepancy in an overlap is acquisition error.
2. **The 2008 data retain per-pulse attribution.** `point_source_id` (flight
   line), `gps_time`, and scan angle all survive — exactly what a per-swath,
   along-track correction needs.

We have no raw trajectory for the delivered clouds (a constraint DeLong et al.,
2022, also faced), so the correction is **data-driven against 3DEP** rather than
system-driven from the navigation logs.

## The workflow, and why each step is what it is

`lidar_diff_icp.pipeline.difference_dem` runs the whole thing. Every choice below
was earned by a measured comparison on the pilot; the non-obvious ones are the
point.

1. **Bare earth: CSF ground classification by default; last-return the fast
   alternative.** A physically based ground filter — PDAL's Cloth Simulation Filter,
   tuned for sparse steep/wooded terrain (`rigidness=1, resolution=1.0,
   threshold=1.5, hdiff=0.5`; ~96% cell / ~94% steep-cell coverage while removing
   ~12% as non-ground) — gives the cleanest, most general bare-earth, removing
   structures and forest understory the heuristic keeps. Opt out with
   `ground_source="last_return"`: last return, `return_number == number_of_returns`,
   *including single returns* (singles dominate flat open ground; dropping them, a
   common `filters.returns groups=last` mistake, empties the agricultural fields).
   On this pilot the two give a near-identical DoD, so `last_return` is the right
   choice to skip CSF's per-tile cost; CSF earns its keep on forest/structure.
2. **Ground = the MEDIAN per cell, taken NORMAL to the local slope.** This is the
   single most important choice. Once the cloth (CSF) has classified ground, the
   returns scatter symmetrically about the true surface, so the unbiased estimate is
   the central tendency — the median (robust to any residual high outlier). A *low*
   percentile instead sits ~1.28σ *below* the surface by an amount set by the cell's
   roughness; because 2021 is ~14× denser and smoother than 2008 that offset differs
   between the epochs and does **not** cancel — it becomes **coherent false change**,
   physically-impossible ridgetop "deposition" that grows with slope. *(History: the
   original heuristic took a low percentile (10th) on RAW last-return points, where
   the true ground sits at the* bottom *of the return distribution and a low pick
   rejects canopy — the right call before classification, dropping ~16–32% of convex
   hillslopes reading as falsely depositional to ~4%. Kept on top of CSF it
   double-counts the cloth — that is the slope-correlated bias above — and the median
   removes it.)* Either way the pick is taken relative to a shared smoothed surface
   (both epochs, `ground="slope_normal"`), which the difference cancels: this also
   removes the downhill bias a *horizontal* pick has on a slope (~35% lower stable σ,
   validated to preserve a known change). `ground="low_q"` is the older horizontal
   pick. Coherent bias, not incoherent noise, fools change detection.
3. **Correct the 2008 points in the acquisition frame, per point, *before*
   gridding** (not post-hoc on the difference) — one **instrumental** term first,
   then the **empirical** ones:
   1. **scanner boresight roll — instrumental, so applied first (opt-in, off by
      default).** A residual scanner-to-IMU mounting-angle error puts an elevation
      bias *proportional to scan angle* into every flight line identically (a sensor
      constant), unlike the per-swath offsets below. With no raw trajectory from the
      vendor, it is self-calibrated from 2008 flight-line **overlap** — where the
      between-line offset difference cancels terrain, so its slope against the
      between-line scan-angle difference *is* the roll — then removed per point as
      `z −= b·scan_angle`. On the pilot the overlap fit gives `b = +2.19 mm/deg`,
      but the ±0.7 spread between flight-line pairs is too wide to call a residual
      resolved, so nothing is applied by default and `boresight_roll_mm_per_deg`
      stays `None` in the delivered corrections. Referencing 2008 against *itself*
      decouples it from the 3DEP lateral tie below, so the chain needs no iteration.
      It removes the cross-track scan-angle asymmetry cleanly (the within-cell tilt
      drops +2.2 → −0.1 mm/deg),
      but its tile-wide DoD footprint is small — much of it self-cancels in the
      per-cell median where swaths overlap, and the per-line-mean part is already
      absorbed by the alignment below — so it is a correctness/consistency fix (a
      tilt attributed to a tilt, reusable across a lift), *not* a scatter reduction.
      Enable with `correct_boresight=True`;
   2. per-swath internal alignment (translation, lowest swath pinned);
   3. lateral **Nuth–Kääb (x,y) registration** to 3DEP — get the horizontal right
      before touching z — then the **deterministic geoid-model vertical offset**
      (`N_gen1 − N_gen2`, GEOID03 → GEOID18), auto-computed per tile from the PROJ
      geoid grids as a constant plus the model's small tilt. Nothing is *fitted* —
      no pad constant, no plane on "stable" surfaces — so the datum cannot absorb
      real hillslope change (the removed reference-plane and order-2 parabola ties
      could; git history keeps them);
   4. **per-swath along-track GNSS-drift spline `f(gps_time)`** — the deterministic,
      physical form of the residual error, and the reusable core: the same failure
      mode statewide, only the coefficients differ per tile.

   The residual warp and real localized change share the same ~100–400 m scale, so
   no data-driven interpolator on the elevation residual can separate them — only
   the acquisition geometry can, which is what the drift uses. A DeLong 400 m
   correction surface (data-driven IDW, TPI floodplain buffer) is available for
   legacy data lacking `gps_time`, but it absorbs localized flat change up to its
   threshold and adds only ~4 mm here, so it is **off by default**.
4. **Difference:** gridded median ground, DoD = 3DEP − 2008 (positive =
   deposition). Cell size (default 5 m) is set by the sparse 2008 density
   (~0.8 pts/m² → ~20 points per 5 m cell for a stable median).

Convention, held everywhere: **DoD is `after − before`; red = erosion, blue =
deposition; standard NW (315°/45°) hillshade.**

## Uncertainty

- **Stable-ground 1σ** (empirical NMAD on low-slope, non-floodplain ground) is the
  trustworthy number — ~0.09 m on the pilot.
- The **per-cell LoD** (`lod.tif`) is a calibrated heteroscedastic error model
  inherited from **xdem** (Hugonnet et al., 2022): the stable-ground DoD
  dispersion is modeled as a function of slope, curvature, and the **ground-estimate
  standard error** — `sqrt(Σ_epoch roughness²/density)`, which combines the two
  distinct within-cell signals: **detrended roughness** (the surface's real internal
  variability, slope removed) as the numerator, and **ground-return density** (how
  much data supports the estimate) as the denominator. These are separate, both
  significant, factors (Aguilar et al. 2005; the Wheaton et al. 2010 covariate set).
  σ is then predicted everywhere, so it honestly rises with slope, roughness, and
  sparse data (~0.04 m flat → ~0.20 m steep on the pilot) rather than being
  relief-inflated. Needs `xdem` (`pip install .[uncertainty]`; its import requires
  `PROJ_DATA` unset); falls back to a within-cell-spread proxy otherwise. The
  slope-dependence is real uncertainty — modeled, not detrended.

## The absolute vertical datum

Swath alignment makes the flight lines mutually consistent. It does not tell you where
the resulting surface sits. `coreg.align_swaths` solves a free network and subtracts the
reference swath's value afterwards, so the gauge touches no swath-to-swath difference –
but the mosaic inherits **the reference line's own vertical error** as its absolute level.
Measured on elbaext, the six per-swath `dz` span

    133  +0.00   134 +22.00   135  +6.20   136  -9.80   137 -18.40   138 -22.60

so re-gauging on a different line moves every elevation by up to **44.60 mm**. An
uncorrected elevation is therefore an arbitrary implementation detail, not a measurement.

**Ground control supplies the one number the network is blind to.** Each epoch is tied to
its own contemporaneous control – the 2008 MnGeo/MnDNR validation checkpoints for gen1,
the 2021 USGS held-out NVA/VVA checkpoints for gen2 – and the correction is applied to
both, so the DoD moves by the difference. Applying it removes the gauge dependence
exactly: with `corrected = z + c` and `c` measured against the same gauged product,
re-gauging by `d` shifts `z` by `+d` and `c` by `-d`, and they cancel.
`ground_control/tests/test_apply_datum.py` demonstrates this rather than asserting it –
uncorrected spread 44.60 mm across the six gauges, corrected spread below 1e-9.

### The relation that governs it

    DoD = c1 - c2 - g

where `c1` and `c2` are each epoch's constant against its own control and `g` is the
geoid-model term the pipeline adds to gen1. **The geoid does not cancel between the two
constants**, which is the trap in this problem: NAVD88 is the datum, whereas GEOID03 and
GEOID18 are models for converting GPS ellipsoidal heights to orthometric ones. Both
control sets publish NAVD88 and are directly comparable, but each epoch's lidar `z` was
converted with a different model, so the two constants reference surfaces in different
frames.

Two independent checks close on this. The measured DoD on stable open ground predicts
`-2.12 mm` against `-2.12 mm` observed over 116,507 cells. Furthermore the control's own
epoch separation, `c1 - c2 = +69.30 mm`, recovers the PROJ geoid difference of
`+67.38 mm` to **1.92 mm** – two survey networks reproducing a geoid model that neither
knows anything about.

### Elba, measured

| quantity | value |
|---|---|
| gen1, delivered surface | **+62.74 ± 23.38 mm** (open ground, 8 marks on 5 lines) |
| bridge, delivered → our reconstruction | **-4.04 ± 11.12 mm** (29 open marks) |
| gen1, our surface | **+58.70 ± 25.89 mm** |
| gen2, its own held-out control | **-2.37 ± 2.37 mm** project-wide; **-6.83 ± 2.96** in the QL1 block |
| geoid term added to gen1 | **+67.38 mm** |
| **DoD shift** | **+2.18 mm**, and it puts stable open ground at **-0.003 mm** |

The DoD shift is small, which is **not** a reason to skip it: the gauge choice it removes
is 21× larger, and its smallness here is a property of line 133 having been a lucky pin.

### Four rules that changed the answer

1. **Epoch-matched control.** 2008 for gen1, 2021 for gen2, never crossed. A 2021 mark on
   a 2008 surface carries thirteen years of real ground change.
2. **Open ground only.** Pooling cover classes bakes canopy response into the datum and
   pre-decides the canopy-versus-erosion question. At Elba this moved the answer 17.17 mm
   and collapsed a disambiguation sensitivity from 8.69 mm to 1.90 mm.
3. **The flight line is the unit of replication.** Marks under one line share that line's
   unknown constant; treating them as independent understates the standard error by a
   measured design effect of 1.40×.
4. **The returns assign the line, never the geometry.** `point_source_id` is reused across
   missions, and a near-north–south line drifts about 1.1 km in easting over 94 km of
   track, so across-track separation is no evidence of a second line. Passes are merged by
   collinearity scaled by the extrapolation's own prediction standard deviation.

### Measuring it at a new site

```bash
# 1. flight-line tracks, once per acquisition (46 tiles, ~60 s, committed thereafter)
./lidar-icp/bin/python ground_control/run_derive_tracks.py \
    --tiles 'data/before/*.laz' --exclude-substring merged \
    --out ground_control/data/gen1_line_tracks.json --chunk-size 2000000

# 2. both epochs' constants at the site
env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python ground_control/run_site_datum.py \
    --easting 578762.8 --northing 4884487.6 --site elbaext \
    --corrections data/derived/elbaext/corrections_geoid.json \
    --tracks ground_control/data/gen1_line_tracks.json \
    --psids 133 134 135 136 137 138 --covers L1O --collinear-sigma 3 \
    --tiles data/before --res 5.0 --gen2-surface ql1_laz \
    --bridge-mm -4.04 --bridge-source "products/bridge_wide_L1O.json" \
    --max-lags-m 20000 40000 80000 160000 --n-lags 25 --n-pairs 800000 \
    --estimators dowd matheron --seed 0 --out SITE_DATUM_elbaext.json

# 3. pass its absolute_datum block straight to the pipeline
#    difference_dem(..., absolute_datum=json.load(open("SITE_DATUM_elbaext.json"))["absolute_datum"])
```

`difference_dem` checks that the constant's `gauge_ref` matches the run's own
`swath_gauge_ref` and raises on a mismatch, because a constant measured against one
reference line belongs to that product and would silently mis-level another. Nothing is
defaulted: `covers`, `gen2_surface`, `collinear_sigma` and the bridge are all required,
and each of them moved the Elba answer by more than the correction itself.

### What is still open

- **Near versus far marks.** Six of Elba's eight open marks sit 14–63 km away and disagree
  with the two near ones by 59.13 mm, which is enough to change the correction's sign.
- **Mechanism, not relation.** The relation is verified; the story that "most of the
  difference was the geoid" is consistent but unproven, and it competes with the
  **unpublished vendor bias adjustments** that both epochs carry and neither publishes.
- **gen2's bridge** is bounded to 0 ± 26 mm by the closure but was never measured
  directly: its checkpoints sit on engineered ground, giving radius spreads of 131–715 mm.
- **The statewide per-line correction** is where control actually pays off. The weighted
  uncertainty falls from 22.75 mm toward 11 mm only as per-line constants improve, which
  needs many marks per line rather than more marks at one site.

## Quick start

```bash
python3 -m venv --system-site-packages lidar-icp   # over apt geospatial libs
pip install -e .
# CSF ground classification (the default ground_source) needs PDAL with filters.csf,
# e.g. `conda create -n lidar-icp pdal` — found automatically on PATH or in conda
# envs, or pass --pdal/csf_pdal. Skip it with --no-csf for a PDAL-free last-return run.

# 1. fetch 3DEP over the tile's bbox (EPT via curl; readers.ept is unreliable here).
#    --auto resolves the covering gen2 project from the bbox and refuses to run
#    unless its boundary fully covers the tile (--base <EPT_URL> pins it instead).
env PROJ_DATA=/usr/share/proj GDAL_DATA=/usr/share/gdal \
  python scripts/fetch_3dep_curl.py --auto --bounds <minx miny maxx maxy> \
  --max-depth 12 --out data/after/3dep_fulldensity.laz

# 2. last-return filter (rn==nr, singles kept)
python scripts/filter_last_return.py data/after/3dep_fulldensity.laz data/after/3dep_last.laz

# 3. difference (PROJ_DATA UNSET so pip rasterio uses its bundled PROJ)
env -u PROJ_DATA -u GDAL_DATA python scripts/gridded_ground_dod.py \
  data/before/<tile>.laz data/after/3dep_last.laz --bounds <minx miny maxx maxy>
# -> data/derived/final/dod.tif, lod.tif, corrections.json ; figures/final_dod.png
```

Or from Python:

```python
from lidar_diff_icp.pipeline import difference_dem
r = difference_dem("before.laz", "3dep_last.laz", bounds, res=5.0)  # ground_q defaults to 0.50 (median)
# r["dod"], r["lod"], r["corrections"], r["stable_sigma"]
```

## Reusing across Minnesota

The correction is deterministic and lives in the acquisition frame, so the same
`difference_dem` runs on any first-generation tile + its 3DEP overlap; only the
fitted coefficients (`corrections.json`: per-swath alignment, tie, and
per-flightline `f(gps_time)` drift) change. That reusability is the goal — the
Elba tile is the pilot.

Both data sources resolve from a coordinate: `tiles.county_for_lonlat` picks the
MnGeo county directory (verified against the live listing), and
`threedep.resolve_reference` picks the covering gen2 3DEP project (most recent,
non-mosaic) and refuses to proceed unless its boundary fully covers the tile bbox.

## Forest structure and large clouds

- **Forest metrics** (`analysis/forest_metrics_pfs.py`) — per-cell canopy cover and
  PAI from the gen2 cloud via **PyForestScan** (plant-area density), a geometry-robust
  land-cover signal that replaces the scan-angle-confounded ground-return "penetration"
  proxy. Runs in the conda `lidar-icp` env, tiles small (400 m) to stay memory-bounded,
  and uses our own `z_after` as the height-above-ground DTM.
- **Large clouds → COPC.** `pdal translate` builds a COPC by holding every point in
  RAM and OOMs on big tiles; **untwine** (conda-forge, isolated env) builds it out-of-core
  (~0.4 GB RAM, external-sorted to disk). The COPC spatial index turns per-tile crops into
  fast indexed seeks — the enabler for forest metrics at statewide scale.
- **Ridgeline tracer** (`analysis/ridgelines/trace_ridgelines.py`) — ridgelines as the
  Scherler & Schwanghart (2020) divide network (via `rivernetworkx.dreich`), generalized
  to run on any tile (grid read from the tile's corrections JSON).

## Data

- **Before (`gen1`)** — MnGeo First-Generation statewide lidar (2008-2012; the SE
  block is Fall 2008); LAZ via MnGeo, organized per county. GPS time present; **no
  CRS embedded** — assign EPSG:26915 (UTM 15N / NAD83). Old-laszip encoding: read
  with **laspy**, not PDAL. `scripts/fetch_tile.py` retrieves a tile by coordinate
  (county resolved automatically) or name.
  https://www.mngeo.state.mn.us/chouse/elevation/lidar_2008-2012.html
- **After (`gen2`)** — USGS Second-Generation 3DEP (pilot: `MN_SEDriftless_2_2021`,
  EPT on AWS `usgs-lidar-public`, stored EPSG:3857). Leaf state differs: gen1 is
  leaf-off (Fall 2008 dormant), gen2 is leaf-on (2021 spring green-up) — the
  mismatch biases the forest DoD (gen2 canopy sits high), which the forest-structure
  tooling above is for.

## Pilot study area

Elba, Minnesota — Whitewater River valley (Winona County): a meander bend
(expected change), an adjacent hillslope (expected stable), away from town.
Reference point 44.101944, −92.004137 (E 579705.72, N 4883677.71, EPSG:26915).

## Repository layout

- `src/lidar_diff_icp/` — the package: `pipeline` (the end-to-end
  `difference_dem`), `coreg` (per-swath alignment, Nuth & Kääb registration, DeLong
  correction surface, along-track drift, `estimate_boresight_roll`), `boresight`
  (scanner-roll self-calibration from flight-line overlap — `estimate_boresight`,
  `apply_boresight`, and the boresight/lateral coupling self-check),
  `references` (deterministic geoid-model datum), `io`, `tiles`
  (county-parametrized gen1 tile discovery + coordinate→county), `threedep`
  (gen2 3DEP project lookup + coverage check), `swathdiff`, `variogram`.
- **Change detection.** `detect.detect_change_standard` is the recommended
  detector: **Wheaton et al. (2010) spatial-coherence Bayesian thresholding**
  (`coherence.py`) + a systematic-error amplitude floor, with an optional
  `wetland.wetland_flag` water mask. It supersedes the earlier two-axis
  `detect.detect_change` (kept for reference). `viz.hillshade` renders shaded
  relief via `gdaldem` (oriented by the geotransform, so it can't be mis-flipped).
- `scripts/` — CLIs: `fetch_3dep_curl`, `filter_last_return`, `gridded_ground_dod`
  (final product), `m3c2_pointcloud` (point-based cross-check),
  `along_track_drift`, `decimation_test`, `fetch_tile`.
- `ground_control/` — the absolute vertical datum subsystem: `control` (epoch-agnostic
  access to both control tables), `lines` (flight-line tracks, one per pass, committed as
  `data/gen1_line_tracks.json`), `same_line` (the site's own lines, marks assigned by
  their returns), `our_surface` (local reconstruction of either epoch's surface anywhere a
  tile is on disk), `site_datum` + `run_site_datum` (both epochs' constants at any site),
  `apply_datum` (gauge-invariant application). `FRAME.md` is the state anchor,
  `REPORT.md` the method record, `INTEGRATION.md` what belongs in `src/` on promotion.
- `analysis/` — documented studies (density decimation, method comparisons).
- `tests/` — regression tests (coreg sign conventions, correction surface,
  synthetic-warp recovery).

## References

DeLong et al. (2022), *Regional-Scale Landscape Response to an Extreme
Precipitation Event From Repeat Lidar*, Earth and Space Science,
[10.1029/2022EA002420](https://doi.org/10.1029/2022EA002420) — the published
analog on this same MN-DNR lidar program (correction surface, uncertainty).

Wheaton, Brasington, Darby & Sear (2010), *Accounting for uncertainty in DEMs from
repeat topographic surveys: improved sediment budgets*, Earth Surface Processes and
Landforms 35(2):136–156, [10.1002/esp.1886](https://doi.org/10.1002/esp.1886) — the
spatial-coherence Bayesian DoD thresholding implemented in `coherence.py`, from
their Geomorphic Change Detection (GCD) software (https://gcd.riverscapes.net).

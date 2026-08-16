# lidar-diff-icp

A reusable workflow to build **bare-earth DEMs of Difference** between the
2008-era Minnesota statewide lidar and modern USGS **3DEP**, correcting the early
survey's navigation error so that what remains is real geomorphic change.

The Elba pilot (below) is the test case; the method is written to generalize to
any 2008 MN tile with 3DEP coverage.

## The problem

The 2008 lidar carries navigation (GNSS/IMU trajectory) error that differs from
one flight line to the next — adjacent swaths over a single tile were flown up to
~1.75 h apart, so each pass has its own, largely independent error state. A single
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

1. **Bare earth = last return, `return_number == number_of_returns`, *including
   single returns*.** Singles dominate flat open ground; dropping them (a common
   `filters.returns groups=last` mistake) empties the agricultural fields.
2. **Ground = a LOW PERCENTILE (10th) of last-return elevation per cell — never
   mean or median.** This is the single most important choice. On rough,
   vegetated, or sloping cells the true ground sits at the *bottom* of the return
   distribution; any central-tendency estimate rides above it, and because 2021 is
   ~14× denser than 2008 that offset becomes **coherent false change** — 16–32% of
   convex hillslopes read as falsely depositional. A low percentile tracks the
   ground and drops that to ~4%. Coherent bias, not incoherent noise, is what
   fools change detection, so this matters more than point-cloud vs raster.
3. **Correct the 2008 points in acquisition-honest order, per point, *before*
   gridding** (not post-hoc on the difference):
   1. per-swath internal alignment (translation, lowest swath pinned);
   2. spatially varying **quadratic tie** to 3DEP on stable ground — removes the
      smooth cross-epoch warp a rigid tie would leave;
   3. **DeLong 400 m correction surface** on flats only, with a **topographic
      position index** floodplain buffer (flow accumulation can't place a channel
      in a flat floodplain);
   4. **per-swath along-track GNSS-drift spline `f(gps_time)`** — the deterministic,
      physical form of the residual error, and the reusable core: the same failure
      mode statewide, only the coefficients differ per tile.
4. **Difference:** gridded low-percentile ground, DoD = 3DEP − 2008 (positive =
   deposition). Cell size (default 5 m) is set by the sparse 2008 density
   (~0.8 pts/m² → ~20 points per 5 m cell for a stable percentile).

Convention, held everywhere: **DoD is `after − before`; red = erosion, blue =
deposition; standard NW (315°/45°) hillshade.**

## Uncertainty

- **Stable-ground 1σ** (empirical NMAD on low-slope, non-floodplain ground) is the
  trustworthy number — ~0.09 m on the pilot.
- A **per-cell LoD** (`lod.tif`) from within-cell spread guides where change is
  detectable, but it is *conservative on slopes* (intra-cell relief inflates it
  yet cancels in the difference), so quote the empirical σ, not the per-cell LoD.

## Quick start

```bash
python3 -m venv --system-site-packages lidar-icp   # over apt geospatial libs
pip install -e .

# 1. fetch 3DEP over the tile's bbox (EPT via curl; readers.ept is unreliable here)
env PROJ_DATA=/usr/share/proj GDAL_DATA=/usr/share/gdal \
  python scripts/fetch_3dep_curl.py --base <EPT_URL> --bounds <minx miny maxx maxy> \
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
r = difference_dem("before.laz", "3dep_last.laz", bounds, res=5.0, ground_q=0.10)
# r["dod"], r["lod"], r["corrections"], r["stable_sigma"]
```

## Reusing across Minnesota

The correction is deterministic and lives in the acquisition frame, so the same
`difference_dem` runs on any 2008 tile + its 3DEP overlap; only the fitted
coefficients (`corrections.json`: per-swath alignment, tie, and per-flightline
`f(gps_time)` drift) change. That reusability is the goal — the Elba tile is the
pilot.

## Data

- **Before** — MN statewide lidar, SE block, Fall 2008; LAZ via MnGeo. GPS time
  present; **no CRS embedded** — assign EPSG:26915 (UTM 15N / NAD83). Old-laszip
  encoding: read with **laspy**, not PDAL. `scripts/fetch_tile.py` retrieves a
  tile by coordinate or name.
  https://www.mngeo.state.mn.us/chouse/elevation/lidar_2008-2012.html
- **After** — USGS 3DEP (pilot: `MN_SEDriftless_2_2021`, EPT on AWS
  `usgs-lidar-public`, stored EPSG:3857; both 2008 and 2021 are leaf-off).

## Pilot study area

Elba, Minnesota — Whitewater River valley (Winona County): a meander bend
(expected change), an adjacent hillslope (expected stable), away from town.
Reference point 44.101944, −92.004137 (E 579705.72, N 4883677.71, EPSG:26915).

## Repository layout

- `src/lidar_diff_icp/` — the package: `pipeline` (the end-to-end
  `difference_dem`), `coreg` (per-swath alignment, quadratic tie, DeLong
  correction surface, along-track drift, Nuth & Kääb), `io`, `tiles`,
  `swathdiff`, `variogram`.
- `scripts/` — CLIs: `fetch_3dep_curl`, `filter_last_return`, `gridded_ground_dod`
  (final product), `m3c2_pointcloud` (point-based cross-check),
  `along_track_drift`, `decimation_test`, `fetch_tile`.
- `analysis/` — documented studies (density decimation, method comparisons).
- `tests/` — regression tests (coreg sign conventions, correction surface,
  synthetic-warp recovery).

## References

DeLong et al. (2022), *Regional-Scale Landscape Response to an Extreme
Precipitation Event From Repeat Lidar*, Earth and Space Science,
[10.1029/2022EA002420](https://doi.org/10.1029/2022EA002420) — the published
analog on this same MN-DNR lidar program (correction surface, uncertainty).

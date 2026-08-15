# lidar-diff-icp

Detecting **real geomorphic change** between two lidar epochs while separating it
from **per-flight-line navigation error** in the earlier survey.

## The problem

The earlier lidar survey carries navigation (GNSS/IMU trajectory) error that
differs from one flight line to the next — adjacent swaths over a single tile
were flown up to ~1.75 h apart, so each pass has its own, largely independent
error state. A single rigid alignment of the whole cloud is therefore wrong: the
error is piecewise per swath. The strategy here is to (1) correct each flight
line individually and (2) certify a difference as real only where it exceeds a
spatially varying detection limit.

Two facts make this tractable:

1. **The overlaps are self-calibrating.** Adjacent swaths were flown minutes to
   hours apart, so there is *no real land-surface change* between them. Every bit
   of elevation discrepancy in a swath overlap is pure acquisition/navigation
   error — the signal we invert to build corrected surfaces.
2. **The earlier data retain per-pulse attribution.** `point_source_id`
   (flight line), GPS time, and scan angle all survive, so each pass can be
   isolated and modeled.

## Study area

Elba, Minnesota — the Whitewater River valley (Winona County). A meander bend
(expected change), an adjacent hillslope (expected stable), sited away from town
to limit anthropogenic change. Reference point: 44.101944, -92.004137
(UTM 15N NAD83: E 579705.72, N 4883677.71).

## Data

- **Before** — Minnesota statewide lidar, SE Minnesota block, flown Fall 2008;
  LAZ point clouds distributed by MnGeo. Point format 1 (GPS time present); no
  CRS is embedded — assign EPSG:26915 (UTM 15N / NAD83) explicitly.
  https://www.mngeo.state.mn.us/chouse/elevation/lidar_2008-2012.html
- **After** — USGS 3DEP, 2021.

Data are not committed. `scripts/fetch_tile.py` retrieves a tile reproducibly by
coordinate or tile name.

## Setup

See [`environment/`](environment/). In brief: a Python venv built with
`--system-site-packages` over apt-installed geospatial libraries.

```bash
python3 -m venv --system-site-packages lidar-icp
# then: sudo apt install $(cat environment/apt-packages.txt)
pip install -e .
```

## Usage

```bash
# fetch the tile covering a coordinate
python scripts/fetch_tile.py --lon -92.004137 --lat 44.101944 --out data/before

# inter-swath consistency (the self-calibration check) for a tile
python scripts/swath_consistency.py data/before/4342-29-64.laz
```

## Repository layout

- `src/lidar_diff_icp/` — the reusable package: `tiles` (tile discovery/download),
  `io` (tile reader), `swathdiff` (density-robust inter-swath difference),
  `coreg` (Nuth & Kaeaeb co-registration + free-network swath alignment),
  `variogram` (robust Dowd variogram + correlated-error detection limit).
- `scripts/` — entry-point pipeline tools: `fetch_tile`, `fetch_naip`,
  `fetch_3dep` / `fetch_3dep_curl`, `swath_consistency`, `swath_coregister`,
  `error_structure`, `naip_cover_error`, `lod`, `stitch_swaths`.
- `analysis/` — documented one-off studies and data-QA (see `analysis/README.md`).
- `tests/` — regression tests (e.g. the coreg sign-convention recovery test).
- `environment/` — apt package list and venv/PROJ setup notes.
- `data/`, `figures/` — inputs and outputs (git-ignored; data is re-fetchable).

## Status

Early. 2008 inter-swath error is resolved and stitched into a common frame
(lowest swath pinned as local reference); the independent-cover error model shows
error is driven by forest/steep terrain, with a long (>=370 m) correlation
length. The 2021 3DEP reference is pulled for the patch and its one tested
flight-line pair shows no detectable internal offset. Next: cross-epoch Nuth &
Kaeaeb tie on stable ground (robust to real change), then the difference map.

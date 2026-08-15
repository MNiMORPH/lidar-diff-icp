# analysis/

Documented one-off studies and data-QA that answer specific questions, kept for
reproducibility and provenance. These reuse the validated core in
`src/lidar_diff_icp/` and the pipeline tools in `scripts/`; they are *not*
themselves core pipeline stages.

- **`dep_internal_check.py`** — check a lidar survey's *own* internal
  flight-line consistency (per `PointSourceId`) on ground returns, before
  trusting it as a reference. Used to verify the 2021 3DEP: the one tested
  overlapping flight-line pair showed no detectable relative offset (dx,dy,dz ~
  0.001,0.000,0.000 m). This measures a relative offset; it does not, by itself,
  establish the vendor processing or the full-tile accuracy.

- **`naip_qa.py`** — quality-checks a fetched NAIP mosaic (band/NDVI
  distributions, nodata) and the k-means cover classification, writing figures
  to `figures/`. In this dissected terrain NAIP carries optical shadow on steep
  slopes, so the "forest" cluster conflates canopy with slope; the QA figures
  make that visible.

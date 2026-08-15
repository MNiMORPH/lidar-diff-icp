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

- **`m3c2_lesson.py`** — M3C2 (point-cloud, along-surface-normal change) vs the
  gridded-vertical difference on bare sloped ground. Finding: M3C2 is ~1.2x
  lower-noise on flat but ~2.8x lower on 15-30 deg slopes, because its local-plane
  fit avoids the gridding-to-vertical artifact (worst on slopes). It does NOT fix
  forest/vegetation, and part of the flat gain is just a larger averaging
  footprint. Needs `py4dgeo`.

- **`naip_qa.py`** — quality-checks a fetched NAIP mosaic (band/NDVI
  distributions, nodata) and the k-means cover classification, writing figures
  to `figures/`. In this dissected terrain NAIP carries optical shadow on steep
  slopes, so the "forest" cluster conflates canopy with slope; the QA figures
  make that visible.

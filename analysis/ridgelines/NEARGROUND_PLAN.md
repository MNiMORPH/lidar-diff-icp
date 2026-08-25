# Near-ground PDF -> gen1/gen2 offset: plan (2026-08-26)

## The physical premise (Andy)
**Plants grow VERTICALLY, even on hillslopes.** Intra-cell terrain relief runs ALONG the
slope. These are different geometries and must not share a coordinate.

- vegetation height: vertical, above the local ground plane
- the offset we want (`d_mm`): slope-PERPENDICULAR, gen1 vs gen2
- the two differ by a known factor: perpendicular = vertical x cos(slope) = vertical / |n|

## What went wrong in the first attempt
1. **Wrong coordinate.** `nearground_cells.py` binned `z - z_cell` (cell-centre elevation,
   no plane correction). Intra-cell relief is +-3.54 m x tan(slope) = +-2.0 m at 30 deg,
   against a +-1 m window at 0.02 m bins. On sloping cells the "near-ground PDF" was mostly
   terrain. This invalidates the slope-dependent conclusions, NOT the flat-cell ones.
2. **Pooling across cells.** Spatial windows mixed cells of different slope, cover and
   vegetation. RMS kept falling to 355 m with no turnover -- averaging away a problem the
   coordinate error had created. Per-cell PDFs were behaving; the pooling broke it.
3. My proposed "fix" (go slope-normal) was wrong the other way: dividing by |n| compresses
   plant heights by cos(slope).

## The rebuild
1. **Recompute the cube in VERTICAL height above the local PLANE**:
   `h = z - (z_cell + gx*(x-xc) + gy*(y-yc))`, no division by |n|.
   Same streaming pattern as `nearground_cells.py`; gx, gy from `np.gradient(z_after, res)`
   exactly as `gen1_save_angles_slope.py` does.
2. **Stay per-cell.** No spatial pooling in the first pass. If returns are too few (gen1
   median ~17/cell), pool only over cells MATCHED in slope and cover, never over raw space.
3. **Characterise the PDFs, both epochs**: quantiles, mode, width, skew, and the
   gen2-minus-gen1 CDF difference (which localises the height where the epochs diverge).
4. **Relate PDF features to the slope-perpendicular offset**, carrying cos(slope) explicitly
   rather than fitting it.

## The test that decides it
Calibrate on divides with slope < 12 deg, apply to held-out divides at 12-18, 18-24,
24-30, 30-36 deg. **If terrain contamination was the culprit, the residual should now stay
flat across slope** instead of degrading and flipping sign at 24 deg (the previous result:
-1.0 / -18.1 / -20.5 / -30.3 / +8.8 / +17.8 mm).

Flat => one calibration covers the whole slope range, and the scope-vs-mass-wasting tension
dissolves. Still degrading => the PDF genuinely does not carry the offset on steep ground.

## Established facts this rests on (verified today, do not re-derive)
- Registration (geoid + lateral tie + per-swath align + along-track drift) is applied and
  tested; open non-eroding ground closes to 2-5 mm. Columns `d_mm_corr`, `dz_*_mm` in
  `beam_offset_table.parquet` on both tiles.
- After registration the offset is FLAT vs incidence 0-45 deg. Slope/incidence is NOT a
  driver; the old ~27 deg knee was per-swath misregistration (+4.29 -> -0.03 mm/deg).
- Cover-dependent offset survives on no-change ground: ~0 open, -20 mm at cover 0.3,
  -60 to -130 mm above 0.5. Two tiles.
- Forest TYPE matters: conifer/brush ~20 mm below deciduous at matched cover (z 2.7).
  Above cover 0.5, 80-94% of cells are conifer-like -- cover and type are confounded here.
- Per-cell noise floor on no-change ground: ~48 mm median, ~108-122 mm RMS.

# ELBAEXT build — expanding the Elba pilot across the Whitewater valley

**Goal.** Extend the Elba pilot tile to a larger study area (`elbaext`) so a set of
vertical tie-point coordinates spread across the Whitewater valley fall inside the
grid. The build mirrors `analysis/slope_bias/fulldensity_regrid.py` exactly in
method — the only changes are the extent, the input clouds, and the CSF cache path.

## Extent

- **elbaext bounds (EPSG:26915 / UTM 15N):** E 575600 – 580050, N 4882200 – 4886250
  (4.45 × 4.05 km; 5 m grid → 890 × 810 cells).
- Extends ~1.85 km **west** and ~0.5 km **south** of the current elba tile
  (E 577493–580035, N 4882738–4886238).
- Data were fetched/clipped to a **150 m buffer** on these bounds to avoid grid-edge
  and CSF-edge effects.

## gen1 (2008 MnGeo lidar) — tiles fetched and merged

The original elba tile `4342-29-64` covers only the EAST part of elbaext. Using the
MnGeo statewide centroid index (`tiles.find_tile` / `tiles.header_bbox`), the 2008
tiles whose footprints intersect the elbaext bounds are (all Winona county):

| tile | bbox E | bbox N | role | status |
|------|--------|--------|------|--------|
| 4342-29-64 | 577493–580035 | 4882738–4886238 | east / main | already present |
| 4342-29-63 | 574993–577533 | 4882709–4886208 | west neighbor | **fetched** |
| 4342-30-64 | 577534–580077 | 4879267–4882767 | south (bottom strip) | **fetched** |
| 4342-30-63 | 575032–577574 | 4879238–4882737 | SW corner | **fetched** |
| 4342-28-64 | 577452–579992 | 4886209–4889709 | top sliver (N≥4886209) | **fetched** |
| 4342-28-63 | 574953–577492 | 4886180–4886250 | NW top sliver | **fetched** |

The six tiles were clipped to the elbaext+150 m buffer and merged into
`data/before/elbaext_gen1_merged.laz`.

**PSID is a GLOBAL swath id here — the merge is safe.** Point-source-id values 133–138
repeat across tiles, but by gps_time each PSID is ONE continuous flight line spanning
tile boundaries (e.g. PSID 138 runs gps_time 237314→237448 across three tiles). So
merging tiles under a shared PSID correctly reconstructs each physical swath, which is
exactly what the per-swath `align_swaths` + along-track drift need. gps_time,
return_number, number_of_returns, and classification are all preserved.

- Merged gen1: **17,354,958 points**, extent E 575450–580043 N 4882050–4886400,
  PSIDs [133,134,135,136,137,138], gps_time 100% present.

## gen2 (2021 3DEP) — fresh full-density EPT pull

The existing `3dep2021_fulldensity.laz` covers only the current elba extent, so a fresh
full-density pull was made over the elbaext+buffer with `scripts/fetch_3dep_curl.py`
(`--auto --max-depth 12`, matching the fulldensity recipe):

- Auto-resolved gen2 reference: **MN_SEDriftless_2_2021** (2021); boundary fully covers
  the elbaext bbox (mandatory coverage check passed).
- 9651 overlapping EPT node tiles at depth ≤ 12; reprojected EPSG:3857 → EPSG:26915,
  clipped, streamed to `data/after/elbaext_3dep_fulldensity.laz`.
- gen2: **415,080,034 points**, extent E 575450–580200 N 4882050–4886400, class-2
  ground present (~28–39%). (Note: as with the earlier merged 3DEP products, gps_time
  is zeroed in this delivery; the `after_ground="class2"` path does not need gen2
  gps_time/scan-angle, so this is not a blocker for the DoD.)

## Blocker hit and fixed: CSF OOM → tiled CSF

The first `elbaext_regrid.py` run **died**: PDAL `filters.csf` was **OOM-killed
(SIGKILL 9)**. `classify_ground_csf` is monolithic — it materialises the whole 17.35M-
point gen1 cloud as an uncompressed LAS plus PDAL's own in-memory copy plus a 1 m cloth
over the full ~4.6×4.35 km footprint. The machine had only ~18 GB free (swap exhausted),
and elbaext at ~4× the original elba point count/area exceeded it. The original elba
tile (7.7M pts) fit; elbaext did not. No cache was written on failure (no corruption).

**Fix — tiled CSF into the cache** (`analysis/slope_bias/elbaext_csf_tiled.py`): CSF is
spatially local (cloth relaxes on a 1 m grid, edge reach a few cells), so the gen1 cloud
was classified in a **2×2 grid (~4.3–5.0M pts/tile) with a 150 m overlap halo**, keeping
each tile's CORE (verified to partition every point exactly once) and concatenating the
class-2 ground into `data/csf_cache/elbaext.las`. Peak RAM stayed ~14 GB free. Because
`difference_dem` loads an existing `csf_cache` and SKIPS CSF, the subsequent regrid runs
the identical fulldensity recipe with no library change.

- CSF ground cache: **15,614,665 class-2 ground points**, PSIDs [133–138], gps_time
  100% present, extent E 575450–580043 N 4882050–4886400, 99.7% of 100 m cells covered
  (only a ~50 m NE-corner sliver uncovered — the gen1 east edge, expected).

## Pipeline settings (identical to fulldensity_regrid.py)

```
difference_dem(BEFORE, AFTER, BOUNDS,
    res=5.0, ground="slope_normal", ground_source="csf",
    after_ground="class2", stream=True, robust_stable=True,
    csf_cache="data/csf_cache/elbaext.las")
```

Defaults in effect (unchanged from fulldensity): `ground_q=0.50`, `along_track_drift=True`,
`tie="reference"` (which, with `stream=True` → no in-memory after cloud, falls back to the
`parabola` datum, exactly as on the current elba products), `datum_tilt=True`,
`geoid_datum=None`, `before_crs=MN_GEN1_CRS`.

## Second blocker and the principled fix: parabola-tie guardrail → reference datum

A second run (now with the CSF cache, so CSF was skipped) then hit a deliberate
guardrail: **`RuntimeError: PARABOLA TIE IS DEACTIVATED`** (pipeline.py:727). With
`stream=True` the after cloud is never in memory, so `references.flat_hard_cells`
cannot find the pavement/pad reference surfaces; the pipeline would fall back to the
order-2 parabola datum, which "warps gen1's z and ABSORBS real hillslope change," and
the current code refuses to do that silently. (The `elba_fulldensity` products predate
this guardrail — their `corrections.json` shows `method: "parabola"`. Reproducing them
byte-for-byte today would require `allow_parabola=True`, a throwaway-test override.)

Rather than override, I took the fix the error message itself prescribes and the
principled path the code now enforces: **run NON-streaming with the gen2 ground in
memory so the reference (const+tilt pad) datum runs.** The full 415M cloud won't fit,
so the gen2 **class-2 ground was extracted to its own file**
(`data/after/elbaext_3dep_fd_class2.laz`, 131,988,280 pts, via
`elbaext_extract_gen2_ground.py`, streamed so the extraction itself is memory-safe),
and the regrid was run with `stream=False`. Non-streaming read loaded 113.7M in-bounds
ground points at ~9.4 GB RSS (14 GB free) — comfortable.

This is the one substantive departure from `fulldensity_regrid.py`: **elbaext uses the
reference-plane datum, not the deprecated parabola.** That is an improvement, not a
regression — it is the datum the pipeline is designed to use.

## Results (completed run)

- Grid: **890 × 810 cells** at 5 m (matches the target extent exactly).
- gen2 ground loaded (non-streaming, class2): **113,696,048 points in-bounds**.
- **Cross-epoch datum = `reference_plane`** (NOT parabola): **382,574 flat-hard
  reference cells**, horizontal shift (−0.750, −0.189) m, constant −84.9 mm, tilt
  7.8 mm/km (residual NMAD 60.7 mm; geoid-tilt cross-check 4 mm/km).
- Per-swath internal alignment + along-track drift solved for all six swaths
  **133–138** (the reference swath is 133, the westmost; swaths 135–138 also appear in
  the original elba pilot). Example dx/dy/dz (m): 138 → (1.389, 0.278, −0.024).
- Stable-ground **1σ = 0.052 m**, robust_stable clip fraction 4.6%; median LoD 0.112 m.
- DoD: **97.8% of cells finite** (704,987 / 720,900), median **+4.7 mm**. The known
  +tan(slope) leaf-on/leaf-off signature is present (median DoD by slope bin: −1.7 mm
  at 0–3°, +3.3 at 3–10°, +7.4 at 10–20°, +13.9 at 20–35°, +29.4 at 35–90°),
  consistent with the Elba pilot — the products behave as expected, not anomalously.

## Tie points in grid

The tie-point set is `data/derived/elba_fulldensity/andy_tie_points.npz` (10 points)
plus the one extra at E 579144.2, N 4883394.2 = **11 points total** (the brief said 12;
the data hold 11 — reported honestly). Only 2 fell inside the original elba tile; **all
11 fall inside elbaext.**

## Products written (`data/derived/elbaext/`)

`z_after.npy`, `dod.npy`, `lod.npy`, `slope.npy`, `stable.npy`, `corrections.json`,
`meta.json`. (Same set `elba_fulldensity` saved — dod/lod/z_after/slope/corrections —
plus `stable.npy` and `meta.json`.)

## Scripts

- `analysis/slope_bias/elbaext_regrid.py` — the regrid driver (mirrors
  `fulldensity_regrid.py`; the differences are documented inline).
- `analysis/slope_bias/elbaext_csf_tiled.py` — tiled CSF that fills the gen1 ground cache.
- `analysis/slope_bias/elbaext_extract_gen2_ground.py` — extracts the gen2 class-2 ground.

## Input artifacts

- `data/before/elbaext_gen1_merged.laz` — 6 merged gen1 tiles (17.35M pts).
- `data/before/4342-29-63.laz`, `4342-30-64.laz`, `4342-30-63.laz`, `4342-28-64.laz`,
  `4342-28-63.laz` — the 5 newly fetched gen1 tiles.
- `data/after/elbaext_3dep_fulldensity.laz` — full-density 3DEP pull (415.08M pts).
- `data/after/elbaext_3dep_fd_class2.laz` — gen2 class-2 ground (131.99M pts).
- `data/csf_cache/elbaext.las` — tiled-CSF gen1 ground cache (15.61M class-2 pts).

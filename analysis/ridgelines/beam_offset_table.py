#!/usr/bin/env python3
"""Assemble a PER-BEAM feature table for gen1 CSF ground returns: for every return, its
beam geometry (scan angle, local slope, incidence angle), its elevation offset vs gen2,
forest cover, return intensity, and other per-return beam characteristics.

No binning, no fitting -- this is the raw per-beam table that Step 2 (relate offset to
incidence) and Step 3 (pure-incidence vs direction-dependent) build on.

Alignment is exact and cheap: gen1_csf_angles.npz was written straight from the cached
CSF LAS in file order with NO masking or reindexing (see gen1_save_angles_slope.py), so
npz row i == LAS point i. We re-read the same cached LAS to pull the fields the npz did
not carry (intensity, return structure, edge/overlap flags, gps_time, elevation) and
column-bind them -- no re-CSF, no fuzzy matching.

Continuous forest cover comes from the PyForestScan canopy raster indexed by the same
per-cell code the producer used for the categorical strata.

    env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python \
        analysis/ridgelines/beam_offset_table.py [tile_dir] [cached_las]

Defaults to the elba pilot. Writes <tile_dir>/beam_offset_table.parquet (the CANONICAL
per-beam file we compare against and fit corrections from) plus a small CSV head sample
for eyeballing, and prints a per-column summary and a few sample rows.

Parquet, not CSV, is the source of truth: columnar + typed (compact dtypes below), ~10x
smaller than CSV at this row count, and read natively by pandas/polars/R/DuckDB/QGIS.
"""
import sys, numpy as np, pandas as pd, laspy

import os, json
TILE = sys.argv[1] if len(sys.argv) > 1 else "data/derived/elba_fulldensity"
LAS  = sys.argv[2] if len(sys.argv) > 2 else "data/csf_cache/elba.las"
def _grid(tile):           # (NY, NX, RES) from tile meta/corrections (cell = iy*NX+ix, C-order)
    for fn in ("meta.json", "corrections_geoid.json", "corrections.json"):
        p = f"{tile}/{fn}"
        if os.path.exists(p):
            j = json.load(open(p)); b = j["bounds"]; r = float(j.get("res") or j.get("res_m"))
            return int(j.get("ny") or round((b[3]-b[1])/r)), int(j.get("nx") or round((b[2]-b[0])/r)), r
    raise SystemExit(f"no grid meta in {tile}")
NY, NX, RES = _grid(TILE)

# --- geometry + offset already computed per return (parabola-free; geoid tie applied later) ---
ang = np.load(f"{TILE}/gen1_csf_angles.npz")
n = ang["d_mm"].shape[0]
cell = ang["cell"].astype(np.int64)
print(f"gen1_csf_angles.npz: {n:,} returns")

# --- continuous forest cover at each return's cell (PyForestScan) ---
cc = np.load(f"{TILE}/canopy_cover_pfs.npy")
assert cc.shape == (NY, NX), f"canopy raster {cc.shape} != {(NY, NX)}"
canopy_cover = cc.ravel()[cell].astype(np.float32)

# --- local surface curvature (Laplacian of gen2 elevation) at each return's cell ---
# per-cell covariate: ~0 on planar ridgetops/hillslopes, signed on convex/concave forms.
lap = np.load(f"{TILE}/curv_laplacian.npy")
assert lap.shape == (NY, NX), f"curvature raster {lap.shape} != {(NY, NX)}"
curv_laplacian = lap.ravel()[cell].astype(np.float32)

# --- surface ASPECT (downslope azimuth) at each return's cell, from gen2 elevation ---
# orientation covariate: which way the terrain faces. Degrees clockwise from North, of the
# DOWNSLOPE direction; NaN on near-flat cells where aspect is undefined. Grid rows increase
# northward (cell = iy*NX+ix, iy=(y-Y0)/RES), so np.gradient gives (d/dNorth, d/dEast) uphill.
Zaf = np.load(f"{TILE}/z_after.npy"); assert Zaf.shape == (NY, NX)
gN, gE = np.gradient(Zaf, RES)                      # uphill gradient components (north, east)
aspect = np.degrees(np.arctan2(-gE, -gN)) % 360.0   # downslope azimuth, CW from North
aspect[np.hypot(gE, gN) < 1e-4] = np.nan            # undefined on flat cells
aspect_deg = aspect.ravel()[cell].astype(np.float32)

# --- per-return fields the npz did not carry: read the SAME cached LAS, in order ---
las = laspy.read(LAS)
assert len(las.x) == n, f"LAS {len(las.x):,} != npz {n:,} -- alignment broken"
_dims = set(las.point_format.dimension_names)
def _opt(name, dt):        # optional LAS dim (PF6+ fields absent in PF<=5) -> zeros
    return np.asarray(getattr(las, name), dt) if name in _dims else np.zeros(n, dt)

cols = {
    # beam geometry (from npz)
    "scan_angle":   ang["scan_angle"],          # deg, signed, off-nadir beam angle
    "slope":        ang["slope"],               # deg, local surface slope
    "incidence":    ang["incidence"],           # deg, beam vs local surface normal
    # response
    "d_mm":         ang["d_mm"],                # mm, slope-normal offset of return vs gen2 surface
    # forest cover
    "canopy_cover": canopy_cover,               # PyForestScan cover fraction at the cell
    # local surface form
    "curv_laplacian": curv_laplacian,           # Laplacian of gen2 elevation at the cell (curvature)
    "aspect_deg":   aspect_deg,                 # downslope azimuth (deg CW from N); NaN on flat cells
    "pfs_forest":   ang["pfs_forest"],          # bool, PyForestScan forest mask -- the cross-tile
    "pfs_open":     ang["pfs_open"],            # bool, PyForestScan open mask    -- cover definition
    "core_forest":  ang["core_forest"],         # bool, forest-core stratum (penetration-derived,
    "core_open":    ang["core_open"],           # bool, open-core stratum   -- elba only, ALL FALSE
    "stratum":      ang["stratum"],             # 1 forest / 2 open / 0 other  elsewhere; see pfs_*
    # per-return beam characteristics (from LAS)
    "intensity":    np.asarray(las.intensity, np.uint16),      # raw, uncalibrated return intensity
    "return_number":     np.asarray(las.return_number, np.uint8),
    "number_of_returns": np.asarray(las.number_of_returns, np.uint8),
    "edge_of_flight_line": np.asarray(las.edge_of_flight_line, np.uint8),  # swath-edge flag
    "overlap":      _opt("overlap", np.uint8),                 # flightline-overlap flag (PF6+)
    "scanner_channel":   _opt("scanner_channel", np.uint8),    # PF6+
    "scan_direction_flag": np.asarray(las.scan_direction_flag, np.uint8),
    "gps_time":     np.asarray(las.gps_time, np.float64),      # along-track time (drift / flight line)
    "z":            np.asarray(las.z, np.float64),             # return elevation (m)
    "point_source_id": ang["point_source_id"],                # flight line id
    "cell":         ang["cell"],
    "in_grid":      ang["in_grid"],                            # bool, return falls in the analysis grid
}

df = pd.DataFrame(cols)

# --- per-cell aggregates (each 5 m cell = one bin; all returns in a cell see the same ground) ---
# cell_mean is the TOTAL offset (real change + mean beam-geometry + datum, inseparable here);
# cell_std is the within-cell scatter (reliability / precision floor); d_resid = the per-beam
# within-cell residual that isolates beam geometry with all per-cell effects differenced out.
ingm = df["in_grid"].to_numpy().astype(bool)
sub = df.loc[ingm]
gd = sub.groupby("cell")["d_mm"]
for name, col in (("cell_n", gd.transform("size")),
                  ("cell_mean_d_mm", gd.transform("mean")),
                  ("cell_std_d_mm", gd.transform("std")),      # NaN where cell n == 1
                  ("cell_inc_std", sub.groupby("cell")["incidence"].transform("std"))):
    df[name] = np.float32("nan")
    df.loc[ingm, name] = col.to_numpy(np.float32)
df["d_resid_mm"] = (df["d_mm"] - df["cell_mean_d_mm"]).astype(np.float32)  # within-cell residual

df.to_parquet(f"{TILE}/beam_offset_table.parquet", index=False, compression="zstd")
df.head(50).to_csv(f"{TILE}/beam_offset_table.head.csv", index=False)   # eyeball preview only
print(f"wrote {TILE}/beam_offset_table.parquet  ({df.shape[1]} cols, {len(df):,} rows)")
print(f"wrote {TILE}/beam_offset_table.head.csv  (first 50 rows, preview)\n")

# --- per-column summary (in-grid returns only, where geometry/offset are meaningful) ---
ing = cols["in_grid"].astype(bool)
print(f"summary over {ing.sum():,} in-grid returns:")
print(f"  {'column':22s} {'min':>12s} {'median':>12s} {'max':>12s}  {'note'}")
notes = {"scan_angle": "deg", "slope": "deg", "incidence": "deg", "d_mm": "mm offset vs gen2",
         "canopy_cover": "fraction", "curv_laplacian": "elev Laplacian (curvature)",
         "intensity": "raw DN", "gps_time": "s", "z": "m elev"}
for k, v in cols.items():
    a = np.asarray(v)[ing]
    if a.dtype == bool:
        print(f"  {k:22s} {'':>12s} {'':>12s} {'':>12s}  {a.mean()*100:5.1f}% true")
    else:
        af = a.astype(float)
        print(f"  {k:22s} {np.nanmin(af):12.3f} {np.nanmedian(af):12.3f} {np.nanmax(af):12.3f}"
              f"  {notes.get(k,'')}{'  ('+str(int(np.isnan(af).sum()))+' nan)' if np.isnan(af).any() else ''}")

# --- a few sample rows ---
print("\nsample rows (in-grid):")
idx = np.where(ing)[0][:: max(1, ing.sum() // 5)][:5]
show = ["scan_angle", "slope", "incidence", "d_mm", "canopy_cover", "intensity",
        "return_number", "number_of_returns", "overlap"]
print("  " + " ".join(f"{c:>10s}" for c in show))
for i in idx:
    print("  " + " ".join(f"{float(np.asarray(cols[c])[i]):10.3f}" for c in show))

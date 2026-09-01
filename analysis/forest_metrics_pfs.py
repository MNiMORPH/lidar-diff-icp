#!/usr/bin/env python3
"""Reusable forest-structure metrics via PyForestScan (canopy cover, PAI) for ANY tile.

Replaces the self-rolled ground-return "penetration" proxy (scan-angle/overlap-confounded)
with PyForestScan's plant-area-density canopy cover -- the field-standard, geometry-robust
metric derived from the vertical distribution of returns above a height threshold.

MUST run in the conda `lidar-icp` env (PDAL + GDAL + pyforestscan), with the env's own PROJ:
    PROJ_DATA=/home/awickert/anaconda3/envs/lidar-icp/share/proj \
    GDAL_DATA=/home/awickert/anaconda3/envs/lidar-icp/share/gdal \
    /home/awickert/anaconda3/envs/lidar-icp/bin/python analysis/forest_metrics_pfs.py <TILE_DIR> <AFTER_COPC>

The full LAZ is too big to read whole (OOM on voxelize), so we crop grid-aligned tiles
(+halo) and mosaic. Give it a **COPC** cloud (`pdal translate in.laz out.copc.laz`): the
octree index makes each tile crop a fast indexed seek instead of a full-file scan (plain
LAZ still works via filters.crop, but re-reads the whole file per tile -> minutes each).
Grid (origin, res) is read from the tile's corrections*/meta JSON; outputs
canopy_cover_pfs.npy + pai_pfs.npy aligned to that grid, plus forest/open masks by threshold.

WHAT THIS MEASURES, AND WHAT IT DOES NOT
----------------------------------------
Canopy cover is the fraction of plant area ABOVE ``--min-height`` (default 2.0 m). It is a
CANOPY metric and it is BLIND TO UNDERGROWTH: measured at Elba, 17.3% of the cells whose
near-ground return fraction (0.15-4 m) exceeds 0.25 -- thick understory -- are labelled
`open` by this filter. If the question is undergrowth, this is the wrong instrument.

THE THRESHOLDS ARE DECLARED, NOT CALIBRATED
-------------------------------------------
``--forest-cover`` and ``--open-cover`` are choices. An attempt to calibrate them against
the only epoch-matched ground truth available -- the 2021 survey's own NVA (non-vegetated)
/ VVA (vegetated) checkpoint classes, 227 and 162 marks -- does NOT support any threshold:

    canopy cover r10                 AUC 0.548   NVA median 0.000   VVA median 0.000
    near-ground fraction 0.15-4.0 m  AUC 0.739   NVA median 0.022   VVA median 0.092
    at the default 0.5: sensitivity 0.006, specificity 0.996

AUC 0.548 is chance. The honest reading is NOT that cover is a bad canopy metric, but that
control marks are sited in the open by design and contain almost no canopy (median 0.000 in
BOTH classes), so this ground truth cannot calibrate a canopy threshold. No threshold is
therefore proposed here, and the defaults are carried forward unchanged and unjustified.
Anything downstream that depends on their exact value should say so.

THE MASKS ARE NOT A PARTITION. ``forest`` and ``open`` leave a gap between them, and at
Elba that gap is most of the tile (63.1% of cells with cover, against forest 1.5% and open
35.4%). The run prints the unclassified fraction so it cannot be dropped silently.
"""
import os, sys, json, argparse, subprocess, tempfile, numpy as np

CENV = "/home/awickert/anaconda3/envs/lidar-icp"


def _grid_from_tile(tile_dir):
    for fn in ("corrections_geoid.json", "corrections.json", "meta.json"):
        p = os.path.join(tile_dir, fn)
        if os.path.exists(p):
            j = json.load(open(p)); b = j.get("bounds"); r = j.get("res_m") or j.get("res")
            if b and r:
                return float(b[0]), float(b[1]), float(b[2]), float(b[3]), float(r)
    raise SystemExit(f"no bounds/res in {tile_dir}")


def _crop(after_cloud, xmin, ymin, xmax, ymax, out_laz):
    """Crop the after cloud to a bbox (keep all returns) -> out_laz.

    A **COPC** input (``*.copc.laz``) is read via ``readers.copc`` with ``bounds`` on the
    reader, which pushes the box down to the octree index -> fast, indexed seek (no full-file
    scan). A plain LAZ falls back to ``filters.crop`` (slow: scans the whole file). Convert
    once with ``pdal translate in.laz out.copc.laz`` to get the fast path.
    """
    bnds = f"([{xmin},{xmax}],[{ymin},{ymax}])"
    if ".copc" in os.path.basename(after_cloud):
        pipe = [{"type": "readers.copc", "filename": after_cloud, "bounds": bnds}, out_laz]
    else:                                                     # slow fallback for un-indexed LAZ
        pipe = [after_cloud, {"type": "filters.crop", "bounds": bnds}, out_laz]
    with tempfile.NamedTemporaryFile("w", suffix=".json", delete=False) as f:
        json.dump(pipe, f); pj = f.name
    env = dict(os.environ, PROJ_DATA=f"{CENV}/share/proj", GDAL_DATA=f"{CENV}/share/gdal")
    subprocess.run([f"{CENV}/bin/pdal", "pipeline", pj], check=True, env=env,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    os.unlink(pj)
    return os.path.exists(out_laz)


def _write_dtm(z, x0, y0, res, crs, out_tif):
    """Write our ground surface (z_after, (ny,nx), row 0 = y0 bottom) as a GeoTIFF DTM for
    fast height-above-ground (raster subtraction, not Delaunay). GeoTIFF is top-down -> flip."""
    import rasterio
    from rasterio.transform import from_origin
    ny, nx = z.shape; nod = -9999.0
    zt = np.where(np.isfinite(z), z, nod).astype("float32")[::-1]     # top-down: row 0 = ymax
    with rasterio.open(out_tif, "w", driver="GTiff", height=ny, width=nx, count=1,
                       dtype="float32", crs=crs, transform=from_origin(x0, y0 + ny * res, res, res),
                       nodata=nod) as dst:
        dst.write(zt, 1)


def canopy_cover_raster(after_laz, bounds, res, *, crs="EPSG:26915", tile_m=400.0,
                        halo_m=25.0, voxel_z=1.0, min_height=2.0, dtm_path=None):
    # tile_m kept small on purpose: per-tile memory scales with tile area (~24M pts at 1km,
    # ~3M at 400m). COPC makes small tiles cheap, so favour many small tiles over OOM risk.
    """Canopy-cover + PAI rasters (ny, nx) aligned to (bounds, res), tiled to fit memory.

    bounds=(X0,Y0,X1,Y1). Returns (canopy_cover, pai) with NaN where no data. The gen2
    cloud is cropped to grid-aligned tiles with a halo (voxels need neighbours); only the
    core of each tile is written into the mosaic.
    """
    import pyforestscan.handlers as H, pyforestscan.calculate as C
    X0, Y0, X1, Y1 = bounds
    nx = int(round((X1 - X0) / res)); ny = int(round((Y1 - Y0) / res))
    cover = np.full((ny, nx), np.nan); pai = np.full((ny, nx), np.nan)
    nxt = int(np.ceil((X1 - X0) / tile_m)); nyt = int(np.ceil((Y1 - Y0) / tile_m))
    tmpdir = tempfile.mkdtemp(prefix="pfs_")
    for j in range(nyt):
        for i in range(nxt):
            tx0 = X0 + i * tile_m; tx1 = min(tx0 + tile_m, X1)
            ty0 = Y0 + j * tile_m; ty1 = min(ty0 + tile_m, Y1)
            lp = os.path.join(tmpdir, f"t_{i}_{j}.laz")
            if not _crop(after_laz, tx0 - halo_m, ty0 - halo_m, tx1 + halo_m, ty1 + halo_m, lp):
                continue
            arr = (H.read_lidar(lp, crs, hag_dtm=True, dtm=dtm_path) if dtm_path
                   else H.read_lidar(lp, crs, hag=True))       # DTM-HAG (fast) vs Delaunay-HAG
            a0 = arr[0] if isinstance(arr, list) else arr
            os.unlink(lp)
            if a0 is None or len(a0) == 0:
                continue
            vox, ext = C.assign_voxels(a0, (res, res, voxel_z))       # ext = [xmin,xmax,ymin,ymax]
            pad = C.calculate_pad(vox, voxel_height=voxel_z)
            cc = C.calculate_canopy_cover(pad, voxel_height=voxel_z, min_height=min_height)
            pi = C.calculate_pai(pad, voxel_height=voxel_z, min_height=min_height)
            # assign_voxels returns hist[x_index, y_index]: x_index from ext xmin (increasing x);
            # y_index=0 is the TOP (ext ymax), increasing DOWNWARD (image convention). extent=[xmin,xmax,ymin,ymax].
            ex0, ex1, ey0, ey1 = ext
            for arr2d, dst in ((cc, cover), (pi, pai)):
                ncol, nrow = arr2d.shape
                for c in range(ncol):
                    xc = ex0 + (c + 0.5) * res
                    if not (tx0 <= xc < tx1):                         # keep CORE only (drop halo)
                        continue
                    gi = int((xc - X0) / res)
                    for r in range(nrow):
                        yc = ey1 - (r + 0.5) * res                    # row 0 = top (ymax), y decreases downward
                        if not (ty0 <= yc < ty1):
                            continue
                        gj = int((yc - Y0) / res)
                        v = arr2d[c, r]
                        if np.isfinite(v) and 0 <= gj < ny and 0 <= gi < nx:
                            dst[gj, gi] = v
    return cover, pai


def classify(cover, forest_cover, open_cover):
    """forest / open masks and a full accounting of every cell.

    The two masks are NOT complementary: cells between the thresholds belong to neither,
    and cells without cover belong to nothing. Both counts are returned so a caller cannot
    lose them by accident.
    """
    if not (open_cover < forest_cover):
        raise ValueError(f"open_cover ({open_cover}) must be below forest_cover "
                         f"({forest_cover}); otherwise the masks overlap and a cell is "
                         f"both forest and open")
    cover = np.asarray(cover, float)
    fin = np.isfinite(cover)
    forest = fin & (cover >= forest_cover)
    openg = fin & (cover <= open_cover)
    acct = dict(n_cells=int(cover.size), n_cover=int(fin.sum()),
                n_forest=int(forest.sum()), n_open=int(openg.sum()),
                n_between=int((fin & ~forest & ~openg).sum()),
                n_nocover=int((~fin).sum()))
    assert acct["n_forest"] + acct["n_open"] + acct["n_between"] == acct["n_cover"], \
        "the accounting does not close: a cell is both forest and open"
    return forest, openg, acct


def report(acct, forest_cover, open_cover, cover=None):
    n = max(acct["n_cover"], 1)
    pct = lambda k: 100.0 * acct[k] / n
    lines = [f"canopy cover: {acct['n_cover']:,} of {acct['n_cells']:,} cells carry a value"
             f"  ({acct['n_nocover']:,} do not)"]
    if cover is not None and acct["n_cover"]:
        lines.append(f"  cover: median {np.nanmedian(cover):.3f}  p90 "
                     f"{np.nanpercentile(cover, 90):.3f}  max {np.nanmax(cover):.3f}")
    lines += [
        f"  forest (cover >= {forest_cover:g}) {acct['n_forest']:8,d}  {pct('n_forest'):5.1f}%",
        f"  open   (cover <= {open_cover:g}) {acct['n_open']:8,d}  {pct('n_open'):5.1f}%",
        f"  NEITHER, between the thresholds {acct['n_between']:8,d}  {pct('n_between'):5.1f}%"
        f"   <- classified by nothing; do not let this vanish",
        "  thresholds are DECLARED, not calibrated -- see the module docstring.",
        "  this is a CANOPY metric and is blind to undergrowth.",
    ]
    return "\n".join(lines)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("tile_dir"); ap.add_argument("after_laz")
    ap.add_argument("--crs", default="EPSG:26915")
    ap.add_argument("--tile-m", type=float, default=400.0)   # small tiles bound per-tile memory
    ap.add_argument("--min-height", type=float, default=2.0)
    ap.add_argument("--forest-cover", type=float, default=0.5, help="cover>=this -> forest")
    ap.add_argument("--open-cover", type=float, default=0.1, help="cover<=this -> open")
    a = ap.parse_args(); d = a.tile_dir.rstrip("/")
    X0, Y0, X1, Y1, res = _grid_from_tile(d)
    dtm_tif = os.path.join(d, "z_after_dtm.tif")           # our ground -> fast DTM-based HAG
    _write_dtm(np.load(os.path.join(d, "z_after.npy")), X0, Y0, res, a.crs, dtm_tif)
    cover, pai = canopy_cover_raster(a.after_laz, (X0, Y0, X1, Y1), res, crs=a.crs,
                                     tile_m=a.tile_m, min_height=a.min_height, dtm_path=dtm_tif)
    np.save(os.path.join(d, "canopy_cover_pfs.npy"), cover)
    np.save(os.path.join(d, "pai_pfs.npy"), pai)
    forest, openg, acct = classify(cover, a.forest_cover, a.open_cover)
    np.save(os.path.join(d, "forest_pfs.npy"), forest); np.save(os.path.join(d, "open_pfs.npy"), openg)
    print(report(acct, a.forest_cover, a.open_cover, cover))
    print(f"saved -> {d}/  (canopy_cover_pfs, pai_pfs, forest_pfs, open_pfs)")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Rasterize a point-based change product into a DEM of Difference + error raster.

Takes any change point cloud carrying ``m3c2`` (change, after - before) and
``lod`` dims -- e.g. the output of scripts/m3c2_pointcloud.py -- and writes a
gridded DoD GeoTIFF plus a rasterized error GeoTIFF (EPSG:26915), so the
point-based path yields the same raster deliverables as the gridded workflow.

Error options:
  * default: per-cell median of the per-point LoD already in the cloud;
  * ``--calibrated-error``: the calibrated heteroscedastic model (xdem / Hugonnet
    2022) fit on stable ground of the rasterized DoD (needs xdem; PROJ_DATA unset).

    env -u PROJ_DATA -u GDAL_DATA python scripts/rasterize_change.py \
      data/derived/change_core2008_robust.laz --res 5 --outdir data/derived/m3c2_raster
"""
import argparse
from pathlib import Path
import numpy as np
import laspy

from lidar_diff_icp.pipeline import rasterize, heteroscedastic_lod
from lidar_diff_icp import coreg, terrain


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("change_laz", help="point cloud with m3c2 + lod dims")
    ap.add_argument("--bounds", nargs=4, type=float, help="minx miny maxx maxy; else from data")
    ap.add_argument("--res", type=float, default=5.0)
    ap.add_argument("--agg", choices=("median", "mean"), default="median")
    ap.add_argument("--calibrated-error", action="store_true",
                    help="use the xdem heteroscedastic model instead of the per-point LoD")
    ap.add_argument("--valley-top", dest="valley_top", default="histogram",
                    help="valley top for the stable mask: an elevation in metres, "
                         "'registry', or 'histogram'. Never chosen for you.")
    ap.add_argument("--tile-dir", dest="tile_dir", default=None,
                    help="tile directory, for --valley-top registry/histogram")
    ap.add_argument("--outdir", default="data/derived/raster")
    ap.add_argument("--figdir", default="figures")
    a = ap.parse_args()

    h = laspy.read(a.change_laz)
    x = np.asarray(h.x); y = np.asarray(h.y); z = np.asarray(h.z)
    change = np.asarray(h.m3c2); lod_pt = np.asarray(h.lod)
    if a.bounds:
        bounds = tuple(a.bounds)
    else:  # snap the data extent out to whole cells
        bounds = (np.floor(x.min() / a.res) * a.res, np.floor(y.min() / a.res) * a.res,
                  np.ceil(x.max() / a.res) * a.res, np.ceil(y.max() / a.res) * a.res)
    res = a.res; X0, Y0, X1, Y1 = bounds
    nx = int(round((X1 - X0) / res)); ny = int(round((Y1 - Y0) / res))

    dod = rasterize(x, y, change, bounds, res, a.agg)
    Zg = rasterize(x, y, z, bounds, res, a.agg)          # surface for hillshade / terrain

    if a.calibrated_error:
        # ONE definition, from terrain.py -- slope, the LoD's curvature covariate and the
        # stable mask together. This hand-rolled all three, and its slope used sigma=1
        # where the pipeline uses 2, so the LoD it reported was never quite the pipeline's.
        _tm = terrain.terrain_masks(Zg, res, valley_top_m=a.valley_top,
                                    tile_dir=a.tile_dir)
        sdeg = _tm["slope_deg"]; curv = _tm["abs_curv"]
        stable = _tm["stable"] & np.isfinite(dod)
        err = heteroscedastic_lod(dod, sdeg, curv, stable)
        emethod = "xdem heteroscedastic (slope,curv), calibrated on stable ground"
        if err is None:
            err = rasterize(x, y, lod_pt, bounds, res, a.agg); emethod = "per-point LoD (xdem unavailable)"
    else:
        err = rasterize(x, y, lod_pt, bounds, res, a.agg)
        emethod = "per-point LoD (rasterized)"

    ex = np.isfinite(dod) & np.isfinite(err)
    print(f"rasterized {a.change_laz} @ {res:.0f} m: {np.isfinite(dod).sum():,} cells; "
          f"median |change| {np.nanmedian(np.abs(dod)):.3f} m; median error {np.nanmedian(err):.3f} m; "
          f"{100*np.mean(np.abs(dod[ex]) > err[ex]):.0f}% exceed error; error = {emethod}", flush=True)

    Path(a.outdir).mkdir(parents=True, exist_ok=True)
    _tif(dod, res, X0, Y0, ny, f"{a.outdir}/dod.tif")
    _tif(err, res, X0, Y0, ny, f"{a.outdir}/error.tif")
    print(f"wrote {a.outdir}/dod.tif, error.tif", flush=True)
    _fig(Zg, dod, err, res, X0, Y0, nx, ny, a.figdir)


def _tif(arr, res, x0, y0, ny, out):
    import rasterio
    from rasterio.transform import from_origin
    with rasterio.open(out, "w", driver="GTiff", height=arr.shape[0], width=arr.shape[1],
                       count=1, dtype="float32", crs="EPSG:26915", nodata=np.nan,
                       transform=from_origin(x0, y0 + ny * res, res, res)) as d:
        d.write(np.flipud(arr).astype("float32"), 1)


def _fig(Z, dod, err, res, X0, Y0, nx, ny, figdir):
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    from matplotlib.colors import LightSource
    Path(figdir).mkdir(exist_ok=True)
    ls = LightSource(azdeg=315, altdeg=45)
    q = np.nan_to_num(Z, nan=np.nanmin(Z))
    hs = np.flipud(ls.hillshade(np.flipud(q), vert_exag=2, dx=res, dy=res))
    ext = (X0, X0 + nx * res, Y0, Y0 + ny * res); v = 0.3
    fig, ax = plt.subplots(1, 2, figsize=(15, 9))
    ax[0].imshow(hs, extent=ext, origin="lower", cmap="gray", alpha=0.6)
    im0 = ax[0].imshow(dod, extent=ext, origin="lower", cmap="RdBu", vmin=-v, vmax=v)
    ax[0].set_title("rasterized DoD, after - before (m)\nred = erosion, blue = deposition")
    fig.colorbar(im0, ax=ax[0], shrink=0.6, extend="both")
    im1 = ax[1].imshow(err, extent=ext, origin="lower", cmap="viridis", vmin=0, vmax=0.2)
    ax[1].set_title("rasterized error (m)")
    fig.colorbar(im1, ax=ax[1], shrink=0.6, extend="max")
    for a in ax: a.set_xlabel("Easting (m)"); a.set_ylabel("Northing (m)")
    fig.savefig(f"{figdir}/rasterized_change.png", dpi=130, bbox_inches="tight")
    print(f"wrote {figdir}/rasterized_change.png", flush=True)


if __name__ == "__main__":
    main()

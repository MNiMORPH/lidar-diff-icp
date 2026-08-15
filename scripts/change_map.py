#!/usr/bin/env python3
"""Cross-epoch change map via M3C2, thresholded by per-point level of detection.

Difference engine = M3C2 (point-based, along local surface normal), which avoids
the gridding-to-vertical artifact that inflates slope noise. Steps:

  1. internally align the 2008 swaths (align_swaths, lowest swath pinned) and
     rigid-tie the mosaic to the 2021 3DEP on stable ground (Nuth & Kaeaeb);
  2. M3C2 between 2008 ground and 2021 ground on a regular core-point grid,
     giving a change distance AND a per-point LoD (from local roughness +
     a registration term) at every cell;
  3. flag change as significant only where |distance| > LoD.

Positive distance = 2021 surface higher than 2008 (deposition); negative =
erosion. Forest cells have sparse ground returns -> high LoD -> not flagged,
so the LoD itself down-weights unreliable (vegetated/steep) ground.

Example:
    env PROJ_DATA=/usr/share/proj GDAL_DATA=/usr/share/gdal \
      python scripts/change_map.py data/before/4342-29-64.laz \
      data/after/3dep2021_fulltile.laz --bounds 577492.8 4882737.6 580035.0 4886238.3 \
      --res 3 --outdir data/derived --figdir figures
"""
import argparse
from pathlib import Path
import numpy as np
import laspy
import py4dgeo

from lidar_diff_icp import io, coreg
from lidar_diff_icp.swathdiff import _median_grid


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("before_laz"); ap.add_argument("after_laz")
    ap.add_argument("--bounds", nargs=4, type=float, required=True,
                    metavar=("MINX", "MINY", "MAXX", "MAXY"))
    ap.add_argument("--res", type=float, default=3.0, help="core-point grid (m)")
    ap.add_argument("--normal-radius", type=float, default=4.0)
    ap.add_argument("--cyl-radius", type=float, default=2.0)
    ap.add_argument("--registration-error", type=float, default=0.02,
                    help="co-registration uncertainty added to the LoD (m)")
    ap.add_argument("--outdir", default="data/derived")
    ap.add_argument("--figdir")
    a = ap.parse_args()
    X0, Y0, X1, Y1 = a.bounds
    nx = int(round((X1 - X0) / a.res)); ny = int(round((Y1 - Y0) / a.res))

    # --- 2008: internal align + rigid tie to 2021 (applied to ALL points) ---
    f = laspy.read(a.before_laz)
    x8 = np.asarray(f.x); y8 = np.asarray(f.y); z8 = np.asarray(f.z)
    ps8 = np.asarray(f.point_source_id); cl8 = np.asarray(f.classification)
    rn8 = np.asarray(f.return_number); nr8 = np.asarray(f.number_of_returns)
    pc = io.PointCloud(x8, y8, z8, ps8, cl8, np.zeros_like(z8),
                       np.zeros_like(ps8), io.MN_2008_CRS)
    corr, _, _ = coreg.align_swaths(pc, ref=int(ps8.min()))
    xc, yc, zc = x8.copy(), y8.copy(), z8.copy()
    for s, (dx, dy, dz) in corr.items():
        m = ps8 == s; xc[m] += dx; yc[m] += dy; zc[m] += dz

    g = laspy.read(a.after_laz)
    x2 = np.asarray(g.x); y2 = np.asarray(g.y); z2 = np.asarray(g.z)
    c2 = np.asarray(g.classification)
    g21 = (c2 == 2) & (x2 >= X0) & (x2 < X1) & (y2 >= Y0) & (y2 < Y1)
    terr08 = (nr8 == 1) & ~np.isin(cl8, [5, 6, 9]) & (xc >= X0) & (xc < X1) & (yc >= Y0) & (yc < Y1)
    Z21 = _median_grid(x2[g21], y2[g21], z2[g21], a.res, X0, Y0, nx, ny)
    Z08 = _median_grid(xc[terr08], yc[terr08], zc[terr08], a.res, X0, Y0, nx, ny)
    tie = coreg.nuth_kaab(Z21, Z08, a.res)
    xc += tie.dx; yc += tie.dy; zc += tie.dz
    print(f"rigid tie 2008->2021: dx={tie.dx:+.3f} dy={tie.dy:+.3f} dz={tie.dz:+.3f} m")

    # --- point clouds: ground class both epochs ---
    grd08 = (cl8 == 2) & (xc >= X0) & (xc < X1) & (yc >= Y0) & (yc < Y1)
    p08 = np.column_stack([xc[grd08], yc[grd08], zc[grd08]]).astype(np.float64)
    p21 = np.column_stack([x2[g21], y2[g21], z2[g21]]).astype(np.float64)

    # --- core points: regular grid, z from the 2008 ground surface ---
    Z08g = _median_grid(xc[grd08], yc[grd08], zc[grd08], a.res, X0, Y0, nx, ny)
    gy_i, gx_i = np.mgrid[0:ny, 0:nx]
    Xc = X0 + (gx_i + 0.5) * a.res; Yc = Y0 + (gy_i + 0.5) * a.res
    valid = np.isfinite(Z08g)
    core = np.column_stack([Xc[valid], Yc[valid], Z08g[valid]]).astype(np.float64)
    print(f"2008 ground: {len(p08):,}  2021 ground: {len(p21):,}  core cells: {len(core):,}")

    m3c2 = py4dgeo.M3C2(epochs=(py4dgeo.Epoch(p08), py4dgeo.Epoch(p21)),
                        corepoints=core, normal_radii=(a.normal_radius,),
                        cyl_radius=a.cyl_radius, max_distance=15.0,
                        registration_error=a.registration_error)
    dist, unc = m3c2.run()
    lod = unc["lodetection"]

    D = np.full((ny, nx), np.nan); L = np.full((ny, nx), np.nan)
    D[valid] = dist; L[valid] = lod
    sig = np.isfinite(D) & (np.abs(D) > L)
    change = np.where(sig, D, np.nan)
    n_ok = np.isfinite(D).sum()
    print(f"M3C2 computed at {n_ok:,} cells; significant change (|d|>LoD): "
          f"{100*sig.sum()/max(n_ok,1):.1f}%  (median LoD {np.nanmedian(L):.3f} m)")

    Path(a.outdir).mkdir(parents=True, exist_ok=True)
    _write(D, a.res, X0, Y0, ny, f"{a.outdir}/m3c2_distance.tif")
    _write(L, a.res, X0, Y0, ny, f"{a.outdir}/m3c2_lod.tif")
    _write(change, a.res, X0, Y0, ny, f"{a.outdir}/change_significant.tif")
    if a.figdir:
        _fig(Z21, change, (X0, X1, Y0, Y1), a.figdir)


def _write(arr, res, x0, y0, ny, out):
    import rasterio
    from rasterio.transform import from_origin
    with rasterio.open(out, "w", driver="GTiff", height=arr.shape[0],
                       width=arr.shape[1], count=1, dtype="float32",
                       crs="EPSG:26915", nodata=np.nan,
                       transform=from_origin(x0, y0 + ny * res, res, res)) as d:
        d.write(np.flipud(arr).astype("float32"), 1)
    print(f"  wrote {out}")


def _fig(Z21, change, ext, figdir):
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LightSource
    Path(figdir).mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(7, 8))
    ls = LightSource(azdeg=315, altdeg=45)
    hs = ls.hillshade(np.nan_to_num(Z21, nan=np.nanmin(Z21)), vert_exag=2)
    ax.imshow(hs, extent=ext, cmap="gray", origin="lower", alpha=0.8)
    im = ax.imshow(change, extent=ext, origin="lower", cmap="RdBu_r",
                   vmin=-0.5, vmax=0.5)
    ax.set_title("Significant change 2008->2021 (M3C2 > LoD)\n"
                 "blue = gain/deposition, red = loss/erosion")
    ax.set_xlabel("Easting (m)"); ax.set_ylabel("Northing (m)")
    fig.colorbar(im, ax=ax, shrink=0.6, label="change along surface normal (m)")
    out = Path(figdir) / "change_map_m3c2.png"
    fig.savefig(out, dpi=120, bbox_inches="tight"); print(f"  wrote {out}")


if __name__ == "__main__":
    main()

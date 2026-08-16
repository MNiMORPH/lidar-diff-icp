#!/usr/bin/env python3
# SUPERSEDED (kept for git history; do not use for products): -> coreg.tie_polynomial
"""Cross-epoch tie and difference: align the internally-consistent 2008 mosaic
to the 2021 3DEP by robust Nuth & Kaeaeb, then difference into a change map.

Steps: (1) internally align the 2008 swaths (align_swaths, lowest swath pinned);
(2) grid 2008 single-return terrain and 2021 ground to bare-earth DEMs over the
common bounds; (3) robust N&K (iterative 3-sigma rejection) ties 2008 -> 2021,
so real change is rejected FROM THE FIT (not from the difference); (4) write the
2021 - aligned2008 difference and a figure.

The tie removes a rigid offset; its dz absorbs any vertical-datum difference, so
the raw dz is not physical. A spatially uniform real change is unrecoverable
(absorbed by dz); the differential pattern is what survives.

Example:
    env PROJ_DATA=/usr/share/proj GDAL_DATA=/usr/share/gdal \
      python scripts/cross_epoch_tie.py data/before/4342-29-64.laz \
      data/after/3dep2021_subpatch.laz \
      --bounds 578014 4883738 579514 4885238 --res 3 \
      --out data/derived/diff_2021_2008.tif --figdir figures
"""
import argparse
from pathlib import Path
import numpy as np
import laspy

from lidar_diff_icp import io, coreg
from lidar_diff_icp.swathdiff import _median_grid


def _nmad(v):
    return 1.4826 * np.median(np.abs(v - np.median(v)))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("before_laz"); ap.add_argument("after_laz")
    ap.add_argument("--bounds", nargs=4, type=float, required=True,
                    metavar=("MINX", "MINY", "MAXX", "MAXY"))
    ap.add_argument("--res", type=float, default=3.0)
    ap.add_argument("--out", help="difference GeoTIFF (2021 - aligned 2008)")
    ap.add_argument("--figdir")
    a = ap.parse_args()
    X0, Y0, X1, Y1 = a.bounds
    nx = int(round((X1 - X0) / a.res)); ny = int(round((Y1 - Y0) / a.res))

    # 2008: build PointCloud from laspy (need return number), internally align
    f = laspy.read(a.before_laz)
    x8 = np.asarray(f.x); y8 = np.asarray(f.y); z8 = np.asarray(f.z)
    ps8 = np.asarray(f.point_source_id); cl8 = np.asarray(f.classification)
    rn8 = np.asarray(f.return_number); nr8 = np.asarray(f.number_of_returns)
    pc = io.PointCloud(x8, y8, z8, ps8, cl8, np.zeros_like(z8),
                       np.zeros_like(ps8), io.MN_2008_CRS)
    corr, _, _ = coreg.align_swaths(pc, ref=int(np.min(ps8)))
    xc, yc, zc = x8.copy(), y8.copy(), z8.copy()
    for s, (dx, dy, dz) in corr.items():
        m = ps8 == s; xc[m] += dx; yc[m] += dy; zc[m] += dz
    t8 = (nr8 == 1) & ~np.isin(cl8, [5, 6, 9]) & \
         (xc >= X0) & (xc < X1) & (yc >= Y0) & (yc < Y1)
    Z08 = _median_grid(xc[t8], yc[t8], zc[t8], a.res, X0, Y0, nx, ny)

    # 2021: ground surface
    g = laspy.read(a.after_laz)
    x2 = np.asarray(g.x); y2 = np.asarray(g.y); z2 = np.asarray(g.z)
    gm = (np.asarray(g.classification) == 2) & \
         (x2 >= X0) & (x2 < X1) & (y2 >= Y0) & (y2 < Y1)
    Z21 = _median_grid(x2[gm], y2[gm], z2[gm], a.res, X0, Y0, nx, ny)

    both = np.isfinite(Z08) & np.isfinite(Z21)
    print(f"grid {ny}x{nx} @ {a.res} m; {both.sum():,} cells with both epochs")
    print(f"before tie: median dz={np.nanmedian((Z21-Z08)[both]):+.3f}  "
          f"NMAD={_nmad((Z21-Z08)[both]):.3f} m")

    c = coreg.nuth_kaab(Z21, Z08, a.res)      # align 2008 onto 2021 (robust)
    print(f"tie 2008->2021: dx={c.dx:+.3f} dy={c.dy:+.3f} dz={c.dz:+.3f} m "
          f"(converged={c.converged}; dz not physical - absorbs datum)")
    diff = Z21 - (coreg._shift_grid(Z08, c.dx, c.dy, a.res) + c.dz)
    m = np.isfinite(diff); dd = diff[m]; nm = _nmad(dd)
    print(f"difference: median={np.median(dd):+.3f} NMAD={nm:.3f} std={dd.std():.3f} m")
    print(f"candidate-change cells |diff|>3*NMAD: {100*(np.abs(dd-np.median(dd))>3*nm).mean():.1f}%")

    if a.out:
        _write(diff, a.res, X0, Y0, ny, a.out)
    if a.figdir:
        _fig(Z21, diff, (X0, X1, Y0, Y1), nm, a.figdir)


def _write(diff, res, x0, y0, ny, out):
    import rasterio
    from rasterio.transform import from_origin
    Path(out).parent.mkdir(parents=True, exist_ok=True)
    with rasterio.open(out, "w", driver="GTiff", height=diff.shape[0],
                       width=diff.shape[1], count=1, dtype="float32",
                       crs="EPSG:26915", nodata=np.nan,
                       transform=from_origin(x0, y0 + ny * res, res, res)) as d:
        d.write(np.flipud(diff).astype("float32"), 1)
    print(f"wrote {out}")


def _fig(Z21, diff, ext, nm, figdir):
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LightSource
    Path(figdir).mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(1, 2, figsize=(13, 6))
    ls = LightSource(azdeg=315, altdeg=45)
    hs = ls.hillshade(np.nan_to_num(Z21, nan=np.nanmin(Z21)), vert_exag=2)
    ax[0].imshow(hs, extent=ext, cmap="gray", origin="lower")
    ax[0].set_title("2021 3DEP shaded relief")
    im = ax[1].imshow(diff, extent=ext, origin="lower", cmap="RdBu_r",
                      vmin=-0.5, vmax=0.5)
    ax[1].set_title("2021 - aligned 2008 (m)\nblue=deposition, red=erosion")
    fig.colorbar(im, ax=ax[1], shrink=0.7, label="elevation change (m)")
    for x in ax:
        x.set_xlabel("Easting (m)"); x.set_ylabel("Northing (m)")
    out = Path(figdir) / "cross_epoch_diff.png"
    fig.savefig(out, dpi=110, bbox_inches="tight"); plt.close(fig)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()

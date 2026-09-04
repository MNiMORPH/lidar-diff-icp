#!/usr/bin/env python3
"""What the floodplain cut actually REMOVES: the TPI mask against the elevation cut.

Using TPI to DETECT a floodplain is one thing; using it to DELETE cells from an analysis is
another, and it has to be looked at. TPI < -2 m over a 600-800 m window is a
topographic-position heuristic: what it removes depends on the window width relative to the
valley, it keeps flat terrace ground sitting at valley level, and on a wide valley the
interior reads as flat-stable rather than as floodplain. The elevation antimode cuts on the
quantity that actually defines a floodplain.

    env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python analysis/plot_floodplain_cuts.py \
        --tiles data/derived/elba data/derived/whitewater
"""
import argparse, json, os, re
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from lidar_diff_icp.refcells import reference_cells
from lidar_diff_icp.figures import grid_of
from lidar_diff_icp.viz import hillshade

ap = argparse.ArgumentParser()
ap.add_argument("--tiles", nargs="+", required=True)
ap.add_argument("--figdir", default="figures/method")
A = ap.parse_args()
os.makedirs(A.figdir, exist_ok=True)

fig, ax = plt.subplots(len(A.tiles), 3, figsize=(16.5, 6.2*len(A.tiles)), squeeze=False)
for row, tile in zip(ax, A.tiles):
    nm = os.path.basename(tile.rstrip("/"))
    X0, Y0, res, nx, ny = grid_of(tile)
    Z = np.load(f"{tile}/z_after.npy")
    ext = (X0, X0+nx*res, Y0, Y0+ny*res)
    hs = hillshade(Z, res, X0, Y0, fill_gaps=True)
    fld = np.load(f"{tile}/floodplain_mask.npy").astype(bool)
    _, rep = reference_cells(tile)                       # elevation cut, as shipped
    m = [re.search(r"antimode ([\d.]+) m", k) for k in rep]
    zthr = float([x.group(1) for x in m if x][0]) if any(m) else np.nan
    low = np.isfinite(Z) & (Z < zthr)
    both = fld & low; only_tpi = fld & ~low; only_z = low & ~fld
    print(f"{nm}: elevation antimode {zthr:.1f} m")
    print(f"   TPI mask {int(fld.sum()):,} cells;  below antimode {int(low.sum()):,};  "
          f"both {int(both.sum()):,};  TPI only {int(only_tpi.sum()):,};  "
          f"elevation only {int(only_z.sum()):,}")
    row[0].imshow(hs, extent=ext, origin="lower", cmap="gray")
    im = row[0].imshow(np.where(np.isfinite(Z), Z, np.nan), extent=ext, origin="lower",
                       cmap="terrain", alpha=0.55)
    row[0].contour(np.where(np.isfinite(Z), Z, np.nan), levels=[zthr], extent=ext,
                   colors="red", linewidths=1.6)
    row[0].set_title(f"{nm}: elevation, antimode {zthr:.1f} m (red)")
    fig.colorbar(im, ax=row[0], shrink=0.6, label="elevation (m)")
    ov = np.full(Z.shape, np.nan)
    ov[only_z] = 0; ov[both] = 1; ov[only_tpi] = 2
    row[1].imshow(hs, extent=ext, origin="lower", cmap="gray")
    # EXPLICIT colours. "brg" runs blue -> RED -> green, so the middle class rendered red
    # and the last green -- the reverse of the legend, and it inverted the reading of this
    # figure (Andy caught it 2026-09-04). Never index a categorical overlay into a named
    # continuous colormap.
    from matplotlib.colors import ListedColormap, BoundaryNorm
    cmap3 = ListedColormap(["tab:blue", "tab:green", "tab:red"])
    row[1].imshow(ov, extent=ext, origin="lower", cmap=cmap3,
                  norm=BoundaryNorm([-0.5, 0.5, 1.5, 2.5], cmap3.N), alpha=0.75)
    row[1].set_title(f"{nm}: what each cut removes\nblue = elevation only "
                     f"({int(only_z.sum()):,}), green = both ({int(both.sum()):,}), "
                     f"red = TPI only ({int(only_tpi.sum()):,})")
    zz = Z[np.isfinite(Z)]
    row[2].hist(zz, bins=200, color="0.75")
    row[2].hist(Z[fld & np.isfinite(Z)], bins=200, color="tab:red", alpha=0.65,
                label="removed by the TPI mask")
    row[2].axvline(zthr, color="k", lw=2, ls="--", label=f"antimode {zthr:.1f} m")
    row[2].set_xlabel("elevation (m)"); row[2].set_ylabel("cells")
    row[2].set_title("where each cut sits in the elevation histogram")
    row[2].legend(fontsize=8)
    for a in row[:2]:
        a.set_xlabel("Easting (m)")
    row[0].set_ylabel("Northing (m)")
fig.suptitle("Floodplain REMOVAL: the TPI mask against the elevation antimode", y=0.995)
fig.tight_layout(rect=(0, 0, 1, 0.985))
out = f"{A.figdir}/floodplain_cut_comparison.png"
fig.savefig(out, dpi=125, bbox_inches="tight"); plt.close(fig)
print("wrote", out)

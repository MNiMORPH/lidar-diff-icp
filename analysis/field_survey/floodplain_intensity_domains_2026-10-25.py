#!/usr/bin/env python
"""Two domains on the floodplain, and whether the southern patchwork is trees only (Andy, 2026-11-20).
Per 1 m cell over the field site (same extent as floodplain_intensity_grid_2026-10-25.py): all 2021 returns, their
height above the pipeline's 5 m 2021 ground grid (z_after.npy of tiles 4342-29-64 / 4358-29-01, the vendor class-2
median per cell, read bottom-up, nearest cell), mean intensity and its coefficient of variation, and the share of
returns more than 2 m up. A cell is TREE when that share exceeds 0.10 or its highest return is above 3 m (both MINE).
Non-tree cells are then split into a northern and a southern domain at N 4885480 (the boundary the 2008-offset
clustering gave; MINE) and restricted to the 10 m floodplain mask of the optical stack (nearest); their intensity and
CV distributions are compared, and non-tree cells are classed by intensity: dark < 2500, mid 2500-3800, bright > 3800
(cuts MINE, from the mean-intensity histogram's valleys reported in the log).
Outputs  data/field/2026-10-25/intensity/intensity_2021_1m_{tree,above2m_share,class}.tif, domains_log.txt
         figures/rtk_2026-10-25_floodplain_intensity_domains.png  class map; intensity and CV histograms of non-tree cells, north vs south
"""
import os, json
import numpy as np, pandas as pd, laspy, rasterio
from rasterio.transform import from_origin
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
ROOT = os.path.dirname(os.path.abspath(__file__)) + "/../.."; os.chdir(ROOT)
OUT = "data/field/2026-10-25/intensity"; LOG = open(f"{OUT}/domains_log.txt", "w")
def say(*a):
    s = " ".join(str(x) for x in a); print(s); LOG.write(s + "\n")
E0, E1, N0, N1 = 579780.0, 580160.0, 4885280.0, 4885760.0; RES = 1.0; nx, ny = int((E1 - E0) / RES), int((N1 - N0) / RES); NSPLIT = 4885480.0
grids = {}
for t in ("4342-29-64", "4358-29-01"):
    m = json.load(open(f"data/derived/{t}/meta.json")); grids[t] = (m["bounds"], m["res"], np.load(f"data/derived/{t}/z_after.npy", mmap_mode="r"))
def ground(x, y):
    z = np.full(len(x), np.nan)
    for t, (b, r, a) in grids.items():
        i = ((y - b[1]) / r).astype(int); j = ((x - b[0]) / r).astype(int); ok = (i >= 0) & (i < a.shape[0]) & (j >= 0) & (j < a.shape[1]) & np.isnan(z)
        if ok.any(): z[ok] = np.asarray(a)[i[ok], j[ok]]
    return z
acc = {k: np.zeros(ny * nx) for k in ("n", "s", "s2", "n2m")}; hmax = np.full(ny * nx, -np.inf)
for t in ("4342-29-64", "4358-29-01"):
    with laspy.open(f"data/after/valley/{t}_3dep_full.laz") as r:
        for ch in r.chunk_iterator(5_000_000):
            x = np.asarray(ch.x); y = np.asarray(ch.y); m = (x >= E0) & (x < E1) & (y >= N0) & (y < N1)
            if t == "4358-29-01": m &= x > 580030.0
            if not m.any(): continue
            x, y, z, inten = x[m], y[m], np.asarray(ch.z)[m], np.asarray(ch.intensity)[m].astype(float); h = z - ground(x, y)
            i = np.clip(((N1 - y) / RES).astype(int), 0, ny - 1); j = np.clip(((x - E0) / RES).astype(int), 0, nx - 1); cell = i * nx + j
            acc["n"] += np.bincount(cell, minlength=ny * nx); acc["s"] += np.bincount(cell, weights=inten, minlength=ny * nx); acc["s2"] += np.bincount(cell, weights=inten ** 2, minlength=ny * nx)
            ok = np.isfinite(h); acc["n2m"] += np.bincount(cell[ok & (h > 2)], minlength=ny * nx); np.maximum.at(hmax, cell[ok], h[ok])
n = acc["n"]
with np.errstate(invalid="ignore", divide="ignore"):
    mean = acc["s"] / n; cv = np.sqrt(np.maximum(acc["s2"] / n - mean ** 2, 0)) / mean; share2 = acc["n2m"] / n
hmax[~np.isfinite(hmax)] = np.nan; tree = (share2 > 0.10) | (hmax > 3.0); tree &= n > 0
z = np.load("data/field/2026-10-25/imagery/floodplain_model_features.npz", allow_pickle=True); xs, ys, fpm = z["xs"], z["ys"], z["fpm"].astype(bool)
cx = E0 + (np.arange(ny * nx) % nx + 0.5) * RES; cy = N1 - (np.arange(ny * nx) // nx + 0.5) * RES
fi = np.clip(np.round((ys[0] - cy) / 10).astype(int), 0, len(ys) - 1); fj = np.clip(np.round((cx - xs[0]) / 10).astype(int), 0, len(xs) - 1); fp = fpm[fi, fj]
cls = np.zeros(ny * nx, int); cls[n > 0] = 1; cls[(n > 0) & (mean < 2500)] = 2; cls[(n > 0) & (mean >= 2500) & (mean < 3800)] = 3; cls[(n > 0) & (mean >= 3800)] = 4; cls[tree] = 5   # 0 no data, 2 dark, 3 mid, 4 bright, 5 tree
say(f"cells {nx * ny}; with returns {int((n > 0).sum())}; tree cells {int(tree.sum())} ({100 * tree.sum() / (n > 0).sum():.1f} %); floodplain-mask cells {int(fp.sum())}")
hist, edges = np.histogram(mean[(n > 0) & ~tree & fp], bins=np.arange(0, 7001, 250)); say("mean-intensity histogram of non-tree floodplain cells (250-unit bins, counts): " + " ".join(f"{int(e)}:{c}" for e, c in zip(edges[:-1], hist)))
north = (cy >= NSPLIT) & fp & (n > 0); south = (cy < NSPLIT) & fp & (n > 0)
for name, dom in (("north (N >= 4885480)", north), ("south (N < 4885480)", south)):
    nt = dom & ~tree; q = lambda a, p: np.nanpercentile(a[nt], p)
    say(f"\n{name}: floodplain cells {int(dom.sum())}; tree {100 * tree[dom].mean():.1f} %; non-tree {int(nt.sum())}: intensity p10/p50/p90 {q(mean, 10):.0f}/{q(mean, 50):.0f}/{q(mean, 90):.0f}, CV p10/p50/p90 {q(cv, 10):.3f}/{q(cv, 50):.3f}/{q(cv, 90):.3f}; "
        f"classes of non-tree cells: dark {100 * (cls[nt] == 2).mean():.1f} %, mid {100 * (cls[nt] == 3).mean():.1f} %, bright {100 * (cls[nt] == 4).mean():.1f} %; share of returns above 2 m in non-tree cells median {np.nanmedian(share2[nt]):.3f}; highest return in non-tree cells median {np.nanmedian(hmax[nt]):.2f} m, p90 {np.nanpercentile(hmax[nt], 90):.2f}")
    say(f"   non-tree cells with highest return between 0.5 and 3 m (shrub or tall herb layer, not tree): {100 * ((hmax[nt] > 0.5) & (hmax[nt] <= 3)).mean():.1f} %")

# the same comparison with a 2 m buffer around tree cells removed (edge mixing inside 1 m cells; buffer MINE)
from scipy.ndimage import binary_dilation
near_tree = binary_dilation(tree.reshape(ny, nx), iterations=2).ravel()
say("\nnon-tree cells more than 2 m from any tree cell:")
for name, dom in (("north", north), ("south", south)):
    nt = dom & ~near_tree; q = lambda a, p: np.nanpercentile(a[nt], p)
    say(f"   {name}: n {int(nt.sum())}; intensity p10/p50/p90 {q(mean, 10):.0f}/{q(mean, 50):.0f}/{q(mean, 90):.0f}; CV p10/p50/p90 {q(cv, 10):.3f}/{q(cv, 50):.3f}/{q(cv, 90):.3f}; dark {100 * (cls[nt] == 2).mean():.1f} %, mid {100 * (cls[nt] == 3).mean():.1f} %, bright {100 * (cls[nt] == 4).mean():.1f} %; top return 0.5-3 m {100 * ((hmax[nt] > 0.5) & (hmax[nt] <= 3)).mean():.1f} %; top return median {np.nanmedian(hmax[nt]):.2f} m, p90 {np.nanpercentile(hmax[nt], 90):.2f}")
tr = from_origin(E0, N1, RES, RES); prof = dict(driver="GTiff", height=ny, width=nx, count=1, crs="EPSG:6344", transform=tr, compress="deflate")
for name, a, dt, nod in (("tree", tree.astype(np.uint8), "uint8", 255), ("above2m_share", np.where(n > 0, share2, -9999.0), "float32", -9999.0), ("class", cls.astype(np.uint8), "uint8", 0)):
    with rasterio.open(f"{OUT}/intensity_2021_1m_{name}.tif", "w", dtype=dt, nodata=nod, **prof) as dst: dst.write(np.asarray(a).reshape(ny, nx).astype(dt), 1)
say("wrote tree, above2m_share and class GeoTIFFs")
P = pd.read_csv("data/field/2026-10-25/rtk_points_all.csv"); F = P[P.group.str.startswith("floodplain") & (P.status == "FIX")]
fig = plt.figure(figsize=(18, 8)); ax = fig.add_subplot(1, 2, 1); cm = ListedColormap(["white", "0.85", "#2c7bb6", "#fdae61", "#ffffbf", "#1a9641"])
img = cls.reshape(ny, nx).astype(float); img[~fp.reshape(ny, nx) & (cls.reshape(ny, nx) != 5)] = 1
ax.imshow(img, extent=(E0, E1, N0, N1), cmap=cm, vmin=-0.5, vmax=5.5, interpolation="nearest"); ax.plot(F.easting_tied, F.northing_tied, ".", color="k", ms=2); ax.axhline(NSPLIT, color="k", lw=0.8, ls="--"); ax.set_aspect("equal"); ax.tick_params(labelsize=7)
from matplotlib.patches import Patch; ax.legend(handles=[Patch(color="#2c7bb6", label="dark non-tree (< 2500)"), Patch(color="#fdae61", label="mid non-tree (2500-3800)"), Patch(color="#ffffbf", label="bright non-tree (> 3800)"), Patch(color="#1a9641", label="tree (> 10 % of returns above 2 m, or top > 3 m)"), Patch(color="0.85", label="outside the floodplain mask")], fontsize=7, loc="lower left")
for k, (a, lab, bins) in enumerate(((mean, "mean intensity, non-tree floodplain cells", np.arange(0, 7001, 100)), (cv, "coefficient of variation, non-tree floodplain cells", np.arange(0, 1.0, 0.01)))):
    axh = fig.add_subplot(2, 2, 2 + 2 * k)
    for dom, c, lab2 in ((north, "#0072b2", "north"), (south, "#d55e00", "south")):
        axh.hist(a[dom & ~tree], bins=bins, density=True, histtype="step", color=c, lw=1.5, label=f"{lab2} (n {int((dom & ~tree).sum())})")
    axh.set_xlabel(lab); axh.set_ylabel("density"); axh.legend(frameon=False, fontsize=8); axh.spines[["top", "right"]].set_visible(False)
fig.tight_layout(); fig.savefig("figures/rtk_2026-10-25_floodplain_intensity_domains.png", dpi=int(os.environ.get("FIGDPI", 170))); say("wrote figures/rtk_2026-10-25_floodplain_intensity_domains.png")
LOG.close()

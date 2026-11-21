#!/usr/bin/env python
"""The forested hillslope: histogram of 2021 return height and brightness, slope-normal (Andy, 2026-11-21).
Area: the bluff west of the road at the hillslope traverses, E 579620-579900, N 4885260-4885540 (MINE), cells of the
official USGS 1 m DEM (tiles x57y489/x58y489, fetched 2026-11-21 into dem3dep/usgs_1m_dem_bluff.tif) with slope >= 10 deg
(MINE) and classed tree by the earlier rule where the 1 m grids exist; here simply slope >= 10 deg west of the road line.
Height of each return = (z - DEM at the return's own x, y by bilinear interpolation) x cos(slope at x, y): the distance
perpendicular to the local DEM surface, so that a return low in a sloping cell is not counted as low. Slope and aspect
from the DEM by central differences at 1 m. Returns: the 2021 cloud of tile 4342-29-64 (all returns). Shown: the full
range on a symlog axis with log density and intensity-colored bars; a near-ground zoom (-0.5..+1.0 m); the near-ground
peak width by slope class (10-20, 20-30, > 30 deg) as the check of the slope-normal referencing; the floodplain sward
region's curve for comparison; the hillslope RTK FIX points against the DEM.
Outputs  data/field/2026-10-25/intensity/hillslope_returns_log.txt; figures/rtk_2026-10-25_hillslope_returns.png
"""
import os
import numpy as np, pandas as pd, laspy, rasterio
from rasterio.windows import from_bounds
from rasterio.merge import merge
from rasterio.io import MemoryFile
from scipy.ndimage import map_coordinates
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
ROOT = os.path.dirname(os.path.abspath(__file__)) + "/../.."; os.chdir(ROOT)
OUT = "data/field/2026-10-25/intensity"; LOG = open(f"{OUT}/hillslope_returns_log.txt", "w")
def say(*a):
    s = " ".join(str(x) for x in a); print(s); LOG.write(s + "\n")
E0, E1, N0, N1 = 579620.0, 579900.0, 4885260.0, 4885540.0; RES = 1.0; nx, ny = int(E1 - E0), int(N1 - N0)
DEMF = "data/field/2026-10-25/dem3dep/usgs_1m_dem_bluff.tif"
if not os.path.exists(DEMF):
    urls = ["https://prd-tnm.s3.amazonaws.com/StagedProducts/Elevation/1m/Projects/MN_SE_Driftless_2021_B21/TIFF/USGS_1M_15_x57y489_MN_SE_Driftless_2021_B21.tif"]
    srcs = []
    for u in urls:
        r = rasterio.open("/vsicurl/" + u); w = from_bounds(E0, N0, E1, N1, transform=r.transform); a = r.read(1, window=w); t = r.window_transform(w)
        m = MemoryFile(); d = m.open(driver="GTiff", height=a.shape[0], width=a.shape[1], count=1, dtype=a.dtype, crs=r.crs, transform=t, nodata=r.nodata); d.write(a, 1); d.close(); srcs.append(m.open())
    mos, tr = merge(srcs)
    with rasterio.open(DEMF, "w", driver="GTiff", height=mos.shape[1], width=mos.shape[2], count=1, dtype="float32", crs=srcs[0].crs, transform=tr, nodata=srcs[0].nodata, compress="deflate") as dst: dst.write(mos[0].astype("float32"), 1)
with rasterio.open(DEMF) as r: DEM = r.read(1).astype(float); DEM[DEM == r.nodata] = np.nan; tr = r.transform
assert abs(tr.c - E0) < 1e-6 and abs(tr.f - N1) < 1e-6, (tr.c, tr.f)
gy, gx = np.gradient(DEM, RES); gy = -gy   # rows run north to south
SLOPE = np.degrees(np.arctan(np.hypot(gx, gy)))
def at(a, x, y): return map_coordinates(a, [((N1 - y) / RES - 0.5), ((x - E0) / RES - 0.5)], order=1, mode="nearest")
# road line (from the intensity map): keep cells west of it by 6 m
road = np.array([(579893, 4885760), (579900, 4885700), (579912, 4885640), (579925, 4885580), (579937, 4885520), (579947, 4885460), (579955, 4885400), (579962, 4885340), (579968, 4885280)], float)
cy = N1 - (np.arange(ny) + 0.5) * RES; cx = E0 + (np.arange(nx) + 0.5) * RES; CX, CY = np.meshgrid(cx, cy)
roadE = np.interp(CY, road[::-1, 1], road[::-1, 0]); BLUFF = (CX < roadE - 6) & (SLOPE >= 10) & np.isfinite(DEM)
say(f"bluff cells: {int(BLUFF.sum())} of {nx * ny}; slope median {np.median(SLOPE[BLUFF]):.1f} deg, p10 {np.percentile(SLOPE[BLUFF], 10):.1f}, p90 {np.percentile(SLOPE[BLUFF], 90):.1f}")
parts = []
with laspy.open("data/after/valley/4342-29-64_3dep_full.laz") as r:
    for ch in r.chunk_iterator(5_000_000):
        x = np.asarray(ch.x); y = np.asarray(ch.y); m = (x >= E0 + 1) & (x < E1 - 1) & (y >= N0 + 1) & (y < N1 - 1)
        if not m.any(): continue
        i = ((N1 - y[m]) / RES).astype(int); j = ((x[m] - E0) / RES).astype(int); k = BLUFF[i, j]
        if not k.any(): continue
        xx, yy, zz = x[m][k], y[m][k], np.asarray(ch.z)[m][k]; d0 = at(DEM, xx, yy); s0 = np.radians(at(SLOPE, xx, yy))
        parts.append(pd.DataFrame(dict(h=(zz - d0) * np.cos(s0), hv=zz - d0, slope=np.degrees(s0), inten=np.asarray(ch.intensity)[m][k].astype(float), nr=np.asarray(ch.number_of_returns)[m][k], rn=np.asarray(ch.return_number)[m][k], c=np.asarray(ch.classification)[m][k])))
R = pd.concat(parts, ignore_index=True); R = R[np.isfinite(R.h)]
say(f"returns on the bluff: {len(R):,}; single-return share {100 * (R.nr == 1).mean():.1f} %; vendor classes {R.c.value_counts().to_dict()}")
def nmad(v): return 1.4826 * np.median(np.abs(v - np.median(v)))
near = R[(R.h > -0.5) & (R.h < 0.5)]
say(f"slope-normal height: median {R.h.median():+.3f}, share within +-0.05 m {100 * (R.h.abs() < 0.05).mean():.1f} %, below -0.05 {100 * (R.h < -0.05).mean():.2f} %, 0.05-0.5 m {100 * ((R.h > 0.05) & (R.h < 0.5)).mean():.1f} %, 0.5-2 m {100 * ((R.h >= 0.5) & (R.h < 2)).mean():.1f} %, above 2 m {100 * (R.h >= 2).mean():.1f} %")
say(f"near-ground returns (|h| < 0.5 m): {len(near):,}; peak NMAD about the median {nmad(near.h):.3f} m; vertical (not slope-normal) NMAD for comparison {nmad(near.hv):.3f} m")
say("intensity medians: ground peak |h| < 0.05: " + f"{R.inten[R.h.abs() < 0.05].median():.0f}; 0.05-0.5 m: {R.inten[(R.h > 0.05) & (R.h < 0.5)].median():.0f}; 0.5-2 m: {R.inten[(R.h >= 0.5) & (R.h < 2)].median():.0f}; canopy > 2 m: {R.inten[R.h >= 2].median():.0f}; below -0.05: {R.inten[R.h < -0.05].median():.0f}")
say("ground-peak returns by return type: single " + f"{100 * ((R.h.abs() < 0.05) & (R.nr == 1)).sum() / (R.h.abs() < 0.05).sum():.0f} %, last-of-several {100 * ((R.h.abs() < 0.05) & (R.nr > 1) & (R.rn == R.nr)).sum() / (R.h.abs() < 0.05).sum():.0f} %")
say("\ncheck of the slope-normal referencing: near-ground peak width (NMAD of |h| < 0.5 m returns) by slope class, slope-normal vs vertical:")
for lo, hi in ((10, 20), (20, 30), (30, 60)):
    q = near[(near.slope >= lo) & (near.slope < hi)]; say(f"   slope {lo}-{hi} deg: n {len(q):7d}  slope-normal NMAD {nmad(q.h):.3f} m, median {q.h.median():+.3f}  |  vertical NMAD {nmad(q.hv):.3f} m, median {q.hv.median():+.3f}")
# RTK hillslope points against the DEM
P = pd.read_csv("data/field/2026-10-25/rtk_points_all.csv"); H = P[P.group.str.startswith("hillslope") & (P.status == "FIX")].copy()
H = H[(H.easting_tied > E0 + 1) & (H.easting_tied < E1 - 1) & (H.northing_tied > N0 + 1) & (H.northing_tied < N1 - 1)]
H["dem_minus_rtk_v"] = at(DEM, H.easting_tied.values, H.northing_tied.values) - H.navd88_tied_m.values; H["slope"] = at(SLOPE, H.easting_tied.values, H.northing_tied.values); H["dem_minus_rtk_n"] = H.dem_minus_rtk_v * np.cos(np.radians(H.slope))
say(f"\nhillslope RTK FIX points in the area: {len(H)}; official DEM minus RTK, vertical: median {H.dem_minus_rtk_v.median():+.3f} NMAD {nmad(H.dem_minus_rtk_v):.3f}; slope-normal: median {H.dem_minus_rtk_n.median():+.3f} NMAD {nmad(H.dem_minus_rtk_n):.3f}; by group " + ", ".join(f"{g} {v:+.3f} (n {len(H[H.group == g])})" for g, v in H.groupby("group").dem_minus_rtk_n.median().items()))
# floodplain sward comparison curve (region returns minus DEM, from the earlier histogram csv)
F = pd.read_csv(f"{OUT}/region_vs_dem_hist.csv")
norm = Normalize(vmin=1000, vmax=4500); cmap = plt.get_cmap("viridis")
fig, axs = plt.subplots(1, 3, figsize=(20, 5.5))
edges = np.concatenate([np.arange(-0.60, 0.30, 0.01), 0.30 * 1.15 ** np.arange(0, 40)]); edges = edges[edges <= R.h.max() * 1.15]; mids = 0.5 * (edges[1:] + edges[:-1])
cnt = np.histogram(R.h, bins=edges)[0]; dens = cnt / (len(R) * np.diff(edges)); b = np.digitize(R.h, edges) - 1; med = np.array([R.inten[b == k].median() if (b == k).sum() >= 20 else np.nan for k in range(len(mids))])
ax = axs[0]; ax.bar(mids, np.where(dens > 0, dens, np.nan), width=np.diff(edges), color=[cmap(norm(m)) if np.isfinite(m) else "0.85" for m in med], edgecolor="none", align="center"); ax.set_xscale("symlog", linthresh=0.3, linscale=2.5); ax.set_yscale("log"); ax.set_xlim(-0.6, edges[-1])
ax.set_xticks([-0.5, -0.2, 0, 0.2, 1, 3, 10, 30]); ax.set_xticklabels(["-0.5", "-0.2", "0", "0.2", "1", "3", "10", "30"]); ax.axvline(0, color="0.4", lw=0.8, ls=":"); ax.set_xlabel("2021 return, slope-normal distance from the official DEM (m)"); ax.set_ylabel("density (1/m), log"); ax.set_title(f"forested bluff, all returns (n {len(R):,})", fontsize=10)
ax = axs[1]; e2 = np.arange(-0.5, 1.0 + 0.005, 0.01); m2 = 0.5 * (e2[1:] + e2[:-1]); c2 = np.histogram(R.h, bins=e2)[0]; d2 = c2 / (len(R) * 0.01); b2 = np.digitize(R.h, e2) - 1; md2 = np.array([R.inten[b2 == k].median() if (b2 == k).sum() >= 20 else np.nan for k in range(len(m2))])
ax.bar(m2, np.where(d2 > 0, d2, np.nan), width=0.01, color=[cmap(norm(m)) if np.isfinite(m) else "0.85" for m in md2], edgecolor="none"); ax.set_yscale("log"); ax.axvline(0, color="0.4", lw=0.8, ls=":"); ax.set_xlim(-0.5, 1.0)
ax.plot(F.height_minus_dem_m, F.density, color="#d55e00", lw=1.3, label="floodplain sward region, returns minus DEM"); ax.legend(frameon=False, fontsize=8); ax.set_xlabel("slope-normal distance from the official DEM (m)"); ax.set_ylabel("density (1/m), log"); ax.set_title("near the ground", fontsize=10)
ax = axs[2]
for (lo, hi), c in (((10, 20), "#0072b2"), ((20, 30), "#e69f00"), ((30, 60), "#d55e00")):
    q = R[(R.slope >= lo) & (R.slope < hi)]; ax.hist(q.h, bins=e2, density=True, histtype="step", color=c, lw=1.4, label=f"slope {lo}-{hi} deg, slope-normal (n {len(q):,})"); ax.hist(q.hv, bins=e2, density=True, histtype="step", color=c, lw=0.8, ls=":", label=f"same, vertical")
ax.set_yscale("log"); ax.axvline(0, color="0.4", lw=0.8, ls=":"); ax.set_xlim(-0.5, 1.0); ax.set_xlabel("distance from the official DEM (m)"); ax.set_ylabel("density (1/m), log"); ax.legend(frameon=False, fontsize=7); ax.set_title("by slope class: solid slope-normal, dotted vertical", fontsize=10)
sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap); sm.set_array([]); fig.colorbar(sm, ax=axs[:2], shrink=0.8, label="median intensity in the bar")
for a in axs: a.spines[["top", "right"]].set_visible(False)
fig.savefig("figures/rtk_2026-10-25_hillslope_returns.png", dpi=int(os.environ.get("FIGDPI", 160)), bbox_inches="tight"); say("wrote figures/rtk_2026-10-25_hillslope_returns.png")

# full-height version (Andy, 2026-11-21): symlog height axis through the canopy top, more room for the log part, plus a linear canopy profile
say(f"canopy: returns above 2 m {100 * (R.h >= 2).mean():.1f} %; their height median {R.h[R.h >= 2].median():.1f} m, p90 {R.h[R.h >= 2].quantile(.9):.1f}, p99 {R.h[R.h >= 2].quantile(.99):.1f}, max {R.h.max():.1f} m")
fig, axs = plt.subplots(1, 2, figsize=(19, 5.5), gridspec_kw=dict(width_ratios=[1.4, 1]))
ax = axs[0]; ax.bar(mids, np.where(dens > 0, dens, np.nan), width=np.diff(edges), color=[cmap(norm(m)) if np.isfinite(m) else "0.85" for m in med], edgecolor="none", align="center")
ax.set_xscale("symlog", linthresh=0.3, linscale=1.0); ax.set_yscale("log"); ax.set_xlim(-0.6, edges[-1]); ax.set_xticks([-0.5, -0.2, 0, 0.2, 0.5, 1, 2, 5, 10, 20, 30]); ax.set_xticklabels(["-0.5", "-0.2", "0", "0.2", "0.5", "1", "2", "5", "10", "20", "30"])
ax.axvline(0, color="0.4", lw=0.8, ls=":"); ax.set_xlabel("2021 return, slope-normal distance from the official DEM (m); symmetric-log axis, linear within 0.3 m of zero"); ax.set_ylabel("density (1/m), log"); ax.set_title(f"forested bluff, ground to canopy top (n {len(R):,})", fontsize=10)
ax = axs[1]; eL = np.arange(0, 32.5, 0.5); mL = 0.5 * (eL[1:] + eL[:-1]); cL = np.histogram(R.h, bins=eL)[0]; dL = cL / (len(R) * 0.5); bL = np.digitize(R.h, eL) - 1; mdL = np.array([R.inten[bL == k].median() if (bL == k).sum() >= 20 else np.nan for k in range(len(mL))])
ax.barh(mL, dL, height=0.5, color=[cmap(norm(m)) if np.isfinite(m) else "0.85" for m in mdL], edgecolor="none"); ax.set_xscale("log"); ax.set_ylim(0, 32); ax.set_xlabel("density (1/m), log"); ax.set_ylabel("slope-normal height above the DEM (m), linear, 0.5 m bins"); ax.set_title("the same record as a vertical profile", fontsize=10)
sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap); sm.set_array([]); fig.colorbar(sm, ax=axs, shrink=0.8, label="median intensity in the bar")
for a in axs: a.spines[["top", "right"]].set_visible(False)
fig.savefig("figures/rtk_2026-10-25_hillslope_returns_canopy.png", dpi=int(os.environ.get("FIGDPI", 160)), bbox_inches="tight"); say("wrote figures/rtk_2026-10-25_hillslope_returns_canopy.png")
LOG.close()

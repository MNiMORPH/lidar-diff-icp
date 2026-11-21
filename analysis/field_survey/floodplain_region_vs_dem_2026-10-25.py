#!/usr/bin/env python
"""A contiguous sward region picked by lidar intensity and CV around the walked loops, and the PDF of ALL its 2021
returns relative to the official USGS 1 m 3DEP DEM (Andy, 2026-11-21: with a consistent offset between mapped and true
ground, the all-data PDF against the official DEM should be our RTK-referenced PDF translated upward).
Region (cuts MINE): seed = the 1 m cells holding the RS+ B and RS+ A FIX points classed low or smooth sward; keep cells
whose mean intensity lies within the seeds' 5th-95th percentile and whose CV is below the seeds' 95th percentile;
close 1-cell holes; take the connected component holding the most seeds.
Official DEM: dem3dep/usgs_1m_dem_field.tif, the USGS 1 m DEM tiles x57y489 and x58y489 of MN_SE_Driftless_2021_B21
windowed to the field (fetched 2026-11-21 from prd-tnm.s3.amazonaws.com; EPSG:26915 label, NAVD88). At the 147
floodplain FIX points it sits +0.130 m above the tied RTK ground (NMAD 0.031) and -0.005 m from our 5 m grid.
Returns: both 2021 clouds, overlap from the first tile; dz_dem = z minus the DEM value of the return's 1 m cell.
Outputs  data/field/2026-10-25/intensity/region_sward_mask.tif; region_vs_dem_hist.csv; region_vs_dem_log.txt
         figures/rtk_2026-10-25_floodplain_region_vs_dem.png  region map; pooled PDFs (region returns minus DEM; RTK-point
                                                               returns minus RTK, same bins) with intensity-colored bars
"""
import os
import numpy as np, pandas as pd, laspy, rasterio
from rasterio.transform import from_origin
from scipy.ndimage import binary_closing, label
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
ROOT = os.path.dirname(os.path.abspath(__file__)) + "/../.."; os.chdir(ROOT)
OUT = "data/field/2026-10-25/intensity"; LOG = open(f"{OUT}/region_vs_dem_log.txt", "w")
def say(*a):
    s = " ".join(str(x) for x in a); print(s); LOG.write(s + "\n")
E0, E1, N0, N1 = 579780.0, 580160.0, 4885280.0, 4885760.0; RES = 1.0; nx, ny = int((E1 - E0) / RES), int((N1 - N0) / RES)
with rasterio.open(f"{OUT}/intensity_2021_1m_mean.tif") as r: MI = r.read(1).astype(float); MI[MI == r.nodata] = np.nan
with rasterio.open(f"{OUT}/intensity_2021_1m_cv.tif") as r: CV = r.read(1).astype(float); CV[CV == r.nodata] = np.nan
with rasterio.open(f"{OUT}/cell_profiles_cluster.tif") as r: CL = r.read(1).astype(float); CL[CL == 255] = np.nan
with rasterio.open("data/field/2026-10-25/dem3dep/usgs_1m_dem_field.tif") as r:
    assert abs(r.transform.c - E0) < 1e-6 and abs(r.transform.f - N1) < 1e-6 and r.width == nx and r.height == ny, "DEM window must match the 1 m grid"; DEM = r.read(1).astype(float); DEM[DEM == r.nodata] = np.nan
P = pd.read_csv("data/field/2026-10-25/rtk_points_all.csv"); F = P[P.group.isin(["floodplain_RSplus_B", "floodplain_RSplus_A"]) & (P.status == "FIX")]
fi = np.clip(((N1 - F.northing_tied.values) / RES).astype(int), 0, ny - 1); fj = np.clip(((F.easting_tied.values - E0) / RES).astype(int), 0, nx - 1)
seed_ok = np.isin(CL[fi, fj], [0, 1]); si, sj = fi[seed_ok], fj[seed_ok]
ilo, ihi = np.nanpercentile(MI[si, sj], [5, 95]); cvhi = np.nanpercentile(CV[si, sj], 95)
say(f"seeds: {int(seed_ok.sum())} of {len(F)} RS+ A/B FIX point cells classed low or smooth sward; intensity window {ilo:.0f}-{ihi:.0f}, CV below {cvhi:.3f} (5th/95th percentiles of the seeds)")
cand = (MI >= ilo) & (MI <= ihi) & (CV <= cvhi); cand = binary_closing(cand, iterations=1) | cand
lab, nlab = label(cand); counts = np.bincount(lab[si, sj][lab[si, sj] > 0], minlength=nlab + 1); k = int(np.argmax(counts)); REG = lab == k
say(f"candidate cells {int(cand.sum())}; connected components {nlab}; chosen component holds {counts[k]} seeds, {int(REG.sum())} cells = {REG.sum() * RES * RES / 1e4:.2f} ha; its extent E {E0 + np.where(REG.any(0))[0].min():.0f}-{E0 + np.where(REG.any(0))[0].max():.0f}, N {N1 - np.where(REG.any(1))[0].max():.0f}-{N1 - np.where(REG.any(1))[0].min():.0f}")
say(f"region: mean intensity median {np.nanmedian(MI[REG]):.0f}, CV median {np.nanmedian(CV[REG]):.3f}; structure classes inside: " + str({int(c): int(n) for c, n in zip(*np.unique(CL[REG][np.isfinite(CL[REG])], return_counts=True))}) + " (0 smooth sward, 1 low sward, 2 tall herb, 3 under crowns, 4 shrub, 5 tree)")
with rasterio.open(f"{OUT}/region_sward_mask.tif", "w", driver="GTiff", height=ny, width=nx, count=1, dtype="uint8", nodata=255, crs="EPSG:6344", transform=from_origin(E0, N1, RES, RES), compress="deflate") as dst: dst.write(REG.astype("uint8"), 1)
# returns in the region against the DEM
parts = []
for t in ("4342-29-64", "4358-29-01"):
    with laspy.open(f"data/after/valley/{t}_3dep_full.laz") as r:
        for ch in r.chunk_iterator(5_000_000):
            x = np.asarray(ch.x); y = np.asarray(ch.y); m = (x >= E0) & (x < E1) & (y >= N0) & (y < N1)
            if t == "4358-29-01": m &= x > 580030.0
            if not m.any(): continue
            i = np.clip(((N1 - y[m]) / RES).astype(int), 0, ny - 1); j = np.clip(((x[m] - E0) / RES).astype(int), 0, nx - 1); inreg = REG[i, j]
            if inreg.any(): parts.append(pd.DataFrame(dict(dz=np.asarray(ch.z)[m][inreg] - DEM[i[inreg], j[inreg]], inten=np.asarray(ch.intensity)[m][inreg].astype(float), nr=np.asarray(ch.number_of_returns)[m][inreg], cell=(i[inreg] * nx + j[inreg]))))
R = pd.concat(parts, ignore_index=True); R = R[np.isfinite(R.dz)]
say(f"\nreturns in the region: {len(R):,} over {R.cell.nunique():,} cells ({len(R) / R.cell.nunique():.1f} per cell); single-return share {100 * (R.nr == 1).mean():.1f} %")
say(f"return minus official DEM: median {R.dz.median():+.3f}, p05 {R.dz.quantile(.05):+.3f}, p10 {R.dz.quantile(.1):+.3f}, p90 {R.dz.quantile(.9):+.3f}, p95 {R.dz.quantile(.95):+.3f}, NMAD {1.4826 * np.median(np.abs(R.dz - R.dz.median())):.3f}; share below -0.08 m (the shelf, 0.13 - 0.05) {100 * (R.dz < -0.08).mean():.2f} %, below -0.11 {100 * (R.dz < -0.11).mean():.2f} %, above +0.20 {100 * (R.dz > 0.20).mean():.2f} %")
# the RTK-referenced comparison set: RS+ B and RS+ A FIX points' 1 m returns, dz vs RTK, from the per-return tables, restricted to points inside the region
parts = []
for tag, g in (("rsplus_b", "floodplain_RSplus_B"), ("rsplus_a", "floodplain_RSplus_A")):
    t = pd.read_csv(f"data/field/2026-10-25/returns/{tag}_returns_2021.csv.gz"); pts = F[F.group == g]
    inside = pts[REG[np.clip(((N1 - pts.northing_tied.values) / RES).astype(int), 0, ny - 1), np.clip(((pts.easting_tied.values - E0) / RES).astype(int), 0, nx - 1)]].point
    parts.append(t[t.radius_1m & t.point.isin(inside)][["dz_m", "intensity"]])
Q = pd.concat(parts, ignore_index=True); off = np.nanmedian(DEM[fi, fj] - F.navd88_tied_m.values)
say(f"RTK-referenced set: {len(Q)} returns within 1 m of the FIX points inside the region; minus RTK median {Q.dz_m.median():+.3f}, p05 {Q.dz_m.quantile(.05):+.3f}, p95 {Q.dz_m.quantile(.95):+.3f}; DEM minus RTK at the RS+ A/B points median {off:+.3f}")
say(f"translation test: region returns minus DEM, shifted by +{off:.3f}: median {R.dz.median() + off:+.3f} vs RTK set {Q.dz_m.median():+.3f}; p05 {R.dz.quantile(.05) + off:+.3f} vs {Q.dz_m.quantile(.05):+.3f}; p95 {R.dz.quantile(.95) + off:+.3f} vs {Q.dz_m.quantile(.95):+.3f}; shelf share (region below -0.08 vs RTK set below +0.05): {100 * (R.dz < -0.08).mean():.2f} % vs {100 * (Q.dz_m < 0.05).mean():.2f} %")
# per-cell medians
cm = R.groupby("cell").dz.median(); say(f"per-cell median of return minus DEM: median {cm.median():+.3f}, p10 {cm.quantile(.1):+.3f}, p90 {cm.quantile(.9):+.3f} over {len(cm)} cells")
edges = np.arange(-0.30, 0.40 + 0.005, 0.01); mids = 0.5 * (edges[1:] + edges[:-1])
h = np.histogram(R.dz, bins=edges)[0] / (len(R) * 0.01); b = np.digitize(R.dz, edges) - 1; med = np.array([R.inten[b == k].median() if (b == k).sum() >= 20 else np.nan for k in range(len(mids))])
pd.DataFrame(dict(height_minus_dem_m=mids, density=h, median_intensity=med, n=np.histogram(R.dz, bins=edges)[0])).to_csv(f"{OUT}/region_vs_dem_hist.csv", index=False, float_format="%.4f")
fig = plt.figure(figsize=(19, 7)); norm = Normalize(vmin=2500, vmax=5000); cmap = plt.get_cmap("viridis")
ax = fig.add_subplot(1, 3, 1); ax.imshow(MI, extent=(E0, E1, N0, N1), cmap="gray", vmin=1500, vmax=4800); ov = np.where(REG, 1.0, np.nan); ax.imshow(ov, extent=(E0, E1, N0, N1), cmap=matplotlib.colors.ListedColormap(["#d55e00"]), alpha=0.35, interpolation="nearest")
ax.plot(F.easting_tied, F.northing_tied, ".", color="#0072b2", ms=3); ax.set_xlim(579850, 580110); ax.set_ylim(4885320, 4885700); ax.set_aspect("equal"); ax.tick_params(labelsize=7); ax.set_title(f"sward region by intensity and CV: {REG.sum() / 1e4:.2f} ha", fontsize=10)
ax = fig.add_subplot(1, 3, 2); ax.bar(mids, h, width=0.01, color=[cmap(norm(m)) if np.isfinite(m) else "0.85" for m in med], edgecolor="none"); ax.axvline(0, color="0.4", lw=0.8, ls=":"); ax.set_xlim(-0.30, 0.40)
ax.set_xlabel("2021 return minus official USGS 1 m DEM (m)"); ax.set_ylabel("probability density (1/m)"); ax.set_title(f"all returns in the region (n {len(R):,})", fontsize=10); ax.spines[["top", "right"]].set_visible(False)
ax = fig.add_subplot(1, 3, 3); e2 = np.arange(-0.30 + off, 0.40 + off + 0.005, 0.01)
ax.hist(R.dz + off, bins=e2, density=True, histtype="step", color="#d55e00", lw=1.6, label=f"region returns minus DEM, shifted +{off:.3f} m")
ax.hist(Q.dz_m, bins=e2, density=True, histtype="step", color="#0072b2", lw=1.6, label=f"returns within 1 m of the RTK points, minus RTK (n {len(Q)})")
ax.axvline(0, color="0.4", lw=0.8, ls=":"); ax.set_xlim(-0.30 + off, 0.40 + off); ax.set_xlabel("height above the surveyed ground (m)"); ax.set_ylabel("probability density (1/m)"); ax.legend(frameon=False, fontsize=8); ax.spines[["top", "right"]].set_visible(False)
sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap); sm.set_array([]); fig.colorbar(sm, ax=fig.axes[1], shrink=0.7, label="median intensity in the bar")
fig.tight_layout(); fig.savefig("figures/rtk_2026-10-25_floodplain_region_vs_dem.png", dpi=int(os.environ.get("FIGDPI", 160))); say("wrote figures/rtk_2026-10-25_floodplain_region_vs_dem.png")

# log-density version (Andy, 2026-11-21): wider range so the canopy tail and the shelf both show
eL = np.arange(-0.50, 1.00 + 0.005, 0.01); mL = 0.5 * (eL[1:] + eL[:-1]); hL = np.histogram(R.dz, bins=eL)[0] / (len(R) * 0.01); bL = np.digitize(R.dz, eL) - 1
medL = np.array([R.inten[bL == k].median() if (bL == k).sum() >= 20 else np.nan for k in range(len(mL))])
fig, axs = plt.subplots(1, 2, figsize=(15, 5.2))
ax = axs[0]; ax.bar(mL, np.where(hL > 0, hL, np.nan), width=0.01, color=[cmap(norm(m)) if np.isfinite(m) else "0.85" for m in medL], edgecolor="none"); ax.set_yscale("log"); ax.axvline(0, color="0.4", lw=0.8, ls=":")
ax.set_xlim(-0.50, 1.00); ax.set_xlabel("2021 return minus official USGS 1 m DEM (m)"); ax.set_ylabel("probability density (1/m), log"); ax.set_title(f"all returns in the region (n {len(R):,}); {100 * ((R.dz < -0.5) | (R.dz > 1.0)).mean():.2f} % outside the axis", fontsize=10); ax.spines[["top", "right"]].set_visible(False)
ax = axs[1]; e3 = np.arange(-0.50 + off, 1.00 + off + 0.005, 0.01)
ax.hist(R.dz + off, bins=e3, density=True, histtype="step", color="#d55e00", lw=1.6, label=f"region returns minus DEM, shifted +{off:.3f} m")
ax.hist(Q.dz_m, bins=e3, density=True, histtype="step", color="#0072b2", lw=1.6, label=f"returns within 1 m of the RTK points, minus RTK (n {len(Q)})")
ax.set_yscale("log"); ax.axvline(0, color="0.4", lw=0.8, ls=":"); ax.set_xlim(-0.50 + off, 1.00 + off); ax.set_xlabel("height above the surveyed ground (m)"); ax.set_ylabel("probability density (1/m), log"); ax.legend(frameon=False, fontsize=8); ax.spines[["top", "right"]].set_visible(False)
fig.colorbar(sm, ax=axs[0], shrink=0.7, label="median intensity in the bar"); fig.tight_layout(); fig.savefig("figures/rtk_2026-10-25_floodplain_region_vs_dem_log.png", dpi=int(os.environ.get("FIGDPI", 160))); say("wrote figures/rtk_2026-10-25_floodplain_region_vs_dem_log.png")
LOG.close()

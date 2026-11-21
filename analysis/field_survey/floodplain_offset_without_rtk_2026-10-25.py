#!/usr/bin/env python
"""Can the sward depth (DEM minus soil, 0.130 m by RTK) be read from the lidar data alone? (Andy, 2026-11-21.)
Three internal estimators, each scored against the RTK value, which is NOT used in building them:
 A. the lower tail of the region PDF (all returns minus the official 1 m DEM, region_sward_mask): the log-density slope
    per 1 cm and the median intensity per 1 cm, looking for a shoulder or an intensity floor where soil hits accumulate;
 B. pits: 1 m DEM cells inside the sward region that lie well below their 7 x 7 m neighborhood median and are dim
    (mean intensity below 3500; both cuts MINE), i.e. bare patches; the distribution of pit depths;
 C. the track near RS+ A 13: cross-track DEM profiles E 579945-579985 for N 4885470-4885530, each aligned on its
    intensity minimum, stacked; the depression depth of the mean profile.
Outputs  data/field/2026-10-25/intensity/offset_without_rtk_log.txt; figures/rtk_2026-10-25_floodplain_offset_without_rtk.png
"""
import os
import numpy as np, pandas as pd, laspy, rasterio
from scipy.ndimage import median_filter
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
ROOT = os.path.dirname(os.path.abspath(__file__)) + "/../.."; os.chdir(ROOT)
OUT = "data/field/2026-10-25/intensity"; LOG = open(f"{OUT}/offset_without_rtk_log.txt", "w")
def say(*a):
    s = " ".join(str(x) for x in a); print(s); LOG.write(s + "\n")
E0, E1, N0, N1 = 579780.0, 580160.0, 4885280.0, 4885760.0; RES = 1.0; nx, ny = 380, 480
with rasterio.open(f"{OUT}/intensity_2021_1m_mean.tif") as r: MI = r.read(1).astype(float); MI[MI == r.nodata] = np.nan
with rasterio.open(f"{OUT}/region_sward_mask.tif") as r: REG = r.read(1) == 1
with rasterio.open("data/field/2026-10-25/dem3dep/usgs_1m_dem_field.tif") as r: DEM = r.read(1).astype(float); DEM[DEM == r.nodata] = np.nan
# ---- A. tail shape and intensity profile from the region returns
parts = []
for t in ("4342-29-64", "4358-29-01"):
    with laspy.open(f"data/after/valley/{t}_3dep_full.laz") as r:
        for ch in r.chunk_iterator(5_000_000):
            x = np.asarray(ch.x); y = np.asarray(ch.y); m = (x >= E0) & (x < E1) & (y >= N0) & (y < N1)
            if t == "4358-29-01": m &= x > 580030.0
            if not m.any(): continue
            i = np.clip(((N1 - y[m]) / RES).astype(int), 0, ny - 1); j = np.clip(((x[m] - E0) / RES).astype(int), 0, nx - 1); k = REG[i, j]
            if k.any(): parts.append(pd.DataFrame(dict(dz=np.asarray(ch.z)[m][k] - DEM[i[k], j[k]], inten=np.asarray(ch.intensity)[m][k].astype(float))))
R = pd.concat(parts, ignore_index=True); R = R[np.isfinite(R.dz)]
edges = np.arange(-0.40, 0.10 + 0.005, 0.01); mids = 0.5 * (edges[1:] + edges[:-1]); cnt = np.histogram(R.dz, bins=edges)[0]; dens = cnt / (len(R) * 0.01); b = np.digitize(R.dz, edges) - 1
med = np.array([R.inten[b == k].median() if (b == k).sum() >= 20 else np.nan for k in range(len(mids))])
logd = np.log10(np.where(dens > 0, dens, np.nan)); slope = np.gradient(logd, 0.01)   # decades per m
say("A. region returns minus DEM, lower tail: bin center, n, density, log10 slope (decades/m, smoothed over 5 cm), median intensity")
sm = pd.Series(slope).rolling(5, center=True).mean().to_numpy()
for k in range(len(mids)):
    if -0.36 <= mids[k] <= 0.0: say(f"   {mids[k]:+.3f}  {cnt[k]:7d}  {dens[k]:8.4f}  {sm[k]:+6.1f}  {med[k]:6.0f}")
# ---- B. pits
nb = median_filter(np.where(np.isfinite(DEM), DEM, np.nanmedian(DEM)), size=7); depth = nb - DEM
cand = REG & np.isfinite(depth) & (MI < 3500); pits = depth[cand & (depth > 0.04)]
say(f"\nB. pits: region cells dim (< 3500) and lower than the 7 x 7 m median by > 0.04 m: {len(pits)} cells; depth median {np.median(pits):.3f}, p25 {np.percentile(pits, 25):.3f}, p75 {np.percentile(pits, 75):.3f}; mode (1 cm bins) {0.01 * np.argmax(np.histogram(pits, bins=np.arange(0.04, 0.40, 0.01))[0]) + 0.045:.3f}")
bright = depth[REG & np.isfinite(depth) & (MI >= 3500) & (depth > 0.04)]; say(f"   for comparison, bright low cells (>= 3500): {len(bright)} cells, depth median {np.median(bright):.3f}")
for lo, hi in ((0.04, 0.08), (0.08, 0.12), (0.12, 0.16), (0.16, 0.20), (0.20, 0.30), (0.30, 1.0)): say(f"   dim pits {lo:.2f}-{hi:.2f} m deep: {int(((pits >= lo) & (pits < hi)).sum())}")
# ---- C. the track transect
rows = []
for n in np.arange(4885470, 4885531, 1.0):
    i = int((N1 - n) / RES); j0, j1 = int((579945 - E0) / RES), int((579985 - E0) / RES); seg_i = MI[i, j0:j1]; seg_d = DEM[i, j0:j1]
    if not np.isfinite(seg_i).all(): continue
    c = int(np.argmin(seg_i)); 
    if 6 <= c <= len(seg_i) - 7: rows.append(seg_d[c - 6:c + 7] - np.nanmedian(np.r_[seg_d[c - 6:c - 3], seg_d[c + 4:c + 7]]))
prof = np.array(rows); mp = np.nanmedian(prof, 0); say(f"\nC. track transect, {len(rows)} rows aligned on the intensity minimum; median cross-track DEM profile relative to the flanks (m), offsets -6..+6 m: " + " ".join(f"{v:+.3f}" for v in mp) + f"; center depression {mp[6]:+.3f} m, mean of center 3 cells {np.nanmean(mp[5:8]):+.3f}")
say(f"   track intensity at the minimum: median {np.median([MI[int((N1 - n) / RES), int((579945 - E0) / RES) + int(np.argmin(MI[int((N1 - n) / RES), int((579945 - E0) / RES):int((579985 - E0) / RES)]))] for n in np.arange(4885470, 4885531, 1.0)]):.0f}")
say("\nRTK answer, for the score only: DEM minus soil 0.130 m (NMAD 0.031).")
fig, axs = plt.subplots(1, 3, figsize=(17, 4.6))
ax = axs[0]; ax.plot(mids, dens, color="#0072b2", lw=1.6); ax.set_yscale("log"); ax.set_xlabel("return minus official DEM (m)"); ax.set_ylabel("density (1/m), log"); ax.axvline(-0.13, color="#d55e00", lw=0.8, ls="--"); ax2 = ax.twinx(); ax2.plot(mids, med, color="0.3", lw=1.2); ax2.set_ylabel("median intensity"); ax.set_xlim(-0.40, 0.10)
ax = axs[1]; ax.hist(pits, bins=np.arange(0.04, 0.40, 0.01), color="#0072b2", alpha=0.8, label="dim pits"); ax.hist(bright, bins=np.arange(0.04, 0.40, 0.01), histtype="step", color="#e69f00", lw=1.4, label="bright low cells"); ax.axvline(0.13, color="#d55e00", lw=0.8, ls="--"); ax.set_xlabel("DEM depth below the 7 x 7 m neighborhood median (m)"); ax.set_ylabel("cells"); ax.legend(frameon=False, fontsize=8)
ax = axs[2]; ax.plot(np.arange(-6, 7), mp, "o-", color="#0072b2"); ax.axhline(-0.13, color="#d55e00", lw=0.8, ls="--"); ax.set_xlabel("distance from the track center (m)"); ax.set_ylabel("DEM relative to the flanks (m)")
for a in axs: a.spines[["top", "right"]].set_visible(False)
fig.tight_layout(); fig.savefig("figures/rtk_2026-10-25_floodplain_offset_without_rtk.png", dpi=int(os.environ.get("FIGDPI", 160))); say("wrote figures/rtk_2026-10-25_floodplain_offset_without_rtk.png (dashed = the RTK value, drawn for the comparison only)")
LOG.close()

#!/usr/bin/env python
"""The shoulder near -13 cm in the stacked return record (Andy, 2026-11-21). For the three distance-from-crown bands
and the sward region, returns minus the official 1 m DEM: (1) the log-density slope per centimeter of the lower tail
and the depth where the steep peak flank hands over to the flat tail (maximum curvature of log density, smoothed over
3 cm; MINE); (2) a mirror test: the excess of the lower flank over the upper flank at the same distance from the peak,
dens(-h) - dens(+h), which removes the symmetric part of the sward peak and leaves the sub-surface population; its
mode and median. Scored afterwards against the RTK offset of each band (DEM minus soil: crown edge +0.120, near +0.145,
open +0.121, all +0.130 m), which is not used in the estimate.
Outputs  data/field/2026-10-25/intensity/shoulder_log.txt; figures/rtk_2026-10-25_floodplain_shoulder.png
"""
import os
import numpy as np, pandas as pd, laspy, rasterio
from scipy.ndimage import distance_transform_edt
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
ROOT = os.path.dirname(os.path.abspath(__file__)) + "/../.."; os.chdir(ROOT)
OUT = "data/field/2026-10-25/intensity"; LOG = open(f"{OUT}/shoulder_log.txt", "w")
def say(*a):
    s = " ".join(str(x) for x in a); print(s); LOG.write(s + "\n")
E0, E1, N0, N1 = 579780.0, 580160.0, 4885280.0, 4885760.0; RES = 1.0; nx, ny = 380, 480
with rasterio.open(f"{OUT}/intensity_2021_1m_tree.tif") as r: TREE = r.read(1) == 1
with rasterio.open(f"{OUT}/region_sward_mask.tif") as r: REG = r.read(1) == 1
with rasterio.open("data/field/2026-10-25/dem3dep/usgs_1m_dem_field.tif") as r: DEM = r.read(1).astype(float); DEM[DEM == r.nodata] = np.nan
z10 = np.load("data/field/2026-10-25/imagery/floodplain_model_features.npz", allow_pickle=True); xs, ys, fpm = z10["xs"], z10["ys"], z10["fpm"].astype(bool)
cy = N1 - (np.arange(ny) + 0.5) * RES; cx = E0 + (np.arange(nx) + 0.5) * RES
FP = fpm[np.clip(np.round((ys[0] - cy) / 10).astype(int), 0, len(ys) - 1)[:, None], np.clip(np.round((cx - xs[0]) / 10).astype(int), 0, len(xs) - 1)[None, :]]
dist = distance_transform_edt(~TREE) * RES; BAND = np.full((ny, nx), -1); BAND[(dist > 0) & (dist <= 3)] = 0; BAND[(dist > 3) & (dist <= 10)] = 1; BAND[dist > 10] = 2; BAND[~FP | TREE] = -1
parts = []
for t in ("4342-29-64", "4358-29-01"):
    with laspy.open(f"data/after/valley/{t}_3dep_full.laz") as r:
        for ch in r.chunk_iterator(5_000_000):
            x = np.asarray(ch.x); y = np.asarray(ch.y); m = (x >= E0) & (x < E1) & (y >= N0) & (y < N1)
            if t == "4358-29-01": m &= x > 580030.0
            if not m.any(): continue
            i = np.clip(((N1 - y[m]) / RES).astype(int), 0, ny - 1); j = np.clip(((x[m] - E0) / RES).astype(int), 0, nx - 1); b = BAND[i, j]; k = (b >= 0) | REG[i, j]
            parts.append(pd.DataFrame(dict(dz=np.asarray(ch.z)[m][k] - DEM[i[k], j[k]], band=b[k], reg=REG[i, j][k])))
R = pd.concat(parts, ignore_index=True); R = R[np.isfinite(R.dz)]
sets = {"crown edge 0-3 m": R[R.band == 0].dz, "near trees 3-10 m": R[R.band == 1].dz, "open > 10 m": R[R.band == 2].dz, "sward region": R[R.reg].dz}
RTK = {"crown edge 0-3 m": 0.120, "near trees 3-10 m": 0.145, "open > 10 m": 0.121, "sward region": 0.130}
W = 0.01; edges = np.arange(-0.45, 0.45 + W / 2, W); mids = 0.5 * (edges[1:] + edges[:-1])
fig, axs = plt.subplots(1, 2, figsize=(13, 5)); COL = ["#d55e00", "#e69f00", "#0072b2", "k"]   # the mirror test is logged but not drawn: the upper flank carries the tall-herb returns and is the heavier side
say(f"{'set':20s} {'n':>8s} {'peak':>6s} {'knee (max curvature of log density)':>36s} {'mirror-excess mode':>19s} {'mirror-excess median':>21s} {'excess share':>12s} {'RTK DEM-soil':>12s}")
for c, (name, v) in zip(COL, sets.items()):
    cnt = np.histogram(v, bins=edges)[0]; dens = cnt / (len(v) * W); pk = mids[np.argmax(dens)]
    ld = np.log10(np.where(dens > 0, dens, np.nan)); ls = pd.Series(ld).rolling(3, center=True).mean().to_numpy(); d1 = np.gradient(ls, W); d2 = np.gradient(d1, W)
    lower = (mids < pk - 0.04) & (mids > pk - 0.35); knee = mids[lower][np.nanargmax(d2[lower])]   # most positive curvature = flank-to-tail handover
    h = np.arange(W, 0.40, W); lo = np.interp(pk - h, mids, dens); hi = np.interp(pk + h, mids, dens); ex = np.clip(lo - hi, 0, None)
    mode = h[np.argmax(ex)]; cdf = np.cumsum(ex) / ex.sum(); med = h[np.searchsorted(cdf, 0.5)]; share = 100 * ex.sum() * W
    # handover depth: where the smoothed slope first comes within 25 % of its tail plateau (median slope over -0.25..-0.15 below the peak); MINE
    rel = mids - pk; plateau = np.nanmedian(d1[(rel > -0.25) & (rel < -0.15)]); cand_ = rel[(rel < -0.02) & (rel > -0.30) & (d1 <= 1.25 * plateau)]; hand = cand_.max() if len(cand_) else np.nan
    say(f"{name:20s} {len(v):8d} {pk:+6.3f} {knee - pk:+36.3f} {-mode:+19.3f} {-med:+21.3f} {share:11.2f}% {-RTK[name]:+12.3f} | tail plateau {plateau:4.1f} decades/m (e-folding {100 / (plateau * np.log(10)):.1f} cm), handover at {hand:+.3f} m below the peak")
    say("     slope (decades/m) at -0.04/-0.06/-0.08/-0.10/-0.12/-0.14/-0.16/-0.20: " + " ".join(f"{np.interp(pk - h, mids, d1):5.1f}" for h in (0.04, 0.06, 0.08, 0.10, 0.12, 0.14, 0.16, 0.20)))
    axs[0].plot(mids - pk, dens, color=c, lw=1.5, label=name); axs[1].plot(mids - pk, d1, color=c, lw=1.2, label=name); axs[1].axvline(-RTK[name], color=c, lw=0.8, ls="--")
axs[0].set_yscale("log"); axs[0].set_xlim(-0.40, 0.20); axs[0].set_xlabel("return minus DEM, relative to the peak (m)"); axs[0].set_ylabel("density (1/m), log")
axs[1].set_xlim(-0.40, -0.02); axs[1].set_ylim(0, 40); axs[1].set_xlabel("return minus DEM, relative to the peak (m)"); axs[1].set_ylabel("slope of log10 density (decades per m)")
for a in axs: a.spines[["top", "right"]].set_visible(False)
axs[0].legend(frameon=False, fontsize=8); fig.tight_layout(); fig.savefig("figures/rtk_2026-10-25_floodplain_shoulder.png", dpi=int(os.environ.get("FIGDPI", 160))); say("wrote figures/rtk_2026-10-25_floodplain_shoulder.png (dashed in the slope panel: the RTK soil depth of each band, for the score only)")
LOG.close()

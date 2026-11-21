#!/usr/bin/env python
"""Isolating the soil component of the stacked return-height distribution without RTK, two constrained ways
(Andy, 2026-11-21), after free Gaussian mixtures failed (floodplain_mixture_fit_2026-10-25.py: the core has exponential
flanks and the second Gaussian copies them instead of finding the soil).
 1. Core from the upper flank: the sward-top core is fitted on the UPPER side only, +0.01 to +0.06 m above the peak
    (soil-free; below the tall-herb layer), as a two-sided exponential (Laplace) in log density; it is reflected to the
    lower side and subtracted; what remains below the peak is the soil-plus-bridge population, whose mode and median
    are read off. Fit window MINE.
 2. Intensity conditioning: height histograms of the region returns by intensity class; the dim classes should be the
    soil hits, so their height mode should sit at the soil.
Sets: sward region and the three distance-from-crown bands, returns minus the official 1 m DEM. RTK soil depths
(0.130 region, 0.120 crown edge, 0.145 near, 0.121 open) are used only to score.
Outputs  data/field/2026-10-25/intensity/soil_component_log.txt; figures/rtk_2026-10-25_floodplain_soil_component.png
"""
import os
import numpy as np, pandas as pd, laspy, rasterio
from scipy.ndimage import distance_transform_edt
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
ROOT = os.path.dirname(os.path.abspath(__file__)) + "/../.."; os.chdir(ROOT)
OUT = "data/field/2026-10-25/intensity"; LOG = open(f"{OUT}/soil_component_log.txt", "w")
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
            parts.append(pd.DataFrame(dict(dz=np.asarray(ch.z)[m][k] - DEM[i[k], j[k]], inten=np.asarray(ch.intensity)[m][k].astype(float), band=b[k], reg=REG[i, j][k])))
R = pd.concat(parts, ignore_index=True); R = R[np.isfinite(R.dz)]
sets = {"sward region": (R[R.reg], 0.130), "crown edge 0-3 m": (R[R.band == 0], 0.120), "near trees 3-10 m": (R[R.band == 1], 0.145), "open > 10 m": (R[R.band == 2], 0.121)}
W = 0.01; edges = np.arange(-0.45, 0.45 + W / 2, W); mids = 0.5 * (edges[1:] + edges[:-1])
fig, axs = plt.subplots(2, 4, figsize=(21, 9))
say("1. core fitted on the upper flank (+0.01..+0.06 m above the peak) as exp(a - |h|/b), reflected, subtracted:")
say(f"{'set':20s} {'peak':>6s} {'core e-fold b':>13s} {'residual below the peak: share':>30s} {'mode':>7s} {'median':>7s} {'p25-p75':>15s} {'RTK soil':>8s}")
for col, (name, (q, rtk)) in enumerate(sets.items()):
    v = q.dz.to_numpy(); cnt = np.histogram(v, bins=edges)[0]; dens = cnt / (len(v) * W); pk = mids[np.argmax(dens)]; h = mids - pk
    up = (h >= 0.01) & (h <= 0.06) & (dens > 0); A = np.polyfit(h[up], np.log(dens[up]), 1); bcore = -1 / A[0]; core = np.exp(A[1] - np.abs(h) / bcore)
    resid = np.clip(dens - core, 0, None); low = h < -0.005; rr = resid[low]; hh = h[low]; share = rr.sum() * W
    cdf = np.cumsum(rr[::-1])[::-1] / rr.sum() if rr.sum() > 0 else None   # cumulative from the peak downward
    mode = hh[np.argmax(rr)]; med = hh[np.argmin(np.abs(np.cumsum(rr) / rr.sum() - 0.5))]; p25 = hh[np.argmin(np.abs(np.cumsum(rr) / rr.sum() - 0.25))]; p75 = hh[np.argmin(np.abs(np.cumsum(rr) / rr.sum() - 0.75))]
    deep = hh < -0.05; rd = rr[deep]; hd = hh[deep]; cd = np.cumsum(rd) / rd.sum()   # the residual beyond the core's own near-peak asymmetry (5 cm cut, MINE)
    say(f"{name:20s} {pk:+6.3f} {bcore:13.4f} {100 * share:29.2f}% {mode:+7.3f} {med:+7.3f} {p25:+7.3f} {p75:+7.3f} {-rtk:+8.3f} | residual below -0.05: share {100 * rd.sum() * W:.2f} %, mode {hd[np.argmax(rd)]:+.3f}, median {hd[np.argmin(np.abs(cd - 0.5))]:+.3f}, p25 {hd[np.argmin(np.abs(cd - 0.25))]:+.3f}, p75 {hd[np.argmin(np.abs(cd - 0.75))]:+.3f}")
    ax = axs[0, col]; ax.plot(h, dens, color="0.3", lw=1.2, label="data"); ax.plot(h, core, color="#0072b2", lw=1.2, ls="--", label="core from the upper flank"); ax.plot(h[low], rr, color="#d55e00", lw=1.6, label="residual below the peak"); ax.axvline(-rtk, color="k", lw=0.8, ls=":", label=f"RTK soil {-rtk:+.3f}")
    ax.set_yscale("log"); ax.set_ylim(5e-3, 40); ax.set_xlim(-0.45, 0.20); ax.set_title(name, fontsize=10); ax.set_xlabel("return minus DEM, relative to the peak (m)"); ax.set_ylabel("density (1/m), log"); ax.legend(frameon=False, fontsize=7); ax.spines[["top", "right"]].set_visible(False)
    # 2. by intensity class
    ax = axs[1, col]; cols = plt.get_cmap("viridis")(np.linspace(0.1, 0.9, 5)); say(f"   {name}: height by intensity class (n, mode, median, p10, share below -0.08):") if col == 0 else None
    for c_, (lo_, hi_) in zip(cols, ((0, 2500), (2500, 3000), (3000, 3500), (3500, 4000), (4000, 9000))):
        s_ = q[(q.inten >= lo_) & (q.inten < hi_)].dz; 
        if len(s_) < 200: continue
        hc = np.histogram(s_, bins=edges)[0] / (len(s_) * W); ax.plot(mids, hc, color=c_, lw=1.4, label=f"intensity {lo_}-{hi_} (n {len(s_):,})")
        say(f"      {name:18s} {lo_:5d}-{hi_:<5d} n {len(s_):7d} mode {mids[np.argmax(hc)]:+.3f} median {s_.median():+.3f} p10 {s_.quantile(.1):+.3f} below -0.08: {100 * (s_ < -0.08).mean():5.2f} %")
    ax.axvline(-rtk, color="k", lw=0.8, ls=":"); ax.set_yscale("log"); ax.set_ylim(5e-3, 60); ax.set_xlim(-0.45, 0.20); ax.set_xlabel("return minus official DEM (m)"); ax.set_ylabel("density (1/m), log"); ax.legend(frameon=False, fontsize=7); ax.spines[["top", "right"]].set_visible(False)
fig.tight_layout(); fig.savefig("figures/rtk_2026-10-25_floodplain_soil_component.png", dpi=int(os.environ.get("FIGDPI", 150))); say("wrote figures/rtk_2026-10-25_floodplain_soil_component.png (dotted: RTK soil, for the score only)")
LOG.close()

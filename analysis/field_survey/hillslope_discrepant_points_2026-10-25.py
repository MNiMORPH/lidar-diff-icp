#!/usr/bin/env python
"""Hillslope points where RTK and lidar disagree most (Andy, 2026-11-21): per point, the near-ground median of the
slope-normal return heights within 1 m (as in hillslope_points_returns_2026-10-25.py); chosen are the three highest,
the three lowest, and the two points with no return within 5 cm of the ground that hold the most returns (MINE).
Their 1 m distributions with intensity-colored bars, and one map: hillshade of the official 1 m bluff DEM with all
hillslope FIX points and the chosen ones circled and labeled.
Outputs  data/field/2026-10-25/intensity/hillslope_discrepant_log.txt; figures/rtk_2026-10-25_hillslope_discrepant_points.png
"""
import os
import numpy as np, pandas as pd, laspy, rasterio
from scipy.ndimage import map_coordinates
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, LightSource
ROOT = os.path.dirname(os.path.abspath(__file__)) + "/../.."; os.chdir(ROOT)
OUT = "data/field/2026-10-25/intensity"; LOG = open(f"{OUT}/hillslope_discrepant_log.txt", "w")
def say(*a):
    s = " ".join(str(x) for x in a); print(s); LOG.write(s + "\n")
with rasterio.open("data/field/2026-10-25/dem3dep/usgs_1m_dem_bluff_wide.tif") as r: DEM = r.read(1).astype(float); DEM[DEM == r.nodata] = np.nan; tr = r.transform; E0, N1 = tr.c, tr.f; E1 = E0 + r.width; N0 = N1 - r.height
gy, gx = np.gradient(DEM, 1.0); gy = -gy
def at(a, x, y): return map_coordinates(a, [((N1 - y) - 0.5), ((x - E0) - 0.5)], order=1, mode="nearest")
P = pd.read_csv("data/field/2026-10-25/rtk_points_all.csv"); H = P[P.group.str.startswith("hillslope") & (P.status == "FIX")].copy()
assert (H.easting_tied.between(E0 + 2, E0 + DEM.shape[1] - 2) & H.northing_tied.between(N1 - DEM.shape[0] + 2, N1 - 2)).all(), "a hillslope point lies outside the DEM window"
H["gx"] = at(gx, H.easting_tied.values, H.northing_tied.values); H["gy"] = at(gy, H.easting_tied.values, H.northing_tied.values); H["slope"] = np.degrees(np.arctan(np.hypot(H.gx, H.gy))); H["dem"] = at(DEM, H.easting_tied.values, H.northing_tied.values)
bbox = (H.easting_tied.min() - 4, H.easting_tied.max() + 4, H.northing_tied.min() - 4, H.northing_tied.max() + 4); parts = []
with laspy.open("data/after/valley/4342-29-64_3dep_full.laz") as r:
    for ch in r.chunk_iterator(5_000_000):
        x = np.asarray(ch.x); y = np.asarray(ch.y); m = (x >= bbox[0]) & (x <= bbox[1]) & (y >= bbox[2]) & (y <= bbox[3])
        if m.any(): parts.append(pd.DataFrame(dict(x=x[m], y=y[m], z=np.asarray(ch.z)[m], inten=np.asarray(ch.intensity)[m].astype(float), nr=np.asarray(ch.number_of_returns)[m], rn=np.asarray(ch.return_number)[m])))
G = pd.concat(parts, ignore_index=True); rows = []; stacks = {}
for _, p in H.iterrows():
    d = np.hypot(G.x - p.easting_tied, G.y - p.northing_tied); w = G[d <= 1].copy()
    plane = p.navd88_tied_m + p.gx * (w.x - p.easting_tied) + p.gy * (w.y - p.northing_tied); w["h"] = (w.z - plane) * np.cos(np.radians(p.slope)); key = (p.group, int(p.point)); stacks[key] = w
    near = w.h[w.h.abs() < 0.5]
    rows.append(dict(group=p.group, point=int(p.point), t=p.t_hhmm, rms_v=p.rms_v, tilt=p.tilt_deg, slope=p.slope, E=p.easting_tied, N=p.northing_tied, n=len(w), med_near=near.median() if len(near) else np.nan, ground_share=(w.h.abs() < 0.05).mean(), canopy_share=(w.h >= 2).mean(), dem_minus_rtk_n=(p.dem - p.navd88_tied_m) * np.cos(np.radians(p.slope)), inten_near=w.inten[w.h.abs() < 0.5].median() if len(near) else np.nan))
S = pd.DataFrame(rows); S.to_csv(f"{OUT}/hillslope_points_summary.csv", index=False, float_format="%.4f")
pd.set_option("display.width", 250); say("per point (1 m, slope-normal to the surveyed ground):"); say(S.round(3).to_string(index=False))
hi = S.dropna(subset=["med_near"]).nlargest(3, "med_near"); lo = S.dropna(subset=["med_near"]).nsmallest(3, "med_near"); nog = S[S.ground_share == 0].nlargest(2, "n")
sel = pd.concat([hi.assign(kind="lidar floor highest above the pole tip"), lo.assign(kind="lidar floor lowest / below the pole tip"), nog.assign(kind="no return within 5 cm")])
say("\nchosen:"); say(sel[["kind", "group", "point", "t", "rms_v", "tilt", "slope", "n", "med_near", "ground_share", "canopy_share", "dem_minus_rtk_n", "inten_near"]].round(3).to_string(index=False))
norm = Normalize(vmin=1000, vmax=4500); cmap = plt.get_cmap("viridis"); e2 = np.arange(-0.6, 1.0 + 0.005, 0.02); m2 = 0.5 * (e2[1:] + e2[:-1])
fig = plt.figure(figsize=(22, 10)); gs = fig.add_gridspec(2, 5, width_ratios=[1, 1, 1, 1, 1.6])
for k, (_, r) in enumerate(sel.iterrows()):
    ax = fig.add_subplot(gs[k // 4, k % 4]); w = stacks[(r.group, r.point)]; c = np.histogram(w.h, bins=e2)[0]; b = np.digitize(w.h, e2) - 1; med = np.array([w.inten[b == i].median() if (b == i).sum() >= 2 else np.nan for i in range(len(m2))])
    ax.bar(m2, c, width=0.02, color=[cmap(norm(v)) if np.isfinite(v) else "0.85" for v in med], edgecolor="none"); ax.axvline(0, color="0.4", lw=0.8, ls=":"); ax.axvline(r.dem_minus_rtk_n, color="#d55e00", lw=1.0, ls="--")
    ax.set_xlim(-0.6, 1.0); short = {"lidar floor highest above the pole tip": "floor high", "lidar floor lowest / below the pole tip": "floor low", "no return within 5 cm": "no ground hit"}[r.kind]
    ax.set_title(f"{r.group.replace('hillslope_', '')} {r.point}: {short}, {r.t}\nslope {r.slope:.0f} deg, rms_v {r.rms_v:.3f}, n {r.n}, canopy {100 * r.canopy_share:.0f} %", fontsize=8); ax.set_xlabel("slope-normal height above the surveyed ground (m)", fontsize=8); ax.set_ylabel("returns per 2 cm", fontsize=8); ax.tick_params(labelsize=7); ax.spines[["top", "right"]].set_visible(False)
ax = fig.add_subplot(gs[:, 4]); ls = LightSource(azdeg=315, altdeg=45); ax.imshow(ls.hillshade(np.nan_to_num(DEM, nan=np.nanmin(DEM)), vert_exag=1, dx=1, dy=1), cmap="gray", extent=(E0, E1, N0, N1))
for g, c in (("hillslope_RS3", "#0072b2"), ("hillslope_T4a", "#d55e00")):
    q = S[S.group == g]; ax.plot(q.E, q.N, ".", color=c, ms=6, label=g.replace("hillslope_", "RS3 "))
for _, r in sel.iterrows():
    ax.plot([r.E], [r.N], "o", ms=13, mfc="none", mec="red", mew=1.8); ax.annotate(f"{r.group.replace('hillslope_', '')} {r.point}", (r.E, r.N), xytext=(8, 4), textcoords="offset points", color="red", fontsize=8, fontweight="bold")
pad = 15; ax.set_xlim(S.E.min() - pad, S.E.max() + pad); ax.set_ylim(S.N.min() - pad, S.N.max() + pad); ax.set_aspect("equal"); ax.legend(fontsize=8, loc="lower right"); ax.tick_params(labelsize=7); ax.set_title("hillshade of the official 1 m DEM; chosen points circled", fontsize=9)
sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap); sm.set_array([]); fig.colorbar(sm, ax=fig.axes[:8], shrink=0.5, label="median intensity in the bar (dashed orange: official DEM minus RTK at the point)")
fig.savefig("figures/rtk_2026-10-25_hillslope_discrepant_points.png", dpi=int(os.environ.get("FIGDPI", 150)), bbox_inches="tight"); say("wrote figures/rtk_2026-10-25_hillslope_discrepant_points.png")
LOG.close()

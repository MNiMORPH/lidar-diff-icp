import numpy as np, pandas as pd, rasterio, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import os; os.chdir("/home/awickert/projects/lidar-diff-icp")
P = pd.read_csv("data/field/2026-10-25/rtk_points_all.csv"); F = P[P.group.str.startswith("floodplain") & (P.status == "FIX")]
with rasterio.open("data/field/2026-10-25/intensity/intensity_2021_1m_mean.tif") as r: INT = r.read(1).astype(float); INT[INT == r.nodata] = np.nan; ext = (r.bounds.left, r.bounds.right, r.bounds.bottom, r.bounds.top)
fig, ax = plt.subplots(figsize=(8.5, 11)); ax.imshow(INT, extent=ext, cmap="gray", vmin=1500, vmax=4800)
COL = {"floodplain_RSplus_B": "#0072b2", "floodplain_RSplus_A": "#e69f00", "floodplain_CORS_Wickert": "#009e73"}; LAB = {"floodplain_RSplus_B": "RS+ B", "floodplain_RSplus_A": "RS+ A", "floodplain_CORS_Wickert": "CORS"}
for g, q in F.groupby("group"): ax.plot(q.easting_tied, q.northing_tied, ".", color=COL[g], ms=5, label=LAB[g])
for g, p, lab in (("floodplain_RSplus_B", 3, "RS+ B 3: crown edge, low intensity overall"), ("floodplain_RSplus_B", 100, "RS+ B 100: canopy gap, high intensity only below the surveyed ground"), ("floodplain_RSplus_A", 13, "RS+ A 13: open sward beside a track, high intensity overall")):
    r = F[(F.group == g) & (F.point == p)].iloc[0]; ax.plot([r.easting_tied], [r.northing_tied], "o", ms=14, mfc="none", mec="red", mew=2); ax.annotate(f"{LAB[g]} {p}", (r.easting_tied, r.northing_tied), xytext=(12, 6), textcoords="offset points", color="red", fontsize=10, fontweight="bold")
ax.plot([579979], [4885340], "w+", ms=12, mew=2); ax.annotate("lot", (579979, 4885340), xytext=(6, -14), textcoords="offset points", color="w", fontsize=9)
ax.set_xlim(579850, 580110); ax.set_ylim(4885320, 4885700); ax.set_aspect("equal"); ax.legend(loc="upper right", fontsize=9); ax.tick_params(labelsize=8)
fig.tight_layout(); fig.savefig("figures/rtk_2026-10-25_floodplain_three_points_overview.png", dpi=160); print("ok")

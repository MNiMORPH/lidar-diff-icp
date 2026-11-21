#!/usr/bin/env python
"""The three floodplain points whose 1 m window holds returns below the surveyed ground, RS+ B 3, RS+ B 100 and
RS+ A 13 (Andy, 2026-11-21: cut grass, or RTK standing high?). Per point: the 1 m and 3 m return-height histograms
of the 2021 lidar with bars colored by median intensity (1 cm bins, MINE), beside a normal neighbor of each for
comparison; and a table of what the other layers say at those spots.
Outputs  figures/rtk_2026-10-25_floodplain_three_points.png; data/field/2026-10-25/returns/floodplain_three_points_log.txt
"""
import os
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
ROOT = os.path.dirname(os.path.abspath(__file__)) + "/../.."; os.chdir(ROOT)
D = "data/field/2026-10-25/returns"; LOG = open(f"{D}/floodplain_three_points_log.txt", "w")
def say(*a):
    s = " ".join(str(x) for x in a); print(s); LOG.write(s + "\n")
P = pd.read_csv("data/field/2026-10-25/rtk_points_all.csv"); C = pd.read_csv("data/field/2026-10-25/intensity/cell_profiles_1m.csv"); OP = pd.read_csv(f"{D}/floodplain_offset_points.csv")
TAG = {"floodplain_RSplus_B": "rsplus_b", "floodplain_RSplus_A": "rsplus_a"}
cases = [("floodplain_RSplus_B", 3, 2), ("floodplain_RSplus_B", 100, 98), ("floodplain_RSplus_A", 13, 12)]   # (group, suspect point, neighbor for comparison)
names = {0: "smooth sward", 1: "low sward", 2: "tall herb", 3: "sward under crowns", 4: "shrub", 5: "tree"}
E0, N1, nx = 579780.0, 4885760.0, 380
def cell_of(e, n): return int((N1 - n)) * nx + int((e - E0))
rows = []
for g, p, nb in cases:
    t21 = pd.read_csv(f"{D}/{TAG[g]}_returns_2021.csv.gz"); t08 = pd.read_csv(f"{D}/{TAG[g]}_returns_2008.csv")
    for q, lab in ((p, "suspect"), (nb, "neighbor")):
        r = P[(P.group == g) & (P.point == q)].iloc[0]; a = t21[(t21.point == q) & t21.radius_1m]; a3 = t21[(t21.point == q) & t21.radius_3m]; b = t08[(t08.point == q) & t08.radius_3m]
        c = C.iloc[cell_of(r.easting_tied, r.northing_tied)]; o = OP[(OP.group == g) & (OP.point == q)].iloc[0]
        rows.append(dict(point=f"{g.replace('floodplain_', '')} {q}", role=lab, t=r.t_hhmm, status=r.status, rms_v=r.rms_v, rtk_z=r.navd88_tied_m, n1m=len(a), med_dz_1m=a.dz_m.median(), p10_dz_1m=a.dz_m.quantile(.1), p90_dz_1m=a.dz_m.quantile(.9), med_int_1m=a.intensity.median(),
                         med_dz_3m=a3.dz_m.median(), n08_3m=len(b), med_dz08_3m=b.dz_star_m.median() if len(b) else np.nan, min_dz08_3m=b.dz_star_m.min() if len(b) else np.nan,
                         cell_cluster=names.get(c.cluster, "-") if np.isfinite(c.cluster) else "-", cell_int=c.intensity, cell_ground_proxy_minus_rtk=c.ground_proxy_m - r.navd88_tied_m, ndvi_apr26=o.metric_ndvi_apr, naip21=o.naip_2021, gen2_product=r.gen2_minus_rtk_m, slope=r.slope_deg))
Tb = pd.DataFrame(rows); pd.set_option("display.width", 300); pd.set_option("display.max_columns", 40); say(Tb.round(3).to_string())
# RTK heights along the line around each suspect point
for g, p, nb in cases:
    s = P[(P.group == g) & (P.point.between(p - 3, p + 3))][["point", "t_hhmm", "status", "rms_v", "navd88_tied_m", "gen2_minus_rtk_m", "gen1_minus_rtk_m"]]
    say(f"\n{g.replace('floodplain_', '')} around point {p}:"); say(s.round(3).to_string(index=False))
fig, axs = plt.subplots(2, 3, figsize=(17, 8)); norm = Normalize(vmin=2500, vmax=5000); cmap = plt.get_cmap("viridis"); edges = np.arange(-0.15, 0.45, 0.01); mids = 0.5 * (edges[1:] + edges[:-1])
for j, (g, p, nb) in enumerate(cases):
    t21 = pd.read_csv(f"{D}/{TAG[g]}_returns_2021.csv.gz")
    for i, (q, lab) in enumerate(((p, "suspect point"), (nb, "neighbor"))):
        ax = axs[i, j]; a = t21[(t21.point == q) & t21.radius_1m]; a3 = t21[(t21.point == q) & t21.radius_3m]
        bidx = np.digitize(a.dz_m, edges) - 1; cnt = np.array([(bidx == k).sum() for k in range(len(mids))]); med = np.array([a.intensity[bidx == k].median() if (bidx == k).sum() >= 2 else np.nan for k in range(len(mids))])
        ax.bar(mids, cnt, width=0.01, color=[cmap(norm(m)) if np.isfinite(m) else "0.8" for m in med], edgecolor="none", label=f"within 1 m (n {len(a)})")
        h3, _ = np.histogram(a3.dz_m, bins=edges); ax.step(mids, h3 * len(a) / max(len(a3), 1), where="mid", color="k", lw=0.9, label=f"within 3 m, scaled (n {len(a3)})")
        ax.axvline(0, color="0.4", lw=0.8, ls=":"); ax.set_xlim(-0.15, 0.45); ax.set_title(f"{g.replace('floodplain_', '')} point {q} ({lab})", fontsize=10); ax.set_xlabel("2021 return minus tied RTK ground (m)"); ax.set_ylabel("returns per 1 cm"); ax.legend(frameon=False, fontsize=8); ax.spines[["top", "right"]].set_visible(False)
sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap); sm.set_array([]); fig.colorbar(sm, ax=axs, shrink=0.6, label="median intensity in the bar")
fig.savefig("figures/rtk_2026-10-25_floodplain_three_points.png", dpi=int(os.environ.get("FIGDPI", 160)), bbox_inches="tight"); say("wrote figures/rtk_2026-10-25_floodplain_three_points.png")
LOG.close()

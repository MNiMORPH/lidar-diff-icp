#!/usr/bin/env python3
"""Figures for the 25 October 2026 RTK survey against the two lidar surfaces.
Reads data/field/2026-10-25/rtk_points_all.csv (compile_rtk_2026-10-25.py). Writes
  figures/rtk_2026-10-25_map.png       points on 2025 FSA imagery, coloured by tied NAVD88 minus gen2
  figures/rtk_2026-10-25_profiles.png  the hillslope traverse (RS3 + T4a) and the floodplain
                                       loop (RS+ A/B + CORS) as RTK, gen2 and gen1 elevations
                                       in survey order, FIX points only
"""
import io
import sys
import urllib.request

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import TwoSlopeNorm
from PIL import Image

P = pd.read_csv("data/field/2026-10-25/rtk_points_all.csv")
SUF = sys.argv[1] if len(sys.argv) > 1 else ""   # output-name suffix, so a replot need not clobber an earlier figure
_b2 = P[P.base == "B"].tie_source.iloc[0]; _pc = P[P.antenna_corr_m != 0].groupby("group").antenna_corr_m.max()
TIE = f"B tied to {_b2.split(' under B (')[0]} ({_b2.split(' under B (')[1].split('; ')[0]}" + ("; pole heights corrected: " + ", ".join(f"{g} {v:+.3f} m" for g, v in _pc.items()) if len(_pc) else "") + ")"
F = P[P.status.str.startswith("FIX") & P.tied_minus_gen2_m.notna()].copy()
bb = (579650, 4885200, 580200, 4885750); X0, Y0, X1, Y1 = bb
url = f"https://imageserver.gisdata.mn.gov/cgi-bin/wms?SERVICE=WMS&VERSION=1.1.1&REQUEST=GetMap&LAYERS=fsa2025&STYLES=&SRS=EPSG:26915&BBOX={X0},{Y0},{X1},{Y1}&WIDTH=1400&HEIGHT=1400&FORMAT=image/jpeg"
try: img = np.asarray(Image.open(io.BytesIO(urllib.request.urlopen(url, timeout=120).read())))
except Exception as e: img = None; print("imagery not fetched:", e)
fig, ax = plt.subplots(figsize=(10, 10.6))
if img is not None: ax.imshow(img, extent=[X0, X1, Y0, Y1])
norm = TwoSlopeNorm(vmin=-0.6, vcenter=0, vmax=0.6); mk = {"hillslope_RS3": "o", "hillslope_T4a": "s", "floodplain_RSplus_A": "^", "floodplain_RSplus_B": "v", "floodplain_CORS_Wickert": "D", "benchmark_DR9925": "*"}
for g, s in F.groupby("group"):
    sc = ax.scatter(s.easting_tied, s.northing_tied, c=s.tied_minus_gen2_m, cmap="RdBu", norm=norm, s=48 if g != "benchmark_DR9925" else 200, marker=mk.get(g, "o"), edgecolors="black", linewidths=0.5, label=g)
for lbl, x, y in (("B", 580148.2, 4885379.6),): ax.plot(x, y, marker="P", ms=14, color="yellow", mec="black"); ax.annotate(lbl, (x, y), xytext=(6, 6), textcoords="offset points", color="yellow", weight="bold")
ax.set_xlim(X0, X1); ax.set_ylim(Y0, Y1); ax.set_xticks([]); ax.set_yticks([]); ax.legend(loc="lower left", fontsize=8)
cb = fig.colorbar(sc, ax=ax, shrink=0.6, pad=0.01); cb.set_label("RTK (tied) minus 2021 lidar ground, m; red = lidar above the ground")
ax.set_title(f"25 Oct 2026 RTK survey, FIX points; B1 (a student receiver) removed; yellow P = base B\n{TIE}", fontsize=8)
fig.tight_layout(); fig.savefig(f"figures/rtk_2026-10-25_map{SUF}.png", dpi=110); print(f"figures/rtk_2026-10-25_map{SUF}.png")

fig, axs = plt.subplots(2, 1, figsize=(12, 9))
for ax, groups, title in ((axs[0], ["hillslope_RS3", "hillslope_T4a"], "hillslope traverse west of the road (RS3 circles, T4a squares), survey order within group"), (axs[1], ["floodplain_RSplus_A", "floodplain_RSplus_B", "floodplain_CORS_Wickert"], "floodplain field east of the road (RS+ A triangles up, RS+ B down, CORS diamonds)")):
    off = 0
    for g in groups:
        s = F[F.group == g].sort_values(["t_start"]); d = np.r_[0, np.cumsum(np.hypot(np.diff(s.easting_tied), np.diff(s.northing_tied)))] + off
        ax.plot(d, s.navd88_tied_m, "k-", lw=1); ax.scatter(d, s.navd88_tied_m, c="black", s=16, marker=mk[g], label=f"{g}: RTK tied" if g == groups[0] else None)
        ax.plot(d, s.gen2_z_m, color="#2a78d6", lw=1.2, label="2021 lidar (gen2)" if g == groups[0] else None); ax.plot(d, s.gen1_z_m, color="#e34948", lw=1.2, label="2008 lidar (gen1)" if g == groups[0] else None)
        off = d[-1] + 30
    ax.set_ylabel("NAVD 88 (m)"); ax.set_title(title, fontsize=9); ax.legend(fontsize=8); ax.grid(alpha=0.3)
axs[1].set_xlabel("distance along the walked path (m), groups laid end to end")
fig.tight_layout(); fig.savefig(f"figures/rtk_2026-10-25_profiles{SUF}.png", dpi=110); print(f"figures/rtk_2026-10-25_profiles{SUF}.png")

#!/usr/bin/env python3
"""Map of RTK (tied) minus the lidar surfaces at every FIX point, 25 October 2026 survey.
Reads data/field/2026-10-25/rtk_points_all.csv (compile_rtk_2026-10-25.py): FIX points, B tied
through the CORS-11 = T1a-63 pair, B1 removed. Writes figures/rtk_2026-10-25_diffmap.png:
two panels on 2025 FSA imagery (MnGeo WMS), left minus the 2021 lidar ground (gen2), right minus
the 2008 lidar ground (gen1); one coloured marker per point, red = ground below the lidar
surface (lidar high), blue = ground above it; a point with no lidar value in that panel is drawn
as an open gray marker. Colour scale +/-0.5 m, values beyond it saturate (the scale is mine).
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

P = pd.read_csv("data/field/2026-10-25/rtk_points_all.csv"); F = P[P.status.str.startswith("FIX")].copy()
SUF = sys.argv[1] if len(sys.argv) > 1 else ""   # output-name suffix, so a replot need not clobber an earlier figure
X0, Y0, X1, Y1 = 579700, 4885300, 580100, 4885720
url = f"https://imageserver.gisdata.mn.gov/cgi-bin/wms?SERVICE=WMS&VERSION=1.1.1&REQUEST=GetMap&LAYERS=fsa2025&STYLES=&SRS=EPSG:26915&BBOX={X0},{Y0},{X1},{Y1}&WIDTH=1200&HEIGHT={int(1200*(Y1-Y0)/(X1-X0))}&FORMAT=image/jpeg"
try: img = np.asarray(Image.open(io.BytesIO(urllib.request.urlopen(url, timeout=120).read())))
except Exception as e: img = None; print("imagery not fetched:", e)
mk = {"hillslope_RS3": "o", "hillslope_T4a": "s", "floodplain_RSplus_A": "^", "floodplain_RSplus_B": "v", "floodplain_CORS_Wickert": "D"}
norm = TwoSlopeNorm(vmin=-0.5, vcenter=0, vmax=0.5)
fig, axs = plt.subplots(1, 2, figsize=(14, 7.2)); fig.subplots_adjust(left=0.01, right=0.99, top=0.9, bottom=0.01, wspace=0.03)
for ax, col, title in ((axs[0], "tied_minus_gen2_m", "2021 lidar"), (axs[1], "tied_minus_gen1_m", "2008 lidar")):
    if img is not None: ax.imshow(img, extent=[X0, X1, Y0, Y1], alpha=0.85)
    nan = F[col].isna()
    print(f"{col}: {int((~nan).sum())} points with a value, {int(nan.sum())} without")
    for g, s in F[nan].groupby("group"):
        if g not in mk: continue
        ax.scatter(s.easting_tied, s.northing_tied, facecolors="none", edgecolors="0.35", s=40, marker=mk[g], linewidths=0.8)
    for g, s in F[~nan].groupby("group"):
        if g not in mk: continue
        sc = ax.scatter(s.easting_tied, s.northing_tied, c=s[col], cmap="RdBu", norm=norm, s=46, marker=mk[g], edgecolors="black", linewidths=0.5, label={'floodplain_RSplus_B': 'RS+ B (T2a)', 'floodplain_RSplus_A': 'RS+ A (T1a)', 'floodplain_CORS_Wickert': 'CORS (Wickert)', 'hillslope_RS3': 'RS3 (T3a)', 'hillslope_T4a': 'T4a', 'benchmark_DR9925': 'DR9925'}[g])
    ax.set_xlim(X0, X1); ax.set_ylim(Y0, Y1); ax.set_xticks([]); ax.set_yticks([]); ax.set_title(title, fontsize=10)
    ax.plot([X0 + 20, X0 + 120], [Y0 + 15, Y0 + 15], "w-", lw=3); ax.text(X0 + 70, Y0 + 22, "100 m", color="white", ha="center", fontsize=9)
axs[0].legend(loc="upper right", fontsize=8); axs[1].scatter([], [], facecolors="none", edgecolors="0.35", s=40, label="no value"); axs[1].legend(loc="upper right", fontsize=8)
cb = fig.colorbar(sc, ax=axs, shrink=0.75, pad=0.01, extend="both"); cb.set_label("RTK minus lidar (m)")
fig.savefig(f"figures/rtk_2026-10-25_diffmap{SUF}.png", dpi=120, bbox_inches="tight"); print(f"figures/rtk_2026-10-25_diffmap{SUF}.png")

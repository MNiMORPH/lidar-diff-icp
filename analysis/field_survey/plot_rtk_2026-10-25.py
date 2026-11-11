#!/usr/bin/env python3
"""Figures for the 25 October 2026 RTK survey against the two lidar surfaces.
Reads data/field/2026-10-25/rtk_points_all.csv (compile_rtk_2026-10-25.py). Writes
  figures/rtk_2026-10-25_profiles{SUF}.png  the hillslope traverse (RS3 + T4a) and the floodplain
                                       loop (RS+ A/B + CORS) as RTK, gen2 and gen1 elevations
                                       in survey order, FIX points only
"""
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

P = pd.read_csv("data/field/2026-10-25/rtk_points_all.csv")
SUF = sys.argv[1] if len(sys.argv) > 1 else ""   # output-name suffix, so a replot need not clobber an earlier figure
F = P[P.status.str.startswith("FIX") & P.tied_minus_gen2_m.notna()].copy()
mk = {"hillslope_RS3": "o", "hillslope_T4a": "s", "floodplain_RSplus_A": "^", "floodplain_RSplus_B": "v", "floodplain_CORS_Wickert": "D", "benchmark_DR9925": "*"}

fig, axs = plt.subplots(2, 1, figsize=(12, 9))
for ax, groups in ((axs[0], ["hillslope_RS3", "hillslope_T4a"]), (axs[1], ["floodplain_RSplus_A", "floodplain_RSplus_B", "floodplain_CORS_Wickert"])):
    off = 0
    for g in groups:
        s = F[F.group == g].sort_values(["t_start"]); d = np.r_[0, np.cumsum(np.hypot(np.diff(s.easting_tied), np.diff(s.northing_tied)))] + off
        ax.plot(d, s.navd88_tied_m, "k-", lw=1); ax.scatter(d, s.navd88_tied_m, c="black", s=16, marker=mk[g], label="RTK" if g == groups[0] else None)
        ax.plot(d, s.gen2_z_m, color="#2a78d6", lw=1.2, label="2021 lidar" if g == groups[0] else None); ax.plot(d, s.gen1_z_m, color="#e34948", lw=1.2, label="2008 lidar" if g == groups[0] else None)
        off = d[-1] + 30
    ax.set_ylabel("NAVD 88 (m)"); ax.legend(fontsize=8); ax.grid(alpha=0.3)
axs[1].set_xlabel("distance along path (m)")
fig.tight_layout(); fig.savefig(f"figures/rtk_2026-10-25_profiles{SUF}.png", dpi=110); print(f"figures/rtk_2026-10-25_profiles{SUF}.png")

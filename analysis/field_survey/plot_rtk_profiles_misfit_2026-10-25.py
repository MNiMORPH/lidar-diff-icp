#!/usr/bin/env python3
"""Lidar-misfit companion to the transect figure (plot_rtk_2026-10-25.py, profiles panel): the same walked-path
axis and markers, but the y-axis is each lidar surface minus the tied RTK in mm (the lidar misfit). Points only: no lines,
no medians, no titles (Andy, 2026-11-11). Reads data/field/2026-10-25/rtk_points_all.csv. Optional argv[1] = output-name suffix.
Writes figures/rtk_2026-10-25_profiles_misfit{SUF}.png
"""
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

P = pd.read_csv(sys.argv[2] if len(sys.argv) > 2 else "data/field/2026-10-25/rtk_points_all.csv")   # argv[2]: another table with the same columns, e.g. rtk_points_all_fitted.csv
SUF = sys.argv[1] if len(sys.argv) > 1 else ""
F = P[P.status.str.startswith("FIX") & P.gen2_minus_rtk_m.notna()].copy()
mk = {"hillslope_RS3": "o", "hillslope_T4a": "s", "floodplain_RSplus_A": "^", "floodplain_RSplus_B": "v", "floodplain_CORS_Wickert": "D"}
C21, C08 = "#2a78d6", "#e34948"

fig, axs = plt.subplots(2, 1, figsize=(12, 9))
for ax, groups in ((axs[0], ["hillslope_RS3", "hillslope_T4a"]), (axs[1], ["floodplain_RSplus_A", "floodplain_RSplus_B", "floodplain_CORS_Wickert"])):
    off = 0
    for g in groups:
        s = F[F.group == g].sort_values(["t_start"]); d = np.r_[0, np.cumsum(np.hypot(np.diff(s.easting_tied), np.diff(s.northing_tied)))] + off
        m21, m08 = 1000 * s.gen2_minus_rtk_m.values, 1000 * s.gen1_minus_rtk_m.values
        ax.scatter(d, m21, c=C21, s=22, marker=mk[g], edgecolors="black", linewidths=0.4, label="2021 lidar" if g == groups[0] else None)
        ok = np.isfinite(m08); ax.scatter(d[ok], m08[ok], c=C08, s=22, marker=mk[g], edgecolors="black", linewidths=0.4, label="2008 lidar" if g == groups[0] else None)
        ax.axvspan(d[0], d[-1], color="0.93" if groups.index(g) % 2 else "white", zorder=0)
        off = d[-1] + 30
    ax.axhline(0, color="black", lw=0.8); ax.set_ylabel("lidar minus RTK (mm)"); ax.legend(fontsize=8, loc="lower right"); ax.grid(alpha=0.3)
axs[1].set_xlabel("distance along path (m)")
fig.tight_layout(); fig.savefig(f"figures/rtk_2026-10-25_profiles_misfit{SUF}.png", dpi=110); print(f"figures/rtk_2026-10-25_profiles_misfit{SUF}.png")

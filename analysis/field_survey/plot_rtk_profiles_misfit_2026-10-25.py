#!/usr/bin/env python3
"""Misfit companion to the transect figure (plot_rtk_2026-10-25.py, profiles panel): the same walked-path
axis and markers, but the y-axis is tied RTK minus each lidar surface in mm, so the offsets can be read
directly. Reads data/field/2026-10-25/rtk_points_all.csv. Optional argv[1] = output-name suffix.
Writes figures/rtk_2026-10-25_profiles_misfit{SUF}.png
"""
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

P = pd.read_csv("data/field/2026-10-25/rtk_points_all.csv")
SUF = sys.argv[1] if len(sys.argv) > 1 else ""
_b = P[P.base == "B"].tie_source.iloc[0]; _pc = P[P.antenna_corr_m != 0].groupby("group").antenna_corr_m.max()
TIE = f"B tied to {_b.split(' under B (')[0]} ({_b.split(' under B (')[1].split('; ')[0]}" + ("; pole heights corrected: " + ", ".join(f"{g} {v:+.3f} m" for g, v in _pc.items()) if len(_pc) else "") + ")"
F = P[P.status.str.startswith("FIX") & P.tied_minus_gen2_m.notna()].copy()
mk = {"hillslope_RS3": "o", "hillslope_T4a": "s", "floodplain_RSplus_A": "^", "floodplain_RSplus_B": "v", "floodplain_CORS_Wickert": "D"}
C21, C08 = "#2a78d6", "#e34948"

fig, axs = plt.subplots(2, 1, figsize=(12, 9))
for ax, groups, title in ((axs[0], ["hillslope_RS3", "hillslope_T4a"], "hillslope traverse west of the road (RS3 circles, T4a squares), survey order within group"),
                          (axs[1], ["floodplain_RSplus_A", "floodplain_RSplus_B", "floodplain_CORS_Wickert"], "floodplain field east of the road (RS+ A triangles up, RS+ B down, CORS diamonds)")):
    off = 0
    for g in groups:
        s = F[F.group == g].sort_values(["t_start"]); d = np.r_[0, np.cumsum(np.hypot(np.diff(s.easting_tied), np.diff(s.northing_tied)))] + off
        m21, m08 = 1000 * s.tied_minus_gen2_m.values, 1000 * s.tied_minus_gen1_m.values
        ax.plot(d, m21, color=C21, lw=1); ax.scatter(d, m21, c=C21, s=22, marker=mk[g], edgecolors="black", linewidths=0.4, label="RTK (tied) minus 2021 lidar (gen2)" if g == groups[0] else None)
        ok = np.isfinite(m08); ax.plot(d[ok], m08[ok], color=C08, lw=1); ax.scatter(d[ok], m08[ok], c=C08, s=22, marker=mk[g], edgecolors="black", linewidths=0.4, label="RTK (tied) minus 2008 lidar (gen1)" if g == groups[0] else None)
        for arr, c in ((m21, C21), (m08, C08)):
            if np.isfinite(arr).sum(): ax.hlines(np.nanmedian(arr), d[0], d[-1], color=c, lw=0.8, ls="--")
        ax.axvspan(d[0], d[-1], color="0.93" if groups.index(g) % 2 else "white", zorder=0)
        ax.annotate(f"{g}: n {int(np.isfinite(m21).sum())} / {int(ok.sum())}", (0.5 * (d[0] + d[-1]), 1.0), xycoords=("data", "axes fraction"), ha="center", va="bottom", fontsize=7)
        off = d[-1] + 30
    ax.axhline(0, color="black", lw=0.8); ax.set_ylabel("RTK (tied) minus lidar ground (mm)\nnegative = ground below the lidar surface"); ax.set_title(title, fontsize=9, pad=14); ax.legend(fontsize=8, loc="lower right"); ax.grid(alpha=0.3)
axs[1].set_xlabel("distance along the walked path (m), groups laid end to end; dashed = group median; FIX points only")
fig.suptitle(f"25 Oct 2026 RTK survey, misfits along the same transects as the elevation figure\n{TIE}", fontsize=9)
fig.tight_layout(); fig.savefig(f"figures/rtk_2026-10-25_profiles_misfit{SUF}.png", dpi=110); print(f"figures/rtk_2026-10-25_profiles_misfit{SUF}.png")

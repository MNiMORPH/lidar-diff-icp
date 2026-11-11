#!/usr/bin/env python3
"""Probability densities of RTK (tied) minus the lidar surfaces, 25 October 2026 survey.
Reads data/field/2026-10-25/rtk_points_all.csv (compile_rtk_2026-10-25.py): FIX points only,
B2 tied through the CORS-11 = T1a-63 pair (+0.473 m vertical, dE +0.802 dN -1.123 m), B1 removed.
Writes figures/rtk_2026-10-25_pdf.png and .pdf: two rows (floodplain field, hillslope bluff), two
columns (minus 2021 lidar gen2, minus 2008 lidar gen1); a Gaussian KDE per group over a stacked
rug of the points, with median and NMAD per group printed in the legend. The KDE bandwidth is
Scott's rule (scipy default); no cut is applied to the residuals.
"""
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde
from matplotlib.transforms import blended_transform_factory

P = pd.read_csv("data/field/2026-10-25/rtk_points_all.csv")
SUF = sys.argv[1] if len(sys.argv) > 1 else ""   # output-name suffix, so a replot need not clobber an earlier figure
_b2 = P[P.base == "B2"].tie_source.iloc[0]; _pc = P[P.antenna_corr_m != 0].groupby("group").antenna_corr_m.max()
TIE = f"B2 tied to {_b2.split(' under B2 (')[0]} ({_b2.split(' under B2 (')[1].split('; ')[0]}" + ("; pole heights corrected: " + ", ".join(f"{g} {v:+.3f} m" for g, v in _pc.items()) if len(_pc) else "") + ")"
F = P[P.status.str.startswith("FIX")].copy()
rows = (("floodplain", ["floodplain_RSplus_B", "floodplain_RSplus_A", "floodplain_CORS_Wickert"], "floodplain: the sward field east of the road (R1-v140)"),
        ("hillslope", ["hillslope_RS3", "hillslope_T4a"], "hillslope: the forested bluff west of the road"))
cols = (("tied_minus_gen2_m", "RTK minus 2021 lidar ground (gen2), mm"), ("tied_minus_gen1_m", "RTK minus 2008 lidar ground (gen1), mm"))
colour = {"floodplain_RSplus_B": "#2a78d6", "floodplain_RSplus_A": "#e3a008", "floodplain_CORS_Wickert": "#1a1a1a", "hillslope_RS3": "#2a78d6", "hillslope_T4a": "#e3a008"}
label = {"floodplain_RSplus_B": "T2a, RS+ (B2)", "floodplain_RSplus_A": "T1a/T1b, RS+ (B2)", "floodplain_CORS_Wickert": "Wickert, RS2+ (CORS, no tie)", "hillslope_RS3": "T3a/T3b/T3c, RS3 (B2)", "hillslope_T4a": "T4a/T4b (B2)"}
fig, axs = plt.subplots(2, 2, figsize=(12, 8), sharex="row")
print("row, group, column, n, median_mm, NMAD_mm")
for i, (rk, groups, title) in enumerate(rows):
    lo = min(F.loc[F.group.isin(groups), c].min() for c, _ in cols) * 1000 - 100
    hi = max(F.loc[F.group.isin(groups), c].max() for c, _ in cols) * 1000 + 100
    x = np.linspace(lo, hi, 800)
    for j, (c, xl) in enumerate(cols):
        ax = axs[i, j]; rug = 0
        for g in groups:
            v = (F.loc[F.group == g, c].dropna() * 1000).to_numpy()
            if len(v) == 0: continue
            med = np.median(v); nmad = 1.4826 * np.median(np.abs(v - med))
            print(f"{rk}, {g}, {c}, {len(v)}, {med:.1f}, {nmad:.1f}")
            lab = f"{label[g]}: n={len(v)}, median {med:+.0f}, NMAD {nmad:.0f} mm"
            if len(v) >= 3 and nmad > 0:
                ax.plot(x, gaussian_kde(v)(x), color=colour[g], lw=2, label=lab)
            else:
                ax.plot([], [], color=colour[g], lw=2, label=lab)
            ax.plot(v, np.full_like(v, 0.02 + rug * 0.035), "|", color=colour[g], ms=9, mew=1.2, alpha=0.8, transform=blended_transform_factory(ax.transData, ax.transAxes), clip_on=False); rug += 1
            ax.axvline(med, color=colour[g], lw=0.9, ls="--", alpha=0.7)
        ax.axvline(0, color="gray", lw=0.8); ax.set_ylim(bottom=0)
        ax.set_xlabel(xl); ax.set_ylabel("density (per mm)"); ax.legend(fontsize=7.5, loc="upper left"); ax.grid(alpha=0.25)
        ax.set_title(title if j == 0 else "", fontsize=9, loc="left")
fig.suptitle(f"25 Oct 2026 RTK survey, FIX points, {TIE}; negative = ground below the lidar surface", fontsize=9)
fig.tight_layout()
for ext in ("png", "pdf"):
    fig.savefig(f"figures/rtk_2026-10-25_pdf{SUF}.{ext}", dpi=120); print(f"figures/rtk_2026-10-25_pdf{SUF}.{ext}")

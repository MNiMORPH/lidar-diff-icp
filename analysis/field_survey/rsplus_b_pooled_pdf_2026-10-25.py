#!/usr/bin/env python
"""Pooled PDF of lidar return heights relative to the RTK ground, every RS+ B point, returns within 1 m.

Reads the per-return tables written by rsplus_b_returns_2026-10-25.py (benchmark-tied frame), keeps the rows with
radius_1m, pools the 100 points, and plots probability densities of return height minus the point's tied RTK
elevation for both flights: "ground" returns (2021 vendor class 2; 2008 vendor classes 2 and 8, 8 being model
key-points, also ground) as the heavy lines and all returns as the thin lines. Bins 0.01 m over -0.20 to +0.80 m (both MINE; the density is normalized over all returns, the clipped counts are in the log); a return that lies
within 1 m of two points counts once per point. Writes the pooled offsets, the statistics and the figure.
  inputs   data/field/2026-10-25/returns/rsplus_b_returns_2021.csv.gz, rsplus_b_returns_2008.csv
  outputs  data/field/2026-10-25/returns/rsplus_b_pooled_1m.csv (point, epoch, classification, dz_m)
           data/field/2026-10-25/returns/rsplus_b_pooled_1m_log.txt
           figures/rtk_2026-10-25_rsplusB_pooled_pdf_1m.png (FIGDPI, default 200)
"""
import os
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
ROOT = os.path.dirname(os.path.abspath(__file__)) + "/../.."; os.chdir(ROOT)
D = "data/field/2026-10-25/returns"; LOG = open(f"{D}/rsplus_b_pooled_1m_log.txt", "w")
def say(*a):
    s = " ".join(str(x) for x in a); print(s); LOG.write(s + "\n")
t21 = pd.read_csv(f"{D}/rsplus_b_returns_2021.csv.gz"); t08 = pd.read_csv(f"{D}/rsplus_b_returns_2008.csv")
a21 = t21[t21.radius_1m][["point", "classification", "dz_m"]].assign(epoch=2021)
a08 = t08[t08.radius_1m][["point", "classification", "dz_star_m"]].rename(columns={"dz_star_m": "dz_m"}).assign(epoch=2008)
pool = pd.concat([a21, a08], ignore_index=True); pool.to_csv(f"{D}/rsplus_b_pooled_1m.csv", index=False, float_format="%.4f")
say("inputs:", f"{D}/rsplus_b_returns_2021.csv.gz", f"{D}/rsplus_b_returns_2008.csv", "(benchmark-tied frame; 2008 in the product's star frame); rows with radius_1m; dz = return z minus the point's tied RTK elevation")
sets = {"2021 ground (class 2)": a21[a21.classification == 2].dz_m, "2021 all returns": a21.dz_m,
        "2008 ground (classes 2, 8)": a08[a08.classification.isin([2, 8])].dz_m, "2008 all returns": a08.dz_m}
def nmad(v): return 1.4826 * np.median(np.abs(v - np.median(v)))
say(f"{'set':28s} {'n':>6s} {'points':>6s} {'min':>7s} {'p05':>7s} {'p10':>7s} {'median':>7s} {'p90':>7s} {'p95':>7s} {'max':>7s} {'NMAD':>6s} {'<0.02':>6s}")
for k, v in sets.items():
    npts = pool[(pool.epoch == int(k[:4])) & (pool.dz_m.isin(v))].point.nunique()
    say(f"{k:28s} {len(v):6d} {npts:6d} {v.min():+7.3f} {v.quantile(.05):+7.3f} {v.quantile(.1):+7.3f} {v.median():+7.3f} {v.quantile(.9):+7.3f} {v.quantile(.95):+7.3f} {v.max():+7.3f} {nmad(v):6.3f} {int((v < 0.02).sum()):6d}")
say("2021 classes within 1 m:", a21.classification.value_counts().to_dict(), "| 2008 classes:", a08.classification.value_counts().to_dict())
say("per-point counts within 1 m, median (min-max): 2021 all", int(a21.groupby('point').size().median()), f"({a21.groupby('point').size().min()}-{a21.groupby('point').size().max()})",
    "| 2008 all", int(a08.groupby('point').size().median()), f"({a08.groupby('point').size().min()}-{a08.groupby('point').size().max()}); points with no 2008 return within 1 m:", 100 - a08.point.nunique())
XLIM = (-0.20, 0.80)   # plotted range, MINE: the 2021 "all returns" tail reaches +14 m (canopy); counts outside are logged
bins = np.arange(XLIM[0], XLIM[1] + 0.005, 0.01)
for k, v in sets.items(): say(f"   outside the plotted range {XLIM}: {k}: {int(((v < XLIM[0]) | (v > XLIM[1])).sum())} of {len(v)} returns (density normalized over all returns)")
fig, ax = plt.subplots(figsize=(8, 4.6))
style = {"2021 ground (class 2)": ("#0072b2", 1.8, "-"), "2021 all returns": ("#0072b2", 0.9, "--"), "2008 ground (classes 2, 8)": ("#d55e00", 1.8, "-"), "2008 all returns": ("#d55e00", 0.9, "--")}
for k, v in sets.items():
    c, lw, ls = style[k]; ax.hist(v, bins=bins, density=True, histtype="step", color=c, lw=lw, ls=ls, label=k)
ax.axvline(0, color="0.4", lw=0.8, ls=":"); ax.set_xlim(XLIM)
ax.set_xlabel("lidar return minus tied RTK ground (m)"); ax.set_ylabel("probability density (1/m)")
ax.legend(frameon=False, fontsize=9); ax.spines[["top", "right"]].set_visible(False); fig.tight_layout()
png = "figures/rtk_2026-10-25_rsplusB_pooled_pdf_1m.png"; fig.savefig(png, dpi=int(os.environ.get("FIGDPI", 200))); say("wrote", png)
LOG.close()

#!/usr/bin/env python
"""Pooled PDF of lidar return heights relative to the RTK ground, every RS+ B point, returns within 1 m.

Reads the per-return tables written by rsplus_b_returns_2026-10-25.py (benchmark-tied frame), keeps the rows with
radius_1m, pools the 100 points, and plots probability densities of return height minus the point's tied RTK
elevation for both flights, ALL returns with no ground-class cut (the log still tabulates the vendor ground classes, 2021
class 2 and 2008 classes 2 and 8, for reference). Bins 0.01 m over -0.20 to +0.80 m (both MINE; the density is normalized over all returns, the clipped counts are in the log); a return that lies
within 1 m of two points counts once per point. Writes the pooled offsets, the statistics and the figure.
  inputs   data/field/2026-10-25/returns/rsplus_b_returns_2021.csv.gz, rsplus_b_returns_2008.csv
  outputs  data/field/2026-10-25/returns/rsplus_b_pooled_1m.csv (point, epoch, classification, dz_m)
           data/field/2026-10-25/returns/rsplus_b_pooled_1m_log.txt
           figures/rtk_2026-10-25_rsplusB_pooled_pdf_1m.png and _1m_2021.png (2021 alone) (FIGDPI, default 200)
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
# all returns only, no ground-class cut (Andy, 2026-11-20: "eschew the ground-return class cutoff")
for k, c in (("2021 all returns", "#0072b2"), ("2008 all returns", "#d55e00")):
    ax.hist(sets[k], bins=bins, density=True, histtype="step", color=c, lw=1.6, label=k[:4] + " lidar, all returns")
ax.axvline(0, color="0.4", lw=0.8, ls=":"); ax.set_xlim(XLIM)
ax.set_xlabel("lidar return minus tied RTK ground (m)"); ax.set_ylabel("probability density (1/m)")
ax.legend(frameon=False, fontsize=9); ax.spines[["top", "right"]].set_visible(False); fig.tight_layout()
png = "figures/rtk_2026-10-25_rsplusB_pooled_pdf_1m.png"; fig.savefig(png, dpi=int(os.environ.get("FIGDPI", 200))); say("wrote", png)
# 2021 alone (Andy, 2026-11-20): same bins, range cut to -0.10..+0.40 m (MINE) so the peak fills the axis
X2 = (-0.10, 0.40); b2 = np.arange(X2[0], X2[1] + 0.005, 0.01)
fig, ax = plt.subplots(figsize=(8, 4.6))
ax.hist(sets["2021 all returns"], bins=b2, density=True, histtype="stepfilled", color="#0072b2", alpha=0.18, lw=0)
ax.hist(sets["2021 all returns"], bins=b2, density=True, histtype="step", color="#0072b2", lw=1.6, label="2021 lidar, all returns")
ax.axvline(0, color="0.4", lw=0.8, ls=":"); ax.set_xlim(X2)
ax.set_xlabel("2021 lidar return minus tied RTK ground (m)"); ax.set_ylabel("probability density (1/m)")
ax.legend(frameon=False, fontsize=9); ax.spines[["top", "right"]].set_visible(False); fig.tight_layout()
png2 = "figures/rtk_2026-10-25_rsplusB_pooled_pdf_1m_2021.png"; fig.savefig(png2, dpi=int(os.environ.get("FIGDPI", 200))); say("wrote", png2, f"(range {X2}; 2021 all returns outside it: {int(((sets['2021 all returns'] < X2[0]) | (sets['2021 all returns'] > X2[1])).sum())} of {len(sets['2021 all returns'])})")

# ---- all 2021 returns, canopy tail included (Andy, 2026-11-20: "nonuniform bins if needed, and/or a method for
# sharing data density more smoothly"). Left: density per meter with bins 0.01 m wide from -0.12 to +0.30 m and then
# geometric (ratio 1.15) to the top return, on a symlog height axis (linear within 0.3 m) and a log density axis.
# Right: the empirical cumulative distribution, no bins. Edges and ratio are MINE.
v_all = sets["2021 all returns"].to_numpy(); v_g = sets["2021 ground (class 2)"].to_numpy()
edges = np.concatenate([np.arange(-0.12, 0.30, 0.01), 0.30 * 1.15 ** np.arange(0, 40)]); edges = edges[edges <= v_all.max() * 1.15]
if edges[-1] < v_all.max(): edges = np.append(edges, edges[-1] * 1.15)
fig, (a1, a2) = plt.subplots(1, 2, figsize=(12, 4.8))
for v, lw, ls, lab in ((v_all, 1.6, "-", "2021 lidar, all returns"),):
    c, _ = np.histogram(v, bins=edges); dens = c / (len(v) * np.diff(edges)); dens = np.where(dens > 0, dens, np.nan)
    a1.stairs(dens, edges, color="#0072b2", lw=lw, ls=ls, label=lab, baseline=None)
    a2.plot(np.sort(v), np.arange(1, len(v) + 1) / len(v), color="#0072b2", lw=lw, ls=ls, label=lab)
for ax in (a1, a2):
    ax.set_xscale("symlog", linthresh=0.3, linscale=2.5); ax.set_xlim(-0.12, edges[-1]); ax.axvline(0, color="0.4", lw=0.8, ls=":")
    ax.set_xlabel("2021 lidar return minus tied RTK ground (m)"); ax.spines[["top", "right"]].set_visible(False)
    ax.set_xticks([-0.1, 0, 0.1, 0.2, 0.3, 1, 3, 10]); ax.set_xticklabels(["-0.1", "0", "0.1", "0.2", "0.3", "1", "3", "10"])
a1.set_yscale("log"); a1.set_ylabel("probability density (1/m)"); a2.set_ylabel("fraction of returns below"); a2.set_ylim(0, 1)
a1.legend(frameon=False, fontsize=9, loc="upper left"); fig.tight_layout()
png3 = "figures/rtk_2026-10-25_rsplusB_pooled_1m_2021_alltail.png"; fig.savefig(png3, dpi=int(os.environ.get("FIGDPI", 200))); say("wrote", png3)
say("2021 within 1 m, fraction of returns below a height (all returns | class 2):")
for h in (0.0, 0.02, 0.05, 0.10, 0.15, 0.20, 0.30, 0.50, 1.0, 2.0, 5.0):
    say(f"   {h:5.2f} m  {100 * (v_all < h).mean():6.2f} %  |  {100 * (v_g < h).mean():6.2f} %")
say(f"   bins: {len(edges) - 1}; widest {np.diff(edges).max():.2f} m, narrowest {np.diff(edges).min():.3f} m; top edge {edges[-1]:.2f} m")
LOG.close()

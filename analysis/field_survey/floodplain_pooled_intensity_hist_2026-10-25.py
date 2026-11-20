#!/usr/bin/env python
"""The pooled 1 m return-height histogram of the 2021 lidar at the RTK points, each bar colored by the median
intensity of the returns in it (Andy, 2026-11-20: do the returns that reach the ground differ in intensity from
those that do not?). Returns within 1 m of the floodplain FIX points (per-return tables, benchmark-tied frame),
all returns, 1 cm bins from -0.10 to +0.40 m (MINE; returns above +0.40 m, the canopy, are counted in the density
normalization and listed in the log). Two panels: all 147 floodplain points pooled, and the 98 RS+ B points alone. A
third panel: median intensity per bin with the interquartile band, for the same bins. Log: intensity of the returns
within 5 cm of the ground against those in the sward peak (0.10-0.20 m), per rover, with a Mann-Whitney U test.
Outputs  figures/rtk_2026-10-25_floodplain_pooled_1m_2021_intensity_hist.png; data/field/2026-10-25/returns/floodplain_pooled_intensity_hist_log.txt
"""
import os
import numpy as np, pandas as pd
from scipy import stats
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
ROOT = os.path.dirname(os.path.abspath(__file__)) + "/../.."; os.chdir(ROOT)
D = "data/field/2026-10-25/returns"; LOG = open(f"{D}/floodplain_pooled_intensity_hist_log.txt", "w")
def say(*a):
    s = " ".join(str(x) for x in a); print(s); LOG.write(s + "\n")
P = pd.read_csv("data/field/2026-10-25/rtk_points_all.csv")
parts = []
for tag, group in (("rsplus_b", "floodplain_RSplus_B"), ("rsplus_a", "floodplain_RSplus_A"), ("cors", "floodplain_CORS_Wickert")):
    t = pd.read_csv(f"{D}/{tag}_returns_2021.csv.gz"); t = t[t.radius_1m].merge(P[(P.group == group) & (P.status == "FIX")][["point"]], on="point"); parts.append(t.assign(group=group))
T = pd.concat(parts, ignore_index=True)
X0, X1, W = -0.10, 0.40, 0.01; edges = np.arange(X0, X1 + W / 2, W); mids = 0.5 * (edges[1:] + edges[:-1])
sets = {"all floodplain points (147)": T, "RS+ B (98)": T[T.group == "floodplain_RSplus_B"]}
say(f"returns within 1 m: {len(T)} over {T.groupby(['group', 'point']).ngroups} FIX points; intensity 14-bit units")
for name, q in sets.items():
    g = q[q.dz_m < 0.05].intensity; s = q[(q.dz_m >= 0.10) & (q.dz_m < 0.20)].intensity; c = q[q.dz_m > 0.40].intensity
    u = stats.mannwhitneyu(g, s, alternative="two-sided")
    say(f"{name}: ground returns (< 0.05 m) n {len(g)} intensity median {g.median():.0f} IQR {g.quantile(.25):.0f}-{g.quantile(.75):.0f} | sward peak (0.10-0.20 m) n {len(s)} median {s.median():.0f} IQR {s.quantile(.25):.0f}-{s.quantile(.75):.0f} | canopy (> 0.40 m) n {len(c)} median {c.median():.0f}; Mann-Whitney ground vs sward p {u.pvalue:.2g}; returns above +0.40 m not drawn: {len(c)}")
say("per rover, ground vs sward-peak median intensity (single returns only, which carry the whole pulse energy):")
for grp, q in T[T.number_of_returns == 1].groupby("group"):
    g = q[q.dz_m < 0.05].intensity; s = q[(q.dz_m >= 0.10) & (q.dz_m < 0.20)].intensity
    say(f"   {grp.replace('floodplain_', ''):16s} ground n {len(g):4d} median {g.median() if len(g) else np.nan:6.0f} | sward n {len(s):5d} median {s.median():6.0f}" + (f" | p {stats.mannwhitneyu(g, s).pvalue:.2g}" if len(g) > 5 else ""))
fig, axs = plt.subplots(1, 3, figsize=(19, 5)); norm = Normalize(vmin=2500, vmax=5000); cmap = plt.get_cmap("viridis")
for ax, (name, q) in zip(axs[:2], sets.items()):
    b = np.digitize(q.dz_m, edges) - 1; dens = np.array([(b == i).sum() for i in range(len(mids))]) / (len(q) * W); med = np.array([q.intensity[b == i].median() if (b == i).sum() >= 5 else np.nan for i in range(len(mids))])
    ax.bar(mids, dens, width=W, color=[cmap(norm(m)) if np.isfinite(m) else "0.85" for m in med], edgecolor="none"); ax.axvline(0, color="0.4", lw=0.8, ls=":")
    ax.set_xlim(X0, X1); ax.set_xlabel("2021 lidar return minus tied RTK ground (m)"); ax.set_ylabel("probability density (1/m)"); ax.set_title(name, fontsize=10); ax.spines[["top", "right"]].set_visible(False)
sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap); sm.set_array([]); fig.colorbar(sm, ax=axs[:2], shrink=0.8, label="median intensity of the returns in the bar (bars with < 5 returns gray)")
ax = axs[2]; q = T; b = np.digitize(q.dz_m, edges) - 1
med = np.array([q.intensity[b == i].median() if (b == i).sum() >= 5 else np.nan for i in range(len(mids))]); lo = np.array([q.intensity[b == i].quantile(.25) if (b == i).sum() >= 5 else np.nan for i in range(len(mids))]); hi = np.array([q.intensity[b == i].quantile(.75) if (b == i).sum() >= 5 else np.nan for i in range(len(mids))])
ax.fill_between(mids, lo, hi, color="#0072b2", alpha=0.2, lw=0); ax.plot(mids, med, color="#0072b2", lw=1.6); ax.axvline(0, color="0.4", lw=0.8, ls=":")
ax.set_xlim(X0, X1); ax.set_xlabel("2021 lidar return minus tied RTK ground (m)"); ax.set_ylabel("intensity, median and interquartile band per 1 cm bin"); ax.set_title("all floodplain points", fontsize=10); ax.spines[["top", "right"]].set_visible(False)
fig.savefig("figures/rtk_2026-10-25_floodplain_pooled_1m_2021_intensity_hist.png", dpi=int(os.environ.get("FIGDPI", 170)), bbox_inches="tight"); say("wrote figures/rtk_2026-10-25_floodplain_pooled_1m_2021_intensity_hist.png")
LOG.close()

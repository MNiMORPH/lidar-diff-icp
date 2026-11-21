#!/usr/bin/env python
"""Is the surveyed ground on the bluff systematically below the 2021 lidar floor? (Andy, 2026-11-21.) Per hillslope FIX
point, from hillslope_points_summary.csv: the lidar floor within 1 m minus RTK (slope-normal) and the official DEM minus
RTK (slope-normal). Median with a 20,000-draw bootstrap 95 % interval, sign test and Wilcoxon test, overall and per
traverse; Spearman with slope, canopy share, ground share (share of the point's 1 m returns within 5 cm of the RTK),
formal RMS, tilt; the outliers in both directions with their raw receiver records.
Outputs  data/field/2026-10-25/intensity/hillslope_rtk_vs_floor_log.txt; figures/rtk_2026-10-25_hillslope_rtk_vs_floor.png
"""
import os
import numpy as np, pandas as pd
from scipy import stats
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
ROOT = os.path.dirname(os.path.abspath(__file__)) + "/../.."; os.chdir(ROOT)
OUT = "data/field/2026-10-25/intensity"; LOG = open(f"{OUT}/hillslope_rtk_vs_floor_log.txt", "w")
def say(*a):
    s = " ".join(str(x) for x in a); print(s); LOG.write(s + "\n")
S = pd.read_csv(f"{OUT}/hillslope_points_summary.csv"); P = pd.read_csv("data/field/2026-10-25/rtk_points_all.csv")
S = S.merge(P[["group", "point", "gen2_minus_rtk_m", "rms_e", "rms_n"]], on=["group", "point"]); S["team"] = S.group.map({"hillslope_RS3": "first RS3 traverse (T3)", "hillslope_T4a": "second RS3 traverse (T4)"})
def boot_ci(v, n=20000, seed=0):
    r = np.random.default_rng(seed); v = np.asarray(v); return np.percentile([np.median(r.choice(v, len(v))) for _ in range(n)], [2.5, 97.5])
for lab, col in (("lidar floor within 1 m minus RTK (slope-normal)", "med_near"), ("official 1 m DEM minus RTK (slope-normal)", "dem_minus_rtk_n")):
    v = S[col].dropna(); ci = boot_ci(v); pos = int((v > 0).sum())
    say(f"{lab}: n {len(v)}, median {v.median():+.3f} m (95 % {ci[0]:+.3f}..{ci[1]:+.3f}), mean {v.mean():+.3f}, NMAD {1.4826 * np.median(np.abs(v - v.median())):.3f}; lidar above RTK at {pos} of {len(v)} points, sign test p {stats.binomtest(pos, len(v), 0.5).pvalue:.3g}, Wilcoxon p {stats.wilcoxon(v).pvalue:.3g}")
    for t, q in S.groupby("team"):
        w = q[col].dropna(); ci = boot_ci(w); pos = int((w > 0).sum()); say(f"   {t}: n {len(w)}, median {w.median():+.3f} (95 % {ci[0]:+.3f}..{ci[1]:+.3f}), above at {pos} of {len(w)}, sign p {stats.binomtest(pos, len(w), 0.5).pvalue:.3g}")
    w = v.drop(v.abs().nlargest(2).index); say(f"   without the two largest |offsets|: median {w.median():+.3f} (95 % {boot_ci(w)[0]:+.3f}..{boot_ci(w)[1]:+.3f})")
say("Spearman of DEM minus RTK with: " + ", ".join(f"{c} {S[['dem_minus_rtk_n', c]].dropna().corr(method='spearman').iloc[0, 1]:+.2f}" for c in ("slope", "canopy_share", "ground_share", "rms_v", "tilt")))
say("points with no 1 m return within 5 cm of the RTK: " + f"{int((S.ground_share == 0).sum())} of {len(S)}; their DEM minus RTK median {S[S.ground_share == 0].dem_minus_rtk_n.median():+.3f}; points with ground share > 0.2: {int((S.ground_share > 0.2).sum())}, median {S[S.ground_share > 0.2].dem_minus_rtk_n.median():+.3f}")
say("\noutliers (see the per-point table in hillslope_discrepant_log.txt): RS3 10 at 14:49, RTK 1.42 m below the DEM floor with no 1 m return lower than +1.12 m and a clean receiver record (26 samples, PDOP 1.9, 27 satellites, RMS 0.019), neighbors 9 and 11 at +0.05 and +0.05: a position displaced ~2.8 m along a 28 deg slope, or a cavity the laser never entered; RS3 24 at 15:12, RTK 0.21 m above the floor, recorded on GPS+GLONASS only (Galileo 0, BeiDou 0) with the base columns lost; T4a 23 at 15:13, RTK 0.18 m above the floor, full constellation, PDOP 1.4, the pole base on a rock, log or root being the plain reading.")
fig, axs = plt.subplots(1, 2, figsize=(14, 5), gridspec_kw=dict(width_ratios=[1.6, 1]))
ax = axs[0]
for t, c in (("first RS3 traverse (T3)", "#0072b2"), ("second RS3 traverse (T4)", "#d55e00")):
    q = S[S.team == t].sort_values("point"); ax.plot(q.point, q.dem_minus_rtk_n, "o-", color=c, ms=5, lw=0.8, label=t); ax.plot(q.point, q.med_near, "s", color=c, ms=4, mfc="none")
ax.axhline(0, color="0.4", lw=0.8, ls=":"); ax.set_ylim(-0.3, 0.3); ax.set_xlabel("point number along the traverse"); ax.set_ylabel("lidar above the surveyed ground, slope-normal (m)"); ax.set_title("filled: official DEM minus RTK; open: lidar floor within 1 m minus RTK (RS3 10 at +1.42 off the axis)", fontsize=9); ax.legend(frameon=False, fontsize=8); ax.spines[["top", "right"]].set_visible(False)
ax = axs[1]
for k, (t, c) in enumerate((("first RS3 traverse (T3)", "#0072b2"), ("second RS3 traverse (T4)", "#d55e00"))):
    v = S[S.team == t].dem_minus_rtk_n.dropna(); v = v[v.abs() < 0.5]; ci = boot_ci(v); ax.plot(np.full(len(v), k) + np.random.default_rng(1).uniform(-0.12, 0.12, len(v)), v, "o", color=c, ms=5, alpha=0.7); ax.plot([k - 0.25, k + 0.25], [v.median()] * 2, color="k", lw=2); ax.plot([k, k], ci, color="k", lw=1)
ax.axhline(0, color="0.4", lw=0.8, ls=":"); ax.set_xticks([0, 1]); ax.set_xticklabels(["first traverse", "second traverse"]); ax.set_ylabel("official DEM minus RTK, slope-normal (m)"); ax.set_title("median and bootstrap 95 % interval", fontsize=9); ax.spines[["top", "right"]].set_visible(False)
fig.tight_layout(); fig.savefig("figures/rtk_2026-10-25_hillslope_rtk_vs_floor.png", dpi=int(os.environ.get("FIGDPI", 150))); say("wrote figures/rtk_2026-10-25_hillslope_rtk_vs_floor.png")
LOG.close()

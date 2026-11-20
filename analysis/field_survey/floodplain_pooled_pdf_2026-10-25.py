#!/usr/bin/env python
"""Pooled densities of lidar return height above the tied RTK ground, returns within 1 m of each floodplain point,
by rover and by distance from the parking lot, ALL returns (no ground-class cut; the vendor classes are tabulated in
the log for reference).

Sets (FIX points only; the 100 m cut is MINE):
  RS+ B                 the 98 FIX points of the RS+ B loop, 66-352 m from the lot center (E 579979, N 4885340)
  RS+ A near lot        RS+ A FIX points within 100 m of the lot center
  RS+ A far             RS+ A FIX points beyond 100 m
  CORS                  the 17 CORS points, 27-68 m from the lot center
Inputs   data/field/2026-10-25/returns/<tag>_returns_{2021.csv.gz,2008.csv} (floodplain_returns_2026-10-25.py; benchmark-tied
         frame, 2008 in the product's star frame), rtk_points_all.csv for status and position
Outputs  data/field/2026-10-25/returns/floodplain_pooled_1m.csv (set, point, epoch, classification, dz_m) and _log.txt
         figures/rtk_2026-10-25_rsplusB_pooled_pdf_1m.png            RS+ B, both flights, 1 cm bins, -0.20..+0.80 m
         figures/rtk_2026-10-25_rsplusB_pooled_pdf_1m_2021.png       RS+ B, 2021 alone, -0.10..+0.40 m
         figures/rtk_2026-10-25_rsplusB_pooled_1m_2021_alltail.png   RS+ B, 2021, nonuniform bins on a symlog axis + ECDF
         figures/rtk_2026-10-25_floodplain_pooled_1m_2021_sets.png   the four sets, 2021, symlog density + ECDF
         figures/rtk_2026-10-25_floodplain_pooled_1m_2021_sets_linear.png  the four sets, 2021, 1 cm bins, -0.10..+0.50 m
         figures/rtk_2026-10-25_floodplain_pooled_1m_2008_sets.png   the four sets, 2008, symlog density + ECDF
         figures/rtk_2026-10-25_floodplain_pooled_1m_{2021,2008}_bins.png and _2021_bins_linear.png
             the same by LOCATION bins, which come from the spatially connected clustering of the per-point 2008
             offset in floodplain_offset_clusters_2026-10-25.py (returns/floodplain_offset_points.csv, cluster_med08_k4):
             "south zone" = clusters 1 and 2 (the lot end and tree line, N < 4885455; cluster 2 is one CORS point under
             a tree, merged by place), "swale stretch" = cluster 3 (14 RS+ B points on the west leg, N 4885496-4885579),
             "open field" = cluster 0 (the rest). The 2021 offset has no spatial structure to cluster on (flat
             semivariogram); the 2008 offset has, with a range near 60 m, and these are its clusters. Adopting them as
             bins is a PROPOSAL.
Bins 0.01 m in the linear figures; in the symlog figures 0.01 m from -0.12 to +0.30 m then geometric (ratio 1.15); the
height axis is linear within 0.3 m of zero and logarithmic beyond; the density is per meter (count / (n * bin width)).
All of these are MINE. FIGDPI sets the raster resolution (default 200).
"""
import os
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
ROOT = os.path.dirname(os.path.abspath(__file__)) + "/../.."; os.chdir(ROOT)
D = "data/field/2026-10-25/returns"; LOG = open(f"{D}/floodplain_pooled_1m_log.txt", "w")
def say(*a):
    s = " ".join(str(x) for x in a); print(s); LOG.write(s + "\n")
DPI = int(os.environ.get("FIGDPI", 200)); LOT = (579979.0, 4885340.0); NEAR = 100.0
P = pd.read_csv("data/field/2026-10-25/rtk_points_all.csv"); P["d_lot"] = np.hypot(P.easting_tied - LOT[0], P.northing_tied - LOT[1])
def load(tag, group):
    t21 = pd.read_csv(f"{D}/{tag}_returns_2021.csv.gz"); t08 = pd.read_csv(f"{D}/{tag}_returns_2008.csv")
    a = t21[t21.radius_1m][["point", "classification", "dz_m"]].assign(epoch=2021)
    b = t08[t08.radius_1m][["point", "classification", "dz_star_m"]].rename(columns={"dz_star_m": "dz_m"}).assign(epoch=2008)
    q = pd.concat([a, b], ignore_index=True).merge(P[P.group == group][["point", "status", "d_lot"]], on="point", how="left")
    return q[q.status == "FIX"].drop(columns="status")
rb = load("rsplus_b", "floodplain_RSplus_B").assign(set="RS+ B")
ra = load("rsplus_a", "floodplain_RSplus_A"); ra["set"] = np.where(ra.d_lot < NEAR, "RS+ A near lot", "RS+ A far")
rc = load("cors", "floodplain_CORS_Wickert").assign(set="CORS")
pool = pd.concat([rb, ra, rc], ignore_index=True); pool.to_csv(f"{D}/floodplain_pooled_1m.csv", index=False, float_format="%.4f")
ORDER = ["RS+ B", "RS+ A near lot", "RS+ A far", "CORS"]; COL = {"RS+ B": "#0072b2", "RS+ A near lot": "#d55e00", "RS+ A far": "#e69f00", "CORS": "k"}
# location bins from the offset clustering
OP = pd.read_csv(f"{D}/floodplain_offset_points.csv")[["group", "point", "cluster_med08_k4"]]
BINNAME = {0: "open field", 1: "south zone (lot, tree line)", 2: "south zone (lot, tree line)", 3: "swale stretch"}
OP["loc_bin"] = OP.cluster_med08_k4.map(BINNAME)
for df, g in ((rb, "floodplain_RSplus_B"), (ra, "floodplain_RSplus_A"), (rc, "floodplain_CORS_Wickert")):
    df["loc_bin"] = df.point.map(OP[OP.group == g].set_index("point").loc_bin)
BORDER = ["south zone (lot, tree line)", "open field", "swale stretch"]; BCOL = {"south zone (lot, tree line)": "#d55e00", "open field": "#0072b2", "swale stretch": "#009e73"}
pool = pd.concat([rb, ra, rc], ignore_index=True); pool.to_csv(f"{D}/floodplain_pooled_1m.csv", index=False, float_format="%.4f")
say("inputs: <tag>_returns tables in", D, "(benchmark-tied frame; 2008 star frame); rows with radius_1m; FIX points only; lot center", LOT, f"; near = within {NEAR:.0f} m (MINE)")
def nmad(v): return 1.4826 * np.median(np.abs(v - np.median(v)))
say(f"{'set':16s} {'epoch':>5s} {'n':>6s} {'points':>6s} {'d_lot m':>12s} {'min':>7s} {'p05':>7s} {'p10':>7s} {'median':>7s} {'p90':>7s} {'p95':>7s} {'max':>7s} {'NMAD':>6s} {'<0.02':>6s} {'<0.05':>6s}")
for st in ORDER:
    for ep in (2021, 2008):
        q = pool[(pool.set == st) & (pool.epoch == ep)]; v = q.dz_m
        if not len(v): say(f"{st:16s} {ep:5d}      0"); continue
        dl = P[(P.status == "FIX")].merge(q[["point"]].drop_duplicates(), on="point"); 
        say(f"{st:16s} {ep:5d} {len(v):6d} {q.point.nunique():6d} {q.d_lot.min():5.0f}-{q.d_lot.max():<5.0f} {v.min():+7.3f} {v.quantile(.05):+7.3f} {v.quantile(.1):+7.3f} {v.median():+7.3f} {v.quantile(.9):+7.3f} {v.quantile(.95):+7.3f} {v.max():+7.3f} {nmad(v):6.3f} {int((v < 0.02).sum()):6d} {int((v < 0.05).sum()):6d}")
say("vendor classes within 1 m, for reference:", {st: {ep: pool[(pool.set == st) & (pool.epoch == ep)].classification.value_counts().to_dict() for ep in (2021, 2008)} for st in ORDER})
say("2021, fraction of returns below a height, by set:")
hs = (0.0, 0.02, 0.05, 0.10, 0.15, 0.20, 0.30, 0.50, 1.0, 2.0)
say("   height " + " ".join(f"{st:>15s}" for st in ORDER))
for h in hs:
    say(f"   {h:5.2f} m " + " ".join(f"{100 * (pool[(pool.set == st) & (pool.epoch == 2021)].dz_m < h).mean():14.2f}%" for st in ORDER))

def symlog_axes(ax, top):
    ax.set_xscale("symlog", linthresh=0.3, linscale=2.5); ax.set_xlim(-0.12, top); ax.axvline(0, color="0.4", lw=0.8, ls=":")
    ax.set_xticks([-0.1, 0, 0.1, 0.2, 0.3, 1, 3, 10]); ax.set_xticklabels(["-0.1", "0", "0.1", "0.2", "0.3", "1", "3", "10"]); ax.spines[["top", "right"]].set_visible(False)
def edges_for(vmax, lin_bin=0.01, ratio=1.15):
    e = np.concatenate([np.arange(-0.12, 0.30, lin_bin), 0.30 * ratio ** np.arange(0, 60)]); e = e[e <= vmax * ratio]
    return np.append(e, e[-1] * ratio) if e[-1] < vmax else e
def density_ecdf(series, labels, colors, lws, xlabel, png, lin_bin=0.01, ratio=1.15):
    vmax = max(v.max() for v in series); edges = edges_for(vmax, lin_bin, ratio)
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(12, 4.8))
    for v, lab, c, lw in zip(series, labels, colors, lws):
        cnt, _ = np.histogram(v, bins=edges); dens = cnt / (len(v) * np.diff(edges)); dens = np.where(dens > 0, dens, np.nan)
        a1.stairs(dens, edges, color=c, lw=lw, label=lab, baseline=None); a2.plot(np.sort(v), np.arange(1, len(v) + 1) / len(v), color=c, lw=lw, label=lab)
    for ax in (a1, a2): symlog_axes(ax, edges[-1]); ax.set_xlabel(xlabel)
    a1.set_yscale("log"); a1.set_ylabel("probability density (1/m)"); a2.set_ylabel("fraction of returns below"); a2.set_ylim(0, 1)
    a1.legend(frameon=False, fontsize=9, loc="upper left"); fig.tight_layout(); fig.savefig(png, dpi=DPI); plt.close(fig); say("wrote", png)

# ---- RS+ B figures (as before, FIX points)
b21 = rb[rb.epoch == 2021].dz_m; b08 = rb[rb.epoch == 2008].dz_m
X1 = (-0.20, 0.80); bins = np.arange(X1[0], X1[1] + 0.005, 0.01)
fig, ax = plt.subplots(figsize=(8, 4.6))
ax.hist(b21, bins=bins, density=True, histtype="step", color="#0072b2", lw=1.6, label="2021 lidar, all returns")
ax.hist(b08, bins=bins, density=True, histtype="step", color="#d55e00", lw=1.6, label="2008 lidar, all returns")
ax.axvline(0, color="0.4", lw=0.8, ls=":"); ax.set_xlim(X1); ax.set_xlabel("lidar return minus tied RTK ground (m)"); ax.set_ylabel("probability density (1/m)")
ax.legend(frameon=False, fontsize=9); ax.spines[["top", "right"]].set_visible(False); fig.tight_layout()
fig.savefig("figures/rtk_2026-10-25_rsplusB_pooled_pdf_1m.png", dpi=DPI); plt.close(fig); say("wrote figures/rtk_2026-10-25_rsplusB_pooled_pdf_1m.png", f"(2021 outside {X1}: {int(((b21 < X1[0]) | (b21 > X1[1])).sum())} of {len(b21)}; 2008: {int(((b08 < X1[0]) | (b08 > X1[1])).sum())} of {len(b08)})")
X2 = (-0.10, 0.40); b2 = np.arange(X2[0], X2[1] + 0.005, 0.01)
fig, ax = plt.subplots(figsize=(8, 4.6))
ax.hist(b21, bins=b2, density=True, histtype="stepfilled", color="#0072b2", alpha=0.18, lw=0); ax.hist(b21, bins=b2, density=True, histtype="step", color="#0072b2", lw=1.6, label="2021 lidar, all returns")
ax.axvline(0, color="0.4", lw=0.8, ls=":"); ax.set_xlim(X2); ax.set_xlabel("2021 lidar return minus tied RTK ground (m)"); ax.set_ylabel("probability density (1/m)")
ax.legend(frameon=False, fontsize=9); ax.spines[["top", "right"]].set_visible(False); fig.tight_layout()
fig.savefig("figures/rtk_2026-10-25_rsplusB_pooled_pdf_1m_2021.png", dpi=DPI); plt.close(fig); say("wrote figures/rtk_2026-10-25_rsplusB_pooled_pdf_1m_2021.png", f"(outside {X2}: {int(((b21 < X2[0]) | (b21 > X2[1])).sum())} of {len(b21)})")
density_ecdf([b21.to_numpy()], ["2021 lidar, all returns"], ["#0072b2"], [1.6], "2021 lidar return minus tied RTK ground (m)", "figures/rtk_2026-10-25_rsplusB_pooled_1m_2021_alltail.png")

# ---- the four sets
for ep, png in ((2021, "figures/rtk_2026-10-25_floodplain_pooled_1m_2021_sets.png"), (2008, "figures/rtk_2026-10-25_floodplain_pooled_1m_2008_sets.png")):
    ser = [pool[(pool.set == st) & (pool.epoch == ep)].dz_m.to_numpy() for st in ORDER]
    # 2008 has 23-220 returns per set within 1 m: 0.05 m linear bins and a 1.3 ratio there (MINE), 0.01 m and 1.15 for 2021
    density_ecdf(ser, [f"{st} (n {len(v)})" for st, v in zip(ORDER, ser)], [COL[st] for st in ORDER], [1.6, 1.6, 1.2, 1.2], f"{ep} lidar return minus tied RTK ground (m)", png, **({"lin_bin": 0.05, "ratio": 1.3} if ep == 2008 else {}))
X3 = (-0.10, 0.50); b3 = np.arange(X3[0], X3[1] + 0.005, 0.01)
fig, ax = plt.subplots(figsize=(8, 4.6))
for st, lw in zip(ORDER, (1.6, 1.6, 1.2, 1.2)):
    v = pool[(pool.set == st) & (pool.epoch == 2021)].dz_m; ax.hist(v, bins=b3, density=True, histtype="step", color=COL[st], lw=lw, label=f"{st} (n {len(v)})")
ax.axvline(0, color="0.4", lw=0.8, ls=":"); ax.set_xlim(X3); ax.set_xlabel("2021 lidar return minus tied RTK ground (m)"); ax.set_ylabel("probability density (1/m)")
ax.legend(frameon=False, fontsize=9); ax.spines[["top", "right"]].set_visible(False); fig.tight_layout()
fig.savefig("figures/rtk_2026-10-25_floodplain_pooled_1m_2021_sets_linear.png", dpi=DPI); plt.close(fig); say("wrote figures/rtk_2026-10-25_floodplain_pooled_1m_2021_sets_linear.png")

# ---- by location bin
say(f"\nLOCATION BINS (from the 2008-offset clustering; MINE as bins): points per bin by rover:")
pts = pool.drop_duplicates(["set", "point", "loc_bin"]).groupby(["loc_bin", "set"]).size().unstack(fill_value=0); say(pts.to_string())
say(f"{'bin':30s} {'epoch':>5s} {'n':>6s} {'points':>6s} {'min':>7s} {'p05':>7s} {'p10':>7s} {'median':>7s} {'p90':>7s} {'p95':>7s} {'max':>7s} {'NMAD':>6s} {'<0.02':>6s} {'<0.05':>6s} {'>0.30':>6s}")
for bn in BORDER:
    for ep in (2021, 2008):
        q = pool[(pool.loc_bin == bn) & (pool.epoch == ep)]; v = q.dz_m
        say(f"{bn:30s} {ep:5d} {len(v):6d} {q.groupby(['set', 'point']).ngroups:6d} {v.min():+7.3f} {v.quantile(.05):+7.3f} {v.quantile(.1):+7.3f} {v.median():+7.3f} {v.quantile(.9):+7.3f} {v.quantile(.95):+7.3f} {v.max():+7.3f} {nmad(v):6.3f} {int((v < 0.02).sum()):6d} {int((v < 0.05).sum()):6d} {int((v > 0.30).sum()):6d}")
say("2021, fraction of returns below a height, by location bin:")
say("   height " + " ".join(f"{bn:>30s}" for bn in BORDER))
for h in hs:
    say(f"   {h:5.2f} m " + " ".join(f"{100 * (pool[(pool.loc_bin == bn) & (pool.epoch == 2021)].dz_m < h).mean():29.2f}%" for bn in BORDER))
for ep, png in ((2021, "figures/rtk_2026-10-25_floodplain_pooled_1m_2021_bins.png"), (2008, "figures/rtk_2026-10-25_floodplain_pooled_1m_2008_bins.png")):
    ser = [pool[(pool.loc_bin == bn) & (pool.epoch == ep)].dz_m.to_numpy() for bn in BORDER]
    density_ecdf(ser, [f"{bn} (n {len(v)})" for bn, v in zip(BORDER, ser)], [BCOL[bn] for bn in BORDER], [1.6, 1.6, 1.6], f"{ep} lidar return minus tied RTK ground (m)", png, **({"lin_bin": 0.05, "ratio": 1.3} if ep == 2008 else {}))
fig, ax = plt.subplots(figsize=(8, 4.6))
for bn in BORDER:
    v = pool[(pool.loc_bin == bn) & (pool.epoch == 2021)].dz_m; ax.hist(v, bins=b3, density=True, histtype="step", color=BCOL[bn], lw=1.6, label=f"{bn} (n {len(v)})")
ax.axvline(0, color="0.4", lw=0.8, ls=":"); ax.set_xlim(X3); ax.set_xlabel("2021 lidar return minus tied RTK ground (m)"); ax.set_ylabel("probability density (1/m)")
ax.legend(frameon=False, fontsize=9); ax.spines[["top", "right"]].set_visible(False); fig.tight_layout()
fig.savefig("figures/rtk_2026-10-25_floodplain_pooled_1m_2021_bins_linear.png", dpi=DPI); plt.close(fig); say("wrote figures/rtk_2026-10-25_floodplain_pooled_1m_2021_bins_linear.png")
fig, ax = plt.subplots(figsize=(8, 4.6)); b4 = np.arange(-0.10, 0.90, 0.02)
for bn in BORDER:
    v = pool[(pool.loc_bin == bn) & (pool.epoch == 2008)].dz_m; ax.hist(v, bins=b4, density=True, histtype="step", color=BCOL[bn], lw=1.6, label=f"{bn} (n {len(v)})")
ax.axvline(0, color="0.4", lw=0.8, ls=":"); ax.set_xlim(-0.10, 0.90); ax.set_xlabel("2008 lidar return minus tied RTK ground (m)"); ax.set_ylabel("probability density (1/m)")
ax.legend(frameon=False, fontsize=9); ax.spines[["top", "right"]].set_visible(False); fig.tight_layout()
fig.savefig("figures/rtk_2026-10-25_floodplain_pooled_1m_2008_bins_linear.png", dpi=DPI); plt.close(fig); say("wrote figures/rtk_2026-10-25_floodplain_pooled_1m_2008_bins_linear.png (2 cm bins)")
LOG.close()

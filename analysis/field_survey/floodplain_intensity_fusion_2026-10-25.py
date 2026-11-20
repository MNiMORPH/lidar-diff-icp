#!/usr/bin/env python
"""Fuse 2021 lidar intensity with the return-shape metrics and the optical stack (Andy, 2026-11-20: "Fuse with lidar
intensity. See what you get.").

Part 1, the 147 floodplain FIX points. Intensity of the 2021 SINGLE returns within 3 m below +0.50 m (the sward; single
returns because a split pulse shares its energy; the cuts are MINE), per swath (point_source_id 3040 and 3041, scan
angle ranks -5..-14 and -13..-20 deg here). A per-swath linear fit of intensity on |scan angle| over the open-field
points gives a normalized intensity = observed / fitted (a ratio, so 1 = typical sward at that angle in that swath).
Then Spearman of raw and normalized intensity with the shape metrics (floodplain_return_shape_points.csv) and the NDVI
metrics, and medians by location bin and imagery type.
Part 2, the whole floodplain at 10 m on the optical feature grid (imagery/floodplain_model_features.npz, 120 x 120
cells, floodplain mask fpm). From the 2021 clouds of tiles 4342-29-64 and 4358-29-01 (chunked, aggregated per cell so
no point table is held): count of returns, count and mean intensity of single returns, class-2 share, single-return
share. Scan angle is not available per cell without the full table, so Part 2 uses raw single-return mean intensity
and reports the swath mix per cell; the within-swath test is Part 1. Correlations (Spearman) of cell intensity with
each NDVI/NDRE/NDMI layer and the NE-like (reed-canary-like) probabilities over floodplain cells with >= 20 single
returns (MINE), split by swath where one swath dominates the cell (>= 90 % of returns).
Outputs  data/field/2026-10-25/returns/floodplain_intensity_points.csv, floodplain_intensity_cells.csv, _log.txt
         figures/rtk_2026-10-25_floodplain_intensity_points.png   normalized intensity vs April NDVI, sward median, ground share
         figures/rtk_2026-10-25_floodplain_intensity_map.png      cell intensity map beside NAIP-2021 NDVI and April-2026 NDVI
"""
import os, json
import numpy as np, pandas as pd, laspy
from scipy import stats
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
ROOT = os.path.dirname(os.path.abspath(__file__)) + "/../.."; os.chdir(ROOT)
D = "data/field/2026-10-25/returns"; LOG = open(f"{D}/floodplain_intensity_log.txt", "w")
def say(*a):
    s = " ".join(str(x) for x in a); print(s); LOG.write(s + "\n")
DPI = int(os.environ.get("FIGDPI", 170)); SUB = 0.50
S = pd.read_csv(f"{D}/floodplain_return_shape_points.csv"); P = pd.read_csv("data/field/2026-10-25/rtk_points_all.csv")
rows = []
for tag, group in (("rsplus_b", "floodplain_RSplus_B"), ("rsplus_a", "floodplain_RSplus_A"), ("cors", "floodplain_CORS_Wickert")):
    t = pd.read_csv(f"{D}/{tag}_returns_2021.csv.gz"); t = t[t.radius_3m & (t.number_of_returns == 1) & (t.dz_m < SUB)]
    for p, q in t.groupby("point"):
        for sw, qq in q.groupby("point_source_id"):
            rows.append(dict(group=group, point=int(p), swath=int(sw), n_single=len(qq), intensity_med=qq.intensity.median(), scan_abs=qq.scan_angle_rank.abs().median()))
I = pd.DataFrame(rows).merge(S[["group", "point", "loc_bin", "cls_naip_2021_k5", "q50", "iqr", "ground_share", "canopy_share", "metric_ndvi_apr", "metric_ndvi_sum", "metric_ndvi_oct", "naip_2021", "med08_3m_m"]], on=["group", "point"])
say(f"Part 1: {len(I)} point-swath rows from {I.point.nunique()} points; intensity range in the tables 471-7403 (14-bit digitizer, values well under 16383); swaths {sorted(I.swath.unique())}")
fits = {}
for sw, g in I[(I.loc_bin == "open field") & (I.n_single >= 30)].groupby("swath"):
    b = np.polyfit(g.scan_abs, g.intensity_med, 1); fits[sw] = b; say(f"   swath {sw}: open-field points {len(g)}, |scan| {g.scan_abs.min():.0f}-{g.scan_abs.max():.0f} deg, intensity median {g.intensity_med.median():.0f}, fit intensity = {b[1]:.0f} {b[0]:+.1f} x |scan| (Spearman {g[['scan_abs', 'intensity_med']].corr(method='spearman').iloc[0, 1]:+.2f})")
I["intensity_norm"] = [r.intensity_med / np.polyval(fits[r.swath], r.scan_abs) if r.swath in fits else np.nan for r in I.itertuples()]
# one row per point: weight swaths by their single-return counts
Ip = I.groupby(["group", "point"]).apply(lambda g: pd.Series(dict(intensity_med=np.average(g.intensity_med, weights=g.n_single), intensity_norm=np.average(g.intensity_norm, weights=g.n_single), n_single=g.n_single.sum(), swaths=len(g)))).reset_index()
Ip = Ip.merge(S[["group", "point", "loc_bin", "cls_naip_2021_k5", "q50", "iqr", "ground_share", "canopy_share", "metric_ndvi_apr", "metric_ndvi_sum", "metric_ndvi_oct", "naip_2021", "med08_3m_m", "E", "N"]], on=["group", "point"])
Ip.to_csv(f"{D}/floodplain_intensity_points.csv", index=False, float_format="%.4f")
pd.set_option("display.width", 250)
say("medians by location bin:"); say(Ip.groupby("loc_bin")[["intensity_med", "intensity_norm", "n_single"]].median().round(3).to_string())
say("medians by NAIP-2021 imagery type:"); say(Ip.groupby("cls_naip_2021_k5")[["intensity_med", "intensity_norm", "canopy_share"]].median().round(3).to_string())
cols = ["q50", "iqr", "ground_share", "canopy_share", "metric_ndvi_apr", "metric_ndvi_sum", "metric_ndvi_oct", "naip_2021", "med08_3m_m"]
def sp(df, a, b): d = df[[a, b]].dropna(); return d.corr(method="spearman").iloc[0, 1], len(d)
say("Spearman with raw and normalized intensity, all 147 points | open field only:")
of = Ip[Ip.loc_bin == "open field"]
for c in cols:
    r1, n1 = sp(Ip, "intensity_med", c); r2, _ = sp(Ip, "intensity_norm", c); r3, n3 = sp(of, "intensity_med", c); r4, _ = sp(of, "intensity_norm", c)
    say(f"   {c:16s} raw {r1:+.2f} norm {r2:+.2f} (n {n1}) | open field raw {r3:+.2f} norm {r4:+.2f} (n {n3})")
fig, axs = plt.subplots(1, 3, figsize=(14, 4.2)); COL = {"south zone": "#d55e00", "open field": "#0072b2", "swale stretch": "#009e73"}
for ax, c, lab in zip(axs, ("metric_ndvi_apr", "q50", "ground_share"), ("April-2026 NDVI", "2021 sward median above ground (m)", "share of returns within 5 cm of ground")):
    for b, col in COL.items():
        q = Ip[Ip.loc_bin == b]; ax.scatter(q[c], q.intensity_norm, color=col, s=16, edgecolors="k", linewidths=0.3, label=b)
    ax.set_xlabel(lab); ax.set_ylabel("normalized single-return intensity"); ax.spines[["top", "right"]].set_visible(False)
axs[0].legend(frameon=False, fontsize=8); fig.tight_layout(); fig.savefig("figures/rtk_2026-10-25_floodplain_intensity_points.png", dpi=DPI); plt.close(fig); say("wrote figures/rtk_2026-10-25_floodplain_intensity_points.png")

# ---------------------------------------------------------------- Part 2: floodplain-wide cells
z = np.load("data/field/2026-10-25/imagery/floodplain_model_features.npz", allow_pickle=True); names = list(z["names"]); stack = z["stack"]; xs = z["xs"]; ys = z["ys"]; fpm = z["fpm"].astype(bool)
res = float(xs[1] - xs[0]); x0 = xs[0] - res / 2; ytop = ys[0] + res / 2; nx, ny = len(xs), len(ys)   # ys descend: row 0 is north
acc = {k: np.zeros(ny * nx) for k in ("n", "n_single", "sum_i_single", "sum_i2_single", "n_cls2", "n_3040")}
for t in ("4342-29-64", "4358-29-01"):
    with laspy.open(f"data/after/valley/{t}_3dep_full.laz") as r:
        for ch in r.chunk_iterator(5_000_000):
            x = np.asarray(ch.x); y = np.asarray(ch.y); j = np.floor((x - x0) / res).astype(int); i = np.floor((ytop - y) / res).astype(int)
            m = (i >= 0) & (i < ny) & (j >= 0) & (j < nx)
            if not m.any(): continue
            cell = (i[m] * nx + j[m]); nr = np.asarray(ch.number_of_returns)[m]; inten = np.asarray(ch.intensity)[m].astype(float); cls = np.asarray(ch.classification)[m]; ps = np.asarray(ch.point_source_id)[m]
            if t == "4358-29-01":   # the two files overlap by 35 m; keep the overlap from the first tile only
                keep = x[m] > 580030.0; cell, nr, inten, cls, ps = cell[keep], nr[keep], inten[keep], cls[keep], ps[keep]
            acc["n"] += np.bincount(cell, minlength=ny * nx); s1 = nr == 1
            acc["n_single"] += np.bincount(cell[s1], minlength=ny * nx); acc["sum_i_single"] += np.bincount(cell[s1], weights=inten[s1], minlength=ny * nx); acc["sum_i2_single"] += np.bincount(cell[s1], weights=inten[s1] ** 2, minlength=ny * nx)
            acc["n_cls2"] += np.bincount(cell[cls == 2], minlength=ny * nx); acc["n_3040"] += np.bincount(cell[ps == 3040], minlength=ny * nx)
n = acc["n"]; ns = acc["n_single"]
with np.errstate(invalid="ignore", divide="ignore"):
    mean_i = acc["sum_i_single"] / ns; sd_i = np.sqrt(acc["sum_i2_single"] / ns - mean_i ** 2); single_share = ns / n; cls2_share = acc["n_cls2"] / n; share_3040 = acc["n_3040"] / n
cells = pd.DataFrame(dict(cell=np.arange(ny * nx), row=np.arange(ny * nx) // nx, col=np.arange(ny * nx) % nx, n=n, n_single=ns, intensity_mean_single=mean_i, intensity_sd_single=sd_i, single_share=single_share, cls2_share=cls2_share, share_swath3040=share_3040, floodplain=fpm.ravel()))
cells["E"] = xs[cells.col]; cells["N"] = ys[cells.row]
lay = [k for k in names if any(s in k for s in ("NDVI", "NDRE", "NDMI")) or k.startswith("nelike")]
for k in lay: cells[k] = stack[names.index(k)].ravel()
cells.to_csv(f"{D}/floodplain_intensity_cells.csv", index=False, float_format="%.4f")
ok = cells[cells.floodplain & (cells.n_single >= 20)]
say(f"\nPart 2: {int(fpm.sum())} floodplain cells at {res:.0f} m; {len(ok)} with >= 20 single returns; returns per cell median {int(ok.n.median())}; single share median {ok.single_share.median():.2f}; cells in the lidar-covered grid {int((cells.n > 0).sum())} of {ny * nx}")
one = ok[(ok.share_swath3040 >= 0.9) | (ok.share_swath3040 <= 0.1)]; say(f"   cells dominated (>= 90 %) by one swath: {len(one)}; swath 3040-dominated {int((one.share_swath3040 >= 0.9).sum())}, 3041-dominated {int((one.share_swath3040 <= 0.1).sum())}")
say("   Spearman of cell mean single-return intensity with the optical layers over floodplain cells (all | 3040-dominated | 3041-dominated | all cells with single_share >= 0.9, i.e. no trees):")
grass = ok[ok.single_share >= 0.9]
for k in lay:
    out = []
    for d in (ok, one[one.share_swath3040 >= 0.9], one[one.share_swath3040 <= 0.1], grass):
        r, m = sp(d, "intensity_mean_single", k); out.append(f"{r:+.2f} (n {m})")
    say(f"   {k:26s} " + " | ".join(out))
say("   intensity by NE-like class (nelike_p_s2_20211018 >= 0.5 vs < 0.5, grass cells): " + ", ".join(f"{lab} median {g.intensity_mean_single.median():.0f} (n {len(g)})" for lab, g in (("NE-like", grass[grass.nelike_p_s2_20211018 >= 0.5]), ("other", grass[grass.nelike_p_s2_20211018 < 0.5]))))
fig, axs = plt.subplots(1, 3, figsize=(16, 5.8)); ext = (xs[0] - res / 2, xs[-1] + res / 2, ys[-1] - res / 2, ys[0] + res / 2)
img = np.where(fpm, mean_i.reshape(ny, nx), np.nan); im = axs[0].imshow(img, extent=ext, cmap="gray", vmin=np.nanpercentile(img, 2), vmax=np.nanpercentile(img, 98)); axs[0].set_title("2021 single-return mean intensity, floodplain cells", fontsize=9); plt.colorbar(im, ax=axs[0], shrink=0.6)
for ax, k in zip(axs[1:], ("naip_2021-08_NDVI", "s2_2026-04-10_NDVI")):
    a = np.where(fpm, stack[names.index(k)], np.nan); im = ax.imshow(a, extent=ext, cmap="viridis", vmin=np.nanpercentile(a, 2), vmax=np.nanpercentile(a, 98)); ax.set_title(k, fontsize=9); plt.colorbar(im, ax=ax, shrink=0.6)
for ax in axs: ax.plot(P[P.group.str.startswith("floodplain") & (P.status == "FIX")].easting_tied, P[P.group.str.startswith("floodplain") & (P.status == "FIX")].northing_tied, ".", color="#d55e00", ms=2); ax.tick_params(labelsize=7)
fig.tight_layout(); fig.savefig("figures/rtk_2026-10-25_floodplain_intensity_map.png", dpi=DPI); plt.close(fig); say("wrote figures/rtk_2026-10-25_floodplain_intensity_map.png")
LOG.close()

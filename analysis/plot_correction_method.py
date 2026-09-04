#!/usr/bin/env python3
"""The correction, drawn: how the curve is measured, and what it does to a tile.

Figure A -- HOW IT IS MADE. Two real control marks' ground-return columns with the surveyed
elevation marked, so the response is visible as a position rather than a number; then all
the marks as rank vs spread with the isotonic fit through them.

Figure B -- WHAT IT DOES. The tile's own ground-return SD, the percentile that earns, and
the resulting shift in millimetres.

    env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python analysis/plot_correction_method.py \
        --tile data/derived/elba
"""
import argparse, os, sys, json
import numpy as np, pandas as pd, laspy, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, "analysis")
from control_mode_shift import STRUCT, BOX, marks
from lidar_diff_icp import groundq
from lidar_diff_icp.groundtruth.tie import _design

ap = argparse.ArgumentParser()
ap.add_argument("--tile", default="data/derived/elba")
ap.add_argument("--set", dest="set_", default="gen2_2021_control")
ap.add_argument("--curve",
                default="data/derived/ground_q_vs_class2sd_gen2_2021_control_LCP-NVA-VVA.npz")
ap.add_argument("--figdir", default="figures/method")
A = ap.parse_args()
os.makedirs(A.figdir, exist_ok=True)
CUR = groundq.load_curve(A.curve)
T = pd.read_csv(f"data/derived/control_marks_by_line_{A.set_}.csv")
M = marks(A.set_).set_index("point_id")


def column(pid):
    z = np.load(f"{STRUCT}/{A.set_}__{pid}.npz"); coef = z["surface_coef"]
    E, N, R = float(z["easting"]), float(z["northing"]), float(z["struct_radius"])
    f = laspy.read(f"{BOX}/{A.set_}__{pid}.laz")
    x, y, zz, cl = np.asarray(f.x), np.asarray(f.y), np.asarray(f.z), np.asarray(f.classification)
    g = (np.hypot(x - E, y - N) <= R) & (cl == 2)
    nn = np.sqrt(1 + coef[1]**2 + coef[2]**2)
    hg = (zz[g] - (_design(x[g]-E, y[g]-N, 2) @ coef)) / nn
    mu = (float(M.loc[pid, "elevation"]) - float(coef[0])) / nn
    return hg * 1000.0, mu * 1000.0


# ---- FIGURE A -----------------------------------------------------------------------
lo = T.iloc[(T.sd*1000 - 20).abs().argsort()].iloc[0]        # a clean column
hi = T.iloc[(T.sd*1000 - 200).abs().argsort()].iloc[0]       # a dirty one
fig, ax = plt.subplots(1, 3, figsize=(16.5, 5.0))
for a, row in zip(ax[:2], (lo, hi)):
    hg, mu = column(row.point_id)
    q = float(groundq.q_from_spread([row.sd*1000], CUR)[0])
    a.hist(hg, bins=60, color="0.75", edgecolor="none")
    a.axvline(mu, color="crimson", lw=2.4, label=f"surveyed ground (rank {row['rank']:.3f})")
    a.axvline(np.median(hg), color="tab:blue", lw=2.0, ls="--", label="class-2 median (q=0.50)")
    a.axvline(np.quantile(hg, q), color="tab:green", lw=2.0, ls="-.",
              label=f"the curve's q = {q:.3f}")
    a.set_title(f"{row.point_id}   SD = {row.sd*1000:.0f} mm   n = {int(row.n_class2)}")
    a.set_xlabel("height above the mark's fitted surface (mm)")
    a.set_ylabel("class-2 returns"); a.legend(fontsize=7.5)
sd = T.sd.to_numpy()*1000; rk = T["rank"].to_numpy()
xs = np.logspace(np.log10(sd.min()), np.log10(sd.max()), 300)
ax[2].scatter(sd, rk, s=9, c="0.72", label=f"all {len(T)} marks")
ax[2].plot(xs, groundq.q_from_spread(xs, CUR), "b-", lw=2.6, label="isotonic fit q(SD)")
ax[2].axhline(0.5, color="0.4", ls="--", lw=1.2, label="q = 0.50, the old default")
for a, row, c in ((ax[2], lo, "tab:orange"), (ax[2], hi, "tab:red")):
    a.scatter([row.sd*1000], [row["rank"]], s=130, c=c, edgecolor="k", zorder=4)
ax[2].set_xscale("log"); ax[2].set_xlabel("ground-return spread, SD (mm)")
ax[2].set_ylabel("percentile at which surveyed ground sits")
ax[2].set_title(f"the relation, fitted on {int(CUR['n_marks'])} marks")
ax[2].legend(fontsize=8, loc="lower left")
fig.suptitle("HOW THE CORRECTION IS MADE: where true ground sits in the lidar's own "
             "ground returns, against how wide that column is", y=0.99)
fig.tight_layout(rect=(0, 0, 1, 0.93))
pa = f"{A.figdir}/correction_method_A.png"
fig.savefig(pa, dpi=130, bbox_inches="tight"); plt.close(fig)
print("wrote", pa)

# ---- FIGURE B -----------------------------------------------------------------------
j = json.load(open(f"{A.tile}/corrections.json")); b = j["bounds"]; RES = float(j["res_m"])
sdg = np.load(f"{A.tile}/class2_sd_mm.npy")
qg = np.load(f"{A.tile}/gen2_q2_used.npy")
shift = (np.load(f"{A.tile}/dod_cover_q2.npy") - np.load(f"{A.tile}/dod_gen2_median.npy"))*1000
ny, nx = sdg.shape
ext = (b[0], b[0]+nx*RES, b[1], b[1]+ny*RES)
fig, ax = plt.subplots(1, 3, figsize=(17, 6.2))
for a, D, ttl, kw in (
        (ax[0], sdg, "the cell's own ground-return SD (mm)",
         dict(cmap="magma", vmin=0, vmax=200)),
        (ax[1], qg, "the percentile it earns, q(SD)",
         dict(cmap="viridis_r", vmin=0.1, vmax=0.51)),
        (ax[2], shift, "the resulting shift (mm), corrected - median",
         dict(cmap="RdBu", vmin=-120, vmax=120))):
    im = a.imshow(D, extent=ext, origin="lower", **kw)
    a.set_title(ttl); a.set_xlabel("Easting (m)")
    fig.colorbar(im, ax=a, shrink=0.62, extend="both")
ax[0].set_ylabel("Northing (m)")
fin = np.isfinite(shift)
fig.suptitle(f"WHAT IT DOES at {os.path.basename(A.tile)}: median shift "
             f"{np.median(shift[fin]):+.1f} mm, p10 {np.percentile(shift[fin],10):+.1f}, "
             f"{int((shift[fin] < -1).sum()):,} cells lowered", y=0.98)
fig.tight_layout(rect=(0, 0, 1, 0.94))
pb = f"{A.figdir}/correction_method_B_{os.path.basename(A.tile)}.png"
fig.savefig(pb, dpi=130, bbox_inches="tight"); plt.close(fig)
print("wrote", pb)

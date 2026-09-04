#!/usr/bin/env python3
"""Plot the ground-q curve and the marks it rests on -- what --diagnostics prints, drawn.

The shape panel is the one Andy asked to see at Elba: every mark shown, the binned medians
with bootstrap 95% CIs on top, and the isotonic curve through them. The bins are the ones
analysis/calibrate_ground_q.py --diagnostics already uses; nothing is re-specified here.

``--lines <tile>`` overlays the marks that sit on that tile's own gen2 flight lines, read
from the tile's gen2_psid_counts.npy (written by analysis/marks_by_flight_line.py). Those
marks are drawn on top of, not instead of, the fitted population -- a site rarely has enough
control of its own to carry a curve, and the figure should show that rather than hide it.

    ./lidar-icp/bin/python analysis/plot_ground_q_curve.py --point-types NVA \
        --lines data/derived/whitewater
"""
import argparse, os
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from lidar_diff_icp import groundq

BINS = [0, 30, 45, 60, 80, 110, 160, 250, 1e9]      # as in calibrate_ground_q --diagnostics

ap = argparse.ArgumentParser()
ap.add_argument("--set", dest="set_", default="gen2_2021_control")
ap.add_argument("--point-types", nargs="+", required=True,
                help="which control types the plotted curve is fitted on. Required for the "
                     "same reason it is required in calibrate_ground_q.py: NVA, VVA and LCP "
                     "are different populations.")
ap.add_argument("--lines", default=None, metavar="TILE",
                help="tile dir whose gen2 flight lines pick out the marks to highlight")
ap.add_argument("--out", default=None)
A = ap.parse_args()
TYPES = [t.upper() for t in A.point_types]
TAG = "-".join(sorted(TYPES))

M = pd.read_csv(f"data/derived/control_marks_by_line_{A.set_}.csv")
ctl = pd.read_csv(f"src/lidar_diff_icp/groundtruth/data/mn_se_driftless_2021_control.csv"
                  if A.set_ == "gen2_2021_control"
                  else "src/lidar_diff_icp/groundtruth/data/mn_dnr_2008_control_semn.csv")
M = M.merge(ctl[["point_id", "point_type"]].drop_duplicates("point_id"), on="point_id",
            how="left")
fit = M[M.point_type.str.upper().isin(TYPES)]
sd = fit.sd.to_numpy() * 1000.0
rk = fit["rank"].to_numpy()
iso = groundq.fit_curve(sd, rk)
print(f"{A.set_}: {len(fit)} marks of point type(s) {'+'.join(sorted(TYPES))}")

site = None
if A.lines:
    psid = np.load(os.path.join(A.lines, "gen2_psid_counts.npy"))[:, 0]
    keep = M.psid_all.apply(lambda s: bool(set(psid.tolist())
                                           & {int(v) for v in str(s).split("|")}))
    site = M[keep]
    print(f"{os.path.basename(A.lines.rstrip('/'))}: gen2 lines {sorted(psid.tolist())}; "
          f"{len(site)} control marks touch them "
          f"({site.point_type.value_counts().to_dict()})")

fig, ax = plt.subplots(1, 2, figsize=(14.5, 5.8))
xs = np.logspace(np.log10(max(sd.min(), 1.0)), np.log10(sd.max()), 300)

ax[0].scatter(sd, rk, s=12, c="0.70", label=f"{'+'.join(sorted(TYPES))} marks (n={len(fit)})")
ax[0].plot(xs, iso.predict(np.log(xs)), "b-", lw=2.6,
           label=f"isotonic fit on {'+'.join(sorted(TYPES))}")
if site is not None and len(site):
    for tp, mk in (("NVA", "o"), ("VVA", "^"), ("LCP", "s")):
        g = site[site.point_type == tp]
        if len(g):
            ax[0].scatter(g.sd*1000, g["rank"], s=95, marker=mk, c="crimson",
                          edgecolor="k", zorder=3, label=f"on the site's lines: {tp} (n={len(g)})")
ax[0].set_xscale("log"); ax[0].set_ylim(-0.04, 1.04)
ax[0].set_xlabel("class-2 spread (mm)")
ax[0].set_ylabel("rank of surveyed ground in the class-2 returns")
ax[0].set_title("every mark, and the curve fitted through them")
ax[0].legend(fontsize=8, loc="lower left")

rng = np.random.default_rng(0)
cx, cy, lo, hi, ns = [], [], [], [], []
for a, b in zip(BINS[:-1], BINS[1:]):
    s = (sd >= a) & (sd < b)
    if s.sum() < 8:
        continue
    v = rk[s]
    bs = np.array([np.median(rng.choice(v, v.size)) for _ in range(2000)])
    cx.append(np.median(sd[s])); cy.append(np.median(v)); ns.append(int(s.sum()))
    lo.append(np.percentile(bs, 2.5)); hi.append(np.percentile(bs, 97.5))
cx, cy = np.array(cx), np.array(cy)
ax[1].errorbar(cx, cy, yerr=[cy-np.array(lo), np.array(hi)-cy], fmt="ko", capsize=4, ms=6,
               label="binned median rank, bootstrap 95% CI")
for x, y, n in zip(cx, cy, ns):
    ax[1].annotate(f"n={n}", (x, y), textcoords="offset points", xytext=(0, 11),
                   ha="center", fontsize=7.5)
ax[1].plot(xs, iso.predict(np.log(xs)), "b-", lw=2.6, label="isotonic fit")
ax[1].axhline(0.50, color="0.5", ls="--", lw=1.2, label="q = 0.50, the pipeline default")
if site is not None and len(site):
    ax[1].scatter(site.sd*1000, site["rank"], s=70, c="crimson", edgecolor="k", zorder=3,
                  label=f"the site's own marks (n={len(site)})")
ax[1].set_xscale("log"); ax[1].set_ylim(-0.04, 1.04)
ax[1].set_xlabel("class-2 spread (mm)"); ax[1].set_ylabel("rank of surveyed ground")
ax[1].set_title("the shape, with its uncertainty")
ax[1].legend(fontsize=8, loc="lower left")

ttl = f"{A.set_}: ground percentile vs class-2 spread, fitted on {'+'.join(sorted(TYPES))}"
if A.lines:
    ttl += f"   |   red = marks on {os.path.basename(A.lines.rstrip('/'))}'s own flight lines"
fig.suptitle(ttl, y=0.98)
fig.tight_layout(rect=(0, 0, 1, 0.95))
out = A.out or (f"figures/sites/ground_q_curve_{TAG}"
                + (f"_{os.path.basename(A.lines.rstrip('/'))}" if A.lines else "") + ".png")
os.makedirs(os.path.dirname(out), exist_ok=True)
fig.savefig(out, dpi=130, bbox_inches="tight"); plt.close(fig)
print("wrote", out)

#!/usr/bin/env python3
"""q2(cover) fits with the floodplain masked, and the DEMs showing what the mask removes.

Top row: each tile's gen2 ground surface, hillshaded, with the floodplain mask drawn over
it -- so the population the fit rests on is visible, not just asserted.
Bottom row: the fit. Points are cover bins, sized by the cells behind them; the anchor at
q2(0) = 0.50 is imposed, not fitted, and is drawn as such.

    env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python analysis/ridgelines/plot_q2_fit.py
"""
import argparse, json, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

from lidar_diff_icp.q2cover import fit_tile, Q2_AT_ZERO
from lidar_diff_icp.viz import hillshade

ap = argparse.ArgumentParser(description=__doc__,
                             formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--tiles", nargs="+",
                default=["data/derived/elba_fulldensity", "data/derived/elbaext"])
ap.add_argument("--out", default="figures/q2_fit_and_floodplain.png")
ap.add_argument("--include-valley", action="store_true",
                help="fit WITHOUT the floodplain cut (the pre-2026-08-26-16:49 population)")
A = ap.parse_args()

n = len(A.tiles)
fig, ax = plt.subplots(2, n, figsize=(7.2 * n, 11.4),
                       gridspec_kw=dict(height_ratios=[1.35, 1.0]))
ax = np.atleast_2d(ax).reshape(2, n)
red = ListedColormap([(0, 0, 0, 0), (0.85, 0.15, 0.15, 0.45)])

for j, T in enumerate(A.tiles):
    T = T.rstrip("/"); site = os.path.basename(T)
    for fn in ("corrections_geoid.json", "corrections.json", "meta.json"):
        p = f"{T}/{fn}"
        if os.path.exists(p):
            cfg = json.load(open(p)); b = cfg["bounds"]
            res = float(cfg.get("res_m") or cfg.get("res")); break
    Z = np.load(f"{T}/z_after.npy")
    fp = f"{T}/floodplain_mask.npy"
    has_fld = os.path.exists(fp)
    # No mask means nothing is DRAWN as floodplain -- which looks identical to a tile that
    # genuinely has none. has_fld is carried into the panel title so the figure says which.
    fld = np.load(fp).astype(bool) if has_fld else np.zeros(Z.shape, bool)
    ny, nx = Z.shape
    ext = (b[0], b[0] + nx * res, b[1], b[1] + ny * res)

    a0 = ax[0, j]
    a0.imshow(hillshade(Z, res, b[0], b[1], fill_gaps=True), extent=ext, origin="lower",
              cmap="gray")
    a0.imshow(fld.astype(int), extent=ext, origin="lower", cmap=red, vmin=0, vmax=1,
              interpolation="nearest")
    if has_fld:
        sub = (f"{int(fld.sum()):,} of {fld.size:,} cells ({100*fld.mean():.1f}%) "
               f"excluded from the reference population")
    else:
        sub = ("NO floodplain_mask.npy FOR THIS TILE -- refcells skips the cut silently,\n"
               "so this tile's reference population is NOT comparable with one that has it")
        a0.text(0.5, 0.06, "no floodplain mask", transform=a0.transAxes, ha="center",
                fontsize=13, color="#b22222", fontweight="bold",
                bbox=dict(fc="white", ec="#b22222", alpha=0.85))
    a0.set_title(f"{site}: gen2 ground, floodplain mask in red\n" + sub, fontsize=10.5)
    a0.set_xlabel("Easting (m)"); a0.set_ylabel("Northing (m)")

    r = fit_tile(T, exclude_valley=not A.include_valley)
    cov = np.array([x["cover"] for x in r["bins"]])
    q2 = np.array([x["q2"] for x in r["bins"]])
    cnt = np.array([x["n"] for x in r["bins"]], float)

    a1 = ax[1, j]
    a1.scatter(cov, q2, s=18 + 90 * cnt / cnt.max(), c="#1f77b4", zorder=3,
               label=f"cover bins (area $\\propto$ cells; {len(cov)} bins)")
    xx = np.linspace(0, max(cov.max() * 1.05, 0.05), 100)
    a1.plot(xx, Q2_AT_ZERO + r["slope"] * xx, "C3", lw=2.2, zorder=4,
            label=f"fit  $q_2 = 0.50 {r['slope']:+.4f}\\,c$   "
                  f"($\\pm${r['slope_se']:.4f})")
    a1.axhline(Q2_AT_ZERO, color="0.5", ls="--", lw=1.1, zorder=1)
    a1.annotate("$q_2(0)=0.50$ imposed, not fitted", (xx[-1], Q2_AT_ZERO),
                xytext=(-6, 5), textcoords="offset points", ha="right", fontsize=8.5,
                color="0.35")
    fi = r["free_intercept"]
    a1.plot(xx, fi["a"] + fi["b"] * xx, ":", color="#7f7f7f", lw=1.6, zorder=2,
            label=f"free-intercept CHECK: a={fi['a']:.3f} "
                  f"({fi['sigma_from_imposed']:+.1f}$\\sigma$ from 0.50)")
    for u in r["unmatchable"]:
        a1.axvspan(u["lo"], u["hi"], color="0.85", zorder=0)
        a1.annotate(f"unmatchable\nn={u['n']}", (0.5 * (u["lo"] + u["hi"]), q2.min()),
                    ha="center", va="bottom", fontsize=7.5, color="0.4")
    a1.set_xlabel("canopy cover fraction (PyForestScan, >2 m, gen2)")
    a1.set_ylabel("$q_2^*$ : gen2 percentile matching gen1's median")
    fl_state = ("INCLUDED" if A.include_valley
                else ("masked" if has_fld else "NO MASK EXISTS -- cut skipped"))
    a1.set_title(f"{site}: {r['cells']['used']:,} of {r['cells']['stable']:,} stable cells, "
                 f"floodplain {fl_state}", fontsize=10)
    a1.grid(alpha=0.25); a1.legend(fontsize=8, loc="lower left")
    print(f"{site:18s} slope {r['slope']:+.4f} +/- {r['slope_se']:.4f}  "
          f"cells {r['cells']['used']:,}  bins {len(r['bins'])}  "
          f"unmatchable {len(r['unmatchable'])}  free a {fi['a']:.4f}")

fig.suptitle("q2(cover): the gen2 percentile that matches gen1's median ground\n"
             "per-site fit -- the relation depends on each pair's phenology and undergrowth",
             fontsize=12)
fig.tight_layout()
os.makedirs(os.path.dirname(A.out) or ".", exist_ok=True)
fig.savefig(A.out, dpi=140)
print("wrote", A.out)

#!/usr/bin/env python3
"""Choose the FORM of the cover-vs-offset relationship from the data, not by assumption.

Fits competing shapes to the binned medians on a declared no-change population (divides,
low curvature, low incidence) and ranks them by AIC. The candidates are chosen to be
distinguishable in kind, not just in flexibility:

    linear            no structure beyond a constant rate
    segmented         two straight lines meeting at a break, continuous (Andy's suggestion)
    segmented-smooth  the same with a softplus splice of width s, so the corner is rounded
    quadratic         smooth curvature, no preferred break
    optical depth     -ln(1-cover): Beer-Lambert if cover is 1 minus a gap fraction
    power             a + b*cover^p

Binning is by QUANTILE of the cover distribution rather than round numbers, because that
distribution is extremely skewed here (most of a divide crest is open ground) and fixed
round bins put nearly all the leverage in bins that hold almost no data. Bins are also
capped at --cover-max, beyond which the population does not constrain anything.

    env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python \
        analysis/ridgelines/cover_curve_select.py --tile data/derived/elbaext
"""
import argparse, os
import numpy as np, pandas as pd
from scipy.optimize import curve_fit

from lidar_diff_icp import binstats as bs
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

ap = argparse.ArgumentParser()
ap.add_argument("--tile", default="data/derived/elbaext")
ap.add_argument("--inc-max", type=float, default=5.0)
ap.add_argument("--curv-max", type=float, default=0.015)
ap.add_argument("--cover-max", type=float, default=None,
                help="optional upper cut. Default None: EVERY observation is binned, with the "
                     "sparse extremes carrying honest cluster-robust error rather than being "
                     "truncated -- in this data they hold the largest effect")
ap.add_argument("--block-m", type=float, default=50.0,
                help="spatial block size for the cluster-robust error (the unit of "
                     "independence: returns inside one woodlot are not independent samples)")
ap.add_argument("--nbins", type=int, default=12, help="quantile bins above the open-ground mass")
ap.add_argument("--open-cut", type=float, default=0.02, help="the near-zero cover mass gets its own bin")
A = ap.parse_args()
TILE = os.path.basename(A.tile.rstrip("/"))
TAG = "" if TILE == "elba_fulldensity" else f"_{TILE}"

df = pd.read_parquet(f"{A.tile}/beam_offset_table.parquet",
                     columns=["cell", "d_mm_corr", "incidence", "canopy_cover",
                              "curv_laplacian", "in_grid"])
df = df[df.in_grid.values]
rm = np.load(f"{A.tile}/ridge_mask.npy").astype(bool).ravel()
sel = (rm[df.cell.to_numpy()] & (df.curv_laplacian.abs().to_numpy() <= A.curv_max)
       & (df.incidence.to_numpy() < A.inc_max)
       & np.isfinite(df.d_mm_corr.to_numpy()) & np.isfinite(df.canopy_cover.to_numpy()))
d = df.d_mm_corr.to_numpy(float)[sel]; c = df.canopy_cover.to_numpy(float)[sel]
nmad = lambda a: 1.4826*np.median(np.abs(a - np.median(a)))
print("=" * 90)
print(f"COVER CURVE SELECTION  [{TILE}]  divides, |Laplacian|<={A.curv_max:g}, incidence<{A.inc_max:g}°")
print(f"{sel.sum():,} returns;  cover p50={np.percentile(c,50):.3f} p90={np.percentile(c,90):.3f} "
      f"p99={np.percentile(c,99):.3f} max={c.max():.3f}")
print(f"cover cut: {A.cover_max if A.cover_max is not None else 'none — all observations included'}")
print("=" * 90)

# one uniform treatment: quantile bins spanning EVERY observation, weighted by the
# independent information they carry (spatial blocks), not by raw return count.
import json as _json
_meta = next(_json.load(open(f"{A.tile}/{fn}")) for fn in
             ("meta.json", "corrections_geoid.json", "corrections.json")
             if os.path.exists(f"{A.tile}/{fn}"))
_res = float(_meta.get("res") or _meta.get("res_m"))
_nx = int(_meta.get("nx") or round((_meta["bounds"][2]-_meta["bounds"][0])/_res))
blk = bs.block_ids(df.cell.to_numpy()[sel], nx=_nx, res=_res, block_m=A.block_m)
if A.cover_max is not None:
    keep = c <= A.cover_max
    c, d, blk = c[keep], d[keep], blk[keep]
edges = bs.quantile_edges(c, A.nbins, first_edge=A.open_cut)
st = bs.binned_stats(c, d, edges, block=blk, min_n=200)
print(f"\n{'cover bin':>16s} {'mean cover':>11s} {'median d':>10s} {'SE_ret':>8s} "
      f"{'SE_block':>9s} {'blocks':>7s} {'n':>10s}")
for i in range(len(st)):
    print(f"{st.lo[i]:7.3f}-{st.hi[i]:<8.3f} {st.x[i]:>11.3f} {st.y[i]:>10.1f} "
          f"{st.se_return[i]:>8.2f} {st.se_block[i]:>9.2f} {st.n_block[i]:>7,} {st.n[i]:>10,}")
print(f"\n  cluster-robust errors are {np.nanmedian(st.se_block/st.se_return):.1f}x the naive ones; "
      f"all {st.n.sum():,} returns are binned, nothing truncated")
x, y, se, nn = st.x, st.y, st.se, st.n
w = st.weights

def seg(c_, a, b1, b2, c0):                      # continuous two-segment
    return a + b1*c_ + (b2-b1)*np.maximum(c_-c0, 0.0)
def segs(c_, a, b1, b2, c0, s):                  # softplus splice, width s
    s = max(abs(s), 1e-4)
    return a + b1*c_ + (b2-b1)*s*np.logaddexp(0.0, (c_-c0)/s)
def quad(c_, a, b, e): return a + b*c_ + e*c_**2
def optd(c_, a, b):    return a + b*(-np.log(1-np.clip(c_, 0, 0.98)))
def lin(c_, a, b):     return a + b*c_
def powr(c_, a, b, p): return a + b*np.power(np.clip(c_, 1e-9, None), p)

CAND = [("linear", lin, [0., -50.]),
        ("segmented", seg, [0., -20., -100., 0.15]),
        ("segmented-smooth", segs, [0., -20., -100., 0.15, 0.03]),
        ("quadratic", quad, [0., -30., -70.]),
        ("optical depth", optd, [0., -45.]),
        ("power", powr, [0., -50., 1.5])]
print(f"\n{'form':>18s} {'k':>3s} {'chi2_red':>9s} {'AIC':>8s}  parameters")
best = None
for nm, fn, p0 in CAND:
    try:
        p, _ = curve_fit(fn, x, y, p0=p0, sigma=se, absolute_sigma=True, maxfev=200000)
    except Exception as e:                                    # noqa: BLE001
        print(f"{nm:>18s}  did not converge ({type(e).__name__})"); continue
    chi2 = float(np.sum(w*(y - fn(x, *p))**2)); k = len(p); aic = chi2 + 2*k
    print(f"{nm:>18s} {k:>3d} {chi2/max(len(x)-k,1):>9.2f} {aic:>8.1f}  "
          + ", ".join(f"{v:+.4g}" for v in p))
    if best is None or aic < best[2]: best = (nm, fn, aic, p)
nm, fn, aic, p = best
print(f"\nSELECTED: {nm}   " + ", ".join(f"{v:+.5g}" for v in p))
if nm.startswith("segmented"):
    print(f"   break at cover = {p[3]:.3f};  slope {p[1]:+.1f} mm/unit below, {p[2]:+.1f} above"
          + (f";  splice width {abs(p[4]):.3f}" if len(p) > 4 else ""))
print("   predictions: " + ", ".join(f"c={cc:.2f}: {fn(np.array([cc]), *p)[0]:+.0f} mm"
                                     for cc in (0.05, 0.10, 0.20, 0.30, 0.40)))

fig, ax = plt.subplots(figsize=(8.8, 5.8), dpi=130)
ax.errorbar(x, y, yerr=se, fmt="o", ms=5, capsize=3, color="k", zorder=5,
            label=f"binned median ± robust SE (n={sel.sum():,})")
xs = np.linspace(0, float(np.max(x))*1.02, 300)
for nmi, fni, p0 in CAND:
    try: pi, _ = curve_fit(fni, x, y, p0=p0, sigma=se, absolute_sigma=True, maxfev=200000)
    except Exception: continue
    ax.plot(xs, fni(xs, *pi), "-" if nmi == nm else "--", lw=2.0 if nmi == nm else 0.9,
            alpha=1.0 if nmi == nm else .55, label=nmi + (" (selected)" if nmi == nm else ""))
ax.axhline(0, color="k", lw=.6)
ax.set_xlabel("PyForestScan canopy cover (fraction)")
ax.set_ylabel("median offset d (mm)   [gen1 − gen2; + = lower in 2021]")
ax.set_title(f"cover–offset curve on declared no-change ground — {TILE}\n"
             f"divides, |Laplacian|≤{A.curv_max:g}, incidence<{A.inc_max:g}°, block-weighted ({A.block_m:g} m)",
             fontsize=10)
ax.legend(fontsize=8); ax.grid(alpha=.3)
out = f"figures/refdatum/cover_curve_select{TAG}.png"
fig.savefig(out, bbox_inches="tight"); plt.close(fig)
print(f"\nwrote {out}")

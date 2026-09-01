#!/usr/bin/env python3
"""q2*(cover): the gen2 percentile whose elevation matches gen1's median ground.

PHYSICAL CONSTRAINT: at zero canopy cover both epochs see the true ground, so the two
medians must agree and q2(0) = 0.50 exactly. Every form here is written to satisfy that
by construction, leaving only the shape of the decline to be fitted.

    q2 = 0.5 + b*c                    linear      (1 parameter)
    q2 = 0.5 + b*c + d*c^2            quadratic   (2)
    q2 = 0.5 - b*c^k                  power       (2)
    q2 = 0.5 - b*(1 - exp(-k*c))      saturating  (2)

Inputs: gen1 = per-cell median of `d_mm_corr` (CSF ground + the four registration terms);
gen2 = the per-cell vendor class-2 near-ground column. Bins carry cluster-robust SEs from
50 m spatial blocks, converted from mm to rank units by the bin's own mm-per-rank slope.

    ./lidar-icp/bin/python analysis/ridgelines/q2_cover_fit.py
"""
import argparse, os
import numpy as np, pyarrow.parquet as pq
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.optimize import brentq, curve_fit

from lidar_diff_icp.binstats import block_ids
from lidar_diff_icp.refcells import reference_cells

_ap = argparse.ArgumentParser()
_ap.add_argument("--tile", default="data/derived/elba_fulldensity",
                 help="tile directory; the fit is PER SITE because the relation depends on "
                      "each pair's phenology, so there is no site-invariant slope to reuse")
_ap.add_argument("--weight", choices=["se", "cells"], default="se",
                 help="se = cluster-robust SE per bin; cells = weight by cell count")
_ap.add_argument("--binw", type=float, default=None,
                 help="uniform bin width in cover units (default: the quantile-ish edges)")
_ap.add_argument("--minn", type=int, default=1,
                 help="minimum cells to keep a bin; 1 = keep everything (default)")
# Per-cell sample-size requirements. Both are DEFINITIONAL at 1 and DISCRETIONARY above it:
#   --min-gen1  a per-cell median of the gen1 offsets needs >=1 return  (definitional: 1)
#   --min-gen2  a quantile of the gen2 near-ground histogram needs >=1 point (definitional: 1)
# Anything higher is a quality judgement, and on this tile it is a biased one: the cells it
# removes are canopy-enriched (median cover 0.35 against 0.19 for those kept), so raising
# either re-selects the sample toward open ground -- the exact mechanism FRAME_2026-08-26
# warns about. They were previously hardcoded at 5 and 10 with no way to see or change them,
# and at values that did not even match percentile_float_fit.py's 3 and 5 on the same cells.
_ap.add_argument("--min-gen1", type=int, default=1,
                 help="gen1 returns needed per cell (1 = definitional: a median needs a point)")
_ap.add_argument("--min-gen2", type=int, default=1,
                 help="gen2 class-2 near-ground points needed per cell (1 = definitional)")
ARGS = _ap.parse_args()

exec(open("analysis/ridgelines/percentile_float_fit.py").read().split("D = A.tile")[0]
     .replace("ap.parse_args()", "ap.parse_args([])"))

D = ARGS.tile.rstrip("/")
SITE = os.path.basename(D)
cube = np.load(f"{D}/nearground_cells_sn.npz"); cells = cube["cells"]
dz = float(cube["dz"]); zlo = float(cube["zlo"])
zf = np.load(f"{D}/z_after.npy"); N = zf.size; NX = zf.shape[1]
cover = np.load(f"{D}/canopy_cover_pfs.npy").ravel()[cells]
t = pq.read_table(f"{D}/beam_offset_table.parquet",
                  columns=["cell", "d_mm_corr", "in_grid"])
g = t["in_grid"].to_numpy().astype(bool)
ce = t["cell"].to_numpy()[g]; dc = t["d_mm_corr"].to_numpy()[g].astype(float)
vs, off, n1 = ragged_sorted(ce, dc, N)
sp = np.load(f"{D}/nearground_gen2_class_split.npz"); Hg = sp["Hg"]
Cg = np.cumsum(Hg, 1).astype(float); ng = Cg[:, -1]
stable, _ = reference_cells(D, cells=cells, slope_max=90.0)
ok = (stable & (n1[cells] >= max(1, ARGS.min_gen1)) & (ng >= max(1, ARGS.min_gen2))
      & np.isfinite(cover))
print(f"cells: {ok.sum():,} of {stable.sum():,} stable "
      f"(gen1 >= {max(1, ARGS.min_gen1)} returns, gen2 class-2 >= {max(1, ARGS.min_gen2)}; "
      f"median cover of those kept {np.nanmedian(cover[ok]):.3f}, "
      f"of those dropped {np.nanmedian(cover[stable & ~ok]):.3f})")
sel = cells[ok]; cv = cover[ok]; Cs = Cg[ok]; ns = ng[ok]
g1 = ragged_quantile(vs, off, n1, 0.50, sel); blk = block_ids(sel, NX, 5.0, 50.0)

if ARGS.binw:
    EDGES = list(np.arange(0.0, cv.max() + ARGS.binw, ARGS.binw))
else:
    EDGES = [0, 0.02, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.50, 0.65, 1.01]
X, Y, S, W, MM, LO, HI, NB = [], [], [], [], [], [], [], []
for lo, hi in zip(EDGES[:-1], EDGES[1:]):
    m = (cv > lo - 1e-9) & (cv <= hi)
    if m.sum() < ARGS.minn:
        continue
    f = lambda q: float(np.median(g1[m] - hist_quantile(Cs[m], ns[m], q, zlo, dz)))
    q = brentq(f, 1e-4, 1 - 1e-4, xtol=1e-6)
    s = (f(max(q - 0.05, 0.01)) - f(min(q + 0.05, 0.99))) / 10.0     # mm per 0.01 rank
    r = g1[m] - hist_quantile(Cs[m], ns[m], q, zlo, dz)
    ub, inv = np.unique(blk[m], return_inverse=True)
    bm = np.array([np.median(r[inv == i]) for i in range(ub.size)])
    se_mm = float(np.std(bm, ddof=1) / np.sqrt(ub.size))
    X.append(cv[m].mean()); Y.append(q); MM.append(s); W.append(int(m.sum()))
    S.append(se_mm / (100 * s))
    NB.append(int(ub.size)); LO.append(lo); HI.append(min(hi, cv[m].max()))
X, Y, S, W, MM, LO, HI, NB = map(np.array, (X, Y, S, W, MM, LO, HI, NB))
FIT = np.isfinite(Y) & np.isfinite(S) & (S > 0)      # bins usable in the fit

# Weighting. "se" uses the cluster-robust SE per bin (n enters once, through the SE).
# "cells" weights by cell count instead, so sigma ∝ 1/sqrt(n).
SIG = S if ARGS.weight == "se" else 1.0 / np.sqrt(W.astype(float))
Xf, Yf, SIGf = X[FIT], Y[FIT], SIG[FIT]

FORMS = {
    "linear      0.5+b·c":        (lambda c, b: 0.5 + b * c, [-0.25]),
    "quadratic   0.5+b·c+d·c²":   (lambda c, b, d: 0.5 + b * c + d * c**2, [-0.1, -0.2]),
    "power       0.5−b·cᵏ":       (lambda c, b, k: 0.5 - b * np.power(np.clip(c, 0, None), k),
                                   [0.34, 1.5]),
    "saturating  0.5−b(1−e^−kc)": (lambda c, b, k: 0.5 - b * (1 - np.exp(-k * c)), [0.4, 2.0]),
}
fits = {}
for name, (fn, p0) in FORMS.items():
    par, _ = curve_fit(fn, Xf, Yf, p0=p0, sigma=SIGf, absolute_sigma=(ARGS.weight=='se'), maxfev=200000)
    pred = fn(X, *par)
    chi2 = float(np.sum(((Yf - fn(Xf, *par)) / SIGf) ** 2) / (FIT.sum() - len(par)))
    fits[name] = (fn, par, pred, chi2)

print(f"q2* per cover bin, fits constrained to q2(0)=0.50, weighting = {ARGS.weight}\n")
print(f"  {'cover bin':13s} {'cells':>7s} {'blk':>5s} {'mean c':>7s} | {'q2*':>6s} {'+-':>5s} | "
      + " ".join(f"{n.split()[0][:9]:>9s}" for n in FORMS))
for i in range(len(X)):
    qs = f"{Y[i]:6.3f}" if np.isfinite(Y[i]) else "  none"
    flag = "" if FIT[i] else ("   <- no q2 in (0,1) matches" if not np.isfinite(Y[i])
                              else "   <- SE unusable")
    print(f"  {LO[i]:.2f}-{HI[i]:<8.2f} {W[i]:7,d} {NB[i]:5d} {X[i]:7.3f} | {qs} {S[i]:5.3f} | "
          + " ".join(f"{fits[n][2][i]:9.3f}" for n in FORMS) + flag)
print(f"\n  {'form':28s} {'chi2/dof':>9s} {'max|resid| mm':>14s}  parameters")
for n, (fn, par, pred, chi2) in fits.items():
    mm = np.abs((Y - pred) * 100 * MM)[FIT]
    print(f"  {n:28s} {chi2:9.2f} {mm.max():14.1f}  " + "  ".join(f"{v:+.4f}" for v in par))
print(f"\n  elevation residual left per bin (mm):")
print(f"  {'cover':>6s} " + " ".join(f"{n.split()[0][:9]:>9s}" for n in FORMS))
for i in range(len(X)):
    print(f"  {X[i]:6.3f} " + " ".join(
        (f"{(Y[i]-fits[n][2][i])*100*MM[i]:9.1f}" if np.isfinite(Y[i]) else f"{'--':>9s}")
        for n in FORMS))

xx = np.linspace(0, 0.8, 400)
fig, ax = plt.subplots(figsize=(9, 5.6), dpi=150)
ax.axhline(0.50, color="0.75", lw=0.8, ls=":", zorder=0)
ax.plot([0], [0.50], marker="*", ms=15, color="0.3", zorder=6,
        label="theory: bare ground → $q_2$ = 0.50")
ax.errorbar(X, Y, yerr=S, xerr=[X - LO, HI - X], fmt="o", ms=6, color="k", lw=1.1,
            capsize=0, zorder=5, label="measured $q_2^*$ per cover bin")
for x, y, w in zip(X, Y, W):
    ax.annotate(f"{w:,}", (x, y), textcoords="offset points", xytext=(0, 9),
                ha="center", fontsize=6.5, color="0.35")
styles = [("-", 2.0, "#1f77b4"), ("--", 1.6, "#d62728"),
          ("-.", 1.6, "#2ca02c"), (":", 2.2, "#9467bd")]
for (n, (fn, par, pred, chi2)), (ls, lw, col) in zip(fits.items(), styles):
    ax.plot(xx, fn(xx, *par), ls, lw=lw, color=col,
            label=f"{n}   ($\\chi^2$/dof {chi2:.2f})")
ax.set_xlabel("canopy cover fraction (PyForestScan, >2 m, gen2)")
ax.set_ylabel("$q_2^*$ : gen2 percentile matching gen1's median")
ax.set_title(f"gen2 percentile vs canopy cover, pinned to the median at bare ground — {SITE} "
             f"(weights: {ARGS.weight})\n"
             "labels = cells per bin; error bars = cluster-robust SE and bin span", fontsize=10)
ax.set_xlim(-0.02, 0.80); ax.set_ylim(0.10, 0.56)
ax.legend(loc="lower left", fontsize=8.5); ax.grid(alpha=0.25)
fig.tight_layout()
tag = f"{SITE}_{ARGS.weight}" + (f"_w{ARGS.binw:g}" if ARGS.binw else "")
out = f"analysis/ridgelines/q2_vs_cover_fits_{tag}.png"
fig.savefig(out)
print("\nwrote", out)

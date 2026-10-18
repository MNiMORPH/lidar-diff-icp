#!/usr/bin/env python3
"""A continuous transfer function: the percentile of gen1's returns at which the gen2
surface sits, as a function of gen1's return span -- fitted at several pooling windows,
and scored to say which window is best.

Andy, 2026-10-18: "What I really want is a continuous function representing that transfer.
Please generate that relationship at a range of different spatial window sizes and see
which one is the best."

FORM. The binned curve runs from a high plateau at low span to a low plateau past ~400 mm,
so a four-parameter logistic, each parameter with a meaning:

    pct0(span) = p_lo + (p_hi - p_lo) / (1 + exp((span - s0) / w))

    p_hi  percentile at zero vegetation      p_lo  percentile under the densest sedge
    s0    span at the midpoint of the fall   w     width of the transition (mm)

Fitted by robust least squares (soft-L1) to every flat-floodplain cell, not to bin medians.

WINDOWS. For each k in 1,3,5,7,9,15 the cell's histogram is pooled over k x k cells and
BOTH quantities are read from that pooled cloud -- pct0 = share of returns at or below the
gen2 surface, span = its p90 - p10. Self-consistent: what the function would be applied
to. (A variant with the cell's OWN span as x is reported alongside, because the choice
matters: pooling both sides smooths both and tightens any relation mechanically.)

CRITERION FOR "BEST" -- MINE, a proposal: out-of-sample skill. 5-fold cross-validation
with folds assigned by 500 m spatial BLOCK (100 x 100 cells), far wider than any window,
so no held-out cell shares pooled returns with a training cell. Reported per window:
CV R^2, residual NMAD (percentile points), and the fitted parameters with their stability
across windows. A window past which the parameters stop moving has converged.

WHAT THIS FUNCTION IS AND IS NOT, stated once. It describes where gen2's surface sits in
gen1's cloud as a function of gen1's scatter. If it were used to CHOOSE gen1's ground
quantile per cell, the floodplain DoD would be zero at every span BY CONSTRUCTION -- it
is a cross-epoch calibration (strategy B of STRATEGY_2026-09-23.md), and it would absorb
real change on the floodplain exactly as the forest-vs-open calibration did. Whether and
how to use it is a decision, not something this script takes.

    ./lidar-icp/bin/python analysis/floodplain_transfer_function.py
"""
import argparse
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pyarrow.parquet as pq
from scipy.ndimage import uniform_filter
from scipy.optimize import least_squares

ap = argparse.ArgumentParser()
ap.add_argument("--tile", default="data/derived/elbaext_delong_star_nosurf")
ap.add_argument("--windows", default="1,3,5,7,9,15")
ap.add_argument("--block-cells", type=int, default=100, help="CV block side, cells (100 = 500 m)")
ap.add_argument("--folds", type=int, default=5)
ap.add_argument("--seed", type=int, default=0)
ap.add_argument("--figdir", default="figures")
ap.add_argument("--out", default="data/derived/floodplain_transfer_function.json")
A = ap.parse_args()
S = A.tile; WINDOWS = [int(v) for v in A.windows.split(",")]
LO, HI, DZ = -1500.0, 1500.0, 20.0
edges = np.arange(LO, HI + DZ, DZ); ctr = 0.5 * (edges[:-1] + edges[1:]); NB = ctr.size

sl = np.load(f"{S}/slope.npy"); fp = np.load(f"{S}/floodplain_mask.npy").astype(bool)
sp_own = np.load(f"{S}/gen1_spread_mm.npy"); ny, nx = sl.shape
t = pq.read_table(f"{S}/beam_offset_table.parquet", columns=["cell", "d_mm_corr", "in_grid"])
g = t["in_grid"].to_numpy().astype(bool)
cell = t["cell"].to_numpy()[g]; d = t["d_mm_corr"].to_numpy()[g].astype(float)
ok = np.isfinite(d) & (d > LO) & (d < HI); cell, d = cell[ok], d[ok]
H = np.zeros((ny * nx, NB), np.int32)
np.add.at(H, (cell, np.clip(((d - LO) / DZ).astype(int), 0, NB - 1)), 1)
H = H.reshape(ny, nx, NB); n_own = H.sum(2)
flat = fp & (sl <= 2.0) & np.isfinite(sp_own) & (n_own >= 3)
iy, ix = np.nonzero(flat)
block = (iy // A.block_cells) * 10_000 + (ix // A.block_cells)
ub, binv = np.unique(block, return_inverse=True)
rng = np.random.default_rng(A.seed)
fold_of_block = rng.permutation(np.arange(ub.size) % A.folds)
fold = fold_of_block[binv]
print(f"INPUTS: {S}/beam_offset_table.parquet (d_mm_corr); floodplain_mask, slope, gen1_spread_mm")
print(f"POPULATION: flat floodplain {flat.sum():,} cells in {ub.size} blocks of {A.block_cells*5} m, {A.folds} folds")
print("FORM: pct0 = p_lo + (p_hi - p_lo) / (1 + exp((span - s0)/w)); robust least squares on all cells")
print("COLUMNS: CV R2 and NMAD = out-of-sample, pooled over folds; parameters = fit on all cells")


def model(p, x):
    return p[1] + (p[0] - p[1]) / (1.0 + np.exp(np.clip((x - p[2]) / p[3], -50, 50)))


def fit(x, y):
    r = least_squares(lambda p: model(p, x) - y, x0=[70, 20, 200, 80], loss="soft_l1", f_scale=10,
                      bounds=([0, 0, 0, 1], [100, 100, 2000, 1000]))
    return r.x


def cv(x, y):
    pred = np.empty_like(y)
    for f in range(A.folds):
        tr = fold != f
        pred[~tr] = model(fit(x[tr], y[tr]), x[~tr])
    res = y - pred
    return 1 - res.var() / y.var(), 1.4826 * np.median(np.abs(res - np.median(res)))


def quant_from_hist(F, q):
    C = np.cumsum(F, 1); tot = C[:, -1]
    return ctr[np.argmax(C >= q * tot[:, None], axis=1)]


results = []
fig, ax = plt.subplots(figsize=(8.5, 5.5))
cmap = plt.get_cmap("Blues")
xs = np.linspace(20, 800, 300)
print(f"\n  {'window':>7}{'ret/cell':>9}{'x = span of':>13}{'CV R2':>7}{'NMAD':>6}{'p_hi':>7}{'p_lo':>7}{'s0':>7}{'w':>7}")
for k_i, k in enumerate(WINDOWS):
    Hp = H.astype(np.float32) if k == 1 else np.rint(uniform_filter(H.astype(np.float32), size=(k, k, 1), mode="nearest") * k * k)
    F = Hp.reshape(-1, NB)[flat.ravel()]
    tot = F.sum(1)
    y = 100.0 * F[:, ctr <= 0].sum(1) / np.maximum(tot, 1)
    x_pool = quant_from_hist(F, 0.90) - quant_from_hist(F, 0.10)
    for xlab, x in (("pooled", x_pool), ("own cell", sp_own[flat])):
        r2, nmad = cv(x, y); p = fit(x, y)
        results.append(dict(window=k, x=xlab, returns_per_cell=float(np.median(tot)), cv_r2=float(r2),
                            cv_nmad=float(nmad), p_hi=float(p[0]), p_lo=float(p[1]), s0=float(p[2]), w=float(p[3])))
        print(f"  {f'{k}x{k}':>7}{np.median(tot):>9.0f}{xlab:>13}{r2:>7.3f}{nmad:>6.1f}{p[0]:>7.1f}{p[1]:>7.1f}{p[2]:>7.0f}{p[3]:>7.0f}", flush=True)
        if xlab == "pooled":
            ax.plot(xs, model(p, xs), color=cmap(0.35 + 0.6 * k_i / (len(WINDOWS) - 1)), lw=2,
                    label=f"{k}x{k} ({k*5} m): CV R² {r2:.2f}")
    del Hp, F
ax.axhline(50, color="#52514e", lw=0.8, ls=":")
ax.set_xlabel("gen1 return span, p90 - p10 of the pooled cloud (mm)")
ax.set_ylabel("percentile of gen1 returns at which the gen2 surface sits")
ax.set_ylim(0, 100); ax.set_xlim(0, 800)
for s_ in ("top", "right"): ax.spines[s_].set_visible(False)
ax.legend(frameon=False, fontsize=9, title="pooling window")
ax.set_title("Flat floodplain: fitted logistic transfer pct0(span) by pooling window", fontsize=10)
fig.tight_layout(); figp = f"{A.figdir}/floodplain_transfer_function.png"; fig.savefig(figp, dpi=140)
json.dump(dict(form="p_lo + (p_hi - p_lo)/(1 + exp((span - s0)/w))", population="flat floodplain, elbaext star-only",
               cv=f"{A.folds}-fold by {A.block_cells*5} m block", results=results), open(A.out, "w"), indent=1)
print(f"\nwrote {A.out}\nfigure: {figp}")

#!/usr/bin/env python3
"""The percentile of gen1's returns at which the gen2 surface sits, as a continuous
function of gen1's return span, on the flat floodplain.

Andy, 2026-10-18: "Are you able to plot the percentile of gen1 returns that correspond to
the gen2 surface as a function of vegetation cover?"

Per cell: pct0 = 100 * (share of gen1 returns with d_mm_corr <= 0), i.e. where the gen2
surface falls in the gen1 cloud. Vegetation = gen1's own p90-p10 span (the accepted proxy,
FRAME decision 0). Shown at 1x1 (the cell's own ~17 returns) and 3x3 pooled (~144), as a
binned median with the 25-75% band, over 20 equal-count bins of span. The 50th percentile
line marks where gen1's median would sit ON the gen2 surface.

Note the circularity, stated once: a cell with a wide span ABOVE the surface has both high
span and low pct0 by construction; a wide span BELOW it would have high span and HIGH pct0.
So the curve's DIRECTION is informative (which side the scatter is on), and its level at
low span is the unforced part.

    ./lidar-icp/bin/python analysis/floodplain_gen2_percentile.py
"""
import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pyarrow.parquet as pq
from scipy.ndimage import uniform_filter

ap = argparse.ArgumentParser()
ap.add_argument("--tile", default="data/derived/elbaext_delong_star_nosurf")
ap.add_argument("--figdir", default="figures")
ap.add_argument("--bins", type=int, default=20)
A = ap.parse_args()
S = A.tile
LO, HI, DZ = -1500.0, 1500.0, 20.0
edges = np.arange(LO, HI + DZ, DZ); ctr = 0.5 * (edges[:-1] + edges[1:]); NB = ctr.size

sl = np.load(f"{S}/slope.npy"); fp = np.load(f"{S}/floodplain_mask.npy").astype(bool)
sp = np.load(f"{S}/gen1_spread_mm.npy"); ny, nx = sl.shape
t = pq.read_table(f"{S}/beam_offset_table.parquet", columns=["cell", "d_mm_corr", "in_grid"])
g = t["in_grid"].to_numpy().astype(bool)
cell = t["cell"].to_numpy()[g]; d = t["d_mm_corr"].to_numpy()[g].astype(float)
ok = np.isfinite(d) & (d > LO) & (d < HI); cell, d = cell[ok], d[ok]
H = np.zeros((ny * nx, NB), np.int32)
np.add.at(H, (cell, np.clip(((d - LO) / DZ).astype(int), 0, NB - 1)), 1)
H = H.reshape(ny, nx, NB); n_own = H.sum(2)
flat = fp & (sl <= 2.0) & np.isfinite(sp) & (n_own >= 3)
below = ctr <= 0.0
print(f"INPUTS: {S}/beam_offset_table.parquet (d_mm_corr), floodplain_mask, slope, gen1_spread_mm")
print(f"POPULATION: flat floodplain (mask AND slope <= 2 deg, >= 3 returns): {flat.sum():,} cells")
print(f"COLUMNS: span = gen1 p90-p10 (mm); pct0 = percentile of gen1 returns at or below the gen2 surface; "
      f"binned median [p25, p75] over {A.bins} equal-count span bins")

x = sp[flat]
qe = np.percentile(x, np.linspace(0, 100, A.bins + 1)); qe[-1] += 1
bi = np.clip(np.searchsorted(qe, x, side="right") - 1, 0, A.bins - 1)
fig, ax = plt.subplots(figsize=(8, 5))
cols = {"1x1": "#2a78d6", "3x3": "#eb6834"}
for k in (1, 3):
    Hp = H.astype(np.float32) if k == 1 else np.rint(uniform_filter(H.astype(np.float32), size=(k, k, 1), mode="nearest") * 9)
    tot = Hp.sum(2); p0 = 100.0 * Hp[..., below].sum(2) / np.maximum(tot, 1)
    y = p0[flat]
    xm = np.array([np.median(x[bi == b]) for b in range(A.bins)])
    med = np.array([np.median(y[bi == b]) for b in range(A.bins)])
    lo = np.array([np.percentile(y[bi == b], 25) for b in range(A.bins)])
    hi = np.array([np.percentile(y[bi == b], 75) for b in range(A.bins)])
    lab = f"{k}x{k}"
    print(f"\n{lab} pooling (~{np.median(tot[flat]):.0f} returns/cell)")
    print(f"  {'span bin':>14}{'cells':>7}{'pct0 med':>10}{'p25':>7}{'p75':>7}")
    for b in range(A.bins):
        print(f"  {f'{qe[b]:.0f}-{qe[b+1]:.0f}':>14}{(bi==b).sum():>7,}{med[b]:>10.1f}{lo[b]:>7.1f}{hi[b]:>7.1f}")
    ax.fill_between(xm, lo, hi, color=cols[lab], alpha=0.15, lw=0)
    ax.plot(xm, med, color=cols[lab], lw=2, marker="o", ms=4, label=f"{lab} pooled, median and 25-75% band")
    del Hp
ax.axhline(50, color="#52514e", lw=0.8, ls=":")
ax.text(np.percentile(x, 99.5) * 0.98, 51.5, "gen1 median on the gen2 surface", ha="right", fontsize=8, color="#52514e")
ax.set_xlabel("gen1 return span in the cell, p90 - p10 (mm)  [vegetation proxy]")
ax.set_ylabel("percentile of gen1 returns at which the gen2 surface sits")
ax.set_ylim(0, 100); ax.set_xlim(0, np.percentile(x, 99.5))
for s_ in ("top", "right"): ax.spines[s_].set_visible(False)
ax.legend(frameon=False, fontsize=9, loc="upper right")
ax.set_title("Flat floodplain: where gen2's surface falls in gen1's return cloud, vs gen1 vegetation span", fontsize=10)
fig.tight_layout(); out = f"{A.figdir}/floodplain_gen2_percentile_vs_span.png"; fig.savefig(out, dpi=140)
print(f"\nfigure: {out}")

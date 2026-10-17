#!/usr/bin/env python3
"""Shape of the gen1 floodplain return PDF vs pooling width: skew, half-widths, tail form,
and the percentile at which the gen2 surface sits inside the gen1 cloud.

Andy, 2026-10-17: "How skewed is the PDF at each pooling width? Do we see differences in
its shape? Can we characterize its tail? At what percentile would the gen2 surface land
within the gen1 cloud at each size of cell?"

Same population, histograms and pooling as floodplain_return_pdf.py (flat floodplain, gen1
spread quartiles, d_mm_corr in 20 mm bins, pooled 1x1..15x15). Two levels:

  PER CELL (median over cells): moment skewness; Bowley quartile skewness
  (q75+q25-2q50)/(q75-q25), which 17 returns can support where a third moment cannot;
  the percentile of d = 0 (the gen2 surface) in the cell's returns; p10/p50/p90.
  QUARTILE-AVERAGE PDF (cells weighted equally): the same, plus the half-widths of the
  peak on each side (distance from the mode to half the peak density, 40 mm smoothing) and
  the UPPER TAIL fitted two ways over mode+60 .. mode+1000 mm -- exponential
  (log density linear in d: e-folding length L) and power law (log density linear in
  log(d - mode): exponent a) -- with each fit's R^2 in log space. A per-cell tail cannot be
  fitted at 17 returns, so the tail is a population statement.

Mine, stated: the tail window (mode+60 to mode+1000 mm; 60 mm clears the peak's own
half-width, 1000 mm is where counts run out) and the 40 mm smoothing for the half-widths.

    ./lidar-icp/bin/python analysis/floodplain_return_shape.py
"""
import argparse
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pyarrow.parquet as pq
from scipy.ndimage import gaussian_filter1d, uniform_filter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from hsm import hsm_from_hist   # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("--tile", default="data/derived/elbaext_delong_star_nosurf")
ap.add_argument("--scales", default="1,3,5,9,15")
ap.add_argument("--figdir", default="figures")
A = ap.parse_args()
S = A.tile; SCALES = [int(v) for v in A.scales.split(",")]
LO, HI, DZ = -1500.0, 1500.0, 20.0
edges = np.arange(LO, HI + DZ, DZ); ctr = 0.5 * (edges[:-1] + edges[1:]); NB = ctr.size
TAIL_LO, TAIL_HI, BW = 60.0, 1000.0, 40.0

print("INPUTS AND ROLES: as floodplain_return_pdf.py (beam table d_mm_corr; floodplain, slope, spread masks)")
print(f"PARAMETERS: bins {DZ:g} mm; scales {SCALES}; tail window mode+{TAIL_LO:g}..mode+{TAIL_HI:g} mm; half-width smoothing {BW:g} mm")
print("COLUMNS: skew = moment skewness; bowley = (q75+q25-2q50)/(q75-q25); pct0 = percentile of the gen2 surface (d=0) in the gen1 cloud;")
print("  hw- / hw+ = lower / upper half-width of the peak (mm); L = exponential e-folding length of the upper tail (mm);")
print("  a = power-law exponent of the upper tail; R2 in log-density space")

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
qs = np.percentile(sp[flat], [25, 50, 75])
qlab = [f"Q1 <{qs[0]:.0f}", f"Q2 {qs[0]:.0f}-{qs[1]:.0f}", f"Q3 {qs[1]:.0f}-{qs[2]:.0f}", f"Q4 >{qs[2]:.0f}"]
qmask = [flat & (sp >= lo) & (sp < hi) for lo, hi in zip([0] + list(qs), list(qs) + [1e9])]
print(f"\nflat floodplain cells: {flat.sum():,}; quartile edges {qs.round(1)} mm")


def cell_shape(F):
    """Vectorised per-cell statistics from histogram rows F (ncell x NB)."""
    tot = F.sum(1); keep = tot >= 3; F = F[keep]; tot = tot[keep]
    P = F / tot[:, None]
    m1 = P @ ctr; m2 = P @ ctr ** 2; m3 = P @ ctr ** 3
    var = m2 - m1 ** 2; sd = np.sqrt(np.maximum(var, 1e-9))
    skew = (m3 - 3 * m1 * m2 + 2 * m1 ** 3) / sd ** 3
    C = np.cumsum(P, 1)
    def q(p):
        return ctr[np.argmax(C >= p, axis=1)]
    q10, q25, q50, q75, q90 = (q(p) for p in (0.10, 0.25, 0.50, 0.75, 0.90))
    bowley = (q75 + q25 - 2 * q50) / np.maximum(q75 - q25, DZ)
    pct0 = 100.0 * P[:, ctr <= 0].sum(1)
    return dict(n=tot, skew=skew, bowley=bowley, pct0=pct0, q10=q10, q50=q50, q90=q90, sd=sd)


def pdf_shape(pdf):
    """Shape of one (quartile-average) density on ctr: mode, skews, pct0, half-widths, tail fits."""
    P = pdf * DZ; P = P / P.sum()
    m1 = P @ ctr; sd = np.sqrt(P @ ctr ** 2 - m1 ** 2)
    skew = (P @ (ctr - m1) ** 3) / sd ** 3
    C = np.cumsum(P); q = lambda p: ctr[np.argmax(C >= p)]
    bowley = (q(.75) + q(.25) - 2 * q(.5)) / (q(.75) - q(.25))
    pct0 = 100.0 * P[ctr <= 0].sum()
    s = gaussian_filter1d(pdf, BW / DZ); i = int(np.argmax(s)); mode = ctr[i]; half = 0.5 * s[i]
    lo_i = np.flatnonzero(s[:i] < half); hi_i = np.flatnonzero(s[i + 1:] < half)
    hwm = mode - ctr[lo_i[-1]] if lo_i.size else np.nan
    hwp = ctr[i + 1 + hi_i[0]] - mode if hi_i.size else np.nan
    # upper tail: fit log density vs (d - mode) [exponential] and vs log(d - mode) [power law]
    x = ctr - mode; w = (x >= TAIL_LO) & (x <= TAIL_HI) & (pdf > 0)
    y = np.log(pdf[w]); xe = x[w]; xp = np.log(xe)
    def fit(xx):
        c = np.polyfit(xx, y, 1); r = y - np.polyval(c, xx)
        return c[0], 1 - r.var() / y.var()
    be, r2e = fit(xe); bp, r2p = fit(xp)
    return dict(mode=mode, skew=skew, bowley=bowley, pct0=pct0, hwm=hwm, hwp=hwp,
                L=-1.0 / be, r2e=r2e, a=bp, r2p=r2p, q10=q(.10), q50=q(.50), q90=q(.90))


blues = plt.get_cmap("Blues")
fig, axes = plt.subplots(1, 4, figsize=(16, 4.2), sharey=True)
print("\nPER CELL (median over cells)")
print(f"  {'quartile':<12}{'scale':>6}{'n/cell':>7}{'skew':>7}{'bowley':>8}{'pct0':>6}{'p10':>7}{'p50':>7}{'p90':>7}{'sd':>7}")
pop_rows = []
for k_i, k in enumerate(SCALES):
    Hp = H.astype(np.float32) if k == 1 else np.rint(uniform_filter(H.astype(np.float32), size=(k, k, 1), mode="nearest") * k * k)
    flatH = Hp.reshape(-1, NB)
    for q_i, (lab, m) in enumerate(zip(qlab, qmask)):
        F = flatH[m.ravel()]
        cs = cell_shape(F)
        print(f"  {lab:<12}{f'{k}x{k}':>6}{np.median(cs['n']):>7.0f}{np.median(cs['skew']):>7.2f}{np.median(cs['bowley']):>8.2f}"
              f"{np.median(cs['pct0']):>6.0f}{np.median(cs['q10']):>7.0f}{np.median(cs['q50']):>7.0f}{np.median(cs['q90']):>7.0f}{np.median(cs['sd']):>7.0f}")
        tot = F.sum(1); okc = tot > 0
        pdf = (F[okc] / tot[okc, None]).mean(0) / DZ
        ps = pdf_shape(pdf); ps.update(q=lab, k=k); pop_rows.append(ps)
        axes[q_i].semilogy(ctr, np.where(pdf > 0, pdf, np.nan), color=blues(0.35 + 0.6 * k_i / (len(SCALES) - 1)),
                           lw=1.5, label=f"{k}x{k} ({k*5} m)")
    del Hp, flatH
print("\nQUARTILE-AVERAGE PDF (cells weighted equally)")
print(f"  {'quartile':<12}{'scale':>6}{'mode':>6}{'skew':>7}{'bowley':>8}{'pct0':>6}{'hw-':>6}{'hw+':>6}{'p10':>6}{'p50':>6}{'p90':>6}"
      f"{'L exp':>8}{'R2':>6}{'a pow':>7}{'R2':>6}")
for r in pop_rows:
    print(f"  {r['q']:<12}{f'{r['k']}x{r['k']}':>6}{r['mode']:>6.0f}{r['skew']:>7.2f}{r['bowley']:>8.2f}{r['pct0']:>6.0f}"
          f"{r['hwm']:>6.0f}{r['hwp']:>6.0f}{r['q10']:>6.0f}{r['q50']:>6.0f}{r['q90']:>6.0f}"
          f"{r['L']:>8.0f}{r['r2e']:>6.3f}{r['a']:>7.2f}{r['r2p']:>6.3f}")
for q_i, lab in enumerate(qlab):
    ax = axes[q_i]; ax.set_title(f"gen1 spread {lab} mm", fontsize=10)
    ax.axvline(0, color="#52514e", lw=0.8, ls=":"); ax.set_xlim(-600, 1200); ax.set_ylim(1e-5, 1e-2)
    ax.set_xlabel("gen1 return height above gen2 surface (mm)")
    for s_ in ("top", "right"): ax.spines[s_].set_visible(False)
axes[0].set_ylabel("density per mm (log)"); axes[0].legend(frameon=False, fontsize=8, title="pooling")
fig.suptitle("Flat floodplain: gen1 return PDF on a log axis, by pooling width -- a straight upper limb is an exponential tail", fontsize=11)
fig.tight_layout(); fig.savefig(f"{A.figdir}/floodplain_return_shape_logpdf.png", dpi=130); plt.close(fig)
print(f"\nfigure: {A.figdir}/floodplain_return_shape_logpdf.png")

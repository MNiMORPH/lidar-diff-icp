#!/usr/bin/env python3
"""Two models of what reed canary grass does to gen1's returns, tested head to head.

Andy, 2026-10-19: "One is that reed-canary grass creates a mode approx. 100 mm above the
ground. Another is that there is a fixed percentile difference between the mode and the
real ground, with actual returns on the real ground occurring. Which of these methods is
effective, or both, in recovering a bare-ground surface in the DEM?"

MODEL 1, MODE SHIFT: the sward moves the densest return cluster up by ~c; the ground is
  (mode - c). Correction = a constant below the mode (or the median), c fitted per class.
MODEL 2, FIXED PERCENTILE: ground returns exist beneath the mat; the ground sits at a fixed
  low percentile q of the cell's returns. Correction = the q-th percentile, q fitted.

POPULATION: the imagery class NE-like 2008 on open flat floodplain, NE patch EXCLUDED (the
independent members, 937 cells), and the 'both years' subset (NE-like in 2008 AND 2021,
where gen2 reads a cut sward and is the cleanest ground reference). Returns pooled 3x3
(~150) so a mode and a 5-15th percentile are estimable; residual = estimate - gen2 surface,
in mm (the gen2 surface is the ground reference; its own offset is ~0 in this class, Q).

SCORING, 5-fold CV by 500 m block: c (model 1) and q (model 2) are fitted on training
blocks to zero the class median residual, then applied to held-out cells. Reported:
held-out bias, NMAD of the residual (the per-cell error the DEM would carry), and the
slope of the residual on gen1 span (does the correction remove the vegetation dependence
or just the mean?). Baselines: the raw median, the raw mode, the minimum.

DIAGNOSTIC for the two models before any fitting: if the sward merely SHIFTS the
distribution, its lower width (median - p10) should equal the open floodplain's; if ground
returns persist beneath a mat, the lower width should be LARGER in the class.

    ./lidar-icp/bin/python analysis/rcg/mode_vs_percentile.py
"""
import sys

import numpy as np
import pandas as pd
import pyarrow.parquet as pq
import rioxarray
from scipy.ndimage import label, uniform_filter

sys.path.insert(0, "analysis"); from hsm import hsm_from_hist  # noqa: E402
S = "data/derived/elbaext_delong_star_nosurf"; R = f"{S}/rcg"; X0, Y0 = 575600.0, 4882200.0; res = 5.0
dod = np.load(f"{S}/dod.npy") * 1e3; new = np.load(f"{S}/dod_floodplain_vegcorr.npy") * 1e3; fp = np.load(f"{S}/floodplain_mask.npy").astype(bool)
sl = np.load(f"{S}/slope.npy"); lod = np.load(f"{S}/lod.npy") * 1e3; sp = np.load(f"{S}/gen1_spread_mm.npy"); ny, nx = dod.shape
flat = fp & (sl <= 2) & np.isfinite(dod); below = fp & np.isfinite(new) & (new < -lod); lab, n = label(below, structure=np.ones((3, 3))); sizes = np.bincount(lab.ravel()); NE = lab == (np.argsort(sizes[1:])[::-1][0] + 1)
yy, xx = np.mgrid[0:ny, 0:nx]; E = X0 + (xx + 0.5) * res; N = Y0 + (yy + 0.5) * res
def sample(tif):
    da = rioxarray.open_rasterio(tif).squeeze(drop=True); xs, ys = da.x.values, da.y.values; r = abs(float(xs[1] - xs[0]))
    ix = np.clip(((E - (xs[0] - r / 2)) / r).astype(int), 0, len(xs) - 1); iy = np.clip(((ys[0] + r / 2 - N) / r).astype(int), 0, len(ys) - 1); return da.values[iy, ix]
C08 = flat & (sample(f"{R}/nelike_p_ls_20081029.tif") > 0.5); C21 = flat & (sample(f"{R}/nelike_p_s2_20211018.tif") > 0.5)
z = np.load(f"{R}/gen2_cellstats.npz"); G2 = pd.DataFrame({k: z[k] for k in z.files if k != "cells"}, index=z["cells"])
open5 = np.zeros(ny * nx, bool); open5[G2.index[(G2.above2m == 0)].to_numpy()] = True; open5 = open5.reshape(ny, nx)
LO, HI, DZ = -1500.0, 1500.0, 20.0; edges = np.arange(LO, HI + DZ, DZ); ctr = 0.5 * (edges[:-1] + edges[1:]); NB = ctr.size
t = pq.read_table(f"{S}/beam_offset_table.parquet", columns=["cell", "d_mm_corr", "in_grid"]); g = t["in_grid"].to_numpy().astype(bool)
cell = t["cell"].to_numpy()[g]; d = t["d_mm_corr"].to_numpy()[g].astype(float); ok = np.isfinite(d) & (d > LO) & (d < HI)
H = np.zeros((ny * nx, NB), np.int32); np.add.at(H, (cell[ok], np.clip(((d[ok] - LO) / DZ).astype(int), 0, NB - 1)), 1); H = H.reshape(ny, nx, NB)
Hp = np.rint(uniform_filter(H.astype(np.float32), size=(3, 3, 1), mode="nearest") * 9)
def qtl(F, q):
    c = np.cumsum(F, 1); tot = c[:, -1]; return ctr[(c >= q * tot[:, None]).argmax(1)]
nmad = lambda v: 1.4826 * np.nanmedian(np.abs(v - np.nanmedian(v)))

pops = {"NE-like 2008, open, no NE patch": C08 & open5 & ~NE, "both years (2008 & 2021), open, no NE patch": C08 & C21 & open5 & ~NE,
        "open floodplain outside both classes": flat & open5 & ~C08 & ~C21, "NE patch itself (closed canopy)": NE & flat}
print("DIAGNOSTIC: lower and upper width of the pooled-3x3 return distribution (mm); a rigid shift keeps the widths, a mat over ground widens the lower one")
print(f"  {'population':<46}{'cells':>7}{'median':>8}{'mode':>7}{'p10':>7}{'med-p10':>9}{'p90-med':>9}{'share<=0':>9}")
F_all = {}
for name, m in pops.items():
    idx = np.flatnonzero(m.ravel()); F = Hp.reshape(-1, NB)[idx]; keep = F.sum(1) >= 30; idx, F = idx[keep], F[keep]; F_all[name] = (idx, F)
    med = qtl(F, 0.5); p10 = qtl(F, 0.1); p90 = qtl(F, 0.9); mode = np.array([hsm_from_hist(F[i], ctr) for i in range(len(F))]); le0 = F[:, ctr <= 0].sum(1) / F.sum(1)
    print(f"  {name:<46}{idx.size:>7,}{np.median(med):>8.0f}{np.median(mode):>7.0f}{np.median(p10):>7.0f}{np.median(med - p10):>9.0f}{np.median(p90 - med):>9.0f}{np.median(le0):>9.2f}")

print("\nHEAD TO HEAD, 5-fold CV by 500 m block; residual = estimate - gen2 surface (mm). bias = held-out median; NMAD = per-cell error; slope = residual vs gen1 span (mm/mm)")
QS = [0.05, 0.10, 0.15, 0.20, 0.25]
for name in ("NE-like 2008, open, no NE patch", "both years (2008 & 2021), open, no NE patch", "NE patch itself (closed canopy)"):
    idx, F = F_all[name]; blk = ((idx // nx) // 100) * 1000 + ((idx % nx) // 100); ub, binv = np.unique(blk, return_inverse=True); fold = np.random.default_rng(0).permutation(np.arange(ub.size) % 5)[binv]
    med = qtl(F, 0.5); mode = np.array([hsm_from_hist(F[i], ctr) for i in range(len(F))]); mn = qtl(F, 0.0001); span = sp.ravel()[idx]
    quant = {q: qtl(F, q) for q in QS}
    print(f"\n  {name} ({idx.size:,} cells, {ub.size} blocks)")
    print(f"  {'estimator':<40}{'bias':>7}{'NMAD':>7}{'slope':>8}")
    def score(est_fn):
        out = np.full(idx.size, np.nan)
        for k in range(5):
            tr = fold != k; out[~tr] = est_fn(tr)[~tr]
        okm = np.isfinite(out) & np.isfinite(span); s_ = np.polyfit(span[okm], out[okm], 1)[0] if okm.sum() > 10 else np.nan
        return np.nanmedian(out), nmad(out), s_
    rows = [("raw median (no correction)", lambda tr: med), ("raw mode", lambda tr: mode), ("minimum", lambda tr: mn),
            ("MODEL 1: median - c (c fitted)", lambda tr: med - np.median(med[tr])), ("MODEL 1: mode - c (c fitted)", lambda tr: mode - np.median(mode[tr]))]
    for q in QS: rows.append((f"MODEL 2: fixed p{int(q*100)} (no fit)", (lambda q_: (lambda tr: quant[q_]))(q)))
    def best_q(tr):
        # q chosen on training blocks to zero the median residual, from a fine ladder
        best, bq = 1e9, 0.5
        for q in np.arange(0.02, 0.6, 0.01):
            v = abs(np.median(qtl(F[tr], q)))
            if v < best: best, bq = v, q
        return bq
    def m2(tr):
        q = best_q(tr); return qtl(F, q)
    rows.append(("MODEL 2: percentile q fitted per fold", m2))
    for lab_, fn in rows:
        b, nm_, s_ = score(fn); print(f"  {lab_:<40}{b:>7.0f}{nm_:>7.0f}{s_:>8.3f}")
    qf = [best_q(fold != k) for k in range(5)]; print(f"  (fitted q per fold: {[f'{q:.2f}' for q in qf]}; fitted c per fold, median: {[f'{np.median(med[fold != k]):.0f}' for k in range(5)]})")

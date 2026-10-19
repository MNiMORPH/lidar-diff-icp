#!/usr/bin/env python3
"""Is there a characteristic SHAPE in the return-height distribution that marks reed canary
grass presence, or its cover intensity -- detectable from the lidar alone?

Andy, 2026-10-19: use the imagery maps of presence/intensity as the independent truth,
then look for coherent structure in the lidar returns, in 2008 and/or 2021.

SHAPE FEATURES, translation-invariant (measured about the cell's OWN mode or median, so
no cross-epoch reference enters; the detector would run on one survey by itself):
  gen1 (3x3-pooled CSF ground returns, ~150/cell, 20 mm bins):
    lo_w = median - p10;  hi_w = p90 - median;  mode_med = mode - median;  bowley skew;
    peak = share within +/-50 mm of the mode;  tail_L = e-folding length of the upper tail
    (mean excess above mode+100 mm);  tail_mass = share above mode+100;  low_mass = share
    below mode-100;  n = returns;  intensity median (lidar-only, not shape).
  gen2 (full-density 2021 cloud streamed once into per-cell histograms of height above
    gen2's own surface, +/-1.5 m, 20 mm bins, flat floodplain cells; cached):
    the same shape set, plus vendor class-2 share and total returns.
TRUTH: the imagery class of the survey's own year (presence, p > 0.5) and the continuous
class score logit(p) (intensity proxy). Domain: open flat floodplain, NE patch EXCLUDED
(pass B) -- and, for reference, with the patch included.
SCORING: presence = LDA AUC, 5-fold CV by 500 m block, (a) shape-only features, (b) shape +
level features (the earlier fingerprint, for comparison); intensity = OLS of logit(p) on
shape features, block-CV R^2. Plus the characteristic shape itself: mode-centred mean PDFs
by quartile of imagery class score.

    ./lidar-icp/bin/python analysis/rcg/return_shape.py
"""
import os
import sys

import laspy
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
import rioxarray
from scipy.ndimage import label, uniform_filter
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.metrics import roc_auc_score

sys.path.insert(0, "analysis"); from hsm import hsm_from_hist  # noqa: E402
S = "data/derived/elbaext_delong_star_nosurf"; R = f"{S}/rcg"; FIG = "figures"; X0, Y0 = 575600.0, 4882200.0; res = 5.0
dod = np.load(f"{S}/dod.npy") * 1e3; new = np.load(f"{S}/dod_floodplain_vegcorr.npy") * 1e3; fp = np.load(f"{S}/floodplain_mask.npy").astype(bool)
sl = np.load(f"{S}/slope.npy"); lod = np.load(f"{S}/lod.npy") * 1e3; z2 = np.load(f"{S}/z_after.npy"); ny, nx = dod.shape
flat = fp & (sl <= 2) & np.isfinite(dod); below = fp & np.isfinite(new) & (new < -lod); lab, n = label(below, structure=np.ones((3, 3))); sizes = np.bincount(lab.ravel()); NE = lab == (np.argsort(sizes[1:])[::-1][0] + 1)
yy, xx = np.mgrid[0:ny, 0:nx]; E = X0 + (xx + 0.5) * res; N = Y0 + (yy + 0.5) * res
def sample(tif):
    da = rioxarray.open_rasterio(tif).squeeze(drop=True); xs, ys = da.x.values, da.y.values; r = abs(float(xs[1] - xs[0]))
    ix = np.clip(((E - (xs[0] - r / 2)) / r).astype(int), 0, len(xs) - 1); iy = np.clip(((ys[0] + r / 2 - N) / r).astype(int), 0, len(ys) - 1); return da.values[iy, ix]
p08 = sample(f"{R}/nelike_p_ls_20081029.tif"); p21 = sample(f"{R}/nelike_p_s2_20211018.tif")
logit = lambda p: np.log(np.clip(p, 1e-4, 1 - 1e-4) / (1 - np.clip(p, 1e-4, 1 - 1e-4)))
zc = np.load(f"{R}/gen2_cellstats.npz"); G2s = pd.DataFrame({k: zc[k] for k in zc.files if k != "cells"}, index=zc["cells"])
open5 = np.zeros(ny * nx, bool); open5[G2s.index[(G2s.above2m == 0)].to_numpy()] = True; open5 = open5.reshape(ny, nx)
LO, HI, DZ = -1500.0, 1500.0, 20.0; edges = np.arange(LO, HI + DZ, DZ); ctr = 0.5 * (edges[:-1] + edges[1:]); NB = ctr.size

# gen1 pooled histograms (+ intensity)
t = pq.read_table(f"{S}/beam_offset_table.parquet", columns=["cell", "d_mm_corr", "in_grid", "intensity"]); g = t["in_grid"].to_numpy().astype(bool)
cell = t["cell"].to_numpy()[g]; d = t["d_mm_corr"].to_numpy()[g].astype(float); it = t["intensity"].to_numpy()[g].astype(float); ok = np.isfinite(d) & (d > LO) & (d < HI)
H1 = np.zeros((ny * nx, NB), np.int32); np.add.at(H1, (cell[ok], np.clip(((d[ok] - LO) / DZ).astype(int), 0, NB - 1)), 1); H1 = H1.reshape(ny, nx, NB)
H1p = np.rint(uniform_filter(H1.astype(np.float32), size=(3, 3, 1), mode="nearest") * 9).reshape(-1, NB)
int1 = pd.Series(it[ok]).groupby(cell[ok]).median()
# gen2 per-cell histograms, streamed once, flat cells only (NOT pooled: ~200 returns per cell already)
cache = f"{R}/gen2_cellhist.npz"
if not os.path.exists(cache):
    flat_idx = np.flatnonzero(flat.ravel()); pos = np.full(ny * nx, -1); pos[flat_idx] = np.arange(flat_idx.size); H2 = np.zeros((flat_idx.size, NB), np.int32)
    with laspy.open("data/after/3dep2021_fulldensity.laz") as f2:
        print(f"streaming gen2: {f2.header.point_count:,} points", flush=True)
        for ch in f2.chunk_iterator(8_000_000):
            cx = np.asarray(ch.x); cy = np.asarray(ch.y); cz = np.asarray(ch.z); ccl = np.asarray(ch.classification)
            jx = ((cx - X0) / res).astype(int); jy = ((cy - Y0) / res).astype(int); okc = (jx >= 0) & (jx < nx) & (jy >= 0) & (jy < ny) & (ccl != 7) & (ccl != 18)
            cc = (jy * nx + jx)[okc]; pp = pos[cc]; keep = pp >= 0; hh = (cz[okc][keep] - z2.ravel()[cc[keep]]) * 1e3; inr = (hh > LO) & (hh < HI)
            np.add.at(H2, (pp[keep][inr], np.clip(((hh[inr] - LO) / DZ).astype(int), 0, NB - 1)), 1)
    np.savez_compressed(cache, cells=flat_idx, H=H2)
zz = np.load(cache); H2full = np.zeros((ny * nx, NB), np.int32); H2full[zz["cells"]] = zz["H"]

def shape_feats(F):
    tot = F.sum(1); P = F / np.maximum(tot, 1)[:, None]; C = np.cumsum(P, 1)
    q = lambda p: ctr[np.argmax(C >= p, axis=1)]
    p10, p25, p50, p75, p90 = (q(p) for p in (0.1, 0.25, 0.5, 0.75, 0.9))
    mode = np.array([hsm_from_hist(F[i], ctr) if tot[i] > 0 else np.nan for i in range(len(F))])
    rel = ctr[None, :] - mode[:, None]
    peak = (P * ((rel >= -50) & (rel <= 50))).sum(1); tail_mass = (P * (rel > 100)).sum(1); low_mass = (P * (rel < -100)).sum(1)
    w = P * (rel > 100); tail_L = (w * (rel - 100)).sum(1) / np.maximum(w.sum(1), 1e-9)
    return pd.DataFrame(dict(lo_w=p50 - p10, hi_w=p90 - p50, mode_med=mode - p50, bowley=(p75 + p25 - 2 * p50) / np.maximum(p75 - p25, DZ), peak=peak, tail_L=tail_L, tail_mass=tail_mass, low_mass=low_mass, n=tot)), mode, p50

nmad = lambda v: 1.4826 * np.nanmedian(np.abs(v - np.nanmedian(v)))
def dprime(a, b): return abs(np.nanmedian(a) - np.nanmedian(b)) / np.sqrt((nmad(a) ** 2 + nmad(b) ** 2) / 2 + 1e-12)
def blocks(idx): blk = ((idx // nx) // 100) * 1000 + ((idx % nx) // 100); ub, binv = np.unique(blk, return_inverse=True); return np.random.default_rng(0).permutation(np.arange(ub.size) % 5)[binv]
def cv_auc(X, y, fold):
    pred = np.zeros(y.size)
    for k in range(5):
        tr = fold != k; m = LinearDiscriminantAnalysis().fit(X[tr], y[tr]); pred[~tr] = m.decision_function(X[~tr])
    return roc_auc_score(y, pred), LinearDiscriminantAnalysis().fit(X, y).coef_[0] * X.std(0)
def cv_r2(X, y, fold):
    Xm = np.column_stack([np.ones(len(X)), X]); pred = np.zeros(y.size)
    for k in range(5):
        tr = fold != k; b, *_ = np.linalg.lstsq(Xm[tr], y[tr], rcond=None); pred[~tr] = Xm[~tr] @ b
    return 1 - (y - pred).var() / y.var()

def run(name, Hall, cls, score, extra, domain, pdf_ax):
    idx = np.flatnonzero(domain.ravel()); F = Hall[idx]; keep = F.sum(1) >= 30; idx, F = idx[keep], F[keep]
    feats, mode, p50 = shape_feats(F); 
    for k, v in extra.items(): feats[k] = v.reindex(idx).to_numpy() if isinstance(v, pd.Series) else v.ravel()[idx]
    y = cls.ravel()[idx].astype(int); s_ = score.ravel()[idx]; fin = np.all(np.isfinite(feats.to_numpy(float)), 1) & np.isfinite(s_)
    feats, y, s_, idx, F, mode = feats[fin], y[fin], s_[fin], idx[fin], F[fin], mode[fin]; fold = blocks(idx)
    print(f"\n{'=' * 100}\n{name}: {idx.size:,} cells, {y.sum():,} in the imagery class")
    print(f"  {'shape feature':<10}{'in class':>18}{'out':>18}{'d':>7}{'corr with class score':>24}")
    rows = []
    for c in feats.columns:
        v = feats[c].to_numpy(float); dd = dprime(v[y == 1], v[y == 0]); r = np.corrcoef(v, s_)[0, 1]
        print(f"  {c:<10}{np.nanmedian(v[y==1]):>10.2f} ({nmad(v[y==1]):5.2f}){np.nanmedian(v[y==0]):>10.2f} ({nmad(v[y==0]):5.2f}){dd:>7.2f}{r:>24.3f}"); rows.append(dict(survey=name, feature=c, in_med=np.nanmedian(v[y == 1]), out_med=np.nanmedian(v[y == 0]), d=dd, r_score=r))
    shape_cols = [c for c in feats.columns if c not in ("n", "intensity", "c2_share")]
    Xs = feats[shape_cols].to_numpy(float); Xa = feats.to_numpy(float)
    a1, w1 = cv_auc(Xs, y, fold); a2, w2 = cv_auc(Xa, y, fold)
    print(f"  PRESENCE, lidar-only LDA, block CV:  shape only AUC {a1:.3f};  shape + n/intensity/class2 AUC {a2:.3f}")
    print("   weights (shape only): " + ", ".join(f"{c} {w:+.2f}" for c, w in sorted(zip(shape_cols, w1), key=lambda t: -abs(t[1]))))
    print(f"  INTENSITY, OLS of imagery class score logit(p) on features, block-CV R2:  shape only {cv_r2(Xs, s_, fold):.3f};  all {cv_r2(Xa, s_, fold):.3f}")
    # characteristic shape: mode-centred mean PDF by quartile of class score (cells weighted equally)
    # groups by imagery class probability (logit(p) is too skewed for quartiles: most cells sit at p ~ 0)
    pp = 1 / (1 + np.exp(-s_)); groups = [("outside class, p < 0.5", pp < 0.5, "#9e9e9e"), ("in class, 0.5 <= p < 0.95", (pp >= 0.5) & (pp < 0.95), "#eb6834"), ("in class, p >= 0.95", pp >= 0.95, "#8b0000")]
    for lab_, m, col in groups:
        if m.sum() == 0: continue
        grid = np.arange(-400, 801, DZ); acc = np.zeros(grid.size)
        for i in np.flatnonzero(m):
            P = F[i] / F[i].sum(); acc += np.interp(grid, ctr - mode[i], P, left=0, right=0)
        pdf_ax.plot(grid, acc / m.sum() / DZ, color=col, lw=1.8, label=f"{lab_} (n={m.sum():,})")
    pdf_ax.axvline(0, color="#52514e", lw=0.8, ls=":"); pdf_ax.set_xlabel("height relative to the cell's own mode (mm)"); pdf_ax.set_ylabel("density per mm"); pdf_ax.set_title(name, fontsize=10); pdf_ax.legend(fontsize=8, frameon=False)
    for sp_ in ("top", "right"): pdf_ax.spines[sp_].set_visible(False)
    return pd.DataFrame(rows)

fig, axes = plt.subplots(2, 2, figsize=(16, 10))
dom_b = flat & open5 & ~NE; dom_a = flat & open5
out = []
out.append(run("gen1 (Nov 2008), open floodplain, NE patch EXCLUDED vs 2008 imagery class", H1p, (p08 > 0.5), logit(p08), {"intensity": int1}, dom_b, axes[0, 0]))
out.append(run("gen2 (May 2021), open floodplain, NE patch EXCLUDED vs 2021 imagery class", H2full, (p21 > 0.5), logit(p21), {"c2_share": G2s.c2_share}, dom_b, axes[1, 0]))
out.append(run("gen1, open floodplain, NE patch INCLUDED", H1p, (p08 > 0.5), logit(p08), {"intensity": int1}, dom_a, axes[0, 1]))
out.append(run("gen2, open floodplain, NE patch INCLUDED", H2full, (p21 > 0.5), logit(p21), {"c2_share": G2s.c2_share}, dom_a, axes[1, 1]))
pd.concat(out).to_csv(f"{R}/return_shape.csv", index=False)
fig.suptitle("Mode-centred return-height distributions, cells weighted equally, grouped by the imagery class probability", fontsize=11)
fig.tight_layout(); fig.savefig(f"{FIG}/rcg_return_shape.png", dpi=110); print(f"\nfigure: {FIG}/rcg_return_shape.png; table: {R}/return_shape.csv")

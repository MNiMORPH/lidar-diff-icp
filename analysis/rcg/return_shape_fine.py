#!/usr/bin/env python3
"""Two refinements of return_shape.py.

(1) gen2 at 5 mm bins. At 20 mm bins every gen2 shape feature was identical inside and
outside the class (lower and upper widths 40 mm, peak share 0.86): the May-2021 surface is
too tight for 20 mm bins (instrument sigma ~11 mm), while the vendor's class-2 share
differed by d' 4.1. So gen2's per-cell histogram is rebuilt at 5 mm over +/-600 mm (one
more streaming pass, cached), and the shape windows scale with it: peak within +/-20 mm of
the mode, tails beyond +/-50 mm.

(2) Cover INTENSITY from a continuous spectral quantity instead of logit(p), which
saturates: the class's own strongest band, NIR reflectance relative to the flat
floodplain's median in the same scene (Landsat-5 2008-10-29 for gen1; Sentinel-2
2021-10-18 B08 for gen2), sampled to 5 m. Reported: block-CV R2 of that index on the
shape features, and the shape features by quartile of the index.

    ./lidar-icp/bin/python analysis/rcg/return_shape_fine.py
"""
import os
import sys

import laspy
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
import rioxarray
from scipy.ndimage import label, uniform_filter
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.metrics import roc_auc_score

sys.path.insert(0, "analysis"); from hsm import hsm_from_hist  # noqa: E402
S = "data/derived/elbaext_delong_star_nosurf"; R = f"{S}/rcg"; X0, Y0 = 575600.0, 4882200.0; res = 5.0
dod = np.load(f"{S}/dod.npy") * 1e3; new = np.load(f"{S}/dod_floodplain_vegcorr.npy") * 1e3; fp = np.load(f"{S}/floodplain_mask.npy").astype(bool)
sl = np.load(f"{S}/slope.npy"); lod = np.load(f"{S}/lod.npy") * 1e3; z2 = np.load(f"{S}/z_after.npy"); ny, nx = dod.shape
flat = fp & (sl <= 2) & np.isfinite(dod); below = fp & np.isfinite(new) & (new < -lod); lab, n = label(below, structure=np.ones((3, 3))); sizes = np.bincount(lab.ravel()); NE = lab == (np.argsort(sizes[1:])[::-1][0] + 1)
yy, xx = np.mgrid[0:ny, 0:nx]; E = X0 + (xx + 0.5) * res; N = Y0 + (yy + 0.5) * res
def sample(tif, band=None):
    da = rioxarray.open_rasterio(tif); da = da.squeeze(drop=True) if band is None else da.isel(band=band); xs, ys = da.x.values, da.y.values; r = abs(float(xs[1] - xs[0]))
    ix = np.clip(((E - (xs[0] - r / 2)) / r).astype(int), 0, len(xs) - 1); iy = np.clip(((ys[0] + r / 2 - N) / r).astype(int), 0, len(ys) - 1); return da.values[iy, ix]
p08 = sample(f"{R}/nelike_p_ls_20081029.tif"); p21 = sample(f"{R}/nelike_p_s2_20211018.tif")
nir08 = sample(f"{R}/ls_20081029_bands.tif", band=3); nir21 = sample(f"{R}/s2_20211018_bands.tif", band=6)      # nir08 is band index 3 of 6; B08 index 6 of 10
zc = np.load(f"{R}/gen2_cellstats.npz"); G2s = pd.DataFrame({k: zc[k] for k in zc.files if k != "cells"}, index=zc["cells"])
open5 = np.zeros(ny * nx, bool); open5[G2s.index[(G2s.above2m == 0)].to_numpy()] = True; open5 = open5.reshape(ny, nx)
dom = flat & open5 & ~NE
relnir08 = nir08 / np.nanmedian(nir08[flat & open5]); relnir21 = nir21 / np.nanmedian(nir21[flat & open5])

# gen1 (20 mm bins, pooled) as before
LO, HI, DZ = -1500.0, 1500.0, 20.0; edges = np.arange(LO, HI + DZ, DZ); ctr = 0.5 * (edges[:-1] + edges[1:]); NB = ctr.size
t = pq.read_table(f"{S}/beam_offset_table.parquet", columns=["cell", "d_mm_corr", "in_grid"]); g = t["in_grid"].to_numpy().astype(bool)
cell = t["cell"].to_numpy()[g]; d = t["d_mm_corr"].to_numpy()[g].astype(float); ok = np.isfinite(d) & (d > LO) & (d < HI)
H1 = np.zeros((ny * nx, NB), np.int32); np.add.at(H1, (cell[ok], np.clip(((d[ok] - LO) / DZ).astype(int), 0, NB - 1)), 1); H1 = H1.reshape(ny, nx, NB)
H1p = np.rint(uniform_filter(H1.astype(np.float32), size=(3, 3, 1), mode="nearest") * 9).reshape(-1, NB)
# gen2 at 5 mm bins
LO2, HI2, DZ2 = -600.0, 600.0, 5.0; edges2 = np.arange(LO2, HI2 + DZ2, DZ2); ctr2 = 0.5 * (edges2[:-1] + edges2[1:]); NB2 = ctr2.size
cache = f"{R}/gen2_cellhist5.npz"
if not os.path.exists(cache):
    flat_idx = np.flatnonzero(flat.ravel()); pos = np.full(ny * nx, -1); pos[flat_idx] = np.arange(flat_idx.size); H2 = np.zeros((flat_idx.size, NB2), np.int32)
    with laspy.open("data/after/3dep2021_fulldensity.laz") as f2:
        print(f"streaming gen2 at 5 mm bins: {f2.header.point_count:,} points", flush=True)
        for ch in f2.chunk_iterator(8_000_000):
            cx = np.asarray(ch.x); cy = np.asarray(ch.y); cz = np.asarray(ch.z); ccl = np.asarray(ch.classification)
            jx = ((cx - X0) / res).astype(int); jy = ((cy - Y0) / res).astype(int); okc = (jx >= 0) & (jx < nx) & (jy >= 0) & (jy < ny) & (ccl != 7) & (ccl != 18)
            cc = (jy * nx + jx)[okc]; pp = pos[cc]; keep = pp >= 0; hh = (cz[okc][keep] - z2.ravel()[cc[keep]]) * 1e3; inr = (hh > LO2) & (hh < HI2)
            np.add.at(H2, (pp[keep][inr], np.clip(((hh[inr] - LO2) / DZ2).astype(int), 0, NB2 - 1)), 1)
    np.savez_compressed(cache, cells=flat_idx, H=H2)
zz = np.load(cache); H2full = np.zeros((ny * nx, NB2), np.int32); H2full[zz["cells"]] = zz["H"]

def shape_feats(F, c, dz, pk, tl):
    tot = F.sum(1); P = F / np.maximum(tot, 1)[:, None]; C = np.cumsum(P, 1); q = lambda p: c[np.argmax(C >= p, axis=1)]
    p10, p25, p50, p75, p90 = (q(p) for p in (0.1, 0.25, 0.5, 0.75, 0.9)); mode = np.array([hsm_from_hist(F[i], c) if tot[i] > 0 else np.nan for i in range(len(F))]); rel = c[None, :] - mode[:, None]
    peak = (P * (np.abs(rel) <= pk)).sum(1); tail_mass = (P * (rel > tl)).sum(1); low_mass = (P * (rel < -tl)).sum(1); w = P * (rel > tl); tail_L = (w * (rel - tl)).sum(1) / np.maximum(w.sum(1), 1e-9)
    return pd.DataFrame(dict(lo_w=p50 - p10, hi_w=p90 - p50, mode_med=mode - p50, bowley=(p75 + p25 - 2 * p50) / np.maximum(p75 - p25, dz), peak=peak, tail_L=tail_L, tail_mass=tail_mass, low_mass=low_mass))
nmad = lambda v: 1.4826 * np.nanmedian(np.abs(v - np.nanmedian(v)))
def dprime(a, b): return abs(np.nanmedian(a) - np.nanmedian(b)) / np.sqrt((nmad(a) ** 2 + nmad(b) ** 2) / 2 + 1e-12)
def blocks(idx): blk = ((idx // nx) // 100) * 1000 + ((idx % nx) // 100); ub, binv = np.unique(blk, return_inverse=True); return np.random.default_rng(0).permutation(np.arange(ub.size) % 5)[binv]
def cv_auc(X, y, fold):
    pred = np.zeros(y.size)
    for k in range(5):
        tr = fold != k; m = LinearDiscriminantAnalysis().fit(X[tr], y[tr]); pred[~tr] = m.decision_function(X[~tr])
    return roc_auc_score(y, pred)
def cv_r2(X, y, fold):
    Xm = np.column_stack([np.ones(len(X)), X]); pred = np.zeros(y.size)
    for k in range(5):
        tr = fold != k; b, *_ = np.linalg.lstsq(Xm[tr], y[tr], rcond=None); pred[~tr] = Xm[~tr] @ b
    return 1 - (y - pred).var() / y.var()
def run(name, Hall, c, dz, pk, tl, cls, inten, extra):
    idx = np.flatnonzero(dom.ravel()); F = Hall[idx]; keep = F.sum(1) >= 30; idx, F = idx[keep], F[keep]; feats = shape_feats(F, c, dz, pk, tl)
    for k, v in extra.items(): feats[k] = v.reindex(idx).to_numpy()
    y = cls.ravel()[idx].astype(int); s_ = inten.ravel()[idx]; fin = np.all(np.isfinite(feats.to_numpy(float)), 1) & np.isfinite(s_); feats, y, s_, idx = feats[fin], y[fin], s_[fin], idx[fin]; fold = blocks(idx)
    print(f"\n{'=' * 100}\n{name}: {idx.size:,} cells, {y.sum():,} in class; intensity index = NIR / floodplain-median NIR, quartile edges {np.percentile(s_, [25, 50, 75]).round(2)}")
    qe = np.percentile(s_, [25, 50, 75]); qm = [(s_ >= lo) & (s_ < hi) for lo, hi in zip([-9] + list(qe), list(qe) + [9])]
    print(f"  {'feature':<10}{'in class':>16}{'out':>16}{'d':>6}{'r(NIR)':>8}   by NIR-index quartile Q1..Q4 (median)")
    for col in feats.columns:
        v = feats[col].to_numpy(float); print(f"  {col:<10}{np.nanmedian(v[y==1]):>9.2f} ({nmad(v[y==1]):4.2f}){np.nanmedian(v[y==0]):>9.2f} ({nmad(v[y==0]):4.2f}){dprime(v[y==1], v[y==0]):>6.2f}{np.corrcoef(v, s_)[0,1]:>8.3f}   " + "  ".join(f"{np.nanmedian(v[m]):7.2f}" for m in qm))
    shape_cols = [col for col in feats.columns if col not in extra]; Xs = feats[shape_cols].to_numpy(float); Xa = feats.to_numpy(float)
    print(f"  PRESENCE AUC (block CV): shape only {cv_auc(Xs, y, fold):.3f}; shape + {list(extra)} {cv_auc(Xa, y, fold):.3f}")
    print(f"  INTENSITY R2 (block CV) of the NIR index: shape only {cv_r2(Xs, s_, fold):.3f}; all {cv_r2(Xa, s_, fold):.3f}")
run("gen1 (Nov 2008), 20 mm bins, pooled 3x3, open floodplain, NE patch excluded", H1p, ctr, DZ, 50, 100, p08 > 0.5, relnir08, {})
run("gen2 (May 2021), 5 mm bins, own cell, open floodplain, NE patch excluded", H2full, ctr2, DZ2, 20, 50, p21 > 0.5, relnir21, {"c2_share": G2s.c2_share})

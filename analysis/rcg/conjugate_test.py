#!/usr/bin/env python3
"""The conjugate test: cells in the reed-canary-grass-type class in 2021 but NOT in 2008
(and the reverse). If the gen1 fingerprint and offset belong to the sward's 2008 STATE and
not to the place, the 2021-only cells must look like ordinary floodplain in gen1, and the
2008-only cells must carry the fingerprint even though the 2021 imagery no longer flags
them.

Four groups on open flat floodplain, NE patch excluded: neither / 2008 only / 2021 only /
both. Per group: gen1 shape about its own mode (peak share, upper width, upper-tail mass),
gen1 level vs gen2 (median, share at/below the surface), the raw and corrected DoD, and
gen2's May-2021 class-2 share. d' of each group against 'neither'. Then the gen1 shape-only
LDA trained on 2008-class vs neither is APPLIED to the 2021-only cells: its score there
should sit with 'neither' if the fingerprint is the 2008 state.

    ./lidar-icp/bin/python analysis/rcg/conjugate_test.py
"""
import sys

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
sl = np.load(f"{S}/slope.npy"); lod = np.load(f"{S}/lod.npy") * 1e3; pf = np.load(f"{S}/penetration_failure.npy").astype(bool); ny, nx = dod.shape
flat = fp & (sl <= 2) & np.isfinite(dod); below = fp & np.isfinite(new) & (new < -lod); lab, n = label(below, structure=np.ones((3, 3))); sizes = np.bincount(lab.ravel()); NE = lab == (np.argsort(sizes[1:])[::-1][0] + 1)
yy, xx = np.mgrid[0:ny, 0:nx]; E = X0 + (xx + 0.5) * res; N = Y0 + (yy + 0.5) * res
def sample(tif):
    da = rioxarray.open_rasterio(tif).squeeze(drop=True); xs, ys = da.x.values, da.y.values; r = abs(float(xs[1] - xs[0]))
    ix = np.clip(((E - (xs[0] - r / 2)) / r).astype(int), 0, len(xs) - 1); iy = np.clip(((ys[0] + r / 2 - N) / r).astype(int), 0, len(ys) - 1); return da.values[iy, ix]
p08 = sample(f"{R}/nelike_p_ls_20081029.tif"); p21 = sample(f"{R}/nelike_p_s2_20211018.tif"); p21L = sample(f"{R}/nelike_p_ls_20211017.tif")
zc = np.load(f"{R}/gen2_cellstats.npz"); G2s = pd.DataFrame({k: zc[k] for k in zc.files if k != "cells"}, index=zc["cells"])
open5 = np.zeros(ny * nx, bool); open5[G2s.index[(G2s.above2m == 0)].to_numpy()] = True; open5 = open5.reshape(ny, nx)
import os
USE_L8 = os.environ.get("CONJ_2021", "s2") == "l8"           # CONJ_2021=l8 -> 2021 class from Landsat-8 30 m (same resolution as 2008)
dom = flat & open5 & ~NE; C08 = p08 > 0.5; C21 = (p21L if USE_L8 else p21) > 0.5
print("2021 class source:", "Landsat-8 2021-10-17, 30 m (like-with-like)" if USE_L8 else "Sentinel-2 2021-10-18, 10 m")
LO, HI, DZ = -1500.0, 1500.0, 20.0; edges = np.arange(LO, HI + DZ, DZ); ctr = 0.5 * (edges[:-1] + edges[1:]); NB = ctr.size
t = pq.read_table(f"{S}/beam_offset_table.parquet", columns=["cell", "d_mm_corr", "in_grid"]); g = t["in_grid"].to_numpy().astype(bool); cell = t["cell"].to_numpy()[g]; d = t["d_mm_corr"].to_numpy()[g].astype(float); ok = np.isfinite(d) & (d > LO) & (d < HI)
H = np.zeros((ny * nx, NB), np.int32); np.add.at(H, (cell[ok], np.clip(((d[ok] - LO) / DZ).astype(int), 0, NB - 1)), 1); Hp = np.rint(uniform_filter(H.reshape(ny, nx, NB).astype(np.float32), size=(3, 3, 1), mode="nearest") * 9).reshape(-1, NB)
idx = np.flatnonzero(dom.ravel()); F = Hp[idx]; keep = F.sum(1) >= 30; idx, F = idx[keep], F[keep]
tot = F.sum(1); P = F / tot[:, None]; C = np.cumsum(P, 1); q = lambda p: ctr[np.argmax(C >= p, axis=1)]
mode = np.array([hsm_from_hist(F[i], ctr) for i in range(len(F))]); rel = ctr[None, :] - mode[:, None]
feat = pd.DataFrame(dict(peak=(P * (np.abs(rel) <= 50)).sum(1), hi_w=q(.9) - q(.5), lo_w=q(.5) - q(.1), tail_mass=(P * (rel > 100)).sum(1), low_mass=(P * (rel < -100)).sum(1), mode_med=mode - q(.5),
                         tail_L=((P * (rel > 100)) * (rel - 100)).sum(1) / np.maximum((P * (rel > 100)).sum(1), 1e-9), bowley=(q(.75) + q(.25) - 2 * q(.5)) / np.maximum(q(.75) - q(.25), DZ),
                         d_med=q(.5), share_le0=(F[:, ctr <= 0].sum(1) / tot), raw_dod=dod.ravel()[idx], corr_dod=new.ravel()[idx], penfail=pf.ravel()[idx].astype(float), g2_c2=G2s.c2_share.reindex(idx).to_numpy(), g2_h50=G2s.h_p50.reindex(idx).to_numpy()))
c08 = C08.ravel()[idx]; c21 = C21.ravel()[idx]
groups = [("neither", ~c08 & ~c21), ("2008 only (lost by 2021)", c08 & ~c21), ("2021 only (gained since 2008)", ~c08 & c21), ("both years", c08 & c21)]
nmad = lambda v: 1.4826 * np.nanmedian(np.abs(v - np.nanmedian(v)))
def dprime(a, b): return abs(np.nanmedian(a) - np.nanmedian(b)) / np.sqrt((nmad(a) ** 2 + nmad(b) ** 2) / 2 + 1e-12)
print(f"open flat floodplain, NE patch excluded: {idx.size:,} cells; groups: " + ", ".join(f"{nm} {m.sum():,}" for nm, m in groups))
print(f"  (2021 class from Sentinel-2 10 m; with Landsat-8 2021-10-17 instead the 2021-only group would be {((~c08) & (p21L.ravel()[idx] > 0.5)).sum():,} cells)")
base = feat[groups[0][1]]
print(f"\n  {'feature':<10}" + "".join(f"{nm[:22]:>24}" for nm, _ in groups) + "   d' vs neither: 2008only  2021only  both")
for k in feat.columns:
    v = feat[k]; line = f"  {k:<10}"
    for nm, m in groups: line += f"{np.nanmedian(v[m]):>14.2f} ({nmad(v[m]):6.2f})"
    line += f"   {dprime(v[groups[1][1]], v[groups[0][1]]):>19.2f}{dprime(v[groups[2][1]], v[groups[0][1]]):>10.2f}{dprime(v[groups[3][1]], v[groups[0][1]]):>6.2f}"
    print(line)
shape_cols = ["peak", "hi_w", "lo_w", "tail_mass", "low_mass", "mode_med", "tail_L", "bowley"]
X = feat[shape_cols].to_numpy(float); fin = np.all(np.isfinite(X), 1)
tr = fin & (groups[0][1] | c08); m = LinearDiscriminantAnalysis().fit(X[tr], c08[tr].astype(int)); sc = m.decision_function(X)
print("\n  gen1 shape-only LDA, trained 2008-class vs neither (resubstitution for the training groups): median score (NMAD) per group")
for nm, mm in groups: print(f"    {nm:<32} {np.nanmedian(sc[mm & fin]):>7.2f} ({nmad(sc[mm & fin]):5.2f})")
a = roc_auc_score(np.r_[np.zeros((groups[0][1] & fin).sum()), np.ones((groups[2][1] & fin).sum())], np.r_[sc[groups[0][1] & fin], sc[groups[2][1] & fin]])
b = roc_auc_score(np.r_[np.zeros((groups[0][1] & fin).sum()), np.ones((groups[1][1] & fin).sum())], np.r_[sc[groups[0][1] & fin], sc[groups[1][1] & fin]])
print(f"  AUC of that score, 2021-only vs neither: {a:.3f}  (0.5 = the 2021-only cells look like 'neither' in 2008);  2008-only vs neither: {b:.3f}")

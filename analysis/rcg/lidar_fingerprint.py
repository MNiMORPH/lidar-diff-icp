#!/usr/bin/env python3
"""The lidar fingerprint of the reed-canary-grass-type class, per survey, from the imagery
classes of each survey's own year.

Andy, 2026-10-19: "find the patches with and without in each time, and then find the
characteristic LiDAR signature from each time ... see if there is a fingerprint of reed
canary grass in each LiDAR survey return set."

DESIGN. The imagery class ("NE-like", map_ne_like.py / map_ne_like_2021.py) is defined
WITHOUT the lidar; the lidar statistics are then read inside and outside it. That is the
right direction: no lidar quantity trains the class.
  gen1 (Nov 2008) is matched to NE-like on Landsat-5 2008-10-29 (30 m, read at 5 m cells).
  gen2 (May 2021) is matched to NE-like on Sentinel-2 2021-10-18 (10 m).
Population: flat floodplain cells (mask AND slope <= 2 deg). Null for scale: the F3 cropped
cells against the rest, and F1 sedge against the rest.

FEATURES, per 5 m cell.
  gen1, from the beam table (CSF ground returns, slope-normal vs gen2): n_csf, d median,
    p10, p90, span, share <= 0, intensity median; and from the RAW merged cloud: n_all,
    class-2 share, class-1 (unclassified) share, single-return share, p90 of height above
    gen2's surface, share above +0.5 m.
  gen2, from the full-density 2021 cloud streamed once (cached to rcg/gen2_cellstats.npz):
    n_all, class-2 share, height above gen2's own surface p10/p50/p90, share <= 0, share in
    0.15-2 m (lowveg), share > 2 m, single-return share, intensity median.
SEPARABILITY d' per feature; then a lidar-ONLY classifier (LDA on the survey's features)
predicting the imagery class, 5-fold CV by 500 m block, reported as AUC -- "can the lidar
see it alone". Then the OFFSETS: raw DoD and corrected DoD inside/outside NE-like 2008.

    ./lidar-icp/bin/python analysis/rcg/lidar_fingerprint.py
"""
import os

import laspy
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
import rioxarray
from scipy.ndimage import label
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.metrics import roc_auc_score

S = "data/derived/elbaext_delong_star_nosurf"; R = f"{S}/rcg"
X0, Y0, X1, Y1 = 575600.0, 4882200.0, 580050.0, 4886250.0; res = 5.0
dod = np.load(f"{S}/dod.npy") * 1e3; new = np.load(f"{S}/dod_floodplain_vegcorr.npy") * 1e3; D = np.load(f"{S}/floodplain_vegcorr_delta.npy") * 1e3
fp = np.load(f"{S}/floodplain_mask.npy").astype(bool); sl = np.load(f"{S}/slope.npy"); lod = np.load(f"{S}/lod.npy") * 1e3
sp = np.load(f"{S}/gen1_spread_mm.npy"); pf = np.load(f"{S}/penetration_failure.npy").astype(bool); z2 = np.load(f"{S}/z_after.npy"); ny, nx = dod.shape
flat = fp & (sl <= 2) & np.isfinite(dod)
below = fp & np.isfinite(new) & (new < -lod); lab, n = label(below, structure=np.ones((3, 3))); sizes = np.bincount(lab.ravel()); NE = lab == (np.argsort(sizes[1:])[::-1][0] + 1)
F3 = flat & (sp < 140) & (np.abs(dod) < 30)
yy, xx = np.mgrid[0:ny, 0:nx]; E = X0 + (xx + 0.5) * res; N = Y0 + (yy + 0.5) * res
def sample(tif):
    da = rioxarray.open_rasterio(tif).squeeze(drop=True); xs, ys = da.x.values, da.y.values; r = abs(float(xs[1] - xs[0]))
    ix = np.clip(((E - (xs[0] - r / 2)) / r).astype(int), 0, len(xs) - 1); iy = np.clip(((ys[0] + r / 2 - N) / r).astype(int), 0, len(ys) - 1)
    return da.values[iy, ix]
p08 = sample(f"{R}/nelike_p_ls_20081029.tif"); p21 = sample(f"{R}/nelike_p_s2_20211018.tif")
C08 = flat & (p08 > 0.5); C21 = flat & (p21 > 0.5)
print(f"flat floodplain cells {flat.sum():,}; NE-like 2008 (Landsat) {C08.sum():,} = {100*C08.sum()/flat.sum():.1f}%; NE-like 2021 (Sentinel-2) {C21.sum():,} = {100*C21.sum()/flat.sum():.1f}%; overlap {(C08 & C21).sum():,}; NE patch cells in 2008 class {100*(C08 & NE).sum()/(NE & flat).sum():.0f}%")

# ---------------- gen1 features
t = pq.read_table(f"{S}/beam_offset_table.parquet", columns=["cell", "d_mm_corr", "in_grid", "intensity"]); g = t["in_grid"].to_numpy().astype(bool)
cell = t["cell"].to_numpy()[g]; d = t["d_mm_corr"].to_numpy()[g].astype(float); it = t["intensity"].to_numpy()[g].astype(float); ok = np.isfinite(d)
df = pd.DataFrame(dict(cell=cell[ok], d=d[ok], it=it[ok])); gb = df.groupby("cell")
G1 = pd.DataFrame(dict(n_csf=gb.size(), d_med=gb.d.median(), d_p10=gb.d.quantile(0.10), d_p90=gb.d.quantile(0.90), share_le0=gb.d.apply(lambda v: np.mean(v <= 0)), csf_int=gb.it.median()))
G1["span"] = G1.d_p90 - G1.d_p10
f = laspy.read("data/before/elbaext_gen1_merged.laz"); x = np.asarray(f.x); y = np.asarray(f.y); z = np.asarray(f.z); cl = np.asarray(f.classification); nr = np.asarray(f.number_of_returns); inten = np.asarray(f.intensity).astype(float)
ix = ((x - X0) / res).astype(int); iy = ((y - Y0) / res).astype(int); ing = (ix >= 0) & (ix < nx) & (iy >= 0) & (iy < ny) & (cl != 7); c1 = (iy * nx + ix)[ing]
h = (z[ing] - z2.ravel()[c1]) * 1e3
raw = pd.DataFrame(dict(cell=c1, h=h, c2=(cl[ing] == 2), c1u=(cl[ing] == 1), single=(nr[ing] == 1), it=inten[ing])); gb = raw.groupby("cell")
G1r = pd.DataFrame(dict(n_all=gb.size(), c2_share=gb.c2.mean(), c1_share=gb.c1u.mean(), single_share=gb.single.mean(), h_p90_all=gb.h.quantile(0.90), above05=gb.h.apply(lambda v: np.mean(v > 500)), all_int=gb.it.median()))
G1 = G1.join(G1r, how="outer"); del f, x, y, z, raw
# ---------------- gen2 features (streamed once, cached)
cache = f"{R}/gen2_cellstats.npz"
if not os.path.exists(cache):
    acc = {}
    with laspy.open("data/after/3dep2021_fulldensity.laz") as f2:
        print(f"streaming gen2 full density: {f2.header.point_count:,} points", flush=True)
        parts = []
        for ch in f2.chunk_iterator(8_000_000):
            cx = np.asarray(ch.x); cy = np.asarray(ch.y); cz = np.asarray(ch.z); ccl = np.asarray(ch.classification); cnr = np.asarray(ch.number_of_returns); cit = np.asarray(ch.intensity).astype(float)
            jx = ((cx - X0) / res).astype(int); jy = ((cy - Y0) / res).astype(int); okc = (jx >= 0) & (jx < nx) & (jy >= 0) & (jy < ny) & (ccl != 7) & (ccl != 18)
            cc = (jy * nx + jx)[okc]; keep = flat.ravel()[cc]
            cc = cc[keep]; hh = (cz[okc][keep] - z2.ravel()[cc]) * 1e3
            parts.append(pd.DataFrame(dict(cell=cc, h=hh, c2=(ccl[okc][keep] == 2), single=(cnr[okc][keep] == 1), it=cit[okc][keep])))
    raw2 = pd.concat(parts); gb = raw2.groupby("cell")
    G2 = pd.DataFrame(dict(n_all=gb.size(), c2_share=gb.c2.mean(), h_p10=gb.h.quantile(0.10), h_p50=gb.h.median(), h_p90=gb.h.quantile(0.90), share_le0=gb.h.apply(lambda v: np.mean(v <= 0)),
                           lowveg=gb.h.apply(lambda v: np.mean((v > 150) & (v <= 2000))), above2m=gb.h.apply(lambda v: np.mean(v > 2000)), single_share=gb.single.mean(), int_med=gb.it.median()))
    np.savez_compressed(cache, cells=G2.index.to_numpy(), **{k: G2[k].to_numpy() for k in G2.columns}); del raw2, parts
z2c = np.load(cache); G2 = pd.DataFrame({k: z2c[k] for k in z2c.files if k != "cells"}, index=z2c["cells"])

nmad = lambda v: 1.4826 * np.nanmedian(np.abs(v - np.nanmedian(v)))
def dprime(a, b): return abs(np.nanmedian(a) - np.nanmedian(b)) / np.sqrt((nmad(a) ** 2 + nmad(b) ** 2) / 2 + 1e-12)
def report(G, name, cls, extra_null, domain=None, exclude_ne=False, drop=()):
    """domain: 5 m boolean restricting BOTH sides (e.g. open floodplain); exclude_ne: drop the NE
    patch from the class, so the class members are the ones the lidar did not help define."""
    dom = flat if domain is None else (flat & domain)
    idx = np.flatnonzero(dom.ravel()); G = G.reindex(idx).drop(columns=[c for c in drop if c in G.columns]); inC = cls.ravel()[idx]; f3 = F3.ravel()[idx]; nep = NE.ravel()[idx]
    if exclude_ne:
        keep = ~nep; idx, inC, f3, nep = idx[keep], inC[keep], f3[keep], nep[keep]; G = G.iloc[keep]
    print(f"\n  [domain cells {idx.size:,}; class members {inC.sum():,}" + ("; NE patch EXCLUDED" if exclude_ne else "") + "]")
    print(f"\n{'=' * 110}\n{name}: per-cell features, flat floodplain; medians (NMAD) inside / outside the imagery class, and d'  [null: F3 cropped vs rest]")
    print(f"  {'feature':<13}{'in class':>18}{'out of class':>18}{'d in/out':>10}{'NE patch':>18}{'F3 cropped':>18}{'d F3/rest':>11}")
    rows = []
    for k in G.columns:
        v = G[k].to_numpy(float); a, b = v[inC], v[~inC]; dd = dprime(a, b); dn = dprime(v[f3], v[~f3])
        print(f"  {k:<13}{np.nanmedian(a):>10.2f} ({nmad(a):5.2f}){np.nanmedian(b):>10.2f} ({nmad(b):5.2f}){dd:>10.2f}{np.nanmedian(v[nep]):>10.2f} ({nmad(v[nep]):5.2f}){np.nanmedian(v[f3]):>10.2f} ({nmad(v[f3]):5.2f}){dn:>11.2f}")
        rows.append(dict(survey=name, feature=k, in_med=np.nanmedian(a), in_nmad=nmad(a), out_med=np.nanmedian(b), out_nmad=nmad(b), d_in_out=dd, ne_med=np.nanmedian(v[nep]), f3_med=np.nanmedian(v[f3]), d_f3=dn))
    # lidar-only classifier, block CV
    Xf = G.to_numpy(float); okr = np.all(np.isfinite(Xf), 1); Xf = Xf[okr]; yv = inC[okr].astype(int)
    blk = ((idx[okr] // nx) // 100) * 1000 + ((idx[okr] % nx) // 100); ub, binv = np.unique(blk, return_inverse=True); fold = np.random.default_rng(0).permutation(np.arange(ub.size) % 5)[binv]
    pred = np.zeros(yv.size)
    for k_ in range(5):
        tr = fold != k_; m = LinearDiscriminantAnalysis().fit(Xf[tr], yv[tr]); pred[~tr] = m.decision_function(Xf[~tr])
    auc = roc_auc_score(yv, pred); m = LinearDiscriminantAnalysis().fit(Xf, yv); coef = m.coef_[0] * Xf.std(0)
    print(f"  LIDAR-ONLY LDA predicting the imagery class, 5-fold CV by 500 m block: AUC {auc:.3f} (0.5 = nothing, 1 = perfect); n {yv.size:,}, positives {yv.sum():,}")
    print("  standardised LDA weights: " + ", ".join(f"{c} {w:+.2f}" for c, w in sorted(zip(G.columns, coef), key=lambda t: -abs(t[1]))))
    return pd.DataFrame(rows), auc
# OPEN floodplain: no gen2 return above 2 m in the cell (gen2 full density sees every tree), so the
# comparison is sward vs sward, not sward vs woods. Defined from gen2 for BOTH surveys (same cells).
open5 = np.zeros(ny * nx, bool); open5[G2.index[(G2.above2m == 0)].to_numpy()] = True; open5 = open5.reshape(ny, nx)
print(f"\nopen floodplain (flat, no gen2 return above 2 m): {(flat & open5).sum():,} of {flat.sum():,} flat cells")
drop1 = ("c1_share", "single_share", "above05"); drop2 = ("int_med",)
print("\n" + "#" * 110 + "\n# A. ALL flat floodplain, NE patch INCLUDED (the first pass)")
r1, auc1 = report(G1, "gen1 (Nov 2008) vs NE-like 2008 (Landsat-5 2008-10-29)", C08, None, drop=drop1)
r2, auc2 = report(G2, "gen2 (May 2021, full density) vs NE-like 2021 (Sentinel-2 2021-10-18)", C21, None, drop=drop2)
print("\n" + "#" * 110 + "\n# B. OPEN flat floodplain only, NE patch EXCLUDED -- the independent fingerprint")
r1b, auc1b = report(G1, "gen1 (Nov 2008) vs NE-like 2008, open, no NE patch", C08, None, domain=open5, exclude_ne=True, drop=drop1)
r2b, auc2b = report(G2, "gen2 (May 2021) vs NE-like 2021, open, no NE patch", C21, None, domain=open5, exclude_ne=True, drop=drop2)
pd.concat([r1.assign(pass_="A"), r2.assign(pass_="A"), r1b.assign(pass_="B"), r2b.assign(pass_="B")]).to_csv(f"{R}/lidar_fingerprint.csv", index=False)
print(f"\n{'=' * 110}\nOFFSETS by imagery class, flat floodplain (mm): median (NMAD)")
for nm, cls in (("NE-like 2008", C08), ("NE-like 2008, no NE patch", C08 & ~NE), ("NE-like 2021", C21), ("NE-like 2021, no NE patch", C21 & ~NE), ("both years", C08 & C21), ("both, no NE patch", C08 & C21 & ~NE), ("neither, open", flat & open5 & ~C08 & ~C21), ("neither", flat & ~C08 & ~C21)):
    m = cls & flat; print(f"  {nm:<14} n {m.sum():>6,}  raw DoD {np.median(dod[m]):>7.1f} ({nmad(dod[m]):5.1f})  corrected DoD {np.median(new[m]):>7.1f} ({nmad(new[m]):5.1f})  gen1 span {np.nanmedian(sp[m]):>5.0f}  penfail {100*pf[m].mean():>4.0f}%")
print(f"table: {R}/lidar_fingerprint.csv; cache: {cache}")

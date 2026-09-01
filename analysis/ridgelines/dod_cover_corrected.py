#!/usr/bin/env python3
"""Rebuild the DEM of difference with the canopy-cover correction applied to gen2.

gen2's ground is taken at the cover-dependent percentile of its own near-ground return
column instead of at the median:

    q2(c) = 0.5 - 0.1922 * c        c = canopy cover fraction  (Q2_COVER_RELATION.md)

gen1's ground is its median, with the four registration terms applied. Both are per-cell
quantiles of the slope-normal residual to the same plane, so

    DoD = (h2(q2) - h1(0.50)) * |n|          positive = elevation ROSE 2008 -> 2021

Requires one streaming pass over the gen2 cloud to build the class-2 near-ground histogram
for EVERY grid cell (the existing class-split cube covers only the divide reference cells).

    env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python \\
        analysis/ridgelines/dod_cover_corrected.py
"""
import argparse, json, os
import numpy as np, laspy, pyarrow.parquet as pq
from lidar_diff_icp import registration as reg
from scipy.ndimage import distance_transform_edt

ap = argparse.ArgumentParser()
ap.add_argument("--tile", default="data/derived/elba_fulldensity")
ap.add_argument("--gen2", default="data/after/3dep2021_fulldensity.laz")
ap.add_argument("--slope", type=float, default=None,
                help="q2 = 0.5 + slope * cover. Default: READ from the tile's own "
                     "q2_cover_fit.json (analysis/ridgelines/q2_cover_fit.py), never typed. "
                     "It used to default to -0.1922, which was Elba's number AND stale -- "
                     "the shipped dod_cover_q2.json records -0.1835 -- so on any other "
                     "region it silently applied Elba's correction. Refuses if absent.")
ap.add_argument("--chunk", type=int, default=3_000_000)
ap.add_argument("--zlo", type=float, default=-1.0)
ap.add_argument("--zhi", type=float, default=2.0)
ap.add_argument("--dz", type=float, default=0.02)
A = ap.parse_args()

D = A.tile


def _q2_slope(tile_dir):
    """The q2 slope from THIS tile's own fit. The relation is per-site -- it depends on each
    pair's phenology and undergrowth -- so there is no site-invariant value to default to."""
    p = os.path.join(tile_dir, "q2_cover_fit.json")
    if not os.path.exists(p):
        raise SystemExit(
            f"no {p}. The q2 slope is NOT defaulted: it is per-site, and a value carried "
            f"from another tile would be applied here without saying so. Produce it:\n"
            f"    env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python "
            f"analysis/ridgelines/q2_cover_fit.py --tile {tile_dir}\n"
            f"or state your own with --slope.")
    j = json.load(open(p))
    st, pop = j["settings"], j["population"]
    print(f"q2 slope {j['linear_slope']:+.4f}, read from {p}")
    print(f"  fitted on {pop['cells_fitted']:,} cells in {pop['bins_used_in_fit']} of "
          f"{pop['bins']} cover bins; weight={st['weight']}, include_valley="
          f"{st['include_valley']}, valley_top_m={st['valley_top_m']}")
    mt = j.get("inputs_mtime", {})
    if "corrections.json" in mt and "beam_offset_table.parquet" in mt:
        if mt["beam_offset_table.parquet"] < mt["corrections.json"]:
            print(f"  WARNING: that fit read a beam_offset_table ({mt['beam_offset_table.parquet']}) "
                  f"OLDER than corrections.json ({mt['corrections.json']}), so its gen1 "
                  f"offsets carry registration terms that are no longer in force.")
    return float(j["linear_slope"])


SLOPE = A.slope if A.slope is not None else _q2_slope(D)
j = reg.read_corrections(D)   # geoid sidecar wins where a tile carries both;
                              # reading "corrections.json" by name picks elbaext's
                              # obsolete reference_plane product
b = j["bounds"]; RES = float(j["res_m"]); X0, Y0 = b[0], b[1]
zf = np.load(f"{D}/z_after.npy"); NY, NX = zf.shape; NC = zf.size
_zf = zf.copy(); _m = ~np.isfinite(_zf)
if _m.any():
    _zf = _zf[tuple(distance_transform_edt(_m, return_distances=False, return_indices=True))]
gy, gx = np.gradient(_zf, RES)
gxf = gx.ravel(); gyf = gy.ravel(); zflat = _zf.ravel()
nnorm = np.sqrt(gxf**2 + gyf**2 + 1.0)
NZ = int(round((A.zhi - A.zlo) / A.dz))

# ---- gen2: class-2 near-ground histogram for every cell -----------------------------
H = np.zeros((NC, NZ), np.int32)
n_in = 0
with laspy.open(A.gen2) as f:
    for pts in f.chunk_iterator(A.chunk):
        cl = np.asarray(pts.classification)
        keep = cl == 2
        x = np.asarray(pts.x)[keep]; y = np.asarray(pts.y)[keep]; z = np.asarray(pts.z)[keep]
        ix = ((x - X0) / RES).astype(np.int64); iy = ((y - Y0) / RES).astype(np.int64)
        ing = (ix >= 0) & (ix < NX) & (iy >= 0) & (iy < NY)
        cc = iy[ing] * NX + ix[ing]
        xc = X0 + ((cc % NX) + 0.5) * RES; yc = Y0 + ((cc // NX) + 0.5) * RES
        h = (z[ing] - (zflat[cc] + gxf[cc] * (x[ing] - xc) + gyf[cc] * (y[ing] - yc))) / nnorm[cc]
        zi = np.floor((h - A.zlo) / A.dz).astype(np.int64)
        m = (zi >= 0) & (zi < NZ)
        np.add.at(H, (cc[m], zi[m]), 1)
        n_in += int(m.sum())
print(f"gen2: {n_in:,} class-2 returns in the window over {int((H.sum(1) > 0).sum()):,} cells",
      flush=True)

# ---- gen1: per-cell median of the registered per-return offsets ----------------------
t = pq.read_table(f"{D}/beam_offset_table.parquet", columns=["cell", "d_mm_corr", "in_grid"])
g = t["in_grid"].to_numpy().astype(bool)
ce = t["cell"].to_numpy()[g]; dc = t["d_mm_corr"].to_numpy()[g].astype(float)
o = np.argsort(ce, kind="stable"); cs = ce[o]; ds = dc[o]
u, st = np.unique(cs, return_index=True); sp = np.r_[st[1:], cs.size]
h1 = np.full(NC, np.nan)
for uu, a, bb in zip(u, st, sp):
    h1[uu] = np.median(ds[a:bb])
print(f"gen1: ground on {int(np.isfinite(h1).sum()):,} cells", flush=True)

# ---- gen2 at the cover-dependent percentile -----------------------------------------
cover = np.load(f"{D}/canopy_cover_pfs.npy").ravel()
q2 = 0.5 + SLOPE * np.where(np.isfinite(cover), cover, 0.0)
C = np.cumsum(H, 1).astype(float); ntot = C[:, -1]
have = ntot > 0
idx = np.arange(NC)
r = q2 * ntot
k = (C >= r[:, None]).argmax(1)
below = np.where(k > 0, C[idx, np.maximum(k - 1, 0)], 0.0)
inbin = C[idx, k] - below
frac = np.where(inbin > 0, (r - below) / np.maximum(inbin, 1e-9), 0.0)
h2 = np.where(have, (A.zlo + (k + np.clip(frac, 0, 1)) * A.dz) * 1000.0, np.nan)
h2_med = np.where(have, (A.zlo + ((C >= 0.5 * ntot[:, None]).argmax(1) + 0.5) * A.dz) * 1000.0,
                  np.nan)

dod_corr = (h2 - h1) / 1000.0 * nnorm            # m, positive = elevation rose
dod_med = (h2_med - h1) / 1000.0 * nnorm         # same but gen2 at its plain median
np.save(f"{D}/dod_cover_q2.npy", dod_corr.reshape(NY, NX))
np.save(f"{D}/dod_gen2_median.npy", dod_med.reshape(NY, NX))
np.save(f"{D}/gen2_q2_used.npy", np.where(have, q2, np.nan).reshape(NY, NX))
json.dump({"relation": "q2 = 0.5 + slope * cover", "slope": SLOPE,
           "source": "analysis/ridgelines/Q2_COVER_RELATION.md",
           "gen1": "beam_offset_table.parquet median of d_mm_corr (4 registration terms)",
           "gen2": f"{A.gen2}, class-2, near-ground column {A.zlo}..{A.zhi} m, {A.dz} m bins",
           "sign": "gen2 - gen1, positive = elevation rose",
           "cells_with_both": int((np.isfinite(dod_corr)).sum())},
          open(f"{D}/dod_cover_q2.json", "w"), indent=2)
ok = np.isfinite(dod_corr) & np.isfinite(dod_med)
print(f"\nwrote {D}/dod_cover_q2.npy ({int(ok.sum()):,} cells)")
print(f"      {D}/dod_gen2_median.npy   {D}/gen2_q2_used.npy   {D}/dod_cover_q2.json")
print(f"\nmedian DoD (mm): gen2 median {np.median(dod_med[ok])*1000:+.1f}  ->  "
      f"cover-corrected {np.median(dod_corr[ok])*1000:+.1f}")
old = np.load(f"{D}/dod.npy").ravel()
oo = ok & np.isfinite(old)
print(f"existing dod.npy on the same cells: {np.median(old[oo])*1000:+.1f} mm  "
      f"(n={int(oo.sum()):,})")

#!/usr/bin/env python3
"""LOOK AT THE DATA: the near-ground return column of each epoch, against canopy.

No classification, no height threshold, no ground/vegetation decision, no binning of the
covariate beyond what a picture needs. Every return within a few metres of the gen2
bare-earth surface is placed by its height, and the columns are stacked against canopy
fraction so the shape of the near-ground distribution can be read directly.

The point is that near-ground foliage CANNOT be separated return-by-return -- a leaf at
0.2 m and the ground beneath it give one return. But it does not vanish: it lifts and skews
the upper part of the near-ground distribution, while its floor (the pulses that got all the
way down) stays put. Range walk on weak returns does the opposite -- it slides the whole
distribution down, floor included. Those are different pictures, so plotting the
distribution answers what counting a height band cannot.

Each column is normalised to its own maximum so the shape is visible at every canopy
fraction rather than being swamped by the open ground that dominates the tile. Percentile
curves (p05/p25/p50) are read off the same columns -- they summarise the picture, they do
not replace it.

    env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python \
        analysis/ridgelines/nearground_profile.py --tile data/derived/elbaext
"""
import argparse, json, os
import numpy as np, laspy
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

ap = argparse.ArgumentParser()
ap.add_argument("--tile", default="data/derived/elbaext")
ap.add_argument("--gen1", default="data/before/elbaext_gen1_merged.laz")
ap.add_argument("--gen2", default="data/after/elbaext_3dep_fulldensity.laz")
ap.add_argument("--zlo", type=float, default=-1.0)
ap.add_argument("--zhi", type=float, default=3.0)
ap.add_argument("--dz", type=float, default=0.02, help="height resolution of the picture (m)")
ap.add_argument("--ncol", type=int, default=80, help="covariate columns (display resolution)")
ap.add_argument("--chunk", type=int, default=3_000_000)
A = ap.parse_args()
TILE = os.path.basename(A.tile.rstrip("/"))


def grid(tile):
    for fn in ("meta.json", "corrections_geoid.json", "corrections.json"):
        p = f"{tile}/{fn}"
        if os.path.exists(p):
            j = json.load(open(p)); b = j["bounds"]; r = float(j.get("res") or j.get("res_m"))
            return b[0], b[1], r, int(j.get("ny") or round((b[3]-b[1])/r)), \
                   int(j.get("nx") or round((b[2]-b[0])/r))
    raise SystemExit(f"no grid meta in {tile}")


X0, Y0, RES, NY, NX = grid(A.tile)
zflat = np.load(f"{A.tile}/z_after.npy").ravel()
f2008 = np.load(f"{A.tile}/gen1_canopy_frac.npz")["frac"].ravel()
cov21 = np.load(f"{A.tile}/canopy_cover_pfs.npy").ravel()
zedges = np.arange(A.zlo, A.zhi + 0.5*A.dz, A.dz); NZ = zedges.size - 1
cedges = np.linspace(0.0, 1.0, A.ncol + 1)
COVS = {"2008 canopy return fraction": f2008, "2021 canopy cover (PFS)": cov21}


def profile(path):
    H = {k: np.zeros((NZ, A.ncol), np.int64) for k in COVS}
    n = 0
    with laspy.open(path) as f:
        for pts in f.chunk_iterator(A.chunk):
            cl = np.asarray(pts.classification)
            keep = cl != 7
            x = np.asarray(pts.x)[keep]; y = np.asarray(pts.y)[keep]; z = np.asarray(pts.z)[keep]
            ix = ((x - X0)/RES).astype(np.int64); iy = ((y - Y0)/RES).astype(np.int64)
            ok = (ix >= 0) & (ix < NX) & (iy >= 0) & (iy < NY)
            cell = iy[ok]*NX + ix[ok]
            h = z[ok] - zflat[cell]
            zi = np.floor((h - A.zlo)/A.dz).astype(np.int64)
            inz = (zi >= 0) & (zi < NZ)
            for k, cv in COVS.items():
                c = cv[cell]
                ci = np.floor(c*A.ncol).astype(np.int64)
                m = inz & np.isfinite(c) & (ci >= 0) & (ci < A.ncol)
                np.add.at(H[k], (zi[m], ci[m]), 1)
            n += int(ok.sum())
    return H, n


print(f"streaming gen1: {A.gen1}", flush=True)
H1, n1 = profile(A.gen1)
print(f"  {n1:,} returns in grid", flush=True)
print(f"streaming gen2: {A.gen2}", flush=True)
H2, n2 = profile(A.gen2)
print(f"  {n2:,} returns in grid", flush=True)

zc = 0.5*(zedges[:-1] + zedges[1:]); cc = 0.5*(cedges[:-1] + cedges[1:])


def pct(col, q):
    """Percentile height read straight off a column of the picture."""
    tot = col.sum()
    if tot < 200: return np.nan
    cdf = np.cumsum(col)/tot
    return float(np.interp(q, cdf, zc))


fig, axes = plt.subplots(2, 2, figsize=(15, 10), dpi=130, sharey=True)
for j, (cname, _) in enumerate(COVS.items()):
    for i, (H, ep, n) in enumerate(((H1, "gen1  2008 leaf-off", n1), (H2, "gen2  2021 leaf-on", n2))):
        ax = axes[i, j]; M = H[cname].astype(float)
        colmax = M.max(axis=0); colmax[colmax == 0] = 1
        ax.imshow(M/colmax, origin="lower", aspect="auto", cmap="magma",
                  extent=(0, 1, A.zlo, A.zhi), vmin=0, vmax=1)
        for q, st, lab in ((0.05, "-", "p05"), (0.25, "--", "p25"), (0.50, ":", "p50")):
            v = [pct(M[:, k], q) for k in range(A.ncol)]
            ax.plot(cc, v, st, lw=1.6, color="cyan", label=lab)
        ax.axhline(0, color="w", lw=0.8, alpha=.7)
        ax.set_xlim(0, min(1.0, float(np.nanmax(cc[colmax > 200])) if (colmax > 200).any() else 1.0))
        ax.set_title(f"{ep}  —  vs {cname}", fontsize=10)
        if j == 0: ax.set_ylabel("height above gen2 bare earth (m)")
        if i == 1: ax.set_xlabel(cname)
        if i == 0 and j == 0: ax.legend(fontsize=8, loc="upper left")
fig.suptitle(f"near-ground return column, every return, no classification — {TILE}\n"
             f"columns normalised to their own maximum; cyan = percentiles of each column", y=1.0)
out = f"figures/refdatum/nearground_profile_{TILE}.png"
fig.savefig(out, bbox_inches="tight"); plt.close(fig)
print(f"wrote {out}")

print(f"\nfloor (p05) and median (p50) of the near-ground column, read off the data:")
print(f"{'covariate':>28s} {'value':>7s} {'gen1 p05':>9s} {'gen1 p50':>9s} {'gen2 p05':>9s} {'gen2 p50':>9s}")
for cname in COVS:
    for k in (2, 8, 16, 24, 32, 40):
        if k >= A.ncol: continue
        a1, b1 = pct(H1[cname][:, k], .05), pct(H1[cname][:, k], .50)
        a2, b2 = pct(H2[cname][:, k], .05), pct(H2[cname][:, k], .50)
        if not np.isfinite(a1) or not np.isfinite(a2): continue
        print(f"{cname:>28s} {cc[k]:>7.3f} {a1:>9.3f} {b1:>9.3f} {a2:>9.3f} {b2:>9.3f}")

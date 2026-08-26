#!/usr/bin/env python3
"""Per-cell CANOPY RETURN FRACTION for one epoch, computed the SAME way for either epoch.

Fraction of that epoch's OWN returns falling more than --hag metres above the gen2
bare-earth surface (z_after.npy), on the tile's 5 m grid. This exists so the 2008 and 2021
canopies can be compared on an IDENTICAL definition: PyForestScan cover (canopy_cover_pfs)
is a plant-area-density product and is not commensurate with a return fraction, so a
2008-vs-2021 contrast built from one of each would confound phenology with method.

Streamed with laspy's chunk_iterator -- never holds more than one chunk, so it runs on a
laptop against a multi-GB cloud. Classification 7 (noise) and withheld points are dropped.

    env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python \
        analysis/ridgelines/epoch_canopy_frac.py CLOUD.laz OUT.npz [--tile DIR] [--hag 2.0]

Writes n_ret / n_above / frac rasters (NY, NX) to OUT.npz.
"""
import argparse, json, os
import numpy as np, laspy

ap = argparse.ArgumentParser()
ap.add_argument("cloud"); ap.add_argument("out")
ap.add_argument("--tile", default="data/derived/elbaext")
ap.add_argument("--hag", type=float, default=2.0, help="height above gen2 bare earth (m)")
ap.add_argument("--chunk", type=int, default=2_000_000)
A = ap.parse_args()


def grid(tile):
    for fn in ("meta.json", "corrections_geoid.json", "corrections.json"):
        p = f"{tile}/{fn}"
        if os.path.exists(p):
            j = json.load(open(p)); b = j["bounds"]; r = float(j.get("res") or j.get("res_m"))
            ny = int(j.get("ny") or round((b[3] - b[1]) / r))
            nx = int(j.get("nx") or round((b[2] - b[0]) / r))
            return b[0], b[1], r, ny, nx
    raise SystemExit(f"no grid meta in {tile}")


X0, Y0, RES, NY, NX = grid(A.tile)
zg = np.load(f"{A.tile}/z_after.npy")
assert zg.shape == (NY, NX), f"{zg.shape} != {(NY, NX)}"
zflat = zg.ravel()

n_ret = np.zeros(NY * NX, np.int64)
n_above = np.zeros(NY * NX, np.int64)
n_bldg = np.zeros(NY * NX, np.int64)   # ASPRS class 6, to flag structures as false "canopy"
seen = kept = 0
with laspy.open(A.cloud) as f:
    print(f"{A.cloud}: {f.header.point_count:,} points, PF{f.header.point_format.id}")
    for pts in f.chunk_iterator(A.chunk):
        seen += len(pts)
        x = np.asarray(pts.x); y = np.asarray(pts.y); z = np.asarray(pts.z)
        cls = np.asarray(pts.classification)
        good = cls != 7
        if "withheld" in pts.point_format.dimension_names:
            good &= ~np.asarray(pts.withheld).astype(bool)
        ix = np.floor((x - X0) / RES).astype(np.int64)
        iy = np.floor((y - Y0) / RES).astype(np.int64)
        good &= (ix >= 0) & (ix < NX) & (iy >= 0) & (iy < NY)
        if not good.any():
            continue
        c = iy[good] * NX + ix[good]
        zz = z[good]; cc6 = cls[good] == 6
        zref = zflat[c]
        ok = np.isfinite(zref)
        c = c[ok]; zz = zz[ok]; zref = zref[ok]; cc6 = cc6[ok]
        kept += c.size
        n_ret += np.bincount(c, minlength=NY * NX)
        hi = c[zz > zref + A.hag]
        if hi.size:
            n_above += np.bincount(hi, minlength=NY * NX)
        if cc6.any():
            n_bldg += np.bincount(c[cc6], minlength=NY * NX)
        print(f"  ...{seen:,} read, {kept:,} kept", end="\r", flush=True)

frac = np.full(NY * NX, np.nan)
m = n_ret > 0
frac[m] = n_above[m] / n_ret[m]
np.savez_compressed(A.out, n_ret=n_ret.reshape(NY, NX), n_above=n_above.reshape(NY, NX),
                    n_bldg=n_bldg.reshape(NY, NX), frac=frac.reshape(NY, NX),
                    hag=A.hag, cloud=A.cloud)
print(f"\n{kept:,} returns gridded over {int(m.sum()):,} cells "
      f"(median {np.median(n_ret[m]):.0f} returns/cell)")
print(f"canopy fraction (cells with >=10 returns): median "
      f"{np.nanmedian(frac[n_ret >= 10]):.3f}, "
      f"P90 {np.nanpercentile(frac[n_ret >= 10], 90):.3f}")
print(f"wrote {A.out}")

#!/usr/bin/env python3
"""Slope-normal low-percentile ground: remove the downhill bias of horizontal
low-percentile ground estimation on sloped terrain.

Horizontal low-percentile ground picks the lowest points in an x,y cell, which on
a slope are the physically *downhill* ones -- a bias of ~(offset * slope). Taking
the low percentile *normal to the local surface* (detrend each point by a common
regional slope before the low-pick) removes that bias. Applied to BOTH epochs
against the SAME regional plane, the tilt cancels in the difference and neither
side carries a residual bias to mismatch.

On the Elba pilot this cut the DoD band-pass 22% and the total stable scatter 31%,
while a synthetic +0.80 m change was recovered identically (0.726 m) -- a targeted
debias, not smoothing. See analysis/banding_source_investigation.md.

Memory note: the full 2021 cloud is ~107 M last-return points; the per-cell
percentile over all of them peaks near 6 GB. A production version should stream
the quantile (chunk by cell block) rather than materialise all residuals.

Usage:
    ./lidar-icp/bin/python analysis/slope_normal_ground.py \
        data/before/4342-29-64.laz data/after/3dep2021_last.laz \
        --bounds 577492.8 4886237.6 --res 5 --shape 700 508
"""
from __future__ import annotations

import argparse

import laspy
import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter


def _last_return(path):
    las = laspy.read(path)
    rn = np.asarray(las.return_number)
    nr = np.asarray(las.number_of_returns)
    last = rn == nr
    return (np.asarray(las.x)[last], np.asarray(las.y)[last],
            np.asarray(las.z)[last].astype(np.float64))


def _cell_index(x, y, X0, Y0, res, H, W):
    cj = ((x - X0) / res).astype(np.int32)     # column (east)
    ci = ((Y0 - y) / res).astype(np.int32)     # row (north at top)
    ok = (cj >= 0) & (cj < W) & (ci >= 0) & (ci < H)
    return ci, cj, ok


def horizontal_low(path, X0, Y0, res, H, W, q=0.10):
    """Current method: low percentile of raw z per horizontal cell."""
    x, y, z = _last_return(path)
    ci, cj, ok = _cell_index(x, y, X0, Y0, res, H, W)
    f = (ci[ok].astype(np.int64) * W + cj[ok])
    g = pd.Series(z[ok]).groupby(f).quantile(q)
    G = np.full(H * W, np.nan)
    G[g.index.values] = g.values
    return G.reshape(H, W)


def slope_normal_low(path, plane, res, X0, Y0, H, W, q=0.10):
    """Low percentile of the surface-normal residual per cell.

    ``plane`` is (Z_reg, dz_deast, dz_dnorth, cos_normal) flattened per cell, from a
    lightly smoothed reference surface shared by both epochs (so it cancels in the
    difference).
    """
    Zr, dz_de, dz_dn, cosn = plane
    x, y, z = _last_return(path)
    ci, cj, ok = _cell_index(x, y, X0, Y0, res, H, W)
    x, y, z, ci, cj = x[ok], y[ok], z[ok], ci[ok], cj[ok]
    f = ci.astype(np.int64) * W + cj
    dxe = x - (X0 + (cj + 0.5) * res)          # offset east of cell centre
    dyn = y - (Y0 - (ci + 0.5) * res)          # offset north of cell centre
    resid = (z - (Zr[f] + dxe * dz_de[f] + dyn * dz_dn[f])) * cosn[f]
    g = pd.Series(resid.astype(np.float32)).groupby(f).quantile(q)
    G = np.full(H * W, np.nan)
    G[g.index.values] = g.values
    return G.reshape(H, W)


def regional_plane(z_ref_grid, res, smooth=1.2):
    """Common regional slope frame from a reference surface (e.g. dense-epoch low-q)."""
    m = np.isfinite(z_ref_grid)
    zf = z_ref_grid.copy()
    zf[~m] = np.nanmean(z_ref_grid[m])
    Zreg = gaussian_filter(zf, smooth)
    gyr, gxr = np.gradient(Zreg, res)          # per-row (south+), per-col (east+)
    cosn = 1.0 / np.sqrt(1.0 + (gxr ** 2 + gyr ** 2))
    return Zreg.ravel(), gxr.ravel(), (-gyr).ravel(), cosn.ravel()


def _bandpass(a, m, s_lo=1.6, s_hi=8.0):
    w = m.astype(float)
    lo = gaussian_filter(np.where(m, a, 0) * w, s_lo) / np.maximum(gaussian_filter(w, s_lo), 1e-6)
    hi = gaussian_filter(np.where(m, a, 0) * w, s_hi) / np.maximum(gaussian_filter(w, s_hi), 1e-6)
    return np.where(m, lo - hi, np.nan)


def _rstd(v):
    v = v[np.isfinite(v)]
    return 1.4826 * np.median(np.abs(v - np.median(v)))


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("before_laz")
    p.add_argument("after_laz")
    p.add_argument("--bounds", type=float, nargs=2, required=True,
                   metavar=("X0", "Y0"), help="upper-left corner (easting, northing)")
    p.add_argument("--res", type=float, default=5.0)
    p.add_argument("--shape", type=int, nargs=2, default=(700, 508), metavar=("H", "W"))
    p.add_argument("--q", type=float, default=0.10)
    p.add_argument("--smooth", type=float, default=1.2, help="regional-slope smoothing (cells)")
    args = p.parse_args()
    X0, Y0 = args.bounds
    res = args.res
    H, W = args.shape

    # dense-epoch horizontal low-q is both the reference surface and the current-method 2021 side
    z21 = horizontal_low(args.after_laz, X0, Y0, res, H, W, args.q)
    g08_h = horizontal_low(args.before_laz, X0, Y0, res, H, W, args.q)

    plane = regional_plane(z21, res, args.smooth)
    g08_n = slope_normal_low(args.before_laz, plane, res, X0, Y0, H, W, args.q)
    g21_n = slope_normal_low(args.after_laz, plane, res, X0, Y0, H, W, args.q)

    dod_h = z21 - g08_h                         # current
    dod_n = g21_n - g08_n                        # slope-normal (common plane cancels)
    m = np.isfinite(dod_h) & np.isfinite(dod_n)

    bh, bn = _rstd(_bandpass(dod_h, m)[m]), _rstd(_bandpass(dod_n, m)[m])
    print(f"band-pass rstd:  horizontal {bh:.4f} m  ->  slope-normal {bn:.4f} m  ({100*(bn/bh-1):+.0f}%)")
    print(f"total DoD rstd:  horizontal {_rstd(dod_h[m]):.4f} m  ->  slope-normal {_rstd(dod_n[m]):.4f} m")


if __name__ == "__main__":
    main()

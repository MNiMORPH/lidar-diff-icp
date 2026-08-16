#!/usr/bin/env python3
"""Test a DeLong et al. (2022) vertical correction surface against the SE
undulation.

Build gridded ground surfaces (low-percentile per cell) for both epochs, form
the raw vertical difference on ground, and fit a stable-area IDW correction
surface (coreg.correction_surface). If that surface -- built ONLY from stable,
low-slope, small-|dz| ground -- reproduces the SE ~500 m undulation, the
undulation is a spurious vertical warp (it appears on ground that should not
have changed) and subtracting the surface flattens it. If the surface is flat
there, the undulation is not captured by stable ground (more likely real).

No M3C2 rerun: the correction surface is vertical and the SE undulation is not
slope-correlated (R^2=0), so the gridded vertical difference carries it.
"""
import argparse
from pathlib import Path
import numpy as np
import laspy
import pandas as pd

from lidar_diff_icp import io, coreg


def ground_grid(x, y, z, q, res, X0, Y0, nx, ny):
    ix = ((x - X0) / res).astype(int); iy = ((y - Y0) / res).astype(int)
    ok = (ix >= 0) & (ix < nx) & (iy >= 0) & (iy < ny)
    s = pd.Series(z[ok]).groupby(iy[ok] * nx + ix[ok]).quantile(q)
    out = np.full(nx * ny, np.nan); out[s.index.values] = s.values
    return out.reshape(ny, nx)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("before_laz"); ap.add_argument("after_last_laz")
    ap.add_argument("--bounds", nargs=4, type=float, required=True)
    ap.add_argument("--res", type=float, default=5.0)
    ap.add_argument("--figdir", default="figures")
    a = ap.parse_args()
    X0, Y0, X1, Y1 = a.bounds
    res = a.res; nx = int(round((X1 - X0) / res)); ny = int(round((Y1 - Y0) / res))

    # 2008 last return, internal align (translation), ground grid (low 10%)
    f = laspy.read(a.before_laz)
    x8 = np.asarray(f.x); y8 = np.asarray(f.y); z8 = np.asarray(f.z)
    ps8 = np.asarray(f.point_source_id); cl8 = np.asarray(f.classification)
    rn8 = np.asarray(f.return_number); nr8 = np.asarray(f.number_of_returns)
    be8 = rn8 == nr8
    pc = io.PointCloud(x8, y8, z8, ps8, cl8, np.zeros_like(z8), np.zeros_like(ps8),
                       io.MN_2008_CRS)
    corr, _, _ = coreg.align_swaths(pc, ref=int(ps8.min()))
    xc, yc, zc = x8.copy(), y8.copy(), z8.copy()
    for s, (dx, dy, dz) in corr.items():
        m = ps8 == s; xc[m] += dx; yc[m] += dy; zc[m] += dz
    Z08 = ground_grid(xc[be8], yc[be8], zc[be8], 0.10, res, X0, Y0, nx, ny)

    g = laspy.read(a.after_last_laz)
    x2 = np.asarray(g.x); y2 = np.asarray(g.y); z2 = np.asarray(g.z)
    m2 = (x2 >= X0) & (x2 < X1) & (y2 >= Y0) & (y2 < Y1)
    Z21 = ground_grid(x2[m2], y2[m2], z2[m2], 0.10, res, X0, Y0, nx, ny)

    cs = coreg.correction_surface(Z21, Z08, res, X0, Y0,
                                  slope_thresh_deg=3.0, dz_thresh=0.7, radius=400.0)
    dz = cs["dz"]; C = cs["C"]; corrected = dz - C
    print(f"stable ground cells: {cs['n_stable']:,} ({100*cs['frac_stable']:.1f}% of grid)")
    print(f"raw dz (2021-2008) NMAD {coreg._nmad(dz[np.isfinite(dz)]):.3f} m")
    print(f"correction surface C: range [{np.nanmin(C):+.3f}, {np.nanmax(C):+.3f}] m")
    print(f"corrected dz NMAD    {coreg._nmad(corrected[np.isfinite(corrected)]):.3f} m")

    # SE transverse profiles (cross axis perpendicular to 30 deg trend)
    XX, YY = np.meshgrid(X0 + (np.arange(nx) + 0.5) * res, Y0 + (np.arange(ny) + 0.5) * res)
    SE = (XX >= 578900) & (YY <= 4884400)
    w = (XX - 579500) * np.cos(np.radians(30)) - (YY - 4883500) * np.sin(np.radians(30))
    bins = np.arange(-400, 600, 20)

    def prof(A):
        m = SE & np.isfinite(A); idx = np.digitize(w[m], bins)
        s = pd.Series(A[m]).groupby(idx).median(); c = pd.Series(A[m]).groupby(idx).size()
        return s, c
    pdz, cnt = prof(dz); pC, _ = prof(C); pcorr, _ = prof(corrected)

    def amp(s, c):
        v = s[c[c > 100].index]; return float(v.max() - v.min())
    print(f"\nSE transverse undulation amplitude (max-min of median):")
    print(f"  raw dz          : {amp(pdz, cnt):.3f} m")
    print(f"  correction C    : {amp(pC, cnt):.3f} m")
    print(f"  corrected dz-C  : {amp(pcorr, cnt):.3f} m")

    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    Path(a.figdir).mkdir(exist_ok=True)
    ext = (X0, X1, Y0, Y1); v = float(np.nanpercentile(np.abs(dz), 98))
    fig, ax = plt.subplots(2, 3, figsize=(18, 11))
    for A, axi, ttl in [(dz, ax[0, 0], "raw dz (2021-2008 ground)"),
                        (C, ax[0, 1], "correction surface C (stable-area IDW)"),
                        (corrected, ax[0, 2], "corrected: dz - C")]:
        im = axi.imshow(A, extent=ext, origin="lower", cmap="RdBu", vmin=-v, vmax=v)
        axi.set_title(ttl); fig.colorbar(im, ax=axi, shrink=0.6)
        axi.axhline(4884400, color="k", lw=0.5); axi.axvline(578900, color="k", lw=0.5)
    st = np.where(cs["stable"], 1.0, np.nan)
    ax[1, 0].imshow(st, extent=ext, origin="lower", cmap="Greys", vmin=0, vmax=1.5)
    ax[1, 0].set_title(f"stable mask ({100*cs['frac_stable']:.0f}%)")
    bx = bins[np.clip(pdz.index[cnt > 100], 0, len(bins) - 1)]
    ax[1, 1].axhline(0, color="k", lw=0.6)
    ax[1, 1].plot(bx, pdz[cnt > 100].values, "-o", ms=3, label="raw dz", color="firebrick")
    ax[1, 1].plot(bx, pC[cnt > 100].values, "-s", ms=3, label="correction C", color="navy")
    ax[1, 1].plot(bx, pcorr[cnt > 100].values, "-^", ms=3, label="corrected dz-C", color="green")
    ax[1, 1].set_xlabel("cross-feature distance (m)"); ax[1, 1].set_ylabel("change (m)")
    ax[1, 1].legend(); ax[1, 1].set_title("SE transverse profiles")
    ax[1, 2].axis("off")
    fig.savefig(f"{a.figdir}/correction_surface_test.png", dpi=110, bbox_inches="tight")
    print(f"wrote {a.figdir}/correction_surface_test.png")


if __name__ == "__main__":
    main()

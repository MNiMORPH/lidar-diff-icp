#!/usr/bin/env python3
"""Local curvature (d2z/dx2, d2z/dy2) at every ridgecrest point + hillslope-diffusion K on
agricultural crests.

(1) Grid-aligned curvature by a windowed 2nd-order polynomial fit (2-D Savitzky-Golay):
    fit z = a + b x + c y + d x^2 + e y^2 + f xy over a +/-L window (x,y in metres from the
    cell centre); d2z/dx2 = 2d, d2z/dy2 = 2e. Convex-up (hilltop) -> NEGATIVE. Computed over
    the ENTIRE map, saved, and sampled at every ridgecrest pixel.

(2) On AGRICULTURAL (open) ridgecrest cells -- where there is no incoming material AND the
    gen2 canopy artifact ~0, so DoD is real erosion -- solve the hillslope-diffusion eq
        dz/dt = K d2z/dx2   ->   K = slope of (dz/dt) vs (d2z/dx2).
    dz/dt = DoD / dt, dt from the flight dates. Reported through-origin (as written),
    with-intercept (isolates any datum offset), and via the full Laplacian (physical form).

    env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python analysis/ridgelines/curvature_diffusion.py
"""
import numpy as np
from datetime import date
from scipy.ndimage import correlate, distance_transform_edt

import sys, os, json
TILE = sys.argv[1] if len(sys.argv) > 1 else "elba_fulldensity"
D = f"data/derived/{TILE}"
def _res(tile):
    for fn in ("meta.json", "corrections_geoid.json", "corrections.json"):
        p = f"data/derived/{tile}/{fn}"
        if os.path.exists(p): j = json.load(open(p)); return float(j.get("res") or j.get("res_m"))
    return 5.0
RES = _res(TILE)
z = np.load(f"{D}/z_after.npy")
ny, nx = z.shape
zf = z.copy(); nm = ~np.isfinite(zf)
if nm.any():
    zf = zf[tuple(distance_transform_edt(nm, return_distances=False, return_indices=True))]

# --- (1) windowed-quadratic curvature kernels (x = east, y = north) ----------------------
L = 3                                                  # +/-15 m window (7x7)
off = RES * np.arange(-L, L + 1)
Xg, Yg = np.meshgrid(off, off)                          # metres from centre
x = Xg.ravel(); y = Yg.ravel()
G = np.column_stack([np.ones_like(x), x, y, x*x, y*y, x*y])
pinv = np.linalg.pinv(G)                                # (6, N)
kxx = 2 * pinv[3].reshape(2*L+1, 2*L+1)                 # d2z/dx2 kernel
kyy = 2 * pinv[4].reshape(2*L+1, 2*L+1)                 # d2z/dy2 kernel
zxx = correlate(zf, kxx, mode="nearest")                # 1/m, convex-up NEGATIVE
zyy = correlate(zf, kyy, mode="nearest")
lap = zxx + zyy                                         # Laplacian (orientation-free)
for nm_, arr in [("curv_xx", zxx), ("curv_yy", zyy), ("curv_laplacian", lap)]:
    np.save(f"{D}/{nm_}.npy", arr)
print(f"saved {D}/curv_xx/curv_yy/curv_laplacian (L=+/-{L*RES:.0f} m)")

# elba-only: add curvature to the ridgecrest table + report crest diffusion (needs crest/dod/pen)
if os.path.exists(f"{D}/ridgecrest_pixels.npz") and os.path.exists(f"{D}/crest_mask.npy"):
    dod = np.load("data/derived/elba_refdatum/dod_geoid.npy"); crest = np.load(f"{D}/crest_mask.npy")
    R = dict(np.load(f"{D}/ridgecrest_pixels.npz", allow_pickle=True))
    rr = R["row"].astype(int); cci = R["col"].astype(int)
    R["curv_xx"] = zxx[rr, cci]; R["curv_yy"] = zyy[rr, cci]; R["curv_laplacian"] = lap[rr, cci]
    np.savez(f"{D}/ridgecrest_pixels.npz", **R)
    cmask = crest & np.isfinite(dod)
    print(f"crest curvature (1/m): d2z/dx2 med={np.median(zxx[cmask]):+.4f}  "
          f"d2z/dy2 med={np.median(zyy[cmask]):+.4f}  Laplacian med={np.median(lap[cmask]):+.4f}")

# --- (2) diffusion K on agricultural (open) crests (elba only; needs crest/pen/dod) ------
if not (os.path.exists(f"{D}/crest_mask.npy") and os.path.exists(f"{D}/penetration.npy")):
    import sys; sys.exit(0)
crest = np.load(f"{D}/crest_mask.npy"); pen = np.load(f"{D}/penetration.npy")
dod = np.load("data/derived/elba_refdatum/dod_geoid.npy")
dt_yr = (date(2021, 5, 1) - date(2008, 11, 21)).days / 365.25   # flight-date span
dzdt = dod / dt_yr                                       # m/yr per cell
ag = crest & (pen >= 0.45) & np.isfinite(dod)
xx = zxx[ag]; rate = dzdt[ag]; lp = lap[ag]
print(f"\n=== diffusion on AGRICULTURAL crests (open, pen>=0.45): dz/dt = K d2z/dx2 ===")
print(f"  dt = {dt_yr:.2f} yr (2008-11-21 -> 2021-05-01)   n = {ag.sum()} ag-crest cells")
print(f"  median dz/dt = {np.median(rate)*1000:+.2f} mm/yr   median d2z/dx2 = {np.median(xx):+.4f} 1/m")

K_orig = np.sum(rate*xx) / np.sum(xx*xx)                 # least-squares through origin
A = np.column_stack([xx, np.ones_like(xx)])
(Kc, c0), *_ = np.linalg.lstsq(A, rate, rcond=None)     # with intercept
K_lap = np.sum(rate*lp) / np.sum(lp*lp)                  # physical Laplacian form
print(f"\n  K (dz/dt = K d2z/dx2, through origin) = {K_orig:.4f} m^2/yr")
print(f"  K (with intercept)                    = {Kc:.4f} m^2/yr   intercept = {c0*1000:+.2f} mm/yr")
print(f"  K (dz/dt = K*Laplacian, physical)     = {K_lap:.4f} m^2/yr")
print(f"\n  (literature: soil creep ~1e-3..1e-2 m^2/yr; TILLAGE diffusion ~1e-2..1e-1 -- ag "
      f"land expected high.)")

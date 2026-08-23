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

RES = 5.0; X0, Y0 = 577492.8, 4882737.6
z = np.load("data/derived/elba_fulldensity/z_after.npy")
dod = np.load("data/derived/elba_refdatum/dod_geoid.npy")
pen = np.load("data/derived/elba_fulldensity/penetration.npy")
crest = np.load("data/derived/elba_fulldensity/crest_mask.npy")
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
    np.save(f"data/derived/elba_fulldensity/{nm_}.npy", arr)

# add curvature to the per-pixel ridgecrest table (re-save)
R = dict(np.load("data/derived/elba_fulldensity/ridgecrest_pixels.npz", allow_pickle=True))
rr = R["row"].astype(int); cci = R["col"].astype(int)
R["curv_xx"] = zxx[rr, cci]; R["curv_yy"] = zyy[rr, cci]; R["curv_laplacian"] = lap[rr, cci]
np.savez("data/derived/elba_fulldensity/ridgecrest_pixels.npz", **R)
print(f"saved curv_xx/curv_yy/curv_laplacian (full map + per-crest-pixel table, L=+/-{L*RES:.0f} m)")
cmask = crest & np.isfinite(dod)
print(f"crest curvature (1/m): d2z/dx2 med={np.median(zxx[cmask]):+.4f}  "
      f"d2z/dy2 med={np.median(zyy[cmask]):+.4f}  Laplacian med={np.median(lap[cmask]):+.4f}")

# --- (2) diffusion K on agricultural (open) crests --------------------------------------
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

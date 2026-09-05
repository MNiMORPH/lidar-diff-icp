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

Part (1) needs only z_after and is written for any tile. Part (2) needs a DoD, a crest mask,
a cover layer and the two flight dates; every one of those is per-tile, so none is defaulted
to Elba's value -- a wrong dt or a wrong DoD would rescale K silently.

NOTE ON ORDER: part (1) AUGMENTS ridgecrest_pixels.npz in place with curv_xx / curv_yy /
curv_laplacian. convexity_dod_landcover.py rewrites that file from scratch, so run this
AFTER it, or those three columns are dropped.

    env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python analysis/ridgelines/curvature_diffusion.py \
        --tile elba_fulldensity --dod data/derived/elba_refdatum/dod_geoid.npy \
        --gen1-date 2008-11-21 --gen2-date 2021-05-01
"""
import numpy as np
from datetime import date
from scipy.ndimage import correlate, distance_transform_edt

import argparse, sys, os, json
_ap = argparse.ArgumentParser()
_ap.add_argument("--tile", default="elba_fulldensity",
                 help="tile name under data/derived/; resolution is read from its own meta")
_ap.add_argument("--dod", default=None,
                 help="DoD grid for parts (1)-(2) (m, + = elevation rose). Not defaulted: "
                      "Elba's shipped run reads data/derived/elba_refdatum/dod_geoid.npy, "
                      "outside the tile directory. Omit to run part (1) only.")
_ap.add_argument("--gen1-date", default=None, help="gen1 flight date, ISO (Elba: 2008-11-21)")
_ap.add_argument("--gen2-date", default=None, help="gen2 flight date, ISO (Elba: 2021-05-01)")
ARGS = _ap.parse_args()
TILE = ARGS.tile
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

# add curvature to the ridgecrest table + report crest curvature (needs the crest table)
if os.path.exists(f"{D}/ridgecrest_pixels.npz") and os.path.exists(f"{D}/crest_mask.npy"):
    crest = np.load(f"{D}/crest_mask.npy")
    R = dict(np.load(f"{D}/ridgecrest_pixels.npz", allow_pickle=True))
    rr = R["row"].astype(int); cci = R["col"].astype(int)
    R["curv_xx"] = zxx[rr, cci]; R["curv_yy"] = zyy[rr, cci]; R["curv_laplacian"] = lap[rr, cci]
    np.savez(f"{D}/ridgecrest_pixels.npz", **R)
    if ARGS.dod is None:
        cmask = crest
        print("  no --dod given: crest curvature is summarised over ALL crest cells, not "
              "over the DoD-finite subset the Elba run reported")
    else:
        dod = np.load(ARGS.dod)
        if dod.shape != z.shape:
            raise SystemExit(f"--dod {ARGS.dod} is {dod.shape} but {TILE}'s grid is "
                             f"{z.shape}; these are different grids")
        cmask = crest & np.isfinite(dod)
    print(f"crest curvature (1/m): d2z/dx2 med={np.median(zxx[cmask]):+.4f}  "
          f"d2z/dy2 med={np.median(zyy[cmask]):+.4f}  Laplacian med={np.median(lap[cmask]):+.4f}")
else:
    print(f"  no ridgecrest table or crest mask under {D}: curvature columns NOT added "
          f"(run convexity_dod_landcover.py --tile {TILE} first)")

# --- (2) diffusion K on agricultural (open) crests (needs crest + cover + DoD + dates) ---
# The OPEN class comes from open_pfs.npy -- the repo's own PyForestScan mask, written by
# the declared pfs_cover step at a stated cut (cover <= 0.1, forest_metrics_pfs.py:201),
# identical on every tile. It replaces `penetration >= 0.45`, retired 2026-09-05: penetration
# was geometry-confounded (ground-return fraction correlates -0.84 with scan angle), so it
# never measured cover in the first place. No threshold is invented here -- if open_pfs.npy
# is absent, Part 2 does not run.
_missing = [n for n in ("crest_mask.npy", "open_pfs.npy") if not os.path.exists(f"{D}/{n}")]
if ARGS.dod is None: _missing.append("--dod")
if not ARGS.gen1_date or not ARGS.gen2_date: _missing.append("--gen1-date/--gen2-date")
if _missing:
    print(f"\n=== PART 2 SKIPPED: {', '.join(_missing)} not available for {TILE} ===")
    print("  K is NOT reported for this tile; it is absent, not zero. The flight dates are "
          "per-acquisition and are never defaulted -- a wrong dt rescales K directly.")
    sys.exit(0)
crest = np.load(f"{D}/crest_mask.npy"); openg = np.load(f"{D}/open_pfs.npy").astype(bool)
dod = np.load(ARGS.dod)
_d1 = date.fromisoformat(ARGS.gen1_date); _d2 = date.fromisoformat(ARGS.gen2_date)
dt_yr = (_d2 - _d1).days / 365.25                               # flight-date span
dzdt = dod / dt_yr                                       # m/yr per cell
ag = crest & openg & np.isfinite(dod)
xx = zxx[ag]; rate = dzdt[ag]; lp = lap[ag]
print(f"\n=== diffusion on AGRICULTURAL crests (open_pfs, cover<=0.1): dz/dt = K d2z/dx2 ===")
print(f"  dt = {dt_yr:.2f} yr ({_d1} -> {_d2})   n = {ag.sum()} ag-crest cells")
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

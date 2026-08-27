#!/usr/bin/env python3
"""The NGV-corrected DoD judged on STEADY-STATE cells, not on convex ones.

Uses this package's own `steady_state_cells` module rather than a new definition:
|grad^2 z| < eps_curv (locally planar -> zero flux divergence), slope < max_slope
(not mass-wasting), restricted to the drainage-divide network and outside the
floodplain. On such cells linear diffusion predicts NO elevation change, so the
DoD there is expected to centre on zero.

This matters because the pipeline's `stable` mask is flat-OR-CONVEX. Convex divides
are precisely where dz/dt = K grad^2 z is negative, so their DoD is expected to be
somewhat negative and cannot be read as a datum check.

eps_curv is chosen by the module's own quantile helper (central band of the
curvature distribution over the base population), not by a number picked here.

    ./lidar-icp/bin/python analysis/steady_state/dod_ngv_on_steady_cells.py
"""
import argparse, os, sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from steady_state_cells import steady_state_mask, eps_curv_from_quantile, nmad

ap = argparse.ArgumentParser(description=__doc__,
                             formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--tile", default="data/derived/elba_fulldensity")
ap.add_argument("--max-slope", type=float, default=15.0)
ap.add_argument("--save", action="store_true", help="write the mask at the default cuts")
A = ap.parse_args()
T = A.tile

curv = np.load(f"{T}/curv_laplacian.npy")
slope = np.load(f"{T}/slope.npy")
dod = np.load(f"{T}/dod.npy")
corr = np.load(f"{T}/dod_ngv.npy")
ngv = np.load(f"{T}/ngv.npy")
flood = np.load(f"{T}/floodplain_mask.npy").astype(bool)
divides = np.isfinite(np.load(f"{T}/kappa_L20.npy"))     # the S&S divide network

cover = divides & ~flood
base = cover & (slope < A.max_slope) & np.isfinite(curv) & np.isfinite(slope) & np.isfinite(dod)
print(f"divide network {int(divides.sum()):,} cells; outside floodplain {int(cover.sum()):,}; "
      f"base (slope < {A.max_slope:g} deg, finite) {int(base.sum()):,}")
print(f"curvature over base: median {np.median(curv[base]):+.5f}  "
      f"NMAD {nmad(curv[base]):.5f} 1/m\n")

print("DoD on locally PLANAR divide cells -- where diffusion predicts zero change.")
print("eps_curv from the module's central-band quantile helper.")
print(f"  {'central':>8} {'eps_curv':>10} {'n':>7} {'med before':>11} {'NMAD before':>12} "
      f"{'med after':>11} {'NMAD after':>12} {'NGV':>6}")
for cf in (0.10, 0.20, 0.30, 0.50):
    eps = eps_curv_from_quantile(curv, base, central_frac=cf)
    m = steady_state_mask(curv, slope, dod, cover_mask=cover,
                          eps_curv=eps, max_slope=A.max_slope) & np.isfinite(corr)
    if m.sum() < 30:
        print(f"  {cf:8.2f} {eps:10.5f} {int(m.sum()):7d}   -- too few")
        continue
    print(f"  {cf:8.2f} {eps:10.5f} {int(m.sum()):7d} {1000*np.median(dod[m]):+11.1f} "
          f"{1000*nmad(dod[m]):12.1f} {1000*np.median(corr[m]):+11.1f} "
          f"{1000*nmad(corr[m]):12.1f} {np.nanmedian(ngv[m]):6.3f}")

eps30 = eps_curv_from_quantile(curv, base, central_frac=0.30)
m30 = steady_state_mask(curv, slope, dod, cover_mask=cover, eps_curv=eps30,
                        max_slope=A.max_slope) & np.isfinite(corr)
print(f"\nFor contrast, the CONVEX tail of the same divides "
      f"(curv < -eps_curv, central_frac 0.30), where diffusion predicts LOWERING:")
mc = (cover & (slope < A.max_slope) & (curv < -eps30) & np.isfinite(dod) & np.isfinite(corr))
print(f"  {'convex divides':22s} n {int(mc.sum()):7d} "
      f"med before {1000*np.median(dod[mc]):+7.1f}  after {1000*np.median(corr[mc]):+7.1f} mm")
print(f"  {'planar divides':22s} n {int(m30.sum()):7d} "
      f"med before {1000*np.median(dod[m30]):+7.1f}  after {1000*np.median(corr[m30]):+7.1f} mm")

if A.save:
    np.save(f"{T}/steady_state_divides.npy", m30)
    print(f"\nwrote {T}/steady_state_divides.npy "
          f"({int(m30.sum()):,} cells, eps_curv {eps30:.5f}, slope < {A.max_slope:g} deg)")

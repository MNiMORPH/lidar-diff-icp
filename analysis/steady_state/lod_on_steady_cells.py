#!/usr/bin/env python3
"""Refit the heteroscedastic LoD on STEADY-STATE planar divide cells instead of the
pipeline's flat-OR-CONVEX stable mask, and compare.

WHAT THIS IS AND IS NOT. The pipeline's LoD uses THREE covariates: slope, |curv|, and
`stderr` = sqrt(sum_epoch roughness^2 / n), built from per-epoch roughness and ground
density grids that difference_dem does NOT save. So this cannot reproduce lod.npy, and
does not try. Both fits here use the SAME two covariates (slope, |curv|) and differ ONLY
in the calibration mask, which is the question. Absolute values are not comparable with
lod.npy; the mask-to-mask difference is.

The concern being tested: the pipeline mask is 78% convex, and convex divides are where
dz/dt = K grad^2 z is negative, so part of the "scatter" it calibrates on is real
diffusive lowering -- which would INFLATE the LoD and make the detector under-sensitive.

    env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python \
        analysis/steady_state/lod_on_steady_cells.py
"""
import argparse, json, os, sys
import numpy as np
from scipy.ndimage import gaussian_filter, distance_transform_edt as edt

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from stable_mask_repro import stable_mask, clip_stable, nmad
from lidar_diff_icp.pipeline import heteroscedastic_lod
from lidar_diff_icp.detect import detect_change_standard

ap = argparse.ArgumentParser(description=__doc__,
                             formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--tile", default="data/derived/elba_fulldensity")
ap.add_argument("--valley-top", dest="valley_top", default="registry",
                help="valley top for the stable mask: an elevation in metres, "
                     "'registry', or 'histogram'. Never chosen for you.")
A = ap.parse_args()
T = A.tile
cfg = json.load(open(f"{T}/corrections.json")); res = float(cfg["res_m"])
dod = np.load(f"{T}/dod.npy")
Z21 = np.load(f"{T}/z_after.npy"); slope = np.load(f"{T}/slope.npy")
lod_pipe = np.load(f"{T}/lod.npy")
steady = np.load(f"{T}/steady_state_divides.npy")

Zf = Z21.copy(); nanm = np.isnan(Zf)
if nanm.any():
    Zf = Zf[tuple(edt(nanm, return_distances=False, return_indices=True))]
Zs = gaussian_filter(Zf, 1.0)                                  # pipeline.py:766
abs_curv = np.abs(np.gradient(np.gradient(Zs, res, axis=0), res, axis=0)
                  + np.gradient(np.gradient(Zs, res, axis=1), res, axis=1))
sdeg = slope

st_pipe, _ = clip_stable(stable_mask(Z21, res, valley_top_m=A.valley_top,
                                     tile_dir=T), dod)
fin = np.isfinite(dod) & np.isfinite(lod_pipe)

print(f"calibration masks:  pipeline {int(st_pipe.sum()):,} cells   "
      f"steady-state divides {int(steady.sum()):,} cells")
print(f"  convex fraction:  pipeline {100*np.mean(np.load(f'{T}/curv_laplacian.npy')[st_pipe]<0):.1f}%"
      f"   steady {100*np.mean(np.load(f'{T}/curv_laplacian.npy')[steady]<0):.1f}%")
print(f"  stable-ground NMAD of the DoD: pipeline {1000*nmad(dod[st_pipe]):.1f} mm   "
      f"steady {1000*nmad(dod[steady]):.1f} mm\n")

out = {}
for nm, mask in (("pipeline mask", st_pipe), ("steady-state", steady)):
    L = heteroscedastic_lod(dod, sdeg, abs_curv, mask & np.isfinite(dod))
    if L is None:
        print(f"  {nm}: model returned None (xdem unavailable or fit degenerate)")
        continue
    out[nm] = L
    print(f"  {nm:14s} LoD (mm)  median {1000*np.median(L[fin]):7.1f}  "
          f"p90 {1000*np.percentile(L[fin],90):7.1f}  max {1000*np.nanmax(L[fin]):7.1f}")
print(f"  {'lod.npy (3-cov, not comparable)':14s}  median "
      f"{1000*np.median(lod_pipe[fin]):7.1f}  p90 {1000*np.percentile(lod_pipe[fin],90):7.1f}")

if len(out) == 2:
    a_, b_ = out["pipeline mask"], out["steady-state"]
    d = (b_ - a_)[fin]
    print(f"\n  steady-state minus pipeline-mask LoD: median {1000*np.median(d):+.1f} mm "
          f"({100*np.median(d)/np.median(a_[fin]):+.1f}%)")

    print(f"\nDETECTION with each LoD (same DoD, same detector, same stable set for tau_sys)")
    print(f"  {'LoD from':16s} {'DoD':14s} {'regions':>8} {'cells':>9} {'%':>7} {'net m3':>13}")
    for lnm, L in out.items():
        for dnm, d_, s_ in (("dod.npy", dod, st_pipe),):
            det = detect_change_standard(d_, L, s_, res)
            ch = det["change"]
            print(f"  {lnm:16s} {dnm:14s} {len(det['regions']):8d} {int(ch.sum()):9,d} "
                  f"{100*ch.sum()/fin.sum():6.2f}% "
                  f"{sum(r['volume_m3'] for r in det['regions']):+13,.0f}")
    np.save(f"{T}/lod_steady.npy", out["steady-state"])
    print(f"\nwrote {T}/lod_steady.npy")

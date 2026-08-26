#!/usr/bin/env python3
"""Recalibrate the level of detection for the cover-corrected DoD.

`lod.npy` was fitted on the UNcorrected product, so comparing detection rates against it
is unfair to the corrected one. This refits `pipeline.heteroscedastic_lod` on the corrected
DoD, and -- for a like-for-like comparison -- refits it on the uncorrected DoD too, using
the SAME stable mask and the SAME covariates in both cases.

Covariates: slope and |curvature| (the pipeline's fallback pair; the stderr covariate needs
per-epoch roughness and counts, which are not saved to disk). Stable mask: the offset-
independent criteria of `refcells` WITHOUT the ridge and curvature cuts -- those cap
|curv| at 0.015, which collapses the curvature axis and makes xdem's Delaunay fit fail --
then the pipeline's iterative 3-NMAD clip of the DoD being fitted.

    env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python analysis/ridgelines/lod_cover_q2.py
"""
import argparse
import numpy as np
from scipy.ndimage import gaussian_filter

from lidar_diff_icp.pipeline import heteroscedastic_lod
from lidar_diff_icp.refcells import reference_cells

ap = argparse.ArgumentParser()
ap.add_argument("--tile", default="data/derived/elba_fulldensity")
A = ap.parse_args()
D = A.tile

sdeg = np.load(f"{D}/slope.npy")
Zf = np.load(f"{D}/z_after.npy")
res = 5.0
_z = np.where(np.isfinite(Zf), Zf, np.nanmedian(Zf))
abs_curv = np.abs(np.gradient(np.gradient(gaussian_filter(_z, 1.0), res, axis=0), res, axis=0)
                  + np.gradient(np.gradient(gaussian_filter(_z, 1.0), res, axis=1), res, axis=1))
stable0, _ = reference_cells(D, slope_max=90.0, gross_change_mm=500.0,
                            curv_max=np.inf, require_ridge=False)
stable0 = stable0.reshape(Zf.shape)


def clip_stable(dod, base):
    """The pipeline's iterative 3-NMAD sigma-clip of the reporting stable mask."""
    s = base & np.isfinite(dod)
    for _ in range(8):
        v = dod[s]; med = np.median(v)
        nm = 1.4826 * np.median(np.abs(v - med))
        keep = s & (np.abs(dod - med) < 3.0 * max(nm, 1e-3))
        if keep.sum() == s.sum():
            break
        s = keep
    return s


print(f"{'product':22s} {'cells':>8s} {'stable':>7s} {'sigma':>7s} {'LoD med':>8s} "
      f"{'LoD p90':>8s} {'detected':>9s}")
out = {}
for nm, fn in (("dod.npy (uncorrected)", "dod.npy"),
               ("dod_cover_q2.npy", "dod_cover_q2.npy")):
    dod = np.load(f"{D}/{fn}")
    st = clip_stable(dod, stable0)
    sig = 1.4826 * np.median(np.abs(dod[st] - np.median(dod[st])))
    lod = heteroscedastic_lod(dod, sdeg, abs_curv, st)
    if lod is None:
        print(f"  {nm:22s}  xdem unavailable or fit failed")
        continue
    ok = np.isfinite(dod) & np.isfinite(lod)
    det = 100 * np.mean(np.abs(dod[ok]) > lod[ok])
    print(f"{nm:22s} {ok.sum():8,d} {st.sum():7,d} {1000*sig:7.1f} {1000*np.median(lod[ok]):8.1f} "
          f"{1000*np.percentile(lod[ok],90):8.1f} {det:8.1f}%")
    out[fn] = lod
np.save(f"{D}/lod_cover_q2.npy", out["dod_cover_q2.npy"])
np.save(f"{D}/lod_refit_uncorrected.npy", out["dod.npy"])
print(f"\nwrote {D}/lod_cover_q2.npy and {D}/lod_refit_uncorrected.npy  (sigma, LoD in mm)")

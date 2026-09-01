#!/usr/bin/env python3
"""Reproduce the pipeline's `stable` mask outside difference_dem, and the LoD it feeds.

`difference_dem` does not save the mask it calibrates the LoD on, so any analysis that
wants to re-run the detector has to rebuild it. This does, with the pipeline's OWN lines,
in BOTH stages -- the geometric mask (pipeline.py 517-526) and the iterative 3-NMAD
sigma-clip against the DoD (pipeline.py 751-760). It is the CLIPPED mask that
difference_dem returns as "stable" and that run_all_sites.py hands the detector.

Verified against the values the tile recorded when it was built:
stable_clip_fraction 0.049444 vs 0.049500 recorded, stable_1sigma_m 0.051572 vs 0.051600.

The clip is DoD-dependent, so it must be recomputed for any alternative DoD rather than
reused -- two DoDs do not share a stable set.
"""
import argparse, json, os
import numpy as np
from scipy.ndimage import uniform_filter, gaussian_filter, distance_transform_edt as edt
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from lidar_diff_icp import coreg
from lidar_diff_icp.detect import detect_change_standard
from lidar_diff_icp.viz import hillshade

Z95 = 1.96


def stable_mask(Z21, res):
    """Verbatim from pipeline.py: terrain masks from the reference ground."""
    Zf = Z21.copy(); nanm = np.isnan(Zf)
    if nanm.any():
        Zf = Zf[tuple(edt(nanm, return_distances=False, return_indices=True))]
    tpi = Z21 - uniform_filter(Zf, size=int(2 * 300 / res), mode="nearest")
    sdeg = np.degrees(coreg.slope_aspect(gaussian_filter(Zf, 2.0), res)[0])
    Zsm = gaussian_filter(Zf, 50 / res / 2)
    lap = (np.gradient(np.gradient(Zsm, res, axis=0), res, axis=0)
           + np.gradient(np.gradient(Zsm, res, axis=1), res, axis=1))
    convex = (sdeg > 5) & (sdeg < 35) & (tpi > -2) & (lap < 0)
    return ((sdeg < 3) & (tpi > -2)) | convex


def clip_stable(stable, dod):
    """pipeline.py 751-760: iterative 3-NMAD sigma-clip, so real change in a floodplain
    wider than the TPI window does not bleed into the stable-ground error."""
    s = stable & np.isfinite(dod)
    n0 = int(s.sum())
    for _ in range(8):
        v = dod[s]; med = np.median(v); nm = 1.4826 * np.median(np.abs(v - med))
        keep = s & (np.abs(dod - med) < 3.0 * max(nm, 1e-3))
        if keep.sum() == s.sum():
            break
        s = keep
    return s, (1.0 - s.sum() / n0) if n0 else 0.0


def nmad(v):
    v = v[np.isfinite(v)]
    return float(1.4826 * np.median(np.abs(v - np.median(v)))) if v.size else np.nan


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tile", default="data/derived/elba_fulldensity")
    a = ap.parse_args()
    T = a.tile
    cfg = json.load(open(os.path.join(T, "corrections.json")))
    res = float(cfg["res_m"])
    dod = np.load(f"{T}/dod.npy"); Z21 = np.load(f"{T}/z_after.npy")
    geom = stable_mask(Z21, res)
    st, cf = clip_stable(geom, dod)
    print(f"geometric {int(geom.sum()):,} cells -> after the 3-NMAD clip {int(st.sum()):,} "
          f"({100*cf:.1f}% removed)")
    print(f"  recorded when the tile was built: stable_clip_fraction "
          f"{cfg['stable_clip_fraction']:.4f}, reproduced {cf:.4f}")
    print(f"  stable_1sigma_m recorded {cfg['stable_1sigma_m']:.6f}, "
          f"reproduced {nmad(dod[st]):.6f}")
    print(f"  DoD on stable ground: median {1000*np.median(dod[st]):+.1f} mm, "
          f"NMAD {1000*nmad(dod[st]):.1f} mm")


if __name__ == "__main__":
    main()

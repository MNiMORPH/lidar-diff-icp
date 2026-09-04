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


def stable_mask(Z21, res, *, valley_top_m, tile_dir=None, curv_max=0.005):
    """The pipeline's stable mask. NO LONGER A COPY.

    This was a verbatim duplicate of a block inside difference_dem, and it is why the
    project had two stable masks that could silently disagree: when the valley cut moved
    from TPI to elevation on 2026-09-04, the pipeline's copy was fixed and this one was not.
    It now delegates to lidar_diff_icp.terrain.terrain_masks, which is the single
    definition, so the two cannot drift again.

    ``valley_top_m`` is required and is never chosen for you: an elevation in metres,
    ``"registry"``, or ``"histogram"``. See terrain.resolve_valley_top.
    """
    from lidar_diff_icp import terrain
    return terrain.terrain_masks(Z21, res, valley_top_m=valley_top_m, tile_dir=tile_dir,
                                 curv_max=curv_max, verbose=False)["stable"]


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

#!/usr/bin/env python3
"""GEN1 GROUND-SELECTION TEST (Task 4) — does the +slope bias come from gen1's ground
POINT selection (CSF cloth) vs using gen1's own vendor class-2?

Both epochs use the identical median estimator; the ONLY asymmetry in the pipeline is
stage-1 selection (gen1=CSF, gen2=class2). Swap gen1 to its own class-2 and compare the
forested-slope DoD curve to the CSF baseline. gen2 (class-2, full density) is the SAME
reference for both, so:  gen1_ground = z_after - dod  =>  gen1_class2 - gen1_CSF = dod_csf - dod_class2.
Andy's physics: CSF cloth bridges roughness -> reads HIGH -> expect gen1_CSF > gen1_class2
(so dod_class2 > dod_csf, i.e. swapping to class2 pushes the bias UP, not down).

    env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python analysis/slope_bias/gen1_class2_test.py
"""
import time, numpy as np
from scipy.ndimage import distance_transform_edt
from lidar_diff_icp.pipeline import difference_dem
from lidar_diff_icp.canopy import ground_penetration

BEFORE = "data/before/4342-29-64.laz"
AFTER  = "data/after/3dep2021_fulldensity.laz"
BOUNDS = (577492.8, 4882737.6, 580032.8, 4886237.6)
RES = 5.0
COMMON = dict(res=RES, ground="slope_normal", after_ground="class2", stream=True,
              robust_stable=True)

runs = {}
for tag, gs, cache in (("gen1=CSF",   "csf",    "data/csf_cache/elba.las"),
                       ("gen1=class2","class2",  None)):
    t = time.time()
    print(f"[{tag}] difference_dem (ground_source={gs}) ...", flush=True)
    kw = dict(COMMON, ground_source=gs)
    if cache: kw["csf_cache"] = cache
    runs[tag] = difference_dem(BEFORE, AFTER, BOUNDS, **kw)
    r = runs[tag]; ex = np.isfinite(r["dod"])
    print(f"[{tag}] {time.time()-t:.0f}s  medDoD={np.nanmedian(r['dod'][ex])*1000:+.1f} mm  "
          f"sigma={r['stable_sigma']:.3f}", flush=True)

base = runs["gen1=CSF"]; nx, ny = base["nx"], base["ny"]; bnds = base["bounds"]
Z = base["z_after"].copy(); nm = ~np.isfinite(Z)
if nm.any():
    Z = Z[tuple(distance_transform_edt(nm, return_distances=False, return_indices=True))]
gy, gx = np.gradient(Z, RES); slope = np.degrees(np.arctan(np.hypot(gx, gy)))
pen = ground_penetration(AFTER, bnds, RES, nx, ny); forest = pen < 0.25

BINS = [0, 5, 10, 15, 20, 25, 30, 90]
dcsf, dc2 = runs["gen1=CSF"]["dod"], runs["gen1=class2"]["dod"]
print(f"\n=== median DoD (mm) vs slope on FORESTED cells ===")
print(f"{'slope':>8} {'n':>6} | {'gen1=CSF':>9} {'gen1=class2':>11} | {'class2-CSF':>10}")
for lo, hi in zip(BINS[:-1], BINS[1:]):
    m = np.isfinite(dcsf) & np.isfinite(dc2) & forest & (slope >= lo) & (slope < hi)
    if not m.any(): print(f"{lo:>3}-{hi:<3}  0"); continue
    a, b = np.median(dcsf[m])*1000, np.median(dc2[m])*1000
    print(f"{lo:>3}-{hi:<3} {m.sum():>6} | {a:>+9.1f} {b:>+11.1f} | {b-a:>+10.1f}")

print("\ngen1_CSF - gen1_class2 ground (= dod_class2 - dod_csf) by slope, forested:")
for lo, hi in zip(BINS[:-1], BINS[1:]):
    m = np.isfinite(dcsf) & np.isfinite(dc2) & forest & (slope >= lo) & (slope < hi)
    if m.any():
        print(f"  {lo:>3}-{hi:<3}: gen1_CSF is {np.median((dc2-dcsf)[m])*1000:+.1f} mm "
              f"{'HIGHER' if np.median((dc2-dcsf)[m])>0 else 'LOWER'} than gen1_class2")
print("\nVERDICT: if swapping gen1 CSF->class2 leaves the +slope bias (or raises it), gen1 "
      "selection is NOT the source -> bias is gen2-side (leaf-on). If it COLLAPSES, gen1 "
      "CSF classification was making gen1 low on forested slopes.")

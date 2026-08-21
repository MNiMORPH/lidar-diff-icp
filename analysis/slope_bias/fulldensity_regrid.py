#!/usr/bin/env python3
"""FULL-DENSITY RE-GRID (Task 4 decisive test) — does differencing gen2 ground from the
FULL 3DEP return set (--max-depth 12, 5.78 ground pts/m2) instead of the decimated
default-depth fetch (0.30 pts/m2) resolve the false hillslope aggradation (+tan(slope)
gen2-high bias)?

Identical methods for both runs -- only the after-cloud changes:
  gen1 CSF ground, slope_normal gridding, tie + drift corrections, robust_stable.
gen1 CSF is cached (data/csf_cache/elba.las) so it runs ONCE and both runs reuse it.

Metric: median DoD (gen2 - gen1) in slope bins on FORESTED cells (where the false
aggradation lives) and on BEDROCK (stable control). If the full cloud flattens the
slope trend, the bias was the decimation.

    env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python analysis/slope_bias/fulldensity_regrid.py
"""
import json, time
import numpy as np
from scipy.ndimage import distance_transform_edt
from lidar_diff_icp.pipeline import difference_dem
from lidar_diff_icp.canopy import ground_penetration

BEFORE = "data/before/4342-29-64.laz"
BOUNDS = (577492.8, 4882737.6, 580032.8, 4886237.6)
RES = 5.0
AFTERS = {"decimated":   "data/after/3dep2021_fulltile.laz",
          "fulldensity": "data/after/3dep2021_fulldensity.laz"}
COMMON = dict(res=RES, ground="slope_normal", ground_source="csf",
              after_ground="class2", stream=True, robust_stable=True,
              csf_cache="data/csf_cache/elba.las")

runs = {}
for tag, after in AFTERS.items():
    t = time.time()
    print(f"[{tag}] difference_dem  after={after}", flush=True)
    runs[tag] = difference_dem(BEFORE, after, BOUNDS, **COMMON)
    r = runs[tag]; dod = r["dod"]; ex = np.isfinite(dod)
    print(f"[{tag}] done {time.time()-t:.0f}s  sigma={r['stable_sigma']:.3f}  "
          f"medLoD={np.nanmedian(r['lod']):.3f}  medDoD={np.nanmedian(dod[ex])*1000:+.1f} mm",
          flush=True)

# --- common slope from the DENSE gen2 surface (clean reference) --------------------------
fd = runs["fulldensity"]; nx, ny = fd["nx"], fd["ny"]; bnds = fd["bounds"]
Zf = fd["z_after"].copy(); nm = ~np.isfinite(Zf)
if nm.any():
    Zf = Zf[tuple(distance_transform_edt(nm, return_distances=False, return_indices=True))]
gy, gx = np.gradient(Zf, RES)
slope = np.degrees(np.arctan(np.hypot(gx, gy)))

# forest via full-density ground penetration; bedrock if the mask exists
pen = ground_penetration(AFTERS["fulldensity"], bnds, RES, nx, ny)
forest = pen < 0.25
try:
    bedrock = np.load("data/derived/elba/bedrock_mask.npy").astype(bool)
    if bedrock.shape != slope.shape: bedrock = None
except Exception:
    bedrock = None

BINS = [0, 5, 10, 15, 20, 25, 30, 90]
def curve(dod, mask):
    ex = np.isfinite(dod) & mask
    out = []
    for lo, hi in zip(BINS[:-1], BINS[1:]):
        m = ex & (slope >= lo) & (slope < hi)
        out.append((f"{lo}-{hi}", int(m.sum()),
                    float(np.median(dod[m])*1000) if m.any() else np.nan))
    return out

def robust_tanfit(dod, mask):
    ex = np.isfinite(dod) & mask & (slope < 60)
    t = np.tan(np.radians(slope[ex])); d = dod[ex]*1000
    # Theil-sen-ish: slope via median of pairwise not needed; use lstsq + report R2
    A = np.vstack([t, np.ones_like(t)]).T
    c, *_ = np.linalg.lstsq(A, d, rcond=None)
    pred = A @ c; ss = 1 - np.sum((d-pred)**2)/np.sum((d-np.mean(d))**2)
    return c[0], c[1], ss, ex.sum()

print("\n" + "="*72)
print("median DoD (mm) vs slope on FORESTED cells  (false aggradation = rising +)")
print(f"{'slope bin':>10} | {'decimated n':>12} {'dec mm':>8} | {'fulldens n':>11} {'FD mm':>8}")
cA = curve(runs["decimated"]["dod"], forest); cB = curve(runs["fulldensity"]["dod"], forest)
for (b, nA, mA), (_, nB, mB) in zip(cA, cB):
    print(f"{b:>10} | {nA:12,} {mA:+8.1f} | {nB:11,} {mB:+8.1f}")

for tag in ("decimated", "fulldensity"):
    s, b, r2, npt = robust_tanfit(runs[tag]["dod"], forest)
    print(f"  [{tag:11s}] DoD = {s:+.1f}*tan(S) {b:+.1f} mm   R2={r2:.2f}  (n={npt:,} forest cells)")

if bedrock is not None:
    print("\nBEDROCK control (stable, should be ~flat vs slope):")
    for tag in ("decimated", "fulldensity"):
        s, b, r2, npt = robust_tanfit(runs[tag]["dod"], bedrock)
        print(f"  [{tag:11s}] DoD = {s:+.1f}*tan(S) {b:+.1f} mm  (n={npt:,} bedrock cells)")
else:
    print("\n(no bedrock_mask.npy compatible with this grid; skipped)")

# persist the full-density products
import os; os.makedirs("data/derived/elba_fulldensity", exist_ok=True)
np.save("data/derived/elba_fulldensity/dod.npy", runs["fulldensity"]["dod"])
np.save("data/derived/elba_fulldensity/lod.npy", runs["fulldensity"]["lod"])
np.save("data/derived/elba_fulldensity/z_after.npy", runs["fulldensity"]["z_after"])
np.save("data/derived/elba_fulldensity/slope.npy", slope)
with open("data/derived/elba_fulldensity/corrections.json", "w") as fh:
    json.dump(runs["fulldensity"]["corrections"], fh, indent=2)
print("\nsaved -> data/derived/elba_fulldensity/  (dod, lod, z_after, slope, corrections)")
print("\nVERDICT: if fulldensity FORESTED median DoD is ~flat vs slope while decimated "
      "rises to +29 mm, the false hillslope aggradation was the DECIMATION -> fix = use "
      "the full ground-return set. Residual rise = a real leaf-on/veg term to handle next.")

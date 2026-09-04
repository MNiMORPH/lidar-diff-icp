#!/usr/bin/env python3
"""Does the DoD depend on ASPECT, and does a correction create or merely unmask it?

A smooth high-on-one-side / low-on-the-other pattern across aspect is the signature of
LATERAL misregistration: shift two surfaces horizontally and every slope reads high on the
side you shifted toward. It is not something a vertical percentile correction can produce,
but a correction that removes a LEVEL offset will make a pre-existing aspect pattern far more
visible, because the field then sits centred on zero.

Reports, per aspect octant on sloping ground: the covariate (class-2 SD), the percentile
applied, and the DoD before and after -- so "created" and "unmasked" can be told apart by
whether the PEAK-TO-PEAK grows or the LEVEL moves.

    env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python analysis/aspect_audit.py \
        --tile data/derived/elba
"""
import argparse, json
import numpy as np
from scipy.ndimage import distance_transform_edt

ap = argparse.ArgumentParser()
ap.add_argument("--tile", required=True)
ap.add_argument("--slope-min", type=float, default=5.0,
                help="aspect is undefined on flat ground and the misregistration signal "
                     "scales with tan(slope), so the test needs slope. MINE.")
A = ap.parse_args()

j = json.load(open(f"{A.tile}/corrections.json")); RES = float(j["res_m"])
z = np.load(f"{A.tile}/z_after.npy"); zf = z.copy(); m = ~np.isfinite(zf)
if m.any():
    zf = zf[tuple(distance_transform_edt(m, return_distances=False, return_indices=True))]
gN, gE = np.gradient(zf, RES)
aspect = np.degrees(np.arctan2(-gE, -gN)) % 360.0
slope = np.load(f"{A.tile}/slope.npy")
fp = np.load(f"{A.tile}/floodplain_mask.npy").astype(bool)
med = np.load(f"{A.tile}/dod_gen2_median.npy") * 1000
corr = np.load(f"{A.tile}/dod_cover_q2.npy") * 1000
sd = np.load(f"{A.tile}/class2_sd_mm.npy")
q = np.load(f"{A.tile}/gen2_q2_used.npy")
ok = (np.isfinite(med) & np.isfinite(corr) & np.isfinite(sd) & ~fp
      & (slope > A.slope_min))
print(f"{A.tile}: sloping ground (>{A.slope_min:g} deg), floodplain excluded: "
      f"{int(ok.sum()):,} cells")
print(f"\n{'aspect':>12s} {'n':>7s} {'SD mm':>7s} {'q':>6s} {'uncorr mm':>10s} "
      f"{'corr mm':>9s} {'shift mm':>9s}")
LAB = ["N", "NNE-NE", "ENE-E", "ESE-SE", "SSE-S", "SSW-SW", "WSW-W", "WNW-NW"]
U, C, S = [], [], []
for i, nm in enumerate(LAB):
    s = ok & (((aspect - (i * 45 - 22.5)) % 360) < 45)
    if s.sum() < 200:
        continue
    U.append(np.median(med[s])); C.append(np.median(corr[s])); S.append(np.median(sd[s]))
    print(f"{nm:>12s} {int(s.sum()):7,} {np.median(sd[s]):7.1f} {np.nanmedian(q[s]):6.3f} "
          f"{np.median(med[s]):10.2f} {np.median(corr[s]):9.2f} "
          f"{np.median(corr[s] - med[s]):9.2f}")
pp = lambda v: float(np.max(v) - np.min(v))
print(f"\n  peak-to-peak across aspect:  uncorrected {pp(U):.2f} mm   "
      f"corrected {pp(C):.2f} mm   class-2 SD {pp(S):.1f} mm")
print(f"  level:  uncorrected median {np.median(med[ok]):+.2f} mm   "
      f"corrected {np.median(corr[ok]):+.2f} mm")
print("\n  peak-to-peak roughly unchanged + level moved => the correction UNMASKED it.")
print("  peak-to-peak grown                            => the correction CREATED it.")

#!/usr/bin/env python3
"""Which cross-epoch datum does each derived tile actually carry, and what does the
retired one cost where it is still in force?

Two datum methods were removed from the pipeline in favour of the auto-computed geoid
difference: a fitted order-2 PARABOLA over "stable" ground, and a fitted reference PLANE.
Removing them from the code does not remove them from products already on disk, and a
product records which it used, so this reads that record rather than assuming.

For a tile still on the parabola it also evaluates BOTH surfaces over the tile grid --
the parabola from its stored coefficients (basis from coreg._poly_basis at 5e756ed:
[1, xn, yn, xn^2, xn*yn, yn^2], xn=(x-xm)/xhr) and the geoid-difference plane -- so the
difference is measured rather than argued.

    env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python analysis/datum_method_audit.py
    env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python analysis/datum_method_audit.py --compare carlton
"""
import argparse, glob, json, os

import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("--compare", default=None,
                help="tile name: evaluate its stored parabola against the geoid plane")
A = ap.parse_args()


def method_of(d):
    if "cross_epoch_datum" in d:
        return "geoid: " + str(d["cross_epoch_datum"].get("method"))
    if "cross_epoch_tie_order2_coef" in d:
        return "order-2 PARABOLA (retired)"
    return "NONE"


print(f"  {'tile':<20}{'cross-epoch tie':<28}{'swath_tie':<16}{'zero_line'}")
for p in sorted(glob.glob("data/derived/*/corrections.json")):
    t = os.path.basename(os.path.dirname(p))
    d = json.load(open(p))
    print(f"  {t:<20}{method_of(d):<28}{str(d.get('swath_tie')):<16}{d.get('zero_line')}")

if A.compare:
    from lidar_diff_icp import references, io
    D = f"data/derived/{A.compare}"
    c = json.load(open(f"{D}/corrections.json"))
    if "cross_epoch_tie_order2_coef" not in c:
        raise SystemExit(f"{A.compare} does not carry a parabola; nothing to compare")
    b = c["bounds"]; res = float(c["res_m"])
    dz = np.asarray(c["cross_epoch_tie_order2_coef"]["dz"], float)
    xm, xhr, ym, yhr = c["cross_epoch_tie_order2_coef"]["norm_xm_xhr_ym_yhr"]

    z = np.load(f"{D}/z_after.npy"); ny, nx = z.shape
    XX, YY = np.meshgrid(b[0] + (np.arange(nx) + 0.5) * res,
                         b[1] + (np.arange(ny) + 0.5) * res)
    xn = (XX - xm) / xhr; yn = (YY - ym) / yhr
    para = sum(dz[k] * v for k, v in enumerate(
        [np.ones_like(xn), xn, yn, xn * xn, xn * yn, yn * yn])) * 1000.0

    a, tb, tc = references.geoid_difference(tuple(b), io.MN_GEN1_CRS)
    cx, cy = 0.5 * (b[0] + b[2]), 0.5 * (b[1] + b[3])
    geo = (a + tb * (XX - cx) / 1000.0 + tc * (YY - cy) / 1000.0) * 1000.0

    fin = np.isfinite(z)
    print(f"\n{A.compare} grid {z.shape} at {res} m, {int(fin.sum()):,} finite cells\n")
    for nm, f in [("parabola (applied)", para), ("geoid plane", geo),
                  ("parabola - geoid", para - geo)]:
        v = f[fin]
        print(f"  {nm:<22} median {np.median(v):+8.2f}  min {v.min():+8.2f}  "
              f"max {v.max():+8.2f}  ptp {np.ptp(v):7.2f}   mm")
    print(f"\n  {A.compare} stable_1sigma_m = {c['stable_1sigma_m'] * 1000:.1f} mm "
          f"(the LoD scale this sits against)")

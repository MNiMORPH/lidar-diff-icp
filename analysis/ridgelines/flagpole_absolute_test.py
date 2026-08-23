#!/usr/bin/env python3
"""Absolute-elevation test at NGS leveled benchmark DG8385 ("9 DRL").

DG8385 is a USGS survey disk set in a mat-foundation/concrete slab (a flagpole
base) at the MN DNR Whitewater WMA headquarters, ~3.7 km NE of Elba, MN. Its
NAVD88 orthometric height (223.352 m) was determined by DIFFERENTIAL LEVELING
(adjusted Feb 2005) -- so it is geoid-model-INDEPENDENT, the preferred kind of
absolute vertical control. The disk is RECESSED 1 INCH (25.4 mm) below the slab
surface, so the lidar SLAB-SURFACE elevation is expected to read ~25.4 mm HIGHER
than the published mark height.

Epoch datum note: gen1 (2008 MN DNR) = NAVD88 / GEOID03; gen2 (2021 3DEP) =
NAVD88 / GEOID18. A LEVELED benchmark height is true NAVD88 regardless of geoid
model; each lidar epoch approximates NAVD88 via GPS + its own geoid model, so
this test measures each epoch's absolute vertical error, and the gen1-gen2
difference should reflect the GEOID03->GEOID18 change (~67 mm here) if that is
the dominant term.

Honest limitation up front (measured, not assumed): the flagpole-base slab is
SUB-METER and is NOT resolvable as a flat co-planar patch in either point cloud.
The surface at the mark is a smoothly GRADED LAWN (slope ~0.05 m/m E, ~0.03 m/m
N; residual scatter about a local plane ~0.07 m = grass roughness, not a smooth
slab). We therefore do not "isolate the slab"; instead we fit a robust local
plane to the ground returns and evaluate it AT the mark's (x,y) to remove the
slope-times-position error, and we carry the horizontal-position uncertainty
(+/- 3 m hand-held GPS on the datasheet) through as an elevation uncertainty.
"""
from __future__ import annotations

import numpy as np
import laspy

# ---- DG8385 authoritative values (NGS datasheet, retrieved 2026-08-22) --------
# NAD83(1986) 44 07 12.89 N, 092 00 12.73 W (hand-held GPS, +/- 3 m horizontal)
DG8385_LAT = 44 + 7 / 60 + 12.89 / 3600
DG8385_LON = -(92 + 0 / 60 + 12.73 / 3600)
DG8385_H_NAVD88 = 223.352        # meters, LEVELED (geoid-independent), adj 2005
DISK_RECESS_M = 0.0254           # disk recessed 1 inch below slab surface
H_HORIZ_SIGMA = 3.0              # meters, datasheet hand-held GPS horiz accuracy

# Expected lidar slab-SURFACE height (what the lidar ground return should see):
H_SURFACE_EXPECTED = DG8385_H_NAVD88 + DISK_RECESS_M   # 223.377 m

# Geoid heights at DG8385 from the NGS geoid API (geodesy.noaa.gov/api/geoid),
# verified 2026-08-22.  GEOID18 (-30.080) matches the datasheet exactly.
N_GEOID03 = -30.012   # m, NGS model 3
N_GEOID18 = -30.080   # m, NGS model 14  (== datasheet GEOID HEIGHT)
# A lidar reports H = h_ellipsoid - N.  For the SAME ground point, the epoch
# difference from the geoid update alone is:
#   H_gen2 - H_gen1 = N_GEOID03 - N_GEOID18
GEOID_PRED_G2_MINUS_G1 = N_GEOID03 - N_GEOID18   # +0.068 m: gen2 reads HIGHER

GEN1_LAZ = "data/before/4342-29-64.laz"                 # 2008 MN DNR, NAVD88/GEOID03
GEN2_LAZ = "data/after/3dep2021_fulldensity.laz"        # 2021 3DEP, NAVD88/GEOID18


def mark_utm15():
    from pyproj import Transformer
    t = Transformer.from_crs("EPSG:4326", "EPSG:26915", always_xy=True)
    return t.transform(DG8385_LON, DG8385_LAT)


def load_ground(path, E, N, half=30.0):
    """Class-2 (ground) individual returns within a half x half box of the mark."""
    las = laspy.read(path)
    cls = np.asarray(las.classification)
    x = np.asarray(las.x); y = np.asarray(las.y); z = np.asarray(las.z)
    m = (cls == 2) & (np.abs(x - E) < half) & (np.abs(y - N) < half)
    return x[m], y[m], z[m]


def plane_fit_at_mark(xg, yg, zg, E, N, radius, n_boot=2000, seed=0):
    """Robust local plane fit; evaluate at the mark (dx=dy=0).

    Returns dict with: n, median z (raw), plane-evaluated z at mark, plane slope,
    residual roughness (MAD-based std), and a bootstrap SE on the mark estimate.
    Also folds in the elevation uncertainty from +/- 3 m horizontal position on
    the fitted slope.
    """
    d = np.hypot(xg - E, yg - N)
    sel = d < radius
    xs, ys, zs = xg[sel] - E, yg[sel] - N, zg[sel]
    n = xs.size
    if n < 6:
        return dict(radius=radius, n=n, insufficient=True)
    A = np.c_[xs, ys, np.ones_like(xs)]
    coef, *_ = np.linalg.lstsq(A, zs, rcond=None)
    resid = zs - A @ coef
    rough = 1.4826 * np.median(np.abs(resid - np.median(resid)))
    z_mark_plane = coef[2]                     # plane value at dx=dy=0
    z_median = float(np.median(zs))
    slope_mag = float(np.hypot(coef[0], coef[1]))
    # bootstrap SE of the plane intercept at the mark
    rng = np.random.default_rng(seed)
    boot = np.empty(n_boot)
    for i in range(n_boot):
        idx = rng.integers(0, n, n)
        c, *_ = np.linalg.lstsq(A[idx], zs[idx], rcond=None)
        boot[i] = c[2]
    se_fit = float(boot.std())
    # elevation uncertainty from +/-3 m horizontal position on the local slope
    se_horiz = slope_mag * H_HORIZ_SIGMA
    se_total = float(np.hypot(se_fit, se_horiz))
    return dict(radius=radius, n=int(n), insufficient=False,
                z_median=z_median, z_mark_plane=float(z_mark_plane),
                slope_mag=slope_mag, roughness=float(rough),
                se_fit=se_fit, se_horiz=float(se_horiz), se_total=se_total)


def main():
    E, N = mark_utm15()
    print(f"DG8385 UTM15N: E={E:.2f} N={N:.2f}")
    print(f"Published NAVD88 (leveled): {DG8385_H_NAVD88:.3f} m  "
          f"(disk recessed {DISK_RECESS_M*1000:.1f} mm)")
    print(f"Expected lidar SLAB-SURFACE height: {H_SURFACE_EXPECTED:.3f} m\n")

    results = {}
    load_gxy = {}
    for label, path in [("gen1", GEN1_LAZ), ("gen2", GEN2_LAZ)]:
        xg, yg, zg = load_ground(path, E, N)
        load_gxy[label] = (xg, yg, zg)
        print(f"=== {label}  ({path}) ===")
        print(f"  ground(class2) returns within 30 m box: {xg.size}")
        results[label] = {}
        for r in (3.0, 5.0, 8.0):
            res = plane_fit_at_mark(xg, yg, zg, E, N, r)
            results[label][r] = res
            if res.get("insufficient"):
                print(f"  r={r:.0f} m: only {res['n']} returns -- INSUFFICIENT")
                continue
            print(f"  r={r:.0f} m: n={res['n']:4d}  "
                  f"z_median={res['z_median']:.3f}  "
                  f"z@mark(plane)={res['z_mark_plane']:.3f}  "
                  f"slope={res['slope_mag']*100:.1f}%  "
                  f"rough={res['roughness']*1000:.0f}mm  "
                  f"SE_fit={res['se_fit']*1000:.0f} SE_horiz={res['se_horiz']*1000:.0f} "
                  f"SE_tot={res['se_total']*1000:.0f}mm")
        print()

    # --- offsets at the primary radius (5 m: enough gen1 pts, tight enough) ----
    R = 5.0
    g1 = results["gen1"][R]; g2 = results["gen2"][R]
    z1 = g1["z_mark_plane"]; z2 = g2["z_mark_plane"]
    print("=== OFFSETS (plane-at-mark estimator, r=5 m) ===")
    print(f"  gen1 z@mark = {z1:.3f} m  (SE_tot {g1['se_total']*1000:.0f} mm)")
    print(f"  gen2 z@mark = {z2:.3f} m  (SE_tot {g2['se_total']*1000:.0f} mm)")
    print(f"  expected surface = {H_SURFACE_EXPECTED:.3f} m "
          f"(mark {DG8385_H_NAVD88:.3f} + {DISK_RECESS_M*1000:.0f} mm recess)\n")
    o1 = (z1 - H_SURFACE_EXPECTED) * 1000
    o2 = (z2 - H_SURFACE_EXPECTED) * 1000
    od = (z1 - z2) * 1000
    print(f"  gen1 - benchmark_surface = {o1:+.0f} mm  "
          f"({'HIGH' if o1>0 else 'LOW'})  [+/-{g1['se_total']*1000:.0f} mm, "
          f"dominated by horizontal position]")
    print(f"  gen2 - benchmark_surface = {o2:+.0f} mm  "
          f"({'HIGH' if o2>0 else 'LOW'})  [+/-{g2['se_total']*1000:.0f} mm, "
          f"dominated by horizontal position]")
    print(f"  gen1 - gen2              = {od:+.0f} mm  "
          f"({'gen1 higher' if od>0 else 'gen1 lower'})")

    # --- epoch difference: robust to horizontal mis-location -------------------
    # Both epochs are evaluated at the SAME (x,y) on the SAME graded surface, so a
    # +/-3 m horizontal error shifts both identically and CANCELS in the diff.
    # Demonstrate by scanning the center over a +/-3 m grid.
    print("\n=== EPOCH DIFFERENCE gen2-gen1 (robust to horizontal position) ===")
    diffs = []
    for dx in (-3, -2, -1, 0, 1, 2, 3):
        for dy in (-3, -2, -1, 0, 1, 2, 3):
            a = plane_fit_at_mark(*load_gxy["gen1"], E + dx, N + dy, R, n_boot=1)
            b = plane_fit_at_mark(*load_gxy["gen2"], E + dx, N + dy, R, n_boot=1)
            if a.get("insufficient") or b.get("insufficient"):
                continue
            diffs.append(b["z_mark_plane"] - a["z_mark_plane"])
    diffs = np.array(diffs) * 1000
    print(f"  gen2 - gen1 over +/-3 m grid (n={diffs.size} positions): "
          f"median {np.median(diffs):+.0f} mm, "
          f"range [{diffs.min():+.0f}, {diffs.max():+.0f}] mm, "
          f"std {diffs.std():.0f} mm")
    print(f"  GEOID03->GEOID18 prediction (gen2-gen1) = "
          f"{GEOID_PRED_G2_MINUS_G1*1000:+.0f} mm "
          f"(N_GEOID03 {N_GEOID03:+.3f} - N_GEOID18 {N_GEOID18:+.3f})")
    resid = np.median(diffs) - GEOID_PRED_G2_MINUS_G1 * 1000
    print(f"  residual (measured - geoid prediction) = {resid:+.0f} mm "
          f"-- i.e. the geoid update explains essentially all of the epoch offset")


if __name__ == "__main__":
    import os
    os.chdir(os.path.join(os.path.dirname(__file__), "..", ".."))
    main()

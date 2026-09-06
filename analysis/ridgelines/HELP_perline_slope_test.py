#!/usr/bin/env python3
"""HELP #3 -- the RIGHT discriminator: remove the per-FLIGHT-LINE gen1 datum offset,
then re-test the slope-deepening.  gen1-only.

WHAT HELP #2/#2b FOUND (and why this is the real test)
------------------------------------------------------
The beam-vs-aspect split (HELP #2, #2b) failed its own flat-ground sanity check:
even after incidence histogram-matching, flat forest showed a +12 mm up-vs-down
split where aspect is meaningless.  Root cause, verified: gen1 has a large
PER-FLIGHT-LINE datum/boresight offset.  On flat forest the floor d by line is:
    line 135: -91.6 mm   line 136: -69.8 mm   line 137: -21.2 mm   line 138: +5.0 mm
-- a ~97 mm spread at the SAME scan angle, unrelated to slope or aspect.  (This is
the 2008 along-track GPS-drift/boresight signature this project exists to correct.)
The "aspect" split was just re-sorting these lines.

CONSEQUENCE FOR THE SLOPE-DEEPENING: if the four flight lines are unevenly
distributed across slope bands (steep bands drawn more from the low-datum lines),
their fixed per-line offsets MASQUERADE as a slope trend.  So before claiming the
-2 mm/deg deepening is real ground, we must remove the per-line offset.

METHOD
------
For EACH flight line separately, anchor its forest floor to ITS OWN open-ground
floor at matched incidence (per-line datum), then read the forest-floor
differential vs slope WITHIN that line.  A per-line offset (constant in slope)
cancels.  If the slope-deepening survives WITHIN individual lines, it is not a
flight-line-composition artifact.  We also report, as the confound's magnitude,
how the line MIX shifts across slope bands.

DATA SOURCES: all gen1 (CSF ground returns; open/forest label from gen1-only
above-ground return fraction; incidence from validated reconstruction).  Aspect/
gen2 canopy NOT used.  d_mm's constant is removed per-line via the open anchor.

    env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python analysis/ridgelines/HELP_perline_slope_test.py
"""
import numpy as np, laspy, math
from scipy.ndimage import distance_transform_edt

NY, NX = 700, 508
X0, Y0 = 577492.8, 4882737.6
RES = 5.0
BUDGET = 20.0
CSF = "data/csf_cache/elba.las"
GEN1 = "data/before/4342-29-64.laz"
D = "data/derived/elba_fulldensity/"

Zg = np.load(D + "z_after.npy"); Zf = Zg.copy(); m = ~np.isfinite(Zf)
if m.any():
    Zf = Zf[tuple(distance_transform_edt(m, return_distances=False, return_indices=True))]
gy, gx = np.gradient(Zf, RES)
slope = np.degrees(np.arctan(np.hypot(gx, gy)))
gxf = gx.ravel(); gyf = gy.ravel(); nnorm = np.sqrt(gxf**2 + gyf**2 + 1.0)
Zflat = Zf.ravel(); slp_grid = slope.ravel()

tot = np.zeros(NY * NX); above = np.zeros(NY * NX)
with laspy.open(GEN1) as f:
    for pts in f.chunk_iterator(5_000_000):
        x = np.asarray(pts.x, np.float64); y = np.asarray(pts.y, np.float64)
        cl = np.asarray(pts.classification)
        ix = ((x - X0) / RES).astype(np.int64); iy = ((y - Y0) / RES).astype(np.int64)
        keep = (ix >= 0) & (ix < NX) & (iy >= 0) & (iy < NY) & (cl != 7)
        cell = iy[keep] * NX + ix[keep]; cl = cl[keep]
        np.add.at(tot, cell, 1)
        np.add.at(above, cell, (cl != 2).astype(float))
g1af = np.where(tot >= 5, above / np.maximum(tot, 1), np.nan)

las = laspy.read(CSF)
x = np.asarray(las.x, np.float64); y = np.asarray(las.y, np.float64); z = np.asarray(las.z, np.float64)
sa = np.asarray(las.scan_angle).astype(float) * 0.006
psid = np.asarray(las.point_source_id); gt = np.asarray(las.gps_time)
ix = ((x - X0) / RES).astype(np.int64); iy = ((y - Y0) / RES).astype(np.int64)
ing = (ix >= 0) & (ix < NX) & (iy >= 0) & (iy < NY)
cell = np.where(ing, iy * NX + ix, 0)

bhx = np.zeros(len(x)); bhy = np.zeros(len(x))
for pl in np.unique(psid):
    mm = psid == pl
    vx = np.polyfit(gt[mm], x[mm], 1)[0]; vy = np.polyfit(gt[mm], y[mm], 1)[0]
    H = math.atan2(vy, vx); cxu = np.array([-math.sin(H), math.cos(H)])
    cross = (x[mm] - x[mm].mean()) * cxu[0] + (y[mm] - y[mm].mean()) * cxu[1]
    sgn = np.sign(np.corrcoef(cross, sa[mm])[0, 1])
    bhx[mm] = -np.sign(sa[mm]) * sgn * cxu[0]; bhy[mm] = -np.sign(sa[mm]) * sgn * cxu[1]
th = np.radians(np.abs(sa)); bx = np.sin(th) * bhx; by = np.sin(th) * bhy; bz = np.cos(th)
inc = np.degrees(np.arccos(np.clip((bx * (-gxf[cell]) + by * (-gyf[cell]) + bz) / nnorm[cell], -1, 1)))
slp = slp_grid[cell]
xc = X0 + ((cell % NX) + 0.5) * RES; yc = Y0 + ((cell // NX) + 0.5) * RES
d = (z - (Zflat[cell] + gxf[cell] * (x - xc) + gyf[cell] * (y - yc))) * (1.0 / nnorm[cell]) * 1000

af = g1af[cell]
m0 = ing & np.isfinite(d) & np.isfinite(inc) & np.isfinite(slp) & np.isfinite(af)
FOREST = m0 & (af >= 0.55)
OPEN = m0 & (af <= 0.15)
lines = np.unique(psid)


def med(mask):
    return np.median(d[mask]) if mask.sum() else np.nan


# ---- 1. the confound: does the flight-line MIX shift across slope bands? ----
print("=" * 78)
print("STEP A -- the confound: per-line datum offset + line-mix vs slope")
print("=" * 78)
print("\nPer-line forest-floor datum (flat forest slope<3, all incidence):")
for pl in lines:
    b = FOREST & (psid == pl) & (slp < 3)
    print(f"  line {pl}: floor {med(b):+7.1f} mm  N={b.sum():,}")

SL_BANDS = [(3, 6), (6, 9), (9, 12), (12, 15), (15, 18), (18, 21), (21, 25)]
print("\nFlight-line MIX (%) of forest returns within each slope band:")
hdr = "  slope  " + "".join(f"line{pl:>4} " for pl in lines)
print(hdr)
for lo, hi in SL_BANDS:
    b = FOREST & (slp >= lo) & (slp < hi)
    n = b.sum()
    frac = [f"{(b & (psid == pl)).sum() / max(n, 1) * 100:5.0f} " for pl in lines]
    print(f"  {lo:2d}-{hi:2d}  " + "".join(frac) + f"  (N={n:,})")
print("  -> if the mix shifts toward low-datum lines as slope rises, the per-line")
print("     offset alone would fake a slope trend.  STEP B removes it.")

# ---- 2. per-line, incidence-matched, datum-removed slope trend ----
INC_BANDS = [(6, 10), (10, 14)]
print()
print("=" * 78)
print("STEP B -- WITHIN each flight line: forest-floor differential vs slope,")
print("          per-line open datum removed, incidence held fixed")
print("=" * 78)
all_line_trends = []
for pl in lines:
    print(f"\n  ---- flight line {pl} ----")
    for ilo, ihi in INC_BANDS:
        ref_b = OPEN & (psid == pl) & (slp < 6) & (inc >= ilo) & (inc < ihi)
        ref = med(ref_b)
        if not np.isfinite(ref) or ref_b.sum() < 300:
            print(f"    inc {ilo}-{ihi}: no per-line open ref (N={ref_b.sum()}), skip")
            continue
        xs, ys = [], []
        for lo, hi in SL_BANDS:
            b = FOREST & (psid == pl) & (slp >= lo) & (slp < hi) & (inc >= ilo) & (inc < ihi)
            if b.sum() < 400:
                continue
            xs.append((lo + hi) / 2); ys.append(med(b) - ref)
        if len(xs) > 2:
            A = np.polyfit(xs, ys, 1)
            span = A[0] * (xs[-1] - xs[0])
            all_line_trends.append(A[0])
            pts = "  ".join(f"{xx:.0f}d:{yy:+.0f}" for xx, yy in zip(xs, ys))
            print(f"    inc {ilo}-{ihi} (open ref {ref:+.0f}): {A[0]:+.2f} mm/deg, "
                  f"{span:+.0f} mm over {xs[0]:.0f}-{xs[-1]:.0f}d ({span/BUDGET:+.1f}x)")
            print(f"        [{pts}]")

if all_line_trends:
    t = np.array(all_line_trends)
    print(f"\n  Per-line slope trends (mm/deg-slope): {np.round(t,2)}")
    print(f"  mean {t.mean():+.2f}  median {np.median(t):+.2f}  min {t.min():+.2f}  max {t.max():+.2f}")

print()
print("=" * 78)
print("VERDICT LOGIC")
print(" - Per-line datum spread (STEP A) confirms a large gen1 boresight/GPS-drift")
print("   offset (~90+ mm across lines) -- this is a real confound.")
print(" - If, WITHIN individual lines (STEP B), the forest floor STILL deepens with")
print("   slope at ~-2 mm/deg after removing that line's own datum, the deepening is")
print("   NOT a flight-line-composition artifact -- it is intrinsic to each line's")
print("   own view of the ground, i.e. real ground OR a per-line slope/footprint bias")
print("   (which the aspect test could not isolate).")
print(" - If the per-line trends SCATTER around 0 / disagree in sign, the pooled")
print("   deepening was largely a line-mix artifact.")
print("=" * 78)

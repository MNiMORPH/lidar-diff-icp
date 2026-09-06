#!/usr/bin/env python3
"""HELP #4 -- the land-cover control: does OPEN (canopy-free) ground ALSO deepen
with slope?  gen1-only.

THE DISCRIMINATOR
-----------------
The slope-deepening survives a gen1-only forest label (HELP #1) and is consistent
in SIGN across all four flight lines when self-anchored (HELP #3).  Remaining
question: is it a FOREST-FLOOR signal (canopy x slope interaction, or real ground
change under trees) or a land-cover-INDEPENDENT slope-correlated geometry/DEM bias?

Clean test: repeat the matched-incidence, self-anchored slope profile on gen1-only
OPEN ground (above-ground return fraction <= 0.15 -- essentially bare farmland,
NO canopy).  If OPEN ground deepens with slope the SAME way forest does, the effect
is NOT canopy-driven: it is a slope-correlated ground/geometry effect common to all
land cover (footprint range-walk on the grazing beam, or a gen1<->gen2 slope-normal
residual).  If ONLY forest deepens, a canopy-under-slope mechanism is implicated.

CAVEAT (stated honestly): steep OPEN ground is SCARCE at Elba (open is mostly the
low-slope valley floor).  We can only populate open slopes up to ~12-15 deg, and
steep-open cells are geographically special (eroding banks, road/quarry cuts) where
REAL change is also plausible.  So this control is directional, not conclusive on
its own; read it together with HELP #3.

    env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python analysis/ridgelines/HELP_open_vs_forest_control.py [--tile elbaext]

The cover proxy stays GEN1-INTERNAL (gen1 above-ground return fraction) on purpose:
this asks a gen1-only mechanism question, so it must not be stratified by a
gen2-derived cover raster.
"""
import argparse, json, os
import numpy as np, laspy, math
from scipy.ndimage import distance_transform_edt

_ap = argparse.ArgumentParser()
_ap.add_argument("--tile", default="elba_fulldensity")
_ap.add_argument("--csf",  default=None, help="CSF ground cache (default data/csf_cache/<tile>.las)")
_ap.add_argument("--gen1", default=None, help="raw gen1 cloud for the ABOVE-GROUND FRACTION cover proxy")
_A = _ap.parse_args()
TILE = os.path.basename(_A.tile.rstrip("/"))
D = f"data/derived/{TILE}/"
_DEF_CSF = {"elba_fulldensity": "data/csf_cache/elba.las"}
_DEF_G1  = {"elba_fulldensity": "data/before/4342-29-64.laz",
            "elbaext":          "data/before/elbaext_gen1_merged.laz"}
CSF  = _A.csf  or _DEF_CSF.get(TILE, f"data/csf_cache/{TILE}.las")
GEN1 = _A.gen1 or _DEF_G1.get(TILE)
if GEN1 is None:
    raise SystemExit(f"no default gen1 cloud for tile {TILE}; pass --gen1")


def _grid(tile):                                        # (X0,Y0,NX,NY,RES) from tile meta/corrections
    for fn in ("meta.json", "corrections_geoid.json", "corrections.json"):
        p = f"data/derived/{tile}/{fn}"
        if os.path.exists(p):
            j = json.load(open(p)); b = j["bounds"]; r = float(j.get("res") or j.get("res_m"))
            nx = int(j.get("nx") or round((b[2]-b[0])/r)); ny = int(j.get("ny") or round((b[3]-b[1])/r))
            return b[0], b[1], nx, ny, r
    raise SystemExit(f"no grid meta for {tile}")


X0, Y0, NX, NY, RES = _grid(TILE)
BUDGET = 20.0

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
try:                                                    # scan angle -> DEGREES, both formats
    sa = np.asarray(las.scan_angle).astype(float) * 0.006      # PF6+ (0.006 deg units)
except Exception:
    sa = np.asarray(las.scan_angle_rank).astype(float)         # PF<=5 (integer degrees)
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
base = ing & np.isfinite(d) & np.isfinite(inc) & np.isfinite(slp) & np.isfinite(af)
FOREST = base & (af >= 0.55)
OPEN = base & (af <= 0.15)


def med(mask):
    return np.median(d[mask]) if mask.sum() else np.nan


SL_BANDS = [(3, 6), (6, 9), (9, 12), (12, 15), (15, 18), (18, 21), (21, 25)]
ILO, IHI = 6, 10  # the incidence band both open and forest populate best on slopes

print("=" * 78)
print(f"HELP #4 -- OPEN vs FOREST slope-deepening at MATCHED incidence 6-10 deg  [{TILE}]")
print("self-anchored to each cover's own flat (<3 deg) floor")
print("=" * 78)
for lbl, sel in [("OPEN (no canopy)", OPEN), ("FOREST", FOREST)]:
    ref = med(sel & (slp < 3) & (inc >= ILO) & (inc < IHI))
    refN = (sel & (slp < 3) & (inc >= ILO) & (inc < IHI)).sum()
    print(f"\n  {lbl}  (flat anchor {ref:+.1f} mm, N={refN:,}):")
    xs, ys = [], []
    for lo, hi in SL_BANDS:
        b = sel & (slp >= lo) & (slp < hi) & (inc >= ILO) & (inc < IHI)
        if b.sum() < 200:
            print(f"      slope {lo:2d}-{hi:2d}: (N={b.sum()} too few)")
            continue
        diff = med(b) - ref
        xs.append((lo + hi) / 2); ys.append(diff)
        print(f"      slope {lo:2d}-{hi:2d}: floor {med(b):+7.1f}  diff {diff:+6.1f} mm ({diff/BUDGET:+.1f}x)  N={b.sum():,}")
    if len(xs) > 2:
        A = np.polyfit(xs, ys, 1)
        span = A[0] * (xs[-1] - xs[0])
        print(f"      => {lbl} slope trend: {A[0]:+.2f} mm/deg, {span:+.0f} mm over {xs[0]:.0f}-{xs[-1]:.0f} deg")

print()
print("=" * 78)
print("INTERPRETATION")
print(" - OPEN ground deepening with slope COMPARABLE to forest => the effect is")
print("   land-cover-INDEPENDENT: a slope-correlated ground/geometry bias, NOT a")
print("   forest-floor / canopy signal.  (Caveat: steep-open N is small & special.)")
print(" - OPEN flat / forest-only deepening => canopy-under-slope is implicated.")
print("=" * 78)

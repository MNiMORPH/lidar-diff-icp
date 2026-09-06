#!/usr/bin/env python3
"""HELP #2 -- discriminate REAL slope-deepening from a gen1 FOOTPRINT/RANGE/ASPECT bias,
using beam-vs-slope-aspect geometry.  gen1-only.

THE LOGIC
---------
The audit established (and HELP_gen1only_strata.py confirms with a gen1-only label)
a slope-dependent DEEPENING of the gen1 forest floor of ~-2 mm/deg-slope at MATCHED
incidence, datum removed.  The open question: is that a REAL lowering of the 2008
ground surface on steep forested slopes, or a gen1 GEOMETRY bias that grows with slope?

A geometry/footprint/range bias has a DIRECTION: on a slope, a finite laser footprint
and a first-return-biased range gate do not sample the slope symmetrically.  A beam
arriving from the DOWNHILL side vs the UPHILL side of the slope aspect illuminates the
footprint's high edge vs low edge differently, and the elongated grazing footprint on a
slope biases the recorded range toward the near (uphill) edge.  So a footprint/aspect
bias must produce a SYSTEMATIC SPLIT in floor d between beams pointing UP-aspect and
DOWN-aspect, at the SAME slope and SAME incidence magnitude.

A REAL ground-surface deepening cannot know which way the beam came from: at matched
slope and matched incidence, up-aspect and down-aspect beams must read the SAME floor.

So the discriminator is:
   Delta_aspect(slope) = median d[beam points UP-aspect] - median d[beam points DOWN-aspect]
   at matched slope band and matched incidence band.
 - |Delta_aspect| small (<< the -2 mm/deg deepening) and slope-flat  => deepening is REAL.
 - |Delta_aspect| large and growing with slope, comparable to the deepening => the
   deepening is (at least partly) a gen1 footprint/range/aspect geometry bias.

We reconstruct the signed horizontal beam direction per flight line EXACTLY as the
validated incidence_angle.py does (heading from x,y vs gps_time; cross-track side from
corr(cross-track pos, scan_angle)), then project it on the slope-aspect (downhill)
direction to classify each return as up-aspect or down-aspect.

DATA SOURCES (labeled):
  - x,y,z,scan_angle,gps_time,point_source_id : gen1 CSF cloth-ground returns (raw las).
  - slope aspect (gx,gy) : from the gen2 reference plane z_after (the common DoD datum).
    NOTE: aspect direction is a property of the terrain, not of gen2 canopy; using the
    gen2 bare-earth gradient for aspect is the same unavoidable reference-plane use as
    inside d_mm, NOT a gen2 canopy magnitude.  Incidence itself is beam-vs-(gen2 normal).
  - forest label : gen1-only above-ground return fraction (as in HELP_gen1only_strata.py),
    recomputed here so this script is standalone and uses NO gen2-cut stratum.

    env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python analysis/ridgelines/HELP_beam_aspect_discriminator.py
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

# --- terrain: gen2 reference plane -> slope, aspect (downhill horizontal dir) ---
Zg = np.load(D + "z_after.npy"); Zf = Zg.copy(); m = ~np.isfinite(Zf)
if m.any():
    Zf = Zf[tuple(distance_transform_edt(m, return_distances=False, return_indices=True))]
gy, gx = np.gradient(Zf, RES)                       # gy=d/dNorth, gx=d/dEast (m/m)
slope = np.degrees(np.arctan(np.hypot(gx, gy)))
gxf = gx.ravel(); gyf = gy.ravel()
nnorm = np.sqrt(gxf**2 + gyf**2 + 1.0)
Zflat = Zf.ravel()
slp_grid = slope.ravel()
# downhill horizontal unit vector = +gradient of -z = (gx,gy) direction of steepest DESCENT
# z increases along (gx,gy)? np.gradient gives d z/dx; steepest ASCENT is (gx,gy); DESCENT is -(gx,gy).
downhill_x = np.where(np.hypot(gxf, gyf) > 0, -gxf / np.maximum(np.hypot(gxf, gyf), 1e-9), 0.0)
downhill_y = np.where(np.hypot(gxf, gyf) > 0, -gyf / np.maximum(np.hypot(gxf, gyf), 1e-9), 0.0)

# --- gen1-only forest label (above-ground return fraction), same as HELP #1 ---
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
g1_abovefrac = np.where(tot >= 5, above / np.maximum(tot, 1), np.nan)
FOREST_MIN = 0.55

# --- gen1 CSF ground returns: beam direction + incidence + d, per return ---
las = laspy.read(CSF)
x = np.asarray(las.x, np.float64); y = np.asarray(las.y, np.float64); z = np.asarray(las.z, np.float64)
sa = np.asarray(las.scan_angle).astype(float) * 0.006
psid = np.asarray(las.point_source_id); gt = np.asarray(las.gps_time)
ix = ((x - X0) / RES).astype(np.int64); iy = ((y - Y0) / RES).astype(np.int64)
ing = (ix >= 0) & (ix < NX) & (iy >= 0) & (iy < NY)
cell = np.where(ing, iy * NX + ix, 0)

bhx = np.zeros(len(x)); bhy = np.zeros(len(x))     # beam horizontal unit (ground->sensor)
for pl in np.unique(psid):
    mm = psid == pl
    vx = np.polyfit(gt[mm], x[mm], 1)[0]; vy = np.polyfit(gt[mm], y[mm], 1)[0]
    H = math.atan2(vy, vx)
    cxu = np.array([-math.sin(H), math.cos(H)])
    cross = (x[mm] - x[mm].mean()) * cxu[0] + (y[mm] - y[mm].mean()) * cxu[1]
    sgn = np.sign(np.corrcoef(cross, sa[mm])[0, 1])
    bhx[mm] = -np.sign(sa[mm]) * sgn * cxu[0]
    bhy[mm] = -np.sign(sa[mm]) * sgn * cxu[1]

th = np.radians(np.abs(sa))
bx = np.sin(th) * bhx; by = np.sin(th) * bhy; bz = np.cos(th)   # beam unit ground->sensor
inc = np.degrees(np.arccos(np.clip((bx * (-gxf[cell]) + by * (-gyf[cell]) + bz) / nnorm[cell], -1, 1)))
slp = slp_grid[cell]

xc = X0 + ((cell % NX) + 0.5) * RES; yc = Y0 + ((cell // NX) + 0.5) * RES
d = (z - (Zflat[cell] + gxf[cell] * (x - xc) + gyf[cell] * (y - yc))) * (1.0 / nnorm[cell]) * 1000  # mm

# --- classify: does the beam HORIZONTAL direction point UP-aspect or DOWN-aspect? ---
# beam horizontal (bhx,bhy) points from ground TOWARD sensor.  Project on DOWNHILL dir.
# proj>0 : sensor is downhill of the point => beam came from downhill => "beam points DOWN-aspect"
# proj<0 : sensor is uphill  of the point => beam came from uphill   => "beam points UP-aspect"
proj = bhx * downhill_x[cell] + bhy * downhill_y[cell]
bh_mag = np.hypot(bhx, bhy)
cos_aspect = np.where(bh_mag > 1e-6, proj / bh_mag, 0.0)   # cos of angle between beam-horiz and downhill
UP = ing & (cos_aspect < -0.5)    # beam clearly points up-aspect (>120deg from downhill)
DOWN = ing & (cos_aspect > 0.5)   # beam clearly points down-aspect (within 60deg of downhill)

forest = np.isfinite(g1_abovefrac[cell]) & (g1_abovefrac[cell] >= FOREST_MIN)
m0 = ing & forest & np.isfinite(d) & np.isfinite(inc) & np.isfinite(slp)


def med(mask):
    return np.median(d[mask]) if mask.sum() else np.nan


print("=" * 78)
print("HELP #2 -- beam-vs-aspect discriminator (gen1-only forest)")
print("real ground deepening: up- and down-aspect beams read the SAME floor.")
print("footprint/range/aspect bias: systematic split growing with slope.")
print("=" * 78)

# sanity: on FLAT ground there is no aspect, so the split should be ~0 by construction.
flat = m0 & (slp < 3)
print(f"\nSanity (flat forest slope<3, no meaningful aspect):")
print(f"  up-aspect med {med(flat & UP):+.1f}  down-aspect med {med(flat & DOWN):+.1f}  "
      f"split {med(flat & UP) - med(flat & DOWN):+.1f} mm  (N up={ (flat&UP).sum():,}, down={(flat&DOWN).sum():,})")

INC_BANDS = [(6, 10), (10, 14), (14, 18)]
SL_BANDS = [(6, 9), (9, 12), (12, 15), (15, 18), (18, 21), (21, 25), (25, 30)]
print("\nDelta_aspect = med(up-aspect) - med(down-aspect), at MATCHED slope & incidence:")
for ilo, ihi in INC_BANDS:
    print(f"\n  === incidence held {ilo}-{ihi} deg ===")
    print(f"    {'slope':>7} {'up med':>8} {'down med':>9} {'Delta':>7}   {'N_up':>7} {'N_down':>7}")
    ds = []
    for lo, hi in SL_BANDS:
        base = m0 & (slp >= lo) & (slp < hi) & (inc >= ilo) & (inc < ihi)
        u = base & UP; dn = base & DOWN
        if u.sum() < 400 or dn.sum() < 400:
            continue
        delta = med(u) - med(dn)
        ds.append(((lo + hi) / 2, delta))
        print(f"    {lo:2d}-{hi:2d}   {med(u):+8.1f} {med(dn):+9.1f} {delta:+7.1f}   {u.sum():>7,} {dn.sum():>7,}")
    if len(ds) > 2:
        sx = np.array([a for a, _ in ds]); sy = np.array([b for _, b in ds])
        A = np.polyfit(sx, sy, 1)
        print(f"    -> Delta_aspect trend vs slope: {A[0]:+.2f} mm/deg-slope "
              f"(mean |Delta| = {np.mean(np.abs(sy)):.1f} mm over {sx[0]:.0f}-{sx[-1]:.0f} deg)")

print()
print("=" * 78)
print("INTERPRETATION")
print(" - If |Delta_aspect| stays SMALL (few mm) and does NOT grow with slope, the")
print("   -2 mm/deg deepening is direction-INDEPENDENT => consistent with REAL ground.")
print(" - If |Delta_aspect| is LARGE and grows with slope (comparable to the deepening),")
print("   a gen1 footprint/range/aspect geometry bias is contributing to the deepening.")
print("=" * 78)

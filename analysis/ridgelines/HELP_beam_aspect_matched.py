#!/usr/bin/env python3
"""HELP #2b -- beam-vs-aspect discriminator, with the incidence confound REMOVED.

WHY THIS SUPERSEDES HELP_beam_aspect_discriminator.py
-----------------------------------------------------
The first version split forest returns into "beam points up-aspect" vs
"down-aspect" inside a 4-deg-wide incidence band and compared floor d.  That test
FAILED a sanity check: on flat forest (slope<3, where aspect is meaningless) it
still showed a +12 mm up-vs-down split.  Diagnosis: the up- and down-aspect groups
sit at DIFFERENT median incidence WITHIN the 4-deg band (11.0 vs 8.9 deg on flat
ground), and the intrinsic incidence effect (~+2 mm/deg) turns that 2.1-deg
difference into several mm of spurious split.  Beam-aspect direction and incidence
are geometrically COUPLED, so a coarse incidence band does not decouple them.

FIX: for each (slope, incidence) cell, resample the up- and down-aspect groups to a
COMMON incidence DISTRIBUTION (histogram matching in 0.5-deg incidence sub-bins,
downsampling the larger group per sub-bin), so the two groups have the SAME
incidence distribution to sub-degree precision.  Then any residual up-vs-down split
is genuinely the aspect/footprint effect, not leftover incidence.

Validated first on FLAT forest: after matching, the flat-ground split MUST collapse
to ~0.  Only then is the sloped-ground split trustworthy.

Same data sources as HELP #2 (all gen1; aspect from gen2 reference plane; forest
label from gen1-only above-ground return fraction; NO gen2 canopy magnitude).

    env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python analysis/ridgelines/HELP_beam_aspect_matched.py
"""
import numpy as np, laspy, math
from scipy.ndimage import distance_transform_edt

NY, NX = 700, 508
X0, Y0 = 577492.8, 4882737.6
RES = 5.0
CSF = "data/csf_cache/elba.las"
GEN1 = "data/before/4342-29-64.laz"
D = "data/derived/elba_fulldensity/"
rng = np.random.default_rng(0)

Zg = np.load(D + "z_after.npy"); Zf = Zg.copy(); m = ~np.isfinite(Zf)
if m.any():
    Zf = Zf[tuple(distance_transform_edt(m, return_distances=False, return_indices=True))]
gy, gx = np.gradient(Zf, RES)
slope = np.degrees(np.arctan(np.hypot(gx, gy)))
gxf = gx.ravel(); gyf = gy.ravel(); nnorm = np.sqrt(gxf**2 + gyf**2 + 1.0)
Zflat = Zf.ravel(); slp_grid = slope.ravel()
hyp = np.hypot(gxf, gyf)
downhill_x = np.where(hyp > 1e-6, -gxf / np.maximum(hyp, 1e-9), 0.0)
downhill_y = np.where(hyp > 1e-6, -gyf / np.maximum(hyp, 1e-9), 0.0)

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

proj = bhx * downhill_x[cell] + bhy * downhill_y[cell]
bh = np.hypot(bhx, bhy)
cosa = np.divide(proj, bh, out=np.zeros_like(proj), where=bh > 1e-6)
forest = np.isfinite(g1af[cell]) & (g1af[cell] >= 0.55)
UP = ing & forest & (cosa < -0.5)
DOWN = ing & forest & (cosa > 0.5)


def matched_split(sel_up, sel_dn, sub=0.5):
    """Histogram-match up/down on incidence in `sub`-deg bins; return matched medians."""
    iu, idn = inc[sel_up], inc[sel_dn]
    du, ddn = d[sel_up], d[sel_dn]
    if len(iu) < 100 or len(idn) < 100:
        return np.nan, np.nan, 0, 0
    lo = min(iu.min(), idn.min()); hi = max(iu.max(), idn.max())
    edges = np.arange(lo, hi + sub, sub)
    ku, kdn = [], []
    for a, b in zip(edges[:-1], edges[1:]):
        mu = np.where((iu >= a) & (iu < b))[0]
        mdn = np.where((idn >= a) & (idn < b))[0]
        k = min(len(mu), len(mdn))
        if k == 0:
            continue
        ku.append(rng.choice(mu, k, replace=False))
        kdn.append(rng.choice(mdn, k, replace=False))
    if not ku:
        return np.nan, np.nan, 0, 0
    ku = np.concatenate(ku); kdn = np.concatenate(kdn)
    return np.median(du[ku]), np.median(ddn[kdn]), len(ku), len(kdn)


print("=" * 78)
print("HELP #2b -- incidence-MATCHED beam-vs-aspect split (gen1-only forest)")
print("=" * 78)

# --- sanity: flat forest, incidence-matched, split MUST be ~0 ---
fu = UP & (slp < 3); fdn = DOWN & (slp < 3)
mu, mdn, nu, ndn = matched_split(fu, fdn)
print(f"\nSANITY flat forest slope<3, incidence-MATCHED:")
print(f"  up med {mu:+.1f}  down med {mdn:+.1f}  split {mu - mdn:+.1f} mm  (matched N up={nu:,}, down={ndn:,})")
print(f"  (unmatched flat split was +12.2 mm -- driven by 11.0 vs 8.9 deg incidence)")

print("\nMATCHED up-vs-down split at each slope band (incidence histogram-matched):")
print(f"  {'slope':>7} {'up med':>8} {'down med':>9} {'split':>7}   {'N(matched each)':>15}")
SL_BANDS = [(3, 6), (6, 9), (9, 12), (12, 15), (15, 18), (18, 21), (21, 25)]
rows = []
for lo, hi in SL_BANDS:
    su = UP & (slp >= lo) & (slp < hi)
    sdn = DOWN & (slp >= lo) & (slp < hi)
    mu, mdn, nu, ndn = matched_split(su, sdn)
    if not np.isfinite(mu):
        continue
    rows.append(((lo + hi) / 2, mu - mdn))
    print(f"  {lo:2d}-{hi:2d}   {mu:+8.1f} {mdn:+9.1f} {mu - mdn:+7.1f}   {min(nu, ndn):>15,}")
if len(rows) > 2:
    sx = np.array([a for a, _ in rows]); sy = np.array([b for _, b in rows])
    A = np.polyfit(sx, sy, 1)
    print(f"\n  matched split trend vs slope: {A[0]:+.2f} mm/deg-slope; "
          f"mean |split| {np.mean(np.abs(sy)):.1f} mm over {sx[0]:.0f}-{sx[-1]:.0f} deg")

print()
print("=" * 78)
print("INTERPRETATION")
print(" - If the matched flat-ground split is ~0 (confound removed) AND the matched")
print("   sloped split stays small (few mm, << the -40 mm slope-deepening), then the")
print("   deepening is direction-INDEPENDENT => consistent with REAL ground lowering.")
print(" - A large, slope-growing matched split => a gen1 footprint/aspect geometry bias")
print("   contributes to the deepening.")
print("=" * 78)

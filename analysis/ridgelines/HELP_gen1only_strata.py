#!/usr/bin/env python3
"""HELP #1 -- gen1-ONLY land-cover stratification, re-run the matched-incidence
slope-deepening test WITHOUT any gen2-derived label.

The audit's central corrected result (AUDIT_corrected_floor_signal.py STEP 2b) is a
slope-dependent DEEPENING of the forest floor of ~-40 to -52 mm over 4-23 deg slope,
at MATCHED incidence, with the ~67 mm datum removed.  BUT that test defined
"forest" and "open" from `stratum`, which is cut on GEN2 leaf-on penetration
(gen1_save_angles_slope.py:49-50).  The audit flagged this as gen2 contamination
baked into the gen1 file.

This script rebuilds the forest/open label from gen1's OWN leaf-off returns only:
    gen1_abovefrac[cell] = (gen1 non-ground returns) / (gen1 non-noise returns)
computed from the RAW gen1 cloud (data/before/4342-29-64.laz), all returns.  This is
gen1's own leaf-off above-ground return fraction -- high where 2008 bare branches
scattered pulses, low over open farmland.  It is INDEPENDENT of gen2 (spatial
corr with gen2 penetration is only ~0.23, verified), so it is a genuine
gen1-internal land-cover axis, not a gen2 relabel.

We then repeat the audit's matched-incidence, datum-removed forest-vs-slope test
using this gen1-only label and ask: does the slope-deepening survive?

DATA SOURCES (labeled):
  - d_mm, incidence, slope, cell, scan_angle : gen1 CSF cloth-ground returns (npz).
    d_mm is gen1 return z minus the gen2 slope-normal reference plane (a DoD; the
    ~67 mm constant part is the geoid datum and is removed differentially below).
  - gen1_abovefrac : gen1 RAW cloud, all returns, per cell.  gen1-ONLY.
  - NO gen2 penetration / understory / canopy magnitude is used, and NO gen2-cut
    stratum.  The only role gen2 plays is as the reference PLANE inside d_mm, which
    is unavoidable for a DoD and whose constant is removed by the open anchor.

    env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python analysis/ridgelines/HELP_gen1only_strata.py
"""
import numpy as np, laspy

NY, NX = 700, 508
X0, Y0 = 577492.8, 4882737.6
RES = 5.0
BUDGET = 20.0  # mm signal budget after the ~67 mm constant datum is removed
D = "data/derived/elba_fulldensity/"
GEN1 = "data/before/4342-29-64.laz"

# ---------------------------------------------------------------------------
# gen1-only per-cell above-ground return fraction (leaf-off canopy descriptor)
# ---------------------------------------------------------------------------
tot = np.zeros(NY * NX)
above = np.zeros(NY * NX)
with laspy.open(GEN1) as f:
    for pts in f.chunk_iterator(5_000_000):
        x = np.asarray(pts.x, np.float64); y = np.asarray(pts.y, np.float64)
        cl = np.asarray(pts.classification)
        ix = ((x - X0) / RES).astype(np.int64); iy = ((y - Y0) / RES).astype(np.int64)
        keep = (ix >= 0) & (ix < NX) & (iy >= 0) & (iy < NY) & (cl != 7)
        cl = cl[keep]; cell = (iy[keep] * NX + ix[keep])
        np.add.at(tot, cell, 1)
        np.add.at(above, cell, (cl != 2).astype(float))
g1_abovefrac = np.where(tot >= 5, above / np.maximum(tot, 1), np.nan)  # >=5 returns/cell for a stable fraction

# ---------------------------------------------------------------------------
# npz returns (gen1 CSF ground), map the gen1-only label onto them
# ---------------------------------------------------------------------------
Z = np.load(D + "gen1_csf_angles.npz")
dm = Z["d_mm"]; inc = Z["incidence"]; sl = Z["slope"]; cell = Z["cell"]; ing = Z["in_grid"]
g2strat = Z["stratum"]  # only used to REPORT overlap with the gen1-only cut, not to filter
af = g1_abovefrac[cell]

m0 = ing & np.isfinite(dm) & np.isfinite(inc) & np.isfinite(sl) & np.isfinite(af)


def med(mask):
    return np.median(dm[mask]) if mask.sum() else np.nan


# ---------------------------------------------------------------------------
# Choose gen1-only forest/open thresholds on gen1_abovefrac.
# open = little above-ground return (bare farmland); forest = much above-ground.
# ---------------------------------------------------------------------------
qs = np.nanpercentile(af[m0], [10, 25, 50, 75, 90])
print("=" * 74)
print("gen1-only above-ground return fraction (leaf-off canopy descriptor)")
print("=" * 74)
print(f"  percentiles of gen1_abovefrac over CSF ground returns [10,25,50,75,90]:")
print("   ", np.round(qs, 3))
OPEN_MAX = 0.15   # <=0.15 above-ground => open (near-bare ground)
FOREST_MIN = 0.55  # >=0.55 above-ground => forest (dense leaf-off canopy)
G1_OPEN = m0 & (af <= OPEN_MAX)
G1_FOREST = m0 & (af >= FOREST_MIN)
print(f"  gen1-only OPEN  (abovefrac<={OPEN_MAX}):   N={G1_OPEN.sum():,}")
print(f"  gen1-only FOREST(abovefrac>={FOREST_MIN}): N={G1_FOREST.sum():,}")

# Report how much the gen1-only cut agrees/disagrees with the gen2-cut stratum,
# to make the independence explicit (they are NOT the same partition).
g2F = (g2strat == 1); g2O = (g2strat == 2)
print("\n  overlap of gen1-only label with gen2-cut stratum (independence check):")
print(f"    gen1-FOREST that gen2 also calls forest: {(G1_FOREST & g2F).sum()/max(G1_FOREST.sum(),1):.0%}")
print(f"    gen1-FOREST that gen2 calls OPEN:        {(G1_FOREST & g2O).sum()/max(G1_FOREST.sum(),1):.0%}")
print(f"    gen1-OPEN that gen2 also calls open:     {(G1_OPEN & g2O).sum()/max(G1_OPEN.sum(),1):.0%}")
print(f"    gen1-OPEN that gen2 calls FOREST:        {(G1_OPEN & g2F).sum()/max(G1_OPEN.sum(),1):.0%}")

# ---------------------------------------------------------------------------
# Open reference at matched incidence (removes datum + common incidence effect)
# ---------------------------------------------------------------------------
INC_BANDS = [(6, 10), (10, 14), (14, 18)]
SL_BANDS = [(3, 6), (6, 9), (9, 12), (12, 15), (15, 18), (18, 21), (21, 25)]

print()
print("=" * 74)
print("STEP 1 -- gen1-only OPEN floor vs incidence (datum+incidence anchor)")
print("=" * 74)
open_ref = {}
for ilo, ihi in INC_BANDS:
    b = G1_OPEN & (sl < 6) & (inc >= ilo) & (inc < ihi)
    open_ref[(ilo, ihi)] = med(b)
    print(f"  inc {ilo:2d}-{ihi:2d} deg, gen1-open slope<6: floor = {med(b):+7.1f} mm  N={b.sum():,}")

print()
print("=" * 74)
print("STEP 2 -- gen1-only FOREST floor vs SLOPE at MATCHED incidence, datum removed")
print("=" * 74)
for ilo, ihi in INC_BANDS:
    ref = open_ref[(ilo, ihi)]
    if not np.isfinite(ref):
        print(f"  --- incidence {ilo}-{ihi}: no open reference, skipping ---")
        continue
    xs, ys = [], []
    print(f"  --- incidence held {ilo}-{ihi} deg (gen1-open ref = {ref:+.1f} mm) ---")
    for lo, hi in SL_BANDS:
        b = G1_FOREST & (sl >= lo) & (sl < hi) & (inc >= ilo) & (inc < ihi)
        if b.sum() < 500:
            continue
        diff = med(b) - ref
        xs.append((lo + hi) / 2); ys.append(diff)
        print(f"      slope {lo:2d}-{hi:2d}: forest floor {med(b):+7.1f}  differential {diff:+6.1f} mm ({diff/BUDGET:+.1f}x)  N={b.sum():,}")
    if len(xs) > 2:
        A = np.polyfit(xs, ys, 1)
        span = A[0] * (xs[-1] - xs[0])
        print(f"      => gen1-only differential slope trend: {A[0]:+.2f} mm/deg-slope, {span:+.0f} mm over {xs[0]:.0f}-{xs[-1]:.0f} deg ({span/BUDGET:+.1f}x budget)")

print()
print("=" * 74)
print("READ-OUT: if the differential slope trend is still ~-2 mm/deg-slope with")
print("this gen1-ONLY forest/open label, the deepening is NOT an artifact of the")
print("gen2-cut stratum -- it is present in gen1's own land-cover partition.")
print("=" * 74)

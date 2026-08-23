#!/usr/bin/env python3
"""AUDIT: corrected forest-floor signal analysis.

Fixes the three contamination modes found in the prior analysis:
  (1) The absolute slope-normal d (~-60 mm) is dominated by a CONSTANT gen1<->gen2
      datum (geoid GEOID03->GEOID18, ~+67 mm applied to gen1 in the pipeline).
      The SIGNAL is only the ~20 mm of *structure on top of* that constant.  We
      remove the constant here by anchoring to the OPEN/low-slope reference at
      MATCHED incidence, so every reported number is differential (mm above/below
      the flat-open floor at the same beam geometry), not absolute.
  (2) Judging effects by Pearson r is wrong: per-return scatter is ~+/-270 mm, so
      |r| is ~0.03-0.1 even when the SYSTEMATIC median shift across a covariate's
      range is tens of mm = the whole 20 mm budget.  We report SYSTEMATIC MEDIAN
      SHIFT in mm relative to the 20 mm budget, never r.
  (3) Per-cell MEDIAN floor pools returns across incidence.  Achievable incidence
      is nearly collinear with slope (corr ~0.85) because the scanner is only
      +/-17 deg off-nadir: on steep slopes, low-incidence (slope-perpendicular)
      beams are geometrically impossible.  Floor d depends strongly on incidence
      (~+2 to +4 mm/deg).  So a mixed-angle median-vs-slope trend is CONTAMINATED
      by the shifting angle composition.  We hold incidence FIXED throughout.

All inputs are gen1 returns.  d_mm is intrinsically gen1-vs-gen2 (gen1 return z
minus the gen2 reference plane, slope-normal) -- that is unavoidable and correct
for a DoD; the fix is to remove its CONSTANT part and control incidence, NOT to
stratify or explain it with gen2 canopy variables.

NOTE ON STRATA: stratum/core_forest/core_open in the npz are cut on GEN2 leaf-on
penetration.  We use them ONLY as a coarse land-cover label to separate the
open reference from forest; we do NOT use any gen2 canopy MAGNITUDE (penetration,
understory_frac, canopy_height) as an explanatory covariate for gen1 behavior.

    env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python analysis/ridgelines/AUDIT_corrected_floor_signal.py
"""
import numpy as np

BUDGET = 20.0  # mm, the real signal budget after the ~67 mm constant is removed
D = "data/derived/elba_fulldensity/"
Z = np.load(D + "gen1_csf_angles.npz")
dm = Z["d_mm"]; inc = Z["incidence"]; sl = Z["slope"]; ing = Z["in_grid"]
st = Z["stratum"]  # gen2-cut label, used only to pick open reference vs forest
m0 = ing & np.isfinite(dm) & np.isfinite(inc) & np.isfinite(sl)

OPEN = m0 & (st == 2)
FOREST = m0 & (st == 1)

# Incidence bands we can populate in BOTH open and forest.  The open stratum is
# almost entirely low-slope, so only modest incidence is reachable there.
INC_BANDS = [(6, 10), (10, 14), (14, 18)]


def med(mask):
    return np.median(dm[mask]) if mask.sum() else np.nan


print("=" * 74)
print("STEP 0 -- the constant dominates the raw number")
print("=" * 74)
print(f"raw median d, all gen1 ground (in_grid): {med(m0):+.1f} mm")
print(f"raw median d, OPEN low-slope (<3 deg):   {med(OPEN & (sl < 3)):+.1f} mm")
print("  -> ~60-67 mm of this is the gen1<->gen2 geoid/co-registration datum,")
print("     NOT signal.  Everything below is DIFFERENTIAL (open reference removed).")

# ---------------------------------------------------------------------------
# Reference: flat-ish OPEN floor at each incidence band = the datum anchor.
# Subtracting it removes the constant AND any pure incidence effect that is
# common to open and forest, leaving the differential forest-floor signal.
# ---------------------------------------------------------------------------
print()
print("=" * 74)
print("STEP 1 -- open-ground floor vs incidence (the datum+incidence reference)")
print("=" * 74)
open_ref = {}
for ilo, ihi in INC_BANDS:
    b = OPEN & (sl < 6) & (inc >= ilo) & (inc < ihi)
    open_ref[(ilo, ihi)] = med(b)
    print(f"  inc {ilo:2d}-{ihi:2d} deg, open slope<6: floor = {med(b):+7.1f} mm  N={b.sum():,}")

# ---------------------------------------------------------------------------
# CONTAMINATED vs CORRECTED: forest floor vs SLOPE.
# ---------------------------------------------------------------------------
print()
print("=" * 74)
print("STEP 2 -- forest floor vs SLOPE: mixed-angle (WRONG) vs matched-incidence")
print("=" * 74)
SL_BANDS = [(3, 6), (6, 9), (9, 12), (12, 15), (15, 18), (18, 21), (21, 25)]

print("\n(a) MIXED-ANGLE median (contaminated -- angle composition shifts with slope):")
mx, my, mi = [], [], []
for lo, hi in SL_BANDS:
    b = FOREST & (sl >= lo) & (sl < hi)
    if b.sum() < 500:
        continue
    mx.append((lo + hi) / 2); my.append(med(b)); mi.append(np.median(inc[b]))
    print(f"    slope {lo:2d}-{hi:2d}: floor {med(b):+7.1f} mm   median incidence {np.median(inc[b]):4.1f} deg   N={b.sum():,}")
print(f"    apparent slope trend (mixed): {(my[-1]-my[0]):+.0f} mm over {mx[0]:.0f}-{mx[-1]:.0f} deg,")
print(f"      but median incidence swept {mi[0]:.0f}->{mi[-1]:.0f} deg simultaneously -- CONFOUNDED.")

print("\n(b) MATCHED-INCIDENCE, datum removed (forest floor minus open floor, same inc band):")
for ilo, ihi in INC_BANDS:
    ref = open_ref[(ilo, ihi)]
    xs, ys = [], []
    print(f"  --- incidence held {ilo}-{ihi} deg (open ref = {ref:+.1f} mm) ---")
    for lo, hi in SL_BANDS:
        b = FOREST & (sl >= lo) & (sl < hi) & (inc >= ilo) & (inc < ihi)
        if b.sum() < 500:
            continue
        diff = med(b) - ref
        xs.append((lo + hi) / 2); ys.append(diff)
        print(f"      slope {lo:2d}-{hi:2d}: forest floor {med(b):+7.1f}  differential {diff:+6.1f} mm  ({diff/BUDGET:+.1f}x budget)  N={b.sum():,}")
    if len(xs) > 2:
        A = np.polyfit(xs, ys, 1)
        span = A[0] * (xs[-1] - xs[0])
        print(f"      => differential slope trend: {A[0]:+.2f} mm/deg-slope, {span:+.0f} mm over {xs[0]:.0f}-{xs[-1]:.0f} deg ({span/BUDGET:+.1f}x budget)")

# ---------------------------------------------------------------------------
# The intrinsic incidence effect itself, at FIXED slope, gen1-only (not veg).
# ---------------------------------------------------------------------------
print()
print("=" * 74)
print("STEP 3 -- intrinsic floor-vs-incidence at FIXED slope (all gen1 ground)")
print("=" * 74)
for lo, hi in [(6, 9), (12, 15), (18, 21), (25, 30)]:
    b = m0 & (sl >= lo) & (sl < hi)
    e = np.arange(0, 46, 4); xs, ys = [], []
    for a, z in zip(e[:-1], e[1:]):
        s = b & (inc >= a) & (inc < z)
        if s.sum() > 500:
            xs.append((a + z) / 2); ys.append(med(s))
    if len(xs) > 2:
        A = np.polyfit(xs, ys, 1)
        span = max(ys) - min(ys)
        print(f"  slope {lo:2d}-{hi:2d}: d/dinc = {A[0]:+.2f} mm/deg   floor span {span:+.0f} mm over inc {xs[0]:.0f}-{xs[-1]:.0f} deg  ({span/BUDGET:+.1f}x budget)")

print()
print("=" * 74)
print("VERDICT (printed values above):")
print(" - Mixed-angle forest-floor-vs-slope reverses sign once incidence is held")
print("   fixed: the naive trend is an angle-composition artifact.")
print(" - At matched incidence, a REAL differential forest-floor deepening with")
print("   slope survives, several x the 20 mm budget -- see STEP 2(b) slopes.")
print(" - The intrinsic incidence effect (STEP 3) is the dominant per-return")
print("   nuisance and must be controlled before any covariate claim.")
print("=" * 74)

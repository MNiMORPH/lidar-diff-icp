#!/usr/bin/env python3
"""RE-EXAMINE THE BEDROCK / CLIFF CONTROL (Task 4, gen1-low reopened).

The bedrock outcrop control "refuted" pulse-broadening: bare rock on slopes did NOT show
the +tan(S) bias (it read -47 mm at steep vs +29 forested). But it lumped ALL 743 outcrop
cells, INCLUDING near-vertical cliff faces where a few-cm horizontal registration error x a
near-vertical slope = huge vertical DoD noise (NMAD 169 mm). That control may be dominated
by cliffs, not clean bare-rock SLOPES -- so it may not validly test broadening.

STRATIFY the bedrock cells by slope AND by per-epoch ground-return count. The clean test of
a cover-independent geometric mechanism (broadening) is on MODERATE bedrock slopes (~15-35
deg) with ADEQUATE returns and MANAGEABLE noise. If those show +tan(S) too, broadening/a
geometric term is back in play; if they're flat/negative even there, broadening stays out.

    env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python analysis/slope_bias/cliff_cells.py
"""
import laspy, numpy as np

RES = 5.0
X0, Y0 = 577492.8, 4882737.6                     # elba AOI grid origin (matches dod.npy)
dod = np.load("data/derived/elba_fulldensity/dod.npy")
slope = np.load("data/derived/elba_fulldensity/slope.npy")
bed = np.load("data/derived/elba/bedrock_mask.npy").astype(bool)
ny, nx = dod.shape

def count_ground(fn, cls=2, chunk=20_000_000):
    """class-2 ground returns per 5 m cell on the AOI grid (row 0 = south)."""
    g = np.zeros(ny * nx, np.int64)
    with laspy.open(fn) as fh:
        for pts in fh.chunk_iterator(chunk):
            c = np.asarray(pts.classification); m = c == cls
            px = np.asarray(pts.x)[m]; py = np.asarray(pts.y)[m]
            ix = ((px - X0) / RES).astype(np.int64); iy = ((py - Y0) / RES).astype(np.int64)
            ok = (ix >= 0) & (ix < nx) & (iy >= 0) & (iy < ny)
            g += np.bincount(iy[ok] * nx + ix[ok], minlength=ny * nx)
    return g.reshape(ny, nx)

print("counting per-cell ground returns (gen1 class-2, gen2 full class-2) ...", flush=True)
g1 = count_ground("data/before/4342-29-64.laz")             # gen1 vendor ground
g2 = count_ground("data/after/3dep2021_fulldensity.laz")    # gen2 full ground

def nmad(a): return 1.4826 * np.median(np.abs(a - np.median(a))) if a.size else np.nan

fin = np.isfinite(dod)
BINS = [0, 15, 25, 35, 45, 60, 90]
print(f"\nbedrock cells: {bed.sum()} total, {int((bed&fin).sum())} with finite DoD")
print("\n=== BEDROCK/CLIFF cells stratified by slope ===")
print(f"{'slope':>8} {'n':>5} {'medDoD':>7} {'NMAD':>6} {'g1/cell':>7} {'g2/cell':>7} "
      f"{'%g1<3':>6} {'%g2<3':>6}")
for lo, hi in zip(BINS[:-1], BINS[1:]):
    m = bed & fin & (slope >= lo) & (slope < hi)
    if m.sum() == 0:
        print(f"{lo:>3}-{hi:<3}  0"); continue
    d = dod[m] * 1000
    print(f"{lo:>3}-{hi:<3} {m.sum():>5} {np.median(d):>7.0f} {nmad(d):>6.0f} "
          f"{np.median(g1[m]):>7.0f} {np.median(g2[m]):>7.0f} "
          f"{100*np.mean(g1[m]<3):>5.0f}% {100*np.mean(g2[m]<3):>5.0f}%")

# CLEAN bedrock test: moderate slope + adequate returns in BOTH epochs + not a noise cell
clean = bed & fin & (slope >= 15) & (slope < 40) & (g1 >= 3) & (g2 >= 6)
noisy = bed & fin & (slope >= 45)
print(f"\nCLEAN bare-rock slopes (15-40 deg, g1>=3 & g2>=6): n={clean.sum()}")
if clean.sum():
    d = dod[clean] * 1000
    print(f"  median DoD = {np.median(d):+.1f} mm   NMAD = {nmad(d):.0f} mm   "
          f"slope med = {np.median(slope[clean]):.0f} deg")
    # tan(S) trend within the clean set
    t = np.tan(np.radians(slope[clean]))
    A = np.vstack([t, np.ones_like(t)]).T
    c, *_ = np.linalg.lstsq(A, d, rcond=None)
    print(f"  DoD ~ {c[0]:+.0f}*tan(S) {c[1]:+.0f} mm  (median-based bins more robust below)")
    for lo, hi in [(15,20),(20,25),(25,30),(30,40)]:
        mm = clean & (slope>=lo) & (slope<hi)
        if mm.sum(): print(f"    {lo}-{hi} deg: n={mm.sum():>4}  medDoD={np.median(dod[mm]*1000):+.1f} mm")
print(f"\nNEAR-VERTICAL cliffs (slope>=45): n={noisy.sum()}  "
      f"medDoD={np.median(dod[noisy]*1000):+.0f} mm  NMAD={nmad(dod[noisy]*1000):.0f} mm  "
      f"(g1/cell med={np.median(g1[noisy]):.0f}) -- registration x slope noise")
print("\nRead: does CLEAN bare-rock slope carry the +tan(S) bias? If yes -> a cover-"
      "independent (geometric/gen1-low) term is real; the old 'bedrock refutes broadening' "
      "rested on cliff noise. If flat/negative with tight NMAD -> broadening stays refuted.")

#!/usr/bin/env python3
"""Does gen1's ground surface float on dense grass?  (Andy's hypothesis, 2026-08-28)

gen1's per-cell ground elevation is the MEDIAN of its CSF-classified ground returns. If a
dense sward puts returns above true ground, that median floats up, and DoD = gen2 - gen1
reads as EROSION. gen2 has ~12x the ground-return density (5.78 vs 0.49 pts/m^2), so far
more of its shots reach the true surface through the same grass -- this is a SAMPLING
argument, not a leaf-state one, and it is independent of the canopy-cover work.

The diagnostic is gen1-INTERNAL, as it must be: p50 - p10 of each cell's OWN gen1 ground
returns (`d_mm` in gen1_csf_angles.npz). A surface read cleanly gives a small lift; returns
sitting on grass give a large one-sided upward spread.

Controls that matter, and are reported: canopy cover (so this is not the leaf-on effect
again) and slope (intra-cell relief also widens the spread -- but the floodplain is FLAT,
so if it still shows a bigger lift than the upland, relief cannot be the cause).

    ./lidar-icp/bin/python analysis/gen1_grass_lift.py
"""
import argparse, os
import numpy as np

ap = argparse.ArgumentParser(description=__doc__,
                             formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--tile", default="data/derived/elba_fulldensity")
ap.add_argument("--min-returns", type=int, default=8)
ap.add_argument("--save", action="store_true", help="write the per-cell lift grid")
A = ap.parse_args()
T = A.tile

d = np.load(f"{T}/gen1_csf_angles.npz", allow_pickle=True)
cell = d["cell"]; dmm = d["d_mm"].astype(float); ing = d["in_grid"]
dod = np.load(f"{T}/dod.npy"); lod = np.load(f"{T}/lod.npy"); ny, nx = dod.shape
flood = np.load(f"{T}/floodplain_mask.npy").astype(bool).ravel()
cov = np.load(f"{T}/canopy_cover_pfs.npy").ravel()
slope = np.load(f"{T}/slope.npy").ravel()

k = ing & (cell >= 0) & (cell < ny * nx) & np.isfinite(dmm)
cell, dmm = cell[k], dmm[k]
o = np.argsort(cell, kind="stable"); c, v = cell[o], dmm[o]
uc, st = np.unique(c, return_index=True); sp = np.r_[st[1:], c.size]
s = (sp - st) >= A.min_returns
p50 = np.array([np.percentile(v[a:b], 50) for a, b in zip(st[s], sp[s])])
p10 = np.array([np.percentile(v[a:b], 10) for a, b in zip(st[s], sp[s])])
cid = uc[s]; lift = p50 - p10
fl = flood[cid]; covg = cov[cid]; slg = slope[cid]
dg = dod.ravel()[cid] * 1000; lg = lod.ravel()[cid] * 1000

print(f"gen1 ground returns in grid {cell.size:,}; cells with >= {A.min_returns} returns "
      f"{int(s.sum()):,} of {uc.size:,}")
print(f"\np50 - p10 of each cell's OWN gen1 ground returns (mm) -- the upward lift")
print(f"  {'group':26s} {'ncell':>7} {'p50-p10':>9} {'slope':>7} {'cover':>7} {'DoD mm':>9}")
for nm, m in (("floodplain", fl), ("upland", ~fl),
              ("floodplain, cover<0.15", fl & (covg < 0.15)),
              ("upland,     cover<0.15", (~fl) & (covg < 0.15)),
              ("floodplain, cover>=0.30", fl & (covg >= 0.30)),
              ("upland,     cover>=0.30", (~fl) & (covg >= 0.30))):
    if m.sum() < 50:
        continue
    print(f"  {nm:26s} {int(m.sum()):7d} {np.nanmedian(lift[m]):9.1f} {np.nanmedian(slg[m]):7.1f} "
          f"{np.nanmedian(covg[m]):7.3f} {np.nanmedian(dg[m]):+9.1f}")

er = fl & np.isfinite(dg) & np.isfinite(lg) & (dg < -lg)
print(f"\n  floodplain EROSION cells n={int(er.sum()):,}  lift {np.nanmedian(lift[er]):.1f} mm  "
      f"slope {np.nanmedian(slg[er]):.1f} deg  cover {np.nanmedian(covg[er]):.3f}  "
      f"DoD {np.nanmedian(dg[er]):+.1f} mm")
print(f"  floodplain, not eroding  n={int((fl & ~er).sum()):,}  lift {np.nanmedian(lift[fl & ~er]):.1f} mm  "
      f"slope {np.nanmedian(slg[fl & ~er]):.1f} deg  cover {np.nanmedian(covg[fl & ~er]):.3f}  "
      f"DoD {np.nanmedian(dg[fl & ~er]):+.1f} mm")
print(f"\n  NOT a correction, a diagnostic: re-gridding gen1 at p10 would raise the DoD by the")
print(f"  lift, {np.nanmedian(lift[fl]):.1f} mm in the floodplain and {np.nanmedian(lift[er]):.1f} mm on")
print(f"  the eroding cells. p10 is a number nobody has justified -- it is not proposed here.")

if A.save:
    G = np.full(ny * nx, np.nan); G[cid] = lift
    np.save(f"{T}/gen1_lift_p50_p10.npy", G.reshape(ny, nx))
    print(f"\nwrote {T}/gen1_lift_p50_p10.npy")

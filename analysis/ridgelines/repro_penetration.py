"""Is data/derived/elba_fulldensity/penetration.npy the output of
canopy.ground_penetration() on the full-density gen2 cloud?

Same arithmetic as canopy.py:19-33, accumulated over chunks instead of one full read.
bincount is additive and the division happens once at the end, so this is the SAME
computation, not an approximation of it.
"""
import numpy as np, laspy

X0, Y0, X1, Y1 = 577492.8, 4882737.6, 580032.8, 4886237.6
RES, NX, NY = 5.0, 508, 700
GROUND_CLASS, NOISE_CLASS = 2, 7          # canopy.ground_penetration defaults

tot = np.zeros(NX * NY, float)
gnd = np.zeros(NX * NY, float)
n_read = 0
with laspy.open("data/after/3dep2021_fulldensity.laz") as fh:
    for pts in fh.chunk_iterator(20_000_000):
        x = np.asarray(pts.x); y = np.asarray(pts.y)
        cl = np.asarray(pts.classification)
        m = (x >= X0) & (x < X1) & (y >= Y0) & (y < Y1) & (cl != NOISE_CLASS)
        cid = ((y[m] - Y0) / RES).astype(int) * NX + ((x[m] - X0) / RES).astype(int)
        tot += np.bincount(cid, minlength=NX * NY)
        gnd += np.bincount(cid[cl[m] == GROUND_CLASS], minlength=NX * NY)
        n_read += len(x)
        print(f"  {n_read:,} points", flush=True)

frac = np.where(tot > 0, gnd / np.maximum(tot, 1), np.nan).reshape(NY, NX)
stored = np.load("data/derived/elba_fulldensity/penetration.npy")

print(f"\nrecomputed: finite {np.isfinite(frac).sum():,}  cells with zero returns "
      f"{int((tot == 0).sum()):,}")
print(f"stored:     finite {np.isfinite(stored).sum():,}")
same = np.array_equal(frac, stored, equal_nan=True)
print(f"\nBIT-IDENTICAL: {same}")
if not same:
    d = frac - stored
    f = np.isfinite(d)
    print(f"  differing cells: {int((d[f] != 0).sum()):,} of {int(f.sum()):,}")
    if f.any() and (d[f] != 0).any():
        print(f"  max |difference|: {np.abs(d[f]).max():.6g}")
    print(f"  NaN disagreement: recomputed-NaN-stored-finite "
          f"{int((~np.isfinite(frac) & np.isfinite(stored)).sum()):,}; "
          f"stored-NaN-recomputed-finite "
          f"{int((np.isfinite(frac) & ~np.isfinite(stored)).sum()):,}")

# what did the stored file put in the cells that have NO returns at all?
empty = (tot == 0).reshape(NY, NX)
v = stored[empty]
print(f"\nthe {int(empty.sum()):,} cells with ZERO gen2 returns hold, in the stored file:")
u, c = np.unique(v, return_counts=True)
for uu, cc in zip(u[:6], c[:6]):
    print(f"    value {uu!r}  in {cc:,} cells")
print(f"  distinct values: {u.size}")
fr = stored[empty] < 0.25
op = stored[empty] >= 0.45
print(f"\n  under the strata cuts used across the repo (forest pen<0.25, open pen>=0.45):")
print(f"    classified FOREST: {int(fr.sum()):,}   OPEN: {int(op.sum()):,}   "
      f"mixed: {int((~fr & ~op).sum()):,}")
cf = np.load("data/derived/elba_fulldensity/core_forest.npy")
print(f"  of those, currently inside core_forest: {int((empty & cf).sum()):,}")

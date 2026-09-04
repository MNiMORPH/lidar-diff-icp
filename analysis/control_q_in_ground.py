#!/usr/bin/env python3
# RECOVERED 2026-09-04 from the session transcript. This was written, run, and its
# output shipped -- but the code itself was never committed, so the CSV in
# data/derived had no producer and the measurement under it was not reproducible.
# Restored verbatim apart from this header and the output path.
"""Where does the SURVEYED ground fall within the lidar's own GROUND-CLASS returns?

Not within the whole column -- that was q_ctrl, and it is dominated by vegetation. This
asks the operational question the pipeline actually faces: it takes the per-cell MEDIAN of
class-2 returns (ground_q = 0.50). Is 0.50 right, and does it drift with cover?

If the rank is 0.50 everywhere, the class-2 median IS the ground and the error lies
elsewhere. If it drifts, that drift is the percentile correction, measured against truth.
"""
import os, sys, numpy as np, pandas as pd, laspy
sys.path.insert(0, "analysis")
from control_lowveg_offset import lowveg
from control_mode_shift import CONTROL, STRUCT, BOX, marks
from lidar_diff_icp.groundtruth.tie import _design
SET = "gen2_2021_control"
rows = []
for t in marks(SET).itertuples():
    sp, bp = f"{STRUCT}/{SET}__{t.point_id}.npz", f"{BOX}/{SET}__{t.point_id}.laz"
    if not (os.path.exists(sp) and os.path.exists(bp)):
        continue
    z = np.load(sp); coef = z["surface_coef"]
    E, N, R = float(z["easting"]), float(z["northing"]), float(z["struct_radius"])
    f = laspy.read(bp)
    x, y, zz, cl = np.asarray(f.x), np.asarray(f.y), np.asarray(f.z), np.asarray(f.classification)
    sel = (np.hypot(x - E, y - N) <= R)
    g = sel & (cl == 2)
    if g.sum() < 20:
        continue
    nn = np.sqrt(1 + coef[1]**2 + coef[2]**2)
    hg = (zz[g] - (_design(x[g]-E, y[g]-N, 2) @ coef)) / nn
    mu = (float(t.elevation) - float(coef[0])) / nn
    rows.append(dict(point_id=t.point_id, n_g=int(g.sum()), n_all=int(sel.sum()),
                     q_in_ground=float(np.mean(hg < mu)),          # rank of truth in class 2
                     ground_med_minus_truth=float(np.median(hg) - mu) * 1000,
                     class2_frac=float(g.sum() / sel.sum())))
F = pd.DataFrame(rows)
F["lowveg"] = [lowveg(p, 0.15, 2.00, setname=SET) for p in F.point_id]
F.to_csv("data/derived/control_q_in_ground.csv", index=False)
print(f"  gen2, {len(F)} marks with >=20 class-2 returns in the 7.5 m disc")
print(f"  RANK OF SURVEYED GROUND WITHIN THE CLASS-2 RETURNS")
print(f"    median {F.q_in_ground.median():.4f}   mean {F.q_in_ground.mean():.4f}   "
      f"p16 {F.q_in_ground.quantile(.16):.4f}   p84 {F.q_in_ground.quantile(.84):.4f}")
print(f"    at exactly 0 (truth below every ground return): {(F.q_in_ground<=0).sum()}   "
      f"at 1 (above all): {(F.q_in_ground>=1).sum()}")
print(f"\n  class-2 MEDIAN minus truth: median {F.ground_med_minus_truth.median():+.1f} mm   "
      f"mean {F.ground_med_minus_truth.mean():+.1f}   sd {F.ground_med_minus_truth.std():.1f}")
q = np.nanpercentile(F.lowveg, np.linspace(0, 100, 6))
print(f"\n  {'lowveg range':>18s} {'n':>4s} {'rank median':>12s} {'class2 med - truth mm':>22s}")
for i in range(5):
    s = (F.lowveg >= q[i]) & (F.lowveg <= q[i+1])
    if s.sum() < 5: continue
    print(f"  {q[i]:8.4f}-{q[i+1]:<8.4f} {int(s.sum()):4d} {F[s].q_in_ground.median():12.4f} "
          f"{F[s].ground_med_minus_truth.median():22.1f}")

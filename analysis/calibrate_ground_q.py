#!/usr/bin/env python3
"""Calibrate ground_q from the class-2 spread, against surveyed control.

WHAT THIS REPLACES. The pipeline takes the per-cell MEDIAN of class-2 returns as the ground
(`ground_q = 0.50`). Measured against 519 surveyed marks, that is right where the ground
class is clean and wrong where it is not -- and the class-2 spread says which is which,
with no cover product, no windows and no external layer.

The finding: the rank of true ground within a cell's class-2 returns is FLAT (~0.57) while
the class-2 spread is no wider than bare-ground noise, and FALLS once it is wider -- i.e.
once something other than ground is in the class. Measured: rho +0.006 (p 0.94) below 60 mm,
rho -0.183 (p 0.001) above, medians 0.571 vs 0.390, Mann-Whitney p 3.5e-06.

No break is imposed. An ISOTONIC (monotone non-increasing) regression of rank on log spread
reproduces the flat-then-falling shape without a threshold, without a functional form, and
without any cutoff to defend.

The test that decides it: held-out marks, spatially blocked. Does taking the q(SD) percentile
of a mark's class-2 returns land closer to surveyed ground than taking the median?
"""
import argparse
import numpy as np, pandas as pd, os, sys, laspy
sys.path.insert(0, "analysis")
from sklearn.isotonic import IsotonicRegression
from control_mode_shift import CONTROL, STRUCT, BOX, marks
from lidar_diff_icp.groundtruth.tie import _design
_ap = argparse.ArgumentParser()
_ap.add_argument("--set", dest="set_", default="gen2_2021_control",
                 choices=["gen1_2008_control", "gen2_2021_control"],
                 help="which epoch's control marks to calibrate on. The curve is NOT "
                      "transferable between epochs: 2008 was flown leaf-off in November and "
                      "2021 at green-up in May, and the classifiers differ too.")
_ap.add_argument("--out", default=None)
_A = _ap.parse_args()
SET = _A.set_
OUTNPZ = _A.out or f"data/derived/ground_q_vs_class2sd_{SET}.npz"

rows = []
for t in marks(SET).itertuples():
    sp, bp = f"{STRUCT}/{SET}__{t.point_id}.npz", f"{BOX}/{SET}__{t.point_id}.laz"
    if not (os.path.exists(sp) and os.path.exists(bp)):
        continue
    z = np.load(sp); coef = z["surface_coef"]
    E, N, R = float(z["easting"]), float(z["northing"]), float(z["struct_radius"])
    f = laspy.read(bp)
    x, y, zz, cl = np.asarray(f.x), np.asarray(f.y), np.asarray(f.z), np.asarray(f.classification)
    g = (np.hypot(x - E, y - N) <= R) & (cl == 2)
    if g.sum() < 20:
        continue
    nn = np.sqrt(1 + coef[1]**2 + coef[2]**2)
    hg = np.sort((zz[g] - (_design(x[g]-E, y[g]-N, 2) @ coef)) / nn)
    mu = (float(t.elevation) - float(coef[0])) / nn
    rows.append(dict(point_id=t.point_id, easting=t.easting, northing=t.northing,
                     sd=float(np.std(hg)), rank=float(np.mean(hg < mu)),
                     mu=mu, hg=hg))
F = pd.DataFrame(rows)
ls = np.log(F.sd.to_numpy() * 1000)
rk = F["rank"].to_numpy()
blk = (F.easting//10000).astype(int).astype(str) + "_" + (F.northing//10000).astype(int).astype(str)
ub = blk.unique(); rng = np.random.default_rng(0); rng.shuffle(ub)
fold = {b: i % 5 for i, b in enumerate(ub)}; f5 = blk.map(fold).to_numpy()

def err_at(q_vec):
    """|estimated ground - truth| in mm, taking the q-th percentile of each mark's class 2."""
    out = np.empty(len(F))
    for i, (hg, mu, q) in enumerate(zip(F.hg, F.mu, q_vec)):
        out[i] = (np.quantile(hg, min(max(q, 0.0), 1.0)) - mu) * 1000
    return out

e50 = err_at(np.full(len(F), 0.50))
e57 = err_at(np.full(len(F), float(np.median(rk))))
pred = np.empty(len(F))
for k in range(5):
    tr, te = f5 != k, f5 == k
    if te.sum() == 0 or tr.sum() < 30: pred[te] = np.median(rk[tr]); continue
    iso = IsotonicRegression(increasing=False, out_of_bounds="clip").fit(ls[tr], rk[tr])
    pred[te] = iso.predict(ls[te])
eiso = err_at(pred)
print(f"  gen2, {len(F)} marks, 5-fold spatially blocked (10 km) CV")
print(f"  {'ground estimator':>34s} {'median err':>11s} {'|median|':>9s} {'RMS':>8s} {'p90|err|':>9s}")
for nm, e in (("q = 0.50 (pipeline default)", e50),
              (f"q = {np.median(rk):.3f} constant (calibrated)", e57),
              ("q = isotonic(log class-2 SD), held out", eiso)):
    print(f"  {nm:>34s} {np.median(e):11.1f} {abs(np.median(e)):9.1f} "
          f"{np.sqrt(np.mean(e**2)):8.1f} {np.percentile(np.abs(e),90):9.1f}")
iso_full = IsotonicRegression(increasing=False, out_of_bounds="clip").fit(ls, rk)
print(f"\n  calibration curve q(class-2 SD), fitted on all {len(F)} marks:")
for v in (20, 40, 60, 80, 120, 200, 400):
    print(f"    SD {v:4d} mm -> q = {float(iso_full.predict([np.log(v)])[0]):.3f}")
# The curve travels with its provenance: a consumer must be able to see what population it
# was fitted on and what it is entitled to be applied to.
np.savez(OUTNPZ,
         log_sd_mm=iso_full.f_.x, q=iso_full.f_.y, n_marks=len(F),
         set=SET,
         fitted_on=f"{SET}: control marks, class-2 returns within 7.5 m",
         response="rank of surveyed ground within the mark's class-2 returns",
         covariate="natural log of the class-2 standard deviation, mm",
         shape="isotonic, monotone non-increasing -- flat then falling, no break imposed",
         cv="5-fold spatially blocked on 10 km blocks",
         known_limits=("calibrated on 7.5 m discs and applied to 5 m cells; per-epoch -- "
                       "this curve is valid only for the epoch named in `set`; corrects a "
                       "population, not a cell"))
print(f"\n  wrote {OUTNPZ}")

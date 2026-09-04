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
from control_mode_shift import CONTROL, STRUCT, BOX, marks
from lidar_diff_icp import groundq
from lidar_diff_icp.groundtruth.tie import _design
_ap = argparse.ArgumentParser()
_ap.add_argument("--set", dest="set_", default="gen2_2021_control",
                 choices=["gen1_2008_control", "gen2_2021_control"],
                 help="which epoch's control marks to calibrate on. The curve is NOT "
                      "transferable between epochs: 2008 was flown leaf-off in November and "
                      "2021 at green-up in May, and the classifiers differ too.")
_ap.add_argument("--out", default=None)
_ap.add_argument("--diagnostics", action="store_true",
                 help="also print every number quoted in "
                      "analysis/GROUND_Q_FROM_CLASS2_SPREAD.md: the rank summary, the "
                      "covariate comparison against lowveg, the spread-bin table with "
                      "bootstrap CIs, and the two-regime test. Without this the doc is only "
                      "half reproducible, which is the gap that keeps recurring on this "
                      "project.")
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
    # The covariate and the response are groundq's, so the statistic fitted here is
    # literally the one groundq.spread_from_histogram measures on a tile.
    rows.append(dict(point_id=t.point_id, easting=t.easting, northing=t.northing,
                     mu=mu, hg=hg, **groundq.mark_statistics(hg, mu)))
F = pd.DataFrame(rows)
SD_MM = F.sd.to_numpy() * 1000
ls = np.log(SD_MM)
rk = F["rank"].to_numpy()
f5, _blocks = groundq.spatial_folds(F.easting.to_numpy(), F.northing.to_numpy())

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
    iso = groundq.fit_curve(SD_MM[tr], rk[tr])
    pred[te] = iso.predict(ls[te])
eiso = err_at(pred)
print(f"  gen2, {len(F)} marks, 5-fold spatially blocked (10 km) CV")
print(f"  {'ground estimator':>34s} {'median err':>11s} {'|median|':>9s} {'RMS':>8s} {'p90|err|':>9s}")
for nm, e in (("q = 0.50 (pipeline default)", e50),
              (f"q = {np.median(rk):.3f} constant (calibrated)", e57),
              ("q = isotonic(log class-2 SD), held out", eiso)):
    print(f"  {nm:>34s} {np.median(e):11.1f} {abs(np.median(e)):9.1f} "
          f"{np.sqrt(np.mean(e**2)):8.1f} {np.percentile(np.abs(e),90):9.1f}")
if _A.diagnostics:
    from scipy.stats import spearmanr, mannwhitneyu
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from control_lowveg_offset import lowveg
    F["lowveg"] = [lowveg(p, 0.15, 2.00, setname=SET) for p in F.point_id]
    print(f"\n  RANK OF SURVEYED GROUND WITHIN THE CLASS-2 RETURNS")
    print(f"    median {F['rank'].median():.4f}   mean {F['rank'].mean():.4f}   "
          f"p16 {F['rank'].quantile(.16):.4f}   p84 {F['rank'].quantile(.84):.4f}")
    print(f"    at exactly 0 (truth below every ground return): {(F['rank']<=0).sum()}   "
          f"at 1 (above all): {(F['rank']>=1).sum()}")
    print(f"    class-2 MEDIAN minus truth: median {F.med_minus_truth.median():+.1f} mm   "
          f"mean {F.med_minus_truth.mean():+.1f}   sd {F.med_minus_truth.std():.1f}")
    print(f"\n  COVARIATE COMPARISON (vs the rank). lowveg is built from TWO chosen bands;")
    print(f"  the spread measures from none.")
    print(f"    {'covariate':>12s} {'windows?':>10s} {'rho':>8s} {'p':>11s}")
    for c, wnd in (("sd","none"), ("nmad","none"), ("iqr","none"), ("skew","none"),
                   ("lowveg","TWO")):
        m = np.isfinite(F[c]) & np.isfinite(F["rank"])
        r = spearmanr(F[c][m], F["rank"][m])
        print(f"    {c:>12s} {wnd:>10s} {r.statistic:+8.3f} {r.pvalue:11.2e}")
    x = F.sd.to_numpy() * 1000; y = F["rank"].to_numpy()
    rng2 = np.random.default_rng(0)
    print(f"\n  SHAPE: median rank per class-2 SD bin, bootstrap SE (2000 draws)")
    print(f"    {'SD bin mm':>14s} {'n':>5s} {'median rank':>12s} {'SE':>7s} {'95% CI':>16s}")
    for a, b in zip([0,30,45,60,80,110,160,250], [30,45,60,80,110,160,250,1e9]):
        sel = (x >= a) & (x < b)
        if sel.sum() < 8:
            continue
        v = y[sel]
        bs = np.array([np.median(rng2.choice(v, v.size)) for _ in range(2000)])
        lab = f"{a:.0f}-{b:.0f}" if b < 1e8 else f">{a:.0f}"
        print(f"    {lab:>14s} {int(sel.sum()):5d} {np.median(v):12.3f} {bs.std():7.3f} "
              f"[{np.percentile(bs,2.5):5.3f},{np.percentile(bs,97.5):5.3f}]")
    lo, hi = x < 60, x >= 60
    print(f"\n  TWO REGIMES, split at the 59.3 mm bare-ground class-2 NMAD:")
    print(f"    within <60 mm:  rho {spearmanr(x[lo],y[lo]).statistic:+.3f}  "
          f"p {spearmanr(x[lo],y[lo]).pvalue:.3f}   n={int(lo.sum())}")
    print(f"    within >=60 mm: rho {spearmanr(x[hi],y[hi]).statistic:+.3f}  "
          f"p {spearmanr(x[hi],y[hi]).pvalue:.3e}   n={int(hi.sum())}")
    print(f"    medians {np.median(y[lo]):.3f} vs {np.median(y[hi]):.3f}   "
          f"Mann-Whitney p {mannwhitneyu(y[lo],y[hi]).pvalue:.3e}")

iso_full = groundq.fit_curve(SD_MM, rk)
print(f"\n  calibration curve q(class-2 SD), fitted on all {len(F)} marks:")
for v in (20, 40, 60, 80, 120, 200, 400):
    print(f"    SD {v:4d} mm -> q = {float(iso_full.predict([np.log(v)])[0]):.3f}")
# The curve travels with its provenance: a consumer must be able to see what population it
# was fitted on and what it is entitled to be applied to.
groundq.save_curve(
    OUTNPZ, iso_full, n_marks=len(F), epoch=SET,
    fitted_on=f"{SET}: control marks, class-2 returns within 7.5 m",
    response="rank of surveyed ground within the mark's class-2 returns",
    covariate="natural log of the class-2 standard deviation, mm",
    shape="isotonic, monotone non-increasing -- flat then falling, no break imposed",
    cv="5-fold spatially blocked on 10 km blocks",
    known_limits=("calibrated on 7.5 m discs and applied to 5 m cells; per-epoch -- "
                  "this curve is valid only for the epoch named in `set`; corrects a "
                  "population, not a cell"))
print(f"\n  wrote {OUTNPZ}")

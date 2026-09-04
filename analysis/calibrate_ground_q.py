#!/usr/bin/env python3
"""Calibrate ground_q from the class-2 spread, against surveyed control.

READ THE RESULT BEFORE USING THE CURVE. On open ground this correction does not help. Held
out on 5 folds of 10 km spatially blocked marks, over the 227 NVA (non-vegetated) marks:

           q = 0.50 (pipeline default)        -3.5      3.5     49.1     73.1
       q = 0.527 constant (calibrated)         0.1      0.1     48.7     73.9
    q = isotonic(log class-2 SD), held out    -5.8      5.8     52.5     76.5

The plain median is already within 3.5 mm of truth and the curve costs more RMS than the bias
it removes. `difference_dem` therefore defaults to ground_q = 0.50 and requires a curve to be
named before it will use one.

WHAT THE CURVE IS. An ISOTONIC (monotone non-increasing) regression of the rank of true ground
within a mark's class-2 returns on log(class-2 spread). No break is imposed, no functional
form, no cutoff. The only constraint is physical: more contamination cannot mean a higher
ground rank.

WHY --point-types IS REQUIRED. The control set is three populations, and which of them you fit
decides what the curve means:

    NVA  n=227   class-2 median  -3.5 mm from truth   non-vegetated, open ground
    VVA  n=162                 +103.3 mm              sited UNDER VEGETATION, by design
    LCP  n=130                  -23.1 mm              the acquisition's calibration points

Pooling all three produced a curve whose falling limb is entirely the VVA marks, and a
headline gain (RMS 124.5 -> 104.6 mm) that vanishes on open ground. This mirrors
ground_control/run_bridge_gen2.py, where --point-types has been required all along:
"--point-types NVA is the consequence, not a preference".

    ./lidar-icp/bin/python analysis/calibrate_ground_q.py --set gen2_2021_control \
        --point-types NVA --diagnostics
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
_ap.add_argument("--point-types", nargs="+", required=True,
                 help="which control point types to fit on, e.g. NVA. REQUIRED and with no "
                      "default, exactly as ground_control/run_bridge_gen2.py requires it: "
                      "'--point-types NVA is the consequence, not a preference'. The three "
                      "types are different populations and pooling them bakes the canopy "
                      "response into the percentile. MEASURED on these marks: NVA (open "
                      "ground, n=227) puts the class-2 median -3.5 mm from truth and needs "
                      "almost no spread dependence; VVA (sited UNDER vegetation by design, "
                      "n=162) puts it +103.3 mm high; LCP (the 143 calibration points, "
                      "n=130) -23.1 mm. The pooled curve's whole falling limb is the VVA "
                      "population, and the +8.1 mm bias that motivated this correction is a "
                      "pooling artifact.")
_ap.add_argument("--shape", default="piecewise", choices=["piecewise", "isotonic"],
                 help="piecewise (default, Andy 2026-09-04): two straight segments in log "
                      "spread with the break scanned, fitted to equal-count bin medians "
                      "weighted by their bootstrap SEs. isotonic: the earlier staircase. "
                      "Both are stored as knots and applied by identical code.")
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
TYPES = [t.upper() for t in _A.point_types]
SHAPE = _A.shape
OUTNPZ = _A.out or (f"data/derived/ground_q_vs_class2sd_{SET}_"
                    f"{'-'.join(sorted(TYPES))}.npz")

_M = marks(SET)
_types = pd.read_csv(CONTROL[SET])[["point_id", "point_type"]].drop_duplicates("point_id")
_M = _M.merge(_types, on="point_id", how="left")
_keep = _M.point_type.str.upper().isin(TYPES)
print(f"point types: fitting on {TYPES}; {int(_keep.sum())} of {len(_M)} marks kept "
      f"({_M.point_type.value_counts().to_dict()})")
_M = _M[_keep]

rows = []
for t in _M.itertuples():
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
                     point_type=t.point_type, n_g=int(g.sum()), mu=mu, hg=hg,
                     **groundq.mark_statistics(hg, mu)))
F = pd.DataFrame(rows)

# THE PER-MARK TABLE. Written because it was NOT: the relationship this whole route rests on
# -- the percentile of true ground, and the class-2 spread it is indexed by, one row per mark
# -- existed only as data/derived/control_q_vs_sigma.csv, which NOTHING in the tree produced.
# It had been made by ad-hoc code and the producer never committed, so the curve was
# reproducible and the measurement under it was not. groundq.mark_statistics reproduces that
# orphan file exactly (519/519 marks, every column, max |diff| 1.4e-14), so this supersedes it.
MARKS_OUT = OUTNPZ.replace(".npz", "_marks.csv").replace("ground_q_vs_class2sd_",
                                                         "control_marks_")
F.drop(columns=["hg"]).to_csv(MARKS_OUT, index=False)
print(f"wrote {MARKS_OUT}  ({len(F)} marks, one row each)")

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
    _lk, _qk = groundq.fit_curve(SD_MM[tr], rk[tr], shape=SHAPE)
    pred[te] = np.interp(ls[te], _lk, _qk)
eiso = err_at(pred)
print(f"  gen2, {len(F)} marks, 5-fold spatially blocked (10 km) CV")
print(f"  {'ground estimator':>34s} {'median err':>11s} {'|median|':>9s} {'RMS':>8s} {'p90|err|':>9s}")
for nm, e in (("q = 0.50 (pipeline default)", e50),
              (f"q = {np.median(rk):.3f} constant (calibrated)", e57),
              (f"q = {SHAPE}(log class-2 SD), held out", eiso)):
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

KNOTS = groundq.fit_curve(SD_MM, rk, shape=SHAPE)
print(f"\n  calibration curve q(class-2 SD), fitted on all {len(F)} marks:")
for v in (20, 40, 60, 80, 120, 200, 400):
    print(f"    SD {v:4d} mm -> q = {float(np.interp(np.log(v), *KNOTS)):.3f}")
# The curve travels with its provenance: a consumer must be able to see what population it
# was fitted on and what it is entitled to be applied to.
groundq.save_curve(
    OUTNPZ, KNOTS, n_marks=len(F), epoch=SET,
    fitted_on=(f"{SET}: point types {'+'.join(sorted(TYPES))}, class-2 returns "
               f"within 7.5 m"),
    point_types="+".join(sorted(TYPES)),
    response="rank of surveyed ground within the mark's class-2 returns",
    covariate="natural log of the class-2 standard deviation, mm",
    shape=("two straight segments in log spread, joined continuously, break SCANNED over "
           "the observed bin centres" if SHAPE == "piecewise"
           else "isotonic, monotone non-increasing -- flat then falling, no break imposed"),
    cv="5-fold spatially blocked on 10 km blocks",
    known_limits=("calibrated on 7.5 m discs and applied to 5 m cells; per-epoch -- "
                  "this curve is valid only for the epoch named in `set`; corrects a "
                  "population, not a cell"))
print(f"\n  wrote {OUTNPZ}")

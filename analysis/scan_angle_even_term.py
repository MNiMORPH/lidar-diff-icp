#!/usr/bin/env python3
"""Measure gen1's EVEN scan-angle elevation term -- the incidence-angle penetration bias.

Andy, 2026-09-20: model it, because it is a SYSTEMATIC and everything else in the budget
is random. Per-cell scatter (~47 mm NMAD) averages down over an area as 47/sqrt(N); an
even scan-angle bias does not. They cross at N ~ 22 cells (~550 m2 at 5 m), so for any
hillslope or catchment the systematic is the error that reaches the answer.

WHAT IT IS. On stable ground the DoD's dependence on scan angle decomposes into an ODD
part (a roll -- a tilt through nadir) and an EVEN part (symmetric about nadir). Measured at
elbaext: RMS odd 2.74 mm, RMS even 2.95 mm, and the even part is a clean Lambda -- gen1
reads progressively HIGHER toward both swath edges, ~10 mm peak-to-trough. That is the
signature of penetration degrading at oblique incidence, not of a geometry error, which is
why no rigid-body registration term touches it.

HOW IT IS MEASURED HERE, and why this way. Fitting f(theta) directly against gen2 is
confounded: within a swath, scan angle is nearly collinear with position, so terrain leaks
in. Instead this uses BETWEEN-LINE differences in SHARED cells. For two lines A and B
covering one cell, the ground cancels exactly:

    d_A - d_B = f(theta_A) - f(theta_B)   + (per-line constant difference)

so regressing on basis DIFFERENCES isolates f's shape, terrain-free. A per-pair intercept
absorbs the per-line constants, and f's own constant is unidentifiable -- which is correct,
because a vertical offset is GAUGE, not error (see level-is-gauge-compare-structure).

TWO BASES, reported side by side, neither preferred by fiat:
  sec(theta) - 1   ONE parameter, physical: path length through near-ground vegetation
                   scales as 1/cos(theta), so a penetration bias should follow it.
  theta^2          one free even parameter, no mechanism. A shape check on the above.
An ODD term (theta) is carried alongside so the roll cannot leak into the even estimate.

VALIDATION IS BUILT IN and is falsifiable: the even term necessarily contributes to the
between-line residual, so removing it must shrink that residual. If it does not, the fit is
describing something else and should be discarded.

    ./lidar-icp/bin/python analysis/scan_angle_even_term.py --tile elbaext_delong
"""
import argparse

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

ap = argparse.ArgumentParser(description=__doc__,
                             formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--tile", required=True)
ap.add_argument("--min-returns", type=int, default=3,
                help="minimum gen1 returns per line per cell before its mean is used. "
                     "MINE, and stated: a cell mean from 1-2 returns is dominated by "
                     "within-cell relief. 3 is the smallest value that is not 1 or 2; the "
                     "sensitivity is printed so the choice is visible, not buried.")
A = ap.parse_args()

D = f"data/derived/{A.tile}"
stable = np.load(f"{D}/stable.npy").astype(bool).ravel()
t = pq.read_table(f"{D}/beam_offset_table.parquet",
                  columns=["cell", "d_mm_corr", "scan_angle", "point_source_id", "in_grid"])
g = t["in_grid"].to_numpy().astype(bool)
cell = t["cell"].to_numpy()[g]
d = t["d_mm_corr"].to_numpy()[g].astype(float)
sa = t["scan_angle"].to_numpy()[g].astype(float)
ps = t["point_source_id"].to_numpy()[g]
ok = stable[cell] & np.isfinite(d) & np.isfinite(sa)
df = pd.DataFrame({"cell": cell[ok], "d": d[ok], "sa": sa[ok], "ps": ps[ok]})
cl = df.groupby(["cell", "ps"]).agg(d=("d", "mean"), sa=("sa", "mean"),
                                    n=("d", "size")).reset_index()
cl = cl[cl.n >= A.min_returns]
print(f"{A.tile}: {len(cl):,} (cell, line) means on stable ground, "
      f">= {A.min_returns} returns each")

lines = sorted(cl.ps.unique())
rows = []
for a_, b_ in zip(lines[:-1], lines[1:]):
    Aa = cl[cl.ps == a_].set_index("cell")
    Bb = cl[cl.ps == b_].set_index("cell")
    j = Aa.join(Bb, how="inner", lsuffix="_a", rsuffix="_b")
    if len(j) < 500:
        continue
    rows.append(pd.DataFrame({"pair": f"{a_}-{b_}", "y": (j.d_a - j.d_b).to_numpy(),
                              "ta": j.sa_a.to_numpy(), "tb": j.sa_b.to_numpy()}))
R = pd.concat(rows, ignore_index=True)
print(f"  {len(R):,} shared-cell line pairs across {R.pair.nunique()} pairs\n")

sec = lambda th: 1.0 / np.cos(np.radians(th)) - 1.0
BASES = {"sec(theta)-1  [physical]": sec, "theta^2       [shape check]": lambda th: th ** 2}
pairs = sorted(R.pair.unique())
P = np.column_stack([(R.pair == p).to_numpy().astype(float) for p in pairs])
odd = (R.ta - R.tb).to_numpy()
y = R.y.to_numpy()

print("FIT: y = per-pair constant + c_odd*(ta-tb) + c_even*(g(ta)-g(tb))")
print(f"  {'basis':<28}{'c_even':>12}{'SE':>9}{'sigma':>8}{'c_odd mm/deg':>14}{'resid sd':>11}")
base_sd = None
fits = {}
for name, gfun in BASES.items():
    X = np.column_stack([P, odd, gfun(R.ta.to_numpy()) - gfun(R.tb.to_numpy())])
    coef, *_ = np.linalg.lstsq(X, y, rcond=None)
    res = y - X @ coef
    s2 = float(res @ res) / (X.shape[0] - X.shape[1])
    cov = s2 * np.linalg.pinv(X.T @ X)
    se = float(np.sqrt(cov[-1, -1]))
    fits[name] = (coef[-1], gfun)
    if base_sd is None:
        base_sd = float(np.std(y - np.column_stack([P, odd]) @
                               np.linalg.lstsq(np.column_stack([P, odd]), y, rcond=None)[0]))
    print(f"  {name:<28}{coef[-1]:>12.3f}{se:>9.3f}{abs(coef[-1])/se:>8.1f}"
          f"{coef[-2]:>14.3f}{float(np.std(res)):>11.2f}")
print(f"  {'(no even term)':<28}{'-':>12}{'-':>9}{'-':>8}{'-':>14}{base_sd:>11.2f}")

print("\nTHE FITTED CURVE, gen1 elevation bias vs scan angle (constant removed -- gauge)")
th = np.array([0, 3, 6, 9, 12, 15, 18], float)
print(f"  {'scan angle':>12}" + "".join(f"{v:>8.0f}" for v in th))
for name, (c, gfun) in fits.items():
    v = c * (gfun(th) - gfun(0.0))
    print(f"  {name.split()[0]:>12}" + "".join(f"{x:>8.2f}" for x in v))
print("  (mm; positive = gen1 reads HIGH, so the DoD reads LOW there)")

print("\nVALIDATION -- removing it must shrink the between-line residual, or it is not real")
for name, (c, gfun) in fits.items():
    corr = y - c * (gfun(R.ta.to_numpy()) - gfun(R.tb.to_numpy()))
    fit0 = np.linalg.lstsq(np.column_stack([P, odd]), corr, rcond=None)[0]
    r = corr - np.column_stack([P, odd]) @ fit0
    print(f"  {name:<28} residual sd {base_sd:>7.2f} -> {float(np.std(r)):>7.2f} mm "
          f"({100*(np.std(r)-base_sd)/base_sd:+.1f}%)")

print("\nSENSITIVITY to --min-returns (the one parameter I chose):")
for mr in (2, 3, 5, 10):
    c2 = cl[cl.n >= mr]
    rr = []
    for a_, b_ in zip(lines[:-1], lines[1:]):
        Aa = c2[c2.ps == a_].set_index("cell"); Bb = c2[c2.ps == b_].set_index("cell")
        j = Aa.join(Bb, how="inner", lsuffix="_a", rsuffix="_b")
        if len(j) >= 500:
            rr.append(pd.DataFrame({"pair": f"{a_}-{b_}", "y": (j.d_a - j.d_b).to_numpy(),
                                    "ta": j.sa_a.to_numpy(), "tb": j.sa_b.to_numpy()}))
    R2 = pd.concat(rr, ignore_index=True)
    P2 = np.column_stack([(R2.pair == p).to_numpy().astype(float)
                          for p in sorted(R2.pair.unique())])
    X2 = np.column_stack([P2, (R2.ta - R2.tb).to_numpy(),
                          sec(R2.ta.to_numpy()) - sec(R2.tb.to_numpy())])
    cf, *_ = np.linalg.lstsq(X2, R2.y.to_numpy(), rcond=None)
    print(f"  min_returns={mr:>3}  n={len(R2):>8,}  c_even(sec) {cf[-1]:>8.3f}")

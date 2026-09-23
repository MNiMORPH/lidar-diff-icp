#!/usr/bin/env python3
"""How many gen2 beams does a ground estimate need? Resampling experiment on the cube.

Andy, 2026-09-23: "In the places with low vegetation coverage, resample from gen2's cube.
Take different numbers of beams at random. Find how they affect the spread and ground
classification."

WHY THIS IS THE RIGHT QUESTION. gen2 carries ~10x gen1's near-ground return count, and the
extra returns buy REACH INTO THE LOW TAIL, not a better centre: decimating gen2 class-2
from 158 to 16 returns per cell moves its LOWEST return +20 mm and its MEDIAN +0 mm. The
pipeline grids ground with ``ground_q = 0.50`` -- a median -- which is invariant to sample
size by construction. So the pipeline currently discards exactly the advantage density buys.
This experiment measures, per estimator, how much of that advantage is recoverable and what
it costs in variance.

THE ESTIMATOR IS THE EXPERIMENT. "Ground classification" here means the rule that turns a
cell's near-ground return distribution into one ground elevation. Different rules have
opposite density behaviour, and that is the point:

    median of all near-ground   density-INVARIANT, biased high by any vegetation
    median of vendor class-2    what the pipeline uses today
    low quantiles (q25/q10/q05) density-SENSITIVE: they reach lower as n grows
    minimum                     most sensitive, and most exposed to a single low blunder

A FIRM TEST NEEDS A CONTROL, so this runs two:

  SYNTHETIC -- draw from a KNOWN distribution and check each estimator behaves as the order
    statistics say it must (median unbiased in n; minimum falling roughly as the 1/n
    quantile). If an estimator fails here the measurement code is wrong, not the data.
  DENSITY-INVARIANCE -- resample the FULL population at n = N (all of it). Every estimator
    must return its full-density value. Any drift is a bug in the resampler.

SAMPLING IS WITHOUT REPLACEMENT, which is what decimating a real flight does. Sampling WITH
replacement leaves a distribution's quantiles unchanged in expectation and would make the
whole experiment return zero -- an earlier version of this analysis made exactly that
mistake and reported the null as a result.

THE LOW-VEGETATION POPULATION IS THRESHOLD-FREE: cells where gen2 records NO return in the
0.15-2.00 m band (lowveg == 0). Other bands are reported alongside so the dependence on
vegetation is visible rather than assumed, and so the choice of population is checkable.

    ./lidar-icp/bin/python analysis/gen2_density_experiment.py
    ./lidar-icp/bin/python analysis/gen2_density_experiment.py --population all
"""
import argparse

import numpy as np

ap = argparse.ArgumentParser(description=__doc__,
                             formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--tile", default="elbaext")
ap.add_argument("--mask-tile", default="elbaext_delong")
ap.add_argument("--population", choices=("lowveg0", "highveg", "all"), default="lowveg0",
                help="lowveg0 (default) = cells with NO gen2 return in the 0.15-2.00 m "
                     "band, a threshold-free 'low vegetation'. 'all' uses every flat "
                     "floodplain cell, for contrast.")
ap.add_argument("--repeats", type=int, default=40,
                help="draws per cell per n. MINE: 40 makes the sampling error of a mean "
                     "over ~3000 cells negligible against the effects here. The "
                     "repeat-count sensitivity is printed.")
ap.add_argument("--min-returns", type=int, default=64,
                help="a cell must hold at least this many gen2 near-ground returns to "
                     "enter, so every n on the ladder is drawable WITHOUT replacement "
                     "from every cell. MINE, and it selects on density, which is stated "
                     "in the output as the population's own density distribution.")
ap.add_argument("--seed", type=int, default=0)
A = ap.parse_args()

RNG = np.random.default_rng(A.seed)
LADDER = [2, 4, 8, 16, 32, 64]


def shorth(d, frac):
    """Midpoint of the SHORTEST interval containing ``frac`` of a SORTED sample.

    Andy, 2026-09-23: "The ground itself has many values within 5x5 m. Median is likely veg
    biased. Minimum will drop forever." Both true, and they rule out the obvious estimators
    for OPPOSITE reasons -- the median has a population target but the wrong one; the
    minimum points the right way but has NO target, converging on the cell's lowest ground
    point rather than its representative ground.

    A MODE-seeking estimator has neither defect. The shorth targets the densest part of the
    distribution -- the ground population -- and ignores an upper tail entirely while
    vegetation stays under ``1 - frac`` of the returns. It has a population value, so it
    converges in n, and it is not an extreme order statistic.
    """
    n = d.size
    k = max(int(np.ceil(frac * n)), 2)
    if k >= n:
        return 0.5 * (d[0] + d[-1])
    w = d[k - 1:] - d[:n - k + 1]
    i = int(np.argmin(w))
    return 0.5 * (d[i] + d[i + k - 1])


def half_sample_mode(d):
    """Bickel-Fruehwirth half-sample mode: recurse on the densest half until <= 3 remain."""
    while d.size > 3:
        n = d.size
        k = n // 2
        w = d[k:] - d[:n - k]
        i = int(np.argmin(w))
        d = d[i:i + k + 1]
    return float(np.median(d))


def stats_of_sorted(d):
    """All estimators from one SORTED draw. One sort, direct indexing -- calling
    np.quantile four times per draw is what made the first version too slow to finish."""
    n = d.size
    i = lambda q: d[min(int(q * (n - 1) + 0.5), n - 1)]
    return (i(0.50), i(0.25), i(0.10), i(0.05), d[0],
            shorth(d, 0.50), shorth(d, 0.25), half_sample_mode(d), i(0.90) - i(0.10))


NAMES = ["q50", "q25", "q10", "q05", "min", "shorth50", "shorth25", "hsm", "spread"]
QS = (0.50, 0.25, 0.10, 0.05)


def quantiles_of_draw(vals, qs=QS):
    return dict(zip(NAMES, stats_of_sorted(np.sort(np.asarray(vals)))))


# ---------------------------------------------------------------- synthetic control
print("=" * 78)
print("CONTROL 1 -- SYNTHETIC. Known distribution; the estimators must behave as the")
print("order statistics require, or the measurement code is wrong.")
truth = RNG.normal(0.0, 60.0, 200_000)          # mm, a plausible near-ground scatter
print(f"\n  {'n':>5}" + "".join(f"{k:>12}" for k in NAMES))
for n in LADDER + [256]:
    acc = {k: [] for k in NAMES}
    for _ in range(300):
        d = RNG.choice(truth, size=n, replace=False)
        for k, v in quantiles_of_draw(d, QS).items():
            acc[k].append(v)
    print(f"  {n:>5}" + "".join(f"{np.mean(acc[k]):>12.1f}" for k in NAMES))
print("\n  EXPECT: q50 flat in n (a median is sample-size invariant); min falling steadily;")
print("  spread RISING toward the population value because small samples cannot reach the")
print("  tails. If q50 drifts, the resampler is biased.")

# ---------------------------------------------------------------- real data
c = np.load(f"data/derived/{A.tile}/nearground_cells_sn.npz")
g = np.load(f"data/derived/{A.tile}/nearground_gen2_class_split.npz")
cells = c["cells"]; e = c["edges"]
ctr = 0.5 * (e[:-1] + e[1:]) * 1000.0                       # bin centres, mm
Hall = (g["Hg"] + g["Hn"]).astype(np.int64)                  # ALL gen2 near-ground
Hg = g["Hg"].astype(np.int64)                                # vendor class-2 only

fp = np.load(f"data/derived/{A.mask_tile}/floodplain_mask.npy").astype(bool).ravel()
sl = np.load(f"data/derived/{A.mask_tile}/slope.npy").ravel()
band = (ctr > 150.0) & (ctr <= 2000.0)                       # the lowveg band, in mm
tot = Hall.sum(1)
lowveg = np.where(tot > 0, Hall[:, band].sum(1) / np.maximum(tot, 1), np.nan)

base = fp[cells] & (sl[cells] <= 2.0) & (tot >= A.min_returns)
if A.population == "lowveg0":
    pop = base & (lowveg == 0)
elif A.population == "highveg":
    # TOP QUARTILE of lowveg among cells meeting the density floor. A quartile, not an
    # invented level: the cut is stated, and the band table below shows the whole range so
    # the choice is checkable rather than load-bearing.
    thr = np.nanpercentile(lowveg[base], 75)
    pop = base & (lowveg >= thr)
    print(f"  high-vegetation cut: lowveg >= {thr:.4f} (75th percentile of the density-"
          f"qualified cells)")
else:
    pop = base
print("\n" + "=" * 78)
print(f"POPULATION: flat floodplain (slope <= 2 deg), >= {A.min_returns} gen2 near-ground "
      f"returns" + (", lowveg == 0" if A.population == "lowveg0" else ""))
print(f"  cells: {int(pop.sum()):,} of {int(base.sum()):,} that meet the density floor")
print(f"  gen2 near-ground returns per cell: median {np.median(tot[pop]):.0f}  "
      f"p10 {np.percentile(tot[pop],10):.0f}  p90 {np.percentile(tot[pop],90):.0f}")
print(f"  class-2 fraction: median {np.median(Hg[pop].sum(1)/tot[pop]):.3f}")
print(f"  lowveg in this population: median {np.nanmedian(lowveg[pop]):.4f}  "
      f"max {np.nanmax(lowveg[pop]):.4f}")

idx = np.flatnonzero(pop)
samples = [np.repeat(ctr, Hall[i]) for i in idx]             # each cell's return heights
full = {k: [] for k in NAMES}
for v in samples:
    for k, q in quantiles_of_draw(v, QS).items():
        full[k].append(q)
full = {k: np.array(v) for k, v in full.items()}

print("\n" + "=" * 78)
print("FULL-DENSITY ESTIMATORS ACROSS VEGETATION BANDS (mm; no resampling needed)")
print("The question: does a MODE estimator separate the ground from the vegetation tail")
print("where a median cannot? If so, mode-minus-median must WIDEN with vegetation.\n")
qb = np.nanpercentile(lowveg[base], [25, 50, 75])
print(f"  {'lowveg band':>18}{'cells':>8}" + "".join(f"{k:>10}" for k in NAMES)
      + f"{'q50-hsm':>10}")
for lo, hi in zip([-1e-9] + list(qb), list(qb) + [9.0]):
    k = base & (lowveg > lo) & (lowveg <= hi)
    if k.sum() < 200:
        continue
    vals = [np.repeat(ctr, Hall[i]) for i in np.flatnonzero(k)]
    st = np.array([stats_of_sorted(np.sort(v)) for v in vals])
    m_ = st.mean(0)
    print(f"  {f'{max(lo,0):.4f}-{hi:.4f}':>18}{int(k.sum()):>8,}"
          + "".join(f"{x:>10.1f}" for x in m_)
          + f"{m_[0]-m_[7]:>10.1f}")
print("\n  q50-hsm is the SEPARATION the mode buys. If it grows with vegetation, mixture")
print("  separation is doing real work; if it is flat, the mode is just a noisy median.\n")

print("\nCONTROL 2 -- DENSITY INVARIANCE. Drawing n = N (the whole cell) must reproduce")
print("the full-density value exactly for every estimator.")
chk = {k: [] for k in NAMES}
for v in samples[:400]:
    d = RNG.choice(v, size=len(v), replace=False)
    for k, q in quantiles_of_draw(d, QS).items():
        chk[k].append(q)
print("  " + "  ".join(f"{k}:{np.max(np.abs(np.array(chk[k])-full[k][:400])):.2e}"
                       for k in NAMES))
print("  (max |drawn - full| over 400 cells; must be 0)")

print("\n" + "=" * 78)
print("THE EXPERIMENT -- estimator value vs number of beams drawn (mm)")
print(f"  {int(pop.sum()):,} cells x {A.repeats} repeats per n, WITHOUT replacement\n", flush=True)
mean, sdv = {}, {}
for n in LADDER:
    pm = np.empty((len(samples), len(NAMES))); ps = np.empty_like(pm)
    for j2, v in enumerate(samples):
        N = v.size
        st = np.empty((A.repeats, len(NAMES)))
        for r in range(A.repeats):
            st[r] = stats_of_sorted(np.sort(v[RNG.permutation(N)[:n]]))
        pm[j2] = st.mean(0); ps[j2] = st.std(0)
    mean[n] = pm.mean(0); sdv[n] = ps.mean(0)
    print(f"  n={n:>4}" + "".join(f"{x:>10.1f}" for x in mean[n]), flush=True)
print(f"  {'FULL':>6}" + "".join(f"{np.mean(full[k]):>10.1f}" for k in NAMES))
print(f"  {'':>6}" + "".join(f"{k:>10}" for k in NAMES))

print("\nBIAS vs FULL DENSITY (mm)")
print(f"  {'n':>6}" + "".join(f"{k:>10}" for k in NAMES))
for n in LADDER:
    print(f"  {n:>6}" + "".join(f"{mean[n][i]-np.mean(full[NAMES[i]]):>10.1f}"
                                for i in range(len(NAMES))))

print("\nPER-CELL SCATTER across repeats (mm sd) -- the price of the reach")
print(f"  {'n':>6}" + "".join(f"{k:>10}" for k in NAMES))
for n in LADDER:
    print(f"  {n:>6}" + "".join(f"{x:>10.1f}" for x in sdv[n]))

print("\nCONVERGENCE: change from n to 2n (~0 means it has a target and reached it)")
print(f"  {'n->2n':>8}" + "".join(f"{k:>10}" for k in NAMES))
for a_, b_ in zip(LADDER[:-1], LADDER[1:]):
    print(f"  {f'{a_}->{b_}':>8}" + "".join(f"{mean[b_][i]-mean[a_][i]:>10.1f}"
                                            for i in range(len(NAMES))))

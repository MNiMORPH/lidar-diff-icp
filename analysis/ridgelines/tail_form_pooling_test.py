#!/usr/bin/env python3
"""Is the open-ground POWER-LAW tail real, or manufactured by pooling cells?

analysis/ridgelines/ground_mixture_fit.py measured the near-ground tail form on histograms
POOLED per land-cover stratum: exponential in forest, power-law on open/bare. A mixture of
exponentials with DIFFERENT rates is heavier-tailed than any of its parts -- with the rates
Gamma-distributed it is exactly Pareto -- so a power law can appear in a pooled sample when
no individual cell has one. Open ground is where the rate should vary most between cells
(bare soil, stubble rows, clods, tussocks); forest cells should be more alike.

THE TEST. Estimate each cell's own exponential rate, group cells by that rate, and pool only
WITHIN a group. If pooling was the cause, the power-law signature must weaken as the group
narrows, while the exponential fit holds. If the tail is genuinely algebraic, stratifying
changes nothing.

Origin is each cell's OWN modal bin -- the ground peak -- not a chosen height, so no cutoff
enters the tail definition. The one discretionary number is how many tail returns a cell
needs before its rate is estimable; results are printed for several so the choice is visible
rather than assumed.

Form is judged by which straight line fits better, on the same pooled counts:
    exponential   log p   vs   d          linear
    power law     log p   vs   log d      linear

    ./lidar-icp/bin/python analysis/ridgelines/tail_form_pooling_test.py --tile data/derived/elba
"""
import argparse
import os

import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("--tile", default="data/derived/elba")
ap.add_argument("--cube", default="nearground_cells_sn.npz")
# WHICH RETURNS. ground_mixture_fit.py measured the tail of the CLASS-2 (ground-class)
# distribution, so that is the default here: a test on all returns would be measuring the
# canopy above the mode, which is a different distribution and not the documented claim.
ap.add_argument("--returns", choices=["class2", "all"], default="class2")
ap.add_argument("--groups", type=int, default=5, help="rate strata, by quantile of lambda")
A = ap.parse_args()

D = A.tile
z = np.load(os.path.join(D, A.cube))
cells, edges = z["cells"], z["edges"]
if A.returns == "class2":
    H = np.load(os.path.join(D, "nearground_gen2_class_split.npz"))["Hg"].astype(float)
else:
    H = z["H2"].astype(float)
ctr = 0.5 * (edges[:-1] + edges[1:])
dz = float(edges[1] - edges[0])

cover = np.load(os.path.join(D, "canopy_cover_pfs.npy")).ravel()[cells]
# Repo conventions, declared not calibrated (analysis/forest_metrics_pfs.py).
FOREST, OPEN = 0.5, 0.1
strata = {"open (cover <= 0.1)": np.isfinite(cover) & (cover <= OPEN),
          "forest (cover >= 0.5)": np.isfinite(cover) & (cover >= FOREST)}


def line_r2(x, y):
    """R^2 of a straight-line fit, the comparison both forms are judged by."""
    ok = np.isfinite(x) & np.isfinite(y)
    if ok.sum() < 3:
        return np.nan
    x, y = x[ok], y[ok]
    p = np.polyfit(x, y, 1)
    r = y - np.polyval(p, x)
    ss = np.sum((y - y.mean()) ** 2)
    return float(1.0 - np.sum(r ** 2) / ss) if ss > 0 else np.nan


def forms(counts, d):
    """(exponential R^2, power-law R^2) for a pooled tail: counts at heights d above mode."""
    m = counts > 0
    if m.sum() < 4:
        return np.nan, np.nan
    lp = np.log(counts[m])
    return line_r2(d[m], lp), line_r2(np.log(d[m]), lp)


def cell_rates(sel, min_tail):
    """Per-cell exponential rate above the cell's own modal bin, and the stacked tails."""
    Hs = H[sel]
    mode = Hs.argmax(1)
    nb = Hs.shape[1]
    off = np.arange(nb)[None, :] - mode[:, None]        # bins above each cell's own mode
    tail_mask = off > 0
    n_tail = np.where(tail_mask, Hs, 0).sum(1)
    keep = n_tail >= min_tail
    # lambda from the MLE for an exponential: 1 / mean height above the mode
    mean_h = np.where(keep, (np.where(tail_mask, Hs * off * dz, 0).sum(1)
                             / np.maximum(n_tail, 1)), np.nan)
    lam = np.where(mean_h > 0, 1.0 / mean_h, np.nan)
    return keep, lam, off, tail_mask, Hs


MAXOFF = 100                      # stack tails to 1 m above the mode (50 bins at 2 cm)
print(f"tile {os.path.basename(D)}   returns={A.returns}   bins {dz*100:.0f} cm   tails stacked to "
      f"{MAXOFF*dz:.2f} m above each cell's own mode")
for label, sel in strata.items():
    print(f"\n{label}: {int(sel.sum()):,} cells")
    for min_tail in (10, 30, 100):
        keep, lam, off, tmask, Hs = cell_rates(sel, min_tail)
        if keep.sum() < 20:
            print(f"  min_tail={min_tail:4d}: only {int(keep.sum())} cells qualify; skipped")
            continue
        lk = lam[keep]
        # ALL cells of this stratum pooled -- the measurement that gave the power law
        acc = np.zeros(MAXOFF)
        idx = np.clip(off, 0, MAXOFF)
        for j in range(1, MAXOFF):
            acc[j] = Hs[keep][(idx[keep] == j) & tmask[keep]].sum()
        d = np.arange(MAXOFF) * dz
        e_all, p_all = forms(acc[1:], d[1:])
        q = np.nanpercentile(lk, np.linspace(0, 100, A.groups + 1))
        print(f"  min_tail={min_tail:4d}   cells {int(keep.sum()):6,d}   "
              f"lambda 1/m: p10 {np.nanpercentile(lk,10):6.2f}  median "
              f"{np.nanmedian(lk):6.2f}  p90 {np.nanpercentile(lk,90):6.2f}  "
              f"spread p90/p10 {np.nanpercentile(lk,90)/max(np.nanpercentile(lk,10),1e-9):5.2f}x")
        print(f"      POOLED all cells      exp R2 {e_all:6.4f}   power R2 {p_all:6.4f}"
              f"   {'power' if p_all > e_all else 'exp'} wins")
        for g in range(A.groups):
            sub = keep.copy()
            sub[keep] = (lk >= q[g]) & (lk <= q[g + 1])
            if sub.sum() < 10:
                continue
            acc = np.zeros(MAXOFF)
            for j in range(1, MAXOFF):
                acc[j] = Hs[sub][(idx[sub] == j) & tmask[sub]].sum()
            e_g, p_g = forms(acc[1:], d[1:])
            print(f"      lambda group {g+1}/{A.groups} [{q[g]:6.2f},{q[g+1]:6.2f}]  "
                  f"cells {int(sub.sum()):6,d}   exp R2 {e_g:6.4f}   power R2 {p_g:6.4f}"
                  f"   {'power' if p_g > e_g else 'exp'} wins")

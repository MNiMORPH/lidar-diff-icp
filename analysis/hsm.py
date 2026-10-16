"""Half-sample mode (Bickel & Fruehwirth 2006), on raw values and on a binned histogram.

THE TIE DEFECT THIS EXISTS TO FIX (found 2026-10-16 by the control-mark null check,
``control_mode_null_check.py``). Reading the mode off a histogram means recursing on a
sample full of REPEATED bin centres. Many candidate half-windows then tie at the minimum
width, and ``np.argmin`` returns the FIRST of them -- the lowest-starting window -- every
time, at every level of the recursion. The result is a mode biased LOW by an amount that
GROWS WITH THE SPREAD of the distribution (synthetic, 20 mm bins, 400 draws each):

    sigma  n     first-of-ties   middle-of-ties   raw values
      60   150        -29.5             +1.1          -0.2
      60    16        -15.7             +7.5          +0.8
     150   150        -58.1             +5.8          -1.4
      30  3000        -18.5             -0.5          -0.3

A spread-dependent downward bias is exactly the signature the floodplain work was looking
for in the data, so every result that used the first-of-ties form on binned values has to
be re-run with this one. Taking the MIDDLE of the tied windows removes the bias at n >= 150;
a residual of order +7 mm remains at n = 16, so prefer ``half_sample_mode`` on the raw
values whenever they are available.
"""
import numpy as np


def half_sample_mode(d):
    """Half-sample mode of raw values. Sorts; recurses on the densest half until <= 3 remain."""
    d = np.sort(np.asarray(d, float))
    while d.size > 3:
        n = d.size
        k = n // 2
        w = d[k:] - d[:n - k]
        ties = np.flatnonzero(w == w.min())
        i = int(ties[ties.size // 2])
        d = d[i:i + k + 1]
    return float(np.median(d)) if d.size else np.nan


def hsm_from_hist(counts, ctr):
    """Half-sample mode read off a histogram: ``counts`` per bin with centres ``ctr``.

    Ties among equal-width windows are broken at the MIDDLE tied window, not the first.
    """
    v = np.repeat(ctr, np.asarray(counts).astype(int))
    if v.size == 0:
        return np.nan
    while v.size > 3:
        k = v.size // 2
        w = v[k:] - v[:v.size - k]
        ties = np.flatnonzero(w == w.min())
        a = int(ties[ties.size // 2])
        v = v[a:a + k + 1]
    return float(np.median(v))

#!/usr/bin/env python3
"""Does the modal-surface Gaussian land on the surveyed ground -- and does occlusion move it?

THE TEST (Andy, 2026-09-03). Ground returns are missing wherever a pulse did not reach the
surface. If they are missing AT RANDOM with respect to ground height within the cell, the
survivors are still an unbiased sample of the ground Gaussian, mu_g needs no shift, and the
component above it is a nuisance to exclude rather than a model to get right. If instead
occlusion is correlated with local ground height -- vegetation collecting in furrows, or a
bump shadowing the ground behind it for an off-nadir beam -- mu_g is biased, and the bias
should grow as the ground component thins.

So: fit the mixture at each surveyed mark, and regress (mu_g - surveyed) on the ground
weight w_g. No correlation => missing-at-random holds and Andy's conclusion is right. A
correlation => that correlation IS the correction, measured against truth.

MODEL, on ALL returns (not class-2: a classifier has already cut that one, and it keeps
misclassified low vegetation -- our class-2 surface sits ~62 mm high in vegetation):

    counts_j ~ Poisson( N * [ w * Normal(d_j; mu, sigma)
                              + (1-w) * ExpNormal(d_j; mu, sigma, tau) ] * dz )

Component 1 is the modal surface: range error plus within-cell roughness. Component 2 is
everything above it -- plants in forest, microtopography on bare ground -- as an exponential
arrival density convolved with the same measurement error, which is why it is an
exponentially-modified Gaussian sharing sigma rather than a free second peak.

Poisson likelihood because bin counts are Poisson: least-squares would assume constant
variance and log-least-squares is undefined at the empty bins that carry the tail shape.

TWO LIMITS OF THE STORED DATA, stated because they bound what this can conclude:
  * bins are 0.02 m. The instrumental floor is sigma ~ 0.011 m (instrumental-precision note),
    so a truly clean ground peak is NARROWER THAN ONE BIN and sigma_g is resolved only where
    roughness or a surface mat widens it.
  * the boxes store returns, not pulses, so w_g here is the ground component's share OF
    RETURNS IN THE WINDOW, not a per-pulse detection probability. It still orders marks by
    how thin the ground component is, which is what the test needs.

    ./lidar-icp/bin/python analysis/control_mixture_validation.py
"""
import os
import sys

import numpy as np
from scipy.ndimage import gaussian_filter1d
from scipy.stats import exponnorm, norm, spearmanr

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from control_lowveg_offset import STRUCT, load

SET = "gen2_2021_control"


MAXOFF_TAU = 3.0        # tau cannot exceed the window it is fitted in

def fit_mark(point_id, surveyed_z, iters=400, tol=1e-9):
    """Generalised EM for the two-component near-ground mixture on binned counts.

    Nelder-Mead on all four parameters at once was the first attempt and it FAILED: the
    likelihood is nearly flat in the weight, so 155 of 389 marks came back at exactly the
    0.500 the optimiser started from, and tau ran past the window edge (p90 2404.5 mm
    against a +2.00 m top). Those were not fits.

    EM does not have that failure mode, because every update is closed form given the
    responsibilities:

        E   r_j = w N(d_j) / [ w N(d_j) + (1-w) EMG(d_j) ]        per bin
        M   w   = sum(r_j c_j) / sum(c_j)                          exact
            mu, sigma = r-weighted Gaussian MLE                    exact
            tau = (1-r)-weighted mean of (d - mu)                  exact, EMG mean = mu + tau

    The weight is never optimised numerically, which is precisely the parameter the
    general-purpose optimiser could not move.
    """
    f = os.path.join(STRUCT, f"{SET}__{point_id}.npz")
    if not os.path.exists(f):
        return None
    z = np.load(f)
    coef = z["surface_coef"]
    nn = np.sqrt(1.0 + coef[1] ** 2 + coef[2] ** 2)
    h_s = (float(surveyed_z) - float(coef[0])) / nn
    e = z["ng_edges"]
    c = z["ng_all"].astype(float)
    n = c.sum()
    if n < 50:
        return None
    d = 0.5 * (e[:-1] + e[1:])

    # Initialise from the data, not from a constant. The LOWER flank is essentially pure
    # ground (nothing lies below the surface), so its spread about the mode is an honest
    # first sigma -- doubled because it is a half-width.
    mu = d[int(np.argmax(c))]
    lo = d <= mu
    sig = max(np.sqrt(np.sum(c[lo] * (d[lo] - mu) ** 2) / max(c[lo].sum(), 1.0)), 0.01)
    hi = d > mu
    tau = max(float(np.sum(c[hi] * (d[hi] - mu)) / max(c[hi].sum(), 1.0)), 0.01)
    w = 0.5
    prev = None
    for it in range(iters):
        g = norm.pdf(d, mu, sig)
        with np.errstate(over="ignore", invalid="ignore"):
            v = exponnorm.pdf(d, max(tau / sig, 1e-6), loc=mu, scale=sig)
        v = np.where(np.isfinite(v), v, 0.0)
        den = w * g + (1 - w) * v
        ok = den > 0
        r = np.zeros_like(d)
        r[ok] = w * g[ok] / den[ok]
        cr = c * r
        cw = c * (1 - r)
        sw = cr.sum()
        w = float(sw / max(c.sum(), 1.0))
        if sw > 0:
            mu = float(np.sum(cr * d) / sw)
            sig = float(max(np.sqrt(np.sum(cr * (d - mu) ** 2) / sw), 0.005))
        if cw.sum() > 0:
            tau = float(np.clip(np.sum(cw * (d - mu)) / cw.sum(), 0.005, MAXOFF_TAU))
        ll = float(np.sum(c[ok] * np.log(den[ok])))
        if prev is not None and abs(ll - prev) < tol * max(abs(prev), 1.0):
            break
        prev = ll
    return dict(point_id=point_id, mu_g=mu, sigma_g=sig, tau=tau, w_g=w, n=n,
                h_survey=h_s, resid_mm=(mu - h_s) * 1000.0, iters=it + 1)


def fit_mark_scalesep(point_id, surveyed_z, bw=0.15, iters=200, tol=1e-9):
    """Ground = the only SHARP feature; vegetation = anything SMOOTH. No form assumed above.

    The exponential was the wrong object. Its rate is a mean free path, 1/tau = G a / cos
    theta, which is a density measurement -- but only if the leaf area density a is CONSTANT
    over the fitted range, and it is not: herb layer, shrubs, trunk space, crown. The
    measured profile shows it, staying near-flat from +50 to +500 mm where an exponential
    decaying from the ground would have fallen away.

    So leave a(h) unknown. The two components are separated by SCALE, not by family:

        ground      the only narrow feature -- instrumental sigma ~11 mm plus within-cell
                    roughness
        vegetation  smooth in height, varying over hundreds of mm

    Iterate as EM, but the M-step for the vegetation component is a SMOOTH of its own
    responsibility-weighted counts rather than a parameter update. The smoother's bandwidth
    is what makes the separation: wide enough that the vegetation component CANNOT represent
    a ground-width peak.

    bw is MINE, and it is the one discretionary number here: the reported runs sweep it.
    """
    f = os.path.join(STRUCT, f"{SET}__{point_id}.npz")
    if not os.path.exists(f):
        return None
    z = np.load(f)
    coef = z["surface_coef"]
    nn = np.sqrt(1.0 + coef[1] ** 2 + coef[2] ** 2)
    h_s = (float(surveyed_z) - float(coef[0])) / nn
    e = z["ng_edges"]; c = z["ng_all"].astype(float)
    if c.sum() < 50:
        return None
    d = 0.5 * (e[:-1] + e[1:]); dz = float(e[1] - e[0])
    mu = d[int(np.argmax(c))]
    lo = d <= mu
    sig = max(np.sqrt(np.sum(c[lo] * (d[lo] - mu) ** 2) / max(c[lo].sum(), 1.0)), 0.01)
    w = 0.5
    r = np.full(d.size, 0.5)
    prev = None
    for it in range(iters):
        veg = gaussian_filter1d(np.maximum(c * (1 - r), 0.0), bw / dz, mode="nearest")
        veg = veg / max(veg.sum() * dz, 1e-12)
        g = norm.pdf(d, mu, sig)
        den = w * g + (1 - w) * veg
        ok = den > 0
        r = np.zeros_like(d); r[ok] = w * g[ok] / den[ok]
        cr = c * r; sw = cr.sum()
        w = float(sw / max(c.sum(), 1.0))
        if sw > 0:
            mu = float(np.sum(cr * d) / sw)
            sig = float(max(np.sqrt(np.sum(cr * (d - mu) ** 2) / sw), 0.005))
        ll = float(np.sum(c[ok] * np.log(den[ok])))
        if prev is not None and abs(ll - prev) < tol * max(abs(prev), 1.0):
            break
        prev = ll
    return dict(point_id=point_id, mu_g=mu, sigma_g=sig, tau=np.nan, w_g=w, n=c.sum(),
                h_survey=h_s, resid_mm=(mu - h_s) * 1000.0, iters=it + 1)


m = load(0.15, 2.00)
MODEL = os.environ.get("MODEL", "emg")
BW = float(os.environ.get("BW", "0.15"))
if MODEL == "scalesep":
    print(f"model: narrow Gaussian + UNKNOWN smooth (bandwidth {BW*1000:.0f} mm)")
    rows = [r for r in (fit_mark_scalesep(t.point_id, t.elevation, bw=BW)
                        for t in m.itertuples()) if r]
else:
    print("model: Gaussian + EMG (superseded: assumes constant leaf density)")
    rows = [r for r in (fit_mark(t.point_id, t.elevation) for t in m.itertuples()) if r]
import pandas as pd
F = pd.DataFrame(rows).merge(m[["point_id", "lowveg"]], on="point_id", how="left")
print(f"marks fitted: {len(F):,} of {len(m):,}   EM iterations: median {F.iters.median():.0f}"
      f"  max {F.iters.max():.0f}")
_pin = int(((F.w_g - 0.5).abs() < 5e-3).sum())
print(f"  w_g within 0.005 of the 0.500 START: {_pin}  "
      f"(Nelder-Mead left 155 of 389 there; EM should not)")
print(f"  sigma_g  median {F.sigma_g.median()*1000:6.1f} mm   p90 {F.sigma_g.quantile(.9)*1000:6.1f}"
      f"   (instrumental floor ~11 mm; excess = roughness or a surface mat)")
print(f"  tau      median {F.tau.median()*1000:6.1f} mm   p90 {F.tau.quantile(.9)*1000:6.1f}")
print(f"  w_g      median {F.w_g.median():6.3f}   p10 {F.w_g.quantile(.1):6.3f}")
print(f"  mu_g - surveyed:  median {F.resid_mm.median():+7.1f} mm   "
      f"mean {F.resid_mm.mean():+7.1f}   sd {F.resid_mm.std():6.1f}")

print("\nTHE TEST -- does the residual move with the ground weight?")
ok = np.isfinite(F.resid_mm) & np.isfinite(F.w_g)
rho = spearmanr(F.w_g[ok], F.resid_mm[ok])
# A null result is worth nothing without a bound on what it could have hidden, so the
# slope carries its standard error and a 95% interval.
_x = F.w_g[ok].to_numpy(); _y = F.resid_mm[ok].to_numpy()
_A = np.vstack([np.ones_like(_x), _x]).T
_p, _res, *_ = np.linalg.lstsq(_A, _y, rcond=None)
_dof = len(_x) - 2
_s2 = float(np.sum((_y - _A @ _p) ** 2) / _dof)
_se = np.sqrt(np.diag(np.linalg.inv(_A.T @ _A) * _s2))
print(f"  Spearman rho(w_g, resid) = {rho.statistic:+.4f}   p = {rho.pvalue:.3e}   n = {int(ok.sum())}")
print(f"  linear: resid_mm = {_p[0]:+.1f} {_p[1]:+.1f} * w_g   "
      f"slope SE {_se[1]:.1f}   95% CI [{_p[1]-1.96*_se[1]:+.1f}, {_p[1]+1.96*_se[1]:+.1f}] mm per unit")
print(f"  {'':4s}{'w_g bin':>14s} {'n':>5s} {'mean w_g':>9s} {'median resid mm':>16s} {'SE':>7s}")
q = np.nanpercentile(F.w_g[ok], np.linspace(0, 100, 6))
for i in range(5):
    s = ok & (F.w_g >= q[i]) & (F.w_g <= q[i + 1])
    if s.sum() < 3:
        continue
    se = F.resid_mm[s].std(ddof=1) / np.sqrt(s.sum())
    print(f"    {q[i]:6.3f}-{q[i+1]:<7.3f} {int(s.sum()):5d} {F.w_g[s].mean():9.3f} "
          f"{F.resid_mm[s].median():16.1f} {se:7.1f}")
print("\n  no trend  => missing-at-random holds; mu_g needs no occlusion shift")
print("  a trend   => that slope IS the correction, measured against surveyed ground")

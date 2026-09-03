#!/usr/bin/env python3
"""Ground from the MODE of the return column, and its measured shift from surveyed truth.

The method this replaces, and why. The cover correction has been a PERCENTILE: pick the
gen2 quantile whose elevation matches gen1's median ground, as a function of cover. A
percentile is a MASS statistic, so it moves whenever the component weights move -- and the
weights are exactly what vegetation changes. Measured at these marks, the ground component's
weight runs from 0.922 on bare ground to 0.000 under heavy cover, and the percentile route
needed corrections of hundreds of mm and drove q2 outside [0,1] on 8-11% of tile cells.

The MODE is a LOCATION statistic. It does not move with the weights at all, as long as the
dominant peak stays the same peak. Measured against surveyed ground at these marks it sits
at +10 mm -- one bin, indistinguishable from zero -- for 353 of 389 marks, stepping to
+70 mm only where cover is heavy enough that the dominant peak stops being the ground and
becomes the vegetation mat.

WHAT THE MAT IS (measured here, referenced to surveyed ground, not to a fitted mode):
mode 69.8 mm above true ground on the barest marks rising to 378.7 mm under heavy cover,
spread 78.8 -> 293.8 mm, and SKEW NEAR ZERO throughout (0.09 overall). It is a symmetric,
compact hump -- a Gaussian-like peak in its own right, better sampled than the broad, sparse
ground beneath it. That is why a free two-component fit prefers it: not a malfunction, but
the correct identification of the wrong peak.

NO WINDOWS (Andy, 2026-09-03). Every return in the mark's box enters, out to whatever height
the cloud reaches -- h_max hits 12.88 m at p90. The -1..+2 m window used earlier was
inherited from the archived histograms, never justified, and made the ground weight
denominator-dependent: forest marks had their canopy excluded, inflating the ground share
(w_g median 0.849 windowed against 0.564 unwindowed).

THE ONE SCALE THAT REMAINS, and it is a window in disguise: locating a mode requires
smoothing. So the shift is reported at SEVERAL smoothing scales rather than one, and the
sensitivity is the answer, not a footnote.

WITHIN-EPOCH ONLY. The shift is (mode - surveyed) for one epoch against its own control,
which are on the same vertical datum -- gen1 against NAVD88(GEOID03) 2008 control, gen2
against its own. Comparing the two epochs' shifts is a further step and needs the geoid
difference handled, which this does not do.

    ./lidar-icp/bin/python analysis/control_mode_shift.py --set gen2_2021_control
    ./lidar-icp/bin/python analysis/control_mode_shift.py --set gen1_2008_control
"""
import argparse
import os
import sys

import laspy
import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter1d
from scipy.stats import norm

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from control_lowveg_offset import STRUCT, lowveg          # same metric, same edges

from lidar_diff_icp.groundtruth.tie import _design

BOX = "data/derived/control_boxes"
CONTROL = {"gen1_2008_control": "src/lidar_diff_icp/groundtruth/data/mn_dnr_2008_control_semn.csv",
           "gen2_2021_control": "src/lidar_diff_icp/groundtruth/data/mn_se_driftless_2021_control.csv"}

ap = argparse.ArgumentParser()
ap.add_argument("--set", dest="set_", choices=sorted(CONTROL), default="gen2_2021_control")
ap.add_argument("--smooth", default="0.005,0.01,0.02,0.05,0.10",
                help="mode-finding smoothing scales in m, comma separated. MINE, and the "
                     "reason they are swept rather than chosen: locating a mode needs a "
                     "scale, which is a window in disguise.")
ap.add_argument("--bw", type=float, default=0.15,
                help="smoother bandwidth for the mixture's vegetation component, m. MINE.")
ap.add_argument("--dz", type=float, default=0.005,
                help="histogram bin, m. Far below the ~50 mm ground-component width measured "
                     "at these marks, so it is not the limiting scale.")
ap.add_argument("--out", default=None)
A = ap.parse_args()
SCALES = [float(s) for s in A.smooth.split(",")]
OUT = A.out or f"data/derived/control_mode_shift_{A.set_}.csv"


def marks(set_):
    """point_id, easting, northing, surveyed elevation -- from that epoch's OWN control."""
    c = pd.read_csv(CONTROL[set_]).drop_duplicates(subset=["easting", "northing", "elevation"])
    need = {"point_id", "easting", "northing", "elevation"}
    if not need <= set(c.columns):
        raise SystemExit(f"{CONTROL[set_]} lacks {sorted(need - set(c.columns))}")
    return c[["point_id", "easting", "northing", "elevation"]].dropna()


def column(pid, set_):
    """Every return in the mark's box, as slope-normal height above the fitted surface.

    The reference surface and its coefficients are the STORED ones, so this height is in the
    same frame as everything else built from these marks -- not a second estimator.
    """
    sp, bp = f"{STRUCT}/{set_}__{pid}.npz", f"{BOX}/{set_}__{pid}.laz"
    if not (os.path.exists(sp) and os.path.exists(bp)):
        return None
    z = np.load(sp)
    coef = z["surface_coef"]
    E, N, R = float(z["easting"]), float(z["northing"]), float(z["struct_radius"])
    f = laspy.read(bp)
    x, y, zz = np.asarray(f.x), np.asarray(f.y), np.asarray(f.z)
    sel = np.hypot(x - E, y - N) <= R
    if sel.sum() < 50:
        return None
    nn = np.sqrt(1.0 + coef[1] ** 2 + coef[2] ** 2)
    h = (zz[sel] - (_design(x[sel] - E, y[sel] - N, 2) @ coef)) / nn
    return h, float(coef[0]), nn


def fit_one(pid, surveyed_z, set_):
    got = column(pid, set_)
    if got is None:
        return None
    h, c0, nn = got
    mu_true = (float(surveyed_z) - c0) / nn        # the surveyed ground, in this frame
    e = np.arange(np.floor(h.min() / A.dz) * A.dz, np.ceil(h.max() / A.dz) * A.dz + A.dz, A.dz)
    c = np.histogram(h, bins=e)[0].astype(float)
    d = 0.5 * (e[:-1] + e[1:])
    out = dict(point_id=pid, n_all=float(h.size), h_min=float(h.min()), h_max=float(h.max()),
               mu_true=mu_true)

    # 1. THE MODE SHIFT, at every smoothing scale. No model, no window: just where the
    #    dominant peak of the observed column sits relative to surveyed ground.
    for s in SCALES:
        k = max(s / A.dz, 1e-6)
        cs = gaussian_filter1d(c, k)
        i = int(np.argmax(cs))
        out[f"mode_shift_{s:g}"] = float(d[i] - mu_true)
        # WIDTH of the dominant peak, from the lidar ALONE -- no surveyed elevation enters,
        # so it can predict the shift rather than encode it. Full width at half maximum:
        # walk out from the peak to the first bin below half its height on each side. Half
        # maximum is a definition, not a tuned parameter; if the peak never falls to half
        # within the column the width is not defined and is left NaN rather than capped.
        half = 0.5 * cs[i]
        li = np.where(cs[:i] < half)[0]
        ri = np.where(cs[i + 1:] < half)[0]
        out[f"fwhm_{s:g}"] = (float(d[i + 1 + ri[0]] - d[li[-1]])
                              if li.size and ri.size else np.nan)
        # Robust spread about the mode, again lidar-only: the median absolute deviation of
        # the returns from the mode, scaled to a Gaussian sigma.
        out[f"nmad_{s:g}"] = float(1.4826 * np.median(np.abs(h - d[i])))

    # 2. The pinned mixture: mu held AT the surveyed ground, so sigma_g and w_g are
    #    measurements of the real ground return rather than of wherever a fit wandered.
    sig, w, r = 0.02, 0.5, np.full(d.size, 0.5)
    for _ in range(150):
        veg = gaussian_filter1d(np.maximum(c * (1 - r), 0.0), A.bw / A.dz, mode="nearest")
        veg = veg / max(veg.sum() * A.dz, 1e-12)
        g = norm.pdf(d, mu_true, sig)
        den = w * g + (1 - w) * veg
        ok = den > 0
        r = np.zeros_like(d)
        r[ok] = w * g[ok] / den[ok]
        cr = c * r
        sw = cr.sum()
        w = float(sw / max(c.sum(), 1.0))
        if sw > 0:
            sig = float(max(np.sqrt(np.sum(cr * (d - mu_true) ** 2) / sw), 0.002))
    out["w_g"] = w
    out["sigma_g"] = sig

    # 3. The mat: the above-ground component, described against SURVEYED ground.
    matc = c * (1 - r)
    tot = matc.sum()
    if tot > 0:
        hh = d - mu_true
        m2 = matc / tot
        mean = float(np.sum(m2 * hh))
        sd = float(np.sqrt(np.sum(m2 * (hh - mean) ** 2)))
        out.update(mat_frac=1 - w,
                   mat_mode=float(hh[int(np.argmax(gaussian_filter1d(matc, 0.02 / A.dz)))]),
                   mat_med=float(hh[np.searchsorted(np.cumsum(m2), 0.5)]),
                   mat_mean=mean, mat_sd=sd,
                   mat_skew=float(np.sum(m2 * (hh - mean) ** 3) / max(sd ** 3, 1e-12)))
    return out


def main():
    M = marks(A.set_)
    print(f"set: {A.set_}   control: {CONTROL[A.set_]}   marks in table: {len(M):,}")
    print(f"  NO vertical window; bins {A.dz*1000:g} mm; mode scales {SCALES} m; "
          f"mixture bandwidth {A.bw*1000:g} mm")
    rows = []
    for t in M.itertuples():
        r = fit_one(t.point_id, t.elevation, A.set_)
        if r:
            rows.append(r)
        if rows and len(rows) % 50 == 0:
            pd.DataFrame(rows).to_csv(OUT, index=False)
            print(f"  ... {len(rows)} fitted, saved", flush=True)
    if not rows:
        raise SystemExit(f"no marks fitted: are the boxes present for {A.set_} in {BOX}?")
    F = pd.DataFrame(rows)
    F["lowveg"] = [lowveg(p, 0.15, 2.00, setname=A.set_) for p in F.point_id]
    F.to_csv(OUT, index=False)
    print(f"\nfitted {len(F):,} marks   h_max median {F.h_max.median():+.2f} m  "
          f"p90 {F.h_max.quantile(.9):+.2f}")

    bands = [("lowveg<=0.02", F.lowveg <= 0.02),
             ("0.02-0.10", (F.lowveg > 0.02) & (F.lowveg <= 0.10)),
             ("0.10-0.25", (F.lowveg > 0.10) & (F.lowveg <= 0.25)),
             ("lowveg>0.25", F.lowveg > 0.25)]
    print("\nMODE SHIFT (mode - surveyed, mm) -- median per cover band, per smoothing scale")
    print(f"  {'cover band':14s} {'n':>4s} " + " ".join(f"{f's={s:g}':>9s}" for s in SCALES))
    for lab, sel in bands:
        g = F[sel]
        if len(g) < 3:
            continue
        print(f"  {lab:14s} {len(g):4d} "
              + " ".join(f"{g[f'mode_shift_{s:g}'].median()*1000:9.1f}" for s in SCALES))
    print(f"  {'ALL':14s} {len(F):4d} "
          + " ".join(f"{F[f'mode_shift_{s:g}'].median()*1000:9.1f}" for s in SCALES))

    print("\nGROUND COMPONENT (mu pinned at surveyed) and THE MAT")
    print(f"  {'cover band':14s} {'n':>4s} {'w_g':>7s} {'sigma_g mm':>11s} {'mat mode mm':>12s} "
          f"{'mat sd mm':>10s} {'skew':>6s}")
    for lab, sel in bands:
        g = F[sel]
        if len(g) < 3:
            continue
        print(f"  {lab:14s} {len(g):4d} {g['w_g'].median():7.3f} "
              f"{g['sigma_g'].median()*1000:11.1f} {g['mat_mode'].median()*1000:12.1f} "
              f"{g['mat_sd'].median()*1000:10.1f} {g['mat_skew'].median():6.2f}")
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()

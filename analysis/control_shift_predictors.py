#!/usr/bin/env python3
"""Which LIDAR-ONLY property of the return column predicts the mode's offset from truth?

The response is `mode - surveyed ground` per control mark (analysis/control_mode_shift.py).
Every candidate below is measured from the lidar alone -- no surveyed elevation enters -- so
each could serve as a correction covariate on a tile where no truth exists. Anything fitted
with the truth in it (e.g. sigma_g from the mu-pinned mixture) is deliberately excluded.

NO WINDOWS: the whole column is used, out to whatever height the cloud reaches. Where a
candidate needs a scale (locating a mode, smoothing) the scale is SWEPT, not chosen.

The candidates, and why each might carry the signal:

  fwhm        width of the dominant peak. A mat both lifts and broadens the peak.
  nmad        robust spread of returns about the mode. Same idea, no peak-shape assumption.
  skew        asymmetry about the mode. A mat adds mass on ONE side only.
  kurt        peakedness. A clean ground return is sharp; a mat is not.
  mode_minus_median   pure shape asymmetry of the whole column, no scale at all.
  mode_minus_p10      how far the column extends BELOW the mode -- ground seen through gaps.
  p90_minus_mode      how far it extends above -- vegetation depth.
  n_peaks     prominent local maxima in the smoothed column. TWO means the ground and the
              mat are separable in this column; ONE means they are not.
  ret_per_pulse   returns divided by first returns. A standard vegetation indicator that
              needs no height coordinate at all -- a pulse that splits met something.
  frac_multi  fraction of returns belonging to multi-return pulses. Same signal, per return.

    ./lidar-icp/bin/python analysis/control_shift_predictors.py --set gen2_2021_control
"""
import argparse
import os
import sys

import laspy
import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter1d
from scipy.signal import find_peaks
from scipy.stats import spearmanr

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from control_mode_shift import CONTROL, STRUCT, BOX, marks
from lidar_diff_icp.groundtruth.tie import _design

ap = argparse.ArgumentParser()
ap.add_argument("--set", dest="set_", choices=sorted(CONTROL), default="gen2_2021_control")
ap.add_argument("--smooth", default="0.005,0.02,0.05")
ap.add_argument("--dz", type=float, default=0.005)
ap.add_argument("--out", default=None)
A = ap.parse_args() if __name__ == "__main__" else ap.parse_args([])
SCALES = [float(s) for s in A.smooth.split(",")]
OUT = A.out or f"data/derived/control_shift_predictors_{A.set_}.csv"


def one(pid, surveyed_z, set_):
    sp, bp = f"{STRUCT}/{set_}__{pid}.npz", f"{BOX}/{set_}__{pid}.laz"
    if not (os.path.exists(sp) and os.path.exists(bp)):
        return None
    z = np.load(sp)
    coef = z["surface_coef"]
    E, N, R = float(z["easting"]), float(z["northing"]), float(z["struct_radius"])
    f = laspy.read(bp)
    x, y, zz = np.asarray(f.x), np.asarray(f.y), np.asarray(f.z)
    rn, nr = np.asarray(f.return_number), np.asarray(f.number_of_returns)
    inten = np.asarray(f.intensity).astype(float)
    sel = np.hypot(x - E, y - N) <= R
    if sel.sum() < 50:
        return None
    nn = np.sqrt(1.0 + coef[1] ** 2 + coef[2] ** 2)
    h = (zz[sel] - (_design(x[sel] - E, y[sel] - N, 2) @ coef)) / nn
    mu_true = (float(surveyed_z) - float(coef[0])) / nn
    rn_s, nr_s, it_s = rn[sel], nr[sel], inten[sel]

    e = np.arange(np.floor(h.min() / A.dz) * A.dz, np.ceil(h.max() / A.dz) * A.dz + A.dz, A.dz)
    c = np.histogram(h, bins=e)[0].astype(float)
    d = 0.5 * (e[:-1] + e[1:])
    out = dict(point_id=pid, n=float(h.size))

    # Pulse-based candidates: no height coordinate at all, so no scale and no window.
    n_first = float(np.sum(rn_s == 1))
    out["ret_per_pulse"] = float(h.size / n_first) if n_first else np.nan
    out["frac_multi"] = float(np.mean(nr_s > 1))

    # INTENSITY. 3DEP intensity is NOT radiometrically calibrated -- it varies with block,
    # AGC, range and incidence -- so absolute values are not comparable BETWEEN marks. The
    # measures here are therefore WITHIN-column: a rank correlation of intensity against
    # height, and contrasts between parts of the same column, both of which are invariant to
    # any per-mark scaling. int_median is kept only as the control that should NOT work.
    out["int_median"] = float(np.median(it_s))
    out["int_iqr"] = float(np.percentile(it_s, 75) - np.percentile(it_s, 25))
    if it_s.size > 20 and np.std(it_s) > 0:
        out["int_h_rho"] = float(spearmanr(it_s, h).statistic)      # scale-free, window-free
    fm = rn_s == 1
    lm = rn_s == nr_s
    if fm.any() and lm.any():
        # ratio, not difference: a per-mark gain cancels
        out["int_last_over_first"] = float(np.median(it_s[lm]) / max(np.median(it_s[fm]), 1e-9))
    out["int_multi_over_single"] = (float(np.median(it_s[nr_s > 1]) /
                                          max(np.median(it_s[nr_s == 1]), 1e-9))
                                    if (nr_s > 1).any() and (nr_s == 1).any() else np.nan)

    for s in SCALES:
        k = max(s / A.dz, 1e-6)
        cs = gaussian_filter1d(c, k)
        i = int(np.argmax(cs))
        mode = d[i]
        out[f"shift_{s:g}"] = float(mode - mu_true)              # THE RESPONSE
        half = 0.5 * cs[i]
        li = np.where(cs[:i] < half)[0]
        ri = np.where(cs[i + 1:] < half)[0]
        out[f"fwhm_{s:g}"] = (float(d[i + 1 + ri[0]] - d[li[-1]]) if li.size and ri.size
                              else np.nan)
        out[f"nmad_{s:g}"] = float(1.4826 * np.median(np.abs(h - mode)))
        dev = h - mode
        sd = float(np.std(dev))
        out[f"skew_{s:g}"] = float(np.mean(dev ** 3) / sd ** 3) if sd > 0 else np.nan
        out[f"kurt_{s:g}"] = float(np.mean(dev ** 4) / sd ** 4 - 3.0) if sd > 0 else np.nan
        out[f"mode_minus_median_{s:g}"] = float(mode - np.median(h))
        out[f"mode_minus_p10_{s:g}"] = float(mode - np.percentile(h, 10))
        out[f"p90_minus_mode_{s:g}"] = float(np.percentile(h, 90) - mode)
        pk, props = find_peaks(cs, prominence=0.10 * cs[i])
        out[f"n_peaks_{s:g}"] = float(pk.size)
        # Intensity CONTRAST across the mode: returns below it against returns at/above it.
        # If the mode is a vegetation mat, the returns beneath it reached the ground through
        # gaps and are attenuated, so the ratio should fall. Ratio again, so gain cancels.
        blo, bhi = h < mode, h >= mode
        out[f"int_below_over_above_{s:g}"] = (
            float(np.median(it_s[blo]) / max(np.median(it_s[bhi]), 1e-9))
            if blo.sum() > 5 and bhi.sum() > 5 else np.nan)
    return out


def main():
    M = marks(A.set_)
    rows = []
    for t in M.itertuples():
        r = one(t.point_id, t.elevation, A.set_)
        if r:
            rows.append(r)
        if rows and len(rows) % 50 == 0:
            pd.DataFrame(rows).to_csv(OUT, index=False)
            print(f"  ... {len(rows)}", flush=True)
    F = pd.DataFrame(rows)
    F.to_csv(OUT, index=False)
    print(f"\n{A.set_}: {len(F):,} marks. Spearman rho of each LIDAR-ONLY candidate against")
    print("the mode offset from surveyed ground. |rho| ranked within each smoothing scale.")
    CANDS = ["fwhm", "nmad", "skew", "kurt", "mode_minus_median", "mode_minus_p10",
             "p90_minus_mode", "n_peaks", "int_below_over_above"]
    for s in SCALES:
        y = F[f"shift_{s:g}"] * 1000
        res = []
        for cnd in CANDS + ["ret_per_pulse", "frac_multi", "int_median", "int_iqr",
                            "int_h_rho", "int_last_over_first", "int_multi_over_single"]:
            col = f"{cnd}_{s:g}" if cnd in CANDS else cnd
            x = F[col]
            m = np.isfinite(x) & np.isfinite(y)
            if m.sum() < 20 or x[m].nunique() < 3:
                continue
            r = spearmanr(x[m], y[m])
            res.append((abs(r.statistic), cnd, r.statistic, r.pvalue, int(m.sum())))
        res.sort(reverse=True)
        print(f"\n  smoothing {s*1000:g} mm   (response median {y.median():+.1f} mm)")
        print(f"    {'candidate':>18s} {'rho':>8s} {'p':>11s} {'n':>5s}")
        for _, cnd, rho, p, n in res:
            print(f"    {cnd:>18s} {rho:+8.3f} {p:11.2e} {n:5d}")
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()

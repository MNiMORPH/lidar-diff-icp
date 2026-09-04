#!/usr/bin/env python3
# RECOVERED 2026-09-04 from the session transcript. This was written, run, and its
# output shipped -- but the code itself was never committed, so the CSV in
# data/derived had no producer and the measurement under it was not reproducible.
# Restored verbatim apart from this header and the output path.
"""Mat characterisation with NO vertical window: every return in the box.

The -1..+2 m window was inherited from the archived histograms and never justified. It
also made w_g window-dependent -- a forest mark's canopy returns were excluded from the
denominator, inflating the ground share. Andy, 2026-09-03: no windows.

The only remaining spatial restriction is the mark's 7.5 m struct_radius, which defines
what "at this mark" means and is not a vertical cut.
"""
import os, sys, numpy as np, pandas as pd, laspy
sys.path.insert(0, "analysis")
from scipy.ndimage import gaussian_filter1d
from scipy.stats import norm
from lidar_diff_icp.groundtruth.tie import _design
from control_lowveg_offset import STRUCT, load
BOX, SET = "data/derived/control_boxes", "gen2_2021_control"
DZ, BW, ITERS = 0.005, 0.15, 150          # BW is MINE; DZ << the ~50 mm ground width
OUT = "data/derived/control_mat_fit_nowindow.csv"
GC = np.arange(-0.20, 1.50, 0.01)

def fit(pid, surveyed_z):
    sp, bp = f"{STRUCT}/{SET}__{pid}.npz", f"{BOX}/{SET}__{pid}.laz"
    if not (os.path.exists(sp) and os.path.exists(bp)): return None
    z = np.load(sp); coef = z["surface_coef"]
    E, N, R = float(z["easting"]), float(z["northing"]), float(z["struct_radius"])
    f = laspy.read(bp); x, y, zz = np.asarray(f.x), np.asarray(f.y), np.asarray(f.z)
    sel = np.hypot(x - E, y - N) <= R
    if sel.sum() < 50: return None
    nn = np.sqrt(1 + coef[1]**2 + coef[2]**2)
    h = (zz[sel] - (_design(x[sel]-E, y[sel]-N, 2) @ coef)) / nn      # NO window
    mu = (float(surveyed_z) - float(coef[0])) / nn
    e = np.arange(np.floor(h.min()/DZ)*DZ, np.ceil(h.max()/DZ)*DZ + DZ, DZ)
    c = np.histogram(h, bins=e)[0].astype(float); d = 0.5*(e[:-1]+e[1:])
    sig, w, r = 0.02, 0.5, np.full(d.size, 0.5)
    for _ in range(ITERS):
        veg = gaussian_filter1d(np.maximum(c*(1-r), 0.0), BW/DZ, mode="nearest")
        veg = veg/max(veg.sum()*DZ, 1e-12)
        g = norm.pdf(d, mu, sig); den = w*g + (1-w)*veg; ok = den > 0
        r = np.zeros_like(d); r[ok] = w*g[ok]/den[ok]
        cr = c*r; sw = cr.sum(); w = float(sw/max(c.sum(), 1.0))
        if sw > 0: sig = float(max(np.sqrt(np.sum(cr*(d-mu)**2)/sw), 0.002))
    matc = c*(1-r); tot = matc.sum()
    if tot <= 0: return None
    hh = d-mu; m2 = matc/tot
    mean = float(np.sum(m2*hh)); sd = float(np.sqrt(np.sum(m2*(hh-mean)**2)))
    return dict(point_id=pid, w_g=w, sigma_g=sig, mat_frac=1-w, n_all=float(h.size),
                h_min=float(h.min()), h_max=float(h.max()),
                mat_mode=float(hh[int(np.argmax(gaussian_filter1d(matc, 0.02/DZ)))]),
                mat_med=float(hh[np.searchsorted(np.cumsum(m2), 0.5)]),
                mat_mean=mean, mat_sd=sd,
                mat_skew=float(np.sum(m2*(hh-mean)**3)/max(sd**3, 1e-12)),
                mat_below_500=float(m2[(hh > 0) & (hh <= 0.5)].sum()))

m = load(0.15, 2.00)
rows = []
for i, t in enumerate(m.itertuples()):
    r = fit(t.point_id, t.elevation)
    if r: rows.append(r)
    if rows and len(rows) % 40 == 0:
        pd.DataFrame(rows).to_csv(OUT, index=False)
        print(f"  ... {len(rows)} fitted, saved", flush=True)
pd.DataFrame(rows).to_csv(OUT, index=False)
print(f"  wrote {OUT}  ({len(rows)} marks)", flush=True)

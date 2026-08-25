#!/usr/bin/env python3
"""Degree-by-degree gen1 residual vs slope for two beam selections, with error bars:
  - low SCAN angle  |scan_angle| < 5 deg  (near-nadir; beam ~vertical)
  - low INCIDENCE   |incidence|  < 5 deg  (near-surface-perpendicular)
Residual r = d_mm - datum, datum = flat-ground (<3 deg) median d_mm, so r<0 = gen1
reads the ground low. Error bar = robust SE of the median = 1.2533 * NMAD / sqrt(n).

    ./lidar-icp/bin/python analysis/ridgelines/nearnadir_vs_perp_slope.py
"""
import argparse, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ap = argparse.ArgumentParser()
ap.add_argument("--tile", default="data/derived/elba_fulldensity")   # any tile with gen1_csf_angles
A_ = ap.parse_args()
TILE = os.path.basename(A_.tile.rstrip("/"))
TAG = "" if TILE == "elba_fulldensity" else f"_{TILE}"
# default tile keeps its original title/filename, so the elba figures stay byte-identical
TITLE_TAG = "" if not TAG else f"  ({TILE})"
A = np.load(f"{A_.tile}/gen1_csf_angles.npz")
d = A["d_mm"].astype(float); sa = np.abs(A["scan_angle"].astype(float))
inc = np.abs(A["incidence"].astype(float)); sl = A["slope"].astype(float)
# in_grid is essential off elba: returns outside the gen2 grid carry a fill slope (~24.5 deg)
# and |d| of ~1 km, which would dump ~1.65M garbage points into one slope bin on elbaext.
ing = A["in_grid"].astype(bool)
ok = ing & np.isfinite(d) & np.isfinite(sa) & np.isfinite(inc) & np.isfinite(sl)
d, sa, inc, sl = d[ok], sa[ok], inc[ok], sl[ok]
datum = np.median(d[sl < 3]); r = d - datum

def nmad(x): return 1.4826 * np.median(np.abs(x - np.median(x)))

def curve(mask, nmin=100):
    edges = np.arange(0, 41, 1.0)
    xc, med, se = [], [], []
    for lo, hi in zip(edges[:-1], edges[1:]):
        v = r[mask & (sl >= lo) & (sl < hi)]
        if v.size < nmin:
            continue
        xc.append((lo + hi) / 2); med.append(np.median(v))
        se.append(1.2533 * nmad(v) / np.sqrt(v.size))
    return np.array(xc), np.array(med), np.array(se)

xs, ms, es = curve(sa < 5)      # near-nadir
xp, mp, ep = curve(inc < 5)     # perpendicular

fig, ax = plt.subplots(figsize=(9, 5.2), dpi=150)
ax.axhline(0, color="0.6", lw=0.8, zorder=0)
ax.errorbar(xs, ms, yerr=es, fmt="o-", ms=4, lw=1.3, capsize=2, color="#1f77b4",
            label=f"low scan angle  |scan|<5°  (near-nadir), datum {datum:+.0f} mm")
ax.errorbar(xp, mp, yerr=ep, fmt="s-", ms=4, lw=1.3, capsize=2, color="#d62728",
            label="low incidence  |inc|<5°  (perpendicular to surface)")
ax.set_xlabel("slope (deg)")
ax.set_ylabel("gen1 depth below gen2 surface  (mm)   [r<0 = gen1 low]")
ax.set_title("gen1 near-nadir vs surface-perpendicular beams, per 1° slope bin" + TITLE_TAG)
ax.legend(loc="lower left", fontsize=9)
ax.set_xlim(0, 40)
ax.grid(True, alpha=0.25)
fig.tight_layout()
out = f"analysis/ridgelines/nearnadir_vs_perp_slope{TAG}.png"
fig.savefig(out)
print("wrote", out, "size", fig.get_size_inches() * fig.dpi)
print("near-nadir  last x:", xs[-1] if len(xs) else None,
      " perpendicular last x:", xp[-1] if len(xp) else None)

#!/usr/bin/env python3
"""Degree-by-degree gen1 residual vs slope for two beam selections, with error bars:
  - large SCAN angle  |scan_angle| > 12 deg  (most oblique; beam closest to horizontal)
  - low   SCAN angle  |scan_angle| <  5 deg  (near-nadir; reference curve)
Residual r = d_mm - datum, datum = flat-ground (<3 deg) median d_mm, so r<0 = gen1
reads the ground low. Error bar = robust SE of the median = 1.2533 * NMAD / sqrt(n).

    ./lidar-icp/bin/python analysis/ridgelines/lowangle_slope_dependence.py
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
sl = A["slope"].astype(float)
ok = np.isfinite(d) & np.isfinite(sa) & np.isfinite(sl)
d, sa, sl = d[ok], sa[ok], sl[ok]
datum = np.median(d[sl < 3]); r = d - datum

def nmad(x): return 1.4826 * np.median(np.abs(x - np.median(x)))

def curve(mask, nmin=100):
    edges = np.arange(0, 41, 1.0)
    xc, med, se, ns = [], [], [], []
    for lo, hi in zip(edges[:-1], edges[1:]):
        v = r[mask & (sl >= lo) & (sl < hi)]
        if v.size < nmin:
            continue
        xc.append((lo + hi) / 2); med.append(np.median(v))
        se.append(1.2533 * nmad(v) / np.sqrt(v.size)); ns.append(v.size)
    return np.array(xc), np.array(med), np.array(se), np.array(ns)

xh, mh, eh, nh = curve(sa > 12)     # closest to horizontal (subject)
xs, ms, es, ns = curve(sa < 5)      # near-nadir (reference)

# --- degree-by-degree table for the |scan|>12 subject subset ---
print(f"datum (flat-ground <3 deg median d_mm) = {datum:+.1f} mm")
print("|scan|>12 deg  subject subset")
print(f"{'slope_bin':>10}  {'median_r_mm':>12}  {'n':>8}")
edges = np.arange(0, 41, 1.0)
for lo, hi in zip(edges[:-1], edges[1:]):
    v = r[(sa > 12) & (sl >= lo) & (sl < hi)]
    if v.size < 100:
        print(f"{lo:4.0f}-{hi:<4.0f}  {'--':>12}  {v.size:>8}")
    else:
        print(f"{lo:4.0f}-{hi:<4.0f}  {np.median(v):>12.1f}  {v.size:>8}")

fig, ax = plt.subplots(figsize=(9, 5.2), dpi=150)
ax.axhline(0, color="0.6", lw=0.8, zorder=0)
ax.errorbar(xh, mh, yerr=eh, fmt="o-", ms=4, lw=1.3, capsize=2, color="#d62728",
            label=f"large scan angle  |scan|>12°  (closest to horizontal), datum {datum:+.0f} mm")
ax.errorbar(xs, ms, yerr=es, fmt="s-", ms=4, lw=1.3, capsize=2, color="#1f77b4",
            label="low scan angle  |scan|<5°  (near-nadir, reference)")
ax.set_xlabel("slope (deg)")
ax.set_ylabel("gen1 depth below gen2 surface  (mm)   [r<0 = gen1 low]")
ax.set_title("gen1 most-oblique vs near-nadir beams, per 1° slope bin" + TITLE_TAG)
ax.legend(loc="lower left", fontsize=9)
ax.set_xlim(0, 40)
ax.grid(True, alpha=0.25)
fig.tight_layout()
out = f"analysis/ridgelines/lowangle_slope_dependence{TAG}.png"
fig.savefig(out)
print("wrote", out, "size", fig.get_size_inches() * fig.dpi)

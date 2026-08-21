#!/usr/bin/env python3
"""Statistical analysis of the ridgeline data + CDF of curvature (transverse convexity kappa
at every ridgecrest pixel), split by land cover. Displayed as a figure.
    env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python analysis/ridgelines/curvature_cdf.py
"""
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

R = np.load("data/derived/elba_fulldensity/ridgecrest_pixels.npz", allow_pickle=True)
kappa = R["curvature_kappa"].astype(float)
slope = R["slope_deg"].astype(float)
dod = R["dod_m"].astype(float)*1000
lc = R["landcover"].astype(str)
forest = lc == "forest"; openc = lc == "open"

def describe(name, k):
    q = np.percentile(k, [5,25,50,75,95])
    print(f"  {name:14s} n={k.size:>5}  min={k.min():+.4f} "
          f"p5={q[0]:+.4f} p25={q[1]:+.4f} med={q[2]:+.4f} p75={q[3]:+.4f} "
          f"p95={q[4]:+.4f} max={k.max():+.4f}")

print("=== curvature kappa (1/m; +=convex-up) statistics on ridgecrest pixels ===")
describe("ALL crests", kappa)
describe("forest", kappa[forest])
describe("open", kappa[openc])
print("\n=== companion stats ===")
print(f"  slope (deg): med={np.median(slope):.1f} p25={np.percentile(slope,25):.1f} "
      f"p75={np.percentile(slope,75):.1f}")
print(f"  DoD (mm):    med={np.median(dod):+.1f} NMAD={1.4826*np.median(np.abs(dod-np.median(dod))):.0f}")

def cdf(a):
    s = np.sort(a); return s, np.arange(1, s.size+1)/s.size

fig, ax = plt.subplots(1, 2, figsize=(15, 6))
for lab, k, col in [("all crests", kappa, "k"), ("forest", kappa[forest], "green"),
                    ("open", kappa[openc], "goldenrod")]:
    x, y = cdf(k); ax[0].plot(x, y, color=col, lw=2, label=f"{lab} (n={k.size})")
ax[0].axvline(0.004, color="red", ls="--", lw=1, label="crest threshold κ=0.004")
ax[0].set_xlabel("transverse curvature κ (1/m, +=convex)"); ax[0].set_ylabel("CDF")
ax[0].set_title("CDF of ridgecrest curvature"); ax[0].set_xlim(0, np.percentile(kappa,99))
ax[0].legend(); ax[0].grid(alpha=0.3)
# zoomed/log view
for lab, k, col in [("forest", kappa[forest], "green"), ("open", kappa[openc], "goldenrod")]:
    x, y = cdf(k); ax[1].plot(x, y, color=col, lw=2, label=lab)
ax[1].set_xscale("log"); ax[1].set_xlabel("κ (1/m, log)"); ax[1].set_ylabel("CDF")
ax[1].set_title("CDF of curvature (log-x), forest vs open"); ax[1].legend(); ax[1].grid(alpha=0.3)
fig.savefig("figures/ridgecrest_curvature_cdf.png", dpi=130, bbox_inches="tight"); plt.close(fig)
print("\nwrote figures/ridgecrest_curvature_cdf.png")

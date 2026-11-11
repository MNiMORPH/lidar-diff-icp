#!/usr/bin/env python3
"""Does a horizontal tie of base B2 tighten the hillslope points against the lidar?
Andy, 2026-11-11: "I am guessing our tie point is 11. Try it and see how well it helps the transform."
Candidate ties: T1a 63 occupied the same spot as CORS 11 (Andy's guess) or CORS 10; the
shift = CORS minus T1a, applied to every B2 point (and the benchmark occupation), then
gen2 re-sampled bilinearly at the shifted positions. Reported per group: median and NMAD of
tied RTK minus gen2, FIX points, for no shift / shift-by-10 / shift-by-11, plus the shift that
MINIMISES the hillslope NMAD on a 0.1 m grid within +/-2.5 m (a free fit, for comparison only).
"""
import glob
import json
import os

import numpy as np
import pandas as pd

P = pd.read_csv("data/field/2026-10-25/rtk_points_all.csv")
tiles = {}
for c in glob.glob("data/derived/4*/corrections.json"):
    t = os.path.basename(os.path.dirname(c)); j = json.load(open(c)); tiles[t] = (j["bounds"], float(j.get("res") or j.get("res_m")), np.load(f"data/derived/{t}/z_after.npy", mmap_mode="r"))
def z2(x, y):
    for t, (b, res, a) in tiles.items():
        if b[0] <= x < b[2] and b[1] <= y < b[3]:
            fy, fx = (y - b[1]) / res - 0.5, (x - b[0]) / res - 0.5; i0, j0 = int(np.floor(fy)), int(np.floor(fx)); wy, wx = fy - i0, fx - j0
            i0 = min(max(i0, 0), a.shape[0] - 2); j0 = min(max(j0, 0), a.shape[1] - 2); q = np.asarray(a[i0:i0 + 2, j0:j0 + 2], float)
            if not np.isfinite(q).all(): return np.nan
            return float((1 - wy) * ((1 - wx) * q[0, 0] + wx * q[0, 1]) + wy * ((1 - wx) * q[1, 0] + wx * q[1, 1]))
    return np.nan
F = P[P.status.str.startswith("FIX") & P.navd88_tied_m.notna()].copy(); isB2 = F.base == "B2"
d63 = P[(P.group == "floodplain_RSplus_A") & (P.point == 63)].iloc[0]
cands = {"no shift (benchmark vertical only)": (0.0, 0.0)}
for n in (10, 11):
    c = P[(P.group == "floodplain_CORS_Wickert") & (P.point == n)].iloc[0]; cands[f"B2 shifted to CORS {n} (dE {c.easting - d63.easting:+.2f}, dN {c.northing - d63.northing:+.2f})"] = (c.easting - d63.easting, c.northing - d63.northing)
groups = {"hillslope B2 (RS3 + T4a)": isB2 & F.group.str.startswith("hillslope"), "floodplain B2 (RS+ A + B)": isB2 & F.group.str.startswith("floodplain"), "benchmark": isB2 & (F.group == "benchmark_DR9925")}
def stats(dx, dy):
    g2 = np.array([z2(x + (dx if b else 0), y + (dy if b else 0)) for x, y, b in zip(F.easting, F.northing, isB2)]); r = (F.navd88_tied_m.values - g2) * 1000
    return {k: (np.nanmedian(r[m]), 1.4826 * np.nanmedian(np.abs(r[m] - np.nanmedian(r[m]))), int(np.isfinite(r[m]).sum())) for k, m in groups.items()}
print("tied RTK minus gen2, mm [median, NMAD, n]; the shift is applied to B2 points only")
for name, (dx, dy) in cands.items():
    st = stats(dx, dy); print(f"\n{name}:"); [print(f"  {k:<30} {v[0]:+7.0f} {v[1]:6.0f} {v[2]:4d}") for k, v in st.items()]
# free fit on the hillslope NMAD, for comparison
best = None
for dx in np.arange(-2.5, 2.51, 0.25):
    for dy in np.arange(-2.5, 2.51, 0.25):
        st = stats(dx, dy)["hillslope B2 (RS3 + T4a)"]
        if best is None or st[1] < best[2]: best = (dx, dy, st[1], st[0])
print(f"\nfree fit (0.25 m grid, +/-2.5 m): hillslope NMAD is smallest at dE {best[0]:+.2f} dN {best[1]:+.2f}: NMAD {best[2]:.0f} mm, median {best[3]:+.0f} mm  [comparison only; not a tie]")

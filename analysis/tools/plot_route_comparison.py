#!/usr/bin/env python3
"""Two routes, side by side, plus the difference between them.

Andy, 2026-09-11: "Build elbaext under both routes ... And build the usual figures so I
can see the differences."

`plot_dod_comparison.py` compares two ARRAYS inside one tile directory. The route fork
puts them in two SEPARATE tile directories built from the same inputs, so this compares
across directories and adds the panel that actually answers the question: what the choice
of route did to the DoD, cell by cell.

Panels, all on one colour scale for the two DoDs so they are visually comparable:
  1. route A's DoD            2. route B's DoD            3. B - A

The third panel has its OWN scale, set from the data rather than inherited, because the
between-route difference is expected to be far smaller than the DoD itself and would be
invisible at the DoD's limits. Its limit is printed in the title so it is never mistaken
for the same scale as the first two.

    ./lidar-icp/bin/python analysis/tools/plot_route_comparison.py \
        --a data/derived/elbaext_independent --b data/derived/elbaext_delong
"""
import argparse
import json
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from lidar_diff_icp.figures import grid_of
from lidar_diff_icp.viz import hillshade

ap = argparse.ArgumentParser()
ap.add_argument("--a", required=True, help="tile dir, route A (left)")
ap.add_argument("--b", required=True, help="tile dir, route B (middle)")
ap.add_argument("--dod", default="dod.npy")
ap.add_argument("--v", type=float, default=0.3,
                help="colour limit for the two DoD panels, m. MINE, matches figures.py")
ap.add_argument("--figdir", default="figures/sites")
ap.add_argument("--out", default=None)
A = ap.parse_args()


def label(tile):
    """Name the route from the product itself, never from the directory name -- a
    directory can be renamed, and a figure that mislabels which pipeline made it is worse
    than no figure."""
    p = os.path.join(tile, "corrections.json")
    c = json.load(open(p)) if os.path.exists(p) else {}
    r = c.get("route", "?")
    bits = []
    if c.get("correction_surface"):
        bits.append("surface")
    if c.get("along_track_drift"):
        bits.append("drift")
    if c.get("geoid_applied"):
        bits.append("geoid")
    return f"{r} [{' + '.join(bits) or 'no corrections'}]", c


X0, Y0, res, nx, ny = grid_of(A.a)
gb = grid_of(A.b)
if (X0, Y0, res, nx, ny) != gb:
    raise SystemExit(f"REFUSED: the two tiles are on different grids, {(X0,Y0,res,nx,ny)} "
                     f"vs {gb}. Differencing them cell-by-cell would compare unlike places.")

la, ca = label(A.a)
lb, cb = label(A.b)
da = np.load(os.path.join(A.a, A.dod))
db = np.load(os.path.join(A.b, A.dod))
Z21 = np.load(os.path.join(A.a, "z_after.npy"))
hs = hillshade(Z21, res, X0, Y0, fill_gaps=False)
ext = (X0, X0 + nx * res, Y0, Y0 + ny * res)

d = db - da
m = np.isfinite(d)
# a symmetric limit from the data, not inherited: the between-route difference is expected
# to be much smaller than the DoD and would be invisible at +/-0.3 m
vd = float(np.nanpercentile(np.abs(d[m]), 99)) if m.any() else 0.0
vd = max(vd, 1e-4)

fig, ax = plt.subplots(1, 3, figsize=(22, 8))
for a_, arr, lab in ((ax[0], da, la), (ax[1], db, lb)):
    a_.imshow(hs, extent=ext, origin="lower", cmap="gray", alpha=0.6)
    im = a_.imshow(arr, extent=ext, origin="lower", cmap="RdBu", vmin=-A.v, vmax=A.v)
    a_.set_title(f"{lab}\nDoD gen2 - gen1 (m); red = erosion, blue = deposition")
    fig.colorbar(im, ax=a_, shrink=0.55, extend="both")
ax[2].imshow(hs, extent=ext, origin="lower", cmap="gray", alpha=0.6)
im2 = ax[2].imshow(d, extent=ext, origin="lower", cmap="PuOr", vmin=-vd, vmax=vd)
ax[2].set_title(f"B - A: what the ROUTE did (m)\nOWN scale, +/-{vd:.3f} m (99th pct |diff|)"
                f" -- NOT the same as the panels left")
fig.colorbar(im2, ax=ax[2], shrink=0.55, extend="both")
for a_ in ax:
    a_.set_xlabel("Easting (m)"); a_.set_ylabel("Northing (m)")

s = (f"median {1000*np.nanmedian(d[m]):+.1f} mm    "
     f"NMAD {1000*1.4826*np.nanmedian(np.abs(d[m]-np.nanmedian(d[m]))):.1f} mm    "
     f"max|diff| {1000*np.nanmax(np.abs(d[m])):.0f} mm    "
     f"cells {int(m.sum()):,}    "
     f"stable 1sigma A {1000*(ca.get('stable_1sigma_m') or float('nan')):.1f} / "
     f"B {1000*(cb.get('stable_1sigma_m') or float('nan')):.1f} mm")
fig.suptitle(f"{os.path.basename(A.a.rstrip('/'))} vs {os.path.basename(A.b.rstrip('/'))}"
             f"\n{s}", y=1.02)
os.makedirs(A.figdir, exist_ok=True)
out = A.out or f"{A.figdir}/route_comparison_{os.path.basename(A.b.rstrip('/'))}.png"
fig.savefig(out, dpi=130, bbox_inches="tight"); plt.close(fig)
print(f"  {s}")
print(f"  wrote {out}")

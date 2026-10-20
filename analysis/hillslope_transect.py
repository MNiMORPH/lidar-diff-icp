#!/usr/bin/env python3
"""A walkable hillslope transect by lot 204 (Andy, 2026-10-20): "a hillslope transect that the
class could walk to look for spatial variability in land-surface elevation to compare
against the lidar data."

Candidate lines start at lot 204 (E 579979 N 4885340) and run uphill to the WNW through the
gen2-high patches and the big hillslope control (FIELD_SITES_RTK.md, "Hillslope gen2-high
patches by lot 204"). Each candidate is sampled every 5 m from the 5 m grids of the star
product (z_after = gen2 surface; dod = gen2 - gen1 after the stable-ground gauge; slope;
canopy cover; crest mask; floodplain mask) and, when it exists, the gen2 lowveg grid
(data/derived/elbaext_independent/lowveg_grid.npy). The line to WALK is the one that stays on
public land, keeps its median slope moderate, and spans the widest DoD range: the choice is
printed with its numbers, and all candidates are listed so Andy can pick another.
Waypoints every 25 m with lat/lon go to data/derived/valley/hillslope_transect.csv; the
profile figure to figures/hillslope_transect.png.

    ./lidar-icp/bin/python analysis/hillslope_transect.py [--length 450] [--az 290,300,...]
"""
import argparse
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyproj
import shapely
from scipy.ndimage import map_coordinates
from shapely.geometry import shape

S = "data/derived/elbaext_delong_star_nosurf"; X0, Y0, RES = 575600.0, 4882200.0, 5.0
LOT = (579979.0, 4885340.0)
ap = argparse.ArgumentParser(); ap.add_argument("--length", type=float, default=450.0); ap.add_argument("--az", default="250,260,270,280,290,300,310,320"); ap.add_argument("--step", type=float, default=5.0)
ap.add_argument("--start", type=float, nargs=2, default=None, metavar=("E", "N"), help="start point (default lot 204)"); ap.add_argument("--tag", default="")
A = ap.parse_args(); AZ = [float(a) for a in A.az.split(",")]
if A.start: LOT = tuple(A.start)
print(f"start E {LOT[0]:.1f} N {LOT[1]:.1f}")
L = {n: np.load(f"{S}/{n}.npy") for n in ("z_after", "dod", "slope", "canopy_cover_pfs", "lod", "stable", "floodplain_mask", "crest_mask")}
base = L["dod"] * 1000 - np.nanmedian((L["dod"] * 1000)[L["stable"].astype(bool) & np.isfinite(L["dod"])])
lv_path = "data/derived/elbaext_independent/lowveg_grid.npy"; lowveg = np.load(lv_path) if os.path.exists(lv_path) else None
print("gen2 lowveg grid:", "present" if lowveg is not None else "NOT built yet (profile omits it)")
feats = lambda n: json.load(open(f"data/derived/valley/public_land/{n}"))["features"]
pub = shapely.union_all([shape(f["geometry"]) for f in feats("wma_units.geojson")] + [shape(f["geometry"]) for f in feats("state_parks.geojson")])
to_ll = pyproj.Transformer.from_crs("EPSG:26915", "EPSG:4326", always_xy=True)

def sample(arr, xs, ys, order=1):
    return map_coordinates(np.nan_to_num(arr, nan=np.nan), [(ys - Y0) / RES - 0.5, (xs - X0) / RES - 0.5], order=order, mode="nearest")

cands = []
for az in AZ:
    s = np.arange(0, A.length + A.step / 2, A.step); xs = LOT[0] + s * np.sin(np.radians(az)); ys = LOT[1] + s * np.cos(np.radians(az))
    on_pub = shapely.contains_xy(pub, xs, ys); z = sample(L["z_after"], xs, ys); d = sample(base, xs, ys); sl = sample(L["slope"], xs, ys, 0); cc = sample(L["canopy_cover_pfs"], xs, ys); fp = sample(L["floodplain_mask"].astype(float), xs, ys, 0) > 0.5
    rec = dict(az=az, public=float(on_pub.mean()), rise_m=float(np.nanmax(z) - z[0]), slope_med=float(np.nanmedian(sl)), slope_max=float(np.nanmax(sl)), dod_min=float(np.nanmin(d)), dod_max=float(np.nanmax(d)), dod_range=float(np.nanmax(d) - np.nanmin(d)), dod_sd=float(np.nanstd(d)), cc_med=float(np.nanmedian(cc)), floodplain_share=float(fp.mean()), n_high=int(np.sum(d > sample(L["lod"] * 1000, xs, ys))))
    cands.append((rec, s, xs, ys, z, d, sl, cc))
C = pd.DataFrame([c[0] for c in cands]); pd.set_option("display.width", 220)
print(f"candidate transects from the start, length {A.length:.0f} m, sampled every {A.step:.0f} m (az = compass bearing):"); print(C.round(2).to_string(index=False))
# choice: fully public, median slope <= 25 deg (walkable; MINE), then widest DoD range
ok = C[(C.public >= 0.999) & (C.slope_med <= 25)]
pick = (ok if len(ok) else C).sort_values("dod_range", ascending=False).iloc[0]; k = int(C.index[C.az == pick.az][0])
rec, s, xs, ys, z, d, sl, cc = cands[k]
print(f"\nCHOSEN: bearing {pick.az:.0f} deg; rise {pick.rise_m:.0f} m over {A.length:.0f} m; median slope {pick.slope_med:.0f} deg (max {pick.slope_max:.0f}); DoD {pick.dod_min:.0f}..{pick.dod_max:.0f} mm; {pick.n_high} of {len(s)} samples above +LoD")
lv = sample(lowveg, xs, ys) if lowveg is not None else np.full_like(s, np.nan)
W = pd.DataFrame(dict(dist_m=s, easting=xs.round(1), northing=ys.round(1), z_gen2_m=z.round(3), dod_mm=d.round(0), slope_deg=sl.round(1), canopy_cover=cc.round(2), lowveg=lv.round(3)))
W["lat"], W["lon"] = zip(*[to_ll.transform(x, y) for x, y in zip(xs, ys)]); W["lat"], W["lon"] = W["lon"].round(6), W["lat"].round(6)
W.to_csv(f"data/derived/valley/hillslope_transect{A.tag}.csv", index=False); wp = W[np.isclose(W.dist_m % 25, 0)]
print("\nWAYPOINTS every 25 m:"); print(wp.to_string(index=False))
fig, ax = plt.subplots(4 if lowveg is not None else 3, 1, figsize=(10, 9), sharex=True)
ax[0].plot(s, z, color="#2a78d6"); ax[0].set_ylabel("gen2 z (m)"); ax[0].set_title(f"Hillslope transect from E {LOT[0]:.0f} N {LOT[1]:.0f}, bearing {pick.az:.0f} deg: gen2 surface, DoD, slope and cover along the line")
ax[1].plot(s, d, color="#e34948"); ax[1].axhline(0, color="k", lw=0.5); ax[1].fill_between(s, -sample(L["lod"] * 1000, xs, ys), sample(L["lod"] * 1000, xs, ys), color="#ccc", alpha=0.5, label="+/- LoD"); ax[1].set_ylabel("DoD gen2-gen1 (mm)"); ax[1].legend(loc="upper left")
ax[2].plot(s, sl, color="#eda100", label="slope (deg)"); ax[2].set_ylabel("slope (deg)"); a2 = ax[2].twinx(); a2.plot(s, cc, color="#1baf7a"); a2.set_ylabel("canopy cover", color="#1baf7a")
if lowveg is not None: ax[3].plot(s, lv, color="#4a3aa7"); ax[3].set_ylabel("gen2 lowveg")
ax[-1].set_xlabel("distance along the line (m)"); fig.tight_layout(); fig.savefig(f"figures/hillslope_transect{A.tag}.png", dpi=110); print(f"figure: figures/hillslope_transect{A.tag}.png; waypoints: data/derived/valley/hillslope_transect{A.tag}.csv")

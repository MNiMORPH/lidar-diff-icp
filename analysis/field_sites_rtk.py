#!/usr/bin/env python3
"""Candidate RTK-GNSS survey sites on PUBLIC land (Whitewater WMA) in the Elba tile, chosen
from the product's own layers, for Andy's geomorphology class.

Andy, 2026-10-18: an RTK survey independent of the reed canary grass; floodplain sites; also
hillslope areas that could use a resurvey; public land (WMA, state park) only.

WHAT THE SURVEY MEASURES, per class of site, and why each class is here:
  F1 sedge-correction test  flat floodplain cells in the most-vegetated pooled-span quartile
                            (>= 400 mm) that read erosion and were corrected: RTK ground
                            vs gen2 tests gen2's floodplain surface; RTK ground vs the mat
                            top (probe through, and the top) in late autumn measures the
                            sedge-mat offset the logistic encodes.
  F2 closed canopy          flat cells with NO gen1 return at the gen2 surface (the
                            penetration-failure flag): the reed-canary-grass-type case.
                            Ground + thatch top gives the species-and-season offset.
  F3 bare control           flat, least-vegetated quartile, |DoD| < 30 mm: the offset must
                            be ~0; this is the zero point of the survey.
  H1 forest hillslope       10-35 deg, forest, DoD > +LoD: gen2 (May green-up) reads high
                            here by hypothesis; RTK ground decides which epoch is right.
                            Leaf-off survey, longer occupations under canopy.
  H2 steep gen1-low         > 25 deg, DoD < -LoD: the near-nadir gen1 reads-low effect.
  H3 open hillslope control 10-35 deg, open cover, |DoD| < LoD/2: the hillslope zero point.

Public land: DNR 'Publicly Accessible State Wildlife Management Areas' (Geospatial Commons
feature service us_mn_state_dnr/bdry_dnr_wildlife_mgmt_areas_pub), fetched 2026-10-18 into
data/derived/public_land/. State park statutory boundaries intersecting the tile: none.
Parking lots and roads: us_mn_state_dnr/struc_dnr_wld_mgmt_areas_pub layers 0 and 1.
Boundaries are 'for general reference only' per DNR; verify at the gate.

MINE: cluster minimum 40 cells (0.1 ha) so a site holds a few dozen RTK points; top 4 per
class by size. Both are stated, not load-bearing.

    ./lidar-icp/bin/python analysis/field_sites_rtk.py
"""
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyproj
import shapely
from matplotlib.image import imread
from scipy.ndimage import binary_dilation, label, uniform_filter
from shapely.geometry import shape

S = "data/derived/elbaext_delong_star_nosurf"; PUB = "data/derived/public_land"; MIN_CELLS = 40; TOP = 4
dod = np.load(f"{S}/dod.npy") * 1e3; new = np.load(f"{S}/dod_floodplain_vegcorr.npy") * 1e3; D = np.load(f"{S}/floodplain_vegcorr_delta.npy") * 1e3
fp = np.load(f"{S}/floodplain_mask.npy").astype(bool); sl = np.load(f"{S}/slope.npy"); lod = np.load(f"{S}/lod.npy") * 1e3
sp = np.load(f"{S}/gen1_spread_mm.npy"); pf = np.load(f"{S}/penetration_failure.npy").astype(bool)
forest = np.load(f"{S}/forest_pfs.npy").astype(bool); openc = np.load(f"{S}/open_pfs.npy").astype(bool); cc = np.load(f"{S}/canopy_cover_pfs.npy")
c = json.load(open(f"{S}/corrections.json")); X0, Y0, X1, Y1 = c["bounds"]; res = c["res_m"]; ny, nx = dod.shape
yy, xx = np.mgrid[0:ny, 0:nx]; E = X0 + (xx + 0.5) * res; N = Y0 + (yy + 0.5) * res

wma = [shape(f["geometry"]) for f in json.load(open(f"{PUB}/whitewater_wma_units.geojson"))["features"]]
units = [f["properties"]["name"] for f in json.load(open(f"{PUB}/whitewater_wma_units.geojson"))["features"]]
wma_u = shapely.union_all(wma)
public = shapely.contains_xy(wma_u, E.ravel(), N.ravel()).reshape(ny, nx)
fac0 = json.load(open(f"{PUB}/wma_fac_0.geojson"))["features"]; fac1 = json.load(open(f"{PUB}/wma_fac_1.geojson"))["features"]
parking = [shape(f["geometry"]) for f in fac0 if f["properties"].get("facility_subgroup") == "Parking Lot"]
roads = [shape(f["geometry"]) for f in fac1]
print(f"public (WMA) cells on the grid: {public.sum():,} of {ny*nx:,} ({100*public.mean():.1f}%); units: {units}")
print(f"WMA parking lots in the search box: {len(parking)}; WMA roads: {len(roads)}")

# pooled 3x3 span, as the correction used it
import pyarrow.parquet as pq
LO, HI, DZ = -5000.0, 5000.0, 20.0; edges = np.arange(LO, HI + DZ, DZ); ctr = 0.5 * (edges[:-1] + edges[1:]); NB = ctr.size
t = pq.read_table(f"{S}/beam_offset_table.parquet", columns=["cell", "d_mm_corr", "in_grid"]); g = t["in_grid"].to_numpy().astype(bool)
cell = t["cell"].to_numpy()[g]; d = t["d_mm_corr"].to_numpy()[g].astype(float); ok = np.isfinite(d) & (d > LO) & (d < HI)
H = np.zeros((ny * nx, NB), np.int32); np.add.at(H, (cell[ok], np.clip(((d[ok] - LO) / DZ).astype(int), 0, NB - 1)), 1); H = H.reshape(ny, nx, NB)
Hp = np.rint(uniform_filter(H.astype(np.float32), size=(3, 3, 1), mode="nearest") * 9); C = np.cumsum(Hp, 2); tot = C[..., -1]
q = lambda p: ctr[np.argmax(C >= p * tot[..., None], axis=2)]; span_pool = np.where(tot > 0, q(0.90) - q(0.10), np.nan); del H, Hp, C

flat = fp & (sl <= 2.0) & np.isfinite(dod)
fin = np.isfinite(dod)
classes = {
    "F1 sedge-correction test": flat & public & (span_pool >= 400) & (dod < -lod) & (D > 0) & ~pf,
    "F2 closed canopy (zero penetration)": flat & public & pf,
    "F3 bare floodplain control": flat & public & (sp < 140) & (np.abs(dod) < 30),
    "H1 forest hillslope, gen2-high": fin & public & (sl >= 10) & (sl <= 35) & forest & (dod > lod),
    "H2 steep slope, gen1-low": fin & public & (sl > 25) & (dod < -lod),
    "H3 open hillslope control": fin & public & (sl >= 10) & (sl <= 35) & openc & (np.abs(dod) < 0.5 * lod),
}
to_ll = pyproj.Transformer.from_crs("EPSG:26915", "EPSG:4326", always_xy=True)
rows = []
for cname, m in classes.items():
    lab, n = label(m, structure=np.ones((3, 3))); sizes = np.bincount(lab.ravel()); order = [k for k in np.argsort(sizes[1:])[::-1] + 1 if sizes[k] >= MIN_CELLS][:TOP]
    print(f"\n{cname}: {m.sum():,} cells; clusters >= {MIN_CELLS} cells: {int((sizes[1:] >= MIN_CELLS).sum())}")
    for r, k in enumerate(order):
        mm = lab == k; e0, n0 = E[mm].mean(), N[mm].mean(); lon, lat = to_ll.transform(e0, n0)
        pt = shapely.Point(e0, n0)
        dpark = min(pt.distance(p) for p in parking) if parking else np.nan; droad = min(pt.distance(p) for p in roads) if roads else np.nan
        edge = (binary_dilation(mm, structure=np.ones((3, 3))) & ~mm & flat & ~pf).sum() if cname.startswith("F2") else np.nan
        rows.append(dict(site=f"{cname.split()[0]}-{r+1}", cls=cname, cells=int(mm.sum()), ha=round(mm.sum() * res * res / 1e4, 2),
                         easting=round(e0, 1), northing=round(n0, 1), lat=round(lat, 6), lon=round(lon, 6),
                         raw_dod_mm=round(float(np.median(dod[mm])), 0), corrected_dod_mm=round(float(np.median(new[mm])), 0) if cname.startswith("F") else np.nan,
                         gen1_span_mm=round(float(np.nanmedian(sp[mm])), 0), slope_deg=round(float(np.median(sl[mm])), 1),
                         canopy_cover=round(float(np.nanmedian(cc[mm])), 2), penfail_pct=round(100 * float(pf[mm].mean()), 0),
                         edge_cells=edge, dist_parking_m=round(dpark, 0), dist_wma_road_m=round(droad, 0)))
F = pd.DataFrame(rows); F.to_csv("data/derived/field_sites_rtk.csv", index=False)
pd.set_option("display.width", 250); print("\n" + F.drop(columns=["cls"]).to_string(index=False))

marks = [("1079_2021_MN (2021 LCP calibration)", 579659.824, 4885857.071, 219.366), ("L1O101 (2008 open mark, L1O)", 578790.555, 4882696.274, 223.291)]
print("\nControl marks inside the tile for the RTK tie:")
for nm, e_, n_, z_ in marks:
    lon, lat = to_ll.transform(e_, n_); print(f"  {nm}: E {e_:.1f} N {n_:.1f}  lat {lat:.6f} lon {lon:.6f}  z {z_:.3f} m  public: {bool(wma_u.contains(shapely.Point(e_, n_)))}")

img = imread("data/derived/ne_patch_imagery/elbaext_fsa2023_2m.jpg")
fig, axes = plt.subplots(1, 2, figsize=(22, 10.5))
cols = {"F1": "#2a78d6", "F2": "#eb6834", "F3": "#1baf7a", "H1": "#e87ba4", "H2": "#e34948", "H3": "#eda100"}
for ax, base, ttl in ((axes[0], None, "2023 air photo"), (axes[1], "dod", "raw DoD (mm), +/-300")):
    if base is None: ax.imshow(img, extent=[X0, X1, Y0, Y1], origin="upper")
    else:
        im = ax.imshow(np.where(fin, dod, np.nan), extent=[X0, X1, Y0, Y1], origin="lower", cmap="RdBu", vmin=-300, vmax=300, interpolation="nearest"); plt.colorbar(im, ax=ax, fraction=0.03, pad=0.01)
    for geom in wma:
        for poly in (geom.geoms if hasattr(geom, "geoms") else [geom]):
            xs, ys = poly.exterior.xy; ax.plot(xs, ys, color="white" if base is None else "black", lw=1.2)
    for p in parking: ax.plot(p.x, p.y, marker="P", color="yellow", ms=11, mec="black")
    for rd in roads:
        xs, ys = rd.xy; ax.plot(xs, ys, color="yellow", lw=1.0, alpha=0.8)
    for _, r in F.iterrows():
        k = r.site[:2]; ax.plot(r.easting, r.northing, marker="o", ms=11, color=cols[k], mec="black")
        ax.annotate(r.site, (r.easting, r.northing), xytext=(6, 6), textcoords="offset points", fontsize=9, color="white" if base is None else "black", weight="bold")
    for nm, e_, n_, _ in marks: ax.plot(e_, n_, marker="^", ms=12, color="magenta", mec="black"); ax.annotate(nm.split()[0], (e_, n_), xytext=(6, -12), textcoords="offset points", fontsize=8, color="magenta", weight="bold")
    ax.set_xlim(X0, X1); ax.set_ylim(Y0, Y1); ax.set_title(ttl, fontsize=11); ax.set_xticks([]); ax.set_yticks([])
from matplotlib.lines import Line2D
h = [Line2D([], [], marker="o", color=v, mec="black", ls="", ms=9, label=k) for k, v in cols.items()] + [Line2D([], [], marker="P", color="yellow", mec="black", ls="", ms=9, label="WMA parking"), Line2D([], [], color="yellow", label="WMA road"), Line2D([], [], color="white", label="WMA boundary"), Line2D([], [], marker="^", color="magenta", mec="black", ls="", ms=9, label="control mark")]
axes[0].legend(handles=h, loc="lower left", fontsize=9, framealpha=0.85)
fig.suptitle("Candidate RTK survey sites on Whitewater WMA land, Elba tile: F = floodplain, H = hillslope (see data/derived/field_sites_rtk.csv)", fontsize=12)
fig.tight_layout(); fig.savefig("figures/field_sites_rtk.png", dpi=100); print("\nfigure: figures/field_sites_rtk.png; table: data/derived/field_sites_rtk.csv")

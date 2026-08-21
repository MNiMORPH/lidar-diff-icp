#!/usr/bin/env python3
"""Label cells by TRUE Macrostrat (scale=large) bedrock unit, then split the FORESTED
hillside cells into Oneota Dolomite (caprock) vs the sandstones (Jordan + Lone Rock) and
compare their dz/dt-vs-curvature behaviour. (K_ag on farmland is loess atop dolostone, not
bedrock — this is the bedrock split, on the forested walls where the units outcrop.)

    PROJ_DATA=/usr/share/proj ./lidar-icp/bin/python analysis/ridgelines/geology_forest_split.py
"""
import requests, numpy as np
from pyproj import Transformer
from matplotlib.path import Path
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.ndimage import distance_transform_edt as edt
from lidar_diff_icp.viz import hillshade

RES = 5.0; X0, Y0 = 577492.8, 4882737.6; DT = 12.44
ll2u = Transformer.from_crs(4326, 26915, always_xy=True)
u2ll = Transformer.from_crs(26915, 4326, always_xy=True)
D = "data/derived/elba_fulldensity/"
z = np.load(D + "z_after.npy"); ny, nx = z.shape
zf = z.copy(); m = ~np.isfinite(zf)
if m.any(): zf = zf[tuple(edt(m, return_distances=False, return_indices=True))]

# --- fetch true Macrostrat polygons over the tile -------------------------------------
polys = {}
for fy in np.linspace(0.05, 0.95, 6):
    for fx in np.linspace(0.1, 0.9, 5):
        lon, lat = u2ll.transform(X0 + fx*nx*RES, Y0 + fy*ny*RES)
        try:
            r = requests.get(f"https://macrostrat.org/api/v2/geologic_units/map?lat={lat:.5f}&lng={lon:.5f}&scale=large&format=geojson_bare", timeout=30).json()
        except Exception: continue
        for f in r.get("features", []):
            mid = f["properties"].get("map_id")
            if mid in polys: continue
            rings = []
            for poly in f["geometry"]["coordinates"]:
                ext = np.asarray(poly[0]); xs, ys = ll2u.transform(ext[:,0], ext[:,1])
                rings.append(np.column_stack([xs, ys]))
            polys[mid] = (str(f["properties"].get("strat_name")), str(f["properties"].get("lith")), rings)
print("Macrostrat units found:", [(v[0]) for v in polys.values()])

# --- rasterize to unit code: 1=dolostone(Oneota) 2=sandstone(Jordan/Lone Rock) ---------
cx = X0 + (np.arange(nx)+0.5)*RES; cy = Y0 + (np.arange(ny)+0.5)*RES
CX, CY = np.meshgrid(cx, cy); pts = np.column_stack([CX.ravel(), CY.ravel()])
unit = np.zeros(ny*nx, np.int8)
for mid, (strat, lith, rings) in polys.items():
    inside = np.zeros(len(pts), bool)
    for ring in rings: inside |= Path(ring).contains_points(pts)
    dolo = ("Oneota" in strat) or ("dolostone" in lith.lower())
    unit[inside] = 1 if dolo else 2
unit = unit.reshape(ny, nx)
np.save(D + "geo_unit.npy", unit)
print(f"cells: dolostone(Oneota)={int((unit==1).sum())}  sandstone(Jordan/LoneRock)={int((unit==2).sum())}  "
      f"unlabeled={int((unit==0).sum())}")

# --- split FORESTED cells by unit, compare dz/dt vs curvature --------------------------
lap = np.load(D + "curv_laplacian.npy"); pen = np.load(D + "penetration.npy")
cover = np.load(D + "canopy_cover.npy")
dod = np.load("data/derived/elba_refdatum/dod_osm.npy"); rate = dod/DT
forest = (pen < 0.25) & np.isfinite(dod) & np.isfinite(lap)
def fit(y, X):
    c,*_=np.linalg.lstsq(X,y,rcond=None); p=X@c; return c,1-np.sum((y-p)**2)/np.sum((y-np.mean(y))**2)
print("\nFORESTED cells, dz/dt = K*Lap + b*cover, by bedrock unit:")
for code, lbl in [(1,"Oneota dolostone (caprock)"), (2,"Jordan/Lone Rock sandstone")]:
    mm = forest & (unit==code)
    if mm.sum() < 100: print(f"  {lbl}: n={mm.sum()} (too few)"); continue
    (K,b), r2 = fit(rate[mm], np.c_[lap[mm], cover[mm]])
    med = np.median(dod[mm])*1000
    print(f"  {lbl:30s}: n={mm.sum():>6}  K={K:+.4f} m2/yr  beta={b*1000:+.2f} mm/yr/cover  "
          f"medDoD={med:+.1f} mm  medElev={np.median(zf[mm]):.0f}m  R2={r2:.3f}")

# --- figures: unit map + dz/dt vs curvature by unit (forested) -------------------------
ext = (X0, X0+nx*RES, Y0, Y0+ny*RES); hs = hillshade(zf, RES, X0, Y0, fill_gaps=True)
fig, ax = plt.subplots(figsize=(10,11)); ax.imshow(hs, extent=ext, origin="lower", cmap="gray")
ov = np.zeros((ny,nx,4)); ov[unit==1]=(0.7,0.2,0.7,0.45); ov[unit==2]=(0.95,0.75,0.25,0.45)
ax.imshow(ov, extent=ext, origin="lower")
ax.set_title("Macrostrat bedrock: Oneota Dolomite (purple) / Jordan+Lone Rock sandstone (tan)")
ax.set_xlabel("Easting (m)"); ax.set_ylabel("Northing (m)")
fig.savefig("figures/refdatum/geo_units.png", dpi=130, bbox_inches="tight"); plt.close(fig)

fig, ax = plt.subplots(figsize=(9,6))
for code,lbl,col in [(1,"Oneota dolostone","purple"),(2,"sandstone","goldenrod")]:
    mm = forest & (unit==code)
    if mm.sum()<100: continue
    L=lap[mm]; y=rate[mm]*1000; e=np.quantile(L,np.linspace(0,1,13)); mid=0.5*(e[:-1]+e[1:])
    md=[np.median(y[(L>=e[i])&(L<e[i+1])]) for i in range(len(e)-1)]
    ax.plot(mid,md,"o-",color=col,label=f"{lbl} (n={mm.sum()})")
ax.axhline(0,color="k",lw=.5); ax.axvline(0,color="k",lw=.5)
ax.set_xlabel("Laplacian curvature ∇²z (1/m)"); ax.set_ylabel("dz/dt (mm/yr)")
ax.set_title("Forested hillslopes: dz/dt vs curvature by bedrock unit"); ax.legend(); ax.grid(alpha=.3)
fig.savefig("figures/refdatum/geo_forest_dzdt.png", dpi=130, bbox_inches="tight"); plt.close(fig)
print("wrote figures/refdatum/geo_units.png, geo_forest_dzdt.png")

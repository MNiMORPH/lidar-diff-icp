#!/usr/bin/env python3
"""Realism of the Clubb floodplain, scored the way Clubb et al. 2017 scored theirs: against the
FEMA National Flood Hazard Layer's Special Flood Hazard Area (zones A/AE, the 1%-annual-chance
floodplain), with reliability r = TP/(TP+FP), sensitivity s = TP/(TP+FN), quality
Q = TP/(TP+FP+FN), on the 10 m grid within 800 m of the DNR Whitewater centrelines (the FEMA
mapping covers only the studied reaches, so the score is computed where SFHA polygons exist:
within 1500 m of any SFHA polygon). The Elba pipeline mask is scored the same way on its tile.

Andy, 2026-10-20: "consider the influence of those parameters and the realism of where the
floodplain is identified." Inputs: data/derived/valley/fema/nfhl_sfha_esrijson.json (97 SFHA
polygons fetched 2026-10-20 from hazards.fema.gov NFHL layer 28), the Clubb runs under
data/derived/valley/clubb/, the Elba mask. Writes data/derived/valley/clubb/realism.csv and
figures/clubb_realism.png (baseline vs FEMA vs pipeline on the lower Main Branch).
    ./lidar-icp/bin/python analysis/clubb_realism.py [--runs corr1500_10m sens_qq0.01 ...]
"""
import argparse
import glob
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import rasterio
import shapely
from shapely.geometry import shape

OUT = "data/derived/valley/clubb"
ap = argparse.ArgumentParser(); ap.add_argument("--runs", nargs="*", default=None); A = ap.parse_args()
runs = A.runs or sorted(os.path.basename(os.path.dirname(p)) for p in glob.glob(f"{OUT}/*/*_valley.bil"))
feats = json.load(open("data/derived/valley/fema/nfhl_sfha_esrijson.json"))
def esri_poly(g): return shapely.MultiPolygon([shapely.Polygon(r) for r in g["rings"]]).buffer(0)
sfha = shapely.union_all([esri_poly(f["geometry"]) for f in feats]); print(f"FEMA SFHA: {len(feats)} polygons, {sfha.area / 1e4:,.0f} ha in the box")
riv = shapely.union_all([shape(f["geometry"]) for f in json.load(open("data/derived/valley/public_land/rivers_whitewater.geojson"))["features"]])
with rasterio.open(f"{OUT}/corr1500_10m/corr1500_10m_valley.bil") as r: tr = r.transform; ny, nx = r.shape
jy, jx = np.mgrid[0:ny, 0:nx]; X, Y = tr * (jx + 0.5, jy + 0.5)
# memory-lean: test only cells inside the corridor's bounding box, in chunks, against PREPARED
# geometries (the first version tested all 15.5 M cells three times and was killed for memory)
def inside(geom, limit_to=None):
    shapely.prepare(geom); out = np.zeros((ny, nx), bool); bx0, by0, bx1, by1 = geom.bounds
    cand = (X >= bx0) & (X <= bx1) & (Y >= by0) & (Y <= by1)
    if limit_to is not None: cand &= limit_to
    idx = np.flatnonzero(cand.ravel()); xs, ys = X.ravel()[idx], Y.ravel()[idx]
    for i in range(0, idx.size, 500_000): out.ravel()[idx[i:i + 500_000]] = shapely.contains_xy(geom, xs[i:i + 500_000], ys[i:i + 500_000])
    return out
corridor = inside(riv.buffer(800)); mapped = inside(sfha.buffer(1500), corridor); fema = inside(sfha, corridor)
dom = corridor & mapped; print(f"scoring domain: within 800 m of the Whitewater AND within 1500 m of an SFHA polygon: {dom.sum() * 0.01:,.0f} ha; FEMA SFHA inside it {np.sum(fema & dom) * 0.01:,.0f} ha")
def score(pred, d):
    tp = np.sum(pred & fema & d); fp = np.sum(pred & ~fema & d); fn = np.sum(~pred & fema & d)
    return dict(r=tp / max(tp + fp, 1), s=tp / max(tp + fn, 1), Q=tp / max(tp + fp + fn, 1), pred_ha=np.sum(pred & d) * 0.01)
rows = []
for run in runs:
    with rasterio.open(f"{OUT}/{run}/{run}_valley.bil") as r: fp = r.read(1) == 1
    if fp.shape != (ny, nx): print(f"  {run}: different grid {fp.shape}, skipped"); continue
    rows.append(dict(run=run, **score(fp, dom)))
# the Elba pipeline mask, on its tile only
m5 = np.load("data/derived/elbaext_delong_star_nosurf/floodplain_mask.npy").astype(bool); EX0, EY0 = 575600.0, 4882200.0
tile = (X >= EX0) & (X < EX0 + m5.shape[1] * 5) & (Y >= EY0) & (Y < EY0 + m5.shape[0] * 5)
pm = np.zeros((ny, nx), bool); c, rr = np.clip(((X[tile] - EX0) / 5).astype(int), 0, m5.shape[1] - 1), np.clip(((Y[tile] - EY0) / 5).astype(int), 0, m5.shape[0] - 1); pm[tile] = m5[rr, c]
rows.append(dict(run="pipeline mask (Elba tile only)", **score(pm, dom & tile)))
with rasterio.open(f"{OUT}/corr1500_10m/corr1500_10m_valley.bil") as r: fpb = r.read(1) == 1
rows.append(dict(run="corr1500_10m (Elba tile only)", **score(fpb, dom & tile)))
D = pd.DataFrame(rows); D.to_csv(f"{OUT}/realism.csv", index=False); pd.set_option("display.width", 200)
print("\nAGREEMENT WITH FEMA SFHA (r reliability, s sensitivity, Q quality; Clubb 2017 Table 4 had Q 0.39-0.68 on its sites):"); print(D.round(3).to_string(index=False))
fig, axs = plt.subplots(1, 3, figsize=(18, 7)); bb = (574000, 4880000, 590000, 4897000)
for ax, (nm, m) in zip(axs, (("Clubb corr1500_10m", fpb), ("FEMA SFHA (A/AE)", fema), ("pipeline convexity mask (Elba tile)", pm))):
    ax.imshow(np.where(m, 1, np.nan), extent=[tr.c, tr.c + nx * tr.a, tr.f + ny * tr.e, tr.f], cmap="Blues", vmin=0, vmax=1.4)
    for g in (riv.geoms if hasattr(riv, "geoms") else [riv]): ax.plot(*g.xy, color="#e34948", lw=0.6)
    ax.set_xlim(bb[0], bb[2]); ax.set_ylim(bb[1], bb[3]); ax.set_xticks([]); ax.set_yticks([]); ax.set_title(nm)
fig.tight_layout(); fig.savefig("figures/clubb_realism.png", dpi=90); print("figures/clubb_realism.png")

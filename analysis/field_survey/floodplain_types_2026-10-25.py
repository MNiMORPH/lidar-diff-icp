"""Plant-type classes from multispectral imagery over the floodplain survey field, and lidar minus RTK by class
(Andy, 2026-11-18: "whether plant types might correlate with the offset"). Two unsupervised classifications,
Wisconsin-style (all bands, one scene, clustering; no NDVI threshold), both MINE in their k:
  A. Sentinel-2 L2A 2026-10-17 (cloud 0.4 %), ten reflective bands at 10 m (20 m bands resampled by odc), k-means
     with k = 6 over the pixels of the product floodplain mask (data/derived/elbaext_delong_star_nosurf/
     floodplain_mask.npy, 5 m, majority at 10 m) within 1.5 km of the field, bands standardized; plus the same on
     a FOUR-DATE stack (2026-04-10, 07-14, 08-08, 10-17) so phenology enters (early green-up, hay cut, late green).
  B. NAIP 2021-08-15 (0.6 m) and 2025-09-29, four bands, k-means k = 5 within the field box, for within-field
     patches (shrub, swale, sward, cut, bare).
For each classification: cluster spectra (median per band), the 147 FIX points per cluster, and the residuals
r21 / r08 (lidar minus RTK) by cluster, raw and with rover medians removed; Kruskal-Wallis across clusters.
OUTPUT: data/field/2026-10-25/imagery/floodplain_types_* (CSV, PNG maps), log on stdout. Written incrementally."""
import pandas as pd, numpy as np, json, os, planetary_computer as pc, odc.stac, rasterio
from rasterio.merge import merge
from pystac_client import Client
from pyproj import Transformer
from scipy.stats import kruskal
from sklearn.cluster import KMeans
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from cmcrameri import cm as cmc
OUT = "data/field/2026-10-25/imagery"; S = "data/derived/elbaext_delong_star_nosurf"
F = pd.read_csv(f"{OUT}/floodplain_points_covariates.csv")
E0, N0 = F.easting_tied.mean(), F.northing_tied.mean()
j = json.load(open(f"{S}/corrections.json")); B5 = j["bounds"]; res5 = float(j.get("res") or j.get("res_m")); fp = np.load(f"{S}/floodplain_mask.npy").astype(bool)
cat = Client.open("https://planetarycomputer.microsoft.com/api/stac/v1", modifier=pc.sign_inplace)
to_ll = Transformer.from_crs("EPSG:6344", "EPSG:4326", always_xy=True)
def bbox_of(half):
    lo0, la0 = to_ll.transform(E0 - half, N0 - half); lo1, la1 = to_ll.transform(E0 + half, N0 + half); return [lo0, la0, lo1, la1]
BANDS = ["B02", "B03", "B04", "B05", "B06", "B07", "B08", "B8A", "B11", "B12"]
def s2_scene(date, half=1500):
    its = list(cat.search(collections=["sentinel-2-l2a"], bbox=bbox_of(half), datetime=f"{date}/{date}T23:59:59Z").items())
    it = sorted(its, key=lambda i: i.properties.get("eo:cloud_cover", 99))[0]
    ds = odc.stac.load([it], bands=BANDS + ["SCL"], bbox=bbox_of(half), resolution=10, chunks={}, crs="EPSG:26915")
    arr = np.stack([ds[b].astype("float32").squeeze().values for b in BANDS]); scl = ds.SCL.squeeze().values
    clear = ~np.isin(scl, [0, 1, 3, 8, 9, 10, 11]); arr[:, ~clear] = np.nan
    return arr, ds.x.values, ds.y.values, it.id
def fp_at(xs, ys):
    X, Y = np.meshgrid(xs, ys); i = ((Y - B5[1]) / res5).astype(int); jj = ((X - B5[0]) / res5).astype(int)
    ok = (i >= 0) & (i < fp.shape[0]) & (jj >= 0) & (jj < fp.shape[1]); m = np.zeros(X.shape, bool); m[ok] = fp[i[ok], jj[ok]]; return m
def idx(xs, ys, x, y):
    jx = np.clip(np.searchsorted(xs, x) - 1, 0, len(xs) - 1); iy = np.clip(len(ys) - np.searchsorted(ys[::-1], y), 0, len(ys) - 1); return iy, jx
def report(tag, labels_at_points, k):
    F[f"cls_{tag}"] = labels_at_points
    print(f"\n{tag}: points per cluster and lidar minus RTK by cluster (mm; raw median, rover-median-removed median, n):")
    rows = []
    for c in range(k):
        s = F[F[f"cls_{tag}"] == c]
        if len(s) == 0: continue
        rows.append(dict(cluster=c, n=len(s), r21=s.r21_mm.median(), r21_rel=s.r21_mm_rel.median(), r08=s.r08_mm.median(), r08_rel=s.r08_mm_rel.median(), n08=int(s.r08_mm.notna().sum()), rovers=dict(s.rover.value_counts())))
        print(f"  cluster {c}: n {len(s):3d}  r21 {s.r21_mm.median():+5.0f} / rel {s.r21_mm_rel.median():+4.0f}   r08 {s.r08_mm.median():+5.0f} / rel {s.r08_mm_rel.median():+4.0f} (n08 {int(s.r08_mm.notna().sum())})  rovers {dict(s.rover.value_counts())}")
    g = [s.r21_mm_rel.values for _, s in F.groupby(f"cls_{tag}") if len(s) >= 5]
    g8 = [s.r08_mm_rel.dropna().values for _, s in F.groupby(f"cls_{tag}") if s.r08_mm_rel.notna().sum() >= 5]
    if len(g) >= 2: print(f"  Kruskal-Wallis across clusters with n >= 5: r21_rel p {kruskal(*g).pvalue:.3f}; r08_rel p {kruskal(*g8).pvalue:.3f}")
    pd.DataFrame(rows).to_csv(f"{OUT}/floodplain_types_{tag}_clusters.csv", index=False); F.to_csv(f"{OUT}/floodplain_types_points.csv", index=False)
def map_fig(tag, lab, xs, ys, k, extent_half=260):
    fig, ax = plt.subplots(figsize=(7, 8)); ex = [xs[0], xs[-1], ys[-1], ys[0]]
    ax.imshow(np.ma.masked_invalid(lab), extent=ex, cmap=plt.get_cmap("tab10", k), vmin=-0.5, vmax=k - 0.5, interpolation="nearest")
    sc = ax.scatter(F.easting_tied, F.northing_tied, c=F.r08_mm, cmap=cmc.bam_r, vmin=-600, vmax=600, s=12, edgecolors="k", linewidths=0.3)
    ax.set_xlim(E0 - extent_half, E0 + extent_half); ax.set_ylim(N0 - extent_half, N0 + extent_half); ax.set_xticks([]); ax.set_yticks([])
    plt.colorbar(sc, ax=ax, shrink=0.6, label="2008 lidar minus RTK (mm)"); fig.savefig(f"{OUT}/floodplain_types_{tag}.png", dpi=140, bbox_inches="tight"); plt.close(fig); print("  wrote", f"{OUT}/floodplain_types_{tag}.png")
# ---- A1: one October scene, ten bands
K = 6
arr, xs, ys, sid = s2_scene("2026-10-17"); m = fp_at(xs, ys) & np.isfinite(arr).all(0)
print("A1 Sentinel-2", sid, "| floodplain-mask pixels in the 3 km box:", int(m.sum()))
X = arr[:, m].T; mu, sd = X.mean(0), X.std(0); km = KMeans(K, n_init=10, random_state=0).fit((X - mu) / sd)
lab = np.full(m.shape, np.nan); lab[m] = km.labels_
print("  cluster median reflectance (x1e-4) per band", BANDS, "and NDVI, NDRE, NDMI; pixel count:")
for c in range(K):
    v = X[km.labels_ == c]; med = np.median(v, 0); ndvi = (med[6] - med[2]) / (med[6] + med[2]); ndre = (med[7] - med[3]) / (med[7] + med[3]); ndmi = (med[6] - med[8]) / (med[6] + med[8])
    print(f"   {c}: {np.round(med).astype(int).tolist()}  NDVI {ndvi:.2f} NDRE {ndre:.2f} NDMI {ndmi:.2f}  n {int((km.labels_ == c).sum())}")
iy, jx = idx(xs, ys, F.easting_tied.values, F.northing_tied.values); report("s2_oct2026_k6", lab[iy, jx], K); map_fig("s2_oct2026_k6", lab, xs, ys, K)
# ---- A2: four dates stacked
stack = []; ids = []
for d in ("2026-04-10", "2026-07-14", "2026-08-08", "2026-10-17"):
    a, xs2, ys2, sid2 = s2_scene(d); assert a.shape == arr.shape; stack.append(a); ids.append(sid2)
A = np.concatenate(stack, 0); m2 = fp_at(xs, ys) & np.isfinite(A).all(0); print("\nA2 four-date stack", ids, "| pixels:", int(m2.sum()))
X = A[:, m2].T; mu, sd = X.mean(0), X.std(0); km = KMeans(K, n_init=10, random_state=0).fit((X - mu) / sd)
lab2 = np.full(m2.shape, np.nan); lab2[m2] = km.labels_
print("  cluster NDVI per date (Apr 10, Jul 14, Aug 8, Oct 17) from median bands; pixel count:")
for c in range(K):
    v = X[km.labels_ == c]; med = np.median(v, 0); nd = [(med[10 * t + 6] - med[10 * t + 2]) / (med[10 * t + 6] + med[10 * t + 2]) for t in range(4)]
    print(f"   {c}: NDVI {[round(x, 2) for x in nd]}  n {int((km.labels_ == c).sum())}")
report("s2_4dates_k6", lab2[iy, jx], K); map_fig("s2_4dates_k6", lab2, xs, ys, K)
# ---- B: NAIP within the field box
for want in ("2021-08-15", "2025-09-29"):
    its = [i for i in cat.search(collections=["naip"], bbox=bbox_of(300)).items() if i.properties["datetime"][:7] == want[:7]]
    srcs = [rasterio.open(i.assets["image"].href) for i in its]; crs = srcs[0].crs; T = Transformer.from_crs("EPSG:6344", crs, always_xy=True)
    x0, y0 = T.transform(E0 - 260, N0 - 260); x1, y1 = T.transform(E0 + 260, N0 + 260); img, tr = merge(srcs, bounds=(x0, y0, x1, y1)); [s.close() for s in srcs]
    Xn = img.reshape(4, -1).T.astype(float); kn = 5; km = KMeans(kn, n_init=5, random_state=0).fit((Xn - Xn.mean(0)) / Xn.std(0)); labn = km.labels_.reshape(img.shape[1:])
    print(f"\nB NAIP {want} ({len(its)} tiles, pixel {tr.a:.1f} m): cluster median DN [R,G,B,NIR], NDVI(DN), share of chip:")
    for c in range(kn):
        v = Xn[km.labels_ == c]; med = np.median(v, 0); print(f"   {c}: {np.round(med).astype(int).tolist()}  NDVI {(med[3] - med[0]) / (med[3] + med[0]):+.2f}  {100 * (km.labels_ == c).mean():.0f} %")
    px, py = T.transform(F.easting_tied.values, F.northing_tied.values); col = ((px - tr.c) / tr.a).astype(int); row = ((py - tr.f) / tr.e).astype(int)
    # majority label within 3 m of each point
    yy, xx = np.indices(labn.shape); gx = tr.c + (xx + 0.5) * tr.a; gy = tr.f + (yy + 0.5) * tr.e
    lab_pt = []
    for a, b in zip(px, py):
        w = labn[(gx - a) ** 2 + (gy - b) ** 2 <= 9.0]; lab_pt.append(np.bincount(w).argmax() if w.size else np.nan)
    report(f"naip_{want[:4]}_k5", np.array(lab_pt, float), kn)
    fig, ax = plt.subplots(figsize=(7, 8)); ex = [tr.c, tr.c + tr.a * img.shape[2], tr.f + tr.e * img.shape[1], tr.f]
    ax.imshow(labn, extent=ex, cmap=plt.get_cmap("tab10", kn), vmin=-0.5, vmax=kn - 0.5, interpolation="nearest")
    sc = ax.scatter(px, py, c=F.r08_mm, cmap=cmc.bam_r, vmin=-600, vmax=600, s=12, edgecolors="k", linewidths=0.3); ax.set_xticks([]); ax.set_yticks([])
    plt.colorbar(sc, ax=ax, shrink=0.6, label="2008 lidar minus RTK (mm)"); fig.savefig(f"{OUT}/floodplain_types_naip_{want[:4]}_k5.png", dpi=140, bbox_inches="tight"); plt.close(fig); print("  wrote", f"{OUT}/floodplain_types_naip_{want[:4]}_k5.png")

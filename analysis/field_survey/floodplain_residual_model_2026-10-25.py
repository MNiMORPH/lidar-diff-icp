"""A data product that best matches the observed spatial spread of lidar minus RTK on the floodplain field
(Andy, 2026-11-18: "generate a data product that best matches the observed spatial spread ... any combination of
spectral bands and/or other data sources"). TARGETS: r08 = 2008 lidar minus tied RTK (126 points) and r21 = 2021
lidar minus tied RTK (147), raw (the rover ties agree to ~20-35 mm on the 2021 surface, so the between-rover part
of the 2008 spread is treated as ground signal). FEATURES, all available as rasters so the fitted model can be
mapped: Sentinel-2 L2A ten reflective bands at 10 m on 2021-05-01 (the lidar day), 2021-08-14, 2026-04-10 and
2026-10-17, plus NDVI, NDRE, NDMI per scene; NAIP NDVI (digital numbers) for 2019-08-27, 2021-08-15, 2025-09-29
averaged to the 10 m grid; the pipeline's 5 m layers averaged to 10 m: z_after (2021 ground elevation),
gen1_spread_mm, lowveg_grid; the reed-canary-grass work's rasters nelike_p (2008 Landsat, 2021 S2, 2025 S2),
metric_retained, metric_senesc_doy, metric_ndvi_oct, metric_ndvi_apr. MODELS: (1) each single feature, Theil-Sen;
(2) ridge on all standardized features (alpha by inner CV); (3) lasso; (4) random forest (300 trees, depth 4);
(5) gradient boosting (shallow). SELECTION: spatially BLOCKED cross-validation, 40 m x 40 m blocks (MINE), repeated
over 20 random block-to-fold assignments with 5 folds; score = CV R^2 and the NMAD of the CV residual against the
NMAD of the target. The best model by CV R^2 is refit on all points and PREDICTED over the floodplain mask within
a 600 m box around the field at 10 m, written as GeoTIFF (EPSG:26915) and PNG with the observed points drawn.
Nothing is chosen by looking at the test folds except the final model; all scores are printed. OUTPUT:
data/field/2026-10-25/imagery/floodplain_model_* ; log on stdout. Written incrementally."""
import pandas as pd, numpy as np, json, os, warnings, planetary_computer as pc, odc.stac, rasterio
from rasterio.merge import merge
from rasterio.transform import from_origin
from pystac_client import Client
from pyproj import Transformer
from scipy.stats import theilslopes, spearmanr
from sklearn.linear_model import RidgeCV, LassoCV
from sklearn.ensemble import RandomForestRegressor, GradientBoostingRegressor
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from cmcrameri import cm as cmc
warnings.filterwarnings("ignore")
OUT = "data/field/2026-10-25/imagery"; S = "data/derived/elbaext_delong_star_nosurf"; R = f"{S}/rcg"
F = pd.read_csv(f"{OUT}/floodplain_points_covariates.csv")
E0, N0 = F.easting_tied.mean(), F.northing_tied.mean(); HALF = 600; RES = 10
xs = np.arange(E0 - HALF, E0 + HALF, RES) + RES / 2; ys = np.arange(N0 + HALF, N0 - HALF, -RES) - RES / 2; X, Y = np.meshgrid(xs, ys)
to_ll = Transformer.from_crs("EPSG:26915", "EPSG:4326", always_xy=True); lo0, la0 = to_ll.transform(E0 - HALF - 50, N0 - HALF - 50); lo1, la1 = to_ll.transform(E0 + HALF + 50, N0 + HALF + 50); BBOX = [lo0, la0, lo1, la1]
cat = Client.open("https://planetarycomputer.microsoft.com/api/stac/v1", modifier=pc.sign_inplace)
feat = {}
def nmad(v): v = np.asarray(v, float); v = v[np.isfinite(v)]; return 1.4826 * np.median(np.abs(v - np.median(v)))
# ---- Sentinel-2 scenes on the common grid
BANDS = ["B02", "B03", "B04", "B05", "B06", "B07", "B08", "B8A", "B11", "B12"]
for d in ("2021-05-01", "2021-08-14", "2026-04-10", "2026-10-17"):
    its = sorted(cat.search(collections=["sentinel-2-l2a"], bbox=BBOX, datetime=f"{d}/{d}T23:59:59Z").items(), key=lambda i: i.properties.get("eo:cloud_cover", 99))
    ds = odc.stac.load([its[0]], bands=BANDS + ["SCL"], bbox=BBOX, resolution=RES, chunks={}, crs="EPSG:26915")
    sx = ds.x.values; sy = ds.y.values; jx = np.clip(np.searchsorted(sx, xs) - 1, 0, len(sx) - 1); iy = np.clip(len(sy) - np.searchsorted(sy[::-1], ys), 0, len(sy) - 1)
    scl = ds.SCL.squeeze().values[np.ix_(iy, jx)]; clear = ~np.isin(scl, [0, 1, 3, 8, 9, 10, 11])
    b = {k: ds[k].astype("float32").squeeze().values[np.ix_(iy, jx)] for k in BANDS}
    for k in BANDS: a = b[k].copy(); a[~clear] = np.nan; feat[f"s2_{d}_{k}"] = a
    feat[f"s2_{d}_NDVI"] = (b["B08"] - b["B04"]) / (b["B08"] + b["B04"]); feat[f"s2_{d}_NDRE"] = (b["B8A"] - b["B05"]) / (b["B8A"] + b["B05"]); feat[f"s2_{d}_NDMI"] = (b["B08"] - b["B11"]) / (b["B08"] + b["B11"])
    for k in ("NDVI", "NDRE", "NDMI"): feat[f"s2_{d}_{k}"][~clear] = np.nan
    print("S2", its[0].id, "clear fraction on grid %.2f" % clear.mean())
# ---- NAIP NDVI averaged to the grid
def grid_mean(vals, gx, gy):
    j = ((gx - (E0 - HALF)) / RES).astype(int); i = ((N0 + HALF - gy) / RES).astype(int); ok = (i >= 0) & (i < len(ys)) & (j >= 0) & (j < len(xs)) & np.isfinite(vals)
    acc = np.zeros(X.shape); cnt = np.zeros(X.shape); np.add.at(acc, (i[ok], j[ok]), vals[ok]); np.add.at(cnt, (i[ok], j[ok]), 1); out = acc / np.where(cnt == 0, np.nan, cnt); return out
for want in ("2019-08", "2021-08", "2025-09"):
    its = [i for i in cat.search(collections=["naip"], bbox=BBOX).items() if i.properties["datetime"][:7] == want]
    srcs = [rasterio.open(i.assets["image"].href) for i in its]; img, tr = merge(srcs, bounds=(E0 - HALF, N0 - HALF, E0 + HALF, N0 + HALF)); [s.close() for s in srcs]
    r, nir = img[0].astype(float), img[3].astype(float); nd = (nir - r) / np.where(nir + r == 0, np.nan, nir + r)
    yy, xx = np.indices(nd.shape); gx = tr.c + (xx + 0.5) * tr.a; gy = tr.f + (yy + 0.5) * tr.e
    feat[f"naip_{want}_NDVI"] = grid_mean(nd.ravel(), gx.ravel(), gy.ravel()); print("NAIP", want, len(its), "tiles")
# ---- pipeline 5 m layers averaged to 10 m, and rcg rasters sampled
import glob
TILES = {}
for c in glob.glob("data/derived/4*/corrections.json"):
    jj = json.load(open(c)); b = jj["bounds"]
    if b[0] < E0 + HALF and b[2] > E0 - HALF and b[1] < N0 + HALF and b[3] > N0 - HALF: TILES[os.path.dirname(c)] = (b, float(jj.get("res") or jj.get("res_m")))
print("5 m tiles covering the box:", {os.path.basename(k): v[0] for k, v in TILES.items()})
def layer5(name, dirs=None):
    """Mean of the 5 m layer over each 10 m cell, mosaicked over the per-tile products that cover the box (the field
    straddles 4342-29-64 and 4358-29-01); a cell covered by two tiles averages both."""
    vals = []; gxs = []; gys = []
    for d, (B5, res5) in (dirs or TILES).items():
        if not os.path.exists(f"{d}/{name}"): continue
        # the tile arrays are stored BOTTOM-UP (row 0 = south, as compile_rtk's sample() and viz.py read them); the first
        # run indexed them top-down and so put the elevation of the mirror location into the feature stack (62 of the 126
        # field points read 234-250 m against an RTK ground of 217-218 m) -- fixed 2026-11-18
        a = np.load(f"{d}/{name}", mmap_mode="r"); i0 = max(int(((N0 - HALF) - B5[1]) / res5), 0); i1 = min(int(((N0 + HALF) - B5[1]) / res5) + 1, a.shape[0]); j0 = max(int(((E0 - HALF) - B5[0]) / res5), 0); j1 = min(int(((E0 + HALF) - B5[0]) / res5) + 1, a.shape[1])
        if i1 <= i0 or j1 <= j0: continue
        sub = np.asarray(a[i0:i1, j0:j1], float); ii, jj2 = np.indices(sub.shape); vals.append(sub.ravel()); gys.append((B5[1] + (i0 + ii + 0.5) * res5).ravel()); gxs.append((B5[0] + (j0 + jj2 + 0.5) * res5).ravel())
    return grid_mean(np.concatenate(vals), np.concatenate(gxs), np.concatenate(gys))
js = json.load(open(f"{S}/corrections.json")); SITE = {S: (js["bounds"], float(js.get("res") or js.get("res_m")))}
feat["z_after_m"] = layer5("z_after.npy"); feat["gen1_spread_mm"] = layer5("gen1_spread_mm.npy"); feat["lowveg_gen2"] = layer5("lowveg_grid.npy"); fpm = (layer5("floodplain_mask.npy", SITE) > 0.5) | (layer5("floodplain_mask_clubb.npy") > 0.5)   # union of the site geometric mask and the tiles Clubb mask (the field is Clubb floodplain, not geometric)
for name in ("nelike_p_ls_20081029", "nelike_p_s2_20211018", "nelike_p_s2_20251002", "metric_retained", "metric_senesc_doy", "metric_ndvi_oct", "metric_ndvi_apr"):
    with rasterio.open(f"{R}/{name}.tif") as src: v = np.array([x[0] for x in src.sample(zip(X.ravel(), Y.ravel()))], float).reshape(X.shape); v[v == src.nodata] = np.nan if src.nodata is not None else v
    feat[name] = v
names = list(feat); stack = np.stack([feat[k] for k in names]); print("features:", len(names), "grid", X.shape, "floodplain pixels", int(fpm.sum()))
np.savez_compressed(f"{OUT}/floodplain_model_features.npz", names=np.array(names), stack=stack.astype("float32"), xs=xs, ys=ys, fpm=fpm)
# ---- sample at the points
jx = np.clip(((F.easting_tied.values - (E0 - HALF)) / RES).astype(int), 0, len(xs) - 1); iy = np.clip(((N0 + HALF - F.northing_tied.values) / RES).astype(int), 0, len(ys) - 1)
P = pd.DataFrame(stack[:, iy, jx].T, columns=names); _dz = P.z_after_m.values - F.gen2_z_m.values; print("CHECK 10 m z_after at the points minus compile gen2_z_m (bilinear 5 m): median %+.2f m, max |.| %.2f m (a mis-indexed mosaic gives tens of m)" % (np.nanmedian(_dz), np.nanmax(np.abs(_dz)))); assert np.nanmax(np.abs(_dz)) < 2.0, "tile layer mis-indexed"
P["r08"] = F.r08_mm.values; P["r21"] = F.r21_mm.values; P["E"] = F.easting_tied.values; P["N"] = F.northing_tied.values; P["rover"] = F.rover.values
P.to_csv(f"{OUT}/floodplain_model_points.csv", index=False)
# ---- blocked CV
BLOCK = 40.0; rng = np.random.default_rng(0)
def block_ids(E, N): return (np.floor((E - E.min()) / BLOCK).astype(int) * 1000 + np.floor((N - N.min()) / BLOCK).astype(int))
def cv_score(model_fn, Xm, y, blocks, reps=20, k=5):
    ub = np.unique(blocks); preds = np.zeros((reps, len(y)))
    for r in range(reps):
        perm = rng.permutation(len(ub)); fold_of = {b: f % k for f, b in zip(range(len(ub)), ub[perm])}; folds = np.array([fold_of[b] for b in blocks])
        for f in range(k):
            te = folds == f; tr = ~te
            if te.sum() == 0: continue
            m = model_fn(); m.fit(Xm[tr], y[tr]); preds[r, te] = m.predict(Xm[te])
    pm = preds.mean(0); ss = np.sum((y - pm) ** 2); st = np.sum((y - y.mean()) ** 2); return 1 - ss / st, nmad(y - pm), pm
results = {}
for target in ("r08", "r21"):
    D = P.dropna(subset=[target]); cols = [c for c in names if D[c].notna().mean() > 0.95]; D = D.dropna(subset=cols); y = D[target].values; Xm = D[cols].values; blocks = block_ids(D.E.values, D.N.values)
    print(f"\n==== TARGET {target}: n {len(D)}, features {len(cols)}, blocks {len(np.unique(blocks))}, target NMAD {nmad(y):.0f} mm, sd {y.std():.0f}")
    print("single features, Spearman and blocked-CV R^2 of a Theil-Sen line, top 12:")
    single = []
    for c in cols:
        rho = spearmanr(D[c], y).correlation
        class TS:
            def fit(s, X_, y_): t = theilslopes(y_, X_[:, 0]); s.a, s.b = t[0], t[1]; return s
            def predict(s, X_): return s.a * X_[:, 0] + s.b
        r2, nm, _ = cv_score(lambda: TS(), D[[c]].values, y, blocks, reps=5); single.append((c, rho, r2, nm))
    single.sort(key=lambda t: -t[2])
    for c, rho, r2, nm in single[:12]: print(f"   {c:28s} rho {rho:+.2f}  CV R2 {r2:+.3f}  CV NMAD {nm:.0f}")
    models = {"ridge": lambda: make_pipeline(StandardScaler(), RidgeCV(alphas=np.logspace(-2, 4, 25))),
              "lasso": lambda: make_pipeline(StandardScaler(), LassoCV(cv=5, max_iter=20000, n_alphas=40)),
              "random_forest": lambda: RandomForestRegressor(300, max_depth=4, min_samples_leaf=5, random_state=0),
              "gradient_boosting": lambda: GradientBoostingRegressor(n_estimators=200, max_depth=2, learning_rate=0.05, subsample=0.8, random_state=0)}
    best = (single[0][0], single[0][2], single[0][3], "single")
    print("multivariate, blocked CV (40 m blocks, 5 folds x 20 repeats):")
    for nm_, fn in models.items():
        r2, nm, pm = cv_score(fn, Xm, y, blocks); print(f"   {nm_:18s} CV R2 {r2:+.3f}  CV NMAD {nm:.0f} mm"); 
        if r2 > best[1]: best = (nm_, r2, nm, "multi")
    print(f"BEST for {target}: {best[0]} (CV R2 {best[1]:+.3f}, CV NMAD {best[2]:.0f} against target NMAD {nmad(y):.0f})")
    # refit best on all points, report coefficients / importances, predict the raster
    if best[3] == "multi":
        m = models[best[0]](); m.fit(Xm, y)
        if best[0] in ("ridge", "lasso"):
            coef = m[-1].coef_; order = np.argsort(-np.abs(coef)); print("   standardized coefficients (mm per sd), top 10:", [(cols[i], round(float(coef[i]), 1)) for i in order[:10]])
        else:
            imp = m.feature_importances_; order = np.argsort(-imp); print("   importances, top 10:", [(cols[i], round(float(imp[i]), 3)) for i in order[:10]])
        full = stack[[names.index(c) for c in cols]].reshape(len(cols), -1).T; ok = np.isfinite(full).all(1); pred = np.full(full.shape[0], np.nan); pred[ok] = m.predict(full[ok]); pred = pred.reshape(X.shape)
    else:
        c = best[0]; t = theilslopes(y, D[c].values); pred = t[0] * feat[c] + t[1]; print(f"   Theil-Sen on {c}: {t[0]:.3g} per unit, intercept {t[1]:.0f}")
    _, _, pm = cv_score(models[best[0]] if best[3] == "multi" else (lambda: TS()), Xm if best[3] == "multi" else D[[best[0]]].values, y, blocks)
    pd.DataFrame(dict(E=D.E, N=D.N, rover=D.rover, observed=y, cv_predicted=pm)).to_csv(f"{OUT}/floodplain_model_{target}_cv.csv", index=False)
    pred_fp = np.where(fpm, pred, np.nan)
    for suffix, arr in (("pred", pred_fp), ("pred_fullbox", pred)):
        with rasterio.open(f"{OUT}/floodplain_model_{target}_{suffix}.tif", "w", driver="GTiff", height=len(ys), width=len(xs), count=1, dtype="float32", crs="EPSG:26915", transform=from_origin(E0 - HALF, N0 + HALF, RES, RES), nodata=np.nan) as dst: dst.write(np.asarray(arr, "float32"), 1)
    vmax = 600 if target == "r08" else 250
    fig, axs = plt.subplots(1, 2, figsize=(13, 6.5), gridspec_kw=dict(width_ratios=[1.15, 1]))
    ax = axs[0]; im = ax.imshow(pred_fp, extent=[xs[0] - RES / 2, xs[-1] + RES / 2, ys[-1] - RES / 2, ys[0] + RES / 2], cmap=cmc.bam_r, vmin=-vmax, vmax=vmax, interpolation="nearest")
    ax.scatter(D.E, D.N, c=y, cmap=cmc.bam_r, vmin=-vmax, vmax=vmax, s=14, edgecolors="k", linewidths=0.3); ax.set_xlim(E0 - 320, E0 + 320); ax.set_ylim(N0 - 320, N0 + 320); ax.set_xticks([]); ax.set_yticks([])
    plt.colorbar(im, ax=ax, shrink=0.7, label=f"{'2008' if target == 'r08' else '2021'} lidar minus RTK (mm)")
    ax = axs[1]; ax.scatter(pm, y, s=12, c=[{"CORS": "k", "RS+ A": "#e69f00", "RS+ B": "#0072b2"}[r] for r in D.rover]); lim = [min(pm.min(), y.min()) - 20, max(pm.max(), y.max()) + 20]; ax.plot(lim, lim, color="#999", lw=0.8); ax.set_xlabel("blocked-CV prediction (mm)"); ax.set_ylabel("observed (mm)"); ax.set_aspect("equal")
    fig.savefig(f"{OUT}/floodplain_model_{target}.png", dpi=140, bbox_inches="tight"); plt.close(fig); print("   wrote", f"{OUT}/floodplain_model_{target}.png", "and _pred.tif")

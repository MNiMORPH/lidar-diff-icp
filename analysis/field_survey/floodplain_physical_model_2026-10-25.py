"""A PHYSICAL model of the 2008 lidar minus RTK offset on the floodplain field (Andy, 2026-11-18: "keep it more physical
than gradient boosting ... find some combination of information from the DEM and Sentinel to explain it").
TARGET r08 = 2008 lidar minus tied RTK at the 126 floodplain points where the 2008 grid has a value (mm; r21 run the
same way as a control). CANDIDATE TERMS, each a named physical quantity on the 10 m grid of floodplain_model_features.npz
(the Sentinel-2 bands and the 5 m pipeline layers averaged to 10 m; the npz is rebuilt by floodplain_residual_model with
the bottom-up tile indexing):
  DEM (2021 ground, z_after): z_plane_m = elevation minus a plane fitted to the floodplain pixels of the box (height above
  the valley-floor trend); z_local_m = elevation minus the median within 150 m over floodplain pixels (swale / ridge);
  slope_deg from the 10 m grid; dist_upland_m = distance to the nearest pixel more than 5 m above the plane (bluff toe).
  Sentinel-2 L2A (surface reflectance / 10000): NDVI on 2026-04-10 (early spring), 2021-05-01 (late spring, the 2021 lidar
  day), 2021-08-14 (summer), 2026-10-17 (autumn); NDRE and NDMI in summer; vis_may = mean of B02-B04 on 2021-05-01 (spring
  brightness, litter); d_ndvi_spring = summer minus late-spring NDVI (late green-up); d_ndvi_autumn = autumn minus summer.
MODEL: ordinary least squares, terms added by FORWARD SELECTION on the spatially blocked CV score (40 m blocks, 5 folds x
20 random assignments, as the data product), stopping when the gain in CV R^2 is below 0.01 or at four terms (both MINE);
every step's CV score is printed. The chosen equation is printed in natural units with OLS standard errors, scored on the
rover-median-removed target too (level is gauge), refit, and mapped over the box (GeoTIFF + PNG) with the points drawn.
OUTPUT data/field/2026-10-25/imagery/floodplain_physical_* . Written incrementally."""
import numpy as np, pandas as pd, rasterio, warnings
from rasterio.transform import from_origin
from scipy import ndimage
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from cmcrameri import cm as cmc
warnings.filterwarnings("ignore")
OUT = "data/field/2026-10-25/imagery"; RES = 10
d = np.load(f"{OUT}/floodplain_model_features.npz"); names = list(d["names"]); st = d["stack"]; xs = d["xs"]; ys = d["ys"]; fpm = d["fpm"]
X, Y = np.meshgrid(xs, ys); g = lambda n: st[names.index(n)].astype(float)
P = pd.read_csv(f"{OUT}/floodplain_model_points.csv"); F = pd.read_csv(f"{OUT}/floodplain_points_covariates.csv")
assert np.nanmax(np.abs(P.z_after_m.values - F.gen2_z_m.values)) < 2.0, "feature stack elevation disagrees with the compile: rebuild the npz"
def nmad(v): v = np.asarray(v, float); v = v[np.isfinite(v)]; return 1.4826 * np.median(np.abs(v - np.median(v)))
# ---- DEM terms
z = g("z_after_m"); fpz = fpm & np.isfinite(z)
A = np.c_[X[fpz], Y[fpz], np.ones(fpz.sum())]; c = np.linalg.lstsq(A, z[fpz], rcond=None)[0]
for _ in range(3):  # robust: refit to pixels within 1 m of the plane (MINE), so the bluff toe and the channel do not tilt it
    r = z[fpz] - A @ c; keep = np.abs(r) < 1.0; c = np.linalg.lstsq(A[keep], z[fpz][keep], rcond=None)[0]
plane = c[0] * X + c[1] * Y + c[2]; z_plane = z - plane
print(f"valley-floor plane over {int(fpz.sum())} floodplain pixels ({int(keep.sum())} within 1 m after 3 robust passes): slope {np.hypot(c[0], c[1])*1000:.2f} m/km towards {np.degrees(np.arctan2(-c[0], -c[1])) % 360:.0f} deg (azimuth of descent)")
zf = np.where(fpz & (np.abs(z_plane) < 5), z, np.nan)
def nanmedian_filter(a, size):
    out = np.full(a.shape, np.nan); h = size // 2
    for i in range(a.shape[0]):
        for j in range(a.shape[1]):
            w = a[max(i - h, 0):i + h + 1, max(j - h, 0):j + h + 1]; w = w[np.isfinite(w)]
            if w.size >= 5: out[i, j] = np.median(w)
    return out
z_local = z - nanmedian_filter(zf, 15)   # 15 cells = 150 m
dzdy, dzdx = np.gradient(z, RES); slope = np.degrees(np.arctan(np.hypot(dzdx, dzdy)))
upland = np.isfinite(z_plane) & (z_plane > 5); dist_upland = ndimage.distance_transform_edt(~upland) * RES
# ---- Sentinel terms
S = {}
S["ndvi_apr26"] = g("s2_2026-04-10_NDVI"); S["ndvi_may21"] = g("s2_2021-05-01_NDVI"); S["ndvi_aug21"] = g("s2_2021-08-14_NDVI"); S["ndvi_oct26"] = g("s2_2026-10-17_NDVI")
S["ndre_aug21"] = g("s2_2021-08-14_NDRE"); S["ndmi_aug21"] = g("s2_2021-08-14_NDMI")
S["vis_may21"] = (g("s2_2021-05-01_B02") + g("s2_2021-05-01_B03") + g("s2_2021-05-01_B04")) / 3 / 10000
S["d_ndvi_spring"] = S["ndvi_aug21"] - S["ndvi_may21"]; S["d_ndvi_autumn"] = S["ndvi_oct26"] - S["ndvi_aug21"]
TERMS = dict(z_plane_m=z_plane, z_local_m=z_local, slope_deg=slope, dist_upland_m=dist_upland, **S)
UNITS = dict(z_plane_m="m", z_local_m="m", slope_deg="deg", dist_upland_m="m", vis_may21="refl"); UNITS.update({k: "index" for k in S if k not in UNITS})
# ---- sample at the points (same cell as the data product)
E0, N0 = F.easting_tied.mean(), F.northing_tied.mean(); HALF = 600
jx = np.clip(((F.easting_tied.values - (E0 - HALF)) / RES).astype(int), 0, len(xs) - 1); iy = np.clip(((N0 + HALF - F.northing_tied.values) / RES).astype(int), 0, len(ys) - 1)
T = pd.DataFrame({k: v[iy, jx] for k, v in TERMS.items()}); T["r08"] = F.r08_mm.values; T["r21"] = F.r21_mm.values; T["E"] = F.easting_tied.values; T["N"] = F.northing_tied.values; T["rover"] = F.rover.values
for c_ in ("r08", "r21"): T[c_ + "_rel"] = T[c_] - T.groupby("rover")[c_].transform("median")
T.to_csv(f"{OUT}/floodplain_physical_points.csv", index=False)
print("COLUMNS floodplain_physical_points.csv:", ", ".join(f"{k} ({UNITS[k]})" for k in TERMS), "; r08/r21 = lidar minus tied RTK (mm); *_rel = minus the rover median")
print("terms at the points (p10 / p50 / p90):"); [print(f"   {k:16s} {np.nanpercentile(T[k], 10):+9.3f} {np.nanpercentile(T[k], 50):+9.3f} {np.nanpercentile(T[k], 90):+9.3f} {UNITS[k]}") for k in TERMS]
# ---- blocked CV OLS forward selection
BLOCK = 40.0; rng = np.random.default_rng(0)
def block_ids(E, N): return (np.floor((E - E.min()) / BLOCK).astype(int) * 1000 + np.floor((N - N.min()) / BLOCK).astype(int))
def ols(Xm, y): A_ = np.c_[Xm, np.ones(len(y))]; b = np.linalg.lstsq(A_, y, rcond=None)[0]; return b
def cv_ols(Xm, y, blocks, reps=20, k=5):
    ub = np.unique(blocks); preds = np.zeros((reps, len(y)))
    for r in range(reps):
        perm = rng.permutation(len(ub)); fold_of = {b: f % k for f, b in zip(range(len(ub)), ub[perm])}; folds = np.array([fold_of[b] for b in blocks])
        for f in range(k):
            te = folds == f; tr = ~te
            if te.sum() == 0: continue
            b = ols(Xm[tr], y[tr]); preds[r, te] = np.c_[Xm[te], np.ones(te.sum())] @ b
    pm = preds.mean(0); return 1 - np.sum((y - pm) ** 2) / np.sum((y - y.mean()) ** 2), nmad(y - pm), pm
for target in ("r08", "r21"):
    D = T.dropna(subset=[target] + list(TERMS)); y = D[target].values; blocks = block_ids(D.E.values, D.N.values); cand = list(TERMS)
    print(f"\n==== TARGET {target}: n {len(D)}, blocks {len(np.unique(blocks))}, target NMAD {nmad(y):.0f} mm, sd {y.std():.0f}")
    print("single terms, Spearman and blocked-CV R^2 of OLS:")
    rows = [(k, D[k].corr(D[target], method="spearman"), *cv_ols(D[[k]].values, y, blocks, reps=5)[:2]) for k in cand]
    print("   within-rover Spearman (T2a RS+ B n %d, T1a RS+ A n %d):" % ((D.rover == "RS+ B").sum(), (D.rover == "RS+ A").sum()), {k: tuple(round(float(D[D.rover == r][k].corr(D[D.rover == r][target], method="spearman")), 2) for r in ("RS+ B", "RS+ A")) for k in cand})
    for k, rho, r2, nm in sorted(rows, key=lambda t: -t[2]): print(f"   {k:16s} rho {rho:+.2f}  CV R2 {r2:+.3f}  CV NMAD {nm:.0f}")
    def forward(y_, label):
        chosen = []; best_r2 = -np.inf
        while len(chosen) < 4:
            trial = [(k, *cv_ols(D[chosen + [k]].values, y_, blocks)[:2]) for k in cand if k not in chosen]; k, r2, nm = max(trial, key=lambda t: t[1])
            print(f"   [{label}] step {len(chosen)+1}: add {k:16s} -> CV R2 {r2:+.3f}  CV NMAD {nm:.0f}" + ("" if r2 - best_r2 >= 0.01 else "   (gain < 0.01, stop)"))
            if r2 - best_r2 < 0.01: break
            chosen.append(k); best_r2 = r2
        return chosen
    yr = D[target + "_rel"].values; chosen_rel = forward(yr, "rover-median removed"); br = ols(D[chosen_rel].values, yr) if chosen_rel else None
    if chosen_rel: print(f"   within-rover structure equation: {target}_rel (mm) = {br[-1]:+.0f}" + "".join(f" {br[i]:+.4g} x {k} [{UNITS[k]}]" for i, k in enumerate(chosen_rel)), "; blocked-CV R2 %+.3f" % cv_ols(D[chosen_rel].values, yr, blocks)[0])
    chosen = forward(y, "raw target")
    Xm = D[chosen].values; b = ols(Xm, y); A_ = np.c_[Xm, np.ones(len(y))]; res = y - A_ @ b; dof = len(y) - A_.shape[1]; cov = np.linalg.inv(A_.T @ A_) * (res @ res) / dof; se = np.sqrt(np.diag(cov))
    r2cv, nmcv, pm = cv_ols(Xm, y, blocks); r2fit = 1 - np.sum(res ** 2) / np.sum((y - y.mean()) ** 2)
    print(f"EQUATION for {target} (OLS, n {len(y)}): {target} (mm) = {b[-1]:+.0f} (se {se[-1]:.0f})" + "".join(f" {b[i]:+.4g} (se {se[i]:.2g}) x {k} [{UNITS[k]}]" for i, k in enumerate(chosen)))
    print(f"   fit R2 {r2fit:.3f}, blocked-CV R2 {r2cv:+.3f}, CV NMAD {nmcv:.0f} mm against target NMAD {nmad(y):.0f}; residual by rover (median, NMAD):", {r: (round(float(np.median(res[D.rover == r]))), round(float(nmad(res[D.rover == r])))) for r in sorted(set(D.rover))})
    print("   effect of each term across its p10-p90 at the points (mm):", {k: round(float(b[i] * (np.percentile(D[k], 90) - np.percentile(D[k], 10)))) for i, k in enumerate(chosen)})
    r2rel, nmrel, _ = cv_ols(Xm, yr, blocks); print(f"   same terms on the rover-median-removed target: blocked-CV R2 {r2rel:+.3f}, CV NMAD {nmrel:.0f} against {nmad(yr):.0f}")
    pd.DataFrame(dict(E=D.E, N=D.N, rover=D.rover, observed=y, fitted=A_ @ b, cv_predicted=pm)).to_csv(f"{OUT}/floodplain_physical_{target}_cv.csv", index=False)
    pred = b[-1] + sum(b[i] * TERMS[k] for i, k in enumerate(chosen)); pred_fp = np.where(fpm, pred, np.nan)
    for suffix, arr in (("pred", pred_fp), ("pred_fullbox", pred)):
        with rasterio.open(f"{OUT}/floodplain_physical_{target}_{suffix}.tif", "w", driver="GTiff", height=len(ys), width=len(xs), count=1, dtype="float32", crs="EPSG:26915", transform=from_origin(E0 - HALF, N0 + HALF, RES, RES), nodata=np.nan) as dst: dst.write(np.asarray(arr, "float32"), 1)
    vmax = 600 if target == "r08" else 250; cols = {"CORS": "k", "RS+ A": "#e69f00", "RS+ B": "#0072b2"}
    fig, axs = plt.subplots(1, 2 + len(chosen), figsize=(5.5 * (2 + len(chosen)), 5.2), gridspec_kw=dict(width_ratios=[1.3] + [1] * (1 + len(chosen))))
    ax = axs[0]; im = ax.imshow(pred_fp, extent=[xs[0] - RES / 2, xs[-1] + RES / 2, ys[-1] - RES / 2, ys[0] + RES / 2], cmap=cmc.bam_r, vmin=-vmax, vmax=vmax, interpolation="nearest")
    ax.scatter(D.E, D.N, c=y, cmap=cmc.bam_r, vmin=-vmax, vmax=vmax, s=14, edgecolors="k", linewidths=0.3); ax.set_xlim(E0 - 320, E0 + 320); ax.set_ylim(N0 - 320, N0 + 320); ax.set_xticks([]); ax.set_yticks([])
    plt.colorbar(im, ax=ax, shrink=0.7, label=f"{'2008' if target == 'r08' else '2021'} lidar minus RTK (mm)")
    ax = axs[1]; ax.scatter(pm, y, s=12, c=[cols[r] for r in D.rover]); lim = [min(pm.min(), y.min()) - 20, max(pm.max(), y.max()) + 20]; ax.plot(lim, lim, color="#999", lw=0.8); ax.set_xlabel("blocked-CV prediction (mm)"); ax.set_ylabel("observed (mm)"); ax.set_aspect("equal")
    for i, k in enumerate(chosen):   # partial residual: y minus the other terms, against this term, with the fitted line
        ax = axs[2 + i]; part = res + b[i] * D[k].values; ax.scatter(D[k], part, s=12, c=[cols[r] for r in D.rover]); xx = np.linspace(D[k].min(), D[k].max(), 2); ax.plot(xx, b[i] * xx, color="#999", lw=0.8); ax.set_xlabel(f"{k} ({UNITS[k]})"); ax.set_ylabel(f"{target} partial residual (mm)")
    for r, c_ in cols.items(): axs[1].scatter([], [], c=c_, s=12, label=r)
    axs[1].legend(frameon=False, fontsize=8)
    fig.savefig(f"{OUT}/floodplain_physical_{target}.png", dpi=140, bbox_inches="tight"); plt.close(fig); print("   wrote", f"{OUT}/floodplain_physical_{target}.png", "_pred.tif, _pred_fullbox.tif, _cv.csv")

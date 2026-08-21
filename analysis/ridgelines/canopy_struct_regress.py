"""Regress the forest-floor measurement offset (dz/dt, mm/yr) on each canopy
STRUCTURE metric, individually and jointly.

Offset = dod_osm / 12.44 yr * 1000 (mm/yr), on forest_region cells with finite
dod_osm. For each metric report per-cell R^2 and binned-median R^2 (~10 quantile
bins, median offset per bin, line fit, R^2 on binned medians). Then a multiple
regression on the best 2-3 metrics (5-fold CV per-cell R^2 and binned R^2).

Run:
  cd /home/awickert/projects/lidar-diff-icp
  env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python analysis/ridgelines/canopy_struct_regress.py
"""
import numpy as np

DERIVED = "data/derived/elba_fulldensity"
DT_YR = 12.44

cs = np.load(f"{DERIVED}/canopy_struct.npz")
forest = np.load(f"{DERIVED}/forest_region.npy")
dod = np.load("data/derived/elba_refdatum/dod_osm.npy")
# reference metrics for comparison
canopy_cover = np.load(f"{DERIVED}/canopy_cover.npy")
penetration = np.load(f"{DERIVED}/penetration.npy")

offset = dod / DT_YR * 1000.0  # mm/yr

metrics = {
    "canopy_height_p95": cs["canopy_height_p95"],
    "understory_frac": cs["understory_frac"],
    "midstory_frac": cs["midstory_frac"],
    "veg_frac": cs["veg_frac"],
    "ground_return_density": cs["ground_return_density"],
    "low_gap": cs["low_gap"],
    "canopy_cover(ref)": canopy_cover,
    "penetration(ref)": penetration,
}

base = forest & np.isfinite(offset)


def r2_line(x, y):
    """R^2 of an OLS line y ~ a + b x."""
    b, a = np.polyfit(x, y, 1)
    yhat = a + b * x
    ss_res = np.sum((y - yhat) ** 2)
    ss_tot = np.sum((y - y.mean()) ** 2)
    return 1 - ss_res / ss_tot, b


def binned_r2(x, y, nbin=10):
    q = np.quantile(x, np.linspace(0, 1, nbin + 1))
    q[-1] += 1e-9
    # unique edges (guard against ties in low-variance predictors)
    q = np.unique(q)
    bidx = np.clip(np.digitize(x, q[1:-1]), 0, len(q) - 2)
    bx, by = [], []
    for k in range(len(q) - 1):
        m = bidx == k
        if m.sum() >= 20:
            bx.append(np.median(x[m]))
            by.append(np.median(y[m]))
    bx = np.array(bx); by = np.array(by)
    if len(bx) < 3:
        return np.nan, np.nan, len(bx)
    r2, b = r2_line(bx, by)
    return r2, b, len(bx)


print(f"forest cells with finite offset: {base.sum():,}")
print(f"offset (mm/yr): mean {offset[base].mean():.2f}  std {offset[base].std():.2f}")
print()
print(f"{'metric':<22} {'n':>7} {'percell_R2':>11} {'binned_R2':>10} {'slope_sign':>11} {'nbins':>6}")
print("-" * 72)

results = {}
for name, arr in metrics.items():
    m = base & np.isfinite(arr)
    x = arr[m].astype(float)
    y = offset[m].astype(float)
    if x.std() == 0:
        continue
    pc_r2, pc_b = r2_line(x, y)
    bn_r2, bn_b, nb = binned_r2(x, y)
    sign = "+" if bn_b > 0 else "-"
    results[name] = dict(percell=pc_r2, binned=bn_r2, slope=bn_b, n=m.sum())
    print(f"{name:<22} {m.sum():>7} {pc_r2:>11.4f} {bn_r2:>10.4f} {sign:>11} {nb:>6}")

# ---- multiple regression on best structure metrics ---------------------------
print()
struct_only = {k: v for k, v in results.items() if "(ref)" not in k}
best = sorted(struct_only, key=lambda k: -struct_only[k]["binned"])[:3]
print("best 3 structure metrics by binned R^2:", best)


def cv_r2(X, y, k=5, seed=0):
    rng = np.random.default_rng(seed)
    idx = rng.permutation(len(y))
    folds = np.array_split(idx, k)
    ss_res = ss_tot = 0.0
    ybar = y.mean()
    for i in range(k):
        te = folds[i]
        tr = np.concatenate([folds[j] for j in range(k) if j != i])
        A = np.column_stack([np.ones(len(tr)), X[tr]])
        coef, *_ = np.linalg.lstsq(A, y[tr], rcond=None)
        Ate = np.column_stack([np.ones(len(te)), X[te]])
        yhat = Ate @ coef
        ss_res += np.sum((y[te] - yhat) ** 2)
        ss_tot += np.sum((y[te] - ybar) ** 2)
    return 1 - ss_res / ss_tot


for combo in (best[:2], best[:3]):
    m = base.copy()
    for c in combo:
        m &= np.isfinite(metrics[c])
    X = np.column_stack([metrics[c][m].astype(float) for c in combo])
    y = offset[m].astype(float)
    # in-sample multiple-R^2
    A = np.column_stack([np.ones(len(y)), X])
    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    yhat = A @ coef
    r2_in = 1 - np.sum((y - yhat) ** 2) / np.sum((y - y.mean()) ** 2)
    r2cv = cv_r2(X, y)
    print(f"multi {combo}: in-sample R2={r2_in:.4f}  5foldCV R2={r2cv:.4f}  n={m.sum():,}")

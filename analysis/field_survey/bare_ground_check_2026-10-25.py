#!/usr/bin/env python
"""Is the 2021 lidar's +0.13 m over the floodplain a registration error of the flight, or vegetation?
(Andy, 2026-11-20: "check the points that we took early on, in the open parking area, compared to the lidar.
Make sure that our points are good since one student set up his RTK GNSS as a false base station.")

What this script does, in order, all from the raw exports and the point clouds:
  A. rebuilds every logged point with the base it was corrected from (B1 = the student's receiver, lat 44.11685493;
     B = the class base, lat 44.11721583; CORS), positions in EPSG:6344, and writes them;
  B. the B1-era points (13:20-13:28) carry B1's unknown autonomous-position error, so they are compared with the
     lidar only INTERNALLY: lot-edge points against the lower-bluff points in the same B1 frame (the error cancels);
  C. every 2021 return within 3 m and 1 m of chosen tied points (benchmark tie, rtk_points_all.csv), with a plane
     fitted to the class-2 returns and evaluated at the point, plus the 2008 returns (star frame) where the
     4342-29-64 tile covers the point: the monument DR9925 (flush disk on a mown roadside), the first points of
     each floodplain rover near the lot, and the B1-era lot points (lidar structure only);
  D. the acquisition's own held-out checkpoints (src/lidar_diff_icp/groundtruth/data/mn_se_driftless_2021_control.csv):
     NVA = bare, VVA = vegetated; error = surveyed minus lidar, so NEGATIVE means the lidar reads HIGH;
  E. a map of the lot area: NAIP 2021 (Planetary Computer, skipped if offline), 2021 intensity, class-2 fraction and
     return-height spread at 1 m, with the points.
Outputs: data/field/2026-10-25/bare_ground/{raw_points_with_base.csv, returns_at_points.csv, lot_map.png, log.txt}.
Cuts that are MINE: 3 m / 1 m radii, lot = within 15 m of E 579979 N 4885340, lower bluff = slope < 10 deg.
"""
import os, re, json, sys
import numpy as np, pandas as pd, laspy
from pyproj import Transformer
ROOT = os.path.dirname(os.path.abspath(__file__)) + "/../.."; os.chdir(ROOT)
F = "data/field/2026-10-25"; RAW = f"{F}/raw"; OUT = f"{F}/bare_ground"; os.makedirs(OUT, exist_ok=True)
LOG = open(f"{OUT}/log.txt", "w")
def say(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); LOG.write(s + "\n"); LOG.flush()
to_utm = Transformer.from_crs("EPSG:6318", "EPSG:6344", always_xy=True)
LOT = (579979.0, 4885340.0); B1_LAT, B_LAT = 44.11685493, 44.11721583
say("inputs: raw Emlid exports in", RAW, "| tied table", f"{F}/rtk_points_all.csv (benchmark tie, no fitted shift)",
    "| 2021 clouds data/after/valley/{4342-29-64,4358-29-01}_3dep_full.laz | 2008 cloud data/before/4342-29-64.laz + star shifts",
    "| checkpoints src/lidar_diff_icp/groundtruth/data/mn_se_driftless_2021_control.csv")

# ---------------------------------------------------------------- A. raw points with their base
def load(path, group):
    d = pd.read_excel(path) if path.endswith(".xlsx") else pd.read_csv(path, encoding="utf-8-sig")
    d = d.dropna(subset=["Longitude"]); rows = []
    for _, r in d.iterrows():
        name = r["Name"] if "Name" in r else r["Point"]; t = re.search(r"\d\d:\d\d:\d\d", str(r["Averaging start"]))
        h = r["Ellipsoidal height"] if "Ellipsoidal height" in r else r["Ellipsoidal height (m)"]
        rows.append(dict(group=group, point=int(float(name)), t=t[0] if t else "", status=r["Solution status"], lon=float(r["Longitude"]), lat=float(r["Latitude"]),
                         h_ell_m=float(h), base_lat=pd.to_numeric(r.get("Base latitude", np.nan), errors="coerce"), baseline_m=pd.to_numeric(r.get("Baseline", np.nan), errors="coerce"),
                         rms_v=pd.to_numeric(r.get("Elevation RMS", np.nan), errors="coerce")))
    return rows
rows = []
rows += load(f"{RAW}/rover_RS3_8243D2AFB94BDAE3_hillslope_T3a_T3b_T3c.csv", "hillslope_RS3")
rows += load(f"{RAW}/rover_RS3_82434D3F6234D057_hillslope_T4a_T4b.csv", "hillslope_T4a")
rows += load(f"{RAW}/rover_RSplus_824350D44A583F7B_floodplain_T2a.csv", "floodplain_RSplus_B")
rows += load(f"{RAW}/rover_RSplus_8243CA7256E53180_floodplain_T1a_T1b.xlsx", "floodplain_RSplus_A")
rows += load(f"{RAW}/rover_CORS_Wickert_Whitewater_floodplain_LiDAR_veg.csv", "floodplain_CORS_Wickert")
A = pd.DataFrame(rows); A["E"], A["N"] = to_utm.transform(A.lon.values, A.lat.values)
A["base"] = np.where(np.isclose(A.base_lat, B1_LAT, rtol=0, atol=1e-6), "B1", np.where(np.isclose(A.base_lat, B_LAT, rtol=0, atol=1e-6), "B", np.where(A.group == "floodplain_CORS_Wickert", "CORS", "none")))
P = pd.read_csv(f"{F}/rtk_points_all.csv"); GEOID_N = float((P.h_ell_m - P.navd88_m).median())
A["navd88_m"] = A.h_ell_m - GEOID_N
m = json.load(open("data/derived/4342-29-64/meta.json")); b = m["bounds"]; res = m["res"]
Z = np.load("data/derived/4342-29-64/z_after.npy", mmap_mode="r"); SL = np.load("data/derived/4342-29-64/slope.npy", mmap_mode="r"); LV = np.load("data/derived/4342-29-64/lowveg_grid.npy", mmap_mode="r")
def samp(a, x, y):
    i = int((y - b[1]) / res); j = int((x - b[0]) / res)     # arrays are bottom-up (row 0 = south)
    return float(a[i, j]) if 0 <= i < a.shape[0] and 0 <= j < a.shape[1] else np.nan
A["z21_grid_m"] = [samp(Z, x, y) for x, y in zip(A.E, A.N)]; A["slope_deg"] = [samp(SL, x, y) for x, y in zip(A.E, A.N)]; A["lowveg_gen2"] = [samp(LV, x, y) for x, y in zip(A.E, A.N)]
A["d_lot_m"] = np.hypot(A.E - LOT[0], A.N - LOT[1])
A.to_csv(f"{OUT}/raw_points_with_base.csv", index=False, float_format="%.4f")
say(f"A. {len(A)} logged points; base counts {A.base.value_counts().to_dict()} (none = export lost the base columns); geoid N = h_ell - navd88 from the compile table, {GEOID_N:+.3f} m")
b1 = A[A.base == "B1"].sort_values("t")
say("   B1-era points (the false base, lat 44.11685493, self-reported E 579986.1 N 4885337.5, i.e. AT the lot):")
say(b1[["group", "point", "t", "status", "E", "N", "baseline_m", "slope_deg", "lowveg_gen2", "d_lot_m", "navd88_m", "z21_grid_m"]].round(3).to_string())

# ---------------------------------------------------------------- B. B1 internal test
b1 = b1[b1.status == "FIX"].copy(); b1["rtk_minus_z21"] = b1.navd88_m - b1.z21_grid_m
lot = b1[b1.d_lot_m < 15]; foot = b1[(b1.d_lot_m >= 15) & (b1.slope_deg < 10)]
say(f"B. B1 frame: RTK minus 2021 grid, lot-edge points (d_lot < 15 m, MINE) median {lot.rtk_minus_z21.median():+.3f} m (n {len(lot)}), lower bluff (slope < 10 deg, MINE) {foot.rtk_minus_z21.median():+.3f} (n {len(foot)}); "
    f"foot minus lot {foot.rtk_minus_z21.median() - lot.rtk_minus_z21.median():+.3f} m (B1's offset cancels). The B1 frame itself sits {b1.rtk_minus_z21.median():+.2f} m from the lidar: B1's autonomous position is wrong by meters, so these points have no absolute use.")

# ---------------------------------------------------------------- C. returns at points
def sel(group, pts, lab):
    out = []
    for q in pts:
        r = P[(P.group == group) & (P.point == q)]
        if len(r): r = r.iloc[0]; out.append(dict(label=f"{lab}{q}", frame="B tie" if r.base == "B" else "CORS", E=r.easting_tied, N=r.northing_tied, z=r.navd88_tied_m, status=r.status, lowveg=r.lowveg_gen2, slope=r.slope_deg, product21=r.gen2_minus_rtk_m))
    return out
pts = sel("benchmark_DR9925", [1, 2], "DR9925_") + sel("floodplain_RSplus_A", [2, 3, 4, 5, 6], "RS+A_") + sel("floodplain_RSplus_B", [1, 2, 3, 17], "RS+B_") + sel("floodplain_CORS_Wickert", [1, 2, 3, 4, 5, 6, 7], "CORS_") + sel("hillslope_T4a", [9, 10, 11], "RS3-T4_")
for _, r in A[(A.base == "B1") & (A.d_lot_m < 15) & (A.status == "FIX")].iterrows():
    pts.append(dict(label=f"B1lot_{r.group[:9]}{r.point}", frame="B1 (meters off)", E=r.E, N=r.N, z=r.navd88_m, status=r.status, lowveg=r.lowveg_gen2, slope=r.slope_deg, product21=np.nan))
pts = pd.DataFrame(pts)
def box(path, bbox, cols):
    parts = []
    with laspy.open(path) as r:
        for ch in r.chunk_iterator(5_000_000):
            x = np.asarray(ch.x); y = np.asarray(ch.y); mm = (x >= bbox[0]) & (x <= bbox[1]) & (y >= bbox[2]) & (y <= bbox[3])
            if mm.any():
                d = dict(x=x[mm], y=y[mm], z=np.asarray(ch.z)[mm]); d.update({k: np.asarray(ch[k])[mm] for k in cols}); parts.append(pd.DataFrame(d))
    return pd.concat(parts) if parts else pd.DataFrame(columns=["x", "y", "z"] + cols)
PAD = 4.0; bbox = (pts.E.min() - PAD, pts.E.max() + PAD, pts.N.min() - PAD, pts.N.max() + PAD)
g2 = pd.concat([box(f"data/after/valley/{t}_3dep_full.laz", bbox, ["classification", "number_of_returns", "intensity"]) for t in ("4342-29-64", "4358-29-01")])
g1 = box("data/before/4342-29-64.laz", bbox, ["point_source_id"]); star = {int(k): v for k, v in json.load(open("data/derived/4342-29-64/corrections.json"))["per_swath_star_frame_dxdydz_m"].items()}
if len(g1): g1["z"] = g1.z + np.array([star[p][2] for p in g1.point_source_id])
def plane_at(s, E, N):
    if len(s) < 6: return np.nan
    M = np.c_[s.x - E, s.y - N, np.ones(len(s))]; c, *_ = np.linalg.lstsq(M, s.z, rcond=None); return float(c[2])
recs = []
for _, p in pts.iterrows():
    d = np.hypot(g2.x - p.E, g2.y - p.N); s3 = g2[d <= 3]; s1 = g2[d <= 1]; c3 = s3[s3.classification == 2]; c1 = s1[s1.classification == 2]
    d1 = np.hypot(g1.x - p.E, g1.y - p.N) if len(g1) else np.array([]); o3 = g1[d1 <= 3] if len(g1) else g1
    recs.append(dict(label=p.label, frame=p.frame, status=p.status, slope_deg=p.slope, lowveg_gen2=p.lowveg, product_gen2_minus_rtk=p.product21,
                     n21_3m=len(s3), min21_3m=s3.z.min() - p.z, p10_21_3m=s3.z.quantile(.1) - p.z, med21_3m=s3.z.median() - p.z, p90_21_3m=s3.z.quantile(.9) - p.z,
                     cls2_share_3m=(s3.classification == 2).mean(), single_share_3m=(s3.number_of_returns == 1).mean(), intensity_med_3m=s3.intensity.median(),
                     n21_1m=len(s1), min21_1m=(s1.z.min() - p.z) if len(s1) else np.nan, med21_1m=(s1.z.median() - p.z) if len(s1) else np.nan,
                     plane_cls2_3m=plane_at(c3, p.E, p.N) - p.z, plane_cls2_1m=plane_at(c1, p.E, p.N) - p.z,
                     n08_3m=len(o3), med08_3m=(o3.z.median() - p.z) if len(o3) else np.nan, min08_3m=(o3.z.min() - p.z) if len(o3) else np.nan))
R = pd.DataFrame(recs); R.to_csv(f"{OUT}/returns_at_points.csv", index=False, float_format="%.4f")
say("C. returns about each point, m above the point's elevation (benchmark tie; the B1 rows are in B1's frame and only their lidar structure means anything):")
pd.set_option("display.width", 300); pd.set_option("display.max_columns", 40)
say(R[["label", "frame", "slope_deg", "lowveg_gen2", "product_gen2_minus_rtk", "n21_3m", "min21_3m", "p10_21_3m", "med21_3m", "p90_21_3m", "cls2_share_3m", "single_share_3m", "n21_1m", "min21_1m", "med21_1m", "plane_cls2_3m", "plane_cls2_1m", "n08_3m", "med08_3m", "min08_3m"]].round(3).to_string())

# ---------------------------------------------------------------- D. the acquisition's own checkpoints
C = pd.read_csv("src/lidar_diff_icp/groundtruth/data/mn_se_driftless_2021_control.csv"); C["dist_km"] = np.hypot(C.easting - 579950, C.northing - 4885400) / 1000
say("D. MN_SE_Driftless_2021 held-out checkpoints, error = surveyed NAVD88 minus delivered lidar (negative = lidar HIGH), whole project:")
for t, g in C[C.point_type.isin(["NVA", "VVA"])].groupby("point_type"):
    say(f"   {t} ({'bare, non-vegetated' if t == 'NVA' else 'vegetated'}): n {g.usgs_ql1_laz_error_m.notna().sum()}, LAZ-surface error median {g.usgs_ql1_laz_error_m.median():+.4f} mean {g.usgs_ql1_laz_error_m.mean():+.4f} m; DEM error median {g.usgs_ql1_dem_error_m.median():+.4f} m; nearest to the site {g.dist_km.min():.1f} km, {int((g.dist_km < 10).sum())} within 10 km")
near = C[(C.point_type != "LCP") & (C.dist_km < 12)].sort_values("dist_km")
say("   within 12 km: " + "; ".join(f"{r.point_id} {r.point_type} {r.dist_km:.1f} km error {r.usgs_ql1_laz_error_m:+.3f}" for _, r in near.iterrows()))
D8 = pd.read_csv("src/lidar_diff_icp/groundtruth/data/mn_dnr_2008_control_semn.csv"); D8["dist_km"] = np.hypot(D8.easting - 579950, D8.northing - 4885400) / 1000; r = D8.sort_values("dist_km").iloc[0]
say(f"   aside, 2008: the DNR validation point {r.point_id} (class {r.point_type}) is {1000 * r.dist_km:.0f} m from the site at E {r.easting:.1f} N {r.northing:.1f}; 2008 surface minus surveyed {-r.dnr_error_m:+.3f} m")

# ---------------------------------------------------------------- E. lot map
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
E0, E1, N0, N1 = 579840, 580110, 4885270, 4885480
g = pd.concat([box(f"data/after/valley/{t}_3dep_full.laz", (E0, E1, N0, N1), ["classification", "intensity"]) for t in ("4342-29-64", "4358-29-01")])
nx, ny = E1 - E0, N1 - N0; cell = ((g.y - N0).astype(int).clip(0, ny - 1)) * nx + (g.x - E0).astype(int).clip(0, nx - 1)
agg = pd.DataFrame(dict(cell=cell.values, z=g.z.values, c=g.classification.values, i=g.intensity.values)).groupby("cell").agg(zmin=("z", "min"), zmax=("z", "max"), cls2=("c", lambda v: (v == 2).mean()), inten=("i", "median"))
def ras(col):
    a = np.full(nx * ny, np.nan); a[agg.index.values] = agg[col].values; return a.reshape(ny, nx)
naip = None
try:
    import planetary_computer as pc, rasterio
    from pystac_client import Client
    from rasterio.warp import transform_bounds
    from rasterio.windows import from_bounds
    tr = Transformer.from_crs("EPSG:6344", "EPSG:4326", always_xy=True); lon0, lat0 = tr.transform(E0, N0); lon1, lat1 = tr.transform(E1, N1)
    cat = Client.open("https://planetarycomputer.microsoft.com/api/stac/v1", modifier=pc.sign_inplace)
    it = [i for i in cat.search(collections=["naip"], bbox=[lon0, lat0, lon1, lat1]).items() if i.datetime.year == 2021][0]
    with rasterio.open(it.assets["image"].href) as rr:
        w = from_bounds(*transform_bounds("EPSG:6344", rr.crs, E0, N0, E1, N1), transform=rr.transform); naip = np.moveaxis(rr.read([1, 2, 3], window=w), 0, -1)
    say(f"E. NAIP {it.id} ({it.datetime.date()}) fetched for the map")
except Exception as e:
    say("E. NAIP not fetched:", repr(e)[:120])
ext = (E0, E1, N0, N1); panels = [("2021 intensity, median per 1 m", ras("inten"), "gray", None), ("2021 class-2 share per 1 m", ras("cls2"), "viridis", (0, 1)), ("2021 return-height spread per 1 m (m, clipped at 1)", np.clip(ras("zmax") - ras("zmin"), 0, 1), "magma", (0, 1))]
fig, axs = plt.subplots(1, 3 + (naip is not None), figsize=(6 * (3 + (naip is not None)), 6.2)); k = 0
if naip is not None: axs[0].imshow(naip, extent=ext); axs[0].set_title("NAIP 2021-08-15"); k = 1
for j, (t, a, cm, cl) in enumerate(panels):
    im = axs[k + j].imshow(a, origin="lower", extent=ext, cmap=cm, vmin=None if cl is None else cl[0], vmax=None if cl is None else cl[1]); axs[k + j].set_title(t, fontsize=10); plt.colorbar(im, ax=axs[k + j], shrink=0.7)
for ax in axs:
    q = A[A.base == "B1"]; ax.plot(q.E, q.N, "s", ms=5, mfc="none", mec="red", label="B1-era points (raw positions)")
    for gname, c, l in [("floodplain_RSplus_A", "#e69f00", "RS+ A"), ("floodplain_RSplus_B", "#0072b2", "RS+ B"), ("floodplain_CORS_Wickert", "k", "CORS"), ("hillslope_T4a", "lime", "RS3 (T4)")]:
        q = P[P.group == gname]; ax.plot(q.easting_tied, q.northing_tied, "o", ms=3, color=c, label=l)
    ax.plot([579986.141], [4885337.535], "r*", ms=12, label="B1 base, self-reported"); ax.plot([LOT[0]], [LOT[1]], "w+", ms=12, mew=2, label="lot center used here")
    ax.set_xlim(E0, E1); ax.set_ylim(N0, N1); ax.set_aspect("equal")
axs[0].legend(fontsize=7, loc="upper left")
for _, r in P[(P.group == "floodplain_RSplus_A") & (P.point <= 6)].iterrows(): axs[0].annotate(f"A{int(r.point)}", (r.easting_tied, r.northing_tied), fontsize=7, color="#e69f00")
for _, r in P[(P.group == "floodplain_CORS_Wickert") & (P.point <= 7)].iterrows(): axs[0].annotate(f"C{int(r.point)}", (r.easting_tied, r.northing_tied), fontsize=7)
fig.tight_layout(); fig.savefig(f"{OUT}/lot_map.png", dpi=int(os.environ.get("FIGDPI", 110))); say(f"   wrote {OUT}/lot_map.png")
LOG.close()

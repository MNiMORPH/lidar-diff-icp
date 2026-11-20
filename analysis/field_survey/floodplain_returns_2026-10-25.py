#!/usr/bin/env python
"""Every lidar return within 1, 3 and 5 m of each point of one floodplain rover, both epochs, as offsets from the
tied RTK elevation of that point; and the PDF of the 3 m returns at one representative point.

Usage  GROUP=floodplain_RSplus_B|floodplain_RSplus_A|floodplain_CORS_Wickert  python floodplain_returns_2026-10-25.py [point]
       (default GROUP floodplain_RSplus_B; the RS+ B outputs keep the names of the first, RS+ B-only version of this
       script, rsplus_b_returns_2026-10-25.py, and hold the same values plus a `tile` column and, in 2008, the 13 returns
       of tile 4358-29-01 within 5 m of points 98-100: the two 2008 LAZ tiles overlap by 42 m in easting, 579993-580035,
       and in that strip 4342-29-64 holds 1,295 returns against 4358-29-01's 5,770, so the first version, reading one
       tile, under-sampled the east edge.)

Inputs
  data/field/2026-10-25/rtk_points_all.csv            benchmark-tied RTK points (env RTK_TABLE overrides; NOT the fitted table:
                                                      it carries the 0.13 m whole-set shift and would put the 2021 lidar on zero)
  data/before/<tile>.laz                              2008 MN DNR, ALL returns, vendor NAVD88 (GEOID03), raw frame
  data/after/valley/<tile>_3dep_full.laz              2021 3DEP, ALL returns (classes 1, 2, 7)
  data/derived/<tile>/corrections.json                per-swath star-frame shifts ADDED to 2008 x, y, z by the product (swath_frame
                                                      "star"), lateral shift ADDED to x, y
  for the tiles 4342-29-64 and 4358-29-01, which between them hold every floodplain point (both clouds are read in the
  points' bounding box; the 2008 star shift is looked up per tile and swath).

Outputs (data/field/2026-10-25/returns/, <tag> = rsplus_b | rsplus_a | cors)
  <tag>_returns_2008.csv, <tag>_returns_2021.csv.gz
      one row per (point, return) pair with dist <= 5 m; a return near two points appears once per point. Columns: point,
      rtk_z_m, x, y, z, dist_m, dz_m (= z - rtk_z_m), return_number, number_of_returns, classification, intensity,
      scan_angle_rank, point_source_id, gps_time, tile; 2008 also x_star, y_star, z_star (the product's frame), dist_star_m,
      dz_star_m. radius_1m/3m/5m booleans from dist_m (dist_star_m for 2008). ALL statuses of the rover are included;
      filter on status downstream (rtk_points_all.csv has it).
  <tag>_returns_summary.csv   per point: status, n and median and lowest dz within each radius, each epoch
  <tag>_returns_log.txt
  figures/rtk_2026-10-25_<TAG>_returns_pdf_p<point>.png   PDF of the 3 m returns, both flights

The RTK "surface" at a point is that point's own tied elevation (floodplain points are meters apart on slopes under 1 deg).
"""
import json, os, sys, time
import numpy as np, pandas as pd, laspy

ROOT = os.path.dirname(os.path.abspath(__file__)) + "/../.."
os.chdir(ROOT)
GROUP = os.environ.get("GROUP", "floodplain_RSplus_B")
TAGS = {"floodplain_RSplus_B": ("rsplus_b", "rsplusB", "RS+ B"), "floodplain_RSplus_A": ("rsplus_a", "rsplusA", "RS+ A"), "floodplain_CORS_Wickert": ("cors", "cors", "CORS")}
tag, TAG, NAME = TAGS[GROUP]
OUT = "data/field/2026-10-25/returns"; os.makedirs(OUT, exist_ok=True); os.makedirs("figures", exist_ok=True)
RTK = os.environ.get("RTK_TABLE", "data/field/2026-10-25/rtk_points_all.csv")
TILES = ["4342-29-64", "4358-29-01"]
RADII = (1.0, 3.0, 5.0); RMAX = max(RADII)
LOG = open(f"{OUT}/{tag}_returns_log.txt", "w")
def say(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); LOG.write(s + "\n"); LOG.flush()

say("inputs:", RTK, "| 2008 data/before/<tile>.laz (all returns, raw vendor frame) | 2021 data/after/valley/<tile>_3dep_full.laz (all returns) | tiles", TILES)
P = pd.read_csv(RTK); B = P[P.group == GROUP].sort_values("point").reset_index(drop=True)
say(f"{NAME} points: {len(B)} (status {B.status.value_counts().to_dict()}); device {B.device.unique().tolist()}; tiles at the points {B.tile.unique().tolist()}")
E = B.easting_tied.to_numpy(float); N = B.northing_tied.to_numpy(float); Z = B.navd88_tied_m.to_numpy(float)
PT = B.point.to_numpy(int)
bbox = (E.min() - RMAX, E.max() + RMAX, N.min() - RMAX, N.max() + RMAX)
say(f"bbox with {RMAX} m pad: E {bbox[0]:.1f}-{bbox[1]:.1f}, N {bbox[2]:.1f}-{bbox[3]:.1f}")

STAR, LAT = {}, {}
for t in TILES:
    corr = json.load(open(f"data/derived/{t}/corrections.json")); assert corr["swath_frame"] == "star", (t, corr["swath_frame"])
    STAR[t] = {int(k): np.array(v, float) for k, v in corr["per_swath_star_frame_dxdydz_m"].items()}; LAT[t] = np.array(corr["cross_epoch_datum"]["horizontal_shift_m"], float)
    say(f"2008 product frame, tile {t}: star shifts ADDED to x,y,z per point_source_id {corr['per_swath_star_frame_dxdydz_m']}; lateral shift ADDED to x,y {LAT[t].tolist()}")

DIMS = ["return_number", "number_of_returns", "classification", "intensity", "scan_angle_rank", "point_source_id", "gps_time"]

def read_box(path, tile, chunk=5_000_000):
    """x, y, z and DIMS for returns inside bbox, chunked so the 1 GB 2021 files never sit in memory."""
    parts = []; n_all = 0; t0 = time.time()
    with laspy.open(path) as r:
        for ch in r.chunk_iterator(chunk):
            x = np.asarray(ch.x, float); y = np.asarray(ch.y, float); n_all += len(x)
            m = (x >= bbox[0]) & (x <= bbox[1]) & (y >= bbox[2]) & (y <= bbox[3])
            if m.any():
                d = {"x": x[m], "y": y[m], "z": np.asarray(ch.z, float)[m]}
                for k in DIMS: d[k] = np.asarray(ch[k])[m]
                d["tile"] = tile; parts.append(pd.DataFrame(d))
    df = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame(columns=["x", "y", "z"] + DIMS + ["tile"])
    say(f"{path}: {len(df):,} returns in bbox of {n_all:,} in file ({time.time()-t0:.0f} s)" + (f"; classification {dict(zip(*np.unique(df.classification, return_counts=True)))}; point_source_id {np.unique(df.point_source_id).tolist()}" if len(df) else ""))
    return df

def pair(df, xcol, ycol, zcol, suffix):
    """One row per (point, return) with distance <= RMAX, using the given coordinate columns."""
    rows = []
    x = df[xcol].to_numpy(); y = df[ycol].to_numpy(); z = df[zcol].to_numpy()
    for i in range(len(B)):
        d = np.hypot(x - E[i], y - N[i]); m = d <= RMAX
        if not m.any(): continue
        sub = df.loc[m].copy(); sub.insert(0, "point", PT[i]); sub.insert(1, "rtk_z_m", Z[i])
        sub[f"dist{suffix}_m"] = d[m]; sub[f"dz{suffix}_m"] = z[m] - Z[i]
        rows.append(sub)
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()

# ---------------------------------------------------------------- 2008
g1 = pd.concat([read_box(f"data/before/{t}.laz", t) for t in TILES], ignore_index=True)
sh = np.stack([STAR[t][p] for t, p in zip(g1.tile, g1.point_source_id)]); lat = np.stack([LAT[t] for t in g1.tile])
g1["x_star"] = g1.x + sh[:, 0] + lat[:, 0]; g1["y_star"] = g1.y + sh[:, 1] + lat[:, 1]; g1["z_star"] = g1.z + sh[:, 2]
say(f"2008 raw: {g1.shape[0]} returns in bbox, rn==nr {100*(g1.return_number==g1.number_of_returns).mean():.1f}%, single {100*(g1.number_of_returns==1).mean():.1f}%")
T1 = pair(g1, "x", "y", "z", "")
S1 = pair(g1, "x_star", "y_star", "z_star", "_star")[["point", "x", "y", "z", "gps_time", "dist_star_m", "dz_star_m"]]
T1 = T1.merge(S1, on=["point", "x", "y", "z", "gps_time"], how="outer")
for r in RADII: T1[f"radius_{int(r)}m"] = T1.dist_star_m <= r
T1 = T1.sort_values(["point", "dist_star_m"]).reset_index(drop=True)
T1.to_csv(f"{OUT}/{tag}_returns_2008.csv", index=False, float_format="%.4f")
say(f"wrote {OUT}/{tag}_returns_2008.csv: {len(T1)} point-return pairs (radius flags from the STAR-frame distance; dz_m is raw vendor z minus RTK, dz_star_m is the product's frame)")

# ---------------------------------------------------------------- 2021
g2 = pd.concat([read_box(f"data/after/valley/{t}_3dep_full.laz", t) for t in TILES], ignore_index=True)
# the two 2021 files were cut from the same 3DEP EPT cloud on tile bounds that overlap by 35 m (E 579995-580030), so a
# return in the overlap is in BOTH files; keep one copy (first tile in TILES). The 2008 MN DNR tiles partition their
# returns (checked: no identical return in both) and are not deduplicated.
_k = ["x", "y", "z", "gps_time", "point_source_id", "return_number"]; _n0 = len(g2); g2 = g2.drop_duplicates(subset=_k, keep="first").reset_index(drop=True)
say(f"2021: {_n0 - len(g2):,} duplicate returns in the tile overlap dropped, {len(g2):,} kept")
say(f"2021: rn==nr {100*(g2.return_number==g2.number_of_returns).mean():.1f}%, single {100*(g2.number_of_returns==1).mean():.1f}%")
T2 = pair(g2, "x", "y", "z", "")
for r in RADII: T2[f"radius_{int(r)}m"] = T2.dist_m <= r
T2 = T2.sort_values(["point", "dist_m"]).reset_index(drop=True)
T2.to_csv(f"{OUT}/{tag}_returns_2021.csv.gz", index=False, float_format="%.4f", compression={"method": "gzip", "mtime": 0})
say(f"wrote {OUT}/{tag}_returns_2021.csv.gz: {len(T2)} point-return pairs (gzipped; pandas reads it directly)")

# ---------------------------------------------------------------- per-point summary
rows = []
for i, p in enumerate(PT):
    row = {"point": p, "rtk_z_m": Z[i], "status": B.status[i], "gen1_minus_rtk_product_m": B.gen1_minus_rtk_m[i], "gen2_minus_rtk_product_m": B.gen2_minus_rtk_m[i]}
    for r in RADII:
        a = T1[(T1.point == p) & (T1[f"radius_{int(r)}m"])] if len(T1) else T1; b = T2[(T2.point == p) & (T2[f"radius_{int(r)}m"])]
        row[f"n2008_{int(r)}m"] = len(a); row[f"med_dz2008_star_{int(r)}m"] = a.dz_star_m.median() if len(a) else np.nan
        row[f"min_dz2008_star_{int(r)}m"] = a.dz_star_m.min() if len(a) else np.nan
        row[f"n2021_{int(r)}m"] = len(b); row[f"med_dz2021_{int(r)}m"] = b.dz_m.median() if len(b) else np.nan
        bg = b[b.classification == 2]; row[f"n2021_class2_{int(r)}m"] = len(bg); row[f"med_dz2021_class2_{int(r)}m"] = bg.dz_m.median() if len(bg) else np.nan
        row[f"min_dz2021_{int(r)}m"] = b.dz_m.min() if len(b) else np.nan
    rows.append(row)
S = pd.DataFrame(rows); S.to_csv(f"{OUT}/{tag}_returns_summary.csv", index=False, float_format="%.4f")
say(f"wrote {OUT}/{tag}_returns_summary.csv")
say(f"counts per point (min / median / max over the {len(S)} points):")
for r in RADII:
    for k in (f"n2008_{int(r)}m", f"n2021_{int(r)}m", f"n2021_class2_{int(r)}m"):
        say(f"   {k:18s} {S[k].min():5d} / {int(S[k].median()):5d} / {S[k].max():5d}")
fx = S[S.status == "FIX"]
say(f"median over the {len(fx)} FIX points of the per-point median offset from RTK (m):")
for r in RADII:
    say(f"   {int(r)} m: 2008 star {fx[f'med_dz2008_star_{int(r)}m'].median():+.3f}  2021 all {fx[f'med_dz2021_{int(r)}m'].median():+.3f}  "
        f"2021 class2 {fx[f'med_dz2021_class2_{int(r)}m'].median():+.3f}")
say(f"check: 2008 product gen1_minus_rtk (CSF ground, 5 m cell, star frame) median {fx.gen1_minus_rtk_product_m.median():+.3f}; "
    f"2021 product gen2_minus_rtk median {fx.gen2_minus_rtk_product_m.median():+.3f}")

# ---------------------------------------------------------------- representative point (MINE)
fix = S[S.status == "FIX"].copy()
def zs(v): return np.abs(v - v.median()) / v.std() if v.std() > 0 else 0 * v
fix["score"] = zs(fix.n2008_3m) + zs(fix.n2021_3m) + zs(fix.gen1_minus_rtk_product_m).fillna(0)
rep = int(sys.argv[1]) if len(sys.argv) > 1 else int(fix.sort_values("score").point.iloc[0])
rr = S[S.point == rep].iloc[0]
say(f"representative point: {rep} (MINE: FIX point nearest the {NAME} medians of product gen1_minus_rtk, n2008_3m and n2021_3m, equally weighted in sd units; override with argv[1]); "
    f"product gen1_minus_rtk {rr.gen1_minus_rtk_product_m:+.3f}, n2008_3m {rr.n2008_3m}, n2021_3m {rr.n2021_3m}, median dz 2008 star {rr.med_dz2008_star_3m:+.3f}, 2021 {rr.med_dz2021_3m:+.3f}")

import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
a = T1[(T1.point == rep) & T1.radius_3m].dz_star_m.to_numpy() if len(T1) else np.array([]); b = T2[(T2.point == rep) & T2.radius_3m].dz_m.to_numpy()
say(f"PDF data: 2008 n {len(a)}" + (f" range {a.min():+.2f}..{a.max():+.2f}" if len(a) else "") + f"; 2021 n {len(b)} range {b.min():+.2f}..{b.max():+.2f} (bins 0.02 m for 2021, 0.05 m for 2008, MINE)")
lo = min(a.min() if len(a) else b.min(), b.min()) - 0.05; hi = max(a.max() if len(a) else b.max(), b.max()) + 0.05
fig, ax = plt.subplots(figsize=(7, 4.2))
ax.hist(b, bins=np.arange(lo, hi, 0.02), density=True, histtype="stepfilled", color="#0072b2", alpha=0.35, lw=0)
ax.hist(b, bins=np.arange(lo, hi, 0.02), density=True, histtype="step", color="#0072b2", lw=1.5, label="2021 lidar")
if len(a):
    ax.hist(a, bins=np.arange(lo, hi, 0.05), density=True, histtype="step", color="#d55e00", lw=1.5, label="2008 lidar")
    ax.plot(a, np.full_like(a, -0.02 * ax.get_ylim()[1]), "|", color="#d55e00", ms=10, mew=1.2, clip_on=False)
ax.axvline(0, color="0.4", lw=0.8, ls="--")
ax.set_xlabel("lidar return minus RTK ground (m)"); ax.set_ylabel("probability density (1/m)")
ax.legend(frameon=False); ax.spines[["top", "right"]].set_visible(False)
fig.tight_layout()
png = f"figures/rtk_2026-10-25_{TAG}_returns_pdf_p{rep}.png"
fig.savefig(png, dpi=int(os.environ.get("FIGDPI", 150))); say(f"wrote {png}")
LOG.close()

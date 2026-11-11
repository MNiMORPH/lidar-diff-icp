#!/usr/bin/env python3
"""Compile the class RTK-GNSS survey of 25 October 2026 (Whitewater WMA, lot 204 and the bluff
north of Elba) into one table, and attach the lidar surfaces at every point.

Andy, 2026-11-11: "Can you build an accurate data set of ground elevations from this?";
"Mine ran CORS. All others used the base station."; "T5a occupied the base station."

SOURCES (data/field/2026-10-25/raw/, copied from the Canvas submissions; student-facing
names kept so the provenance is traceable):
  rover_RS3_8243D2AFB94BDAE3_hillslope_T3a_T3b_T3c.csv   Emlid Flow CSV, 24 pts
  rover_hillslope_T4a_T4b_from_docx_table.csv                   parsed from T4a's
        F07 .docx table (E, N, ellipsoid h, base, time, RMS); the raw export was not uploaded
  rover_RSplus_8243CA7256E53180_floodplain_T1a_T1b.xlsx     Emlid Flow export, 63 pts
  rover_RSplus_824350D44A583F7B_floodplain_T2a.csv              Emlid Flow CSV from the
        shared Google Sheet, 100 pts
  (T3a's other CSV is a 9 October practice session in St Paul and is excluded.)
  benchmark_DR9925_occupation_T5a_flow_export.csv   T5a's unrounded Flow export
        (received 2026-11-11): TWO 60-second FIX occupations of monument DR9925 ("8508 Y") under
        base B2 at 15:03 and 15:37 CDT, Reach RS2, ellipsoid h 189.202 and 189.210 m. The PDF
        transcription (benchmark_DR9925_occupation_T5a.csv, positions rounded) is superseded.
        "Purpose: To reference other teams' GNSS data to a known location."
  rover_CORS_Wickert_Whitewater_floodplain_LiDAR_veg.csv   Andy's Reach RS2+ on MnDOT CORS
        (mount point RTCM_34_NAD83(2011), a network solution; CS "NAD83(2011) / UTM 15N +
        NAVD88(GEOID18)"): 17 FIX points 15:50-15:56 CDT in the SW corner of the field east
        of the road. These are in NAD83(2011) directly and need no tie; where a student B2
        point lies within 3 m of one of them the two are compared, which checks the
        benchmark tie against an independent frame.

THE TIE (Andy, 2026-11-11: "Let's trust the reference point"). DR9925's published NAD83(2011)
position is 44 07 18.13863 N 092 00 16.79318 W, ellipsoid height 188.669 m (NGS datasheet,
network accuracy 0.61 cm horizontal, 1.49 cm height). T5a's two occupations (unrounded
Flow export, raw/benchmark_DR9925_occupation_T5a_flow_export.csv) read 189.202 and
189.210 m, mean 189.206, so everything measured under B2 is HIGH by +0.537 m, SUBTRACTED
from every B2 point: column `h_ell_tied_m`. Horizontally the occupations sit at E 579636.167
/ .176, N 4885873.225 / .203 against the published E 579636.942 N 4885872.102, so every B2
point is SHIFTED by their difference: columns `easting_tied`, `northing_tied`, and the lidar
is read there. The CORS 11 = T1a 63 pair, the tie applied before this export arrived, is
printed as the check (horizontal within 0.03 m; vertical 64 mm short, see the antenna note).

BASES. The student rovers show two base positions. "B1" (lon -92.000382, lat 44.116855,
h 187.878, ~13:20-13:30 CDT) was NOT the class base: a student's receiver was accidentally
configured as a base and confused the radio link (Andy, 2026-11-11: "Remove it and all
points that count on it"). Every point whose corrections came from B1 is DROPPED from the
table (listed in the log). B2 (lon -91.998351, lat 44.117216, h 191.032, ~14:40 on) is the
class base; its averaged position carries an unknown offset, removed by the ties below.

HEIGHTS. Ellipsoid heights are as exported (NAD83(2011) ellipsoid, the base's frame).
NAVD 88 is computed with GEOID18 through pyproj (EPSG:6319 -> 6318+5703) when the grid is
available; the geoid source is printed. The Elevation column in the Flow exports (where
present) is the app's own NAVD88 value and is kept as `elev_app_m` for comparison.

LIDAR at each point: gen2 z_after (2021 ground, NAVD88) and the 2021-2008 difference from
the Phase B tiles (data/derived/<tile>/, star frame, no correction surface, each tile's own
stable-ground zero), plus the Clubb floodplain flag, slope and the gen2 lowveg share.

Writes data/field/2026-10-25/rtk_points_all.csv and prints the per-base / per-status counts
and the bare-ground check. Nothing is filtered out: `solution_status` is a column.
"""
import glob
import json
import os
import re

import numpy as np
import pandas as pd
import pyproj

RAW = "data/field/2026-10-25/raw"; OUT = "data/field/2026-10-25/rtk_points_all.csv"
# PROJECTION ONLY, no datum shift: the lidar grids are NAD83(2011) UTM 15N (3DEP; the pipeline's "EPSG:26915"
# is a label on the same meters). EPSG:6318 -> EPSG:26915 made pyproj insert a NAD83(2011) -> NAD83 datum
# shift of +0.012 m E, -0.242 m N here, which put every surveyed point 0.24 m south of the lidar (fixed 2026-11-11).
to_utm = pyproj.Transformer.from_crs("EPSG:6318", "EPSG:6344", always_xy=True)   # NAD83(2011) -> NAD83(2011) / UTM 15N
to_ll = pyproj.Transformer.from_crs("EPSG:6344", "EPSG:6318", always_xy=True)
B1 = dict(label="B1", lon=-92.00038215, lat=44.11685493, h=187.878, where="lot 204 (gravel lot)")
B2 = dict(label="B2", lon=-91.99835123, lat=44.11721583, h=191.032, where="field east of the WMA road")

def base_label(blon):
    if pd.isna(blon): return "none recorded"
    return "B1" if abs(blon - B1["lon"]) < 1e-4 else ("B2" if abs(blon - B2["lon"]) < 1e-4 else f"other({blon:.5f})")

rows = []
# 1. T3a / T3b / T3c, Reach RS3, hillslope west of the road
d = pd.read_csv(f"{RAW}/rover_RS3_8243D2AFB94BDAE3_hillslope_T3a_T3b_T3c.csv")
for _, r in d.iterrows():
    rows.append(dict(group="hillslope_RS3", people="T3a, T3b, T3c", device=r["Device type"], serial=r["Device serial number"], point=int(r["Name"]),
                     lon=r["Longitude"], lat=r["Latitude"], h_ell_m=r["Ellipsoidal height"], elev_app_m=r["Elevation"], antenna_m=r["Antenna height"], status=r["Solution status"],
                     rms_e=r["Easting RMS"], rms_n=r["Northing RMS"], rms_v=r["Elevation RMS"], t_start=r["Averaging start"], base=base_label(r["Base longitude"]), baseline_m=r["Baseline"], tilt_deg=r.get("Tilt angle", np.nan), source="rover_RS3_8243D2AFB94BDAE3_hillslope_T3a_T3b_T3c.csv"))
# 2. T4a / T4b, from the docx table (UTM, ellipsoid h, base label, HH:MM CDT)
d = pd.read_csv(f"{RAW}/rover_hillslope_T4a_T4b_from_docx_table.csv")
for _, r in d.iterrows():
    lon, lat = to_ll.transform(r.Easting, r.Northing)
    rows.append(dict(group="hillslope_T4a", people="T4a, T4b", device="Reach (model not in the table)", serial="", point=int(r.Pt), lon=lon, lat=lat, h_ell_m=r.Ellipsoidal_height, elev_app_m=np.nan, antenna_m=np.nan,
                     status="FIX (table lists RMS only; status not exported)", rms_e=r.Lateral_RMS, rms_n=np.nan, rms_v=r.Vertical_RMS, t_start=f"2026-10-25 {r.Time_CDT} CDT", base={"B1": "B1", "B2": "B2"}.get(str(r.Base).strip(), "none recorded"), baseline_m=np.nan, tilt_deg=np.nan, source="rover_hillslope_T4a_T4b_from_docx_table.csv"))
# 3. T1a / T1b, Reach RS+, floodplain field
d = pd.read_excel(f"{RAW}/rover_RSplus_8243CA7256E53180_floodplain_T1a_T1b.xlsx").dropna(subset=["Point"])
for _, r in d.iterrows():
    rows.append(dict(group="floodplain_RSplus_A", people="T1a, T1b", device="Reach RS+", serial="8243CA7256E53180", point=int(r["Point"]), lon=r["Longitude"], lat=r["Latitude"], h_ell_m=r["Ellipsoidal height (m)"], elev_app_m=np.nan, antenna_m=r["Antenna height (m)"],
                     status=r["Solution status"], rms_e=r["Easting RMS"], rms_n=r["Northing RMS"], rms_v=r["Elevation RMS"], t_start=r["Averaging start"], base=base_label(r["Base longitude"]), baseline_m=r["Baseline"], tilt_deg=np.nan, source="rover_RSplus_8243CA7256E53180_floodplain_T1a_T1b.xlsx"))
# 4. T2a (+ T5a?), Reach RS+, floodplain field
d = pd.read_csv(f"{RAW}/rover_RSplus_824350D44A583F7B_floodplain_T2a.csv")
for _, r in d.iterrows():
    rows.append(dict(group="floodplain_RSplus_B", people="T2a (+ group)", device=r["Device type"], serial=r["Device serial number"], point=int(r["Name"]), lon=r["Longitude"], lat=r["Latitude"], h_ell_m=r["Ellipsoidal height"], elev_app_m=r["Elevation"], antenna_m=r["Antenna height"],
                     status=r["Solution status"], rms_e=r["Easting RMS"], rms_n=r["Northing RMS"], rms_v=r["Elevation RMS"], t_start=r["Averaging start"], base=base_label(r["Base longitude"]), baseline_m=r["Baseline"], tilt_deg=np.nan, source="rover_RSplus_824350D44A583F7B_floodplain_T2a.csv"))
# 5. Andy's CORS rover
d = pd.read_csv(f"{RAW}/rover_CORS_Wickert_Whitewater_floodplain_LiDAR_veg.csv")
for _, r in d.iterrows():
    rows.append(dict(group="floodplain_CORS_Wickert", people="Wickert", device=r["Device type"], serial=r["Device serial number"], point=int(r["Name"]), lon=r["Longitude"], lat=r["Latitude"], h_ell_m=r["Ellipsoidal height"], elev_app_m=r["Elevation"], antenna_m=r["Antenna height"],
                     status=r["Solution status"], rms_e=r["Easting RMS"], rms_n=r["Northing RMS"], rms_v=r["Elevation RMS"], t_start=r["Averaging start"], base="CORS", baseline_m=r["Baseline"], tilt_deg=np.nan, source="rover_CORS_Wickert_Whitewater_floodplain_LiDAR_veg.csv"))
P = pd.DataFrame(rows)
# the benchmark occupation itself, as a fifth record set (so it is in the table and the check is visible)
# T5a's unrounded Flow export (2026-11-11) replaces the PDF transcription: the MEASURED position, under B2
# like every other student point. Andy (2026-11-11): "The reference point was made on the metal disk, properly."
bm = pd.read_csv(f"{RAW}/benchmark_DR9925_occupation_T5a_flow_export.csv")
for _, r in bm.iterrows():
    rows.append(dict(group="benchmark_DR9925", people="T5a, T5b, T5c", device=r["Device type"], serial=r["Device serial number"], point=int(r["Name"]), lon=r["Longitude"], lat=r["Latitude"], h_ell_m=r["Ellipsoidal height"], elev_app_m=r["Elevation"], antenna_m=r["Antenna height"],
                     status=r["Solution status"], rms_e=r["Easting RMS"], rms_n=r["Northing RMS"], rms_v=r["Elevation RMS"], t_start=r["Averaging start"], base="B2", baseline_m=r["Baseline"], tilt_deg=np.nan, source="benchmark_DR9925_occupation_T5a_flow_export.csv (measured position under B2; tied position = the published monument)"))
P = pd.DataFrame(rows)
# a few exports lost the base columns on individual points (the link dropped and came back); the base is
# assigned by TIME: B1 ran ~13:20-13:30 CDT, B2 from ~14:40 (the earliest B2 point is 14:40:14). Stated per point.
P["t_hhmm"] = P.t_start.str.extract(r"(\d\d:\d\d)")[0]
lost = P.base == "none recorded"; P.loc[lost & (P.t_hhmm >= "14:30"), "base"] = "B2"; P.loc[lost & (P.t_hhmm < "14:30"), "base"] = "B1"
print(f"base assigned by time on {int(lost.sum())} points whose export had no base columns: {P.loc[lost, ['group', 'point', 't_hhmm', 'base']].to_dict('records')}")
# DROP B1: a student's receiver acting as a base; everything that used its corrections goes
b1 = P.base == "B1"; print(f"DROPPED {int(b1.sum())} points corrected from B1 (a student's receiver, not the class base): {P.loc[b1, ['group', 'point', 't_hhmm', 'status']].to_dict('records')}")
P = P[~b1].reset_index(drop=True)
# DR9925, NGS datasheet retrieved 2026-11-11 (raw/DR9925_NGS_datasheet_2026-11-11.txt): NAD 83(2011) 44 07 18.13863 N
# 092 00 16.79318 W, ellipsoid 188.669 m (network accuracy 0.61 cm horizontal, 1.49 cm ellipsoid height, 95 %),
# UTM 15 N 4,885,872.102 E 579,636.942 (computed from the position, no datum shift; to_utm reproduces it to 1 mm).
DR9925_H_ELL = 188.669; DR9925_LON, DR9925_LAT = -(92 + 0 / 60 + 16.79318 / 3600), 44 + 7 / 60 + 18.13863 / 3600
B2_OFFSET_BM = float(bm["Ellipsoidal height"].mean() - DR9925_H_ELL)
print(f"B2 vertical, benchmark: DR9925 measured {bm['Ellipsoidal height'].tolist()} m vs published {DR9925_H_ELL} m -> B2 reads {B2_OFFSET_BM:+.3f} m  [ADOPTED, Andy 2026-11-11: trust the reference point]")
P["easting"], P["northing"] = to_utm.transform(P.lon.values, P.lat.values)
P["easting"] = P.easting.round(3); P["northing"] = P.northing.round(3)

# B2 HORIZONTAL tie: the published monument minus the measured occupation (mean of the two FIX points),
# applied to every B2 point as easting_tied / northing_tied. ADOPTED, Andy 2026-11-11: "Let's trust the
# reference point and not on my comparison with CORS."
_pe, _pn = to_utm.transform(DR9925_LON, DR9925_LAT); _bmrows = P[P.group == "benchmark_DR9925"]
B2_DX, B2_DY = float(_pe - _bmrows.easting.mean()), float(_pn - _bmrows.northing.mean())
print(f"B2 HORIZONTAL TIE, benchmark: published DR9925 E {_pe:.3f} N {_pn:.3f} (datasheet UTM 579636.942 / 4885872.102) minus measured E {_bmrows.easting.mean():.3f} N {_bmrows.northing.mean():.3f} -> B2 shifted by dE {B2_DX:+.3f} m, dN {B2_DY:+.3f} m  [ADOPTED]; CORS unshifted")
# THE CHECK: T1a point 63 (15:53:33) and Andy's CORS point 11 (15:53:52), the same spot one after the other
# (Andy: "accurate to within a centimeter"). CORS minus T1a is an independent tie in both axes; it was the
# applied tie until the unrounded benchmark export arrived. Its vertical value sits 64 mm below the benchmark's;
# T1a's Reach RS+ was exported with antenna height 2.134 m (the RS2 value; Emlid's RS+ offset is 65 mm, and
# T2a's RS+ carries 2.065), which would put his points 69 mm low and the pair tie 69 mm short. NOT applied.
_d63 = P[(P.group == "floodplain_RSplus_A") & (P.point == 63)].iloc[0]; _c11 = P[(P.group == "floodplain_CORS_Wickert") & (P.point == 11)].iloc[0]
B2_OFFSET_C11 = float(_d63.h_ell_m - _c11.h_ell_m); B2_OFFSET = B2_OFFSET_BM
print(f"CHECK, CORS 11 pair: horizontal CORS minus T1a dE {float(_c11.easting - _d63.easting):+.3f} dN {float(_c11.northing - _d63.northing):+.3f} (benchmark minus pair dE {B2_DX - float(_c11.easting - _d63.easting):+.3f} dN {B2_DY - float(_c11.northing - _d63.northing):+.3f}); vertical T1a 63 ({_d63.h_ell_m:.3f}) minus CORS 11 ({_c11.h_ell_m:.3f}) = {B2_OFFSET_C11:+.3f} m, benchmark minus pair {1000 * (B2_OFFSET_BM - B2_OFFSET_C11):+.0f} mm; T1a RS+ antenna 2.134 vs the RS+ offset 2.065 = 69 mm")
P["easting_tied"] = np.where(P.base == "B2", P.easting + B2_DX, P.easting).round(3); P["northing_tied"] = np.where(P.base == "B2", P.northing + B2_DY, P.northing).round(3)
# tied ellipsoid heights: B2 by the benchmark; CORS direct
P["h_ell_tied_m"] = np.where(P.base == "B2", P.h_ell_m - B2_OFFSET, np.where(P.base == "CORS", P.h_ell_m, np.nan)); P["tie_source"] = np.where(P.base == "B2", f"DR9925 benchmark under B2 ({B2_OFFSET:+.3f} m removed, dE {B2_DX:+.3f} dN {B2_DY:+.3f} m shifted; CORS 11 = T1a 63 pair gives {B2_OFFSET_C11:+.3f} m)", np.where(P.base == "CORS", "MnDOT CORS network RTK, NAD83(2011) direct", ""))
# NAVD 88 via GEOID18 (PROJ network fetches the grid)
geoid_source = "none"
try:
    pyproj.network.set_network_enabled(True)
    tg = pyproj.Transformer.from_crs("EPSG:6319", "EPSG:6318+5703", always_xy=True)
    x, y, H = tg.transform(-92.00466477, 44.12170518, 188.669)
    if abs(H - 218.75) < 0.05:
        _, _, P["navd88_m"] = tg.transform(P.lon.values, P.lat.values, P.h_ell_m.values); geoid_source = f"pyproj EPSG:6319->6318+5703 (GEOID18); check at DR9925: {H:.3f} vs published 218.75"
    else: geoid_source = f"pyproj transform returned {H:.3f} at DR9925 (published 218.75): NOT GEOID18, not used"
except Exception as e: geoid_source = f"pyproj geoid transform unavailable: {e}"
if "navd88_m" not in P: P["navd88_m"] = P.h_ell_m + 30.081; geoid_source += "; FALLBACK: constant N = -30.081 m from DR9925 (ellipsoid 188.669, NAVD88 218.75), i.e. the geoid slope across the 1 km site is neglected"
print("geoid:", geoid_source)

# lidar at each point, from the Phase B tiles
tiles = {}
for c in glob.glob("data/derived/4*/corrections.json"):
    t = os.path.basename(os.path.dirname(c)); j = json.load(open(c)); tiles[t] = (j["bounds"], float(j.get("res") or j.get("res_m")))
_zero = {}
def sample(name, x, y, kind="float"):
    """z_after and dod are read BILINEARLY from the 5 m grid (cell centers at +2.5 m): on a 30 deg
    slope the nearest-cell value is off by up to 3.5 m x tan(30) = 2 m, which the first run showed
    as an NMAD over 1 m on the hillslope groups. Masks and slope are nearest-cell."""
    for t, (b, res) in tiles.items():
        if b[0] <= x < b[2] and b[1] <= y < b[3]:
            p = f"data/derived/{t}/{name}"
            if not os.path.exists(p): return np.nan, t
            a = np.load(p, mmap_mode="r")
            if kind == "bilinear":
                fy, fx = (y - b[1]) / res - 0.5, (x - b[0]) / res - 0.5; i0, j0 = int(np.floor(fy)), int(np.floor(fx)); wy, wx = fy - i0, fx - j0
                i0 = min(max(i0, 0), a.shape[0] - 2); j0 = min(max(j0, 0), a.shape[1] - 2); q = np.asarray(a[i0:i0 + 2, j0:j0 + 2], float)
                if not np.isfinite(q).all(): v = float(a[int((y - b[1]) / res), int((x - b[0]) / res)])
                else: v = float((1 - wy) * ((1 - wx) * q[0, 0] + wx * q[0, 1]) + wy * ((1 - wx) * q[1, 0] + wx * q[1, 1]))
            else: v = float(a[int((y - b[1]) / res), int((x - b[0]) / res)])
            if name == "dod.npy":
                if t not in _zero:
                    st = np.load(f"data/derived/{t}/stable.npy", mmap_mode="r"); d = np.asarray(np.load(p, mmap_mode="r")); _zero[t] = float(np.nanmedian(d[np.asarray(st).astype(bool) & np.isfinite(d)]))
                v = v - _zero[t]
            return v, t
    return np.nan, ""
for col, lay, kind in (("gen2_z_m", "z_after.npy", "bilinear"), ("gen2_zref_m", "z_after_differenced.npy", "bilinear"), ("dod_2021_2008_m", "dod.npy", "bilinear"), ("slope_deg", "slope.npy", "nearest"), ("lowveg_gen2", "lowveg_grid.npy", "nearest"), ("clubb_floodplain", "floodplain_mask_clubb.npy", "nearest"), ("penetration_failure_2008", "penetration_failure.npy", "nearest")):
    vals = [sample(lay, x, y, kind) for x, y in zip(P.easting_tied, P.northing_tied)]; P[col] = [v for v, _ in vals]   # lidar read at the TIED position
P["tile"] = [sample("z_after.npy", x, y)[1] for x, y in zip(P.easting_tied, P.northing_tied)]
# gen1 = the 2021 surface the difference was built on (z_after_differenced = slope-normal median at the
# cell center, Zref) minus the difference. z_after.npy is the HORIZONTAL per-cell median (Z21), which the
# point-cloud check reproduces at the survey points (-8 mm median) and is kept as gen2_z_m; on the bluff
# the two differ by up to 0.13 m at single points, so mixing them put that error into gen1 (fixed 2026-11-11).
P["gen1_z_m"] = P.gen2_zref_m - P.dod_2021_2008_m
P["rtk_minus_gen2_m"] = P.navd88_m - P.gen2_z_m; P["rtk_minus_gen1_m"] = P.navd88_m - P.gen1_z_m
P["navd88_tied_m"] = P.navd88_m - (P.h_ell_m - P.h_ell_tied_m)
# CORS vs B2 cross-check: nearest B2 FIX point within 3 m of each CORS point, raw ellipsoid-height difference (B2 minus CORS); the benchmark says +B2_OFFSET (the B2 points here are T1a's RS+, see the antenna note above)
cors = P[P.base == "CORS"]; b2f = P[(P.base == "B2") & P.status.str.startswith("FIX") & (P.group != "benchmark_DR9925")]
pairs = []
for _, c in cors.iterrows():
    d_ = np.hypot(b2f.easting - c.easting, b2f.northing - c.northing); k = d_.idxmin()
    if d_[k] <= 3.0: pairs.append(dict(cors_pt=c.point, b2_group=b2f.loc[k, "group"], b2_pt=b2f.loc[k, "point"], sep_m=round(float(d_[k]), 2), slope_deg=round(float(c.slope_deg), 1) if np.isfinite(c.slope_deg) else np.nan, b2_minus_cors_m=round(float(b2f.loc[k, "h_ell_m"] - c.h_ell_m), 3)))
PAIRS = pd.DataFrame(pairs)
print(f"\nCROSS-CHECK, B2 against CORS on the same ground (pairs within 3 m; the 3 m is mine): {len(PAIRS)} pairs")
if len(PAIRS): print(PAIRS.to_string(index=False)); print(f"  median B2 minus CORS {PAIRS.b2_minus_cors_m.median():+.3f} m, NMAD {1.4826 * np.median(np.abs(PAIRS.b2_minus_cors_m - PAIRS.b2_minus_cors_m.median())):.3f} m; benchmark tie said {B2_OFFSET:+.3f} m")
P["tied_minus_gen2_m"] = P.navd88_tied_m - P.gen2_z_m; P["tied_minus_gen1_m"] = P.navd88_tied_m - P.gen1_z_m
P = P[["group", "people", "device", "serial", "point", "t_start", "t_hhmm", "base", "baseline_m", "status", "lon", "lat", "easting", "northing", "easting_tied", "northing_tied", "h_ell_m", "h_ell_tied_m", "tie_source", "navd88_m", "navd88_tied_m", "elev_app_m", "antenna_m", "rms_e", "rms_n", "rms_v", "tilt_deg", "tile", "gen2_z_m", "gen2_zref_m", "gen1_z_m", "dod_2021_2008_m", "rtk_minus_gen2_m", "rtk_minus_gen1_m", "tied_minus_gen2_m", "tied_minus_gen1_m", "slope_deg", "lowveg_gen2", "clubb_floodplain", "penetration_failure_2008", "source"]]
os.makedirs(os.path.dirname(OUT), exist_ok=True); P.to_csv(OUT, index=False)
pd.set_option("display.width", 250); pd.set_option("display.max_columns", 40)
print(f"\n{len(P)} points written to {OUT}")
print("\nby group / base / status:"); print(P.groupby(["group", "base", P.status.str[:5]]).size().to_string())
fix = P[P.status.str.startswith("FIX")]
print("\nFIX points, RTK NAVD88 minus gen2 (2021) ground, mm, by group and base [median, NMAD, n]:")
for (g, b), s in fix.groupby(["group", "base"]):
    v = s.rtk_minus_gen2_m.dropna() * 1000; print(f"  {g:<22} {b:<14} {np.median(v):+8.0f}  {1.4826 * np.median(np.abs(v - np.median(v))):6.0f}  {len(v):4d}")
print("\nAFTER THE TIE (FIX points): tied NAVD88 minus gen2, mm, by group and base [median, NMAD, n]:")
for (g, b), s_ in fix.groupby(["group", "base"]):
    v = P.loc[s_.index, "tied_minus_gen2_m"].dropna() * 1000
    if len(v): print(f"  {g:<22} {b:<14} {np.median(v):+8.0f}  {1.4826 * np.median(np.abs(v - np.median(v))):6.0f}  {len(v):4d}")
print("\nAFTER THE TIE (FIX points): tied NAVD88 minus gen1 (2008) ground, mm, by group and base [median, NMAD, n]:")
for (g, b), s_ in fix.groupby(["group", "base"]):
    v = P.loc[s_.index, "tied_minus_gen1_m"].dropna() * 1000
    if len(v): print(f"  {g:<22} {b:<14} {np.median(v):+8.0f}  {1.4826 * np.median(np.abs(v - np.median(v))):6.0f}  {len(v):4d}")
print("\nBARE-GROUND CHECK (slope <= 3 deg, lowveg <= 0.05, FIX): median tied RTK minus gen2 per base (note: on the floodplain the gen2 lowveg filter does not separate cut sward from bare ground):")
bare = fix[(fix.slope_deg <= 3) & (fix.lowveg_gen2 <= 0.05)]
for b, s in bare.groupby("base"):
    v = s.tied_minus_gen2_m.dropna() * 1000; print(f"  {b}: {np.median(v):+.0f} mm (NMAD {1.4826 * np.median(np.abs(v - np.median(v))):.0f}, n {len(v)}; groups {sorted(s.group.unique())})")
print("  (the 0.05 lowveg and 3 deg cut-offs are MINE, for the check only)")

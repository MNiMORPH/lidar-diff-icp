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
  benchmark_DR9925_occupation_T5a.csv   transcribed from the PDF embedded in T5a
        T5a's .docx: TWO 60-second FIX occupations of monument DR9925 ("8508 Y") under
        base B2 at 15:03 and 15:37 CDT, Reach RS2, ellipsoid h 189.202 and 189.210 m.
        "Purpose: To reference other teams' GNSS data to a known location."
  rover_CORS_Wickert_Whitewater_floodplain_LiDAR_veg.csv   Andy's Reach RS2+ on MnDOT CORS
        (mount point RTCM_34_NAD83(2011), a network solution; CS "NAD83(2011) / UTM 15N +
        NAVD88(GEOID18)"): 17 FIX points 15:50-15:56 CDT in the SW corner of the field east
        of the road. These are in NAD83(2011) directly and need no tie; where a student B2
        point lies within 3 m of one of them the two are compared, which checks the
        benchmark tie against an independent frame.

THE TIE. DR9925's published NAD83(2011) ellipsoid height is 188.669 m (NGS datasheet,
network accuracy 1.49 cm). T5a's mean occupation reads 189.206 m, so everything
measured under B2 is HIGH by +0.537 m (the two occupations agree to 8 mm, 34 minutes
apart), and that constant is SUBTRACTED from every B2 point here: column `h_ell_tied_m`.
Horizontally the printed occupation (E 579636.2, N 4885873) sits within the PDF's own
rounding (0.1 m, 1 m) of the published E 579636.94 N 4885872.10, so no horizontal shift is
applied and B2's horizontal error is bounded only at the metre level by this record. B1 has
no benchmark occupation: its points are tied PROVISIONALLY by the bare-ground check below
(their offset from the 2021 lidar on the gravel lot), stated as such in `tie_source`.

BASES. Every student rover worked from the class base station, which was set TWICE: B1 at
lot 204 (lon -92.000382, lat 44.116855, ellipsoid h 187.878) from ~13:20 to ~13:30 CDT, and
B2 in the field east of the road (lon -91.998351, lat 44.117216, h 191.032) from ~14:40 on.
Both base positions came from the base's own averaging (no known point was entered), so
each carries an unknown absolute error of order a metre and the two differ by an unknown
constant. Rover coordinates are therefore RELATIVE to their base. This script keeps the
base label on every point and reports, for each base, the vertical offset of its bare-ground
points from the 2021 lidar (gen2, which reads bare ground to -3.5 mm against surveyed
control) as a provisional datum check; the real tie is the benchmark on Andy's rover.

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
to_utm = pyproj.Transformer.from_crs("EPSG:6318", "EPSG:26915", always_xy=True)   # NAD83(2011) -> UTM 15N (NAD83)
to_ll = pyproj.Transformer.from_crs("EPSG:26915", "EPSG:6318", always_xy=True)
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
bm = pd.read_csv(f"{RAW}/benchmark_DR9925_occupation_T5a.csv", comment="#")
for _, r in bm.iterrows():
    rows.append(dict(group="benchmark_DR9925", people="T5a, T5b, T5c", device="Reach RS2", serial="824377FA5B116E0F", point=int(r.point), lon=-92.00466477, lat=44.12170518, h_ell_m=r.h_ell_m, elev_app_m=r.elev_app_m, antenna_m=2.134,
                     status=r.status, rms_e=r.rms_e, rms_n=r.rms_n, rms_v=r.rms_v, t_start=f"{r.t_start_cdt} UTC-05:00", base="B2", baseline_m=r.baseline_m, tilt_deg=np.nan, source="benchmark_DR9925_occupation_T5a.csv (position = the published monument position; the PDF prints the measured one rounded)"))
P = pd.DataFrame(rows)
# a few exports lost the base columns on individual points (the link dropped and came back); the base is
# assigned by TIME: B1 ran ~13:20-13:30 CDT, B2 from ~14:40 (the earliest B2 point is 14:40:14). Stated per point.
P["t_hhmm"] = P.t_start.str.extract(r"(\d\d:\d\d)")[0]
lost = P.base == "none recorded"; P.loc[lost & (P.t_hhmm >= "14:30"), "base"] = "B2"; P.loc[lost & (P.t_hhmm < "14:30"), "base"] = "B1"
print(f"base assigned by time on {int(lost.sum())} points whose export had no base columns: {P.loc[lost, ['group', 'point', 't_hhmm', 'base']].to_dict('records')}")
DR9925_H_ELL = 188.669; B2_OFFSET_BM = float(bm.h_ell_m.mean() - DR9925_H_ELL)
print(f"B2 vertical, benchmark: DR9925 measured {bm.h_ell_m.tolist()} m vs published {DR9925_H_ELL} m -> B2 reads {B2_OFFSET_BM:+.3f} m")
P["easting"], P["northing"] = to_utm.transform(P.lon.values, P.lat.values)
P["easting"] = P.easting.round(3); P["northing"] = P.northing.round(3)

# B2 HORIZONTAL tie (provisional, Andy 2026-11-11: "I am guessing our tie point is 11"): T1a 63
# (15:53:33) and Andy's CORS point 11 (15:53:52) were the same spot, one after the other. The shift
# CORS minus T1a is applied to every B2 point: easting_tied / northing_tied. Test
# (test_b2_horizontal_tie.py): hillslope NMAD vs gen2 540 -> 139 mm, and a free fit over +/-2.5 m
# lands within 0.1 m of the same shift. T5a's unrounded benchmark export replaces it when it comes.
_d63 = P[(P.group == "floodplain_RSplus_A") & (P.point == 63)].iloc[0]; _c11 = P[(P.group == "floodplain_CORS_Wickert") & (P.point == 11)].iloc[0]
B2_DX, B2_DY = float(_c11.easting - _d63.easting), float(_c11.northing - _d63.northing)
# VERTICAL tie from the same pair (Andy, 2026-11-11: "And vertical tie too"): B2 reads T1a 63 minus CORS 11.
# ADOPTED for B2 at Andy's word; the benchmark value is printed beside it and the two disagree by their difference.
B2_OFFSET_C11 = float(_d63.h_ell_m - _c11.h_ell_m); B2_OFFSET = B2_OFFSET_C11
print(f"B2 vertical, CORS 11 pair: T1a 63 ({_d63.h_ell_m:.3f}) minus CORS 11 ({_c11.h_ell_m:.3f}) -> B2 reads {B2_OFFSET_C11:+.3f} m  [ADOPTED]; benchmark said {B2_OFFSET_BM:+.3f}: the two ties differ by {1000 * (B2_OFFSET_BM - B2_OFFSET_C11):+.0f} mm")
P["easting_tied"] = np.where(P.base == "B2", P.easting + B2_DX, P.easting).round(3); P["northing_tied"] = np.where(P.base == "B2", P.northing + B2_DY, P.northing).round(3)
print(f"B2 HORIZONTAL TIE (provisional): T1a 63 = CORS 11 -> B2 shifted by dE {B2_DX:+.3f} m, dN {B2_DY:+.3f} m; B1 and CORS unshifted")
# tied ellipsoid heights: B2 by the benchmark; B1 provisional (bare-ground check, filled in below after the lidar sampling)
P["h_ell_tied_m"] = np.where(P.base == "B2", P.h_ell_m - B2_OFFSET, np.where(P.base == "CORS", P.h_ell_m, np.nan)); P["tie_source"] = np.where(P.base == "B2", f"CORS 11 = T1a 63 pair under B2 ({B2_OFFSET:+.3f} m removed; benchmark gives {B2_OFFSET_BM:+.3f})", np.where(P.base == "CORS", "MnDOT CORS network RTK, NAD83(2011) direct", ""))
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
    """z_after and dod are read BILINEARLY from the 5 m grid (cell centres at +2.5 m): on a 30 deg
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
for col, lay, kind in (("gen2_z_m", "z_after.npy", "bilinear"), ("dod_2021_2008_m", "dod.npy", "bilinear"), ("slope_deg", "slope.npy", "nearest"), ("lowveg_gen2", "lowveg_grid.npy", "nearest"), ("clubb_floodplain", "floodplain_mask_clubb.npy", "nearest"), ("penetration_failure_2008", "penetration_failure.npy", "nearest")):
    vals = [sample(lay, x, y, kind) for x, y in zip(P.easting_tied, P.northing_tied)]; P[col] = [v for v, _ in vals]   # lidar read at the TIED position
P["tile"] = [sample("z_after.npy", x, y)[1] for x, y in zip(P.easting_tied, P.northing_tied)]
P["gen1_z_m"] = P.gen2_z_m - P.dod_2021_2008_m
P["rtk_minus_gen2_m"] = P.navd88_m - P.gen2_z_m; P["rtk_minus_gen1_m"] = P.navd88_m - P.gen1_z_m
# PROVISIONAL B1 tie: median offset of B1's bare, flat, FIX points from gen2 (gen2 reads bare ground to -3.5 mm against control). MINE, replaced when Andy's CORS rover gives a B1-era common point.
fixm = P.status.str.startswith("FIX"); bare1 = fixm & (P.base == "B1") & (P.slope_deg <= 3) & (P.lowveg_gen2 <= 0.05) & np.isfinite(P.rtk_minus_gen2_m)
B1_OFFSET = float(np.median(P.loc[bare1, "rtk_minus_gen2_m"])) if bare1.sum() else np.nan
print(f"B1 PROVISIONAL TIE: {int(bare1.sum())} bare flat FIX points under B1 read {B1_OFFSET:+.3f} m against gen2 (median); removed from every B1 point as a provisional constant")
P.loc[P.base == "B1", "h_ell_tied_m"] = P.loc[P.base == "B1", "h_ell_m"] - B1_OFFSET; P.loc[P.base == "B1", "tie_source"] = f"PROVISIONAL: bare-ground offset from gen2 under B1 ({B1_OFFSET:+.3f} m removed)"
P["navd88_tied_m"] = P.navd88_m - (P.h_ell_m - P.h_ell_tied_m)
# CORS vs B2 cross-check: nearest B2 FIX point within 3 m of each CORS point, raw ellipsoid-height difference (B2 minus CORS); the benchmark says +B2_OFFSET
cors = P[P.base == "CORS"]; b2f = P[(P.base == "B2") & P.status.str.startswith("FIX") & (P.group != "benchmark_DR9925")]
pairs = []
for _, c in cors.iterrows():
    d_ = np.hypot(b2f.easting - c.easting, b2f.northing - c.northing); k = d_.idxmin()
    if d_[k] <= 3.0: pairs.append(dict(cors_pt=c.point, b2_group=b2f.loc[k, "group"], b2_pt=b2f.loc[k, "point"], sep_m=round(float(d_[k]), 2), slope_deg=round(float(c.slope_deg), 1) if np.isfinite(c.slope_deg) else np.nan, b2_minus_cors_m=round(float(b2f.loc[k, "h_ell_m"] - c.h_ell_m), 3)))
PAIRS = pd.DataFrame(pairs)
print(f"\nCROSS-CHECK, B2 against CORS on the same ground (pairs within 3 m; the 3 m is mine): {len(PAIRS)} pairs")
if len(PAIRS): print(PAIRS.to_string(index=False)); print(f"  median B2 minus CORS {PAIRS.b2_minus_cors_m.median():+.3f} m, NMAD {1.4826 * np.median(np.abs(PAIRS.b2_minus_cors_m - PAIRS.b2_minus_cors_m.median())):.3f} m; benchmark tie said {B2_OFFSET:+.3f} m")
P["tied_minus_gen2_m"] = P.navd88_tied_m - P.gen2_z_m; P["tied_minus_gen1_m"] = P.navd88_tied_m - P.gen1_z_m
P = P[["group", "people", "device", "serial", "point", "t_start", "t_hhmm", "base", "baseline_m", "status", "lon", "lat", "easting", "northing", "easting_tied", "northing_tied", "h_ell_m", "h_ell_tied_m", "tie_source", "navd88_m", "navd88_tied_m", "elev_app_m", "antenna_m", "rms_e", "rms_n", "rms_v", "tilt_deg", "tile", "gen2_z_m", "gen1_z_m", "dod_2021_2008_m", "rtk_minus_gen2_m", "rtk_minus_gen1_m", "tied_minus_gen2_m", "tied_minus_gen1_m", "slope_deg", "lowveg_gen2", "clubb_floodplain", "penetration_failure_2008", "source"]]
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
print("\nBARE-GROUND CHECK (slope <= 3 deg, lowveg <= 0.05, FIX): median RTK minus gen2 per base, the provisional vertical tie of each base to the 2021 lidar:")
bare = fix[(fix.slope_deg <= 3) & (fix.lowveg_gen2 <= 0.05)]
for b, s in bare.groupby("base"):
    v = s.rtk_minus_gen2_m.dropna() * 1000; print(f"  {b}: {np.median(v):+.0f} mm (NMAD {1.4826 * np.median(np.abs(v - np.median(v))):.0f}, n {len(v)}; groups {sorted(s.group.unique())})")
print("  (the 0.05 lowveg and 3 deg cut-offs are MINE, for the check only)")

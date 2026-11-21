#!/usr/bin/env python
"""Error budget of the tied RTK ground heights (Andy, 2026-11-21: propagate the receiver error, the tilt compensation,
the benchmark tie and the slope). Per FIX point, one-sigma terms in meters:
  receiver       sv = Elevation RMS, sh = hypot(Easting RMS, Northing RMS), the Emlid solution's own formal precision
  tilt           RS3 (IMU tilt compensation): residual angle error d = 0.5 deg (MINE; Emlid publishes no figure we hold):
                 sh_tilt = L sin d, sv_tilt = L (1 - cos d), L = 2.134 m; the "compensation off" bound uses the logged tilt
                 itself: sh = L sin(tilt), sv = L (1 - cos tilt). RS+/RS2+ (no IMU, bubble-leveled): d = 0.5 deg (MINE).
  tie            B-based points only, common to all of them (a shift, not scatter): vertical from the datasheet network
                 accuracy 1.49 cm at 95 % (-> 0.0076 sigma), the two 60 s occupations (8 mm apart; 0.013 RMS each),
                 combined 0.011 m; horizontal 0.03 m (MINE: the independent CORS pair agreed with the benchmark tie
                 to +0.031 E / -0.011 N). CORS points: no tie term.
  slope          a horizontal error e_h on a slope theta displaces the point along the surface: vertical error at the
                 recorded position e_h tan(theta) (if along the fall line; / sqrt 2 for an isotropic error, both given)
Combined: vertical at the recorded position  sz = sqrt(sv_tot^2 + (sh_tot tan theta / sqrt 2)^2)  and its fall-line bound;
          slope-normal distance to a surface  sn = sqrt((sv_tot cos theta)^2 + (sh_tot sin theta / sqrt 2)^2).
Slope per point from the compiled table (5 m 2021 grid). Observed scatter for comparison: the NMAD of lidar floor minus
RTK, 0.085 m on the bluff (1 m stacks), 0.031 m on the floodplain sward (DEM minus RTK).
Outputs  data/field/2026-10-25/rtk_error_budget.csv (per point), rtk_error_budget_log.txt
"""
import os
import numpy as np, pandas as pd
ROOT = os.path.dirname(os.path.abspath(__file__)) + "/../.."; os.chdir(ROOT)
LOG = open("data/field/2026-10-25/rtk_error_budget_log.txt", "w")
def say(*a):
    s = " ".join(str(x) for x in a); print(s); LOG.write(s + "\n")
L = 2.134; D_RES = np.radians(0.5); TIE_V = float(np.sqrt((0.0149 / 1.96) ** 2 + (0.013 / np.sqrt(2)) ** 2 + (0.008 / 2) ** 2)); TIE_H = 0.03
P = pd.read_csv("data/field/2026-10-25/rtk_points_all.csv"); F = P[(P.status == "FIX") & P.group.str.contains("floodplain|hillslope")].copy()
F["sv_rec"] = F.rms_v; F["sh_rec"] = np.hypot(F.rms_e, F.rms_n); F["theta"] = np.radians(F.slope_deg.fillna(0))
rs3 = F.device.str.contains("RS3"); tilt = np.radians(F.tilt_deg.fillna(0))
F["sh_tilt"] = L * np.sin(D_RES); F["sv_tilt"] = L * (1 - np.cos(D_RES))
F["sh_tilt_off"] = np.where(rs3, L * np.sin(tilt), L * np.sin(D_RES)); F["sv_tilt_off"] = np.where(rs3, L * (1 - np.cos(tilt)), L * (1 - np.cos(D_RES)))
F["sv_tie"] = np.where(F.base == "B", TIE_V, 0.0); F["sh_tie"] = np.where(F.base == "B", TIE_H, 0.0)
def combine(sv_cols, sh_cols, tag):
    sv = np.sqrt(sum(F[c] ** 2 for c in sv_cols)); sh = np.sqrt(sum(F[c] ** 2 for c in sh_cols))
    F[f"sv_{tag}"] = sv; F[f"sh_{tag}"] = sh
    F[f"sz_{tag}"] = np.sqrt(sv ** 2 + (sh * np.tan(F.theta) / np.sqrt(2)) ** 2); F[f"sz_{tag}_fallline"] = np.sqrt(sv ** 2 + (sh * np.tan(F.theta)) ** 2)
    F[f"sn_{tag}"] = np.sqrt((sv * np.cos(F.theta)) ** 2 + (sh * np.sin(F.theta) / np.sqrt(2)) ** 2)
combine(["sv_rec"], ["sh_rec"], "rec"); combine(["sv_rec", "sv_tilt"], ["sh_rec", "sh_tilt"], "rec_tilt"); combine(["sv_rec", "sv_tilt", "sv_tie"], ["sh_rec", "sh_tilt", "sh_tie"], "all"); combine(["sv_rec", "sv_tilt_off", "sv_tie"], ["sh_rec", "sh_tilt_off", "sh_tie"], "all_tiltoff")
F.to_csv("data/field/2026-10-25/rtk_error_budget.csv", index=False, float_format="%.4f")
say(f"terms (1 sigma, m): tilt residual angle {np.degrees(D_RES):.1f} deg -> horizontal {L * np.sin(D_RES):.4f}, vertical {L * (1 - np.cos(D_RES)):.5f}; tie vertical {TIE_V:.4f}, tie horizontal {TIE_H:.3f} (common mode for B points)")
F["where"] = np.where(F.group.str.startswith("hillslope"), "hillslope", "floodplain")
pd.set_option("display.width", 250)
cols = ["sv_rec", "sh_rec", "sz_rec", "sn_rec", "sz_rec_tilt", "sz_all", "sz_all_fallline", "sn_all", "sz_all_tiltoff", "sn_all_tiltoff"]
say("\nmedians by group (m): sv_rec/sh_rec = receiver; sz = vertical error at the recorded position; sn = slope-normal distance error; _rec receiver only, _rec_tilt + tilt residual, _all + tie, _fallline = horizontal error along the fall line, _tiltoff = RS3 compensation assumed OFF (full lean enters)")
say(F.groupby("group")[["slope_deg", "tilt_deg"] + cols].median().round(3).to_string())
say("\n90th percentiles by place:"); say(F.groupby("where")[cols].quantile(0.9).round(3).to_string())
say("\nshare of the combined vertical variance (sz_all^2) by term, median over points: receiver vertical, tilt vertical, tie vertical, slope-converted horizontal (receiver + tilt + tie):")
for w, q in F.groupby("where"):
    tot = q.sz_all ** 2; hs = (q.sh_all * np.tan(q.theta) / np.sqrt(2)) ** 2
    say(f"   {w:10s}: receiver v {100 * (q.sv_rec ** 2 / tot).median():4.0f} %, tilt v {100 * (q.sv_tilt ** 2 / tot).median():4.0f} %, tie v {100 * (q.sv_tie ** 2 / tot).median():4.0f} %, horizontal via slope {100 * (hs / tot).median():4.0f} %  (of which tie horizontal {100 * ((q.sh_tie * np.tan(q.theta) / np.sqrt(2)) ** 2 / tot).median():.0f} %)")
say("\nobserved for comparison: lidar floor minus RTK on the bluff NMAD 0.085 m (1 m stacks), 0.103 slope-normal DEM minus RTK; floodplain DEM minus RTK NMAD 0.031 m (which includes the sward's own 2-3 cm texture)")
LOG.close()

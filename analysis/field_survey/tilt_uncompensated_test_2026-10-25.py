"""TEST (Andy, 2026-11-17): treat the two RS3 exports as if tilt compensation had been OFF, i.e. as if the
exported ground were (antenna phase center) minus a VERTICAL 2.134 m at the antenna's horizontal position.
The tip would then be 2.134*(1 - cos(tilt)) HIGHER than the exported ground and 2.134*sin(tilt) away
horizontally in the lean direction, which the export does not record. The vertical part is applied to
every RS3 point (bluff lines from the compiled table, tied; B1 toe points from the raw exports) and lidar
minus RTK is shown before and after; the horizontal part is bounded by |lateral| * tan(slope). A regression
of lidar-minus-RTK on the predicted vertical correction within each team has slope ~1 if compensation was
off and ~0 if it was on. No parameter is chosen here."""
import pandas as pd, numpy as np, sys
AH = 2.134  # antenna height entered on both RS3s (2.000 m pole + 0.134 phase-center offset)
def vcorr(t): return AH * (1 - np.cos(np.radians(t)))
def lat(t): return AH * np.sin(np.radians(t))
P = pd.read_csv("data/field/2026-10-25/rtk_points_all.csv")
H = P[P.group.str.startswith("hillslope") & (P.status == "FIX")].copy()
H["vcorr_m"] = vcorr(H.tilt_deg); H["lateral_m"] = lat(H.tilt_deg)
H["l21mr_mm"] = 1000 * H.gen2_minus_rtk_m                       # 2021 lidar minus tied RTK, as compiled
H["l21mr_tiltoff_mm"] = H.l21mr_mm - 1000 * H.vcorr_m          # tip higher -> lidar minus RTK smaller
H["hbound_mm"] = 1000 * H.lateral_m * np.tan(np.radians(H.slope_deg))
pd.set_option("display.width", 250)
print("INPUT: data/field/2026-10-25/rtk_points_all.csv (tied, antenna-checked), hillslope FIX points; tilt_deg from the RS3 exports.")
print("COLUMNS: vcorr_m = 2.134*(1-cos tilt), the tip's height above the exported ground IF uncompensated; lateral_m = 2.134*sin tilt;")
print("         l21mr_mm = 2021 lidar minus RTK as compiled; l21mr_tiltoff_mm = the same after raising the RTK ground by vcorr;")
print("         hbound_mm = lateral_m * tan(slope at the point), the largest vertical effect the unrecorded lean direction could have.\n")
print(H[["group", "point", "tilt_deg", "vcorr_m", "lateral_m", "slope_deg", "hbound_mm", "l21mr_mm", "l21mr_tiltoff_mm"]].round(3).to_string(index=False))
print("\nBLUFF LINES, 2021 lidar minus RTK, FIX, mm [median, NMAD, n]:")
def nmad(x): x = np.asarray(x); return 1.4826 * np.median(np.abs(x - np.median(x)))
for g, s in H.groupby("group"):
    print(f"  {g:16s} as compiled {s.l21mr_mm.median():+6.0f}, {nmad(s.l21mr_mm):4.0f}, {len(s)} | tilt-off hypothesis {s.l21mr_tiltoff_mm.median():+6.0f}, {nmad(s.l21mr_tiltoff_mm):4.0f} | tilt median {s.tilt_deg.median():.1f} max {s.tilt_deg.max():.1f} deg | vcorr median {1000*s.vcorr_m.median():.0f} max {1000*s.vcorr_m.max():.0f} mm | hbound median {s.hbound_mm.median():.0f} max {s.hbound_mm.max():.0f} mm")
b = H[H.group == "hillslope_T4a"]; m = H[H.group == "hillslope_RS3"]
print(f"  RS3 (T3a) minus T4a, medians: as compiled {m.l21mr_mm.median() - b.l21mr_mm.median():+.0f} mm; tilt-off hypothesis {m.l21mr_tiltoff_mm.median() - b.l21mr_tiltoff_mm.median():+.0f} mm")
print("\nREGRESSION of lidar-minus-RTK (mm) on the predicted vertical correction (mm), per team (slope ~ +1 if uncompensated, ~0 if compensated; OLS and Theil-Sen):")
from scipy.stats import theilslopes, linregress
for g, s in H.groupby("group"):
    x = 1000 * s.vcorr_m.values; y = s.l21mr_mm.values
    o = linregress(x, y); t = theilslopes(y, x)
    print(f"  {g:16s} OLS slope {o.slope:+.2f} (se {o.stderr:.2f}, r {o.rvalue:+.2f}, n {len(s)}); Theil-Sen {t[0]:+.2f} (95 % {t[2]:+.2f} to {t[3]:+.2f}); vcorr range {x.min():.0f}-{x.max():.0f} mm")
# B1 pair from the raw exports (both rovers under B1; B1's error is common)
R = pd.read_csv("data/field/2026-10-25/b1_rover_pair.csv")
R["vcorr_m"] = vcorr(R["Tilt angle"]); R["lateral_m"] = lat(R["Tilt angle"])
print("\nB1 PAIR (both rovers under the same base; from b1_rover_pair.csv):")
for tm, nm in (("T4a", 4), ("T3a", 1)):
    r = R[(R.team == tm) & (R.Name == nm)].iloc[0]
    print(f"  {tm} {nm}: tilt {r['Tilt angle']:.1f} deg, vcorr {1000*r.vcorr_m:+.0f} mm, lateral {r.lateral_m:.2f} m on {r.slope:.1f} deg ground (hbound {1000*r.lateral_m*np.tan(np.radians(r.slope)):.0f} mm); exported ell. h {r['Ellipsoidal height']:.3f} -> tilt-off tip {r['Ellipsoidal height'] + r.vcorr_m:.3f}")
b4 = R[(R.team == "T4a") & (R.Name == 4)].iloc[0]; m1 = R[(R.team == "T3a") & (R.Name == 1)].iloc[0]
dz = b4.z21_raw - m1.z21_raw
before = (b4["Ellipsoidal height"] - m1["Ellipsoidal height"]) - dz
after = before + (b4.vcorr_m - m1.vcorr_m)
print(f"  T4a minus T3a on the same ground (RTK dh minus lidar dz {1000*dz:+.0f} mm): as exported {1000*before:+.0f} mm; tilt-off hypothesis {1000*after:+.0f} mm")
print("\nREADING: under the tilt-off hypothesis the B1 pair closes if 'after' is near 0 while the bluff gap changes by the difference of the teams' median vcorr; the regression slopes say which hypothesis the bluff data themselves prefer.")
print("\nSCATTER vs the lean bound: if uncompensated, |residual from the team median| should grow with hbound (Spearman, per team):")
from scipy.stats import spearmanr
for g, s in H.groupby("group"):
    r = spearmanr(s.hbound_mm, np.abs(s.l21mr_mm - s.l21mr_mm.median()))
    print(f"  {g:16s} rho {r.correlation:+.2f}, p {r.pvalue:.2f}, n {len(s)}")

#!/usr/bin/env python3
"""Save per-return ANGLE information for gen1 CSF cloth ground returns (incidence to local
surface normal, signed scan angle, local slope, forest-floor elevation d, cell, flight line,
stratum), and (elba only) plot the SLOPE DEPENDENCY of the forest-floor elevation.

Incidence reconstruction validated in incidence_angle.py (flat ground -> |scan angle|).
Generalized over tiles: grid (origin/res) read from the tile's meta/corrections JSON; the
old penetration/floodplain/core strata are OPTIONAL (zeros if absent). The incidence and
d_mm math is identical across tiles. Writes data/derived/<tile>/gen1_csf_angles.npz.

    env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python analysis/ridgelines/gen1_save_angles_slope.py \
        [tile=elba_fulldensity] [csf=data/csf_cache/elba.las]
"""
import sys, os, json, numpy as np, laspy, math
from scipy.ndimage import distance_transform_edt
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

TILE = sys.argv[1] if len(sys.argv) > 1 else "elba_fulldensity"
CSF  = sys.argv[2] if len(sys.argv) > 2 else "data/csf_cache/elba.las"
D = f"data/derived/{TILE}"

def _grid(tile):                                        # (X0,Y0,NX,NY,RES) from tile meta/corrections
    for fn in ("meta.json", "corrections_geoid.json", "corrections.json"):
        p = f"data/derived/{tile}/{fn}"
        if os.path.exists(p):
            j = json.load(open(p)); b = j["bounds"]; r = float(j.get("res") or j.get("res_m"))
            nx = int(j.get("nx") or round((b[2]-b[0])/r)); ny = int(j.get("ny") or round((b[3]-b[1])/r))
            return b[0], b[1], nx, ny, r
    raise SystemExit(f"no grid meta for {tile}")

X0, Y0, NX, NY, RES = _grid(TILE)
Zg = np.load(f"{D}/z_after.npy"); Zf = Zg.copy(); m = ~np.isfinite(Zf)
if m.any(): Zf = Zf[tuple(distance_transform_edt(m, return_distances=False, return_indices=True))]
gy, gx = np.gradient(Zf, RES); slope_deg = np.degrees(np.arctan(np.hypot(gx, gy)))
gxf = gx.ravel(); gyf = gy.ravel(); nnorm = np.sqrt(gxf**2 + gyf**2 + 1.0); Zflat = Zf.ravel()

def _opt(name, default):                                # optional strata inputs (elba has them; elbaext may not)
    p = f"{D}/{name}.npy"; return np.load(p).ravel() if os.path.exists(p) else default
g2pen = _opt("penetration", None)
fld   = _opt("floodplain_mask", None); fld = fld.astype(bool) if fld is not None else None
core  = _opt("core_forest", None); copen = _opt("core_open", None)

las = laspy.read(CSF)
x = np.asarray(las.x, np.float64); y = np.asarray(las.y, np.float64); z = np.asarray(las.z, np.float64)
try:                                                    # scan angle -> DEGREES, both formats
    sa = np.asarray(las.scan_angle).astype(float) * 0.006      # PF6+ (0.006 deg units)
except Exception:
    sa = np.asarray(las.scan_angle_rank).astype(float)         # PF<=5 (integer degrees, 0=nadir)
psid = np.asarray(las.point_source_id); gt = np.asarray(las.gps_time)
ix = ((x - X0) / RES).astype(np.int64); iy = ((y - Y0) / RES).astype(np.int64)
ing = (ix >= 0) & (ix < NX) & (iy >= 0) & (iy < NY); cell = np.where(ing, iy*NX + ix, 0)

bhx = np.zeros(len(x)); bhy = np.zeros(len(x))
for pl in np.unique(psid):
    m = psid == pl
    vx = np.polyfit(gt[m], x[m], 1)[0]; vy = np.polyfit(gt[m], y[m], 1)[0]; H = math.atan2(vy, vx)
    cxu = np.array([-math.sin(H), math.cos(H)])
    cross = (x[m]-x[m].mean())*cxu[0] + (y[m]-y[m].mean())*cxu[1]
    sgn = np.sign(np.corrcoef(cross, sa[m])[0, 1])
    bhx[m] = -np.sign(sa[m])*sgn*cxu[0]; bhy[m] = -np.sign(sa[m])*sgn*cxu[1]
th = np.radians(np.abs(sa)); bx = np.sin(th)*bhx; by = np.sin(th)*bhy; bz = np.cos(th)
inc = np.degrees(np.arccos(np.clip((bx*(-gxf[cell]) + by*(-gyf[cell]) + bz)/nnorm[cell], -1, 1)))
slp = slope_deg.ravel()[cell]
xc = X0 + ((cell % NX)+0.5)*RES; yc = Y0 + ((cell // NX)+0.5)*RES
d = (z - (Zflat[cell] + gxf[cell]*(x-xc) + gyf[cell]*(y-yc))) * (1.0/nnorm[cell]) * 1000  # mm

# stratum code: 1 forest, 2 farmland(open), 0 other ; plus core flags (all 0/False if inputs absent)
strat = np.zeros(len(x), np.int8)
if g2pen is not None and fld is not None:
    strat[ing & ((g2pen[cell] < 0.25) & ~fld[cell])] = 1
    strat[ing & ((g2pen[cell] >= 0.45) & ~fld[cell])] = 2
cf = (core[cell] & ing) if core is not None else np.zeros(len(x), bool)
co = (copen[cell] & ing) if copen is not None else np.zeros(len(x), bool)
np.savez_compressed(f"{D}/gen1_csf_angles.npz",
    incidence=inc.astype(np.float32), scan_angle=sa.astype(np.float32), slope=slp.astype(np.float32),
    d_mm=d.astype(np.float32), cell=cell.astype(np.int32), point_source_id=psid.astype(np.int32),
    stratum=strat, core_forest=cf, core_open=co, in_grid=ing)
print(f"saved {D}/gen1_csf_angles.npz  (n=%d returns, grid {NX}x{NY})" % len(x))

# --- SLOPE DEPENDENCY plot (elba only; needs the penetration strata) ---
if TILE == "elba_fulldensity" and g2pen is not None:
    F = ing & (strat == 1); O = ing & (strat == 2)
    def bin_med(mask, xv, nb=10, lo=0, hi=40):
        e = np.linspace(lo, hi, nb+1); mx=[];my=[];mc=[]
        for i in range(nb):
            b = mask & (xv >= e[i]) & (xv < e[i+1])
            if b.sum() < 300: continue
            mx.append((e[i]+e[i+1])/2); my.append(np.median(d[b])); mc.append(b.sum())
        return np.array(mx), np.array(my), mc
    sx, sy, sc = bin_med(F, slp)
    print("\nforest-floor d vs SLOPE (forest CSF ground):")
    for a, b, c in zip(sx, sy, sc): print(f"  slope {a:4.0f}deg: d {b:+7.1f} mm  n={c:,}")
    print(f"corr(d, slope) forest = {np.corrcoef(d[F], slp[F])[0,1]:+.3f}")
    fig, ax = plt.subplots(1, 2, figsize=(14, 6))
    sxo, syo, _ = bin_med(O, slp, nb=8, lo=0, hi=20)
    ax[0].plot(sx, sy, "C0o-", label="forest"); ax[0].plot(sxo, syo, "C2s-", label="farmland")
    ax[0].set_xlabel("slope (deg)"); ax[0].set_ylabel("forest-floor elevation d (mm)")
    ax[0].set_title("floor elevation vs SLOPE"); ax[0].legend(); ax[0].grid(alpha=.3)
    ix2, iy2, _ = bin_med(F, inc, nb=12, lo=0, hi=45)
    ax[1].plot(ix2, iy2, "C3o-", label="forest")
    ax[1].set_xlabel("incidence angle to surface (deg)"); ax[1].set_ylabel("forest-floor elevation d (mm)")
    ax[1].set_title("floor elevation vs INCIDENCE (beam vs surface)"); ax[1].legend(); ax[1].grid(alpha=.3)
    fig.suptitle("gen1 CSF forest-floor elevation: slope and incidence dependency", y=1.0)
    fig.savefig("figures/refdatum/gen1_floor_vs_slope_incidence.png", dpi=130, bbox_inches="tight"); plt.close(fig)
    print("\nwrote figures/refdatum/gen1_floor_vs_slope_incidence.png")

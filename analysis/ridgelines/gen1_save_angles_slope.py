#!/usr/bin/env python3
"""Save per-return ANGLE information for gen1 CSF cloth ground returns (incidence to local
surface normal, signed scan angle, local slope, forest-floor elevation d, cell, flight line,
stratum), and (elba only) plot the SLOPE DEPENDENCY of the forest-floor elevation.

Incidence reconstruction validated in incidence_angle.py (flat ground -> |scan angle|).
Generalized over tiles: grid (origin/res) read from the tile's meta/corrections JSON; the
cover and floodplain strata are REQUIRED unless you say otherwise. The incidence and d_mm
math is identical across tiles. Writes data/derived/<tile>/gen1_csf_angles.npz.

The `stratum` column (1 forest, 2 farmland/open) came from `penetration.npy` until
2026-09-05. That layer is RETIRED -- it correlates -0.84 with SCAN ANGLE, which made it
unsafe as a stratum inside a beam-angle analysis, the exact use it had here. The column is
now built from the PyForestScan masks forest_pfs/open_pfs, which this script already carried
alongside it for cross-tile comparison. core_forest/core_open went with strata_core, whose
two classes penetration defined.

A MISSING STRATUM REFUSES; it is not silently substituted. These layers used to default to
zeros when absent, which turns "this tile has no forest_core layer" into a stratum that
reads as MEASURED AND EMPTY -- a table that looks like a finding of "no forest". Running
without them is allowed and often right, but it has to be stated:

    --without forest_pfs,open_pfs                    run without the cover strata
    --without all                                    run with none of the optional strata

    env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python analysis/ridgelines/gen1_save_angles_slope.py \
        [tile=elba_fulldensity] [csf=data/csf_cache/elba.las] [--without ...]
"""
import sys, os, json, numpy as np, laspy, math
from scipy.ndimage import distance_transform_edt
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

# --without <list> may appear anywhere; strip it before reading positionals.
_argv, _wo_vals = [], []
_it = iter(sys.argv[1:])
for _a in _it:
    if _a == "--without":
        _wo_vals.append(next(_it, ""))
    elif _a.startswith("--without="):
        _wo_vals.append(_a.split("=", 1)[1])
    else:
        _argv.append(_a)
TILE = _argv[0] if len(_argv) > 0 else "elba_fulldensity"
CSF  = _argv[1] if len(_argv) > 1 else "data/csf_cache/elba.las"
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

OPTIONAL_STRATA = ("floodplain_mask", "forest_pfs", "open_pfs", "canopy_cover_pfs")
_ABSENT = []
WITHOUT = set()
for _w in _wo_vals:
    WITHOUT |= set(OPTIONAL_STRATA) if _w.strip() == "all" else {t.strip() for t in _w.split(",") if t.strip()}
_bad = WITHOUT - set(OPTIONAL_STRATA)
if _bad:
    raise SystemExit(f"--without names layers that are not optional strata: {sorted(_bad)}; "
                     f"choose from {list(OPTIONAL_STRATA)} or 'all'")
if WITHOUT:
    print(f"  running WITHOUT, as stated: {sorted(WITHOUT)}", flush=True)


def _opt(name, default):
    """An absent stratum yields None, and the arrays that depend on it are OMITTED from the
    archive. It used to REFUSE unless named in --without, because absence was zero-filled --
    the requirement was compensating for the defect rather than fixing it. With the zero-fill
    gone, absence cannot manufacture a measurement, so it no longer needs announcing: a tile
    without cover simply produces a file without the cover arrays.

    --without still means something: exclude a layer that IS present, deliberately.
    """
    if name in WITHOUT:
        return default
    p = f"{D}/{name}.npy"
    if not os.path.exists(p):
        _ABSENT.append(name)
        return default
    return np.load(p).ravel()
fld   = _opt("floodplain_mask", None); fld = fld.astype(bool) if fld is not None else None
# PyForestScan cover: the cover definition used for cross-tile comparison, and now the only
# one. It replaced penetration, which correlates -0.84 with SCAN ANGLE and was therefore
# unsafe as a stratum inside a beam-angle analysis -- the very thing this archive feeds.
pfsf  = _opt("forest_pfs", None); pfso = _opt("open_pfs", None)
pfsc  = _opt("canopy_cover_pfs", None)

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

# stratum code: 1 forest, 2 farmland(open), 0 other. When a cover mask or the floodplain
# mask is absent this used to be written as ALL ZEROS -- i.e. "every return is 'other'", a
# claim, where the truth is "not computed". The array is OMITTED instead, so a consumer
# raises KeyError naming the key.
strat = None
if pfsf is not None and pfso is not None and fld is not None:
    strat = np.zeros(len(x), np.int8)
    strat[ing & (pfsf[cell].astype(bool) & ~fld[cell])] = 1
    strat[ing & (pfso[cell].astype(bool) & ~fld[cell])] = 2
# A boolean has no "unmeasured" value, so an ABSENT pfs mask must not be written as all
# False -- that is a claim of "no forest here" where the truth is "not measured", the same
# zero-fill defect as penetration.npy's 677 no-return cells. Absent masks are OMITTED from
# the archive, so a consumer raises KeyError naming the key instead of reading False.
# canopy_cover_pfs stays as NaN: for a float, NaN IS the unmeasured value.
_out = dict(
    incidence=inc.astype(np.float32), scan_angle=sa.astype(np.float32), slope=slp.astype(np.float32),
    d_mm=d.astype(np.float32), cell=cell.astype(np.int32), point_source_id=psid.astype(np.int32),
    in_grid=ing,
    canopy_cover_pfs=(np.where(ing, pfsc[cell], np.nan) if pfsc is not None
                      else np.full(len(x), np.nan)).astype(np.float32))
for _k, _v in (("stratum", strat),
               ("pfs_forest", pfsf[cell].astype(bool) & ing if pfsf is not None else None),
               ("pfs_open", pfso[cell].astype(bool) & ing if pfso is not None else None)):
    if _v is not None:
        _out[_k] = _v
_omitted = [k for k in ("stratum", "pfs_forest", "pfs_open") if k not in _out]
if _ABSENT:
    print(f"  layers absent in this tile: {sorted(set(_ABSENT))}", flush=True)
if _omitted:
    print(f"  OMITTED from the archive (not computed, NOT empty): {_omitted}", flush=True)
np.savez_compressed(f"{D}/gen1_csf_angles.npz", **_out)
print(f"saved {D}/gen1_csf_angles.npz  (n=%d returns, grid {NX}x{NY})" % len(x))

# --- SLOPE DEPENDENCY plot (elba only; needs the cover strata) ---
if TILE == "elba_fulldensity" and strat is not None:
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

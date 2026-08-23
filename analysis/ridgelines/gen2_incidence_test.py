#!/usr/bin/env python3
"""GEN2 (2021 USGS 3DEP, leaf-on) test of the beam-to-surface INCIDENCE-ANGLE control on
apparent forest-floor elevation. Reconstructs per-return incidence the SAME validated way as
gen1 (analysis/ridgelines/incidence_angle.py, gen1_save_angles_slope.py) and asks whether gen2
shows the same sign/magnitude of dependence: near-surface-perpendicular (low incidence) reads the
ground DEEPER, oblique (high incidence) reads it HIGHER.

DATA SOURCE (gen2) -- IMPORTANT:
  The MERGED products data/after/3dep2021_*.laz have scan_angle_rank AND gps_time ZEROED and
  point_source_id collapsed to 5 values (a PDAL/entwine merge stripped them). They are UNUSABLE
  for incidence reconstruction. The per-return angle/time/flight-line fields survive ONLY in the
  raw EPT tiles data/after/_ept_tiles/*.laz (3DEP EPT pull, EPSG:3857 web-mercator). This script
  therefore streams those EPT tiles, reprojects 3857 -> 26915 (UTM 15N, the grid CRS), and uses
  scan_angle_rank (point format 1, 1 deg/unit, signed, range ~+/-24 deg here).

  classification==2 is ground; 7 = noise (excluded). Reference bare earth = z_after.npy (this IS
  gen2's own bare earth, so d_mm centers near 0; the TEST is the RELATIVE dependence of d_mm on
  incidence, ideally at fixed slope).

METHOD (per gen1):
  1. Per flight line (point_source_id): heading H from fitting (x,y) vs gps_time; cross-track unit
     c=(-sinH,cosH); side of +scan_angle from sign of corr(cross-track pos, scan_angle) (|corr|
     ~0.9-0.98 here). Beam horizontal unit (ground->sensor) = -sign(scan_angle)*sgn*c.
     Beam vector b = sin(theta)*h + cos(theta)*z, theta=|scan_angle|.
  2. Surface normal n = (-gx,-gy,1)/|.| from z_after gradient.
  3. incidence = arccos(b.n), degrees. VALIDATE: on flat open ground incidence -> |scan angle|.
  d_mm = slope-normal distance of each ground return to the z_after plane, x1000 (as in gen1).

Run (needs PROJ_DATA for pyproj -> use conda lidar-icp env):
  PROJ_DATA=/home/awickert/anaconda3/envs/lidar-icp/share/proj \
    /home/awickert/anaconda3/envs/lidar-icp/bin/python analysis/ridgelines/gen2_incidence_test.py
"""
import os, glob, math
import numpy as np, laspy
from scipy.ndimage import distance_transform_edt
from pyproj import Transformer
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

# ---------------- grid + reference surfaces ----------------
NY, NX = 700, 508; X0, Y0 = 577492.8, 4882737.6; RES = 5.0
X1, Y1 = X0 + NX*RES, Y0 + NY*RES
DER = "data/derived/elba_fulldensity"
EPT = sorted(glob.glob("data/after/_ept_tiles/*.laz"))

Zg = np.load(f"{DER}/z_after.npy"); Zf = Zg.copy(); m = ~np.isfinite(Zf)
if m.any():
    Zf = Zf[tuple(distance_transform_edt(m, return_distances=False, return_indices=True))]
gy, gx = np.gradient(Zf, RES)                       # gy=d/dNorth, gx=d/dEast (m/m)
slope_deg = np.degrees(np.arctan(np.hypot(gx, gy)))
gxf = gx.ravel(); gyf = gy.ravel(); nnorm = np.sqrt(gxf**2 + gyf**2 + 1.0)
Zflat = Zf.ravel(); slpflat = slope_deg.ravel()
g2pen  = np.load(f"{DER}/penetration.npy").ravel()
fld    = np.load(f"{DER}/floodplain_mask.npy").astype(bool).ravel()
core_f = np.load(f"{DER}/core_forest.npy").ravel()   # core forest
core_o = np.load(f"{DER}/core_open.npy").ravel()     # core farmland (open)

trans = Transformer.from_crs("EPSG:3857", "EPSG:26915", always_xy=True)

# ---------------- pass over EPT tiles: load class-2 ground, reproject ----------------
def overlapping_tiles():
    keep = []
    for fn in EPT:
        with laspy.open(fn) as f: h = f.header
        cx = [h.x_min, h.x_max, h.x_min, h.x_max]; cy = [h.y_min, h.y_min, h.y_max, h.y_max]
        ux, uy = trans.transform(cx, cy)
        if max(ux) < X0 or min(ux) > X1 or max(uy) < Y0 or min(uy) > Y1: continue
        keep.append(fn)
    return keep

TILES = overlapping_tiles()
print(f"gen2 EPT tiles overlapping grid: {len(TILES)} of {len(EPT)}")

xs=[]; ys=[]; zs=[]; gts=[]; sas=[]; pss=[]
n_noise = 0
for k, fn in enumerate(TILES):
    las = laspy.read(fn)
    cl = np.asarray(las.classification)
    g = cl == 2
    n_noise += int((cl == 7).sum())
    if not g.any(): continue
    ux, uy = trans.transform(np.asarray(las.x)[g], np.asarray(las.y)[g])
    ux = np.asarray(ux); uy = np.asarray(uy)
    ing = (ux >= X0) & (ux < X1) & (uy >= Y0) & (uy < Y1)
    if not ing.any(): continue
    xs.append(ux[ing]); ys.append(uy[ing]); zs.append(np.asarray(las.z)[g][ing])
    gts.append(np.asarray(las.gps_time)[g][ing])
    sas.append(np.asarray(las.scan_angle_rank).astype(np.float64)[g][ing])
    pss.append(np.asarray(las.point_source_id)[g][ing])
    if (k+1) % 500 == 0: print(f"  ...{k+1}/{len(TILES)} tiles")

x = np.concatenate(xs); y = np.concatenate(ys); z = np.concatenate(zs)
gt = np.concatenate(gts); sa = np.concatenate(sas); psid = np.concatenate(pss)
del xs, ys, zs, gts, sas, pss
print(f"gen2 class-2 ground returns inside grid: {len(x):,}  (excluded class-7 noise seen: {n_noise:,})")
print(f"gen2 scan_angle_rank range: {sa.min():.0f} .. {sa.max():.0f} deg  (n flight lines: {len(np.unique(psid))})")

ix = ((x - X0)/RES).astype(np.int64); iy = ((y - Y0)/RES).astype(np.int64)
cell = iy*NX + ix                                    # all in-grid by construction

# ---------------- per-flight-line heading + side sign ----------------
bhx = np.zeros(len(x)); bhy = np.zeros(len(x))
print("\nper-flight-line reconstruction (heading, side-sign):")
lines = np.unique(psid); nfit = 0
for pl in lines:
    m = psid == pl
    if m.sum() < 200:  # too few to fit a trajectory reliably
        continue
    vx = np.polyfit(gt[m], x[m], 1)[0]; vy = np.polyfit(gt[m], y[m], 1)[0]
    H = math.atan2(vy, vx); cxu = np.array([-math.sin(H), math.cos(H)])
    cross = (x[m]-x[m].mean())*cxu[0] + (y[m]-y[m].mean())*cxu[1]
    cc = np.corrcoef(cross, sa[m])[0, 1]
    sgn = np.sign(cc) if np.isfinite(cc) and cc != 0 else 1.0
    bhx[m] = -np.sign(sa[m])*sgn*cxu[0]; bhy[m] = -np.sign(sa[m])*sgn*cxu[1]
    nfit += 1
    if m.sum() > 100000:
        print(f"  line {pl:>5}: n={m.sum():>9,} heading={math.degrees(H):+4.0f}deg  corr(cross,scan)={cc:+.2f}")
print(f"  ({nfit} flight lines fit; {len(lines)-nfit} skipped for <200 ground pts)")

# ---------------- incidence + slope-normal d ----------------
th = np.radians(np.abs(sa))
bx = np.sin(th)*bhx; by = np.sin(th)*bhy; bz = np.cos(th)      # beam unit (ground->sensor)
cosi = (bx*(-gxf[cell]) + by*(-gyf[cell]) + bz*1.0)/nnorm[cell]
inc = np.degrees(np.arccos(np.clip(cosi, -1, 1)))
slp = slpflat[cell]
xc = X0 + ((cell % NX) + 0.5)*RES; yc = Y0 + ((cell // NX) + 0.5)*RES
d = (z - (Zflat[cell] + gxf[cell]*(x-xc) + gyf[cell]*(y-yc)))*(1.0/nnorm[cell])*1000.0  # mm

# ---------------- VALIDATION: flat open ground incidence == |scan angle| ----------------
flat = (g2pen[cell] >= 0.45) & (~fld[cell]) & (slp < 2.0)
print("\nVALIDATION (gen2) on flat OPEN ground (slope<2deg, farmland): incidence vs |scan angle|")
for lo, hi in [(0, 2), (4, 6), (8, 10), (12, 16), (16, 24)]:
    b = flat & (np.abs(sa) >= lo) & (np.abs(sa) < hi)
    if b.sum() < 500: continue
    print(f"  |scan| {lo:2d}-{hi:2d} deg: median incidence {np.median(inc[b]):5.1f} deg  (should ~{(lo+hi)/2:.0f})  n={b.sum():,}")
steep = (g2pen[cell] < 0.25) & (~fld[cell]) & (slp > 20)
if steep.sum() > 100:
    print(f"  STEEP forest (slope>20): median |scan| {np.median(np.abs(sa)[steep]):.1f} deg vs "
          f"median incidence {np.median(inc[steep]):.1f} deg (differ on slopes)  n={steep.sum():,}")

# ---------------- SAVE gen2 angles ----------------
strat = np.zeros(len(x), np.int8)                    # 1 forest, 2 farmland(open), 0 other
strat[(g2pen[cell] < 0.25) & ~fld[cell]] = 1
strat[(g2pen[cell] >= 0.45) & ~fld[cell]] = 2
outnpz = f"{DER}/gen2_csf_angles.npz"
np.savez_compressed(outnpz,
    incidence=inc.astype(np.float32), scan_angle=sa.astype(np.float32), slope=slp.astype(np.float32),
    d_mm=d.astype(np.float32), cell=cell.astype(np.int32), point_source_id=psid.astype(np.int32),
    stratum=strat, core_forest=core_f[cell], core_open=core_o[cell])
print(f"\nsaved {outnpz}  (n={len(x):,} gen2 ground returns)")

# ================= incidence dependence: tables + figure =================
INC_BANDS = [(0, 5), (5, 10), (10, 15), (15, 25), (25, 40)]
BAND_COL  = ["C0", "C1", "C2", "C3", "C4"]

def band_table(sel, title):
    print(f"\n{title}  (median gen2 d_mm per incidence band):")
    rows = []
    for lo, hi in INC_BANDS:
        b = sel & (inc >= lo) & (inc < hi)
        if b.sum() < 50:
            print(f"  inc {lo:2d}-{hi:2d} deg: (n={b.sum()}, too few)"); rows.append((lo, hi, np.nan, b.sum())); continue
        md = np.median(d[b]); rows.append((lo, hi, md, b.sum()))
        print(f"  inc {lo:2d}-{hi:2d} deg: median d {md:+7.1f} mm   n={b.sum():,}")
    # mm/deg slope of median-d vs band-center (finite bands only)
    valid = [(0.5*(lo+hi), md) for lo, hi, md, n in rows if np.isfinite(md)]
    if len(valid) >= 2:
        vc = np.array([v[0] for v in valid]); vm = np.array([v[1] for v in valid])
        s = np.polyfit(vc, vm, 1)[0]
        print(f"  --> gen2 d vs incidence: {s:+.2f} mm/deg (linear fit through band medians)")
        # also report the low->high span across the 0-5 vs 15-25 (surface-perp vs oblique) bands
        return rows, s
    return rows, np.nan

# stratifications requested
sel_combined = (core_o[cell] | core_f[cell])          # CORE farmland + CORE forest
sel_forest   =  core_f[cell]                          # CORE forest alone
sel_open     =  core_o[cell]                          # CORE farmland alone (context)

rows_comb, slope_comb = band_table(sel_combined, "STRATIFICATION 1: CORE farmland + CORE forest (combined)")
rows_for,  slope_for  = band_table(sel_forest,   "STRATIFICATION 2: CORE forest alone")
_,         slope_open = band_table(sel_open,     "(context) CORE farmland alone")

# fixed slope band 16-24 deg (the gen1 headline band), within core forest
sb = sel_forest & (slp >= 16) & (slp < 24)
print("\ngen2 CORE FOREST at FIXED slope 16-24 deg (median d_mm per incidence band):")
for lo, hi in INC_BANDS:
    b = sb & (inc >= lo) & (inc < hi)
    if b.sum() < 30:
        print(f"  inc {lo:2d}-{hi:2d} deg: (n={b.sum()}, too few)"); continue
    print(f"  inc {lo:2d}-{hi:2d} deg: median d {np.median(d[b]):+7.1f} mm   n={b.sum():,}")

# ---------------- FIGURE: PDFs of gen2 d_mm by incidence band, for both stratifications ----------------
bins = np.arange(-350, 251, 10)
def pdf_panel(ax, sel, title, slo=None, shi=None):
    s = sel if slo is None else sel & (slp >= slo) & (slp < shi)
    for (lo, hi), col in zip(INC_BANDS, BAND_COL):
        b = s & (inc >= lo) & (inc < hi)
        if b.sum() < 50: continue
        md = np.median(d[b])
        ax.hist(d[b], bins=bins, density=True, histtype="step", lw=1.8, color=col,
                label=f"inc {lo}-{hi}° (n={b.sum():,}, med {md:+.0f})")
        ax.axvline(md, color=col, ls=":", lw=1.0)
    ax.set_xlabel("gen2 ground-return elevation d (mm)"); ax.set_ylabel("density")
    ax.set_title(title, fontsize=10); ax.legend(fontsize=7); ax.grid(alpha=.3)

fig, ax = plt.subplots(2, 2, figsize=(13, 9))
pdf_panel(ax[0, 0], sel_combined, "GEN2  CORE farmland+forest (all slopes)")
pdf_panel(ax[0, 1], sel_forest,   "GEN2  CORE forest (all slopes)")
pdf_panel(ax[1, 0], sel_forest,   "GEN2  CORE forest, slope 8-16°", 8, 16)
pdf_panel(ax[1, 1], sel_forest,   "GEN2  CORE forest, slope 16-24°", 16, 24)
fig.suptitle("GEN2 (2021 3DEP, leaf-on) ground-return elevation PDFs vs beam-to-surface INCIDENCE angle",
             y=0.995, fontsize=12)
fig.tight_layout(rect=[0, 0, 1, 0.98])
out = "figures/refdatum/gen2_ground_pdf_vs_incidence.png"
fig.savefig(out, dpi=100, bbox_inches="tight"); plt.close(fig)
print(f"\nwrote {out}  ({fig.get_size_inches()[0]*100:.0f} px wide nominal)")

# stash key numbers for the markdown
np.savez(f"{DER}/gen2_incidence_summary.npz",
         inc_bands=np.array(INC_BANDS), rows_comb=np.array(rows_comb, float),
         rows_for=np.array(rows_for, float), slope_comb=slope_comb, slope_for=slope_for,
         slope_open=slope_open)
print("done.")

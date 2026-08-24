#!/usr/bin/env python3
"""Separate the gen1 cross-track offset signal into a per-flight-line VERTICAL REGISTRATION
term and a common scanner BORESIGHT (roll) term, and test for a terrain-ORIENTATION offset.

Uses the flight-line overlap (~49% of cells here). In each shared cell the between-line
offset DIFFERENCE cancels terrain, real change, and datum (all shared by the cell), leaving:
    d_lineA - d_lineB  =  (reg_A - reg_B)  +  b * (scan_A - scan_B)
so regressing the per-cell between-line offset difference on the between-line scan-angle
difference gives:  slope b = BORESIGHT roll (mm/deg, common to the sensor);
                   intercept = pure REGISTRATION offset between the two lines.
A least-squares solve over all pairwise intercepts recovers per-line registration offsets.

Part B tests whether the offset depends on terrain ASPECT (downslope azimuth) within a
fixed slope band on near-planar cells -- an orientation-based offset independent of flight line.

    env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python \
        analysis/ridgelines/per_flightline_offset.py [tile_dir]
"""
import sys, math, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

TILE = sys.argv[1] if len(sys.argv) > 1 else "data/derived/elba_fulldensity"
NX, X0, Y0, RES = 508, 577492.8, 4882737.6, 5.0
MIN_CELL_LINE = 3        # min returns for a stable per-(cell,line) mean
CURV_MAX = 0.002         # near-planar cells for the aspect test

df = pd.read_parquet(f"{TILE}/beam_offset_table.parquet")
df = df[df.in_grid.values].copy()
lines = sorted(df.point_source_id.unique())

# ---- per-line overview: heading, scan sampling, raw offset ----
print("PER-LINE OVERVIEW:")
print(f"  {'psid':>5s} {'n':>10s} {'heading':>8s} {'mean_scan':>9s} {'median_d':>9s}")
head = {}
for pid, g in df.groupby("point_source_id"):
    x = X0 + ((g.cell % NX) + 0.5) * RES; y = Y0 + ((g.cell // NX) + 0.5) * RES
    vx = np.polyfit(g.gps_time, x, 1)[0]; vy = np.polyfit(g.gps_time, y, 1)[0]
    head[pid] = math.degrees(math.atan2(vy, vx))
    print(f"  {pid:5d} {len(g):10,d} {head[pid]:+8.1f} {g.scan_angle.mean():+9.2f} {g.d_mm.median():+9.1f}")

# ---- A. overlap-cell between-line differences -> boresight + registration ----
cl = (df.groupby(["cell", "point_source_id"])
        .agg(d=("d_mm", "mean"), sc=("scan_angle", "mean"), nn=("d_mm", "size")).reset_index())
cl = cl[cl.nn >= MIN_CELL_LINE]
m = cl.merge(cl, on="cell", suffixes=("_a", "_b"))
m = m[m.point_source_id_a < m.point_source_id_b].copy()
m["dd"] = m.d_a - m.d_b; m["dsc"] = m.sc_a - m.sc_b

print(f"\nA. PAIRWISE between-line (overlap cells, per-cell means), Dd = reg_diff + boresight*Dscan:")
print(f"  {'pair':>9s} {'n_cells':>8s} {'median_Dd':>10s} {'boresight':>10s} {'registration':>13s} {'headings'}")
pair_int = []
for (a, b), g in m.groupby(["point_source_id_a", "point_source_id_b"]):
    if len(g) < 50: continue
    sl, ic = np.polyfit(g.dsc, g.dd, 1)
    rel = "opp" if abs((head[a] - head[b] + 180) % 360 - 180) > 90 else "same"
    print(f"  {a}-{b} {len(g):8,d} {g.dd.median():+10.1f} {sl:+9.2f}/d {ic:+12.1f}   {rel} ({head[a]:+.0f},{head[b]:+.0f})")
    pair_int.append((a, b, ic))
b_pool, ic_pool = np.polyfit(m.dsc, m.dd, 1)
print(f"  POOLED boresight roll b = {b_pool:+.2f} mm/deg   (pooled intercept {ic_pool:+.1f} mm)")

# per-line registration offsets: least squares on pairwise intercepts, mean-zero constraint
idx = {p: i for i, p in enumerate(lines)}; rows, rhs = [], []
for a, b, ic in pair_int:
    r = np.zeros(len(lines)); r[idx[a]] = 1; r[idx[b]] = -1; rows.append(r); rhs.append(ic)
rows.append(np.ones(len(lines))); rhs.append(0.0)          # sum(reg) = 0 gauge
reg, *_ = np.linalg.lstsq(np.array(rows), np.array(rhs), rcond=None)
print("  per-line REGISTRATION offset (mean-zero gauge):  " +
      "  ".join(f"{p}:{reg[idx[p]]:+.1f}mm" for p in lines))

# ---- B. terrain-ASPECT offset (near-planar cells, fixed slope band) ----
stab = df[df.curv_laplacian.abs() <= CURV_MAX]
pc = (stab.groupby("cell")
          .agg(d=("d_mm", "mean"), asp=("aspect_deg", "first"), slp=("slope", "first")).dropna())
band = pc[(pc.slp >= 8) & (pc.slp <= 15)]
print(f"\nB. OFFSET vs terrain ASPECT (near-planar |curv|<={CURV_MAX}, slope 8-15 deg, {len(band):,} cells):")
edges = np.arange(0, 361, 30)
print(f"  {'aspect(deg)':>12s} {'median_d(mm)':>12s} {'n':>8s}")
a = band.asp.to_numpy(); dd = band.d.to_numpy()
for i in range(len(edges) - 1):
    bb = (a >= edges[i]) & (a < edges[i + 1])
    if bb.sum() < 30: continue
    print(f"  {(edges[i]+edges[i+1])/2:12.0f} {np.median(dd[bb]):+12.1f} {bb.sum():8,d}")
# directional (cosine) fit: d = c0 + A cos(aspect - phi)
ar = np.radians(a); M = np.c_[np.ones_like(ar), np.cos(ar), np.sin(ar)]
c0, ca, cb = np.linalg.lstsq(M, dd, rcond=None)[0]
amp = math.hypot(ca, cb); phi = math.degrees(math.atan2(cb, ca)) % 360
print(f"  cosine fit: amplitude {amp:.1f} mm, peak at aspect {phi:.0f} deg, mean {c0:+.1f} mm "
      f"(amplitude = orientation-based offset)")

# ---- figure ----
fig, ax = plt.subplots(1, 2, figsize=(14, 6))
sedg = np.arange(-24, 26, 3)
sc_c, dd_m = [], []
xs = m.dsc.to_numpy(); ys = m.dd.to_numpy()
for i in range(len(sedg) - 1):
    bb = (xs >= sedg[i]) & (xs < sedg[i + 1])
    if bb.sum() < 100: continue
    sc_c.append((sedg[i] + sedg[i + 1]) / 2); dd_m.append(np.median(ys[bb]))
ax[0].hexbin(xs, ys, gridsize=50, bins="log", cmap="viridis", mincnt=1, extent=(-24, 24, -200, 200))
ax[0].plot(sc_c, dd_m, "C3o-", label="binned median")
xx = np.array([-24, 24]); ax[0].plot(xx, b_pool * xx + ic_pool, "w--", lw=2, label=f"fit b={b_pool:+.2f} mm/deg")
ax[0].axhline(0, color="k", lw=.6); ax[0].axvline(0, color="k", lw=.6)
ax[0].set_xlim(-24, 24); ax[0].set_ylim(-200, 200)
ax[0].set_xlabel("between-line scan-angle diff (deg)"); ax[0].set_ylabel("between-line offset diff (mm)")
ax[0].set_title("A. boresight (slope) vs registration (intercept)"); ax[0].legend()
ac, am = [], []
for i in range(len(edges) - 1):
    bb = (a >= edges[i]) & (a < edges[i + 1])
    if bb.sum() < 30: continue
    ac.append((edges[i] + edges[i + 1]) / 2); am.append(np.median(dd[bb]))
ax[1].plot(ac, am, "C0o-", label="median offset")
xf = np.linspace(0, 360, 200)
ax[1].plot(xf, c0 + amp * np.cos(np.radians(xf) - math.radians(phi)), "C3-",
           label=f"cos fit A={amp:.1f}mm @ {phi:.0f}deg")
ax[1].axhline(c0, color="k", lw=.6); ax[1].set_xlabel("terrain aspect (deg CW from N)")
ax[1].set_ylabel("median offset d (mm)"); ax[1].set_title("B. offset vs terrain aspect (slope 8-15 deg)")
ax[1].legend(); ax[1].grid(alpha=.3)
fig.suptitle("gen1 offset: per-flight-line registration + boresight, and terrain-aspect dependence (elba)", y=1.0)
fig.savefig("figures/refdatum/per_flightline_offset.png", dpi=130, bbox_inches="tight"); plt.close(fig)
print("\nwrote figures/refdatum/per_flightline_offset.png")

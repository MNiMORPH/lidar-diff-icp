#!/usr/bin/env python3
"""Does the correction chain converge in ONE STEP or need ITERATION, and is the terrain-aspect
offset really the gen1-vs-gen2 lateral misregistration?

Two corrections act on the raw gen1 offset table:
  - BORESIGHT roll b (mm/deg) -- from flight-line overlap (gen2 cancels; module).
  - LATERAL shift (sE, sN) -- estimated here by a robust plane fit of the per-cell offset in
    gen2-gradient space: d ~ 1000*(sE*gE + sN*gN), the signature of a horizontal misregistration.

Coupling test (both directions):
  * re-estimate the lateral shift AFTER removing boresight  -> does it move?
  * re-estimate boresight AFTER removing the lateral shift  -> does it move?
If both cross-effects are below their uncertainties, the chain is one-pass (triangular coupling);
otherwise it needs iteration.

Confirmation that aspect == lateral (non-circular):
  * the fitted shift magnitude should match the independently-known ~0.66 m elba lateral;
  * the aspect cosine amplitude should scale ~tan(slope) across slope bands (misregistration),
    not be slope-independent (which a sensor/illumination directional effect would be).

    env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python \
        analysis/ridgelines/boresight_lateral_coupling.py [tile_dir]
"""
import sys, math, numpy as np, pandas as pd
from lidar_diff_icp.boresight import estimate_boresight
from lidar_diff_icp.swathdiff import _robust_plane_fit          # reuse the package's Huber plane fit

TILE = sys.argv[1] if len(sys.argv) > 1 else "data/derived/elba_fulldensity"
NY, NX, RES = 700, 508, 5.0
CURV_MAX = 0.002

df = pd.read_parquet(f"{TILE}/beam_offset_table.parquet")
df = df[df.in_grid.values].copy()
Z = np.load(f"{TILE}/z_after.npy"); gN, gE = np.gradient(Z, RES)   # uphill grad: (north, east)
df["gE"] = gE.ravel()[df.cell.to_numpy()]; df["gN"] = gN.ravel()[df.cell.to_numpy()]

# ---- boresight (all in-grid returns, via module) ----
sol = estimate_boresight(df.cell.values, df.point_source_id.values,
                         df.scan_angle.values, df.d_mm.values, min_cell_line=3, min_pair_cells=50)
b0 = sol.b
print(f"boresight roll b0 = {b0:+.2f} +/- {sol.b_pair_std:.2f} mm/deg")

# ---- lateral shift from the offset's gradient dependence (robust plane fit, per-cell) ----
stab = df[df.curv_laplacian.abs() <= CURV_MAX].copy()

def lateral_fit(dcol, data):
    g = data.groupby("cell").agg(d=(dcol, "mean"), gE=("gE", "first"), gN=("gN", "first"),
                                 slp=("slope", "first"), asp=("aspect_deg", "first")).dropna()
    a, b, c = _robust_plane_fit(g.gE.to_numpy(), g.gN.to_numpy(), g.d.to_numpy())
    return a, b, c, g                                    # a,b in mm per unit gradient; shift = (a,b)/1000 m

a0, b0l, c0, graw = lateral_fit("d_mm", stab)
sh0 = math.hypot(a0, b0l) / 1000.0
print(f"lateral shift (from offset-gradient fit) = {sh0:.3f} m, "
      f"direction {math.degrees(math.atan2(b0l, a0)) % 360:.0f} deg  "
      f"(independent elba Nuth-Kaeaeb was ~0.66 m)")

# ---- COUPLING (both directions) ----
stab["d_b"] = stab.d_mm - b0 * stab.scan_angle
a1, b1, c1, _ = lateral_fit("d_b", stab)
sh1 = math.hypot(a1, b1) / 1000.0
d_lat_full = df.d_mm.to_numpy() - (a0 * df.gE.to_numpy() + b0l * df.gN.to_numpy())
b_after = estimate_boresight(df.cell.values, df.point_source_id.values,
                             df.scan_angle.values, d_lat_full, min_cell_line=3, min_pair_cells=50).b
print("\nCOUPLING:")
print(f"  lateral shift: raw {sh0:.3f} m  ->  after boresight removal {sh1:.3f} m  "
      f"(delta {abs(sh1-sh0)*1000:.1f} mm)")
print(f"  boresight    : raw {b0:+.2f}  ->  after lateral removal   {b_after:+.2f} mm/deg  "
      f"(delta {abs(b_after-b0):.2f} vs uncertainty {sol.b_pair_std:.2f})")
one_step = (abs(b_after - b0) < sol.b_pair_std) and (abs(sh1 - sh0) * 1000 < 50)
print(f"  => {'ONE STEP (coupling below uncertainty)' if one_step else 'ITERATION NEEDED'}")

# ---- confirm aspect == lateral: cosine amplitude scales ~tan(slope) ----
print("\nASPECT == LATERAL check (amplitude should scale ~tan(slope)):")
print(f"  {'slope band':>12s} {'tan(mid)':>8s} {'cos_amp(mm)':>12s} {'amp/tan':>8s} {'n':>7s}")
def cosamp(g):
    ar = np.radians(g.asp.to_numpy()); M = np.c_[np.ones_like(ar), np.cos(ar), np.sin(ar)]
    _, ca, cb = np.linalg.lstsq(M, g.d.to_numpy(), rcond=None)[0]
    return math.hypot(ca, cb)
for lo, hi in [(3, 8), (8, 15), (15, 25), (25, 40)]:
    gb = graw[(graw.slp >= lo) & (graw.slp < hi)]
    if len(gb) < 100: continue
    mid = math.tan(math.radians((lo + hi) / 2)); amp = cosamp(gb)
    print(f"  {lo:4d}-{hi:<4d}    {mid:8.2f} {amp:12.1f} {amp/mid:8.0f} {len(gb):7,d}")

# aspect amplitude before vs after removing the lateral model (should collapse)
graw2 = graw.copy(); graw2["d"] = graw2.d - (a0 * graw2.gE + b0l * graw2.gN)
band = lambda g: g[(g.slp >= 8) & (g.slp < 15)]
print(f"\n  aspect cos amplitude (slope 8-15): raw {cosamp(band(graw)):.1f} mm  ->  "
      f"after lateral removal {cosamp(band(graw2)):.1f} mm")

#!/usr/bin/env python3
"""Solve each gen1 flight line's HORIZONTAL shift against gen2 topography (a STAR),
and compare with the line-to-line chain that ``coreg.align_swaths`` solves.

Andy, 2026-09-18: "The horizontal shift is an interesting idea because we previously did
this as a later step comparing to topography. Try it."

WHY. On elbaext the swath network is a bare CHAIN -- 6 lines, 5 edges, zero loops, because
non-adjacent lines have disjoint bounding boxes. A chain has no redundancy: every line's
constants are the running sum of the edge observations, and `max |misclosure|` is 0.0000 mm
by construction, so the method cannot report its own error. Its `dx` accumulates
monotonically, every edge the same sign, to 1.3889 m.

Comparing each line DIRECTLY to gen2's topography makes the network a star instead: every
line gets an independent determination, so they can disagree, which is the whole point.
It reuses the estimator the pipeline already trusts -- `coreg.tie_polynomial(order=0)`,
the Nuth & Kaeaeb lateral shift -- applied per line instead of to the whole cloud.

THE CONTROL IS NOT OPTIONAL and runs first. Two shortcuts were tried and BOTH failed to
reproduce the shipped whole-cloud shift, which is why this script reads the CSF cache and
rebuilds the pipeline's own estimator:

    per-cell median of raw z, tie_polynomial   dx -0.1322   (shipped -0.7519)  FAILED
    plain least-squares NK on the beam table   dx +0.0860                      FAILED
    pipeline estimator, points NOT pre-aligned dx -0.1342                      FAILED
    pipeline estimator, align_swaths APPLIED   dx -0.7519   dy -0.1897         PASSES

The last line matches the shipped [-0.7519, -0.1896] to four decimals. The reason the
others fail is the finding: ``align_swaths`` runs INSIDE ``register_gen1``, BEFORE the
correction chain, so its per-line dx/dy are already in the points when ``LateralShift``
sees them.

The script REPORTS. It changes no product and no pipeline code.

    ./lidar-icp/bin/python scripts/perline_lateral_shift.py --tile elbaext_delong \
        --csf data/csf_cache/elbaext.las
"""
import argparse
import json

import numpy as np
from scipy.ndimage import gaussian_filter

from lidar_diff_icp import coreg, groundest, terrain

ap = argparse.ArgumentParser(description=__doc__,
                             formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--tile", required=True)
ap.add_argument("--csf", required=True, help="CSF-classified gen1 cache (pre-alignment)")
ap.add_argument("--valley-top-m", type=float, required=True,
                help="CITED per tile, never chosen here; elbaext = 230.0 (refcells.py)")
ap.add_argument("--ground-q", type=float, default=0.50)
ap.add_argument("--sn-smooth-cells", type=float, default=1.2)
A = ap.parse_args()

import laspy

D = f"data/derived/{A.tile}"
c = json.load(open(f"{D}/corrections.json"))
b = c["bounds"]; res = float(c["res_m"]); X0, Y0 = b[0], b[1]
Zref = np.load(f"{D}/z_after_differenced.npy"); ny, nx = Zref.shape
grid = (X0, Y0, res, nx, ny)

Zf = terrain.terrain_masks(Zref, res, valley_top_m=A.valley_top_m)["filled"]
Zreg = gaussian_filter(Zf, A.sn_smooth_cells)
plane = (Zreg.ravel(), np.gradient(Zreg, res, axis=1).ravel(),
         np.gradient(Zreg, res, axis=0).ravel())

f = laspy.read(A.csf)
ps = np.asarray(f.point_source_id)
x = np.asarray(f.x); y = np.asarray(f.y); z = np.asarray(f.z)
SW = c["per_swath_internal_alignment_dxdydz_m"]


def ground_of(xx, yy, zz):
    return groundest.estimate_ground(xx, yy, zz, grid, A.ground_q, "slope_normal",
                                     plane=plane)


def nk(xx, yy, zz):
    h = coreg.tie_polynomial(Zref, ground_of(xx, yy, zz), res, X0, Y0, order=0)
    return float(h["a"][0]), float(h["b"][0]), float(h["c"][0])


print("CONTROL -- reproduce the shipped whole-cloud lateral shift, or stop here")
X, Y, Z = x.copy(), y.copy(), z.copy()
for k, v in SW.items():
    m = ps == int(k)
    X[m] += v[0]; Y[m] += v[1]; Z[m] += v[2]
a1, b1, _ = nk(X, Y, Z)
want = c["cross_epoch_datum"]["horizontal_shift_m"]
print(f"  align_swaths applied, then NK : dx {a1:+.4f}  dy {b1:+.4f}")
print(f"  shipped                       : dx {want[0]:+.4f}  dy {want[1]:+.4f}")
ok = abs(a1 - want[0]) < 5e-4 and abs(b1 - want[1]) < 5e-4
print(f"  -> {'PASS' if ok else 'FAIL'}")
if not ok:
    raise SystemExit("control failed; per-line numbers would not be comparable")

a0, b0, _ = nk(x, y, z)
print(f"\nSame NK WITHOUT align_swaths' horizontal: dx {a0:+.4f}  dy {b0:+.4f}")
print(f"  -> align_swaths moves the cloud {a1 - a0:+.4f} m in x, which NK then removes")

print(f"\nPER-LINE, star against gen2 topography (raw points), vs the chain")
print(f"{'psid':>6}{'points':>11}{'star dx':>9}{'star dy':>9}{'star dz':>9}"
      f"{'chain dx':>10}{'chain dy':>9}{'chain dz':>9}")
star = {}
for p in sorted(set(ps.tolist())):
    m = ps == p
    star[p] = nk(x[m], y[m], z[m])
    v = SW[str(p)]
    print(f"{p:>6}{int(m.sum()):>11,}{star[p][0]:>9.4f}{star[p][1]:>9.4f}{star[p][2]:>9.4f}"
          f"{v[0]:>10.4f}{v[1]:>9.4f}{v[2]:>9.4f}")
sx = [v[0] for v in star.values()]; sy = [v[1] for v in star.values()]
cx = [SW[str(p)][0] for p in sorted(star)]; cy = [SW[str(p)][1] for p in sorted(star)]
print(f"\n  dx spread: star {max(sx)-min(sx):.4f} m   chain {max(cx)-min(cx):.4f} m")
print(f"  dy spread: star {max(sy)-min(sy):.4f} m   chain {max(cy)-min(cy):.4f} m")
print("\nREPORTED, NOT ADOPTED. Neither solution is checked against truth here, and the")
print("two disagree by up to ~1.5 m in dx. The star's value is that its lines CAN")
print("disagree; the chain's misclosure is 0.0000 mm by construction and cannot.")

# ---------------------------------------------------------------------------
# STAGE 2. Andy, 2026-09-18: "I was thinking of the entire set of swaths being
# horizontally co-registered but vertically mis-registered. But this was an assumption.
# It seems more principled to run Nuth & Kaeaeb by swath (against gen2), and then to see
# vertical misalignment afterwards. Indeed, some vertical misalignment could actually just
# be the result of horizontal misalignment."
#
# THE TEST. Correct each line HORIZONTALLY ONLY (three ways), then measure what vertical
# offset each line still has against gen2 on STABLE cells. If a horizontal fix shrinks the
# per-line vertical spread, then part of what align_swaths calls vertical misregistration
# was horizontal misregistration read through the surface gradient.
# ---------------------------------------------------------------------------
stable = np.load(f"{D}/stable.npy").astype(bool)


def vertical_by_line(dxdy):
    """Per-line median (gen1 - gen2) on stable cells, after a horizontal-only shift."""
    out = {}
    for p in sorted(set(ps.tolist())):
        m = ps == p
        ddx, ddy = dxdy(p)
        G = ground_of(x[m] + ddx, y[m] + ddy, z[m])
        ok = stable & np.isfinite(G) & np.isfinite(Zref)
        out[p] = float(np.median((G[ok] - Zref[ok]))) * 1000.0 if ok.sum() else float("nan")
    return out


print("\nSTAGE 2 -- per-line VERTICAL offset (mm, gen1 - gen2, median on stable cells)")
print("           after a HORIZONTAL-ONLY correction. No vertical term is applied.")
cases = [("no horizontal correction", lambda p: (0.0, 0.0)),
         ("align_swaths horizontal", lambda p: (SW[str(p)][0], SW[str(p)][1])),
         ("per-line star NK horizontal", lambda p: (star[p][0], star[p][1]))]
res_rows = {}
for lab, fn in cases:
    v = vertical_by_line(fn)
    res_rows[lab] = v
    vals = [v[p] for p in sorted(v)]
    print(f"\n  {lab}")
    print("   " + "".join(f"{p:>10}" for p in sorted(v)))
    print("   " + "".join(f"{v[p]:>10.1f}" for p in sorted(v)))
    print(f"    spread {max(vals)-min(vals):>8.1f} mm     sd {np.std(vals):>7.1f} mm")
print("\n  If the star row is TIGHTER, some 'vertical' misregistration was horizontal.")

# The FAIR method-vs-method comparison: each method's FULL solution (dx, dy AND dz), which
# is what the pipeline actually applies. The horizontal-only rows above isolate one
# component; no product is ever in that state.
print("\nSTAGE 3 -- each method's FULL per-line solution applied (dx, dy AND dz)")


def vertical_full(dxdydz):
    out = {}
    for p in sorted(set(ps.tolist())):
        m = ps == p
        ddx, ddy, ddz = dxdydz(p)
        G = ground_of(x[m] + ddx, y[m] + ddy, z[m] + ddz)
        ok = stable & np.isfinite(G) & np.isfinite(Zref)
        out[p] = float(np.median(G[ok] - Zref[ok])) * 1000.0 if ok.sum() else float("nan")
    return out


for lab, fn in (("align_swaths full (shipped)",
                 lambda p: (SW[str(p)][0], SW[str(p)][1], SW[str(p)][2])),
                ("per-line star NK full",
                 lambda p: (star[p][0], star[p][1], star[p][2]))):
    v = vertical_full(fn); vals = [v[p] for p in sorted(v)]
    print(f"\n  {lab}")
    print("   " + "".join(f"{p:>10}" for p in sorted(v)))
    print("   " + "".join(f"{v[p]:>10.1f}" for p in sorted(v)))
    print(f"    spread {max(vals)-min(vals):>8.1f} mm     sd {np.std(vals):>7.1f} mm")
print("\n  This is the comparison that decides which method registers better VERTICALLY.")
print("\n  ⚠️ CONFOUND, and it bounds what these spreads mean: each line covers a DIFFERENT")
print("  part of the tile, so its median offset against gen2 carries real spatial variation")
print("  in gen1-gen2 (terrain, vegetation, the correction surface) and not registration")
print("  alone. align_swaths ties lines on SHARED cells, which is apples-to-apples; this")
print("  metric is not. Footprints are identical across the cases above, so COMPARING the")
print("  cases is fair -- but the residual spread must NOT be read as pure misregistration.")

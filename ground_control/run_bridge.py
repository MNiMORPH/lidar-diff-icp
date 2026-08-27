"""The BRIDGE: from the surface the control measured to the surface our DoD runs on.

The control's residual describes the **delivered** product -- the vendor's own
``Surface Z`` at the mark.  Our DoD runs on a different surface: CSF-reclassified ground,
swath-aligned, geoid-shifted, gridded at 5 m by the slope-normal median.  The constant
that puts the delivered surface on NAVD88 is not the constant that puts ours there.

    bridge = z_delivered - z_ours          so   constant_ours = constant_delivered + bridge

Positive bridge = our surface sits BELOW the delivered one.

Two different gaps, which must not be pooled
--------------------------------------------
**A. estimator gap.**  Our tie estimator (radius ladder, order-2 local surface, quantile)
reading the *delivered* cloud, against the vendor's own ``Surface Z`` at the same mark.
Same points, different reduction.  Measurable at any mark whose tile is on disk.

**B. processing gap.**  Our gridded DoD surface -- CSF ground, swath alignment, geoid
shift -- against the vendor's ``Surface Z``.  This is the one the products actually need.
Measurable ONLY where our reconstructed grid exists, i.e. inside elba/elbaext.

They are different quantities and this script reports them apart.  A is not a proxy for B.

The geoid trap, and why gap B must not be read raw
--------------------------------------------------
Our gridded gen1 surface has been carried onto gen2's frame, which is **NAVD88(GEOID18)**.
The 2008 control is **NAVD88(GEOID03)**.  Differencing them raw compares two geoids and
reports the conversion as if it were a processing difference: at Elba that conversion is
``references.geoid_difference`` = +67.38 mm constant, and the raw gap B at the one open
mark inside our grid is -59.7 mm, almost all of which IS the conversion.  Gap B is
therefore reported both raw and geoid-corrected, and the corrected column is the one that
means anything.  Gap A needs no such term: the raw gen1 cloud and the 2008 control are
both on GEOID03, so it is one frame.

Recovering gen1's grid
----------------------
No ``z_before.npy`` is stored.  gen1's surface is recovered as ``z_after - dod``, verified
on ``elba_fulldensity`` where an independent gen1 grid also exists
(``z_before_absolute.npy`` minus its own datum constant): median +0.489 mm, NMAD 0.391 mm
over 341,239 cells.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

_HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(_HERE))
sys.path.insert(0, str(_HERE.parent))
sys.path.insert(0, str(_HERE.parent / "src"))

import lines as L  # noqa: E402
import same_line as S  # noqa: E402
from lidar_diff_icp.groundtruth import gen1_datum as G  # noqa: E402
from trust.provenance import Run  # noqa: E402


def read_grid_at(site_dir, dod_name, easting, northing, *, radius_m):
    """gen1's gridded surface AT a coordinate, by the repo's own tie estimator.

    A window median drifts on sloping terrain -- swept over 5 to 55 m it moved the answer
    at L1O101 by 72.3 mm and at L2T51 by 224.8 mm, so it was measuring the terrain, not
    the surface.  Instead the finite grid cells within ``radius_m`` are handed to
    ``groundtruth.tie.ground_elevation_at`` as points: the same order-2 local surface fit
    plus residual-quantile the tie machinery uses at a control mark, so grid and cloud are
    read the SAME WAY and the bridge differences like with like.
    """
    from lidar_diff_icp.groundtruth import tie
    d = json.loads((Path(site_dir) / "corrections.json").read_text())
    x0, y0, x1, y1 = d["bounds"]
    res = float(d["res_m"])
    za = np.load(Path(site_dir) / "z_after.npy")
    dod = np.load(Path(site_dir) / f"{dod_name}.npy")
    if za.shape != dod.shape:
        raise ValueError(f"{site_dir}: z_after {za.shape} != {dod_name} {dod.shape}")
    zb = za - dod
    ny, nx = zb.shape
    ix = int((easting - x0) // res)
    iy = int((northing - y0) // res)      # row 0 = south; pipeline uses +row = north
    if not (0 <= ix < nx and 0 <= iy < ny):
        return None, None, res, zb.shape
    h = int(np.ceil(radius_m / res)) + 1
    i0, i1 = max(ix - h, 0), min(ix + h + 1, nx)
    j0, j1 = max(iy - h, 0), min(iy + h + 1, ny)
    sub = zb[j0:j1, i0:i1]
    gx = x0 + (np.arange(i0, i1) + 0.5) * res
    gy = y0 + (np.arange(j0, j1) + 0.5) * res
    GX, GY = np.meshgrid(gx, gy)
    fin = np.isfinite(sub)
    if fin.sum() < 6:
        return None, None, res, zb.shape
    zhat, info = tie.ground_elevation_at(GX[fin], GY[fin], sub[fin],
                                         easting, northing, radius_m,
                                         surface_order=2, quantile=0.50)
    if not np.isfinite(zhat):
        return None, None, res, zb.shape
    return float(info["n"]), float(zhat), res, zb.shape


def geoid_shift_mm_at(site_dir, easting, northing):
    """The geoid conversion the pipeline added to gen1 at this coordinate, mm.

    ``references.geoid_difference`` returns ``(a, b, c)`` to ADD to gen1, as a plane about
    the bounds centroid -- the same centroid the pipeline uses.
    """
    from lidar_diff_icp import references
    d = json.loads((Path(site_dir) / "corrections.json").read_text())
    b = d["bounds"]
    a, bx, cy = references.geoid_difference(b, d["crs"])
    cx, cyy = (b[0] + b[2]) / 2.0, (b[1] + b[3]) / 2.0
    return (a + bx * (easting - cx) / 1000.0 + cy * (northing - cyy) / 1000.0) * 1000.0


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--easting", type=float, required=True)
    p.add_argument("--northing", type=float, required=True)
    p.add_argument("--tracks", required=True)
    p.add_argument("--psids", type=int, nargs="+", required=True)
    p.add_argument("--covers", nargs="+", required=True)
    p.add_argument("--tiles", nargs="+", required=True)
    p.add_argument("--res", type=float, required=True)
    p.add_argument("--grid-sites", nargs="+", required=True,
                   help="data/derived/<site> dirs holding z_after.npy and a dod")
    p.add_argument("--dod-name", nargs="+", required=True,
                   help="one per --grid-sites; which DoD recovers gen1 there")
    p.add_argument("--radii-m", type=float, nargs="+", required=True,
                   help="swept; a single radius would hide the sensitivity")
    a = p.parse_args(argv)

    ts = L.load_tracks(a.tracks)
    control = G.load_control()

    R = Run("How far is the surface our DoD runs on from the surface the 2008 control "
            "actually measured?")
    R.input(control.origin, role="gen1's own 2008 control; dnr_surface_z_m is the "
                                 "vendor's DELIVERED surface elevation at the mark")
    R.input(a.tracks, role="flight-line tracks, to find the site's own-line marks")
    for s, dn in zip(a.grid_sites, a.dod_name):
        R.input(str(Path(s) / "z_after.npy"), role=f"gen2 gridded surface at {Path(s).name}")
        R.input(str(Path(s) / f"{dn}.npy"),
                role=f"DoD at {Path(s).name}; gen1's grid is recovered as z_after - this")
    R.param("site", (a.easting, a.northing), src="andy")
    R.param("psids", tuple(a.psids), src="andy")
    R.param("covers", tuple(a.covers), src="andy")
    R.param("dod_name", tuple(a.dod_name), src="MINE",
            why="which DoD product recovers gen1 at each site; verified on "
                "elba_fulldensity against an independent gen1 grid to median +0.489 mm, "
                "NMAD 0.391 mm over 341,239 cells")
    R.param("radii_m", tuple(a.radii_m), src="MINE",
            why="radius of grid cells handed to tie.ground_elevation_at; SWEPT, because "
                "a single radius would hide the sensitivity. The predecessor window "
                "MEDIAN moved 72.3 mm over 5-55 m, which is why an order-2 surface fit "
                "replaced it")
    R.column("point_id", "control mark id")
    R.column("cover", "MnDNR land-cover class of the mark, unitless code: L1O open, "
                      "L2T tall weeds/crops, L3B brush, L4F forest, L5U urban")
    R.column("gap", "A = estimator gap (our reduction of the DELIVERED cloud vs the "
                    "vendor's Surface Z). B = processing gap (our gridded CSF/aligned/"
                    "geoid-shifted surface vs the vendor's Surface Z)")
    R.column("delivered_z_m", "the vendor's own dnr_surface_z_m at the mark")
    R.column("ours_z_m", "our surface AT the mark, m: for A the cloud reduction, for B "
                         "the same order-2 surface fit + residual median applied to the "
                         "grid cell centres within radius")
    R.column("geoid_mm", "geoid conversion applied to gen1 at this mark to carry it onto "
                         "gen2's GEOID18 frame, mm (references.geoid_difference); 0 for "
                         "gap A, which is entirely within GEOID03")
    R.column("bridge_raw_mm", "z_delivered - z_ours, mm, WITHOUT undoing the geoid "
                              "conversion; for gap B this differences two geoids and is "
                              "reported only so the conversion's size is visible")
    R.column("bridge_mm", "z_delivered - z_ours with both on NAVD88(GEOID03), mm; ADD to "
                          "a delivered-surface constant to carry it onto our surface")
    R.notes.append("A and B are DIFFERENT quantities. A is not a proxy for B: it changes "
                   "only the reduction, while B also changes the ground classification, "
                   "the swath alignment and the vertical datum.")
    R.notes.append("Gap B's raw value differences NAVD88(GEOID03) control against a "
                   "GEOID18-framed surface. The geoid column undoes that; read bridge_mm, "
                   "not bridge_raw_mm.")
    R.banner()

    rows = []

    # ---- A: estimator gap, at the site's own-line marks
    sc, sites, meas, skipped, est = S.estimate(
        ts, psids=a.psids, easting=a.easting, northing=a.northing, scope="pass",
        half_width_m=S.SEAM_HALF_SPACING_M, covers=a.covers, tile_dirs=a.tiles,
        res=a.res, control=control)
    A = []
    for m in S.marks_on_scope_psids(meas, a.psids):
        pub = m.site.mark.dnr_error_m
        if pub is None:
            continue
        br = m.tie_mm - pub * 1000.0
        A.append(br)
        rows.append([m.point_id, m.site.mark.cover_class, "A",
                     f"{m.site.mark.dnr_surface_z_m:.3f}",
                     f"{m.site.mark.checkpoint.elevation - m.tie_mm/1000.0:.3f}",
                     "0.0", f"{br:+.1f}", f"{br:+.1f}"])

    # ---- B: processing gap, wherever our grid actually covers a mark
    B = []
    for win in a.radii_m:
      for sdir, dn in zip(a.grid_sites, a.dod_name):
        for mk in control:
            if mk.dnr_surface_z_m is None:
                continue
            cell, wmed, res, shape = read_grid_at(sdir, dn, mk.checkpoint.easting,
                                                 mk.checkpoint.northing, radius_m=win)
            if wmed is None:
                continue
            raw = (mk.dnr_surface_z_m - wmed) * 1000.0
            g = geoid_shift_mm_at(sdir, mk.checkpoint.easting, mk.checkpoint.northing)
            br = raw + g          # undo the GEOID03 -> GEOID18 carry
            B.append((mk.aliases[0], mk.cover_class, raw, g, br, win))
            rows.append([f"{mk.aliases[0]} r{win:g}", mk.cover_class, "B",
                         f"{mk.dnr_surface_z_m:.3f}", f"{wmed:.3f}",
                         f"{g:+.1f}", f"{raw:+.1f}", f"{br:+.1f}"])
    R.table(["point_id", "cover", "gap", "delivered_z_m", "ours_z_m", "geoid_mm",
             "bridge_raw_mm", "bridge_mm"], rows)

    print()
    A = np.array(A)
    if A.size:
        print(f"  A  estimator gap : n={A.size}  mean {A.mean():+.2f}  median "
              f"{np.median(A):+.2f}  sd {A.std(ddof=1):.2f}  "
              f"SE {A.std(ddof=1)/np.sqrt(A.size):.2f} mm")
    print(f"  B  processing gap: n={len(B)}"
          + ("  -- too few to carry an uncertainty; listed above, not summarised"
             if len(B) < 3 else ""))
    import collections
    per = collections.defaultdict(list)
    for pid, cv, raw, g, br, win in B:
        per[(pid, cv)].append((win, br))
    for (pid, cv), v in per.items():
        v.sort()
        vals = np.array([b for _, b in v])
        print(f"       {pid:<24s} {cv}  " +
              "  ".join(f"r{w:g}:{b:+7.1f}" for w, b in v) +
              f"   |  spread {np.ptp(vals):.1f} mm")
    print()
    print("  The products need B. A is reported because it is measurable and because")
    print("  conflating the two would let an estimator difference stand in for a")
    print("  processing difference.")
    R.done(headline=f"A n={A.size} mean {A.mean():+.2f} mm; B on {len(per)} marks swept over {len(a.radii_m)} radii")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

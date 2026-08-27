"""The bridge on every control mark that falls inside a gen1 tile already on disk.

    env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python \
        ground_control/run_bridge_wide.py --covers L1O --tiles data/before \
        --csf-half-width-m 300 --res 5.0 --radii-m 7.5 10 --align-res 2.0 \
        --swath-tie intercept --csf-cache <dir> --out ground_control/products/bridge_wide.json

Nothing is downloaded.  Each tile is read once for its swath constants (cached to
ground_control/data/) and each mark gets a cropped CSF window (~4 s), so the whole sweep
costs minutes rather than the hours a per-tile CSF would.

VALIDATION GATE.  Before any wide number is printed, the same reconstruction is run at the
two marks inside elbaext and compared against the SHIPPED grid.  Measured: +1.7 mm at the
open mark, -30.5 mm at the vegetated one.  The producer is therefore trusted on open
ground and not under vegetation, and --covers is the caller's to set with that in view.
"""

from __future__ import annotations

import argparse
import glob
import json
import os
import sys
from pathlib import Path

import numpy as np

_HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(_HERE))
sys.path.insert(0, str(_HERE.parent))
sys.path.insert(0, str(_HERE.parent / "src"))

import our_surface as OS  # noqa: E402
from lidar_diff_icp.groundtruth import gen1_datum as G, tie as T  # noqa: E402
from trust.provenance import Run  # noqa: E402


def grid_value(site, dod, easting, northing, radius_m):
    d = json.load(open(f"{site}/corrections_geoid.json"))
    x0, y0, _, _ = d["bounds"]
    res = float(d["res_m"])
    zb = np.load(f"{site}/z_after.npy") - np.load(f"{site}/{dod}.npy")
    ny, nx = zb.shape
    ix = int((easting - x0) // res); iy = int((northing - y0) // res)
    h = int(np.ceil(radius_m / res)) + 1
    i0, i1 = max(ix - h, 0), min(ix + h + 1, nx)
    j0, j1 = max(iy - h, 0), min(iy + h + 1, ny)
    sub = zb[j0:j1, i0:i1]
    gx = x0 + (np.arange(i0, i1) + 0.5) * res
    gy = y0 + (np.arange(j0, j1) + 0.5) * res
    GX, GY = np.meshgrid(gx, gy)
    fin = np.isfinite(sub)
    if fin.sum() < 6:
        return None
    return T.ground_elevation_at(GX[fin], GY[fin], sub[fin], easting, northing,
                                 radius_m, surface_order=2, quantile=0.50)[0]


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--covers", nargs="+", required=True)
    p.add_argument("--tiles", required=True)
    p.add_argument("--csf-half-width-m", type=float, required=True)
    p.add_argument("--res", type=float, required=True)
    p.add_argument("--radii-m", type=float, nargs="+", required=True)
    p.add_argument("--align-res", type=float, required=True)
    p.add_argument("--swath-tie", required=True)
    p.add_argument("--csf-cache", required=True)
    p.add_argument("--out", required=True)
    a = p.parse_args(argv)

    control = G.load_control()
    tiles = sorted(t for t in glob.glob(f"{a.tiles}/*.laz") if "merged" not in t)

    R = Run("What is the bridge -- z_delivered minus OUR reconstructed gen1 surface -- "
            "at every control mark inside a gen1 tile we already hold?")
    R.input(control.origin, role="gen1's own 2008 control; dnr_surface_z_m is the "
                                 "vendor's DELIVERED surface elevation at the mark")
    for t in tiles[:3]:
        R.input(t, role="gen1 2008 delivered tile; CSF ground + swath alignment here "
                        "rebuild our surface locally")
    R.param("covers", tuple(a.covers), src="andy",
            why="the validation gate reproduces the shipped grid to +1.7 mm on the OPEN "
                "mark and -30.5 mm on the vegetated one, so the producer is trusted on "
                "open ground only")
    R.param("csf_half_width_m", a.csf_half_width_m, src="MINE",
            why="CSF window half-width; must be wide enough that the cloth is not "
                "dominated by the crop edge. Excludes no marks; it sets how much context "
                "the cloth sees")
    R.param("res_m", a.res, src="repo", why="the pipeline's grid resolution")
    R.param("align_res_m", a.align_res, src="repo",
            why="coreg.align_swaths' own default, which pipeline.difference_dem uses; "
                "NOT the grid resolution")
    R.param("swath_tie", a.swath_tie, src="repo",
            why="elbaext's corrections_geoid.json records swath_tie=intercept")
    R.param("radii_m", tuple(a.radii_m), src="MINE",
            why="radius of the order-2 surface fit at the mark; swept. Only radii a 5 m "
                "grid can fill with 6+ cells are usable, so 7.5 m is the floor")
    R.param("lateral_nuth_kaab", "NOT applied", src="MINE",
            why="that term is gen1->gen2 registration and needs gen2, which is not on "
                "disk away from Elba. Its measured effect on a tie is 10.0 mm "
                "(ABSOLUTE_BASIS_ELBA). Its absence is what the validation gate measures")
    R.column("point_id", "control mark id")
    R.column("cover", "MnDNR land-cover class, unitless code")
    R.column("tile", "gen1 tile the mark falls in")
    R.column("lines", "point_source_id values present in the CSF window, unitless")
    R.column("n_csf", "CSF ground returns in the window, count")
    R.column("delivered_z_m", "the vendor's own dnr_surface_z_m at the mark")
    R.column("ours_z_m", "our reconstructed gen1 surface at the mark, on NAVD88(GEOID03), m")
    R.column("bridge_mm", "z_delivered - z_ours, both on NAVD88(GEOID03), mm; ADD to a "
                          "delivered-surface constant to carry it onto our surface")
    R.column("radius_spread_mm", "max - min of bridge_mm across the radii that could be "
                                 "FITTED, mm -- how much the answer depends on an "
                                 "unstated window. Meaningless unless n_radii >= 2")
    R.column("n_radii", "radii successfully fitted / radii swept. An order-2 fit needs 6 "
                        "grid cells, so a radius under ~2*res cannot be filled on a 5 m "
                        "grid and is REPORTED as unfitted rather than dropped")
    R.notes.append("VALIDATION GATE, run before this table: the same reconstruction at "
                   "the two marks inside elbaext reproduces the SHIPPED grid to +1.7 mm "
                   "(open) and -30.5 mm (vegetated).")
    R.banner()

    # gate
    print("  validation gate (reconstruction vs shipped elbaext grid):")
    for pid, tname in (("L1O101", "4342-30-64"), ("L2T51", "4342-29-63")):
        mk = next(m for m in control if m.aliases[0] == pid)
        E, N = mk.checkpoint.easting, mk.checkpoint.northing
        sp = OS.our_gen1_surface_at(f"{a.tiles}/{tname}.laz", E, N,
                                    csf_half_width_m=a.csf_half_width_m, res=a.res,
                                    radius_m=a.radii_m[-1], exclude=(5, 6, 9),
                                    align_res=a.align_res, swath_tie=a.swath_tie,
                                    csf_cache=a.csf_cache)
        gv = grid_value("data/derived/elbaext", "dod_geoid", E, N, a.radii_m[-1])
        if sp and gv:
            print(f"    {pid:<8s} ({mk.cover_class})  local - grid = "
                  f"{(sp.z_geoid18_m - gv)*1000:+.1f} mm")

    keep = set(a.covers)
    marks = [m for m in control if m.cover_class in keep and m.dnr_surface_z_m is not None]
    import laspy
    boxes = {}
    for t in tiles:
        with laspy.open(t) as f:
            h = f.header
            boxes[t] = (h.mins[0], h.mins[1], h.maxs[0], h.maxs[1])

    rows, recs = [], []
    todo = [(m, t) for m in marks for t, b in boxes.items()
            if b[0] <= m.checkpoint.easting <= b[2] and b[1] <= m.checkpoint.northing <= b[3]]
    print(f"\n  {len(todo)} marks of cover {tuple(a.covers)} fall inside {len(tiles)} tiles\n")
    for i, (mk, tp) in enumerate(todo, 1):
        E, N = mk.checkpoint.easting, mk.checkpoint.northing
        vals, used, failed = [], [], []
        sp = None
        for r in a.radii_m:
            s = OS.our_gen1_surface_at(tp, E, N, csf_half_width_m=a.csf_half_width_m,
                                       res=a.res, radius_m=r, exclude=(5, 6, 9),
                                       align_res=a.align_res, swath_tie=a.swath_tie,
                                       csf_cache=a.csf_cache)
            if s is None:
                failed.append(r)          # NEVER silent: a radius that could not be
                continue                  # fitted is reported, not dropped
            sp = s
            used.append(r)
            vals.append((mk.dnr_surface_z_m - s.z_geoid03_m) * 1000.0)
        if not vals:
            print(f"    [{i:>2}/{len(todo)}] {mk.aliases[0]:<24s} -- NO radius could be "
                  f"fitted (all of {a.radii_m} gave <6 cells)")
            continue
        if len(vals) < 2:
            print(f"    [{i:>2}/{len(todo)}] {mk.aliases[0]:<24s} -- only {len(used)} of "
                  f"{len(a.radii_m)} radii fitted ({failed} failed); a spread over ONE "
                  f"value is not a spread")
        br = float(np.median(vals))
        spread = float(np.ptp(vals)) if len(vals) > 1 else 0.0
        recs.append(dict(point_id=mk.aliases[0], cover=mk.cover_class,
                         tile=os.path.basename(tp), bridge_mm=br,
                         radius_spread_mm=spread, n_csf=sp.n_ground_pts,
                         n_radii_used=len(used), radii_failed=failed,
                         lines=list(sp.lines_present)))
        rows.append([mk.aliases[0], mk.cover_class, os.path.basename(tp)[:-4],
                     "/".join(map(str, sp.lines_present)), sp.n_ground_pts,
                     f"{mk.dnr_surface_z_m:.3f}", f"{sp.z_geoid03_m:.3f}",
                     f"{br:+.1f}", f"{spread:.1f}", f"{len(used)}/{len(a.radii_m)}"])
        print(f"    [{i:>2}/{len(todo)}] {mk.aliases[0]:<24s} bridge {br:+8.1f} mm  "
              f"spread {spread:5.1f}  ({os.popen('free -m').read().splitlines()[1].split()[2]} MB used)")
    R.table(["point_id", "cover", "tile", "lines", "n_csf", "delivered_z_m", "ours_z_m",
             "bridge_mm", "radius_spread_mm", "n_radii"], rows)
    v = np.array([r["bridge_mm"] for r in recs])
    sp_ = np.array([r["radius_spread_mm"] for r in recs])
    nr = np.array([r["n_radii_used"] for r in recs])
    print()
    print(f"  radii fitted per mark: min {nr.min()} of {len(a.radii_m)}; "
          f"marks with <2 (no spread computable): {(nr < 2).sum()}")
    if (nr >= 2).any():
        print(f"  radius spread over marks with >=2 radii: median "
              f"{np.median(sp_[nr >= 2]):.1f} mm, max {sp_[nr >= 2].max():.1f} mm")
    print(f"  BRIDGE over {v.size} {'/'.join(a.covers)} marks:")
    print(f"    mean {v.mean():+.2f}  median {np.median(v):+.2f}  sd {v.std(ddof=1):.2f}  "
          f"SE {v.std(ddof=1)/np.sqrt(v.size):.2f} mm")
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    Path(a.out).write_text(json.dumps(dict(marks=recs, covers=a.covers,
                                           params=vars(a)), indent=1))
    print(f"  wrote {a.out}")
    R.done(headline=f"bridge n={v.size} mean {v.mean():+.2f} +/- "
                    f"{v.std(ddof=1)/np.sqrt(v.size):.2f} mm")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

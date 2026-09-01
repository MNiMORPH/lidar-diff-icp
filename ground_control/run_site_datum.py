"""Measure both epochs' datum constants at a site and emit the pipeline's absolute_datum.

    env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python ground_control/run_site_datum.py \
        --easting 578762.8 --northing 4884487.6 --site elbaext \
        --corrections data/derived/elbaext/corrections_geoid.json \
        --tracks ground_control/data/gen1_line_tracks.json \
        --psids 133 134 135 136 137 138 --covers L1O --collinear-sigma 3 \
        --tiles data/before --res 5.0 --gen2-surface ql1_laz \
        --bridge-mm -4.04 --bridge-source "ground_control/products/bridge_wide_L1O.json" \
        --max-lags-m 20000 40000 80000 160000 --n-lags 25 --n-pairs 800000 \
        --estimators dowd matheron --seed 0 --out ground_control/products/SITE_DATUM_elbaext.json

The output JSON's ``absolute_datum`` block is passed straight to
``pipeline.difference_dem(..., absolute_datum=...)``.
"""
from __future__ import annotations
import argparse, json, sys
from pathlib import Path
_HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(_HERE)); sys.path.insert(0, str(_HERE.parent))
sys.path.insert(0, str(_HERE.parent / "src"))
import site_datum as SD  # noqa: E402
from trust.provenance import Run  # noqa: E402


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    for a_ in ("--easting", "--northing", "--res", "--collinear-sigma", "--bridge-mm"):
        p.add_argument(a_, type=float, required=True)
    p.add_argument("--site", required=True)
    p.add_argument("--corrections", required=True)
    p.add_argument("--tracks", required=True)
    p.add_argument("--psids", type=int, nargs="+", required=True)
    p.add_argument("--covers", nargs="+", required=True)
    p.add_argument("--tiles", nargs="+", required=True)
    p.add_argument("--gen2-surface", required=True)
    p.add_argument("--bridge-source", required=True)
    p.add_argument("--max-lags-m", type=float, nargs="+", required=True)
    p.add_argument("--n-lags", type=int, required=True)
    p.add_argument("--n-pairs", type=int, required=True)
    p.add_argument("--estimators", nargs="+", required=True)
    p.add_argument("--seed", type=int, required=True)
    p.add_argument("--out", required=True)
    a = p.parse_args(argv)

    corr = json.loads(Path(a.corrections).read_text())
    gauge = corr.get("zero_line")
    if gauge is None:
        gauge = min(int(k) for k in corr["per_swath_internal_alignment_dxdydz_m"])

    R = Run("What datum constants place BOTH epochs on surveyed NAVD88 at this site, and "
            "what does that do to the DoD?")
    R.input(a.corrections, role="the site's pipeline corrections; supplies bounds, crs and "
                                "the swath gauge the gen1 constant must be tied to")
    R.input(a.tracks, role="gen1 flight-line tracks, for reused-psid disambiguation only")
    R.param("site", (a.easting, a.northing), src="andy")
    R.param("psids", tuple(a.psids), src="repo", why="the lines this tile is built from")
    R.param("covers", tuple(a.covers), src="andy",
            why="OPEN GROUND ONLY -- pooling covers bakes canopy response into the datum")
    R.param("gen2_surface", a.gen2_surface, src="andy",
            why="gen2 publishes four; the choice moved the Elba answer by 12.55 mm")
    R.param("collinear_sigma", a.collinear_sigma, src="andy")
    R.param("bridge_mm", a.bridge_mm, src="andy", why=f"supplied: {a.bridge_source}")
    R.param("zero_line", gauge, src="repo",
            why="read from the corrections file; the gen1 constant is TIED to it and the "
                "pipeline refuses a mismatch")
    R.column("quantity", "what is being reported")
    R.column("value_mm", "millimetres; positive = the surface reads LOW, so ADD it")
    R.banner()

    sd = SD.measure_site_datum(
        easting=a.easting, northing=a.northing, bounds=corr["bounds"], crs=corr["crs"],
        psids=a.psids, tile_dirs=a.tiles, tracks_path=a.tracks, covers=a.covers,
        collinear_sigma=a.collinear_sigma, res=a.res, gen2_surface=a.gen2_surface,
        bridge_mm=a.bridge_mm, bridge_source=a.bridge_source, gauge_ref=gauge,
        max_lags_m=a.max_lags_m, n_lags=a.n_lags, n_pairs=a.n_pairs,
        estimators=a.estimators, seed=a.seed)

    pipe = sd.to_pipeline(source=f"ground_control/run_site_datum.py @ {a.site}")
    R.table(["quantity", "value_mm"], [
        [f"gen1 delivered ({sd.gen1_n_marks} marks, {sd.gen1_n_lines} lines)",
         f"{sd.gen1_delivered_mm:+.2f} +/- {sd.gen1_delivered_se_mm:.2f}"],
        ["bridge, delivered -> ours", f"{sd.bridge_mm:+.2f}"],
        ["gen1 on OUR surface", f"{sd.gen1_our_surface_mm:+.2f}"],
        [f"gen2 delivered ({sd.gen2_surface})", f"{sd.gen2_delivered_mm:+.2f}"],
        ["geoid term added to gen1", f"{sd.geoid_mm:+.2f}"],
        ["gen1 constant IN the DoD frame", f"{pipe['gen1_mm']:+.2f}"],
        ["DoD shift (gen2_mm - gen1_mm)", f"{sd.dod_shift_mm:+.2f}"]])
    for w in sd.warnings:
        print(f"  WARNING: {w}")
    out = sd.to_dict(); out["absolute_datum"] = pipe
    Path(a.out).write_text(json.dumps(out, indent=1))
    print(f"\n  wrote {a.out}")
    print("  pass its 'absolute_datum' block to pipeline.difference_dem(absolute_datum=...)")
    R.done(headline=f"gen1 {sd.gen1_our_surface_mm:+.2f}, gen2 {sd.gen2_delivered_mm:+.2f}, "
                    f"DoD shift {sd.dod_shift_mm:+.2f} mm")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

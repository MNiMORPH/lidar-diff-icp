#!/usr/bin/env python3
"""Trace ridgelines on ANY tile with the Scherler & Schwanghart (2020) divide network,
via Andy's verified reimplementation (rivernetworkx.dreich.drainage_divides).

Divides are basin BOUNDARIES of the channel network (the dual of the drainage net), so
sub-channelization-threshold farm-furrow micro-topography cannot create spurious ridges --
the failure mode of the naive inverted-flow-accumulation approach (which traced tillage
furrows in the fields). QC each threshold by furrow contamination (open+flat field cells).

Generic over tiles: grid (origin, res) is read from the tile's corrections*/meta JSON, and
the DEM/slope/penetration are loaded from the tile dir. slope+penetration are OPTIONAL
(QC only); a tile with just z_after + a corrections json still traces.

    env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python analysis/ridgelines/trace_ridgelines.py [TILE_DIR]
        TILE_DIR default data/derived/elba_fulldensity; e.g. data/derived/elbaext
"""
import sys, os, json, argparse, numpy as np
sys.path.insert(0, "/home/awickert/dataanalysis/r.fluvial")
from rivernetworkx import dreich as D
from scipy.ndimage import uniform_filter, distance_transform_edt, gaussian_filter


def _fill_nan(z):
    zf = np.asarray(z, float).copy(); nm = ~np.isfinite(zf)
    if nm.any():
        zf = zf[tuple(distance_transform_edt(nm, return_distances=False, return_indices=True))]
    return zf


def trace_ridgelines(z, res, *, threshold=200, smooth_sigma_m=12.5, tpi_window_m=305.0,
                     slope=None, penetration=None, sweep=(50, 100, 200, 500, 1000)):
    """Ridge mask (bool, z.shape) from S&S divides on the gen2 bare earth.

    Returns (ridge, qc, sweep_stats). ``z`` NaNs are filled internally. Tracing is done on
    a Gaussian-smoothed DEM (kill sub-~15 m tillage furrows so D8 does not route along them);
    convexity/analysis downstream is measured on the ORIGINAL DEM. ``slope``+``penetration``
    (gen2) are optional and used only for the field/forest QC stats.
    """
    zf = _fill_nan(z)
    ztrace = gaussian_filter(zf.astype(np.float64), sigma=smooth_sigma_m / res)
    filled = D.fill(ztrace, nodata=-9999.0, cellsize=res)
    fi = D.build_flowinfo(filled, nodata=-9999.0, cellsize=res)
    tpi_large = zf - uniform_filter(zf, size=max(3, int(round(tpi_window_m / res))))
    field = None
    if slope is not None and penetration is not None:
        field = (np.asarray(penetration) >= 0.45) & (np.asarray(slope) < 5)   # open + flat = cultivated
    stats = {}
    for T in sorted(set(sweep) | {threshold}):
        divide, _ = D.drainage_divides(fi, threshold=T)
        s = {"cells": int(divide.sum())}
        if divide.sum():
            s["pct_on_highs"] = 100.0 * float(np.mean(tpi_large[divide] > 0))
            if field is not None:
                s["pct_field_furrow"] = 100.0 * float(np.mean(field[divide]))
                s["pct_forest"] = 100.0 * float(np.mean(np.asarray(penetration)[divide] < 0.25))
        stats[T] = (divide, s)
    ridge, qc = stats[threshold]
    return ridge, qc, {T: v[1] for T, v in stats.items()}


def _grid_from_tile(tile_dir):
    """(x0, y0, res) from the tile's corrections*/meta JSON (bounds + res)."""
    for fn in ("corrections_geoid.json", "corrections.json", "meta.json"):
        p = os.path.join(tile_dir, fn)
        if os.path.exists(p):
            j = json.load(open(p)); b = j.get("bounds"); r = j.get("res_m") or j.get("res")
            if b and r:
                return float(b[0]), float(b[1]), float(r)
    raise SystemExit(f"no bounds/res found in {tile_dir} (need a corrections*.json or meta.json)")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("tile_dir", nargs="?", default="data/derived/elba_fulldensity")
    ap.add_argument("--threshold", type=int, default=200)
    ap.add_argument("--out", default="ridge_mask.npy")
    ap.add_argument("--no-fig", action="store_true")
    a = ap.parse_args(); d = a.tile_dir.rstrip("/")
    x0, y0, res = _grid_from_tile(d)
    z = np.load(os.path.join(d, "z_after.npy"))
    slope = np.load(os.path.join(d, "slope.npy")) if os.path.exists(os.path.join(d, "slope.npy")) else None
    penp = os.path.join(d, "penetration.npy")
    pen = np.load(penp) if os.path.exists(penp) else None
    ridge, qc, sweep = trace_ridgelines(z, res, threshold=a.threshold, slope=slope, penetration=pen)

    print(f"{'threshold':>10} {'cells':>8} {'%highs':>7} {'%furrow':>8} {'%forest':>8}")
    for T, s in sorted(sweep.items()):
        f = lambda k: s.get(k, float('nan'))
        print(f"{T:>10} {s.get('cells', 0):>8} {f('pct_on_highs'):>6.0f}% "
              f"{f('pct_field_furrow'):>7.0f}% {f('pct_forest'):>7.0f}%")
    outp = os.path.join(d, a.out); np.save(outp, ridge)
    print(f"\nthreshold={a.threshold}: {int(ridge.sum())} ridge cells, "
          f"{qc.get('pct_on_highs', float('nan')):.0f}% on highs -> {outp}")
    if pen is None:
        print("  (no penetration.npy in tile -> field-furrow/forest QC skipped)")

    if not a.no_fig:
        import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
        from lidar_diff_icp.viz import hillshade
        zf = _fill_nan(z); ny, nx = z.shape
        hs = hillshade(zf, res, x0, y0, fill_gaps=True); ext = (x0, x0 + nx * res, y0, y0 + ny * res)
        fig, ax = plt.subplots(figsize=(11, 13))
        ax.imshow(hs, extent=ext, origin="lower", cmap="gray")
        ov = np.zeros((ny, nx, 4)); ov[ridge] = (0.9, 0.1, 0.1, 1.0)
        ax.imshow(ov, extent=ext, origin="lower")
        ax.set_title(f"S&S divide ridgelines (threshold={a.threshold}): {int(ridge.sum())} cells, "
                     f"{qc.get('pct_on_highs', float('nan')):.0f}% on highs")
        ax.set_xlabel("Easting (m)"); ax.set_ylabel("Northing (m)")
        os.makedirs("figures", exist_ok=True)
        figp = f"figures/ridgelines_{os.path.basename(d)}.png"
        fig.savefig(figp, dpi=130, bbox_inches="tight"); plt.close(fig)
        print("wrote", figp)


if __name__ == "__main__":
    main()

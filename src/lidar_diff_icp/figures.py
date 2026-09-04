"""The two standard per-site figures, rebuildable from the saved products.

These lived inside ``scripts/run_all_sites.py`` and took in-memory arrays, so regenerating
either one meant re-running the whole site pipeline -- minutes of point-cloud work to redraw
a picture whose every input was already on disk. They read a tile directory now, so a figure
is reproducible on its own and can be re-made after a colour or label change without
touching the data.

    A. ``<name>_dod_lod.png``  the DoD beside the level of detection
    B. ``<name>_change.png``   change above the detection limit, 70% opaque over hillshade

    lidar-diff-figures --tile data/derived/carlton
"""
from __future__ import annotations

import argparse
import json
import os

import numpy as np

DEFAULT_FIGDIR = "figures/sites"


def grid_of(tile_dir):
    """(X0, Y0, res, nx, ny) from the tile's own metadata, not from an argument."""
    for fn in ("meta.json", "corrections_geoid.json", "corrections.json"):
        p = os.path.join(tile_dir, fn)
        if os.path.exists(p):
            j = json.load(open(p))
            b = j["bounds"]
            res = float(j.get("res") or j.get("res_m"))
            return b[0], b[1], res, int(round((b[2] - b[0]) / res)), int(round((b[3] - b[1]) / res))
    raise SystemExit(f"no grid metadata in {tile_dir}")


def _need(tile_dir, *names):
    out = []
    for n in names:
        p = os.path.join(tile_dir, n)
        if not os.path.exists(p):
            raise SystemExit(f"{p} is missing; this figure cannot be drawn without it")
        out.append(np.load(p))
    return out


def dod_lod_figure(tile_dir, figdir=DEFAULT_FIGDIR, name=None,
                   dod_name="dod.npy", lod_name="lod.npy", suffix=""):
    """A: the DoD raster beside the level of detection.

    ``dod_name``/``lod_name`` choose WHICH DoD is drawn; the defaults are the base products,
    so every existing caller is unaffected. Pass ``dod_cover_q2.npy`` / ``lod_cover_q2.npy``
    with a ``suffix`` to draw the vegetation-corrected pair instead -- which had no figure
    producer at all, so the corrected DoD could only be seen by loading the array.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from lidar_diff_icp.viz import hillshade

    name = name or os.path.basename(str(tile_dir).rstrip("/"))
    X0, Y0, res, nx, ny = grid_of(tile_dir)
    Z21, dod, lod = _need(tile_dir, "z_after.npy", dod_name, lod_name)

    hs = hillshade(Z21, res, X0, Y0, fill_gaps=False)  # nodata -> white, as in the LoD panel
    ext = (X0, X0 + nx * res, Y0, Y0 + ny * res); v = 0.3
    fig, ax = plt.subplots(1, 2, figsize=(15, 9))
    ax[0].imshow(hs, extent=ext, origin="lower", cmap="gray", alpha=0.6)
    im0 = ax[0].imshow(dod, extent=ext, origin="lower", cmap="RdBu", vmin=-v, vmax=v)
    ax[0].set_title(f"{name}: DEM of Difference (gridded ground): gen2 - gen1 (m)\n"
                    f"red = erosion, blue = deposition   [{dod_name}]")
    fig.colorbar(im0, ax=ax[0], shrink=0.6, extend="both")
    im1 = ax[1].imshow(lod, extent=ext, origin="lower", cmap="viridis", vmin=0, vmax=0.2)
    ax[1].set_title(f"level of detection (m)   [{lod_name}]")
    fig.colorbar(im1, ax=ax[1], shrink=0.6, extend="max")
    for a in ax:
        a.set_xlabel("Easting (m)"); a.set_ylabel("Northing (m)")
    os.makedirs(figdir, exist_ok=True)
    out = f"{figdir}/{name}_dod_lod{suffix}.png"
    fig.savefig(out, dpi=130, bbox_inches="tight"); plt.close(fig)
    return out


def change_figure(tile_dir, figdir=DEFAULT_FIGDIR, name=None,
                  dod_name="dod.npy", change_name="change.npy", suffix=""):
    """B: the hillshade with robustly-detected DoD cells at 70% opacity."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from lidar_diff_icp.viz import hillshade

    name = name or os.path.basename(str(tile_dir).rstrip("/"))
    X0, Y0, res, nx, ny = grid_of(tile_dir)
    Z21, dod, change = _need(tile_dir, "z_after.npy", dod_name, change_name)

    hs = hillshade(Z21, res, X0, Y0, fill_gaps=True)  # gap-filled backdrop, no white holes
    ext = (X0, X0 + nx * res, Y0, Y0 + ny * res); v = 0.3
    over = np.where(change, dod, np.nan)
    fig, ax = plt.subplots(figsize=(10, 11))
    ax.imshow(hs, extent=ext, origin="lower", cmap="gray")
    im = ax.imshow(over, extent=ext, origin="lower", cmap="RdBu", vmin=-v, vmax=v, alpha=0.7)
    ax.set_title(f"{name}: topographic change above detection limits")
    ax.set_xlabel("Easting (m)"); ax.set_ylabel("Northing (m)")
    fig.colorbar(im, ax=ax, shrink=0.6, extend="both", label="detected DoD (m)")
    os.makedirs(figdir, exist_ok=True)
    out = f"{figdir}/{name}_change{suffix}.png"
    fig.savefig(out, dpi=130, bbox_inches="tight"); plt.close(fig)
    return out


FIGURES = {"dod_lod": dod_lod_figure, "change": change_figure}


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--tile", required=True, help="tile directory under data/derived/")
    ap.add_argument("--figdir", default=DEFAULT_FIGDIR)
    ap.add_argument("--which", nargs="*", choices=sorted(FIGURES), default=sorted(FIGURES))
    ap.add_argument("--name", default=None, help="title/filename stem (default: tile name)")
    a = ap.parse_args(argv)
    for k in a.which:
        print(f"wrote {FIGURES[k](a.tile, a.figdir, a.name)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

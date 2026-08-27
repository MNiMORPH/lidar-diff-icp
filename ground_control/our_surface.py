"""Reconstruct OUR gen1 surface at a point, anywhere a gen1 tile is on disk.

Why
---
The bridge -- ``z_delivered - z_ours`` -- can only be measured where our reconstructed
surface exists, and it existed on two tiles holding two control marks.  This rebuilds the
same surface locally, around any mark, so the bridge can be measured on the 89 marks that
fall inside the 46 gen1 tiles already on disk.

The recipe, matching ``pipeline.difference_dem``'s gen1 side
------------------------------------------------------------
1. **CSF ground** in a window about the mark (``ground_source="csf"``), via
   ``groundtruth.tie.csf_ground_near`` -> ``ground.classify_ground_csf`` at repo defaults.
   Cropping first is what makes this affordable: ~4 s on a 600 m box against minutes for
   a whole tile.
2. **Swath alignment**: the tile's own ``coreg.align_swaths`` constants, added per
   ``point_source_id``.  Computed once per tile and cached, because a tile read is ~400 MB
   and this is a shared laptop.
3. **Geoid**: ``references.geoid_difference`` over the window, added to gen1 to carry it
   onto gen2's NAVD88(GEOID18) frame -- then UNDONE again when comparing against the
   GEOID03 control.  It is applied and removed explicitly rather than skipped, so the
   recipe stays the product's recipe.
4. **Grid** at the pipeline's ``res`` with ``ground="slope_normal"``, ``ground_q=0.50``,
   then read AT the mark with ``tie.ground_elevation_at`` -- the same way the elbaext grid
   is read, so the two are compared like with like.

What is deliberately NOT included
---------------------------------
The **Nuth-Kaeaeb lateral (x, y) shift**.  That term is gen1->gen2 registration and needs
gen2, which is not on disk away from Elba.  ``ABSOLUTE_BASIS_ELBA.md`` measures its effect
on a tie at 10.0 mm.  Its absence is the known difference between this reconstruction and
the shipped grid, and :func:`validate_against_grid` measures it rather than assuming it.
"""

from __future__ import annotations

import json
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np

_HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(_HERE))
sys.path.insert(0, str(_HERE.parent / "src"))

from lidar_diff_icp import coreg, io, references  # noqa: E402
from lidar_diff_icp.groundtruth import tie as T  # noqa: E402

CACHE = _HERE / "data" / "swath_constants_cache.json"


@dataclass(frozen=True)
class SurfacePoint:
    """Our gen1 surface at one coordinate, and everything that went into it."""

    z_geoid18_m: float          # on gen2's frame, as the product carries it
    z_geoid03_m: float          # geoid term removed, comparable to the 2008 control
    geoid_mm: float
    n_ground_pts: int
    n_cells: int
    radius_m: float
    csf_half_width_m: float
    swath_source: str
    lines_present: tuple
    note: str = ""


def swath_constants(tile_path, *, res, exclude, swath_tie, cache_path=CACHE):
    """``{line: (dx, dy, dz)}`` from ``coreg.align_swaths`` on the tile, cached on disk.

    Matches ``pipeline.difference_dem``: ``align_swaths(pc, ref=int(ps.min()),
    tie=swath_tie)`` on the tile's cloud, with ``res`` at align_swaths' own default of
    2.0 -- NOT the 5 m grid resolution, which is a different parameter.

    Read once per tile and cached: a gen1 tile is ~7 M points and ~400 MB in memory.
    """
    key = str(Path(tile_path).name)
    cache = {}
    p = Path(cache_path)
    if p.exists():
        cache = json.loads(p.read_text())
    ck = f"{key}|res={res}|tie={swath_tie}|exclude={','.join(map(str, exclude))}"
    if ck in cache:
        return {int(k): tuple(v) for k, v in cache[ck].items()}
    pc = io.read_tile(tile_path)
    corr, _edges, _mis = coreg.align_swaths(
        pc, res=res, exclude=exclude, ref=int(np.min(pc.point_source_id)),
        tie=swath_tie)
    del pc
    out = {int(k): tuple(float(z) for z in v) for k, v in corr.items()}
    cache[ck] = {str(k): list(v) for k, v in out.items()}
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(cache, indent=1))
    return out


def our_gen1_surface_at(tile_path, easting, northing, *, csf_half_width_m, res,
                        radius_m, exclude, align_res, swath_tie, pdal=None,
                        csf_cache=None, apply_swath=True, geoid_bounds=None,
                        crs="EPSG:26915"):
    """Our gen1 surface at ``(easting, northing)``, both frames, from the tile on disk."""
    g = T.csf_ground_near(tile_path, easting, northing, csf_half_width_m,
                          pdal=pdal, cache_dir=csf_cache)
    x, y, z = np.asarray(g.x), np.asarray(g.y), np.asarray(g.z)
    ps = np.asarray(g.point_source_id)
    if x.size == 0:
        return None
    src = "none"
    if apply_swath:
        sc = swath_constants(tile_path, res=align_res, exclude=exclude,
                             swath_tie=swath_tie)
        src = f"align_swaths(tile, res={align_res}, tie={swath_tie}, ref=min)"
        for ln, (dx, dy, dz) in sc.items():
            m = ps == ln
            if m.any():
                x[m] += dx; y[m] += dy; z[m] += dz

    b = geoid_bounds or (easting - csf_half_width_m, northing - csf_half_width_m,
                         easting + csf_half_width_m, northing + csf_half_width_m)
    a0, bx, cy = references.geoid_difference(list(b), crs)
    cx, cyy = (b[0] + b[2]) / 2.0, (b[1] + b[3]) / 2.0
    gshift = a0 + bx * (easting - cx) / 1000.0 + cy * (northing - cyy) / 1000.0

    # grid at the pipeline resolution, slope-normal median, then read AT the mark
    x0 = easting - csf_half_width_m; y0 = northing - csf_half_width_m
    nx = int(np.ceil(2 * csf_half_width_m / res)); ny = nx
    ix = ((x - x0) / res).astype(int); iy = ((y - y0) / res).astype(int)
    ok = (ix >= 0) & (ix < nx) & (iy >= 0) & (iy < ny)
    if ok.sum() < 6:
        return None
    fid = iy[ok] * nx + ix[ok]
    import pandas as pd
    s = pd.Series(z[ok]).groupby(fid).quantile(0.50)
    cx_ = x0 + (s.index.values % nx + 0.5) * res
    cy_ = y0 + (s.index.values // nx + 0.5) * res
    zc = s.values
    zhat, info = T.ground_elevation_at(cx_, cy_, zc, easting, northing, radius_m,
                                       surface_order=2, quantile=0.50)
    if not np.isfinite(zhat):
        return None
    return SurfacePoint(z_geoid18_m=float(zhat + gshift), z_geoid03_m=float(zhat),
                        geoid_mm=float(gshift * 1000.0), n_ground_pts=int(x.size),
                        n_cells=int(zc.size), radius_m=float(radius_m),
                        csf_half_width_m=float(csf_half_width_m), swath_source=src,
                        lines_present=tuple(sorted(int(v) for v in np.unique(ps))))

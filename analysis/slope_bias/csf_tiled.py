#!/usr/bin/env python3
"""Pre-classify a gen1 cloud's ground with CSF, TILED, to fit in RAM.

The elbaext gen1 merge (17.35M pts over ~4.6x4.35 km) OOM-killed a single PDAL CSF
(machine had ~18 GB free; CSF materialises an uncompressed LAS copy plus its own
in-memory cloud + cloth). CSF is spatially LOCAL -- the cloth relaxes on a 1 m grid
and edge effects reach only a few cells -- so classifying in spatial tiles with a
generous overlap buffer and keeping each tile's CORE gives the same ground as one
run, at a fraction of peak RAM.

Output = a class-2-ground-only LAS (all attributes preserved) under
``data/csf_cache/``, which difference_dem loads and SKIPS CSF for (pipeline.py: if
csf_cache exists -> reuse). So elbaext_regrid.py then runs identically to the
fulldensity recipe.

The CSF parameters themselves are NOT set here: ``ground.classify_ground_csf`` is
called at its repo defaults (rigidness 1, everything else PDAL's own), so every cloud
classified through this script gets the identical ground filter. Only the TILING is a
parameter, and tiling is a memory device: the halo makes it a partition of the same
computation, not a different one.

Defaults reproduce the elbaext build exactly (3x3, 150 m halo). Any other cloud is a
matter of ``--src``/``--out``; ``--nx``/``--ny`` size the memory peak only.

    env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python analysis/slope_bias/csf_tiled.py
    env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python analysis/slope_bias/csf_tiled.py \
        --src data/before/4342-28-64.laz --out data/csf_cache/4342-28-64.las --nx 2 --ny 2
"""
import argparse
import os
import shutil

import numpy as np
import laspy

from lidar_diff_icp.ground import classify_ground_csf


def copy_all_dims(src, out_hdr, mask):
    """Copy EVERY dimension of ``src`` (masked) into a new LasData on ``out_hdr``.

    Never a hand-picked field list: that is what silently dropped scan_angle_rank
    from the first elbaext build.  x/y/z go through the scaled accessors so a tile
    whose header carries different scales/offsets still lands in the right place.
    """
    out = laspy.LasData(out_hdr)
    out.x = np.asarray(src.x)[mask]
    out.y = np.asarray(src.y)[mask]
    out.z = np.asarray(src.z)[mask]
    for d in src.point_format.dimension_names:
        if d in ("X", "Y", "Z"):
            continue
        setattr(out, d, np.asarray(getattr(src, d))[mask])
    return out


def core_mask(gx, gy, i, j, nx, ny, cx0, cx1, cy0, cy1):
    """Which of a halo tile's points belong to tile (i, j)'s CORE.

    Half-open [cx0, cx1) x [cy0, cy1), widened to swallow the outer domain edges so
    the NX*NY cores partition the domain exactly once -- no seam gap, no double count.
    """
    core = (gx >= cx0) & (gx < cx1) & (gy >= cy0) & (gy < cy1)
    if i == 0:
        core |= (gx < cx0)
    if i == nx - 1:
        core |= (gx >= cx1)
    if j == 0:
        core |= (gy < cy0)
    if j == ny - 1:
        core |= (gy >= cy1)
    return core


def run(src, out, nx, ny, overlap, tmpdir="data/derived/_csf_tiles_tmp"):
    os.makedirs(os.path.dirname(out) or ".", exist_ok=True)
    os.makedirs(tmpdir, exist_ok=True)

    f = laspy.read(src)
    x = np.asarray(f.x)
    y = np.asarray(f.y)
    X0, X1 = x.min(), x.max()
    Y0, Y1 = y.min(), y.max()
    dx = (X1 - X0) / nx
    dy = (Y1 - Y0) / ny
    print(f"src N={len(x):,}  extent E{X0:.0f}-{X1:.0f} N{Y0:.0f}-{Y1:.0f}", flush=True)

    hdr = laspy.LasHeader(point_format=f.header.point_format.id,
                          version=str(f.header.version))
    hdr.scales = f.header.scales
    hdr.offsets = f.header.offsets

    # The cache writer's header is taken from PDAL's OUTPUT (opened on the first tile),
    # not from the source: PDAL's writers.las promotes LAS 1.1/PF1 -> 1.4/PF7, which
    # RENAMES scan_angle_rank (int deg) -> scan_angle (int16, 0.006 deg).  Re-imposing
    # the PF1 source header here is what silently zeroed the elbaext scan angles.
    # gen1_save_angles_slope.py reads either name, and the elba cache is PF7 already,
    # so keeping PDAL's format makes both tiles identical in method.
    total = 0
    w = None
    ghdr = None
    try:
        for j in range(ny):
            for i in range(nx):
                cx0, cx1 = X0 + i * dx, X0 + (i + 1) * dx      # tile CORE (no halo)
                cy0, cy1 = Y0 + j * dy, Y0 + (j + 1) * dy
                # read points in core + halo, CSF them, then KEEP only core points
                hm = ((x >= cx0 - overlap) & (x <= cx1 + overlap) &
                      (y >= cy0 - overlap) & (y <= cy1 + overlap))
                tile_in = os.path.join(tmpdir, f"tile_{i}_{j}.laz")
                copy_all_dims(f, hdr, hm).write(tile_in)
                print(f"tile[{i},{j}] halo N={hm.sum():,} -> CSF ...", flush=True)
                gpath = classify_ground_csf(tile_in)         # class-2 ground LAS
                g = laspy.read(gpath)
                core = core_mask(np.asarray(g.x), np.asarray(g.y),
                                 i, j, nx, ny, cx0, cx1, cy0, cy1)
                if w is None:                                # first tile fixes the format
                    ghdr = laspy.LasHeader(point_format=g.header.point_format.id,
                                           version=str(g.header.version))
                    ghdr.scales = g.header.scales
                    ghdr.offsets = g.header.offsets
                    w = laspy.open(out, mode="w", header=ghdr)
                    print(f"cache format: PF{ghdr.point_format.id} v{ghdr.version} "
                          f"dims={len(list(ghdr.point_format.dimension_names))}", flush=True)
                w.write_points(copy_all_dims(g, ghdr, core).points)
                total += int(core.sum())
                print(f"tile[{i},{j}] ground-in-core N={core.sum():,}  (cum {total:,})",
                      flush=True)
                os.remove(tile_in)
                shutil.rmtree(os.path.dirname(gpath), ignore_errors=True)
    finally:
        if w is not None:
            w.close()

    # fail loudly rather than leave a silently angle-less cache behind again
    chk = laspy.read(out)
    sa = (np.asarray(chk.scan_angle) * 0.006
          if "scan_angle" in chk.point_format.dimension_names
          else np.asarray(chk.scan_angle_rank).astype(float))
    nz = 100 * (sa != 0).mean()
    print(f"wrote {out}  N={total:,} class-2 ground pts  "
          f"scan angle: nonzero%={nz:.1f} range[{sa.min():.1f},{sa.max():.1f}] deg", flush=True)
    assert nz > 50, "scan angle lost in the CSF round-trip -- cache is unusable, do NOT ship it"
    return total


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--src", default="data/before/elbaext_gen1_merged.laz")
    ap.add_argument("--out", default="data/csf_cache/elbaext.las")
    ap.add_argument("--nx", type=int, default=3,
                    help="tile columns. 3x3 -> ~2.5M pts/tile incl. halo on elbaext. 2x2 "
                         "(~5.4M) ran to completion but PDAL CSF costs ~1 kB/point, and its "
                         "peak on top of a loaded desktop drove this shared laptop into swap.")
    ap.add_argument("--ny", type=int, default=3, help="tile rows")
    ap.add_argument("--overlap", type=float, default=150.0,
                    help="m halo; CSF edge reach is a few cloth cells, 150 m is ample")
    A = ap.parse_args()
    run(A.src, A.out, A.nx, A.ny, A.overlap)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Pre-classify the elbaext gen1 ground with CSF, TILED, to fit in RAM.

The elbaext gen1 merge (17.35M pts over ~4.6x4.35 km) OOM-killed a single PDAL CSF
(machine had ~18 GB free; CSF materialises an uncompressed LAS copy plus its own
in-memory cloud + cloth). CSF is spatially LOCAL -- the cloth relaxes on a 1 m grid
and edge effects reach only a few cells -- so classifying in spatial tiles with a
generous overlap buffer and keeping each tile's CORE gives the same ground as one
run, at a fraction of peak RAM.

Output = data/csf_cache/elbaext.las (class-2 ground only, all attributes preserved),
which difference_dem loads and SKIPS CSF for (pipeline.py: if csf_cache exists ->
reuse). So elbaext_regrid.py then runs identically to the fulldensity recipe.

    env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python analysis/slope_bias/elbaext_csf_tiled.py
"""
import os, numpy as np, laspy
from lidar_diff_icp.ground import classify_ground_csf

SRC   = "data/before/elbaext_gen1_merged.laz"
CACHE = "data/csf_cache/elbaext.las"
NX, NY = 2, 2          # 2x2 tiles -> ~4.3M pts/tile, well within RAM
OVERLAP = 150.0        # m halo; CSF edge reach is a few cloth cells, 150 m is ample

os.makedirs("data/csf_cache", exist_ok=True)
os.makedirs("data/derived/_csf_tiles_tmp", exist_ok=True)

f = laspy.read(SRC)
x = np.asarray(f.x); y = np.asarray(f.y)
X0, X1 = x.min(), x.max(); Y0, Y1 = y.min(), y.max()
dx = (X1 - X0) / NX; dy = (Y1 - Y0) / NY
print(f"src N={len(x):,}  extent E{X0:.0f}-{X1:.0f} N{Y0:.0f}-{Y1:.0f}", flush=True)

hdr = laspy.LasHeader(point_format=f.header.point_format.id, version=str(f.header.version))
hdr.scales = f.header.scales; hdr.offsets = f.header.offsets


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


# The cache writer's header is taken from PDAL's OUTPUT (opened on the first tile),
# not from the source: PDAL's writers.las promotes LAS 1.1/PF1 -> 1.4/PF7, which
# RENAMES scan_angle_rank (int deg) -> scan_angle (int16, 0.006 deg).  Re-imposing
# the PF1 source header here is what silently zeroed the elbaext scan angles.
# gen1_save_angles_slope.py reads either name, and the elba cache is PF7 already,
# so keeping PDAL's format makes both tiles identical in method.
total = 0
w = None
try:
    for j in range(NY):
        for i in range(NX):
            cx0, cx1 = X0 + i*dx, X0 + (i+1)*dx      # tile CORE (no halo)
            cy0, cy1 = Y0 + j*dy, Y0 + (j+1)*dy
            # read points in core + halo, CSF them, then KEEP only core points
            hm = ((x >= cx0-OVERLAP) & (x <= cx1+OVERLAP) &
                  (y >= cy0-OVERLAP) & (y <= cy1+OVERLAP))
            tile_in = f"data/derived/_csf_tiles_tmp/tile_{i}_{j}.laz"
            copy_all_dims(f, hdr, hm).write(tile_in)
            print(f"tile[{i},{j}] halo N={hm.sum():,} -> CSF ...", flush=True)
            gpath = classify_ground_csf(tile_in)         # class-2 ground LAS
            g = laspy.read(gpath)
            gx = np.asarray(g.x); gy = np.asarray(g.y)
            core = (gx >= cx0) & (gx < cx1) & (gy >= cy0) & (gy < cy1)
            # include the outer domain edges fully (>= / <) so no seam gap at X0/Y0/X1/Y1
            if i == 0:      core |= (gx < cx0)
            if i == NX-1:   core |= (gx >= cx1)
            if j == 0:      core |= (gy < cy0)
            if j == NY-1:   core |= (gy >= cy1)
            if w is None:                                # first tile fixes the output format
                ghdr = laspy.LasHeader(point_format=g.header.point_format.id,
                                       version=str(g.header.version))
                ghdr.scales = g.header.scales; ghdr.offsets = g.header.offsets
                w = laspy.open(CACHE, mode="w", header=ghdr)
                print(f"cache format: PF{ghdr.point_format.id} v{ghdr.version} "
                      f"dims={len(list(ghdr.point_format.dimension_names))}", flush=True)
            w.write_points(copy_all_dims(g, ghdr, core).points)
            total += int(core.sum())
            print(f"tile[{i},{j}] ground-in-core N={core.sum():,}  (cum {total:,})", flush=True)
            os.remove(tile_in)
            import shutil; shutil.rmtree(os.path.dirname(gpath), ignore_errors=True)
finally:
    if w is not None:
        w.close()

# fail loudly rather than leave a silently angle-less cache behind again
chk = laspy.read(CACHE)
sa = (np.asarray(chk.scan_angle) * 0.006 if "scan_angle" in chk.point_format.dimension_names
      else np.asarray(chk.scan_angle_rank).astype(float))
nz = 100 * (sa != 0).mean()
print(f"wrote {CACHE}  N={total:,} class-2 ground pts  "
      f"scan angle: nonzero%={nz:.1f} range[{sa.min():.1f},{sa.max():.1f}] deg", flush=True)
assert nz > 50, "scan angle lost in the CSF round-trip -- cache is unusable, do NOT ship it"

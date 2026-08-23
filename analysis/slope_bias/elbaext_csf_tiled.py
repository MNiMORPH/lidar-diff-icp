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

total = 0
with laspy.open(CACHE, mode="w", header=hdr) as w:
    for j in range(NY):
        for i in range(NX):
            cx0, cx1 = X0 + i*dx, X0 + (i+1)*dx      # tile CORE (no halo)
            cy0, cy1 = Y0 + j*dy, Y0 + (j+1)*dy
            # read points in core + halo, CSF them, then KEEP only core points
            hm = ((x >= cx0-OVERLAP) & (x <= cx1+OVERLAP) &
                  (y >= cy0-OVERLAP) & (y <= cy1+OVERLAP))
            tile_in = f"data/derived/_csf_tiles_tmp/tile_{i}_{j}.laz"
            sub = laspy.LasData(hdr)
            for d in ("x","y","z","intensity","return_number","number_of_returns",
                      "classification","point_source_id","gps_time","scan_angle_rank"):
                try: setattr(sub, d, np.asarray(getattr(f, d))[hm])
                except Exception: pass
            sub.write(tile_in)
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
            gg = laspy.LasData(hdr)
            for d in ("x","y","z","intensity","return_number","number_of_returns",
                      "classification","point_source_id","gps_time","scan_angle_rank"):
                try: setattr(gg, d, np.asarray(getattr(g, d))[core])
                except Exception: pass
            w.write_points(gg.points)
            total += int(core.sum())
            print(f"tile[{i},{j}] ground-in-core N={core.sum():,}  (cum {total:,})", flush=True)
            os.remove(tile_in)
            import shutil; shutil.rmtree(os.path.dirname(gpath), ignore_errors=True)

print(f"wrote {CACHE}  N={total:,} class-2 ground pts", flush=True)

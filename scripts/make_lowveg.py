#!/usr/bin/env python3
"""Per-cell ``lowveg`` for a tile: the near-ground vegetation COMPOSITION.

    lowveg = (returns with lo < h <= hi) / (returns in the whole near-ground window)

This is the metric of ``analysis/CONTROL_LOWVEG_OFFSET.md``, which fitted it against
SURVEYED CONTROL and found it outperforms canopy cover -- "It is a composition, not a
canopy density, which is why it outperforms canopy cover here." Nothing here is invented:
the band edges, the population and the bin-centre test are all read off that document and
``analysis/control_lowveg_offset.py``, so a coefficient fitted at the control marks and a
value computed here refer to the same quantity.

Three things make the tile cube and the control-mark boxes commensurate, and all three
were checked rather than assumed:

  window      ``nearground_cells_sn.npz`` stores -1.00..+2.00 m at 0.02 m; the control
              boxes store 150 bins over -1.00..+2.00 m. Same window, same bin width.
  population  ``H2`` is ALL gen2 returns -- nearground_cells.py: "no classification, no
              ground/vegetation decision, no height threshold beyond the window" -- which
              is the ``ng_all`` array ``control_lowveg_offset.lowveg()`` reads. The
              class-2 split lives in nearground_gen2_class_split.npz and is NOT used here.
  bin-centre  the band test is on bin CENTRES, so with 20 mm bins the (0.14, 0.16] bin is
              counted whole and the metric is really "fraction above 0.14 m". The document
              flags this explicitly and states that the -290 mm/unit coefficient is correct
              FOR THE BINNED METRIC and must not be applied to an exactly-computed one.
              Reproducing the binning is therefore required, not incidental.

TWO COVERAGES, TWO FILES, because they are not interchangeable and a single name would let
one be used where the other belongs:

  lowveg.npy       from the near-ground cube. Only the divide cells nearground_cells.py
                   selected; NaN elsewhere. Enough to FIT q2(lowveg), which is evaluated on
                   exactly those cells.
  lowveg_grid.npy  --gen2: one streaming pass over the cloud, every cell in the tile. What
                   an APPLIED correction needs. Same definition, same bin-centre test; only
                   the set of cells differs. Nothing is interpolated into cells the cloud
                   does not cover -- they stay NaN.

dod_cover_corrected.py builds the grid-wide layer itself while it streams (it needs the
histogram anyway, so there is no pass to save). This producer exists so the layer can be
inspected, mapped and tracked for staleness WITHOUT running a correction to see it.

    ./lidar-icp/bin/python scripts/make_lowveg.py --tile elba
    ./lidar-icp/bin/python scripts/make_lowveg.py --tile elba --check
    env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python scripts/make_lowveg.py \
        --tile elba --gen2 data/after/3dep2021_fulldensity.laz
"""
import argparse
import os

import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("--tile", required=True, help="tile name under data/derived/, or a path")
ap.add_argument("--cube", default="nearground_cells_sn.npz")
# Band edges are CITED, not chosen here. CONTROL_LOWVEG_OFFSET.md, "lowveg, EXACT
# DEFINITION -- it must travel with any coefficient fitted from it": (0.15, 2.00] m.
ap.add_argument("--lo", type=float, default=0.15, help="band lower edge, m (doc: 0.15)")
ap.add_argument("--hi", type=float, default=2.00, help="band upper edge, m (doc: 2.00)")
ap.add_argument("--gen2", default=None,
                help="gen2 cloud. With it, lowveg is built for EVERY cell by streaming, and "
                     "written to lowveg_grid.npy. Without it, only the cube's cells are "
                     "covered and the output is lowveg.npy.")
ap.add_argument("--chunk", type=int, default=3_000_000)
ap.add_argument("--out", default=None)
ap.add_argument("--check", action="store_true",
                help="compare against the existing file and write nothing")
A = ap.parse_args()

D = A.tile if os.path.sep in A.tile else os.path.join("data", "derived", A.tile)
OUT = A.out or ("lowveg_grid.npy" if A.gen2 else "lowveg.npy")

if A.gen2:
    # Grid-wide: the same definition, evaluated on every cell instead of the cube's.
    import laspy
    from lidar_diff_icp import registration as reg
    from scipy.ndimage import distance_transform_edt
    j = reg.read_corrections(D)
    b = j["bounds"]; RES = float(j["res_m"]); X0, Y0 = b[0], b[1]
    zf = np.load(os.path.join(D, "z_after.npy")); NY, NX = zf.shape; NC = zf.size
    _zf = zf.copy(); _m = ~np.isfinite(_zf)
    if _m.any():
        _zf = _zf[tuple(distance_transform_edt(_m, return_distances=False, return_indices=True))]
    gy, gx = np.gradient(_zf, RES)
    gxf = gx.ravel(); gyf = gy.ravel(); zflat = _zf.ravel()
    nn = np.sqrt(gxf ** 2 + gyf ** 2 + 1.0)
    DZ = 0.02; ZLO, ZHI = -1.0, 2.0          # the cube's window, so the two files agree
    edges = np.arange(ZLO, ZHI + 0.5 * DZ, DZ); NZ = edges.size - 1
    H = np.zeros((NC, NZ), np.int32); n_in = 0
    with laspy.open(A.gen2) as f:
        for pts in f.chunk_iterator(A.chunk):
            x = np.asarray(pts.x); y = np.asarray(pts.y); z = np.asarray(pts.z)
            ix = ((x - X0) / RES).astype(np.int64); iy = ((y - Y0) / RES).astype(np.int64)
            ing = (ix >= 0) & (ix < NX) & (iy >= 0) & (iy < NY)
            cc = iy[ing] * NX + ix[ing]
            xc = X0 + ((cc % NX) + 0.5) * RES; yc = Y0 + ((cc // NX) + 0.5) * RES
            h = (z[ing] - (zflat[cc] + gxf[cc] * (x[ing] - xc) + gyf[cc] * (y[ing] - yc))) / nn[cc]
            zi = np.floor((h - ZLO) / DZ).astype(np.int64)
            mm = (zi >= 0) & (zi < NZ)
            np.add.at(H, (cc[mm], zi[mm]), 1)
            n_in += int(mm.sum())
    mid = 0.5 * (edges[:-1] + edges[1:])
    band = (mid > A.lo) & (mid <= A.hi)
    tot = H.sum(1).astype(float)
    with np.errstate(invalid="ignore", divide="ignore"):
        out = np.where(tot > 0, H[:, band].sum(1) / tot, np.nan).reshape(NY, NX)
    print(f"lowveg = (returns with {A.lo:g} < h <= {A.hi:g} m) / (returns in "
          f"{ZLO:+.2f}..{ZHI:+.2f} m), bin-centre test on {DZ:.3f} m bins")
    print(f"  population: ALL gen2 returns in the window (not class-2)")
    print(f"  {n_in:,} returns placed; cells with a value: {int(np.isfinite(out).sum()):,} "
          f"of {out.size:,}")
    q = np.nanpercentile(out, [10, 50, 90])
    print(f"  lowveg  p10 {q[0]:.4f}  median {q[1]:.4f}  p90 {q[2]:.4f}  max {np.nanmax(out):.4f}")
    dst = os.path.join(D, OUT)
    if A.check:
        old = np.load(dst)
        same = np.array_equal(np.nan_to_num(old, nan=-9e9), np.nan_to_num(out, nan=-9e9))
        print(f"  --check: {'IDENTICAL' if same else 'DIFFERS'} from {dst}; nothing written")
        raise SystemExit(0 if same else 1)
    np.save(dst, out)
    print(f"wrote {dst}")
    raise SystemExit(0)

cube_path = os.path.join(D, A.cube)
if not os.path.exists(cube_path):
    raise SystemExit(f"{cube_path} is missing. It is the product of "
                     f"analysis/ridgelines/nearground_cells.py --tile {D} "
                     f"--out {A.cube}; lowveg is a slice of that cube and cannot be "
                     f"computed without it.")

z = np.load(cube_path)
cells, edges, H2 = z["cells"], z["edges"], z["H2"].astype(float)
mid = 0.5 * (edges[:-1] + edges[1:])
band = (mid > A.lo) & (mid <= A.hi)              # bin-CENTRE test, as at the control marks

tot = H2.sum(1)
num = H2[:, band].sum(1)
with np.errstate(invalid="ignore", divide="ignore"):
    lv_cells = np.where(tot > 0, num / tot, np.nan)

ny, nx = np.load(os.path.join(D, "z_after.npy")).shape
out = np.full(ny * nx, np.nan)
out[cells] = lv_cells
out = out.reshape(ny, nx)

fin = np.isfinite(out)
print(f"lowveg = (returns with {A.lo:g} < h <= {A.hi:g} m) / (returns in "
      f"{edges[0]:+.2f}..{edges[-1]:+.2f} m), bin-centre test on {edges[1]-edges[0]:.3f} m bins")
print(f"  population: H2 = ALL gen2 returns in the window (not class-2)")
print(f"  band covers {int(band.sum())} of {band.size} bins, "
      f"centres {mid[band][0]:+.3f}..{mid[band][-1]:+.3f} m")
print(f"  cells with a value: {int(fin.sum()):,} of {out.size:,} grid cells "
      f"({int((tot == 0).sum()):,} cube cells had no returns)")
if fin.any():
    q = np.nanpercentile(out, [10, 50, 90])
    print(f"  lowveg  p10 {q[0]:.4f}  median {q[1]:.4f}  p90 {q[2]:.4f}  max {np.nanmax(out):.4f}")

dst = os.path.join(D, OUT)
if A.check:
    if not os.path.exists(dst):
        raise SystemExit(f"--check: {dst} does not exist")
    old = np.load(dst)
    same = np.array_equal(np.nan_to_num(old, nan=-9e9), np.nan_to_num(out, nan=-9e9))
    print(f"  --check: {'IDENTICAL' if same else 'DIFFERS'} from {dst}; nothing written")
    raise SystemExit(0 if same else 1)

np.save(dst, out)
print(f"wrote {dst}")

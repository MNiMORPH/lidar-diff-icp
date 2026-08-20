"""Correct shaded relief via GDAL.

Hand-rolling a hillshade from ``matplotlib.colors.LightSource`` is easy to get
wrong: the illumination azimuth is applied in the array's row/column frame, so on
an ``origin='lower'`` array (row 0 = south, as our grids are) the north-south
gradient is inverted and NW light silently becomes SW light -- the terrain reads
pseudoscopically (valleys look like ridges). Rather than track flips by hand, we
let **GDAL** (``gdaldem hillshade``) compute it from the raster's geotransform
using Horn's algorithm, which is the standard and cannot be mis-oriented.
"""
from __future__ import annotations

import glob
import os
import shutil
import subprocess
import tempfile

import numpy as np
from scipy.ndimage import distance_transform_edt as edt


def _find_gdaldem(gdaldem=None):
    if gdaldem and os.path.exists(gdaldem):
        return gdaldem
    found = shutil.which(gdaldem or "gdaldem")
    if found:
        return found
    for pat in ("~/anaconda3/bin/gdaldem", "~/miniconda3/bin/gdaldem",
                "/opt/conda/bin/gdaldem", "/usr/bin/gdaldem"):
        for p in glob.glob(os.path.expanduser(pat)):
            return p
    raise FileNotFoundError("gdaldem not found (install GDAL, or pass gdaldem=<path>)")


def hillshade(dem, res, x0, y0, *, az=315, alt=45, vert_exag=2.0, crs="EPSG:26915",
              gdaldem=None, fill_gaps=False):
    """NW (315/45) shaded relief in [0, 1] via ``gdaldem hillshade``.

    ``dem``: elevation grid, ``origin='lower'`` (row 0 = south) to match how our
    DoDs are displayed. ``res``: cell size (m); ``x0, y0``: SW corner (the grid
    origin). Returns an array in the SAME orientation as ``dem`` (so use it with
    ``imshow(..., origin='lower')``); NaNs where ``dem`` is NaN. Illumination is
    fixed by the written geotransform, so it is always geographically NW.

    ``fill_gaps``: for a VISUALIZATION backdrop only. Gaps in ``dem`` (water,
    dropouts) are always nearest-filled before the shading is computed (to avoid
    edge artifacts); by default they are re-masked to NaN afterwards. Set
    ``fill_gaps=True`` to KEEP the filled shading, so isolated no-data cells do not
    show as distracting white squares behind a semi-transparent data overlay. This
    only fills the relief backdrop -- never fabricate data in the overlay itself.
    """
    import rasterio
    from rasterio.transform import from_origin
    gdaldem = _find_gdaldem(gdaldem)
    ny, nx = dem.shape
    m = np.isfinite(dem)
    fill = dem.copy()
    if (~m).any():                                   # fill gaps by nearest so the
        fill = fill[tuple(edt(~m, return_distances=False, return_indices=True))]
    d = tempfile.mkdtemp(prefix="hs_")
    try:
        demtif = os.path.join(d, "dem.tif"); hstif = os.path.join(d, "hs.tif")
        # write north-up (flipud): the geotransform's top-left is (x0, y0 + ny*res)
        with rasterio.open(demtif, "w", driver="GTiff", height=ny, width=nx, count=1,
                           dtype="float32", crs=crs,
                           transform=from_origin(x0, y0 + ny * res, res, res)) as w:
            w.write(np.flipud(fill).astype("float32"), 1)
        subprocess.run([gdaldem, "hillshade", "-az", str(az), "-alt", str(alt),
                        "-z", str(vert_exag), "-compute_edges", demtif, hstif],
                       check=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
        with rasterio.open(hstif) as r:
            hs = np.flipud(r.read(1).astype(float))  # back to origin='lower'
    finally:
        shutil.rmtree(d, ignore_errors=True)
    return np.where(m, hs / 255.0, np.nan)

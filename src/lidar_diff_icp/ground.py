"""Bare-earth ground classification for the before-epoch cloud (CSF via PDAL).

Why this exists: our default bare-earth heuristic (last return + low percentile per
cell) is fast and, on open/gentle terrain, as good as anything -- but a physically
based ground filter gives a cleaner, more general bare-earth (removes buildings,
structures, and forest understory the heuristic can keep). We use PDAL's Cloth
Simulation Filter (CSF; Zhang et al., 2016).

Two practical wrinkles this module handles:
1. **Old-LAZ compatibility.** The 2008-era MN lidar is delivered in an old LAZ
   chunk-table format that current PDAL reads the *header* of but cannot
   decompress ("Invalid version ... in LAZ chunk table"). laspy reads it fine, so
   we first rewrite the cloud to an uncompressed LAS that PDAL can read.
2. **We DEPART FROM PDAL's DEFAULTS -- read this before using the output.**
   PDAL's ``filters.csf`` defaults are ``rigidness=3``, ``threshold=0.5``,
   ``hdiff=0.3``, ``resolution=1``. This module's defaults are ``rigidness=1``,
   ``threshold=1.5``, ``hdiff=0.5``, ``resolution=1.0`` -- i.e. the SOFTEST cloth,
   a classification distance 3x looser, and a looser height difference. All three
   departures push the same way: keep more points and follow lower ones.

   They were chosen for sparse 2008 data on steep, wooded terrain, where they
   retain ~96% of 5 m cells including ~94% of steep (>=20 deg) cells while still
   removing ~12% of points (canopy/structure/blunders). A coarser cloth (>=1.5 m)
   or tighter thresholds delete steep wooded hillslopes -- do not change these
   without re-checking steep-cell coverage.

   **Known consequence.** At ``threshold=1.5`` against a near-ground return column
   whose IQR is 69 mm (open) to 159 mm (dense canopy), the threshold is ~10 IQRs
   out and filters essentially nothing: 96.5% of gen1 returns in a -1..+2 m window
   are retained as ground. The per-cell median then does all the work, and it is
   computed over a population that still contains sub-ground returns (25% of gen1's
   window lies >0.1 m below its own ground, vs 10% for gen2). Since nothing
   scatters from below the ground, those are error, and retaining them biases the
   gen1 ground LOW -- more so under canopy, where the column is ~2x wider. This
   parameter choice therefore propagates to every tile processed with this module.
   See analysis/ridgelines/FRAME_2026-08-26.md.
"""
from __future__ import annotations

import glob
import json
import os
import shutil
import subprocess
import tempfile

import laspy


def find_pdal(pdal=None):
    """Locate a PDAL binary. Checks an explicit path, then PATH, then common
    conda env locations. Raises FileNotFoundError with guidance if not found."""
    if pdal and os.path.exists(pdal):
        return pdal
    found = shutil.which(pdal or "pdal")
    if found:
        return found
    for pat in ("~/anaconda3/envs/*/bin/pdal", "~/miniconda3/envs/*/bin/pdal",
                "~/mambaforge/envs/*/bin/pdal", "/opt/conda/envs/*/bin/pdal"):
        for p in glob.glob(os.path.expanduser(pat)):
            if os.path.exists(p):
                return p
    raise FileNotFoundError(
        "PDAL (with filters.csf) not found. Install PDAL or pass pdal=<path>. "
        "The conda 'lidar-icp' env has it (e.g. ~/anaconda3/envs/lidar-icp/bin/pdal).")


def classify_ground_csf(in_path, out_path=None, *, pdal=None, resolution=1.0,
                        rigidness=1, threshold=1.5, hdiff=0.5, smooth=True,
                        iterations=500, elm=True, outlier=False):
    """Classify ground with PDAL CSF; return a path to a LAS of ground points.

    Rewrites ``in_path`` to uncompressed LAS via laspy first (old-LAZ
    compatibility), then runs, in order:

    * ``filters.elm`` (default ``elm=True``) -- Extended Local Minimum, marks LOW
      blunder points (multipath, sub-canopy noise, negative-noise below true
      ground) as noise. CSF is vulnerable to exactly these: one low point drags
      the cloth down and makes it miss ground or carve a pit. Standard, protective,
      cheap; disable with ``elm=False``.
    * ``filters.outlier`` (optional, ``outlier=True``) -- statistical isolated
      high/low points.
    * a drop of the noise class before the ground filter, then ``filters.csf``,
      then keep ``Classification == 2``.

    All point attributes (gps_time, point_source_id, return numbers) are preserved.
    If ``out_path`` is None a temp file is created; the caller owns it (and the
    intermediate LAS is removed automatically). CSF at 1 m resolution is slow
    (minutes per tile) -- that is the cost of the default.
    """
    pdal = find_pdal(pdal)
    tmpdir = tempfile.mkdtemp(prefix="csf_")
    las_in = os.path.join(tmpdir, "in.las")
    laspy.read(str(in_path)).write(las_in)            # old-LAZ -> LAS PDAL can read
    if out_path is None:
        out_path = os.path.join(tmpdir, "ground.las")
    stages = [{"type": "readers.las", "filename": las_in}]
    if elm:                                           # low blunders -> noise (class 7)
        stages.append({"type": "filters.elm"})
    if outlier:                                       # isolated high/low -> class 7
        stages.append({"type": "filters.outlier"})
    if elm or outlier:                                # drop noise before ground filter
        stages.append({"type": "filters.expression", "expression": "Classification != 7"})
    stages += [
        {"type": "filters.csf", "resolution": resolution, "rigidness": rigidness,
         "threshold": threshold, "hdiff": hdiff, "smooth": smooth,
         "iterations": iterations},
        {"type": "filters.expression", "expression": "Classification == 2"},
        {"type": "writers.las", "filename": out_path},
    ]
    pipe = {"pipeline": stages}
    pj = os.path.join(tmpdir, "csf.json")
    with open(pj, "w") as fh:
        json.dump(pipe, fh)
    try:
        subprocess.run([pdal, "pipeline", pj], check=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
    except subprocess.CalledProcessError as e:
        shutil.rmtree(tmpdir, ignore_errors=True)
        raise RuntimeError("PDAL CSF failed: " + e.stderr.decode(errors="replace")) from e
    os.remove(las_in)                                 # drop the large intermediate
    return out_path

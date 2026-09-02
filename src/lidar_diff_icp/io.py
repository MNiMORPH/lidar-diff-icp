"""Reading MN 2008 lidar tiles into plain arrays, with the CRS assigned.

The 2008 LAZ files are point format 1 (GPS time present) but carry no
georeferencing VLR, so the CRS is attached here explicitly. laspy reads these
old-laszip files via either the ``laszip`` (apt) or ``lazrs`` backend; PDAL's
lazperf backend does not (transcode to LAS first if PDAL is needed).
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import laspy

from . import MN_GEN1_CRS


@dataclass
class PointCloud:
    """Per-point arrays for one tile. Coordinates are in ``crs`` (UTM 15N)."""

    x: np.ndarray
    y: np.ndarray
    z: np.ndarray
    point_source_id: np.ndarray  # flight line
    classification: np.ndarray
    gps_time: np.ndarray
    scan_angle: np.ndarray
    crs: str

    def __len__(self) -> int:
        return self.x.size

    @property
    def swaths(self) -> np.ndarray:
        """Distinct flight-line ids present, sorted."""
        return np.unique(self.point_source_id)


def read_tile(path: str | Path, crs: str = MN_GEN1_CRS) -> PointCloud:
    """Read a LAZ/LAS tile into a :class:`PointCloud`.

    ``gps_time`` and ``scan_angle`` are NaN when the file does not record them, never
    zero. Zero is a measurement here -- scan angle 0 is nadir, gps_time 0 is a timestamp --
    and the previous fallback wrote exactly those, so a file lacking the dimension read as
    "every return was taken at nadir, at time zero". That is not hypothetical: point format
    6 and later drop ``scan_angle_rank`` in favour of a scaled ``scan_angle``, and the CSF
    caches this project writes are PF7, so every one of them hit the fallback. Consumers of
    ``pc.scan_angle`` include the across-track intercept tie in
    :func:`coreg.coregister_swaths`, which fits against ``tan(scan_angle)``.

    PF6+ ``scan_angle`` is stored in 0.006-degree units; ``scan_angle_rank`` (PF<=5) is
    already in degrees. Both are returned in DEGREES.
    """
    f = laspy.read(str(path))
    dims = set(f.point_format.dimension_names)
    n = len(f.x)

    if "scan_angle" in dims:                       # PF6+, 0.006-deg units
        scan_angle = np.asarray(f.scan_angle).astype(float) * 0.006
    elif "scan_angle_rank" in dims:                # PF<=5, already degrees
        scan_angle = np.asarray(f.scan_angle_rank).astype(float)
    else:
        scan_angle = np.full(n, np.nan)

    return PointCloud(
        x=np.asarray(f.x),
        y=np.asarray(f.y),
        z=np.asarray(f.z),
        point_source_id=np.asarray(f.point_source_id),
        classification=np.asarray(f.classification),
        gps_time=np.asarray(f.gps_time) if "gps_time" in dims else np.full(n, np.nan),
        scan_angle=scan_angle,
        crs=crs,
    )

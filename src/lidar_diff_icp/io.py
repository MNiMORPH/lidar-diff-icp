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
    """Read a LAZ/LAS tile into a :class:`PointCloud`."""
    f = laspy.read(str(path))
    dims = set(f.point_format.dimension_names)
    return PointCloud(
        x=np.asarray(f.x),
        y=np.asarray(f.y),
        z=np.asarray(f.z),
        point_source_id=np.asarray(f.point_source_id),
        classification=np.asarray(f.classification),
        gps_time=np.asarray(f.gps_time) if "gps_time" in dims
        else np.zeros(len(f.x)),
        scan_angle=np.asarray(f.scan_angle_rank) if "scan_angle_rank" in dims
        else np.zeros(len(f.x), dtype=np.int8),
        crs=crs,
    )

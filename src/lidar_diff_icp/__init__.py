"""Change detection between lidar epochs with per-flight-line correction.

The earlier (2008) survey carries per-swath navigation error. Adjacent swaths
were flown minutes to hours apart, so overlap discrepancies are pure acquisition
error (no real change) and self-calibrate the correction. See the package README.
"""

__version__ = "0.0.1"

# Coordinate reference system of the MN 2008 SE lidar (not embedded in the files).
MN_2008_CRS = "EPSG:26915"  # UTM zone 15N, NAD83

from . import tiles, io, swathdiff  # noqa: E402,F401

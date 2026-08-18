"""Change detection between lidar epochs with per-flight-line correction.

The MnGeo **First-Generation** statewide lidar (``gen1``, flown 2008-2012) carries
per-swath navigation error. Adjacent swaths were flown minutes to hours apart, so
overlap discrepancies are pure acquisition error (no real change) and
self-calibrate the correction against modern **Second-Generation** 3DEP
(``gen2``, 2020s). See the package README.
"""

__version__ = "0.0.1"

# CRS assigned to the First-Generation MN lidar tiles (not embedded in the files):
# UTM zone 15N, NAD83, the MnGeo statewide standard.
MN_GEN1_CRS = "EPSG:26915"
MN_2008_CRS = MN_GEN1_CRS       # backward-compatible alias (deprecated: 'gen1' spans 2008-2012)

from . import tiles, io, swathdiff  # noqa: E402,F401

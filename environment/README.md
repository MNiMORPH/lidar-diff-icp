# Environment

The processing environment is a Python **venv built with
`--system-site-packages`** on the system interpreter, layered over
apt-installed geospatial libraries. This reuses the distribution's compiled
GDAL/PROJ/laspy stack rather than re-downloading it, and keeps the project
isolated for anything installed on top.

## 1. System packages (apt)

```bash
sudo apt install $(grep -v '^#' apt-packages.txt)
```

`cloudcompare` and `qgis` are standalone GUI apps (ICP/M3C2 inspection, DEM and
tile-index viewing). The `python3-*` packages are the libraries the venv sees.
Note: PDAL is **not** in apt — see the conda note below if you need it.

## 2. Project venv

```bash
python3 -m venv --system-site-packages lidar-icp
source lidar-icp/bin/activate
pip install -e .          # installs this package; deps are satisfied by apt
```

## Gotchas

- **PROJ/GDAL data path.** If your shell auto-activates conda `base`, it exports
  `PROJ_DATA`/`GDAL_DATA` pointing at Anaconda's copies, which breaks the venv's
  system PROJ/GDAL ("Internal Proj Error: ... no database context"). Fix once
  with `conda config --set auto_activate_base false` (affects new shells). As a
  per-call fallback: `env PROJ_DATA=/usr/share/proj GDAL_DATA=/usr/share/gdal`.

- **PDAL (optional).** PDAL is not available via apt and cannot be pip-installed
  without `libpdal`. If you need it (e.g. `readers.ept` for 3DEP, `filters.icp`),
  create a conda env instead: `conda create -n lidar-icp-pdal -c conda-forge
  python pdal python-pdal gdal rasterio laspy lazrs-python`. Note that PDAL's
  lazperf backend cannot read the old-laszip 2008 tiles directly — transcode
  LAZ→LAS with laspy first.

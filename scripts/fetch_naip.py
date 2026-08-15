#!/usr/bin/env python3
"""Fetch a NAIP 4-band mosaic for a bounding box from Microsoft Planetary
Computer, as an *independent* land-cover source (not derived from the lidar).

Planetary Computer serves NAIP as signed 4-band (R,G,B,NIR) COGs in EPSG:26915
for Minnesota, so this reads only the windows overlapping the bbox and mosaics
them at the requested resolution. Saves an .npz with the bands and NDVI.

Example (2010, ~contemporaneous with the fall-2008 lidar; cover is stable):
    python scripts/fetch_naip.py --bounds 577492.8 4882737.6 580035.0 4886238.3 \
        --year 2010 --res 2 --out data/naip/naip2010_2m.npz
"""
import os
# If a conda base leaks PROJ_DATA/GDAL_DATA, run with:
#   env PROJ_DATA=/usr/share/proj GDAL_DATA=/usr/share/gdal python scripts/fetch_naip.py ...
os.environ.update(GDAL_HTTP_TIMEOUT="60", GDAL_HTTP_CONNECTTIMEOUT="15",
                  GDAL_HTTP_MAX_RETRY="5", GDAL_HTTP_RETRY_DELAY="2",
                  GDAL_DISABLE_READDIR_ON_OPEN="EMPTY_DIR",
                  CPL_VSIL_CURL_ALLOWED_EXTENSIONS=".tif", VSI_CACHE="TRUE")
import argparse, json, urllib.request, warnings
from pathlib import Path
import numpy as np
warnings.filterwarnings("ignore")

STAC = "https://planetarycomputer.microsoft.com/api/stac/v1/search"
SAS = "https://planetarycomputer.microsoft.com/api/sas/v1/sign?href="


def _get(url):
    return json.load(urllib.request.urlopen(url, timeout=60))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bounds", nargs=4, type=float, required=True,
                    metavar=("MINX", "MINY", "MAXX", "MAXY"), help="EPSG:26915")
    ap.add_argument("--year", type=int, default=2010)
    ap.add_argument("--res", type=float, default=2.0)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    minx, miny, maxx, maxy = a.bounds

    from pyproj import Transformer
    t = Transformer.from_crs("EPSG:26915", "EPSG:4326", always_xy=True)
    lo = t.transform(minx, miny); hi = t.transform(maxx, maxy)
    bbox = [min(lo[0], hi[0]), min(lo[1], hi[1]), max(lo[0], hi[0]), max(lo[1], hi[1])]
    feats = _get(f"{STAC}?collections=naip&bbox={','.join(map(str,bbox))}"
                 f"&datetime={a.year}-01-01/{a.year}-12-31&limit=50")["features"]
    print(f"{len(feats)} NAIP {a.year} items over bbox", flush=True)

    import rasterio
    from rasterio.windows import from_bounds
    from rasterio.transform import from_origin
    nx = int(round((maxx - minx) / a.res))
    ny = int(round((maxy - miny) / a.res))
    rgbn = np.zeros((4, ny, nx), np.uint8)
    # Per-tile boundless windowed read onto the common target grid, then
    # composite (quads are non-overlapping, so max keeps the valid pixel).
    for f in feats:
        href = f["assets"]["image"]["href"]
        a4 = None
        for attempt in range(1, 4):                      # retry flaky /vsicurl reads
            try:
                signed = _get(SAS + urllib.request.quote(href, safe=''))["href"]
                with rasterio.open(signed) as src:
                    win = from_bounds(minx, miny, maxx, maxy, transform=src.transform)
                    a4 = src.read(range(1, 5), window=win, boundless=True,
                                  fill_value=0, out_shape=(4, ny, nx))
                break
            except Exception as exc:
                print(f"  {f['id']} attempt {attempt} failed: {exc}", flush=True)
        if a4 is None:
            raise RuntimeError(f"could not read {f['id']} after retries")
        rgbn = np.maximum(rgbn, a4.astype(np.uint8))
        print(f"  read {f['id']}", flush=True)
    transform = from_origin(minx, maxy, a.res, a.res)
    mosaic = rgbn
    R = mosaic[0].astype(np.float32); NIR = mosaic[3].astype(np.float32)
    ndvi = (NIR - R) / (NIR + R + 1e-6)
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(a.out, rgbn=rgbn, ndvi=ndvi.astype(np.float32),
                        transform=np.array(transform).reshape(-1)[:6],
                        bounds=np.array([minx, miny, maxx, maxy]),
                        crs="EPSG:26915", res=a.res)
    print(f"wrote {a.out}  shape={rgbn.shape}  NDVI range [{ndvi.min():.2f},{ndvi.max():.2f}]",
          flush=True)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Fetch the public-land and access layers used by field_sites_rtk.py (needs the network).

Source: Minnesota Geospatial Commons feature services (ArcGIS REST, EPSG:26915 in/out):
  WMA boundaries  us_mn_state_dnr/bdry_dnr_wildlife_mgmt_areas_pub   (Publicly Accessible State WMAs)
  WMA facilities  us_mn_state_dnr/struc_dnr_wld_mgmt_areas_pub        layers 0 points (parking), 1 lines (roads)
  State parks     us_mn_state_dnr/bdry_dnr_lrs_prk                    (none intersect the Elba box)
The GDRS zip downloads for these datasets 404 as of 2026-10-18; the feature services work.

    ./lidar-icp/bin/python analysis/fetch_public_land.py
"""
import os
import subprocess

OUT = "data/derived/public_land"; os.makedirs(OUT, exist_ok=True)
BB = "574000,4880000,582000,4888000"
Q = (f"where=1%3D1&geometry={BB}&geometryType=esriGeometryEnvelope&inSR=26915&spatialRel=esriSpatialRelIntersects"
     "&outFields=*&outSR=26915&f=geojson")
BASE = "https://enterprise.gisdata.mn.gov/aghost/rest/services/us_mn_state_dnr"
for name, url in (("whitewater_wma_units.geojson", f"{BASE}/bdry_dnr_wildlife_mgmt_areas_pub/FeatureServer/0/query"),
                  ("wma_fac_0.geojson", f"{BASE}/struc_dnr_wld_mgmt_areas_pub/FeatureServer/0/query"),
                  ("wma_fac_1.geojson", f"{BASE}/struc_dnr_wld_mgmt_areas_pub/FeatureServer/1/query"),
                  ("state_parks.geojson", f"{BASE}/bdry_dnr_lrs_prk/FeatureServer/2/query")   # layer 2 = management-unit POLYGONS; layer 0 is reference points, layer 1 statutory boundaries):
    subprocess.run(["curl", "-s", "-m", "180", f"{url}?{Q}", "-o", f"{OUT}/{name}"], check=True)
    print(name, os.path.getsize(f"{OUT}/{name}"), "bytes")

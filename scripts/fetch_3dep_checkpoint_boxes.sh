#!/usr/bin/env bash
# Fetch one full-density gen2 (2021 3DEP) box around each surveyed QA checkpoint near
# Elba, so gen2 can be tested DIRECTLY against the marks that certify it.
#
#   bash scripts/fetch_3dep_checkpoint_boxes.sh [half_width_m]
#
# Coordinates come from the bundled checkpoint CSV
# (src/lidar_diff_icp/groundtruth/data/mn_se_driftless_2021_ql1_near_elba.csv); keep the
# two in step if a mark is added. --max-depth 12 is the FULL-density fetch the products
# use (README quick start step 1); the default depth 9 thins ground returns by ~19x and
# would measure a different surface (see the decimation note in
# analysis/decimation_result.md).
#
# Six 400 m boxes are ~57 MB and ~1 minute total on a warm connection.
set -euo pipefail
cd "$(dirname "$0")/.."

HW="${1:-200}"
OUT=data/after/checkpoints
mkdir -p "$OUT"

CSV=src/lidar_diff_icp/groundtruth/data/mn_se_driftless_2021_ql1_near_elba.csv

tail -n +2 "$CSV" | while IFS=, read -r pid ptype e n rest; do
    short="${pid%%_*}"
    dst="$OUT/cp${short}_gen2.laz"
    if [ -f "$dst" ]; then
        echo "== $pid: $dst exists, skipping"
        continue
    fi
    echo "== $pid ($ptype) at $e $n -> $dst"
    env PROJ_DATA=/usr/share/proj GDAL_DATA=/usr/share/gdal \
        ./lidar-icp/bin/python scripts/fetch_3dep_curl.py --auto \
        --bounds "$(echo "$e - $HW" | bc)" "$(echo "$n - $HW" | bc)" \
                 "$(echo "$e + $HW" | bc)" "$(echo "$n + $HW" | bc)" \
        --max-depth 12 --workers 8 \
        --ept-cache data/derived/groundtruth/ept_boundaries.json \
        --out "$dst"
done

echo
echo "boxes on disk:"
ls -la "$OUT"/*.laz

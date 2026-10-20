#!/usr/bin/env bash
# Phase B loop: the 29 public-land Whitewater tiles, one at a time, gently (D6, D7).
# Order: the WMA Main Branch reach first (where the field sites are), then the rest.
# Skips a tile whose valley_run.log already ends in DONE; stops if free memory is short.
set -u
cd "$(dirname "$0")/.."
TILES="4342-29-64 4342-28-64 4358-28-01 4358-29-01 4358-27-01 4342-27-64 4358-26-01 4358-26-02 4358-26-03 4358-27-02 4358-29-02 4358-30-01 4358-30-02 4342-30-64 4342-30-63 4342-29-63 4342-29-62 4342-30-62 4342-29-61 4342-29-60 4342-29-59 4342-30-60 4342-30-61 4342-31-62 4342-31-63 4358-31-01 4358-31-02 4358-32-01 4358-32-02"
for t in $TILES; do
  if grep -qs "DONE; stages" data/derived/$t/valley_run.log; then echo "== $t done"; continue; fi
  free_mb=$(free -m | awk '/Mem:/{print $7}')
  if [ "$free_mb" -lt 4000 ]; then echo "== $t: only ${free_mb} MB available; stopping the loop here ($(date))"; exit 2; fi
  echo "== $t start $(date) (available ${free_mb} MB)"
  mkdir -p data/derived/$t
  nice -n 10 ./lidar-icp/bin/python scripts/run_valley_tile.py --tile $t >> data/derived/$t/valley_run.log 2>&1
  rc=$?; echo "== $t exit $rc $(date)"
  if [ $rc -ne 0 ]; then echo "== $t FAILED; continuing with the next tile"; fi
done
echo "== LOOP DONE $(date)"

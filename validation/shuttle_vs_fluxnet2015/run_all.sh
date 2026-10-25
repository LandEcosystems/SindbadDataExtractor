#!/bin/bash
# Compare the FLUXNET Shuttle extraction with the FLUXNET2015 SINDBAD site files for
# every site on both lists, then write the dashboard csv files.
#
# usage: run_all.sh <work folder> [python] [parallel jobs]
# Set SINDBAD_VALIDATION_DATA when the data are not in ~/research/RnD/data.
set -e
HERE="$(cd "$(dirname "$0")" && pwd)"
WORK="${1:?usage: run_all.sh <work folder> [python] [parallel jobs]}"
PY="${2:-$HOME/miniforge3/envs/sinproc-prefetch/bin/python}"
JOBS="${3:-4}"
DATA="${SINDBAD_VALIDATION_DATA:-$HOME/research/RnD/data}"
IDX="$DATA/sindbad_ancillary_data/prefetch/fluxnet_shuttle/index_fluxnet_shuttle.csv"
mkdir -p "$WORK"

# Sites on the shuttle and in the daily or hourly FLUXNET2015 files
tail -n +2 "$IDX" | cut -d, -f1 | sort -u > "$WORK/shuttle_sites.txt"
(ls "$DATA/fn"; ls "$DATA/hourly") | sed 's/\..*//' | sort -u > "$WORK/fn_sites.txt"
comm -12 "$WORK/shuttle_sites.txt" "$WORK/fn_sites.txt" > "$WORK/common.txt"
echo "$(wc -l < "$WORK/common.txt") common sites"

cd "$HERE"
xargs -P "$JOBS" -I{} "$PY" -W ignore compare_site.py {} "$WORK/out" \
    < "$WORK/common.txt" | grep -c "^done"
"$PY" aggregate.py "$WORK"

ls "$WORK/out" | grep -E '\.csv$' | sed -E 's/_(daily|hourly)\.csv$//' | sort -u \
    | xargs -P "$JOBS" -I{} "$PY" -W ignore era_check.py {} "$WORK/era"

"$PY" make_dashboard_data.py "$WORK"

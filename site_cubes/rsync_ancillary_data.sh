#!/bin/bash
# Copy the local ancillary data to ancillary_data of the site cubes on the cluster.
# Run it on the local machine. Files already copied are skipped, and nothing is
# deleted on either side. Extra arguments go to rsync, such as -n for a dry run or
# -P to resume large files and show progress.

SRC=~/research/RnD/data/sindbad_ancillary_data
DST=skoirala@dialog.bgc-jena.mpg.de:/Net/Groups/BGI/work_3/sindbad/site_cubes/ancillary_data
rsync -av --exclude '*.lock' --exclude 'caches_before_*' \
    --exclude era5_API --exclude era5land_API "$@" "$SRC/" "$DST/"

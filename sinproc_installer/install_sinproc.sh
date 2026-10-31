#!/bin/bash
# The one conda environment of the site cubes: Snakemake, staging, the cubes, the
# remote products of prefetch.py, fluxcom for the v0.07 releases and PyCLIFF.
#
#   ./install_sinproc.sh [prefix]
#
# The prefix defaults to tools/env of the site cubes, which is what envs of
# configuration/pipeline.json and pycliff_python of the CLIFF configs point to.
# fluxcom and PyCLIFF are installed in development mode from their clones next to
# this repository, tools/fluxcom and tools/PyCLIFF, which FLUXCOM_DIR and PYCLIFF_DIR
# override. The environment is made again when it exists.
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
tools=$(dirname "$(dirname "$here")")
prefix=${1:-/Net/Groups/BGI/work_3/sindbad/site_cubes/tools/env}
fluxcom=${FLUXCOM_DIR:-$tools/fluxcom}
pycliff=${PYCLIFF_DIR:-$tools/PyCLIFF}

mamba env create -y -p "$prefix" -f "$here/sinproc.yml"
"$prefix/bin/python" -m pip install --no-deps -e "$fluxcom"
"$prefix/bin/python" -m pip install --no-deps -e "$pycliff"
"$prefix/bin/python" - <<'EOF'
import fluxcom.providers, fluxnet_shuttle.shuttle, missingpy, pycliff, rasterio
import snakemake, xarray, zarr
print("environment ready")
EOF

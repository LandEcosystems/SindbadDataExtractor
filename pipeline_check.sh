#!/bin/bash
# What each step of the site cubes has done, from pipeline_check.py. Nothing is
# fetched or extracted.
#
#   ./pipeline_check.sh                    every step
#   ./pipeline_check.sh setup registry     only these steps
#   ./pipeline_check.sh store --hash       also the sites the next run extracts again
#   ./pipeline_check.sh store --products shuttle_obs ndvi_landsat
#
# The steps are setup, registry, prefetch, store and cubes, as in README.md. PYTHON
# overrides the python of the environment. Run it from anywhere.
cd "$(dirname "$0")"
python=${PYTHON:-/Net/Groups/BGI/work_3/sindbad/site_cubes/tools/env/bin/python}
exec "$python" pipeline_check.py "$@"

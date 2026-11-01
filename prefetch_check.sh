#!/bin/bash
# Show the state of the prefetch: the prefetch jobs still in SLURM, then the sites
# cached and missing for each dataset, from prefetch.py --check.
#
#   ./prefetch_check.sh                   every dataset with a cache
#   ./prefetch_check.sh fa_gami vf_s2     only these
#
# Nothing is fetched. PYTHON and REGISTRY override the defaults, as in
# prefetch_array.sh. Run it from anywhere.
cd "$(dirname "$0")"
python=${PYTHON:-/Net/Groups/BGI/work_3/sindbad/site_cubes/tools/env/bin/python}
registry=${REGISTRY:-/Net/Groups/BGI/work_3/sindbad/site_cubes/ancillary_data/insitu/sites_consolidated.csv}

if command -v squeue >/dev/null 2>&1; then
    echo "== prefetch jobs in SLURM"
    squeue -u "$USER" -o "%.12i %.30j %.9T %.10M %R" | grep -E 'JOBID|prefetch_' \
        || echo "none"
    echo
fi

echo "== sites in the caches"
"$python" prefetch.py configuration/sindbad/sindbad_extraction_main.json "$@" \
    --sites-csv "$registry" --all-sites --check

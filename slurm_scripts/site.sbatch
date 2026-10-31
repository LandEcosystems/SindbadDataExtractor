#!/bin/bash
#SBATCH --job-name=sdx_site
#SBATCH --ntasks=1
#SBATCH --mail-type=FAIL

# The workflow for one site with run_site.py. Submitted by sites_as_jobs.sh, which
# passes the resources and the settings of env.sh. A site that ends without an
# error gets a marker in $SDX_OUTPUT/done. See README_SLURM.md.

# the proxy module gives the job internet access
if ! type module >/dev/null 2>&1; then source /etc/profile; fi
module load proxy

set -euo pipefail

: "${SDX_REPO:?submit through the scripts in slurm_scripts/}"
site="${1:?give the site}"
cd "$SDX_REPO"

# the forcing is prefetched for all sites at once, never by a site job
args=(--workers "$SLURM_CPUS_PER_TASK" --no-prefetch)
if [ -n "$SDX_COLLECTIONS" ]; then args+=(--collections "$SDX_COLLECTIONS"); fi
if [ -n "$SDX_FREQS" ]; then args+=(--freqs "$SDX_FREQS"); fi
if [ "$SDX_OVERWRITE" = "1" ]; then args+=(--force); fi

echo "${site} on ${SLURMD_NODENAME:-$(hostname -s)} with ${SLURM_CPUS_PER_TASK} cpus," \
    "started $(date)"
status=0
"$SDX_PYTHON" run_site.py "$site" "${args[@]}" || status=$?
if [ "$status" -eq 0 ]; then
    mkdir -p "${SDX_OUTPUT}/done"
    date > "${SDX_OUTPUT}/done/${site}"
fi
echo "${site} finished with status ${status} $(date)"
exit "$status"

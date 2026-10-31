#!/bin/bash
# The workflow for some sites on one node, with run_site.py. The sites run
# --workers at once, one per cpu of the job, and the arguments after the sites go
# to run_site.py. Jobs of different sites may run at the same time.
#
#   sbatch run_site.sh DE-Hai
#   sbatch run_site.sh DE-Hai --collections Shuttle --freqs daily
#   sbatch --cpus-per-task=32 run_site.sh @sites_01.txt
#
# Run it from the repository. The proxy module gives the jobs internet access.
# PYTHON overrides the python of the environment.
#SBATCH --job-name=run_site
#SBATCH --partition=big
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=4G
#SBATCH --time=24:00:00
#SBATCH --output=run_site_%j.out

if ! type module >/dev/null 2>&1; then source /etc/profile; fi
module load proxy

python=${PYTHON:-/Net/Groups/BGI/work_3/sindbad/site_cubes/tools/env/bin/python}
cd "${SLURM_SUBMIT_DIR:-.}"
"$python" run_site.py "$@" --workers "${SLURM_CPUS_PER_TASK:-1}"

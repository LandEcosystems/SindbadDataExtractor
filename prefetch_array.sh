#!/bin/bash
# Prefetch one dataset as a SLURM array. Task i of n fetches part i/n of the work
# into the one cache of the dataset: the years of the Copernicus NDVI and LAI
# products, the sites of Landsat or the variables of the CLIFF forcing, so the tasks
# never write the same file. Work already in the cache is skipped, so a task that
# stopped can simply be submitted again.
#
#   sbatch --array=1-13 prefetch_array.sh ndvi_clms300     # 2014 to 2026
#   sbatch --array=1-22 prefetch_array.sh ndvi_clms1km     # 1999 to 2020
#   sbatch --array=1-13 prefetch_array.sh lai_clms300
#   sbatch --array=1-15 prefetch_array.sh lai_clms1km      # 1999 to 2013
#   sbatch --array=1-40 prefetch_array.sh ndvi_landsat     # about 33 sites a task
#   sbatch --array=1-9 prefetch_array.sh cliff_forcing_era5  # a variable each
#   sbatch --array=1-10 prefetch_array.sh cliff_forcing_crujra  # a variable each
#   sbatch prefetch_array.sh fa_gami                         # no array: all of it
#
# Submitted without --array, the job fetches all the work of the dataset, which suits
# the datasets that cannot split it. prefetch_all.sh submits every dataset.
#
# The logs go to output/slurm_logs of the site cubes as <job name>_<job>_<task>.log.
# SLURM does not make that folder, and prefetch_all.sh does. Run it from the
# repository. The Copernicus products read the Copernicus Data Space
# with the cdse profile of ~/.aws/credentials, and the first task to start lists the
# products while the others wait for it. PYTHON and REGISTRY override the defaults.
# WORKERS sets the threads or processes of prefetch.py, by default the cpus of the
# job. Datasets that mostly wait on the network run more threads than cpus.
#SBATCH --job-name=prefetch
#SBATCH --partition=big
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=24:00:00
#SBATCH --output=/Net/Groups/BGI/work_3/sindbad/site_cubes/output/slurm_logs/%x_%A_%a.log

# the proxy module gives the jobs internet access
if ! type module >/dev/null 2>&1; then source /etc/profile; fi
module load proxy

dataset=$1
shift
python=${PYTHON:-/Net/Groups/BGI/work_3/sindbad/site_cubes/tools/env/bin/python}
registry=${REGISTRY:-/Net/Groups/BGI/work_3/sindbad/site_cubes/ancillary_data/insitu/sites_consolidated.csv}
part=()
if [ -n "${SLURM_ARRAY_TASK_ID:-}" ]; then
    n_tasks=${SLURM_ARRAY_TASK_COUNT:-$SLURM_ARRAY_TASK_MAX}
    part=(--part "${SLURM_ARRAY_TASK_ID}/${n_tasks}")
fi

cd "${SLURM_SUBMIT_DIR:-.}"
"$python" prefetch.py configuration/sindbad/sindbad_extraction_main.json "$dataset" \
    --sites-csv "$registry" --all-sites --workers "${WORKERS:-${SLURM_CPUS_PER_TASK:-8}}" \
    "${part[@]}" "$@"

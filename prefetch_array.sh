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
#
# Run it from the repository. The Copernicus products read the Copernicus Data Space
# with the cdse profile of ~/.aws/credentials, and the first task to start lists the
# products while the others wait for it. PYTHON and REGISTRY override the defaults.
#SBATCH --job-name=prefetch
#SBATCH --partition=big
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=24:00:00
#SBATCH --output=prefetch_%x_%A_%a.out

dataset=$1
shift
python=${PYTHON:-/Net/Groups/BGI/work_3/sindbad/site_cubes/tools/env/bin/python}
registry=${REGISTRY:-/Net/Groups/BGI/work_3/sindbad/site_cubes/ancillary_data/insitu/sites_consolidated.csv}
n_tasks=${SLURM_ARRAY_TASK_COUNT:-$SLURM_ARRAY_TASK_MAX}

cd "${SLURM_SUBMIT_DIR:-.}"
"$python" prefetch.py configuration/sindbad/sindbad_extraction_main.json "$dataset" \
    --sites-csv "$registry" --all-sites --workers "${SLURM_CPUS_PER_TASK:-8}" \
    --part "${SLURM_ARRAY_TASK_ID}/${n_tasks}" "$@"

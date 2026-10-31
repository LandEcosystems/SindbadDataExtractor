#!/bin/bash
# The whole workflow through SLURM with a share of the cluster cpus.
#
#   ./run_workflow.sh 50              half of TOTAL_CPUS, 400 of 800
#   ./run_workflow.sh 25 -n           a dry run, at a quarter
#   ./run_workflow.sh 100 store       stage every product, no cubes
#   TOTAL_CPUS=1200 ./run_workflow.sh 50
#
# The first argument is the percent of TOTAL_CPUS, 800 by default, that the running
# jobs may use together, and no single job asks for more. The other arguments go to
# snakemake. Snakemake itself runs here until the workflow ends, so start it in
# screen or tmux. The proxy module gives the jobs internet access. SNAKEMAKE
# overrides the snakemake of the environment.
set -e
percent=${1:?give the percent of the cpus to use, such as 50}
shift
total=${TOTAL_CPUS:-800}
cpus=$(( total * percent / 100 ))
if [ "$cpus" -lt 1 ]; then cpus=1; fi

if ! type module >/dev/null 2>&1; then source /etc/profile; fi
module load proxy

snakemake=${SNAKEMAKE:-/Net/Groups/BGI/work_3/sindbad/site_cubes/tools/env/bin/snakemake}
cd "$(dirname "$0")"
echo "running with $cpus of $total cpus ($percent%)"
export SINDBAD_MAX_CPUS=$cpus
# targets and options first: --resources takes several values and would swallow them
exec "$snakemake" "$@" --profile profiles/slurm --resources cpus="$cpus"

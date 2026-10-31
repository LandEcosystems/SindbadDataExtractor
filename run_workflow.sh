#!/bin/bash
# The whole workflow through SLURM with a share of the free cpus of the partition.
#
#   ./run_workflow.sh                 75% of the idle cpus of the partition
#   ./run_workflow.sh 50              50% of them
#   ./run_workflow.sh 50 -n           a dry run
#   ./run_workflow.sh 75 store        extract every product, no cubes
#   TOTAL_CPUS=800 ./run_workflow.sh  75% of a fixed pool of 800 cpus
#
# The share is a percent of the cpus that sinfo reports idle on PARTITION, big by
# default, when the run starts, or of TOTAL_CPUS when that is set. It is fixed for
# the run, since Snakemake does not change it later. The running jobs use at most
# that many cpus together, and no single job asks for more. A first argument that is
# a number is the percent; the other arguments go to snakemake. Snakemake itself
# runs here until the workflow ends, so start it in screen or tmux. The proxy module
# gives the jobs internet access. SNAKEMAKE overrides the snakemake of the
# environment.
set -e
percent=75
if [[ "${1:-}" =~ ^[0-9]+$ ]]; then
    percent=$1
    shift
fi
partition=${PARTITION:-big}

if [ -n "${TOTAL_CPUS:-}" ]; then
    pool=$TOTAL_CPUS
    from="a pool of $pool cpus"
else
    # sinfo %C gives allocated/idle/other/total cpus of the partition
    pool=$(sinfo -h -p "$partition" -o %C | awk -F/ '{s += $2} END {print s + 0}')
    from="$pool idle cpus of $partition"
fi
cpus=$(( pool * percent / 100 ))
if [ "$cpus" -lt 1 ]; then
    echo "no cpus to use: $percent% of $from" >&2
    exit 1
fi

if ! type module >/dev/null 2>&1; then source /etc/profile; fi
module load proxy

snakemake=${SNAKEMAKE:-/Net/Groups/BGI/work_3/sindbad/site_cubes/tools/env/bin/snakemake}
cd "$(dirname "$0")"
echo "running with $cpus cpus, $percent% of $from"
export SINDBAD_MAX_CPUS=$cpus
# targets and options first: --resources takes several values and would swallow them
exec "$snakemake" "$@" --profile profiles/slurm --resources cpus="$cpus"

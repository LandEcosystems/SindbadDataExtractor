#!/bin/bash
# The whole workflow through SLURM with a share of the free cpus of the partition,
# or on this machine alone with LOCAL=1.
#
#   ./run_workflow.sh                 75% of the idle cpus of the partition
#   ./run_workflow.sh 50              50% of them
#   ./run_workflow.sh 50 -n           a dry run
#   ./run_workflow.sh 75 store        extract every product, no cubes
#   ./run_workflow.sh registry        only the site registry
#   TOTAL_CPUS=800 ./run_workflow.sh  75% of a fixed pool of 800 cpus
#   LOCAL=1 ./run_workflow.sh 50      no SLURM: the jobs run here on 50% of the cpus
#                                     of this machine
#
# The partition is SDX_PARTITION of slurm_scripts/env.sh, which every way to run
# reads, and the jobs go to it through profiles/slurm. The share is a percent of the
# cpus that sinfo reports idle on that partition when the run starts, or of
# TOTAL_CPUS when that is set. With LOCAL=1 it is a percent of the cpus of this
# machine, or of TOTAL_CPUS. It is fixed for the run, since Snakemake does not change
# it later. The running jobs use at most that many
# cpus together, and no single job asks for more. A first argument that is a number
# is the percent; the other arguments go to snakemake. Snakemake itself runs here
# until the workflow ends, so start it in screen or tmux. The proxy module gives the
# jobs internet access where the module exists. SNAKEMAKE overrides the snakemake of
# the environment.
set -e
percent=75
if [[ "${1:-}" =~ ^[0-9]+$ ]]; then
    percent=$1
    shift
fi
# only the partition of env.sh, in a subshell, so its other settings stay out
partition=$(source "$(dirname "$0")/slurm_scripts/env.sh" && echo "$SDX_PARTITION")
if [ -z "$partition" ]; then
    echo "no SDX_PARTITION in slurm_scripts/env.sh" >&2
    exit 1
fi
# the submit command of profiles/slurm reads it
export SDX_PARTITION=$partition

if [ -n "${TOTAL_CPUS:-}" ]; then
    pool=$TOTAL_CPUS
    from="a pool of $pool cpus"
elif [ -n "${LOCAL:-}" ]; then
    pool=$(getconf _NPROCESSORS_ONLN)
    from="the $pool cpus of $(hostname)"
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

if ! type module >/dev/null 2>&1 && [ -f /etc/profile ]; then source /etc/profile; fi
if type module >/dev/null 2>&1; then
    module load proxy
else
    echo "no module command here, so no proxy is loaded"
fi

snakemake=${SNAKEMAKE:-/Net/Groups/BGI/work_3/sindbad/site_cubes/tools/env/bin/snakemake}
cd "$(dirname "$0")"
export SINDBAD_MAX_CPUS=$cpus
# targets and options first: --resources takes several values and would swallow them
if [ -n "${LOCAL:-}" ]; then
    echo "running here with $cpus cpus, $percent% of $from"
    exec "$snakemake" "$@" --cores "$cpus" --resources cpus="$cpus" \
        --rerun-incomplete
fi
echo "running with $cpus cpus, $percent% of $from"
exec "$snakemake" "$@" --profile profiles/slurm --resources cpus="$cpus"

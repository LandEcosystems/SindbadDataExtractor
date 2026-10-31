#!/bin/bash
# Helpers sourced by the submit scripts to size the jobs. See README_SLURM.md.

# convert a SLURM memory value such as 12G or 500M to MB
memory_in_mb() {
    local value="$1"
    case "$value" in
        *T) echo $(( ${value%T} * 1024 * 1024 )) ;;
        *G) echo $(( ${value%G} * 1024 )) ;;
        *M) echo "${value%M}" ;;
        *) echo "$value" ;;
    esac
}

# jobs of the run are named <prefix><site>
job_prefix="${SDX_RUN_NAME}_"

# Prints "<cpus> <memory in MB>" held by the jobs of the run in the given states, such
# as RUNNING or RUNNING,PENDING.
run_usage() {
    squeue -h -u "$USER" -t "$1" -o "%j %C %m" | awk -v prefix="$job_prefix" '
        index($1, prefix) == 1 {
            cpus += $2
            mem = $3
            if (mem ~ /T$/) mem = mem * 1024 * 1024
            else if (mem ~ /G$/) mem = mem * 1024
            mem_total += mem + 0
        }
        END { printf "%d %d\n", cpus, mem_total }'
}

# Prints "<cpus> <memory in MB>" this run may hold in total: the share of the idle
# cpus and unallocated memory of the partition, plus what its running jobs hold, so
# that a submitter that runs again does not grow its share.
resource_budget() {
    local held_cpus held_mem
    read -r held_cpus held_mem <<< "$(run_usage RUNNING)"

    # CPUsState is allocated/idle/other/total; memory counts only on usable nodes
    sinfo -h -N -p "$SDX_PARTITION" -O "CPUsState:40,Memory:20,AllocMem:20,StateCompact:20" |
        awk -v partition="$SDX_PARTITION" -v percent="$SDX_USE_PERCENT" \
            -v held_cpus="$held_cpus" -v held_mem="$held_mem" '
            {
                split($1, cpus, "/")
                idle_cpus += cpus[2]
                if ($4 ~ /^(idle|mix)/) idle_mem += $2 - $3
            }
            END {
                avail_cpus = idle_cpus + held_cpus
                avail_mem = idle_mem + held_mem
                budget_cpus = int(avail_cpus * percent / 100)
                budget_mem = int(avail_mem * percent / 100)
                printf "%d cpus and %d MB memory available in %s, %d%% is %d cpus and %d MB\n", avail_cpus, avail_mem, partition, percent, budget_cpus, budget_mem > "/dev/stderr"
                printf "%d %d\n", budget_cpus, budget_mem
            }'
}

# resources of one site job, for sbatch
site_cpus="$SDX_THREADS"
site_mem=$(memory_in_mb "$SDX_MEM")
sbatch_resources=(--partition="$SDX_PARTITION" --cpus-per-task="$site_cpus"
    --mem="${site_mem}M" --time="$SDX_TIME")
if [ -n "$SDX_EXCLUDE" ]; then
    sbatch_resources+=(--exclude="$SDX_EXCLUDE")
fi

#!/bin/bash
# Submits one SLURM job per site within a share of the resources.
# See README_SLURM.md.
# Usage: ./sites_as_jobs.sh [--watch] [--dry-run] [sites]

set -euo pipefail

source "$(dirname "$0")/env.sh"
source "$slurm_dir/job_limit.sh"

watch=0
dry_run=0
while [[ "${1:-}" == --* ]]; do
    case "$1" in
        --watch) watch=1 ;;
        --dry-run) dry_run=1 ;;
        *) echo "unknown option $1" >&2; exit 1 ;;
    esac
    shift
done

log_dir="${SDX_OUTPUT}/logs"
mkdir -p "$log_dir"

# the sites of the run that have a cube, one per line
path_sites="${SDX_OUTPUT}/sites.txt"
list_args=("${1:-$SDX_SITES}" --list)
if [ -n "$SDX_COLLECTIONS" ]; then list_args+=(--collections "$SDX_COLLECTIONS"); fi
if [ -n "$SDX_FREQS" ]; then list_args+=(--freqs "$SDX_FREQS"); fi
(cd "$SDX_REPO" && "$SDX_PYTHON" run_site.py "${list_args[@]}") > "$path_sites"
n_sites=$(grep -c . "$path_sites" || true)
if [ "$n_sites" -eq 0 ]; then
    echo "no sites with a cube in ${SDX_SITES} ${SDX_COLLECTIONS}" >&2
    exit 1
fi

# sites submitted by this run, so a failed site is not resubmitted while watching
submitted=""

# sets n_left to the sites still waiting for a free place
submit_round() {
    local queued_jobs budget_cpus budget_mem used_cpus used_mem site
    queued_jobs=$(squeue -h -u "$USER" -o "%j")
    read -r budget_cpus budget_mem <<< "$(resource_budget)"
    read -r used_cpus used_mem <<< "$(run_usage RUNNING,PENDING)"

    local n_submitted=0 n_done=0 n_queued=0
    n_left=0
    while read -r site; do
        [ -n "$site" ] || continue
        if [ -f "${SDX_OUTPUT}/done/${site}" ] && [ "$SDX_OVERWRITE" != "1" ]; then
            n_done=$((n_done + 1))
            continue
        fi
        if grep -qx "${job_prefix}${site}" <<< "$queued_jobs"; then
            n_queued=$((n_queued + 1))
            submitted+="${site}"$'\n'
            continue
        fi
        if grep -qxF "$site" <<< "$submitted"; then
            continue
        fi
        # one job always runs
        if [ "$used_cpus" -gt 0 ] && { [ $((used_cpus + site_cpus)) -gt "$budget_cpus" ] ||
            [ $((used_mem + site_mem)) -gt "$budget_mem" ]; }; then
            n_left=$((n_left + 1))
            continue
        fi
        local sbatch_args=(--job-name="${job_prefix}${site}" "${sbatch_resources[@]}"
            --output="${log_dir}/${site}_%j.log" --export=ALL
            "$slurm_dir/site.sbatch" "$site")
        if [ "$dry_run" -eq 1 ]; then
            echo "would run: sbatch ${sbatch_args[*]}"
            # SLURM checks the first request and estimates its start, not queuing it
            if [ "$n_submitted" -eq 0 ]; then
                sbatch --test-only "${sbatch_args[@]}" || true
            fi
        else
            sbatch "${sbatch_args[@]}"
        fi
        submitted+="${site}"$'\n'
        n_submitted=$((n_submitted + 1))
        used_cpus=$((used_cpus + site_cpus))
        used_mem=$((used_mem + site_mem))
    done < "$path_sites"

    local verb="submitted"
    [ "$dry_run" -eq 1 ] && verb="would submit"
    echo "$(date '+%F %T') ${verb} ${n_submitted}, ${n_queued} in the queue," \
        "${n_done} done, ${n_left} left to submit, using ${used_cpus} of ${budget_cpus} cpus" \
        "and ${used_mem} of ${budget_mem} MB"
}

if [ "$dry_run" -eq 1 ]; then
    watch=0
    echo "dry run, nothing is submitted"
    echo "  repo:        ${SDX_REPO}"
    echo "  python:      ${SDX_PYTHON}"
    echo "  output:      ${SDX_OUTPUT}"
    echo "  sites:       ${n_sites} of ${1:-$SDX_SITES}," \
        "collections ${SDX_COLLECTIONS:-all}, freqs ${SDX_FREQS:-all}"
    echo "  per job:     ${SDX_THREADS} cpus, ${SDX_MEM} memory, ${SDX_TIME} on ${SDX_PARTITION}"
    echo "  share:       ${SDX_USE_PERCENT}% of the available resources"
fi

submit_round
while [ "$watch" -eq 1 ] && [ "$n_left" -gt 0 ]; do
    sleep $(( SDX_WATCH_MINUTES * 60 ))
    submit_round
done

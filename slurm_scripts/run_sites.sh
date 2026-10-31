#!/bin/bash
# Submits a job per site of the run, refilled by a watcher in a detached tmux
# session. See README_SLURM.md.
# Usage: ./run_sites.sh [--dry-run] [sites]

set -euo pipefail

source "$(dirname "$0")/env.sh"

# show what would be submitted, in the foreground, no tmux
if [ "${1:-}" = "--dry-run" ]; then
    exec "$slurm_dir/sites_as_jobs.sh" "$@"
fi
session="$SDX_RUN_NAME"
path_log="${SDX_OUTPUT}/logs/submit_watch.log"
mkdir -p "$(dirname "$path_log")"

if ! command -v tmux > /dev/null; then
    echo "tmux is not available on $(hostname)" >&2
    exit 1
fi
if tmux has-session -t "$session" 2> /dev/null; then
    echo "a tmux session ${session} already runs on $(hostname)," \
        "attach with: tmux attach -t ${session}" >&2
    exit 1
fi

command=$(printf "%q " cd "$slurm_dir")
command+="&& "
command+=$(printf "%q " ./sites_as_jobs.sh --watch "$@")
command+="2>&1 | tee -a $(printf "%q" "$path_log")"

tmux new-session -d -s "$session" "$command"
tmux set-option -t "$session" remain-on-exit on > /dev/null

echo "watcher started in tmux session ${session} on $(hostname)"
echo "  attach:   tmux attach -t ${session}   (detach with Ctrl-b d)"
echo "  log:      tail -f ${path_log}"
echo "  stop:     tmux kill-session -t ${session}"

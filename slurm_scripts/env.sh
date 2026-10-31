#!/bin/bash
# All settings of the site jobs. Every script in slurm_scripts/ reads this file. See
# README_SLURM.md.

# paths
slurm_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export SDX_REPO="$(dirname "$slurm_dir")"
export SDX_ROOT="/Net/Groups/BGI/work_3/sindbad/site_cubes"
export SDX_PYTHON="${SDX_ROOT}/tools/env/bin/python"

# the run: its name sets the output folder, the job names and the tmux session, so
# runs with different names go side by side without mixing
export SDX_RUN_NAME="site_cubes"
export SDX_SITES="all"                      # all, a list such as DE-Hai,IT-Lav, or @file
export SDX_COLLECTIONS=""                   # such as Shuttle; empty for every collection
export SDX_FREQS=""                         # daily or hourly; empty for both
export SDX_OVERWRITE=0                      # 1 runs sites that are done again, extracted anew

# one job, one site
export SDX_PARTITION="big"
export SDX_THREADS=2                        # cpus, the sites run_site.py runs at once
export SDX_MEM="8G"
export SDX_TIME="1-00:00:00"
export SDX_EXCLUDE=""                       # nodes to avoid, such as node-r5-he34

# submitting: the share of the available resources for all jobs, and the time between
# submit rounds of the watcher, which reads this file at start, so restart the watcher
# after edits
export SDX_USE_PERCENT=75
export SDX_WATCH_MINUTES=10

export SDX_OUTPUT="${SDX_ROOT}/input/store/_runs/${SDX_RUN_NAME}"

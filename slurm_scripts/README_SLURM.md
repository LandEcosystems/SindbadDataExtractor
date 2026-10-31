# Site jobs on SLURM

Runs the workflow of the site cubes on the cluster, one site per job, with
`run_site.py`. A job extracts every store product of the site and builds its cubes, the
cliffin cube that PyCLIFF reads first. It is the per-site alternative to running
Snakemake with `run_workflow.sh`: the watcher keeps a share of the free cpus and
memory busy, and a site that fails does not stop the others.

## Files

All scripts are in `slurm_scripts/`.

| File | Role |
|---|---|
| `env.sh` | **All settings**: paths, sites, collections, job size, share |
| `run_sites.sh` | Entry point: the watcher in a detached tmux session, or `--dry-run` |
| `sites_as_jobs.sh` | One job per site, `--watch` keeps filling free places |
| `job_limit.sh` | Job size and the resource share |
| `site.sbatch` | The job: loads the proxy, runs `run_site.py` for one site, marks it done |

## Before the jobs

The jobs never prefetch. The CLIFF forcing of all sites comes first, a variable per
array task, so that no job reads the ERA5 files for a single site:

```bash
sbatch --array=1-10 prefetch_array.sh cliff_forcing_crujra
sbatch --array=1-9  prefetch_array.sh cliff_forcing_era5
```

The site registry, `insitu/sites_consolidated.csv`, and the other prefetch caches must
be there as well. A site whose forcing is missing still runs, without its PyCLIFF
products, and the job log says so.

## Running

1. **Pick the sites** in `env.sh`: `SDX_SITES` (`all`, a list or `@file`),
   `SDX_COLLECTIONS` (such as `Shuttle`) and `SDX_FREQS`. Only sites with a cube of the
   collections get a job.

2. Dry run, in the foreground, to check the setup without submitting anything:

   ```bash
   ./slurm_scripts/run_sites.sh --dry-run          # optionally with sites
   ```

   It lists the sites into `$SDX_OUTPUT/sites.txt`, prints the paths, the job size and
   the share, the resource budget and the `sbatch` command of each job it would
   submit. The first request goes to `sbatch --test-only`, so SLURM validates it and
   estimates a start time without queuing it.

3. Start it on the login node:

   ```bash
   ./slurm_scripts/run_sites.sh                   # the sites of env.sh
   ./slurm_scripts/run_sites.sh DE-Hai,IT-Lav     # or others
   ```

   It prints the host name. The tmux session lives only on that login node.

| Action | Command |
|---|---|
| Follow progress | `tail -f $SDX_OUTPUT/logs/submit_watch.log` |
| A site's log | `$SDX_OUTPUT/logs/<site>_<jobid>.log` |
| Attach / detach | `tmux attach -t $SDX_RUN_NAME`, then `Ctrl-b d` |
| Stop submitting | `tmux kill-session -t $SDX_RUN_NAME` (running jobs continue) |
| See jobs | `squeue -u $USER` |

## How it works

**Job size.** Every job has `SDX_THREADS` cpus and `SDX_MEM`. `run_site.py` runs that
many sites at once within a step, so with one site per job the cpus mostly serve the
cube step; the PyCLIFF fit of a site runs in one process.

**Resource share.** Each round, `job_limit.sh` adds the idle cpus and unallocated
memory of the partition (only nodes in `idle` or `mix` state) to what the running jobs
of the run hold, and takes `SDX_USE_PERCENT` % of both as the budget. Jobs are
submitted while the cpus and memory of the running and pending jobs of the run stay
within it. One job is always allowed.

**Watching.** Every `SDX_WATCH_MINUTES` the submitter fills free places and exits once
every site is done or has been submitted. Skipped:

- sites with a marker in `$SDX_OUTPUT/done/`, unless `SDX_OVERWRITE=1`,
- sites with a job `<SDX_RUN_NAME>_<site>` of yours already in the queue,
- sites already submitted by this watch, so a failing site is not retried in a loop.

To retry failed sites, start the watcher again after it has finished.

**Runs.** `SDX_RUN_NAME` names the output folder, the job names and the tmux session, so
runs with different names go side by side without mixing. Jobs of different sites
share the store safely: `extract_product.py` merges the manifest of a product under a
lock.

**After the jobs.** Nothing here writes the markers of Snakemake. A later
`run_workflow.sh` therefore still extracts and assembles everything, skips the
sites whose stage hash is unchanged, and writes the zarr stores, which combine all
sites of a cube.

## Settings

All settings are in `slurm_scripts/env.sh`, with a comment per setting, and every
script reads it, replacing any of these variables set in the shell. The sbatch file
holds no resources; the submit scripts pass partition, cpus, memory and time.

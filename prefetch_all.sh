#!/bin/bash
# Submit the prefetch of every dataset with a cache, for the sites missing from it.
#
#   ./prefetch_all.sh                     every dataset below
#   ./prefetch_all.sh fa_gami vf_s2       only these
#   DRY=1 ./prefetch_all.sh               print the sbatch commands only
#
# Each dataset gets one job of prefetch_array.sh, or an array whose tasks split its
# work. Sites already in a cache are skipped, so the script can simply be run again
# after a failure. jung_C2F waits for shuttle_obs, whose zips it reads. Left out:
# era5_API, paused for its request limits; vf_modis, whose AppEEARS request runs in
# the steps of prefetch.py --step; shuttle_clim, which shares the cache of
# shuttle_obs. The logs are prefetch_<dataset>_<job>[_<task>].log in
# output/slurm_logs of the site cubes. The CLMS and Landsat tasks run 8 threads on 2
# cpus, since their threads wait on the network, which leaves room under the cpu
# limit per user.
cd "$(dirname "$0")"
logs=/Net/Groups/BGI/work_3/sindbad/site_cubes/output/slurm_logs
[ -n "${DRY:-}" ] || mkdir -p "$logs"

# dataset, then the sbatch options of its job
jobs=(
    "cliff_forcing_era5    --array=1-9 --cpus-per-task=8 --mem=32G"
    "cliff_forcing_crujra  --array=1-10 --cpus-per-task=1 --mem=8G"
    "fa_gami               --cpus-per-task=8 --mem=16G"
    "agb_esacci            --cpus-per-task=8 --mem=16G"
    "pft_cgls              --cpus-per-task=8 --mem=16G"
    "pft_worldcover        --cpus-per-task=8 --mem=16G"
    "vf_esacci             --cpus-per-task=8 --mem=16G"
    "vf_cgls               --cpus-per-task=8 --mem=16G"
    "vf_worldcover         --cpus-per-task=8 --mem=16G"
    "vf_fcover             --cpus-per-task=8 --mem=16G"
    "wf_chm                --cpus-per-task=8 --mem=16G"
    "vf_s2                 --cpus-per-task=8 --mem=16G"
    "vf_s2sl2p             --cpus-per-task=8 --mem=16G"
    "ndvi_clms300          --array=1-13 --cpus-per-task=2 --export=ALL,WORKERS=8"
    "ndvi_clms1km          --array=1-22 --cpus-per-task=2 --export=ALL,WORKERS=8"
    "lai_clms300           --array=1-13 --cpus-per-task=2 --export=ALL,WORKERS=8"
    "lai_clms1km           --array=1-15 --cpus-per-task=2 --export=ALL,WORKERS=8"
    "ndvi_landsat          --array=1-40 --cpus-per-task=2 --export=ALL,WORKERS=8"
    "shuttle_obs           --cpus-per-task=4 --mem=16G"
    "jung_C2F              --cpus-per-task=8 --mem=32G"
)

shuttle_id=""
for line in "${jobs[@]}"; do
    read -r name opts <<< "$line"
    if [ $# -gt 0 ] && [[ ! " $* " =~ " $name " ]]; then
        continue
    fi
    after=()
    # %a of a job outside an array is no number, so such a job is named by its id
    if [[ "$opts" != *--array* ]]; then
        after=(--output="$logs/%x_%j.log")
    fi
    if [ "$name" = jung_C2F ] && [ -n "$shuttle_id" ]; then
        after+=(--dependency="afterany:$shuttle_id")
    fi
    # shellcheck disable=SC2086
    cmd=(sbatch --parsable --job-name="prefetch_$name" $opts "${after[@]}"
         prefetch_array.sh "$name")
    if [ -n "${DRY:-}" ]; then
        echo "${cmd[*]}"
        id=DRY
    else
        id=$("${cmd[@]}") || exit 1
        echo "$name: job $id"
    fi
    if [ "$name" = shuttle_obs ]; then
        shuttle_id=$id
    fi
done

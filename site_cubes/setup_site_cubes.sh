#!/bin/bash
# Folders of the SINDBAD site cubes, described in README.md next to this script.
#
#   ./setup_site_cubes.sh dirs     make the folders
#   ./setup_site_cubes.sh tools    clone SindbadDataExtractor, PyCLIFF and fluxcom
#                                  into tools, or show their last commit
#   ./setup_site_cubes.sh legacy   merge the LaThuile2007 site files of the CLIFF
#                                  input and of the fluxnetBGI v0.07 extraction
#                                  into ancillary_data/legacy/LaThuile2007 with
#                                  merge_lathuile.py, with SOURCE.md and MD5SUMS,
#                                  read only
#   ./setup_site_cubes.sh check    list what the pipeline reads, with what exists
#
# The root is the folder of this script, or ROOT when it is set. EXTRACTOR_URL,
# PYCLIFF_URL and FLUXCOM_URL override where the repositories are cloned from.
# PyCLIFF comes by default from the remote of the checkout in toolDev.
set -euo pipefail

root=${ROOT:-$(cd "$(dirname "$0")" && pwd)}
anc=$root/ancillary_data
legacy_ltl07=/Net/Groups/BGI/work_3/cliff/input/fluxnetBGI2021_siteCube-v0.03_iter-5
legacy_v007=/Net/Groups/BGI/scratch/skoirala/fluxnetBGI_v0.07
github=https://github.com/LandEcosystems
extractor_url=${EXTRACTOR_URL:-$github/SindbadDataExtractor.git}
pycliff_url=${PYCLIFF_URL:-}
fluxcom_url=${FLUXCOM_URL:-https://git.bgc-jena.mpg.de/fluxcom/fluxcom}
pycliff_dev=/Net/Groups/BGI/scratch/skoirala/toolDev/PyCLIFF

make_dirs() {
    mkdir -p "$anc"/insitu/site_lists "$anc"/reference "$anc"/prefetch \
        "$anc"/legacy "$root"/input/store "$root"/input/pycliff "$root"/output \
        "$root"/tools
    echo "folders ready below $root"
}

clone_tools() {
    mkdir -p "$root"/tools
    if [ -z "$pycliff_url" ]; then
        pycliff_url=$(git -C "$pycliff_dev" remote get-url origin)
    fi
    for repo in "SindbadDataExtractor $extractor_url" "PyCLIFF $pycliff_url" \
        "fluxcom $fluxcom_url"; do
        set -- $repo
        if [ -d "$root/tools/$1/.git" ]; then
            echo "$1: $(git -C "$root/tools/$1" log -1 --format='%h %ad %s' \
                --date=short)"
        else
            git clone "$2" "$root/tools/$1"
        fi
    done
}

# The LaThuile2007 site files of the CLIFF input, with the meteorology, and of the
# fluxnetBGI v0.07 extraction, with the fluxes, merged into one read-only copy.
merge_legacy() {
    if [ -e "$anc/legacy/LaThuile2007/SOURCE.md" ]; then
        echo "$anc/legacy/LaThuile2007 already holds a merged copy, see its SOURCE.md"
        return
    fi
    "${PYTHON:-$root/tools/env/bin/python}" "$(dirname "$0")/merge_lathuile.py" \
        "$anc/legacy" --ltl07 "$legacy_ltl07" --v007 "$legacy_v007" \
        --workers "${WORKERS:-4}"
}

# Kind and path, relative to the root, of what the pipeline reads. copy is copied by
# hand, clone is a repository of tools and made is written by a tool of
# SindbadDataExtractor or PyCLIFF.
expected() {
    cat <<EOF
clone tools/SindbadDataExtractor
clone tools/PyCLIFF
clone tools/fluxcom
made tools/env/bin/python
copy ancillary_data/all_sites_multi_with_ru_flux_and_oneflux.csv
copy ancillary_data/reference/ne_10m_admin_0_map_units.zip
copy ancillary_data/insitu/PFT_consolidated.csv
copy ancillary_data/insitu/PFT.csv
copy ancillary_data/insitu/PFT_overrides.csv
copy ancillary_data/insitu/site_coords_overrides.csv
copy ancillary_data/insitu/FLUXNET_veg_frac.csv
copy ancillary_data/insitu/agb_De2026_consolidated.csv
copy ancillary_data/insitu/forest_age.csv
copy ancillary_data/insitu/globbiomass_agb_100m.csv
copy ancillary_data/insitu/soilDepth_Shan2026_consolidated.csv
copy ancillary_data/insitu/treeFrac_De2026_consolidated.csv
copy ancillary_data/insitu/vegFrac_De2026_consolidated.csv
copy ancillary_data/prefetch/fluxnet_shuttle/index_fluxnet_shuttle.csv
copy ancillary_data/prefetch/jung_C2F/index_jung_C2F_loose.csv
copy ancillary_data/prefetch/jung_C2F/code/c2f_fluxnet-v1.0
copy ancillary_data/prefetch/cgls_lc100_v3_1km_sites.csv
copy ancillary_data/prefetch/esacci_agb_v7_1km_sites.csv
copy ancillary_data/prefetch/esacci_lc_pft_v2.0.81_1km_sites.csv
copy ancillary_data/prefetch/fcover300_2015/site_means.csv
copy ancillary_data/prefetch/meta_chm_1km_sites.csv
copy ancillary_data/prefetch/meta_chm_tiles.geojson
copy ancillary_data/prefetch/modis_vcf_2010/site_means.csv
copy ancillary_data/prefetch/pft_cgls_lc100_v3_100m_sites.csv
copy ancillary_data/prefetch/pft_worldcover_v200_10m_sites.csv
copy ancillary_data/prefetch/s2_fvc_2019_1km_sites_v2.csv
copy ancillary_data/prefetch/s2_peak_ndvi_2019
copy ancillary_data/prefetch/s2_sl2p_fcover_2019_1km_sites.csv
copy ancillary_data/prefetch/s2_sl2p_fcover_2019_peak
copy ancillary_data/prefetch/worldcover_2021_v200_1km_sites.csv
made ancillary_data/legacy/LaThuile2007/SOURCE.md
made ancillary_data/prefetch/gami_v2.1_1km_sites.csv
made ancillary_data/prefetch/clms_ndvi_300m
made ancillary_data/prefetch/clms_ndvi_1km
made ancillary_data/prefetch/clms_lai_300m
made ancillary_data/prefetch/clms_lai_1km
made ancillary_data/prefetch/landsat_ndvi
made ancillary_data/prefetch/cliff_forcing/ERA5.e1
made ancillary_data/prefetch/cliff_forcing/CRUJRA.v4.10
made ancillary_data/insitu/site_lists
made ancillary_data/insitu/sites_consolidated.csv
made input/store
made input/pycliff
made output
EOF
}

check() {
    n_missing=0
    while read -r kind path; do
        if [ -e "$root/$path" ]; then
            state=ok
        else
            state=missing
            [ "$kind" != made ] && n_missing=$((n_missing + 1))
        fi
        printf "%-8s %-5s %s\n" "$state" "$kind" "$path"
    done < <(expected)
    echo "$n_missing files to copy or clone are missing below $root"
}

case "${1:-}" in
    dirs) make_dirs ;;
    tools) clone_tools ;;
    legacy) merge_legacy ;;
    check) check ;;
    *) sed -n '2,11p' "$0" | sed 's/^# \{0,1\}//'; exit 1 ;;
esac

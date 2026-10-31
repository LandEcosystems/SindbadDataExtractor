#!/bin/bash
# Folders of the SINDBAD site cubes, described in README.md next to this script.
#
#   ./setup_site_cubes.sh dirs     make the folders
#   ./setup_site_cubes.sh tools    clone SindbadDataExtractor, PyCLIFF and fluxcom
#                                  into tools, or show their last commit
#   ./setup_site_cubes.sh legacy   copy the LaThuile2007 site files of the CLIFF
#                                  input into ancillary_data/legacy, write
#                                  SOURCE.md and MD5SUMS and make it read only
#   ./setup_site_cubes.sh check    list what the pipeline reads, with what exists
#
# The root is the folder of this script, or ROOT when it is set. EXTRACTOR_URL,
# PYCLIFF_URL and FLUXCOM_URL override where the repositories are cloned from.
# PyCLIFF comes by default from the remote of the checkout in toolDev.
set -euo pipefail

root=${ROOT:-$(cd "$(dirname "$0")" && pwd)}
anc=$root/ancillary_data
legacy_src=/Net/Groups/BGI/work_3/cliff/input/fluxnetBGI2021_siteCube-v0.03_iter-5
legacy_dst=$anc/legacy/fluxnetBGI2021_siteCube-v0.03_iter-5
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

copy_legacy() {
    if [ -e "$legacy_dst/SOURCE.md" ]; then
        echo "$legacy_dst already holds a copy, see its SOURCE.md" >&2
        exit 1
    fi
    for res in DD HR; do
        for sub in data config_info; do
            src=$legacy_src/fluxnetBGI2021.LTL07.$res/$sub
            if [ -d "$src" ]; then
                mkdir -p "$legacy_dst/fluxnetBGI2021.LTL07.$res"
                rsync -a "$src" "$legacy_dst/fluxnetBGI2021.LTL07.$res/"
            else
                echo "not found: $src"
            fi
        done
    done
    (cd "$legacy_dst" && find . -type f -name '*.nc' | sort | xargs md5sum) \
        > "$legacy_dst/MD5SUMS"
    n_dd=$(find "$legacy_dst"/fluxnetBGI2021.LTL07.DD -name '*.nc' | wc -l)
    n_hr=$(find "$legacy_dst"/fluxnetBGI2021.LTL07.HR -name '*.nc' | wc -l)
    cat > "$legacy_dst/SOURCE.md" <<EOF
# LaThuile2007 site files

Read-only copy of the LTL07 site files of the CLIFF input cube, used as the
LaThuile2007 source of the site cubes. Never edit these files.

- Origin: $legacy_src/fluxnetBGI2021.LTL07.{DD,HR}
- Copied: $(date -u +%Y-%m-%dT%H:%M:%SZ) by $(whoami)
- Files: $n_dd daily, $n_hr hourly
- Checksums: MD5SUMS
EOF
    chmod -R a-w "$legacy_dst"
    echo "copied $n_dd daily and $n_hr hourly files to $legacy_dst"
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
copy ancillary_data/legacy/fluxnetBGI2021_siteCube-v0.03_iter-5/SOURCE.md
made ancillary_data/prefetch/gami_v2.1_1km_sites.csv
made ancillary_data/prefetch/clms_ndvi_300m
made ancillary_data/prefetch/clms_ndvi_1km
made ancillary_data/prefetch/clms_lai_300m
made ancillary_data/prefetch/clms_lai_1km
made ancillary_data/prefetch/landsat_ndvi
made ancillary_data/prefetch/cliff_forcing/ERA5.e1
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
    legacy) copy_legacy ;;
    check) check ;;
    *) sed -n '2,11p' "$0" | sed 's/^# \{0,1\}//'; exit 1 ;;
esac

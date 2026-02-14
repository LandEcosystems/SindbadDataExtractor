#!/bin/bash
cube_version='v0.07'
run_name='fluxnetBGI'
output_path="/Net/Groups/BGI/work_3/cliff/input/${run_name}_${cube_version}/"
pythonpath='/Net/Groups/BGI/scratch/skoirala/sksrc/miniforge3/envs/sinproc007/bin/python'
configuration_file='/Net/Groups/BGI/scratch/skoirala/SindbadSiteExtractor/configuration/cliff/cliffin_exp_config.json'

fluxcom_cube_path_base='/Net/Groups/BGI/work_4/scratch/fluxcom/sitecube_proc/'
fn_version='FLUXNET2015'

# fn_version='LaThuile'
time_scale='daily'
declare -a time_scales=("daily" "hourly")
declare -a time_scales=("hourly")

## declare an array variable
declare -a fn_versions=("LaThuile")
declare -a fn_versions=("AmeriFlux" "FLUXNET2015" "ICOS" "ICOSDrought2018" "OzFLUX" "WarmWinter2020" "CCBY")
declare -a fn_versions=("FLUXNET2015" "ICOS" "ICOSDrought2018" "OzFLUX" "WarmWinter2020" "CCBY")
## now loop through the above array
for time_scale in "${time_scales[@]}"
do
    for fn_version in "${fn_versions[@]}"
    do
        if [ "$fn_version" == "LaThuile" ]; then
            fluxcom_cube_path="/Net/Groups/BGI/data/FLUXCOM/v0.04/site_cube/"
            # when running for LaThuile, we need to use the older cube version (v0.04) since the newer cube versions do not have the LaThuile data. Also, need to use older version of fluxcom. Copy fluxcom.v6 to fluxcom in sinproc_installer directory and do that before running the extraction for LaThuile.
        else
            fluxcom_cube_path="${fluxcom_cube_path_base}${cube_version}/site_cube/"
        fi
        # echo ${pythonpath} ${PWD}/run_extraction.py ${configuration_file} ${fn_version} ${time_scale} ${output_path} ${fluxcom_cube_path} ${run_name}
        ${pythonpath} ${PWD}/run_extraction.py ${configuration_file} ${fn_version} ${time_scale} ${output_path} ${fluxcom_cube_path} ${run_name}
    done
done
exit 0

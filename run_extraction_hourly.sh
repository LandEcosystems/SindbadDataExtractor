#!/bin/bash
output_path='/Net/Groups/BGI/scratch/skoirala/fluxNet_0.07/'
pythonpath='/Net/Groups/BGI/scratch/skoirala/sksrc/miniforge3/envs/sinproc007/bin/python'
configuration_file='/Net/Groups/BGI/scratch/skoirala/SindbadSiteExtractor/configuration/hourly/hourly_extraction_main.json'
# fn_version='FLUXNET2015'
fn_version='LaThuile'
time_scale='hourly'

## declare an array variable
declare -a fn_versions=("LaThuile" "WarmWinter2020" "AmerifluxFLUXNET" "ICOSFLUXNET" "ICOSDrought2018" "FLUXNET2015" "CCBY" "CCBY_v0.04" "CCBY_v0.05")

## now loop through the above array
for fn_version in "${fn_versions[@]}"
do
    echo ${pythonpath} ${PWD}/run_extraction.py ${configuration_file} ${fn_version} ${time_scale} ${output_path}
    ${pythonpath} ${PWD}/run_extraction.py ${configuration_file} ${fn_version} ${time_scale} ${output_path}
done

# You can access them using echo "${arr[0]}", "${arr[1]}" also
exit 0
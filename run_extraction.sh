#!/bin/bash
output_path='/Net/Groups/BGI/scratch/skoirala/fluxNet_0.07/'
pythonpath='/Net/Groups/BGI/scratch/skoirala/sksrc/miniforge3/envs/sinproc007/bin/python'
configuration_file='/Net/Groups/BGI/scratch/skoirala/SindbadSiteExtractor/configuration/sindbad/sindbad_extraction_main.json'
fn_version='FLUXNET2015'
# fn_version='LaThuile'
time_scale='daily'

## declare an array variable
# declare -a fn_versions=("LaThuile")
# declare -a fn_versions=("WarmWinter2020" "AmeriFlux" "ICOS" "ICOSDrought2018" "FLUXNET2015" "CCBY" "CCBY_v0.04" "CCBY_v0.05")
declare -a fn_versions=("AmeriFlux" "FLUXNET2015" "ICOS" "ICOSDrought2018" "OzFLUX" "WarmWinter2020" "CCBY")
## now loop through the above array
for fn_version in "${fn_versions[@]}"
do
    echo ${pythonpath} ${PWD}/run_extraction.py ${configuration_file} ${fn_version} ${time_scale} ${output_path}
    ${pythonpath} ${PWD}/run_extraction.py ${configuration_file} ${fn_version} ${time_scale} ${output_path}
done
exit 0


# /Net/Groups/BGI/scratch/skoirala/sksrc/miniforge3/envs/sinproc007/bin/python /Net/Groups/BGI/scratch/skoirala/SindbadSiteExtractor/run_extraction.py /Net/Groups/BGI/scratch/skoirala/SindbadSiteExtractor/configuration/sindbad/sindbad_extraction_main.json AmeriFlux daily /Net/Groups/BGI/scratch/skoirala/fluxNet_0.07/
# /Net/Groups/BGI/scratch/skoirala/sksrc/miniforge3/envs/sinproc007/bin/python /Net/Groups/BGI/scratch/skoirala/SindbadSiteExtractor/run_extraction.py /Net/Groups/BGI/scratch/skoirala/SindbadSiteExtractor/configuration/sindbad/sindbad_extraction_main.json FLUXNET2015 daily /Net/Groups/BGI/scratch/skoirala/fluxNet_0.07/
# /Net/Groups/BGI/scratch/skoirala/sksrc/miniforge3/envs/sinproc007/bin/python /Net/Groups/BGI/scratch/skoirala/SindbadSiteExtractor/run_extraction.py /Net/Groups/BGI/scratch/skoirala/SindbadSiteExtractor/configuration/sindbad/sindbad_extraction_main.json ICOS daily /Net/Groups/BGI/scratch/skoirala/fluxNet_0.07/
# /Net/Groups/BGI/scratch/skoirala/sksrc/miniforge3/envs/sinproc007/bin/python /Net/Groups/BGI/scratch/skoirala/SindbadSiteExtractor/run_extraction.py /Net/Groups/BGI/scratch/skoirala/SindbadSiteExtractor/configuration/sindbad/sindbad_extraction_main.json ICOSDrought2018 daily /Net/Groups/BGI/scratch/skoirala/fluxNet_0.07/
# /Net/Groups/BGI/scratch/skoirala/sksrc/miniforge3/envs/sinproc007/bin/python /Net/Groups/BGI/scratch/skoirala/SindbadSiteExtractor/run_extraction.py /Net/Groups/BGI/scratch/skoirala/SindbadSiteExtractor/configuration/sindbad/sindbad_extraction_main.json OzFLUX daily /Net/Groups/BGI/scratch/skoirala/fluxNet_0.07/
# /Net/Groups/BGI/scratch/skoirala/sksrc/miniforge3/envs/sinproc007/bin/python /Net/Groups/BGI/scratch/skoirala/SindbadSiteExtractor/run_extraction.py /Net/Groups/BGI/scratch/skoirala/SindbadSiteExtractor/configuration/sindbad/sindbad_extraction_main.json WarmWinter2020 daily /Net/Groups/BGI/scratch/skoirala/fluxNet_0.07/
# /Net/Groups/BGI/scratch/skoirala/sksrc/miniforge3/envs/sinproc007/bin/python /Net/Groups/BGI/scratch/skoirala/SindbadSiteExtractor/run_extraction.py /Net/Groups/BGI/scratch/skoirala/SindbadSiteExtractor/configuration/sindbad/sindbad_extraction_main.json CCBY_v0.07 daily /Net/Groups/BGI/scratch/skoirala/fluxNet_0.07/

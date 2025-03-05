#!/bin/bash
output_path='/Net/Groups/BGI/scratch/skoirala/fluxNet_0.06/'
pythonpath='/Net/Groups/BGI/scratch/skoirala/sksrc/miniforge3/envs/sinproc006/bin/python'
configuration_file='/Net/Groups/BGI/scratch/skoirala/SindbadSiteExtractor/configuration/sindbad/sindbad_extraction_main.json'
# fn_version='FLUXNET2015'
fn_version='LaThuile'
time_scale='daily'
echo ${pythonpath} ${PWD}/run_extraction.py ${configuration_file} ${fn_version} ${time_scale} ${output_path}
${pythonpath} ${PWD}/run_extraction.py ${configuration_file} ${fn_version} ${time_scale} ${output_path}
exit 0
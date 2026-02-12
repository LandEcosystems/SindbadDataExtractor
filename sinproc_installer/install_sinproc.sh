## Install conda environment
# Before running the installer, download the latest version of the fluxcom git repository and save it as a fluxcom folder here
# available at: https://git.bgc-jena.mpg.de/fluxcom/fluxcom
conda env create -f sinproc.yml
conda activate sinproc007
pip install -e fluxcom


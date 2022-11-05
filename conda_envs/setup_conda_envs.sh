#!/usr/bin/env bash

set -e

# NOTE: Please replace `conda` with `mamba` if it is installed for faster installs.
resolverCondaBinary="mamba" # OR mamba

# NOTE: By default, the conda environments are expected by the `conda_local` profile to be created within `xbs-nf/conda_envs` directory

$resolverCondaBinary env create -p xbs-nf-env-1 --file xbs-nf-env-1.yml 

#$resolverCondaBinary env create -p xbs-nf-env-2 --file xbs-nf-env-2.yml

echo "INFO: Activate conda env with tb-profiler and setup the WHO database within the xbs-nf-env-1"
eval "$(conda shell.bash hook)"
conda activate "./xbs-nf-env-1"

echo "INFO: Make a local copy and cd inside it"
cp -r ../resources/resistance_db_who ./
cd resistance_db_who

echo "INFO: Load the database within tb-profiler"
tb-profiler load_library resistance_db_who

echo "INFO: Remove the local copy of the database"
cd ..
rm -rf resistance_db_who

echo "INFO: Deactivate the xbs-nf-env-1 env"
conda deactivate

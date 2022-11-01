#!/usr/bin/env bash

set -e

# NOTE: Please replace `conda` with `mamba` if it is installed for faster installs.
condaBinary="conda" # OR mamba

# NOTE: By default, the conda environments are expected by the `conda_local` profile to be created within `xbs-nf/conda_envs` directory

$condaBinary env create -p xbs-nf-env-1 --file xbs-nf-env-1.yml

$condaBinary env create -p xbs-nf-env-2 --file xbs-nf-env-2.yml

#NOTE: Activate conda env with tb-profiler
eval "$(conda shell.bash hook)"
$condaBinary activate "./xbs-nf-env-1"

#NOTE: Setup the WHO database
cp -r ../resources/resistance_db_who ./
cd resistance_db_who
tb-profiler load_library resistance_db_who
rm -rf resistance_db_who ./
cd ..
$condaBinary deactivate

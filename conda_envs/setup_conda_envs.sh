#!/usr/bin/env bash

set -e

# NOTE: Please replace `conda` with `mamba` if it is installed for faster installs.
resolverCondaBinary="conda" # pick either conda OR mamba

# NOTE: By default, the conda environments are expected by the `conda_local` profile to be created within `magma/conda_envs` directory

$resolverCondaBinary env create -p magma-env-1 --file magma-env-1.yml 

$resolverCondaBinary env create -p magma-env-2 --file magma-env-2.yml

#NOTE: Setup the tbprofiler env with WHO v2 Database

$resolverCondaBinary env create -p magma-tbprofiler-env --file magma-tbprofiler-env.yml

echo "INFO: Activate conda env with tb-profiler and setup the WHO database"
eval "$(conda shell.bash hook)"
conda activate "./magma-tbprofiler-env"

echo "INFO: Use WHO-v2 database in tb-profiler"
tb-profiler update_tbdb --commit 22e9ceb6433d02071222af8dca506c6348604877 --logging DEBUG

echo "INFO: Deactivate the magma-tbprofiler-env "
conda deactivate

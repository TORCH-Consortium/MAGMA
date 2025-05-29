#!/usr/bin/env bash

set -e

# NOTE: Please replace `conda` with `mamba` if it is installed for faster installs.
resolverCondaBinary="conda" # pick either conda OR mamba


#===========================================================
#
# NOTE: By default, the conda environments are expected by the `conda_local` profile to be created within `magma/conda_envs` directory

$resolverCondaBinary env create -p magma-env-1 --file magma-env-1.yml 

$resolverCondaBinary env create -p magma-env-2 --file magma-env-2.yml 

$resolverCondaBinary env create -p magma-ntmprofiler-env --file magma-ntmprofiler-env.yml

$resolverCondaBinary env create -p magma-tbprofiler-env --file magma-tbprofiler-env.yml

#===========================================================

#NOTE: Setup the tbprofiler env with WHO v2 Database

$resolverCondaBinary env create -p magma-tbprofiler-env --file magma-tbprofiler-env.yml

echo "INFO: Activate conda env with tb-profiler and setup the WHO database"
eval "$(conda shell.bash hook)"
conda activate "./magma-tbprofiler-env"

#echo "INFO: Use WHO-v2 database in magma-tbprofiler-env"
#tb-profiler update_tbdb --commit bdace1f82d948ce0001e1dade6eb93d2da9c47e5 --logging DEBUG

#echo "INFO: Use MAGMA branch from tbdb database in magma-tbprofiler-env"
tb-profiler update_tbdb --commit 30f8bc37df15affa378ebbfbd3e1eb4c5903056e --logging DEBUG


echo "INFO: Deactivate the magma-tbprofiler-env "
conda deactivate

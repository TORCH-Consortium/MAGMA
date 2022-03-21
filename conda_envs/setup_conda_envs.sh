#!/usr/bin/env bash

set -xue

# NOTE: Please replace `conda` with `mamba` if it is installed for faster installs.

# NOTE: The conda environments are expected by the `conda_local` profile to be created within `xbs-nf/conda_envs` directory

conda env create -p xbs-nf-env-1 --file xbs-nf-env-1.yml

conda env create -p xbs-nf-env-2 --file xbs-nf-env-2.yml

## Conda based execution

You can run the MAGMA pipeline using the Conda or (Micro)mamba package manager to install all the prerequisite softwares from popular repositories such as bioconda and conda-forge.

You can use the `conda` based setup for the pipeline for running MAGMA 
- On a local linux machine(e.g. your laptop or a university server)
- On an HPC cluster (e.g. SLURM, PBS) in case you don't have access to container systems like Singularity, Podman or Docker 


All the requisite softwares have been provided as a `conda` recipe (i.e. `yml` files) 
- [magma-env-1.yml](./conda_envs/magma-env-1.yml)
- [magma-env-2.yml](./conda_envs/magma-env-2.yml)

These files can be downloaded using the following commands

```console
wget https://raw.githubusercontent.com/TORCH-Consortium/MAGMA/master/conda_envs/magma-env-2.yml
wget https://raw.githubusercontent.com/TORCH-Consortium/MAGMA/master/conda_envs/magma-env-1.yml

```

The `conda` environments are expected by the `conda_local` profile of the pipeline, it is recommended that it should be created **prior** to the use of the pipeline, using the following commands. Note that if you have `mamba` (or `micromamba`) available you can rely upon those as well.


```sh
$ conda env create -n magma-env-1 --file magma-env-1.yml

$ conda env create -n magma-env-2 --file magma-env-2.yml
```

Once the environments are created, you can make use of the `conda_envs_location` parameter to inform the pipeline of the names and location of the conda envs.

TIP: You can find out the location of conda environments using `conda env list`. [Here's](https://docs.conda.io/projects/conda/en/4.6.0/_downloads/52a95608c49671267e40c689e0bc00ca/conda-cheatsheet.pdf) a useful cheatsheet for conda operations.

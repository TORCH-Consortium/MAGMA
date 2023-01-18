## Conda based execution

You can run the MAGMA pipeline using the Conda based package manager to install all the prerequisite softwares.

The `conda` environments are expected by the `conda_local` profile of the pipeline, to be created within `MAGMA/conda_envs` directory

> **NOTE**
> If you do have access to Singularity or Podman, then owing to their compatibility with Docker, you can still use the MAGMA Docker containers mentioned [docker.config](../conf/docker.config).


You can use the `conda` based setup for the pipeline for running MAGMA 
- On a local linux machine (e.g. your laptop or university server)
- On an HPC cluster in case you don't have access to container systems like Singularity, Podman or Docker


### Steps to setup the pipeline locally

> **NOTE**
> These steps are only necessary if you don't have access to any container system, then therefore you'd need to install all softwares using the `conda` package manager.


1. Copy the environment files from [conda_envs](../conda_envs) folder locally

```sh
$ git clone https://github.com/TORCH-Consortium/MAGMA

$ cd MAGMA

```

2. After `cd` in the `conda_envs` folder and execute the following commands to create the env 
 
> **TIP**
> 1. For faster installation process, please download [mamba](https://github.com/mamba-org/mamba) tool and replace `conda` with `mamba` in the above commands.
> 2. The path `-p` should be customized as per you setup


```sh
$ conda env create -p magma-env-1 --file magma-env-1.yml

$ conda env create -p magma-env-2 --file magma-env-2.yml
```


### Run the pipeline

1. Customize the pipeline and process level settings in the [default_params](../default_params.config) file

2. From inside the `MAGMA` folder, invoke the pipeline

```sh
$ nextflow run main.nf -profile conda_local
```
3. Use the `-resume` flag to continue from previously generated output files, rather than restarting the entire analysis.

```sh
$ nextflow run main.nf -profile conda_local -resume
```

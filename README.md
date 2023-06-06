# MAGMA

MAGMA (**M**aximum **A**ccessible **G**enome for **M**tb **A**nalysis) is a pipeline for comprehensive genomic analyses of Mycobacterium tuberculosis with a focus on clinical decision making as well as research.

# Salient features of the implementation

- Fine-grained control over resource allocation (CPU/Memory/Storage)
- Reliance of bioconda for installing packages for reproducibility
- Ease of use on a range of infrastructure (cloud/on-prem HPC clusters/ servers (or local machines))
- Resumability for failed processes
- Centralized locations for specifying analysis parameters and hardware requirements
  - MAGMA parameters (`default_parameters.config`)
  - Hardware requirements (`conf/standard.config`)
  - Execution (software) requirements (`conf/docker.config` or `conf/conda.config`)
- An (optional) GVCF reference dataset for ~600 samples is provided for augmenting smaller datasets

> **Note**
Downloading the reference EXIT_RIF GVCF files from FIXME

## Tutorials and Presentations

For the tutorials(./docs/tutorials.md) and [presentations](./docs/presentations.md) please refer the [docs](./docs) folder.

## Prerequisites

### Nextflow

- `git` : The version control in the pipeline.
- `Java-11` or `Java-17` (preferred)

> :warning: **Check `java` version!**:
The `java` version should NOT be an `internal jdk` release! You can check the release via `java -version`

- Download Nextflow

```bash
$ curl -s https://get.nextflow.io | bash
```

- Make Nextflow executable

```sh
$ chmod +x nextflow
```

- Add `nextflow` to your `path` (for example `/usr/local/bin/`)

```sh
$ mv nextflow /usr/local/bin

```

- Sanity check for `nextflow` installation

```console
$ nextflow info

  Version: 23.04.1 build 5866
  Created: 15-04-2023 06:51 UTC (08:51 SAST)
  System: Mac OS X 12.6.5
  Runtime: Groovy 3.0.16 on OpenJDK 64-Bit Server VM 17.0.7+7-LTS
  Encoding: UTF-8 (UTF-8)

```

> :heavy_check_mark: **With this you're all set with Nextflow. Next stop, conda or docker - pick one!**: <br>

## Customizing pipeline parameters for your dataset

The pipeline parameters are distinct from Nextflow parameters, and therefore it is recommended that they are provided using a `yml` file as shown below


```yml

# Sample contents of my_parameters_1.yml file

input_samplesheet: /path/to/your_samplesheet.csv
only_validate_fastqs: true
conda_envs_location: /path/to/both/conda_envs
```


> **Note**
The `-profile` mechanism is used to enable infrastructure specific settings of the pipeline. The example below, assumes you are using `conda` based setup.


Which could be provided to the pipeline using `-params-file` parameter as shown below

```console
nextflow run 'https://github.com/TORCH-Consortium/MAGMA' \
		 -profile conda_local \ 
		 -r v1.0.1 \
		 -params-file  my_parameters_1.yml

```

## Running MAGMA using conda

You can run the pipeline using Conda, Mamba or Micromamba package managers to install all the prerequisite softwares from popular repositories such as bioconda and conda-forge.

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

The `conda` environments are expected by the `conda_local` profile of the pipeline, it is recommended that it should be created **prior** to the use of the pipeline, using the following commands. Note that if you have `mamba` (or `micromamba`) available you can rely upon that instead of `conda`.


```sh
$ conda env create -n magma-env-1 --file magma-env-1.yml

$ conda env create -n magma-env-2 --file magma-env-2.yml
```

Once the environments are created, you can make use of the pipeline parameter `conda_envs_location` to inform the pipeline of the names and location of the conda envs.

> :information_source: **Conda environments and cheatsheet**: <br>
You can find out the location of conda environments using `conda env list`. [Here's](https://docs.conda.io/projects/conda/en/4.6.0/_downloads/52a95608c49671267e40c689e0bc00ca/conda-cheatsheet.pdf) a useful cheatsheet for conda operations.

## Running MAGMA using docker

We provide [two docker containers](https://github.com/orgs/TORCH-Consortium/packages?repo_name=MAGMA) with the pipeline so that you could just download and run the pipeline with them. There is **NO** need to create any docker containers, just download and enable the `docker` profile.

> :memo: **Have singularity or podman instead?**: <br>
If you do have access to Singularity or Podman, then owing to their compatibility with Docker, you can still use the provided docker containers.

Here's the command which should be used 

```console
nextflow run 'https://github.com/torch-consortium/magma' \
		 -params-file my_parameters_2.yml \
		 -profile docker \
		 -r v1.0.1 
```

> :bulb: **Hint**: <br>
You could use `-r` option of Nextflow for working with any specific version/branch of the pipeline.

## Customizing the pipeline configuration for your infrastructure

There might be cases when you need to customize the default configuration such as `cpus` and `memory` etc. For these cases, it is recommended you refer Nextflow configuration docs as well as the [default_params.config](./default_params.config) file.

Shown below is one sample configuration

- `custom.config` => Ideally this file should only contain hardware level configurations such as 

```nextflow

process {
    errorStrategy = { task.attempt < 3 ? 'retry' : 'ignore' }

    time = '1h'
    cpus = 8
    memory = 8.GB

   withName:FASTQ_VALIDATOR {
      cpus = 2
      memory = 4.GB
   }
}
```

You can then include this configuration as part of the pipeline invocation command 

```console
nextflow run 'https://github.com/torch-consortium/magma' \
		 -profile docker \
		 -r v1.0.1 \
         -c custom.config \
		 -params-file my_parameters_2.yml
```

## Running MAGMA on HPC and cloud executors

1. HPC based execution for MAGMA, please refer [this doc](./docs/hpc_execution.md).
2. Cloud batch (AWS/Google/Azure) based execution for MAGMA, please refer [this doc](./docs/cloud_batch_execution.md)

# Citation 

The MAGMA pipeline paper has been submitted.

The XBS variant calling core was published here: https://doi.org/10.1099%2Fmgen.0.000689

TODO: Update this section and add a citation.cff file 

# Contributions and Interactions with us

Contributions are warmly accepted! We encourage you to interact with us using `Discussions` and `Issues` feature of Github.

# License

Please refer the [GPL 3.0 LICENSE](./LICENSE) file.

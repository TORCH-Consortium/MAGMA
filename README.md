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

# Tutorials and Presentations

For the tutorials(./docs/tutorials.md) and [presentations](./docs/presentations.md) please refer the [docs](./docs) folder.

## Prerequisites

### Nextflow

- `git` : The version control in the pipeline.
- `Java-11` or `Java-17` (preferred)

> **Warning**
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

  Version: 22.04.5 build 5708
  Created: 15-07-2022 16:09 UTC (18:09 SAST)
  System: Mac OS X 10.15.6
  Runtime: Groovy 3.0.10 on OpenJDK 64-Bit Server VM 11.0.9.1+1-LTS
  Encoding: UTF-8 (UTF-8)

```

### Customizing pipeline parameters for your dataset

The pipeline parameters are distinct from Nextflow parameters, and therefore it is recommended that they are provided using a `yml` file as shown below


```yml

# Sample contents of my_parameters_1.yml file

input_samplesheet: /path/to/your_samplesheet.csv
only_validate_fastqs: true
conda_envs_location: /path/to/both/conda_envs
```

Which could be provided to the pipeline using `-params-file` parameter as shown below

```console
nextflow run 'https://github.com/TORCH-Consortium/MAGMA' \
		 -profile conda_local \ 
		 -params-file  my_parameters_1.yml

```

### Running MAGMA using conda

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

> **Tip**
> You can find out the location of conda environments using `conda env list`. [Here's](https://docs.conda.io/projects/conda/en/4.6.0/_downloads/52a95608c49671267e40c689e0bc00ca/conda-cheatsheet.pdf) a useful cheatsheet for conda operations.

### Running MAGMA using docker

> **NOTE**
> If you do have access to Singularity or Podman, then owing to their compatibility with Docker, you can still use the MAGMA Docker containers mentioned [docker.config](../conf/docker.config).

Here's the command which should be used 

```console
nextflow run 'https://github.com/TORCH-Consortium/MAGMA' \
		 -name experiment-1 \
		 -params-file experiment-1.yml \
		 -profile conda_local \
		 -c custom.config \
		 -r v1.1.0 
```

You could use `-r` option of Nextflow for working with any specific version/branch of the pipeline.

---------

And here are the contents of the following files

- `experiment-1.yml` => You could name it as per your convenience. Here's a sample params yaml file

```yaml
input_samplesheet: "/full/path/to/samplesheet.csv"
outdir :  "/full/path/to/magma-results"
optimize_variant_recalibration :  false
compute_minor_variants :  true
dataset_is_not_contaminated :  true
conda_envs_location :  "/home/magma-runs/magma/conda_envs"
```

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

### Running MAGMA on HPC and cloud executors

1. HPC based execution for MAGMA
2. Cloud batch (AWS/Google/Azure) based execution for MAGMA

# Citation 

The MAGMA pipeline paper has been submitted.

The XBS variant calling core was published here: https://doi.org/10.1099%2Fmgen.0.000689

TODO: Update this section and add a citation.cff file 

# Contributions

Contributions are warmly accepted!

# License

Please refer the [GPL 3.0 LICENSE](./LICENSE) file.

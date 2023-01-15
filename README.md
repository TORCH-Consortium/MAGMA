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
> Downloading the reference EXIT_RIF GVCF files from FIXME

# Usage and Tutorial

For the usage and tutorials please refer the [docs](./docs) folder.

## Prerequisites

### Nextflow

- `git` : The version control in the pipeline.
- `Java-11` or `Java-17` (preferred)

**NOTE**: The `java` version should NOT be an `internal jdk` release! You can check the release via `java -version`

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


### Running MAGMA on different environments

1. Local Conda environments for MAGMA
2. Docker based execution for MAGMA
3. HPC based execution for MAGMA
4. Cloud batch (AWS/Google/Azure) based execution for MAGMA

# Citation 

TODO: Update this section and add a citation.cff file 

# Contributions

Contributions are warmly accepted!

# License

Please refer the [GPL 3.0 LICENSE](./LICENSE) file.

# XBS-nf

XBS-nf (compleX Bacterial Samples) is a pipeline for comprehensive genomic analyses of Mycobacterium tuberculosis with a focus on clinical decision making as well as research.

# Salient features of the implementation

- Fine-grained control over resource allocation (CPU/Memory/Storage)
- Reliance of bioconda for installing packages for reproducibility
- Ease of use on a range of infrastructure (cloud/on-prem HPC clusters/ servers (or local machines))
- Resumability for failed processes
- Centralized locations for specifying analysis parameters and hardware requirements
  - XBS-nf parameters (`default_parameters.config`)
  - Hardware requirements (`conf/standard.config`)
  - Execution (software) requirements (`conf/docker.config` or `conf/conda.config`)
- A GVCF reference dataset for ~600 samples

# Usage and Tutorial

For the usage and tutorials please refer the XBS-nf website

## Prerequisites

### Nextflow

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

- Add `nextflow` to your `path` (perhaps `/usr/local/bin/`)

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

### Local Conda environments for XBS-nf

> **NOTE**: The conda environments are expected by the `conda_local` profile to be created within `xbs-nf/conda_envs` directory

- Clone the pipeline locally and `cd` into it

```sh
$ git clone https://github.com/TORCH-Consortium/xbs-nf

$ cd xbs-nf

```

- `cd` in the `conda_envs` folder and execute the following commands

```sh
$ conda env create -p xbs-nf-env-1 --file xbs-nf-env-1.yml

$ conda env create -p xbs-nf-env-2 --file xbs-nf-env-2.yml
```

> TIP: For faster installation process, please download [mamba](https://github.com/mamba-org/mamba) tool and replace `conda` with `mamba` in the above commands.

### Run the pipeline

- Customize the pipeline and process level settings in the `default_params.config` file

- From inside the `xbs-nf` folder, invoke the pipeline

```sh
$ nextflow run main.nf -profile conda
```

<!-- # Citation -->

<!-- TODO: Update this section and add a citation.cff file -->

# Contributions

Contributions are warmly accepted!

# License

Please refer the (LICENSE)[./LICENSE] file.

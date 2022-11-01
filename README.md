# magma

magma (compleX Bacterial Samples) is a pipeline for comprehensive genomic analyses of Mycobacterium tuberculosis with a focus on clinical decision making as well as research.

# Salient features of the implementation

- Fine-grained control over resource allocation (CPU/Memory/Storage)
- Reliance of bioconda for installing packages for reproducibility
- Ease of use on a range of infrastructure (cloud/on-prem HPC clusters/ servers (or local machines))
- Resumability for failed processes
- Centralized locations for specifying analysis parameters and hardware requirements
  - magma parameters (`default_parameters.config`)
  - Hardware requirements (`conf/standard.config`)
  - Execution (software) requirements (`conf/docker.config` or `conf/conda.config`)
- A GVCF reference dataset for ~600 samples

# Usage and Tutorial

For the usage and tutorials please refer the magma website

## Prerequisites

### Git tooling

- `git` and `git-lfs`

> NOTE: Without the `git-lfs` tool the optional bundled wouldn't be downloaded correctly.

### Nextflow

- `Java-11` or `Java-17` (preferred)

**NOTE**: The `java` version should NOT be an `internal jdk` release! You can check the release via `java -version`

- Download Nextflow

```bash
curl -s https://get.nextflow.io | bash
```

- Make Nextflow executable

```sh
chmod +x nextflow
```

- Add `nextflow` to your `path` (perhaps `/usr/local/bin/`)

```sh
mv nextflow /usr/local/bin

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

### Local Conda environments for magma

> **NOTE**: The conda environments are expected by the `conda_local` profile to be created within `magma/conda_envs` directory

- Clone the pipeline locally and `cd` into it

```sh
git clone https://github.com/TORCH-Consortium/magma

cd magma

```

- `cd` in the `conda_envs` folder and execute the following commands

```sh
conda env create -p magma-env-1 --file magma-env-1.yml

conda env create -p magma-env-2 --file magma-env-2.yml
```

> TIP: For faster installation process, please download [mamba](https://github.com/mamba-org/mamba) tool and replace `conda` with `mamba` in the above commands.

### Run the pipeline

- Customize the pipeline and process level settings in the [default_params](./default_params.config) file

- From inside the `magma` folder, invoke the pipeline

```sh
nextflow run main.nf -profile conda
```

- use the ```-resume``` flag to continue from previously generated output files, rather than starting from scratch.

```sh
nextflow run main.nf -profile conda -resume
```

<!-- # Citation -->

<!-- TODO: Update this section and add a citation.cff file -->

# Contributions

Contributions are warmly accepted!

# License

Please refer the [LICENSE](./LICENSE) file.

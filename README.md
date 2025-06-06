# MAGMA

MAGMA (**M**aximum **A**ccessible **G**enome for **M**tb **A**nalysis) is a pipeline for comprehensive genomic analyses of Mycobacterium tuberculosis with a focus on clinical decision making as well as research.


# DOCUMENTATION 

We are actively working on improving the documentation based on user-feedback. Please refer the following links for principal pages.

> [!IMPORTANT]  
> Usage: https://torch-consortium.github.io/MAGMA/usage.html
> 
> Customizable parameters: https://torch-consortium.github.io/MAGMA/customizable-parameters.html
> 
> Output: https://torch-consortium.github.io/MAGMA/output.html


## Go to

- [Prerequisites](#Prerequisites)
- [Customization](#Customization)
- [Analysis](#Analysis)
- [Interpretation](#Interpretation)
- [Citation](#Citation)

# Salient features of the implementation

- Fine-grained control over resource allocation (CPU/Memory/Storage)
- Reliance of bioconda for installing packages for reproducibility
- Ease of use on a range of infrastructure (cloud/on-prem HPC clusters/ servers (or local machines))
- Resumability for failed processes
- Centralized locations for specifying analysis parameters and hardware requirements
  - MAGMA parameters (`default_parameters.config` which can overridden using a params.yaml file)
  - Hardware requirements (`conf/server.config` or `conf/pbs.config` or `conf/low_memory.config`)
  - Execution (software) requirements (`conf/docker.config` or `conf/conda_local.config` or `conf/podman.config`)

## Prerequisites

### Nextflow

- `git` : The version control in the pipeline.
- `Java-11` or `Java-17` LTS release (preferred)

> :warning: **Check `java` version!**:
The `java` version should NOT be an `internal jdk` release! You can check the release via `java --version`
Notice the `LTS` next to `OpenJDK` line.


```bash

$ java -version
openjdk version "17.0.7" 2023-04-18 LTS
OpenJDK Runtime Environment (build 17.0.7+7-LTS)
OpenJDK 64-Bit Server VM (build 17.0.7+7-LTS, mixed mode, sharing)

```


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

## Samplesheet

A dummy `samplesheet` is provided [here](./samplesheet/example_MAGMA_samplesheet.csv)

The minimal samplesheet structure should have the following fields.

```csv
Sample,R1,R2
S0001,/full_path_to_directory_of_fastq_files/S0001_01_R1.fastq.gz,full_path_to_directory_of_fastq_files/S0001_01_R1.fastq.gz
S0002,/full_path_to_directory_of_fastq_files/S0002_01_R1.fastq.gz,full_path_to_directory_of_fastq_files/S0002_01_R2.fastq.gz
S0003,/full_path_to_directory_of_fastq_files/S0003_01_R1.fastq.gz,
```

If you have the metadata from sequencing instrument, you can specify further information in the samplesheet

```csv
Study,Sample,Library,Attempt,R1,R2,Flowcell,Lane,Index Sequence
Study_Name,S0001,1,1,full_path_to_directory_of_fastq_files/S0001_01_R1.fastq.gz,full_path_to_directory_of_fastq_files/S0001_01_R1.fastq.gz,1,1,1
Study_Name,S0002,1,1,full_path_to_directory_of_fastq_files/S0002_01_R1.fastq.gz,full_path_to_directory_of_fastq_files/S0002_01_R2.fastq.gz,1,1,1
Study_Name,S0003,1,1,full_path_to_directory_of_fastq_files/S0003_01_R1.fastq.gz,full_path_to_directory_of_fastq_files/S0003_01_R2.fastq.gz,1,1,1
Study_Name,S0004,1,1,full_path_to_directory_of_fastq_files/S0004_01_R1.fastq.gz,full_path_to_directory_of_fastq_files/S0004_01_R2.fastq.gz,1,1,1
```

Here's a formatted version of the CSV above, including all optional fields


|Study     |Sample|Library|Attempt|R1                                                        |R2                                                        |Flowcell|Lane|Index Sequence|
|----------|------|-------|-------|----------------------------------------------------------|----------------------------------------------------------|--------|----|--------------|
|Study_Name|S0001 |1      |1      |full_path_to_directory_of_fastq_files/S0001_01_R1.fastq.gz|full_path_to_directory_of_fastq_files/S0001_01_R1.fastq.gz|1       |1   |1             |
|Study_Name|S0002 |1      |1      |full_path_to_directory_of_fastq_files/S0002_01_R1.fastq.gz|full_path_to_directory_of_fastq_files/S0002_01_R2.fastq.gz|1       |1   |1             |
|Study_Name|S0003 |1      |1      |full_path_to_directory_of_fastq_files/S0003_01_R1.fastq.gz|full_path_to_directory_of_fastq_files/S0003_01_R2.fastq.gz|1       |1   |1             |
|Study_Name|S0004 |1      |1      |full_path_to_directory_of_fastq_files/S0004_01_R1.fastq.gz|full_path_to_directory_of_fastq_files/S0004_01_R2.fastq.gz|1       |1   |1             |



## Customization


> **Note**
We are currently working on the transition to nf-core standard (see https://github.com/TORCH-Consortium/MAGMA/issues/188), which would add standardized configurations and pipeline structure to benefit from the nf-core [nf-core/modules](https://github.com/nf-core/modules) and [nf-core/configs](https://github.com/nf-core/configs) projects.


The pipeline parameters are distinct from Nextflow parameters, and therefore it is recommended that they are provided using a `yml` file as shown below


```yml

# Sample contents of my_parameters_1.yml file

input_samplesheet: /path/to/your_samplesheet.csv
only_validate_fastqs: true
conda_envs_location: /path/to/folder/with/conda_envs
```

When running the pipeline, use profiles to ensure smooth execution on your computing system. The two types of profiles employed by the pipeline are: execution environment + memory/computing requirements

*Execution environment profiles:*

- **conda_local**
- **docker**
- **podman**

*Memory/computing profiles:*

- **pbs** (good for high performance computing clusters)
- **server** (good for local servers)
- **low_memory** (this can be run on a laptop, even limited to 8 cores and 8 GB of RAM)


> **Advanced Users**
The MAGMA pipeline has default parameters related to minimum QC thresholds that must be reached for samples to be included in the cohort analysis. These default parameters are listed in default_params.config. Users wishing to adjust these parameters should specify these adjustments in the params.yml file supplied when launching the pipeline. An example of adjusted parameters is shown below:
> **Note**
The `-profile` mechanism is used to enable infrastructure specific settings of the pipeline. The example below, assumes you are using `conda` based setup.


Which could be provided to the pipeline using `-params-file` parameter as shown below

```console
nextflow run 'https://github.com/TORCH-Consortium/MAGMA' \
         -profile conda_local, server \
         -r v1.1.1 \
         -params-file  my_parameters_1.yml

```

# Analysis

## Running MAGMA using Nextflow Tower

You can also use Seqera Platform (aka Nextflow Tower) to run the pipeline on any of the supported cloud platforms and monitoring the pipeline execution.

Please refer the [Tower docs](https://help.tower.nf/) for further information.


## Running MAGMA using conda

> :warning::warning::warning: **We discourage running MAGMA via conda, it is prone to challenging-to-reproduce errors**

You can run the pipeline using Conda, Mamba or Micromamba package managers to install all the prerequisite softwares from popular repositories such as bioconda and conda-forge.

> :information_source: **Conda environments and cheatsheet**: <br>
You can find out the location of conda environments using `conda env list`. [Here's](https://docs.conda.io/projects/conda/en/4.6.0/_downloads/52a95608c49671267e40c689e0bc00ca/conda-cheatsheet.pdf) a useful cheatsheet for conda operations.


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
$ conda env create -n magma-tbprofiler-env --file magma-tbprofiler-env.yaml
$ conda env create -n magma-ntmprofiler-env --file magma-ntmprofiler-env.yaml
```
Optionally, you can run `bash ./conda/setup_conda_envs.sh` to build all the necessary conda environments.

Once the environments are created, you can make use of the pipeline parameter `conda_envs_location` to inform the pipeline of the names and location of the conda envs.

Next, you need to load the WHO Resistance Catalog within `tb-profiler`; basically the [instructions](https://github.com/TORCH-Consortium/MAGMA/blob/master/conda_envs/setup_conda_envs.sh#L20-L23), which are used to build the necessary containers.

1. Download [magma_resistance_db_who_v1.zip](https://github.com/TORCH-Consortium/MAGMA/files/14559680/resistance_db_who_v1.zip)  and unzip it

```console
wget https://github.com/TORCH-Consortium/MAGMA/files/14559680/resistance_db_who_v1.zip

unzip resistance_db_who

```

2. Activate `magma-env-1`, which has `tb-profiler`

```console
conda activate magma-env-1

```

3. Move inside that folder and use `tb-profiler load_library` functionality to load the database


```console

cd resistance_db_who

tb-profiler load_library ./resistance_db_who

```

Success, would look like this
<img width="1060" alt="image" src="https://github.com/TORCH-Consortium/MAGMA/assets/12799326/de5eb0fc-c636-44f6-a787-39bbbf8bc8c7">



## Running MAGMA using docker

> :heavy_check_mark::heavy_check_mark::heavy_check_mark:**This is the recommended execution strategy**

We provide [two docker containers](https://github.com/orgs/TORCH-Consortium/packages?repo_name=MAGMA) with the pipeline so that you could just download and run the pipeline with them. There is **NO** need to create any docker containers, just download and enable the `docker` profile.

> 🚧 **Container build script**: The script used to build these containers is provided [here](https://github.com/TORCH-Consortium/MAGMA/tree/master/containers).

Although, you don't need to pull the containers manually, but should you need to, you could use the following commands to pull the pre-built and provided containers

```console
docker pull ghcr.io/torch-consortium/magma/magma-container-1:1.1.1

docker pull ghcr.io/torch-consortium/magma/magma-container-2:1.1.1
```


> :memo: **Have singularity or podman instead?**: <br>
If you do have access to Singularity or Podman, then owing to their compatibility with Docker, you can still use the provided docker containers.

Here's the command which should be used

```console
nextflow run 'https://github.com/torch-consortium/magma' \
         -params-file my_parameters_2.yml \
         -profile docker,pbs \
         -r v1.1.1
```

> :bulb: **Hint**: <br>
You could use `-r` option of Nextflow for working with any specific version/branch of the pipeline.

## Running MAGMA on HPC and cloud executors

1. HPC based execution for MAGMA, please refer [this doc](./docs/hpc_execution.md).
2. Cloud batch (AWS/Google/Azure) based execution for MAGMA, please refer [this doc](./docs/cloud_batch_execution.md)

## MAGMA samplesheets

In order to run the MAGMA pipeline, you must provide a samplesheet as input. The structure of the samplesheet should be that located in [samplesheet](./samplesheet/example_MAGMA_samplesheet.csv)

> :warning: **Make sure to use full paths!!!**:

- Library
```
Certain samples may have had multiple libraries prepared.
This row allows the pipeline to distinguish between
different libraries of the same sample.
```

- Attempt
```
Certain libraries may need to be sequenced multiple times.
This row allows the pipeline to distinguish between
different attempts of the same library.
```

- Flowcell/Lane/Index Sequence
```
Providing this information may allow the VQSR filtering step
to better distinguish between true variants and sequencing
errors. Including these is optional, if unknown or irrelevant,
just fill in with a '1' as shown in example_MAGMA_samplesheet.csv)
```

## (Optional) GVCF datasets

We also provide some reference GVCF files which you could use for specific use-cases.

- For **small datasets (20 samples or less)**, we recommend that you download the `EXIT_RIF GVCF` files from https://zenodo.org/record/8054182
containing GVCF reference dataset for ~600 samples is provided for augmenting smaller datasets

- For including **Mtb lineages and outgroup (M. canettii) in the phylogenetic tree**, you can download the `LineagesAndOutgroup` files from https://zenodo.org/record/8233518


```
use_ref_gvcf = false
ref_gvcf =  "/path/to/FILE.g.vcf.gz"
ref_gvcf_tbi =  "/path/to/FILE.g.vcf.gz.tbi"
```

> :bulb: **Custom GVCF dataset**: <br>
For creating a custom GVCF dataset, you can refer the discussion [here](https://github.com/TORCH-Consortium/MAGMA/issues/162).

## Tutorials and Presentations

Tim Huepink and Lennert Verboven created an in-depth tutorial of the features of the variant calling in MAGMA:

[![Video](https://img.youtube.com/vi/Kic2ItrJHj0/maxresdefault.jpg)](https://www.youtube.com/watch?v=Kic2ItrJHj0)


We have also included a presentation (in PDF format) of the logic and workflow of the MAGMA pipeline as well as posters that have been presented at conferences. Please refer the Zenodo record, https://zenodo.org/records/12633898.

# Interpretation

The results directory produced by MAGMA is as follows:

```bash
/path/to/results_dir/
.
├── QC_statistics
├── analyses
└── vcf_files
```

## QC Statistics Directory

In this directory you will find files related to the quality control carried out by the MAGMA pipeline. The structure is as follows:

```bash
/path/to/results_dir/QC_statistics
├── cohort
│   └── multiqc
│       └── multiqc_data
└── per_sample
    ├── coverage
    ├── fastqc
    └── mapping

```

- **cohort**

Here you will find the `joint.merged_cohort_stats.tsv` which contains the QC statistics for all samples in the samplesheet and allows users to determine why certain samples failed to be incorporated in the cohort analysis steps

In addition, you'll find the cohort-level MultiQC report generated by `per_sample/fastqc` analysis.

- **per_sample/coverage**

Contains the GATK WGSMetrics outputs for each of the samples in the samplesheet

- **per_sample/mapping**

> Contains the FlagStat and samtools stats for each of the samples in the samplesheet

## Analysis Directory

```bash
/path/to/results_dir/analysis
├── cluster_analysis
├── drug_resistance
├── phylogeny
└── snp_distances
```

- **Cluster Analysis**

> Contains files related to clustering based on 5SNP and 12SNP cutoffs
> **.figtree files**: These can be imported directly into Figtree for visualisation

- **Drug Resistance**

Organised based on the different types of variants as well as combined results:

```bash
/path/to/results_dir/analysis/drug_resistance
├── combined_results
├── major_variants
├── minor_variants
└── structural_variants
```

Each of the directories containing results related to the different variants (major | minor | structural) have text files that can be used to annotate the .treefiles produced by MAGMA in iToL (https://itol.embl.de)

The combined resistance results file contains a per-sample drug resistance summary based on the WHO Catalogue of *Mtb* mutations (https://www.who.int/publications/i/item/9789240082410)

MAGMA also notes the presence of all variants in in tier 1 and tier 2 drug resistance genes.

- **Phylogeny**

Contains the outputs of the IQTree phylogenetic tree construction.

> :memo: By default we recommend that you use the **ExDRIncComplex** files as MAGMA was optimized to be able to accurately call positions on the edges of complex regions in the *Mtb* genome

- **SNP distances**

Contains the SNP distance tables.

> :memo: By default we recommend that you use the **ExDRIncComplex** files as MAGMA was optimized to be able to accurately call positions on the edges of complex regions in the *Mtb* genome

## `vcf_files` Directory

```bash
/path/to/results_dir/vcf_files
├── cohort
│   ├── combined_variant_files
│   ├── minor_variants
│   ├── multiple_alignment_files
│   ├── raw_variant_files
│   ├── snp_variant_files
│   └── structural_variants
└── per_sample
    ├── minor_variants
    ├── raw_variant_files
    └── structural_variants
```

- **Combined variant files**

> Contains the cohort gvcfs based on major variants detected by the MAGMA pipeline

- **Minor variants**

> Merged vcfs of all samples, generated by LoFreq

- **Multiple alignment files**

> FASTA files for the generation of phylogenetic trees by IQTree

- **Raw variant files**

> Unfiltered indel and SNPs detected by the MAGMA pipeline

- **SNP variant files**

> Filtered SNPs detected by the MAGMA pipeline

- **Structural variant files**

> Unfiltered structural variants detected by the MAGMA pipeline

## Libraries Directory

> Contains files related to FASTQ validation and FASTQC analysis

## Samples Directory

> Contains vcf files for major|minor|structural variants for each individual samples


# Citations

The MAGMA paper has been published here: https://doi.org/10.1371/journal.pcbi.1011648

The XBS variant calling core was published here: https://doi.org/10.1099%2Fmgen.0.000689

# Contributions and interactions with us

Contributions are warmly accepted! We encourage you to interact with us using [Discussions](https://github.com/TORCH-Consortium/MAGMA/discussions) and [Issues](https://github.com/TORCH-Consortium/MAGMA/issues) feature of Github.

# License

Please refer the [GPL 3.0 LICENSE](./LICENSE) file.

[Here's](https://www.tldrlegal.com/license/gnu-general-public-license-v3-gpl-3#:~:text=You%20may%20copy%2C%20distribute%20and,along%20with%20build%20%26%20install%20instructions.) a quick TLDR for the license terms.

# XBS-nf

# Benefits of the Nextflow wrapper

- Fine-grained control over resource allocation (CPU/Memory/Storage)
- Reliance of bioconda and biocontainers for installing packages for reproducibility
- Ease of use on a range of infrastructure (cloud/on-prem clusters/local machine)
- Resumability for failed processes
- Centralized locations for specifying analysis parameters and hardware requirements
    - XBS-nf parameters (`conf/global_parameters.config`)
    - Hardware requirements (`conf/standard.config`)
    - Software requirements (`conf/docker.config` or `conf/conda.config`)

# Quickstart for a server/laptop

**NOTE**: The instructions for a cluster system like SLURM/PBS are slightly different!

The simplest use case is to analyze a few genomes on a single machine environment. Almost all aspects are customizable but for the sake of brevity, a bare bones guide for any beginner user is as shown below

- [ ] Clone the project 

```shell
git clone https://github.com/abhi18av/xbs-nf
cd xbs-nf
```

- [ ] Move your genomes (`fastq.gz files`) to a specific folder. For example `xbs-nf/data/full_data` folder

- [ ] Prepare a samplesheet using `xbs-nf/resources/reference_set/xbs-nf.test.csv` as a reference for the format.

You can optionally put your sample samplsheet in `xbs-nf/resources/reference_set/` folder.

- [ ] Update the `xbs-nf/conf/server.config` file to point to the reference sheet

- [ ] To run the pipeline, make sure you have `conda` installed. Moreover, if you don't already have `nextflow` installed, you can use the following commands to install it 

```shell
conda create -n xbs-nf-env -c bioconda -c conda-forge nextflow openjdk=11
```


You can confirm the setup by activating that environment and using the `nextflow info`  command

```
conda activate -n xbs-nf-env

nextflow info 
```

- [ ] Then simply issue the following command on the command line 

```
nextflow run main.nf -profile conda,server
```


# Contributions

Contributions are warmly accepted!


# License

TODO

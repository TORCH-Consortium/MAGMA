## Conda based execution

You can run the MAGMA pipeline using the Conda based package manager to install all the prerequisite softwares.

The `conda` environments are expected by the `conda_local` profile of the pipeline, to be created within `MAGMA/conda_envs` directory

> **NOTE**
> If you do have access to Singularity or Podman, then owing to their compatibility with Docker, you can still use the MAGMA Docker containers mentioned [docker.config](../conf/docker.config).


Here's the command which should be used 

```console
nextflow run 'https://github.com/TORCH-Consortium/MAGMA' \
		 -name experiment-1 \
		 -params-file experiment-1.yml \
		 -profile conda_local \
		 -c custom.config \
		 -r v1.0.0 
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

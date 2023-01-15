NOTE: Since we want to test more from the user-perspective, I've started to use the `-params-file` to capture pipeline parameters and the `-c` for a custom config file that anyone could rely upon.


Here's the command which should be used 

```console
nextflow run 'https://github.com/TORCH-Consortium/MAGMA' \
		 -name experiment-analysis-1 \
		 -params-file params.yml \
		 -profile conda_local \
		 -c custom.config \
		 -r v1.0.0 
```

You could use `-r` option of Nextflow for working with any specific version/branch of the pipeline.

---------

And here are the contents of the following files

- `experiment-name.yml` => You could name it as per your convenience etc

```yaml
input_samplesheet: "/full/path/to/samplesheet.csv"
outdir :  "/full/path/to/magma-results"
optimize_variant_recalibration :  false
compute_minor_variants :  true
dataset_is_not_contaminated :  true
conda_envs_location :  "/home/magma-runs/magma/conda_envs"
```

- `custom.config` => Ideally this file should only contain hardware level configurations

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


- `custom.config` => Ideally should only contain hardware level configurations

```nextflow
executor {
    //queueSize = 1
    pollInterval = '5sec'
}

process {

    executor = "slurm"
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

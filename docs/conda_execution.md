
> **NOTE**: The conda environments are expected by the `conda_local` profile of the pipeline, to be created within `MAGMA/conda_envs` directory

- Clone the pipeline locally and `cd` into it

```sh
$ git clone https://github.com/TORCH-Consortium/MAGMA

$ cd MAGMA

```

- `cd` in the `conda_envs` folder and execute the following commands

```sh
$ conda env create -p magma-env-1 --file magma-env-1.yml

$ conda env create -p magma-env-2 --file magma-env-2.yml
```

> TIP: For faster installation process, please download [mamba](https://github.com/mamba-org/mamba) tool and replace `conda` with `mamba` in the above commands.

### Run the pipeline

- Customize the pipeline and process level settings in the [default_params](./default_params.config) file

- From inside the `magma` folder, invoke the pipeline

```sh
$ nextflow run main.nf -profile conda
```
- use the ```-resume``` flag to continue from previously generated output files, rather than starting from scratch.

```sh
$ nextflow run main.nf -profile conda -resume
```



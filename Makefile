# https://makefiletutorial.com/

run_dev:
	nextflow run main.nf -profile conda,dev -entry TEST  -resume -with-tower

run_devslurm:
	nextflow run main.nf -profile conda,dev_slurm -entry TEST  -resume

run_test:
	nextflow run main.nf -params-file params/standard.yml -entry test -resume

sync:
	bash _resources/sync.sh

sync_devslurm:
	bash _resources/sync_slurmhpc.sh

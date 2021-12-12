# https://makefiletutorial.com/

run_stub:
	bash ./data/mock_data/generate_mock_files.sh && nextflow run main.nf -entry TEST -profile dev -stub-run -process.cpus 1 -process.memory 1.GB -resume

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

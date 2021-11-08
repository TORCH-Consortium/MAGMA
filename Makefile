# https://makefiletutorial.com/

run_stub:
	bash ./data/mock_data/generate_mock_files.sh && nextflow run main.nf -entry TEST -profile dev -stub-run -process.cpus 1 -process.memory 1.GB -resume

run_dev:
	nextflow run main.nf -profile standard,conda,dev -entry TEST -queue-size 1  -process.cpus 8 -process.memory 14.GB -resume -with-tower

run_test:
	nextflow run main.nf -params-file params/standard.yml -entry test -resume

sync:
	bash _resources/sync.sh

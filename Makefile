# https://makefiletutorial.com/

run_stub:
	bash ./data/mock_data/generate_mock_files.sh && nextflow run main.nf -entry TEST -profile dev -stub-run

run_dev:
	nextflow run main.nf -profile standard,conda,dev -entry TEST -queueSize 1 -resume

run_test:
	nextflow run main.nf -params-file params/standard.yml -entry test -resume

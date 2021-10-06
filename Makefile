# https://makefiletutorial.com/

run_stub:
	bash ./data/mock_data/generate_mock_files.sh && nextflow run main.nf -profile dev -stub-run

run_dev:
	nextflow run main.nf -profile dev -resume -with-tower

run_test:
	nextflow run main.nf -params-file params/standard.yml -entry test -resume
